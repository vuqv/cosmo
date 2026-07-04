"""Co-translational synthesis through an analytic exit tunnel (``cosmo.csp.cylinder``).

A parallel of :mod:`cosmo.csp.protocol` for the **cylinder ribosome model**. Instead
of explicit rigid CG ribosome beads, the ribosome is a pure boundary condition -- an
analytic exit tunnel (a bore of radius ``r`` in an otherwise-infinite wall,
:func:`add_tunnel_cylinder`). There are no ribosome beads, no tRNA tether and no A/P
translocation, so **each residue is a single MD segment** (not the three sub-stages of
the explicit protocol).

Timing is the **same O'Brien kinetics** as :mod:`cosmo.csp.protocol`: each residue's
MD length comes from its codon mean-first-passage time (:mod:`cosmo.csp.kinetics`),
mapped to integration steps by the ``scale_factor`` / ``dt`` compression. The
C-terminus is held on the tunnel axis at the PTC by a harmonic position restraint; the
bore keeps the in-tunnel segment extended and lets the finished chain leave only
through the exit.

Mirrors the sibling ``topo`` project's ``topo/csp/cylinder.py`` but on cosmo's IDP
force field: a length-``L`` model is :func:`cosmo.models.buildCoarseGrainModel` on the
first ``L`` residues (no STRIDE / Gō contacts). Reused from :mod:`cosmo.csp.core`: the
coordinate seeding, C-terminus restraint, config helpers and :class:`RunParams`. New
here: the analytic tunnel force and the single-segment kinetic driver.

Drive it with an INI control file (see :func:`read_cylinder_config`)::

    cosmo-cylinder -f cylinder.ini
    python -m cosmo.csp.cylinder -f cylinder.ini

**Units:** OpenMM defaults -- length nm, time ps, energy kJ/mol, temperature K,
force constants kJ/mol/nm^2. In-vivo dwell times in the kinetics are **seconds**.
"""
from __future__ import annotations

import argparse
import configparser
import random
import time
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit

import MDAnalysis as mda

from cosmo import engine
from cosmo.core import models
from cosmo.parameters import model_parameters
# Reuse the cosmo.csp.core per-length building blocks -- the cylinder driver swaps only
# the confinement (analytic tunnel) and the timing (single segment / residue).
from cosmo.csp.core import (TUNNEL_AXIS, RunParams, add_cterm_restraint,
                            cold_start_positions, write_subset_structure,
                            _make_cfg, _write_pdb, _quiet, _vprint)
from cosmo.utils.config import strtobool
from cosmo.csp import kinetics

# Default tunnel wall stiffness (kJ/mol/nm^2 = 20 kcal/mol/A^2): the radial +
# exit-face + PTC-end "infinite wall" is a stiff finite harmonic.
TUNNEL_CYL_K = 8368.0


# --------------------------------------------------------------------------
# The analytic tunnel force
# --------------------------------------------------------------------------
def add_tunnel_cylinder(system, nascent_indices, r_nm: float,
                        x_lo_nm: float, x_exit_nm: float,
                        k: float = TUNNEL_CYL_K,
                        y0_nm: float = 0.0, z0_nm: float = 0.0,
                        mouth_round_nm: float = 0.2) -> mm.Force:
    """Add the analytic exit-tunnel restraint (a hole in an infinite wall).

    One ``CustomExternalForce`` over every nascent bead penalising the **penetration
    depth into the solid ribosome** ``S`` -- everything outside the bore up to the exit
    face, plus the closed PTC end::

        S = { x < x_exit AND d > r } cup { x < x_lo },   d = |(y,z) - (y0,z0)|

    The bead escapes ``S`` via whichever face is nearer -- the bore wall (``d - r``, a
    radial inward push that keeps the in-tunnel chain extended) or the exit face
    (``x_exit - x``, a +x push so a cytosol bead can only re-enter through the bore,
    never off-axis). The 90 deg inner corner at the mouth ``(x_exit, r)`` is rounded by
    a fillet of radius ``rho`` so the potential is continuous (-> stable MD)::

        U   = k * max(0, pen)^2 + k * min(0, x - x_lo)^2
        pen = (rounded) min( x_exit - x , d - r )      # 0 outside S; > 0 inside S

    Returns the added ``CustomExternalForce`` (already in ``system``).
    """
    # The wall stiffness global is named `ktun` (not `k`) so it never collides with the
    # C-terminus position restraint's own global `k` (OpenMM shares global parameters by
    # name across all forces and rejects conflicting defaults).
    energy = (
        "ktun*max(0, pen)^2 + ktun*min(0, x - xlo)^2;"
        "pen = select(corner, pround, psharp);"
        "corner = step(rho - qx)*step(rho - qd);"
        "pround = rho - sqrt((rho - qx)^2 + (rho - qd)^2);"
        "psharp = min(qx, qd);"
        "qx = xexit - x;"
        "qd = sqrt((y - y0)^2 + (z - z0)^2) - r"
    )
    force = mm.CustomExternalForce(energy)
    force.addGlobalParameter("ktun", k)
    force.addGlobalParameter("xlo", x_lo_nm)
    force.addGlobalParameter("xexit", x_exit_nm)
    force.addGlobalParameter("r", r_nm)
    force.addGlobalParameter("y0", y0_nm)
    force.addGlobalParameter("z0", z0_nm)
    force.addGlobalParameter("rho", mouth_round_nm)
    for i in nascent_indices:
        force.addParticle(int(i), [])
    system.addForce(force)
    return force


# --------------------------------------------------------------------------
# Per-run parameters (cylinder geometry on top of the shared run params)
# --------------------------------------------------------------------------
@dataclass
class CylinderParams(RunParams):
    """:class:`~cosmo.csp.core.RunParams` (MD + O'Brien kinetics) + the tunnel geometry.

    Subclasses :class:`cosmo.csp.core.RunParams`, so every MD knob (timestep,
    temperature, restraint constant, output, ...) **and** every kinetic field
    (``scale_factor``, ``time_stage_1``/``time_stage_2``, ``uniform_codon_time``,
    ``max_steps_per_stage``, ...) is inherited unchanged; only the analytic-tunnel
    geometry and the post-elongation phase fields are new.
    """
    tunnel_radius_nm: float = 0.9          # bore radius r (~3 CG beads wide)
    tunnel_length_nm: float = 10.0         # bore length; x_exit = x_lo + length
    tunnel_x_lo_nm: float = 0.0            # PTC / closed end
    tunnel_center_nm: Tuple[float, float] = (0.0, 0.0)   # (y0, z0): tunnel axis
    tunnel_k: float = TUNNEL_CYL_K         # wall stiffness (kJ/mol/nm^2)
    tunnel_mouth_round_nm: float = 0.2     # mouth-corner fillet radius rho
    # Post-synthesis free runs (once the chain reaches its final length; 0 = skip).
    # Both drop the C-terminus restraint: 'ejection' lets the finished protein diffuse
    # out the exit, then 'dissociation' continues it as a second free run. The analytic
    # tunnel stays on throughout (the only way out is the exit).
    ejection_steps: int = 0
    dissociation_steps: int = 0


# --------------------------------------------------------------------------
# Single length step (nascent-only; one MD segment)
# --------------------------------------------------------------------------
def run_length(L: int, *, full_pdb: str,
               prev_final: Optional[np.ndarray], out_root: Path,
               params: CylinderParams, cterm_seed: np.ndarray,
               x_lo: float, x_exit: float,
               seed_override: Optional[np.ndarray] = None,
               restrain: bool = True, out_subdir: Optional[str] = None,
               n_steps_override: Optional[int] = None,
               label: Optional[str] = None) -> np.ndarray:
    """Build, seed, (restrain,) minimize and run one length-``L`` nascent System.

    The System is the nascent chain only (no ribosome beads); the analytic tunnel
    (:func:`add_tunnel_cylinder`) supplies all ribosome confinement.

    The same routine drives an elongation step and the post-synthesis free runs
    (ejection / dissociation); these arguments tailor it:

    - ``seed_override`` : use these ``(L, 3)`` nm coordinates directly (the fully
      synthesized structure) instead of cold-start / new-residue placement.
    - ``restrain`` : if False, drop the C-terminus restraint (ejection/dissociation --
      the finished protein is released and free to diffuse out the exit).
    - ``out_subdir`` : output folder under ``out_root`` (default ``L_<L>``); e.g.
      ``ejection`` / ``dissociation``.
    - ``n_steps_override`` : run this many steps instead of ``params.n_steps`` (the
      kinetic driver passes the per-residue codon-dwell step count here).
    - ``label`` : console-summary text.

    Returns the final nascent ``(L, 3)`` nm coordinate array (seeds length L+1).
    """
    _vprint()
    _vprint("#" * 66)
    _vprint("# " + (label or f"Nascent length L = {L}  (analytic tunnel)"))
    _vprint("#" * 66)

    out_dir = out_root / (out_subdir or f"L_{L:03d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    cg_bond_length_nm = float(
        model_parameters.parameters[params.model]["bond_length_protein"])

    # 1. length-L native structure (bonded terms + per-residue properties).
    sub_pdb = str(out_dir / f"native_1_{L}.pdb")
    write_subset_structure(full_pdb, L, sub_pdb)

    # 2. build the length-L cosmo model (sequence-based: first L residues).
    with _quiet():
        cgModel = models.buildCoarseGrainModel(sub_pdb, model=params.model,
                                               minimize=False, check_forces=False)

    # 3. seed coordinates. Post-elongation (seed_override): the finished structure as-is.
    #    Cold start (L == L0): lay residues 1..L0 extended along +x from the PTC
    #    (C-terminus at cterm_seed = (x_lo, y0, z0), N-terminus toward the exit). L > L0:
    #    keep 1..L-1 from the previous final structure and seed the new C-terminal residue
    #    at its rest point cterm_seed.
    if seed_override is not None:
        if seed_override.shape[0] != L:
            raise ValueError(
                f"seed_override has {seed_override.shape[0]} residues but L = {L}.")
        nascent_pos = seed_override
    elif prev_final is None:
        nascent_pos = cold_start_positions(L, cterm_seed, cg_bond_length_nm)
    else:
        if prev_final.shape[0] != L - 1:
            raise ValueError(
                f"prev_final has {prev_final.shape[0]} residues but L-1 = {L - 1}.")
        nascent_pos = np.vstack([prev_final, cterm_seed[None, :]])

    # 4. nascent-only output path: write a small seed PDB and feed it via init_position
    #    (coords used as-is so the absolute tunnel frame is preserved).
    cfg = _make_cfg(out_dir, sub_pdb, params)
    if n_steps_override is not None:
        cfg.md_steps = n_steps_override
    seed_pdb = str(out_dir / "seed.pdb")
    _write_pdb(cgModel.topology, nascent_pos, seed_pdb)
    cfg.init_position = seed_pdb
    with _quiet():
        cgModel.dumpTopology(cfg.output_path(".psf"))

    # 5. hold the current C-terminus (residue L) on the tunnel axis at the PTC with a
    #    harmonic position restraint (no tRNA tether -- no bead in cylinder mode).
    #    Skipped for ejection (restrain=False -> the protein is released).
    if restrain:
        add_cterm_restraint(cgModel.system, L - 1, cterm_seed, params.restraint_k)

    # 5b. the analytic exit tunnel over every nascent bead (includes the closed-PTC-end
    #     term). Kept on during ejection too, so the released protein can only leave via
    #     the exit.
    y0, z0 = params.tunnel_center_nm
    add_tunnel_cylinder(cgModel.system, range(L), r_nm=params.tunnel_radius_nm,
                        x_lo_nm=x_lo, x_exit_nm=x_exit, k=params.tunnel_k,
                        y0_nm=y0, z0_nm=z0, mouth_round_nm=params.tunnel_mouth_round_nm)

    # 6. set up, minimize, run, finalize (reuse cosmo.engine).
    start = time.time()
    with _quiet():
        ctx = engine.setup_simulation(cfg, cgModel, control_file=None,
                                      shift_positions=False)
        if params.minimize:
            ctx.simulation.minimizeEnergy()
            ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        engine.attach_reporters(cfg, ctx.simulation, total_steps=cfg.md_steps)
        ctx.simulation.step(cfg.md_steps)
        final_pe = ctx.simulation.context.getState(getEnergy=True).getPotentialEnergy(
            ).value_in_unit(unit.kilojoule_per_mole)
        engine.finalize_simulation(cfg, ctx, cgModel.topology, start)

    print(f"  L={L:>3d}  {(label or 'run'):<28s}  {cfg.md_steps:>6d} steps  "
          f"{time.time() - start:>6.2f} s  PE={final_pe:>+13.4e} kJ/mol")

    # 7. final nascent coordinates seed the next length.
    final = mm.app.PDBFile(cfg.output_path("_final.pdb")).getPositions(
        asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]


# --------------------------------------------------------------------------
# The continuous-synthesis loop (single stage per residue; O'Brien kinetics)
# --------------------------------------------------------------------------
def run_cylinder_synthesis(full_pdb: str, *, L0: int = 1, L_max: Optional[int] = None,
                           out_root: str = "synth_out",
                           mrna: Optional[str] = None,
                           codon_time_table_path: Optional[str] = None,
                           params: Optional[CylinderParams] = None) -> None:
    """Synthesize ``L = L0 .. L_max`` through the analytic tunnel, one segment per residue.

    Mirrors :func:`cosmo.csp.protocol.run_continuous_synthesis`, but with the cylinder
    confinement and **a single MD segment per residue** (no A/P sub-stages). Each
    residue's segment length is its codon dwell time mapped to integration steps -- the
    same O'Brien kinetics (:mod:`cosmo.csp.kinetics`) the explicit protocol uses, just
    not split three ways. Writes per-length trajectories under ``out_root/L_<L>/`` and a
    per-residue ``dwell_times.dat`` table.

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein (the nascent chain at full length).
    L0, L_max : int
        First / final nascent length (``L_max=None`` -> the full residue count).
    out_root : str
        Root output directory; each length writes to ``<out_root>/L_<L>/``.
    mrna : str, optional
        mRNA sequence file (one codon per residue) for the codon-resolved kinetics.
        Required for per-codon timing (unless ``params.uniform_codon_time`` is set).
    codon_time_table_path : str, optional
        Per-codon mean-time table; ``None`` -> the bundled E. coli 310 K table.
    params : CylinderParams, optional
        Per-length run parameters + tunnel geometry (defaults to the dataclass defaults).

    Raises
    ------
    ValueError
        If the length schedule is invalid, or if non-uniform kinetics are requested
        without ``mrna`` (propagated from
        :func:`cosmo.csp.kinetics.build_codon_time_lists`).
    """
    if params is None:
        params = CylinderParams()

    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # Analytic tunnel geometry (no ribosome PDB). x_exit = the infinite exit wall; the
    # C-terminus rests on the axis at the PTC (cold start lays the chain on-axis from
    # here, and each new residue is seeded here).
    x_lo = params.tunnel_x_lo_nm
    x_exit = x_lo + params.tunnel_length_nm
    y0, z0 = params.tunnel_center_nm
    cterm_seed = np.array([x_lo, y0, z0], dtype=float)

    print(f"Analytic exit tunnel: bore radius r = {params.tunnel_radius_nm} nm, "
          f"length {params.tunnel_length_nm} nm (x in [{x_lo:.2f}, {x_exit:.2f}]), "
          f"axis (y,z) = ({y0}, {z0}) nm, mouth fillet {params.tunnel_mouth_round_nm} nm, "
          f"k = {params.tunnel_k} kJ/mol/nm^2.")
    print(f"C-terminus seeded + restrained on-axis at the PTC {np.round(cterm_seed, 3)} nm "
          f"(k = {params.restraint_k} kJ/mol/nm^2; no tRNA tether).")

    N_full = len(mda.Universe(full_pdb).select_atoms("protein and name CA"))
    if L_max is None:
        L_max = N_full
    if not (1 <= L0 <= L_max <= N_full):
        raise ValueError(f"require 1 <= L0 <= L_max <= N_full; got L0={L0}, "
                         f"L_max={L_max}, N_full={N_full}.")

    # Kinetics: per-codon mean-first-passage times. Need intrinsic[L_max] valid ->
    # at least L_max + 1 codons. Single stage -> each residue's whole codon dwell is one
    # MD segment.
    intrinsic, real, codons = kinetics.build_codon_time_lists(
        L_max + 1, uniform_codon_time=params.uniform_codon_time,
        mrna_path=mrna, codon_time_table_path=codon_time_table_path,
        ribosome_traffic=params.ribosome_traffic, initiation_rate=params.initiation_rate)
    rng = random.Random(params.random_seed)

    print()
    print("=" * 66)
    print("[ cylinder continuous synthesis -- kinetic schedule (single stage/residue) ]")
    print("=" * 66)
    print(f"  timing mode: {'uniform' if params.uniform_codon_time is not None else 'per-codon (mRNA)'}; "
          f"scale_factor={params.scale_factor:g}; dt={params.dt_ps} ps")
    if params.max_steps_per_stage is not None:
        print(f"  TEST CLAMP: <= {params.max_steps_per_stage} steps/residue. "
              f"Remove for production.")
    print(f"Synthesizing {full_pdb}: L = {L0} .. {L_max} (N_full = {N_full}), analytic tunnel.")

    dwell_log = out_path / "dwell_times.dat"
    dwell_fh = open(dwell_log, "w")
    dwell_fh.write(
        "# cylinder continuous-synthesis per-residue dwell times (cosmo.csp.cylinder)\n"
        f"#   scale_factor={params.scale_factor:g}  dt={params.dt_ps} ps  "
        f"timing={'uniform' if params.uniform_codon_time is not None else 'per-codon'}  "
        f"random_seed={params.random_seed}\n"
        "#   t_dwell = sampled codon dwell (s); ns = in-silico ns; steps = integration "
        "steps actually run (single MD segment)\n"
        "# L  codon  t_dwell_s  ns  steps\n")
    dwell_fh.flush()

    prev_final: Optional[np.ndarray] = None
    for L in range(L0, L_max + 1):
        # Single stage: sample the codon's total dwell, map to integration steps, clamp.
        dwell_s = kinetics.sample_fpt(intrinsic[L], rng)
        n_steps = kinetics.seconds_to_steps(dwell_s, params.scale_factor, params.dt_ps)
        if params.max_steps_per_stage is not None:
            n_steps = min(n_steps, int(params.max_steps_per_stage))
        n_steps = max(n_steps, int(params.min_steps_per_stage))

        codon = codons[L - 1] if codons is not None else "uniform"
        print(f"L={L:>3d}  {codon:>5s}  dwell {dwell_s:>9.4g} s  steps {n_steps:>7d}")

        ns = dwell_s * 1e9 / params.scale_factor
        dwell_fh.write(f"{L:4d}  {codon:>5s}  {dwell_s:.6e}  {ns:.6e}  {n_steps:8d}\n")
        dwell_fh.flush()

        prev_final = run_length(
            L, full_pdb=full_pdb, prev_final=prev_final, out_root=out_path,
            params=params, cterm_seed=cterm_seed, x_lo=x_lo, x_exit=x_exit,
            n_steps_override=n_steps,
            label=f"L={L}  {codon}  ({n_steps} steps, analytic tunnel)")

    dwell_fh.close()
    print()
    print(f"Done. Synthesized {L0} -> {L_max}. Per-length outputs under {out_path}/")
    print(f"Per-residue dwell-time table: {dwell_log}")

    # Post-synthesis free runs: once the chain reaches its final length, release the
    # C-terminus restraint and let the finished protein diffuse out the exit (ejection),
    # then continue as a second free run so it drifts fully clear (dissociation). Both
    # continue the length-L_max system from the previous final structure; the analytic
    # tunnel stays on (only way out is the exit).
    if params.ejection_steps > 0:
        print()
        print(f"=== Ejection (L = {L_max}, {params.ejection_steps} steps, "
              f"C-terminus restraint OFF -> free diffusion) -> {out_path / 'ejection'}/ ===")
        prev_final = run_length(
            L_max, full_pdb=full_pdb, prev_final=None, out_root=out_path, params=params,
            cterm_seed=cterm_seed, x_lo=x_lo, x_exit=x_exit, seed_override=prev_final,
            restrain=False, out_subdir="ejection",
            n_steps_override=params.ejection_steps,
            label=f"Ejection (L = {L_max})")
        print(f"Done. Ejection written to {out_path / 'ejection'}/")

    if params.dissociation_steps > 0:
        print()
        print(f"=== Dissociation (L = {L_max}, {params.dissociation_steps} steps, "
              f"C-terminus restraint OFF) -> {out_path / 'dissociation'}/ ===")
        run_length(
            L_max, full_pdb=full_pdb, prev_final=None, out_root=out_path, params=params,
            cterm_seed=cterm_seed, x_lo=x_lo, x_exit=x_exit, seed_override=prev_final,
            restrain=False, out_subdir="dissociation",
            n_steps_override=params.dissociation_steps,
            label=f"Dissociation (L = {L_max})")
        print(f"Done. Dissociation written to {out_path / 'dissociation'}/")


# --------------------------------------------------------------------------
# INI control file
# --------------------------------------------------------------------------
@dataclass
class CylinderConfig:
    """Parsed contents of a cylinder synthesis control file (``cylinder.ini``)."""
    pdb_file: str
    L0: int
    L_max: Optional[int] = None
    outdir: str = "synth_out"
    mrna: Optional[str] = None
    codon_time_table_path: Optional[str] = None
    params: Optional[CylinderParams] = None
    config_file: Optional[str] = None


def read_cylinder_config(config_file: str, verbose: bool = True) -> CylinderConfig:
    """Parse a cylinder synthesis control file (INI) into a :class:`CylinderConfig`.

    Single ``[OPTIONS]`` section. Required: ``pdb_file``, ``L0``, and -- unless
    ``codon_times`` is a positive number (uniform timing) -- ``mrna``. No ribosome PDB
    (the tunnel geometry comes from the params). Recognised keys (optional ones fall
    back to defaults):

    - ``pdb_file`` -- full native PDB of the target protein (the nascent chain).
    - ``L0`` / ``L_max`` -- start / final nascent length (blank ``L_max`` -> full).
    - ``outdir`` -- root output directory (per-length subfolders ``L_<L>/``).
    - ``model`` -- nascent-chain force field (default ``hps_kr``).
    - **Kinetics** (same as CSP): ``mrna``, ``codon_times`` (a codon table path, or a
      positive number of seconds for uniform timing; blank -> bundled E. coli 310 K),
      ``scale_factor``, ``time_stage_1``, ``time_stage_2``, ``random_seed``,
      ``max_steps_per_stage``, ``min_steps_per_stage``.
    - **Integrator / MD**: ``dt``, ``ref_t``, ``tau_t``, ``nstout``, ``device``,
      ``ppn``, ``minimize`` (yes/no), ``constraints`` ('None' | 'AllBonds'),
      ``restraint_k`` (C-terminus position-restraint constant, kJ/mol/nm^2).
    - **Tunnel geometry**: ``tunnel_radius`` (nm), ``tunnel_length`` (nm),
      ``tunnel_x_lo`` (nm), ``tunnel_center`` (``"y0,z0"`` nm), ``tunnel_k``
      (kJ/mol/nm^2), ``tunnel_mouth_round`` (nm).
    - **Post-synthesis**: ``ejection_steps`` / ``dissociation_steps`` -- free runs
      with the C-terminus restraint released (0 = skip).

    Inline ``#``/``;`` comments are ignored. **Units:** OpenMM defaults.
    """
    def log(msg: str) -> None:
        if verbose:
            print(msg)

    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read cylinder config file: {config_file!r}")
    if "OPTIONS" not in cp:
        raise ValueError(f"{config_file}: missing required [OPTIONS] section.")
    o = cp["OPTIONS"]

    def opt(key: str) -> Optional[str]:
        v = o.get(key, None)
        if v is None:
            return None
        v = v.strip()
        return v if v != "" else None

    def req(key: str) -> str:
        v = opt(key)
        if v is None:
            raise ValueError(f"{config_file}: required option '{key}' is missing or blank.")
        return v

    def as_int(s: str) -> int:
        return int(str(s).replace("_", ""))

    log(f"Reading cylinder synthesis parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    L0 = int(req("L0"))
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None
    outdir = opt("outdir") or "synth_out"
    mrna = opt("mrna")
    _uniform_codon_time, codon_time_table_path = kinetics.parse_codon_times(opt("codon_times"))

    p = CylinderParams()
    # --- integrator / MD ---
    if opt("model") is not None:
        p.model = opt("model")
    if opt("dt") is not None:
        p.dt_ps = float(opt("dt"))
    if opt("ref_t") is not None:
        p.ref_t = float(opt("ref_t"))
    if opt("tau_t") is not None:
        p.tau_t = float(opt("tau_t"))
    if opt("nstout") is not None:
        p.nstout = int(opt("nstout"))
    if opt("device") is not None:
        p.device = opt("device")
    if opt("ppn") is not None:
        p.ppn = int(opt("ppn"))
    cons = opt("constraints")
    if cons is not None:
        p.constraints = None if cons.strip().lower() == "none" else cons
    if opt("restraint_k") is not None:
        p.restraint_k = float(opt("restraint_k"))
    if opt("minimize") is not None:
        p.minimize = bool(strtobool(opt("minimize")))
    # --- O'Brien kinetics ---
    if opt("scale_factor") is not None:
        p.scale_factor = float(opt("scale_factor"))
    if opt("time_stage_1") is not None:
        p.time_stage_1 = float(opt("time_stage_1"))
    if opt("time_stage_2") is not None:
        p.time_stage_2 = float(opt("time_stage_2"))
    p.uniform_codon_time = _uniform_codon_time
    if opt("ribosome_traffic") is not None:
        p.ribosome_traffic = bool(strtobool(opt("ribosome_traffic")))
    if opt("initiation_rate") is not None:
        p.initiation_rate = float(opt("initiation_rate"))
    if opt("random_seed") is not None:
        p.random_seed = as_int(opt("random_seed"))
    if opt("max_steps_per_stage") is not None:
        p.max_steps_per_stage = as_int(opt("max_steps_per_stage"))
    if opt("min_steps_per_stage") is not None:
        p.min_steps_per_stage = as_int(opt("min_steps_per_stage"))
    # --- tunnel geometry ---
    if opt("tunnel_radius") is not None:
        p.tunnel_radius_nm = float(opt("tunnel_radius"))
    if opt("tunnel_length") is not None:
        p.tunnel_length_nm = float(opt("tunnel_length"))
    if opt("tunnel_x_lo") is not None:
        p.tunnel_x_lo_nm = float(opt("tunnel_x_lo"))
    if opt("tunnel_center") is not None:
        parts = opt("tunnel_center").replace(",", " ").split()
        if len(parts) != 2:
            raise ValueError(f"{config_file}: tunnel_center must be 'y0,z0'; got "
                             f"{opt('tunnel_center')!r}.")
        p.tunnel_center_nm = (float(parts[0]), float(parts[1]))
    if opt("tunnel_k") is not None:
        p.tunnel_k = float(opt("tunnel_k"))
    if opt("tunnel_mouth_round") is not None:
        p.tunnel_mouth_round_nm = float(opt("tunnel_mouth_round"))
    # --- post-elongation ---
    if opt("ejection_steps") is not None:
        p.ejection_steps = as_int(opt("ejection_steps"))
    if opt("dissociation_steps") is not None:
        p.dissociation_steps = as_int(opt("dissociation_steps"))

    if p.uniform_codon_time is None and mrna is None:
        raise ValueError(
            f"{config_file}: per-codon kinetics need an 'mrna' file (or set "
            f"'codon_times' to a positive number of seconds for a uniform codon time). "
            f"A 'codon_times' table path is optional (defaults to the bundled "
            f"E. coli 310 K table).")

    log(f"  inputs: pdb_file={pdb_file} (ribosome: analytic tunnel, no PDB)")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}, "
        f"model={p.model}, constraints={p.constraints}")
    _timing = (f"uniform (codon_time={p.uniform_codon_time:g} s)" if p.uniform_codon_time is not None
               else f"per-codon (mrna={mrna}, codon_times={codon_time_table_path or 'bundled E. coli 310 K'})")
    log(f"  timing: {_timing}; scale_factor={p.scale_factor:g}, "
        f"time_stage_1={p.time_stage_1:g} s, time_stage_2={p.time_stage_2:g} s")
    log(f"  tunnel: r={p.tunnel_radius_nm} nm, length={p.tunnel_length_nm} nm, "
        f"x_lo={p.tunnel_x_lo_nm} nm, center={p.tunnel_center_nm} nm, "
        f"k={p.tunnel_k} kJ/mol/nm^2, mouth_round={p.tunnel_mouth_round_nm} nm")
    log(f"  mechanics: restraint_k={p.restraint_k} kJ/mol/nm^2, minimize={p.minimize}")
    if p.ejection_steps or p.dissociation_steps:
        log(f"  post-synthesis: ejection={p.ejection_steps} steps, "
            f"dissociation={p.dissociation_steps} steps")
    else:
        log("  post-synthesis: off")
    log(f"  integrator: dt={p.dt_ps} ps, ref_t={p.ref_t} K, tau_t={p.tau_t} /ps, nstout={p.nstout}")
    log(f"  hardware/output: device={p.device}, ppn={p.ppn}, outdir={outdir}")

    return CylinderConfig(pdb_file=pdb_file, L0=L0, L_max=L_max, outdir=outdir,
                          mrna=mrna, codon_time_table_path=codon_time_table_path,
                          params=p, config_file=config_file)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def cylinder(argv: Optional[List[str]] = None) -> None:
    """Console entry point: ``cosmo-cylinder -f cylinder.ini``."""
    import warnings
    warnings.filterwarnings("ignore", category=Warning, module=r"MDAnalysis")

    parser = argparse.ArgumentParser(
        prog="cosmo-cylinder",
        description="Co-translational synthesis through an analytic exit tunnel (the "
                    "cylinder ribosome model) on cosmo's IDP force field. Grows the "
                    "nascent chain N->C, restraining the C-terminus on the tunnel axis "
                    "at the PTC; an analytic cylindrical bore (hole in an infinite wall) "
                    "confines the in-tunnel segment. One MD segment per residue, timed "
                    "by the O'Brien codon kinetics. Controlled by an INI file: "
                    "cosmo-cylinder -f cylinder.ini",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input", "-f", dest="config", type=str,
                        help="synthesis control file (INI, [OPTIONS] section).")
    parser.add_argument("-o", "--outdir", default=None,
                        help="override the output directory from the config file.")
    parser.add_argument("--device", default=None, choices=["CPU", "GPU"],
                        help="override the compute device from the config file.")

    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(argv)
    if not args.config:
        parser.error("a synthesis control file is required: -f cylinder.ini")

    print(f"OpenMM version: {mm.__version__}")

    cfg = read_cylinder_config(args.config)
    if args.outdir:
        cfg.outdir = args.outdir
    if args.device:
        cfg.params.device = args.device

    run_cylinder_synthesis(cfg.pdb_file, L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
                           mrna=cfg.mrna, codon_time_table_path=cfg.codon_time_table_path,
                           params=cfg.params)


if __name__ == "__main__":
    cylinder()
