"""Tutorial 09 -- co-translational synthesis through an **analytic exit tunnel**.

This is a self-contained elongation runner that replaces the explicit
ribosome beads of tutorials 07/08 with an **analytic tunnel** (a hole drilled
through an infinite wall): a cylindrical bore of radius ``r`` along the X-axis
through an infinite wall at ``x_exit``. See ``PLAN.md`` in this folder.

Why a separate script: the analytic ("cylinder") ribosome is a different
*physics* of confinement than the explicit-bead ribosome in
:mod:`cosmo.translation` (build steps v1/v2), so it lives here as a tutorial
variant and does **not** modify the shipped ``cosmo.translation`` package. It
**reuses** that package's tested, unchanged building blocks (the per-length
build/seed/restrain/output machinery) and only adds the one new force --
:func:`add_tunnel_cylinder` -- plus a nascent-only elongation loop around it.

The simulated System is the **nascent chain only** (no ribosome beads -> fast,
no bead clashes/jamming). The C-terminus is seeded and position-restrained on
the tunnel axis at the PTC ``(x_lo, 0, 0)``; the analytic tunnel confines the
in-tunnel segment so the chain stays extended and threads out the exit instead
of collapsing into a globule at the PTC (the 07/08 failure mode).

Run from this folder::

    python cylinder.py -f elongate.ini

then visualise with the shipped movie tool (no ``--ribosome`` -- the tunnel is
analytic, there are no beads to draw)::

    cosmo-elongate-movie -o synth_out
    vmd -e synth_out/movie.tcl
"""
from __future__ import annotations

import argparse
import configparser
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit

import cosmo
from cosmo import engine
from cosmo.core import models
from cosmo.parameters import model_parameters
from cosmo.utils.config import strtobool

# Reuse the *unchanged* cosmo.translation building blocks (PLAN.md §5: "reuse the
# v1 output path"). Importing them keeps this tutorial structurally identical to
# the shipped runner without copying code or touching the package.
from cosmo.translation.elongate import (
    TUNNEL_AXIS,
    ElongationParams,
    _make_cfg,
    _quiet,
    _write_pdb,
    add_cterm_restraint,
    cold_start_positions,
    write_subset_structure,
)

# Default tunnel wall stiffness (kJ/mol/nm^2 = 20 kcal/mol/A^2): the radial +
# exit-face + PTC-end "infinite wall" is a stiff finite harmonic (PLAN.md §3, §4).
TUNNEL_CYL_K = 8368.0


# --------------------------------------------------------------------------
# The analytic tunnel force (PLAN.md §3)
# --------------------------------------------------------------------------
def add_tunnel_cylinder(system, nascent_indices, r_nm: float,
                        x_lo_nm: float, x_exit_nm: float,
                        k: float = TUNNEL_CYL_K,
                        y0_nm: float = 0.0, z0_nm: float = 0.0,
                        mouth_round_nm: float = 0.2) -> mm.Force:
    """Add the analytic exit-tunnel restraint (a hole in an infinite wall).

    One ``CustomExternalForce`` over every nascent bead penalising the
    **penetration depth into the solid ribosome** ``S`` -- everything outside the
    bore up to the exit face, plus the closed PTC end (PLAN.md §2):

    ``S = { x < x_exit AND d > r } ∪ { x < x_lo }``,  ``d = |(y,z) - (y0,z0)|``.

    The bead escapes ``S`` via whichever face is nearer -- the bore wall
    (``d - r``, a radial inward push that keeps the in-tunnel chain extended) or
    the exit face (``x_exit - x``, a ``+x`` push so a cytosol bead can only
    re-enter the tunnel through the bore, never off-axis). The 90° inner corner
    at the mouth ``(x_exit, r)`` is rounded by a fillet of radius ``rho`` so the
    potential is continuous (-> stable MD)::

        U   = k * max(0, pen)^2 + k * min(0, x - x_lo)^2
        pen = (rounded) min( x_exit - x , d - r )      # 0 outside S; > 0 inside S

    Parameters
    ----------
    system : openmm.System
        System to add the force to (the nascent-only System).
    nascent_indices : iterable of int
        Particle indices the tunnel acts on (every nascent bead).
    r_nm : float
        Bore radius (nm).
    x_lo_nm, x_exit_nm : float
        Closed PTC end and the exit-face / infinite wall (nm); ``x_exit`` is the
        far end of the bore (= ``x_lo + tunnel_length``).
    k : float
        Wall stiffness (kJ/mol/nm^2) -- radial + exit-face + PTC-end.
    y0_nm, z0_nm : float
        Tunnel axis position in the y-z plane (nm); default the X-axis.
    mouth_round_nm : float
        Fillet radius ``rho`` rounding the mouth corner (nm).

    Returns
    -------
    openmm.Force
        The added ``CustomExternalForce`` (already in ``system``).
    """
    # Globals only -- no per-particle parameters (PLAN.md §3). `select(c, a, b)`
    # returns a when c != 0 else b; `step(x)` is 1 for x >= 0 else 0.
    # The wall stiffness global is named `ktun` (not `k`) so it never collides with
    # the C-terminus position restraint's own global `k` -- OpenMM shares global
    # parameters by name across all forces and rejects conflicting default values.
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
# Per-run parameters (cylinder additions on top of the shared run params)
# --------------------------------------------------------------------------
@dataclass
class CylinderParams(ElongationParams):
    """Elongation params + the analytic tunnel geometry (PLAN.md §4 defaults).

    Subclasses :class:`cosmo.translation.elongate.ElongationParams` so the shared
    per-length config helper (``_make_cfg``) sees the same base fields (model,
    schedule, integrator, device); only the tunnel-geometry fields are new.
    """
    tunnel_radius_nm: float = 0.9          # bore radius r (~3 CG beads wide)
    tunnel_length_nm: float = 10.0         # bore length; x_exit = x_lo + length
    tunnel_x_lo_nm: float = 0.0            # PTC / closed end
    tunnel_center_nm: Tuple[float, float] = (0.0, 0.0)   # (y0, z0): tunnel axis
    tunnel_k: float = TUNNEL_CYL_K         # wall stiffness (kJ/mol/nm^2)
    tunnel_mouth_round_nm: float = 0.2     # mouth-corner fillet radius rho
    # Post-elongation phase (runs after the chain reaches its final length, for
    # post_elongation_steps steps; 0 = skip). 'ejection' releases the C-terminus
    # restraint and lets the finished protein diffuse out the tunnel exit (-> the
    # 'ejection/' folder); 'stallation' keeps the restraint so the chain stays
    # threaded/stalled in the tunnel (-> 'stallation/'). The analytic tunnel (bore +
    # closed PTC end + exit wall) stays on in both, so the only way out is the exit.
    post_elongation: str = "stallation"
    post_elongation_steps: int = 0


# --------------------------------------------------------------------------
# Single length step (nascent-only; reuses the v1 output path)
# --------------------------------------------------------------------------
def run_length(L: int, *, full_pdb: str, prev_final: Optional[np.ndarray],
               out_root: Path, params: CylinderParams,
               cg_bond_length_nm: float, cterm_seed: np.ndarray,
               x_lo: float, x_exit: float,
               seed_override: Optional[np.ndarray] = None,
               restrain: bool = True, out_subdir: Optional[str] = None,
               n_steps_override: Optional[int] = None,
               label: Optional[str] = None) -> np.ndarray:
    """Build, seed, (restrain,) minimize and run one length-``L`` nascent System.

    The System is the nascent chain only (no ribosome beads); the analytic tunnel
    (:func:`add_tunnel_cylinder`) supplies all ribosome confinement. Mirrors the
    nascent-only branch of :func:`cosmo.translation.elongate.run_length` (the "v1
    output path"), substituting the cylinder for the planar tunnel wall.

    The same routine drives both an elongation step and the post-elongation phase
    (ejection / stallation); these arguments tailor it to the latter:

    - ``seed_override`` : use these ``(L, 3)`` nm coordinates directly (the
      fully-synthesized structure) instead of cold-start / new-residue placement.
    - ``restrain`` : if False, drop the C-terminus restraint (ejection -- the
      finished protein is released and free to diffuse out the exit).
    - ``out_subdir`` : output folder under ``out_root`` (default ``L_<L>``); e.g.
      ``ejection`` / ``stallation``.
    - ``n_steps_override`` : run this many steps instead of ``params.n_steps``.
    - ``label`` : console-banner text.

    Returns the final nascent ``(L, 3)`` nm coordinate array (seeds length L+1).
    """
    print()
    print("#" * 66)
    print("# " + (label or f"Nascent length L = {L}"))
    print("#" * 66)

    out_dir = out_root / (out_subdir or f"L_{L:03d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. length-L native structure (residue identities + connectivity).
    sub_pdb = str(out_dir / f"native_1_{L}.pdb")
    write_subset_structure(full_pdb, L, sub_pdb)

    # 2. build the length-L cosmo model (sequence-based; chatty -> silence it).
    with _quiet():
        cgModel = models.buildCoarseGrainModel(sub_pdb, model=params.model,
                                               minimize=False, check_forces=False)
    print(f"[ build ] {sub_pdb}")
    print(f"  chains={cgModel.n_chains}  CA atoms={cgModel.n_atoms}  "
          f"bonds={cgModel.n_bonds}  model={params.model}  (analytic tunnel)")

    # 3. seed coordinates. Post-elongation (seed_override): use the finished
    #    structure as-is. Cold start (L == L0): lay residues 1..L0 extended along
    #    +x from the PTC (C-terminus at cterm_seed = (x_lo, y0, z0), N-terminus
    #    toward the exit). L > L0: keep 1..L-1 from the previous final structure and
    #    seed the new C-terminal residue at its rest point cterm_seed.
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

    cfg = _make_cfg(out_dir, sub_pdb, params)
    if n_steps_override is not None:
        cfg.md_steps = n_steps_override

    # 4. hold the current C-terminus (residue L) on the tunnel axis at the PTC with
    #    a harmonic position restraint (no tRNA tether -- no bead in cylinder mode).
    #    Skipped for ejection (restrain=False -> the protein is released).
    if restrain:
        add_cterm_restraint(cgModel.system, L - 1, cterm_seed, params.restraint_k)

    # 4b. the analytic exit tunnel over every nascent bead (replaces the planar
    #     tunnel_wall; the cylinder already includes the closed-PTC-end term). Kept
    #     on during ejection too, so the released protein can only leave via the exit.
    y0, z0 = params.tunnel_center_nm
    add_tunnel_cylinder(cgModel.system, range(L), r_nm=params.tunnel_radius_nm,
                        x_lo_nm=x_lo, x_exit_nm=x_exit, k=params.tunnel_k,
                        y0_nm=y0, z0_nm=z0,
                        mouth_round_nm=params.tunnel_mouth_round_nm)

    # 5. v1 output path: seed.pdb via init_position (coords used as-is so the
    #    absolute tunnel frame is preserved), dump the nascent PSF, run, finalize.
    seed_pdb = str(out_dir / "seed.pdb")
    _write_pdb(cgModel.topology, nascent_pos, seed_pdb)
    cfg.init_position = seed_pdb
    with _quiet():
        cgModel.dumpTopology(cfg.output_path(".psf"))

    start = time.time()
    with _quiet():
        ctx = engine.setup_simulation(cfg, cgModel, control_file=None,
                                      shift_positions=False)
        if params.minimize:
            ctx.simulation.minimizeEnergy()
            ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        engine.attach_reporters(cfg, ctx.simulation, total_steps=cfg.md_steps)
        ctx.simulation.step(cfg.md_steps)
        engine.finalize_simulation(cfg, ctx, cgModel.topology, start)
    print(f"  ran {cfg.md_steps} steps ({params.device})"
          f"{' + minimize' if params.minimize else ''} (analytic tunnel"
          f"{', C-terminus RELEASED' if not restrain else ''}) "
          f"-> {cfg.output_path('_final.pdb')}  ({time.time() - start:.1f} s)")

    # 6. final nascent coordinates seed the next length.
    final = mm.app.PDBFile(cfg.output_path("_final.pdb")).getPositions(
        asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]


# --------------------------------------------------------------------------
# Elongation loop
# --------------------------------------------------------------------------
def run_elongation(full_pdb: str, *, L0: int, L_max: Optional[int] = None,
                   out_root: str = "synth_out",
                   params: Optional[CylinderParams] = None) -> None:
    """Run the nascent-chain elongation loop ``L = L0 .. L_max`` through the tunnel.

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein (the nascent chain at full length).
    L0 : int
        Starting nascent-chain length (cold-start layout).
    L_max : int, optional
        Final length; defaults to the full residue count ``N_full``.
    out_root : str
        Root output directory; each length writes to ``<out_root>/L_<L>/``.
    params : CylinderParams, optional
        Per-length run parameters + tunnel geometry (defaults to the §4 settings).
    """
    if params is None:
        params = CylinderParams()

    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # Tunnel geometry (analytic; no ribosome PDB). x_exit = the infinite exit wall.
    x_lo = params.tunnel_x_lo_nm
    x_exit = x_lo + params.tunnel_length_nm
    y0, z0 = params.tunnel_center_nm
    # C-terminus rest point: on the axis at the PTC (the cold start lays the chain
    # on-axis from here, and each new residue is seeded here).
    cterm_seed = np.array([x_lo, y0, z0])

    print(f"Analytic exit tunnel: bore radius r = {params.tunnel_radius_nm} nm, "
          f"length {params.tunnel_length_nm} nm "
          f"(x in [{x_lo:.2f}, {x_exit:.2f}]), axis (y,z) = ({y0}, {z0}) nm, "
          f"mouth fillet {params.tunnel_mouth_round_nm} nm, k = {params.tunnel_k} kJ/mol/nm^2.")
    print(f"C-terminus seeded + restrained on-axis at the PTC {np.round(cterm_seed, 3)} nm "
          f"(k = {params.restraint_k} kJ/mol/nm^2; no tRNA tether).")

    # CG protein bond length for the cold-start spacing (model-dependent).
    cg_bond_length_nm = float(
        model_parameters.parameters[params.model]["bond_length_protein"])

    # Full-length residue count (from the nascent PDB's CA atoms).
    import MDAnalysis as mda
    N_full = len(mda.Universe(full_pdb).select_atoms("protein and name CA"))
    if L_max is None:
        L_max = N_full
    if not (1 <= L0 <= L_max <= N_full):
        raise ValueError(f"require 1 <= L0 <= L_max <= N_full; got L0={L0}, "
                         f"L_max={L_max}, N_full={N_full}.")

    print()
    print(f"Elongating {full_pdb}: L = {L0} .. {L_max} (N_full = {N_full}), "
          f"model={params.model}, {params.n_steps} steps/residue, analytic tunnel.")

    prev_final: Optional[np.ndarray] = None
    for L in range(L0, L_max + 1):
        prev_final = run_length(
            L, full_pdb=full_pdb, prev_final=prev_final, out_root=out_path,
            params=params, cg_bond_length_nm=cg_bond_length_nm,
            cterm_seed=cterm_seed, x_lo=x_lo, x_exit=x_exit)

    print()
    print(f"Done. Elongated {L0} -> {L_max}. Per-length outputs under {out_path}/")

    # Post-elongation phase: once the chain reaches its final length, either release
    # the C-terminus restraint and let the finished protein diffuse out the tunnel
    # exit (ejection) or keep it threaded/stalled (stallation). Continues the same
    # length-L_max system from the final synthesized structure; the analytic tunnel
    # stays on (the only way out is the exit).
    if params.post_elongation_steps > 0:
        phase = params.post_elongation.strip().lower()
        if phase not in ("ejection", "stallation"):
            raise ValueError(f"post_elongation must be 'ejection' or 'stallation', "
                             f"got {params.post_elongation!r}.")
        restrain = phase == "stallation"
        print()
        print(f"=== Post-elongation: {phase} (L = {L_max}, "
              f"{params.post_elongation_steps} steps, C-terminus restraint "
              f"{'ON' if restrain else 'OFF -> free diffusion'}) "
              f"-> {out_path / phase}/ ===")
        run_length(
            L_max, full_pdb=full_pdb, prev_final=None, out_root=out_path,
            params=params, cg_bond_length_nm=cg_bond_length_nm,
            cterm_seed=cterm_seed, x_lo=x_lo, x_exit=x_exit,
            seed_override=prev_final, restrain=restrain, out_subdir=phase,
            n_steps_override=params.post_elongation_steps,
            label=f"Post-elongation: {phase} (L = {L_max})")
        print(f"Done. {phase.capitalize()} written to {out_path / phase}/")


# --------------------------------------------------------------------------
# INI control file
# --------------------------------------------------------------------------
@dataclass
class CylinderConfig:
    """Parsed contents of a cylinder elongation control file (``elongate.ini``)."""
    pdb_file: str
    L0: int
    L_max: Optional[int] = None
    outdir: str = "synth_out"
    params: CylinderParams = None
    config_file: Optional[str] = None


def read_elongate_config(config_file: str, verbose: bool = True) -> CylinderConfig:
    """Parse a cylinder elongation control file (INI) into a :class:`CylinderConfig`.

    Single ``[OPTIONS]`` section. Required: ``pdb_file``, ``L0``. The ``ribosome``
    PDB is **optional** here (the tunnel geometry comes from the params, not a
    structure). Recognised keys (all optional fall back to the §4 defaults):

    - ``pdb_file`` -- full native PDB of the target protein (the nascent chain).
    - ``L0`` / ``L_max`` -- start / final nascent length (blank L_max -> full).
    - ``outdir`` -- root output directory (per-length subfolders ``L_<L>/``).
    - ``model`` -- nascent force field (``hps_urry`` default, ``hps_kr``,
      ``hps_ss``, ``mpipi``).
    - ``n_steps``, ``dt``, ``ref_t``, ``tau_t``, ``nstout`` -- schedule / integrator.
    - ``device`` ('CPU'/'GPU'), ``ppn`` (CPU threads), ``minimize`` (yes/no).
    - ``restraint_k`` -- C-terminus position-restraint constant (kJ/mol/nm^2).
    - ``ribosome_model`` -- must be ``cylinder`` (the only mode this script runs).
    - ``tunnel_radius`` -- bore radius r (nm; default 0.9).
    - ``tunnel_length`` -- bore length (nm; default 10.0); ``x_exit = x_lo + length``.
    - ``tunnel_x_lo`` -- PTC / closed end (nm; default 0.0).
    - ``tunnel_center`` -- tunnel axis ``"y0,z0"`` (nm; default ``0,0``).
    - ``tunnel_k`` -- wall stiffness (kJ/mol/nm^2; default 8368 = 20 kcal/mol/A^2).
    - ``tunnel_mouth_round`` -- mouth-corner fillet radius rho (nm; default 0.2).
    - ``post_elongation`` -- post-synthesis phase at the final length: ``ejection``
      (release the C-terminus restraint -> the protein diffuses out the exit) or
      ``stallation`` (keep it threaded). Written to ``<outdir>/<phase>/``.
    - ``post_elongation_steps`` -- steps for that phase (0 = skip; default 0). Use a
      long run to watch the protein actually clear the tunnel.

    Inline ``#``/``;`` comments are ignored; underscores in ``n_steps`` are allowed.
    **Units:** OpenMM defaults (nm, ps, kJ/mol, K, kJ/mol/nm^2).
    """
    def log(msg: str) -> None:
        if verbose:
            print(msg)

    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read elongation config file: {config_file!r}")
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

    log(f"Reading cylinder elongation parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    L0 = int(req("L0"))
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None
    outdir = opt("outdir") or "synth_out"

    ribosome_model = (opt("ribosome_model") or "cylinder").lower()
    if ribosome_model != "cylinder":
        raise ValueError(
            f"{config_file}: ribosome_model = {ribosome_model!r}; this script only "
            f"runs the analytic 'cylinder' tunnel. For the explicit-bead ribosome "
            f"('beads'), use cosmo-elongate (tutorials 07/08).")

    p = CylinderParams()
    if opt("model") is not None:
        p.model = opt("model")
    if opt("n_steps") is not None:
        p.n_steps = int(str(opt("n_steps")).replace("_", ""))
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
    if opt("restraint_k") is not None:
        p.restraint_k = float(opt("restraint_k"))
    if opt("minimize") is not None:
        p.minimize = bool(strtobool(opt("minimize")))
    # --- tunnel geometry ---
    if opt("tunnel_radius") is not None:
        p.tunnel_radius_nm = float(opt("tunnel_radius"))
    if opt("tunnel_length") is not None:
        p.tunnel_length_nm = float(opt("tunnel_length"))
    if opt("tunnel_x_lo") is not None:
        p.tunnel_x_lo_nm = float(opt("tunnel_x_lo"))
    if opt("tunnel_center") is not None:
        parts = [s for s in opt("tunnel_center").replace(",", " ").split()]
        if len(parts) != 2:
            raise ValueError(f"{config_file}: tunnel_center must be 'y0,z0'; got "
                             f"{opt('tunnel_center')!r}.")
        p.tunnel_center_nm = (float(parts[0]), float(parts[1]))
    if opt("tunnel_k") is not None:
        p.tunnel_k = float(opt("tunnel_k"))
    if opt("tunnel_mouth_round") is not None:
        p.tunnel_mouth_round_nm = float(opt("tunnel_mouth_round"))
    # --- post-elongation ---
    if opt("post_elongation") is not None:
        p.post_elongation = opt("post_elongation")
    if opt("post_elongation_steps") is not None:
        p.post_elongation_steps = int(str(opt("post_elongation_steps")).replace("_", ""))

    log(f"  inputs: pdb_file={pdb_file} (ribosome: analytic tunnel, no PDB)")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}, "
        f"model={p.model}, n_steps={p.n_steps}")
    log(f"  tunnel: r={p.tunnel_radius_nm} nm, length={p.tunnel_length_nm} nm, "
        f"x_lo={p.tunnel_x_lo_nm} nm, center={p.tunnel_center_nm} nm, "
        f"k={p.tunnel_k} kJ/mol/nm^2, mouth_round={p.tunnel_mouth_round_nm} nm")
    log(f"  mechanics: restraint_k={p.restraint_k} kJ/mol/nm^2, minimize={p.minimize}")
    log(f"  post-elongation: {p.post_elongation if p.post_elongation_steps > 0 else 'off'}"
        + (f" ({p.post_elongation_steps} steps)" if p.post_elongation_steps > 0 else ""))
    log(f"  integrator: dt={p.dt_ps} ps, ref_t={p.ref_t} K, tau_t={p.tau_t} /ps, nstout={p.nstout}")
    log(f"  hardware/output: device={p.device}, ppn={p.ppn}, outdir={outdir}")

    return CylinderConfig(pdb_file=pdb_file, L0=L0, L_max=L_max, outdir=outdir,
                          params=p, config_file=config_file)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> None:
    """Console entry point: ``python cylinder.py -f elongate.ini``."""
    parser = argparse.ArgumentParser(
        prog="cylinder.py",
        description="Co-translational synthesis through an analytic exit tunnel "
                    "(tutorial 09). Grows the nascent chain N->C one residue per "
                    "step, restraining the C-terminus on the tunnel axis at the "
                    "PTC; an analytic cylindrical bore (hole in an infinite wall) "
                    "confines the in-tunnel segment. Controlled by an INI file: "
                    "python cylinder.py -f elongate.ini",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input", "-f", dest="config", type=str,
                        help="elongation control file (INI, [OPTIONS] section).")
    parser.add_argument("-o", "--outdir", default=None,
                        help="override the output directory from the config file.")
    parser.add_argument("--device", default=None, choices=["CPU", "GPU"],
                        help="override the compute device from the config file.")

    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(argv)
    if not args.config:
        parser.error("an elongation control file is required: -f elongate.ini")

    print(f"OpenMM version: {mm.__version__}")

    cfg = read_elongate_config(args.config)
    if args.outdir:
        cfg.outdir = args.outdir
    if args.device:
        cfg.params.device = args.device

    run_elongation(cfg.pdb_file, L0=cfg.L0, L_max=cfg.L_max,
                   out_root=cfg.outdir, params=cfg.params)


if __name__ == "__main__":
    main()
