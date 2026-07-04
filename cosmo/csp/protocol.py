"""Continuous Synthesis Protocol (O'Brien), ported to cosmo (``cosmo.csp``).

The **full O'Brien continuous-synthesis protocol** -- the per-codon, 3-stage
elongation cycle of ``continuous_synthesis_v6.py`` -- expressed in cosmo style, on
cosmo's sequence-based IDP force field (HPS / mpipi). It mirrors the sibling ``topo``
project's ``topo/csp/protocol.py`` but drops topo's structure-based Gō machinery
(STRIDE, native-contact precompute, ``domain.yaml`` nscales): a length-``L`` model is
just :func:`cosmo.models.buildCoarseGrainModel` on the first ``L`` residues (see
``cosmo/translation/PLAN.md`` §2).

What is reused vs. new:

- **Reused** from :mod:`cosmo.csp.core`: the per-length MD machinery
  (:func:`run_length` -- build, coordinate seeding, rigid-ribosome scenery, tunnel
  wall, minimize/run/finalize), :func:`read_anchor` and :class:`RunParams`.
- **Reused** from :mod:`cosmo.csp.ribosome`: :func:`load_ribosome` (the rigid scenery).
- **New** (this module): the O'Brien *kinetics* (:mod:`cosmo.csp.kinetics`) and the
  outer loop that calls :func:`run_length` **three times per residue**, switching the
  C-terminus restraint target A->P to reproduce translocation.

The 3-stage mapping onto :func:`run_length` (``L`` = nascent length):

==========  =====================================  ==========================
stage       biology                                restraint target / seed
==========  =====================================  ==========================
1           peptidyl transfer / A-site delivery    A-anchor; new residue placed
2           translocation begins                   A-anchor; continue stage 1
3           translocation to P-site / wait         P-anchor; continue stage 2
==========  =====================================  ==========================

Stage 3's final structure seeds the next residue's stage 1. The cold-start segment
(``L == L0``) is laid down the tunnel from the P-anchor (no A-site delivery yet).
Because CSP needs the restraint target to switch A<->P, it drives the **position
restraint** path (``trna_tether`` is forced off): the supplied ribosome is always
rigid scenery and the tunnel wall + excluded volume + electrostatics are on.

Drive it with an INI control file (see :func:`read_csp_config`)::

    cosmo-csp -f csp.ini
    python -m cosmo.csp -f csp.ini

**Units:** OpenMM defaults -- length nm, time ps, energy kJ/mol, temperature K,
force constants kJ/mol/nm^2. In-vivo dwell times in the kinetics are **seconds**.
"""
from __future__ import annotations

import argparse
import configparser
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np
import openmm as mm

import MDAnalysis as mda

from cosmo.csp.core import (RunParams, TUNNEL_AXIS, read_anchor, run_length,
                            optimal_ptc_targets)
from cosmo.csp.ribosome import load_ribosome, TRNA_TETHER_BOND_NM
from cosmo.parameters import model_parameters
from cosmo.utils.config import strtobool
from cosmo.csp import kinetics


# --------------------------------------------------------------------------
# The continuous-synthesis loop
# --------------------------------------------------------------------------
def run_continuous_synthesis(full_pdb: str, ribosome_pdb: str, *,
                             L0: int = 1, L_max: Optional[int] = None,
                             out_root: str = "synth_out",
                             mrna: Optional[str] = None,
                             codon_time_table_path: Optional[str] = None,
                             params: Optional[RunParams] = None) -> None:
    """Run the full O'Brien continuous synthesis ``L = L0 .. L_max``.

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein (the nascent chain at full length).
    ribosome_pdb : str
        Truncated CG ribosome PDB -- source of the P-/A-anchors and the rigid
        (mass-0) scenery (always loaded; providing it is the signal to use it).
    L0 : int, optional
        First nascent length to synthesize (default ``1``).
    L_max : int or None, optional
        Final nascent length (default ``None`` -> the full residue count).
    out_root : str
        Root output directory; each residue writes ``L_<L>/stage_<1,2,3>/``.
    mrna : str, optional
        mRNA sequence file (one codon per residue) for the codon-resolved kinetics.
        Required for per-codon timing; not needed for uniform timing
        (``params.uniform_codon_time`` set).
    codon_time_table_path : str, optional
        Per-codon mean-time table. ``None`` -> the bundled E. coli 310 K table
        (organism-universal; see
        :func:`cosmo.csp.kinetics.default_codon_time_table_path`).
    params : RunParams, optional
        Kinetic + MD/ribosome run parameters (defaults to the dataclass defaults).

    Raises
    ------
    ValueError
        If the length schedule is invalid (``1 <= L0 <= L_max <= N_full`` fails), or
        if non-uniform kinetics are requested without ``mrna`` (propagated from
        :func:`cosmo.csp.kinetics.build_codon_time_lists`).
    """
    if params is None:
        params = RunParams()
    ep = params
    # CSP switches the C-terminus restraint target A->P across the three stages; that
    # is the position-restraint path. The O'Brien tRNA tether (in cosmo) always targets
    # the P-site bead and does not switch A/P, so it is incompatible with the 3-stage
    # translocation -- force it off here (mirrors cosmo/translation/PLAN.md's CSP note).
    if ep.trna_tether:
        print("[csp] trna_tether is forced off: the 3-stage protocol uses the "
              "position-restraint path so the C-terminus target can switch A->P.")
        ep.trna_tether = False
    print("C-terminus restraint: position restraint to the A/P target point "
          "(a->a->p migration).")

    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # --- rigid ribosome (always: the supplied file is rigid scenery) --------
    # Loaded once; identical at every length. A cosmo CG ribosome PDB (O'Brien 3/4-bead
    # P/R/BR rep, built by cosmo.csp.cg_ribosome + truncate_ribosome); per-bead Rmin/2 +
    # charge from model_parameters. Loaded before the anchors, which are its beads.
    ribo = load_ribosome(ribosome_pdb, model=ep.model)
    print(f"Rigid ribosome: {ribo.n} beads from {ribosome_pdb} "
          f"(mass-0 scenery; ribosome<->nascent 12-10-6 excluded volume + Yukawa on).")

    # --- anchors (fixed points from the truncated ribosome) -----------------
    p_anchor = read_anchor(ribosome_pdb, "PtR", resid=76, bead="R")
    a_anchor = read_anchor(ribosome_pdb, "AtR", resid=76, bead="R")
    print(f"P-anchor (PtR 76 R): {p_anchor} nm")
    print(f"A-anchor (AtR 76 R): {a_anchor} nm")

    # Hold/seed targets (fixed points the C-terminus is position-restrained to).
    if ep.optimize_ptc_geometry:
        # Optimize the PTC geometry: place the A-site and P-site target points exactly
        # one peptide bond apart (the model's bond_length_protein -- 0.380 nm for hps_kr,
        # not topo's 0.381) and clear of the ribosome excluded volume. Seeding the new
        # residue at a_target while the previous C-terminus rests at p_target makes the
        # always-present peptide bond start at its equilibrium length.
        pep_nm = float(model_parameters.parameters[ep.model]["bond_length_protein"])
        a_target, p_target = optimal_ptc_targets(ribo, peptide_nm=pep_nm)
        print(f"[optimize_ptc_geometry] optimal PTC targets "
              f"(|A-P| = {np.linalg.norm(a_target - p_target):.4f} nm = 1 peptide bond "
              f"@ {pep_nm} nm; fixed points, not bonds):")
        print(f"  A-site target (new-AA seed + stage-1/2 restraint): {a_target} nm")
        print(f"  P-site target (prev-AA / stage-3 restraint)       : {p_target} nm")
    else:
        # Default: offset into the tunnel (+x) from each raw tRNA anchor bead so the
        # C-terminus does not sit on top of a ribosome bead (a near-coincident clash).
        # Offset defaults to the tether bond length (override via ptc_offset).
        offset = ep.ptc_offset_nm
        if offset is None:
            offset = TRNA_TETHER_BOND_NM
        p_target = p_anchor + offset * TUNNEL_AXIS
        a_target = a_anchor + offset * TUNNEL_AXIS

    # Tunnel-wall plane (auto-derived -- never a stale user knob): the lower
    # (deeper-in-tunnel, smaller-x) C-terminus hold plane, so the held C-terminus sits
    # at (just on) the plane while the chain cannot slip below the synthesis point.
    wall_x0 = float(min(a_target[0], p_target[0]))
    if ep.tunnel_wall:
        print(f"Tunnel wall plane: x >= {wall_x0:.4f} nm (auto: lower C-terminus hold plane).")

    # --- schedule bounds ----------------------------------------------------
    N_full = len(mda.Universe(full_pdb).select_atoms("protein and name CA"))
    if L_max is None:
        L_max = N_full
    if not (1 <= L0 <= L_max <= N_full):
        raise ValueError(f"require 1 <= L0 <= L_max <= N_full; got L0={L0}, "
                         f"L_max={L_max}, N_full={N_full}.")

    # --- kinetics: intrinsic / real per-codon mFPT lists --------------------
    # Need intrinsic[L_max] valid -> at least L_max + 1 codons.
    intrinsic, real, codons = kinetics.build_codon_time_lists(
        L_max + 1, uniform_codon_time=ep.uniform_codon_time,
        mrna_path=mrna, codon_time_table_path=codon_time_table_path,
        ribosome_traffic=ep.ribosome_traffic, initiation_rate=ep.initiation_rate)
    rng = random.Random(ep.random_seed)

    print()
    print("=" * 66)
    print("[ O'Brien continuous synthesis -- kinetic schedule ]")
    print("=" * 66)
    print(f"  timing mode: {'uniform' if ep.uniform_codon_time is not None else 'per-codon (mRNA)'}; "
          f"scale_factor={ep.scale_factor:g}; dt={ep.dt_ps} ps")
    print(f"  stage means (s): peptidyl-transfer={ep.time_stage_1:g}, "
          f"translocation={ep.time_stage_2:g}, tRNA-binding=remainder")
    if ep.max_steps_per_stage is not None:
        print(f"  TEST CLAMP: <= {ep.max_steps_per_stage} steps/stage "
              f"(~{3 * ep.max_steps_per_stage} steps/residue). Remove for production.")
    print(f"Synthesizing {full_pdb}: L = {L0} .. {L_max} (N_full = {N_full}), "
          f"model={ep.model}.")

    # --- per-residue dwell-time log (one machine-parsable row per residue) --
    dwell_log = out_path / "dwell_times.dat"
    dwell_fh = open(dwell_log, "w")
    dwell_fh.write(
        "# O'Brien continuous-synthesis per-residue dwell times (cosmo.csp)\n"
        f"#   scale_factor={ep.scale_factor:g}  dt={ep.dt_ps} ps  "
        f"time_stage_1={ep.time_stage_1:g} s  time_stage_2={ep.time_stage_2:g} s\n"
        f"#   timing={'uniform' if ep.uniform_codon_time is not None else 'per-codon'}  "
        f"random_seed={ep.random_seed}\n"
        "#   t1/t2/t3 = sampled peptidyl-transfer / translocation / tRNA-binding "
        "dwell (s); steps = clamped integration steps actually run\n"
        "# L  codon  t_invivo_total_s  t1_s  t2_s  t3_s  "
        "ns1  ns2  ns3  steps1  steps2  steps3\n")
    dwell_fh.flush()

    # --- main loop: one residue = three sub-stages --------------------------
    prev_final: Optional[np.ndarray] = None
    for L in range(L0, L_max + 1):
        (s1, s2, s3), (t1, t2, t3) = kinetics.stage_steps(
            L, intrinsic, real, time_stage_1=ep.time_stage_1,
            time_stage_2=ep.time_stage_2, scale_factor=ep.scale_factor,
            dt_ps=ep.dt_ps, rng=rng, max_steps_per_stage=ep.max_steps_per_stage,
            min_steps_per_stage=ep.min_steps_per_stage)

        codon = codons[L - 1] if codons is not None else "uniform"
        print(f"L={L:>3d}  {codon:>5s}  dwell {intrinsic[L]:>9.4g} s  "
              f"steps {s1:>4d}/{s2:>4d}/{s3:>4d}")

        dwell_fh.write(
            f"{L:4d}  {codon:>5s}  {intrinsic[L]:.6e}  "
            f"{t1:.6e}  {t2:.6e}  {t3:.6e}  "
            f"{t1 * 1e9 / ep.scale_factor:.6e}  {t2 * 1e9 / ep.scale_factor:.6e}  "
            f"{t3 * 1e9 / ep.scale_factor:.6e}  {s1:8d}  {s2:8d}  {s3:8d}\n")
        dwell_fh.flush()

        ldir = f"L_{L:03d}"
        cold = prev_final is None
        # Stage 1: deliver the new residue at the A-site and restrain there (cold start
        # lays the initial segment from the P-anchor instead).
        stage1_anchor = p_target if cold else a_target
        f1 = run_length(L, full_pdb=full_pdb, p_anchor=stage1_anchor, a_anchor=a_anchor,
                        prev_final=prev_final, out_root=out_path, params=ep, ribo=ribo,
                        wall_x0_nm=wall_x0,
                        cterm_seed=stage1_anchor, restrain=True,
                        out_subdir=f"{ldir}/stage_1", n_steps_override=s1,
                        label="stage 1 peptidyl-transfer")

        # Stage 2: continue from stage 1, still held at the A-site. Skip the redundant
        # minimization (the seeded structure is stage 1's already-relaxed final).
        f2 = run_length(L, full_pdb=full_pdb, p_anchor=stage1_anchor, a_anchor=a_anchor,
                        prev_final=None, seed_override=f1, out_root=out_path, params=ep,
                        ribo=ribo, wall_x0_nm=wall_x0,
                        cterm_seed=stage1_anchor, restrain=True,
                        out_subdir=f"{ldir}/stage_2", n_steps_override=s2,
                        minimize_override=False, label="stage 2 translocation")

        # Stage 3: translocate A->P (restrain the C-terminus to the P-target).
        f3 = run_length(L, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
                        prev_final=None, seed_override=f2, out_root=out_path, params=ep,
                        ribo=ribo, wall_x0_nm=wall_x0,
                        cterm_seed=p_target, restrain=True,
                        out_subdir=f"{ldir}/stage_3", n_steps_override=s3,
                        label="stage 3 tRNA-binding")
        prev_final = f3

    dwell_fh.close()
    print()
    print(f"Done. Synthesized {L0} -> {L_max}. Per-residue/-stage outputs under {out_path}/")
    print(f"Per-residue dwell-time table: {dwell_log}")

    # --- post-synthesis: ejection then dissociation (both free runs) --------
    if ep.ejection_steps > 0:
        print()
        print(f"=== Ejection (L = {L_max}, {ep.ejection_steps} steps, restraint OFF) "
              f"-> {out_path / 'ejection'}/ ===")
        prev_final = run_length(
            L_max, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
            prev_final=None, seed_override=prev_final, out_root=out_path, params=ep,
            ribo=ribo, wall_x0_nm=wall_x0, restrain=False,
            out_subdir="ejection", n_steps_override=ep.ejection_steps, label="ejection")

    if ep.dissociation_steps > 0:
        print()
        print(f"=== Dissociation (L = {L_max}, {ep.dissociation_steps} steps, restraint OFF) "
              f"-> {out_path / 'dissociation'}/ ===")
        run_length(
            L_max, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
            prev_final=None, seed_override=prev_final, out_root=out_path, params=ep,
            ribo=ribo, wall_x0_nm=wall_x0, restrain=False,
            out_subdir="dissociation", n_steps_override=ep.dissociation_steps,
            label="dissociation")


# --------------------------------------------------------------------------
# INI control file
# --------------------------------------------------------------------------
@dataclass
class CSPConfig:
    """Parsed contents of a CSP control file (``csp.ini``).

    Bundles the run inputs (structures, the ``L0..L_max`` schedule, output directory,
    mRNA) with the kinetic + MD :class:`RunParams`. Produced by :func:`read_csp_config`
    and consumed by :func:`csp` / passed to :func:`run_continuous_synthesis`.
    """
    pdb_file: str
    ribosome: str
    L0: int = 1
    L_max: Optional[int] = None
    outdir: str = "synth_out"
    mrna: Optional[str] = None
    codon_time_table_path: Optional[str] = None
    params: RunParams = field(default_factory=RunParams)
    config_file: Optional[str] = None


def read_csp_config(config_file: str, verbose: bool = True) -> CSPConfig:
    """Parse a CSP control file (INI ``[OPTIONS]``) into a :class:`CSPConfig`.

    The MD / ribosome keys configure the shared :class:`cosmo.csp.core.RunParams`
    machinery; the O'Brien **kinetic keys** are added on top. Required: ``pdb_file``,
    ``ribosome``. ``L0`` (default ``1``) and ``L_max`` (default = full residue count)
    are optional. Per-codon timing additionally requires ``mrna`` (``codon_times`` is
    optional -- defaults to the bundled E. coli 310 K table).

    Kinetic keys
    ------------
    - ``mrna`` -- mRNA sequence file (one codon per residue + 1 stop). Required for
      per-codon timing.
    - ``codon_times`` -- either a **path** to a per-codon mean-time table
      (``CODON  seconds``) **or** a **positive number of seconds** (uniform codon time,
      no ``mrna`` needed). Optional; omitting it uses the bundled E. coli 310 K table.
    - ``scale_factor`` -- in-vivo seconds -> in-silico ns compressor.
    - ``time_stage_1`` / ``time_stage_2`` -- mean peptidyl-transfer / translocation
      dwell (s); stage 3 = codon total minus these.
    - ``random_seed`` -- seed for the FPT sampler (reproducible schedules).
    - ``max_steps_per_stage`` / ``min_steps_per_stage`` -- clamp each stage's step count.
    - ``ejection_steps`` / ``dissociation_steps`` -- post-synthesis free runs (0 = skip).

    MD / ribosome keys: ``model`` (nascent force field, default ``hps_kr`` -- it carries
    the O'Brien Rmin/2 table for the ribosome-NC 12-10-6 excluded volume), ``dt``,
    ``ref_t``, ``tau_t``, ``nstout``, ``device``, ``ppn``, ``constraints``,
    ``restraint_k``, ``buffer``, ``minimize``, ``tunnel_wall``, ``ptc_offset``,
    ``optimize_ptc_geometry`` (yes/no; place the A/P restraint targets one peptide bond
    apart and clear of the ribosome excluded volume -- when on, ``ptc_offset`` is unused).
    (``trna_tether`` is not honoured -- CSP forces the position restraint; the wall
    plane is auto-derived; output is always nascent-only.)

    Inline ``#``/``;`` comments are ignored. **Units:** OpenMM defaults.
    """
    def log(msg: str) -> None:
        if verbose:
            print(msg)

    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read CSP config file: {config_file!r}")
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

    log(f"Reading CSP parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    ribosome = req("ribosome")
    L0 = int(opt("L0")) if opt("L0") is not None else 1
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None
    outdir = opt("outdir") or "synth_out"
    mrna = opt("mrna")
    _uniform_codon_time, codon_time_table_path = kinetics.parse_codon_times(opt("codon_times"))

    p = RunParams()
    # --- MD / ribosome knobs ---
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
        p.constraints = None if cons.lower() == "none" else cons
    if opt("restraint_k") is not None:
        p.restraint_k = float(opt("restraint_k"))
    if opt("buffer") is not None:
        p.buffer_nm = float(opt("buffer"))
    if opt("minimize") is not None:
        p.minimize = bool(strtobool(opt("minimize")))
    if opt("tunnel_wall") is not None:
        p.tunnel_wall = bool(strtobool(opt("tunnel_wall")))
    if opt("ptc_offset") is not None:
        p.ptc_offset_nm = float(opt("ptc_offset"))
    if opt("optimize_ptc_geometry") is not None:
        p.optimize_ptc_geometry = bool(strtobool(opt("optimize_ptc_geometry")))

    # --- O'Brien kinetic knobs ---
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
    if opt("ejection_steps") is not None:
        p.ejection_steps = as_int(opt("ejection_steps"))
    if opt("dissociation_steps") is not None:
        p.dissociation_steps = as_int(opt("dissociation_steps"))

    if p.uniform_codon_time is None and mrna is None:
        raise ValueError(f"{config_file}: per-codon timing needs an 'mrna' file "
                         f"(or set 'codon_times' to a positive number of seconds for a "
                         f"uniform codon time). A 'codon_times' table path is optional "
                         f"(defaults to the bundled E. coli 310 K table).")

    log(f"  inputs: pdb_file={pdb_file}, ribosome={ribosome}")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}, model={p.model}")
    if p.uniform_codon_time is not None:
        log(f"  timing: uniform (codon_time={p.uniform_codon_time:g} s)")
    else:
        _table = codon_time_table_path or "bundled E. coli 310 K"
        log(f"  timing: per-codon (mrna={mrna}, codon_times={_table})")
    log(f"          scale_factor={p.scale_factor:g}, time_stage_1={p.time_stage_1:g} s, "
        f"time_stage_2={p.time_stage_2:g} s")
    if p.max_steps_per_stage is not None:
        log(f"          TEST CLAMP max_steps_per_stage={p.max_steps_per_stage} "
            f"(~{3 * p.max_steps_per_stage} steps/residue)")
    log(f"  ribosome: rigid scenery (O'Brien 12-10-6 excluded volume + Yukawa)"
        f"; tunnel wall: {'on (plane auto-derived)' if p.tunnel_wall else 'off'}")
    log(f"  PTC geometry: {'optimized (A/P targets one peptide bond apart, EV-clear)' if p.optimize_ptc_geometry else 'raw tRNA anchors + ptc_offset'}")
    log(f"  integrator: dt={p.dt_ps} ps, ref_t={p.ref_t} K, tau_t={p.tau_t} /ps, nstout={p.nstout}")
    if p.ejection_steps or p.dissociation_steps:
        log(f"  post-synthesis: ejection={p.ejection_steps} steps, "
            f"dissociation={p.dissociation_steps} steps")
    log(f"  hardware/output: device={p.device}, ppn={p.ppn}, outdir={outdir}")

    return CSPConfig(pdb_file=pdb_file, ribosome=ribosome, L0=L0, L_max=L_max,
                     outdir=outdir, mrna=mrna, codon_time_table_path=codon_time_table_path,
                     params=p, config_file=config_file)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def csp(argv: Optional[List[str]] = None) -> None:
    """Console entry point: ``cosmo-csp -f csp.ini``.

    The O'Brien continuous synthesis protocol (per-codon, 3-stage elongation), driven
    by an INI control file (see :func:`read_csp_config`). ``-o`` / ``--device``
    override the output directory / compute device for sweeps.
    """
    import warnings
    warnings.filterwarnings("ignore", category=Warning, module=r"MDAnalysis")

    parser = argparse.ArgumentParser(
        prog="cosmo-csp",
        description="O'Brien Continuous Synthesis Protocol (per-codon, 3-stage "
                    "elongation) in cosmo style, on cosmo's IDP force field. Times "
                    "every residue from its codon, splits it into peptidyl-transfer / "
                    "translocation / tRNA-binding sub-stages, and grows the nascent "
                    "chain N->C on the rigid ribosome. Controlled by an INI file: "
                    "cosmo-csp -f csp.ini",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input", "-f", dest="config", type=str,
                        help="CSP control file (INI, [OPTIONS] section).")
    parser.add_argument("-o", "--outdir", default=None,
                        help="override the output directory from the config file.")
    parser.add_argument("--device", default=None, choices=["CPU", "GPU"],
                        help="override the compute device from the config file.")

    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(argv)
    if not args.config:
        parser.error("a CSP control file is required: -f csp.ini")

    print(f"OpenMM version: {mm.__version__}")

    cfg = read_csp_config(args.config)
    if args.outdir:
        cfg.outdir = args.outdir
    if args.device:
        cfg.params.device = args.device

    run_continuous_synthesis(
        cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
        mrna=cfg.mrna, codon_time_table_path=cfg.codon_time_table_path, params=cfg.params)


if __name__ == "__main__":
    csp()
