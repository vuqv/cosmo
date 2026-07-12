"""Continuous Synthesis Protocol (O'Brien), ported to cosmo (``cosmo.csp``).

The **full O'Brien continuous-synthesis protocol** -- the per-codon, 3-stage
elongation cycle of ``continuous_synthesis_v6.py`` -- expressed in cosmo style, on
cosmo's sequence-based IDP force field (HPS / mpipi). It mirrors the sibling ``topo``
project's ``topo/csp/protocol.py`` but drops topo's structure-based Gō machinery
(STRIDE, native-contact precompute, ``domain.yaml`` nscales): a length-``L`` model is
just :func:`cosmo.models.buildCoarseGrainModel` on the first ``L`` residues (all cosmo
forces are sequence-local or pairwise-by-type, so the restriction is exact).

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
The C-terminus is held either by a **position restraint** to the A/P target point
(default) or, with ``trna_tether = yes``, by the **O'Brien tRNA tether** (bond +
orienting angles + improper to the A-site tRNA beads in stages 1-2, the P-site beads
in stage 3). Either way the supplied ribosome is rigid scenery and the tunnel wall +
excluded volume + electrostatics are on.

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

from cosmo.csp.core import (RunParams, read_anchor, run_length,
                            optimal_ptc_targets)
from cosmo.csp.ribosome import load_ribosome
from cosmo.parameters import model_parameters
from cosmo.utils.config import strtobool
from cosmo.csp import kinetics
from cosmo.csp import resume as resume_mod
from cosmo.csp import synth_mrna


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
        (``params.uniform_codon_time`` set). (The INI also accepts
        ``fastest``/``slowest``/``median`` to auto-build a synonymous-codon mRNA; that is
        resolved to a written file in :func:`read_csp_config` before this function is called.)
    codon_time_table_path : str
        Per-codon mean-time table (required for per-codon timing; pick one under
        ``assets/csp/codon_dwell_times/``). There is no bundled default.
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
    # Two C-terminus restraint paths (selected by ep.trna_tether, from the INI):
    #  - position restraint (default): a moving harmonic spring to the A/P target POINT,
    #    switched A->A->P over the 3 stages;
    #  - O'Brien tRNA tether (trna_tether = yes): a bond + orienting angles + improper to
    #    the A-site tRNA beads (stages 1-2) then the P-site tRNA beads (stage 3), which
    #    reproduces O'Brien's orientation control. The A<->P switch is by tether site
    #    (tether_segid below), so the tether supports translocation too.
    if ep.trna_tether:
        print("C-terminus restraint: O'Brien tRNA tether (bond + 2 angles + improper; "
              "A-site stages 1-2, P-site stage 3).")
    else:
        print("C-terminus restraint: position restraint to the A/P target point "
              "(a->a->p migration).")

    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # --- rigid ribosome (always: the supplied file is rigid scenery) --------
    # Loaded once; identical at every length. A CG ribosome PDB (O'Brien 3/4-bead
    # P/R/BR rep, prepared with the sibling topo package); per-bead Rmin/2 + charge from
    # model_parameters. Loaded before the anchors, which are its beads.
    ribo = load_ribosome(ribosome_pdb, model=ep.model)
    print(f"Rigid ribosome: {ribo.n} beads from {ribosome_pdb} "
          f"(mass-0 scenery; ribosome<->nascent 12-10-6 excluded volume + Yukawa on).")

    # --- anchors (fixed points from the truncated ribosome) -----------------
    p_anchor = read_anchor(ribosome_pdb, "PtR", resid=76, bead="R")
    a_anchor = read_anchor(ribosome_pdb, "AtR", resid=76, bead="R")
    print(f"P-anchor (PtR 76 R): {p_anchor} nm")
    print(f"A-anchor (AtR 76 R): {a_anchor} nm")

    # The PTC restraint targets (a_target/p_target) and the tunnel-wall plane are set
    # below, in the fresh/resume branch: on a fresh start they are solved once
    # (optimal_ptc_targets, SLSQP) and written into the schedule-file header; on resume
    # they are re-read from that header so the restraint geometry is pinned identical to
    # the residues already on disk (see cosmo.csp.resume).

    # --- schedule bounds ----------------------------------------------------
    N_full = len(mda.Universe(full_pdb).select_atoms("protein and name CA"))
    if L_max is None:
        L_max = N_full
    if not (1 <= L0 <= L_max <= N_full):
        raise ValueError(f"require 1 <= L0 <= L_max <= N_full; got L0={L0}, "
                         f"L_max={L_max}, N_full={N_full}.")

    dwell_log = out_path / "dwell_times.dat"

    # --- decide fresh vs. resume --------------------------------------------
    # The schedule (per-residue step counts) and the PTC restraint geometry are the
    # immutable plan: drawn/solved once at a fresh start, persisted to dwell_times.dat,
    # then re-read verbatim on resume (no RNG, no SLSQP) so an interrupted run continues
    # with a kinetic schedule + restraint geometry identical to the uninterrupted run.
    # See cosmo.csp.resume.
    if ep.resume not in ("auto", "yes", "no"):
        raise ValueError(f"resume must be 'auto', 'yes' or 'no'; got {ep.resume!r}.")
    if ep.resume == "yes":
        if not resume_mod.progress_exists(out_path):
            raise SystemExit(
                f"[resume] resume = yes but no {resume_mod.PROGRESS_FILENAME} in "
                f"{out_path}/ (nothing to resume). Use resume = no to start fresh.")
        do_resume = True
    elif ep.resume == "auto":
        do_resume = resume_mod.progress_exists(out_path) and dwell_log.is_file()
    else:   # "no"
        do_resume = False

    prog = None
    resume_L = L0
    prev_final: Optional[np.ndarray] = None
    if do_resume:
        # Re-read the immutable plan; recover progress; drop the <=1 in-flight unit.
        schedule, a_target, p_target, wall_x0 = resume_mod.read_schedule(dwell_log)
        resume_mod.schedule_covers(schedule, L0, L_max)
        prog = resume_mod.read_progress(out_path)
        resume_mod.verify_completed_units(out_path, prog, L0)
        dropped = resume_mod.drop_running_units(out_path, prog)
        resume_L = max(prog.last_done_residue + 1, L0)
        if prog.last_done_residue >= L0:
            prev_final = resume_mod.load_final_pdb(
                resume_mod.residue_final_path(out_path, prog.last_done_residue))
        print(f"[resume] {prog.last_done_residue} length(s) complete on disk"
              + (f"; dropped in-flight {dropped}" if dropped else "")
              + f"; continuing from L={resume_L}.")
    else:
        # Fresh start: solve the PTC geometry once, draw the whole schedule from the
        # seeded RNG, then persist both (schedule rows + PTC header) and open progress.
        pep_nm = float(model_parameters.parameters[ep.model]["bond_length_protein"])
        a_target, p_target = optimal_ptc_targets(ribo, peptide_nm=pep_nm)
        wall_x0 = float(min(a_target[0], p_target[0]))
        intrinsic, real, codons = kinetics.build_codon_time_lists(
            L_max + 1, uniform_codon_time=ep.uniform_codon_time,
            mrna_path=mrna, codon_time_table_path=codon_time_table_path,
            ribosome_traffic=ep.ribosome_traffic, initiation_rate=ep.initiation_rate)
        rng = random.Random(ep.random_seed)
        schedule = []
        for L in range(L0, L_max + 1):
            (s1, s2, s3), (t1, t2, t3) = kinetics.stage_steps(
                L, intrinsic, real, time_stage_1=ep.time_stage_1,
                time_stage_2=ep.time_stage_2, scale_factor=ep.scale_factor,
                dt_ps=ep.dt_ps, rng=rng, max_steps_per_stage=ep.max_steps_per_stage,
                min_steps_per_stage=ep.min_steps_per_stage)
            codon = codons[L - 1] if codons is not None else "uniform"
            schedule.append(resume_mod.SchedRow(
                L=L, codon=codon, t_total=intrinsic[L],
                times=(t1, t2, t3), steps=(s1, s2, s3)))
        resume_mod.write_schedule(dwell_log, schedule, ep, a_target, p_target, wall_x0)
        resume_mod.write_progress_header(out_path)

    # --- PTC restraint targets + tunnel wall (from whichever branch above) --
    print(f"[PTC geometry] optimal PTC targets "
          f"(|A-P| = {np.linalg.norm(a_target - p_target):.4f} nm; fixed points):")
    print(f"  A-site target (new-AA seed + stage-1/2 restraint): {a_target} nm")
    print(f"  P-site target (prev-AA / stage-3 restraint)       : {p_target} nm")
    if ep.tunnel_wall:
        print(f"Tunnel wall plane: x >= {wall_x0:.4f} nm (auto: lower C-terminus hold plane).")

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

    # --- up-front cost report: exact step total, nominal wall-time ----------
    total_steps = (sum(sum(r.steps) for r in schedule)
                   + max(ep.ejection_steps, 0) + max(ep.dissociation_steps, 0))
    print(f"[schedule] {L_max - L0 + 1} residues, {total_steps:,} planned MD steps"
          f"{resume_mod.est_walltime(total_steps, ep)}")
    print(f"Per-residue dwell-time table: {dwell_log}")
    print(f"Synthesizing {full_pdb}: L = {L0} .. {L_max} (N_full = {N_full}), "
          f"model={ep.model}.")

    # --- main loop: one residue = three sub-stages, one L_<L>/ dir each ------
    for L in range(resume_L, L_max + 1):
        (s1, s2, s3) = schedule[L - L0].steps
        codon = schedule[L - L0].codon
        print(f"L={L:>3d}  {codon:>5s}  dwell {schedule[L - L0].t_total:>9.4g} s  "
              f"steps {s1:>4d}/{s2:>4d}/{s3:>4d}")

        resume_mod.append_progress(out_path, f"L_{L:03d}", "RUNNING")

        ldir = f"L_{L:03d}"
        cold = prev_final is None
        # Per-stage tRNA-tether site (used only when ep.trna_tether): the new residue is
        # held at the A-site in stages 1-2 and translocated to the P-site in stage 3. At
        # cold start (no previous residue) the first residue is held at the P-site. In
        # stage 1 the previous residue (L-1) is additionally tethered to the P-site so
        # both ends of the new peptide bond sit at the equilibrium PTC.
        stage1_segid = "PtR" if cold else "AtR"
        stage1_prev_segid = None if cold else "PtR"
        # Consolidated layout: all three stages share one L_<L>/ directory -- shared
        # traj.psf + native_1_L.pdb, per-stage traj_s{1,2,3}.dcd, one folded
        # traj_runinfo.log, and a single traj_final.pdb (stage 3 only).
        # Stage 1: deliver the new residue at the A-site and restrain there (cold start
        # lays the initial segment from the P-anchor instead).
        stage1_anchor = p_target if cold else a_target
        f1 = run_length(L, full_pdb=full_pdb, p_anchor=stage1_anchor, a_anchor=a_anchor,
                        prev_final=prev_final, out_root=out_path, params=ep, ribo=ribo,
                        wall_x0_nm=wall_x0,
                        cterm_seed=stage1_anchor, restrain=True,
                        out_subdir=ldir, outname="traj_s1", persist_final=False,
                        n_steps_override=s1, tether_segid=stage1_segid,
                        tether_prev_segid=stage1_prev_segid,
                        label="stage 1 peptidyl-transfer")

        # Stage 2: continue from stage 1, still held at the A-site. Skip the redundant
        # minimization (the seeded structure is stage 1's already-relaxed final).
        f2 = run_length(L, full_pdb=full_pdb, p_anchor=stage1_anchor, a_anchor=a_anchor,
                        prev_final=None, seed_override=f1, out_root=out_path, params=ep,
                        ribo=ribo, wall_x0_nm=wall_x0,
                        cterm_seed=stage1_anchor, restrain=True,
                        out_subdir=ldir, outname="traj_s2", persist_final=False,
                        n_steps_override=s2, tether_segid=stage1_segid,
                        minimize_override=False,
                        label="stage 2 translocation")

        # Stage 3: translocate A->P (restrain the C-terminus to the P-target). Its final
        # is the one persisted traj_final.pdb -- the seed for L+1 and the resume target.
        f3 = run_length(L, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
                        prev_final=None, seed_override=f2, out_root=out_path, params=ep,
                        ribo=ribo, wall_x0_nm=wall_x0,
                        cterm_seed=p_target, restrain=True,
                        out_subdir=ldir, outname="traj_s3", persist_final=True,
                        n_steps_override=s3, tether_segid="PtR",
                        label="stage 3 tRNA-binding")
        prev_final = f3
        # The DONE line is the commit point: all three stages of L are on disk.
        resume_mod.append_progress(out_path, f"L_{L:03d}", "DONE")

    print()
    print(f"Done. Synthesized {L0} -> {L_max}. Per-residue outputs under {out_path}/L_<L>/")
    print(f"Per-residue dwell-time table: {dwell_log}")

    # --- post-synthesis: ejection then dissociation (both free runs) --------
    # Each phase is its own progress unit; on resume a completed phase is skipped and its
    # final structure reloaded to seed the next phase.
    if ep.ejection_steps > 0:
        if do_resume and prog.is_done("ejection"):
            prev_final = resume_mod.load_final_pdb(
                resume_mod.phase_final_path(out_path, "ejection"))
            print("[resume] ejection already complete; skipping.")
        else:
            print()
            print(f"=== Ejection (L = {L_max}, {ep.ejection_steps} steps, restraint OFF) "
                  f"-> {out_path / 'ejection'}/ ===")
            resume_mod.append_progress(out_path, "ejection", "RUNNING")
            prev_final = run_length(
                L_max, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
                prev_final=None, seed_override=prev_final, out_root=out_path, params=ep,
                ribo=ribo, wall_x0_nm=wall_x0, restrain=False,
                out_subdir="ejection", n_steps_override=ep.ejection_steps, label="ejection")
            resume_mod.append_progress(out_path, "ejection", "DONE")

    if ep.dissociation_steps > 0:
        if do_resume and prog.is_done("dissociation"):
            print("[resume] dissociation already complete; skipping.")
        else:
            print()
            print(f"=== Dissociation (L = {L_max}, {ep.dissociation_steps} steps, restraint OFF) "
                  f"-> {out_path / 'dissociation'}/ ===")
            resume_mod.append_progress(out_path, "dissociation", "RUNNING")
            run_length(
                L_max, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
                prev_final=None, seed_override=prev_final, out_root=out_path, params=ep,
                ribo=ribo, wall_x0_nm=wall_x0, restrain=False,
                out_subdir="dissociation", n_steps_override=ep.dissociation_steps,
                label="dissociation")
            resume_mod.append_progress(out_path, "dissociation", "DONE")


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
    are optional. Per-codon timing requires both ``mrna`` and a ``codon_times`` table
    path (there is no bundled default -- pick one under ``assets/csp/codon_dwell_times/``).

    Kinetic keys
    ------------
    - ``mrna`` -- mRNA sequence file (one codon per residue + 1 stop). Required for
      per-codon timing. May also be the keyword ``fastest`` / ``slowest`` / ``median`` to
      auto-build a synonymous-codon mRNA (each residue's fastest/slowest/median-dwell-time
      codon per the ``codon_times`` table, written next to the PDB); a real filename must
      not be ``fastest``/``slowest``/``median``.
    - ``codon_times`` -- either a **path** to a per-codon mean-time table
      (``CODON  seconds  amino_acid``; required for per-codon, no bundled default -- pick
      one under ``assets/csp/codon_dwell_times/``) **or** a **positive number of seconds**
      (uniform codon time, no ``mrna`` needed). A table filename must not be a bare number.
    - ``scale_factor`` -- in-vivo seconds -> in-silico ns compressor.
    - ``time_stage_1`` / ``time_stage_2`` -- mean peptidyl-transfer / translocation
      dwell (s); stage 3 = codon total minus these.
    - ``random_seed`` -- seed for the FPT sampler (reproducible schedules).
    - ``max_steps_per_stage`` / ``min_steps_per_stage`` -- clamp each stage's step count.
    - ``ejection_steps`` / ``dissociation_steps`` -- post-synthesis free runs (0 = skip).
    - ``resume`` -- ``auto`` (default) / ``yes`` / ``no``: continue an interrupted run.

    MD / ribosome keys: ``model`` (nascent IDP force field, default ``hps_kr``; any IDP
    model works -- the ribosome-NC 12-10-6 excluded volume is decoupled from it, always
    using the O'Brien Rmin/2 collision radii, and a model lacking its own ``Rmin_2``
    (e.g. ``hps_urry``, ``mpipi``) transparently falls back to the shared hps_kr table
    for the steric radius only), ``dt``,
    ``ref_t``, ``tau_t``, ``nstout``, ``device``, ``ppn``, ``constraints`` (default
    ``None`` -- flexible bonds), ``restraint_k``, ``minimize``, ``tunnel_wall``,
    ``trna_tether`` (default off -- the O'Brien tRNA tether vs. the position restraint).
    (The wall plane is auto-derived; the PTC geometry is always optimized -- each new
    residue is seeded/held one peptide bond from the previous C-terminus, EV-clear;
    output is always nascent-only.)

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
    # mrna = "fastest"/"slowest"/"median": one-shot prep -- write the synonymous-codon
    # mRNA next to the PDB and hand that file to the normal per-codon path. Reserved
    # keywords, so a real mRNA filename must not be "fastest"/"slowest"/"median".
    if mrna is not None and mrna.strip().lower() in synth_mrna.SYNTHETIC_MRNA_MODES:
        _mode = mrna.strip().lower()
        if _uniform_codon_time is not None:
            raise ValueError(
                f"{config_file}: mrna={_mode} is incompatible with a uniform "
                f"'codon_times' (a number, {_uniform_codon_time:g} s): {_mode} "
                f"picks a per-amino-acid synonymous codon, which needs a codon-time "
                f"*table* (e.g. one under assets/csp/codon_dwell_times/), not a single "
                f"uniform time.")
        if codon_time_table_path is None:
            raise ValueError(
                f"{config_file}: mrna={_mode} needs a 'codon_times' table path -- it "
                f"defines which synonymous codon is {_mode}.")
        mrna = synth_mrna.write_synthetic_mrna(pdb_file, codon_time_table_path, _mode)
        log(f"  mrna={_mode}: wrote synonymous-codon mRNA -> {mrna}")

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
    if opt("minimize") is not None:
        p.minimize = bool(strtobool(opt("minimize")))
    if opt("tunnel_wall") is not None:
        p.tunnel_wall = bool(strtobool(opt("tunnel_wall")))
    # C-terminus restraint mode: position restraint (default) vs. O'Brien tRNA tether.
    if opt("trna_tether") is not None:
        p.trna_tether = bool(strtobool(opt("trna_tether")))

    # --- O'Brien kinetic knobs ---
    if opt("scale_factor") is not None:
        p.scale_factor = float(opt("scale_factor"))
    if opt("time_stage_1") is not None:
        p.time_stage_1 = float(opt("time_stage_1"))
    if opt("time_stage_2") is not None:
        p.time_stage_2 = float(opt("time_stage_2"))
    p.uniform_codon_time = _uniform_codon_time
    # Retired keys -> point users at the single codon_times key.
    for _legacy in ("trans_times", "uniform_ta", "uniform_mfpt"):
        if opt(_legacy) is not None:
            raise ValueError(
                f"{config_file}: '{_legacy}' has been replaced by the single "
                f"'codon_times' key -- set it to a codon-time table path (per-codon "
                f"timing) or to a positive number of seconds (uniform codon time).")
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
    # Resume policy: auto (default; resume iff an interrupted run is present), yes
    # (require a resumable run), no (always fresh). See cosmo.csp.resume.
    if opt("resume") is not None:
        r = opt("resume").strip().lower()
        if r not in ("auto", "yes", "no"):
            raise ValueError(f"{config_file}: resume must be 'auto', 'yes' or 'no', "
                             f"got {opt('resume')!r}.")
        p.resume = r

    # Per-codon timing needs both the (protein-specific) mRNA and an explicit codon-time
    # table -- there is no bundled default (pick one under assets/csp/codon_dwell_times/).
    if p.uniform_codon_time is None and (mrna is None or codon_time_table_path is None):
        raise ValueError(f"{config_file}: per-codon timing needs both an 'mrna' file and "
                         f"a 'codon_times' table path (e.g. one under "
                         f"assets/csp/codon_dwell_times/), or set 'codon_times' to a "
                         f"positive number of seconds for a uniform codon time.")

    log(f"  inputs: pdb_file={pdb_file}, ribosome={ribosome}")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}, model={p.model}")
    if p.uniform_codon_time is not None:
        log(f"  timing: uniform (codon_time={p.uniform_codon_time:g} s)")
    else:
        log(f"  timing: per-codon (mrna={mrna}, codon_times={codon_time_table_path})")
    log(f"          scale_factor={p.scale_factor:g}, time_stage_1={p.time_stage_1:g} s, "
        f"time_stage_2={p.time_stage_2:g} s")
    if p.max_steps_per_stage is not None:
        log(f"          TEST CLAMP max_steps_per_stage={p.max_steps_per_stage} "
            f"(~{3 * p.max_steps_per_stage} steps/residue)")
    log(f"  ribosome: rigid scenery (O'Brien 12-10-6 excluded volume + Yukawa)"
        f"; tunnel wall: {'on (plane auto-derived)' if p.tunnel_wall else 'off'}")
    log(f"  PTC geometry: optimized (A/P targets one peptide bond apart, EV-clear)")
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
    parser.add_argument("--no-resume", "--fresh", dest="no_resume", action="store_true",
                        help="ignore any interrupted run and start fresh "
                             "(overrides resume in the config file).")

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
    if args.no_resume:
        cfg.params.resume = "no"

    run_continuous_synthesis(
        cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
        mrna=cfg.mrna, codon_time_table_path=cfg.codon_time_table_path, params=cfg.params)


if __name__ == "__main__":
    csp()
