"""Protein synthesis: the elongation loop (``cosmo.translation``).

The nascent-chain elongation driver, mirroring the sibling ``topo`` project's
``topo/translation/elongate.py`` but built on cosmo's **sequence-based IDP force
field** (HPS / mpipi) instead of topo's structure-based Gō model. See
``cosmo/translation/PLAN.md``.

**Build step v1 — mechanics only, no ribosome forces.** The simulated System is
the *nascent chain only*. The truncated ribosome is used purely as a source of
two fixed points — the **P-anchor** (P-site tRNA residue-76 ``R`` bead) and the
**A-anchor** (A-site tRNA residue-76 ``R`` bead) — for new-residue placement and
the C-terminus restraint target. (Build step v2 — the rigid ribosome scenery —
is deferred; see PLAN.md §6.)

Because cosmo's force field is sequence-based, there is **no STRIDE, no native
contact map and no build-once-subset machinery** (PLAN.md §2): a length-``L`` model
is simply :func:`cosmo.models.buildCoarseGrainModel` on the **first ``L`` residues**
of the nascent PDB.

What the loop does, for ``L = L0 .. N_full`` (``L`` = current nascent length):

1. **Build** the length-``L`` cosmo model on the first ``L`` residues (bonds,
   Yukawa, the short-range non-bonded; ``hps_ss`` adds the local angle/torsion).
2. **Seed coordinates.** ``L == L0``: lay residues ``1..L0`` extended along the
   tunnel axis (+x) from the P-anchor (C-terminus at the P-anchor, N-terminus
   toward +x), one CG bond length apart. ``L > L0``: residues ``1..L-1`` from the
   previous step's final structure; the new residue ``L`` at the A-anchor + buffer.
3. **Restrain only residue ``L``** (the current C-terminus) to the P-anchor with a
   harmonic ``CustomExternalForce`` (``k = 83680 kJ/mol/nm^2`` = 200 kcal/mol/A^2).
   The hand-off is automatic: each rebuilt step restrains only its own C-terminus.
4. **Minimize** (relax the placement / the stretched new bond), draw Boltzmann
   velocities at ``ref_t``, **run ``n_steps`` steps**, write the per-length outputs,
   and seed ``L+1`` from this step's final structure.

Use it as a CLI, driven by an INI control file (see :func:`read_elongate_config`)::

    cosmo-elongate -f elongate.ini
    python -m cosmo.translation -f elongate.ini

or call :func:`run_elongation` from your own script. Setup / reporters / finalize
are reused from :mod:`cosmo.engine`; this module is the outer loop over length ``L``.
"""
from __future__ import annotations

import argparse
import configparser
import contextlib
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np
import openmm as mm
from openmm import unit

import cosmo
from cosmo import engine
from cosmo.core import models
from cosmo.parameters import model_parameters
from cosmo.utils import runinfo
from cosmo.utils.config import strtobool
from cosmo.translation.ribosome import (Ribosome, load_ribosome, append_ribosome,
                                        add_trna_tether, add_tunnel_wall,
                                        TRNA_TETHER_BOND_NM, TUNNEL_WALL_X0_NM,
                                        TUNNEL_WALL_K, RIBO_NC_EPS_KJ)

# --- physical constants / conversions -------------------------------------
# Tunnel / exit axis: the working ribosome is X-aligned (PLAN.md §5), so the chain
# extrudes toward +x.
TUNNEL_AXIS = np.array([1.0, 0.0, 0.0])
# Default C-terminus position-restraint force constant (kJ/mol/nm^2 = OpenMM units;
# 83680 = 200 kcal/mol/A^2).
RESTRAINT_K_KJ = 83680.0
# Small transverse zig-zag amplitude (nm) for the cold-start layout. A perfectly
# collinear chain makes every bond angle exactly 180 degrees, where an angle-force
# gradient would be singular (-> NaN on the first minimization step); the beads are
# offset by +/- this amount perpendicular to the tunnel axis to break collinearity.
# Tiny relative to the ~0.38 nm CG bond, so the chain stays "extended along the
# axis". (Only the hps_ss model has an angle term, but the offset is harmless for
# the others and keeps the layout model-independent.)
COLD_START_ZIGZAG_NM = 0.03


@contextlib.contextmanager
def _quiet():
    """Suppress stdout from cosmo's chatty per-length build / engine calls.

    cosmo's :func:`buildCoarseGrainModel` and :mod:`cosmo.engine` print a detailed
    force-by-force setup log on every call. Repeating that for each nascent length
    would bury the elongation progress, so we silence it inside the loop and emit a
    concise per-length summary instead (matching the sibling ``topo`` runner, whose
    build path is quiet by design).
    """
    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
        yield


# --------------------------------------------------------------------------
# Anchors
# --------------------------------------------------------------------------
def read_anchor(pdb_file: str, segid: str, resid: int = 76,
                bead: str = "R") -> np.ndarray:
    """Return the coordinate (nm) of one ribosome bead, e.g. a tRNA anchor.

    Parses the (CG) ribosome PDB for the single atom whose **segID** (columns
    73-76), **residue number** and **atom name** match ``segid`` / ``resid`` /
    ``bead``. Used to pick the P-anchor (``segid='PtR'``) and A-anchor
    (``segid='AtR'``) from the truncated ribosome (both keep residue 76). PDB
    coordinates are in angstrom; the returned vector is in nm.

    Raises
    ------
    ValueError
        If no matching atom (or more than one) is found.
    """
    matches = []
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            res_seq = line[22:26].strip()
            seg = line[72:76].strip()
            if name == bead and seg == segid and res_seq == str(resid):
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                matches.append(np.array([x, y, z]) / 10.0)  # angstrom -> nm
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one atom (segid={segid!r}, resid={resid}, "
            f"name={bead!r}) in {pdb_file}, found {len(matches)}.")
    return matches[0]


# --------------------------------------------------------------------------
# Subset native structure (first L residues)
# --------------------------------------------------------------------------
def write_subset_structure(full_pdb: str, L: int, out_pdb: str) -> None:
    """Write a CA-only PDB of the first ``L`` residues of ``full_pdb``.

    This length-``L`` native structure supplies the per-residue identities (and so
    the mass/charge/radius/hydropathy parameters) and the chain connectivity for
    the length-``L`` build; the *simulated* coordinates come from the seeding
    scheme, not from here. Residues are taken in file order so particle ``i``
    (``0..L-1``) corresponds to native residue ``i+1`` (the C-terminus is ``L-1``).
    """
    import MDAnalysis as mda  # local import: heavy, only needed for slicing

    u = mda.Universe(full_pdb)
    ca = u.select_atoms("protein and name CA")
    if len(ca) < L:
        raise ValueError(f"requested L={L} but {full_pdb} has only {len(ca)} CA atoms.")
    ca[:L].write(out_pdb)


# --------------------------------------------------------------------------
# Coordinate seeding
# --------------------------------------------------------------------------
def cold_start_positions(L0: int, p_anchor: np.ndarray,
                         bond_length_nm: float) -> np.ndarray:
    """Extended cold-start layout for the first length ``L0`` (PLAN.md §6).

    Residues ``1..L0`` are laid along the tunnel axis (+x) from the P-anchor: the
    C-terminus (residue ``L0``) sits *at* the P-anchor and the N-terminus (residue
    1) points toward the exit (+x), one CG bond length apart. A small alternating
    transverse zig-zag (``COLD_START_ZIGZAG_NM``) breaks the exact collinearity.
    Returns an ``(L0, 3)`` array in nm (row ``i`` = residue ``i+1``).
    """
    positions = np.empty((L0, 3))
    for i in range(L0):  # i = 0..L0-1  -> native residue i+1
        # residue L0 (i = L0-1) at offset 0 (the P-anchor); residue 1 furthest +x.
        offset = (L0 - 1 - i) * bond_length_nm
        pos = p_anchor + offset * TUNNEL_AXIS
        pos = pos + ((-1) ** i) * COLD_START_ZIGZAG_NM * np.array([0.0, 1.0, 0.0])
        positions[i] = pos
    return positions


def seed_positions(prev_final: np.ndarray, a_anchor: np.ndarray,
                   buffer_nm: float) -> np.ndarray:
    """Seed length ``L`` from the previous final structure + the new residue.

    Residues ``1..L-1`` keep their coordinates from step ``L-1``'s final structure
    (``prev_final``, shape ``(L-1, 3)`` nm); the new C-terminal residue ``L`` is
    placed at the A-anchor offset by ``buffer_nm`` along +x (the buffer clears
    excluded volume so the new bead does not get a huge non-bonded kick). Returns
    an ``(L, 3)`` array in nm.
    """
    new_residue = a_anchor + buffer_nm * TUNNEL_AXIS
    return np.vstack([prev_final, new_residue[None, :]])


# --------------------------------------------------------------------------
# C-terminus restraint
# --------------------------------------------------------------------------
def add_cterm_restraint(system: mm.System, particle_index: int,
                        anchor_nm: np.ndarray, k: float) -> mm.Force:
    """Add a harmonic restraint pulling one particle toward a fixed anchor.

    ``U = k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)`` via a ``CustomExternalForce``.
    Only the current C-terminus (``particle_index = L-1``) is restrained -- to the
    P-anchor -- so the tether hand-off between lengths is automatic (each rebuilt
    step restrains only its own C-terminus). ``k`` is in OpenMM units (kJ/mol/nm^2).
    """
    restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    restraint.addParticle(int(particle_index),
                          [float(anchor_nm[0]), float(anchor_nm[1]), float(anchor_nm[2])])
    system.addForce(restraint)
    return restraint


# --------------------------------------------------------------------------
# PDB writing helper
# --------------------------------------------------------------------------
def _write_pdb(topology, positions_nm: np.ndarray, path: str) -> None:
    """Write a PDB from a topology and an ``(N, 3)`` nm coordinate array."""
    coords = [mm.Vec3(float(x), float(y), float(z)) for x, y, z in positions_nm] * unit.nanometer
    with open(path, "w") as fh:
        mm.app.PDBFile.writeFile(topology, coords, fh)


class NascentDCDReporter:
    """A DCD reporter that writes only the first ``n_keep`` atoms each frame.

    Used in build step v2 so the (large, static) ribosome beads are **not** written
    to the trajectory every frame -- only the nascent chain (indices ``0..n_keep-1``)
    is saved. Mirrors :class:`openmm.app.DCDReporter` but slices the positions and
    uses a fixed ``n_keep``-atom topology, so the DCD pairs with the nascent PSF.
    """

    def __init__(self, file, reportInterval, nascent_topology, n_keep, append=False):
        self._reportInterval = reportInterval
        self._topology = nascent_topology
        self._n = n_keep
        self._append = append
        self._out = open(file, "r+b" if append else "wb")
        self._dcd = None

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {"steps": steps, "periodic": False, "include": ["positions"]}

    def report(self, simulation, state):
        if self._dcd is None:
            self._dcd = mm.app.DCDFile(
                self._out, self._topology, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, self._append)
        positions = state.getPositions(asNumpy=True)[:self._n]
        self._dcd.writeModel(positions)

    def __del__(self):
        try:
            self._out.close()
        except Exception:
            pass


def _dump_topology_psf(cgModel, path: str) -> None:
    """Write a PSF for the full (nascent + ribosome) system via parmed.

    The model's own ``dumpTopology`` keys per-atom mass/charge off its nascent-only
    lists, so it cannot describe the v2 system. parmed reads masses/bonds straight
    from the OpenMM topology + System instead.
    """
    import parmed as pmd
    pmd.openmm.load_topology(cgModel.topology, system=cgModel.system).save(
        path, overwrite=True)


def _finalize_nascent(cfg, ctx, nascent_topology, n_keep: int,
                      start_epoch: float) -> None:
    """Finalize a v2 length writing a **nascent-only** final structure.

    Like :func:`cosmo.engine.finalize_simulation` but the written ``_final.pdb`` is
    only the first ``n_keep`` (nascent) atoms -- the rigid ribosome is dropped. The
    saved checkpoint still holds the **full** system (needed for a correct restart).
    """
    sim = ctx.simulation
    sim.saveCheckpoint(ctx.checkpoint)
    final_pdb = cfg.output_path("_final.pdb")
    pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
    pos = pos[:n_keep].value_in_unit(unit.nanometer)
    _write_pdb(nascent_topology, pos, final_pdb)
    if ctx.runinfo_path is not None:
        runinfo.write_run_end(ctx.runinfo_path, simulation=sim,
                              start_epoch=start_epoch, final_structure=final_pdb)


# --------------------------------------------------------------------------
# Per-length configuration
# --------------------------------------------------------------------------
@dataclass
class ElongationParams:
    """Run parameters shared by every length (set once from the CLI)."""
    model: str = "hps_urry"            # nascent-chain force field (PLAN.md §7)
    n_steps: int = 1000
    dt_ps: float = 0.01
    ref_t: float = 300.0
    tau_t: float = 0.01
    nstout: int = 50
    device: str = "GPU"            # default to GPU; override with device=CPU in the INI
    ppn: int = 1
    restraint_k: float = RESTRAINT_K_KJ   # kJ/mol/nm^2 (= 200 kcal/mol/A^2)
    buffer_nm: float = 0.4
    minimize: bool = True
    # How far into the tunnel (+x) from the P-anchor bead to hold/seed the
    # C-terminus (nm). None -> auto: 0 in v1 (no ribosome bead), 0.4 nm in v2 (so
    # the C-terminus does not sit on the P-tRNA bead and explode).
    ptc_offset_nm: Optional[float] = None
    # Build step v2: append the truncated ribosome as rigid (mass-0) scenery and
    # wire the ribosome<->nascent excluded-volume + electrostatics.
    rigid_ribosome: bool = False
    # Ribosome RNA representation for v2 (PLAN.md §5): 'topo' (3/4-bead P/R/BR,
    # default) or 'cosmo' (1 bead/nucleotide). Drives which CG ribosome to use.
    rna_model: str = "topo"
    # v2: write only the nascent chain to the trajectory / PSF / final structure
    # (the rigid ribosome is static). The checkpoint still holds the full system.
    nascent_only_output: bool = True
    # v2: ribosome<->nascent excluded-volume strength (kJ/mol). Default = O'Brien's
    # deliberately-soft value (see ribosome.RIBO_NC_EPS_KJ): the CG tunnel bore is too
    # tight for a hard wall, so the soft wall lets the chain thread it. Raise only for
    # experimentation (a hard wall jams the chain at the PTC).
    ribo_eps_kj: float = RIBO_NC_EPS_KJ
    # v2: override the excluded-volume radius (nm) of the ribosome RNA beads. The topo
    # P/R/BR value (0.71 nm) is from a different model and is large vs the CG tunnel
    # bore; reducing it shrinks the contact distance so a harder wall (`ribo_eps`) can
    # fit the bore without jamming. None -> use the per-bead defaults.
    ribo_rna_radii_nm: Optional[float] = None
    # v2: make the ribosome participate in the FULL Ashbaugh-Hatch HPS potential (real
    # per-bead sigma/hydropathy, ε=0.8368 kJ/mol) instead of the soft, separate
    # pure-repulsion excluded volume -- a hard, full-size HPS wall giving real radial
    # confinement (the chain stays extended instead of collapsing at the PTC). HPS
    # models only (not mpipi); use a model that parameterises RNA too (e.g. hps_kr).
    # When on, `ribo_eps` / `ribo_rna_radii` are unused. Off -> the soft-wall default.
    ribo_full_hps: bool = False
    # v2: tether the C-terminus to the P-site tRNA the O'Brien way (bond + orienting
    # angle) instead of a plain position restraint. Only used in v2 (needs the tRNA
    # bead); ignored in v1.
    trna_tether: bool = True
    # v2: one-sided planar tunnel wall keeping nascent beads at x >= tunnel_wall_x0
    # (forward-only extrusion). x0 None -> auto = P-anchor x + tether bond length.
    tunnel_wall: bool = True
    tunnel_wall_x0_nm: Optional[float] = None
    tunnel_wall_k: float = TUNNEL_WALL_K


def _make_cfg(out_dir: Path, sub_pdb: str,
              params: ElongationParams) -> cosmo.SimulationConfig:
    """Build a per-length :class:`SimulationConfig` for the engine helpers.

    Each length is a self-contained standalone run (its own output folder), so this
    mirrors a single ``cosmo-mdrun`` invocation: a constant-temperature production
    run of ``n_steps`` at ``ref_t``. The caller sets ``init_position`` (v1, seeded
    coordinates as-is) or ``cgModel.positions`` directly (v2). Coordinates are never
    shifted -- the absolute tunnel/anchor frame is deliberate (PLAN.md §4).
    """
    cfg = cosmo.SimulationConfig()
    cfg.model = params.model
    cfg.md_steps = params.n_steps
    cfg.dt = params.dt_ps * unit.picoseconds
    cfg.nstxout = params.nstout
    cfg.nstlog = params.nstout
    cfg.nstchk = params.nstout
    cfg.nstcomm = None             # no COM removal (a fixed restraint sets the frame)
    cfg.tcoupl = True
    cfg.ref_t = params.ref_t * unit.kelvin
    cfg.tau_t = params.tau_t / unit.picoseconds
    cfg.pbc = False
    cfg.box_dimension = None
    cfg.pdb_file = sub_pdb
    cfg.output_dir = str(out_dir)
    cfg.outname = "traj"
    cfg.device = params.device
    cfg.ppn = params.ppn
    cfg.restart = False
    cfg.minimize = False  # we minimize the seeded structure explicitly below
    return cfg


# --------------------------------------------------------------------------
# Single length step
# --------------------------------------------------------------------------
def run_length(L: int, *, full_pdb: str,
               p_anchor: np.ndarray, a_anchor: np.ndarray,
               prev_final: Optional[np.ndarray], out_root: Path,
               params: ElongationParams, cg_bond_length_nm: float,
               ribo: Optional[Ribosome] = None, anchor_bead: str = "R",
               wall_x0_nm: float = TUNNEL_WALL_X0_NM,
               cterm_seed: Optional[np.ndarray] = None,
               restrain: bool = True) -> np.ndarray:
    """Build, seed, restrain, minimize and run one length-``L`` system.

    Build step v1 (``ribo=None``): the System is the nascent chain only. Build step
    v2 (``ribo`` given): the rigid (mass-0) ribosome is appended with the
    ribosome<->nascent excluded-volume + electrostatics.

    Returns the final nascent ``(L, 3)`` nm coordinate array (seeds the next length).
    """
    print()
    print("#" * 66)
    print(f"# Nascent length L = {L}")
    print("#" * 66)

    out_dir = out_root / f"L_{L:03d}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. length-L native structure (residue identities + connectivity) -------
    sub_pdb = str(out_dir / f"native_1_{L}.pdb")
    write_subset_structure(full_pdb, L, sub_pdb)

    # 2. build the length-L cosmo model (sequence-based: no contact precompute).
    #    check_forces=False -- the native-subset PDB geometry is not what we
    #    simulate (the seeded coordinates are, and they are minimized explicitly).
    #    The builder is chatty; silence it and print a one-line summary instead.
    with _quiet():
        cgModel = models.buildCoarseGrainModel(sub_pdb, model=params.model,
                                               minimize=False, check_forces=False)
    print(f"[ build ] {sub_pdb}")
    print(f"  chains={cgModel.n_chains}  CA atoms={cgModel.n_atoms}  "
          f"bonds={cgModel.n_bonds}  model={params.model}"
          + (f"  (+ {ribo.n} rigid ribosome beads)" if ribo is not None else ""))

    # 3. seed nascent coordinates --------------------------------------------
    # cterm_seed = the C-terminus' equilibrium hold point (so its tether/restraint
    # bond starts at rest length, not stretched ~1 nm from the A-anchor). Falls back
    # to p_anchor for the legacy v1 path.
    if cterm_seed is None:
        cterm_seed = p_anchor
    if prev_final is None:          # cold start (L == L0)
        nascent_pos = cold_start_positions(L, cterm_seed, cg_bond_length_nm)
    else:                           # continue from previous length + new residue
        if prev_final.shape[0] != L - 1:
            raise ValueError(
                f"prev_final has {prev_final.shape[0]} residues but L-1 = {L - 1}.")
        if ribo is not None:
            # v2: seed the new C-terminus at its equilibrium hold point so the tether
            # (or restraint) bond starts at rest length -- avoids the large placement
            # strain (and heating) of dropping it ~1 nm away at the A-anchor. The
            # backbone bond to L-1 + the wall then drive forward extrusion.
            nascent_pos = np.vstack([prev_final, cterm_seed[None, :]])
        else:
            # v1 (no ribosome): A-anchor + buffer, restraint/tether migrates it A->P.
            nascent_pos = seed_positions(prev_final, a_anchor, params.buffer_nm)

    cfg = _make_cfg(out_dir, sub_pdb, params)

    # 3b. v2 nascent-only output: capture the nascent (L-atom) topology and write
    #     the nascent PSF BEFORE appending the ribosome (append mutates the topology;
    #     the model's dumpTopology keys off its nascent-only per-atom lists).
    nascent_only_v2 = ribo is not None and params.nascent_only_output
    nascent_topology = None
    if nascent_only_v2:
        nascent_topology = mm.app.PDBFile(sub_pdb).topology
        with _quiet():
            cgModel.dumpTopology(cfg.output_path(".psf"))

    # 3c. v2: append the rigid ribosome (mass-0 scenery + cross-interactions).
    if ribo is not None:
        with _quiet():
            append_ribosome(cgModel, ribo, nc_eps_kj=params.ribo_eps_kj,
                            full_hps=params.ribo_full_hps)
        positions = np.vstack([nascent_pos, ribo.coords_nm])
    else:
        positions = nascent_pos

    # 4. hold the current C-terminus (residue L). v2 + trna_tether: O'Brien
    #    peptidyl-tRNA linkage (bond + CA-CA-tRNA orienting angle). Otherwise a
    #    generic harmonic position restraint to the P-anchor target.
    if restrain:
        if ribo is not None and params.trna_tether:
            prev_index = (L - 2) if L >= 2 else None
            add_trna_tether(cgModel, L - 1, prev_index, ribo, L,
                            segid="PtR", resid=76, bead=anchor_bead)
        else:
            add_cterm_restraint(cgModel.system, L - 1, p_anchor, params.restraint_k)

    # 4b. v2 tunnel wall: keep nascent beads at x >= x0 (forward-only extrusion).
    if ribo is not None and params.tunnel_wall:
        add_tunnel_wall(cgModel.system, range(L), x0_nm=wall_x0_nm,
                        k=params.tunnel_wall_k)

    # 5. coordinates: v1 via init_position seed.pdb; v2 set positions directly.
    #    Either way coordinates are used as-is (shift_positions=False) so the
    #    absolute tunnel/anchor frame is preserved.
    if ribo is None:
        seed_pdb = str(out_dir / "seed.pdb")
        _write_pdb(cgModel.topology, positions, seed_pdb)
        cfg.init_position = seed_pdb
        with _quiet():
            cgModel.dumpTopology(cfg.output_path(".psf"))
    else:
        cfg.init_position = None
        cgModel.positions = [mm.Vec3(*r) for r in positions] * unit.nanometer
        if not nascent_only_v2:
            with _quiet():
                _dump_topology_psf(cgModel, cfg.output_path(".psf"))

    start = time.time()
    with _quiet():
        ctx = engine.setup_simulation(cfg, cgModel, control_file=None,
                                      shift_positions=False)
        if params.minimize:
            ctx.simulation.minimizeEnergy()
            ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        engine.attach_reporters(cfg, ctx.simulation, total_steps=cfg.md_steps)
        if nascent_only_v2:
            # Swap the full-system DCD reporter for a nascent-only one (the rigid
            # ribosome is static -- no need to write its beads every frame).
            ctx.simulation.reporters[1] = NascentDCDReporter(
                cfg.output_path(".dcd"), cfg.nstxout, nascent_topology, L)
        ctx.simulation.step(cfg.md_steps)
        if nascent_only_v2:
            _finalize_nascent(cfg, ctx, nascent_topology, L, start)
        else:
            engine.finalize_simulation(cfg, ctx, cgModel.topology, start)
    print(f"  ran {cfg.md_steps} steps ({params.device})"
          f"{' + minimize' if params.minimize else ''}"
          f"{' (+ rigid ribosome)' if ribo is not None else ''} "
          f"-> {cfg.output_path('_final.pdb')}  ({time.time() - start:.1f} s)")

    # 6. final nascent coordinates seed the next length ----------------------
    final = mm.app.PDBFile(cfg.output_path("_final.pdb")).getPositions(
        asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]


# --------------------------------------------------------------------------
# Elongation loop
# --------------------------------------------------------------------------
def run_elongation(full_pdb: str, ribosome_pdb: str, *,
                   L0: int, L_max: Optional[int] = None,
                   out_root: str = "synth_out",
                   params: Optional[ElongationParams] = None) -> None:
    """Run the full nascent-chain elongation loop ``L = L0 .. L_max``.

    Build step v1 (``params.rigid_ribosome = False``): nascent chain only. Build
    step v2 (``True``): the truncated ribosome is appended as rigid (mass-0) scenery
    with the ribosome<->nascent excluded-volume + electrostatics.

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein (the nascent chain at full length).
    ribosome_pdb : str
        Truncated CG ribosome PDB -- the source of the P-/A-anchor coordinates
        (PtR/AtR residue-76 ``R`` beads), and in v2 the rigid scenery itself.
    L0 : int
        Starting nascent-chain length (cold-start layout).
    L_max : int, optional
        Final length; defaults to the full residue count ``N_full``.
    out_root : str
        Root output directory; each length writes to ``<out_root>/L_<L>/``.
    params : ElongationParams, optional
        Per-length run parameters (defaults to the test settings).
    """
    if params is None:
        params = ElongationParams()

    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # Anchors (fixed points from the truncated ribosome). The tRNA-76 attachment
    # bead is the ribose ``R`` in the topo representation, but the single ``P`` bead
    # in the cosmo 1-bead representation.
    anchor_bead = "R" if params.rna_model == "topo" else "P"
    p_anchor = read_anchor(ribosome_pdb, "PtR", resid=76, bead=anchor_bead)
    a_anchor = read_anchor(ribosome_pdb, "AtR", resid=76, bead=anchor_bead)
    print(f"P-anchor (PtR 76 {anchor_bead}): {p_anchor} nm")
    print(f"A-anchor (AtR 76 {anchor_bead}): {a_anchor} nm")

    # Where the C-terminus is held / seeded, measured into the tunnel (+x) from the
    # P-anchor bead. Auto: 0 in v1 (no ribosome bead to clash with); 0.4 nm in v2 so
    # the restraint does not pin the C-terminus on top of the P-tRNA bead (-> blowup).
    offset = params.ptc_offset_nm
    if offset is None:
        offset = 0.4 if params.rigid_ribosome else 0.0
    p_target = p_anchor + offset * TUNNEL_AXIS
    if offset:
        print(f"C-terminus target: P-anchor + {offset} nm (+x) = {p_target} nm")

    # The C-terminal-AA addition plane (PTC) = P-anchor x + tether bond length: the
    # default tunnel-wall plane.
    ptc_plane_x = float(p_anchor[0]) + TRNA_TETHER_BOND_NM
    wall_x0 = (params.tunnel_wall_x0_nm if params.tunnel_wall_x0_nm is not None
               else ptc_plane_x)

    # Equilibrium seed point for the C-terminus: where its hold force is at rest, so
    # placement injects no large strain (topo's planned fix). With the tRNA tether the
    # rest length is the tether bond from the P-site bead (down the tunnel, +x); with
    # a position restraint it is the restraint target p_target itself.
    if params.rigid_ribosome and params.trna_tether:
        cterm_seed = p_anchor + TRNA_TETHER_BOND_NM * TUNNEL_AXIS
    else:
        cterm_seed = p_target

    # Build step v2: load the rigid ribosome once (identical in every length).
    ribo = None
    if params.rigid_ribosome:
        ribo = load_ribosome(ribosome_pdb, model=params.model,
                             rna_model=params.rna_model,
                             rna_radii=params.ribo_rna_radii_nm)
        _sr = ("full Ashbaugh-Hatch HPS (hard wall + hydropathy)"
               if params.ribo_full_hps else "soft separate excluded volume")
        print(f"Rigid ribosome: {ribo.n} beads from {ribosome_pdb} "
              f"(mass-0 scenery; ribosome<->nascent short-range = {_sr}; Yukawa on; "
              f"rna_model={params.rna_model}).")
        print(f"  C-terminus: {'tRNA tether (bond + angle)' if params.trna_tether else 'position restraint'}"
              f"; seeded at its rest position {np.round(cterm_seed, 3)} nm")
        if params.ribo_full_hps:
            print("  ribosome short-range: full HPS (epsilon = 0.8368 kJ/mol, real per-bead "
                  "sigma/hydropathy); ribo_eps / ribo_rna_radii ignored")
        else:
            print(f"  excluded-volume eps = {params.ribo_eps_kj} kJ/mol"
                  f"{', RNA radii override = %.3f nm' % params.ribo_rna_radii_nm if params.ribo_rna_radii_nm is not None else ''} "
                  f"({'soft wall: chain threads the tight CG bore' if params.ribo_eps_kj < 0.01 else 'HARD wall: needs RNA radii < bore or it jams'})")
        if params.tunnel_wall:
            print(f"  tunnel wall: x >= {wall_x0:.3f} nm (k={params.tunnel_wall_k} kJ/mol/nm^2)")

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
          f"model={params.model}, {params.n_steps} steps/residue, "
          f"ribosome={'on (v2)' if ribo is not None else 'off (v1)'}.")

    prev_final: Optional[np.ndarray] = None
    for L in range(L0, L_max + 1):
        prev_final = run_length(
            L, full_pdb=full_pdb, p_anchor=p_target, a_anchor=a_anchor,
            prev_final=prev_final, out_root=out_path, params=params,
            cg_bond_length_nm=cg_bond_length_nm, ribo=ribo,
            anchor_bead=anchor_bead, wall_x0_nm=wall_x0, cterm_seed=cterm_seed)

    print()
    print(f"Done. Elongated {L0} -> {L_max}. Per-length outputs under {out_path}/")


# --------------------------------------------------------------------------
# INI control file
# --------------------------------------------------------------------------
@dataclass
class ElongateConfig:
    """Parsed contents of an elongation control file (``elongate.ini``)."""
    pdb_file: str
    ribosome: str
    L0: int
    L_max: Optional[int] = None
    outdir: str = "synth_out"
    params: ElongationParams = field(default_factory=ElongationParams)
    config_file: Optional[str] = None


def read_elongate_config(config_file: str, verbose: bool = True) -> ElongateConfig:
    """Parse an elongation control file (INI) into an :class:`ElongateConfig`.

    The file has a single ``[OPTIONS]`` section. Required keys: ``pdb_file``,
    ``ribosome``, ``L0``. All other keys are optional and fall back to the defaults
    in :class:`ElongationParams` / :class:`ElongateConfig`:

    - ``pdb_file`` -- full native PDB of the target protein (the nascent chain).
    - ``ribosome`` -- truncated CG ribosome PDB (source of the P-/A-anchors).
    - ``L0`` -- starting nascent-chain length (cold-start layout).
    - ``L_max`` -- final length (blank -> full residue count).
    - ``outdir`` -- root output directory (per-length subfolders ``L_<L>/``).
    - ``model`` -- nascent-chain force field: ``hps_urry`` (default), ``hps_kr``,
      ``hps_ss`` or ``mpipi``.
    - ``n_steps`` -- integration steps per residue (constant schedule).
    - ``dt`` -- time step (ps); ``ref_t`` -- temperature (K); ``tau_t`` -- Langevin
      friction (1/ps); ``nstout`` -- trajectory/log/checkpoint write frequency.
    - ``device`` -- 'CPU' or 'GPU'; ``ppn`` -- CPU threads (device = CPU).
    - ``restraint_k`` -- C-terminus position-restraint constant (kJ/mol/nm^2).
    - ``buffer`` -- new-residue placement buffer beyond the A-anchor (nm).
    - ``ptc_offset`` -- hold/seed the C-terminus this far (+x) from the P-anchor (nm).
    - ``minimize`` -- yes/no, per-step energy minimization.
    - ``rigid_ribosome`` -- yes/no, build step v2: append the rigid ribosome + the
      ribosome<->nascent excluded-volume + electrostatics.
    - ``rna_model`` -- ribosome RNA representation for v2: 'topo' (3/4-bead P/R/BR,
      default) or 'cosmo' (1 bead/nucleotide).
    - ``nascent_only_output`` -- v2: write only the nascent chain to traj/PSF/final
      (default yes; the rigid ribosome is static).
    - ``trna_tether`` -- v2: hold the C-terminus with the O'Brien tRNA tether (bond +
      orienting angle) vs. a position restraint (default yes).
    - ``tunnel_wall`` -- v2: one-sided planar wall keeping nascent beads at
      ``x >= tunnel_wall_x0`` (forward-only extrusion; default yes).
    - ``tunnel_wall_x0`` -- wall plane (nm). **Leave blank for the auto value**, the
      C-terminal-AA addition plane (PTC) = ``P-anchor x + tether bond length``
      (``TRNA_TETHER_BOND_NM`` = 0.476 nm); set only to override.
    - ``tunnel_wall_k`` -- wall force constant (kJ/mol/nm^2; default 8368 = 20 kcal/mol/A^2).
    - ``ribo_eps`` -- ribosome<->nascent excluded-volume strength (kJ/mol). Default is
      O'Brien's deliberately-soft value (the CG tunnel bore is too tight for a hard
      wall -- it would jam the chain); raise only for experimentation.
    - ``ribo_full_hps`` -- yes/no: make the ribosome participate in the **full
      Ashbaugh-Hatch HPS** potential (real per-bead sigma/hydropathy, hard wall +
      attraction) instead of the soft separate excluded volume. HPS models only;
      use a model parameterising RNA too (e.g. ``hps_kr``). When on, ``ribo_eps`` /
      ``ribo_rna_radii`` are unused. Default no.

    Inline comments starting with ``#`` or ``;`` are ignored; underscores in
    ``n_steps`` (e.g. ``1_000``) are allowed. Paths are resolved relative to the
    current working directory (as for ``md.ini``).

    **Units:** OpenMM defaults throughout -- length nm, time ps, energy kJ/mol,
    temperature K, force constants kJ/mol/nm^2.
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

    log(f"Reading elongation parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    ribosome = req("ribosome")
    L0 = int(req("L0"))
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None
    outdir = opt("outdir") or "synth_out"

    p = ElongationParams()
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
    if opt("buffer") is not None:
        p.buffer_nm = float(opt("buffer"))
    if opt("ptc_offset") is not None:
        p.ptc_offset_nm = float(opt("ptc_offset"))
    if opt("minimize") is not None:
        p.minimize = bool(strtobool(opt("minimize")))
    if opt("rigid_ribosome") is not None:
        p.rigid_ribosome = bool(strtobool(opt("rigid_ribosome")))
    if opt("rna_model") is not None:
        p.rna_model = opt("rna_model")
    if opt("nascent_only_output") is not None:
        p.nascent_only_output = bool(strtobool(opt("nascent_only_output")))
    if opt("trna_tether") is not None:
        p.trna_tether = bool(strtobool(opt("trna_tether")))
    if opt("tunnel_wall") is not None:
        p.tunnel_wall = bool(strtobool(opt("tunnel_wall")))
    if opt("tunnel_wall_x0") is not None:
        p.tunnel_wall_x0_nm = float(opt("tunnel_wall_x0"))
    if opt("tunnel_wall_k") is not None:
        p.tunnel_wall_k = float(opt("tunnel_wall_k"))
    if opt("ribo_eps") is not None:
        p.ribo_eps_kj = float(opt("ribo_eps"))
    if opt("ribo_rna_radii") is not None:
        p.ribo_rna_radii_nm = float(opt("ribo_rna_radii"))
    if opt("ribo_full_hps") is not None:
        p.ribo_full_hps = bool(strtobool(opt("ribo_full_hps")))

    log(f"  inputs: pdb_file={pdb_file}, ribosome={ribosome}")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}, "
        f"model={p.model}, n_steps={p.n_steps}")
    log(f"  mechanics: restraint_k={p.restraint_k} kJ/mol/nm^2, buffer={p.buffer_nm} nm, "
        f"ptc_offset={p.ptc_offset_nm} nm, minimize={p.minimize}")
    log(f"  ribosome forces (v2): {'on (rigid)' if p.rigid_ribosome else 'off (v1)'}"
        + (f", rna_model={p.rna_model}, "
           f"output={'nascent-only' if p.nascent_only_output else 'full system'}, "
           f"tether={'on' if p.trna_tether else 'off'}, "
           f"wall={('x>=%.2f nm' % p.tunnel_wall_x0_nm) if (p.tunnel_wall and p.tunnel_wall_x0_nm is not None) else ('auto' if p.tunnel_wall else 'off')}"
           if p.rigid_ribosome else ""))
    log(f"  integrator: dt={p.dt_ps} ps, ref_t={p.ref_t} K, tau_t={p.tau_t} /ps, nstout={p.nstout}")
    log(f"  hardware/output: device={p.device}, ppn={p.ppn}, outdir={outdir}")

    return ElongateConfig(pdb_file=pdb_file, ribosome=ribosome, L0=L0, L_max=L_max,
                          outdir=outdir, params=p, config_file=config_file)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def elongate(argv: Optional[List[str]] = None) -> None:
    """Console entry point: ``cosmo-elongate -f elongate.ini``.

    The simulation is controlled by an INI file (see :func:`read_elongate_config`).
    ``-o`` / ``--device`` are optional overrides handy for sweeps; everything else
    lives in the control file.
    """
    parser = argparse.ArgumentParser(
        prog="cosmo-elongate",
        description="Protein synthesis elongation loop (build step v1: nascent "
                    "chain only, no ribosome forces). Grows the nascent chain N->C "
                    "one residue per step, restraining the current C-terminus to "
                    "the ribosome P-anchor. Controlled by an INI file: "
                    "cosmo-elongate -f elongate.ini",
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

    run_elongation(cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max,
                   out_root=cfg.outdir, params=cfg.params)


if __name__ == "__main__":
    elongate()
