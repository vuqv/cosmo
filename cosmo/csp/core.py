"""Core nascent-chain MD engine for protein synthesis (a library, not a runner).

The shared per-length machinery of :mod:`cosmo.csp`, mirroring the sibling ``topo``
project's ``topo/csp/core.py`` but built on cosmo's **sequence-based IDP force field**
(HPS / mpipi) instead of topo's structure-based Gō model. Both the O'Brien 3-stage
:mod:`cosmo.csp.protocol` runner and the analytic-tunnel :mod:`cosmo.csp.cylinder`
runner call into here; the module itself has no loop and no CLI.

What it provides:

- :class:`RunParams` -- the single flat parameter block (MD / rigid-ribosome knobs
  **and** the O'Brien kinetic knobs; the per-length helpers ignore the kinetic fields).
- :func:`run_length` -- build + seed + (restrain) + minimize + run + finalize **one**
  length-``L`` system (nascent chain, optionally + the rigid ribosome scenery). Called
  once per residue by the cylinder runner and three times per residue (stage 1/2/3) by
  the CSP protocol runner.
- :func:`read_anchor` -- read a fixed ribosome bead coordinate (the P-/A-anchors).
- coordinate seeding (:func:`cold_start_positions`, :func:`seed_positions`), the
  C-terminus position restraint (:func:`add_cterm_restraint`), and small I/O helpers.

Because cosmo's force field is sequence-based there is **no STRIDE, no native contact
map and no build-once-subset machinery**: a length-``L`` model is simply
:func:`cosmo.models.buildCoarseGrainModel` on the **first ``L`` residues** of the
nascent PDB (all cosmo forces are sequence-local or pairwise-by-type, so the
restriction is exact).

Like topo, :func:`run_length` wraps each stage in a **per-stage dt-halving stability
guard** (:data:`STABILITY_POTE_LIMIT_KJ` / :data:`STABILITY_MAX_ATTEMPTS`). The
divergence it guards is **non-native excluded volume**, not a native-contact well: the
new residue is seeded at the *fixed* A-site target, and a recently-added residue that has
not yet cleared that region can sit inside the stiff repulsive core of the
ribosome<->nascent 12-10-6 EV or the nascent Ashbaugh-Hatch potential. That stiff EV mode
is unstable at the configured dt, the stage diverges (PotE -> ~1e13 kJ/mol), and -- with
no guard -- the chain explodes and never recovers (observed: sandbox/Yeast asyn, L47
stage 1, new bead 0.22 nm from residue 45). Both the 12-10-6 EV and the AH core are
present in cosmo *and* topo, so this guard is genuinely shared architecture, not
structure-based (Gō) physics. It re-runs a diverging stage at ``dt/2`` with ``2x`` steps
(dwell ``= n_steps * dt`` preserved exactly), which stabilises the integration.

**Units:** OpenMM defaults throughout -- length nm, time ps, energy kJ/mol,
temperature K, force constants kJ/mol/nm^2. In-vivo dwell times (kinetics) are seconds.
"""
from __future__ import annotations

import contextlib
import math
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit

import cosmo
from cosmo import engine
from cosmo.core import models
from cosmo.parameters import model_parameters
from cosmo.utils import runinfo
from cosmo.csp.ribosome import (Ribosome, append_ribosome, add_trna_tether,
                                add_tunnel_wall, TRNA_TETHER_BOND_NM,
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
COLD_START_ZIGZAG_NM = 0.03

# --- per-stage stability guard --------------------------------------------
# The seeded structure is integrated at the configured timestep with flexible (un-
# constrained) bonds. For a few configurations the dynamics diverge (potential energy
# -> ~1e13 kJ/mol) and corrupt that stage's frames. The stiff mode responsible is
# **non-native excluded volume**: the new residue is seeded at the fixed A-site target,
# and a recently-added residue that has not cleared that region can land inside the
# repulsive core of the ribosome<->nascent 12-10-6 EV or the nascent Ashbaugh-Hatch
# potential (observed: sandbox/Yeast asyn L47, new bead 0.22 nm from residue 45). Both
# potentials exist in cosmo and topo, so this is shared architecture (mirrors
# topo.csp.core). Rather than rigid AllBonds constraints (topo's Go path; cosmo keeps
# flexible IDP bonds), a diverging stage is detected and re-run with a halved timestep
# and proportionally more steps -- which preserves the physical in-vivo dwell time
# (dwell = n_steps * dt) exactly while stabilising the integration. Up to
# STABILITY_MAX_ATTEMPTS halvings are tried (dt -> dt/2^k). A finite, physically sane CG
# stage has |PotE| of order 1e1-1e4 kJ/mol, so the 1e9 limit cleanly separates
# "diverged" from "fine" with margin.
STABILITY_POTE_LIMIT_KJ = 1.0e9
STABILITY_MAX_ATTEMPTS = 6

# --- PTC-geometry optimization (always on) ---------------------------------
# O'Brien tRNA resting bond lengths (nm; 4.27 / 4.76 A) and orienting-angle stiffness.
_PTC_D_A_NM, _PTC_D_P_NM = 0.427, 0.476
_PTC_ANGLE_K_KJ = 25.0 * 4.184                  # = 104.6 kJ/mol/rad^2 (O'Brien)

# Verbose banners are off by default (CSP runs many stages; the concise per-stage
# summary line is enough). Set COSMO_CSP_VERBOSE=1 to restore the full build banners.
VERBOSE = os.environ.get("COSMO_CSP_VERBOSE", "").strip() not in ("", "0", "no", "false")


def _vprint(*args, **kwargs) -> None:
    """``print`` only when :data:`VERBOSE` (``COSMO_CSP_VERBOSE``) is set."""
    if VERBOSE:
        print(*args, **kwargs)


@contextlib.contextmanager
def _quiet():
    """Suppress stdout from cosmo's chatty per-length build / engine calls.

    :func:`buildCoarseGrainModel` and :mod:`cosmo.engine` print a detailed
    force-by-force setup log on every call; repeating that for every stage of every
    residue would bury the synthesis progress. Silence it inside the loop and emit a
    concise per-stage summary instead (matching the sibling ``topo`` runner). Honoured
    only when not :data:`VERBOSE`.
    """
    if VERBOSE:
        yield
        return
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


def _ptc_bead_index(ribo, segid: str, resid: int, name: str) -> int:
    """Index of one named ribosome bead in a :class:`Ribosome` (or raise)."""
    for i in range(ribo.n):
        if ribo.segids[i] == segid and ribo.resids[i] == resid and ribo.names[i] == name:
            return i
    raise ValueError(f"ribosome bead {segid}:{resid}@{name} not found.")


def optimal_ptc_targets(ribo, *, peptide_nm: float, aa_rmin_2_nm: float = 0.5,
                        n_starts: int = 60
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Optimal A-site / P-site C-terminus restraint **target points** (nm).

    Ported from ``topo.csp.core.optimal_ptc_targets`` (the always-on PTC-geometry
    path). Returns ``(a_target, p_target)`` -- two **fixed points in space** (nm; NOT
    bonds; the CSP path restrains the C-terminus to a fixed point via
    :func:`add_cterm_restraint`) that sit exactly one **peptide bond** (``peptide_nm``,
    the active model's ``bond_length_protein``) apart and clear of the ribosome excluded
    volume. Seeding
    the new residue at ``a_target`` while the previous residue rests at ``p_target``
    makes the always-present peptide bond start at its equilibrium length -- so no
    residue is delivered ~0.9 nm (the raw tRNA-anchor separation) off equilibrium.

    The points minimize the soft O'Brien restraint energy -- the **A/P tRNA bonds**
    (``_PTC_D_A_NM`` / ``_PTC_D_P_NM``, harmonic ``kb = restraint_k``) and orienting
    **angles/improper** (``_PTC_ANGLE_K_KJ``) -- plus the ribosome ``eps*(Rmin/r)^12``
    excluded volume (``eps = RIBO_NC_EPS_KJ``, ``Rmin_ij = (Rmin/2)_AA + (Rmin/2)_bead``
    sum rule, the SAME EV the simulation applies), subject to:

    - the **peptide bond** as the only **hard** distance constraint (``peptide_nm``);
    - an **exit-side inequality** keeping each point at ``x >= x`` of its tRNA R bead
      (between the tRNA and the +x exit port); and
    - O'Brien's two excluded-volume exclusions (new AA<->AtR:76@R, prev AA<->PtR:76@R),

    over ``n_starts`` deterministic full-sphere (Fibonacci) starts. The tRNA bond
    lengths are soft (finite-k), so ``|A-RA|`` / ``|P-RP|`` come out *near* 0.427 /
    0.476 nm rather than exactly, keeping the solve feasible on any geometry.

    Parameters
    ----------
    ribo : Ribosome
        The parsed rigid CG ribosome (supplies the AtR/PtR 76 R/P/BR2 beads, all bead
        coordinates and the excluded-volume Rmin/2).
    aa_rmin_2_nm : float, optional
        Nascent-bead Rmin/2 (nm) for the seed excluded-volume term (O'Brien sum rule).
        Default 0.5 (conservative, so the seed clears the wall for essentially every
        residue).
    n_starts : int, optional
        Number of deterministic multistart initial orientations (default 60).
    peptide_nm : float
        Equilibrium peptide-bond length (nm) held as the hard constraint (**required**;
        no default, so the model's value is never silently wrong). The CSP runner passes
        the active model's ``bond_length_protein`` -- 0.380 nm for ``hps_kr``, 0.382 for
        ``hps_urry``/``hps_ss``, 0.381 for ``mpipi``.

    Returns
    -------
    (numpy.ndarray, numpy.ndarray)
        ``a_target`` and ``p_target``, each a ``(3,)`` coordinate in nm.
    """
    from scipy.optimize import minimize  # lazy: only needed for this geometry mode

    RB = ribo.coords_nm                                     # all bead coords (nm)
    RBr = np.asarray(ribo.Rmin_2_nm)                        # bead Rmin/2 (nm)
    iRA = _ptc_bead_index(ribo, "AtR", 76, "R"); RA = RB[iRA]
    iRP = _ptc_bead_index(ribo, "PtR", 76, "R"); RP = RB[iRP]
    PA = RB[_ptc_bead_index(ribo, "AtR", 76, "P")]
    PP = RB[_ptc_bead_index(ribo, "PtR", 76, "P")]
    U2A = RB[_ptc_bead_index(ribo, "AtR", 76, "BR2")]       # BR2 == O'Brien PU2
    U2P = RB[_ptc_bead_index(ribo, "PtR", 76, "BR2")]

    R_pair = aa_rmin_2_nm + RBr                             # pair contact dist (nm), sum rule
    mA = np.ones(ribo.n, bool); mA[iRA] = False             # O'Brien exclusions
    mP = np.ones(ribo.n, bool); mP[iRP] = False
    ka, eps = _PTC_ANGLE_K_KJ, RIBO_NC_EPS_KJ
    kb = RESTRAINT_K_KJ                                      # tRNA bond k (kJ/mol/nm^2)

    def _ang(p, q, r):
        u = p - q; v = r - q
        return np.arccos(np.clip(u.dot(v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))

    def _dih(p, q, r, s):
        b1, b2, b3 = q - p, r - q, s - r
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        m = np.cross(n1, b2 / np.linalg.norm(b2))
        return np.arctan2(m.dot(n2), n1.dot(n2))

    def _ev(pt, mask):                                      # O'Brien 12-10-6 EV (kJ/mol)
        d = np.linalg.norm(RB[mask] - pt, axis=1)
        x = R_pair[mask] / d
        return np.sum(eps * (13.0 * x ** 12 - 18.0 * x ** 10 + 4.0 * x ** 6))

    d2r = np.deg2rad
    # Only the peptide bond is a HARD equality (a_target/p_target exactly peptide_nm
    # apart). The two tRNA bonds are SOFT penalties in the objective. Inequalities keep
    # each point on the EXIT side (x >= x_tRNA) so a buried minimum cannot win.
    cons = [{"type": "eq", "fun": lambda x: np.linalg.norm(x[:3] - x[3:]) - peptide_nm},
            {"type": "ineq", "fun": lambda x: x[0] - RA[0]},   # A.x >= A-tRNA.x (exit side)
            {"type": "ineq", "fun": lambda x: x[3] - RP[0]}]   # P.x >= P-tRNA.x (exit side)

    def _obj(x):                                            # kJ/mol
        A, P = x[:3], x[3:]
        E = kb * (np.linalg.norm(A - RA) - _PTC_D_A_NM) ** 2     # soft A-tRNA bond
        E += kb * (np.linalg.norm(P - RP) - _PTC_D_P_NM) ** 2    # soft P-tRNA bond
        E += ka * (_ang(A, RA, PA) - d2r(106)) ** 2 + ka * (_ang(A, RA, U2A) - d2r(127)) ** 2
        E += ka * (_ang(P, RP, PP) - d2r(117)) ** 2 + ka * (_ang(P, RP, U2P) - d2r(130)) ** 2
        E += ka * (_dih(A, RA, PA, U2A) - d2r(128)) ** 2 + ka * (_dih(P, RP, PP, U2P) - d2r(-161)) ** 2
        return E + _ev(A, mA) + _ev(P, mP)

    best = None
    golden = np.pi * (1.0 + 5.0 ** 0.5)
    for t in range(n_starts):
        z = 1.0 - 2.0 * (t + 0.5) / n_starts
        rho = np.sqrt(max(0.0, 1.0 - z * z))
        dirn = np.array([rho * np.cos(golden * t), rho * np.sin(golden * t), z])
        sA = RA + _PTC_D_A_NM * dirn
        sP = sA + peptide_nm * (RP - sA) / np.linalg.norm(RP - sA)
        r = minimize(_obj, np.concatenate([sA, sP]), method="SLSQP",
                     constraints=cons, options={"maxiter": 300, "ftol": 1e-12})
        if r.success and (best is None or r.fun < best.fun):
            best = r
    if best is None:
        raise RuntimeError("optimal_ptc_targets: no feasible restraint geometry found.")
    return best.x[:3], best.x[3:]                           # nm (already)


# --------------------------------------------------------------------------
# Subset native structure (first L residues)
# --------------------------------------------------------------------------
def write_subset_structure(full_pdb: str, L: int, out_pdb: str) -> None:
    """Write a CA-only PDB of the first ``L`` residues of ``full_pdb``.

    This length-``L`` native structure supplies the per-residue identities (and so
    the mass/charge/radius/hydropathy parameters) and the chain connectivity for the
    length-``L`` build; the *simulated* coordinates come from the seeding scheme, not
    from here. Residues are taken in file order so particle ``i`` (``0..L-1``)
    corresponds to native residue ``i+1`` (the C-terminus is ``L-1``).
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
def cold_start_positions(L0: int, anchor: np.ndarray,
                         bond_length_nm: float) -> np.ndarray:
    """Extended cold-start layout for the first length ``L0`` (PLAN.md §6).

    Residues ``1..L0`` are laid along the tunnel axis (+x) from ``anchor``: the
    C-terminus (residue ``L0``) sits *at* the anchor and the N-terminus (residue 1)
    points toward the exit (+x), one CG bond length apart. A small alternating
    transverse zig-zag (:data:`COLD_START_ZIGZAG_NM`) breaks the exact collinearity.
    Returns an ``(L0, 3)`` array in nm (row ``i`` = residue ``i+1``).
    """
    positions = np.empty((L0, 3))
    for i in range(L0):  # i = 0..L0-1  -> native residue i+1
        # residue L0 (i = L0-1) at offset 0 (the anchor); residue 1 furthest +x.
        offset = (L0 - 1 - i) * bond_length_nm
        pos = anchor + offset * TUNNEL_AXIS
        pos = pos + ((-1) ** i) * COLD_START_ZIGZAG_NM * np.array([0.0, 1.0, 0.0])
        positions[i] = pos
    return positions


def seed_positions(prev_final: np.ndarray, seed_point: np.ndarray) -> np.ndarray:
    """Seed length ``L`` from the previous final structure + the new residue.

    Residues ``1..L-1`` keep their coordinates from step ``L-1``'s final structure
    (``prev_final``, shape ``(L-1, 3)`` nm); the new C-terminal residue ``L`` is
    placed at ``seed_point`` -- the equilibrium-bond A-site target from
    :func:`optimal_ptc_targets` -- so the always-present peptide bond ``L-1<->L``
    starts at its equilibrium length and the seeded structure minimizes cleanly.
    Returns an ``(L, 3)`` array in nm.
    """
    new_residue = np.asarray(seed_point, dtype=float)
    return np.vstack([prev_final, new_residue[None, :]])


# --------------------------------------------------------------------------
# C-terminus restraint
# --------------------------------------------------------------------------
def add_cterm_restraint(system: mm.System, particle_index: int,
                        anchor_nm: np.ndarray, k: float) -> mm.Force:
    """Add a harmonic restraint pulling one particle toward a fixed anchor.

    ``U = k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)`` via a ``CustomExternalForce``.
    Only the current C-terminus (``particle_index = L-1``) is restrained -- to the
    stage's target point -- so the tether hand-off between stages/lengths is automatic
    (each rebuilt step restrains only its own C-terminus). The CSP protocol switches
    the target A->P by passing a different ``anchor_nm`` per stage. ``k`` is in OpenMM
    units (kJ/mol/nm^2).

    ``k`` is a **per-particle** parameter (not a global) so this force can coexist with
    the tunnel wall (:func:`cosmo.csp.ribosome.add_tunnel_wall`), whose global constant
    is also called ``k``: two forces sharing a global name with different defaults is an
    OpenMM error, which per-particle ``k`` avoids (position restraint + wall is exactly
    what CSP uses together).
    """
    restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    for p in ("k", "x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    restraint.addParticle(int(particle_index),
                          [float(k), float(anchor_nm[0]), float(anchor_nm[1]),
                           float(anchor_nm[2])])
    system.addForce(restraint)
    return restraint


# --------------------------------------------------------------------------
# I/O helpers
# --------------------------------------------------------------------------
def _write_pdb(topology, positions_nm: np.ndarray, path: str) -> None:
    """Write a PDB from a topology and an ``(N, 3)`` nm coordinate array."""
    coords = [mm.Vec3(float(x), float(y), float(z)) for x, y, z in positions_nm] * unit.nanometer
    with open(path, "w") as fh:
        mm.app.PDBFile.writeFile(topology, coords, fh)


class NascentDCDReporter:
    """A DCD reporter that writes only the first ``n_keep`` atoms each frame.

    Used with the rigid ribosome so the (large, static) ribosome beads are **not**
    written to the trajectory every frame -- only the nascent chain (indices
    ``0..n_keep-1``) is saved. Mirrors :class:`openmm.app.DCDReporter` but slices the
    positions and uses a fixed ``n_keep``-atom topology, so the DCD pairs with the
    nascent PSF.
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
    lists, so it cannot describe the combined system. parmed reads masses/bonds
    straight from the OpenMM topology + System instead.
    """
    import parmed as pmd
    pmd.openmm.load_topology(cgModel.topology, system=cgModel.system).save(
        path, overwrite=True)


def _finalize_nascent(cfg, ctx, nascent_topology, n_keep: int,
                      start_epoch: float) -> None:
    """Finalize a length writing a **nascent-only** final structure.

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
# Run parameters (MD + rigid ribosome + O'Brien kinetics -- one flat block)
# --------------------------------------------------------------------------
@dataclass
class RunParams:
    """All run parameters, set once from the CLI/INI -- the single CSP + elongation
    parameter block.

    Holds both the **MD / ribosome** knobs consumed by the per-length helper
    (:func:`run_length`: timestep, temperature, rigid-ribosome tunnel wall, C-terminus
    restraint, ...) **and** the **O'Brien kinetic** knobs consumed by the high-level
    :func:`cosmo.csp.run_continuous_synthesis` / :func:`cosmo.csp.run_cylinder_synthesis`
    (per-codon timing, stage dwell means, post-synthesis phases). The low-level
    per-length helper simply ignores the kinetic fields.
    """
    # --- MD / integrator ---
    # Nascent-chain force field. Any IDP model works (hps_kr / hps_urry / mpipi): the
    # ribosome<->nascent 12-10-6 excluded volume uses the model-independent O'Brien
    # Rmin/2 collision-radius tables (OBRIEN_RMIN_2_NM / OBRIEN_RNA_RMIN_2_BEADS in
    # model_parameters), decoupled from the force field, while the nascent IDP<->IDP
    # interaction is whatever the selected model provides. hps_kr is just the default.
    model: str = "hps_kr"
    n_steps: int = 1000                # fallback per-segment steps (kinetics override it)
    dt_ps: float = 0.01
    ref_t: float = 300.0
    tau_t: float = 0.01
    nstout: int = 50
    device: str = "CPU"
    ppn: int = 1
    # Bond treatment (forwarded to buildCoarseGrainModel). Default None = flexible
    # (harmonic) bonds -- the physically appropriate choice for this IDP/condensate force
    # field, where backbone flexibility matters (cosmo defaults to flexible where topo's
    # Go model defaults to 'AllBonds'). Set 'AllBonds' for rigid bonds (larger-timestep
    # path); constraints act only on the CA/P pseudo-bonds, so the stiff non-native EV is
    # unaffected -- the per-stage dt-halving guard, not AllBonds, is what prevents EV
    # blow-ups. The always-on PTC optimization still seeds each new residue one peptide
    # bond from the previous C-terminus, so the bond starts at equilibrium either way.
    constraints: object = None
    restraint_k: float = RESTRAINT_K_KJ   # kJ/mol/nm^2 (= 200 kcal/mol/A^2)
    minimize: bool = True
    # --- rigid ribosome scenery (the supplied PDB is always rigid for CSP) ---
    # When a ribosome is present the output is always nascent-only (the rigid scenery is
    # static). ribosome<->nascent excluded volume is the O'Brien 12-10-6 (fixed eps /
    # Rmin/2 from model_parameters) -- no soft-wall knobs.
    trna_tether: bool = True           # O'Brien tRNA tether vs. plain position restraint
    tunnel_wall: bool = True           # one-sided planar wall x >= x0 (forward extrusion)
    tunnel_wall_k: float = TUNNEL_WALL_K  # x0 is always auto-derived per-structure (protocol)

    # ------------------------------------------------------------------
    # O'Brien continuous-synthesis kinetics (used only by the high-level runners).
    # ------------------------------------------------------------------
    scale_factor: float = 4331293.0     # in-vivo seconds -> in-silico ns compressor
    time_stage_1: float = 0.00034       # mean peptidyl-transfer dwell (s)
    time_stage_2: float = 0.004201      # mean translocation dwell (s)
    # Uniform-timing mean codon time (s), or None for per-codon timing from the mRNA.
    uniform_codon_time: Optional[float] = None
    ribosome_traffic: bool = False      # apply the external traffic correction if available
    initiation_rate: float = 0.083333   # translation initiation rate (1/s), traffic only
    random_seed: Optional[int] = None   # seed for the FPT sampler (reproducibility)
    # --- test clamps (production: leave both at their defaults / None) ---
    max_steps_per_stage: Optional[int] = None  # cap each stage (tutorial: small)
    min_steps_per_stage: int = 1               # floor each stage
    # --- post-synthesis phases (steps; 0 = skip) ---
    ejection_steps: int = 0             # release the restraint; let the chain leave
    dissociation_steps: int = 0         # free run; protein drifts off the ribosome


def _make_cfg(out_dir: Path, sub_pdb: str, params: RunParams) -> cosmo.SimulationConfig:
    """Build a per-length/-stage :class:`SimulationConfig` for the engine helpers.

    Each stage is a self-contained standalone run (its own output folder), so this
    mirrors a single ``cosmo-mdrun`` invocation: a constant-temperature production run
    of ``n_steps`` at ``ref_t``. The caller sets ``init_position`` (seed coordinates
    fed as-is) or ``cgModel.positions`` directly. Coordinates are never shifted -- the
    absolute tunnel/anchor frame is deliberate (PLAN.md §4).
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
    cfg.constraints = params.constraints
    cfg.output_dir = str(out_dir)
    cfg.outname = "traj"
    cfg.device = params.device
    cfg.ppn = params.ppn
    cfg.restart = False
    cfg.minimize = False  # we minimize the seeded structure explicitly below
    return cfg


# --------------------------------------------------------------------------
# Single length step (nascent chain; optional rigid ribosome scenery)
# --------------------------------------------------------------------------
def run_length(L: int, *, full_pdb: str,
               p_anchor: np.ndarray, a_anchor: np.ndarray,
               prev_final: Optional[np.ndarray], out_root: Path,
               params: RunParams,
               ribo: Optional[Ribosome] = None,
               wall_x0_nm: float,
               cterm_seed: Optional[np.ndarray] = None,
               seed_override: Optional[np.ndarray] = None,
               restrain: bool = True, out_subdir: Optional[str] = None,
               n_steps_override: Optional[int] = None,
               minimize_override: Optional[bool] = None,
               label: Optional[str] = None) -> np.ndarray:
    """Build, seed, (restrain,) minimize and run one length-``L`` system.

    The System is the nascent chain only when ``ribo`` is ``None``; when ``ribo`` is
    given the rigid (mass-0) ribosome is appended with the ribosome<->nascent
    excluded-volume + electrostatics (:func:`cosmo.csp.ribosome.append_ribosome`).

    The same routine drives an elongation stage and the post-synthesis phase; these
    arguments tailor it:

    - ``p_anchor`` : the C-terminus restraint **target** for this call. The CSP
      protocol switches it A->P across stages (stages 1-2 to the A-target, stage 3 to
      the P-target); the hand-off is automatic since each build restrains only its own
      C-terminus.
    - ``cterm_seed`` : the C-terminus' equilibrium hold point (where its restraint bond
      starts at rest). Falls back to ``p_anchor``.
    - ``seed_override`` : use these ``(L, 3)`` nm coordinates directly (continue a stage
      from the previous one's final, or seed a post-synthesis phase from the finished
      structure) instead of cold-start / new-residue placement.
    - ``restrain`` : if ``False``, drop the C-terminus restraint (ejection /
      dissociation -- the chain is released).
    - ``out_subdir`` : output folder under ``out_root`` (default ``L_<L>``); e.g.
      ``L_007/stage_2`` or ``ejection``.
    - ``n_steps_override`` : run this many steps instead of ``params.n_steps`` (the
      kinetic driver passes the per-stage codon-dwell step count here).
    - ``minimize_override`` : per-call minimize gate (default ``params.minimize``);
      e.g. CSP stage 2 continues from an already-relaxed stage-1 final, so it skips it.
    - ``label`` : short stage tag for the concise per-stage summary line.

    Returns the final nascent ``(L, 3)`` nm coordinate array (seeds the next stage).
    """
    _vprint()
    _vprint("#" * 66)
    _vprint("# " + (f"L={L}  {label}" if label else
                    (f"Nascent length L = {L}"
                     + ("  (+ rigid ribosome)" if ribo is not None else ""))))
    _vprint("#" * 66)

    out_dir = out_root / (out_subdir or f"L_{L:03d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    if cterm_seed is None:
        cterm_seed = p_anchor
    cg_bond_length_nm = float(
        model_parameters.parameters[params.model]["bond_length_protein"])

    # 1. length-L native structure (residue identities + connectivity) -------
    sub_pdb = str(out_dir / f"native_1_{L}.pdb")
    write_subset_structure(full_pdb, L, sub_pdb)

    # 2. build the length-L cosmo model (sequence-based: no contact precompute).
    #    check_forces=False -- the native-subset PDB geometry is not what we simulate
    #    (the seeded coordinates are, and they are minimized explicitly).
    with _quiet():
        cgModel = models.buildCoarseGrainModel(sub_pdb, model=params.model,
                                               minimize=False, check_forces=False,
                                               constraints=params.constraints)
    _vprint(f"[ build ] {sub_pdb}")
    _vprint(f"  chains={cgModel.n_chains}  CA atoms={cgModel.n_atoms}  "
            f"bonds={cgModel.n_bonds}  model={params.model}"
            + (f"  (+ {ribo.n} rigid ribosome beads)" if ribo is not None else ""))

    # 3. seed nascent coordinates --------------------------------------------
    if seed_override is not None:   # continue a stage / seed a post-synthesis phase
        if seed_override.shape[0] != L:
            raise ValueError(
                f"seed_override has {seed_override.shape[0]} residues but L = {L}.")
        nascent_pos = seed_override
    elif prev_final is None:        # cold start (L == L0)
        nascent_pos = cold_start_positions(L, cterm_seed, cg_bond_length_nm)
    else:                           # continue from previous length + new residue
        if prev_final.shape[0] != L - 1:
            raise ValueError(
                f"prev_final has {prev_final.shape[0]} residues but L-1 = {L - 1}.")
        # Seed the new C-terminus at its equilibrium hold point (the A-site target one
        # peptide bond from the previous C-terminus), so the backbone bond to L-1 starts
        # at equilibrium and the seeded structure minimizes cleanly; the restraint/wall
        # then drive forward extrusion.
        if cterm_seed is None:
            raise ValueError("cterm_seed is required to seed a continued length "
                             "(L > L0); the CSP runner always supplies it.")
        nascent_pos = seed_positions(prev_final, cterm_seed)

    cfg = _make_cfg(out_dir, sub_pdb, params)
    if n_steps_override is not None:
        cfg.md_steps = n_steps_override

    # 3b. nascent-only output (ribosome present): capture the nascent (L-atom) topology
    #     and write the nascent PSF BEFORE appending the ribosome (append mutates the
    #     topology; dumpTopology keys off the model's nascent-only per-atom lists).
    nascent_only = ribo is not None
    nascent_topology = None
    if nascent_only:
        nascent_topology = mm.app.PDBFile(sub_pdb).topology
        with _quiet():
            cgModel.dumpTopology(cfg.output_path(".psf"))

    # 3c. append the rigid ribosome (mass-0 scenery + cross-interactions).
    if ribo is not None:
        with _quiet():
            append_ribosome(cgModel, ribo)
        positions = np.vstack([nascent_pos, ribo.coords_nm])
    else:
        positions = nascent_pos

    # 4. hold the current C-terminus (residue L). ribo + trna_tether: O'Brien
    #    peptidyl-tRNA linkage. Otherwise a harmonic position restraint to p_anchor
    #    (the CSP protocol uses the position-restraint path so the target can switch A->P).
    if restrain:
        if ribo is not None and params.trna_tether:
            prev_index = (L - 2) if L >= 2 else None
            add_trna_tether(cgModel, L - 1, prev_index, ribo, L,
                            segid="PtR", resid=76)
        else:
            add_cterm_restraint(cgModel.system, L - 1, p_anchor, params.restraint_k)

    # 4b. tunnel wall: keep nascent beads at x >= x0 (forward-only extrusion).
    if ribo is not None and params.tunnel_wall:
        add_tunnel_wall(cgModel.system, range(L), x0_nm=wall_x0_nm,
                        k=params.tunnel_wall_k)

    # 5. coordinates: nascent-only path via init_position seed.pdb; ribosome path sets
    #    positions directly. Either way coordinates are used as-is (shift_positions=False)
    #    so the absolute tunnel/anchor frame is preserved.
    if ribo is None:
        seed_pdb = str(out_dir / "seed.pdb")
        _write_pdb(cgModel.topology, positions, seed_pdb)
        cfg.init_position = seed_pdb
        with _quiet():
            cgModel.dumpTopology(cfg.output_path(".psf"))
    else:
        cfg.init_position = None
        cgModel.positions = [mm.Vec3(*r) for r in positions] * unit.nanometer
        if not nascent_only:
            with _quiet():
                _dump_topology_psf(cgModel, cfg.output_path(".psf"))

    do_minimize = params.minimize if minimize_override is None else bool(minimize_override)

    # 5b. Stability-guarded stage run (mirrors topo.csp.core). The common case runs once
    #     at the configured dt; a diverging stage -- a non-native EV blow-up, see the
    #     module docstring / STABILITY_* constants -- is re-run at dt/2 with 2x steps
    #     (identical dwell = n_steps * dt) until it integrates cleanly. Divergence is
    #     judged from the **maximum** |PotE| over the stage (a mid-run blow-up can cool
    #     back below the limit by the final frame yet still have ruined those frames), and
    #     the chunked stepping aborts a diverging stage early. Each retry's fresh
    #     setup_simulation + attach_reporters truncates the per-stage output, so a
    #     successful attempt cleanly overwrites the aborted one. The `[stability]` notices
    #     print outside _quiet() so a divergence/retry is visible in the run log.
    base_dt = cfg.dt
    base_steps = cfg.md_steps
    ctx = None
    max_pe = float("nan")
    start = time.time()
    for attempt in range(STABILITY_MAX_ATTEMPTS):
        cfg.dt = base_dt / (2 ** attempt)
        cfg.md_steps = base_steps * (2 ** attempt)
        if attempt > 0:
            print(f"[stability] stage diverged (max|PotE| = {max_pe:.3g} kJ/mol > "
                  f"{STABILITY_POTE_LIMIT_KJ:g}); re-running with dt = "
                  f"{cfg.dt.value_in_unit(unit.picoseconds):g} ps x {cfg.md_steps} "
                  f"steps (identical dwell time; attempt {attempt + 1}/"
                  f"{STABILITY_MAX_ATTEMPTS}).")
        with _quiet():
            ctx = engine.setup_simulation(cfg, cgModel, control_file=None,
                                          shift_positions=False)
        diverged = False
        max_pe = 0.0
        if do_minimize:
            try:
                with _quiet():
                    ctx.simulation.minimizeEnergy()
                    ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)
            except Exception as exc:
                # A NaN during minimization (e.g. a bead seeded inside the stiff
                # ribosome 12-10-6 wall) -> treat as divergence and retry with a halved
                # timestep (same seed; smaller dt lets the wall push it out cleanly).
                print(f"[stability] minimization failed ({type(exc).__name__}: "
                      f"{str(exc).splitlines()[0][:80]}); treating as divergence.")
                diverged = True

        if not diverged:
            with _quiet():
                engine.attach_reporters(cfg, ctx.simulation, total_steps=cfg.md_steps)
                if nascent_only:
                    # Swap the full-system DCD reporter for a nascent-only one (the rigid
                    # ribosome is static -- no need to write its beads every frame).
                    ctx.simulation.reporters[1] = NascentDCDReporter(
                        cfg.output_path(".dcd"), cfg.nstxout, nascent_topology, L)
            # Step in chunks so a divergence is caught (and the stage aborted) mid-run.
            chunk = max(cfg.nstxout, cfg.md_steps // 20, 1)
            done = 0
            while done < cfg.md_steps:
                n = min(chunk, cfg.md_steps - done)
                try:
                    ctx.simulation.step(n)
                except Exception as exc:
                    # OpenMM raises "Particle coordinate is NaN" when a stage blows up
                    # outright; catch it so the guard can retry instead of aborting the run.
                    print(f"[stability] integration blew up ({type(exc).__name__}: "
                          f"{str(exc).splitlines()[0][:80]}).")
                    diverged = True
                    break
                done += n
                pe = abs(ctx.simulation.context.getState(getEnergy=True)
                         .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
                max_pe = max(max_pe, pe)
                if not math.isfinite(pe) or pe > STABILITY_POTE_LIMIT_KJ:
                    diverged = True
                    break
        if not diverged:
            break
    else:
        print(f"[stability][warning] stage '{out_subdir or ('L_%03d' % L)}' still "
              f"diverges after {STABILITY_MAX_ATTEMPTS} attempts (max|PotE| = "
              f"{max_pe:.3g} kJ/mol). Continuing, but this stage's frames are suspect.")
    # Restore the configured values (per-call cfg, but keep it tidy for the summary
    # line + finalize; report the physical base step count, not a retry's inflated one).
    cfg.dt = base_dt
    cfg.md_steps = base_steps

    final_pe = ctx.simulation.context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)
    with _quiet():
        if nascent_only:
            _finalize_nascent(cfg, ctx, nascent_topology, L, start)
        else:
            engine.finalize_simulation(cfg, ctx, cgModel.topology, start)

    # One concise, column-aligned line per stage.
    print(f"  L={L:>3d}  {(label or 'run'):<26s}  {cfg.md_steps:>6d} steps  "
          f"{time.time() - start:>6.2f} s  PE={final_pe:>+13.4e} kJ/mol")

    # 6. final nascent coordinates seed the next stage/length ----------------
    final = mm.app.PDBFile(cfg.output_path("_final.pdb")).getPositions(
        asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]
