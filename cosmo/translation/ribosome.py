"""Rigid ribosome scenery + cross-interactions for elongation build step v2.

Build step v1 (:mod:`cosmo.translation.elongate`) simulates the **nascent chain
only** and uses the truncated ribosome merely as two fixed anchor *coordinates*.
**Build step v2** appends the truncated ribosome to the System as **rigid (mass-0)
scenery** and wires the two ribosome ↔ nascent-chain interactions (PLAN.md §6):

1. **Append** the ribosome beads at indices ``L..N-1`` with **mass = 0** (frozen;
   not integrated), coordinates as-is. The P-/A-anchors become real beads.
2. **Nascent short-range force** (Ashbaugh–Hatch for HPS / Wang–Frenkel for
   mpipi): give the ribosome beads a dummy ``addParticle`` entry and restrict the
   force to the interaction group ``{nascent}×{nascent}`` so the ribosome never
   participates. For ``mpipi`` the dummy is an **in-range ``id = 0``** (its tabulated
   potential is id-indexed; PLAN.md §2.1).
3. **Ribosome–NC excluded volume:** a separate ``CustomNonbondedForce`` with the
   pure repulsion ``ε·(σ_ij/r)¹²`` (``ε = 0.000132 kcal/mol``,
   ``σ_ij = ½(σ_i+σ_j)``), cutoff 2.0 nm / switch 1.8 nm, interaction group
   ``{nascent}×{ribosome}``. Nascent σ = the residue collision diameters
   (``rf_sigma``); ribosome σ = the per-bead ``radii`` (§ local table for P/R/BR).
4. **Electrostatics:** extend the existing Yukawa force with the ribosome charges
   (rRNA phosphate −1e, charged residues) restricted to ``{nascent}×{nascent}`` +
   ``{nascent}×{ribosome}`` (no intra-ribosome electrostatics).

The ribosome is held rigid, so **no intra-ribosome forces are ever computed**. The
nascent chain uses flexible bonds and the ribosome has none, so no constraint ever
involves a mass-0 particle.

Mirrors the sibling ``topo`` project's ``topo/translation/ribosome.py``, adapted to
cosmo's HPS/mpipi forces and its **dual RNA representation** (PLAN.md §5): the
P/R/BR bead parameters are kept **local** here (scenery-only), while ribosomal
protein Cα beads and any single-bead nucleotides look up cosmo's own
``model_parameters``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit
from openmm.app.element import carbon as _CARBON

from cosmo.parameters.model_parameters import parameters as MODEL_PARAMS

# Ribosome-NC excluded-volume strength (kJ/mol) for the pure ε·(σ/r)¹² repulsion.
#
# O'Brien's value (= topo): 0.000132 kcal/mol = 0.00055 kJ/mol. This is *deliberately*
# ~4500x smaller than kT, i.e. a very SOFT wall: two beads feel ~kT only after
# overlapping to ~0.5 σ. The reason is geometric -- in the truncated CG ribosome the
# exit-tunnel bore is only ~0.4-0.6 nm in radius while the nascent<->ribosome contact
# distance is sigma ~ 0.66 nm, so a nascent bead does NOT fit through the bore without
# overlapping the wall beads. The soft epsilon lets the chain *thread* the
# geometrically-tight CG tunnel (feeling weak steric guidance) instead of jamming.
# Raising epsilon toward cosmo's own steric scale (~0.8 kJ/mol) turns the ribosome
# into a hard wall that the chain cannot pass -- it jams at the PTC and the run heats
# up / blows up. So keep this small by default; it is exposed as the `ribo_eps` INI
# key only for experimentation (e.g. with a wider-bore / higher-resolution ribosome).
RIBO_NC_EPS_KJ = 0.000132 * 4.184
# NC excluded-volume cutoff / switch (nm).
RIBO_NC_CUTOFF_NM = 2.0
RIBO_NC_SWITCH_NM = 1.8

# O'Brien peptidyl-tRNA tether (P-site resting geometry): a harmonic bond
# CA(L)--tRNA(P-site) plus an orienting angle CA(L-1)--CA(L)--tRNA. All force
# constants in OpenMM units (kJ/mol/nm^2); kcal/mol/A^2 -> kJ/mol/nm^2 is x418.4.
TRNA_TETHER_BOND_NM = 0.476
TRNA_TETHER_BOND_K = 83680.0          # kJ/mol/nm^2 (= 200 kcal/mol/A^2)
# Orienting-angle = the double-Gaussian backbone-angle form (identical to cosmo's
# hps_ss `addGaussianAngleForces`); its basins theta_alpha = 1.6 rad (~91.7 deg)
# and theta_beta = 2.27 rad (~130 deg) ARE O'Brien's tether-angle parameters. A
# dedicated one-angle force is used so the tether works for every nascent model
# (the hps_ss backbone force is only present for that model).
_TETHER_ANGLE_ENERGY = ("(-1/gamma)*log("
                        "exp(-gamma*(k_alpha*(theta-theta_alpha)^2+eps_alpha))"
                        "+exp(-gamma*k_beta*(theta-theta_beta)^2))")
_TETHER_ANGLE_PARAMS = {            # values in kJ/mol-based default units, rad
    "gamma": 0.0239, "eps_alpha": 17.9912, "theta_alpha": 1.6,
    "theta_beta": 2.27, "k_alpha": 445.1776, "k_beta": 110.0392,
}

# Planar tunnel wall (O'Brien): a one-sided restraint U = k*min(x-x0, 0)^2 keeping
# every nascent bead at x >= x0, so the chain can only extrude forward (+x, toward
# the exit). x0 is the C-terminal-AA addition plane (PTC / P-site) ~ P-anchor x +
# tether bond length; the runner sets it from the P-anchor when not given.
TUNNEL_WALL_X0_NM = 1.05
TUNNEL_WALL_K = 8368.0                # kJ/mol/nm^2 (= 20 kcal/mol/A^2)

# Local, scenery-only RNA bead parameters for the `rna_model = topo` representation
# (3/4 beads per nucleotide). Kept here, NOT in cosmo.parameters.model_parameters,
# so the IDP force-field tables stay clean (PLAN.md §5): these beads only ever
# contribute excluded volume + charge to the rigid ribosome, never the HPS/mpipi
# pair potentials. Values follow O'Brien (radii = collision diameter 0.71 nm;
# phosphate carries −1e, ribose/base rings are neutral).
RIBO_RNA_BEADS = {
    "P":  {"radii": 0.71, "charge": -1.0},   # phosphate
    "R":  {"radii": 0.71, "charge": 0.0},    # ribose-ring centroid
    "BR": {"radii": 0.71, "charge": 0.0},    # base-ring centroid (BR1/BR2 -> BR)
}

# Where to source a residue collision diameter when the nascent model has none
# (mpipi residues carry no per-residue `radii`): use the HPS radii as a universal
# steric table -- the bead's physical size is independent of its IDP force field.
_RADII_FALLBACK_MODEL = "hps_urry"


@dataclass
class Ribosome:
    """Parsed rigid CG ribosome: per-bead coordinates (nm) and force parameters."""
    coords_nm: np.ndarray       # (M, 3)
    radii_nm: List[float]       # excluded-volume sigma per bead (nm)
    charges: List[float]
    hps: List[float]            # Ashbaugh-Hatch hydropathy (lambda) per bead
    names: List[str]
    resnames: List[str]
    resids: List[int]
    segids: List[str]

    @property
    def n(self) -> int:
        return len(self.radii_nm)


def _lookup_residue(model: str, key: str) -> Tuple[float, float, float]:
    """Return ``(radii_nm, charge, hps)`` for a protein/nucleotide bead by key.

    Charge / hydropathy come from the active ``model`` when it defines the residue;
    the collision diameter (radii) from the active model if present, else from a
    steric fallback table (``hps_kr``/``hps_urry``) -- needed because ``mpipi``
    carries no per-residue ``radii`` and ``hps_urry`` carries no RNA residues. ``hps``
    (the Ashbaugh-Hatch lambda) is 0.0 when the table has none (e.g. mpipi).
    """
    e = MODEL_PARAMS[model].get(key)
    if isinstance(e, dict) and "radii" in e:
        return e["radii"], e["charge"], e.get("hps", 0.0)
    for m in ("hps_kr", "hps_urry"):
        fe = MODEL_PARAMS[m].get(key)
        if isinstance(fe, dict) and "radii" in fe:
            charge = e["charge"] if isinstance(e, dict) else fe["charge"]
            hps = (e.get("hps") if isinstance(e, dict) and "hps" in e
                   else fe.get("hps", 0.0))
            return fe["radii"], charge, hps
    raise ValueError(f"no parameters for residue/bead {key!r} in model {model!r} "
                     f"or the steric fallback tables.")


def _bead_params(model: str, rna_model: str, name: str, resname: str,
                 rna_radii: Optional[float] = None) -> Tuple[float, float, float]:
    """Return ``(radii_nm, charge, hps)`` for one CG ribosome bead.

    Protein Cα beads and (``rna_model='cosmo'``) single-bead nucleotides look up
    cosmo's ``model_parameters`` by residue name; (``rna_model='topo'``) P/R/BR
    beads use the local scenery-only :data:`RIBO_RNA_BEADS` table (hydropathy 0 --
    they are scenery-only excluded volume). ``rna_radii``, if given, overrides the
    excluded-volume radius of every **RNA** bead (the topo P/R/BR value of 0.71 nm
    comes from a different model and is geometrically large for the CG tunnel bore;
    reducing it lets a harder wall fit without jamming). Protein Cα radii unaffected.
    """
    if name == "CA":                                   # protein (not RNA)
        return _lookup_residue(model, resname)
    # --- RNA bead ---
    if rna_model == "cosmo":                            # 1 bead/nucleotide
        r, q, h = _lookup_residue(model, resname)
    else:                                               # topo P/R/BR
        btype = name.rstrip("0123456789") or resname
        if btype in RIBO_RNA_BEADS:
            d = RIBO_RNA_BEADS[btype]
            r, q, h = d["radii"], d["charge"], 0.0
        else:
            r, q, h = _lookup_residue(model, btype)
    if rna_radii is not None:
        r = rna_radii
    return r, q, h


def load_ribosome(pdb_file: str, model: str = "hps_urry",
                  rna_model: str = "topo",
                  rna_radii: Optional[float] = None) -> Ribosome:
    """Parse a (truncated) CG ribosome PDB into a :class:`Ribosome`.

    Reads each ATOM/HETATM record and looks up the excluded-volume radius (σ) and
    charge (:func:`_bead_params`). ``model`` is the nascent-chain force field (used
    for protein-Cα / single-nucleotide lookups); ``rna_model`` selects how RNA beads
    are interpreted (``topo`` = P/R/BR via the local table; ``cosmo`` = 1 bead per
    nucleotide via the model's nucleic params). ``rna_radii`` (nm), if given,
    overrides every RNA bead's excluded-volume radius. Coordinates -> nm.
    """
    coords, radii, charges, hps = [], [], [], []
    names, resnames, resids, segids = [], [], [], []
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            resid = int(line[22:26])
            seg = line[72:76].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            r, q, h = _bead_params(model, rna_model, name, resname, rna_radii=rna_radii)
            coords.append((x / 10.0, y / 10.0, z / 10.0))
            radii.append(r); charges.append(q); hps.append(h)
            names.append(name); resnames.append(resname)
            resids.append(resid); segids.append(seg)
    if not coords:
        raise ValueError(f"no ATOM records parsed from ribosome file {pdb_file!r}.")
    return Ribosome(np.asarray(coords, dtype=float), radii, charges, hps,
                    names, resnames, resids, segids)


def _nascent_sigmas(cgModel) -> List[float]:
    """Per-nascent-bead collision diameter (nm) for the ribosome-NC excluded volume.

    HPS models carry these as ``rf_sigma``; mpipi has none, so fall back to the HPS
    per-residue radii (steric size is independent of the IDP force field).
    """
    if isinstance(cgModel.rf_sigma, list):
        return list(cgModel.rf_sigma)
    tbl = MODEL_PARAMS[_RADII_FALLBACK_MODEL]
    resnames = [a.residue.name for a in cgModel.topology.atoms()]
    return [tbl[rn]["radii"] for rn in resnames[:cgModel.n_atoms]]


def append_ribosome(cgModel, ribo: Ribosome,
                    nc_eps_kj: float = RIBO_NC_EPS_KJ,
                    full_hps: bool = False) -> Tuple[List[int], List[int]]:
    """Append the rigid ribosome to a built nascent cosmo model (system + topology).

    Mutates ``cgModel`` in place. Must be called **after**
    :func:`cosmo.models.buildCoarseGrainModel` (the nascent forces must exist).
    Returns ``(nascent_indices, ribosome_indices)``.

    The ribosome<->nascent **short-range** interaction has two modes:

    - ``full_hps=False`` (default): the ribosome does **not** participate in the IDP
      pair potential (dummy Ashbaugh-Hatch / Wang-Frenkel entries restricted to
      nascent x nascent); a **separate, soft** pure-repulsion ``eps*(sigma/r)^12``
      (``nc_eps_kj``) provides the ribosome<->nascent excluded volume. This soft wall
      lets the chain thread the tight CG bore without jamming.
    - ``full_hps=True`` (Ashbaugh-Hatch / HPS models only): the ribosome beads are
      given their **real** ``(sigma, hps)`` HPS parameters and the Ashbaugh-Hatch
      force is extended to nascent x ribosome -- so the ribosome is a **hard, full
      HPS wall** (the standard ε=0.8368 kJ/mol, real collision radii, hydropathy
      attraction). No separate soft excluded-volume force is added. Use a model that
      parameterises both protein and RNA (e.g. ``hps_kr``), so the rRNA beads carry
      proper radii/hydropathy. ``nc_eps_kj`` is ignored in this mode.
    """
    system = cgModel.system
    topology = cgModel.topology
    L = cgModel.n_atoms
    M = ribo.n
    N = L + M
    nascent_idx = list(range(L))
    ribo_idx = list(range(L, N))

    # 1. Rigid (mass-0) ribosome particles. OpenMM does not integrate mass-0
    #    particles, so they stay fixed; no constraints involve them.
    for _ in range(M):
        system.addParticle(0.0)

    # 2. Nascent short-range force. Default: dummy ribosome entries restricted to
    #    nascent x nascent (the ribosome stays out of the IDP pair potential).
    #    full_hps: give the ribosome its real (sigma, hps) and extend the
    #    Ashbaugh-Hatch force to nascent x ribosome -> a hard full-HPS wall.
    if cgModel.ashbaugh_HatchForce is not None:
        sr = cgModel.ashbaugh_HatchForce
        if full_hps:
            for i in range(M):
                sr.addParticle((ribo.radii_nm[i], ribo.hps[i]))   # real HPS params
        else:
            for _ in range(M):
                sr.addParticle((0.1, 0.0))      # dummy (sigma, hps); never evaluated
    elif cgModel.wang_Frenkel_Force is not None:
        if full_hps:
            raise ValueError("full_hps is only supported for the Ashbaugh-Hatch (HPS) "
                             "models, not the Wang-Frenkel (mpipi) potential.")
        sr = cgModel.wang_Frenkel_Force
        for _ in range(M):
            sr.addParticle((0,))                # dummy in-range id=0 (PLAN.md §2.1)
    else:
        raise RuntimeError("nascent model has no Ashbaugh-Hatch or Wang-Frenkel force.")
    sr.addInteractionGroup(nascent_idx, nascent_idx)
    if full_hps:
        sr.addInteractionGroup(nascent_idx, ribo_idx)   # ribosome<->nascent via real HPS

    # 3. Yukawa electrostatics: add ribosome charges; restrict to nascent-nascent +
    #    nascent-ribosome (no intra-ribosome electrostatics).
    yf = cgModel.yukawaForce
    for q in ribo.charges:
        yf.addParticle((q,))
    yf.addInteractionGroup(nascent_idx, nascent_idx)
    yf.addInteractionGroup(nascent_idx, ribo_idx)

    # 4. Ribosome-NC excluded volume: a SEPARATE soft pure (sigma/r)^12 force, only
    #    when not full_hps (in full_hps mode the Ashbaugh-Hatch force above already
    #    provides the ribosome<->nascent steric, so this is skipped).
    if not full_hps:
        nc = mm.CustomNonbondedForce("eps*(sigma/r)^12; sigma=0.5*(sigma1+sigma2)")
        nc.addGlobalParameter("eps", nc_eps_kj)
        nc.addPerParticleParameter("sigma")
        for s in _nascent_sigmas(cgModel):          # nascent collision diameters (nm)
            nc.addParticle((s,))
        for s in ribo.radii_nm:                      # ribosome collision diameters (nm)
            nc.addParticle((s,))
        nc.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        nc.setUseSwitchingFunction(True)
        nc.setSwitchingDistance(RIBO_NC_SWITCH_NM * unit.nanometer)
        nc.setCutoffDistance(RIBO_NC_CUTOFF_NM * unit.nanometer)
        nc.addInteractionGroup(nascent_idx, ribo_idx)
        # OpenMM (CPU platform) requires every CustomNonbondedForce to share an
        # identical exclusion list. The nascent Yukawa carries the bonded (1-2,1-3)
        # exclusions; copy them here so the lists match (they are nascent-nascent
        # pairs, irrelevant to this force's {nascent}x{ribosome} group).
        for k in range(yf.getNumExclusions()):
            i, j = yf.getExclusionParticles(k)
            nc.addExclusion(i, j)
        system.addForce(nc)

    # 5. Topology: append a ribosome chain per segID, grouping beads into residues.
    _append_topology(topology, ribo)

    return nascent_idx, ribo_idx


def _append_topology(topology, ribo: Ribosome) -> None:
    """Append the ribosome beads to an OpenMM topology (one chain per segID)."""
    chains: dict = {}
    residues: dict = {}
    for i in range(ribo.n):
        seg = ribo.segids[i] or "RB"
        chain = chains.get(seg)
        if chain is None:
            chain = topology.addChain(id=seg)
            chains[seg] = chain
        rkey = (seg, ribo.resids[i])
        res = residues.get(rkey)
        if res is None:
            res = topology.addResidue(ribo.resnames[i], chain, id=str(ribo.resids[i]))
            residues[rkey] = res
        topology.addAtom(ribo.names[i], _CARBON, res)


def bead_system_index(ribo: Ribosome, L_nascent: int, segid: str,
                      resid: int, bead: str = "R") -> Optional[int]:
    """System index of a named ribosome bead (e.g. the P-site tRNA anchor).

    Ribosome beads are appended after the ``L_nascent`` nascent particles in load
    order, so the system index is ``L_nascent + (position in the Ribosome arrays)``.
    Returns None if no such bead exists.
    """
    for i in range(ribo.n):
        if ribo.segids[i] == segid and ribo.resids[i] == resid and ribo.names[i] == bead:
            return L_nascent + i
    return None


def add_trna_tether(cgModel, cterm_index: int, prev_index: Optional[int],
                    ribo: Ribosome, L_nascent: int, *, segid: str = "PtR",
                    resid: int = 76, bead: str = "R",
                    bond_length_nm: float = TRNA_TETHER_BOND_NM,
                    bond_k: float = TRNA_TETHER_BOND_K) -> int:
    """Tether the nascent C-terminus to the P-site tRNA, O'Brien-style.

    Replaces the generic position restraint with the peptidyl-tRNA linkage:

    - a harmonic **bond** ``CA(L) -- tRNA`` (the frozen tRNA bead holds the
      C-terminus at the PTC), and
    - an **orienting angle** ``CA(L-1) -- CA(L) -- tRNA`` (double-Gaussian basins
      91.7 deg / 130 deg) that aims the chain down the tunnel.

    ``bead`` is the P-site tRNA attachment bead -- the ribose ``R`` in the topo RNA
    representation, the single ``P`` bead in the cosmo one. ``prev_index`` is the
    CA(L-1) particle index, or None for L == 1 (the angle is then skipped). Returns
    the tRNA bead's system index.
    """
    system = cgModel.system
    anchor_idx = bead_system_index(ribo, L_nascent, segid, resid, bead)
    if anchor_idx is None:
        raise ValueError(f"tRNA tether: {segid} resid {resid} bead {bead!r} not found "
                         f"in the ribosome.")

    bond = mm.HarmonicBondForce()
    bond.addBond(int(cterm_index), anchor_idx, bond_length_nm, bond_k)
    system.addForce(bond)

    if prev_index is not None:
        angle = mm.CustomAngleForce(_TETHER_ANGLE_ENERGY)
        for name, value in _TETHER_ANGLE_PARAMS.items():
            angle.addGlobalParameter(name, value)
        angle.addAngle(int(prev_index), int(cterm_index), anchor_idx)
        system.addForce(angle)
    return anchor_idx


def add_tunnel_wall(system, nascent_indices, x0_nm: float = TUNNEL_WALL_X0_NM,
                    k: float = TUNNEL_WALL_K) -> mm.Force:
    """Add O'Brien's one-sided planar tunnel wall on the nascent chain.

    ``U = k * min(x - x0, 0)^2`` per nascent bead -- penalizes ``x < x0`` only, so
    the chain is kept at ``x >= x0`` and can only extrude forward (+x). ``x0_nm`` is
    the plane position (the C-terminal-AA addition plane / PTC) and ``k`` the force
    constant (kJ/mol/nm^2).
    """
    force = mm.CustomExternalForce("k*r^2; r=min(x-x0, 0)")
    force.addGlobalParameter("k", k)
    force.addGlobalParameter("x0", x0_nm)
    for i in nascent_indices:
        force.addParticle(int(i), [])
    system.addForce(force)
    return force
