"""Rigid ribosome scenery + cross-interactions for co-translational synthesis.

The nascent chain is simulated with cosmo's sequence-based IDP force field (the
**hps_kr** Ashbaugh-Hatch potential for the CSP path); the truncated ribosome is added
to the System as **rigid (mass-0) scenery**. This mirrors the sibling ``topo`` project's
``topo/csp/ribosome.py`` -- same excluded-volume form, tRNA tether and tunnel wall --
adapting only the one thing that is genuinely different in cosmo: the *nascent*
short-range interaction is the HPS Ashbaugh-Hatch / mpipi Wang-Frenkel pair potential
(not topo's Gō contact table).

Following O'Brien et al. (2011/2012):

1. **Append** the ribosome beads at indices ``L..N-1`` with **mass = 0** (frozen; not
   integrated), coordinates as-is. The P-/A-anchors are now real beads.
2. **Nascent short-range force** (Ashbaugh-Hatch for HPS / Wang-Frenkel for mpipi):
   give the ribosome beads a dummy ``addParticle`` entry and restrict the force to the
   interaction group ``{nascent}x{nascent}`` -- so the ribosome never participates in
   the IDP pair potential. For ``mpipi`` the dummy is an in-range ``id = 0``. This is
   cosmo's analog of topo restricting its ``L x L`` Gō contact table to
   ``{nascent}x{nascent}``.
3. **Ribosome-NC excluded volume:** a separate ``CustomNonbondedForce`` reproducing
   O'Brien's NC<->ribosome interaction -- the **12-10-6** form
   ``eps[13(R/r)^12 - 18(R/r)^10 + 4(R/r)^6]`` (``eps = 0.000132 kcal/mol``) with the
   **sum** combination rule ``R_ij = Rmin/2_i + Rmin/2_j`` (O'Brien's convention),
   cutoff 2.0 nm / switch 1.8 nm, on ``{nascent}x{ribosome}`` only.
4. **Electrostatics:** extend the existing Yukawa force with the ribosome charges
   (rRNA phosphate -1e, charged residues) restricted to ``{nascent}x{nascent}`` +
   ``{nascent}x{ribosome}`` (no intra-ribosome electrostatics).

The ribosome is held rigid, so **no intra-ribosome forces are ever computed**; the
nascent chain uses flexible bonds and the ribosome has none, so no constraint ever
involves a mass-0 particle.

**Rmin/2 (collision-radius) source.** The per-bead ``Rmin_2`` values are inherited
verbatim from topo (O'Brien's structure-based CG collision radii) and now live in
:mod:`cosmo.parameters.model_parameters` under the **hps_kr** model (per-AA + the rRNA
``P``/``R``/``BR`` beads). This module reads them exactly as topo does -- so the
ribosome<->nascent 12-10-6 excluded volume is identical to topo's, while the nascent
IDP<->IDP interaction stays the hps_kr Ashbaugh-Hatch potential.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit
from openmm.app.element import carbon as _CARBON

from cosmo.parameters.model_parameters import parameters as MODEL_PARAMS

# Ribosome-NC excluded-volume well depth: 0.000132 kcal/mol -> kJ/mol (O'Brien).
RIBO_NC_EPS_KJ = 0.000132 * 4.184

# All force constants below are in OpenMM units (kJ/mol/nm^2); kcal/mol/A^2 ->
# kJ/mol/nm^2 is 4.184 * 100 = 418.4.
# O'Brien peptidyl-tRNA tether (from continuous_synthesis_v6.py):
#   bond   C-term(CA) -- tRNA:76 R      k = 200 kcal/mol/A^2 = 83680 kJ/mol/nm^2
TRNA_TETHER_BOND_NM = 0.476
TRNA_TETHER_BOND_K = 83680.0            # kJ/mol/nm^2 (= 200 kcal/mol/A^2)
TRNA_TETHER_ANGLE_K = 25.0 * 4.184      # kJ/mol/rad^2 (= 104.6)
TRNA_TETHER_IMPROPER_K = 25.0 * 4.184   # kJ/mol/rad^2
# Per-site resting geometry of the aminoacyl-/peptidyl-tRNA tether (O'Brien), tuple
#   (bond N--R nm, angle N-R-P deg, angle N-R-BR2 deg, improper N-R-P-BR2 deg).
# Bead names: R == tRNA:76 R, P == O'Brien's P (R-1), BR2 == O'Brien's PU2 (R+2).
_TRNA_SITE_GEOM = {
    "AtR": (0.427, 106.0, 127.0, 128.0),    # A-site (aminoacyl-tRNA): stages 1-2
    "PtR": (0.476, 117.0, 130.0, -161.0),   # P-site (peptidyl-tRNA):  stage 3 / prev-AA
}
# O'Brien improper-dihedral energy form: periodic harmonic on |theta - theta0|.
_IMPROPER_ENERGY = ("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); "
                    "pi = 3.1415926535")

# O'Brien's 12-10-6 excluded-volume form (U = eps[13(R/r)^12 - 18(R/r)^10 + 4(R/r)^6];
# well minimum -eps at r = R), with the SUM combining rule R = Rmin/2_i + Rmin/2_j.
_NC_126_ENERGY = ("eps*(13*(R/r)^12 - 18*(R/r)^10 + 4*(R/r)^6); R = rm1 + rm2")

# Planar tunnel wall (O'Brien): a one-sided restraint keeping the nascent chain at
# x >= x0, so it can only extrude forward (+x). x0 is the C-terminal-AA addition plane
# (PTC / P-site) ~ P-anchor x + tether bond length; the runner sets it from the anchors.
TUNNEL_WALL_X0_NM = 1.05
TUNNEL_WALL_K = 8368.0                  # kJ/mol/nm^2 (= 20 kcal/mol/A^2)

# The model carrying the O'Brien Rmin/2 table (inherited from topo). The CSP nascent
# chain runs on hps_kr; a model without Rmin_2 for a bead falls back to this one for the
# *excluded-volume radius only* (steric size is force-field-independent).
_RMIN2_MODEL = "hps_kr"

# O'Brien per-AA sidechain Rmin/2 (nm), read from the hps_kr table (== topo's values).
# Used for the nascent side of the NC<->ribosome excluded volume when no explicit
# per-residue array is supplied.
_AA20 = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
OBRIEN_SC_RMIN_2_NM = {aa: MODEL_PARAMS[_RMIN2_MODEL][aa]["Rmin_2"]
                       for aa in _AA20 if "Rmin_2" in MODEL_PARAMS[_RMIN2_MODEL].get(aa, {})}
OBRIEN_SC_RMIN_2_NM.update({t: OBRIEN_SC_RMIN_2_NM["HIS"] for t in ("HSD", "HSE", "HSP")
                            if "HIS" in OBRIEN_SC_RMIN_2_NM})


@dataclass
class Ribosome:
    """Parsed rigid CG ribosome: per-bead coordinates (nm) and force parameters."""
    coords_nm: np.ndarray       # (M, 3)
    Rmin_2_nm: List[float]      # per-bead collision radius Rmin/2 (nm; sum-rule EV)
    charges: List[float]
    names: List[str]
    resnames: List[str]
    resids: List[int]
    segids: List[str]

    @property
    def n(self) -> int:
        return len(self.Rmin_2_nm)


def _bead_type(name: str, resname: str) -> str:
    """Parameter-lookup key for a CG bead.

    Protein Ca beads (atom name ``CA``) look up by residue name; rRNA beads by atom
    name with trailing digits stripped (``P``, ``R``, ``BR1``/``BR2`` -> ``BR``).
    """
    if name == "CA":
        return resname
    return name.rstrip("0123456789")


def _rmin2_charge(model: str, btype: str) -> Tuple[float, float]:
    """Return ``(Rmin/2 nm, charge)`` for a bead type from ``model_parameters``.

    Reads ``Rmin_2`` + ``charge`` from ``model`` if present, else from the
    :data:`_RMIN2_MODEL` (``hps_kr``) O'Brien table -- so any nascent force field can
    still supply the excluded-volume radius (which is force-field-independent), while
    only hps_kr needs to carry ``Rmin_2``.
    """
    e = MODEL_PARAMS[model].get(btype)
    if isinstance(e, dict) and "Rmin_2" in e:
        return float(e["Rmin_2"]), float(e.get("charge", 0.0))
    fe = MODEL_PARAMS[_RMIN2_MODEL].get(btype)
    if isinstance(fe, dict) and "Rmin_2" in fe:
        charge = float(e["charge"]) if isinstance(e, dict) and "charge" in e else float(fe.get("charge", 0.0))
        return float(fe["Rmin_2"]), charge
    raise ValueError(f"bead type {btype!r} has no Rmin_2 in model {model!r} or the "
                     f"{_RMIN2_MODEL!r} O'Brien table (needed for the ribosome-NC "
                     f"12-10-6 excluded volume).")


def load_ribosome(pdb_file: str, model: str = "hps_kr") -> Ribosome:
    """Parse a (truncated) CG ribosome PDB into a :class:`Ribosome`.

    Reads each ATOM/HETATM record, derives its bead type (:func:`_bead_type`) and looks
    up its Rmin/2 (collision radius, nm) and charge from ``model_parameters``
    (:func:`_rmin2_charge`) -- **protein Ca beads** by residue name, **rRNA P/R/BR
    beads** by bead type. These are O'Brien's values (inherited from topo, carried on
    the hps_kr model). Coordinates are converted from angstrom to nm.

    The ribosome CG PDB is the O'Brien 3/4-bead representation (protein Ca; rRNA
    P/R/BR), produced by :mod:`cosmo.csp.cg_ribosome` + :mod:`cosmo.csp.truncate_ribosome`.

    Raises
    ------
    ValueError
        If a parsed bead type has no Rmin_2, or no ATOM records are found.
    """
    coords, radii, charges = [], [], []
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
            btype = _bead_type(name, resname)
            rm, q = _rmin2_charge(model, btype)
            coords.append((x / 10.0, y / 10.0, z / 10.0))
            radii.append(rm); charges.append(q)
            names.append(name); resnames.append(resname)
            resids.append(resid); segids.append(seg)
    if not coords:
        raise ValueError(f"no ATOM records parsed from ribosome file {pdb_file!r}.")
    return Ribosome(np.asarray(coords, dtype=float), radii, charges,
                    names, resnames, resids, segids)


def append_ribosome(nascent_model, ribo: Ribosome,
                    nascent_rmin_2: Optional[np.ndarray] = None
                    ) -> Tuple[List[int], List[int]]:
    """Append the rigid ribosome to a built nascent cosmo model (system + topology).

    Mutates ``nascent_model`` in place. Must be called **after**
    :func:`cosmo.models.buildCoarseGrainModel` (the nascent forces must exist).
    Appends mass-0 ribosome particles; restricts the nascent short-range (HPS/mpipi)
    and Yukawa forces to their nascent interaction groups; adds the O'Brien 12-10-6
    ribosome-NC excluded volume (sum combining rule, ``{nascent}x{ribosome}`` only);
    and extends the topology with one ribosome chain per segID.

    Returns ``(nascent_indices, ribosome_indices)``.

    Parameters
    ----------
    nascent_model
        The built nascent model whose ``.system`` and ``.topology`` are mutated.
    ribo : Ribosome
        The parsed rigid ribosome (coords, per-bead Rmin/2, charges).
    nascent_rmin_2 : np.ndarray, optional
        Per-nascent-bead Rmin/2 (nm) for the NC-ribosome excluded volume (length must
        equal the nascent atom count). ``None`` -> per-AA O'Brien sidechain Rmin/2
        (:data:`OBRIEN_SC_RMIN_2_NM`), keyed by each nascent residue's name.
    """
    system = nascent_model.system
    topology = nascent_model.topology
    L = nascent_model.n_atoms
    M = ribo.n
    N = L + M
    nascent_idx = list(range(L))
    ribo_idx = list(range(L, N))

    # 1. Rigid (mass-0) ribosome particles.
    for _ in range(M):
        system.addParticle(0.0)

    # 2. Nascent short-range force (cosmo analog of topo's Gō-contact restriction):
    #    dummy ribosome entries + restrict to {nascent}x{nascent} so the ribosome never
    #    participates in the HPS/mpipi pair potential.
    if nascent_model.ashbaugh_HatchForce is not None:
        sr = nascent_model.ashbaugh_HatchForce
        for _ in range(M):
            sr.addParticle((0.1, 0.0))      # dummy (sigma, hps); never evaluated
    elif nascent_model.wang_Frenkel_Force is not None:
        sr = nascent_model.wang_Frenkel_Force
        for _ in range(M):
            sr.addParticle((0,))            # dummy in-range id=0 (mpipi WF is id-indexed)
    else:
        raise RuntimeError("nascent model has no Ashbaugh-Hatch or Wang-Frenkel force.")
    sr.addInteractionGroup(nascent_idx, nascent_idx)

    # 3. Yukawa electrostatics: add ribosome charges; restrict to nascent-nascent +
    #    nascent-ribosome (no intra-ribosome electrostatics).
    yf = nascent_model.yukawaForce
    for q in ribo.charges:
        yf.addParticle((q,))
    yf.addInteractionGroup(nascent_idx, nascent_idx)
    yf.addInteractionGroup(nascent_idx, ribo_idx)

    # 4. Ribosome-NC excluded volume: O'Brien 12-10-6, sum rule R = rm_i + rm_j
    #    (nascent x ribosome only). Nascent Rmin/2 = the given per-residue array, else
    #    per-AA O'Brien sidechain values; ribosome Rmin/2 from the Ribosome.
    nc = mm.CustomNonbondedForce(_NC_126_ENERGY)
    nc.addGlobalParameter("eps", RIBO_NC_EPS_KJ)
    nc.addPerParticleParameter("rm")        # per-bead Rmin/2 (nm); pair R = rm1 + rm2
    nascent_atoms = list(topology.atoms())[:L]   # nascent CA beads (ribosome not appended yet)
    if nascent_rmin_2 is not None:
        if len(nascent_rmin_2) != L:
            raise ValueError(f"nascent_rmin_2 has {len(nascent_rmin_2)} entries but L={L}.")
        for rm in nascent_rmin_2:
            nc.addParticle((float(rm),))
    else:
        for atom in nascent_atoms:
            rn = atom.residue.name
            if rn not in OBRIEN_SC_RMIN_2_NM:
                raise ValueError(f"NC-ribosome EV: residue {rn!r} has no O'Brien Rmin/2 "
                                 f"(OBRIEN_SC_RMIN_2_NM); cannot set its excluded-volume radius.")
            nc.addParticle((OBRIEN_SC_RMIN_2_NM[rn],))
    for rm in ribo.Rmin_2_nm:
        nc.addParticle((float(rm),))
    nc.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    nc.setUseSwitchingFunction(True)
    nc.setSwitchingDistance(1.8 * unit.nanometer)
    nc.setCutoffDistance(2.0 * unit.nanometer)
    nc.addInteractionGroup(nascent_idx, ribo_idx)
    # OpenMM (CPU platform) requires every CustomNonbondedForce to share an identical
    # exclusion list. The nascent short-range force carries the bonded (1-2, 1-3)
    # exclusions; copy them (nascent-nascent pairs, irrelevant to this force's
    # {nascent}x{ribosome} group, but the lists must match).
    for k in range(sr.getNumExclusions()):
        i, j = sr.getExclusionParticles(k)
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


def anchor_coord(ribo: Ribosome, segid: str, resid: int = 76,
                 bead: str = "R") -> np.ndarray:
    """Coordinate (nm) of a named ribosome bead from a loaded :class:`Ribosome`.

    The :class:`Ribosome`-object analog of :func:`cosmo.csp.core.read_anchor` (which
    parses a PDB): picks the P-/A-anchors (``segid='PtR'/'AtR'``, resid 76, ``R`` bead)
    directly from the loaded arrays. Raises if the bead is absent or non-unique.
    """
    matches = [i for i in range(ribo.n)
               if ribo.segids[i] == segid and ribo.resids[i] == resid
               and ribo.names[i] == bead]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one bead (segid={segid!r}, resid={resid}, "
                         f"name={bead!r}) in the ribosome, found {len(matches)}.")
    return ribo.coords_nm[matches[0]]


def bead_system_index(ribo: Ribosome, L_nascent: int, segid: str,
                      resid: int, bead: str = "R") -> Optional[int]:
    """System index of a named ribosome bead (e.g. the P-site tRNA R anchor).

    Ribosome beads are appended after the ``L_nascent`` nascent particles in load
    order, so the system index is ``L_nascent + (position in the arrays)``. Returns
    None if no such bead exists (e.g. a pyrimidine has no ``BR2``).
    """
    for i in range(ribo.n):
        if ribo.segids[i] == segid and ribo.resids[i] == resid and ribo.names[i] == bead:
            return L_nascent + i
    return None


def add_trna_tether(nascent_model, cterm_index: int, prev_index,
                    ribo: Ribosome, L_nascent: int,
                    segid: str = "PtR", resid: int = 76) -> None:
    """Tether a nascent residue to a tRNA, the full O'Brien orienting way.

    Reproduces the aminoacyl-/peptidyl-tRNA linkage of ``continuous_synthesis_v6.py``
    for the resting geometry of the given site (``segid``: ``"AtR"`` A-site /
    ``"PtR"`` P-site). For the restrained residue N (= ``cterm_index``):

    - **bond** ``N -- tRNA:R`` (length from :data:`_TRNA_SITE_GEOM`, k
      :data:`TRNA_TETHER_BOND_K`);
    - **orienting harmonic angles** ``N -- R -- P`` and ``N -- R -- BR2``
      (:data:`TRNA_TETHER_ANGLE_K`) -- fix the residue's bearing in the tRNA frame;
    - **improper** ``N -- R -- P -- BR2`` (periodic-harmonic on ``abs(theta-theta0)``,
      :data:`TRNA_TETHER_IMPROPER_K`) -- fixes the out-of-plane sense;
    - **backbone orienting angle** ``prev -- N -- R`` (added to the model's
      double-Gaussian angle force) -- aims the chain down the tunnel. Only present when
      the nascent model has a Gaussian angle force (``hps_ss``); skipped otherwise.

    ``prev_index`` is the CA(N-1) particle index, or None (skips the backbone angle).
    ``P`` / ``BR2`` are resolved by bead name; a site missing either skips that
    angle/improper.

    Raises
    ------
    ValueError
        If ``segid`` is not a known tRNA site, or its ``R`` bead is not found.
    """
    if segid not in _TRNA_SITE_GEOM:
        raise ValueError(f"tRNA tether: unknown site {segid!r} (expected one of "
                         f"{sorted(_TRNA_SITE_GEOM)}).")
    system = nascent_model.system
    R_idx = bead_system_index(ribo, L_nascent, segid, resid, "R")
    if R_idx is None:
        raise ValueError(f"tRNA tether: {segid} resid {resid} R bead not found.")
    P_idx = bead_system_index(ribo, L_nascent, segid, resid, "P")
    U2_idx = bead_system_index(ribo, L_nascent, segid, resid, "BR2")
    bond_nm, ang_P_deg, ang_U2_deg, imp_deg = _TRNA_SITE_GEOM[segid]
    d2r = math.radians

    # 1. bond N -- R (holds the residue at its site's resting length).
    bond = mm.HarmonicBondForce()
    bond.addBond(int(cterm_index), R_idx, bond_nm, TRNA_TETHER_BOND_K)
    system.addForce(bond)

    # 2. orienting harmonic angles N -- R -- P and N -- R -- BR2 (tRNA frame).
    haf = mm.HarmonicAngleForce()
    if P_idx is not None:
        haf.addAngle(int(cterm_index), R_idx, P_idx, d2r(ang_P_deg), TRNA_TETHER_ANGLE_K)
    if U2_idx is not None:
        haf.addAngle(int(cterm_index), R_idx, U2_idx, d2r(ang_U2_deg), TRNA_TETHER_ANGLE_K)
    if haf.getNumAngles() > 0:
        system.addForce(haf)

    # 3. improper N -- R -- P -- BR2 (out-of-plane sense; O'Brien periodic-harmonic form).
    if P_idx is not None and U2_idx is not None:
        imf = mm.CustomTorsionForce(_IMPROPER_ENERGY)
        imf.addPerTorsionParameter("k")
        imf.addPerTorsionParameter("theta0")
        imf.addTorsion(int(cterm_index), R_idx, P_idx, U2_idx,
                       [TRNA_TETHER_IMPROPER_K, d2r(imp_deg)])
        system.addForce(imf)

    # 4. backbone orienting angle prev -- N -- R (aims the chain down the tunnel).
    #    Only if the nascent model has a Gaussian angle force (hps_ss); else skipped.
    if prev_index is not None and getattr(nascent_model, "gaussianAngleForce", None) is not None:
        nascent_model.gaussianAngleForce.addAngle(int(prev_index), int(cterm_index), R_idx)


def add_tunnel_wall(system, nascent_indices, x0_nm: float = TUNNEL_WALL_X0_NM,
                    k: float = TUNNEL_WALL_K) -> mm.Force:
    """Add O'Brien's one-sided planar tunnel wall on the nascent chain.

    ``U = k * min(x - x0, 0)^2`` per nascent bead -- penalizes ``x < x0`` only, so the
    chain is kept at ``x >= x0`` and can only extrude forward (+x). ``k`` is a **global**
    parameter here; the coexisting C-terminus position restraint
    (:func:`cosmo.csp.core.add_cterm_restraint`) makes *its* ``k`` per-particle so the
    two forces do not clash on the shared global name.
    """
    force = mm.CustomExternalForce("k*r^2; r=min(x-x0, 0)")
    force.addGlobalParameter("k", k)
    force.addGlobalParameter("x0", x0_nm)
    for i in nascent_indices:
        force.addParticle(int(i), [])
    system.addForce(force)
    return force
