#!/usr/bin/env python
"""Tutorial 10 -- co-translational synthesis by **iterative growth**.

Port of ``tutorials/growing_reference/`` to the current cosmo codebase. The
nascent chain is grown one residue at a time inside a **frozen explicit
ribosome**: for each residue, a new CA bead is appended to the complex at the
origin, the whole complex is rebuilt, the newest residue is held at the PTC
``(0.38, 0, 0)`` nm and earlier residues are walled forward (``x > 0.38 nm``), and
a short MD run relaxes/extrudes the chain. The final structure seeds the next
residue.

Differences from the legacy reference (all behaviour-preserving):

- ``buildHPSModel`` -> :func:`cosmo.models.buildCoarseGrainModel` (the alias was
  removed; the ``frozen_indices`` / ``except_chains`` arguments still exist).
- ribosome--ribosome nonbonded is excluded with **interaction groups** (post-build,
  as in :func:`cosmo.translation.ribosome.append_ribosome`) instead of the O(n^2)
  ``nb_exclusions`` pair list.
- per-iteration setup/reporters/finalize reuse :mod:`cosmo.engine` (standard
  ``cosmoReporter`` log + ``_runinfo.log`` provenance).
- the **trajectory writes the nascent chain only** (nascent-only PSF + DCD); the
  simulation and each iteration's final PDB still cover the whole complex.
- outputs use the ``<outdir>/L_<L>/traj.*`` layout, so the shipped
  ``cosmo-elongate-movie`` stitches them directly.

This test uses a **constant** ``md_steps`` per residue. Realistic per-codon
translation kinetics (variable steps from ``ctf_utils.codon_translation_time``) is
a later addition -- see ``PLAN.md`` Phase 5.

Run from this folder::

    python grow_nascent.py -f md.ini
    cosmo-elongate-movie -o synth_out      # then: vmd -e synth_out/movie.tcl
"""
from __future__ import annotations

import argparse
import configparser
import os
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import openmm as mm
from openmm import unit

import cosmo
from cosmo import engine
from cosmo.core import models
from cosmo.utils import parse_nascent_sequence, write_pdb_with_chain_ids
from cosmo.utils.config import strtobool

try:                                            # parmed emits a noisy OpenMM warning
    from parmed.exceptions import OpenMMWarning
    warnings.filterwarnings("ignore", category=OpenMMWarning)
except Exception:
    pass

# PTC geometry (kept from the reference). The newest residue is restrained at
# x = PTC_X_NM on the axis; earlier residues are kept at x > PTC_X_NM by a wall.
PTC_X_NM = 0.38
PTC_K = 10000.0          # kJ/mol/nm^2, restraint stiffness
WALL_K = 10000.0         # kJ/mol/nm^2, planar-wall stiffness


# --------------------------------------------------------------------------
# Nascent-only trajectory output (subset PSF + DCD by atom index)
# --------------------------------------------------------------------------
# The complex is [ribosome ... | nascent ...] in atom order, so the nascent beads
# are a (trailing) subset. OpenMM 8.2's DCDReporter has no atomSubset option, so we
# write the nascent subset ourselves -- a CHARMM PSF + a matching DCD over the same
# index list, both in ascending index order, so they pair for MDAnalysis / VMD.
def _particle_property_by_atom_index(atoms, values):
    if isinstance(values, list):
        return {atom.index: values[i] for i, atom in enumerate(atoms)}
    return {atom.index: values for atom in atoms}


def write_subset_psf(cosmo_model, atom_indices, output_file):
    """Write a PSF containing only ``atom_indices`` (in the given order)."""
    atom_indices = list(atom_indices)
    atoms_by_index = {atom.index: atom for atom in cosmo_model.atoms}
    selected_atoms = [atoms_by_index[i] for i in atom_indices]
    index_map = {atom.index: i + 1 for i, atom in enumerate(selected_atoms)}
    masses = _particle_property_by_atom_index(cosmo_model.atoms, cosmo_model.particles_mass)
    charges = _particle_property_by_atom_index(cosmo_model.atoms, cosmo_model.particles_charge)

    bonds = []
    for atom1, atom2 in cosmo_model.topology.bonds():
        if atom1.index in index_map and atom2.index in index_map:
            bonds.append((index_map[atom1.index], index_map[atom2.index]))

    with open(output_file, 'w') as f:
        f.write("PSF CHEQ EXT XPLOR\n\n")
        f.write("         1 !NTITLE\n\n\n")
        f.write(f"{len(selected_atoms):10d} !NATOM\n")
        for psf_index, atom in enumerate(selected_atoms, start=1):
            residue = atom.residue
            segid = residue.chain.id or 'A'
            resid = residue.id or str(psf_index)
            f.write(
                f"{psf_index:10d} {segid:<8s} {resid:<8s} {residue.name:<8s} "
                f"{atom.name:<8s} {atom.name:<8s} {float(charges[atom.index]):14.6f} "
                f"{float(masses[atom.index]):13.4f}           \n"
            )
        f.write(f"\n{len(bonds):10d} !NBOND: bonds\n")
        bond_numbers = [n for bond in bonds for n in bond]
        for i in range(0, len(bond_numbers), 8):
            f.write(''.join(f"{n:10d}" for n in bond_numbers[i:i + 8]) + "\n")
        f.write("\n         0 !NTHETA: angles\n\n\n")
        f.write("         0 !NPHI: dihedrals\n\n\n")
        f.write("         0 !NIMPHI: impropers\n\n\n")
        f.write("         0 !NDON: donors\n\n\n")
        f.write("         0 !NACC: acceptors\n\n\n")
        f.write("         0 !NNB\n\n")
        for i in range(0, len(selected_atoms), 8):
            f.write(''.join(f"{0:10d}" for _ in selected_atoms[i:i + 8]) + "\n")
        f.write("\n         0         0 !NGRP NST2\n\n")


def create_subset_topology(source_topology, atom_indices):
    """An OpenMM Topology with only ``atom_indices`` and their internal bonds."""
    selected = set(atom_indices)
    sub = mm.app.Topology()
    sub.setPeriodicBoxVectors(source_topology.getPeriodicBoxVectors())
    atom_map, chain_map, residue_map = {}, {}, {}
    for chain in source_topology.chains():
        for residue in chain.residues():
            for atom in residue.atoms():
                if atom.index not in selected:
                    continue
                if chain not in chain_map:
                    chain_map[chain] = sub.addChain(id=chain.id)
                if residue not in residue_map:
                    residue_map[residue] = sub.addResidue(
                        residue.name, chain_map[chain], id=residue.id,
                        insertionCode=residue.insertionCode)
                atom_map[atom] = sub.addAtom(atom.name, atom.element, residue_map[residue],
                                             id=atom.id)
    for atom1, atom2 in source_topology.bonds():
        if atom1 in atom_map and atom2 in atom_map:
            sub.addBond(atom_map[atom1], atom_map[atom2])
    return sub


class SubsetDCDReporter:
    """A DCD reporter that writes only ``atom_indices`` each frame (OpenMM 8.2)."""

    def __init__(self, file, report_interval, source_topology, atom_indices,
                 append=False, enforcePeriodicBox=None):
        self._reportInterval = report_interval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        self._atomIndices = list(atom_indices)
        self._topology = create_subset_topology(source_topology, self._atomIndices)
        self._out = open(file, 'r+b' if append else 'wb')
        self._dcd = None

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {'steps': steps, 'periodic': self._enforcePeriodicBox,
                'include': ['positions']}

    def report(self, simulation, state):
        if self._dcd is None:
            self._dcd = mm.app.DCDFile(self._out, self._topology,
                                       simulation.integrator.getStepSize(),
                                       self._reportInterval, self._reportInterval,
                                       self._append)
        positions = state.getPositions(asNumpy=True)
        subset = [positions[i] for i in self._atomIndices]
        self._dcd.writeModel(subset, periodicBoxVectors=state.getPeriodicBoxVectors())

    def __del__(self):
        try:
            self._out.close()
        except Exception:
            pass


# --------------------------------------------------------------------------
# Append one CA bead to the complex (kept from the reference)
# --------------------------------------------------------------------------
def append_amino_acid_to_complex(input_pdb, output_pdb, res_name, res_num,
                                 chain_id='8', position=(0.0, 0.0, 0.0)):
    """Append a new CA atom for one amino acid to an existing PDB file.

    The new bead is written at ``position`` (angstrom) on chain ``chain_id``; the
    ribosome / previously-grown beads are copied through unchanged.
    """
    header_lines, atom_lines, max_serial = [], [], 0
    with open(input_pdb, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    max_serial = max(max_serial, int(line[6:11].strip()))
                except ValueError:
                    pass
                atom_lines.append(line)
            elif line.startswith(('TER', 'END')):
                continue
            else:
                header_lines.append(line)

    with open(output_pdb, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        for line in atom_lines:
            f.write(line + '\n')
        new_serial = max_serial + 1
        x, y, z = position
        f.write(
            f"ATOM  {new_serial:5d}  CA  {res_name:>3s} {chain_id:1s}{res_num:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n")
        f.write(f"TER   {new_serial+1:5d}      {res_name:>3s} {chain_id:1s}{res_num:4d}\n")
        f.write("END\n")


# --------------------------------------------------------------------------
# Ribosome scenery: frozen beads + interaction groups
# --------------------------------------------------------------------------
def extract_ribosome(rib_structure: str, model: str) -> Tuple[List[int], List[str]]:
    """Return ``(atom_indices, chain_ids)`` of the CG ribosome (frozen scenery)."""
    rib = cosmo.system(rib_structure, model)
    rib.coarseGrainingStructure()
    rib.getAtoms()
    atom_indices = [atom.index for atom in rib.atoms]
    chain_ids = [chain.id for chain in rib.topology.chains()]
    return atom_indices, chain_ids


def add_ribosome_interaction_groups(cgModel, nascent_idx, ribo_idx) -> None:
    """Restrict the nonbonded forces to nascent-nascent + nascent-ribosome.

    Excludes ribosome-ribosome nonbonded (never computed; the ribosome is frozen)
    via interaction groups -- the efficient equivalent of the reference's O(n^2)
    intra-ribosome ``nb_exclusions`` (mirrors ``append_ribosome``).
    """
    for force in (cgModel.ashbaugh_HatchForce, cgModel.wang_Frenkel_Force,
                  cgModel.yukawaForce):
        if force is None:
            continue
        force.addInteractionGroup(nascent_idx, nascent_idx)
        force.addInteractionGroup(nascent_idx, ribo_idx)


def add_ptc_restraint_and_wall(system, nascent_indices, restrain_newest=True) -> None:
    """Hold the newest residue at the PTC; wall the rest forward (+x).

    - ``restrain_newest`` (growth/stallation): harmonically restrain the newest
      nascent bead to ``(PTC_X_NM, 0, 0)`` and wall every *earlier* bead at
      ``x > PTC_X_NM``.
    - ``restrain_newest=False`` (ejection): no PTC restraint -- the finished chain
      is released; wall **all** nascent beads at ``x > PTC_X_NM`` so they can only
      diffuse forward out of the ribosome, not back through the closed PTC end.
    """
    if restrain_newest:
        newest = nascent_indices[-1]
        restraint = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint.addGlobalParameter('k', PTC_K)
        restraint.addGlobalParameter('x0', PTC_X_NM)
        restraint.addGlobalParameter('y0', 0.0)
        restraint.addGlobalParameter('z0', 0.0)
        restraint.addParticle(newest, [])
        system.addForce(restraint)
        walled = nascent_indices[:-1]
    else:
        walled = nascent_indices

    if walled:
        wall = mm.CustomExternalForce(
            f'k_wall * step({PTC_X_NM} - x) * ({PTC_X_NM} - x)^2')
        wall.addGlobalParameter('k_wall', WALL_K)
        for idx in walled:
            wall.addParticle(idx, [])
        system.addForce(wall)


# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
@dataclass
class GrowConfig:
    """Parsed contents of the growth control file (``md.ini``)."""
    model: str = "hps_kr"
    rib_structure: str = "rib.pdb"
    nascent_structure: str = "nascent.pdb"
    nascent_chain_id: str = "8"
    md_steps: int = 1000
    organism: str = "ecoli"
    translation_type: str = "fast"
    # Post-elongation phase at full length (0 = skip). 'ejection' releases the PTC
    # restraint and lets the finished chain diffuse out of the ribosome (-> the
    # 'ejection/' folder); 'stallation' keeps it pinned at the PTC.
    post_elongation: str = "ejection"
    post_elongation_steps: int = 0
    dt_ps: float = 0.01
    ref_t: float = 310.0
    tau_t: float = 0.01
    nstxout: int = 100
    nstlog: int = 100
    pbc: bool = False
    box_dimension: Optional[float] = None
    outdir: str = "synth_out"
    device: str = "GPU"
    ppn: int = 4


def read_grow_config(config_file: str) -> GrowConfig:
    """Parse ``md.ini`` ([OPTIONS]) into a :class:`GrowConfig`."""
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read config file: {config_file!r}")
    o = cp["OPTIONS"]
    pbc = bool(strtobool(o.get("pbc", "no")))
    box = o.get("box_dimension", "").strip()
    return GrowConfig(
        model=o.get("model", "hps_kr"),
        rib_structure=o.get("rib_structure", "rib.pdb"),
        nascent_structure=o.get("nascent_structure", "nascent.pdb"),
        nascent_chain_id=o.get("nascent_chain_id", "8"),
        md_steps=int(str(o.get("md_steps", "1000")).replace("_", "")),
        organism=o.get("organism", "ecoli"),
        translation_type=o.get("translation_type", "fast"),
        post_elongation=o.get("post_elongation", "ejection"),
        post_elongation_steps=int(str(o.get("post_elongation_steps", "0")).replace("_", "")),
        dt_ps=float(o.get("dt", 0.01)),
        ref_t=float(o.get("ref_t", 310.0)),
        tau_t=float(o.get("tau_t", 0.01)),
        nstxout=int(o.get("nstxout", 100)),
        nstlog=int(o.get("nstlog", 100)),
        pbc=pbc,
        box_dimension=(float(box) if (pbc and box) else None),
        outdir=o.get("outdir", "synth_out"),
        device=o.get("device", "GPU"),
        ppn=int(o.get("ppn", 4)),
    )


def _make_cfg(out_dir: Path, complex_pdb: str, gc: GrowConfig, md_steps: int):
    """Build a per-iteration :class:`cosmo.SimulationConfig` for the engine."""
    cfg = cosmo.SimulationConfig()
    cfg.model = gc.model
    cfg.md_steps = md_steps
    cfg.dt = gc.dt_ps * unit.picoseconds
    cfg.nstxout = gc.nstxout
    cfg.nstlog = gc.nstlog
    cfg.nstchk = gc.nstxout
    cfg.nstcomm = None             # no COM removal (frozen ribosome + restraint set the frame)
    cfg.tcoupl = True
    cfg.ref_t = gc.ref_t * unit.kelvin
    cfg.tau_t = gc.tau_t / unit.picoseconds
    cfg.pbc = gc.pbc
    cfg.box_dimension = gc.box_dimension
    cfg.pdb_file = complex_pdb
    cfg.init_position = None       # use the built complex coordinates as-is
    cfg.output_dir = str(out_dir)
    cfg.outname = "traj"
    cfg.device = gc.device
    cfg.ppn = gc.ppn
    cfg.restart = False
    cfg.minimize = False
    return cfg


# --------------------------------------------------------------------------
# Single iteration
# --------------------------------------------------------------------------
def run_single_iteration(complex_pdb: str, L: int, gc: GrowConfig,
                         rib_indices: List[int], rib_chains: List[str],
                         out_root: Path, md_steps: int, restrain: bool = True,
                         out_subdir: Optional[str] = None,
                         label: Optional[str] = None) -> str:
    """Build, restrain, run one length-``L`` complex; return the full final PDB.

    The System is the full complex (frozen ribosome + nascent chain). Outputs are
    written nascent-only (PSF + DCD); the ``_final.pdb`` is the full complex (it
    seeds the next residue's append).

    - ``restrain`` (growth/stallation): pin the newest residue at the PTC.
      ``restrain=False`` (ejection): release it and wall all nascent beads forward.
    - ``out_subdir`` overrides the ``L_<L>`` output folder (e.g. ``ejection``).
    - ``label`` overrides the console banner.
    """
    out_dir = out_root / (out_subdir or f"L_{L:03d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'#'*66}\n# {label or f'Nascent length L = {L}'}\n{'#'*66}")

    # 1. build the full complex with the ribosome frozen (mass-0) and excluded from
    #    bonding. No nb_exclusions -- ribosome-ribosome nonbonded is handled by the
    #    interaction groups below.
    cgModel = models.buildCoarseGrainModel(
        complex_pdb, model=gc.model, minimize=False,
        box_dimension=gc.box_dimension, frozen_indices=rib_indices,
        except_chains=rib_chains, check_forces=False)

    rib_set = set(rib_indices)
    nascent_indices = [a.index for a in cgModel.atoms if a.index not in rib_set]
    if not nascent_indices:
        raise ValueError("no nascent atoms found (all atoms are in the frozen ribosome).")
    print(f"  complex atoms={cgModel.n_atoms}  ribosome={len(rib_indices)} (frozen)  "
          f"nascent={len(nascent_indices)}  model={gc.model}")

    # 2. exclude ribosome-ribosome nonbonded via interaction groups (not O(n^2)).
    add_ribosome_interaction_groups(cgModel, nascent_indices, rib_indices)

    # 3. PTC restraint on the newest residue + forward wall on the earlier ones
    #    (ejection: restrain=False -> released chain, all beads walled forward).
    add_ptc_restraint_and_wall(cgModel.system, nascent_indices, restrain_newest=restrain)

    # 4. nascent-only PSF (matches the nascent-only DCD written below).
    cfg = _make_cfg(out_dir, complex_pdb, gc, md_steps)
    write_subset_psf(cgModel, nascent_indices, cfg.output_path('.psf'))

    # 5. run via the engine; swap in the nascent-only DCD reporter; finalize with
    #    the full-complex final PDB (preserving chain IDs for the next append).
    start = time.time()
    ctx = engine.setup_simulation(cfg, cgModel, control_file=None,
                                  shift_positions=False)
    engine.attach_reporters(cfg, ctx.simulation, total_steps=cfg.md_steps)
    ctx.simulation.reporters[1] = SubsetDCDReporter(
        cfg.output_path('.dcd'), cfg.nstxout, cgModel.topology, nascent_indices,
        enforcePeriodicBox=bool(cfg.pbc))
    print(f"  running {md_steps} steps ({gc.device}) ...")
    ctx.simulation.step(md_steps)
    engine.finalize_simulation(cfg, ctx, cgModel.topology, start,
                               write_pdb=write_pdb_with_chain_ids)

    final_pdb = cfg.output_path('_final.pdb')
    print(f"  -> {final_pdb}  ({time.time() - start:.1f} s)")
    return final_pdb


# --------------------------------------------------------------------------
# Growth loop / CLI
# --------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Iterative co-translational growth of a nascent chain in a "
                    "frozen ribosome (tutorial 10). Constant steps/residue; codon "
                    "kinetics is a later addition. python grow_nascent.py -f md.ini")
    parser.add_argument('-f', '--config', required=True, help='control file (md.ini)')
    parser.add_argument('--rib', default=None, help='override ribosome PDB')
    parser.add_argument('--nascent', default=None, help='override nascent PDB')
    parser.add_argument('--start', type=int, default=1, help='first residue (1-indexed)')
    parser.add_argument('--end', type=int, default=None, help='last residue (default: all)')
    parser.add_argument('--steps', type=int, default=None,
                        help='override md_steps per residue (constant schedule)')
    parser.add_argument('-o', '--outdir', default=None, help='override output dir')
    args = parser.parse_args(argv)

    print(f"OpenMM version: {mm.__version__}")
    gc = read_grow_config(args.config)
    if args.rib:
        gc.rib_structure = args.rib
    if args.nascent:
        gc.nascent_structure = args.nascent
    if args.outdir:
        gc.outdir = args.outdir
    md_steps = args.steps if args.steps is not None else gc.md_steps

    print("=" * 66)
    print("Iterative nascent-chain growth (tutorial 10)")
    print("=" * 66)
    print(f"  ribosome (frozen): {gc.rib_structure}")
    print(f"  nascent (sequence): {gc.nascent_structure}  chain id={gc.nascent_chain_id}")
    print(f"  model={gc.model}  organism={gc.organism}  translation_type={gc.translation_type}")
    print(f"  schedule: CONSTANT {md_steps} steps/residue (codon kinetics: Phase 5)")
    print(f"  device={gc.device}  outdir={gc.outdir}")

    # Sequence to grow (CA-only nascent PDB -> [(res_name, res_num), ...]).
    sequence = parse_nascent_sequence(gc.nascent_structure)
    print(f"  nascent sequence: {len(sequence)} residues")

    # Ribosome scenery (frozen beads + chains to exclude from bonding).
    rib_indices, rib_chains = extract_ribosome(gc.rib_structure, gc.model)
    print(f"  ribosome: {len(rib_indices)} beads, chains {rib_chains}")
    print("=" * 66)

    out_root = Path(gc.outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    start_idx = args.start - 1
    end_idx = min(args.end if args.end is not None else len(sequence), len(sequence))

    current_complex = gc.rib_structure
    for i in range(start_idx, end_idx):
        res_name, res_num = sequence[i]
        L = i + 1
        out_dir = out_root / f"L_{L:03d}"
        out_dir.mkdir(parents=True, exist_ok=True)
        new_complex = str(out_dir / "complex.pdb")
        append_amino_acid_to_complex(current_complex, new_complex, res_name, res_num,
                                     chain_id=gc.nascent_chain_id,
                                     position=(0.0, 0.0, 0.0))
        current_complex = run_single_iteration(
            new_complex, L, gc, rib_indices, rib_chains, out_root, md_steps)

    print(f"\n{'='*66}\nDone. Grew residues {args.start}..{end_idx}. "
          f"Per-length outputs under {out_root}/\n{'='*66}")

    # Post-elongation phase: once the chain is full length, either release the PTC
    # restraint and let it diffuse out of the ribosome (ejection) or keep it pinned
    # (stallation). Continues from the final full-length complex; the ribosome stays
    # frozen and the forward wall stays on (only +x escape).
    if gc.post_elongation_steps > 0:
        phase = gc.post_elongation.strip().lower()
        if phase not in ("ejection", "stallation"):
            raise ValueError(f"post_elongation must be 'ejection' or 'stallation', "
                             f"got {gc.post_elongation!r}.")
        restrain = phase == "stallation"
        L = end_idx
        print(f"\n=== Post-elongation: {phase} (L = {L}, {gc.post_elongation_steps} "
              f"steps, PTC restraint {'ON' if restrain else 'OFF -> diffusion'}) "
              f"-> {out_root / phase}/ ===")
        run_single_iteration(
            current_complex, L, gc, rib_indices, rib_chains, out_root,
            gc.post_elongation_steps, restrain=restrain, out_subdir=phase,
            label=f"Post-elongation: {phase} (L = {L})")
        print(f"Done. {phase.capitalize()} written to {out_root / phase}/")


if __name__ == '__main__':
    main()
