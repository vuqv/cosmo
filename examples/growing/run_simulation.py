#!/usr/bin/env python
"""Co-translational ("growing" nascent chain) simulation example.

A specialized variant of the canonical runner. It keeps the bespoke physics of
this example -- a harmonic position restraint on the residue currently in the
ribosome exit tunnel and a planar wall that keeps already-translated residues at
x >= 0 -- while reusing the shared machinery for everything else:

    * config parsing      -> cosmo.read_simulation_config
    * model build + dumps  -> cosmo.engine.build_model
    * simulation setup     -> cosmo.engine.setup_simulation (shift_positions=False;
                              the absolute frame matters here because the restraint
                              is anchored to a fixed point)
    * reporters / finalize -> cosmo.engine.attach_reporters / finalize_simulation
                              (final frame written with chain IDs preserved)

The only extra control-file key beyond the standard set is ``nascent_chain_id``.

Usage:

    python run_simulation.py -f md.ini
"""
import argparse
import configparser
import time
import warnings

import openmm as mm
from openmm import unit

try:
    from parmed.exceptions import OpenMMWarning

    warnings.filterwarnings("ignore", category=OpenMMWarning)
except Exception:
    pass

import cosmo
from cosmo import engine
from cosmo.utils import write_pdb_with_chain_ids


def main():
    parser = argparse.ArgumentParser(description="Co-translational nascent-chain "
                                                 "growing simulation.")
    parser.add_argument('-input', '-f', type=str, required=True,
                        help='simulation config file')
    args = parser.parse_args()

    print(f"OpenMM version: {mm.__version__}")

    # Standard parsing for all the usual settings...
    cfg = cosmo.read_simulation_config(args.input)

    # ...plus this example's extra key (the nascent chain's chain ID).
    extra = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    extra.read(args.input)
    nascent_chain_id = extra['OPTIONS'].get('nascent_chain_id')
    print(f"  nascent_chain_id: {nascent_chain_id}")
    print('-' * 70)

    # Identify the nascent-chain atoms with a lightweight throwaway model.
    dummy_model = cosmo.system(cfg.pdb_file, cfg.model)
    dummy_model.coarseGrainingStructure()
    dummy_model.getAtoms()
    nascent_atom_indices = [atom.index for atom in dummy_model.atoms
                            if atom.residue.chain.id == nascent_chain_id]
    print(f"There are {len(nascent_atom_indices)} atoms in the nascent chain")

    # The residue currently being synthesized: restrained at a fixed point near
    # the exit tunnel. This index also tells the model to ignore interactions of
    # not-yet-translated atoms.
    restraint_atom_index = nascent_atom_indices[29]  # topology atom index
    print('-' * 70)

    # Build the model and write provenance files.
    built = engine.build_model(cfg)

    # --- bespoke forces -------------------------------------------------------
    # Map the topology atom index to its particle index in the System.
    particle_index = None
    for i, atom in enumerate(built.atoms):
        if atom.index == restraint_atom_index:
            particle_index = i
            break

    if particle_index is not None:
        restraint_force = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint_force.addGlobalParameter(
            'k', 10000.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
        restraint_force.addGlobalParameter('x0', 0.0 * unit.nanometer)
        restraint_force.addGlobalParameter('y0', 0.4 * unit.nanometer)
        restraint_force.addGlobalParameter('z0', 0.3 * unit.nanometer)
        restraint_force.addParticle(particle_index, [])
        built.system.addForce(restraint_force)
        print(f"Added position restraint to atom {restraint_atom_index} "
              f"(particle {particle_index})")
    else:
        print(f"Warning: atom index {restraint_atom_index} not found in system")

    # Planar wall keeping already-translated residues (index < restraint atom) at x >= 0.
    translated_atom_indices = [a for a in nascent_atom_indices if a < restraint_atom_index]
    translated_particle_indices = [i for i, atom in enumerate(built.atoms)
                                   if atom.index in translated_atom_indices]
    if translated_particle_indices:
        print(f"Adding planar wall force: {len(translated_particle_indices)} "
              f"translated residues restricted to x >= 0")
        wall_force = mm.CustomExternalForce('k_wall * step(-x) * x^2')
        wall_force.addGlobalParameter(
            'k_wall', 10000.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
        for particle_idx in translated_particle_indices:
            wall_force.addParticle(particle_idx, [])
        built.system.addForce(wall_force)
    # -------------------------------------------------------------------------

    # No COM remover here: the restraint/wall anchor the system to a fixed frame,
    # so positions must NOT be shifted (shift_positions=False).
    print('Simulation started')
    start_time = time.time()
    ctx = engine.setup_simulation(cfg, built, control_file=args.input,
                                  shift_positions=False)
    engine.attach_reporters(cfg, ctx.simulation, append=ctx.restart_active,
                            total_steps=cfg.md_steps)
    ctx.simulation.step(ctx.nsteps_remain)

    # Preserve chain IDs in the final structure (multi-chain ribosome + nascent).
    engine.finalize_simulation(cfg, ctx, built.topology, start_time,
                               write_pdb=write_pdb_with_chain_ids)


if __name__ == '__main__':
    main()
