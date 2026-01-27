#!/usr/bin/env python
"""
Setup simulation for ribosome-nascent protein complex.

The ribosome atoms:
- Only interact with nascent protein via nonbonded forces
- Do NOT interact with each other (no bonds, no nonbonded interactions)

The nascent protein has full interactions (bonds, angles, etc.)
"""

import sys
import numpy as np
import openmm as mm
from openmm import unit
from openmm.app import PDBFile, Modeller

import cosmo
from cosmo.core import system


def combine_structures(ribosome_pdb, nascent_pdb, output_pdb):
    """
    Combine ribosome and nascent protein PDB files into one structure.
    
    Parameters
    ----------
    ribosome_pdb : str
        Path to ribosome PDB file
    nascent_pdb : str
        Path to nascent protein PDB file
    output_pdb : str
        Path to output combined PDB file
    """
    # Load both structures
    ribosome = PDBFile(ribosome_pdb)
    nascent = PDBFile(nascent_pdb)
    
    # Combine using Modeller
    modeller = Modeller(ribosome.topology, ribosome.positions)
    modeller.add(nascent.topology, nascent.positions)
    
    # Write combined structure
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    print(f"Combined structures written to: {output_pdb}")
    return modeller.topology, modeller.positions


def setup_ribosome_nascent_system(combined_pdb, nascent_chain_id='8', model='hps_urry', 
                                   box_dimension=None, minimize=True):
    """
    Setup COSMO system for ribosome-nascent complex with modified interactions.
    
    Parameters
    ----------
    combined_pdb : str
        Path to combined PDB file
    nascent_chain_id : str
        Chain ID of nascent protein (default: '8')
    model : str
        COSMO model to use (default: 'hps_urry')
    box_dimension : float or list, optional
        Box dimensions for PBC
    minimize : bool
        Whether to minimize energy (default: True)
    
    Returns
    -------
    cosmo_model : system
        COSMO system object
    ribosome_atom_indices : list
        List of atom indices belonging to ribosome
    nascent_atom_indices : list
        List of atom indices belonging to nascent protein
    """
    # Build standard COSMO model
    print("Building COSMO model...")
    cosmo_model = system(combined_pdb, model)
    
    # Build the system as usual first
    print("Setting up system...")
    cosmo_model.coarseGrainingStructure()
    cosmo_model.getAtoms()
    
    # Identify ribosome vs nascent atoms by chain ID (after coarse graining)
    ribosome_atom_indices = []
    nascent_atom_indices = []
    
    for i, atom in enumerate(cosmo_model.topology.atoms()):
        chain_id = atom.residue.chain.id
        if chain_id == nascent_chain_id:
            nascent_atom_indices.append(i)
        else:
            ribosome_atom_indices.append(i)
    
    print(f"Ribosome atoms: {len(ribosome_atom_indices)}")
    print(f"Nascent protein atoms: {len(nascent_atom_indices)}")
    
    # Get bonds (will be filtered later in bond force)
    cosmo_model.getBonds()
    cosmo_model.setMassPerResidueType()
    cosmo_model.setChargePerResidueType()
    
    if model in ['hps_kr', 'synthesis_kr', 'hps_urry', 'hps_ss']:
        cosmo_model.setRadiusPerResidueType()
        cosmo_model.setHPSPerResidueType()
    
    # Add forces
    cosmo_model.setBondForceConstants()
    cosmo_model.addHarmonicBondForces()
    
    # Modify bond force to remove ribosome-ribosome interactions
    print("Modifying bond forces (disabling ribosome-ribosome bonds)...")
    bond_force = None
    for force in cosmo_model.system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            bond_force = force
            break
    
    if bond_force:
        ribosome_set = set(ribosome_atom_indices)
        disabled_count = 0
        # Set force constant to 0 for ribosome-ribosome bonds
        for i in range(bond_force.getNumBonds()):
            atom1, atom2, length, k = bond_force.getBondParameters(i)
            if atom1 in ribosome_set and atom2 in ribosome_set:
                bond_force.setBondParameters(i, atom1, atom2, length, 0.0 * unit.kilojoule_per_mole / unit.nanometer**2)
                disabled_count += 1
        print(f"Disabled {disabled_count} ribosome-ribosome bonds")
        # Note: updateParametersInContext will be called when context is created
    
    # Add nonbonded forces
    if box_dimension:
        use_pbc = True
        if isinstance(box_dimension, list):
            cosmo_model.topology.setPeriodicBoxVectors(
                ((box_dimension[0], 0, 0), (0, box_dimension[1], 0), (0, 0, box_dimension[2])))
        else:
            cosmo_model.topology.setPeriodicBoxVectors(
                ((box_dimension, 0, 0), (0, box_dimension, 0), (0, 0, box_dimension)))
        unit_cell = cosmo_model.topology.getPeriodicBoxVectors()
        cosmo_model.system.setDefaultPeriodicBoxVectors(*unit_cell)
    else:
        use_pbc = False
    
    cosmo_model.addYukawaForces(use_pbc)
    cosmo_model.addAshbaughHatchForces(use_pbc)
    
    # Modify nonbonded forces to exclude ribosome-ribosome interactions
    print("Modifying nonbonded forces to exclude ribosome-ribosome interactions...")
    
    ribosome_set = set(ribosome_atom_indices)
    exclusion_count = 0
    
    # Get all nonbonded forces and add exclusions
    for force in cosmo_model.system.getForces():
        if isinstance(force, mm.CustomNonbondedForce):
            # Add exclusions for all ribosome-ribosome pairs
            for i in ribosome_atom_indices:
                for j in ribosome_atom_indices:
                    if i < j:  # Only add each pair once
                        try:
                            force.addExclusion(i, j)
                            exclusion_count += 1
                        except:
                            pass  # Exclusion might already exist
        elif isinstance(force, mm.NonbondedForce):
            # For standard NonbondedForce, we need to use a different approach
            # We can create a CustomNonbondedForce wrapper or modify interaction groups
            # For now, we'll add exclusions if possible
            for i in ribosome_atom_indices:
                for j in ribosome_atom_indices:
                    if i < j:
                        try:
                            force.addExclusion(i, j)
                            exclusion_count += 1
                        except:
                            pass
    
    print(f"Added {exclusion_count} exclusions for ribosome-ribosome pairs")
    
    # Optional: Freeze ribosome atoms by setting mass to 0
    # Uncomment the following lines if you want ribosome atoms to be frozen:
    # print("Freezing ribosome atoms (setting mass to 0)...")
    # for idx in ribosome_atom_indices:
    #     cosmo_model.system.setParticleMass(idx, 0.0 * unit.dalton)
    
    # Create system object
    print("Creating system object...")
    cosmo_model.createSystemObject(minimize=minimize, check_bond_distances=True)
    
    return cosmo_model, ribosome_atom_indices, nascent_atom_indices


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python setup_ribosome_nascent.py <ribosome.pdb> <nascent.pdb> [combined.pdb]")
        print("  This script combines the structures and sets up the simulation system.")
        sys.exit(1)
    
    ribosome_pdb = sys.argv[1]
    nascent_pdb = sys.argv[2]
    combined_pdb = sys.argv[3] if len(sys.argv) > 3 else 'ribosome_nascent_combined.pdb'
    
    # Combine structures
    topology, positions = combine_structures(ribosome_pdb, nascent_pdb, combined_pdb)
    
    # Setup system
    cosmo_model, ribosome_indices, nascent_indices = setup_ribosome_nascent_system(
        combined_pdb, nascent_chain_id='8', model='hps_urry', minimize=True
    )
    
    print("\nSetup complete!")
    print(f"Combined structure: {combined_pdb}")
    print(f"Ribosome atoms: {len(ribosome_indices)}")
    print(f"Nascent protein atoms: {len(nascent_indices)}")
    print("\nNote: Ribosome atoms only interact with nascent protein via nonbonded forces.")
    print("Ribosome-ribosome interactions are excluded.")
