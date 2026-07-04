#!/usr/bin/env python
"""
Iteratively grow nascent chain by adding one amino acid at a time.

This script:
1. Reads nascent.pdb to get the amino acid sequence
2. Starts with rib.pdb
3. For each amino acid:
   - Appends the new amino acid to the current complex
   - Creates complex_l{N}.pdb
   - Runs simulation for 1000 steps
   - Saves outputs (complex_l{N}_final.pdb, complex_l{N}.dcd, etc.)
"""

import os
import sys
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import openmm as mm
from openmm import unit

from parmed.exceptions import OpenMMWarning

import cosmo
from cosmo.utils import write_pdb_with_chain_ids, parse_nascent_sequence, generate_mRNA_sequence
from cosmo.utils.ctf_utils import codon_translation_time

# Suppress OpenMM warnings
warnings.filterwarnings("ignore", category=OpenMMWarning)


def append_amino_acid_to_complex(input_pdb, output_pdb, res_name, res_num, chain_id='8', position=(0.0, 0.0, 0.0)):
    """
    Append a new CA atom for an amino acid to an existing PDB file.
    
    Parameters
    ----------
    input_pdb : str
        Path to input PDB file (ribosome or previous complex)
    output_pdb : str
        Path to output PDB file
    res_name : str
        Three-letter residue name (e.g., 'MET', 'ALA')
    res_num : int
        Residue number
    chain_id : str
        Chain ID (default '8')
    position : tuple
        (x, y, z) position in Angstroms (default (0, 0, 0))
    """
    # Read input PDB and find the last atom serial number
    header_lines = []  # REMARK, CRYST1, etc.
    atom_lines = []     # ATOM and HETATM lines
    max_serial = 0
    
    with open(input_pdb, 'r') as f:
        for line in f:
            line_stripped = line.rstrip('\n')
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract serial number (columns 6-11)
                try:
                    serial = int(line[6:11].strip())
                    max_serial = max(max_serial, serial)
                except ValueError:
                    pass
                atom_lines.append(line_stripped)
            elif line.startswith('TER'):
                # Skip TER records - we'll add new ones at the end
                continue
            elif line.startswith('END'):
                # Skip END - we'll add it at the end
                continue
            else:
                # Keep header lines (REMARK, CRYST1, etc.)
                header_lines.append(line_stripped)
    
    # Write output PDB
    with open(output_pdb, 'w') as f:
        # Write header lines
        for line in header_lines:
            f.write(line + '\n')
        
        # Write all existing atom lines
        for line in atom_lines:
            f.write(line + '\n')
        
        # Append new CA atom
        new_serial = max_serial + 1
        x, y, z = position
        atom_line = (
            f"ATOM  {new_serial:5d}  CA  {res_name:>3s} {chain_id:1s}{res_num:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
        )
        f.write(atom_line)
        
        # Write TER and END
        f.write(f"TER   {new_serial+1:5d}      {res_name:>3s} {chain_id:1s}{res_num:4d}\n")
        f.write("END\n")


def run_single_iteration(complex_pdb, output_prefix, frozen_indices, except_chains, nb_exclusions, sim_params, md_steps=1000):
    """
    Run a single iteration of the simulation.

    Parameters
    ----------
    complex_pdb : str
        Path to input complex PDB file
    output_prefix : str
        Prefix for output files (e.g., 'complex_l1')
    frozen_indices : list
        List of atom indices to freeze
    except_chains : list
        List of chain IDs to exclude from bonding
    nb_exclusions : list
        List of pairs of atom indices to exclude from nonbonded forces
    sim_params : dict
        Simulation parameters: model, box_dimension, device, ppn, traj_dir
    md_steps : int
        Number of MD steps
    """
    model = sim_params['model']
    dt = sim_params['dt']
    ref_t = sim_params['ref_t']
    tau_t = sim_params['tau_t']
    pbc = sim_params['pbc']
    nstxout = sim_params['nstxout']
    nstlog = sim_params['nstlog']
    box_dimension = sim_params['box_dimension']
    device = sim_params['device']
    ppn = sim_params['ppn']
    traj_dir = sim_params['traj_dir']

    print(f"\n{'='*70}")
    print(f"Running simulation for {output_prefix}")
    print(f"{'='*70}")


    # Build HPS model
    # All nascent atoms are translated and interact normally
    # Only ribosome-ribosome interactions are excluded
    hps_model = cosmo.models.buildHPSModel(
        complex_pdb,
        minimize=False,  # Don't minimize for iterative growth
        model=model,
        box_dimension=box_dimension,
        except_chains=except_chains,
        frozen_indices=frozen_indices,
        nb_exclusions=nb_exclusions,
    )

    #convert ribosome atom indices to set for faster lookup
    rib_atom_indices_set = set(frozen_indices)
    # Get nascent atom indices directly from the built model
    nascent_atom_indices = [atom.index for atom in hps_model.atoms
                           if atom.index not in rib_atom_indices_set]  # if not in ribosome index then it is nascent

    if not nascent_atom_indices:
        raise ValueError("No nascent atoms found (all atoms are in the frozen ribosome set)")
    
    # The last atom index is the newest amino acid to be restrained
    # In iterative growth, all nascent atoms are "translated" - no untranslated regions
    growing_atom_index = nascent_atom_indices[-1]
    print(f"Restraining growing amino acid at atom index {growing_atom_index}")    
    if growing_atom_index is not None:
        # Create harmonic position restraint at origin
        restraint_force = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint_force.addGlobalParameter('k', 10000.0 * unit.kilojoule_per_mole / unit.nanometer**2)
        restraint_force.addGlobalParameter('x0', 0.38 * unit.nanometer)
        restraint_force.addGlobalParameter('y0', 0.0 * unit.nanometer)
        restraint_force.addGlobalParameter('z0', 0.0 * unit.nanometer)
        restraint_force.addParticle(growing_atom_index, [])
        hps_model.system.addForce(restraint_force)
        print(f"Added position restraint to atom {growing_atom_index} at (0.38, 0.0, 0.0)")
    
    # Add planar wall force for all translated residues except the newest one (which is restrained)
    # In iterative growth, all nascent atoms are "translated", so we apply wall force to all except the newest
    if len(nascent_atom_indices) > 1:
        # Wall at x = 0.0 Å (0.0 nm) - translated atoms must have x > 0.0 Å
        wall_position = 0.38  # in nanometers (3.8 Angstroms)
        print(f"Adding planar wall force: {len(nascent_atom_indices)-1} nascent residues (all except newest) restricted to x > {wall_position*10:.1f} Å")
        # step(wall_position - x) is 1 when x < wall_position, 0 when x >= wall_position
        # Penalty is (wall_position - x)^2 when x < wall_position
        wall_force = mm.CustomExternalForce(f'k_wall * step({wall_position} - x) * ({wall_position} - x)^2')
        wall_force.addGlobalParameter('k_wall', 10000.0 * unit.kilojoule_per_mole / unit.nanometer**2)
        for atom_idx in nascent_atom_indices[:-1]:
            wall_force.addParticle(atom_idx, [])
        hps_model.system.addForce(wall_force)
    
    # Create trajectory directory if it doesn't exist
    os.makedirs(traj_dir, exist_ok=True)
    
    # Write topology PSF file to traj folder
    psf_file = os.path.join(traj_dir, f'{output_prefix}.psf')
    hps_model.dumpTopology(psf_file)
    print(f"Topology written to {psf_file}")
    
    # Setup platform
    if device == 'GPU':
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
    else:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {'Threads': str(ppn)}
    
    # Create simulation
    integrator = mm.LangevinIntegrator(ref_t, tau_t, dt)
    simulation = mm.app.Simulation(hps_model.topology, hps_model.system, integrator, platform, properties)
    
    # Set initial positions and velocities
    simulation.context.setPositions(hps_model.positions)
    simulation.context.setVelocitiesToTemperature(ref_t)
    
    # Setup reporters - all files go to traj folder
    checkpoint_file = os.path.join(traj_dir, f'{output_prefix}.chk')
    simulation.reporters.append(mm.app.CheckpointReporter(checkpoint_file, nstxout))
    # DCD file goes to traj folder
    dcd_file = os.path.join(traj_dir, f'{output_prefix}.dcd')
    simulation.reporters.append(
        mm.app.DCDReporter(dcd_file, nstxout, enforcePeriodicBox=bool(pbc), append=False)
    )
    # Log file goes to traj folder
    log_file = os.path.join(traj_dir, f'{output_prefix}.log')
    simulation.reporters.append(
        mm.app.StateDataReporter(
            log_file, nstlog, step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, remainingTime=True, speed=True,
            totalSteps=md_steps, separator='\t', append=False
        )
    )
    
    # Run simulation
    print(f"Running {md_steps} steps...")
    start_time = time.time()
    simulation.step(md_steps)
    elapsed = time.time() - start_time
    print(f"Simulation completed in {elapsed:.3f} seconds")
    
    # Write final structure to traj folder
    last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(pbc)).getPositions()
    final_pdb = os.path.join(traj_dir, f'{output_prefix}_final.pdb')
    write_pdb_with_chain_ids(hps_model.topology, last_frame, final_pdb)
    print(f"Final structure written to {final_pdb}")
    
    return final_pdb


def main():
    parser = argparse.ArgumentParser(description='Iteratively grow nascent chain')
    parser.add_argument('-f', '--config', type=str, required=True, help='Config file (md.ini)')
    parser.add_argument('--rib', type=str, default='rib.pdb', help='Ribosome PDB file')
    parser.add_argument('--nascent', type=str, default='nascent.pdb', help='Nascent chain PDB file')
    parser.add_argument('--start', type=int, default=1, help='Start from amino acid number (1-indexed, default: 1)')
    parser.add_argument('--end', type=int, default=None, help='End at amino acid number (default: all)')
    parser.add_argument('--steps', type=int, default=1000, help='MD steps per iteration (default: 1000)')
    
    args = parser.parse_args()
    
    # Read config file
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(args.config)
    params = config['OPTIONS']
    
    model = params.get('model', 'hps_kr')

    rib_structure = params.get('rib_structure', 'rib.pdb')
    nascent_structure = params.get('nascent_structure', 'nascent.pdb')
    scale_factor = float(params.get('scale_factor', 1.0))
    dt = float(params.get('dt', 0.01)) * unit.picoseconds  # timestep in picoseconds
    ref_t = float(params.get('ref_t', 310.0)) * unit.kelvin
    tau_t = float(params.get('tau_t', 0.01)) / unit.picoseconds
    nstxout = int(params.get('nstxout', 100))
    nstlog = int(params.get('nstlog', 100))

    pbc = bool(strtobool(params.get('pbc', 'no')))
    box_dimension = loads(params['box_dimension']) if pbc else None  # only read when pbc=yes

    pcoupl = bool(strtobool(params.get('pcoupl', 'no')))
    if pcoupl:
        ref_p = float(params.get('ref_p', 1.0)) * unit.bar
        frequency_p = int(params.get('frequency_p', 25))
        assert box_dimension is not None, "Pressure coupling requires box dimensions"
        assert pbc, "Pressure coupling requires periodic boundary condition"
        print(f"Pressure coupling (pcoupl): on (ref_p={ref_p}, frequency={frequency_p})")
    else:
        print('Pressure coupling (pcoupl): off')

    nascent_chain_id = params.get('nascent_chain_id', '8')
    organism = params.get('organism', 'ecoli')
    translation_type = params.get('translation_type', 'fast')
    device = params.get('device', 'CPU')
    ppn = int(params.get('ppn', 1))
    
    print("="*70)
    print("Iterative Nascent Chain Growth Simulation")
    print("="*70)
    print(f"Ribosome PDB: {rib_structure}")
    print(f"Nascent PDB: {nascent_structure}")
    print(f"Organism: {organism}")
    print(f"Model: {model}")
    print(f"Translation type: {translation_type}")
    print(f"Scale factor: {scale_factor}")
    print(f"Timestep (dt): {dt} ps")
    print(f"Reference temperature (ref_t): {ref_t} K")
    print(f"Time constant for temperature coupling (tau_t): {tau_t} ps")
    print(f"Nascent chain ID: {nascent_chain_id}")
    print(f"Device: {device}")
    print(f"Steps per iteration: variable (translation_time / dt / scale_factor)")
    print("="*70)
    
    # Parse nascent sequence
    print(f"\nParsing nascent sequence from {nascent_structure}...")
    sequence = parse_nascent_sequence(nascent_structure)
    aa_three = [res_name for res_name, _ in sequence]
    mRNA_sequence, translation_time_ms = generate_mRNA_sequence(aa_three, translation_type, organism=organism)
    print(f"Found {len(sequence)} amino acids in nascent chain")
    print(f"mRNA sequence: {mRNA_sequence}")
    print(f"Estimated translation time: {translation_time_ms:.1f} ms (real time)")
    print(f"Estimated translation time in simulation: {(10**3)*translation_time_ms/scale_factor:.6f} microseconds")
    print(f"Translation type: {translation_type}")
    print(f"Organism: {organism}")
    print("="*70)
    
    
    # Determine range
    start_idx = args.start - 1  # Convert to 0-indexed
    end_idx = args.end if args.end is not None else len(sequence)
    end_idx = min(end_idx, len(sequence))
    
    print(f"\nWill grow from amino acid {args.start} to {end_idx}")
    print("="*70)
    
    # Create trajectory directory and package simulation parameters
    traj_dir = 'traj'
    os.makedirs(traj_dir, exist_ok=True)
    print(f"Output files (PDB, DCD, PSF) will be written to '{traj_dir}/' directory")

    sim_params = {
        'model': model,
        'dt': dt,
        'ref_t': ref_t,
        'tau_t': tau_t,
        'pbc': pbc,
        'nstxout': nstxout,
        'nstlog': nstlog,
        'box_dimension': box_dimension,
        'device': device,
        'ppn': ppn,
        'traj_dir': traj_dir,
    }

    # Start with ribosome
    current_complex = rib_structure

    # get atom indices of ribosome
    rib_model = cosmo.system(rib_structure, model)
    rib_model.coarseGrainingStructure()
    rib_model.getAtoms()
    # Get all atom indices of ribosome to freeze during simulation
    rib_atom_indices = [atom.index for atom in rib_model.atoms]
    print(f"There are {len(rib_atom_indices)} atoms in the ribosome")
    # print(f"Ribosome atom indices: {rib_atom_indices}")
    print('-' * 70)
    # get all ribosomal chain to exclude from bonding (except_chains)
    except_chains = [chain.id for chain in rib_model.topology.chains()]
    print(f"Ribosome chains to exclude from bonding: {except_chains}")
    print('-' * 70)
    # Construct nb_exclusions for nonbonded forces (all pairs of ribosome atoms)
    n_rib = len(rib_atom_indices)
    nb_exclusions = [(rib_atom_indices[i], rib_atom_indices[j])
                     for i in range(n_rib) for j in range(i + 1, n_rib)]
    print(f"Nonbonded exclusions: {len(nb_exclusions)}")
    print('-' * 70)

    # Precompute (nsteps, codon, codon_time_ms) per residue: nsteps = time_ps / (dt_ps * scale_factor)
    # 1 ms = 1e9 ps, so time_ps = codon_time_ms * 1e9
    dt_ps = dt.value_in_unit(unit.picosecond)  # scalar in picoseconds
    residue_steps = []
    for k in range(len(sequence)):
        codon = mRNA_sequence[3 * k : 3 * k + 3]
        codon_time_ms = codon_translation_time[organism][codon]['Median_time']
        time_ps = codon_time_ms * 1e9  # ms -> ps
        nsteps = max(1, int(time_ps / (dt_ps * scale_factor)))
        residue_steps.append((nsteps, codon, codon_time_ms))

    # Iterate through amino acids
    for i in range(start_idx, end_idx):
        res_name, res_num = sequence[i]
        iteration = i + 1  # 1-indexed iteration number
        nsteps, codon, codon_time_ms = residue_steps[i]

        print(f"\n{'#'*70}")
        print(f"# Iteration {iteration}: Adding {res_name} {res_num}")
        print(f"{'#'*70}")
        
        # Create new complex PDB in traj folder
        output_prefix = f'complex_l{iteration}'
        new_complex = os.path.join(traj_dir, f'{output_prefix}.pdb')
        
        print(f"Creating {new_complex} by appending {res_name} {res_num} to {current_complex}...")
        append_amino_acid_to_complex(
            current_complex, 
            new_complex, 
            res_name, 
            res_num, 
            chain_id=nascent_chain_id,
            position=(0.0, 0.0, 0.0)  # Position in Angstroms (x = 0.0 nm)
        )

        print(f"Codon {codon} translation time: {codon_time_ms:.1f} ms -> {nsteps} steps")
        # Run simulation
        final_pdb = run_single_iteration(
            new_complex,
            output_prefix,
            rib_atom_indices,
            except_chains,
            nb_exclusions,
            sim_params,
            nsteps,
        )
        
        # Update current complex for next iteration
        current_complex = final_pdb
        print(f"\nIteration {iteration} complete. Final structure: {final_pdb}")
    
    print(f"\n{'='*70}")
    print("All iterations complete!")
    print(f"{'='*70}")


if __name__ == '__main__':
    print(f"OpenMM version: {mm.__version__}")
    main()
