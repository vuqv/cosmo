# Import OpenMM library
# Import COSMO library
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import numpy as np
import openmm as mm
from openmm import unit

from parmed.exceptions import OpenMMWarning

import cosmo
from cosmo.utils import write_pdb_with_chain_ids

# Suppress OpenMM warnings
warnings.filterwarnings("ignore", category=OpenMMWarning)


def main():
    """
        Run a simulation using the COSMO library and parameters specified in a config file.

        Usage: python run_simulation.py -f md.ini
        or cosmo-simulation -f md.ini
        """

    # Default values:
    # Set default values
    md_steps = 1000
    dt = 0.01
    nstxout = 10
    nstlog = 10
    nstcomm = 100
    model = 'hps_urry'
    tcoupl = True
    ref_t = 300.0
    pcoupl = False
    ref_p = 1.0
    frequency_p = 25
    pbc = False
    device = 'CPU'
    ppn = 1
    restart = False
    minimize = True

    # Other attributes
    tau_t = None
    box_dimension = None
    protein_code = None
    checkpoint = None
    pdb_file = None
    nascent_chain_id = None
    # Parse config file:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-input', '-f', type=str, help='simulation config file', required=True)

    args = parser.parse_args()
    config_file = args.input

    # Read config file
    print(f"Reading simulation parameters from {config_file} file...")
    print('-' * 70)
    print('Simulation parameters:')
    config = configparser.ConfigParser()
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(config_file)
    params = config['OPTIONS']

    # Update simulation parameters
    md_steps = int(params.get('md_steps', md_steps))
    print(f'  md_steps: {md_steps}')
    dt = float(params.get('dt', dt)) * unit.picoseconds
    print(f'  dt: {dt}')
    nstxout = int(params['nstxout'])
    print(f'  nstxout: {nstxout}')
    nstlog = int(params['nstlog'])
    print(f'  nstlog: {nstlog}')
    nstcomm = int(params.get('nstcomm', nstcomm))
    print(f'  nstcomm: {nstcomm}')
    model = params.get('model', model)
    print(f'  model: {model}')
    nascent_chain_id = params.get('nascent_chain_id', nascent_chain_id)
    print(f'  nascent_chain_id: {nascent_chain_id}')
    tcoupl = bool(strtobool(params.get('tcoupl', tcoupl)))
    if tcoupl:
        ref_t = float(params['ref_t']) * unit.kelvin
        tau_t = float(params['tau_t']) / unit.picoseconds
        print(f'  tcoupl: on (ref_t={ref_t}, tau_t={tau_t})')
    else:
        print('  tcoupl: off')
    pbc = bool(strtobool(params.get('pbc', pbc)))
    if pbc:
        box_dimension = loads(params['box_dimension'])
        print(f'  pbc: on (box_dimension={box_dimension} nm)')
    else:
        box_dimension = None
        print('  pbc: off')
    pcoupl = bool(strtobool(params.get('pcoupl', pcoupl)))
    if pcoupl:
        assert pbc, "Pressure coupling requires box dimensions and periodic boundary condition is on"
        ref_p = float(params['ref_p']) * unit.bar
        frequency_p = int(params['frequency_p'])
        print(f'  pcoupl: on (ref_p={ref_p}, frequency={frequency_p})')
    else:
        print('  pcoupl: off')
    pdb_file = params['pdb_file']
    print(f'  pdb_file: {pdb_file}')
    protein_code = params['protein_code']
    print(f'  output_prefix: {protein_code}')
    checkpoint = params.get('checkpoint', checkpoint)
    device = params.get('device', device)
    print(f'  device: {device}')
    if device == "CPU":
        ppn = int(params.get('ppn', ppn))
        print(f'  threads: {ppn}')
    restart = bool(strtobool(params.get('restart', restart)))
    print(f'  restart: {restart}')
    if restart:
        minimize = False
    else:
        minimize = bool(strtobool(params.get('minimize', minimize)))
        print(f'  minimize: {minimize}')
    print('-' * 70)

    """
    End of reading parameters
    """
    #get the number of atoms in the nascent chain
    # build a dummy model to get the nascent atom indices, lightweight and fast
    dummy_model = cosmo.system(pdb_file, model)
    dummy_model.coarseGrainingStructure()
    dummy_model.getAtoms()
    nascent_atom_indices = [atom.index for atom in dummy_model.atoms if atom.residue.chain.id == nascent_chain_id]
    print(f"There are {len(nascent_atom_indices)} atoms in the nascent chain")
    print(f"Nascent atom indices: {nascent_atom_indices}")
    for idx, atom_id in enumerate(nascent_atom_indices, start=1):
        print(f"Atom {idx} has index {atom_id} and chain id {dummy_model.atoms[atom_id].residue.chain.id}")

    print('-' * 70)

    # Particle index to indicate it is synthesizing, then restraint at the origin. This also tells buildHPSModel to ignore interaction of untranslated atoms.
    restraint_atom_index = nascent_atom_indices[29]  # Topology atom index

    # build the HPS model
    hps_model = cosmo.models.buildHPSModel(pdb_file, minimize=minimize, model=model, box_dimension=box_dimension)

    # # Remove center of mass motion
    # hps_model.system.addForce(mm.CMMotionRemover(nstcomm))
    

    # Find particle index in System (particles are added in same order as atoms after getAtoms)
    particle_index = None
    for i, atom in enumerate(hps_model.atoms):
        if atom.index == restraint_atom_index:
            particle_index = i
            break
    
    if particle_index is not None:
        # Create harmonic position restraint
        restraint_force = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint_force.addGlobalParameter('k', 10000.0 * unit.kilojoule_per_mole / unit.nanometer**2)  # Force constant
        restraint_force.addGlobalParameter('x0', 0.0 * unit.nanometer)
        restraint_force.addGlobalParameter('y0', 0.4 * unit.nanometer)
        restraint_force.addGlobalParameter('z0', 0.3 * unit.nanometer)
        restraint_force.addParticle(particle_index, [])
        hps_model.system.addForce(restraint_force)
        print(f"Added position restraint to atom {restraint_atom_index} (particle {particle_index}) at (0,0,0)")
    else:
        print(f"Warning: Atom index {restraint_atom_index} not found in system")

    # Add planar wall force to prevent translated residues from moving into x < 0
    # Translated residues are those in nascent chain with index <= restraint_atom_index
    translated_atom_indices = [atom_idx for atom_idx in nascent_atom_indices if atom_idx < restraint_atom_index]
    translated_particle_indices = []
    for i, atom in enumerate(hps_model.atoms):
        if atom.index in translated_atom_indices:
            translated_particle_indices.append(i)
    
    if translated_particle_indices:
        print(f"Adding planar wall force: {len(translated_particle_indices)} translated residues restricted to x >= 0")
        wall_force = mm.CustomExternalForce('k_wall * step(-x) * x^2')
        wall_force.addGlobalParameter('k_wall', 10000.0 * unit.kilojoule_per_mole / unit.nanometer**2)
        for particle_idx in translated_particle_indices:
            wall_force.addParticle(particle_idx, [])
        hps_model.system.addForce(wall_force)
        print(f"Added planar wall force to particles: {translated_particle_indices}")

    # dump Forcefield File
    if model in ['hps_kr', 'synthesis_kr', 'hps_urry', 'hps_ss']:
        """ current dumpForceFieldData function can only write the standard format of forcefield which require
         sigma, epsilon for each residue.
        """
        hps_model.dumpForceFieldData(f'{protein_code}_ff.dat')
    # dump Structure into PDB file for visualize
    hps_model.dumpStructure(f'{protein_code}_init.pdb')
    hps_model.dumpTopology(f'{protein_code}.psf')

    if device == 'GPU':
        # Run simulation on CUDA
        print("Running on GPU (CUDA)")
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
        # in case of many GPUs present, we can select which one to use

    elif device == 'CPU':
        print(f"Running on CPU with {ppn} cores")
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {'Threads': str(ppn)}

    print('Simulation started')
    start_time = time.time()

    integrator = mm.LangevinIntegrator(ref_t, tau_t, dt)
    simulation = mm.app.Simulation(hps_model.topology, hps_model.system, integrator, platform,
                                   properties)
    if restart:
        simulation.loadCheckpoint(checkpoint)
        print(f"Restarting from step: {simulation.context.getState().getStepCount()}")
        nsteps_remain = md_steps - simulation.context.getState().getStepCount()
    else:
        # xyz = np.array(hps_model.positions / unit.nanometer)
        # xyz[:, 0] -= np.amin(xyz[:, 0])
        # xyz[:, 1] -= np.amin(xyz[:, 1])
        # xyz[:, 2] -= np.amin(xyz[:, 2])
        # hps_model.positions = xyz * unit.nanometer
        simulation.context.setPositions(hps_model.positions)
        simulation.context.setVelocitiesToTemperature(ref_t)
        nsteps_remain = md_steps

    simulation.reporters = []
    simulation.reporters.append(mm.app.CheckpointReporter(checkpoint, nstxout))
    simulation.reporters.append(
        mm.app.DCDReporter(f'{protein_code}.dcd', nstxout, enforcePeriodicBox=bool(pbc),
                           append=restart))
    simulation.reporters.append(
        mm.app.StateDataReporter(f'{protein_code}.log', nstlog, step=True, time=True,
                                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                 temperature=True, remainingTime=True, speed=True,
                                 totalSteps=md_steps, separator='\t', append=restart))
    simulation.step(nsteps_remain)

    # write the last frame (using custom writer to preserve chain IDs)
    last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(pbc)).getPositions()
    write_pdb_with_chain_ids(hps_model.topology, last_frame, f'{protein_code}_final.pdb')
    simulation.saveCheckpoint(checkpoint)

    print(f"Finished in {(time.time() - start_time):.3f} seconds")


if __name__ == '__main__':
    print(f"OpenMM version: {mm.__version__}")
    main()
