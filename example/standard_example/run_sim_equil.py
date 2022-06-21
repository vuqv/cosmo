# Import OpenMM library
# Import sbmOpenMM library
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from sys import stdout

import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *
from parmed.exceptions import OpenMMWarning

import hps

warnings.filterwarnings("ignore", category=OpenMMWarning)

# Parse config file:
parser = argparse.ArgumentParser(description="\n Usage: python run_simulation.py -f config.ini ")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()
# Reading parameters

config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
config.read(args.input)
params = config['DEFAULT']
md_steps = int(params['md_steps'])
equil_steps = int(params['equil_steps'])
dt = float(params['dt']) * picosecond
nstxout = int(params['nstxout'])
nstlog = int(params['nstlog'])

ref_t = float(params['ref_t']) * kelvin
tau_t = float(params['tau_t']) / picosecond
pbc = strtobool(params['pbc'])

protein_code = params['protein_code']
checkpoint_file = protein_code + '.chk'
checkpoint_file_equil = protein_code + '_equil.chk'
pdb_file = params['pdb_file']
device = params['device']
ppn = params['ppn']
gen_vel = strtobool(params['gen_vel'])
restart = strtobool(params['restart'])
"""
End of reading parameters
"""
# @TODO: check if use pbc then initialize hps model with pbc
cgModel = hps.models.getCAModel(pdb_file, hps_scale='kr')

# dump Forcefield File
cgModel.dumpForceFieldData('forcefield.dat')
# dump Structure into PDB file for visualize
cgModel.dumpStructure(f'{protein_code}_init.pdb')
cgModel.dumpTopology(f'{protein_code}.psf')

if device == 'GPU':
    # Run simulation on CUDA
    print(f"Running simulation on GPU CUDA")
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
    # in case of many GPUs present, we can select which one to use


elif device == 'CPU':
    print(f"Running simulation on CPU using {ppn} cores")
    platform = Platform.getPlatformByName('CPU')
    properties = {'Threads': str(ppn)}

print('Simulation started')
start_time = time.time()

if restart:
    # simulation reporter
    integrator = LangevinIntegrator(ref_t, tau_t, dt)
    simulation = Simulation(cgModel.topology, cgModel.system, integrator, platform, properties)
    simulation.loadCheckpoint(checkpoint_file)
    print(
        f"Restart from checkpoint, Time = {simulation.context.getState().getTime()}, Step= {simulation.context.getState().getStepCount()}")
    # number of steps remain to run
    nsteps_remain = md_steps - simulation.context.getState().getStepCount()
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(checkpoint_file, nstxout))
    simulation.reporters.append(DCDReporter(f'{protein_code}.dcd', nstxout, append=True))
    simulation.reporters.append(
        StateDataReporter(stdout, nstlog, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t'))
    simulation.reporters.append(
        StateDataReporter(f'{protein_code}_prod.log', nstlog, step=True, time=True, potentialEnergy=True,
                          kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t', append=True))
    simulation.step(nsteps_remain)
else:
    # Production phase, create new integrator and simulation context to reset number of steps
    integrator = LangevinIntegrator(ref_t, tau_t, dt)
    simulation = Simulation(cgModel.topology, cgModel.system, integrator, platform, properties)
    # Set initial positions: translate input coordinate
    xyz = np.array(cgModel.positions / nanometer)
    xyz[:, 0] -= np.amin(xyz[:, 0])
    xyz[:, 1] -= np.amin(xyz[:, 1])
    xyz[:, 2] -= np.amin(xyz[:, 2])
    cgModel.positions = xyz * nanometer
    simulation.context.setPositions(cgModel.positions)
    simulation.context.setVelocitiesToTemperature(ref_t)
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(checkpoint_file, nstxout))
    simulation.reporters.append(DCDReporter(f'{protein_code}.dcd', nstxout, append=False))
    simulation.reporters.append(
        StateDataReporter(stdout, nstlog, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t'))
    simulation.reporters.append(
        StateDataReporter(f'{protein_code}.log', nstlog, step=True, time=True, potentialEnergy=True,
                          kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t', append=False))
    simulation.step(md_steps)


# write the last frame
lastframe = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(cgModel.topology, lastframe, open(f'{protein_code}_final.pdb', 'w'))
simulation.saveCheckpoint(checkpoint_file)

print("--- Finished in %s seconds ---" % (time.time() - start_time))
