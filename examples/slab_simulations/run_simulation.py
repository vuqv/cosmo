# Import OpenMM library
# Import hpsOpenMM library
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import numpy as np
from openmm import *
from openmm.app import *
from parmed.exceptions import OpenMMWarning

import cosmo

warnings.filterwarnings("ignore", category=OpenMMWarning)

# Parse config file:
parser = argparse.ArgumentParser(description="\n Usage: python run_simulation.py -f md.ini ")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
config.read(args.input)
params = config['DEFAULT']
# Reading parameters
md_steps = int(params['md_steps'])
dt = float(params['dt']) * unit.picoseconds
nstxout = int(params['nstxout'])
nstlog = int(params['nstlog'])
model = params['model']
ref_t = float(params['ref_t']) * unit.kelvin
tau_t = float(params['tau_t']) / unit.picoseconds
ref_p = float(params['ref_p']) * unit.bar
nstpcouple = int(params['nstpcouple'])
pbc = strtobool(params['pbc'])
if pbc:
    box_dimension = loads(params['box_dimension'])
else:
    box_dimension = None

protein_code = params['protein_code']
checkpoint = params['checkpoint']
pdb_file = params['pdb_file']
device = params['device']
ppn = params['ppn']
restart = strtobool(params['restart'])
if not restart:
    minimize = strtobool(params['minimize'])
else:
    minimize = False
"""
End of reading parameters
"""

cgModel = cosmo.models.buildHPSModel(pdb_file, minimize=minimize, model=model, box_dimension=box_dimension)

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
    # Thermostats
    integrator = LangevinIntegrator(ref_t, tau_t, dt)
    # barostats
    barostat = MonteCarloBarostat(ref_p, ref_t, nstpcouple)
    # add barostat to system
    cgModel.system.addForce(barostat)
    simulation = Simulation(cgModel.topology, cgModel.system, integrator, platform, properties)
    simulation.loadCheckpoint(checkpoint)
    print(
        f"Restart from checkpoint, Time = {simulation.context.getState().getTime()}, Step= {simulation.context.getState().getStepCount()}")
    # number of steps remain to run
    nsteps_remain = md_steps - simulation.context.getState().getStepCount()
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(checkpoint, nstxout))
    simulation.reporters.append(DCDReporter(f'{protein_code}.dcd', nstxout, enforcePeriodicBox=bool(pbc), append=True))
    simulation.reporters.append(
        StateDataReporter(f'{protein_code}.log', nstlog, step=True, time=True, potentialEnergy=True,
                          kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t', append=True))
    simulation.step(nsteps_remain)
else:
    # Production phase, create new integrator and simulation context to reset number of steps
    # Thermostats
    integrator = LangevinIntegrator(ref_t, tau_t, dt)
    # barostat, using Monte Carlo algorithm to scale box dimension
    barostat = MonteCarloBarostat(ref_p, ref_t, nstpcouple)
    """
        add barostat to system. We don't need to call reinitialize() since the simulation is not created yet.
        Otherwise, if the simulation has been created and run, e.g in NVT, then we add barostat to run at NPT, we
        need to call reinitialize(), e.g: simulation.context.reinitialize(preserveState=True) such that simulation
        will see the changes. For details, check this: 
        https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/Changing%20Temperature%20and%20Pressure.html
    """
    cgModel.system.addForce(barostat)
    simulation = Simulation(cgModel.topology, cgModel.system, integrator, platform, properties)
    # Set initial positions: translate input coordinate, the coordinate is >=0
    xyz = np.array(cgModel.positions / unit.nanometer)
    xyz[:, 0] -= np.amin(xyz[:, 0])
    xyz[:, 1] -= np.amin(xyz[:, 1])
    xyz[:, 2] -= np.amin(xyz[:, 2])
    cgModel.positions = xyz * unit.nanometer
    simulation.context.setPositions(cgModel.positions)
    simulation.context.setVelocitiesToTemperature(ref_t)
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(checkpoint, nstxout))
    simulation.reporters.append(DCDReporter(f'{protein_code}.dcd', nstxout, enforcePeriodicBox=bool(pbc), append=False))
    simulation.reporters.append(
        StateDataReporter(f'{protein_code}.log', nstlog, step=True, time=True, potentialEnergy=True,
                          kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t', append=False))
    simulation.step(md_steps)

# write the last frame
last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(pbc)).getPositions()
PDBFile.writeFile(cgModel.topology, last_frame, open(f'{protein_code}_final.pdb', 'w'))
simulation.saveCheckpoint(checkpoint)

print("--- Finished in %s seconds ---" % (time.time() - start_time))
