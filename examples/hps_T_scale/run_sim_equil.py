# Import OpenMM library
# Import hpsOpenMM library
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads
from sys import stdout

import numpy as np
from openmm import *
from openmm.app import *
from parmed.exceptions import OpenMMWarning

import hps

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
lj_eps = float(params['LJ_eps'])
ref_t = float(params['ref_t']) * unit.kelvin
tau_t = float(params['tau_t']) / unit.picoseconds
pbc = strtobool(params['pbc'])
if pbc:
    box_dimension = loads(params['box_dimension'])
else:
    box_dimension = None

protein_code = params['protein_code']
checkpoint = params['checkpoint']
pdb_file = params['pdb_file']
machine = params['machine']
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

cgModel = hps.models.buildHPSModel(pdb_file, minimize=minimize, hps_scale=model, box_dimension=box_dimension)

# ashbaugh force has 1 global parameter: epsilon
"""
The eq. 7d in paper has typo mistake in 'T' term, should be minus

lambda(H) = lambda_HPS + -25.47406054 + 0.14536949*T - 0.00020058*T^2
lambda(A) = lambda_HPS + -26.18813544 + 0.15033133*T - 0.00020919*T^2
lambda(O) = lambda_HPS + 2.45798364 - 0.01432925*T + 0.00002037*T^2
lambda(P) = lambda_HPS + 11.795 - 0.067679*T + 0.00009411*T^2
lambda(C) = lambda_HPS + 9.66118972 - 0.05425780*T + 0.00007312*T^2

change hps parameter for every atoms:
using: 
lj=simulation.system.getForce(2) # force 2 in the system is ashbaugh force
lj.getParticleParameters(138) # get Particle parameter of atom 138 (0 based)

lj force has 2 per particle parameter: (sigma, hps)
lj.setParticleParameters(138, (1,2)): change params for particle 138: sigma=1, hps=2. (1,2) can be [1,2]

loop through all particle in cgModel:
for i, atom in enumerate(cgModel.atoms):
    if atom.residue.name in [Hydrophobic/hydrophilic/...]
        new_hps = hps(T)
        old_params = lj.getParticleParameters(i)
        lj.setParticleParameters(i, (old_params[0], new_hps))

"""

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
    simulation.loadCheckpoint(checkpoint)
    print(
        f"Restart from checkpoint, Time = {simulation.context.getState().getTime()}, Step= {simulation.context.getState().getStepCount()}")
    # number of steps remain to run
    nsteps_remain = md_steps - simulation.context.getState().getStepCount()
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(checkpoint, nstxout))
    simulation.reporters.append(DCDReporter(f'{protein_code}.dcd', nstxout, append=True))
    if machine == 'local':
        simulation.reporters.append(
            StateDataReporter(stdout, nstlog, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                              totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                              totalSteps=md_steps, separator='\t'))
    simulation.reporters.append(
        StateDataReporter(f'{protein_code}.log', nstlog, step=True, time=True, potentialEnergy=True,
                          kineticEnergy=True,
                          totalEnergy=True, temperature=True, remainingTime=True, speed=True,
                          totalSteps=md_steps, separator='\t', append=True))
    simulation.step(nsteps_remain)
else:
    # Production phase, create new integrator and simulation context to reset number of steps
    integrator = LangevinIntegrator(ref_t, tau_t, dt)
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
    simulation.reporters.append(DCDReporter(f'{protein_code}.dcd', nstxout, append=False))
    if machine == 'local':
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
simulation.saveCheckpoint(checkpoint)

print("--- Finished in %s seconds ---" % (time.time() - start_time))
