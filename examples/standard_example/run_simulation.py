# Import OpenMM library
# Import hpsOpenMM library
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import numpy as np
import openmm as mm
from openmm import unit
# from openmm.app import *
from parmed.exceptions import OpenMMWarning

import hps

warnings.filterwarnings("ignore", category=OpenMMWarning)

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

# Parse config file:
parser = argparse.ArgumentParser(
    description="\n Usage: python run_simulation.py -f md.ini \n or hps-simulation -f md.ini")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()
config_file = args.input
print(f"Reading simulation parameters from {config_file} file...")
config = configparser.ConfigParser()
config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
config.read(config_file)
params = config['OPTIONS']

md_steps = int(params.get('md_steps', md_steps))
print(f'Setting number of simulation steps to: {md_steps}')
dt = float(params.get('dt', dt)) * unit.picoseconds
print(f'Setting timestep for integration of equations of motion to: {dt}')
nstxout = int(params['nstxout'])
print(f'Setting number of steps to write checkpoint and coordinate: {nstxout}')
nstlog = int(params['nstlog'])
print(f'Setting number of steps to write logfile: {nstlog}')
nstcomm = int(params.get('nstcomm', nstcomm))
print(f'Setting frequency of center of mass motion removal to every {nstcomm} steps')
model = params.get('model', model)
print(f'Setting hydropathy scale to: {model}')
tcoupl = bool(strtobool(params.get('tcoupl', tcoupl)))
if tcoupl:
    ref_t = float(params['ref_t']) * unit.kelvin
    tau_t = float(params['tau_t']) / unit.picoseconds
    print(
        f'Turning on temperature coupling with reference temperature: {ref_t} and time constant: {tau_t}')
else:
    print("Temperature coupling is off")
pbc = bool(strtobool(params.get('pbc', pbc)))
if pbc:
    box_dimension = loads(params['box_dimension'])
    print(f'Turning on periodic boundary conditions with box dimension: {box_dimension} nm')
else:
    box_dimension = None
    print('Periodic boundary conditions are off')
pcoupl = bool(strtobool(params.get('pcoupl', pcoupl)))
if pcoupl:
    assert pbc, "Pressure coupling requires box dimensions and periodic boundary condition is on"
    ref_p = float(params['ref_p']) * unit.bar
    frequency_p = int(params['frequency_p'])
    print(f'Pressure is set to reference of {ref_p} with frequency of coupling {frequency_p}')
else:
    print("Pressure coupling is off")
pdb_file = params['pdb_file']
print(f'Input structure: {pdb_file}')
protein_code = params['protein_code']
print(f'Prefix use to write file: {protein_code}')
checkpoint = params.get('checkpoint', checkpoint)
device = params.get('device', device)
print(f'Running simulation on {device}')
if device == "CPU":
    ppn = int(params.get('ppn', ppn))
    print(f'Using {ppn} threads')
restart = bool(strtobool(params.get('restart', restart)))
print(f'Restart simulation: {restart}')
if restart:
    minimize = False
else:
    minimize = bool(strtobool(params.get('minimize', minimize)))
    print(f'Perform Energy minimization of input structure: {minimize}')

"""
End of reading parameters
"""

hps_model = hps.models.buildHPSModel(pdb_file, minimize=minimize, model=model, box_dimension=box_dimension)

# Remove center of mass motion
hps_model.system.addForce(mm.CMMotionRemover(nstcomm))

# dump Forcefield File
if model in ['hps_kr', 'hps_urry', 'hps_ss']:
    """ current dumpForceFieldData function can only write the standard format of forcefield which require
     sigma, epsilon for each residue.
    """
    hps_model.dumpForceFieldData(f'{protein_code}_ff.dat')
# dump Structure into PDB file for visualize
hps_model.dumpStructure(f'{protein_code}_init.pdb')
hps_model.dumpTopology(f'{protein_code}.psf')

if device == 'GPU':
    # Run simulation on CUDA
    print(f"Running simulation on GPU CUDA")
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
    # in case of many GPUs present, we can select which one to use


elif device == 'CPU':
    print(f"Running simulation on CPU using {ppn} cores")
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {'Threads': str(ppn)}

print('Simulation started')
start_time = time.time()

integrator = mm.LangevinIntegrator(ref_t, tau_t, dt)
simulation = mm.app.Simulation(hps_model.topology, hps_model.system, integrator, platform,
                               properties)
if restart:
    simulation.loadCheckpoint(checkpoint)
    print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
    nsteps_remain = md_steps - simulation.context.getState().getStepCount()
else:
    xyz = np.array(hps_model.positions / unit.nanometer)
    xyz[:, 0] -= np.amin(xyz[:, 0])
    xyz[:, 1] -= np.amin(xyz[:, 1])
    xyz[:, 2] -= np.amin(xyz[:, 2])
    hps_model.positions = xyz * unit.nanometer
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

# write the last frame
last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(pbc)).getPositions()
mm.app.PDBFile.writeFile(hps_model.topology, last_frame, open(f'{protein_code}_final.pdb', 'w'))
simulation.saveCheckpoint(checkpoint)

print("--- Finished in %s seconds ---" % (time.time() - start_time))
