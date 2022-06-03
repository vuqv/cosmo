# Import OpenMM library
from os import system
import time
from sys import stdout
import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *

import hps
import warnings
from parmed.exceptions import OpenMMWarning
warnings.filterwarnings("ignore", category=OpenMMWarning)
# MD parameter
# let's decide here now as long as we want to run the simulation and the file writing period
mdsteps = 200_000
dcdperiod = 100
logperiod = 100
# stage of simulation equil, prod
stage = 'equil'
# which platform to run simulation: CPU/GPU
device = 'GPU'

pdbname = 'asyn'
protein_code = f'{pdbname}'
pdb_file = f'{pdbname}.pdb'

# Create an sbmOpenMM.system() object and store it in "sbmCAModelModel" variable.
# No PBC
# cgModel = hps.models.getCAModel(pdb_file, hps_scale='kr')
# With PBC, or can call: box_dimension= [Lx, Ly, Lz]
cgModel = hps.models.getCAModel(pdb_file, hps_scale='kr', box_dimension= 100)
if device == 'GPU':
    # Run simulation on CUDA
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    # in case of many GPUs present, we can select which one to use
    # properties["DeviceIndex"] = "0"

elif device == 'CPU':
    platform = Platform.getPlatformByName('CPU')
    properties = {'Threads': str(8)}

integrator = LangevinIntegrator(500 * kelvin, 0.01 / picosecond, 10 * femtoseconds)
simulation = Simulation(cgModel.topology, cgModel.system, integrator, platform, properties)

# Set initial positions
xyz = np.array(cgModel.positions/nanometer)
xyz[:,0] -= np.amin(xyz[:,0])
xyz[:,1] -= np.amin(xyz[:,1])
xyz[:,2] -= np.amin(xyz[:,2])
cgModel.positions = xyz*nanometer
simulation.context.setPositions(cgModel.positions)
# set velocity by temperature
simulation.context.setVelocitiesToTemperature(298 * kelvin)


# # dump Forcefield File
cgModel.dumpForceFieldData('forcefield.dat')
# dump Structure into PDB file for visualize
cgModel.dumpStructure(f'{protein_code}_{stage}_init.pdb')
cgModel.dumpTopology(f'{protein_code}.psf')


# Add a DCD reporter that writes coordinates every 100 steps.
# simulation.reporters.append(DCDReporter(f'{protein_code}_{stage}.dcd', dcdperiod, enforcePeriodicBox=True))
simulation.reporters.append(DCDReporter(f'{protein_code}_{stage}.dcd', dcdperiod))

# If we don't need too details, we can use the default reporter from OpenMM
simulation.reporters.append(
    StateDataReporter(stdout, logperiod, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                      totalEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True,
                      totalSteps=mdsteps, separator='\t'))
simulation.reporters.append(
    StateDataReporter(f'{protein_code}_{stage}.log', logperiod, step=True, time=True, potentialEnergy=True,
                      kineticEnergy=True,
                      totalEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True,
                      totalSteps=mdsteps, separator='\t'))

print('Simulation started')
start_time = time.time()
simulation.step(mdsteps)

# write the last frame
lastframe = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(cgModel.topology, lastframe, open(f'{protein_code}_{stage}_final.pdb', 'w'))
simulation.saveCheckpoint(f'checkpoint_{stage}.chk')

print("--- Finished in %s seconds ---" % (time.time() - start_time))
