# Import OpenMM library
# Import sbmOpenMM library
import os
import sys
import time
from sys import stdout

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

import hps

# MD parameter
# let's decide here now long we want to run the simulation and the file writing period
mdsteps = 1000000
dcdperiod = 1000
logperiod = 100

# stage of simulation equil, prod
stage='equil'
# which platform to run simulation: CPU/GPU
device ='GPU'

pdbname = 'asyn'
protein_code = 'asyn'
pdb_file = f'{pdbname}.pdb'

# Create an sbmOpenMM.system() object and store it in "sbmCAModelModel" variable.
cgModel = hps.models.getCAModel(pdb_file)

# dump Forcefield File
cgModel.dumpForceFieldData('forcefield.dat')
# dump Structure into PDB file for visualize
cgModel.dumpStructure(f'{protein_code}_{stage}_init.pdb')
cgModel.dumpTopology(f'{protein_code}.psf')


if device == 'GPU':

    # Run simulation on CUDA
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    ## incase many GPUs present, we can select which one to use
    properties["DeviceIndex"] = "0"

    integrator = LangevinIntegrator(298 * kelvin, 5 / picosecond, 10 * femtoseconds)
    simulation = Simulation(cgModel.topology, cgModel.system, integrator, platform, properties)

elif device == 'CPU':
    integrator = LangevinIntegrator(298*kelvin, 5/picosecond, 10*femtoseconds)
    simulation = Simulation(cgModel.topology, cgModel.system, integrator)

# Set initial positions
simulation.context.setPositions(cgModel.positions)
# set velocity by temperature
# simulation.context.setVelocitiesToTemperature(298 * kelvin)

# Add a DCD reporter that writes coordinates every 100 steps.
simulation.reporters.append(DCDReporter(f'{protein_code}_{stage}.dcd', dcdperiod))

## Reporter using sbmOpenMM reporter, which is customize to print energy group
## Add a SBM reporter that writes energies every 100 steps.
## This is useful for debugging when we know every force component, we can know which one is unrealistic
# simulation.reporters.append(sbmOpenMM.sbmReporter('energy.data', logperiod, sbmObject=cgModel,
#                                                   step=True, potentialEnergy=True, temperature=True))


# #Add a SBM reporter that prints energies every 20000 steps.
# simulation.reporters.append(sbmOpenMM.sbmReporter(stdout, logperiod, sbmObject=cgModel,
#                                                   step=True, potentialEnergy=True, temperature=True))

# If we don't need too details, we can use the default reporter from OpenMM
simulation.reporters.append(
    StateDataReporter(stdout, logperiod, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                          totalEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True,
                          totalSteps=mdsteps, separator='\t'))
simulation.reporters.append(
    StateDataReporter(f'{protein_code}_{stage}.log', logperiod, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
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
