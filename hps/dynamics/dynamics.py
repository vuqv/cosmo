import configparser
from distutils.util import strtobool
from json import loads
from typing import Any

import numpy as np
import openmm
from openmm import unit

from ..core import models


# TODO: check input file and variables that user provide


class Dynamics:
    """
    Dynamics class contains two main functions: read config file and run simulation.
    User only need to provide config file, e.g md.ini and specify parameters control simulation there.
    """

    def __init__(self):
        self.md_steps: int = 1
        self.dt: float = 0.01
        self.nstxout: int = 1
        self.nstlog: int = 1
        self.model: str = 'hps_urry'

        # temperature coupling
        self.tcoupl = None
        self.ref_t = None
        self.tau_t = None

        # pressure coupling
        self.pcoupl = None
        self.ref_p = None
        self.frequency_p = None

        # periodic boundary condition
        self.pbc = None
        self.box_dimension: Any = None

        # prefix and io file name
        self.protein_code = None
        self.checkpoint = None
        self.pdb_file = None

        # simulation platform
        self.device: str = 'CPU'
        self.ppn: int = 1

        # restart simulation or run from beginning
        self.restart: bool = False
        self.minimize = None

    def read_config(self, config_file):
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        config.read(config_file)
        params = config['OPTIONS']
        # Reading parameters
        """
        This can be check like:
        if 'md_steps' in params:
            do something
        else:
            raise error of just ignore this variable if it is not necessary.
        """
        # must have parameters
        self.md_steps = int(params['md_steps'])
        self.dt = float(params['dt']) * unit.picoseconds
        self.nstxout = int(params['nstxout'])
        self.nstlog = int(params['nstlog'])
        self.model = params['model']

        # temperature coupling. Pretty sure it is always on.
        self.tcoupl = strtobool(params['tcoupl'])
        if self.tcoupl:
            self.ref_t = float(params['ref_t']) * unit.kelvin
            self.tau_t = float(params['tau_t']) / unit.picoseconds

        # Periodic Boundary condition
        self.pbc = strtobool(params['pbc'])
        if self.pbc:
            self.box_dimension = loads(params['box_dimension'])
        else:
            self.box_dimension = None

        # Pressure coupling
        self.pcoupl = strtobool(params['pcoupl'])
        if self.pcoupl:
            assert self.pbc, f"Pressure coupling requires box dimensions and periodic boundary condition is on"
            self.ref_p = float(params['ref_p']) * unit.bar
            self.frequency_p = int(params['frequency_p'])

        # prefix IO files
        self.protein_code = params['protein_code']
        self.checkpoint = params['checkpoint']
        self.pdb_file = params['pdb_file']

        # Platform to run simulation
        self.device = params['device']
        self.ppn = params['ppn']

        # Restart simulation or run from beginning
        self.restart = strtobool(params['restart'])
        if not self.restart:
            self.minimize = strtobool(params['minimize'])
        else:
            self.minimize = False
        """
        End of reading parameters
        """

    def dynamics(self):
        """
        Run simulation
        Returns
        -------

        """
        hps_model = models.buildHPSModel(self.pdb_file, minimize=self.minimize, hps_scale=self.model,
                                         box_dimension=self.box_dimension)
        # dump topology PSF file and initial coordinate pdb file
        hps_model.dumpStructure(f'{self.protein_code}_init.pdb')
        hps_model.dumpTopology(f'{self.protein_code}.psf')

        if self.device == 'GPU':
            # Run simulation on CUDA
            print(f"Running simulation on GPU CUDA")
            platform = openmm.Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
            # in case of many GPUs present, we can select which one to use

        elif self.device == 'CPU':
            print(f"Running simulation on CPU using {self.ppn} cores")
            platform = openmm.Platform.getPlatformByName('CPU')
            properties = {'Threads': str(self.ppn)}

        if self.restart:
            # simulation reporter
            integrator = openmm.LangevinIntegrator(self.ref_t, self.tau_t, self.dt)
            if self.pcoupl:
                # if pressure coupling is on, add barostat force to the system.
                barostat = openmm.MonteCarloBarostat(self.ref_p, self.ref_t, self.frequency_p)
                hps_model.system.addForce(barostat)

            simulation = openmm.app.Simulation(hps_model.topology, hps_model.system, integrator, platform, properties)
            simulation.loadCheckpoint(self.checkpoint)
            print(
                f"Restart from checkpoint, Time = {simulation.context.getState().getTime()}, Step= {simulation.context.getState().getStepCount()}")
            # number of steps remain to run
            nsteps_remain = self.md_steps - simulation.context.getState().getStepCount()
            simulation.reporters = []
            simulation.reporters.append(openmm.app.CheckpointReporter(self.checkpoint, self.nstxout))
            simulation.reporters.append(
                openmm.app.DCDReporter(f'{self.protein_code}.dcd', self.nstxout, enforcePeriodicBox=bool(self.pbc),
                                       append=True))
            simulation.reporters.append(
                openmm.app.StateDataReporter(f'{self.protein_code}.log', self.nstlog, step=True, time=True,
                                             potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                             temperature=True, remainingTime=True, speed=True,
                                             totalSteps=self.md_steps, separator='\t', append=True))
            simulation.step(nsteps_remain)
        else:
            # Production phase, create new integrator and simulation context to reset number of steps
            integrator = openmm.LangevinIntegrator(self.ref_t, self.tau_t, self.dt)
            if self.pcoupl:
                # if pressure coupling is on, add barostat force to the system.
                barostat = openmm.MonteCarloBarostat(self.ref_p, self.ref_t, self.frequency_p)
                hps_model.system.addForce(barostat)

            simulation = openmm.app.Simulation(hps_model.topology, hps_model.system, integrator, platform, properties)
            # Set initial positions: translate input coordinate, the coordinate is >=0
            xyz = np.array(hps_model.positions / unit.nanometer)
            xyz[:, 0] -= np.amin(xyz[:, 0])
            xyz[:, 1] -= np.amin(xyz[:, 1])
            xyz[:, 2] -= np.amin(xyz[:, 2])
            hps_model.positions = xyz * unit.nanometer
            simulation.context.setPositions(hps_model.positions)
            simulation.context.setVelocitiesToTemperature(self.ref_t)
            simulation.reporters = []
            simulation.reporters.append(openmm.app.CheckpointReporter(self.checkpoint, self.nstxout))
            simulation.reporters.append(
                openmm.app.DCDReporter(f'{self.protein_code}.dcd', self.nstxout, enforcePeriodicBox=bool(self.pbc),
                                       append=False))
            simulation.reporters.append(
                openmm.app.StateDataReporter(f'{self.protein_code}.log', self.nstlog, step=True, time=True,
                                             potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                             temperature=True, remainingTime=True, speed=True,
                                             totalSteps=self.md_steps, separator='\t', append=False))
            simulation.step(self.md_steps)

        # write the last frame
        last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(self.pbc)).getPositions()
        openmm.app.PDBFile.writeFile(hps_model.topology, last_frame, open(f'{self.protein_code}_final.pdb', 'w'))
        simulation.saveCheckpoint(self.checkpoint)
