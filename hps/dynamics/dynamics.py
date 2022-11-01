import configparser
from dataclasses import dataclass
from distutils.util import strtobool
from json import loads
from typing import Any

import numpy as np
import openmm
from openmm import unit

from ..core import models

# TODO: check input file and variables that user provide


@dataclass
class Dynamics:
    """
    Dynamics class contains two main functions: read config file and run simulation.
    User only need to provide config file, e.g md.ini and specify parameters control simulation there.
    """

    def __init__(self):
        self.md_steps: int = None
        self.dt: float = None
        self.nstxout: int = None
        self.nstlog: int = None
        self.model: str = None
        self.ref_t: float = None
        self.tau_t: float = None
        self.pbc: bool = None
        self.box_dimension: Any = None
        self.protein_code: str = None
        self.checkpoint: str = None
        self.pdb_file: str = None
        self.device: str = None
        self.ppn: int = None
        self.restart: bool = None
        self.minimize: bool = None

    def read_config(self, config_file):
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        config.read(config_file)
        params = config['DEFAULT']
        # Reading parameters
        self.md_steps = int(params['md_steps'])
        self.dt = float(params['dt']) * unit.picoseconds
        self.nstxout = int(params['nstxout'])
        self.nstlog = int(params['nstlog'])
        self.model = params['model']
        self.ref_t = float(params['ref_t']) * unit.kelvin
        self.tau_t = float(params['tau_t']) / unit.picoseconds
        self.pbc = strtobool(params['pbc'])
        if self.pbc:
            self.box_dimension = loads(params['box_dimension'])
        else:
            self.box_dimension = None

        self.protein_code = params['protein_code']
        self.checkpoint = params['checkpoint']
        self.pdb_file = params['pdb_file']
        self.device = params['device']
        self.ppn = params['ppn']
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
