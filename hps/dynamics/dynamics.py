import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads
from typing import Any

import numpy as np
import openmm
from openmm import unit
from parmed.exceptions import OpenMMWarning

from ..core import models

warnings.filterwarnings("ignore", category=OpenMMWarning)


class Dynamics:
    """
    Dynamics class contains two main functions: read config file and run simulation.
    User only need to provide config file, e.g md.ini and specify parameters control simulation there.

    Parameters
    ----------
    config_file: str
        control parameters for simulation

    Attributes
    ----------
    md_steps: int [1, steps]
        Number of steps to perform molecular dynamics simulation
    dt: float [0.01, ps]
        time step for integration
    nstxout: int [1, steps]
        number of steps that elapse between writing coordinates to output trajectory file,
        the last coordinates are always written
    nstlog: int [1, steps]
        number of steps that elapse between writing energies to the log file, the last energies are always written
    nstcomm: int [100, steps]
        frequency for center of mass motion removal
    model: str ['hps_urry']
        Hydropathy scale
    tcoupl: bool
        Using temperature coupling.
    ref_t: float [Kelvin]
        reference temperature for coupling
    tau_t: float [ps]
        ime constant for temperature coupling
    pcoupl: bool
        Pressure coupling
    ref_p: float [bar]
        The reference pressure for coupling.
    frequency_p: int [25, steps]
        The frequency for coupling the pressure.
    pbc: bool
        Use periodic boundary conditions.
    box_dimension: float or list of float
        Box dimension defined the unit cell, better to use rectangular for simplicity
    protein_code: str
        Prefix to write output file based on this parameter
    checkpoint: str
        Checkpoint file name
    pdb_file: str
        Input structure read to generate model.
    device: str
        Device to perform simulation [GPU/CPU] if CPU is used, then need to provide number of threads to run simulation.
    ppn: int [1, cores]
        In case simulation is run on CPU, use this parameter to control the number of threads to run simulation.
    restart: bool [No]
        If simulation run from beginning or restart from checkpoint.
    minimize: bool
        If simulation run from beginning then need to perform energy minimization. If simulation restarted, this
        parameters will be override to False.

    Returns
    -------

    """

    def __init__(self, config_file):
        self.md_steps: int = 1
        self.dt: float = 0.01
        self.nstxout: int = 1
        self.nstlog: int = 1
        self.nstcomm: int = 100
        self.model: str = 'hps_urry'

        # temperature coupling
        self.tcoupl: bool = True
        self.ref_t: float = 300.0
        self.tau_t = None

        # pressure coupling
        self.pcoupl: bool = False
        self.ref_p: float = 1.0
        self.frequency_p: int = 25

        # periodic boundary condition
        self.pbc: bool = False
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
        self.minimize: bool = True

        # call read_config function to initialize these attribute
        self.read_config(config_file)

        # store hps object to inspect
        self.hps_model = None

    def read_config(self, config_file):
        """
        Read simulation control parameters from config file *.ini into class attributes.

        TODO: check parameters in control file more carefully.
                Raise error and exit immediately if something wrong.
        """
        print(f"Reading simulation parameters from {config_file} file...")
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        config.read(config_file)
        params = config['OPTIONS']
        # Reading parameters
        # must have parameters
        if 'md_steps' in params:
            self.md_steps = int(params['md_steps'])
            print(f'Running simulation for {self.md_steps} steps')
        else:
            print(f'not found md_steps keyword in control file, use default md_steps = {self.md_steps}')
        if 'dt' in params:
            self.dt = float(params['dt']) * unit.picoseconds
            print(f'Timestep to integrate equation of motion: dt = {self.dt}')
        else:
            print(f'use default dt={self.dt} ps')
        self.nstxout = int(params['nstxout'])
        self.nstlog = int(params['nstlog'])
        if 'nstcomm' in params:
            self.nstcomm = int(params['nstcomm'])
            print(f'Remove center of mass motion every {self.nstcomm} step')
        else:
            print(f'Use default nstcomm = {self.nstcomm} steps for center of mass remover')
        if 'model' in params:
            self.model = params['model']
            print(f'Hydropathy scale is set to: {self.model}')
        else:
            print(f'Use default model={self.model}')

        # temperature coupling. Pretty sure it is always on.
        if 'tcoupl' in params:
            self.tcoupl = bool(strtobool(params['tcoupl']))
            print('Turn on temperature coupling')
        if self.tcoupl:
            self.ref_t = float(params['ref_t']) * unit.kelvin
            self.tau_t = float(params['tau_t']) / unit.picoseconds
            print(f'Temperature is set to: {self.ref_t} with time constant = {self.tau_t}')

        # Periodic Boundary condition
        if 'pbc' in params:
            self.pbc = bool(strtobool(params['pbc']))
        if self.pbc:
            self.box_dimension = loads(params['box_dimension'])
            print(f'Using Periodic boundary condition with box dimension: {self.box_dimension} nm')
        else:
            self.box_dimension = None
            print('Not using Periodic boundary condition')

        # Pressure coupling
        if 'pcoupl' in params:
            self.pcoupl = bool(strtobool(params['pcoupl']))
            # print('Pressure coupling')
        if self.pcoupl:
            assert self.pbc, f"Pressure coupling requires box dimensions and periodic boundary condition is on"
            self.ref_p = float(params['ref_p']) * unit.bar
            self.frequency_p = int(params['frequency_p'])
            print(f'Pressure is set to reference of {self.ref_p} with frequency of coupling {self.frequency_p}')

        # prefix IO files
        if 'pdb_file' in params:
            self.pdb_file = params['pdb_file']
            print(f'Input structure: {self.pdb_file}')
        else:
            raise EOFError('Not found pdb_file keyword in control file')

        self.protein_code = params['protein_code']
        print(f'Prefix use to write file: {self.protein_code}')

        if 'checkpoint' in params:
            self.checkpoint = params['checkpoint']
            print(f'Checkpoint will be writen to: {self.checkpoint}')
        else:
            self.checkpoint = self.pdb_file.split('.')[0] + '.chk'
            print(f'not found checkpoint keyword in control file. Set based on input structure: {self.checkpoint}')

        # Platform to run simulation
        if 'device' in params:
            self.device = params['device']
            print(f'Platform used to perform simulation: {self.device}')
        else:
            print(f'not found device keyword in control file. Set to {self.device}')

        self.ppn = params['ppn']

        # Restart simulation or run from beginning
        if 'restart' in params:
            self.restart = bool(strtobool(params['restart']))
            print(f'Restart simulation: {self.restart}')
        if not self.restart:
            self.minimize = bool(strtobool(params['minimize']))
            print(f'Perform Energy minimization of input structure: {self.minimize}')
        else:
            self.minimize = False
        print('__________________________________________________________________')
        """
        End of reading parameters
        """

    def run(self):
        """
        Run simulation
        Returns
        -------

        """

        # Initialize model
        self.hps_model = models.buildHPSModel(self.pdb_file, minimize=self.minimize, model=self.model,
                                              box_dimension=self.box_dimension)
        # dump topology PSF file and initial coordinate pdb file
        self.hps_model.dumpStructure(f'{self.protein_code}_init.pdb')
        self.hps_model.dumpTopology(f'{self.protein_code}.psf')
        self.hps_model.system.addForce(openmm.CMMotionRemover(self.nstcomm))

        # setup integrator and simulation object
        integrator = openmm.LangevinIntegrator(self.ref_t, self.tau_t, self.dt)
        if self.pcoupl:
            # if pressure coupling is on, add barostat force to the system.
            barostat = openmm.MonteCarloBarostat(self.ref_p, self.ref_t, self.frequency_p)
            self.hps_model.system.addForce(barostat)

        # Setup platform to run simulation
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

        simulation = openmm.app.Simulation(self.hps_model.topology, self.hps_model.system, integrator, platform,
                                           properties)
        start_time = time.time()
        if self.restart:
            simulation.loadCheckpoint(self.checkpoint)
            prev_time, prev_steps = simulation.context.getState().getTime(), simulation.context.getState().getStepCount()
            print(f"Restart from checkpoint, Time = {prev_time}, Step= {prev_steps}")
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
            # Set initial positions: translate input coordinate, the coordinate is >=0
            xyz = np.array(self.hps_model.positions / unit.nanometer)
            xyz[:, 0] -= np.amin(xyz[:, 0])
            xyz[:, 1] -= np.amin(xyz[:, 1])
            xyz[:, 2] -= np.amin(xyz[:, 2])
            self.hps_model.positions = xyz * unit.nanometer
            simulation.context.setPositions(self.hps_model.positions)
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
        openmm.app.PDBFile.writeFile(self.hps_model.topology, last_frame, open(f'{self.protein_code}_final.pdb', 'w'))
        simulation.saveCheckpoint(self.checkpoint)
        print("--- Finished in %s seconds ---" % (time.time() - start_time))
