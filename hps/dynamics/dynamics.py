import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import numpy as np
import openmm as mm
import openmm.unit as unit
from parmed.exceptions import OpenMMWarning

from ..core import models

warnings.filterwarnings("ignore", category=OpenMMWarning)


class Dynamics:
    """
    The Dynamics class is used to perform molecular dynamics simulations.
    It contains two main functions: reading a configuration file and running a simulation.
    To use the class, a user needs to provide a configuration file (e.g md.ini)
    and specify the parameters that control the simulation.


    Parameters
    ----------
    config_file : str
        The path to the configuration file that contains the parameters for the simulation.

    Attributes
    ----------
    md_steps : int, optional (default: 1000)
        The number of steps to perform in the molecular dynamics simulation.
    dt : float, optional (default: 0.01) [ps]
        The time step for integration in picoseconds.
    nstxout : int, optional (default: 10)
        The number of steps between writing coordinates to the output trajectory file. The last coordinates are always written.
    nstlog : int, optional (default: 10)
        The number of steps between writing energies to the log file. The last energies are always written.
    nstcomm : int, optional (default: 100)
        The frequency for center of mass motion removal.
    model : str, optional (default: 'hps_urry')
        Hydropathy scale used in the simulation.
    tcoupl : bool, optional (default: False)
        Indicates whether temperature coupling is used in the simulation.
    ref_t : float, optional (default: 300.0) [K]
        The reference temperature for coupling in Kelvin.
    tau_t : float, optional (default: 0.1) [ps]
        The time constant for temperature coupling in picoseconds.
    pcoupl : bool, optional (default: False)
        Indicates whether pressure coupling is used in the simulation.
    ref_p : float, optional (default: 1.0) [bar]
        The reference pressure for coupling in bar.
    frequency_p : int, optional (default: 25)
        The frequency for coupling the pressure.
    pbc : bool, optional (default: True)
        Indicates whether periodic boundary conditions are used in the simulation.
    box_dimension : float or list of float, optional (default: None)
        The dimension of the box used in the simulation. It is better to use rectangular for simplicity.
    protein_code : str, optional (default: None)
        A prefix to write output files based on this parameter.
    checkpoint : str, optional (default: None)
        The name of the checkpoint file.
    pdb_file : str, optional (default: None)
        The input structure used to generate the model.
    device : str, optional (default: 'GPU')
        The device used to perform the simulation. Options are 'GPU' or 'CPU'.
    ppn : int, optional (default: 1)
        In case the simulation is run on a CPU, this parameter controls the number of threads used to run the simulation.
    restart : bool, optional (default: False)
        Indicates whether the simulation should be run from the beginning or restarted from a checkpoint.
    minimize : bool, optional (default: True)
        Indicates whether energy minimization should be performed at the start of the simulation


    Returns
    -------

    """

    def __init__(self, config_file):
        # Set default values
        self.md_steps = 1000
        self.dt = 0.01
        self.nstxout = 10
        self.nstlog = 10
        self.nstcomm = 100
        self.model = 'hps_urry'
        self.tcoupl = True
        self.ref_t = 300.0
        self.pcoupl = False
        self.ref_p = 1.0
        self.frequency_p = 25
        self.pbc = False
        self.device = 'CPU'
        self.ppn = 1
        self.restart = False
        self.minimize = True

        # Other attributes
        self.tau_t = None
        self.box_dimension = None
        self.protein_code = None
        self.checkpoint = None
        self.pdb_file = None

        # Initialize attributes with config file
        self.read_config(config_file)

        self.hps_model = None

    def read_config(self, config_file):
        """
           Read simulation control parameters from config file *.ini into class attributes.

           Parameters
           ----------
           config_file: string, default=None
               Control file store simulation parameters.

           TODO: check parameters in control file more carefully.
                   Raise error and exit immediately if something wrong.
           """
        print(f"Reading simulation parameters from {config_file} file...")
        config = configparser.ConfigParser()
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        config.read(config_file)
        params = config['OPTIONS']

        self.md_steps = int(params.get('md_steps', self.md_steps))
        print(f'Setting number of simulation steps to: {self.md_steps}')
        self.dt = float(params.get('dt', self.dt)) * unit.picoseconds
        print(f'Setting timestep for integration of equations of motion to: {self.dt}')
        self.nstxout = int(params['nstxout'])
        self.nstlog = int(params['nstlog'])
        self.nstcomm = int(params.get('nstcomm', self.nstcomm))
        print(f'Setting frequency of center of mass motion removal to every {self.nstcomm} steps')
        self.model = params.get('model', self.model)
        print(f'Setting hydropathy scale to: {self.model}')
        self.tcoupl = bool(strtobool(params.get('tcoupl', self.tcoupl)))
        if self.tcoupl:
            self.ref_t = float(params['ref_t']) * unit.kelvin
            self.tau_t = float(params['tau_t']) / unit.picoseconds
            print(
                f'Turning on temperature coupling with reference temperature: {self.ref_t} and time constant: {self.tau_t}')
        else:
            print("Temperature coupling is off")
        self.pbc = bool(strtobool(params.get('pbc', self.pbc)))
        if self.pbc:
            self.box_dimension = loads(params['box_dimension'])
            print(f'Turning on periodic boundary conditions with box dimension: {self.box_dimension} nm')
        else:
            self.box_dimension = None
            print('Periodic boundary conditions are off')
        self.pcoupl = bool(strtobool(params.get('pcoupl', self.pcoupl)))
        if self.pcoupl:
            assert self.pbc, "Pressure coupling requires box dimensions and periodic boundary condition is on"
            self.ref_p = float(params['ref_p']) * unit.bar
            self.frequency_p = int(params['frequency_p'])
            print(f'Pressure is set to reference of {self.ref_p} with frequency of coupling {self.frequency_p}')
        else:
            print("Pressure coupling is off")
        self.pdb_file = params['pdb_file']
        print(f'Input structure: {self.pdb_file}')
        self.protein_code = params['protein_code']
        print(f'Prefix use to write file: {self.protein_code}')
        self.checkpoint = params.get('checkpoint', self.checkpoint)
        self.device = params.get('device', self.device)
        print(f'Running simulation on {self.device}')
        if self.device == "CPU":
            self.ppn = int(params.get('ppn', self.ppn))
            print(f'Using {self.ppn} threads')
        self.restart = bool(strtobool(params.get('restart', self.restart)))
        print(f'Restart simulation: {self.restart}')
        if self.restart:
            self.minimize = False
        else:
            self.minimize = bool(strtobool(params.get('minimize', self.minimize)))
            print(f'Perform Energy minimization of input structure: {self.minimize}')

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
        if self.model in ['hps_kr', 'hps_urry', 'hps_ss']:
            """ current dumpForceFieldData function can only write the standard format of forcefield which require
             sigma, epsilon for each residue.
            """
            self.hps_model.dumpForceFieldData(f'{self.protein_code}_ff.dat')
        self.hps_model.dumpStructure(f'{self.protein_code}_init.pdb')
        self.hps_model.dumpTopology(f'{self.protein_code}.psf')
        self.hps_model.system.addForce(mm.CMMotionRemover(self.nstcomm))

        # setup integrator and simulation object
        integrator = mm.LangevinIntegrator(self.ref_t, self.tau_t, self.dt)
        if self.pcoupl:
            # if pressure coupling is on, add barostat force to the system.
            barostat = mm.MonteCarloBarostat(self.ref_p, self.ref_t, self.frequency_p)
            self.hps_model.system.addForce(barostat)

        # Setup platform to run simulation
        if self.device == 'GPU':
            # Run simulation on CUDA
            print(f"Running simulation on GPU CUDA")
            platform = mm.Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
            # in case of many GPUs present, we can select which one to use

        elif self.device == 'CPU':
            print(f"Running simulation on CPU using {self.ppn} cores")
            platform = mm.Platform.getPlatformByName('CPU')
            properties = {'Threads': str(self.ppn)}

        simulation = mm.app.Simulation(self.hps_model.topology, self.hps_model.system, integrator, platform,
                                       properties)
        start_time = time.time()
        if self.restart:
            simulation.loadCheckpoint(self.checkpoint)
            print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
            nsteps_remain = self.md_steps - simulation.context.getState().getStepCount()
        else:
            xyz = np.array(self.hps_model.positions / unit.nanometer)
            xyz[:, 0] -= np.amin(xyz[:, 0])
            xyz[:, 1] -= np.amin(xyz[:, 1])
            xyz[:, 2] -= np.amin(xyz[:, 2])
            self.hps_model.positions = xyz * unit.nanometer
            simulation.context.setPositions(self.hps_model.positions)
            simulation.context.setVelocitiesToTemperature(self.ref_t)
            nsteps_remain = self.md_steps

        simulation.reporters = []
        simulation.reporters.append(mm.app.CheckpointReporter(self.checkpoint, self.nstxout))
        simulation.reporters.append(
            mm.app.DCDReporter(f'{self.protein_code}.dcd', self.nstxout, enforcePeriodicBox=bool(self.pbc),
                               append=self.restart))
        simulation.reporters.append(
            mm.app.StateDataReporter(f'{self.protein_code}.log', self.nstlog, step=True, time=True,
                                     potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                     temperature=True, remainingTime=True, speed=True,
                                     totalSteps=self.md_steps, separator='\t', append=self.restart))
        simulation.step(nsteps_remain)

        # write the last frame
        last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(self.pbc)).getPositions()
        mm.app.PDBFile.writeFile(self.hps_model.topology, last_frame, open(f'{self.protein_code}_final.pdb', 'w'))
        simulation.saveCheckpoint(self.checkpoint)
        print("--- Finished in %s seconds ---" % (time.time() - start_time))
