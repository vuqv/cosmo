from dataclasses import dataclass
import configparser
from distutils.util import strtobool
from json import loads
from typing import Any
from openmm import unit


@dataclass
class Dynamics:
    """
    Dynamics class contains two main functions: read config file and run simulation.
    User only need to provide config file, e.g md.ini and specify parameters control simulation there.
    """
    md_steps: int
    dt: float
    nstxout: int
    nstlog: int
    model: str
    ref_t: float
    tau_t: float
    pbc: bool
    box_dimension: Any
    protein_code: str
    checkpoint: str
    pdb_file: str
    device: str
    ppn: int
    restart: bool
    minimize: bool

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
