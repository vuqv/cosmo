# Import OpenMM library
# Import hpsOpenMM library
import argparse
import warnings

# from openmm.app import *
from parmed.exceptions import OpenMMWarning

import hps

warnings.filterwarnings("ignore", category=OpenMMWarning)

# Parse config file:
parser = argparse.ArgumentParser(description="\n Usage: python run_simulation.py -f md.ini ")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

hps_sim = hps.dynamics.Dynamics()

hps_sim.read_config(args.input)
print(vars(hps_sim))
hps_sim.dynamics()
