# Import OpenMM library
# Import COSMO library
import argparse

import cosmo

# Parse config file:
parser = argparse.ArgumentParser(
    description="\n Usage: python run_simulation.py -f md.ini \n or cosmo-simulation -f md.ini")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

hps_sim = cosmo.dynamics.Dynamics(args.input)
hps_sim.run()
