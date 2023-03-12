#!/home/qvv5013/anaconda3/envs/hpsopenmm/bin/python

"""
This script is a command-line program for running a simulation of a protein system using the
HPS (hydrophobic-polar scale) force field. The script uses the argparse library to parse a simulation config file
specified by the user as a command-line argument. The script then imports the hps library,
and creates an instance of the Dynamics class, passing in the parsed config file.
Finally, the script runs the simulation by calling the run() method on the Dynamics class instance.
"""

# Import OpenMM library
# Import hpsOpenMM library
import argparse

import hps

# Parse config file:
parser = argparse.ArgumentParser(
    description="\n Usage: python run_simulation.py -f md.ini \n or hps-simulation -f md.ini")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

hps_sim = hps.dynamics.Dynamics(args.input)
hps_sim.run()
