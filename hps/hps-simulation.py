#!/home/qvv5013/anaconda3/envs/py310/bin/python
# Import OpenMM library
# Import hpsOpenMM library
import argparse
import hps

# Parse config file:
parser = argparse.ArgumentParser(description="\n Usage: python run_simulation.py -f md.ini ")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

hps_sim = hps.dynamics.Dynamics(args.input)
hps_sim.dynamics()
 
