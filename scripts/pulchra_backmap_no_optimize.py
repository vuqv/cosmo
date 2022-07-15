import argparse
import os

import mdtraj as md
import numpy as np

# parser args
parser = argparse.ArgumentParser(description="Backmap CG to All-atom")
parser.add_argument('-top', '-p', type=str, help='Topology')
parser.add_argument('-traj', '-f', type=str, help='Trajectory')
parser.add_argument('-output', '-o', type=str, help='Output Trajectory', default='output')
parser.add_argument('-begin', '-b', type=int, help='Starting frame (default: 0)', default=0)
parser.add_argument('-end', '-e', type=int, help='End frame (default: last)', default=-1)

args = parser.parse_args()

t = md.load_dcd(args.traj, args.top)

if args.begin != 0:
    frame_start = args.begin
else:
    frame_start = 0

if args.end == -1:
    frame_end = t.n_frames
else:
    frame_end = args.end

# save template to get number of atoms in backmap structure
t[0].save('frame.pdb')
os.system('pulchra frame.pdb')
traj_0 = md.load('frame.rebuilt.pdb')

traj_coor = np.empty((0, traj_0.n_atoms, 3))

for i in range(frame_start, frame_end):
    print(f'backmapping frame: {i}')
    t[i].save('frame.pdb')
    os.system('pulchra frame.pdb')
    traj_temp = md.load('frame.rebuilt.pdb')
    traj_coor = np.append(traj_coor, traj_temp.xyz)

n_frames = frame_end - frame_start
new_traj_coor = traj_coor.reshape(n_frames, traj_0.n_atoms, 3)
coor_2save = md.Trajectory(new_traj_coor, traj_0.top)
coor_2save.save_dcd(f'{args.output}.dcd')
