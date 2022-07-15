#!/usr/bin/env python
# coding: utf-8

import mdtraj as md
import numpy as np

traj = md.load('../hcyp_prod.dcd', top='../hcyp.psf')

aa_masses = {'ALA': 71.08, 'ARG': 156.20, 'ASN': 114.10,
             'ASP': 115.10, 'CYS': 103.10, 'GLU': 129.10,
             'GLN': 128.10, 'GLY': 57.05, 'HIS': 137.10,
             'ILE': 113.20, 'LEU': 113.20, 'LYS': 128.20,
             'MET': 131.20, 'PHE': 147.20, 'PRO': 97.12,
             'SER': 87.08, 'THR': 101.10, 'TRP': 186.20,
             'TYR': 163.20, 'VAL': 99.07}

masses = np.zeros(len(list(traj.topology.residues)))

for i, r in enumerate(traj.topology.residues):
    #     print(i,r.name)
    masses[i] = aa_masses[r.name]

# unit in nm
Rg = md.compute_rg(traj, masses)

print(Rg.mean(), Rg.std())
