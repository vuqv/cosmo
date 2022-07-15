#!/usr/bin/env python
# coding: utf-8
r"""
Calculate SAXS profile from simulation.
CG parameters for form factor taken from paper:
1. D. Tong, S. Yang, L. Lu,
Accurate optimization of amino acid form factors for computing small-angle X-ray scattering intensity of atomistic protein structures.
J. Appl. Crystallogr. 49, 1148–1161 (2016).
"""

import itertools
import time
from math import sin

import MDAnalysis as mda
import numpy as np
import pandas as pd
from numba import njit

# In[2]:


u = mda.Universe('nup.psf', 'nup_equil.dcd')
protein = u.select_atoms('protein')
n_atoms = len(protein)
resnames = protein.resnames

n_frames = u.trajectory.n_frames
print(n_frames)
print(n_atoms)

aa_order = ['q', 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN',
            'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
form_factor_data = np.loadtxt('../table_edm_COE_1bead.dat')
form_factor = pd.DataFrame(form_factor_data, columns=aa_order)
# Note that the form-factor are given with the wawelength in unit of A^-1
# convert form_factor from pandas to numpy array
form_factor_np = np.zeros((51, n_atoms), dtype=float)
for i in range(form_factor_np.shape[0]):
    for j, r in enumerate(resnames):
        form_factor_np[i, j] = form_factor.loc[i, r]


@njit(fastmath=True)
def get_I(idx, qs, distance_pair):
    I = 0.0
    if idx == 0:
        for i in range(n_atoms):
            for j in range(n_atoms):
                if i != j:
                    I += form_factor_np[idx, i] * form_factor_np[idx, j]
    # Iq, otherwise
    else:
        for i in range(n_atoms):
            for j in range(n_atoms):
                if i != j:
                    I += form_factor_np[idx, i] * form_factor_np[idx, j] * sin(qs[idx] * distance_pair[i, j]) / (
                            qs[idx] * distance_pair[i, j])
    return I


def numba_saxs_single_frame(positions):
    n_atoms = len(positions)
    pair_array = np.asarray(list(itertools.product(positions, positions)))
    distance_diff = pair_array[:, 0, :] - pair_array[:, 1, :]
    distance_pair = np.linalg.norm(distance_diff, axis=1).reshape(n_atoms, n_atoms)

    qs = np.arange(0, 0.501, 0.01)
    Iq = np.zeros(len(qs))
    for idx in range(len(qs)):
        Iq[idx] = get_I(idx, qs, distance_pair)

    return Iq / Iq[0]


# calculate SAXS along trajectory
saxs_total = np.empty((51, 0), dtype=float)
begin_time = time.time()
for ts in u.trajectory[1000:]:
    positions = protein.positions
    #     res=saxs_single_frame(positions)
    res = numba_saxs_single_frame(positions)
    saxs_total = np.append(saxs_total, np.array([res]).transpose(), axis=1)

# get mean value 
saxs_mean = np.mean(saxs_total, axis=1)
end_time = time.time()
total_run_time = end_time - begin_time
print(f'#REPORT: Total execution time: {total_run_time / 60.0:.3f} mins\n')

qs = np.arange(0, 0.501, 0.01)
new_arr = np.array([qs, saxs_mean])
np.savetxt('saxs_sim.dat', new_arr.T)
