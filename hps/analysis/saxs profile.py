#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
from math import sin

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# from numba import njit


# In[2]:


u = mda.Universe('asyn.psf', 'asyn_equil.dcd')
# u=mda.Universe('asyn.psf', 'asyn_equil_final.pdb')
protein = u.select_atoms('protein')
n_atoms = len(protein)
resnames = protein.resnames

# In[3]:


n_frames = u.trajectory.n_frames
print(n_frames)
print(n_atoms)

# In[4]:


aa_order = ['q', 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN',
            'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
form_factor_data = np.loadtxt('table_edm_CA_1bead.dat')
# form_factor_data = np.loadtxt('table_edm_COE_1bead.dat')
form_factor = pd.DataFrame(form_factor_data, columns=aa_order)


# Note that the form-factor are given with the wavelength in unit of A^-1


# In[5]:


# @njit(fastmath=True)
def saxs_single_frame(positions):
    n_atoms = len(positions)
    pair_array = np.asarray(list(itertools.product(positions, positions)))
    distance_diff = pair_array[:, 0, :] - pair_array[:, 1, :]
    distance_pair = np.linalg.norm(distance_diff, axis=1).reshape(n_atoms, n_atoms)

    qs = np.arange(0, 0.501, 0.01)
    Iq = np.zeros(len(qs))
    for idx, q in enumerate(qs):
        I = 0.0
        # I0
        if idx == 0:

            for i in range(n_atoms):
                for j in range(n_atoms):
                    if i != j:
                        I += form_factor.loc[0, resnames[i]] * form_factor.loc[0, resnames[j]]
                        # Iq, otherwise
        else:
            for i in range(n_atoms):
                for j in range(n_atoms):
                    if i != j:
                        I += form_factor.loc[1, resnames[i]] * form_factor.loc[1, resnames[j]] * sin(
                            q * distance_pair[i, j]) / (q * distance_pair[i, j])
        Iq[idx] = I

    return Iq / Iq[0]


# In[6]:


# calculate SAXS along trajectory
saxs_total = np.empty((51, 0), dtype=float)
for ts in u.trajectory[:50]:
    positions = protein.positions
    res = saxs_single_frame(positions)
    saxs_total = np.append(saxs_total, np.array([res]).transpose(), axis=1)

# In[9]:


# get mean value 
saxs_mean = np.mean(saxs_total, axis=1)

# In[11]:


# Lindorff Data
lindof_dat = np.loadtxt('async.dat')

# In[12]:


qs = np.arange(0, 0.501, 0.01)
plt.figure(figsize=(8, 6), dpi=80)
plt.plot(qs, saxs_mean, '-')
plt.plot(lindof_dat[:, 0], lindof_dat[:, 1] / lindof_dat[0, 1])
plt.yscale('log')
plt.xlabel('q(A^-1)')
plt.ylabel('Iq/I0')
