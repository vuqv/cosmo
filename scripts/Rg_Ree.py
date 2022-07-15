#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.colors as mcolors
import mdtraj as md
import numpy as np
from matplotlib import pyplot as plt


# In[9]:


def compute_Ree(traj):
    n_frames = traj.xyz.shape[0]
    Ree = np.zeros(n_frames)
    for f in range(n_frames):
        coord1 = traj.xyz[f, 0]
        coord2 = traj.xyz[f, -1]
        Ree[f] = np.linalg.norm(coord2 - coord1)

    return Ree


# In[2]:


traj = md.load('csptm_equil.dcd', top='csptm.psf')

# In[3]:


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

# In[12]:


Rg = md.compute_rg(traj[1000:], masses)
print(Rg.mean(), Rg.std())

# In[14]:


Ree = compute_Ree(traj[1000:])

# In[21]:


plt.figure(figsize=(14, 8))
plt.hist2d(Rg, Ree, bins=100, density=True)
plt.xlabel('Rg(nm)')
plt.ylabel('Ree(nm)')
plt.colorbar(norm=mcolors.Colormap('gray'))
plt.savefig('Rg_Ree.png')

# In[ ]:
