import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np

r"""
        Script to calculate mean square displacement of protein.
        Export two plots: log-log and normal scale

"""

u = mda.Universe('asyn.psf', 'asyn.dcd')
MSD = msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
MSD.run(start=10000)

font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}

mtl.rc('font', **font)

msd = MSD.results.timeseries
nframes = MSD.n_frames
timestep = 10  # this needs to be the actual time between frames
lagtimes = np.arange(nframes) * timestep  # make the lag-time axis

fig1 = plt.figure(figsize=(16, 8))
ax = plt.axes()
ax.loglog(lagtimes, msd, "-", label=r'average over all monomer')
plt.xlabel(r'$\tau$ $[ps]$')
plt.ylabel(r'$MSD$ $[\AA^2]$')
plt.legend()
plt.show()
fig1.savefig('log_log_msd_vs_time.png', dpi=600)

fig2 = plt.figure(figsize=(16, 8))
ax2 = plt.axes()
ax2.plot(lagtimes, msd, "-", label=r'average over all monomer')
plt.xlabel(r'$\tau$ $[ps]$')
plt.ylabel(r'$MSD$ $[\AA^2]$')
plt.legend()
plt.show()
fig2.savefig('normal_msd_vs_time.png', dpi=600)
