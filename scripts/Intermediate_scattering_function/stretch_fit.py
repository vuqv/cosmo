"""
Extrapolate unfolding time at 310K using data at 500, 550, 600, 650 and 700K
This script fitting with no restriction on a,b,c parameters
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

"""
Matplotlib control parameters:
"""
######################################
matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['xtick.labelsize'] = 'large'
matplotlib.rcParams['ytick.labelsize'] = 'large'
matplotlib.rcParams['legend.fontsize'] = 'large'
plt.rcParams.update({'font.size': 16})
##############################################################

"""
Functions definition
"""


def stretch_func(t, B):
    """
    Function to fitting Su, get k and t0
    """
    return np.exp(-(t / 7) ** B)


##############################################################
# first col is inverse temperature
# second row is real data, and following rows are bootstrapping data
data = np.loadtxt('result.dat')
t = data[:, 0]
# T_range = np.arange(T.min(), T.max(), 1)
fkt = data[:, 1]
# T_range = np.arange(100, 1000, 1)

kopt, kcov = curve_fit(stretch_func, t, fkt, maxfev=5000)
cc = np.corrcoef([fkt, stretch_func(t, *kopt)])[0, 1] ** 2

fig = plt.figure(figsize=(12, 8))
plt.subplots_adjust(bottom=0.15, wspace=0.3, hspace=0.3)
ax = fig.add_subplot(1, 1, 1)
ax.plot(data[:, 0], data[:, 1], 'o', label='Data points', markersize=10, alpha=0.5)

ax.plot(t, stretch_func(t, *kopt),
        label=(r'Stretch function: $F_k(t)/F_k(0) = \exp^{-(t/t_\alpha)^\beta}$' +
               '\n' + r'$R^2$' + f' = {cc:.4f}\n' +
               r'$\beta$ = ' + f'{kopt[0]:.1f}\n'))

ax.set_xlabel('t(step)')
ax.set_ylabel(r'$\frac{F_k(t)}{F_k(0)}$')
plt.title(r'$\alpha$-synuclein')
ax.legend()
plt.xscale('log')
plt.xlim(1, len(t))
fig.savefig('asyn_FIT_fkt.png', dpi=600)
