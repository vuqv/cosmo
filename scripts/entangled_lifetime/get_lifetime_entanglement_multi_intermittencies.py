#!/usr/bin/env python
# coding: utf-8

import sys

import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.lib.correlations import correct_intermittency, autocorrelation
from scipy.optimize import curve_fit

# setup font for plotting
font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}

mtl.rc('font', **font)

if len(sys.argv) != 2:
    print(f'[0] = script')
    print(f'[1] = Input data (entanglement analysis)')
    quit()

result_file = sys.argv[1]
logfile = result_file.split('.')[0] + '.log'
raw_img_file = result_file.split('.')[0] + '_raw.png'
# system specific

data = np.loadtxt(result_file, usecols=5)
threshold = 0.6
timestep = 10  # ps
tau_max = 500  # maximum number of frames to calculate autocorrelation


################################################
def map_binary(x):
    """
        This function return the binary state of frame. If current frame is entangled, return {1}
        Otherwise, return an empty list then length of it is 0 {}
    """
    if x:
        return set([1])
    else:
        return set([])


def map_binary_unentanled(x):
    """
        Do every that inverse to map_binary function.
        The usage of this function is for finding short disentanglement time
    """
    if x:
        return set([])
    else:
        return set([1])


def fit_biexponential(tau_timeseries, ac_timeseries, initial_guess=[1, 0.5, 1, 2]):
    """Fit a biexponential function to a hydrogen bond time autocorrelation function

    Return the two time constants
    """
    print(f"initial guess: {initial_guess}")

    def model(t, A, tau1, B, tau2):
        """Fit data to a biexponential function.
        """
        return A * np.exp(-t / tau1) + B * np.exp(-t / tau2)

    params, params_covariance = curve_fit(model, tau_timeseries, ac_timeseries, initial_guess,
                                          bounds=(0, [1, np.inf, 1, np.inf]))

    fit_t = np.linspace(tau_timeseries[0], tau_timeseries[-1], 1000)
    fit_ac = model(fit_t, *params)

    return params, fit_t, fit_ac


# Main program

f = open(logfile, 'w')

states = [map_binary(x >= threshold) for x in data]

# use intermittency that allow entangled break for lagframe and form again.
fig = plt.figure(figsize=(16, 9))
ax = plt.axes()
params = [1, 0.5, 1, 2]
intermittencies = [0, 1, 5, 10, 30, 50, 100]
for intermittency in intermittencies:
    print(f"intermittency: {intermittency}")
    state_intermittency = correct_intermittency(states, intermittency=intermittency)
    tau_frames, timeseries, timeseries_data = autocorrelation(state_intermittency, tau_max)

    print(f"Data with lag-frames (using intermittency= {intermittency} frames)", file=f)
    params, fit_t, fit_ac = fit_biexponential(tau_frames, timeseries, initial_guess=params)
    A, tau1, B, tau2 = params
    # this is equivalent to integrate.
    time_constant = A * tau1 + B * tau2
    print(
        f"corrected_intermittency data with intermittency= {intermittency} : time_constant = {time_constant * timestep:.2f} (ps)",
        file=f)
    print(
        f"A = {params[0]:.2f}, tau1 = {params[1] * timestep:.2f} (ps), B = {params[2]:.2f}, tau2 = {params[3] * timestep:.2f} (ps)",
        file=f)

    # plot for data with intermittency
    tau_times = np.array(tau_frames) * timestep
    timeseries = np.array(timeseries)

    ax.plot(tau_times, timeseries, lw=2, label=intermittency)
    plt.title(f"Entanglement lifetime ", weight="bold")
    plt.xlabel(r"$\tau\ \rm (ps)$")
    plt.ylabel(r"$C(\tau)$")
    plt.ylim(0, 1)
    ax.legend()
    # plt.show()

intermit_img_file = result_file.split('.')[0] + '_intermittencies.png'
fig.savefig(intermit_img_file, dpi=600)
f.close()
