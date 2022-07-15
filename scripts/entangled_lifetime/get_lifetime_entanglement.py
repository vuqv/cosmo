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
tau_frames, timeseries, timeseries_data = autocorrelation(states, tau_max)

# raw data, we don't know about the optimal values so just use default [1, 0.5, 1, 2]
params, fit_t, fit_ac = fit_biexponential(tau_frames, timeseries)
A, tau1, B, tau2 = params
# this is equivalent to integrate.
time_constant = A * tau1 + B * tau2
print(f"Raw data: time_constant = {time_constant * timestep:.2f} (ps)", file=f)
print(
    f"A = {params[0]:.2f}, tau1 = {params[1] * timestep:.2f} (ps), B = {params[2]:.2f}, tau2 = {params[3] * timestep:.2f} (ps)",
    file=f)
allow_lag = min(int(params[1]) + 1, int(params[3]) + 1)
print(f"Lagframe to use to get rid of fast event: {allow_lag} frames = {allow_lag * timestep} (ps)", file=f)

# plot for raw data
tau_times = np.array(tau_frames) * timestep  # ps
timeseries = np.array(timeseries)
fig = plt.figure(figsize=(16, 9))
ax = plt.axes()
ax.plot(tau_times, timeseries, lw=2)
ax.plot(fit_t * timestep, fit_ac, label=(
        f"lifetime={time_constant * timestep:.2f} (ps)" + "\n" + r"Fitting function:" +
        "\n" +
        r"$C(\tau)=A\exp(-t/\tau_1)+B\exp(-t/\tau_2)$" +
        "\n" +
        r"$A = %.2f, \tau_1 = %.2f(ps), B = %.2f, \tau_2 = %.2f(ps)$" % (
            params[0], params[1] * timestep, params[2], params[3] * timestep)))
plt.title(r"Entanglement lifetime (raw data)", weight="bold")
plt.xlabel(r"$\tau\ \rm (ps)$")
plt.ylabel(r"$C(\tau)$")
ax.legend()
# plt.show()
fig.savefig(raw_img_file, dpi=600)
# End of raw data

# use intermittency that allow entangled break for lagframe and form again.

state_intermittency = correct_intermittency(states, intermittency=allow_lag)
tau_frames, timeseries, timeseries_data = autocorrelation(state_intermittency, tau_max)

print(f"Data with lag-frames (using intermittency= {allow_lag} frames)", file=f)
params, fit_t, fit_ac = fit_biexponential(tau_frames, timeseries, initial_guess=params)
A, tau1, B, tau2 = params
# this is equivalent to integrate.
time_constant = A * tau1 + B * tau2
print(f"corrected_intermittency data: time_constant = {time_constant * timestep:.2f} (ps)", file=f)
print(
    f"A = {params[0]:.2f}, tau1 = {params[1] * timestep:.2f} (ps), B = {params[2]:.2f}, tau2 = {params[3] * timestep:.2f} (ps)",
    file=f)

# plot for data with intermittency
tau_times = np.array(tau_frames) * timestep
timeseries = np.array(timeseries)

fig = plt.figure(figsize=(16, 9))
ax = plt.axes()
ax.plot(tau_times, timeseries, lw=2, label='data')
ax.plot(fit_t * timestep, fit_ac, label=(
        f"lifetime = {time_constant * timestep:.2f} (ps)" + "\n" + r"Fitting function:" +
        "\n" +
        r"$C(\tau)=A\exp(-t/\tau_1)+B\exp(-t/\tau_2)$" +
        "\n" +
        r"$A = %.2f, \tau_1 = %.2f(ps), B = %.2f, \tau_2 = %.2f(ps)$" % (
            params[0], params[1] * timestep, params[2], params[3] * timestep)))
plt.title(f"Entanglement lifetime (intermittency = {allow_lag} frames)", weight="bold")
plt.xlabel(r"$\tau\ \rm (ps)$")
plt.ylabel(r"$C(\tau)$")
plt.ylim(0, 1)
ax.legend()
# plt.show()
intermit_img_file = result_file.split('.')[0] + '_intermittency_' + f'{allow_lag}_frames.png'
fig.savefig(intermit_img_file, dpi=600)
f.close()
