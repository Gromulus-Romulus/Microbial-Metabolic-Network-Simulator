# Code for generating an Ornstein-Uhlenbeck process.
# This process will be used to model stochastic kinetic rates.
#
# Adapted from C++ code provided by Dr. Louca.
# Author: Nathan Malamud
#

import numpy as np
from scipy import interpolate

def spline(times, means, sigmas, decays, starts, dt=7.0) -> callable:
    """ Interpolates Ornstein-Uhlenbeck nodes. """
    assert len(means) == len(sigmas) == len(decays) == len(starts), "Parameter arrays must be the same size!"

    N = len(means)
    T = len(times)
    
    time_series = np.ndarray(shape=(N, T), dtype=np.double)

    for n in range(N):
        time_series[n] = knots(times, means[n], sigmas[n], decays[n], starts[n], dt)

    vector_spline = interpolate.interp1d(times, time_series, kind='linear')

    return vector_spline

def knots(times, mu, sigma, decay, start, dt) -> np.array:
    """ Generated by Ornstein-Uhlenbeck process """

    NT = len(times)

    knots = np.zeros(NT, dtype=np.double)
    knots[0] = start

    for t in range(NT - 1):
        previous = knots[t]

        std = sigma * np.sqrt(1 - np.exp(-2 * dt * decay))
        expectation = mu + (previous - mu) * np.exp(-dt * decay)
        new_value = expectation + std * np.random.randn()

        if new_value < 0.0:
            new_value = 0.0

        knots[t + 1] = new_value

    return knots

def calibrate(stoich_mat_lim, typical_rates, typical_decay, typical_std, typical_con=1.0, random_seed=None) -> list:
    """
    Calibrates the Ornstein-Uhlenbeck process based on
    simulation parameters from the configuration files.

    The value of random_seed will determine the behavior of the Ornstein-Uhlenbeck processes.

    stoich_mat_lim - stoichiometry for limiting substrates (M x N matrix)
    typical_rates - typical range for reaction rates (N x 2 matrix)
    typical_decay - typical range for decay rates (N x 2 matrix)
    typical_std - typical range for standard deviation values (N x 2 matrix)
    typical_con - typical metabolite concentration value (set to 1.0 uM by default)
    """

    if not random_seed is None:
        np.random.seed(random_seed)

    M, N = stoich_mat_lim.shape

    means = np.zeros(N)
    sigmas = np.zeros(N)
    decays = np.zeros(N)
    starts = np.zeros(N)

    for n in range(N):
        # count number of limiting substrates
        tot_limiting_substrates = 0
        for m in range(M):
            if stoich_mat_lim[m, n] < 0:
                tot_limiting_substrates += 1

        # formulas provided by Dr. Louca
        means[n] = np.random.uniform(typical_rates[n][0], typical_rates[n][1]) / (typical_con ** tot_limiting_substrates)
        sigmas[n] = np.random.uniform(typical_std[n][0]*means[n], typical_std[n][1]*means[n]) 
        decays[n] = abs(pow(10, np.random.uniform(np.log10(typical_decay[n][0]), np.log10(typical_decay[n][1]))))
        starts[n] = abs(np.random.normal(means[n], sigmas[n]))

    return [means, sigmas, decays, starts]
