# Numerical integration code.
# Author: Nathan Malamud

import numpy as np

from . import ornbeck
from .model import ode_model

from .config.general import RUNTIME, TEMPERATURE, MAX_ITERATIONS, E_MIN, E_MAX, DELTA_G_BOUND

def calculate_gibbs(C, S, F, T):
    """
    Calculate Gibbs free energy vector ∆G.
    """
    R = 8.3e-3 # universal gas constant
    molar_C = C * 10e-6

    # This helps me avoid Python's runtime warnings.
    # This doesn't affect the simulation results
    # since I have an adaptive step size scheme in execute
    # which checks for negative metabolite concentrations.
    if (molar_C <= 0).any():
        molar_C[molar_C <= 0] = 1e-37

    return S.transpose() @ (F + R * T * np.log(molar_C))

def calculate_flux(t, h, C, f) -> tuple:
    """
    Uses a 5th order Runge-Kutta-Fehlberg method to determine
    the net metabolite flux at time t. Also returns the error bound.

    Source: [Cheney & Kincaid (2008) - Numerical Mathematics and Computing Sixth Edition pp. 450-453]
    """
    k1 = h * f(t, C)
    k2 = h * f(t + 0.25 * h, C + 0.25 * k1)
    k3 = h * f(t + 0.375 * h, C + 0.09375 * k1 + 0.28125 * k2)
    k4 = h * f(t + 12./13. * h, C + 1932./2197. * k1 + (-7200./2197.) * k2 + 7296./2197 * k3) 
    k5 = h * f(t + h, C + 439./216. * k1 + (-8.) * k2 + (3680./513.) * k3 + (-845./4104.) * k4)
    k6 = h * f(t + 0.5 * h, C + (-8./27.) * k1 + 2. * k2 + (-3544./2565.) * k3 + (-845./4104.) * k4 + (-0.275) * k5)

    C_RK4_flux = (25./216.) * k1 + (1408./2565.) * k3 + (2197./4104.) * k4 + (-0.2) * k5
    C_RK5_flux = (16./135.) * k1 + (6656./12825.) * k3 + (28561./56430.) * k4 + (-0.18) * k5 + (2./55.) * k6

    # Error is determined as the euclidean distance between estimates.
    error = np.linalg.norm(C_RK5_flux - C_RK4_flux)

    return C_RK5_flux, error

def execute(initialC, deltaGf0, stoich_mats, ou_parameters, default_timestep=0.01, random_seed=None) -> object:
    """ 
    Executes a simulation to
    to find the Dead-End of a given ODE
    under the conditions given.

    Returns a dictionary (sol) containing all calculations
    for all variables of interest.

    input:
    - - - - - - - -
    initialC - starting concentrations for all metabolites.
    deltaGf0 - vector of free energies.
    stoich_mats - three stoichiometric matrices from sim_config file.
    ou_parameters - list of parameter vectors for reaction kinetics

    output:
    - - - - - - - -
    sol - dictionary containing time-series data for multiple variables.
    """

    stoich_mat_full, stoich_mat_lim, stoich_mat_nconst = stoich_mats
    M, N = stoich_mat_full.shape

    sol = {'time' : [],
           'composition' : [],
           'deltaG' : [],
           'stochastic process value' : [],
           'net metabolite flux' : [],
           'messages' : []}
    
    composition = initialC.copy()

    time = 0

    # Initiate the Ornstein-Uhlenbeck process
    weeks = int(RUNTIME // 7)
    ornbeck_times = np.linspace(0, (weeks + 1) * 7.0, weeks + 1)

    typical_rates, typical_decay, typical_std = ou_parameters
    [means, sigmas, decays, starts] = ornbeck.calibrate(stoich_mat_lim, typical_rates, typical_decay, typical_std, typical_con=1.0, random_seed=random_seed)
    ornbeck_spline = ornbeck.spline(ornbeck_times, means, sigmas, decays, starts)

    # Boolean flag: indicates whether a dead-end state has been reached
    success = False

    # Set timestep to default
    timestep = default_timestep

    ode_function = lambda t, C : ode_model(t, C, stoich_mats, calculate_gibbs(C, stoich_mat_full, deltaGf0, TEMPERATURE), ornbeck_vector)

    sol['messages'].append('# Time (in days) : message.')
    sol['messages'].append(f'{time:.4f}: simulation starts with timestep {timestep}.')

    # Counts how many times we have changed the timestep without moving forward
    loop_count = 0

    # Counts total iterations
    iter = 0

    while time <= RUNTIME:

        iter += 1

        deltaG = calculate_gibbs(composition, stoich_mat_full, deltaGf0, TEMPERATURE)
        ornbeck_vector = ornbeck_spline(time)

        # Fehlberg scheme
        flux, error = calculate_flux(time, timestep, composition, ode_function)

        # Adaptive timestep
        new_composition = composition + flux

        if iter > MAX_ITERATIONS:
            sol['messages'].append(f'{time:.4f}: simulation terminated after {iter} iterations.')
            break

        # This should never happen
        if np.isnan(new_composition).any():
            timestep /= 2.0
            sol['messages'].append(f'{time:.4f}: timestep halved to {timestep} due to NaN concentrations.')
            continue

        if (new_composition <= 0).any():
            timestep /= 2.0
            sol['messages'].append(f'{time:.4f}: timestep halved to {timestep} due to negative concentrations.')
            continue

        # If this occurs, we are stuck in a loop
        # where the simulation is halving and doubling
        # the time step over and over again.
        #
        # To remedy this, we simply move forward without adjusting the timestep.
        if loop_count > 5:
            sol['messages'].append(f'{time:.4f} : excessive loop count detected. Continuing with timestep {timestep}.')

        # These error bounds can be interpreted as:
        # E_MIN to E_MAX uM uncertainty in predictions.
        elif error > E_MAX:
            timestep /= 2.0
            loop_count += 1
            sol['messages'].append(f'{time:.4f} : timestep halved to {timestep} due to error over {E_MAX}.')
            continue
        elif (error < E_MIN and (new_composition > 1.0).all()):
            timestep *= 2.0
            loop_count += 1
            sol['messages'].append(f'{time:.4f} : timestep doubled to {timestep} due to error under {E_MIN}.')
            continue

        # Abort - avoids an infinite loop
        if timestep == 0:
            sol['messages'].append(f'{time:.4f}: simulation terminated early with timestep {timestep}.')
            break

        # Record current values
        sol['time'].append(time)
        sol['composition'].append(composition)
        sol['deltaG'].append(deltaG)
        sol['stochastic process value'].append(ornbeck_vector)
        sol['net metabolite flux'].append(flux)

        # Dead-end state condition:
        if (deltaG >= DELTA_G_BOUND).all():
            sol['messages'].append(f'{time:.4f}: simulation terminated successfully with timestep {timestep}.')
            success = True
            break

        composition = new_composition
        time += timestep
        loop_count = 0
    
    if success == False:
        sol['messages'].append(f'{time:.4f}: simulation terminated unsuccessfully with timestep {timestep}.')

    sol['messages'].append(f'Simulation terminated after {iter} loop iterations (see sim.execute function).')

    for data in sol:
        sol[data] = np.array(sol[data])
    
    return sol, success
