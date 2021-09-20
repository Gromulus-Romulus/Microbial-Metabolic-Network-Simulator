"""
Based on info provided in reactions.py, metabolites.py,
and stoichiometry.txt, setup.py selects K reactions
and builds the required data structures for a reaction-centric model.
"""

import numpy as np
import numpy.random as np_random

import sys

from . import parser

from .config.reactions import *     # reactions, typical_rates, typical_decay, typical_std
from .config.metabolites import *   # metabolites, initialC, deltaGf0

def build_network(K: int) -> list:
    """
    Via a two-phase process, this function uses
    all of the information provided in the config folder
    to return all required parameters for the simulation.

    To whom it may concern:

        You know what would be pretty cool? If we had the option
        to have K be a list of reactions as opposed to a single integer.

        Just to give this function more flexibility and to give the user
        more options for simulation configs.

        Alas, I do not possess the time to implement this right now.
        I leave it up to you! The far more capable coder!

        - Nathan

    Only argument is K - number of reactions needed in the network.

    Order of returns:
        metabolites - names of all species
        reactions - names of all K reactions
        deltaGf0 - gibbs free energy values
        initalC - concentration values
        stoich_mats - stoichiometry of the system [full, lim, nconst]
        ou_parameters - needed for kinetic rates [means, sigmas, decays, starts]
    """

    global metabolites
    global reactions

    global typical_decay, typical_rates, typical_std
    global deltaGf0, initialC

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    # I. NETWORK MODEL SELECTION and TRIMMING:
    #
    # In order to generate random networks, we need
    # to be able to access a master library (e.g. in reaction config)
    # and filter away what is not needed.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    M = len(metabolites)
    N = len(reactions)

    # select K random reactions, and extract all parameters needed for them
    selection = np.zeros(shape=N, dtype='bool')

    # Choose a random assortment of reaction indices (without replacement - no repeats)
    selection_indices = np_random.choice(np.array([i for i in range(0, N)]), size=K, replace=False)

    for i in selection_indices:
        selection[i] = True

    reactions = reactions[selection].copy()
    typical_rates = typical_rates[selection].copy()
    typical_decay = typical_decay[selection].copy()
    typical_std = typical_std[selection].copy()

    stoichiometry = parser.parse_file('./sim/config/stoichiometry.txt')

    all_reactions = list(stoichiometry.keys())

    for rxn in all_reactions:
        # that is, if the current reaction is not in our selection
        if not rxn in reactions:
            del stoichiometry[rxn]

    # select all metabolites needed for these reactions

    metabolites_needed = []
    for rxn in stoichiometry:
        for met in stoichiometry[rxn]['species']:
            metabolites_needed.append(met)

    selection = []

    for met in metabolites:
        if met in metabolites_needed:
            selection.append(True)
        else:
            selection.append(False)

    selection = np.array(selection)

    metabolites = metabolites[selection].copy()
    deltaGf0 = deltaGf0[selection].copy()
    initialC = initialC[selection].copy()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    # II. BUILDING STOICHIOMETRIC MATRICES
    #
    # After filtering through the configuration library,
    # we need to use the data provided in stoichiometry.txt
    # to build three stoichiometric matices.
    #
    #   * stoich_mat_full   - full stoichiometry used for calculating deltaG values
    #   * stoich_mat_lim    - exclusively "limiting" stoichiometry for kinetic rates
    #   * stoich_mat_nconst - non-constant stoichiometry used for updating concentrations
    #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

    M = len(metabolites)
    N = len(reactions)

    stoich_mat_full = np.zeros(shape=(M, N), dtype=np.double)
    stoich_mat_lim = np.zeros(shape=(M, N), dtype=np.double)
    stoich_mat_nconst = np.zeros(shape=(M, N), dtype=np.double)

    for i in range(N):
        rxn = reactions[i]

        try:
            stoich_info = stoichiometry[rxn]['species']
        except:
            print(f"\nERROR! reaction '{rxn}' not found in stoichiometry.txt file. Aborting!")
            sys.exit(1)

        for j in range(M):
            met = metabolites[j]

            if met in stoich_info:

                coeff = stoich_info[met]['coeff']
                flags = stoich_info[met]['flags']
                stoich_mat_full[j, i] = coeff

                # check for non-limiting and constant concentration flags
                if not ('CONST' in flags or 'C' in flags):
                    stoich_mat_nconst[j, i] = coeff

                if not 'NL' in flags:
                    stoich_mat_lim[j, i] = coeff
                else:
                    stoich_mat_lim[j, i] = 0.0
                    # non-limiting substrates must also be constant by default
                    if coeff < 0:
                        stoich_mat_nconst[j, i] = 0.0

    stoich_mats = [stoich_mat_full, stoich_mat_lim, stoich_mat_nconst]
    ou_parameters = [typical_rates, typical_decay, typical_std]

    return [metabolites, reactions, deltaGf0, initialC, stoich_mats, ou_parameters]