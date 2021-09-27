
# Main simulation script - acts as a command-line interface to the sim package.
# See README.md for usage examples.
#
# Author: Nathan Malamud
#

import sim
import numpy as np

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
# I. READ COMMAND-LINE ARGUMENTS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

import argparse

parser = argparse.ArgumentParser(description='Microbe Metabolism Simulator')

parser.add_argument('-m', '--mode', help='Determines whether simulation mode is \'random\' or \'fixed\'',
                const='random', default='random', nargs='?', type=str)

parser.add_argument('-r', '--runs', help='Set number of simulation runs.',
                    const=50, default=50, nargs='?', type=int)

parser.add_argument('-o', '--out', help='Specify target output directory.',
                    const='data', default='data', nargs='?', type=str)

parser.add_argument('-t', '--timestep', help="Default timestep for ODE solver.",
                    const=5.0, default=5.0, nargs='?', type=float)

parser.add_argument('-d', '--debug', help="Indicates whether to log all sim data (True) or just end states (False).",
                    const=False, default=False, nargs='?', type=bool)

parser.add_argument('-s', '--seed', help="Input random seed for Ornstein-Uhlenbeck Process.",
                    const=None, default=None, nargs='?', type=int)

args = parser.parse_args()

reaction_argument = None # depending on the mode, this is either an integer or a list

if args.mode == 'random':
    network_size = 3
    while True:
        network_size = input("\nCurrently running simulation in 'random' mode.\nPlease specify size of random network: ")
        try:
            network_size = int(network_size)
        except:
            print("\nNot valid input for size of random network. Please try again.")
            continue
        break

    reaction_argument = network_size

elif args.mode == 'fixed':
    network_list = []

    network_list = input("\nCurrently running simulation in 'fixed' mode.\nPlease specify list of reactions (separated by spaces): ")
    network_list = network_list.strip().split(' ')

    reaction_argument = network_list

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
# II. READ SIULATION PARAMETERS from SETUP.py
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
 
[
    metabolites,
    reactions,
    deltaGf0,
    initialC,
    stoich_mats,
    ou_parameters
] = sim.setup.build_network(reaction_argument)

stoich_mat_full, stoich_mat_lim, stoich_mat_nconst = stoich_mats

# - - - - - - - - - - - - - - - - - - - //
# III. BUILD OUTPUT DIRECTORY 
# - - - - - - - - - - - - - - - - - - - //

import os, shutil, csv

if os.path.exists(args.out):
    print("""\nWarning: specified output directory already exists.\nOutput directory will be overwritten.""")
    shutil.rmtree(args.out)

os.makedirs(args.out)

with open(f'{args.out}/dead_ends.tsv', 'a+') as dead_end_file:
    dead_ends = csv.writer(dead_end_file, delimiter='\t')
    dead_ends.writerow(list(metabolites))

    # log all reactions and metabolites
    with open(f'{args.out}/network_desc.txt', 'w') as names:
        names.write(f'metabolites ({len(metabolites)}): {metabolites}\n')
        names.write(f'reactions ({len(reactions)}): {reactions}\n')

    # write all stoichiometric matrices to txt files
    np.savetxt(f'{args.out}/stoich_mat_full.txt', stoich_mat_full)
    np.savetxt(f'{args.out}/stoich_mat_lim.txt', stoich_mat_lim)
    np.savetxt(f'{args.out}/stoich_mat_nconst.txt', stoich_mat_nconst)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    # IV. DETERMINE WHETHER TO VARY METABOLITE CONCENTRATIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

    # if mode == 'fixed', specify variable to be changed and what interval to change it on.
    VARY_METABOLITE = False 
    if args.mode == 'fixed':
        answer = input("Would you like to vary the concentration of a certain metabolite? (y or n) : ")
        if answer.lower() in ['yes', 'y', 'sure']:
            while True:
                metabolite_of_interest = input(f"""
                    Please select from the following list of metabolites:

                    {"n".join(metabolites)}

                    Your choice: """)
                
                try:
                    assert metabolite_of_interest.strip() in metabolites
                except AssertionError:
                    print("Invalid choice: {metabolite_of_interest}")
                    continue

                break

            while True:
                input_range = input("\nPlease input concentration range as a tuple of the form (min, max): ")

                try:
                    input_range = input_range.strip()
                    assert input_range[0] == '(' and input_range[-1] == ')'
                    integer_tuple = input_range[1:-1].split(',')

                    try:
                        integer_tuple = tuple([float(i) for i in integer_tuple])
                    except:
                        raise AssertionError
                
                    assert len(input_range) == 2
                    assert input_range[0] > 0 and input_range[1] > 0
                    assert input_range[1] > input_range[0]

                except AssertionError:
                    print(f"""Invalid input: {input_range}. The rules: max > min, and must be a tuple of two strictly positive floats.
                         \nTry again, please.""") 
                    continue

                break

    if VARY_METABOLITE:
        var_met_index = metabolites.find(metabolite_of_interest)

        a, b = input_range
        var_range = np.linspace(a, b, num=args.runs)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    # V. EXECUTE ALL SIMULATIONS (In Serial)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    for i in range(1, args.runs + 1):
        del sim
        import sim

        if VARY_METABOLITE:
            initialC[var_met_index] = var_range[i]

        sol, success = sim.execute(initialC, deltaGf0, stoich_mats, ou_parameters, default_timestep=args.timestep, random_seed=args.seed)

        # if dead-end state condition met
        if success:
            dead_ends.writerow(sol['composition'][-1, :])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
        # VI. OUTPUT DATA to TSV FILES if DEBUG == TRUE
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
        if args.debug:
            os.mkdir(f'{args.out}/sim_{i:0>2}')

            metabolite_header = '\t'.join(['time'] + list(metabolites))
            reaction_header = '\t'.join(['time'] + list(reactions))

            deltaG_data = np.column_stack([sol['time'], sol['deltaG']])
            np.savetxt(f'{args.out}/sim_{i:0>2}/deltaG_{i:0>2}.tsv', deltaG_data, delimiter='\t',
                    header=reaction_header, comments='')

            composition_data = np.column_stack([sol['time'], sol['composition']])
            np.savetxt(f'{args.out}/sim_{i:0>2}/composition_{i:0>2}.tsv', composition_data, delimiter='\t',
                    header=metabolite_header, comments='')

            ornbeck_data = np.column_stack([sol['time'], sol['stochastic process value']])
            np.savetxt(f'{args.out}/sim_{i:0>2}/ornbeck_{i:0>2}.tsv', ornbeck_data, delimiter='\t',
                    header=reaction_header, comments='')

            flux_data = np.column_stack([sol['time'], sol['net metabolite flux']])
            np.savetxt(f'{args.out}/sim_{i:0>2}/flux_{i:0>2}.tsv', flux_data, delimiter='\t',
                    header=metabolite_header, comments='')
        
            # report whether simulation was successful
            with open(f'{args.out}/sim_{i:0>2}/report.txt', 'w') as report:
                    report.write(f'success: {success}')

            # write all error messages to a txt file
            with open(f'{args.out}/sim_{i:0>2}/messages.txt', 'w') as messages:
                    for line in sol['messages']:
                            messages.write(line + '\n')
