# Main simulation script - acts as a command-line interface to the sim package.
# See README.md for usage examples.
#
# Author: Nathan Malamud

import sim
import numpy as np

import sys
import os, shutil
import csv

import argparse
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
# I. READ COMMAND-LINE ARGUMENTS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

parser = argparse.ArgumentParser(description='Microbe Metabolism Simulator')

parser.add_argument('-r', '--runs', help='Set number of simulation runs.',
                    const=50, default=50, nargs='?', type=int)

parser.add_argument('-o', '--out', help='Specify target output directory.',
                    const='data', default='data', nargs='?', type=str)

parser.add_argument('-t', '--timestep', help="Default timestep for ODE solver.",
                    const=5.0, default=5.0, nargs='?', type=float)

parser.add_argument('-d', '--debug', help="Indicates whether to log all sim data (True) or just end states (False).",
                    const=True, default=False, nargs='?', type=bool)

parser.add_argument('-s', '--seed', help="Input random seed for Ornstein-Uhlenbeck Process.",
                    const=None, default=None, nargs='?', type=int)

args = parser.parse_args()

RUNS = args.runs
OUT = args.out
TIMESTEP = args.timestep
DEBUG = args.debug
SEED = args.seed

if sys.stdin is None:
    sys.stderr.write("Cannot execute program without input file.")
    sys.stderr.flush()
    sys.exit(1)

FILE = sys.stdin

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
# II. READ SIMULATION PARAMETERS from SETUP.py
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

reaction_list = FILE.readline().strip().split(' ')
 
[
    metabolites,
    reactions,
    deltaGf0,
    initialC,
    stoich_mats,
    ou_parameters
] = sim.setup.build_network(reaction_list)

stoich_mat_full, stoich_mat_lim, stoich_mat_nconst = stoich_mats

# - - - - - - - - - - - - - - - - - - - //
# III. BUILD OUTPUT DIRECTORY 
# - - - - - - - - - - - - - - - - - - - //

if os.path.exists(OUT):
    print("""\nWarning: specified output directory already exists.\nOutput directory will be overwritten.""")
    shutil.rmtree(OUT)

os.makedirs(OUT)

with open(f'{OUT}/dead_ends.tsv', 'a+') as dead_end_file:
    # log all reactions and metabolites
    with open(f'{OUT}/network_desc.txt', 'w') as names:
        names.write(f'metabolites ({len(metabolites)}): {metabolites}\n')
        names.write(f'reactions ({len(reactions)}): {reactions}\n')

    # write all stoichiometric matrices to txt files
    np.savetxt(f'{OUT}/stoich_mat_full.txt', stoich_mat_full)
    np.savetxt(f'{OUT}/stoich_mat_lim.txt', stoich_mat_lim)
    np.savetxt(f'{OUT}/stoich_mat_nconst.txt', stoich_mat_nconst)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    # IV. DETERMINE WHETHER TO VARY METABOLITE CONCENTRATIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

    VARY_METABOLITE = False 
    
    answer = FILE.readline().strip()
    if answer.lower() in ['yes', 'y', 'sure']:
        var_met_name = FILE.readline().strip()

        try:
            assert var_met_name in metabolites
        except AssertionError:
            sys.stderr.write("Invalid choice: {var_met_name} not found.")
            sys.stderr.flush()
            sys.exit(1)

        try:
            input_range = FILE.readline().strip()
            assert input_range[0] == '(' and input_range[-1] == ')'
            float_tuple = input_range[1:-1].split(',')

            try:
                float_tuple = tuple([float(i) for i in float_tuple])
            except:
                raise AssertionError
                
                assert len(float_tuple) == 2
                assert float_tuple[0] > 0 and float_tuple[1] > 0
                assert float_tuple[1] > float_tuple[0]

        except AssertionError:
            sys.stderr.write(f"Invalid input: {input_range}. The rules: max > min, and must be a tuple of two strictly positive floats.")
            sys.stderr.flush()
            sys.exit(1)

        input_range = float_tuple
        VARY_METABOLITE = True

    if VARY_METABOLITE:
        var_met_index = list(metabolites).index(var_met_name)

        a, b = input_range
        var_range = np.linspace(a, b, num=RUNS)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    # V. EXECUTE ALL SIMULATIONS (In Serial)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    dead_ends = csv.writer(dead_end_file, delimiter='\t')

    if VARY_METABOLITE:
        dead_ends.writerow(list(metabolites) + [f'{var_met_name}_INIT'])
    else:
        dead_ends.writerow(list(metabolites))

    for i in range(1, RUNS + 1):

        # refresh
        del sim
        import sim

        if VARY_METABOLITE:
            # variable metabolite initial concentration
            var_met_init_con = var_range[i - 1]
            initialC[var_met_index] = var_met_init_con

        sol, success = sim.execute(initialC, deltaGf0, stoich_mats, ou_parameters, default_timestep=TIMESTEP, random_seed=SEED)

        # if dead-end state condition met
        if success:
            if VARY_METABOLITE:
                dead_ends.writerow(np.append(sol['composition'][-1, :], var_met_init_con))
            else:
                dead_ends.writerow(sol['composition'][-1, :])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
        # VI. OUTPUT DATA TO TSV FILES IF DEBUG == TRUE
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
        if DEBUG:
            os.mkdir(f'{OUT}/sim_{i:0>2}')

            if VARY_METABOLITE:
                with open(f'{OUT}/sim_{i:0>2}/initial_condition.txt', 'w') as f:
                    f.write(f'[{var_met_name}]  # in micro-molars\n{var_range[i - 1]}')

            metabolite_header = '\t'.join(['time'] + list(metabolites))
            reaction_header = '\t'.join(['time'] + list(reactions))

            deltaG_data = np.column_stack([sol['time'], sol['deltaG']])
            np.savetxt(f'{OUT}/sim_{i:0>2}/deltaG_{i:0>2}.tsv', deltaG_data, delimiter='\t',
                    header=reaction_header, comments='')

            composition_data = np.column_stack([sol['time'], sol['composition']])
            np.savetxt(f'{OUT}/sim_{i:0>2}/composition_{i:0>2}.tsv', composition_data, delimiter='\t',
                    header=metabolite_header, comments='')

            ornbeck_data = np.column_stack([sol['time'], sol['stochastic process value']])
            np.savetxt(f'{OUT}/sim_{i:0>2}/ornbeck_{i:0>2}.tsv', ornbeck_data, delimiter='\t',
                    header=reaction_header, comments='')

            flux_data = np.column_stack([sol['time'], sol['net metabolite flux']])
            np.savetxt(f'{OUT}/sim_{i:0>2}/flux_{i:0>2}.tsv', flux_data, delimiter='\t',
                    header=metabolite_header, comments='')
        
            # report whether simulation was successful
            with open(f'{OUT}/sim_{i:0>2}/report.txt', 'w') as report:
                    report.write(f'success: {success}')

            # write all error messages to a txt file
            with open(f'{OUT}/sim_{i:0>2}/messages.txt', 'w') as messages:
                    for line in sol['messages']:
                            messages.write(line + '\n')
