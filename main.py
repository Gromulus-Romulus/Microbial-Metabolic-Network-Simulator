# main simulation script
# TYPICAL USAGE: python3 main.py -k=[network size] --runs=[# simulation runs] --out=[output directory]

import sim
import numpy as np

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
# I. READ COMMAND-LINE ARGUMENTS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

import argparse

parser = argparse.ArgumentParser(description='Microbe Metabolism Simulator')

parser.add_argument('-k', '--reactions', help='Set number of reactions.',
                    const=5, default=5, nargs='?', type=int)

parser.add_argument('-r', '--runs', help='Set number of simulation runs.',
                    const=1, default=1, nargs='?', type=int)

parser.add_argument('-o', '--out', help='Specify target output directory.',
                    const='data', default='data', nargs='?', type=str)

parser.add_argument('-t', '--timestep', help="Default timestep for ODE solver.",
                    const=5.0, default=5.0, nargs='?', type=float)

parser.add_argument('-d', '--debug', help="Indicates whether to log all sim data (True) or just end states (False).",
                    const=False, default=False, nargs='?', type=bool)

# Note: before Nathan was playing around with random network structures,
# he simply fixed the structure and varied the concentrations of certain metabolites.
# parser.add_argument('-v', '--variable', help='Specify variable to alter.')
# parser.add_argument('-i', '--interval', help='Range of values for altered variable (-v).', nargs='+', type=float)

args = parser.parse_args()

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
] = sim.setup.build_network(args.reactions)

stoich_mat_full, stoich_mat_lim, stoich_mat_nconst = stoich_mats

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
# III. BUILD OUTPUT DIRECTORY and EXECUTE SIMULATION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

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

    for i in range(1, args.runs + 1):
        del sim
        import sim

        sol, success = sim.execute(initialC, deltaGf0, stoich_mats, ou_parameters, default_timestep=args.timestep)

        # if dead-end state condition met
        if success:
            dead_ends.writerow(sol['composition'][-1, :])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
        # IV. OUTPUT DATA to TSV FILES if DEBUG == TRUE
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
