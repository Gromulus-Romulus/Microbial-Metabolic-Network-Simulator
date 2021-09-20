"""
File parser for reaction.txt
"""

from .patterns import *

import sys

patterns = [('comment', CommentPat), ('reaction', ReactionPat)]

def parse_file(fname: str) -> dict:
    """
    Reads configuration file to extract a directory
    of all reactions used in the simulation.

    Calls match_line_to_pattern and parse_reaction.

    Side effect: opens and closes .txt file.
    """

    rxn_directory = {}

    # parse entire text file - remove comments
    with open(fname) as f:
        for line in f:
            try:
                extracted_fields = match_line_to_pattern(line)

                if extracted_fields['kind'] == 'reaction':
                    rxn_name, rxn_entry = parse_reaction(extracted_fields)

                    # make sure we haven't parsed this reaction before, either
                    already_exists = rxn_name in rxn_directory

                    if already_exists:
                        print(
                        f"""
                            \nWarning: reaction {rxn} has been declared more than once!
                            Previous coefficients, species, and names will be overwritten!
                        """)

                    rxn_directory[rxn_name] = rxn_entry

            except SyntaxError:

                 print("\nUnrecognized line in reactions.txt file!")
                 print(f">>> {line}")
                 sys.exit(1)

    return rxn_directory

def match_line_to_pattern(line: str) -> dict:

    for kind, pattern in patterns:
        match = pattern.fullmatch(line)

        if match:
            fields = match.groupdict()
            fields['kind'] = kind
            return fields

    raise SyntaxError(f"Line in configuration file not recognized! >>> {line[:-1]}")

def parse_reaction(line: str) -> tuple:
    """
    Parses one reaction in the configuration file.
    Calls parse_species.
    """

    rxn_name = line['name']
    rxn_entry = {'species': {}}

    reactants = line['reactants'].split(',')
    reactants = [species.strip() for species in reactants]

    products = line['products'].split(',')
    products = [species.strip() for species in products]

    all_species = reactants + products

    for species in all_species:
        try:
            coeff, metabolite_name, flags = parse_species(species).values()
        except SyntaxError:
            # cowardly approach: take no risks - stop program immediately
            print(f'Species {species} in reaction {rxn_name} not recognized!')
            sys.exit(1)

        # reactants are consumed, products are produced
        species_entry = {'flags': flags, 'coeff': -coeff if species in reactants else +coeff}

        # make sure the same metabolite isn't parsed twice
        if metabolite_name in rxn_entry['species']:
            print(f"\nError: metabolite {met} appears more than once in {rxn}.")
            print("Be sure to combine like terms!")
            sys.exit(1)

        rxn_entry['species'][metabolite_name] = species_entry
    
    return rxn_name, rxn_entry

def parse_species(species: str) -> dict:
    """
    Parses one species found in a line in the configuration file.

    Examples of Use:

    >> parse_species('[2]H20[CONST|NL]')
    {'coeff' : 2.0, 'name' : 'H20', 'flags' : ['CONST', 'NL']}

    >> parse_species('H20')
    {'coeff' : 1.0, 'name' : 'H20', 'flags' : []}

    >> parse_species('')
    ...SyntaxError: Metabolite pattern not recognized!

    """

    # required: must import regex patterns from patterns.py
    match = SpeciesPat.fullmatch(species)

    if match:
        fields = match.groupdict()

        coeff = fields['coeff']
        name = fields['name']
        flags = fields['flags']

        if coeff is None:
            coeff = float(1)
        else:
            coeff = float(coeff)
            if coeff < 0:
                print(f"\nWarning: negative coefficient for {species} is not recognized.")
                print(f"Absolute value of {coeff} will be used instead.")
                coeff = abs(coeff)

        if flags is None:
            flags = []
        else:
            flags = flags.split('|')

        return {'coeff': coeff, 'name': name, 'flags': flags}

    print(f"\nWarning: Metabolite pattern {species} is not recognized. Will treat as null species.")
    raise SyntaxError(f"Metabolite pattern not recognized! >>> {species}")



