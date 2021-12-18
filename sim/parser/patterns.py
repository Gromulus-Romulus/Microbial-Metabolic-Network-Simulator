# Regex patterns for file parsing.
# Author: Nathan Malamud

import re

CommentPat = re.compile(r"""
              \s*

              # note that empty lines are also
              # recognized as comments.

              (
                  [\#;].*
              )?
              \s*
              """, re.VERBOSE | re.MULTILINE)

ReactionPat = re.compile(r"""
              \s*
              (?P<name> [\S(?!:)]+)    # Reaction Name
              \s* : \s*

              (?P<reactants>
                  ([\S(?!,)]+ \s*, \s*)* ([\S(?!,)]+)
              )

              \s* -> \s*

             (?P<products>
                 ([\S(?!,)]+ \s*, \s*)* ([\S(?!,)]+)
             ) \s*

             # optional comment
             (
                 [\#;].*
             )?
             \s*
             """, re.VERBOSE | re.MULTILINE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  SpeciesPat is used exclusively in parse_species.   #
# - - - - - - - - - - - - - - - - - - - - - - - - - - #

SpeciesPat = re.compile(r"""
             (\[ (?P<coeff>
                 [-]?[0-9]+(\.[0-9]+)?
                 ([eE][-]?[0-9]+)?
             ) \])?

             (?P<name> [A-Za-z0-9-\{\}\-\+]+)

             (\[(?P<flags> \s* ([A-Za-z0-9-\{\}]+ \s* [\|] )*
               \s* [A-Za-z0-9-\{\}]+ \s*)\])?
             """, re.VERBOSE | re.MULTILINE)

