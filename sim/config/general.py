# Author: Nathan Malamud

# - - - - - - - - - - - - #
#    GENERAL PARAMETERS   #
# - - - - - - - - - - - - #

# Max simulation runtime (in days)
RUNTIME = 100 * 365 
# Simulation temperature (in Kelvin)
TEMPERATURE = 298

# Optional: tell the simulation to stop after a certain number of iterations
# Ad-Hoc way to avoid infinite loops (or really long simulations).
MAX_ITERATIONS = 1_000_000

# Error bounds for Fehlberg integration scheme (units in uM)
E_MIN = 0.1 
E_MAX = 1 

# Boundary condition - simulation will terminate early
# if ∆G for all reactions is above this value
DELTA_G_BOUND = -1
