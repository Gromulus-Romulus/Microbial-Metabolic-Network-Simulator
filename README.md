# Microbial Metabolic Network Simulator

## Introduction:

This code was written to simulate the long-term biochemical state of a set of microbially catalyzed reactions.

Given a system of `M` metabolites and `N` reactions, the net rates of metabolite consumption and production may be expressed as the following ordinary differential equation.<sup>1, 2</sup>

```dCdt = S * H```

Where `C` is the M-dimensional vector of metabolite concentrations (in μM),  `S` is the `M` by `N` matrix of stoichiometric coefficients, and `H` is the `N`-dimensional vector of reaction rates (in μM / day).

Given a vector of initial metabolite concentrations, the numerical algorithm described in the following pages will integrate Eq. 1 until the system arrives at a dead-end state at which the energetic yield of all N reactions is endergonic (∆G<sub>n</sub> ≥ 0, 1 ≤ n ≤ `N`).

### 1. Executing from the command-line:
`main.py` is the command-line interface to the entire simulation. This is where the user is able to specify what kind of simulation they would like to run, as well as how many runs of that simulation to run and where to output the data.

In order to execute this simulation, the user would type...

```
$ python3 main.py [flag1=val1, flag2=val2, ...]
```

...followed by one or more of the following flags:

- `--mode` (`-m`) : specifies whether to execute a simulation with a specified set of reactions (`--mode=fixed`) or to execute a simulation with a randomized set of reactions (`--mode=random`). There are only two options for `--mode`.
- `--runs` (`-r`) : specifies how many times to execute the simulation in serial. For parallel execution, a user could run `main.py` multiple times via an external shell script (e.g., `script.sh`). 
- `--out` (`-o`) : specifies output directory. Will create a new folder to write all dead-end states to.
- `--timestep` (`-t`) : supplies `simulation.py` with the timestep for the ODE integration algorithm. The timestep must be in units of days.
- `--debug` (`-d`) : this is a boolean flag. If `--debug=True`, then `main.py` will output additional simulation data along with the dead-end states.
- `--seed` (`-s`) : supplies a random seed for the Ornstein-Uhlenbeck processes. This is important for reproducing previous simulation results.

All of these arguments have default settings if nothing is passed to them:

- `--mode` is `random` by default.
- `--runs` is `50` by default.
- `--out` is `data` by default.
- `--timestep` is `5.0` (days) by default.
- `--debug` is `False` by default.
- `--seed` is `None` by default.

**Important things to Note**:

- Regardless of what `mode` the simulation is executed in, there will be a series of input prompts for the user to fill in. For `--mode=random`, the user will need to specify the number of reactions needed in the simulation. For `--mode=fixed`, the user will need to specify whether or not they would like to vary the concentration of a certain metabolite. This may affect how the simulation is being used in shell scripts.

- Whatever output directory specified to `--out` will be overwritten if it already exists.

- Running `main.py` with `--debug` will use up significantly more storage space. A comprehensive test suite has not been done to determine exactly how much more memory would be occupied.

### 2. Modifying the configuration files:

TODO

## Examples of Use:

To further illustrate how this code is meant to be used, here are some concrete examples.

## 1. 

TODO

## 2. 

TODO

## References:

1. Louca, Scranton, Taylor., Astor., Crowe, & Doebeli (2019). Circumventing kinetics in biogeochemical modeling. PNAS 116: 11329-11338 
1. Louca, & Doebeli (2016). Reaction-centric modeling of microbial ecosystems. Ecological Modelling  335: 74-86
1. 
1.
1.
1. Press, Teukolsky, Vetterling &, Flannery (2007). Integration of Ordinary Differential Equations: Runge-Kutta Method. Numerical Recipes: The Art of Scientific Computing (pp. 907-910). Cambridge University Press
1. 
