# Microbial Metabolic Network Simulator

## Introduction:

This code was written to simulate the long-term biochemical state of a set of microbially catalyzed reactions.

Given a system of `M` metabolites and `N` reactions, the net rates of metabolite consumption and production may be expressed as the following ordinary differential equation.<sup>1, 2</sup>

```dCdt = S * H```

Where `C` is the M-dimensional vector of metabolite concentrations (in μM),  `S` is the `M` by `N` matrix of stoichiometric coefficients, and `H` is the `N`-dimensional vector of reaction rates (in μM / day). In this equation, `H` is modeled by a vector of stochastic Ornstein-Uhlenbeck processes.

Given a vector of initial metabolite concentrations, the simulation will numerically integrate the above equation until the system arrives at a dead-end state at which the energetic yield of all N reactions is endergonic (∆G<sub>n</sub> ≥ 0, 1 ≤ n ≤ `N`).

## How to Use:
### 1. Executing from the command-line:
`main.py` is the command-line interface to the entire simulation. This is where the user is able to specify what kind of simulation they would like to run, as well as how many runs of that simulation and where to output the data.

In order to execute this simulation, the user would type the following:

```
$ python3 main.py [flags] < input.txt
```

Where `input.txt` is a text file that describes the chemical reactions that will be simulated. This input file is essential and the program will not run if it is not provided. An example of an input file provided with the code is `cariaco.txt`:

```
asos oxH2SrNO3 oxH2SrNO2 amoA nxr anammox
yes
O2
(5, 200)

# input file for simulating Cariaco
# use: python3 main.py < cariaco.txt
```

The first line of `cariaco.txt` is a list of chemical reactions separated by spaces. The next line indicates whether or not we'd like to vary the concentration of a particular metabolite (`yes` or `no`). Assuming the second line is `yes`, the third line `O2` is the name of the metabolite that we would like to vary. Finally, the fourth line is a tuple showing the concentration range for `O2`. 

Any input file written for `main.py` is expected to follow the format of `cariaco.txt`. If the user does not want to vary the concentration of any metabolites, they should put `no` on the second line.

Along with the input file, there are several program flags:
- `--runs` (`-r`) : specifies how many times to execute the simulation in serial. For parallel execution, a user could run `main.py` multiple times via an external shell script (`script.sh`). 
- `--out` (`-o`) : specifies the output directory to store all simulation data.
- `--timestep` (`-t`) : supplies `simulation.py` with the timestep for the integration algorithm (`sim/simulation.py`). The timestep must be in units of days.
- `--debug` (`-d`) : this is a boolean flag. If `--debug=True`, then `main.py` will output additional simulation data along with the dead-end states. See the third section of "How to Use" for more details.
- `--seed` (`-s`) : supplies a random seed for the Ornstein-Uhlenbeck processes. This is important for reproducing previous simulation results.

All of these arguments have default settings if nothing is passed to them:

- `--runs` is `50` by default.
- `--out` is `data` by default.
- `--timestep` is `5.0` (days) by default.
- `--debug` is `False` by default.
- `--seed` is `None` by default.

### 2. Modifying the configuration files:

In the `sim/config` directory there are 4 files the user can modify:

#### `general.py`

This file describes general simulation parameters such as `RUNTIME` and `TEMPERATURE`. These parameters are put in `general.py` because they are not metabolite- or reaction-specific.

#### `metabolites.py` 

This file includes all information about metabolites, including concentrations and free energies of formation values.


#### `reactions.py` 
Includes all information about reactions, including typical rates and other kinetic parameters. All kinetic parameters are used in the generation of stochastic Ornstein-Uhlenbeck processes in order to model random reaction rates (`sim/ornbeck.py`).

#### `stoichiometry.txt`

Describes how each reaction in `reactions.py` alters the concentrations of the metabolites in `metabolites.py`. Since this is a `.txt` file instead of a `.py` file, it needs to be parsed via regular expressions in order to be used in the simulation. This is done using the parsing script in the `sim/parser` directory.

It is very important to keep in mind that any changes made to these 4 files will change the behavior of *any* simulation being run, regardless of what input file is being fed to `main.py`.

### 3. Parsing and visualizing output:

After the simulation has finished executing, there will be an output directory containing all dead-end states in a file called `dead_ends.tsv`. There are also other files that are useful:

- `network_desc.txt` provides a summary of all chemical reactions and metabolites in the simulation.
- `stoich_mat_full.txt`, `stoich_mat_lim.txt`, `stoich_mat_nconst.txt` are the three different stoichiometric matrices used for model calculations (full, limiting, and non-constant respectively). Full stoichiometry is used for calculating energetic values ($∆G$), limiting stoichiometry is used for calculating Ornstein-Uhlenbeck processes, while non-constant stoichiometry is used for numerical integration. Feel free to view the code for more details.

If the simulation is run with the `--debug` flag, an additional directory for each simulation run will be made. In each directory, the following files will be created:

- `composition.tsv` records the metabolite concentrations at each timestep in micromolar units.
- `deltaG.tsv` records $∆G$ values at each timestep in units of kJ per mole.
- `flux.tsv` records net change in metabolite concentrations at each timestep in micromolar units.
- `ornbeck.tsv` records the value of the vector of Ornstein-Uhlenbeck processes for all reactions.
- `initial_condition.txt` records the initial concentration of the metabolite we are varying (if we specify to vary a metabolite in the input file).
- `messages.txt` - records timesteps for any notable issues during simulation execution, such as negative concentrations, changes in step size, and early termination time.
- `report.txt` - records whether or not the simulation has actually reached a dead-end state.

The data in `dead_ends.tsv` can be used to create a series of bifurcation plots for visualizing the distribution of end states for our differential equation model. To see how this can be produced, look at the code in `visualization.ipynb.`

## Author:
- Nathan Malamud, undergraduate student at the University of Oregon

## Supervisor:
- Dr. Stilianos Louca, professor at the Institute of Ecology and Evolution

## References:

1. Louca, Scranton, Taylor., Astor., Crowe, & Doebeli (2019). Circumventing kinetics in biogeochemical modeling. PNAS 116: 11329-11338 
1. Louca, & Doebeli (2016). Reaction-centric modeling of microbial ecosystems. Ecological Modelling  335: 74-86
1. Press, Teukolsky, Vetterling &, Flannery (2007). Integration of Ordinary Differential Equations: Runge-Kutta Method. Numerical Recipes: The Art of Scientific Computing (pp. 907-910). Cambridge University Press
