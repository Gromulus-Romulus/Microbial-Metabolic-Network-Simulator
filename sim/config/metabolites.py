import numpy as np

metabolites = np.array([
    'C6H12O6',  # glucose
    'CH3COOH',  # acetic acid
    'CH4',      # methane
    'CO',       # carbon monoxide
    'CO2',      # carbon dioxide
    'HCO3-',    # bicarbonate
    'O2',       # oxygen
    'H2O',      # water
    'H+',       # hydrogen ions
    'H2',       # hydrogen gas
    'NO2-',     # nitrite
    'NO3-',     # nitrate
    'N2',       # nitrogen gas
    'NH4+',     # ammonium
    'S',        # sulfur
    'SO4-2',    # sulfate
    'H2S'       # sulfide
])

# units - kJ / mol
# source: CHNOSZ package from CRAN library (R programming lang)
# ask Dr. Louca for the source code for finding these values.

deltaGf0 = np.array([
    -915.288,   # C6H12O6
    -534.715,   # CH3COOH
    -34.0578,   # CH4
    -120.005,   # CO
    -385.974,   # CO2
    -586.94,    # HCO3-
    16.5435,    # O2
    0.00000,    # H2O
    0.00000,    # H+
    17.7234,    # H2
    -32.2168,   # NO2-
    -110.905,   # NO3-
    18.1878,    # N2
    -79.4542,   # NH4+
    0.00000,    # S
    -744.459,   # SO4-2
    -27.9198    # H2S
])

# units - uM (micro-molar)
# source: email exchange w/ Dr. Louca
# it is expected that these concentrations
# somewhat reflect what is typical
# in an anoxic ocean basin (e.g. Cariaco, Venezuela; Saanich Inlet, Canada).

# NOTE: Do not set these concentrations equal to 0.
# This will cause complications in the calculate_gibbs function
# in simulation.py, since the function takes the natural logarithim
# of initalC. Instead, set them to some value close to 0.

initialC = np.array([
    5.0,      # C6H12O6 
    5.0,      # CH3COOH
    5.0,      # CH4
    5.0,      # CO
    5.0,      # CO2
    5.0,      # HCO3-
    20.0,     # O2
    5.5e7,    # H2O
    2.399e-2, # H+
    5.0,      # H2
    0.100,    # NO2-
    10.00,    # NO3-
    470.0,    # N2
    25.00,    # NH4+
    5.0,      # S
    28000,    # SO4-2
    50.00     # H2S
])
