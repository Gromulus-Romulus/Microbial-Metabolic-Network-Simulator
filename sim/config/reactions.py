import numpy as np
import numpy.random as np_random

reactions = np.array([
    ## Reed et al. (2013)
    'cox',      # aerobic respiration via cytochrome-c oxidase
    'narG',     # nitrate reduction
    'nirK',     # dissimilatory nitrtie reduction to ammonium
    'nrf',      # sulfate reduction
    'dsr',      # sulfide oxidation coupled to nitrate reduction
    'nap',      # sulfide oxidation coupled to nitrate reduction
    'sox',      # aerobic sulfide oxidation w/ bicarbonate, water, and CO2
    'asos',     # aerobic sulfide oxidation
    'amoA',     # aerobic ammonium oxidation
    'hzo',      # anaerobic ammonium oxidation
    'anammox',  # anaerobic ammonium oxidation using nitrite
    'nxr',      # aerobic nitrtie oxidation

    ## Louca et al. (2019)
    'oxH2SrNO3',    # oxidation of sulfide to sulfate w/ reduction of nitrate to nitrite
    'oxH2SrNO2',    # oxidation of sulfide to sulfate w/ reduction of nitrite to N2

    ## Lu et al. (2020)
    'B1',   # hydrogen oxidation
    'I5',   # hydrogen oxidation
    'B3',   # methane oxidation 
    'I14',  # methane oxidation
    'M6',   # methane oxidation
    'B9',   # ammonium oxidation
    'B14',  # sulfide oxidation
    'B15',  # sulfide oxidation
    'B16',  # sulfur oxidation
    'I29',  # sulfide oxidation
    'D3',   # methanogenesis
    'M2',   # sulfate reduction

    # Louca & Doebeli. (2017)
    'C',    # glucose fermentation
    'D',    # glucose fermentation
    'K',    # methanogenesis
    'L'     # methanogenesis
])

# - - - - - - - - - - - - - - - - - - - - - - - //
# PARAMETERS FOR ORNSTEIN-UHLENBECK PROCESSES
# - - - - - - - - - - - - - - - - - - - - - - - //

# units - uM / day
# format - (min, max) # reaction name
typical_rates = np.array([
    (0.1, 50), # cox
    (0.1, 50), # narG
    (0.1, 50), # nirK
    (0.1, 50), # nrf
    (0.1, 50), # dsr
    (0.1, 50), # nap
    (0.1, 50), # sox
    (0.1, 50), # asos
    (0.1, 50), # amoA
    (0.1, 50), # hzo
    (0.1, 50), # anammox
    (0.1, 50), # nxr
    (0.1, 50), # oxH2SrNO4
    (0.1, 50), # oxH2SrNO2
    (0.1, 50), # B1
    (0.1, 50), # I5
    (0.1, 50), # B3
    (0.1, 50), # I14
    (0.1, 50), # M6
    (0.1, 50), # B9
    (0.1, 50), # B14
    (0.1, 50), # B15
    (0.1, 50), # B16
    (0.1, 50), # I29
    (0.1, 50), # D3
    (0.1, 50), # M2
    (0.1, 50), # C
    (0.1, 50), # D
    (0.1, 50), # K
    (0.1, 50)  # L
])

# units - uM / day
# format - (min, max) # reaction name
typical_decay = np.array([
    (0.1, 20), # cox
    (0.1, 20), # narG
    (0.1, 20), # nirK
    (0.1, 20), # nrf
    (0.1, 20), # dsr
    (0.1, 20), # nap
    (0.1, 20), # sox
    (0.1, 20), # asos
    (0.1, 20), # amoA
    (0.1, 20), # hzo
    (0.1, 20), # anammox
    (0.1, 20), # nxr
    (0.1, 20), # oxH2SrNO4
    (0.1, 20), # oxH2SrNO2
    (0.1, 20), # B1
    (0.1, 20), # I5
    (0.1, 20), # B3
    (0.1, 20), # I14
    (0.1, 20), # M6
    (0.1, 20), # B9
    (0.1, 20), # B14
    (0.1, 20), # B15
    (0.1, 20), # B16
    (0.1, 20), # I29
    (0.1, 20), # D3
    (0.1, 20), # M2
    (0.1, 20), # C
    (0.1, 20), # D
    (0.1, 20), # K
    (0.1, 20)  # L
])

# units - % value (number between 0 and 1)
# format - (min, max) # reaction name
typical_std = np.array([
    (0.1, 1.0), # cox
    (0.1, 1.0), # narG
    (0.1, 1.0), # nirK
    (0.1, 1.0), # nrf
    (0.1, 1.0), # dsr
    (0.1, 1.0), # nap
    (0.1, 1.0), # sox
    (0.1, 1.0), # asos
    (0.1, 1.0), # amoA
    (0.1, 1.0), # hzo
    (0.1, 1.0), # anammox
    (0.1, 1.0), # nxr
    (0.1, 1.0), # oxH2SrNO4
    (0.1, 1.0), # oxH2SrNO2
    (0.1, 1.0), # B1
    (0.1, 1.0), # I5
    (0.1, 1.0), # B3
    (0.1, 1.0), # I14
    (0.1, 1.0), # M6
    (0.1, 1.0), # B9
    (0.1, 1.0), # B14
    (0.1, 1.0), # B15
    (0.1, 1.0), # B16
    (0.1, 1.0), # I29
    (0.1, 1.0), # D3
    (0.1, 1.0), # M2
    (0.1, 1.0), # C
    (0.1, 1.0), # D
    (0.1, 1.0), # K
    (0.1, 1.0)  # L
])
