# ODE model of microbial metabolism.
#
# Source: Project prompt from Dr. Stilianos Louca.
# Author: Nathan Malamud

import numpy as np

def ode_model(t, C, S_mats, G, X):
    """
    Non-autonomous ODE: returns dCdt.
    """
    S_full, S_lim, S_nconst = S_mats
    M, N = S_full.shape

    H = np.zeros(N)

    for n in range(N):
        if G[n] >= 0:
            continue

        H[n] = X[n]

        for m in range(M):
            if S_lim[m, n] < 0:
                H[n] *= C[m]
        
    return S_nconst @ H
