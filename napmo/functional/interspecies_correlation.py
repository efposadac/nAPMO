# file: interspecies_correlation.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np


def isc_functional_selector(name):
    """
    Returns the function for a given functional name

    Args:
        name (str): Functional name
    """
    if name is not None:
        assert isinstance(name, str)

        database = {
            'epc17-2': napmo.epc17_2
        }

        functional = database.get(name)

        if functional is None:
            raise NotImplementedError(name + " Functional NOT available!")
        else:
            return functional
    else:
        return None


def epc17_2(rhoE, rhoN):
    """
    Calculates the epc17-2 interspecies functional by Hammes-Schiffer

    references:
        http://dx.doi.org/10.1021/acs.jpclett.7b01442

    Args:
        rhoE (ndarray): Density on the grid for electrons
        rhoN (ndarray): Density on the grid for nuclei

    Notes:
        the dimesions of ``rhoE`` and ``rhoN`` must be the same, you have to choose the common points
        between the grids or calculate the extrapolation between them.
    """
    a = 2.35
    b = 2.4
    c = 6.6

    denominator = a - b * np.sqrt(rhoE * rhoN) + c * rhoE * rhoN

    # c_ene = -rhoE * rhoN / denominator
    c_ene = -rhoN / denominator

    if np.any(np.abs(rhoE) < 1.0e-200) or np.any(np.abs(rhoN) < 1.0e-200):
        print("USING OLD FUNCTIONAL")
        c_pot_A = (b * np.sqrt(rhoE) * rhoN**(1.5) - 2.0 * a * rhoN) / denominator**2 / 2.0
        c_pot_B = (b * np.sqrt(rhoN) * rhoE**(1.5) - 2.0 * a * rhoE) / denominator**2 / 2.0
    else:
        c_pot_A = (rhoE * rhoN * (c * rhoN - b * rhoN / (2.0 * np.sqrt(rhoE * rhoN))) - rhoN * denominator) / denominator**2
        c_pot_B = (rhoN * rhoE * (c * rhoE - b * rhoE / (2.0 * np.sqrt(rhoE * rhoN))) - rhoE * denominator) / denominator**2

    return c_ene, c_pot_A, c_pot_B
