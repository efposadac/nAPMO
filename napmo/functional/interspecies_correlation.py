# file: interspecies_correlation.py
# nAPMO package
# Copyright (c) 2020, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

import napmo
import numpy as np


def isc_functional_selector(name):
    """
    Returns the function for a given functional name
    """
    assert isinstance(name, str)

    database = {
        'epc17-2': napmo.epc17_2
    }

    functional = database.get(name)

    if functional is None:
        raise NotImplementedError(name+" Functional NOT available!")
    else:
        return functional


def epc17_2(rhoE, rhoN):
    """
    Calculates the epc17-2 interspecies functional by Hammes-Schiffer
    see http://dx.doi.org/10.1021/acs.jpclett.7b01442
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
