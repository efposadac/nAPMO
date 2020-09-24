# file: interspecies_correlation.py
# nAPMO package
# Copyright (c) 2020, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
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


def epc17_2(rho, other_rho):
    """
    Calculates the epc17-2 interspecies functional by Hammes-Schiffer
    see http://dx.doi.org/10.1021/acs.jpclett.7b01442
    """
    a = 2.35
    b = 2.4
    c = 6.6

    denominator = a - b * np.sqrt(rho * other_rho) + (c * rho * other_rho)
    c_energy = (-rho * other_rho) / denominator
    c_potential = ((b * np.sqrt(rho) * other_rho**1.5 ) - (2.0 * a * other_rho)) / denominator**2 / 2.0
    c_other_potential = ((b * np.sqrt(other_rho) * rho**1.5 ) - (2.0 * a * rho)) / denominator**2 / 2.0

    return c_energy, c_potential, c_other_potential
