# file: test_napmo.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from napmo.system.input_parser import InputParser
from napmo.system.napmo import NAPMO
import numpy as np


def test_napmo_single():
    file = ("""
# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
    H   0.000000     -0.754174      0.528382
}

# basis-set specification
basis {
    e- 6-311G
}

# scf options
scf {
    direct
    hf
    analytic
}
""")

    data = InputParser(file)
    system = NAPMO(data)
    pe = -76.008843007747
    e = system.solve()
    print(e)
    np.testing.assert_allclose(pe, e)

    file = ("""
# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
    H   0.000000     -0.754174      0.528382
}

# basis-set specification
basis {
    e- 6-311G
}

# scf options
scf {
    hf
    analytic
}
""")

    data = InputParser(file)
    system = NAPMO(data)
    pe = -76.008843007747
    e = system.solve()

    np.testing.assert_allclose(pe, e)


def test_napmo_multi():

    file = ("""
# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
    H_1   0.000000     -0.754174      0.528382
}

# basis-set specification
basis {
    e- 6-311G
   H_1 NAKAI-5-SP
}

# scf options
scf {
    maxiter = 1000
    # direct
    hf
    analytic
}
""")

    data = InputParser(file)
    system = NAPMO(data)
    pe = -75.969790197894
    e = system.solve()

    np.testing.assert_allclose(pe, e, rtol=10e-6)

    file = ("""
# Molecule definition
molecule [e-:1] {
    H  0.00 0.00 0.00
    E+ 0.00 0.00 0.00
}

# basis-set specification
basis {
   e- CC-PVTZ
   E+ E+-H-7SP-AUG-CC-PVTZ
}

# scf options
scf {
    maxiter = 1000
    # direct
    uhf
    analytic
}
""")

    data = InputParser(file)
    system = NAPMO(data)
    pe = -0.66217971
    e = system.solve()

    np.testing.assert_allclose(pe, e, rtol=10e-6)

# test_napmo_single()
# test_napmo_multi()
