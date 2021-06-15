# file: test_napmo_system.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from napmo.system.input_parser import InputParser
from napmo.system.napmo_system import NAPMO
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
    system.solve()
    e = system._energy
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
    system.solve()
    e = system._energy

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
    system.solve()
    e = system._energy

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

code {
ref_energy = -0.66217971
diff_energy = ref_energy - Energy
}
""")

    data = InputParser(file)
    system = NAPMO(data)
    pe = -0.66217971
    system.solve()
    e = system._energy
    system.exec_code()

    np.testing.assert_allclose(pe, e, rtol=10e-6)


test_napmo_single()
# test_napmo_multi()
