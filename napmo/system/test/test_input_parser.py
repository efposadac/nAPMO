# file: test_input_parser.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from napmo.system.input_parser import *


def test_input_parser_interface():
    data = ("""# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
    H_1   0.000000     -0.754174      0.528382
}

# basis-set specification
basis {
    e- 6-311G
   H_1 DZSNB
}

# scf options
scf {
    direct
    hf
    analytic
}

# energy of reference
code {
ref_energy = -76.0001
diff = ref_energy - Energy
}

# basis STO-3G
""")

    test = InputParser(data)
    print(test.data)
    print(test.scf)
    print(test.code)

    assert fix_casting('5') == 5
    assert fix_casting('5.0') == 5.0
    assert fix_casting('hello') == 'hello'

    try:
        raise_exception(ValueError, 'test', 'test')
        assert False, 'Expecting Failure!'
    except ValueError:
        assert True

    data = ("""# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
}

# basis-set specification
basis = 6-311G

# scf options
scf {
    hf
    analytic = False
}

# energy of reference
code {
ref_energy = -76.0001
diff = ref_energy - Energy
}

""")

    test = InputParser(data)
    print(test.data)
    print(test.scf)
    print(test.code)

    data = ("""# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
}

# basis-set specification
basis {
    e- 6-311G /usr/local/lib
}

# scf options
scf {
    hf
    analytic = False
}

# energy of reference
code {
ref_energy = -76.0001
diff = ref_energy - Energy
}
""")

    test = InputParser(data)
    print(test.data)
    print(test.scf)
    print(test.code)

    data = ("""# Molecule definition
molecule {
    O   0.000000      0.000000     -0.066575
    H_1   0.000000      0.754175      0.528381
}

# basis-set specification
basis = 6-311G

# scf options
scf {
    hf
    analytic
    direct = False
}

# energy of reference
ref_energy = -76.0001
diff = ref_energy - Energy

""")
    try:
        test = InputParser(data)
        assert False, "Expecting Failure"
    except ValueError:
        assert True

    data = ("""# Molecule definition
# basis-set specification
basis = 6-311G

# scf options
scf {
    hf
    analytic
    direct = False
}

# energy of reference
ref_energy = -76.0001
diff = ref_energy - Energy

""")
    try:
        test = InputParser(data)
        assert False, "Expecting Failure"
    except ValueError:
        assert True

    data = ("""# Molecule

molecule [e-:0] {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
}

# basis-set specification
basis {
    e- 6-311G
    H_1 Nakai-5-sp
}

# scf options
scf {
    hf
    analytic
    direct = False
}

# energy of reference
ref_energy = -76.0001
diff = ref_energy - Energy

""")
    try:
        test = InputParser(data)
        assert False, "Expecting Failure"
    except ValueError:
        assert True

    data = ("""# Molecule definition
molecule [e-:0] {
    O   0.000000      0.000000     -0.066575
    H   0.000000      0.754175      0.528381
}

# basis-set specification
basis = 6-311G

# scf options
scf {
    hf
    analytic
    debug = True
}

# energy of reference
code {
ref_energy = -76.0001
diff = ref_energy - Energy
}
""")

    test = InputParser(data)
    assert test.scf.get('debug') is True

# test_input_parser_interface()
