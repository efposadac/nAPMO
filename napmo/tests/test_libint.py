from __future__ import division
from __future__ import print_function
from napmo.system.molecular_system import MolecularSystem
from napmo.utilities.libint_iface import Libint


def printMatrix(a, name):
    print(name, " [" + ("%d" % a.shape[0]) + "][" + ("%d" % a.shape[1]) + "]")
    rows = a.shape[0]
    cols = a.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            print("%6.3f" % a[i, j], end=" ")
        print()
    print()

basis_file = "results/Coulomb/TEST.json"
molecule = MolecularSystem()
molecule.add_atom("H", [0.000000, 0.754175, 0.528381],
                  basis_kind="GTO", basis_file=basis_file)
molecule.add_atom("H", [0.000000, -0.754175, 0.528381],
                  basis_kind="GTO", basis_file=basis_file)
molecule.show()

integrals = Libint(molecule)
printMatrix(integrals.get_overlap_matrix(), 'Overlap')
printMatrix(integrals.get_kinetic_matrix(), 'Kinetic')
printMatrix(integrals.get_nuclear_matrix(), 'Nuclear')
