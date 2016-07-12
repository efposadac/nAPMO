import os

from napmo.system.timer import timeblock
from napmo.system.atomic_element import AtomicElement
from napmo.system.elementary_particle import ElementaryParticle
from napmo.system.basis_set import BasisSet, BasisSet_C
from napmo.system.contracted_gaussian import ContractedGaussian
from napmo.system.primitive_gaussian import PrimitiveGaussian
from napmo.system.input_parser import InputParser
from napmo.system.napmo import NAPMO
from napmo.system.molecular_system import MolecularSystem
from napmo.system.cext import napmo_library as cext

from napmo.data.databases import AtomicElementsDatabase
from napmo.data.constants import ANGSTROM_TO_BOHR

from napmo.hf.wavefunction import WaveFunction
from napmo.hf.scf import SCF


# Basis-set path
basis_dir = os.path.join(os.path.dirname(__file__), 'data/basis')
