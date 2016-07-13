import os
import numpy as np
import numpy.ctypeslib as npct

from ctypes import *

from napmo.system.cext import napmo_library as cext
from napmo.system.timer import timeblock
from napmo.system.input_parser import InputParser
from napmo.system.atomic_element import AtomicElement
from napmo.system.elementary_particle import ElementaryParticle
from napmo.system.primitive_gaussian import PrimitiveGaussian
from napmo.system.contracted_gaussian import ContractedGaussian
from napmo.system.basis_set import BasisSet, BasisSet_C
from napmo.system.molecular_system import MolecularSystem
from napmo.system.napmo import NAPMO

from napmo.grids.radial import RadialGrid
from napmo.grids.radial_cheb import RadialGridCheb
from napmo.grids.radial_transform import RadialTransform
from napmo.grids.radial_transform import IndentityRadialTransform
from napmo.grids.radial_transform import PowerRadialTransform
from napmo.grids.radial_transform import ChebyshevRadialTransform
from napmo.grids.angular import AngularGrid
from napmo.grids.atomic import AtomicGrid
from napmo.grids.becke import BeckeGrid
from napmo.grids.cubic_spline import CubicSpline
from napmo.grids.extrapolation import Extrapolation
from napmo.grids.extrapolation import CuspExtrapolation
from napmo.grids.extrapolation import PowerExtrapolation

from napmo.utilities.cell import Cell
from napmo.utilities.ode2 import solve_ode2
from napmo.utilities.int1d import Integrator1D
from napmo.utilities.int1d import StubIntegrator1D
from napmo.utilities.int1d import TrapezoidIntegrator1D
from napmo.utilities.int1d import CubicIntegrator1D
from napmo.utilities.int1d import SimpsonIntegrator1D

from napmo.data.databases import AtomicElementsDatabase
from napmo.data.databases import ElementaryParticlesDatabase
from napmo.data.databases import CouplingConstantsDatabase
from napmo.data.constants import ANGSTROM_TO_BOHR
from napmo.data.constants import PROTON_MASS
from napmo.data.constants import NEUTRON_MASS
from napmo.data.constants import SPIN_ELECTRON

from napmo.hf.wavefunction import WaveFunction
from napmo.hf.scf import SCF


# Basis-set path
basis_dir = os.path.join(os.path.dirname(__file__), 'data/basis')

# Types for ctypes
a1df = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
a2df = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

# C functions

# Libint

cext.LibintInterface_new.restype = c_void_p
cext.LibintInterface_new.argtypes = [c_int]

cext.LibintInterface_del.restype = None

cext.LibintInterface_add_basis.restype = None
cext.LibintInterface_add_basis.argtypes = [
    c_void_p, POINTER(BasisSet_C)]

cext.LibintInterface_add_pointcharges.restype = None
cext.LibintInterface_add_pointcharges.argtypes = [
    c_void_p, c_int, a1df]

cext.LibintInterface_get_nbasis.restype = c_int
cext.LibintInterface_get_nbasis.argtypes = [c_void_p]

cext.LibintInterface_compute_1body_ints.restype = None
cext.LibintInterface_compute_1body_ints.argtypes = [
    c_void_p, c_int, a2df]

cext.LibintInterface_init_2body_ints.restype = None
cext.LibintInterface_init_2body_ints.argtypes = [c_void_p]

cext.LibintInterface_compute_2body_ints.restype = c_void_p
cext.LibintInterface_compute_2body_ints.argtypes = [
    c_void_p, a2df]

cext.LibintInterface_compute_2body_direct.restype = None
cext.LibintInterface_compute_2body_direct.argtypes = [
    c_void_p, a2df, a2df]

cext.LibintInterface_init_2body_ints.restype = None
cext.LibintInterface_init_2body_ints.argtypes = [c_void_p]

cext.LibintInterface_compute_coupling_direct.restype = None
cext.LibintInterface_compute_coupling_direct.argtypes = [
    c_void_p, c_void_p, a2df, a2df]

# Wavefunction

cext.wavefunction_guess_hcore.restype = None
cext.wavefunction_guess_hcore.argtypes = [POINTER(WaveFunction)]

cext.wavefunction_iterate.restype = None
cext.wavefunction_iterate.argtypes = [
    POINTER(WaveFunction)]

cext.wavefunction_compute_energy.restype = None
cext.wavefunction_compute_energy.argtypes = [
    POINTER(WaveFunction)]

cext.wavefunction_compute_density.restype = None
cext.wavefunction_compute_density.argtypes = [
    POINTER(WaveFunction)]

cext.wavefunction_compute_2body_matrix.restype = None
cext.wavefunction_compute_2body_matrix.argtypes = [
    POINTER(WaveFunction), c_void_p]

# Primitive Gaussian
cext.gto_normalize_primitive.restype = c_double
cext.gto_normalize_primitive.argtypes = [
    POINTER(PrimitiveGaussian)]

cext.gto_compute_primitive.restype = None
cext.gto_compute_primitive.argtypes = [
    POINTER(PrimitiveGaussian),
    a2df, a1df, c_int]

cext.gto_overlap_primitive.restype = c_double
cext.gto_overlap_primitive.argtypes = [
    POINTER(PrimitiveGaussian),
    POINTER(PrimitiveGaussian)]
