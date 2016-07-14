import os
import numpy as np
import numpy.ctypeslib as npct

from ctypes import *

from napmo.system.cext import napmo_library as cext
from napmo.system.timer import Timer
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

# OMP Threads
threads = int(os.environ.get('OMP_NUM_THREADS', 1))

# Timer singleton
runtime = Timer()

# Basis-set path
basis_dir = os.path.join(os.path.dirname(__file__), 'data/basis')

# Types for ctypes
a1df = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
a2df = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
a1di = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
aptr = npct.ndpointer(c_void_p, flags="C_CONTIGUOUS")

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

# Angular grid
cext.angular_cartesian.restype = None
cext.angular_cartesian.argtypes = [POINTER(AngularGrid)]

cext.angular_to_spherical.restype = None
cext.angular_to_spherical.argtypes = [POINTER(AngularGrid)]

cext.angular_integrate.restype = c_double
cext.angular_integrate.argtypes = [
    POINTER(AngularGrid), c_int, a1df]

# Atomic grid
cext.atomic_grid_init.restype = None
cext.atomic_grid_init.argtypes = [
    POINTER(AtomicGrid), POINTER(AngularGrid),
    POINTER(RadialGrid)]

cext.angular_spherical_expansion.restype = None
cext.angular_spherical_expansion.argtypes = [
    POINTER(AngularGrid), c_int, c_int, a1df, a2df]

cext.angular_eval_expansion.restype = None
cext.angular_eval_expansion.argtypes = [
    POINTER(AngularGrid), c_int, c_int, a2df, a1df]

cext.atomic_grid_integrate.restype = None
cext.atomic_grid_integrate.argtypes = [
    POINTER(AtomicGrid), c_int, c_int, c_int, a1df, a1df]

# Becke grid
cext.becke_weights.restype = None
cext.becke_weights.argtypes = [POINTER(BeckeGrid), a1df]

cext.eval_decomposition_grid.restype = None
cext.eval_decomposition_grid.argtypes = [
    aptr, a1df, a1df, a2df, c_void_p, c_long, c_long]

# CubicSpline
cext.CubicSpline_new.restype = c_void_p
cext.CubicSpline_new.argtypes = [
    a1df, a1df, c_void_p, c_void_p, c_int]

cext.CubicSpline_del.restype = None
cext.CubicSpline_del.argtypes = [c_void_p]

cext.CubicSpline_eval.restype = None
cext.CubicSpline_eval.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.CubicSpline_eval_deriv.restype = None
cext.CubicSpline_eval_deriv.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.CubicSpline_get_first_x.restype = c_double
cext.CubicSpline_get_first_x.argtypes = [c_void_p]

cext.CubicSpline_get_last_x.restype = c_double
cext.CubicSpline_get_last_x.argtypes = [c_void_p]

cext.CubicSpline_get_rtransform.restype = c_void_p
cext.CubicSpline_get_rtransform.argtypes = [c_void_p]

cext.CubicSpline_get_extrapolation.restype = c_void_p
cext.CubicSpline_get_extrapolation.argtypes = [c_void_p]

cext.solve_cubic_spline_system.restype = None
cext.solve_cubic_spline_system.argtypes = [
    a1df, a1df, c_int]


# Extrapolation
cext.Extrapolation_prepare.restype = None
cext.Extrapolation_prepare.argtypes = [c_void_p, c_void_p]

cext.Extrapolation_del.restype = None
cext.Extrapolation_del.argtypes = [c_void_p]

cext.Extrapolation_eval_left.restype = c_double
cext.Extrapolation_eval_left.argtypes = [c_void_p, c_double]

cext.Extrapolation_eval_right.restype = c_double
cext.Extrapolation_eval_right.argtypes = [c_void_p, c_double]

cext.Extrapolation_deriv_left.restype = c_double
cext.Extrapolation_deriv_left.argtypes = [c_void_p, c_double]

cext.Extrapolation_deriv_right.restype = c_double
cext.Extrapolation_deriv_right.argtypes = [c_void_p, c_double]

cext.Extrapolation_has_tail.restype = c_bool
cext.Extrapolation_has_tail.argtypes = [c_void_p]

cext.ZeroExtrapolation_new.restype = c_void_p

cext.CuspExtrapolation_new.restype = c_void_p

cext.PowerExtrapolation_new.restype = c_void_p
cext.PowerExtrapolation_new.argtypes = [c_double]

cext.PowerExtrapolation_get_power.restype = c_double
cext.PowerExtrapolation_get_power.argtypes = [c_void_p]

# Radial grid
cext.radial_integrate.restype = c_double
cext.radial_integrate.argtypes = [
    POINTER(RadialGrid), c_int, a1df]

# RadialCheb
cext.radial_init.restype = None
cext.radial_init.argtypes = [POINTER(RadialGridCheb)]

cext.radial_get_z.restype = None
cext.radial_get_z.argtypes = [POINTER(RadialGridCheb)]

cext.radial_deriv_z.restype = None
cext.radial_deriv_z.argtypes = [POINTER(RadialGridCheb)]

cext.radial_deriv2_z.restype = None
cext.radial_deriv2_z.argtypes = [POINTER(RadialGridCheb)]

# Radial Transform
cext.RTransform_get_npoint.restype = c_int
cext.RTransform_get_npoint.argtypes = [c_void_p]

cext.RTransform_del.restype = None
cext.RTransform_del.argtypes = [c_void_p]

cext.RTransform_radius.restype = c_double
cext.RTransform_radius.argtypes = [c_void_p, c_double]

cext.RTransform_deriv.restype = c_double
cext.RTransform_deriv.argtypes = [c_void_p, c_double]

cext.RTransform_deriv2.restype = c_double
cext.RTransform_deriv2.argtypes = [c_void_p, c_double]

cext.RTransform_deriv3.restype = c_double
cext.RTransform_deriv3.argtypes = [c_void_p, c_double]

cext.RTransform_inv.restype = c_double
cext.RTransform_inv.argtypes = [c_void_p, c_double]

cext.RTransform_radius_array.restype = None
cext.RTransform_radius_array.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.RTransform_deriv_array.restype = None
cext.RTransform_deriv_array.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.RTransform_deriv2_array.restype = None
cext.RTransform_deriv2_array.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.RTransform_deriv3_array.restype = None
cext.RTransform_deriv3_array.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.RTransform_inv_array.restype = None
cext.RTransform_inv_array.argtypes = [
    c_void_p, a1df, a1df, c_int]

cext.IdentityRTransform_new.restype = c_void_p
cext.IdentityRTransform_new.argtypes = [c_int]

cext.PowerRTransform_new.restype = c_void_p
cext.PowerRTransform_new.argtypes = [c_double, c_double, c_int]

cext.PowerRTransform_get_rmin.restype = c_double
cext.PowerRTransform_get_rmin.argtypes = [c_void_p]

cext.PowerRTransform_get_rmax.restype = c_double
cext.PowerRTransform_get_rmax.argtypes = [c_void_p]

cext.PowerRTransform_get_power.restype = c_double
cext.PowerRTransform_get_power.argtypes = [c_void_p]

cext.ChebyshevRTransform_new.restype = c_void_p
cext.ChebyshevRTransform_new.argtypes = [c_double, c_int]

cext.ChebyshevRTransform_get_radii.restype = c_double
cext.ChebyshevRTransform_get_radii.argtypes = [c_void_p]

# Cell
cext.Cell_new.restype = c_void_p
cext.Cell_new.argtypes = [c_void_p, c_int]

cext.Cell_del.restype = None
cext.Cell_del.argtypes = [c_void_p]

# int1d
cext.compute_cubic_spline_int_weights.restype = None
cext.compute_cubic_spline_int_weights.argtypes = [a1df, c_int]

# ODE solver
cext.build_ode2.restype = None
cext.build_ode2.argtypes = [
    a1df, a1df, a1df, a1df, a2df, a1df, c_long]
