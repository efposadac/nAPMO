import os
import numpy as np
import numpy.ctypeslib as npct

from ctypes import *

from napmo.system.cext import napmo_library as cext
from napmo.system.timer import Timer
from napmo.system.input_parser import InputParser
from napmo.system.input_parser import raise_exception
from napmo.system.atomic_element import AtomicElement
from napmo.system.elementary_particle import ElementaryParticle
from napmo.system.molecular_system import MolecularSystem
from napmo.system.napmo_system import NAPMO

from napmo.gto.basis_set import BasisSet
from napmo.gto.primitive_gaussian import PrimitiveGaussian
from napmo.gto.contracted_gaussian import ContractedGaussian

from napmo.grids.radial import RadialGrid
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
from napmo.grids.extrapolation import PotentialExtrapolation
from napmo.grids.lebedev import lebedev_get_order
from napmo.grids.poisson_solver import poisson_solver

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
from napmo.data.constants import BOHR_TO_ANGSTROM
from napmo.data.constants import PROTON_MASS
from napmo.data.constants import NEUTRON_MASS
from napmo.data.constants import SPIN_ELECTRON

from napmo.scf.convergence import Convergence
from napmo.scf.scf import SCF

from napmo.solver.hf_solver import HF
from napmo.solver.dft_solver import DFT

from napmo.wavefunction.psi_base import PSIB
from napmo.wavefunction.psi_analytic import PSIA
from napmo.wavefunction.psi_numeric import PSIN
from napmo.wavefunction.psi_hybrid import PSIH
from napmo.wavefunction.psi_optimization import PSIO
from napmo.wavefunction.nkinetic import compute_kinetic
from napmo.wavefunction.nnuclear import compute_nuclear
from napmo.wavefunction.ntwobody import compute_coulomb
from napmo.wavefunction.ndpsi import compute_dpsi
# from napmo.wavefunction.nexccor import compute_exccor


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
a1dli = npct.ndpointer(dtype=np.int64, ndim=1, flags='CONTIGUOUS')
aptr = npct.ndpointer(c_void_p, flags="C_CONTIGUOUS")
d1pp = npct.ndpointer(dtype=np.uintp, ndim=1, flags='C')
# C functions

# Libint
cext.LibintInterface_new.restype = c_void_p
cext.LibintInterface_new.argtypes = [c_int]

cext.LibintInterface_del.restype = None

cext.LibintInterface_add_basis.restype = None
cext.LibintInterface_add_basis.argtypes = [
    c_void_p, c_void_p]

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

cext.LibintInterface_diis_new.restype = c_void_p
cext.LibintInterface_diis_new.argtypes = [c_int]

cext.LibintInterface_diis.restype = None
cext.LibintInterface_diis.argtypes = [c_void_p, c_void_p]

# Wavefunction
cext.wavefunction_guess_hcore.restype = None
cext.wavefunction_guess_hcore.argtypes = [c_void_p]

cext.wavefunction_transformation_matrix.restype = None
cext.wavefunction_transformation_matrix.argtypes = [c_void_p]

cext.wavefunction_compute_coefficients.restype = None
cext.wavefunction_compute_coefficients.argtypes = [
    c_void_p]

cext.wavefunction_compute_density.restype = None
cext.wavefunction_compute_density.argtypes = [
    c_void_p]

cext.wavefunction_compute_energy.restype = None
cext.wavefunction_compute_energy.argtypes = [
    c_void_p]

cext.wavefunction_iterate.restype = None
cext.wavefunction_iterate.argtypes = [
    c_void_p]

cext.wavefunction_compute_2body_matrix.restype = None
cext.wavefunction_compute_2body_matrix.argtypes = [
    c_void_p, c_void_p]

# cext.wavefunction_compute_exccor_matrix.restype = None
# cext.wavefunction_compute_exccor_matrix.argtypes = [
#     c_void_p]

# NWavefunction
cext.nwavefunction_compute_density_from_dm.restype = None
cext.nwavefunction_compute_density_from_dm.argtypes = [
    c_void_p, c_void_p, a2df, a1df, c_double, a1df]

cext.nwavefunction_compute_2body_matrix_atm.restype = None
cext.nwavefunction_compute_2body_matrix_atm.argtypes = [
    c_void_p, c_void_p, a2df, a1df, a2df]

cext.nwavefunction_compute_2body_matrix_mol.restype = None
cext.nwavefunction_compute_2body_matrix_mol.argtypes = [
    c_void_p, c_void_p, a2df, a1df, a2df]

cext.nwavefunction_compute_coupling.restype = None
cext.nwavefunction_compute_coupling.argtypes = [
    c_void_p, c_void_p, a2df, a1df, a2df]

cext.nwavefunction_compute_exccor_matrix.restype = None
cext.nwavefunction_compute_exccor_matrix.argtypes = [
    c_void_p, c_void_p, a2df, a1df, a1df]

cext.nwavefunction_compute_cor2species_matrix.restype = None
cext.nwavefunction_compute_cor2species_matrix.argtypes = [
    c_void_p, c_void_p, c_void_p, a2df, a1df, a1df, a1df]

# PrimitiveGaussian
cext.PrimitiveGaussian_new.restype = c_void_p
cext.PrimitiveGaussian_new.argtypes = [
    a1di, a1df, c_double, c_double]

cext.PrimitiveGaussian_compute.restype = c_double
cext.PrimitiveGaussian_compute.argtypes = [
    c_void_p, a2df, a1df, c_int]

cext.PrimitiveGaussian_overlap.restype = c_double
cext.PrimitiveGaussian_overlap.argtypes = [c_void_p, c_void_p]

cext.PrimitiveGaussian_get_l.restype = None
cext.PrimitiveGaussian_get_l.argtypes = [c_void_p, a1di]

cext.PrimitiveGaussian_get_origin.restype = None
cext.PrimitiveGaussian_get_origin.argtypes = [c_void_p, a1df]

cext.PrimitiveGaussian_get_zeta.restype = c_double
cext.PrimitiveGaussian_get_zeta.argtypes = [c_void_p]

cext.PrimitiveGaussian_get_coeff.restype = c_double
cext.PrimitiveGaussian_get_coeff.argtypes = [c_void_p]

cext.PrimitiveGaussian_get_norma.restype = c_double
cext.PrimitiveGaussian_get_norma.argtypes = [c_void_p]

# ContractedGaussian
cext.ContractedGaussian_new.restype = c_void_p
cext.ContractedGaussian_new.argtypes = [aptr, c_int]

cext.ContractedGaussian_get_nprim.restype = c_int
cext.ContractedGaussian_get_nprim.argtypes = [c_void_p]

cext.ContractedGaussian_get_l.restype = None
cext.ContractedGaussian_get_l.argtypes = [c_void_p, a1di]

cext.ContractedGaussian_get_origin.restype = None
cext.ContractedGaussian_get_origin.argtypes = [c_void_p, a1df]

cext.ContractedGaussian_compute.restype = None
cext.ContractedGaussian_compute.argtypes = [
    c_void_p, a2df, a1df, c_int]

cext.ContractedGaussian_overlap.restype = c_double
cext.ContractedGaussian_overlap.argtypes = [c_void_p, c_void_p]

cext.ContractedGaussian_get_norma.restype = c_double
cext.ContractedGaussian_get_norma.argtypes = [c_void_p]


# BasisSet
cext.BasisSet_new_empty.restype = c_void_p

cext.BasisSet_new.restype = c_void_p
cext.BasisSet_new.argtypes = [aptr, c_int]

cext.BasisSet_compute.restype = None
cext.BasisSet_compute.argtypes = [
    c_void_p, a2df, a2df, c_int]

cext.BasisSet_update.restype = None
cext.BasisSet_update.argtypes = [
    c_void_p, c_void_p]

cext.BasisSet_get_nbasis.restype = c_int
cext.BasisSet_get_nbasis.argtypes = [c_void_p]

cext.BasisSet_get_max_l.restype = c_int
cext.BasisSet_get_max_l.argtypes = [c_void_p]

cext.BasisSet_get_max_nprim.restype = c_int
cext.BasisSet_get_max_nprim.argtypes = [c_void_p]


# Angular grid
cext.AngularGrid_new.restype = c_void_p
cext.AngularGrid_new.argtypes = [c_int]

cext.AngularGrid_del.restype = None
cext.AngularGrid_del.argtypes = [c_void_p]

cext.AngularGrid_spherical.restype = None
cext.AngularGrid_spherical.argtypes = [c_void_p]

cext.AngularGrid_spherical_expansion.restype = None
cext.AngularGrid_spherical_expansion.argtypes = [
    c_void_p, c_int, c_int, a1df, a2df]

cext.AngularGrid_eval_expansion.restype = None
cext.AngularGrid_eval_expansion.argtypes = [
    c_void_p, c_int, c_int, a2df, a1df]

cext.AngularGrid_integrate.restype = c_double
cext.AngularGrid_integrate.argtypes = [c_void_p, c_int, a1df]

cext.AngularGrid_get_lorder.restype = c_int
cext.AngularGrid_get_lorder.argtypes = [c_void_p]

cext.AngularGrid_get_points.restype = POINTER(c_double)
cext.AngularGrid_get_points.argtypes = [c_void_p]

cext.AngularGrid_get_weights.restype = POINTER(c_double)
cext.AngularGrid_get_weights.argtypes = [c_void_p]

# Radial grid
cext.RadialGrid_new.restype = c_void_p
cext.RadialGrid_new.argtypes = [c_void_p, c_double]

cext.RadialGrid_del.restype = None
cext.RadialGrid_del.argtypes = [c_void_p]

cext.RadialGrid_integrate.restype = c_double
cext.RadialGrid_integrate.argtypes = [c_void_p, c_int, a1df]

cext.RadialGrid_get_size.restype = c_int
cext.RadialGrid_get_size.argtypes = [c_void_p]

cext.RadialGrid_get_radii.restype = c_double
cext.RadialGrid_get_radii.argtypes = [c_void_p]

cext.RadialGrid_get_points.restype = POINTER(c_double)
cext.RadialGrid_get_points.argtypes = [c_void_p]

cext.RadialGrid_get_weights.restype = POINTER(c_double)
cext.RadialGrid_get_weights.argtypes = [c_void_p]


# Atomic grid
cext.AtomicGrid_new.restype = c_void_p
cext.AtomicGrid_new.argtypes = [c_void_p, c_void_p, a1df]

cext.AtomicGrid_del.restype = None
cext.AtomicGrid_del.argtypes = [c_void_p]

cext.AtomicGrid_integrate.restype = POINTER(c_double)
cext.AtomicGrid_integrate.argtypes = [c_void_p, c_int, c_int, c_int, a1df]

cext.AtomicGrid_get_size.restype = c_int
cext.AtomicGrid_get_size.argtypes = [c_void_p]

cext.AtomicGrid_get_radii.restype = c_double
cext.AtomicGrid_get_radii.argtypes = [c_void_p]

cext.AtomicGrid_get_origin.restype = POINTER(c_double)
cext.AtomicGrid_get_origin.argtypes = [c_void_p]

cext.AtomicGrid_get_points.restype = POINTER(c_double)
cext.AtomicGrid_get_points.argtypes = [c_void_p]

cext.AtomicGrid_get_weights.restype = POINTER(c_double)
cext.AtomicGrid_get_weights.argtypes = [c_void_p]

# Becke grid
cext.BeckeGrid_new.restype = c_void_p
cext.BeckeGrid_new.argtypes = [aptr]

cext.BeckeGrid_del.restype = None
cext.BeckeGrid_del.argtypes = [c_void_p]

cext.BeckeGrid_integrate.restype = c_double
cext.BeckeGrid_integrate.argtypes = [c_void_p, a1df]

cext.BeckeGrid_get_ncenter.restype = c_int
cext.BeckeGrid_get_ncenter.argtypes = [c_void_p]

cext.BeckeGrid_get_size.restype = c_int
cext.BeckeGrid_get_size.argtypes = [c_void_p]

cext.BeckeGrid_get_radii.restype = POINTER(c_double)
cext.BeckeGrid_get_radii.argtypes = [c_void_p]

cext.BeckeGrid_get_points.restype = POINTER(c_double)
cext.BeckeGrid_get_points.argtypes = [c_void_p]

cext.BeckeGrid_get_weights.restype = POINTER(c_double)
cext.BeckeGrid_get_weights.argtypes = [c_void_p]

cext.BeckeGrid_get_origin.restype = POINTER(c_double)
cext.BeckeGrid_get_origin.argtypes = [c_void_p]

cext.BeckeGrid_get_becke_weights.restype = POINTER(c_double)
cext.BeckeGrid_get_becke_weights.argtypes = [c_void_p]

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

cext.PotentialExtrapolation_new.restype = c_void_p
cext.PotentialExtrapolation_new.argtypes = [c_int]


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

# Multipolar expansion
cext.dot_multi_moments.restype = None
cext.dot_multi_moments.argtypes = [
    c_long, c_long, d1pp, a2df, a1df, c_long, c_long, a1dli, a2df, c_long]
