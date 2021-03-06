# Copyright (C) 2015, 2016, 2018 Jan Blechta
#
# This file is part of dolfin-tape.
#
# dolfin-tape is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dolfin-tape is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with dolfin-tape. If not, see <http://www.gnu.org/licenses/>.

"""This script finds an approximation of p-Laplace problem and then uses
its residual in W^{-1,q} to demonstrate localization result of

    J. Blechta, J. M\'alek, M. Vohral\'ik, Localization of the $W^{-1,q}$
    norm for local a posteriori efficiency. Submitted, 2016.
    URL https://hal.inria.fr/hal-01332481.
"""

from __future__ import print_function

from dolfin import *
import ufl
import numpy as np
import os

from dolfintape.poincare import poincare_friedrichs_cutoff
from dolfintape.plotting import plot_alongside, pyplot
from dolfintape.utils import mkdir_p, list_timings
from dolfintape.demo_problems import solve_p_laplace_adaptive
from dolfintape.sobolev_norm import sobolev_norm


not_working_in_parallel('This')

# UFLACS issue #49
#parameters['form_compiler']['representation'] = 'uflacs'

parameters['form_compiler']['optimize'] = True


def compute_liftings(name, p, mesh, f, exact_solution=None):
    r"""Find approximation to p-Laplace problem with rhs f,
    and compute global and local liftings of the residual.
    Return tuple (
        u_h,
        \sum_a ||\nabla r     ||_{p,\omega_a}^p \psi_a/|\omega_a|,
        \sum_a ||\nabla r^a   ||_{p,\omega_a}^p \psi_a/|\omega_a|,
        \sum_a ||\nabla(u-u_h)||_{p,\omega_a}^p \psi_a/|\omega_a|,
        \sum_a ||\sigma(\nabla u)-\sigma(\nabla u_h)||_{q,\omega_a}^q \psi_a/|\omega_a|,
        C_{cont,PF},
        ||\nabla r||_p^{p-1},
        ( 1/N \sum_a ||\nabla r_a||_p^p )^{1/q},
        ||\sigma(\nabla u) - \sigma(\nabla u_h)||_q,
        N * C_PF * r_norm_loc / r_norm_glob,
        N * r_norm_loc_PF / r_norm_glob,
        r_norm_glob / r_norm_loc,
        flux_err / r_norm_glob,
    ). First five are P1 functions, the rest are numbers.
    """
    q = p/(p-1) # Dual Lebesgue exponent
    N = mesh.topology().dim() + 1 # Vertices per cell

    # Check that mesh is the coarsest one
    assert mesh.id() == mesh.root_node().id()

    # Get Galerkin approximation of p-Laplace problem -\Delta_p u = f
    log(25, 'Computing residual of p-Laplace problem')
    V = FunctionSpace(mesh, 'Lagrange', 1)
    criterion = lambda u_h, Est_h, Est_eps, Est_tot, Est_up: Est_eps <= 1e-6*Est_tot
    u = solve_p_laplace_adaptive(p, criterion, V, f,
                                 zero(mesh.geometry().dim()), exact_solution)

    # Plot exact solution, approximation and error
    plot_error(exact_solution, u, name)

    # p-Laplacian flux of u
    S = inner(grad(u), grad(u))**(0.5*Constant(p)-1.0) * grad(u)

    # Compute cell-wise norm of flux
    if exact_solution is not None:
        u_ex = exact_solution
    else:
        warning("Don't have exact solution for computation of flux error. Assuming zero.")
        u_ex = Constant(0)
    S_ex = ufl.replace(S, {u: u_ex})
    S_ex = inner(grad(exact_solution), grad(exact_solution))**(0.5*Constant(p)-1.0) * grad(exact_solution)
    mesh_fine = u.function_space().mesh()
    flux_err = ((S - S_ex)**2)**Constant(0.5*q)
    flux_err_fine, flux_err_coarse = compute_cellwise_norm(flux_err, mesh_fine)

    # Distribute cell-wise flux error to patches
    flux_err_p1 = distribute_p0_to_p1(flux_err_coarse, Function(V))

    # Sanity check and logging
    flux_err = sobolev_norm(S_ex-S, q, k=0)
    assert np.isclose(assemble(flux_err_fine  *dx), flux_err**q)
    assert np.isclose(assemble(flux_err_coarse*dx), flux_err**q)
    assert np.isclose(assemble(flux_err_p1    *dx), flux_err**q)
    info_blue(r"||\sigma(u)-\sigma(u_h)||_q^q = %g" % flux_err**q)

    # Global lifting of W^{-1, p'} functional R = f + div(S)
    u.set_allow_extrapolation(True) # Needed hack
    r_glob = compute_global_lifting(p, mesh, f, S)
    u.set_allow_extrapolation(False)

    # Compute cell-wise norm of global lifting
    dr_glob_fine, dr_glob_coarse = compute_cellwise_grad(r_glob, p)

    # Distribute cell-wise norms of global lifting to patches
    dr_glob_p1 = distribute_p0_to_p1(dr_glob_coarse, Function(V))

    # Norm of global lifting, equal to norm of residual
    r_norm_glob = sobolev_norm(r_glob, p)**(p/q)

    # Sanity check
    assert np.isclose(assemble(dr_glob_fine  *dx), r_norm_glob**q)
    assert np.isclose(assemble(dr_glob_coarse*dx), r_norm_glob**q)
    assert np.isclose(assemble(dr_glob_p1    *dx), r_norm_glob**q)

    # Compute local liftings
    r_norm_loc, r_norm_loc_PF, r_loc_p1 = compute_local_liftings(p, V, f, S)

    # Compute energy error
    if exact_solution:
        ee_fine, ee_coarse = compute_cellwise_grad(exact_solution-u, p,
                                 mesh_fine=u.function_space().mesh())
        ee_p1 = distribute_p0_to_p1(ee_coarse, Function(V))
        ee = sobolev_norm(exact_solution-u, p)
        assert np.isclose(assemble(ee_fine  *dx), ee**p)
        assert np.isclose(assemble(ee_coarse*dx), ee**p)
        assert np.isclose(assemble(ee_p1    *dx), ee**p)
        info_blue(r"||\nabla(u-u_h)||_p^p = %g" % ee**p)
    else:
        ee_p1 = None


    # Check effectivity of localization estimates
    C_PF = poincare_friedrichs_cutoff(mesh, p)
    ratio_a = ( N * C_PF * r_norm_loc ) / r_norm_glob
    ratio_b = r_norm_glob / r_norm_loc
    ratio_a_PF = ( N * r_norm_loc_PF ) / r_norm_glob
    ratio_c = flux_err / r_norm_glob
    assert ratio_a >= 1.0 and ratio_b >= 1.0
    assert ratio_a_PF >= 1.0
    assert ratio_c >= 1.0

    # Report
    info_blue(r"||\nabla r||_p^{p-1} = %g, ( 1/N \sum_a ||\nabla r_a||_p^p )^{1/q} = %g"
              % (r_norm_glob, r_norm_loc))
    info_blue("C_{cont,PF} = %g" %  C_PF)
    info_green("(4.8a) ok: rhs/lhs = %g >= 1" % ratio_a)
    info_green("(4.8b) ok: rhs/lhs = %g >= 1" % ratio_b)
    info_green("ratio_c = %g >= 1" % ratio_c)

    return u, dr_glob_p1, r_loc_p1, ee_p1, flux_err_p1, \
        C_PF, r_norm_glob, r_norm_loc, flux_err, ratio_a, ratio_a_PF, ratio_b, ratio_c


def compute_global_lifting(p, mesh, f, S):
    """Return global lifting of
        R = f + div S
    """
    log(25, 'Computing global lifting of the resiual')

    # Use higher order space and better quadrature
    V_high = FunctionSpace(mesh, 'Lagrange', 2)
    parameters['form_compiler']['quadrature_degree'] = 8

    # Compute lifting adaptively
    criterion = lambda u_h, Est_h, Est_eps, Est_tot, Est_up: \
        Est_eps <= 1e-2*Est_tot and Est_tot <= 1e-2*sobolev_norm(u_h, p)**(p-1.0)
    r_glob = solve_p_laplace_adaptive(p, criterion, V_high, f, S,
            solver_parameters={"newton_solver": {"linear_solver": "mumps"}})

    # Rollback side-effect
    parameters['form_compiler']['quadrature_degree'] = -1

    return r_glob


def compute_local_liftings(p, P1, f, S):
    """Compute local liftings of
        R = f + div S
    and return global norm of liftings and P1 representation
    of patch-wise norms."
    """
    assert P1.ufl_element().family() == 'Lagrange'
    assert P1.ufl_element().degree() == 1
    assert P1.ufl_element().value_shape() == ()
    mesh = P1.mesh()

    r_loc_p1 = Function(P1)
    r_loc_p1_dofs = r_loc_p1.vector()
    v2d = vertex_to_dof_map(P1)
    r_norm_loc = 0.0
    r_norm_loc_PF = 0.0

    # Adjust verbosity
    old_log_level = get_log_level()
    set_log_level(WARNING)
    prg = Progress('Solving local liftings on patches', mesh.num_vertices())

    cf = CellFunction('size_t', mesh)
    for v in vertices(mesh):
        # Prepare submesh covering a patch
        cf.set_all(0)
        for c in cells(v):
            cf[c] = 1
        submesh = SubMesh(mesh, cf, 1)

        # Compute p-Laplace lifting on the patch using higher order element
        V_loc = FunctionSpace(submesh, 'Lagrange', 4)
        criterion = lambda u_h, Est_h, Est_eps, Est_tot, Est_up: \
            Est_eps <= 1e-2*Est_tot and Est_tot <= 1e-2*sobolev_norm(u_h, p)**(p-1.0)
        parameters['form_compiler']['quadrature_degree'] = 8
        r = solve_p_laplace_adaptive(p, criterion, V_loc, f, S)
        parameters['form_compiler']['quadrature_degree'] = -1

        # Compute local norm of residual
        r_norm_loc_a = sobolev_norm(r, p)**p
        r_norm_loc += r_norm_loc_a
        r_norm_loc_PF += r_norm_loc_a * poincare_friedrichs_cutoff(v, p)**(p/(p-1))
        scale = 1.0 / sum(c.volume() for c in cells(v))
        r_loc_p1_dofs[v2d[v.index()]] = r_norm_loc_a * scale
        log(18, r"||\nabla r_a||_p = %g" % r_norm_loc_a**(1.0/p))

        # Advance progress bar
        set_log_level(PROGRESS)
        prg += 1
        set_log_level(WARNING)

    # Recover original verbosity
    set_log_level(old_log_level)

    # Scale by 1/N and take q-root of sum finally
    r_norm_loc /= mesh.topology().dim() + 1
    r_norm_loc_PF /= mesh.topology().dim() + 1
    q = p/(p-1)
    r_norm_loc **= 1.0/q
    r_norm_loc_PF **= 1.0/q

    # Sanity check
    e_norm = assemble(r_loc_p1*dx)
    assert np.isclose(e_norm, r_norm_loc**q)

    return r_norm_loc, r_norm_loc_PF, r_loc_p1


def compute_cellwise_grad(r, p, mesh_fine=None):
    r"""Return fine and coarse P0 functions representing cell-wise
    L^p norm of grad(r), i.e. functions having values

        ||\nabla r||_{p, K} \frac{1}{|K|}

    on cell K. First (fine) function is defined on (fine) cells of r;
    second (coarse) function is reduction to coarse mesh (obtained from
    root node of hierarchical chain of r, if any).

    Scaling is chosen such that

        \int D = ||\nabla r||_p^p

    for both returned functions D.
    """
    # Compute desired quantity accurately on fine mesh
    mesh_fine = mesh_fine or r.function_space().mesh()
    dr_fine, dr_coarse = compute_cellwise_norm((grad(r)**2)**Constant(0.5*p), mesh_fine)
    return dr_fine, dr_coarse


def compute_cellwise_norm(expression, mesh_fine):
    # Compute desired quantity accurately on fine mesh
    P0_fine = FunctionSpace(mesh_fine, 'Discontinuous Lagrange', 0)
    dr_fine = project(expression, P0_fine)

    # Special case
    mesh_coarse = mesh_fine.root_node()
    if mesh_fine.id() == mesh_coarse.id():
        return dr_fine, dr_fine.copy(deepcopy=True)

    # Compute parent cells from finest to coarsest
    mesh = mesh_fine
    tdim = mesh.topology().dim()
    parent_cells = slice(None)
    while mesh.parent():
        parent_cells = mesh.data().array('parent_cell', tdim)[parent_cells]
        mesh = mesh.parent()

    # Sanity check
    assert parent_cells.shape == (mesh_fine.num_cells(),)
    assert parent_cells.ptp() + 1 == mesh_coarse.num_cells()

    # Init coarse quantity
    P0_coarse = FunctionSpace(mesh_coarse, 'Discontinuous Lagrange', 0)
    dr_coarse = Function(P0_coarse)

    # Fetch needed objects to speed-up the hot loop
    dofs_fine = P0_fine.dofmap().cell_dofs
    dofs_coarse = P0_coarse.dofmap().cell_dofs
    x_fine = dr_fine.vector()
    x_coarse = dr_coarse.vector()

    # Reduce fine to coarse
    for c in cells(mesh_fine):
        i = c.index()
        scale = c.volume()/Cell(mesh_coarse, parent_cells[i]).volume()
        x_coarse[dofs_coarse(parent_cells[i])] += scale * x_fine[dofs_fine(i)]

    return dr_fine, dr_coarse


def distribute_p0_to_p1(f, out=None):
    r"""Distribute P0 function to P1 function s.t.

        g = \sum_{a \in vertices} \sum_{K \ni a} f_K |K|/|\omega_a|

    Returns P1 function g.
    """
    P0 = f.function_space()
    mesh = P0.mesh()
    if out is None:
        P1 = FunctionSpace(mesh, 'Lagrange', 1)
        out = Function(P1)
    else:
        P1 = out.function_space()

    # Sanity check
    assert P0.ufl_element().family() == 'Discontinuous Lagrange'
    assert P0.ufl_element().degree() == 0
    assert P0.ufl_element().value_shape() == ()
    assert P1.ufl_element().family() == 'Lagrange'
    assert P1.ufl_element().degree() == 1
    assert P1.ufl_element().value_shape() == ()

    vec = out.vector()
    x = np.ndarray(mesh.geometry().dim()) # Dummy coord
    val = np.ndarray(1)
    v2d = vertex_to_dof_map(P1)

    # Collect contribution from cells
    for c in cells(mesh):
        f.eval(val, x, c, c)
        for v in vertices(c):
            # FIXME: it would be faster to assemble vector just once
            vol_cell = c.volume()
            vol_patch = sum(c.volume() for c in cells(v))
            dof = v2d[v.index()]
            vec[dof] = vec[dof][0] + val[0]*vol_cell/vol_patch

    return out


def save_functions(uh, glob, loc, ee, fe, prefix):
    path = "./"
    mkdir_p(path)

    naming = {"uh": uh, "glob": glob, "loc": loc, "ee": ee, "fe": fe}

    for name, func in naming.iteritems():
        filename = os.path.join(path, prefix + "_" + name + ".h5")
        mesh = func.function_space().mesh()
        with HDF5File(mesh.mpi_comm(), filename, 'w') as f:
            f.write(mesh, "mesh")
            f.write(func, "func")


def plot_liftings(glob, loc, ee, fe, prefix):
    path = "./"
    mkdir_p(path)

    # Plot global and local lifting norms on patches
    plot_alongside(glob, loc, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_r_f.pdf"))
    plot_alongside(glob, loc, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_r_g.pdf"))
    plot_alongside(glob, loc, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_r_w.pdf"))
    plot_alongside(glob, loc, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_r_c.pdf"))

    # Plot global lifting norm anf flux error on patches
    plot_alongside(glob, fe, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_gf_f.pdf"))
    plot_alongside(glob, fe, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_gf_g.pdf"))
    plot_alongside(glob, fe, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_gf_w.pdf"))
    plot_alongside(glob, fe, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_gf_c.pdf"))

    # Plot local lifting norm anf flux error on patches
    plot_alongside(loc, fe, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_lf_f.pdf"))
    plot_alongside(loc, fe, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_lf_g.pdf"))
    plot_alongside(loc, fe, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_lf_w.pdf"))
    plot_alongside(loc, fe, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_lf_c.pdf"))

    # Plot global and local lifting and flux error norms on patches
    plot_alongside(glob, loc, fe, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_glf_f.pdf"))
    plot_alongside(glob, loc, fe, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_glf_g.pdf"))
    plot_alongside(glob, loc, fe, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_glf_w.pdf"))
    plot_alongside(glob, loc, fe, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_glf_c.pdf"))

    # Plot global lifting and energy error norms on patches
    plot_alongside(glob, ee, common_cbar=False, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_ee_f.pdf"))
    plot_alongside(glob, ee, common_cbar=False, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_ee_g.pdf"))
    plot_alongside(glob, ee, common_cbar=False, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_ee_w.pdf"))
    plot_alongside(glob, ee, common_cbar=False, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_ee_c.pdf"))


def plot_error(u, uh, prefix):
    path = "./"
    mkdir_p(path)

    # Plot exact solution and approximation
    u = project(u, uh.function_space())
    plot_alongside(u, uh, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_u_f.pdf"))
    plot_alongside(u, uh, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_u_g.pdf"))
    plot_alongside(u, uh, mode="warp")
    pyplot.savefig(os.path.join(path, prefix+"_u_w.pdf"))
    plot_alongside(u, uh, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_u_c.pdf"))

    # Plot error
    e = project(u-uh, uh.function_space())
    plot_alongside(e, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_e_f.pdf"))
    plot_alongside(e, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_e_g.pdf"))
    plot_alongside(e, mode="warp")
    pyplot.savefig(os.path.join(path, prefix+"_e_w.pdf"))
    plot_alongside(e, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_e_c.pdf"))


def plot_cutoff_distribution(p, mesh, prefix):
    """Plot distribution of Poincate-Friedrichs cutoff constant
    across mesh patches
    """
    # Compute cutoff constants and assemble P1 function
    P1 = FunctionSpace(mesh, 'Lagrange', 1)
    dist = Function(P1)
    vec = dist.vector()
    v2d = vertex_to_dof_map(P1)
    for v in vertices(mesh):
        vec[v2d[v.index()]] = poincare_friedrichs_cutoff(v, p)

    # Plot
    path = "./"
    mkdir_p(path)
    plot_alongside(dist, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_cutoff_f.pdf"))
    plot_alongside(dist, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_cutoff_g.pdf"))
    plot_alongside(dist, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_cutoff_w.pdf"))


def format_result(*args):
    assert len(args) == 11
    assert isinstance(args[0], str) and len(args[0].split()) == 1
    assert isinstance(args[2], int)
    assert all(isinstance(args[i], float) for i in [1]+range(3, 11))
    print("#RESULT name, p, num_cells, C_{cont,PF}, " \
          "||E_glob||_q, ||E_loc||_q, ||E_flux||_q, Eff_(4.8a), Eff_(4.8b), Eff_(flux/glob)")
    print("RESULT %s %s %4d %.3f %.4f %.4f %.4f %.1f %.1f %.2f %.2f " % args)


def test_ChaillouSuri(p, N):
    from dolfintape.demo_problems.exact_solutions import pLaplace_ChaillouSuri

    label = 'ChaillouSuri_%s_%02d' % (p, N)

    # Fetch exact solution and rhs of p-Laplacian
    mesh = UnitSquareMesh(N, N, 'crossed')
    print("num cells %s" % mesh.num_cells())
    plot_cutoff_distribution(p, mesh, label)
    u, f = pLaplace_ChaillouSuri(p, domain=mesh, degree=4)

    # Now the heavy lifting
    result = compute_liftings(label, p, mesh, f, u)
    uh, glob, loc, ee, fe = result[0:5]

    # Report
    format_result('Chaillou--Suri', p, mesh.num_cells(), *result[5:])
    plot_liftings(glob, loc, ee, fe, label)
    save_functions(uh, glob, loc, ee, fe, label)
    list_timings(TimingClear_clear, [TimingType_wall])


def test_CarstensenKlose(p, N):
    from dolfintape.demo_problems.exact_solutions import pLaplace_CarstensenKlose
    from dolfintape.mesh_fixup import mesh_fixup
    import mshr

    label = 'CarstensenKlose_%s_%02d' % (p, N)

    # Build mesh on L-shaped domain (-1, 1)^2 \ (0, 1)*(-1, 0)
    b0 = mshr.Rectangle(Point(-1.0, -1.0), Point(1.0, 1.0))
    b1 = mshr.Rectangle(Point(0.0, -1.0), Point(1.0, 0.0))
    mesh = mshr.generate_mesh(b0 - b1, N)
    mesh = mesh_fixup(mesh)
    print("num cells %s" % mesh.num_cells())
    plot_cutoff_distribution(p, mesh, label)

    # Fetch exact solution and rhs of p-Laplacian
    u, f = pLaplace_CarstensenKlose(p=p, eps=0.0, delta=7.0/8,
                                    domain=mesh, degree=4)
    # There are some problems with quadrature element,
    # see https://bitbucket.org/fenics-project/ffc/issues/84,
    # so precompute (f, vh) for vh from P1
    f = project(f, FunctionSpace(mesh, 'Lagrange', 1),
                form_compiler_parameters=
                    {'quadrature_degree': f.ufl_element().degree()})
    f.set_allow_extrapolation(True)

    # Now the heavy lifting
    result = compute_liftings(label, p, mesh, f, u)
    uh, glob, loc, ee, fe = result[0:5]

    # Report
    format_result('Carstensen--Klose', p, mesh.num_cells(), *result[5:])
    plot_liftings(glob, loc, ee, fe, label)
    save_functions(uh, glob, loc, ee, fe, label)
    list_timings(TimingClear_clear, [TimingType_wall])


def main(argv):
    default_tests = [
            ('ChaillouSuri',   10.0,  5),
            ('ChaillouSuri',   10.0, 10),
            ('ChaillouSuri',   10.0, 15),
            ('ChaillouSuri',   10.0, 20),
            ('ChaillouSuri',    1.5,  5),
            ('ChaillouSuri',    1.5, 10),
            ('ChaillouSuri',    1.5, 15),
            ('ChaillouSuri',    1.5, 20),
            ('CarstensenKlose', 4.0,  5),
            ('CarstensenKlose', 4.0, 10),
            ('CarstensenKlose', 4.0, 15),
            ('CarstensenKlose', 4.0, 20),
        ]

    usage = """%s

usage: python %s [-h|--help] [test-name p N]

Without arguments run default test cases. Or run test case with
given value of p when given on command-line.

Default test cases:

%s
""" % (__doc__, __file__, default_tests)

    # Decrease verbosity of DOLFIN
    set_log_level(25)

    # Run all tests
    if len(argv) == 1:
        for test in default_tests:
            exec('test_%s(%s)' % test)
            return

    # Print help
    if argv[1] in ['-h', '--help']:
        print(usage)
        return

    # Now expecting 3 arguments
    if len(argv) != 4:
        print("Command-line arguments not understood!")
        print()
        print(usage)
        return 1

    # Run the selected test
    try:
        exec("tester = test_%s" % argv[1])
    except NameError:
        print ("'Test %s' does not exist!" % argv[1])
        return 1
    else:
        tester(float(argv[2]), int(argv[3]))
        return


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
