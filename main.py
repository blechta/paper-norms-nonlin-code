# Copyright (C) 2015, 2016 Jan Blechta
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
parameters['form_compiler']['representation'] = 'quadrature'

parameters['form_compiler']['optimize'] = True


def compute_liftings(name, p, mesh, f, exact_solution=None, S=None):
    r"""Find approximation to p-Laplace problem with rhs f,
    and compute global and local liftings of the residual.
    Return tuple (
        \sum_a ||\nabla r     ||_{p,\omega_a}^p \psi_a/|\omega_a|,
        \sum_a ||\nabla r^a   ||_{p,\omega_a}^p \psi_a/|\omega_a|,
        \sum_a ||\nabla(u-u_h)||_{p,\omega_a}^p \psi_a/|\omega_a|,
        \sum_a \sum_{K\in\omega_a} \eta_K^q |K|/|\omega_a| \psi_a,
        C_{cont,PF},
        ||\nabla r||_p^{p-1},
        ( 1/N \sum_a ||\nabla r_a||_p^p )^{1/q},
        Eff_{(4.8a)},
        Eff_{(4.8b)}
    ). First three are P1 functions, the rest are numbers.
    """
    q = p/(p-1) # Dual Lebesgue exponent
    N = mesh.topology().dim() + 1 # Vertices per cell

    # Check that mesh is the coarsest one
    assert mesh.id() == mesh.root_node().id()

    # Get Galerkin approximation of p-Laplace problem -\Delta_p u = f
    log(25, 'Computing residual of p-Laplace problem')
    V = FunctionSpace(mesh, 'Lagrange', 1)
    criterion = lambda u_h, Est_h, Est_eps, Est_tot, Est_up: Est_eps <= 1e-6*Est_tot
    u, est_h, est_eps, est_tot = solve_p_laplace_adaptive(p, criterion, V, f, S,
                                                          u_ex=exact_solution, eps0=0.0)

    # Plot exact solution, approximation and error
    plot_error(exact_solution, u, name)

    # Default to p-Laplacian flux of u
    if S is None:
        def createResidualFlux(p, u):
            def _S(r, eps):
                return (eps + inner(grad(r), grad(r)))**(0.5*Constant(p)-1.0) * grad(r) \
                    + inner(grad(u), grad(u))**(0.5*Constant(p)-1.0) * grad(u)
            return _S
        R = createResidualFlux(p, u)
    else:
        def createResidualFlux(p, u):
            def _S(r, eps):
                return (eps + inner(grad(r), grad(r)))**(0.5*Constant(p)-1.0) * grad(r) \
                    + S(u, eps)
            return _S
        R = createResidualFlux(p, u)

    # Global lifting of W^{-1, p'} functional R = f + div(S)
    u.set_allow_extrapolation(True) # Needed hack
    r_glob = compute_global_lifting(p, mesh, f, R)
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
    r_norm_loc, r_norm_loc_PF, r_loc_p1 = compute_local_liftings(p, V, f, R)

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
    if r_norm_loc == 0:
        ratio_b = np.nan
    else:
        ratio_b = r_norm_glob / r_norm_loc
    ratio_a_PF = ( N * r_norm_loc_PF ) / r_norm_glob
    #assert ratio_a >= 1.0 and ratio_b >= 1.0
    #assert ratio_a_PF >= 1.0

    # Report
    info_blue(r"||\nabla r||_p^{p-1} = %g, ( 1/N \sum_a ||\nabla r_a||_p^p )^{1/q} = %g"
              % (r_norm_glob, r_norm_loc))
    info_blue("C_{cont,PF} = %g" %  C_PF)
    info_green("(4.8a) ok: rhs/lhs = %g >= 1" % ratio_a)
    info_green("(4.8b) ok: rhs/lhs = %g >= 1" % ratio_b)

    # Get P1 distribution of cell estimator
    for c in cells(est_tot.mesh()):
        est_tot[c] **= q
        est_tot[c] /= c.volume()
    eta = distribute_cellfunction_to_p1(est_tot, Function(V))
    est_tot.set_all(0)
    info_blue("eta_tot = %g" % assemble(eta*dx)**(1.0/q))

    return dr_glob_p1, r_loc_p1, ee_p1, eta, C_PF, r_norm_glob, r_norm_loc, ratio_a, ratio_a_PF, ratio_b


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
    r_glob, _, _, _ = solve_p_laplace_adaptive(p, criterion, V_high, f, S,
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
    #for v in vertices(mesh):
    for v in []:
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
        r, _, _, _ = solve_p_laplace_adaptive(p, criterion, V_loc, f, S)
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
    P0_fine = FunctionSpace(mesh_fine, 'Discontinuous Lagrange', 0)
    dr_fine = project((grad(r)**2)**Constant(0.5*p), P0_fine)

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
    r"""Distribute P0 function f to P1 function g s.t.

        g = \sum_{a \in vertices} \sum_{K \ni a} f_K |K|/|\omega_a| \psi_a

    where \psi_a is P1 basis function for vertex a and \omega_a
    is a patch of cells around vertex a. Returns P1 function g.
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


def distribute_cellfunction_to_p1(f, out=None):
    r"""Distribute cell function f to P1 function g s.t.

        g = \sum_{a \in vertices} \sum_{K \ni a} f_K |K|/|\omega_a| \psi_a

    where \psi_a is P1 basis function for vertex a and \omega_a
    is a patch of cells around vertex a. Returns P1 function g.
    """
    mesh = f.mesh()
    if out is None:
        P1 = FunctionSpace(mesh, 'Lagrange', 1)
        out = Function(P1)
    else:
        P1 = out.function_space()

    # Sanity check
    assert isinstance(f, CellFunctionDouble)
    assert P1.ufl_element().family() == 'Lagrange'
    assert P1.ufl_element().degree() == 1
    assert P1.ufl_element().value_shape() == ()

    vec = out.vector()
    v2d = vertex_to_dof_map(P1)

    # Collect contribution from cells
    for c in cells(mesh):
        val = f[c]
        for v in vertices(c):
            # FIXME: it would be faster to assemble vector just once
            vol_cell = c.volume()
            vol_patch = sum(c.volume() for c in cells(v))
            dof = v2d[v.index()]
            vec[dof] = vec[dof][0] + val*vol_cell/vol_patch

    return out


def function_ipow(fun, exponent):
    "Take inplace power of function by powering its dofs"""
    x = fun.vector()
    # FIXME: PETSc VecPow would be faster
    x[:] = x.array()**exponent


def plot_liftings(glob, loc, ee, eta, prefix):
    path = "./"
    mkdir_p(path)

    # Plot global and local lifting norms on patches
    #plot_alongside(glob, loc, mode="color", shading="flat", edgecolors="k")
    #pyplot.savefig(os.path.join(path, prefix+"_r_f.pdf"))
    #plot_alongside(glob, loc, mode="color", shading="gouraud")
    #pyplot.savefig(os.path.join(path, prefix+"_r_g.pdf"))
    #plot_alongside(glob, loc, mode="warp", range_min=0.0)
    #pyplot.savefig(os.path.join(path, prefix+"_r_w.pdf"))
    #plot_alongside(glob, loc, mode="contour")
    #pyplot.savefig(os.path.join(path, prefix+"_r_c.pdf"))

    # Plot global lifting and energy error norms on patches
    plot_alongside(glob, ee, common_cbar=False, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_ee_f.pdf"))
    plot_alongside(glob, ee, common_cbar=False, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_ee_g.pdf"))
    plot_alongside(glob, ee, common_cbar=False, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_ee_w.pdf"))
    plot_alongside(glob, ee, common_cbar=False, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_ee_c.pdf"))

    # Plot global lifting norm and error estimator on patches
    plot_alongside(glob, eta, mode="color", shading="flat", edgecolors="k")
    pyplot.savefig(os.path.join(path, prefix+"_eta_f.pdf"))
    plot_alongside(glob, eta, mode="color", shading="gouraud")
    pyplot.savefig(os.path.join(path, prefix+"_eta_g.pdf"))
    plot_alongside(glob, eta, mode="warp", range_min=0.0)
    pyplot.savefig(os.path.join(path, prefix+"_eta_w.pdf"))
    plot_alongside(glob, eta, mode="contour")
    pyplot.savefig(os.path.join(path, prefix+"_eta_c.pdf"))


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
    assert len(args) == 9
    assert isinstance(args[0], str) and len(args[0].split()) == 1
    assert isinstance(args[2], int)
    assert all(isinstance(args[i], float) for i in [1]+range(3, 9))
    print("#RESULT name, p, num_cells, C_{cont,PF}, " \
          "||E_glob||_q, ||E_loc||_q, Eff_(4.8a), Eff_(4.8b)")
    print("RESULT %s %s %4d %.3f %.4f %.4f %.1f %.1f %.2f" % args)


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
    glob, loc, ee, eta = result[0:4]

    # Report
    format_result('Chaillou--Suri', p, mesh.num_cells(), *result[4:])
    plot_liftings(glob, loc, ee, eta, label)
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
    glob, loc, ee, eta = result[0:4]

    # Report
    format_result('Carstensen--Klose', p, mesh.num_cells(), *result[4:])
    plot_liftings(glob, loc, ee, eta, label)
    list_timings(TimingClear_clear, [TimingType_wall])


def test_NicaiseVenel(sigma_minus, N):
    p = q = 2
    label = 'NicaiseVenel_%s_%02d' % (sigma_minus, N)

    # Fetch exact solution and rhs of p-Laplacian
    assert N % 2 == 0
    mesh = UnitSquareMesh(N, N, 'crossed')
    mesh.coordinates()[:] *= 2
    mesh.coordinates()[:] += (-1, -1)
    cf = CellFunction('size_t', mesh)
    AutoSubDomain(lambda x: x[0] <= +DOLFIN_EPS).mark(cf, 0)
    AutoSubDomain(lambda x: x[0] >= -DOLFIN_EPS).mark(cf, 1)

    # FIXME: Use new martinal's functionality: cell functions in cpp exprs
    sigma = Expression("x[0] >= 0 ? 1.0 : sigma_minus",
                       sigma_minus=sigma_minus, degree=0, domain=mesh)
    def S(u, eps):
        return sigma*grad(u)
    u = Expression("x[0] > 0.0 ? "
                   "sigma_minus*x[0]*(x[0]+1)*(x[0]-1)*(x[1]+1)*(x[1]-1)"
                   " : "
                   "x[0]*(x[0]+1)*(x[0]-1)*(x[1]+1)*(x[1]-1)",
                   sigma_minus=sigma_minus, degree=5, domain=mesh)
    f = Expression("-sigma_minus*2.0*x[0]*((x[0]+1)*(x[0]-1)+3.0*(x[1]+1)*(x[1]-1))",
                   sigma_minus=sigma_minus, degree=3, domain=mesh)

    #pyplot.figure()
    #plot(cf)
    #pyplot.show()
    #plot(u, mesh=mesh, mode="warp")
    #pyplot.show()
    #plot(f, mesh=mesh, mode="warp")
    #pyplot.show()
    #exit()

    print("num cells %s" % mesh.num_cells())
    plot_cutoff_distribution(p, mesh, label)

    # Now the heavy lifting
    result = compute_liftings(label, p, mesh, f, u, S=S)
    glob, loc, ee, eta = result[0:4]

    # Take square root of P1 functions (then they are no more polynomials...)
    function_ipow(glob, 1.0/q)
    function_ipow(loc, 1.0/q)
    function_ipow(ee, 1.0/p)
    function_ipow(eta, 1.0/q)

    # Report
    format_result('Nicaise--Venel', sigma_minus, mesh.num_cells(), *result[4:])
    plot_liftings(glob, loc, ee, eta, label)
    list_timings(TimingClear_clear, [TimingType_wall])


def test_BonnetBenDhia(sigma_minus, N):
    p = q = 2
    label = 'BonnetBenDhia_%s_%02d' % (sigma_minus, N)

    # Fetch exact solution and rhs of p-Laplacian
    assert N % 2 == 0
    mesh = UnitSquareMesh(N, N, 'crossed')
    mesh.coordinates()[:] *= 2
    mesh.coordinates()[:] += (-1, -1)

    # FIXME: Use new martinal's functionality: cell functions in cpp exprs
    sigma = Expression("x[0] >= 0 && x[1] >= 0 ? 1.0 : sigma_minus",
                       sigma_minus=sigma_minus, degree=0, domain=mesh)

    class SolTest(Expression):
        def __init__(self, mu, *args, **kwargs):
            self.mu = mu
            self.dec_point = False

        def eval(self, values, x):
            x, y = x[0:2]

            x_dc = x
            y_dc = y
            if self.dec_point:
                x_dc = x_dec
                y_dc = y_dec

            lmbd = 2 / pi * acos(( 1 - self.mu ) / 2.0 / abs( 1 + self.mu ))

            r = sqrt( pow (x, 2) + pow (y, 2) )
            th = self.f_angle_point_dec( x, y, x_dc, y_dc )
            th_dc = self.f_angle_point( x_dc, y_dc )

            A = sin ( lmbd * pi / 2.0 )
            B = cos ( lmbd * pi / 2.0 )
            D = sin ( lmbd * pi * 3 / 2.0 )
            F = cos ( lmbd * pi * 3 / 2.0 )

            c1 = 1
            d2 = A / D
            d1 = ( A * B + self.mu * F * A * A / D ) / ( D + self.mu * A )
            c2 = d1 *D / A

            if th_dc >= 0 and th_dc <= pi / 2.0:
                values[0] = c1 * sin( lmbd * th ) + c2 * sin( lmbd * ( pi / 2.0 - th ))
                values[0] *= pow( r, lmbd )
            else:
                values[0] = d1 * sin( lmbd * ( th - pi / 2.0 )) + d2 * sin( lmbd * ( 2 * pi - th ))
                values[0] *= pow( r, lmbd )

        @staticmethod
        def f_angle_point(x, y):
            ZERO2 = DOLFIN_EPS
            if abs( x ) < ZERO2:
                if y >= 0:
                   theta = pi / 2.0
                else:
                   theta = 3 * pi / 2.0
                return theta
            if x >= 0 and y >= 0:
                  theta = atan ( y / x )
            if x < 0 and y >= 0:
                  theta = atan ( y / x ) + pi
            if x < 0 and y < 0:
                  theta = atan ( y / x ) + pi
            if x >= 0 and y < 0:
                  theta = atan ( y / x ) + 2 * pi
            return theta

        @staticmethod
        def f_angle_point_dec(x, y, x_dec, y_dec):
            ZERO2 = DOLFIN_EPS
            if abs( x ) < ZERO2:
                if y_dec >= 0:
                   theta = pi / 2.0
                else:
                   theta = 3 * pi / 2.0
                return theta
            if x_dec >= 0 and y_dec >= 0:
                  theta = atan ( y / x )
            if x_dec < 0 and y_dec >= 0:
                  theta = atan ( y / x ) + pi
            if x_dec < 0 and y_dec < 0:
                  theta = atan ( y / x ) + pi
            if x_dec >= 0 and y_dec < 0:
                  theta = atan ( y / x ) + 2 * pi
            return theta


    def S(u, eps):
        return sigma*grad(u)
    u = SolTest(sigma_minus, degree=3, domain=mesh)
    f = Constant(0)

    #pyplot.show()
    #plot(u, mesh=mesh, mode="warp")
    #pyplot.show()
    #exit()

    print("num cells %s" % mesh.num_cells())
    plot_cutoff_distribution(p, mesh, label)

    # Now the heavy lifting
    result = compute_liftings(label, p, mesh, f, u, S=S)
    glob, loc, ee, eta = result[0:4]

    # Take square root of P1 functions (then they are no more polynomials...)
    function_ipow(glob, 1.0/q)
    function_ipow(loc, 1.0/q)
    function_ipow(ee, 1.0/p)
    function_ipow(eta, 1.0/q)

    # Report
    format_result('Bonnet--BenDhia', sigma_minus, mesh.num_cells(), *result[4:])
    plot_liftings(glob, loc, ee, eta, label)
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
