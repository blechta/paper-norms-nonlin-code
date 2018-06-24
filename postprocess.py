# Copyright (C) 2018 Jan Blechta
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

from __future__ import print_function

import matplotlib
#matplotlib.rc('font', family='serif', serif=['computer modern roman'])
#matplotlib.rc('font', family='serif')
matplotlib.rc('mathtext', fontset='cm')
matplotlib.use("agg")

from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d

from dolfin import *

import os


def read_p1_function(prefix, name):
    path = "./"
    filename = os.path.join(path, prefix + "_" + name + ".h5")

    mesh = Mesh(mpi_comm_self())
    with HDF5File(mesh.mpi_comm(), filename, 'r') as f:
        f.read(mesh, "mesh", True)

    V = FunctionSpace(mesh, "P", 1)

    func = Function(V)
    with HDF5File(mesh.mpi_comm(), filename, 'r') as f:
        f.read(func, "func")

    return func


def plot_errors(g, l, f, e, label):
    path = "./"

    # Prepare figure
    fig = pyplot.figure(figsize=(10, 8))

    # Plot l, g, f
    functions = [
        (1, g, r"$\epsilon_\mathrm{glob}^q$"),
        (2, l, r"$\epsilon_\mathrm{loc}^q$"),
        (3, f, r"$\epsilon_\mathrm{flux}^q$"),
    ]
    _plot_subplots(2, 2, functions, [0.85, 0.55, 0.02, 0.40])
    frame = pyplot.Polygon((
        ( 92,  16),
        (332,  16),
        (332, 296),
        (720, 296),
        (720, 564),
        ( 92, 564),
    ), fill=False)
    fig.patches.append(frame)

    # Plot e
    functions = [
        (4, e, r"$\epsilon_\mathrm{en}^p$"),
    ]
    _plot_subplots(2, 2, functions, [0.85, 0.05, 0.02, 0.40])
    frame = pyplot.Rectangle((356, 16), 364, 264, fill=False)
    fig.patches.append(frame)

    # Save as PDF
    pyplot.savefig(os.path.join(path, label+"_glfe_w.pdf"))


def plot_effectivities(g, l, f, label):
    path = "./"

    # Prepare figure
    fig = pyplot.figure(figsize=(6, 12))

    # Plot effectivities g/l, f/l, f/g
    functions = [
        (1, g, l, r"\frac{\epsilon_\mathrm{glob}^q}{\epsilon_\mathrm{loc }^q}", 1.0),
        (2, f, l, r"\frac{\epsilon_\mathrm{flux}^q}{\epsilon_\mathrm{loc }^q}", 1.0),
        (3, f, g, r"\frac{\epsilon_\mathrm{flux}^q}{\epsilon_\mathrm{glob}^q}", 0.0),
    ]

    for i, f1, f2, title, range_min in functions:
        #assert(f1.function_space() == f2.function_space())
        assert(f1.function_space().ufl_element() == f2.function_space().ufl_element())
        assert(f1.function_space().mesh().hash() == f2.function_space().mesh().hash())
        eff = Function(f1.function_space())
        eff.vector()[:] = f1.vector().array() / f2.vector().array()
        print(r"\min_\Omega {} = {}".format(title, eff.vector().array().min()))
        _plot_subplots(3, 1, [(i, eff, r"${}$".format(title))],
                       [0.85, 0.27*(3-i)+0.10, 0.02, 0.20],
                       range_min=range_min)

    # Save as PDF
    pyplot.savefig(os.path.join(path, label+"_eff_w.pdf"))


def _plot_subplots(nrows, ncols, functions, cbar_rect, range_min=None):
    # Extract common range
    range_min = range_min or 0.0
    range_max = max(item[1].vector().max() for item in functions)

    # Create subplots
    for i, func, t in functions:
        sp = pyplot.subplot(nrows, ncols, i, projection="3d")
        p = plot(func, backend="matplotlib", mode="warp", title=t,
                 range_min=range_min, range_max=range_max)
        ax = pyplot.gca(projection="3d")
        ax.set_zlim(range_min, range_max)

    # Create common colorbar
    pyplot.subplots_adjust(right=0.8)
    cbar_ax = pyplot.gcf().add_axes(cbar_rect)
    pyplot.colorbar(p, cax=cbar_ax)


def postprocess(name, p, N):
    label = '%s_%s_%02d' % (name, p, N)

    # Read DOLFIN functions
    g = read_p1_function(label, "glob")
    l = read_p1_function(label, "loc")
    f = read_p1_function(label, "fe")
    e = read_p1_function(label, "ee")

    # Create and save plots
    plot_errors(g, l, f, e, label)
    plot_effectivities(g, l, f, label)


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

    # Run all tests
    if len(argv) == 1:
        for test in default_tests:
            postprocess(*test)
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
    postprocess(argv[1], float(argv[2]), int(argv[3]))
    return


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
