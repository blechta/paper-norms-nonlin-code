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
#matplotlib.rc('font', family='serif', serif=['Latin Modern Roman'])
#matplotlib.rc('font', family='serif', serif=['computer modern roman'])
#matplotlib.rc('font', family='serif')
#matplotlib.rc('mathtext', fontset='cm')
matplotlib.rc('mathtext', fontset='dejavuserif')
matplotlib.use("agg")

from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d

import numpy as np

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
    _plot_subplots(2, 2, functions, [0.85, 0.55, 0.02, 0.40], 12, 8, zmin=0.0)
    frame = pyplot.Polygon((
        ( 92,  16),
        (340,  16),
        (340, 296),
        (672, 296),
        (672, 564),
        ( 92, 564),
    ), fill=False)
    fig.patches.append(frame)

    # Plot e
    functions = [
        (4, e, r"$\epsilon_\mathrm{en}^p$"),
    ]
    _plot_subplots(2, 2, functions, [0.85, 0.05, 0.02, 0.40], 12, 8)
    frame = pyplot.Rectangle((364, 16), 308, 264, fill=False)
    fig.patches.append(frame)

    # Save as PDF
    pyplot.savefig(os.path.join(path, label+"_glfe_w.pdf"))


def plot_effectivities(g, l, f, label):
    path = "./"

    # Prepare figure
    fig = pyplot.figure(figsize=(6, 12))

    # Treat CK differently (cut-off effectivity plots)
    if "CarstensenKlose" in label:
        functions = [
            (1, g, l, r"\frac{\epsilon_\mathrm{glob}^q}{\epsilon_\mathrm{loc }^q}", 1.0, 5.0,  5.0,  'max'),
            (2, f, l, r"\frac{\epsilon_\mathrm{flux}^q}{\epsilon_\mathrm{loc }^q}", 1.0, 5.0,  5.0,  'max'),
            (3, f, g, r"\frac{\epsilon_\mathrm{flux}^q}{\epsilon_\mathrm{glob}^q}", 1.0, None, None, 'neither'),
        ]
    else:
        functions = [
            (1, g, l, r"\frac{\epsilon_\mathrm{glob}^q}{\epsilon_\mathrm{loc }^q}", 1.0, None, None, 'neither'),
            (2, f, l, r"\frac{\epsilon_\mathrm{flux}^q}{\epsilon_\mathrm{loc }^q}", 1.0, None, None, 'neither'),
            (3, f, g, r"\frac{\epsilon_\mathrm{flux}^q}{\epsilon_\mathrm{glob}^q}", 1.0, None, None, 'neither'),
        ]

    for i, f1, f2, title, zmin, zmax, cmax, cext in functions:
        #assert(f1.function_space() == f2.function_space())
        assert(f1.function_space().ufl_element() == f2.function_space().ufl_element())
        assert(f1.function_space().mesh().hash() == f2.function_space().mesh().hash())
        eff = Function(f1.function_space())
        eff.vector()[:] = f1.vector().array() / f2.vector().array()

        print(r"\max_\Omega {} = {}".format(title, eff.vector().array().max()))
        print(r"\min_\Omega {} = {}".format(title, eff.vector().array().min()))

        # Trim few very large overshoots to prevent undershoots in plots which
        # are rendering artifacts. Note that this does not change qualitative
        # look of overshoots in the figures, just prevents undershoots.
        eff_arr = eff.vector().array()
        print('Trimming values {} in {} to 80.0'.format(eff_arr[eff_arr > 80.0], title))
        eff_arr[eff_arr > 80.0] = 80.0
        eff.vector()[:] = eff_arr

        _plot_subplots(3, 1, [(i, eff, r"${}$".format(title))],
                       [0.75, 0.32*(3-i)+0.05, 0.02, 0.24], 21, 10,
                       zmin=zmin, zmax=zmax, cmax=cmax, cbar_extend=cext)

    # Save as PDF
    pyplot.savefig(os.path.join(path, label+"_eff_w.pdf"))


def _plot_subplots(nrows, ncols, functions, cbar_rect, tfs, lfs,
                   cmin=None, cmax=None, zmin=None, zmax=None,
                   cbar_extend='neither'):
    # Extract common range
    if cmin is None:
        cmin = min(item[1].vector().min() for item in functions)
    if cmax is None:
        cmax = max(item[1].vector().max() for item in functions)
    if zmin is None:
        zmin = min(item[1].vector().min() for item in functions)
    if zmax is None:
        zmax = max(item[1].vector().max() for item in functions)

    # Create subplots
    for i, func, t in functions:
        sp = pyplot.subplot(nrows, ncols, i, projection="3d")
        p = plot(func, backend="matplotlib", mode="warp",
                 vmin=cmin, vmax=cmax)
        ax = pyplot.gca(projection="3d")
        ax.set_zlim(zmin, zmax)
        ax.set_title(t, fontdict={'fontsize': tfs}, x=0.2)
        ax.ticklabel_format(style='sci', axis='z', scilimits=(-3, 6))
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_zlabel("", labelpad=4)
        ax.zaxis.get_offset_text().set_fontsize(lfs)
        ax.tick_params(axis='both', which='major', labelsize=lfs, pad=4)

    # Create common colorbar
    pyplot.subplots_adjust(right=0.8,wspace=0.15,top=0.92,bottom=0.04,hspace=0.315)
    cbar_ax = pyplot.gcf().add_axes(cbar_rect)
    cbar_ax.tick_params(axis='both', which='major', labelsize=lfs)
    cbar_ax.yaxis.get_offset_text().set_fontsize(lfs)
    cbar = pyplot.colorbar(p, cax=cbar_ax, extend=cbar_extend)
    cbar.formatter.set_powerlimits((-3, 6))
    cbar.update_ticks()


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
