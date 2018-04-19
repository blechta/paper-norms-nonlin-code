===========================================================================
Localization of the :math:`W^{-1,q}` norm for local a posteriori efficiency
===========================================================================

.. image:: https://zenodo.org/badge/21218/blechta/paper-norms-nonlin-code.svg
   :target: https://zenodo.org/badge/latestdoi/21218/blechta/paper-norms-nonlin-code

This repository contains a code to reproduce numerical simulations
of the paper

  Jan Blechta, Josef Málek, and Martin Vohralík.
  Localization of the :math:`W^{-1,q}` norm for local a posteriori efficiency.
  Submitted, 2016. URL https://hal.inria.fr/hal-01332481.


Usage
=====

On Linux machine with Docker installed, typing

.. code:: sh

  make

will

* download Docker image
  https://quay.io/blechta/dolfin-tape@sha256:2425f7973b3ee5d667a1176d7a51f0f334b311298d935b5bdb2296fee7189a3a
  containing necessary software stack (PETSc, FEniCS, dolfin-tape);
  the image can be rebuilt using the files in ``docker`` directory,

* run simulations in the Docker containers to produce log files
  and pdf figures,

* postprocess log files and write out ``tabular.tex``.

Type :code:`make -j2` to run 2 jobs in parallel, or :code:`make -j12` to run
all twelve test cases in parallel. Note that the biggest test case needs around
30GB of RAM and 30 core-hours. Completely parallel run ``-j12`` should be
possible with cca 60GB in cca 60 core-hours. Note that most of the cases can
be computed in few hours and less than 8GB of memory.

Run

.. code:: sh

  make INIT= RUNNER=$\(CMD\)

to avoid running in Docker container and rather use current environment.
You will need

|  PETSc         3.6.3
|  petsc4py      3.6.0
|  mpi4py        2.0.0
|  FIAT          c36b6d7a988a211b04048f64c7155b0c25ed5a52
|  UFL           0c5b1b90498aa4f9a25fb1999463d3c1c010199a
|  Instant       2f355dec4142c56eb4d464e5975a1e6ea3eac493
|  FFC           4dc648a466ad087448a41921ade005f114e41268
|  DOLFIN        fcf70d934d63168e5ed037678e22ac66fb2b3474
|  mshr          7447149c972977ff8ce2e89283d1ed0525fa2bc6
|  dolfin-tape   240e47f193d1f01b52f78effb29b56535274bd1b
|  matplotlib    1.3.1


Copyright
=========

Copyright 2015-2016, 2018 Jan Blechta


License
=======

GNU LGPLv3 except the contents of ``docker`` directory. See ``COPYING``,
``COPYING.LESSER`` and ``docker/LICENSE``.
