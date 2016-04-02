# This dockerfile builds particular version of FEniCS stack and
# dolfin-tape to reproduce numerical experiments of the paper
#
#   Jan Blechta, Josef Malek, and Martin Vohralik.
#   Localization of $W^{-1,q}$ norm for local a posteriori efficiency.
#   In preparation, 2016.
#
# Prebuilt image is stored at:
#
#    https://quay.io/blechta/dolfin-tape:paper0
#
# Authors:
# Jack S. Hale <jack.hale@uni.lu>
# Lizao Li <lixx1445@umn.edu>
# Garth N. Wells <gnw20@cam.ac.uk>
#
# Modified by Jan Blechta <blechta@karlin.mff.cuni.cz>

FROM phusion/baseimage:0.9.18
MAINTAINER Jan Blechta <blechta@karlin.mff.cuni.cz>

# Get Ubuntu updates
USER root
RUN apt-get update && \
    apt-get upgrade -y -o Dpkg::Options::="--force-confold" && \
    apt-get -y install locales && \
    echo "C.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set locale environment
ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    LANGUAGE=C.UTF-8

# Set up user so that we do not run as root
RUN useradd -m -s /bin/bash -G sudo,docker_env fenics && \
    echo "fenics:docker" | chpasswd && \
    echo "fenics ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# OpenBLAS threads should be 1 to ensure performance
RUN echo 1 > /etc/container_environment/OPENBLAS_NUM_THREADS && \
    echo 0 > /etc/container_environment/OPENBLAS_VERBOSE

# See https://github.com/phusion/baseimage-docker/issues/186
RUN touch /etc/service/syslog-forwarder/down
RUN echo "cat /home/fenics/WELCOME" >> /home/fenics/.bashrc

USER fenics
ENV HOME /home/fenics
WORKDIR /home/fenics
ENTRYPOINT ["sudo","/sbin/my_init","--quiet","--","sudo","-u","fenics","/bin/bash","-l","-c"]
CMD ["/bin/bash","-i"]

# Non-Python utilities and libraries
USER root
RUN apt-get -qq update && \
    apt-get -y --with-new-pkgs upgrade && \
    apt-get -y install \
        bison cmake flex g++ gfortran git ipython openmpi-bin pkg-config wget \
        libopenblas-dev libboost-filesystem-dev libboost-system-dev \
        libboost-program-options-dev libboost-thread-dev libboost-timer-dev \
        libboost-iostreams-dev libcln-dev libcppunit-dev libeigen3-dev \ 
        liblapack-dev libmpfr-dev libopenmpi-dev libhwloc-dev libpcre3-dev \
        libhdf5-openmpi-dev libcgal-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Our helper scripts
COPY fenics.conf $HOME/fenics.conf
COPY fenics.env.conf $HOME/fenics.env.conf
COPY bin $HOME/bin
RUN chown -R fenics:fenics $HOME

USER fenics

RUN echo "source /home/fenics/fenics.env.conf" >> /home/fenics/.bash_profile && \
    echo "source /home/fenics/fenics.conf" >> /home/fenics/.bashrc

# Version numbers
ENV PETSC_VERSION=3.6.3 \
    SLEPC_VERSION=3.6.2 \
    SWIG_VERSION=3.0.8 \
    PETSC4PY_VERSION=3.6.0 \
    SLEPC4PY_VERSION=3.6.0 \
    TRILINOS_VERSION=12.6.1

USER root
RUN echo "2" > /etc/container_environment/FENICS_PYTHON_MAJOR_VERSION && \
    echo "7" > /etc/container_environment/FENICS_PYTHON_MINOR_VERSION

WORKDIR /tmp

# Install Python2 based environment
RUN apt-get -qq update && \
    apt-get -y --with-new-pkgs upgrade && \
    apt-get -y install \
        python-dev python-numpy python-six python-ply python-pytest python-h5py \
        python-urllib3 python-setuptools python-pip python-matplotlib && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install PETSc from source
RUN wget --quiet -nc http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${PETSC_VERSION}.tar.gz && \
    tar -xf petsc-lite-${PETSC_VERSION}.tar.gz && \
    cd petsc-${PETSC_VERSION} && \
    ./configure --COPTFLAGS="-O2" \
                --CXXOPTFLAGS="-O2" \
                --FOPTFLAGS="-O2" \
                --with-blas-lib=/usr/lib/libopenblas.a --with-lapack-lib=/usr/lib/liblapack.a \
                --with-c-support \
                --with-debugging=0 \
                --with-shared-libraries \
                --download-suitesparse \
                --download-scalapack \
                --download-metis \
                --download-parmetis \
                --download-ptscotch \
                --download-hypre \
                --download-mumps \
                --download-blacs \
                --download-spai \
                --download-ml \
                --prefix=/usr/local && \
     make && \
     make install && \
     rm -rf /tmp/*

# Install SLEPc from source
RUN export PETSC_DIR=/usr/local && \
    wget -nc --quiet http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-${SLEPC_VERSION}.tar.gz -O slepc-${SLEPC_VERSION}.tar.gz && \
    tar -xf slepc-${SLEPC_VERSION}.tar.gz && \
    cd slepc-${SLEPC_VERSION} && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    rm -rf /tmp/*
ENV SLEPC_DIR=/usr/local \
    PETSC_DIR=/usr/local

RUN pip install --upgrade pip && \
    pip install --upgrade setuptools

# Install sympy, petsc4py, slepc4py and swig from source
RUN pip install sympy && \
    easy_install \
     https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-2.0.0.tar.gz && \
    easy_install \
     https://bitbucket.org/petsc/petsc4py/downloads/petsc4py-${PETSC4PY_VERSION}.tar.gz && \
    easy_install \
     https://bitbucket.org/slepc/slepc4py/downloads/slepc4py-${SLEPC4PY_VERSION}.tar.gz && \
    cd /tmp && \
    wget -nc --quiet http://downloads.sourceforge.net/swig/swig-${SWIG_VERSION}.tar.gz && \
    tar xf swig-${SWIG_VERSION}.tar.gz && \
    cd swig-${SWIG_VERSION} && \
    ./configure && \
    make && \
    make install && \
    rm -rf /tmp/*

USER fenics
WORKDIR $HOME

ENV FENICS_PYTHON_MAJOR_VERSION=2 \
    FENICS_PYTHON_MINOR_VERSION=7 \
    FENICS_BUILD_TYPE=Release

COPY bin $HOME/bin

# Build needed versions of FEniCS and dolfin-tape
RUN bash -l -c "source $HOME/fenics.conf && \
                pull_fenics_paper0 && \
                build_fenics && \
                build_fenics dolfin-tape"

# NFS-safe file-locking in Instant
RUN bash -l -c "pip install --prefix=/home/fenics/build/ flufl.lock"
