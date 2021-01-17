.. _install:

Installation
===============

The DESI spectroscopic pipeline requires and interfaces with other external software packages.  This document assumes that you are setting up a "development" software stack from scratch based on the latest versions of the DESI tools.

External Dependencies
------------------------

In order to set up a working DESI pipeline, there are several external software packages that must be installed on your system.  There are several ways of obtaining these dependencies:  using your OS package manager, using a third-party package manager (e.g. macports, etc), or using one of the conda-based Python distributions (e.g. Anaconda from Continuum Analytics).

The list of general software outside the DESI project which is needed for the spectroscopic pipeline is:

    * BLAS / LAPACK (OpenBLAS, MKL, Accelerate framework, etc)
    * CFITSIO
    * BOOST
    * requests
    * matplotlib
    * scipy
    * astropy
    * fitsio
    * pyyaml
    * speclite

If you wish to run the pipeline on a cluster, also ensure that you have installed mpi4py (which obviously requires a working MPI installation).

Installing Dependencies on a Linux Workstation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For development, debugging, and small tests it is convenient to install the pipeline tools on a laptop or workstation.  If you are using a Linux machine, you can get all the dependencies (except fitsio and speclite, which are pip-installable) from your distribution's package manager.  Alternatively, you can install just the non-python dependencies with your package manager and then use Anaconda for your python stack.  On OS X, I recommend using macports to get at least the non-python dependencies, and perhaps all of the python tools as well.

**Example:  Ubuntu GNU/Linux**

On an Ubuntu machine, you could install all the dependencies with::

    %> sudo apt-get install libboost-all-dev libcfitsio-dev \
        libopenblas-dev liblapack-dev python3-matplotlib \
        python3-scipy python3-astropy python3-requests \
        python3-yaml python3-mpi4py

    %> pip install --no-binary :all: fitsio speclite iniparser


Installing Dependencies on an OS X System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing scientific software on OS X is often more difficult than Linux, since Apple is primarily concerned with development of apps using Xcode.  The approach described here for installing desispec dependencies seems to get the job done with the fewest steps.

First, install homebrew::

    %> /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Now use homebrew to install CFITSIO, BOOST, and OpenMPI::

    %> brew install cfitsio
    %> brew install boost
    %> brew install openmpi

We are going to use homebrew to install our python stack.  Some people prefer Anaconda, but that distribution has several problems on OS X.  When installing python, we add the flag to indicate we want the python3 versions as well::

    %> brew install homebrew/python/mpi4py --with-python3
    %> brew install homebrew/python/scipy --with-python3
    %> brew install homebrew/python/matplotlib --with-python3

The rest of the dependencies we can install with pip::

    %> pip3 install requests pyyaml iniparser speclite astropy
    %> pip3 install --no-binary :all: fitsio


Dependencies at NERSC
~~~~~~~~~~~~~~~~~~~~~~~~~

At NERSC there is already a conda-based python stack and a version of the non-python dependencies installed.  You can add the necessary module search path by doing (you can safely add this to your ~/.bashrc.ext)::

    %> module use /global/common/${NERSC_HOST}/contrib/desi/modulefiles

and then whenever you want to load the software::

    %> module load desi-conda

