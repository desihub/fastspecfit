.. _install:

Installation & Setup
====================

``FastSpecFit`` is simple to install and has a relatively small number of
standard dependencies. Here, we describe three different ways of setting up
``FastSpecFit``:

  1. :ref:`at NERSC (for DESI collaborators)<nersc installation>`;
  2. :ref:`on a laptop<laptop installation>`; or
  3. :ref:`using the Docker container<docker installation>`.

.. _nersc installation:

1. NERSC Installation
---------------------

At `NERSC`_, ``FastSpecFit`` can be loaded trivially on top of the standard DESI
software stack. In a login or interactive `Perlmutter
<https://docs.nersc.gov/systems/perlmutter>`_ node, run the following commands
to load a stable version of ``FastSpecFit`` and its dependencies::

  source /global/cfs/cdirs/desi/software/desi_environment.sh 23.1
  module swap desispec/0.57.0
  module load fastspecfit/2.1.1

Alternatively, the following commands will load the development version of
``FastSpecFit``, which is updated nightly and not guaranteed to be stable::

  source /global/cfs/cdirs/desi/software/desi_environment.sh main
  module load fastspecfit/main

Finally, some users may want to access ``FastSpecFit`` through `NERSC`_'s
`JupyterHub`_ notebook server. To set up the kernel first do::

  mkdir -p ${HOME}/.local/share/jupyter/kernels/fastspecfit
  wget -O ${HOME}/.local/share/jupyter/kernels/fastspecfit/kernel.json \
    https://raw.githubusercontent.com/desihub/fastspecfit/main/etc/jupyter-kernel.json

Then, in `JupyterHub`_, simply select the *FastSpecFit* kernel and you are all
set!

.. _laptop installation:

2. Laptop Installation
----------------------

To install ``FastSpecFit`` and all its dependencies on a laptop we recommend a
dedicated `Miniforge`_ environment. For example, to install the latest stable
version transparently into an environment called, e.g., *fastspec* one would
do::

  conda create -y --name fastspec python=3.10
  conda activate fastspec
  pip install fastspecfit

Note that if you are planning to write documentation you will also need to
install the following dependencies::

  pip install sphinx sphinx-toolbox sphinx-rtd-theme sphinxcontrib-napoleon

Alternatively, some users may want ``FastSpecFit`` and its dependencies to be
installed in a more accessible location (e.g., */path/to/desi/code*), in which
case one would do::
  
  conda create -y --name fastspec python=3.10 numpy scipy numba astropy matplotlib seaborn
  conda activate fastspec
  pip install fitsio healpy speclite

  export DESI_CODE=/path/to/desi/code
  mkdir -p $DESI_CODE
  
  pushd $DESI_CODE 
  for package in desiutil desimodel desitarget desispec fastspecfit; do
    git clone https://github.com/desihub/$package.git
    export PATH=$DESI_CODE/$package/bin:$PATH
    export PYTHONPATH=$DESI_CODE/$package/py:$PYTHONPATH
  done
  popd

Finally, ``FastSpecFit`` has four more data dependencies, each specified with
their own environment variable:

  * ``DESI_ROOT``, which specifies the top-level location of the DESI data;
  * ``DUST_DIR``, which specifies the location of the `Schlegel, Finkbeiner, &
    Davis dust maps`_; 
  * ``DR9_DIR``, which specifies the location of the `DESI Legacy Imaging
    Surveys Data Release 9 (DR9)`_ data; and
  * ``FTEMPLATES_DIR``, which indicates the location of the stellar population
    synthesis models used by ``FastSpecFit``.

.. note::
   
  Currently, the DESI data are only available to DESI collaborators; however,
  the `Early Data Release (EDR)`_ is expected to be publicly available in Spring
  2023 and other data releases will be announced in the `DESI Data Release`_
  page, after which point the instructions here will be updated.

With the preceding caveat in mind, one can set up the remaining dependencies
with the following commands::

  export DESI_ROOT=/path/to/desi/data
  export DUST_DIR=/path/to/dustmaps
  export DR9_DIR=/path/to/dr9/data
  export FTEMPLATES_DIR=/path/to/templates/fastspecfit

  wget -r -np -nH --cut-dirs 5 -A fits -P $DUST_DIR \
    https://portal.nersc.gov/project/cosmo/data/dust/v0_1/maps
  wget -O $FTEMPLATES_DIR/ftemplates-chabrier-1.0.0.fits \
    https://data.desi.lbl.gov/public/external/templates/fastspecfit/1.0.0/ftemplates-chabrier-1.0.0.fits
  
.. _docker installation:

3. Using Docker
---------------

Finally, for production runs and for expert users, ``FastSpecFit`` is also
available as a Docker container which is served publicly in the
`DockerHub/desihub`_ repository.

For example, on a laptop one would retrieve (or update) and enter the *2.1.1*
version of the container with::
  
  docker pull desihub/fastspecfit:2.1.1
  docker run -it desihub/fastspecfit:2.1.1

Alternatively, at `NERSC`_ one would need to use `shifter`_::

  shifterimg pull docker:desihub/fastspecfit:2.1.1
  shifter --image docker:desihub/fastspecfit:2.1.1 bash

However, neither of the preceding commands define the required environment
variables, although we provide a simple setup script which does. For simple
interactive work at `NERSC`_ (e.g., in a login node) do::

  mkdir -p /path/to/fastspecfit/setup/script
  wget https://raw.githubusercontent.com/desihub/fastspecfit/main/bin/fastspecfit-setup.sh \
    -O /path/to/fastspecfit/setup/script/fastspecfit-setup.sh

  /path/to/fastspecfit/setup/script/fastspecfit-setup.sh shifter
  source /path/to/fastspecfit/setup/script/fastspecfit-setup.sh env

.. note::
  To run ``FastSpecFit`` on a large sample of objects (or for a full production
  or data release), please do not use a login node; instead, see the
  :ref:`running_fastspecfit` documentation for instructions and best practices.

.. _`Miniforge`: https://github.com/conda-forge/miniforge

.. _`Schlegel, Finkbeiner, & Davis dust maps`: https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract

.. _`DESI Legacy Imaging Surveys Data Release 9 (DR9)`: https://www.legacysurvey.org/dr9

.. _`NERSC`: https://www.nersc.gov/

.. _`JupyterHub`: https://jupyter.nersc.gov/ 

.. _`DockerHub/desihub`: https://hub.docker.com/u/desihub

.. _`shifter`: https://docs.nersc.gov/development/shifter/

.. _`Early Data Release (EDR)`: https://data.desi.lbl.gov/doc/releases/edr/

.. _`Data Release 1 (DR1)`: https://data.desi.lbl.gov/doc/releases/dr1

.. _`DESI Data Release`: https://data.desi.lbl.gov
