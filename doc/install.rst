.. _install:

Installation & Setup
====================

``FastSpecFit`` is simple to install and has a relatively small number of
standard dependencies. Here, we describe two different ways of setting up
``FastSpecFit``:

  1. :ref:`at NERSC (for DESI collaborators)<nersc installation>`; or
  2. :ref:`on a laptop<laptop installation>`.

.. _nersc installation:

1. NERSC Installation
---------------------

At `NERSC`_, ``FastSpecFit`` can be loaded trivially on top of the standard DESI
software stack. In a login or interactive `Perlmutter
<https://docs.nersc.gov/systems/perlmutter>`_ node, simply run the following
commands::

  source /global/cfs/cdirs/desi/software/desi_environment.sh main
  module load fastspecfit/main

Or, to load a specific version or tag of the code simply replace ``main`` with
the desired version, for example::

  source /global/cfs/cdirs/desi/software/desi_environment.sh 24.8
  module load fastspecfit/2.5.2

JupyterHub
##########

Alternatively, some users may want to access ``FastSpecFit`` through `NERSC`_'s
`JupyterHub`_ notebook server. In order to use this service, however, we first
have to install the appropriate *kernel*. To set up the ``main`` kernel do (just
once!)::

  mkdir -p ${HOME}/.local/share/jupyter/kernels/fastspecfit-main
  wget -O ${HOME}/.local/share/jupyter/kernels/fastspecfit-main/kernel.json \
    https://raw.githubusercontent.com/desihub/fastspecfit/main/etc/jupyter-kernel.json

Or you can install a kernal with the latest version of the code by executing the
following commands (as before, just one time)::

  mkdir -p ${HOME}/.local/share/jupyter/kernels/fastspecfit-latest
  wget -O ${HOME}/.local/share/jupyter/kernels/fastspecfit-latest/kernel.json \
    https://raw.githubusercontent.com/desihub/fastspecfit/main/etc/jupyter-kernel-latest.json

Then, in `JupyterHub`_, select either the *FastSpecFit Main* or *FastSpecFit
Latest* kernel and then have fun coding!


.. _laptop installation:

2. Laptop Installation
----------------------

To install ``FastSpecFit`` and all its dependencies on a laptop we recommend a
dedicated `Miniforge`_ environment. (Unfortunately, ``FastSpecFit`` has not yet
been uploaded to `PyPi`_, although there is an `open ticket`_ to do so.)
Therefore, for the time being, the code and its dependencies must be installed
"by hand" in an accessible location (e.g., */path/to/desi/code*) with the
following commands::

  conda create -y --name fastspec python numpy scipy numba astropy matplotlib seaborn
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

Note that if you are planning to write documentation you will also need to
install the following dependencies::

  pip install sphinx sphinx-toolbox sphinx-rtd-theme sphinxcontrib-napoleon

Finally, ``FastSpecFit`` has four more data dependencies, each specified with
their own environment variable:

  * ``DESI_ROOT``, which specifies the top-level location of the DESI data;
  * ``DUST_DIR``, which specifies the location of the `Schlegel, Finkbeiner, &
    Davis dust maps`_;
  * ``FPHOTO_DIR``, which specifies the location of the `DESI Legacy Imaging
    Surveys Data Release 9 (DR9)`_ data; and
  * ``FTEMPLATES_DIR``, which indicates the location of the stellar population
    synthesis models used by ``FastSpecFit``.

These environment variables can be set with the following commands::

  export DESI_ROOT=/path/to/desi/data
  export DUST_DIR=/path/to/dustmaps
  export FPHOTO_DIR=/path/to/dr9/data
  export FTEMPLATES_DIR=/path/to/templates/fastspecfit

  wget -r -np -nH --cut-dirs 5 -A fits -P $DUST_DIR \
    https://portal.nersc.gov/project/cosmo/data/dust/v0_1/maps
  wget -O $FTEMPLATES_DIR/ftemplates-chabrier-2.0.0.fits \
    https://data.desi.lbl.gov/public/external/templates/fastspecfit/2.0.0/ftemplates-chabrier-2.0.0.fits

.. note::

  The DESI `Early Data Release (EDR)`_ became publicly available in
  June 2023!

.. _`Miniforge`: https://github.com/conda-forge/miniforge

.. _`PyPi`: https://packaging.python.org/en/latest

.. _`open ticket`: https://github.com/desihub/fastspecfit/issues/83

.. _`Schlegel, Finkbeiner, & Davis dust maps`: https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract

.. _`DESI Legacy Imaging Surveys Data Release 9 (DR9)`: https://www.legacysurvey.org/dr9

.. _`NERSC`: https://www.nersc.gov/

.. _`JupyterHub`: https://jupyter.nersc.gov/

.. _`DockerHub/desihub`: https://hub.docker.com/u/desihub

.. _`shifter`: https://docs.nersc.gov/development/shifter/

.. _`Early Data Release (EDR)`: https://data.desi.lbl.gov/doc/releases/edr/

.. _`Data Release 1 (DR1)`: https://data.desi.lbl.gov/doc/releases/dr1

.. _`DESI Data Release`: https://data.desi.lbl.gov
