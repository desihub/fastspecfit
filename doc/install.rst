.. _install:

Installation & Setup
====================

``FastSpecFit`` is simple to install and has a relatively small number of
standard dependencies. Here, we describe three different ways of setting up
``FastSpecFit``:

  1. :ref:`on a laptop<laptop installation>`; 
  2. :ref:`at NERSC (for DESI collaborators)<nersc installation>`; or
  3. :ref:`mounting the Docker container<docker installation>`.

.. _laptop installation:

1. Laptop Installation
----------------------

To install ``FastSpecFit`` and all its dependencies on a laptop we recommend a
dedicated `conda`_ environment. For example, to install everything transparently
into an environment called *fastspec* one would do::

  conda create -y --name fastspec python numpy scipy astropy matplotlib seaborn
  conda activate fastspec
  pip install fitsio speclite
  
  for package in desiutil desitarget desispec redrock fastspecfit; do
    python -m pip install git+https://github.com/desihub/$package.git@main#egg=$package
  done

Alternatively, some users may want the DESI software to be installed in a more
accessible location (e.g., */path/to/desi/code*), in which case one would do::
  
  conda create -y --name fastspec python numpy scipy astropy matplotlib seaborn
  conda activate fastspec
  pip install fitsio speclite

  export DESI_CODE=/path/to/desi/code
  mkdir -p $DESI_CODE
  
  pushd $DESI_CODE 
  for package in desiutil desitarget desispec redrock fastspecfit; do
    git clone https://github.com/desihub/$package.git
    export PATH=$DESI_CODE/$package/bin:$PATH
    export PYTHONPATH=$DESI_CODE/$package/py:$PYTHONPATH
  done
  popd

Finally, ``FastSpecFit`` has three more data dependencies, each specified with
their own environment variable:

  * ``DESI_ROOT``, which specifies the top-level location of the DESI data;
  * ``DUST_DIR``, which specifies the location of the `Schlegel, Finkbeiner, &
    Davis dust maps`_; and
  * ``FASTSPECFIT_TEMPLATES``, which specifies the location of the simple
    stellar population (SSP) templates used to model the stellar continuum.

.. note::
  Currently, the DESI data and spectral templates are only available to DESI
  collaborators; however, once the data are publicly released these instructions
  will be updated.

With the preceding caveat in mind, one can set up the remaining dependencies
with::
  
  export DESI_ROOT=/path/to/desi/data
  export DUST_DIR=/path/to/dust/maps
  export FASTSPECFIT_TEMPLATES=/path/to/fastspecfit/templates

**Finish writing these instructions.**

.. _nersc installation:

2. NERSC Installation
---------------------

At `NERSC`_ ``FastSpecFit`` is part of the standard DESI software stack and can
be loaded trivially. In a login or interactive node simply run the following
commands and you are read to go!::

  source /global/cfs/cdirs/desi/software/desi_environment.sh main
  module load fastspecfit/main
  
  export DESI_ROOT=/global/cfs/cdirs/desi
  export FASTSPECFIT_TEMPLATES=$DESI_ROOT/science/gqp/templates/SSP-CKC14z
  export DUST_DIR=/global/cfs/cdirs/cosmo/data/dust/v0_1

Alternatively, some users may want to access ``FastSpecFit`` through NERSC's
`JupyterHub notebook server`_. To set up the kernel first do::

  wget 




.. _docker installation:

3. Docker Installation
----------------------

Docker.

and all its dependencies are installed inside a Docker container, making it easy
to run either on a personal laptop (if you have the necessary data) or at NERSC
using *shifter*.

.. _`conda`: https://anaconda.org/

.. _`Schlegel, Finkbeiner, & Davis dust maps`: https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract

.. _`NERSC`: https://www.nersc.gov/

.. _`JupyterHub notebook server`: https://jupyter.nersc.gov/ 
