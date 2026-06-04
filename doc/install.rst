.. _install:

Installation & Setup
====================

.. contents:: Contents
    :depth: 3

``FastSpecFit`` can be set up in two ways:

- :ref:`At NERSC<nersc installation>` — for DESI collaborators, using the shared software stack.
- :ref:`On a laptop<laptop installation>` — for interactive work or running on smaller datasets.

.. _nersc installation:

NERSC Installation
------------------

To run ``FastSpecFit`` at `NERSC`_, first request an interactive CPU node::

  salloc -N 1 -C cpu -A desi -t 01:00:00 --qos interactive


Then, load the package on top of the standard DESI software stack::

  source /dvs_ro/cfs/cdirs/desi/software/desi_environment.sh main
  module load fastspecfit/main

To load specific tagged versions, replace ``main`` with the desired
release, for example::

  source /dvs_ro/cfs/cdirs/desi/software/desi_environment.sh 26.3
  module load fastspecfit/3.4.3

Development Checkout
~~~~~~~~~~~~~~~~~~~~

To work with a local git clone of ``FastSpecFit`` at NERSC, first load the
DESI software stack as above (without loading the ``fastspecfit`` module), then
run the following **once** after cloning or after each ``git pull``::

  pip install --no-deps -e /path/to/fastspecfit

This registers the CLI entry points (``fastspec``, ``fastphot``, etc.) in the
active conda environment and generates the ``_version.py`` file that records the
correct package version in output file headers. The ``--no-deps`` flag skips
reinstalling dependencies already provided by the DESI software stack.

Alternatively, in some cases it may be necessary to bypass the *pip
install* step and manually add ``FastSpecFit`` to your ``PYTHONPATH``
and ``PATH`` via::

  export PYTHONPATH=/path/to/fastspecfit/py:${PYTHONPATH}
  export PATH=/path/to/fastspecfit/bin:${PATH}


JupyterHub
~~~~~~~~~~

To use ``FastSpecFit`` within NERSC's `JupyterHub`_, first install the
appropriate Jupyter kernel (one time). Then, to access the ``main``
branch, do::

  mkdir -p ${HOME}/.local/share/jupyter/kernels/fastspecfit-main
  wget -O ${HOME}/.local/share/jupyter/kernels/fastspecfit-main/kernel.json \
    https://raw.githubusercontent.com/desihub/fastspecfit/main/etc/jupyter-kernel.json

Or, to use the latest tagged release, do::

  mkdir -p ${HOME}/.local/share/jupyter/kernels/fastspecfit-latest
  wget -O ${HOME}/.local/share/jupyter/kernels/fastspecfit-latest/kernel.json \
    https://raw.githubusercontent.com/desihub/fastspecfit/main/etc/jupyter-kernel-latest.json

Then select the *FastSpecFit Main* or *FastSpecFit Latest* kernel from the
JupyterHub launcher.

.. _laptop installation:

Laptop Installation
-------------------

Create a Virtual Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend `micromamba`_ for creating isolated Python environments, though
`Miniforge`_ and `Anaconda`_ work equally well.

.. code-block:: bash

   micromamba create -n fastspec python
   micromamba activate fastspec

.. dropdown:: Using Miniforge or Anaconda instead

   .. code-block:: bash

      conda create -n fastspec python
      conda activate fastspec

Install FastSpecFit
~~~~~~~~~~~~~~~~~~~

``FastSpecFit`` is available on `PyPI`_ and installs with a single command::

  pip install fastspecfit

.. dropdown:: Developer installation

   To modify the source code, clone the repository and install in editable mode::

     git clone https://github.com/desihub/fastspecfit.git
     pip install -e fastspecfit

   To also install the test and documentation dependencies::

     pip install -e "fastspecfit[test,doc]"

Data Dependencies
~~~~~~~~~~~~~~~~~

Running ``FastSpecFit`` on DESI data requires four additional data
products, each pointed to by an environment variable.

``DESI_SPECTRO_REDUX``
   Root directory of the DESI spectroscopic reduction outputs. The
   DESI products are organized according to the `DESI data model`_;
   the `DESI data organization`_ page describes the full on-disk
   layout. The DESI `Early Data Release (EDR)`_ and `Data Release 1
   (DR1)`_ are already publicly available, with Data Release 2 (DR2)
   expected to be public in early 2027.

``DUST_DIR``
   Location of the `Schlegel, Finkbeiner, & Davis (1998)`_ dust reddening
   maps, used to correct spectra and photometry for Milky Way dust extinction.
   The maps can be downloaded from the NERSC public portal (see commands
   below).

``FPHOTO_DIR``
   Location of the `DESI Legacy Imaging Surveys DR9`_ broadband photometry,
   which provides optical and near-infrared fluxes used for joint
   spectrophotometric SED fitting. See the `DR9 description`_ for an overview
   of the imaging data and data products.

``FTEMPLATES_DIR``
   Location of the stellar population synthesis template files used by
   ``FastSpecFit`` to model the stellar continuum. Templates are distributed
   via the DESI public data portal and can be downloaded with the command
   below.

Set these variables and download the required external files::

  export DESI_SPECTRO_REDUX=/path/to/spectro/redux
  export DUST_DIR=/path/to/dustmaps
  export FPHOTO_DIR=/path/to/dr9/data
  export FTEMPLATES_DIR=/path/to/templates/fastspecfit

  wget -r -np -nH --cut-dirs 5 -A fits -P $DUST_DIR \
    https://portal.nersc.gov/project/cosmo/data/dust/v0_1/maps
  wget -O $FTEMPLATES_DIR/ftemplates-chabrier-2.0.0.fits \
    https://data.desi.lbl.gov/public/external/templates/fastspecfit/2.0.0/ftemplates-chabrier-2.0.0.fits

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the documentation locally, install the optional documentation
dependencies and run Sphinx::

  pip install "fastspecfit[doc]"
  sphinx-build -W --keep-going -b html doc doc/_build/html

.. _`micromamba`: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
.. _`Miniforge`: https://github.com/conda-forge/miniforge
.. _`Anaconda`: https://www.anaconda.com
.. _`PyPI`: https://pypi.org/project/fastspecfit
.. _`Schlegel, Finkbeiner, & Davis (1998)`: https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract
.. _`DESI data model`: https://desidatamodel.readthedocs.io/en/latest/
.. _`DESI data organization`: https://data.desi.lbl.gov/doc/organization/
.. _`DESI Legacy Imaging Surveys DR9`: https://www.legacysurvey.org/dr9
.. _`DR9 description`: https://www.legacysurvey.org/dr9/description/
.. _`NERSC`: https://www.nersc.gov/
.. _`JupyterHub`: https://jupyter.nersc.gov/
.. _`Early Data Release (EDR)`: https://data.desi.lbl.gov/doc/releases/edr/
.. _`Data Release 1 (DR1)`: https://data.desi.lbl.gov/doc/releases/dr1
