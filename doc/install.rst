.. _install:

Installation & Setup
====================

``FastSpecFit`` can be set up in two ways:

- :ref:`At NERSC<nersc installation>` — for DESI collaborators using the shared software stack.
- :ref:`On a laptop<laptop installation>` — for interactive work or running on smaller datasets.

.. _nersc installation:

NERSC Installation
------------------

On a `NERSC`_ `Perlmutter
<https://docs.nersc.gov/systems/perlmutter/architecture>`_ login or
interactive node, load ``FastSpecFit`` on top of the standard DESI
software stack::

  source /global/cfs/cdirs/desi/software/desi_environment.sh main
  module load fastspecfit/main

To load a specific tagged version, replace ``main`` with the desired
release, for example::

  source /global/cfs/cdirs/desi/software/desi_environment.sh 26.3
  module load fastspecfit/3.2.1

JupyterHub
~~~~~~~~~~

To use ``FastSpecFit`` in `JupyterHub`_, install its Jupyter kernel
once. For the ``main`` branch, do::

  mkdir -p ${HOME}/.local/share/jupyter/kernels/fastspecfit-main
  wget -O ${HOME}/.local/share/jupyter/kernels/fastspecfit-main/kernel.json \
    https://raw.githubusercontent.com/desihub/fastspecfit/main/etc/jupyter-kernel.json

Or, to always track the latest tagged release::

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

Running ``FastSpecFit`` on DESI spectra requires four additional data products,
each specified by an environment variable:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Purpose
   * - ``DESI_SPECTRO_REDUX``
     - Root directory of DESI spectroscopic reduction outputs
   * - ``DUST_DIR``
     - `Schlegel, Finkbeiner, & Davis dust maps`_
   * - ``FPHOTO_DIR``
     - `DESI Legacy Imaging Surveys DR9`_ broadband photometry
   * - ``FTEMPLATES_DIR``
     - Stellar population synthesis template files

Set these variables and download the required external files::

  export DESI_SPECTRO_REDUX=/path/to/spectro/redux
  export DUST_DIR=/path/to/dustmaps
  export FPHOTO_DIR=/path/to/dr9/data
  export FTEMPLATES_DIR=/path/to/templates/fastspecfit

  wget -r -np -nH --cut-dirs 5 -A fits -P $DUST_DIR \
    https://portal.nersc.gov/project/cosmo/data/dust/v0_1/maps
  wget -O $FTEMPLATES_DIR/ftemplates-chabrier-2.0.0.fits \
    https://data.desi.lbl.gov/public/external/templates/fastspecfit/2.0.0/ftemplates-chabrier-2.0.0.fits

.. note::

   The DESI `Early Data Release (EDR)`_ and `Data Release 1 (DR1)`_ are
   publicly available.

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
.. _`Schlegel, Finkbeiner, & Davis dust maps`: https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract
.. _`DESI Legacy Imaging Surveys DR9`: https://www.legacysurvey.org/dr9
.. _`NERSC`: https://www.nersc.gov/
.. _`JupyterHub`: https://jupyter.nersc.gov/
.. _`Early Data Release (EDR)`: https://data.desi.lbl.gov/doc/releases/edr/
.. _`Data Release 1 (DR1)`: https://data.desi.lbl.gov/doc/releases/dr1
