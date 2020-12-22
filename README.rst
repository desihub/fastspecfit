===========
fastspecfit
===========

|Actions Status| |Coveralls Status| |Documentation Status|

.. |Actions Status| image:: https://github.com/desihub/fastspecfit/workflows/CI/badge.svg
    :target: https://github.com/desihub/fastspecfit/actions
    :alt: GitHub Actions CI Status

.. |Coveralls Status| image:: https://coveralls.io/repos/desihub/fastspecfit/badge.svg
    :target: https://coveralls.io/github/desihub/fastspecfit
    :alt: Test Coverage Status
.. image:: https://readthedocs.org/projects/fastspecfit/badge/?version=latest
    :target: http://fastspecfit.readthedocs.org/en/latest/
    :alt: Documentation Status

Introduction
============

This repository contains code and documentation to perform fast, simple spectral
synthesis and emission-line fitting of DESI spectra. For full documentation
please visit `fastspecfit on Read the Docs`_.

.. _DESI: https://desi.lbl.gov
.. _`fastspecfit on Read the Docs`: http://fastspecfit.readthedocs.org/en/latest/

Installation
============

This product is installable using pip_ or it can cloned directly from `Github`_.
For example, to install the "0.1.0" tag (or release) do

.. code-block:: bash

  pip install git+https://github.com/desihub/fastspecfit.git@0.1.0
  
Alternatively, you can clone the master branch from `Github`_ 
  
.. code-block:: bash

  git clone git@github.com:desihub/fastspecfit.git

and then either install the package to an installation directory of your
choice

.. code-block:: bash

  python setup.py install --prefix=$INSTALL_DIR  

or explicitly add the code directory to your ``$PYTHONPATH`` environment, e.g. 

.. code-block:: bash

  export PYTHONPATH=/path/to/fastspecfit/py:$PYTHONPATH
  export PATH=/path/to/fastspecfit/bin:${PATH}

.. _pip: http://pip.readthedocs.org
.. _Github: https://github.com

Product Contents
================

bin/
    Executable scripts.
doc/
    High-level documentation (.rst files).
etc/
    Small data and configuration files.
py/
    Python code.

License
=======

`fastspecfit`_ is free software licensed under a 3-clause BSD-style license. For
details see the ``LICENSE.rst`` file.

.. _`fastspecfit`: https://github.com/desihub/fastspecfit
