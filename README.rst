=======
desigal
=======

.. image:: https://img.shields.io/travis/desihub/desigal.svg
    :target: https://travis-ci.org/desihub/desigal
    :alt: Travis Build Status
.. image:: https://coveralls.io/repos/desihub/desigal/badge.svg?service=github
    :target: https://coveralls.io/github/desihub/desigal
    :alt: Test Coverage Status
.. image:: https://readthedocs.org/projects/desigal/badge/?version=latest
    :target: http://desigal.readthedocs.org/en/latest/
    :alt: Documentation Status

Introduction
============

This repository contains code and documentation for the DESI_ Galaxy and Quasar
Physics (GQP) Working Group.  For full documentation please visit `desigal on
Read the Docs`_.

.. _DESI: https://desi.lbl.gov
.. _`desigal on Read the Docs`: http://desigal.readthedocs.org/en/latest/

Installation
============

This product is installable using pip_ or it can cloned directly from `Github`_.
For example, to install the "0.1.0" tag (or release) do

.. code-block:: bash

  pip install git+https://github.com/desihub/desigal.git@0.1.0
  
Alternatively, you can clone the master branch from `Github`_ 
  
.. code-block:: bash

  git clone git@github.com:desihub/desigal.git

and then either install the package to an installation directory of your
choice

.. code-block:: bash

  python setup.py install --prefix=$INSTALL_DIR  

or explicitly add the code directory to your ``$PYTHONPATH`` environment, e.g. 

.. code-block:: bash

  export PYTHONPATH=/path/to/desigal/py:$PYTHONPATH
  export PATH=/path/to/desigal/bin:${PATH}

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

`desigal`_ is free software licensed under a 3-clause BSD-style license. For
details see the ``LICENSE.rst`` file.

.. _`desigal`: https://github.com/desihub/desigal
