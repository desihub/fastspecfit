============================================
Welcome to the Documentation for FastSpecFit
============================================

.. image:: _static/fastspecfit-logo.png
   :scale: 30%

Overview
--------

``FastSpecFit`` is a stellar continuum and emission-line modeling code for `DESI
<https://desi.lbl.gov>`_ which is optimized for speed and simplicity; it can be
called independently to fit the broadband photometry (``fastphot``) or the
three-camera optical spectrophotometry from DESI (``fastspec``), although in
either case DESI redshifts are required.

Contents
--------

.. toctree::
   :maxdepth: 1

   install.rst
   running.rst
   algorithms.rst
   api.rst
   changes.rst
  

Value-Added Catalogs
--------------------

.. toctree::
   :maxdepth: 1

   fuji.rst
   iron.rst

Data Model
----------

.. toctree::
   :maxdepth: 1

   fastspec.rst
   fastphot.rst

Indices & Tables
----------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
