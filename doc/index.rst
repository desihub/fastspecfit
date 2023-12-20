============================================
Welcome to the Documentation for FastSpecFit
============================================

.. image:: _static/fastspecfit-logo.png
   :scale: 30%

Overview
--------

``FastSpecFit`` is a stellar continuum and emission-line modeling code for `DESI
<https://desi.lbl.gov>`_ which is optimized for speed and simplicity; it uses
physically motivated stellar population synthesis and emission-line templates to
jointly model the three-camera optical spectrophotometry from DESI and the
ultraviolet through infrared broadband photometry (``fastspec``). Alternatively,
it can be used to model just the broadband photometry (``fastphot``), although
in either case DESI redshifts are required.

Contents
--------

.. toctree::
   :maxdepth: 1

   install.rst
   running.rst
   algorithms.rst
   api.rst
   changes.rst

Value Added Catalogs (VACs)
---------------------------

.. toctree::
   :maxdepth: 1

   vacs.rst
   fuji.rst
   iron.rst
   guadalupe.rst
   acknowledgments.rst

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
