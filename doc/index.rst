============================================
Welcome to the Documentation for FastSpecFit
============================================

.. image:: _static/fastspecfit-logo.png
   :scale: 30%

Overview
--------

``FastSpecFit`` is a stellar continuum and emission-line modeling code which is
optimized for speed and simplicity; it uses physically motivated stellar
population synthesis and emission-line templates to jointly model UV/optical
spectrophotometry and ultraviolet through infrared broadband photometry
(``fastspec``). Alternatively, it can be used to model just the broadband
photometry (``fastphot``).

``FastSpecFit`` currently only supports `DESI <https://desi.lbl.gov>`_
spectrophotometric data, although contributions to develop the software for
other instruments and surveys are more than welcome.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   install.rst
   running.rst
   data.rst
   tutorials.rst
   algorithms.rst
   api.rst
   changes.rst

.. toctree::
   :maxdepth: 1
   :caption: Data Model

   fastspec.rst
   fastphot.rst

.. _vacs index:

.. toctree::
   :maxdepth: 1
   :caption: Value Added Catalogs (VACs)

   vacs.rst
   acknowledgments.rst

.. toctree::
   :maxdepth: 1
   :caption: DR2 (Loa)

   loa.rst
   loa-fastphot.rst

.. toctree::
   :maxdepth: 1
   :caption: DR1 (Iron)

   iron.rst
   iron-fastphot.rst

.. toctree::
   :maxdepth: 1
   :caption: DR1 Supplement (Guadalupe)

   guadalupe.rst

.. toctree::
   :maxdepth: 1
   :caption: EDR (Fuji)

   fuji.rst

Indices & Tables
----------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
