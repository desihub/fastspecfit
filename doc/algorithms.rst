.. _algorithms:

Algorithms
==========

.. contents:: Contents
    :depth: 3

Overview
--------

``FastSpecFit`` leverages our understanding of galaxies and quasars to build a
physical model of their observed-frame optical spectra and UV/optical/IR
spectral energy distributions (SEDs). Specifically, ``FastSpecFit`` uses 



Simple Stellar Population Templates
-----------------------------------

Document the SSPs.


Modeling the Stellar Continuum
------------------------------

* Initial line-sigma estimation.
* Finding significant lines; masking.  
* Fast NNLS.
* Reddening. If CONTINUUM_AV_IVAR is zero it means that fitted for the
  (intrinsic) dust extinction failed.
* Velocity dispersion (with criteria).
* Smooth continuum correction.  


Modeling the Emission Lines
---------------------------

* Line-fitting (broad, narrow).
* Continuum, EWs, upper limits, etc.

K-Corrections and Rest-Frame Photometry
---------------------------------------

* K-corrections.
* Default cosmology.
  
If the inverse variance on a given absolutely magnitude is zero it means that
the absolute magnitude was derived from *synthesized* photometry based on the
best-fitting model (i.e., use with care).
  
.. _`planned improvements`:

Planned Improvements
--------------------
  
* Joint fitting of the broadband photometry and spectroscopy.
* Adding dust emission to the templates.
* Additional QSO templates.
* Broadband photometry is missing emission lines.
