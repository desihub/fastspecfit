.. _algorithms:

Algorithms
==========

.. note::  

  This page is still under construction!

.. contents:: Contents
    :depth: 3

Overview
--------

``FastSpecFit`` leverages our understanding of galaxies and quasars to build a
physical model of their observed-frame optical spectra and UV/optical/IR
spectral energy distributions (SEDs). 

Simple Stellar Population Templates
-----------------------------------

Document the templates, physical properties, etc.

Modeling the Stellar Continuum
------------------------------

* Initial line-sigma estimation.
* Finding significant lines; masking.  
* Velocity dispersion (with criteria).
* Aperture correction.  
* Smooth continuum correction.
* No-photometry and photometry-only modes.

K-Corrections and Rest-Frame Photometry
---------------------------------------

* K-corrections.
* Default cosmology.
  
If the inverse variance on a given absolutely magnitude is zero it means that
the absolute magnitude was derived from *synthesized* photometry based on the
best-fitting model (i.e., use with care).
  
Modeling the Emission Lines
---------------------------

* Line-fitting (broad, narrow).
* Continuum, EWs, upper limits, etc.

.. _`planned improvements`:

Planned Improvements
--------------------
  
* Additional QSO templates.
