==========
Change Log
==========

3.2.0 (not released yet)
------------------------

* 

3.1.1 (2025-01-05)
------------------

* Miscellaneous bug fixes [`PR #205`_].
* Pure-MPI implementation; new Podman container; bug fixes [`PR #203`_].
* Updated algorithm for updating QSO redshifts [`PR #201`_].
* Progress toward pure-MPI production code [`PR #200`_].
* Fix <1% bias in fluxes and EWs of tied and free doublet ratios [`PR #198`_].
* Backwards incompatible update to the data model; expanded unit tests [`PR #197`_].

.. _`PR #197`: https://github.com/desihub/fastspecfit/pull/197
.. _`PR #198`: https://github.com/desihub/fastspecfit/pull/198
.. _`PR #200`: https://github.com/desihub/fastspecfit/pull/200
.. _`PR #201`: https://github.com/desihub/fastspecfit/pull/201
.. _`PR #203`: https://github.com/desihub/fastspecfit/pull/203
.. _`PR #205`: https://github.com/desihub/fastspecfit/pull/205

3.1.0 (2024-11-21)
------------------

* Update theoretical [NII] 6548,84 and [OIII] 4959,5007 doublet ratios [`PR #195`_].
* Additional cleanups and speedups to the new Monte Carlo code [`PR #194`_].
* Implement a simplified idiom for Monte Carlo iteration [`PR #193`_].
* Merge two multiprocessing loops into one [`PR #192`_].
* Monte Carlo to estimate model parameters [`PR #189`_].
* Fix issues with log level retention [`PR #187`_].

.. _`PR #187`: https://github.com/desihub/fastspecfit/pull/187
.. _`PR #189`: https://github.com/desihub/fastspecfit/pull/189
.. _`PR #192`: https://github.com/desihub/fastspecfit/pull/192
.. _`PR #193`: https://github.com/desihub/fastspecfit/pull/193
.. _`PR #194`: https://github.com/desihub/fastspecfit/pull/194
.. _`PR #195`: https://github.com/desihub/fastspecfit/pull/195

3.0.0 (2024-10-30)
------------------

* Rewrite of the smooth continuum and continuum-flux algorithms [`PR #186`_].
* Miscellaneous bug fixes [`PR #185`_].
* Use Numba-based Resolution implementation with faster multiply.
  Also tweak parameters of continuum least_squares optimizer for
  speed. [`PR #181`_]
* Cache Fourier transforms of the templates for speed [`PR #180`_].
* Update line-list (drop broad HeI lines) and rest wavelengths [`PR #179`_].
* Near-total rewrite of both the continuum and emission-line fitting engines and
  major updates to the organizational infrastructure of the code, all with the
  goal of maximizing speed and accuracy [`PR #177`_].

.. _`PR #177`: https://github.com/desihub/fastspecfit/pull/177
.. _`PR #179`: https://github.com/desihub/fastspecfit/pull/179
.. _`PR #180`: https://github.com/desihub/fastspecfit/pull/180
.. _`PR #181`: https://github.com/desihub/fastspecfit/pull/181
.. _`PR #185`: https://github.com/desihub/fastspecfit/pull/185
.. _`PR #186`: https://github.com/desihub/fastspecfit/pull/186


2.5.2 (2024-04-28)
------------------

* Add support for processing a custom sample with full MPI [`PR #168`_].
* New optional inputs to ``mpi-fastspecfit``: ``--no-smooth-continuum``,
  ``--fphotodir``, and ``--fphotofile``.
* Documentation bug in ``LOGLNU`` definitions (reported by R. Hada).
* Add ``--zmin`` optional argument to ``fastspecfit``.
* Other miscellaneous changes committed directly to `main`.

.. _`PR #168`: https://github.com/desihub/fastspecfit/pull/168

2.5.1 (2024-01-18)
------------------

* Bug fix: handle Lya falling in a masked camera (only for objects at z>5ish).

2.5.0 (2024-01-13)
------------------

* Address stellar mass bias; more constrained broad+narrow fitting; and max
  velocity dispersion increased to 475 km/s (bump template version to 1.3.0)
  [`PR #166`_].

.. _`PR #166`: https://github.com/desihub/fastspecfit/pull/166

2.4.3 (2023-12-03)
------------------

* Expand velocity dispersion grid to 425 km/s and bump the template version to
  1.2.0 [`PR #158`_].
* Handle variable pixel size in camera-overlap region bug [`PR #157`_].
* Fix minor bug that broke the stackfit functionality [`PR #155`_].

.. _`PR #155`: https://github.com/desihub/fastspecfit/pull/155
.. _`PR #157`: https://github.com/desihub/fastspecfit/pull/157
.. _`PR #158`: https://github.com/desihub/fastspecfit/pull/158

2.4.2 (2023-08-30)
------------------

* Fix incorrect syntax synthesizing photometry for highest-redshift targets bug
  [`PR #152`_].

.. _`PR #152`: https://github.com/desihub/fastspecfit/pull/152

2.4.1 (2023-08-23)
------------------

* Just two rounds of emission-line fitting [`PR #151`_].

.. _`PR #151`: https://github.com/desihub/fastspecfit/pull/151

2.4.0 (2023-08-19)
------------------

* Bug fixes and miscellaneous feature requests for next VACs, including modified
  SPS templates and a user-friendly refactor of the K-correction code [`PR #148`_].

.. _`PR #148`: https://github.com/desihub/fastspecfit/pull/148

2.3.0 (2023-08-07)
------------------

* Support non-DR9 photometry via a new photometric configuration file [`PR #133`_].
* Miscellaneous updates and bug fixes ahead of next version of VACs [`PR #145`_].

.. _`PR #133`: https://github.com/desihub/fastspecfit/pull/133
.. _`PR #145`: https://github.com/desihub/fastspecfit/pull/145

2.2.0 (2023-08-02)
------------------

* Allow the Redrock redshift to be overridden [`PR #115`_].
* Code to support fitting stacked spectra [`PR #116`_].
* Bug fix of reversed tied flux ratio of [OII]7320,7330 doublet [`PR #120`_].
* Do not constrain the SPS age by default [`PR #132`_].
* Bug fix of emission-line subtracted Dn(4000) measurement [`PR #135`_].
* Update IGM attenuation coefficients [`PR #136`_].
* Several significant changes [`PR #137`_]:

  * Record the observed-space emission-line amplitude in ``_AMP`` and move the
    model-space amplitude to ``_MODELAMP``.
  * Demand at least 12 pixels to measure the scatter in the pixels under the
    line (therefore ``_AMP_IVAR`` should be more reliable for narrow lines).
  * Major bug fix whereby the model emission-line spectra were not being
    convolved with the resolution matrix.
  * Redefine ``_CHI2`` for an emission line as the observed not reduced chi2.
  * Switch from (deprecated) ``pkg_resources`` to ``importlib``.
  * Updated documentation (data model) and several non-negligible speed-ups.

* Improved modeling of galaxies with broad+narrow line-emission [`PR #142`_]:

.. _`PR #115`: https://github.com/desihub/fastspecfit/pull/115
.. _`PR #116`: https://github.com/desihub/fastspecfit/pull/116
.. _`PR #120`: https://github.com/desihub/fastspecfit/pull/120
.. _`PR #132`: https://github.com/desihub/fastspecfit/pull/132
.. _`PR #135`: https://github.com/desihub/fastspecfit/pull/135
.. _`PR #136`: https://github.com/desihub/fastspecfit/pull/136
.. _`PR #137`: https://github.com/desihub/fastspecfit/pull/137
.. _`PR #142`: https://github.com/desihub/fastspecfit/pull/142

2.1.2 (2023-04-01)
------------------

* Web-app updates needed for Fuji/v2.0 database load [`PR #107`_].
* Get target cutouts using image coadds on-disk [`PR #108`_].
* Initial hooks to fit stacked spectra [`PR #113`_].
* Updated documentation for v2.0 Fujilupe and v1.0 Iron VACs [`PR #114`_].

.. _`PR #107`: https://github.com/desihub/fastspecfit/pull/107
.. _`PR #108`: https://github.com/desihub/fastspecfit/pull/108
.. _`PR #113`: https://github.com/desihub/fastspecfit/pull/113
.. _`PR #114`: https://github.com/desihub/fastspecfit/pull/114

2.1.1 (2023-02-22)
------------------

* Be robust to synthesizing photometry of the highest-redshift targets [`PR #101`_].
* Fix another corner-case crash to the highest-redshift targets [`PR #102`_].
* Do not crash if there are no lines to optimize [`PR #104`_].

.. _`PR #101`: https://github.com/desihub/fastspecfit/pull/101
.. _`PR #102`: https://github.com/desihub/fastspecfit/pull/102
.. _`PR #104`: https://github.com/desihub/fastspecfit/pull/104

2.1.0 (2023-02-17)
------------------

* Tests, bug fixes, and speed-ups of version 2.0.0 [`PR #99`_].

.. _`PR #99`: https://github.com/desihub/fastspecfit/pull/99

2.0.0 (2023-01-23)
------------------

* Support custom coadds, update laboratory line-wavelengths, and fix major EW
  bug [`PR #87`_].
* Refactor fitting engine to not use fnnls or astropy.modeling [`PR #92`_].
* Additional Fujilupe documentation [`PR #93`_].
* Webapp updates to support latest data model [`PR #94`_].
* Joint spectrophotometric fitting and much more [`PR #95`_].
* Additional fujilupe v2.0 updates [`PR #96`_].

.. _`PR #87`: https://github.com/desihub/fastspecfit/pull/87
.. _`PR #92`: https://github.com/desihub/fastspecfit/pull/92
.. _`PR #93`: https://github.com/desihub/fastspecfit/pull/93
.. _`PR #94`: https://github.com/desihub/fastspecfit/pull/94
.. _`PR #95`: https://github.com/desihub/fastspecfit/pull/95
.. _`PR #96`: https://github.com/desihub/fastspecfit/pull/96

1.0.1 (2022-08-11)
------------------

* Additional cleanup needed to finish Fujilupe processing [`PR #78`_].
* Initial Fuji and Guadalupe VAC documentation [`PR #77`_].

.. _`PR #77`: https://github.com/desihub/fastspecfit/pull/77
.. _`PR #78`: https://github.com/desihub/fastspecfit/pull/78

1.0.0 (2022-08-01)
------------------

* Update Docker container and tag all dependencies [`PR #76`_].
* Numerous backwards-incompatible improvements and changes to the code engine
  and data model in preparation for processing Fuji (EDR)+Guadalupe [`PR #69`_].
* Initial set-up of GitHub Actions and unit tests [`PR #61`_].
* Initial version of the web-application [`PR #60`_].
* First round of development work in preparation for Fuji [`PR #55`_].

.. _`PR #55`: https://github.com/desihub/fastspecfit/pull/55
.. _`PR #60`: https://github.com/desihub/fastspecfit/pull/60
.. _`PR #61`: https://github.com/desihub/fastspecfit/pull/61
.. _`PR #69`: https://github.com/desihub/fastspecfit/pull/69
.. _`PR #76`: https://github.com/desihub/fastspecfit/pull/76

0.3 (2022-01-19)
----------------

* Additional updates needed to complete Everest release [`PR #44`_].

.. _`PR #44`: https://github.com/desihub/fastspecfit/pull/44

0.2 (2021-09-04)
----------------

* Major update to support Everest data release [`PR #40`_].

.. _`PR #40`: https://github.com/desihub/fastspecfit/pull/40

0.1 (2021-07-29
----------------

* Fix spectroscopic Dn(4000) calculation bug [`PR #35`_].
* Add UBV rest-frame photometry [`PR #34`_].
* Additional template work [`PR #24`_].
* Initial code to build spectrophotometric templates [`PR #20`_].
* Additional updates needed to finish fitting all of Denali [`PR #18`_].
* First set of updates for Denali data release [`PR #16`_].

.. _`PR #16`: https://github.com/desihub/fastspecfit/pull/16
.. _`PR #18`: https://github.com/desihub/fastspecfit/pull/18
.. _`PR #20`: https://github.com/desihub/fastspecfit/pull/20
.. _`PR #24`: https://github.com/desihub/fastspecfit/pull/24
.. _`PR #34`: https://github.com/desihub/fastspecfit/pull/34
.. _`PR #35`: https://github.com/desihub/fastspecfit/pull/35

0.0.2 (2021-04-10)
------------------

* More flexible line-fitting and data model updates to handle the Cascades data
  release [`PR #15`_].

.. _`PR #15`: https://github.com/desihub/fastspecfit/pull/15

