==========
Change Log
==========

2.2.0 (not released yet)
------------------------

* Allow the Redrock redshift to be overridden [`PR #115`_].
* Code to support fitting stacked spectra [`PR #116`_].
* Bug fix of reversed tied flux ratio of [OII]7320,7330 doublet [`PR #120`_].
* Do not constrain the SPS age by default [`PR #132`_].
* Bug fix of emission-line subtracted Dn(4000) measurement [`PR #135`_].
* Update IGM attenuation coefficients [`PR #136`_].

.. _`PR #115`: https://github.com/desihub/fastspecfit/pull/115
.. _`PR #116`: https://github.com/desihub/fastspecfit/pull/116
.. _`PR #120`: https://github.com/desihub/fastspecfit/pull/120
.. _`PR #132`: https://github.com/desihub/fastspecfit/pull/132
.. _`PR #135`: https://github.com/desihub/fastspecfit/pull/135
.. _`PR #136`: https://github.com/desihub/fastspecfit/pull/136

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

