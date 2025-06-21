.. _iron vac:

Iron VAC (DR1)
==============

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Iron`` value-added catalog, which was
publicly released as part of the `DESI Data Release 1 (DESI/DR1)`_
sometime in 2025 March.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - iron>`.

Data Content & Access
---------------------

Data from the ``Iron`` VAC can be accessed at any of the following links:

============================ ==================================================================
Data url                     https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v3.0``
============================ ==================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts>`

Summary Statistics
------------------
  
The following four tables summarize the size and number of objects in
each merged catalog. The first table lists the objects in the
``main/backup`` program and in the ``cmx``, ``sv1``, ``sv2``, ``sv3``,
and ``special`` surveys; and the second and third tables list the
number of objects in the ``main/bright`` and ``main/dark`` program,
respectively.

.. note::

   The ``main/bright`` and ``main/dark`` samples have been divided
   into ``nside=`` healpixels in order to keep each individual catalog
   to a reasonable size.

.. rst-class:: columns

================================= ========= =================
File Name                         File Size Number of Targets
================================= ========= =================
fastspec-iron-cmx-other.fits      12 MB     2,762
fastspec-iron-sv1-backup.fits     14.5 MB   3,331
fastspec-iron-sv1-bright.fits     542 MB    126,650
fastspec-iron-sv1-dark.fits       997 MB    233,202
fastspec-iron-sv1-other.fits      146 MB    34,112
fastspec-iron-sv2-backup.fits     647 KB    105
fastspec-iron-sv2-bright.fits     199 MB    46,486
fastspec-iron-sv2-dark.fits       225 MB    52,690
fastspec-iron-sv3-backup.fits     6.7 MB    1,524
fastspec-iron-sv3-bright.fits     1.11 GB   265,293
fastspec-iron-sv3-dark.fits       2.48 GB   592,191
fastspec-iron-special-backup.fits 2.49 MB   552
fastspec-iron-special-bright.fits 181 MB    43,261
fastspec-iron-special-dark.fits   62.6 MB   14,954
fastspec-iron-special-other.fits  178 MB    42,064
fastspec-iron-main-backup.fits    63.4 MB   15,163
================================= ========= =================

========================================== ========= =================
File Name                                  File Size Number of Targets
========================================== ========= =================
fastspec-iron-main-bright-nside1-hp00.fits 427 MB    101,838
fastspec-iron-main-bright-nside1-hp01.fits 4.16 GB   1,017,041
fastspec-iron-main-bright-nside1-hp02.fits 5.56 GB   1,358,627
fastspec-iron-main-bright-nside1-hp03.fits 1.65 GB   403,581
fastspec-iron-main-bright-nside1-hp04.fits 4.02 GB   981,600
fastspec-iron-main-bright-nside1-hp05.fits 1.23 GB   301,057
fastspec-iron-main-bright-nside1-hp06.fits 5.51 GB   1,347,464
fastspec-iron-main-bright-nside1-hp07.fits 2.76 GB   673,711
fastspec-iron-main-bright-nside1-hp08.fits 343 MB    81,734
fastspec-iron-main-bright-nside1-hp09.fits 280 MB    66,856
fastspec-iron-main-bright-nside1-hp10.fits 191 MB    45,570
fastspec-iron-main-bright-nside1-hp11.fits 280 MB    66,848
========================================== ========= =================

======================================== ========= =================
File Name                                File Size Number of Targets
======================================== ========= =================
fastspec-iron-main-dark-nside1-hp00.fits 1.44 GB   352,447
fastspec-iron-main-dark-nside1-hp01.fits 4.58 GB   1,118,746
fastspec-iron-main-dark-nside1-hp02.fits 6.96 GB   1,699,122
fastspec-iron-main-dark-nside1-hp03.fits 901 MB    214,658
fastspec-iron-main-dark-nside1-hp04.fits 7.13 GB   1,739,317
fastspec-iron-main-dark-nside1-hp05.fits 2.37 GB   579,026
fastspec-iron-main-dark-nside1-hp06.fits 11.7 GB   2,851,879
fastspec-iron-main-dark-nside1-hp07.fits 4.23 GB   1,032,151
fastspec-iron-main-dark-nside1-hp08.fits 910 MB    216,757
fastspec-iron-main-dark-nside1-hp09.fits 657 MB    156,454
fastspec-iron-main-dark-nside1-hp10.fits 313 MB    74,463
fastspec-iron-main-dark-nside1-hp11.fits 171 MB    40,609
======================================== ========= =================

Code & Template Versions
------------------------

The following tables document the code versions and environment variables used
to produce this VAC. For details regarding the revision history of
``FastSpecFit``, please see the `change log`_.

Note that the tagged dependencies can be retrieve from any FITS file with the
following bit of code::

  import fitsio
  from desiutil.depend import Dependencies
  codever = Dependencies(fitsio.read_header('/path/to/fastspecfit/file.fits', ext=0))
  for codename, version in codever.items():
      print(codename, version)

.. rst-class:: columns

================ ============
Software Package Version(s)
================ ============
python           3.10.14
numpy            1.22.4
scipy            1.8.1
astropy          6.0.1
yaml             6.0.1
matplotlib       3.8.4
fitsio           1.2.1
desiutil         3.4.3
desispec         0.68.1
desitarget       2.8.0
speclite         0.20
fastspecfit      3.1.5, 3.2.0
================ ============

.. rst-class:: columns

==================== =====
Environment Variable Value
==================== =====
DESI_SPECTRO_REDUX   /global/cfs/cdirs/desi/spectro/redux
DUST_DIR             /global/cfs/cdirs/cosmo/data/dust/v0_1
FPHOTO_DIR           /global/cfs/cdirs/desi/external/legacysurvey/dr9
FTEMPLATES_DIR       /global/cfs/cdirs/desi/public/external/templates/fastspecfit
FTEMPLATES_FILE      ftemplates-chabrier-2.0.0.fits (see `README.txt`_)
FPHOTO_FILE          /global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/code/fastspecfit/3.1.5/lib/python3.10/site-packages/fastspecfit/data/legacysurvey-dr9.yaml
EMLINES_FILE         /global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/code/fastspecfit/3.1.5/lib/python3.10/site-packages/fastspecfit/data/emlines.ecsv
==================== =====

.. _previous versions - iron:

Notes & Known Issues
--------------------

v3.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~

* Release date: June 2025
* ``FastSpecFit`` version: ``3.1.5``, ``3.2.0``
* Templates: ``ftemplates-chabrier-2.0.0.fits``  (see `README.txt`_).
* Notes:

  * Several updates to the spectrophotometric templates aimed at addressing the
    stellar mass bias identified in `issue/#159`_ (see `PR/#166`_):
    
* Known issues:
  
  * None at this time.

v2.1
~~~~

* Release date: January 2024
* ``FastSpecFit`` version: ``2.5.0``, ``2.5.1``
* Templates: ``ftemplates-chabrier-1.3.0.fits``  (see `README.txt`_).
* Notes:

  * Several updates to the spectrophotometric templates aimed at addressing the
    stellar mass bias identified in `issue/#159`_ (see `PR/#166`_):
    
    * Templates are now just solar metallicity (previously 0.1, 1, and 1.6 times
      solar).
    * Five age bins now (vs 8 previously).
    * Expanded velocity dispersion grid (new measurable maximum value is now 475
      km/s).
  * Correction to how the light-weighted ages, dust attenuations, and SFRs were
    being computed.
  * When fitting the broad+narrow emission-line model, [OIII] 4959,5007 is now
    fitted separately and the narrow Balmer+helium+forbidden line-widths and
    velocity shifts are all tied together.
  * All known bugs fixed.
* Known issues:
  
  * None at this time.

v2.0
~~~~

* Release date: August 2023
* ``FastSpecFit`` versions: ``2.4.1``, ``2.4.2``
* Templates: ``ftemplates-chabrier-1.1.0.fits``  (see `README.txt`_).
* Notes:

  * Minor updates to spectrophotometric templates.
  * Just two rounds of emission-line fitting, not three (see `PR/#151`_).
  * Updated IGM attenuation coefficients (see `PR/#136`_).
  * Major algorithmic updates related to how emission-line amplitudes, fluxes,
    and inverse variances are computed, including a bug fix which the
    emission-line model spectra were not being convolved with the resolution
    matrix (see `PR/#137`_). 
* Known Issues:
  
  * **Warning**: Stellar masses are systematically higher (by 0.2-0.5 dex)
    compared to other methods, so they should be used with caution (see
    `issue/#159`_). Similarly, star-formation rates and other SPS model
    parameters have not been fully validated.
  * **Bug**: Fluxes (and EWs) of lines which lie in the camera-overlap region
    are overestimated by a factor of 2 due to a bug handling the different pixel
    scale (fixed in `PR/#157`_).

v1.0
~~~~

* Release date: February 2023
* ``FastSpecFit`` versions: ``2.1.0``, ``2.1.1``
* Templates: ``ftemplates-chabrier-1.0.0.fits``  (see `README.txt`_).
* Known Issues:
  
  * **Bug**: [OII] 7320,7330 doublet amplitude ratio incorrectly inverted (fixed
    in `PR/#120`_).
  * **Bug**: Artificial redshift dependence in derived stellar masses due to age
    prior (fixed in `PR/#132`_). 
  * **Bug**: Emission-line subtracted Dn(4000) values incorrectly computed
    (fixed in `PR/#135`_). 

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/public/external/templates/fastspecfit/README.txt
.. _`issue/#159`: https://github.com/desihub/fastspecfit/issues/159
.. _`PR/#120`: https://github.com/desihub/fastspecfit/pull/120
.. _`PR/#132`: https://github.com/desihub/fastspecfit/pull/132
.. _`PR/#135`: https://github.com/desihub/fastspecfit/pull/135
.. _`PR/#136`: https://github.com/desihub/fastspecfit/pull/136
.. _`PR/#137`: https://github.com/desihub/fastspecfit/pull/137
.. _`PR/#151`: https://github.com/desihub/fastspecfit/pull/151
.. _`PR/#157`: https://github.com/desihub/fastspecfit/pull/157
.. _`PR/#158`: https://github.com/desihub/fastspecfit/pull/158
.. _`PR/#166`: https://github.com/desihub/fastspecfit/pull/166
