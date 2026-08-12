.. _iron vac:
.. _iron fastspec vac:

Iron fastspec VAC (DR1)
=======================

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Iron`` fastspec value-added catalog, which
contains full spectrophotometric fitting results (continuum + emission
lines) and was publicly released in March 2025 as part of `DESI Data
Release 1 (DESI/DR1)`_. For the companion photometry-only catalog, see
the :ref:`Iron fastphot VAC<iron fastphot vac>`.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - iron>`.

Data Content & Access
---------------------

Data from the ``Iron`` VAC can be accessed at any of the following links:

============================ ==================================================================
Data url                     https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v4.0
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v4.0``
============================ ==================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts>`

Summary Statistics
------------------

The following tables summarize the size and number of targets in each merged
catalog, organized by survey and program.

.. rst-class:: columns

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastspec-iron-cmx-other.fits 12 MB     2,762
============================ ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastspec-iron-sv1-backup.fits 15 MB     3,331
fastspec-iron-sv1-bright.fits 542 MB    126,650
fastspec-iron-sv1-dark.fits   999 MB    233,202
fastspec-iron-sv1-other.fits  146 MB    34,112
Total (sv1)                   1.7 GB    397,295
============================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastspec-iron-sv2-backup.fits 647 KB    105
fastspec-iron-sv2-bright.fits 200 MB    46,486
fastspec-iron-sv2-dark.fits   226 MB    52,690
Total (sv2)                   426 MB    99,281
============================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastspec-iron-sv3-backup.fits 6.7 MB    1,524
fastspec-iron-sv3-bright.fits 1.1 GB    265,293
fastspec-iron-sv3-dark.fits   2.5 GB    592,191
Total (sv3)                   3.6 GB    859,008
============================= ========= =================

================================= ========= =================
File Name                         File Size Number of Targets
================================= ========= =================
fastspec-iron-special-backup.fits 2.5 MB    552
fastspec-iron-special-bright.fits 181 MB    43,261
fastspec-iron-special-dark.fits   63 MB     14,954
fastspec-iron-special-other.fits  178 MB    42,064
Total (special)                   425 MB    100,831
================================= ========= =================

========================================== ========= =================
File Name                                  File Size Number of Targets
========================================== ========= =================
fastspec-iron-main-backup.fits             64 MB     15,163
fastspec-iron-main-bright-nside1-hp00.fits 428 MB    101,838
fastspec-iron-main-bright-nside1-hp01.fits 4.2 GB    1,017,041
fastspec-iron-main-bright-nside1-hp02.fits 5.6 GB    1,358,627
fastspec-iron-main-bright-nside1-hp03.fits 1.7 GB    403,581
fastspec-iron-main-bright-nside1-hp04.fits 4.0 GB    981,600
fastspec-iron-main-bright-nside1-hp05.fits 1.2 GB    301,057
fastspec-iron-main-bright-nside1-hp06.fits 5.5 GB    1,347,464
fastspec-iron-main-bright-nside1-hp07.fits 2.8 GB    673,711
fastspec-iron-main-bright-nside1-hp08.fits 343 MB    81,734
fastspec-iron-main-bright-nside1-hp09.fits 281 MB    66,856
fastspec-iron-main-bright-nside1-hp10.fits 191 MB    45,570
fastspec-iron-main-bright-nside1-hp11.fits 281 MB    66,848
fastspec-iron-main-dark-nside1-hp00.fits   1.4 GB    352,447
fastspec-iron-main-dark-nside1-hp01.fits   4.6 GB    1,118,746
fastspec-iron-main-dark-nside1-hp02.fits   7.0 GB    1,699,122
fastspec-iron-main-dark-nside1-hp03.fits   903 MB    214,658
fastspec-iron-main-dark-nside1-hp04.fits   7.1 GB    1,739,317
fastspec-iron-main-dark-nside1-hp05.fits   2.4 GB    579,026
fastspec-iron-main-dark-nside1-hp06.fits   12 GB     2,851,879
fastspec-iron-main-dark-nside1-hp07.fits   4.2 GB    1,032,151
fastspec-iron-main-dark-nside1-hp08.fits   911 MB    216,757
fastspec-iron-main-dark-nside1-hp09.fits   658 MB    156,454
fastspec-iron-main-dark-nside1-hp10.fits   313 MB    74,463
fastspec-iron-main-dark-nside1-hp11.fits   171 MB    40,609
Total (main)                               68 GB     16,536,719
========================================== ========= =================

The following tables summarize the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

.. rst-class:: columns

============================ ================= ===============================
Catalog                      Number of Objects Number with Corrected Redshifts
============================ ================= ===============================
fastspec-iron-cmx-other.fits 2,762             33
============================ ================= ===============================

============================= ================= ===============================
Catalog                       Number of Objects Number with Corrected Redshifts
============================= ================= ===============================
fastspec-iron-sv1-backup.fits 3,331             40
fastspec-iron-sv1-bright.fits 126,650           27
fastspec-iron-sv1-dark.fits   233,202           2,220
fastspec-iron-sv1-other.fits  34,112            68
Total (sv1)                   397,295           2,355
============================= ================= ===============================

============================= ================= ===============================
Catalog                       Number of Objects Number with Corrected Redshifts
============================= ================= ===============================
fastspec-iron-sv2-backup.fits 105               0
fastspec-iron-sv2-bright.fits 46,486            1
fastspec-iron-sv2-dark.fits   52,690            485
Total (sv2)                   99,281            486
============================= ================= ===============================

============================= ================= ===============================
Catalog                       Number of Objects Number with Corrected Redshifts
============================= ================= ===============================
fastspec-iron-sv3-backup.fits 1,524             0
fastspec-iron-sv3-bright.fits 265,293           33
fastspec-iron-sv3-dark.fits   592,191           1,801
Total (sv3)                   859,008           1,834
============================= ================= ===============================

================================= ================= ===============================
Catalog                           Number of Objects Number with Corrected Redshifts
================================= ================= ===============================
fastspec-iron-special-backup.fits 552               0
fastspec-iron-special-bright.fits 43,261            2
fastspec-iron-special-dark.fits   14,954            149
fastspec-iron-special-other.fits  42,064            0
Total (special)                   100,831           151
================================= ================= ===============================

========================================== ================= ===============================
Catalog                                    Number of Objects Number with Corrected Redshifts
========================================== ================= ===============================
fastspec-iron-main-backup.fits             15,163            0
fastspec-iron-main-bright-nside1-hp00.fits 101,838           6
fastspec-iron-main-bright-nside1-hp01.fits 1,017,041         59
fastspec-iron-main-bright-nside1-hp02.fits 1,358,627         95
fastspec-iron-main-bright-nside1-hp03.fits 403,581           21
fastspec-iron-main-bright-nside1-hp04.fits 981,600           48
fastspec-iron-main-bright-nside1-hp05.fits 301,057           18
fastspec-iron-main-bright-nside1-hp06.fits 1,347,464         96
fastspec-iron-main-bright-nside1-hp07.fits 673,711           26
fastspec-iron-main-bright-nside1-hp08.fits 81,734            6
fastspec-iron-main-bright-nside1-hp09.fits 66,856            4
fastspec-iron-main-bright-nside1-hp10.fits 45,570            2
fastspec-iron-main-bright-nside1-hp11.fits 66,848            2
fastspec-iron-main-dark-nside1-hp00.fits   352,447           2,952
fastspec-iron-main-dark-nside1-hp01.fits   1,118,746         8,760
fastspec-iron-main-dark-nside1-hp02.fits   1,699,122         10,383
fastspec-iron-main-dark-nside1-hp03.fits   214,658           2,181
fastspec-iron-main-dark-nside1-hp04.fits   1,739,317         14,399
fastspec-iron-main-dark-nside1-hp05.fits   579,026           3,615
fastspec-iron-main-dark-nside1-hp06.fits   2,851,879         13,275
fastspec-iron-main-dark-nside1-hp07.fits   1,032,151         5,949
fastspec-iron-main-dark-nside1-hp08.fits   216,757           1,969
fastspec-iron-main-dark-nside1-hp09.fits   156,454           644
fastspec-iron-main-dark-nside1-hp10.fits   74,463            431
fastspec-iron-main-dark-nside1-hp11.fits   40,609            433
Total (main)                               16,536,719        65,374
========================================== ================= ===============================

Code & Template Versions
------------------------

The following tables document the code versions and environment variables used
to produce this VAC. For details regarding the revision history of
``FastSpecFit``, please see the `change log`_.

.. rst-class:: columns

================ ==========
Software Package Version(s)
================ ==========
python           3.13.12
numpy            2.3.5
scipy            1.16.3
astropy          7.2.0
yaml             6.0.3
matplotlib       3.10.8
fitsio           1.3.0
mpi4py           4.1.1
healpy           1.19.0
desiutil         3.6.1
desispec         0.71.0
desitarget       4.5.0
speclite         1.0.0
fastspecfit      3.6.1
================ ==========

.. rst-class:: columns

==================== ============================================================
Environment Variable Value
==================== ============================================================
DESI_SPECTRO_REDUX   /dvs_ro/cfs/cdirs/desi/spectro/redux
DUST_DIR             /dvs_ro/cfs/cdirs/cosmo/data/dust/v0_1
FPHOTO_DIR           /dvs_ro/cfs/cdirs/desi/external/legacysurvey/dr9
FTEMPLATES_DIR       /dvs_ro/cfs/cdirs/desi/public/external/templates/fastspecfit
FTEMPLATES_FILE      ftemplates-chabrier-2.2.0.fits
FPHOTO_FILE          legacysurvey-dr9.yaml
EMLINES_FILE         emlines.ecsv
CONSTRAINTS_FILE     emline-constraints.yaml
==================== ============================================================

.. _previous versions - iron:

Notes & Known Issues
--------------------

v4.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~

* Release date: 2027 (Exact date TBD)
* ``FastSpecFit`` version: ``3.6.1``
* Templates: ``ftemplates-chabrier-2.2.0.fits``  (see `README.txt`_).
* Known issues:

  * None at this time.

v3.0
~~~~

* Release date: June 2025
* ``FastSpecFit`` version: ``3.1.5``
* Templates: ``ftemplates-chabrier-2.0.0.fits``  (see `README.txt`_).
* Notes:

  * Near-total rewrite of the stellar continuum and emission-line fitting
    engines, resulting in improved accuracy and significantly faster execution
    [`PR/#177`_, `PR/#186`_].
  * Updated rest wavelengths and emission-line list; broad He I lines removed
    [`PR/#179`_].
  * Monte Carlo uncertainties added for all fitted model parameters (e.g.,
    line fluxes, equivalent widths, stellar masses) [`PR/#189`_].
  * Updated theoretical [N II] 6548, 6584 and [O III] 4959, 5007 doublet
    ratios [`PR/#195`_].
  * Updated algorithm for correcting QSO redshifts [`PR/#201`_].
  * Backwards-incompatible data model update; large merged catalogs now split
    into separate ``SPECPHOT`` and ``FASTSPEC`` extensions [`PR/#197`_,
    `PR/#211`_].
  * Bug fix: <1% bias in fluxes and equivalent widths of tied and free doublet
    ratios [`PR/#198`_].
  * Bug fix: significant errors in flux and equivalent-width inverse-variance
    estimates [`PR/#213`_].

* Known issues:

  * **Bug**: Poor line fits for a small fraction of strong emission lines due
    to a clipping error in the ``ivar2var`` function (fixed in ``3.2.0``; see
    `PR/#228`_).

v2.1
~~~~

* Release date: January 2024
* ``FastSpecFit`` version: ``2.5.0``
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

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/doc/releases/dr1
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
.. _`PR/#177`: https://github.com/desihub/fastspecfit/pull/177
.. _`PR/#179`: https://github.com/desihub/fastspecfit/pull/179
.. _`PR/#186`: https://github.com/desihub/fastspecfit/pull/186
.. _`PR/#189`: https://github.com/desihub/fastspecfit/pull/189
.. _`PR/#195`: https://github.com/desihub/fastspecfit/pull/195
.. _`PR/#197`: https://github.com/desihub/fastspecfit/pull/197
.. _`PR/#198`: https://github.com/desihub/fastspecfit/pull/198
.. _`PR/#201`: https://github.com/desihub/fastspecfit/pull/201
.. _`PR/#211`: https://github.com/desihub/fastspecfit/pull/211
.. _`PR/#213`: https://github.com/desihub/fastspecfit/pull/213
.. _`PR/#227`: https://github.com/desihub/fastspecfit/pull/227
.. _`PR/#228`: https://github.com/desihub/fastspecfit/pull/228
