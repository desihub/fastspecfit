.. _iron vac:

Iron VAC (DR1)
==============

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Iron`` value-added catalog, which will be publicly
released as part of the `DESI Data Release 1 (DESI/DR1)`_ sometime in 2024
(exact date TBD).

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - iron>`.

Data Content & Access
---------------------

Data from the ``Iron`` VAC can be accessed at any of the following links:

============================ ==================================================================
Data url                     https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v2.1
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1``
============================ ==================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts>`

Summary Statistics
------------------
  
The next two tables summarize the size and number of targets in each merged
catalog. The first table gives the sample in each survey and program, while the
second table combines all the individual programs into separate ``main``,
``sv``, and ``special`` catalogs.

.. rst-class:: columns

================================= ========= =================
File Name                         File Size Number of Targets
================================= ========= =================
fastspec-iron-cmx-other.fits      11.2 MB   2,762
fastspec-iron-main-backup.fits    59.1 MB   15,163
fastspec-iron-main-bright.fits    24.6 GB   6,445,926
fastspec-iron-main-dark.fits      38.5 GB   10,075,557
fastspec-iron-special-backup.fits 2.31 MB   552
fastspec-iron-special-bright.fits 169 MB    43,261
fastspec-iron-special-dark.fits   58.3 MB   14,954
fastspec-iron-special-other.fits  166 MB    42,064
fastspec-iron-sv1-backup.fits     13.6 MB   3,331
fastspec-iron-sv1-bright.fits     506 MB    126,650
fastspec-iron-sv1-dark.fits       931 MB    233,201
fastspec-iron-sv1-other.fits      136 MB    34,112
fastspec-iron-sv2-backup.fits     593 KB    105
fastspec-iron-sv2-bright.fits     186 MB    46,486
fastspec-iron-sv2-dark.fits       211 MB    52,690
fastspec-iron-sv3-backup.fits     6.25 MB   1,524
fastspec-iron-sv3-bright.fits     1.04 GB   265,292
fastspec-iron-sv3-dark.fits       2.32 GB   592,190
fastspec-iron.fits                70.9 GB   17,995,820
================================= ========= =================

The following table summarizes the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

.. rst-class:: columns

================================= ================= ===============================
Catalog                           Number of Objects Number with Corrected Redshifts
================================= ================= ===============================
fastspec-iron-cmx-other.fits      2,762             34
fastspec-iron-main-backup.fits    15,163            0
fastspec-iron-main-bright.fits    6,445,926         348
fastspec-iron-main-dark.fits      10,075,557        59,907
fastspec-iron-special-backup.fits 552               0
fastspec-iron-special-bright.fits 43,261            2
fastspec-iron-special-dark.fits   14,954            137
fastspec-iron-special-other.fits  42,064            0
fastspec-iron-sv1-backup.fits     3,331             40
fastspec-iron-sv1-bright.fits     126,650           27
fastspec-iron-sv1-dark.fits       233,201           2,199
fastspec-iron-sv1-other.fits      34,112            67
fastspec-iron-sv2-backup.fits     105               0
fastspec-iron-sv2-bright.fits     46,486            1
fastspec-iron-sv2-dark.fits       52,690            485
fastspec-iron-sv3-backup.fits     1,524             0
fastspec-iron-sv3-bright.fits     265,292           28
fastspec-iron-sv3-dark.fits       592,190           1,668
fastspec-iron.fits                17,995,820        64,943
================================= ================= ===============================

Code & Template Versions
------------------------

The following tables document the code versions and environment variables used
to produce this VAC. For details regarding the revision history of
``FastSpecFit``, please see the `change log`_.

Note that the tagged dependencies can be retrieve from any FITS file with the
following bit of code::

  import fitsio
  from desiutil.depend import Dependencies
  codever = Dependencies(fitsio.read_header('/path/to/fastspecfit/file.fits, ext=0))
  for codename, version in codever.items():
      print(codename, version)

.. rst-class:: columns

================ ==========
Software Package Version(s)
================ ==========
python           3.10.8
numpy            1.22.4
scipy            1.8.1
astropy          5.2.1
yaml             6.0
matplotlib       3.6.2
fitsio           1.1.8
desiutil         3.4.2
desispec         0.60.2
desitarget       2.7.0
desimodel        0.19.0
speclite         0.17
fastspecfit      2.5.0, 2.5.1
================ ==========

.. rst-class:: columns

==================== =====
Environment Variable Value
==================== =====
DESI_ROOT            /dvs_ro/cfs/cdirs/desi
DUST_DIR             /dvs_ro/cfs/cdirs/cosmo/data/dust/v0_1
FPHOTO_DIR           /dvs_ro/cfs/cdirs/desi/external/legacysurvey/dr9
FTEMPLATES_DIR       /dvs_ro/cfs/cdirs/desi/science/gqp/templates/fastspecfit
FTEMPLATES_FILE      ftemplates-chabrier-1.3.0.fits (see `README.txt`_)
FPHOTO_FILE          /global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/code/fastspecfit/2.5.1/lib/python3.10/site-packages/fastspecfit/data/legacysurvey-dr9.yaml
EMLINES_FILE         /global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/code/fastspecfit/2.5.1/lib/python3.10/site-packages/fastspecfit/data/emlines.ecsv
==================== =====

.. _previous versions - iron:

Notes & Known Issues
--------------------

v2.1 (latest release)
~~~~~~~~~~~~~~~~~~~~~

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
