.. _iron fastphot vac:

Iron fastphot VAC (DR1)
=======================

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Iron`` fastphot value-added catalog, which contains
photometry-only continuum fitting results and was publicly released as part of
`DESI Data Release 1 (DESI/DR1)`_. For the companion full spectrophotometric
catalog, see the :ref:`Iron fastspec VAC<iron fastspec vac>`.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - iron fastphot>`.

Data Content & Access
---------------------

Data from the ``Iron`` fastphot VAC can be accessed at any of the following links:

============================ ==================================================================
Data url                     TBD
`NERSC`_ (for collaborators) TBD
============================ ==================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts iron fastphot>`

.. _`qso redshifts iron fastphot`:

Updated QSO Redshifts
---------------------

The QSO redshift corrections applied to this VAC are identical to those used in
the :ref:`Iron fastspec VAC<iron fastspec vac>`; refer to that page for the
algorithm description and per-catalog correction counts.

Summary Statistics
------------------

The following tables summarize the size and number of targets in each merged
catalog, organized by survey and program.

.. rst-class:: columns

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastphot-iron-cmx-other.fits 2.34 MB   2,762
============================ ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-sv1-backup.fits 2.95 MB   3,331
fastphot-iron-sv1-bright.fits 106 MB    126,650
fastphot-iron-sv1-dark.fits   196 MB    233,202
fastphot-iron-sv1-other.fits  28.6 MB   34,112
Total (sv1)                   334 MB    397,295
============================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-sv2-backup.fits 132 KB    105
fastphot-iron-sv2-bright.fits 39.3 MB   46,486
fastphot-iron-sv2-dark.fits   44.4 MB   52,690
Total (sv2)                   83.8 MB   99,281
============================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-sv3-backup.fits 1.33 MB   1,524
fastphot-iron-sv3-bright.fits 229 MB    265,293
fastphot-iron-sv3-dark.fits   511 MB    592,191
Total (sv3)                   741 MB    859,008
============================= ========= =================

================================= ========= =================
File Name                         File Size Number of Targets
================================= ========= =================
fastphot-iron-special-backup.fits 458 KB    552
fastphot-iron-special-bright.fits 32.1 MB   43,261
fastphot-iron-special-dark.fits   11.1 MB   14,954
fastphot-iron-special-other.fits  32.9 MB   42,064
Total (special)                   76.5 MB   100,831
================================= ========= =================

========================================== =========== =================
File Name                                  File Size   Number of Targets
========================================== =========== =================
fastphot-iron-main-backup.fits             11.2 MB     15,163
fastphot-iron-main-bright-nside1-hp00.fits 76.8 MB     101,838
fastphot-iron-main-bright-nside1-hp01.fits 766 MB      1,017,041
fastphot-iron-main-bright-nside1-hp02.fits 1.02e+03 MB 1,358,627
fastphot-iron-main-bright-nside1-hp03.fits 304 MB      403,581
fastphot-iron-main-bright-nside1-hp04.fits 740 MB      981,600
fastphot-iron-main-bright-nside1-hp05.fits 227 MB      301,057
fastphot-iron-main-bright-nside1-hp06.fits 1.02e+03 MB 1,347,464
fastphot-iron-main-bright-nside1-hp07.fits 508 MB      673,711
fastphot-iron-main-bright-nside1-hp08.fits 61.6 MB     81,734
fastphot-iron-main-bright-nside1-hp09.fits 50.4 MB     66,856
fastphot-iron-main-bright-nside1-hp10.fits 34.4 MB     45,570
fastphot-iron-main-bright-nside1-hp11.fits 50.4 MB     66,848
fastphot-iron-main-dark-nside1-hp00.fits   269 MB      352,447
fastphot-iron-main-dark-nside1-hp01.fits   853 MB      1,118,746
fastphot-iron-main-dark-nside1-hp02.fits   1.26 GB     1,699,122
fastphot-iron-main-dark-nside1-hp03.fits   164 MB      214,658
fastphot-iron-main-dark-nside1-hp04.fits   1.29 GB     1,739,317
fastphot-iron-main-dark-nside1-hp05.fits   441 MB      579,026
fastphot-iron-main-dark-nside1-hp06.fits   2.12 GB     2,851,879
fastphot-iron-main-dark-nside1-hp07.fits   787 MB      1,032,151
fastphot-iron-main-dark-nside1-hp08.fits   165 MB      216,757
fastphot-iron-main-dark-nside1-hp09.fits   119 MB      156,454
fastphot-iron-main-dark-nside1-hp10.fits   56.8 MB     74,463
fastphot-iron-main-dark-nside1-hp11.fits   31 MB       40,609
Total (main)                               12.3 GB     16,536,719
========================================== =========== =================

Code & Template Versions
------------------------

The following tables document the code versions and environment variables used
to produce this VAC. For details regarding the revision history of
``FastSpecFit``, please see the `change log`_.

Note that the tagged dependencies can be retrieved from any FITS file with the
following bit of code::

  import fitsio
  from desiutil.depend import Dependencies
  codever = Dependencies(fitsio.read_header('/path/to/fastspecfit/file.fits', ext=0))
  for codename, version in codever.items():
      print(codename, version)

.. rst-class:: columns

================ ==========
Software Package Version(s)
================ ==========
python           3.10.14
numpy            1.22.4
scipy            1.8.1
astropy          6.0.1
yaml             6.0.1
matplotlib       3.8.4
fitsio           1.2.1
mpi4py           3.1.6
healpy           1.16.6
desiutil         3.4.3
desispec         0.68.1
desitarget       2.8.0
speclite         0.20
fastspecfit      3.1.5, 3.2.0
================ ==========

.. rst-class:: columns

==================== ============================================================
Environment Variable Value
==================== ============================================================
DESI_ROOT            /global/cfs/cdirs/desi
DUST_DIR             /dvs_ro/cfs/cdirs/desi/external/dust/v0_1
FPHOTO_DIR           /dvs_ro/cfs/cdirs/desi/external/legacysurvey/dr9
FTEMPLATES_DIR       /dvs_ro/cfs/cdirs/desi/public/external/templates/fastspecfit
FTEMPLATES_FILE      ftemplates-chabrier-2.0.0.fits (see `README.txt`_)
==================== ============================================================

.. _previous versions - iron fastphot:

Notes & Known Issues
--------------------

v1.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~~

* Release date: June 2026
* ``FastSpecFit`` version: ``3.1.5``
* Templates: ``ftemplates-chabrier-2.0.0.fits``  (see `README.txt`_).
* Notes:

  * Initial public release of the ``Iron`` fastphot VAC.

* Known issues:

  * None at this time.

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/doc/releases/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/public/external/templates/fastspecfit/README.txt
