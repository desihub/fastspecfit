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
   information regarding :ref:`Previous Versions<previous versions - iron>`.

Data Content & Access
---------------------

Data from the ``Iron`` VAC can be accessed at any of the following links:

============================ ============================================================
Data url                     https://data.desi.lbl.gov/desi/spectro/fastspecfit/iron/v2.0
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/spectro/fastspecfit/iron/v2.0``
============================ ============================================================

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
fastspec-iron-cmx-other.fits      11.5 MB   2,761
fastspec-iron-main-backup.fits    61.3 MB   15,163
fastspec-iron-main-bright.fits    25.5 GB   6,445,915
fastspec-iron-main-dark.fits      39.9 GB   10,074,406
fastspec-iron-special-backup.fits 2.38 MB   552
fastspec-iron-special-bright.fits 175 MB    43,261
fastspec-iron-special-dark.fits   60.5 MB   14,953
fastspec-iron-special-other.fits  172 MB    42,064
fastspec-iron-sv1-backup.fits     14 MB     3,331
fastspec-iron-sv1-bright.fits     524 MB    126,649
fastspec-iron-sv1-dark.fits       965 MB    233,158
fastspec-iron-sv1-other.fits      141 MB    34,112
fastspec-iron-sv2-backup.fits     593 KB    105
fastspec-iron-sv2-bright.fits     193 MB    46,486
fastspec-iron-sv2-dark.fits       218 MB    52,682
fastspec-iron-sv3-backup.fits     6.45 MB   1,524
fastspec-iron-sv3-bright.fits     1.08 GB   265,291
fastspec-iron-sv3-dark.fits       2.4 GB    592,165
fastspec-iron.fits                73.4 GB   17,994,578
================================= ========= =================

The following table summarizes the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

.. rst-class:: columns

=============================== ================= ===============================
Catalog                         Number of Targets Number with Corrected Redshifts
=============================== ================= ===============================
fastspec-iron-cmx-other.fits    2,771             63
fastspec-iron-special-dark.fits 35,647            389
fastspec-iron-sv1-backup.fits   3,683             119
fastspec-iron-sv1-bright.fits   126,677           402
fastspec-iron-sv1-dark.fits     235,881           4,656
fastspec-iron-sv1-other.fits    34,150            372
fastspec-iron-sv2-backup.fits   107               0
fastspec-iron-sv2-bright.fits   46,510            151
fastspec-iron-sv2-dark.fits     52,771            1,185
fastspec-iron-sv3-backup.fits   1,564             32
fastspec-iron-sv3-bright.fits   265,324           649
fastspec-iron-sv3-dark.fits     592,394           5,973
fastspec-iron.fits              1,397,479         13,991
=============================== ================= ===============================

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
desiutil         3.3.1
desispec         0.59.2
desitarget       2.6.0
desimodel        0.18.0
speclite         0.16
fastspecfit      2.4.1, 2.4.2
================ ==========

.. rst-class:: columns

==================== =====
Environment Variable Value
==================== =====
DESI_ROOT            /dvs_ro/cfs/cdirs/desi
DUST_DIR             /dvs_ro/cfs/cdirs/cosmo/data/dust/v0_1
FPHOTO_DIR           /dvs_ro/cfs/cdirs/desi/external/legacysurvey/dr9
FTEMPLATES_DIR       /dvs_ro/cfs/cdirs/desi/science/gqp/templates/fastspecfit
FTEMPLATES_FILE      ftemplates-chabrier-1.1.0.fits (see `README.txt`_)
FPHOTO_FILE          /global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/code/fastspecfit/2.4.1/lib/python3.10/site-packages/fastspecfit/data/legacysurvey-dr9.yaml
EMLINES_FILE         /global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/code/fastspecfit/2.4.1/lib/python3.10/site-packages/fastspecfit/data/emlines.ecsv
==================== =====

Known Issues
------------

This section documents any issues or problems which were identified with the VAC
after its final release. To report additional problems or to request new
features please `open a ticket`_. 

* Fluxes (and EWs) of lines which lie in the camera-overlap region are
  overestimated by a factor of 2 due to a bug handling the different pixel scale
  (see `issue/#157`_).
* Stellar masses are systematically higher (by 0.2-0.5 dex) compared to other
  methods, so they should be used with caution; see `issue/#159`_. Similarly,
  star-formation rates have not been fully validated.

.. _`issue/#157`: https://github.com/desihub/fastspecfit/issues/157
.. _`issue/#159`: https://github.com/desihub/fastspecfit/issues/159

.. _previous versions - iron:

Previous Versions
-----------------

In this section we document the version of ``FastSpecFit`` used to generate
previous, earlier versions of this VAC. Please see the `change log`_ for a
record of what code and data model changes have occurred since these previous
versions were released.

.. rst-class:: columns

=========== ======================
VAC Version FastSpecFit Version(s)
=========== ======================
v1.0        2.1.0, 2.1.1
=========== ======================

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`issue/#159`: https://github.com/desihub/fastspecfit/issues/159
.. _`README.txt`: https://data.desi.lbl.gov/desi/public/external/templates/fastspecfit/README.txt

