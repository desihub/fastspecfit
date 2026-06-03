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
fastphot-iron-cmx-other.fits TBD       TBD
============================ ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-sv1-backup.fits TBD       TBD
fastphot-iron-sv1-bright.fits TBD       TBD
fastphot-iron-sv1-dark.fits   TBD       TBD
fastphot-iron-sv1-other.fits  TBD       TBD
Total (sv1)                   TBD       TBD
============================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-sv2-backup.fits TBD       TBD
fastphot-iron-sv2-bright.fits TBD       TBD
fastphot-iron-sv2-dark.fits   TBD       TBD
Total (sv2)                   TBD       TBD
============================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-sv3-backup.fits TBD       TBD
fastphot-iron-sv3-bright.fits TBD       TBD
fastphot-iron-sv3-dark.fits   TBD       TBD
Total (sv3)                   TBD       TBD
============================= ========= =================

================================= ========= =================
File Name                         File Size Number of Targets
================================= ========= =================
fastphot-iron-special-backup.fits TBD       TBD
fastphot-iron-special-bright.fits TBD       TBD
fastphot-iron-special-dark.fits   TBD       TBD
fastphot-iron-special-other.fits  TBD       TBD
Total (special)                   TBD       TBD
================================= ========= =================

============================= ========= =================
File Name                     File Size Number of Targets
============================= ========= =================
fastphot-iron-main-backup.fits TBD      TBD
fastphot-iron-main-bright.fits TBD      TBD
fastphot-iron-main-dark.fits   TBD      TBD
Total (main)                   TBD      TBD
============================= ========= =================

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
python           TBD
numpy            TBD
scipy            TBD
astropy          TBD
yaml             TBD
matplotlib       TBD
fitsio           TBD
mpi4py           TBD
healpy           TBD
desiutil         TBD
desispec         TBD
desitarget       TBD
speclite         TBD
fastspecfit      TBD
================ ==========

.. rst-class:: columns

==================== ============================================================
Environment Variable Value
==================== ============================================================
DESI_ROOT            /global/cfs/cdirs/desi
DUST_DIR             /dvs_ro/cfs/cdirs/desi/external/dust/v0_1
FPHOTO_DIR           /dvs_ro/cfs/cdirs/desi/external/legacysurvey/dr9
FTEMPLATES_DIR       /dvs_ro/cfs/cdirs/desi/public/external/templates/fastspecfit
FTEMPLATES_FILE      TBD (see `README.txt`_)
==================== ============================================================

.. _previous versions - iron fastphot:

Notes & Known Issues
--------------------

v1.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~~

* Release date: TBD
* ``FastSpecFit`` version: TBD
* Templates: TBD (see `README.txt`_).
* Notes:

  * Initial public release of the ``Iron`` fastphot VAC.

* Known issues:

  * None at this time.

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/doc/releases/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/public/external/templates/fastspecfit/README.txt
