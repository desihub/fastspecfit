.. _guadalupe vac:

Guadalupe VAC (DR1 Supplement)
==============================

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Guadalupe`` value-added catalog, which will be
publicly released as part of the `DESI Data Release 1 (DESI/DR1)`_ sometime in
2024 (exact date TBD).

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`Previous Versions<previous versions - guadalupe>`.

Data Content & Access
---------------------

Data from the ``Guadalupe`` VAC can be accessed at any of the following links:

============================ =================================================================
Data url                     https://data.desi.lbl.gov/desi/spectro/fastspecfit/guadalupe/v3.0
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/spectro/fastspecfit/guadalupe/v3.0``
============================ =================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts>`

Summary Statistics
------------------
  
The following table summarizes the size and number of targets in each merged
catalog.

.. note::

   The last catalog listed in the table is a super-merge of all the preceding
   catalogs, i.e., a merge over all surveys and programs.

.. rst-class:: columns

====================================== ========= =================
File Name                              File Size Number of Targets
====================================== ========= =================
fastspec-guadalupe-special-dark.fits   17.4 MB   3,848
fastspec-guadalupe-special-bright.fits 43.1 MB   9,598
fastspec-guadalupe-main-bright.fits    4.77 GB   1,092,041
fastspec-guadalupe-main-dark.fits      4.94 GB   1,131,858
fastspec-guadalupe.fits                9.79 GB   2,237,345
====================================== ========= =================

The following table summarizes the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

.. rst-class:: columns

======================================== ================= ===============================
Catalog                                  Number of Objects Number with Corrected Redshifts
======================================== ================= ===============================
{fastspec}-guadalupe-main-bright.fits    1,092,041         156
{fastspec}-guadalupe-main-dark.fits      1,131,858         22,445
{fastspec}-guadalupe-special-bright.fits 9,598             0
{fastspec}-guadalupe-special-dark.fits   3,848             111
{fastspec}-guadalupe.fits                2,237,345         22,712
======================================== ================= ===============================

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

.. _previous versions - guadalupe:

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
v2.0        2.1.0, 2.1.1
v1.0        1.0.0, 1.0.1
=========== ======================

To report projects or to request new features please `open a ticket`_.

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/desi/public/external/templates/fastspecfit/README.txt

