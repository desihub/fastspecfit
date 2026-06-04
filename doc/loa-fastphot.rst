.. _loa fastphot vac:

Loa fastphot VAC (DR2)
======================

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Loa`` fastphot value-added catalog, which contains
photometry-only continuum fitting results and will be publicly released as part
of `DESI Data Release 2 (DESI/DR2)`_. For the companion full spectrophotometric
catalog, see the :ref:`Loa fastspec VAC<loa fastspec vac>`.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - loa fastphot>`.

Data Content & Access
---------------------

Data from the ``Loa`` fastphot VAC can be accessed at any of the following links:

============================ ========================================================
Data url                     https://data.desi.lbl.gov/desi/vac/dr2/fastphot/loa/v1.0
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/vac/dr2/fastphot/loa/v1.0``
============================ ========================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts loa fastphot>`

.. _`qso redshifts loa fastphot`:

Updated QSO Redshifts
---------------------

The QSO redshift corrections applied to this VAC are identical to those used in
the :ref:`Loa fastspec VAC<loa fastspec vac>`; refer to that page for the
algorithm description and per-catalog correction counts.

Summary Statistics
------------------

The following tables summarize the size and number of targets in each merged
catalog, organized by survey and program.

=========================== ========= =================
File Name                   File Size Number of Targets
=========================== ========= =================
fastphot-loa-cmx-other.fits 2.34 MB   2,765
=========================== ========= =================

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastphot-loa-sv1-backup.fits 1.04 MB   1,145
fastphot-loa-sv1-bright.fits 106 MB    126,658
fastphot-loa-sv1-dark.fits   198 MB    235,821
fastphot-loa-sv1-other.fits  28.5 MB   34,004
Total (sv1)                  334 MB    397,628
============================ ========= =================

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastphot-loa-sv2-backup.fits 132 KB    107
fastphot-loa-sv2-bright.fits 39.3 MB   46,491
fastphot-loa-sv2-dark.fits   44.4 MB   52,693
Total (sv2)                  83.8 MB   99,291
============================ ========= =================

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastphot-loa-sv3-backup.fits 1.37 MB   1,573
fastphot-loa-sv3-bright.fits 229 MB    265,307
fastphot-loa-sv3-dark.fits   511 MB    592,128
Total (sv3)                  741 MB    859,008
============================ ========= =================

================================ ========= =================
File Name                        File Size Number of Targets
================================ ========= =================
fastphot-loa-special-backup.fits 430 KB    512
fastphot-loa-special-bright.fits 87.7 MB   105,172
fastphot-loa-special-dark.fits   283 MB    261,706
Total (special)                  371 MB    367,390
================================ ========= =================

========================================= ========= =================
File Name                                 File Size Number of Targets
========================================= ========= =================
fastphot-loa-main-backup.fits             44.5 MB   60,372
fastphot-loa-main-bright-nside1-hp00.fits 247 MB    325,494
fastphot-loa-main-bright-nside1-hp01.fits 2 GB      2,693,426
fastphot-loa-main-bright-nside1-hp02.fits 2.33 GB   3,141,532
fastphot-loa-main-bright-nside1-hp03.fits 487 MB    641,198
fastphot-loa-main-bright-nside1-hp04.fits 1.69 GB   2,285,745
fastphot-loa-main-bright-nside1-hp05.fits 310 MB    408,754
fastphot-loa-main-bright-nside1-hp06.fits 1.75 GB   2,359,953
fastphot-loa-main-bright-nside1-hp07.fits 569 MB    749,725
fastphot-loa-main-bright-nside1-hp08.fits 234 MB    307,648
fastphot-loa-main-bright-nside1-hp09.fits 52.9 MB   69,672
fastphot-loa-main-bright-nside1-hp10.fits 45.9 MB   60,466
fastphot-loa-main-bright-nside1-hp11.fits 108 MB    141,886
fastphot-loa-main-dark-nside1-hp00.fits   645 MB    846,678
fastphot-loa-main-dark-nside1-hp01.fits   2.93 GB   3,937,260
fastphot-loa-main-dark-nside1-hp02.fits   2.78 GB   3,730,928
fastphot-loa-main-dark-nside1-hp03.fits   441 MB    578,709
fastphot-loa-main-dark-nside1-hp04.fits   3.8 GB    5,112,856
fastphot-loa-main-dark-nside1-hp05.fits   659 MB    864,168
fastphot-loa-main-dark-nside1-hp06.fits   4.04 GB   5,426,022
fastphot-loa-main-dark-nside1-hp07.fits   1.2 GB    1,616,210
fastphot-loa-main-dark-nside1-hp08.fits   591 MB    775,136
fastphot-loa-main-dark-nside1-hp09.fits   119 MB    156,454
fastphot-loa-main-dark-nside1-hp10.fits   107 MB    140,162
fastphot-loa-main-dark-nside1-hp11.fits   67.6 MB   88,647
Total (main)                              27.1 GB   36,519,101
========================================= ========= =================

Code & Template Versions
------------------------

The following tables document the code versions and environment variables used
to produce this VAC. For details regarding the revision history of
``FastSpecFit``, please see the `change log`_.

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
fastspecfit      3.2.0
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

.. _previous versions - loa fastphot:

Notes & Known Issues
--------------------

v1.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~~

* Release date: 2027 (Exact date TBD)
* ``FastSpecFit`` version: ``3.2.0``
* Templates: ``ftemplates-chabrier-2.0.0.fits``  (see `README.txt`_).
* Notes:

  * Initial release of the ``Loa`` fastphot VAC.

* Known issues:

  * None at this time.

.. _`DESI Data Release 2 (DESI/DR2)`: https://data.desi.lbl.gov/doc/releases/dr2
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/public/external/templates/fastspecfit/README.txt
