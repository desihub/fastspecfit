.. _loa vac:
.. _loa fastspec vac:

Loa fastspec VAC (DR2)
======================

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Loa`` fastspec value-added catalog, which contains
full spectrophotometric fitting results (continuum + emission lines) and will
be publicly released as part of `DESI Data Release 2 (DESI/DR2)`_. For the
companion photometry-only catalog, see the :ref:`Loa fastphot VAC<loa fastphot vac>`.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - loa>`.

Data Content & Access
---------------------

Data from the ``Loa`` fastspec VAC can be accessed at any of the following links:

============================ ===========================================================
Data url                     https://data.desi.lbl.gov/desi/vac/dr2/fastspecfit/loa/v1.0
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0``
============================ ===========================================================

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

=========================== ========= =================
File Name                   File Size Number of Targets
=========================== ========= =================
fastspec-loa-cmx-other.fits 12 MB     2,765
=========================== ========= =================

============================ =========== =================
File Name                    File Size   Number of Targets
============================ =========== =================
fastspec-loa-sv1-backup.fits 5.11 MB     1,145
fastspec-loa-sv1-bright.fits 542 MB      126,658
fastspec-loa-sv1-dark.fits   1.01e+03 MB 235,821
fastspec-loa-sv1-other.fits  145 MB      34,004
Total (sv1)                  1.66 GB     397,628
============================ =========== =================

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastspec-loa-sv2-backup.fits 652 KB    107
fastspec-loa-sv2-bright.fits 199 MB    46,491
fastspec-loa-sv2-dark.fits   225 MB    52,693
Total (sv2)                  425 MB    99,291
============================ ========= =================

============================ ========= =================
File Name                    File Size Number of Targets
============================ ========= =================
fastspec-loa-sv3-backup.fits 6.91 MB   1,573
fastspec-loa-sv3-bright.fits 1.11 GB   265,307
fastspec-loa-sv3-dark.fits   2.48 GB   592,128
Total (sv3)                  3.6 GB    859,008
============================ ========= =================

================================ ========= =================
File Name                        File Size Number of Targets
================================ ========= =================
fastspec-loa-special-backup.fits 2.32 MB   512
fastspec-loa-special-bright.fits 450 MB    105,172
fastspec-loa-special-dark.fits   1.15 GB   261,706
Total (special)                  1.6 GB    367,390
================================ ========= =================

========================================= ========= =================
File Name                                 File Size Number of Targets
========================================= ========= =================
fastspec-loa-main-backup.fits             252 MB    60,372
fastspec-loa-main-bright-nside1-hp00.fits 1.33 GB   325,494
fastspec-loa-main-bright-nside1-hp01.fits 11 GB     2,693,426
fastspec-loa-main-bright-nside1-hp02.fits 12.9 GB   3,141,532
fastspec-loa-main-bright-nside1-hp03.fits 2.63 GB   641,198
fastspec-loa-main-bright-nside1-hp04.fits 9.36 GB   2,285,745
fastspec-loa-main-bright-nside1-hp05.fits 1.67 GB   408,754
fastspec-loa-main-bright-nside1-hp06.fits 9.67 GB   2,359,953
fastspec-loa-main-bright-nside1-hp07.fits 3.07 GB   749,725
fastspec-loa-main-bright-nside1-hp08.fits 1.26 GB   307,648
fastspec-loa-main-bright-nside1-hp09.fits 292 MB    69,672
fastspec-loa-main-bright-nside1-hp10.fits 254 MB    60,466
fastspec-loa-main-bright-nside1-hp11.fits 595 MB    141,886
fastspec-loa-main-dark-nside1-hp00.fits   3.47 GB   846,678
fastspec-loa-main-dark-nside1-hp01.fits   16.1 GB   3,937,260
fastspec-loa-main-dark-nside1-hp02.fits   15.3 GB   3,730,928
fastspec-loa-main-dark-nside1-hp03.fits   2.37 GB   578,709
fastspec-loa-main-dark-nside1-hp04.fits   21 GB     5,112,856
fastspec-loa-main-dark-nside1-hp05.fits   3.54 GB   864,168
fastspec-loa-main-dark-nside1-hp06.fits   22.2 GB   5,426,022
fastspec-loa-main-dark-nside1-hp07.fits   6.62 GB   1,616,210
fastspec-loa-main-dark-nside1-hp08.fits   3.18 GB   775,136
fastspec-loa-main-dark-nside1-hp09.fits   657 MB    156,454
fastspec-loa-main-dark-nside1-hp10.fits   588 MB    140,162
fastspec-loa-main-dark-nside1-hp11.fits   372 MB    88,647
Total (main)                              150 GB    36,519,101
========================================= ========= =================

The following tables summarize the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

=========================== ================= ===============================
Catalog                     Number of Objects Number with Corrected Redshifts
=========================== ================= ===============================
fastspec-loa-cmx-other.fits 2,765             35
=========================== ================= ===============================

============================ ================= ===============================
Catalog                      Number of Objects Number with Corrected Redshifts
============================ ================= ===============================
fastspec-loa-sv1-backup.fits 1,145             17
fastspec-loa-sv1-bright.fits 126,658           34
fastspec-loa-sv1-dark.fits   235,821           2,958
fastspec-loa-sv1-other.fits  34,004            71
Total (sv1)                  397,628           3,080
============================ ================= ===============================

============================ ================= ===============================
Catalog                      Number of Objects Number with Corrected Redshifts
============================ ================= ===============================
fastspec-loa-sv2-backup.fits 107               0
fastspec-loa-sv2-bright.fits 46,491            2
fastspec-loa-sv2-dark.fits   52,693            674
Total (sv2)                  99,291            676
============================ ================= ===============================

============================ ================= ===============================
Catalog                      Number of Objects Number with Corrected Redshifts
============================ ================= ===============================
fastspec-loa-sv3-backup.fits 1,573             0
fastspec-loa-sv3-bright.fits 265,307           42
fastspec-loa-sv3-dark.fits   592,128           2,465
Total (sv3)                  859,008           2,507
============================ ================= ===============================

================================ ================= ===============================
Catalog                          Number of Objects Number with Corrected Redshifts
================================ ================= ===============================
fastspec-loa-special-backup.fits 512               0
fastspec-loa-special-bright.fits 105,172           1
fastspec-loa-special-dark.fits   261,706           242
Total (special)                  367,390           243
================================ ================= ===============================

========================================= ================= ===============================
Catalog                                   Number of Objects Number with Corrected Redshifts
========================================= ================= ===============================
fastspec-loa-main-backup.fits             60,372            0
fastspec-loa-main-bright-nside1-hp00.fits 325,494           27
fastspec-loa-main-bright-nside1-hp01.fits 2,693,426         260
fastspec-loa-main-bright-nside1-hp02.fits 3,141,532         270
fastspec-loa-main-bright-nside1-hp03.fits 641,198           60
fastspec-loa-main-bright-nside1-hp04.fits 2,285,745         206
fastspec-loa-main-bright-nside1-hp05.fits 408,754           35
fastspec-loa-main-bright-nside1-hp06.fits 2,359,953         220
fastspec-loa-main-bright-nside1-hp07.fits 749,725           54
fastspec-loa-main-bright-nside1-hp08.fits 307,648           25
fastspec-loa-main-bright-nside1-hp09.fits 69,672            3
fastspec-loa-main-bright-nside1-hp10.fits 60,466            5
fastspec-loa-main-bright-nside1-hp11.fits 141,886           11
fastspec-loa-main-dark-nside1-hp00.fits   846,678           6,781
fastspec-loa-main-dark-nside1-hp01.fits   3,937,260         30,565
fastspec-loa-main-dark-nside1-hp02.fits   3,730,928         27,213
fastspec-loa-main-dark-nside1-hp03.fits   578,709           6,712
fastspec-loa-main-dark-nside1-hp04.fits   5,112,856         34,425
fastspec-loa-main-dark-nside1-hp05.fits   864,168           5,805
fastspec-loa-main-dark-nside1-hp06.fits   5,426,022         30,824
fastspec-loa-main-dark-nside1-hp07.fits   1,616,210         10,651
fastspec-loa-main-dark-nside1-hp08.fits   775,136           5,529
fastspec-loa-main-dark-nside1-hp09.fits   156,454           932
fastspec-loa-main-dark-nside1-hp10.fits   140,162           924
fastspec-loa-main-dark-nside1-hp11.fits   88,647            1,027
Total (main)                              36,519,101        162,564
========================================= ================= ===============================


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
desiutil         3.4.3
desispec         0.68.1
desitarget       2.8.0
speclite         0.20
fastspecfit      3.1.4, 3.1.5
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

.. _previous versions - loa:

Notes & Known Issues
--------------------

v1.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~~

* Release date: 2027 (Exact date TBD)
* ``FastSpecFit`` version: ``3.1.4``, ``3.1.5``
* Templates: ``ftemplates-chabrier-2.0.0.fits``  (see `README.txt`_).
* Notes:

  * Initial release of the ``Loa`` fastspec VAC.

* Known issues:

  * None at this time.

.. _`DESI Data Release 2 (DESI/DR2)`: https://data.desi.lbl.gov/doc/releases/dr2
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/public/external/templates/fastspecfit/README.txt
