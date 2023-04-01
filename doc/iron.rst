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

Data Content & Access
---------------------

Data from the ``Iron`` VAC can be accessed at any of the following links:

============================ ===================================================================
Data url                     https://data.desi.lbl.gov/public/edr/vac/fastspecfit/iron/v1.0
`NERSC`_ (for collaborators) /global/cfs/cdirs/desi/public/dr1/vac/fastspecfit/iron/v1.0
============================ ===================================================================

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

Known Issues
------------

This section documents any issues or problems which were identified with the VAC
after its final release. To date, no major issues have been identified!

To report projects or to request new features please `open a ticket`_.

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues

