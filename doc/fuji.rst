.. _fuji vac:

Fuji VAC (EDR)
==============

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Fuji`` value-added catalog, which will be publicly
released as part of the `DESI Early Data Release (DESI/EDR)`_ in late-2023. 

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *v3.0* version of this VAC; the *v1.0* and *v2.0*
   versions have been deprecated and are not documented further.

Data Content & Access
---------------------

Data from the ``Fuji`` VAC can be accessed at any of the following links:

============================ ============================================================
Data url                     https://data.desi.lbl.gov/desi/spectro/fastspecfit/fuji/v3.0
Interactive web-app          https://fastspecfit.desi.lbl.gov
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/spectro/fastspecfit/fuji/v3.0``
============================ ============================================================

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

=============================== ========= =================
File Name                       File Size Number of Targets
=============================== ========= =================
fastspec-fuji-cmx-other.fits    12.8 MB   2,771
fastspec-fuji-special-dark.fits 161 MB    35,649
fastspec-fuji-sv1-backup.fits   17.1 MB   3,683
fastspec-fuji-sv1-bright.fits   579 MB    126,678
fastspec-fuji-sv1-dark.fits     1.05 GB   235,942
fastspec-fuji-sv1-other.fits    156 MB    34,152
fastspec-fuji-sv2-backup.fits   664 KB    107
fastspec-fuji-sv2-bright.fits   213 MB    46,510
fastspec-fuji-sv2-dark.fits     241 MB    52,781
fastspec-fuji-sv3-backup.fits   7.31 MB   1,564
fastspec-fuji-sv3-bright.fits   1.19 GB   265,325
fastspec-fuji-sv3-dark.fits     2.66 GB   592,441
fastspec-fuji.fits              6.29 GB   1,397,603
=============================== ========= =================

The following table summarizes the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

.. rst-class:: columns

================================= ================= ===============================
Catalog                           Number of Objects Number with Corrected Redshifts
================================= ================= ===============================
{fastspec}-fuji-cmx-other.fits    2,771             56
{fastspec}-fuji-special-dark.fits 35,649            313
{fastspec}-fuji-sv1-backup.fits   3,683             100
{fastspec}-fuji-sv1-bright.fits   126,678           65
{fastspec}-fuji-sv1-dark.fits     235,942           3,810
{fastspec}-fuji-sv1-other.fits    34,152            170
{fastspec}-fuji-sv2-backup.fits   107               0
{fastspec}-fuji-sv2-bright.fits   46,510            8
{fastspec}-fuji-sv2-dark.fits     52,781            1,029
{fastspec}-fuji-sv3-backup.fits   1,564             0
{fastspec}-fuji-sv3-bright.fits   265,325           133
{fastspec}-fuji-sv3-dark.fits     592,441           3,444
{fastspec}-fuji.fits              1,397,603         9,128
================================= ================= ===============================

.. _known issues:

Known Issues
------------

This section documents any issues or problems which were identified with the VAC
after its final release. So far, none have been identified!

To report projects or to request new features please `open a ticket`_.

.. _`DESI Early Data Release (DESI/EDR)`: https://data.desi.lbl.gov/public/edr
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
