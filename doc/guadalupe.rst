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

   We only document the *v3.0* version of this VAC; the *v1.0* and *v2.0*
   versions have been deprecated and are not documented further.

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

Known Issues
------------

This section documents any issues or problems which were identified with the VAC
after its final release.

* The tied amplitudes of [OII]7320,7330 doublet were reversed in the line fitting [`#119`_].

.. _`#119`: https://github.com/desihub/fastspecfit/issues/119

To report projects or to request new features please `open a ticket`_.

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues

