.. _fuji vac:

Fuji (Early Data Release)
=========================

.. contents:: Contents
    :depth: 4

Overview
--------

This page describes the `Fuji` value-added catalogs, which will be publicly
released as part of the `DESI Early Data Release (DESI/EDR)`_ in
late-spring 2023.

.. note::

   We only document the **v2.0** version of this VAC; the **v1.0** version is
   deprecated and has not been publicly released.

Data Access
-----------

The Fuji VACs can be accessed at the following url:

========================== ===================================================================
Value-Added Catalog        URL
========================== ===================================================================
Fuji (EDR)                 https://data.desi.lbl.gov/public/edr/vac/fastspecfit/fuji/v2.0
========================== ===================================================================

.. highlight:: none

.. note::

   DESI Collaborators may access the catalogs directly at `NERSC`_ at the
   following directory::
  
     /global/cfs/cdirs/desi/public/edr/vac/fastspecfit/fuji/v2.0

.. highlight:: default



Within the data release directory, there are two key subdirectories, `healpix`
and `catalogs`, which we now describe in more detail.

Healpix Catalogs
~~~~~~~~~~~~~~~~

We run ``FastSpecFit`` on the healpix-coadded DESI spectra, and organize the
files identically to how the spectra, redshift catalogs, and other data products
are organized in the DESI data releases (as documented `here`_). In other words,
for a given spectroscopic production ``SPECPROD=fuji``), the individual
``fastspec`` files (see :ref:`algorithms`) can be found at the following
location::

  VACDIR=/global/cfs/cdirs/desi/public/edr/vac/fastspecfit/fuji/v2.0
  $VACDIR/healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/fastspec-SURVEY-PROGRAM-HEALPIX.fits.gz

where ``SURVEY``, ``PROGRAM``, ``HPIXGROUP``, and ``HEALPIX`` are fully
documented `here`_.

.. note::

   The ``fastspec`` catalogs are gzipped because they contain the fitting
   results as well as the best-fitting model spectra; see the :ref:`fastspec
   data model<fastspec datamodel>` pages for a full description of the contents
   of these files.

Merged Catalogs
~~~~~~~~~~~~~~~

Most users will be interested in the merged ``FastSpecFit`` catalogs, which we
summarize in the tables below. Note that the last row of each table is a
super-merge of all the preceding catalogs (i.e., a merge over all possible
surveys and programs) listed in the table.

Fuji
""""

.. rst-class:: columns

=============================== ========= =================
File Name                       File Size Number of Targets
=============================== ========= =================
fastspec-fuji-cmx-other.fits    11.6 MB   2,771
fastspec-fuji-special-dark.fits 145 MB    35,647
fastspec-fuji-sv1-backup.fits   15.5 MB   3,683
fastspec-fuji-sv1-bright.fits   524 MB    126,677
fastspec-fuji-sv1-dark.fits     976 MB    235,881
fastspec-fuji-sv1-other.fits    141 MB    34,150
fastspec-fuji-sv2-backup.fits   602 KB    107
fastspec-fuji-sv2-bright.fits   193 MB    46,510
fastspec-fuji-sv2-dark.fits     218 MB    52,771
fastspec-fuji-sv3-backup.fits   6.62 MB   1,564
fastspec-fuji-sv3-bright.fits   1.08 GB   265,324
fastspec-fuji-sv3-dark.fits     2.41 GB   592,394
fastspec-fuji.fits              5.7 GB    1,397,479
=============================== ========= =================

.. note::

   In order to keep the size of the files reasonable, the `fastspec` files do
   not contain the ``MODELS`` FITS extension (see the :ref:`fastspec data
   model<fastspec datamodel>` page for a description of this FITS extension).

Updated QSO Redshifts
---------------------

.. rst-class:: columns

=============================== ================= ===============================
Catalog                         Number of Targets Number with Corrected Redshifts
=============================== ================= ===============================
fastspec-fuji-cmx-other.fits    2,771             56
fastspec-fuji-special-dark.fits 35,647            311
fastspec-fuji-sv1-backup.fits   3,683             100
fastspec-fuji-sv1-bright.fits   126,677           64
fastspec-fuji-sv1-dark.fits     235,881           3,749
fastspec-fuji-sv1-other.fits    34,150            168
fastspec-fuji-sv2-backup.fits   107               0
fastspec-fuji-sv2-bright.fits   46,510            8
fastspec-fuji-sv2-dark.fits     52,771            1,019
fastspec-fuji-sv3-backup.fits   1,564             0
fastspec-fuji-sv3-bright.fits   265,324           132
fastspec-fuji-sv3-dark.fits     592,394           3,397
fastspec-fuji.fits              1,397,479         9,004
=============================== ================= ===============================

Known Issues
------------

This section documents any issues or problems which were identified with these
VACs after their final release. To date, no major issues have been identified!
To report projects or to request new features please `open a ticket`_.

.. _`DESI Early Data Release (DESI/EDR)`: https://data.desi.lbl.gov/public/edr
.. _`DESI/EDR`: https://data.desi.lbl.gov/public/edr
.. _`NERSC`: https://nersc.gov
.. _`here`: https://data.desi.lbl.gov/doc/organization/
.. _`Redrock catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/redrock-SURVEY-PROGRAM-PIXNUM.html
.. _`quasarnet catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_qn-SURVEY-PROGRAM-PIXNUM.html
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`DESI Member Institutions`: https://www.desi.lbl.gov/collaborating-institutions
