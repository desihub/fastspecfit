.. _fujilupe vac:

Fuji and Guadalupe VACs
=======================

.. contents:: Contents
    :depth: 4

Overview
--------

This page describes the content and construction of the **version 1.0**
``FastSpecFit`` value-added catalogs (VACs) which were generated from the DESI
**Fuji** and **Guadalupe** spectroscopic productions. Both VACs will be publicly
released in the near future but are available to all DESI collaborators now:

  * Fuji will be publicly released as part of the `DESI Early Data Release
    (DESI/EDR)`_ in early 2023 (exact date TBD); and
  * Guadalupe will be released as a supplemental dataset to `DESI Data Release 1
    (DESI/DR1)`_ sometime in 2024 (exact date TBD).
    
Data Access and Organization
----------------------------

The Fuji and Guadalupe VACs can be accessed at the following urls:

========================== ===================================================================
Value-Added Catalog        URL
========================== ===================================================================
Fuji (EDR)                 https://data.desi.lbl.gov/public/edr/vac/fastspecfit/fuji/v1.0
Guadalupe (DR1 Supplement) https://data.desi.lbl.gov/public/dr1/vac/fastspecfit/guadalupe/v1.0
========================== ===================================================================

.. highlight:: none

.. note::

   DESI Collaborators may access the catalogs directly at `NERSC`_ at the
   following directories::
  
     /global/cfs/cdirs/desi/public/edr/vac/fastspecfit/fuji/v1.0
     /global/cfs/cdirs/desi/public/dr1/vac/fastspecfit/guadalupe/v1.0

.. highlight:: default

Within each data release directory, there are two key subdirectories, `healpix`
and `catalogs`, which we now describe in more detail. 

Healpix Catalogs
~~~~~~~~~~~~~~~~

We run ``FastSpecFit`` on the healpix-coadded DESI spectra, and organize the
files identically to how the spectra, redshift catalogs, and other data products
are organized in the DESI data releases (as documented `here`_). In other words,
for a given spectroscopic production ``SPECPROD={fuji, guadalupe}``), the
individual ``fastspec`` and ``fastphot`` files (see :ref:`algorithms`) can be
found at the following locations::

  healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/fastphot-SURVEY-PROGRAM-HEALPIX.fits
  healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/fastspec-SURVEY-PROGRAM-HEALPIX.fits.gz

where ``SURVEY``, ``PROGRAM``, ``HPIXGROUP``, and ``HEALPIX`` are fully
documented `here`_.

.. note::

   The ``fastspec`` catalogs are gzipped because they contain the fitting
   results as well as the best-fitting model spectra, whereas the ``fastphot``
   files only contain fitting results; see the :ref:`fastspec data
   model<fastspec datamodel>` and :ref:`fastphot data model<fastphot datamodel>`
   pages for a full description of the contents of these files.

.. _`merged catalogs`:

Merged Catalogs
~~~~~~~~~~~~~~~

Most users will be interested in the merged ``FastSpecFit`` catalogs, which we
summarize in the tables below, separately for the Fuji and Guadalupe
productions. Note that the last row of each table is a super-merge of all the
preceding catalogs (i.e., a merge over all possible surveys and programs) listed
in the table.

Fuji
""""

.. rst-class:: columns

=============================== ========= ================= =============================== ========= =================
File Name                       File Size Number of Targets File Name                       File Size Number of Targets
=============================== ========= ================= =============================== ========= =================
fastspec-fuji-cmx-other.fits    9.27 MB   2,771             fastphot-fuji-cmx-other.fits    1.82 MB   2,771
fastspec-fuji-special-dark.fits 119 MB    35,647            fastphot-fuji-special-dark.fits 24.6 MB   35,647
fastspec-fuji-sv1-backup.fits   12.4 MB   3,683             fastphot-fuji-sv1-backup.fits   2.56 MB   3,683
fastspec-fuji-sv1-bright.fits   419 MB    126,677           fastphot-fuji-sv1-bright.fits   82.7 MB   126,677
fastspec-fuji-sv1-dark.fits     780 MB    235,881           fastphot-fuji-sv1-dark.fits     154 MB    235,881
fastspec-fuji-sv1-other.fits    113 MB    34,150            fastphot-fuji-sv1-other.fits    22.2 MB   34,150
fastspec-fuji-sv2-backup.fits   498 KB    107               fastphot-fuji-sv2-backup.fits   101 KB    107
fastspec-fuji-sv2-bright.fits   154 MB    46,510            fastphot-fuji-sv2-bright.fits   30.6 MB   46,510
fastspec-fuji-sv2-dark.fits     175 MB    52,771            fastphot-fuji-sv2-dark.fits     34.6 MB   52,771
fastspec-fuji-sv3-backup.fits   5.31 MB   1,564             fastphot-fuji-sv3-backup.fits   1.06 MB   1,564
fastspec-fuji-sv3-bright.fits   883 MB    265,324           fastphot-fuji-sv3-bright.fits   179 MB    265,324
fastspec-fuji-sv3-dark.fits     1.92 GB   592,394           fastphot-fuji-sv3-dark.fits     400 MB    592,394
fastspec-fuji.fits              4.57 GB   1,397,479         fastphot-fuji.fits              970 MB    1,397,479
=============================== ========= ================= =============================== ========= =================

Guadalupe
"""""""""

.. rst-class:: columns

====================================== ========= ================= ====================================== ========= =================
File Name                              File Size Number of Targets File Name                              File Size Number of Targets
====================================== ========= ================= ====================================== ========= =================
fastspec-guadalupe-special-dark.fits   12.5 MB   3,847             fastphot-guadalupe-special-dark.fits   2.15 MB   3,847
fastspec-guadalupe-special-bright.fits 30.9 MB   9,598             fastphot-guadalupe-special-bright.fits 5.36 MB   9,598
fastspec-guadalupe-main-bright.fits    3.42 GB   1,092,038         fastphot-guadalupe-main-bright.fits    606 MB    1,092,038
fastspec-guadalupe-main-dark.fits      3.54 GB   1,131,601         fastphot-guadalupe-main-dark.fits      622 MB    1,131,601
fastspec-guadalupe.fits                7.02 GB   2,237,084         fastphot-guadalupe.fits                1.23 GB   2,237,084
====================================== ========= ================= ====================================== ========= =================

.. note::

   In order to keep the size of the files reasonable, the `fastspec` files do
   not contain the ``MODELS`` FITS extension (see the :ref:`fastspec data
   model<fastspec datamodel>` page for a description of this FITS extension).

Sample Selection
----------------

The sample selection---in other words, the criteria used the choose which DESI
targets to fit---were chosen to be very inclusive so that modeling results would
be available for as many objects as possible. In brief, we fit *all*
extragalactic (redshift greater than 0.001) non-sky (i.e., object) targets in
both Fuji and Guadalupe, with no cuts on targeting bits, redshift or
fiber-assignment warning bits, or other quality cuts. 

Specifically, let ``redrockfile`` be the full pathname to a given `redrock
catalog`_. The following bit of Python code illustrates which targets we fit:

.. code-block:: python

  import fitsio
  import numpy as np
  from fastspecfit.io import ZWarningMask

  zb = fitsio.read(redrockfile, 'REDSHIFTS')
  fm = fitsio.read(redrockfile, 'FIBERMAP')

  I = np.where((zb['Z'] > 0.001) * (fm['OBJTYPE'] == 'TGT') *
               (zb['ZWARN'] & ZWarningMask.NODATA == 0))[0]

where the ``ZWarningMask.NODATA`` bit indicates a spectrum which contains no
data (all inverse variance pixel values in the extracted spectrum are zero).

Updated QSO Redshifts
~~~~~~~~~~~~~~~~~~~~~

For a small but important fraction of quasar (QSO) targets, the redshift
determined by Redrock is incorrect. To mitigate this issue, the DESI team has
developed an approach to rectify the redshift nominally measured by Redrock
using the machine-learning algorithm ``QuasarNet``. In the Fuji and Guadalupe
``FastSpecFit`` VACs we adopt the same algorithm. 

Specifically, let ``redrockfile`` and ``qnfile`` be the full pathname to a given
`redrock catalog`_ and `QuasarNet catalog`_, respectively. We update the Redrock
redshift ``Z`` (and store the original Redrock redshift in ``Z_RR``; see the
:ref:`fastspec data model<fastspec datamodel>` and :ref:`fastphot data
model<fastphot datamodel>`) using the following bit of code:

.. code-block:: python

  import fitsio
  import numpy as np
  from astropy.table import Table

  zb = Table(fitsio.read(redrockfile, 'REDSHIFTS'))
  qn = Table(fitsio.read(qnfile, 'QN_RR'))

  QNLINES = ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']

  qn['IS_QSO_QN'] = np.max(np.array([qn[name] for name in linecols]), axis=0) > 0.95
  qn['IS_QSO_QN_NEW_RR'] &= qn['IS_QSO_QN']
  if np.count_nonzero(qn['IS_QSO_QN_NEW_RR']) > 0:
      zb['Z'][qn['IS_QSO_QN_NEW_RR']] = qn['Z_NEW'][qn['IS_QSO_QN_NEW_RR']]

For reference, the table below summarizes the number of objects with updated
redshifts in each of the Fuji and Guadalupe :ref:`merged catalogs`:

.. rst-class:: columns

========================================== ================= ===============================
Catalog                                    Number of Targets Number with Corrected Redshifts
========================================== ================= ===============================
{fastspec,fastphot}-fuji-cmx-other.fits    2,771             63
{fastspec,fastphot}-fuji-special-dark.fits 35,647            389
{fastspec,fastphot}-fuji-sv1-backup.fits   3,683             119
{fastspec,fastphot}-fuji-sv1-bright.fits   126,677           402
{fastspec,fastphot}-fuji-sv1-dark.fits     235,881           4,656
{fastspec,fastphot}-fuji-sv1-other.fits    34,150            372
{fastspec,fastphot}-fuji-sv2-backup.fits   107               0
{fastspec,fastphot}-fuji-sv2-bright.fits   46,510            151
{fastspec,fastphot}-fuji-sv2-dark.fits     52,771            1,185
{fastspec,fastphot}-fuji-sv3-backup.fits   1,564             32
{fastspec,fastphot}-fuji-sv3-bright.fits   265,324           649
{fastspec,fastphot}-fuji-sv3-dark.fits     592,394           5,973
{fastspec,fastphot}-fuji.fits              1,397,479         13,991
========================================== ================= ===============================

.. rst-class:: columns

================================================= ================= ===============================
Catalog                                           Number of Targets Number with Corrected Redshifts
================================================= ================= ===============================
{fastspec,fastphot}-guadalupe-main-bright.fits    1,092,038         2,080
{fastspec,fastphot}-guadalupe-main-dark.fits      1,131,601         26,741
{fastspec,fastphot}-guadalupe-special-bright.fits 9,598             13
{fastspec,fastphot}-guadalupe-special-dark.fits   3,847             121
{fastspec,fastphot}-guadalupe.fits                2,237,084         28,955
================================================= ================= ===============================

Known Issues
------------

This section documents any issues or problems which were identified with these
VACs after their final release. To date, no issues have been identified!


.. _`DESI Early Data Release (DESI/EDR)`: https://data.desi.lbl.gov/public/edr
.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`DESI/EDR`: https://data.desi.lbl.gov/public/edr
.. _`DESI/DR1`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`here`: https://data.desi.lbl.gov/doc/organization/
.. _`redrock catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/redrock-SURVEY-PROGRAM-PIXNUM.html
.. _`quasarnet catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_qn-SURVEY-PROGRAM-PIXNUM.html


