.. _fuji vac:

Fuji VAC
========

.. contents:: Contents
    :depth: 4

Overview
--------

This page describes the content and construction of the ``FastSpecFit``
value-added catalogs (VACs) which were generated from the DESI **Fuji**
spectroscopic production. These catalogs will be publicly released as part of
the `DESI Early Data Release (DESI/EDR)`_ in April 2023.
    
Data Access and Organization
----------------------------

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

  healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/fastspec-SURVEY-PROGRAM-HEALPIX.fits.gz

where ``SURVEY``, ``PROGRAM``, ``HPIXGROUP``, and ``HEALPIX`` are fully
documented `here`_.

.. note::

   The ``fastspec`` catalogs are gzipped because they contain the fitting
   results as well as the best-fitting model spectra; see the :ref:`fastspec
   data model<fastspec datamodel>` pages for a full description of the contents
   of these files.

.. _`fuji merged catalogs`:

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
redshifts in the :ref:`fuji merged catalogs`:

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

Acknowledgements
----------------

For questions (or problems) regarding these catalogs or their construction,
please `open a ticket`_ and/or contact `John Moustakas`_. 

JM gratefully acknowledges funding support for this work from the
U.S. Department of Energy, Office of Science, Office of High Energy Physics
under Award Number DE-SC0020086. We also gratefully acknowledge important
contributions to the VACs presented herein from the following individuals:

* Arjun Dey (NSF's NOIRLab)
* Stephen Bailey (Lawrence Berkeley National Lab)
* Rebecca Canning (University of Portsmouth)
* Victoria Fawcett (Durham University)  
* Stephanie Juneau (NSF's NOIRLab)
* Ashod Khederlarian (University of Pittsburgh)
* Dustin Lang (Perimeter Institute of Theoretical Physics)
* Adam Myers (University of Wyoming)
* Jeffrey Newman (University of Pittsburgh)
* Ragadeepika Pucha (University of Arizona)
* Anand Raichoor (Lawrence Berkeley National Lab)
* Khaled Said (Australian National University)  
* David Setton (University of Pittsburgh)
* Benjamin Weaver (NSF's NOIRLab)

DESI research is supported by the Director, Office of Science, Office of High
Energy Physics of the U.S. Department of Energy under Contract
No. DE–AC02–05CH11231, and by the National Energy Research Scientific Computing
Center, a DOE Office of Science User Facility under the same contract;
additional support for DESI is provided by the U.S. National Science Foundation,
Division of Astronomical Sciences under Contract No. AST-0950945 to the NSF’s
National Optical-Infrared Astronomy Research Laboratory; the Science and
Technologies Facilities Council of the United Kingdom; the Gordon and Betty
Moore Foundation; the Heising-Simons Foundation; the French Alternative Energies
and Atomic Energy Commission (CEA); the National Council of Science and
Technology of Mexico (CONACYT); the Ministry of Science and Innovation of Spain
(MICINN), and by the `DESI Member Institutions`_.

.. _`DESI Early Data Release (DESI/EDR)`: https://data.desi.lbl.gov/public/edr
.. _`DESI/EDR`: https://data.desi.lbl.gov/public/edr
.. _`NERSC`: https://nersc.gov
.. _`here`: https://data.desi.lbl.gov/doc/organization/
.. _`redrock catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/redrock-SURVEY-PROGRAM-PIXNUM.html
.. _`quasarnet catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_qn-SURVEY-PROGRAM-PIXNUM.html
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`John Moustakas`: mailto:jmoustakas@siena.edu
.. _`DESI Member Institutions`: https://www.desi.lbl.gov/collaborating-institutions
