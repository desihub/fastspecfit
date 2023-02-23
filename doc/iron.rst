.. _iron vac:

Iron and Guadalupe VACs
=======================

.. contents:: Contents
    :depth: 4

Overview
--------

This page describes the content and construction of the ``FastSpecFit``
value-added catalogs (VACs) which were generated from the DESI **Iron** and
**Guadalupe** spectroscopic productions. Both VACs will be publicly released as
part of the `DESI Data Release 1 (DESI/DR1)`_ sometime in 2024 (exact date TBD).

Data Access and Organization
----------------------------

The Iron and Guadalupe VACs can be accessed at the following urls:

========================== ===================================================================
Value-Added Catalog        URL
========================== ===================================================================
Iron (DR1)                 https://data.desi.lbl.gov/public/edr/vac/fastspecfit/iron/v1.0
Guadalupe (DR1 Supplement) https://data.desi.lbl.gov/public/dr1/vac/fastspecfit/guadalupe/v2.0
========================== ===================================================================

.. highlight:: none

.. note::

   DESI Collaborators may access the catalogs directly at `NERSC`_ at the
   following directories::
  
     /global/cfs/cdirs/desi/public/dr1/vac/fastspecfit/iron/v1.0
     /global/cfs/cdirs/desi/public/dr1/vac/fastspecfit/guadalupe/v2.0

.. highlight:: default

Within each data release directory, there are two key subdirectories, `healpix`
and `catalogs`, which we now describe in more detail. 

Healpix Catalogs
~~~~~~~~~~~~~~~~

We run ``FastSpecFit`` on the healpix-coadded DESI spectra, and organize the
files identically to how the spectra, redshift catalogs, and other data products
are organized in the DESI data releases (as documented `here`_). In other words,
for a given spectroscopic production ``SPECPROD={iron, guadalupe}``), the
individual ``fastspec`` files (see :ref:`algorithms`) can be
found at the following location::

  healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/fastspec-SURVEY-PROGRAM-HEALPIX.fits.gz

where ``SURVEY``, ``PROGRAM``, ``HPIXGROUP``, and ``HEALPIX`` are fully
documented `here`_.

.. note::

   The ``fastspec`` catalogs are gzipped because they contain the fitting
   results as well as the best-fitting model spectra; see the :ref:`fastspec
   data model<fastspec datamodel>` pages for a full description of the contents
   of these files.

.. _`iron merged catalogs`:

Merged Catalogs
~~~~~~~~~~~~~~~

Most users will be interested in the merged ``FastSpecFit`` catalogs, which we
summarize in the tables below, separately for the Iron and Guadalupe
productions. Note that the last row of each table is a super-merge of all the
preceding catalogs (i.e., a merge over all possible surveys and programs) listed
in the table.

Iron
""""

.. rst-class:: columns

=============================== ========= =================
File Name                       File Size Number of Targets
=============================== ========= =================
fastspec-iron-cmx-other.fits    9.27 MB    2,771            
fastspec-iron-special-dark.fits 119 MB     35,647           
fastspec-iron-sv1-backup.fits   12.4 MB    3,683            
fastspec-iron-sv1-bright.fits   419 MB     126,677          
fastspec-iron-sv1-dark.fits     780 MB     235,881          
fastspec-iron-sv1-other.fits    113 MB     34,150           
fastspec-iron-sv2-backup.fits   498 KB     107              
fastspec-iron-sv2-bright.fits   154 MB     46,510           
fastspec-iron-sv2-dark.fits     175 MB     52,771           
fastspec-iron-sv3-backup.fits   5.31 MB    1,564            
fastspec-iron-sv3-bright.fits   883 MB     265,324          
fastspec-iron-sv3-dark.fits     1.92 GB    592,394          
fastspec-iron-main-bright.fits  3.42 GB    1,092,038        
fastspec-iron-main-dark.fits    3.54 GB    1,131,601        
fastspec-iron.fits              4.57 GB    1,397,479        
=============================== ========= =================

Guadalupe
"""""""""

.. rst-class:: columns

====================================== ========= =================
File Name                              File Size Number of Targets
====================================== ========= =================
fastspec-guadalupe-special-dark.fits   15.7 MB    3,847            
fastspec-guadalupe-special-bright.fits 38.9 MB    9,598            
fastspec-guadalupe-main-bright.fits    4.31 GB    1,092,038        
fastspec-guadalupe-main-dark.fits      4.46 GB    1,131,601        
fastspec-guadalupe.fits                8.83 GB    2,237,084        
====================================== ========= =================

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
both Iron and Guadalupe, with no cuts on targeting bits, redshift or
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
using the machine-learning algorithm ``QuasarNet``. In the Iron and Guadalupe
``FastSpecFit`` VACs we adopt the same algorithm. 

Specifically, let ``redrockfile`` and ``qnfile`` be the full pathname to a given
`redrock catalog`_ and `QuasarNet catalog`_, respectively. We update the Redrock
redshift ``Z`` (and store the original Redrock redshift in ``Z_RR``; see the
:ref:`fastspec data model<fastspec datamodel>`) for all QSO targets using the
following bit of code:

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
redshifts in each of the Iron and Guadalupe :ref:`iron merged catalogs`:

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

.. rst-class:: columns

====================================== ================= ===============================
Catalog                                Number of Targets Number with Corrected Redshifts
====================================== ================= ===============================
fastspec-guadalupe-main-bright.fits    1,092,038         153
fastspec-guadalupe-main-dark.fits      1,131,601         26,741
fastspec-guadalupe-special-bright.fits 9,598             13
fastspec-guadalupe-special-dark.fits   3,847             121
fastspec-guadalupe.fits                2,237,084         28,955
====================================== ================= ===============================

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

.. _`DESI Data Release 1 (DESI/DR1)`: https://data.desi.lbl.gov/public/dr1
.. _`DESI/DR1`: https://data.desi.lbl.gov/public/dr1
.. _`NERSC`: https://nersc.gov
.. _`here`: https://data.desi.lbl.gov/doc/organization/
.. _`redrock catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/redrock-SURVEY-PROGRAM-PIXNUM.html
.. _`quasarnet catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_qn-SURVEY-PROGRAM-PIXNUM.html
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`John Moustakas`: mailto:jmoustakas@siena.edu
.. _`DESI Member Institutions`: https://www.desi.lbl.gov/collaborating-institutions
