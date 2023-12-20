.. _vacs:

VAC Overview
============

.. contents:: Contents
    :depth: 3

Description
-----------

The following pages describe the content and construction of the ``FastSpecFit``
value-added catalogs (VACs) for the following DESI Data Releases and Data
Assemblies:

* :ref:`Fuji (Early Data Release)<fuji vac>`
* :ref:`Iron (Data Release 1)<iron vac>`
* :ref:`Guadalupe (Data Release 1 Supplement)<guadalupe vac>`

.. note::

   A DESI *data release* has been *publicly* released whereas a *data assembly*
   is only available internally to DESI collaboration members; not all data
   assemblies become data releases.

.. _`data organization`:

Data Organization
-----------------

Within the data release directory of each VAC, there are two key subdirectories,
``healpix`` and ``catalogs``, which we now describe in more detail.

.. _`healpix catalogs`:

Healpix Catalogs
~~~~~~~~~~~~~~~~

We run ``FastSpecFit`` on the healpix-coadded DESI spectra, and organize the
files identically to how the spectra, redshift catalogs, and other data products
are organized in the DESI data releases and assemblies (as documented
`here`_). In other words, relative to the VAC data release directory ``VACDIR``,
the individual ``fastspec`` files (see :ref:`algorithms`) can be found at the
following location::

  $VACDIR/healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/fastspec-SURVEY-PROGRAM-HEALPIX.fits.gz

where ``SURVEY``, ``PROGRAM``, ``HPIXGROUP``, and ``HEALPIX`` are fully
documented `here`_.

.. note::

   The ``fastspec`` catalogs are gzipped because they contain the fitting
   results as well as the best-fitting model spectra; see the :ref:`fastspec
   data model<fastspec datamodel>` pages for a full description of the contents
   of these files.

.. _`merged catalogs`:

Merged Catalogs
~~~~~~~~~~~~~~~

Most users will be interested in the merged ``FastSpecFit`` catalogs, which
combine all the individual `healpix` catalogs for a given ``SURVEY`` and
``PROGRAM`` into a single file. Relative to the top-level data release directory
``VACDIR``, the merged catalog filenames have the following syntax::

  $VACDIR/catalogs/fastspec-SURVEY-PROGRAM.fits

.. note::

   In order to keep the size of the files reasonable, these files do not contain
   the ``MODELS`` FITS extension (see the :ref:`fastspec data model<fastspec
   datamodel>` page for a description of this FITS extension). 

.. _`sample selection`:

Sample Selection
----------------

The sample selection criteria were chosen to be very inclusive so that modeling
results would be available for as many objects as possible.

In brief, we fit all extragalactic (redshift greater than :math:`10^{-3}`),
non-sky targets with no cuts on targeting bits, redshift or fiber-assignment
warning bits other than a simple cut which throws out spectra with no data.

Specifically, let ``redrockfile`` be the full pathname to a given `Redrock
catalog`_. The following bit of Python code illustrates which targets we fit:

.. code-block:: python

  import fitsio
  import numpy as np
  from fastspecfit.io import ZWarningMask

  zb = fitsio.read(redrockfile, 'REDSHIFTS')
  fm = fitsio.read(redrockfile, 'FIBERMAP')

  I = np.where((zb['Z'] > 0.001) * (fm['OBJTYPE'] == 'TGT') *
               (zb['ZWARN'] & ZWarningMask.NODATA == 0))[0]

Here, the ``ZWarningMask.NODATA`` bit indicates a spectrum which contains no
data (all inverse variance pixel values in the extracted spectrum are zero).

.. _`qso redshifts`:

Updated QSO Redshifts
---------------------

For a small but important fraction of quasar (QSO) targets, the redshift
determined by Redrock is incorrect. To mitigate this issue, the DESI team has
developed an approach to rectify the redshift nominally measured by Redrock
using the machine-learning algorithm ``QuasarNet``.

Let ``redrockfile`` and ``qnfile`` be the full pathname to a given `Redrock
catalog`_ and `QuasarNet catalog`_, respectively. We update the Redrock redshift
``Z`` (and store the original Redrock redshift in ``Z_RR``; see the
:ref:`fastspec data model<fastspec datamodel>`) for all QSO targets using the
following bit of code:

.. code-block:: python

  import fitsio
  import numpy as np
  from astropy.table import Table
  from desitarget.targets import main_cmx_or_sv

  QNLINES = ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']
  QNCOLS = ['TARGETID', 'Z_NEW', 'IS_QSO_QN_NEW_RR', 'C_LYA', 'C_CIV',
            'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']

  zb = Table(fitsio.read(redrockfile, 'REDSHIFTS'))

  # find QSO targets
  surv_target, surv_mask, surv = main_cmx_or_sv(meta)
  if surv == 'cmx':
      desi_target = surv_target[0]
      desi_mask = surv_mask[0]
      # need to check multiple QSO masks
      IQSO = []
      for bitname in desi_mask.names():
          if 'QSO' in bitname:
              IQSO.append(np.where(meta[desi_target] & desi_mask[bitname] != 0)[0])
      if len(IQSO) > 0:
          IQSO = np.sort(np.unique(np.hstack(IQSO)))
  else:
      desi_target, bgs_target, mws_target = surv_target
      desi_mask, bgs_mask, mws_mask = surv_mask
      IQSO = np.where(meta[desi_target] & desi_mask['QSO'] != 0)[0]

  if len(IQSO) > 0:
      qn = Table(fitsio.read(qnfile, 'QN_RR', columns=QNCOLS))
      assert(np.all(qn['TARGETID'] == zb['TARGETID'][IQSO]))
      print('Updating QSO redshifts using a QN threshold of 0.95.')
      qn['IS_QSO_QN'] = np.max(np.array([qn[name] for name in QNLINES]), axis=0) > 0.95
      qn['IS_QSO_QN_NEW_RR'] &= qn['IS_QSO_QN']
      if np.count_nonzero(qn['IS_QSO_QN_NEW_RR']) > 0:
          zb['Z'][IQSO[qn['IS_QSO_QN_NEW_RR']]] = qn['Z_NEW'][qn['IS_QSO_QN_NEW_RR']]

Issues
------

For questions (or problems) regarding any of the VACs catalogs or their
construction, please `open a ticket`_ and/or contact `John Moustakas (Siena
College)`_.

.. _`here`: https://data.desi.lbl.gov/doc/organization/
.. _`Redrock catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/redrock-SURVEY-PROGRAM-PIXNUM.html
.. _`quasarnet catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_qn-SURVEY-PROGRAM-PIXNUM.html
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`John Moustakas (Siena College)`: mailto:jmoustakas@siena.edu
