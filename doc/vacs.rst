.. _vacs:

VAC Overview
============

.. contents:: Contents
    :depth: 3

Description
-----------

The following pages describe the content and construction of the ``FastSpecFit``
value-added catalogs (VACs). For a complete list of available VACs organized by
data release, please see the :ref:`VAC index<vacs index>`.

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
the individual per-healpix files can be found at the following location::

  $VACDIR/healpix/SURVEY/PROGRAM/HPIXGROUP/HEALPIX/{fastspec,fastphot}-SURVEY-PROGRAM-HEALPIX.fits{.gz}

where ``SURVEY``, ``PROGRAM``, ``HPIXGROUP``, and ``HEALPIX`` are fully
documented `here`_.

.. note::

   The ``fastspec`` healpix catalogs are gzipped because they contain the
   best-fitting model spectra in addition to the fitting results (see the
   :ref:`fastspec data model<fastspec datamodel>`). The ``fastphot`` healpix
   catalogs are not gzipped because they do not include a ``MODELS`` extension
   (see the :ref:`fastphot data model<fastphot datamodel>`).

.. _`merged catalogs`:

Merged Catalogs
~~~~~~~~~~~~~~~

Most users will be interested in the merged ``FastSpecFit`` catalogs, which
combine all the individual `healpix` catalogs for a given ``SURVEY`` and
``PROGRAM`` into a single file. Relative to the top-level data release directory
``VACDIR``, the merged catalog filenames have the following syntax::

  $VACDIR/catalogs/{fastspec,fastphot}-SURVEY-PROGRAM.fits

.. note::

   In order to keep the size of the files reasonable, the merged ``fastspec``
   catalogs do not contain the ``MODELS`` FITS extension (see the :ref:`fastspec
   data model<fastspec datamodel>` page for a description of this FITS
   extension). The ``fastphot`` merged catalogs also omit the ``MODELS``
   extension, which is not produced by ``fastphot``.

.. note::

   For large survey-program combinations (e.g., ``main-bright`` and
   ``main-dark``), the merged catalogs are further subdivided by nside=1
   healpixel to keep individual file sizes manageable. In that case the
   filenames follow the pattern::

     {fastspec,fastphot}-SURVEY-PROGRAM-nside1-hpNN.fits

   where ``NN`` runs from ``00`` to ``11``, corresponding to the 12 nside=1
   healpixels covering the full sky. Individual VAC pages document which
   survey-program combinations are split in this way.

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
  from fastspecfit.util import ZWarningMask

  zb = fitsio.read(redrockfile, 'REDSHIFTS')
  fm = fitsio.read(redrockfile, 'FIBERMAP')

  I = (zb['Z'] > 1e-3) & (fm['OBJTYPE'] == 'TGT') & (zb['ZWARN'] & ZWarningMask.NODATA == 0)

Here, the ``ZWarningMask.NODATA`` bit indicates a spectrum which
contains no data (all inverse variance pixel values in the extracted
spectrum are zero). Navigate to the DESI `ZWARN bit definitions`_
documentation for more details.

.. _`qso redshifts`:

Updated QSO Redshifts
---------------------

For a small but important fraction of quasar (QSO) targets, the redshift
determined by Redrock is incorrect. To mitigate this issue, the DESI team has
developed an approach to rectify the redshift nominally measured by Redrock
using the machine-learning algorithm ``QuasarNet``.

We update the Redrock redshift ``Z`` for affected QSO targets using the
`QuasarNet catalog`_ and `MgII catalog`_ afterburners distributed alongside
each Redrock output file. The original Redrock redshift and ``ZWARNING``
bitmask are preserved in the ``Z_RR`` and ``ZWARN_RR`` output columns (see the
:ref:`fastspec data model<fastspec datamodel>`); note that only ``Z`` is
modified, not ``ZWARN``.

The update is applied separately to two classes of targets:

**Primary QSO targets** — objects assigned to a QSO targeting bit (or the
equivalent ``SV0_QSO`` / ``MINI_SV_QSO`` commissioning bits for CMX data) —
have their redshift corrected when two conditions are met: the QuasarNet
afterburner independently classifies the spectrum as a QSO
(``IS_QSO_QN_NEW_RR = True``), and the maximum QuasarNet line confidence
across the six lines Lyα, C IV, C III], Mg II, Hβ, and Hα exceeds a
threshold. For DR2 (Loa) this threshold is 0.99; for DR1 (Iron) it was 0.95.

**Variable WISE QSO secondary targets** (``WISE_VAR_QSO`` secondary targeting
bit, present in main-survey and some SV programs but not CMX data) are updated
under a broader criterion: ``IS_QSO_QN_NEW_RR`` must be set and at least one
of the following must hold — the Redrock spectral type is ``QSO``, the MgII
afterburner identifies the object as a QSO (``IS_QSO_MGII``), or the maximum
QuasarNet line confidence exceeds the threshold above.


Issues
------

For questions (or problems) regarding any of the VACs catalogs or their
construction, please `open a ticket`_ and/or contact `John Moustakas (Siena
University)`_.

.. _`ZWARN bit definitions`: https://desidatamodel.readthedocs.io/en/latest/bitmasks.html#zwarn-bit-definitions
.. _`here`: https://data.desi.lbl.gov/doc/organization/
.. _`Redrock catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/redrock-SURVEY-PROGRAM-PIXNUM.html
.. _`QuasarNet catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_qn-SURVEY-PROGRAM-PIXNUM.html
.. _`MgII catalog`: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/qso_mgii-SURVEY-PROGRAM-PIXNUM.html
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`John Moustakas (Siena University)`: mailto:jmoustakas@siena.edu
