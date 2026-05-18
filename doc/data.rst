.. _data:

Working with FastSpecFit Data
==============================

.. contents:: Contents
    :depth: 2

This page is for users working with existing ``FastSpecFit`` catalogs — either
from a public DESI data release or produced locally. See the :ref:`VACs<vacs>`
page for links to the public catalogs, and the :ref:`fastspec<fastspec
datamodel>` and :ref:`fastphot<fastphot datamodel>` data model pages for
complete column descriptions and units. For worked examples, see the
:ref:`Tutorials<tutorials>` page.

Reading a Catalog
-----------------

The recommended way to open a ``FastSpecFit`` output file is
:func:`fastspecfit.io.read_fastspecfit`, which returns the catalog tables as
:class:`astropy.table.Table` objects::

   from fastspecfit.io import read_fastspecfit

   fastfile = '/path/to/fastspec-iron-main-dark.fits'
   meta, specphot, fastfit, coadd_type, fastphot = read_fastspecfit(fastfile)

The five return values are:

- ``meta`` — targeting metadata and photometry (the ``METADATA`` extension).
- ``specphot`` — spectrophotometric measurements such as synthesized magnitudes
  (the ``SPECPHOT`` extension).
- ``fastfit`` — continuum and emission-line fitting results, or ``None`` in
  ``fastphot`` mode (the ``FASTSPEC`` extension).
- ``coadd_type`` — coadd type string from the file header (e.g., ``healpix``).
- ``fastphot`` — ``True`` if the file was produced in photometry-only mode.

The two tables are row-matched.

.. note::

   The merged ``FastSpecFit`` catalogs distributed as part of the DESI VACs do
   not include the ``MODELS`` extension in order to keep file sizes manageable.
   Model spectra are available in the per-healpix output files produced during
   a production run. To load them::

      meta, specphot, fastfit, coadd_type, fastphot, models = \
          read_fastspecfit(fastfile, read_models=True)

   The ``models`` array has shape ``(N, 3, Nwave)``, where the three components
   are the stellar continuum, the smooth continuum correction, and the
   emission-line model. The corresponding wavelength array is encoded in the
   ``MODELS`` extension header::

      import numpy as np, fitsio
      models, hdr = fitsio.read(fastfile, 'MODELS', header=True)
      modelwave = hdr['CRVAL1'] + np.arange(hdr['NAXIS1']) * hdr['CDELT1']

Accessing Key Quantities
------------------------

Both ``fastfit`` and ``meta`` are standard :class:`astropy.table.Table` objects
that support column access, boolean masking, and cross-matching. A few
commonly used quantities:

.. code-block:: python

   import numpy as np
   from fastspecfit.io import read_fastspecfit

   meta, specphot, fastfit, coadd_type, fastphot = read_fastspecfit(fastfile)

   # Redshifts and stellar masses
   redshift = fastfit['Z']
   logmstar  = fastfit['LOGMSTAR']

   # H-alpha flux and signal-to-noise
   halpha_flux = fastfit['HALPHA_FLUX']
   halpha_snr  = halpha_flux * np.sqrt(fastfit['HALPHA_FLUX_IVAR'])

   # Select star-forming galaxies with significant H-alpha detections
   sf = (halpha_snr > 3) & (logmstar > 9)
   print(f'{sf.sum()} star-forming galaxies selected')

See the :ref:`fastspec data model<fastspec datamodel>` for a complete
description of every column.
