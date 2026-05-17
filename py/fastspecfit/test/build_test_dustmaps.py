"""Build miniature SFD dust maps for use in fastspecfit unit tests.

This script extracts small cutouts from the full-sky Schlegel, Finkbeiner &
Davis (1998) dust maps and writes them to the ``test/data/`` directory. The
cutouts cover a generous border around DESI tile 80613, spanning roughly
RA 104–107 deg, Dec 56–58 deg. Both galactic-pole hemispheres (NGP and SGP)
are written because the SFD map is split into two FITS files; the NGP file
covers the region of interest for tile 80613, while the SGP file is included
so the test data directory mirrors the expected directory structure.

Requirements
------------
The ``DUST_DIR`` environment variable must point to a directory containing the
full-resolution SFD maps at ``$DUST_DIR/maps/SFD_dust_4096_{ngp,sgp}.fits``.

Outputs
-------
``test/data/SFD_dust_4096_ngp.fits`` and ``test/data/SFD_dust_4096_sgp.fits``:
    Cropped FITS images with WCS headers updated to reflect the subarray
    extent. These files are checked into the repository and consumed by
    the photometry unit tests via the ``DUST_DIR`` test fixture.

Notes
-----
Run this script once from the ``test/`` directory whenever the test target or
sky footprint changes. It must not be run as part of the automated test suite.
The pixel-coordinate arithmetic follows the same Lambert equal-area projection
used internally by ``desiutil.dust.SFDMap``.
"""
import os
import numpy as np

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

# Large border around tile 80613.
ra1, dec1 = 104.0, 56.0
ra2, dec2 = 107.0, 58.0

c1 = SkyCoord(ra1, dec1, unit='degree', frame='icrs')
c2 = SkyCoord(ra2, dec2, unit='degree', frame='icrs')

l1, b1 = c1.galactic.l.radian, c1.galactic.b.radian
l2, b2 = c2.galactic.l.radian, c2.galactic.b.radian

for pole in ['ngp', 'sgp']:
    dustmap = os.path.join(os.environ.get('DUST_DIR'), 'maps', 'SFD_dust_4096_{}.fits'.format(pole))
    data, hdr = fits.getdata(dustmap, header=True)

    wcs = WCS(hdr)

    crpix1 = hdr['CRPIX1']
    crpix2 = hdr['CRPIX2']
    lam_scal = hdr['LAM_SCAL']
    sign = hdr['LAM_NSGP']  # north = 1, south = -1

    # Get the pixel indices corresponding to the RA, Dec boundaries listed
    # above. Code taken from desiutil.dust.SFDMap.ebv and
    # desiutil.dust._Hemisphere.
    x1 = np.round((crpix1 - 1.0 + lam_scal * np.cos(l1) * np.sqrt(1.0 - sign * np.sin(b1)))).astype(np.int32)
    x2 = np.round((crpix1 - 1.0 + lam_scal * np.cos(l2) * np.sqrt(1.0 - sign * np.sin(b2)))).astype(np.int32)

    y1 = np.round((crpix2 - 1.0 - sign * lam_scal * np.sin(l1) * np.sqrt(1.0 - sign * np.sin(b1)))).astype(np.int32)
    y2 = np.round((crpix2 - 1.0 - sign * lam_scal * np.sin(l2) * np.sqrt(1.0 - sign * np.sin(b2)))).astype(np.int32)

    x1 = np.clip(x1, 0, data.shape[1]-1)
    x2 = np.clip(x2, 0, data.shape[1]-1)
    y1 = np.clip(y1, 0, data.shape[0]-1)
    y2 = np.clip(y2, 0, data.shape[0]-1)

    # slice
    wcs = wcs[y2:y1, x1:2]
    data = data[y2:y1, x1:x2]

    wcshdr = wcs.to_header()
    for key in wcshdr.keys():
        hdr[key] = wcshdr[key]

    outfile = os.path.join('data', os.path.basename(dustmap))
    print('Writing {}'.format(outfile))
    fits.writeto(outfile, data, header=hdr, overwrite=True)
