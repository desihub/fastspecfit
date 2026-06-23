"""Build miniature SFD dust maps for use in fastspecfit unit tests.

This script extracts small cutouts from the full-sky Schlegel, Finkbeiner &
Davis (1998) dust maps and writes them to the ``test/data/`` directory. The
output is the pixel-space union of bounding boxes around every entry in
:data:`REGIONS`, so all test objects are covered by a single pair of FITS
files. To add coverage for a new sky position, append a
``(ra_center, dec_center, border_deg)`` tuple to :data:`REGIONS`.

Both galactic-pole hemispheres (NGP and SGP) are always written because the
SFD map is split into two FITS files. Regions whose objects lie in the wrong
hemisphere receive a clipped, near-zero-size cutout that satisfies the
expected directory structure without consuming significant disk space.

Requirements
------------
``DUST_DIR`` must point to a directory containing the full-resolution SFD
maps at ``$DUST_DIR/maps/SFD_dust_4096_{ngp,sgp}.fits``.

Run once from the ``test/`` directory whenever the sky footprint changes;
not part of the automated test suite::

    python build_test_dustmaps.py

"""
import os
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


# ---------------------------------------------------------------------------
# Sky regions to cover. Each tuple is (ra_center, dec_center, border_deg).
# The output cutout is the pixel-space union of all bounding boxes.
# ---------------------------------------------------------------------------

REGIONS = [
    (105.5,         57.0,      2.0),   # iron tile 80613 fixture
    (149.4654042,   1.793647,  0.5),   # loa healpix / Suprime fixture
]


def _region_pixels(ra_cen, dec_cen, border_deg, crpix1, crpix2, lam_scal, sign):
    """Return pixel (x, y) arrays for the corners and center of a region.

    Samples five representative points (four corners plus the center) of the
    RA/Dec bounding box and maps them to Lambert equal-area pixel coordinates
    using the SFD projection formula.
    """
    ras  = [ra_cen - border_deg, ra_cen + border_deg,
            ra_cen - border_deg, ra_cen + border_deg,
            ra_cen]
    decs = [dec_cen - border_deg, dec_cen - border_deg,
            dec_cen + border_deg, dec_cen + border_deg,
            dec_cen]
    coords = SkyCoord(ras, decs, unit='degree', frame='icrs')
    ls = coords.galactic.l.radian
    bs = coords.galactic.b.radian
    r  = np.sqrt(1.0 - sign * np.sin(bs))
    xs = np.round(crpix1 - 1.0 + lam_scal * np.cos(ls) * r).astype(np.int32)
    ys = np.round(crpix2 - 1.0 - sign * lam_scal * np.sin(ls) * r).astype(np.int32)
    return xs, ys


for pole in ['ngp', 'sgp']:
    sign = 1 if pole == 'ngp' else -1

    dustmap = os.path.join(os.environ['DUST_DIR'], 'maps', f'SFD_dust_4096_{pole}.fits')
    data, hdr = fits.getdata(dustmap, header=True)
    wcs = WCS(hdr)

    crpix1   = hdr['CRPIX1']
    crpix2   = hdr['CRPIX2']
    lam_scal = hdr['LAM_SCAL']

    # Compute the pixel-space union of all region bounding boxes.
    all_xs, all_ys = [], []
    for ra_cen, dec_cen, border in REGIONS:
        xs, ys = _region_pixels(ra_cen, dec_cen, border, crpix1, crpix2, lam_scal, sign)
        all_xs.append(xs)
        all_ys.append(ys)

    all_xs = np.concatenate(all_xs)
    all_ys = np.concatenate(all_ys)

    x_lo = int(np.clip(all_xs.min(), 0, data.shape[1] - 1))
    x_hi = int(np.clip(all_xs.max(), 0, data.shape[1] - 1))
    y_lo = int(np.clip(all_ys.min(), 0, data.shape[0] - 1))
    y_hi = int(np.clip(all_ys.max(), 0, data.shape[0] - 1))

    # Guarantee a non-empty slice even when all regions fall in the other
    # hemisphere and their pixel coordinates collapse to a single point.
    if x_lo >= x_hi:
        x_hi = x_lo + 1
    if y_lo >= y_hi:
        y_hi = y_lo + 1

    wcs_cut  = wcs[y_lo:y_hi, x_lo:x_hi]
    data_cut = data[y_lo:y_hi, x_lo:x_hi]

    wcshdr = wcs_cut.to_header()
    for key in wcshdr.keys():
        hdr[key] = wcshdr[key]

    outfile = os.path.join('data', os.path.basename(dustmap))
    print(f'Writing {outfile}  [{y_lo}:{y_hi}, {x_lo}:{x_hi}]')
    fits.writeto(outfile, data_cut, header=hdr, overwrite=True)
