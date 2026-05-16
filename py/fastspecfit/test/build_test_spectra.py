"""Build spectral and photometric test fixtures for fastspecfit unit tests.

This script extracts single-target cutouts from a DESI tile reduction and
writes them to the ``test/data/`` directory. The target is drawn from tile
80613 (petal 4, cumulative through night 20210324) of the ``iron``
spectroscopic production.

The following files are produced:

``test/data/redrock-4-80613-thru20210324.fits``
    One-row subset of the Redrock redshift catalog for
    TARGETID 39633345008634465.

``test/data/coadd-4-80613-thru20210324.fits``
    Corresponding one-fiber coadded spectrum.

``test/data/tiles-iron.csv``
    One-row subset of the tile summary CSV.

``test/data/north/tractor/105/tractor-1054p567.fits``
    Legacy Survey DR9 Tractor photometry row for the same target.

Requirements
------------
``DESI_ROOT`` must point to the DESI collaboration filesystem root (e.g.,
``/global/cfs/cdirs/desi`` on NERSC). This script is intended to be run once
by collaborators with access to the full DESI data when the test target or
data model changes; it must not be run as part of the automated test suite.

Run from the ``test/`` directory::

    python build_test_spectra.py
"""
import os
import numpy as np
import fitsio
from astropy.table import Table
from desispec.io import write_spectra, read_spectra

SPECPROD = 'iron'
NIGHT = 20210324
TILE = 80613
PETAL = 4
TARGETID = 39633345008634465


def build_spectra(outdir='data'):
    """Write one-target redrock, coadd, and tile-summary files to *outdir*.

    Returns the subset fibermap for the test target, which can be passed
    directly to :func:`build_photometry`.
    """
    os.environ['SPECPROD'] = SPECPROD

    reduxdir = os.path.join(os.environ['DESI_ROOT'], 'spectro', 'redux', SPECPROD)
    datadir = os.path.join(reduxdir, 'tiles', 'cumulative', str(TILE), str(NIGHT))

    coaddfile = os.path.join(datadir, f'coadd-{PETAL}-{TILE}-thru{NIGHT}.fits')
    redrockfile = os.path.join(datadir, f'redrock-{PETAL}-{TILE}-thru{NIGHT}.fits')

    redhdr = fitsio.read_header(redrockfile)
    zbest = Table.read(redrockfile, 'REDSHIFTS')
    fibermap = Table.read(redrockfile, 'FIBERMAP')
    expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
    tsnr2 = Table.read(redrockfile, 'TSNR2')

    zbest = zbest[np.isin(zbest['TARGETID'], TARGETID)]
    fibermap = fibermap[np.isin(fibermap['TARGETID'], TARGETID)]
    expfibermap = expfibermap[np.isin(expfibermap['TARGETID'], TARGETID)]
    tsnr2 = tsnr2[np.isin(tsnr2['TARGETID'], TARGETID)]

    out_redrockfile = os.path.join(outdir, os.path.basename(redrockfile))
    print(f'Writing {out_redrockfile}')
    with fitsio.FITS(out_redrockfile, 'rw', clobber=True) as ff:
        ff.write(None, header=redhdr)
        ff.write(zbest.as_array(), extname='REDSHIFTS')
        ff.write(fibermap.as_array(), extname='FIBERMAP')
        ff.write(expfibermap.as_array(), extname='EXP_FIBERMAP')
        ff.write(tsnr2.as_array(), extname='TSNR2')

    spec = read_spectra(coaddfile).select(targets=TARGETID)
    out_coaddfile = os.path.join(outdir, os.path.basename(coaddfile))
    print(f'Writing {out_coaddfile}')
    write_spectra(out_coaddfile, spec)

    tilesfile = os.path.join(reduxdir, f'tiles-{SPECPROD}.csv')
    tiles = Table.read(tilesfile)
    tiles = tiles[tiles['TILEID'] == TILE]
    out_tilesfile = os.path.join(outdir, os.path.basename(tilesfile))
    print(f'Writing {out_tilesfile}')
    tiles.write(filename=out_tilesfile, format='csv', overwrite=True)

    return fibermap


def build_photometry(fibermap, outdir='data'):
    """Write the Legacy Survey DR9 Tractor photometry row for the test target.

    Parameters
    ----------
    fibermap : :class:`astropy.table.Table`
        Subset fibermap for the test target, as returned by
        :func:`build_spectra`. Used to look up the Tractor brick.
    outdir : str, optional
        Output directory. Defaults to ``'data'``.
    """
    from desispec.io import photo

    dr9dir = os.path.join(os.environ['DESI_ROOT'], 'external', 'legacysurvey', 'dr9')

    tractor = photo.gather_tractorphot(fibermap, dr9dir=dr9dir)
    tractor.remove_columns(('TARGETID', 'LS_ID'))

    region = 'north'
    brick = tractor['BRICKNAME'][0]
    tractorfile = os.path.join(dr9dir, region, 'tractor', brick[:3], f'tractor-{brick}.fits')
    tractorhdr = fitsio.read_header(tractorfile)
    out_tractorfile = os.path.join(outdir, region, 'tractor', brick[:3], os.path.basename(tractorfile))

    print(f'Writing {out_tractorfile}')
    os.makedirs(os.path.dirname(out_tractorfile), exist_ok=True)
    fitsio.write(out_tractorfile, tractor.as_array(), header=tractorhdr, clobber=True)


if __name__ == '__main__':
    fibermap = build_spectra()
    build_photometry(fibermap)
