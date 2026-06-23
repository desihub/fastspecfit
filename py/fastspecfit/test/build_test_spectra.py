"""Build spectral and photometric test fixtures for fastspecfit unit tests.

This script extracts single-target cutouts from DESI spectroscopic reductions
and writes them to the ``test/data/`` directory. To add a new fixture, append
a dict to :data:`TARGETS` — no other code changes are required.

Current fixtures
----------------

**iron cumulative tile (TARGETID 39633345008634465)**

``test/data/redrock-4-80613-thru20210324.fits``
    One-row Redrock redshift catalog.
``test/data/coadd-4-80613-thru20210324.fits``
    One-fiber coadded spectrum.
``test/data/tiles-iron.csv``
    One-row tile summary CSV.
``test/data/north/tractor/105/tractor-1054p567.fits``
    Legacy Survey DR9 Tractor photometry row.

**loa healpix special/dark (TARGETID 39089837499744649)**

``test/data/redrock-special-dark-27247.fits``
    One-row Redrock redshift catalog.
``test/data/coadd-special-dark-27247.fits``
    One-fiber coadded spectrum.
``test/data/suprime-photinfo.fits``
    One-row external photometry catalog (PHOTINFO extension).

Requirements
------------
``DESI_ROOT`` must point to the DESI collaboration filesystem root (e.g.
``/global/cfs/cdirs/desi`` on NERSC). Run once from the ``test/`` directory
whenever test targets or the data model changes; not part of the automated
test suite::

    python build_test_spectra.py

"""
import os
import numpy as np
import fitsio
from astropy.table import Table
from desispec.io import write_spectra, read_spectra


# ---------------------------------------------------------------------------
# Target registry — add a dict here for each new test fixture.
#
# Required keys for all targets:
#   label       short identifier used in console output
#   specprod    DESI spectroscopic production name (e.g. 'iron', 'loa')
#   coadd_type  'cumulative' | 'healpix'  (extend _build_spectra as needed)
#   targetid    DESI TARGETID (integer)
#   phot_type   'legacysurvey' | 'external'
#
# Additional keys for coadd_type == 'cumulative':
#   tile, petal, night
#
# Additional keys for coadd_type == 'healpix':
#   survey, program, healpix
#
# Additional keys for phot_type == 'external':
#   phot_infile   absolute path to the source photometric catalog
#   phot_ext      FITS extension name containing the catalog
#   phot_outfile  basename of the output file written to outdir
# ---------------------------------------------------------------------------

TARGETS = [
    dict(
        label      = 'iron-tile',
        specprod   = 'iron',
        coadd_type = 'cumulative',
        tile       = 80613,
        petal      = 4,
        night      = 20210324,
        targetid   = 39633345008634465,
        phot_type  = 'legacysurvey',
    ),
    dict(
        label        = 'loa-suprime',
        specprod     = 'loa',
        coadd_type   = 'healpix',
        survey       = 'special',
        program      = 'dark',
        healpix      = 27247,
        targetid     = 39089837499744649,
        phot_type    = 'external',
        phot_infile  = '/global/cfs/cdirs/desi/users/ioannis/desihiz/suprime/v20260616/desi-suprime.fits',
        phot_ext     = 'PHOTINFO',
        phot_outfile = 'suprime-photinfo.fits',
    ),
]


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def build_spectra(target, outdir='data'):
    """Extract spectra and redshift catalog for *target*, writing to *outdir*.

    Returns the fibermap subset for the target (passed to
    :func:`build_photometry` when ``phot_type == 'legacysurvey'``).
    """
    coadd_type = target['coadd_type']
    if coadd_type == 'cumulative':
        return _build_tile_spectra(target, outdir)
    elif coadd_type == 'healpix':
        return _build_healpix_spectra(target, outdir)
    else:
        raise ValueError(f"Unsupported coadd_type '{coadd_type}' for target '{target['label']}'")


def build_photometry(target, fibermap=None, outdir='data'):
    """Write photometry for *target* to *outdir*.

    Parameters
    ----------
    target : dict
        Entry from :data:`TARGETS`.
    fibermap : :class:`astropy.table.Table` or None
        Fibermap subset returned by :func:`build_spectra`.  Required when
        ``phot_type == 'legacysurvey'``.
    outdir : str, optional
        Output directory.  Defaults to ``'data'``.
    """
    phot_type = target['phot_type']
    if phot_type == 'legacysurvey':
        _build_legacysurvey_photometry(target, fibermap, outdir)
    elif phot_type == 'external':
        _build_external_photometry(target, outdir)
    else:
        raise ValueError(f"Unsupported phot_type '{phot_type}' for target '{target['label']}'")


# ---------------------------------------------------------------------------
# Private helpers — one per (coadd_type, phot_type) combination
# ---------------------------------------------------------------------------

def _build_tile_spectra(target, outdir):
    specprod = target['specprod']
    tile     = target['tile']
    petal    = target['petal']
    night    = target['night']
    targetid = target['targetid']

    os.environ['SPECPROD'] = specprod

    reduxdir = os.path.join(os.environ['DESI_ROOT'], 'spectro', 'redux', specprod)
    datadir  = os.path.join(reduxdir, 'tiles', 'cumulative', str(tile), str(night))

    coaddfile   = os.path.join(datadir, f'coadd-{petal}-{tile}-thru{night}.fits')
    redrockfile = os.path.join(datadir, f'redrock-{petal}-{tile}-thru{night}.fits')

    redhdr      = fitsio.read_header(redrockfile)
    zbest       = Table.read(redrockfile, 'REDSHIFTS')
    fibermap    = Table.read(redrockfile, 'FIBERMAP')
    expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
    tsnr2       = Table.read(redrockfile, 'TSNR2')

    zbest       = zbest[np.isin(zbest['TARGETID'], targetid)]
    fibermap    = fibermap[np.isin(fibermap['TARGETID'], targetid)]
    expfibermap = expfibermap[np.isin(expfibermap['TARGETID'], targetid)]
    tsnr2       = tsnr2[np.isin(tsnr2['TARGETID'], targetid)]

    out_redrockfile = os.path.join(outdir, os.path.basename(redrockfile))
    print(f'Writing {out_redrockfile}')
    with fitsio.FITS(out_redrockfile, 'rw', clobber=True) as ff:
        ff.write(None, header=redhdr)
        ff.write(zbest.as_array(), extname='REDSHIFTS')
        ff.write(fibermap.as_array(), extname='FIBERMAP')
        ff.write(expfibermap.as_array(), extname='EXP_FIBERMAP')
        ff.write(tsnr2.as_array(), extname='TSNR2')

    spec = read_spectra(coaddfile).select(targets=targetid)
    out_coaddfile = os.path.join(outdir, os.path.basename(coaddfile))
    print(f'Writing {out_coaddfile}')
    write_spectra(out_coaddfile, spec)

    tilesfile     = os.path.join(reduxdir, f'tiles-{specprod}.csv')
    tiles         = Table.read(tilesfile)
    tiles         = tiles[tiles['TILEID'] == tile]
    out_tilesfile = os.path.join(outdir, os.path.basename(tilesfile))
    print(f'Writing {out_tilesfile}')
    tiles.write(filename=out_tilesfile, format='csv', overwrite=True)

    return fibermap


def _build_healpix_spectra(target, outdir):
    specprod = target['specprod']
    survey   = target['survey']
    program  = target['program']
    healpix  = target['healpix']
    targetid = target['targetid']

    reduxdir = os.path.join(os.environ['DESI_ROOT'], 'spectro', 'redux', specprod)
    datadir  = os.path.join(reduxdir, 'healpix', survey, program,
                            str(healpix // 100), str(healpix))

    stem        = f'{survey}-{program}-{healpix}'
    coaddfile   = os.path.join(datadir, f'coadd-{stem}.fits')
    redrockfile = os.path.join(datadir, f'redrock-{stem}.fits')

    redhdr      = fitsio.read_header(redrockfile)
    zbest       = Table.read(redrockfile, 'REDSHIFTS')
    fibermap    = Table.read(redrockfile, 'FIBERMAP')
    expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
    tsnr2       = Table.read(redrockfile, 'TSNR2')

    zbest       = zbest[np.isin(zbest['TARGETID'], targetid)]
    fibermap    = fibermap[np.isin(fibermap['TARGETID'], targetid)]
    expfibermap = expfibermap[np.isin(expfibermap['TARGETID'], targetid)]
    tsnr2       = tsnr2[np.isin(tsnr2['TARGETID'], targetid)]

    out_redrockfile = os.path.join(outdir, os.path.basename(redrockfile))
    print(f'Writing {out_redrockfile}')
    with fitsio.FITS(out_redrockfile, 'rw', clobber=True) as ff:
        ff.write(None, header=redhdr)
        ff.write(zbest.as_array(), extname='REDSHIFTS')
        ff.write(fibermap.as_array(), extname='FIBERMAP')
        ff.write(expfibermap.as_array(), extname='EXP_FIBERMAP')
        ff.write(tsnr2.as_array(), extname='TSNR2')

    spec = read_spectra(coaddfile).select(targets=targetid)
    out_coaddfile = os.path.join(outdir, os.path.basename(coaddfile))
    print(f'Writing {out_coaddfile}')
    write_spectra(out_coaddfile, spec)

    return fibermap


def _build_legacysurvey_photometry(target, fibermap, outdir):
    from desispec.io import photo

    dr9dir  = os.path.join(os.environ['DESI_ROOT'], 'external', 'legacysurvey', 'dr9')
    tractor = photo.gather_tractorphot(fibermap, dr9dir=dr9dir)
    tractor.remove_columns(('TARGETID', 'LS_ID'))

    region      = 'north'
    brick       = tractor['BRICKNAME'][0]
    tractorfile = os.path.join(dr9dir, region, 'tractor', brick[:3], f'tractor-{brick}.fits')
    tractorhdr  = fitsio.read_header(tractorfile)

    out_tractorfile = os.path.join(outdir, region, 'tractor', brick[:3],
                                   os.path.basename(tractorfile))
    print(f'Writing {out_tractorfile}')
    os.makedirs(os.path.dirname(out_tractorfile), exist_ok=True)
    fitsio.write(out_tractorfile, tractor.as_array(), header=tractorhdr, clobber=True)


def _build_external_photometry(target, outdir):
    targetid     = target['targetid']
    phot_infile  = target['phot_infile']
    phot_ext     = target['phot_ext']
    phot_outfile = os.path.join(outdir, target['phot_outfile'])

    phot = fitsio.read(phot_infile, ext=phot_ext)
    row  = phot[phot['TARGETID'] == targetid]
    if len(row) == 0:
        raise ValueError(f'TARGETID {targetid} not found in {phot_infile}[{phot_ext}]')

    print(f'Writing {phot_outfile}')
    fitsio.write(phot_outfile, row, extname=phot_ext, clobber=True)


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    for target in TARGETS:
        print(f"\n=== {target['label']} ===")
        fibermap = build_spectra(target)
        build_photometry(target, fibermap=fibermap)
