"""Build a stacked LRG spectrum for use as a fastspecfit unit-test fixture.

This script stacks five sv3/dark LRG spectra from the ``loa`` spectroscopic
production into a single rest-frame spectrum and writes the result to
``test/data/``. The five targets span z ~ 0.41–0.45 and were selected by
hand to be representative LRGs with reliable fastspecfit solutions. Spectra
are shifted to the rest frame, normalized in a 50 Å window centered on
4500 Å (rest), and combined with inverse-variance weighting via
``desigal.specutils.stack.stack_spectra``.

Requirements
------------
``DESI_ROOT`` (or a valid ``SPECPROD`` environment) must be set so that
``desispec`` can locate the ``loa`` coadd files. ``desigal`` must be
installed (``pip install fastspecfit[testdata]``). This script is intended
to be run once by collaborators with NERSC access when the test target or
stacking parameters change; it must not be run as part of the automated
test suite.

Run from the ``test/`` directory::

    python build_test_stack.py

Outputs
-------
``test/data/stack-LRG.fits``
    Single-bin stacked spectrum in the format expected by ``stackfit``.

"""
import os
import warnings
import numpy as np
from astropy.table import Table
from desispec.io.spectra import read_spectra_parallel
from desigal.specutils.stack import stack_spectra, write_binned_stacks

SPECPROD = 'loa'
NORM_WINDOW = [4475.0, 4525.0]   # centered on ~4500 Å rest (LRG continuum)
DWAVE = 0.4                       # Å/pixel
OBSWAVE_MIN = 3600.0
OBSWAVE_MAX = 9800.0

TARGETS = Table({
    'TARGETID': np.array([
        39627763635720044,
        39627829507262695,
        39627921001810574,
        39633290050669693,
        39627775727900698,
    ], dtype=np.int64),
    'SURVEY':  ['sv3'] * 5,
    'PROGRAM': ['dark'] * 5,
    'HEALPIX': np.array([26275, 27247, 26091, 11425, 26276], dtype=np.int32),
    'Z': np.array([
        0.4052904043347254,
        0.411587106145884,
        0.4462853971996436,
        0.44740186572631313,
        0.4263606540076314,
    ]),
})


def build_stack(outdir='data'):
    """Stack five LRG spectra into a rest-frame composite and write to *outdir*."""
    zmin = TARGETS['Z'].min()
    zmax = TARGETS['Z'].max()
    restwave = np.arange(OBSWAVE_MIN / (1.0 + zmax), OBSWAVE_MAX / (1.0 + zmin), DWAVE)

    spectra = read_spectra_parallel(
        TARGETS,
        prefix='coadd',
        specprod=SPECPROD,
        nproc=1,
        rdspec_kwargs={
            'skip_hdus': ['EXP_FIBERMAP', 'SCORES', 'EXTRA_CATALOG', 'RESOLUTION'],
        },
    )
    assert np.all(spectra.target_ids() == TARGETS['TARGETID']), \
        'Spectra order does not match targets table.'

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        # stack_spectra returns ((flux, ivar), per-object weights); weights unused here
        (flux, ivar), _ = stack_spectra(
            spectra=spectra,
            redshift=TARGETS['Z'].data,
            stack_redshift=0.0,
            norm_flux_window=NORM_WINDOW,
            norm_method='flux-window',
            resample_method='linear',
            stack_method='ivar-weighted-mean',
            output_wave_grid=restwave,
            n_workers=1,
        )

    igood = np.where(np.isfinite(flux) & np.isfinite(ivar))[0]
    snr = float(np.median(flux[igood] * np.sqrt(ivar[igood])))

    outmeta = Table({
        'IBIN':    np.array([0], dtype=np.int32),
        'NALLOBJ': np.array([len(TARGETS)], dtype=np.int32),
        'NOBJ':    np.array([len(TARGETS)], dtype=np.int32),
        'STACKID': np.array([0], dtype=np.int32),
        'Z':       np.array([0.0]),
        'ZMED':    np.array([np.median(TARGETS['Z'])], dtype=np.float32),
        'SNR':     np.array([snr], dtype=np.float32),
    })

    outfile = os.path.join(outdir, 'stack-LRG.fits')
    write_binned_stacks(
        outfile,
        restwave,
        flux[np.newaxis, :],   # write_binned_stacks expects (nbins, npix)
        ivar[np.newaxis, :],
        resolution=None,
        stackinfo=outmeta,
    )
    print(f'Wrote {outfile}  (SNR={snr:.2f}, {len(igood)} good pixels)')


if __name__ == '__main__':
    build_stack()
