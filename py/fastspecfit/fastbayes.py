"""
fastspecfit.fastbayes
=====================

Grid-based Bayesian broadband-photometric SED fitting.

Given a pre-synthesized photometry grid (see ``bin/build-bayesian-templates``
and ``bin/build-bayesian-photometry``), fit each object's observed broadband
photometry by interpolating the grid to the object's (Redrock) redshift,
solving for the chi2-minimizing stellar-mass amplitude of every grid template
in closed form, and building up the posterior probability distribution of
every grid parameter by weighting each template by its likelihood. Point
estimates (maximum-likelihood, mean, mode, and 16/50/84th percentiles) are
reported for every parameter; a sparse top-K representation of the joint
posterior can optionally be written out as well.

IGM attenuation and cosmological dimming are already baked into the
photometry grid (see ``bin/build-bayesian-photometry``), so this module never
needs the stellar template basis, the IGM model, or the cosmology directly.

"""
import os, sys, time, logging
import numpy as np
import fitsio
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.util import MPPool, fsftime
from fastspecfit.photometry import Photometry
from fastspecfit.cosmo import TabulatedDESI
from fastspecfit.singlecopy import sc_data

# Grid parameters for which point-estimate statistics are reported.
# LOGMASS is derived per-template from the closed-form amplitude solve and
# treated identically to the native grid parameters.
PARAM_NAMES = ('AGE', 'LOGMET', 'DUST', 'UVB', 'UMIN', 'GAMMA', 'QPAH', 'LOGMASS')

# Grid-column name (lower-case, as written by bin/build-bayesian-templates)
# for each reported parameter; LOGMASS has no grid column (computed per-object).
_GRID_COLUMN = {
    'AGE': 'age', 'LOGMET': 'logmet', 'DUST': 'dust', 'UVB': 'uvb',
    'UMIN': 'umin', 'GAMMA': 'gamma', 'QPAH': 'qpah',
}


class BayesianGrid(object):
    """Pre-synthesized Bayesian grid photometry, shared read-only across MPPool workers.

    Populated by :func:`_initialize_fastbayes_workers` once per worker
    process, mirroring the ``sc_data`` singleton pattern in
    ``fastspecfit.singlecopy``.

    """
    def __init__(self):
        self.file = None


    def load(self, gridfile):
        """Load a ``bayesian-photometry-*.fits`` grid file (idempotent)."""
        if self.file == gridfile:
            return

        T = fitsio.FITS(gridfile)
        prihdr = T[0].read_header()

        self.gridnumber = prihdr.get('GRIDNUM')
        self.imf = prihdr.get('IMF')
        self.fphotofile = prihdr.get('FPHOTO')
        self.logspace = bool(prihdr.get('LOGSPACE'))

        photsys_hdr = prihdr.get('PHOTSYS')
        self.photsys_keys = [''] if photsys_hdr in (None, 'NONE') else photsys_hdr.split(',')

        self.redshift = T['REDSHIFT'].read().astype('f8')
        self.zgrid_interp = np.log10(1. + self.redshift) if self.logspace else self.redshift

        self.meta = Table(T['METADATA'].read())
        self.ntemplate = len(self.meta)

        filt = T['FILTERS'].read()
        self.bands = np.array([b.strip() for b in filt['band']])
        self.nband = len(self.bands)

        self.maggies = {}
        for key in self.photsys_keys:
            extname = 'MAGGIES' + ('' if key == '' else f'_{key}')
            self.maggies[key] = T[extname].read() # [ntemplate, nz, nband]

        # File is set last so a failed/partial load is retried rather than
        # treated as already-cached.
        self.file = gridfile
        log.debug(f'Cached Bayesian grid {self.file}')


    def interpolate_at_z(self, photsys, redshift):
        """Linearly interpolate ``MAGGIES(template, z, band)`` to a single redshift.

        Parameters
        ----------
        photsys : :class:`str`
            Photometric-system key (e.g. ``'N'``, ``'S'``, or ``''``).
        redshift : :class:`float`
            Object redshift.

        Returns
        -------
        :class:`numpy.ndarray`
            Model maggies per solar mass formed, shape (ntemplate, nband).

        """
        if photsys not in self.maggies:
            errmsg = (f"photsys '{photsys}' not present in grid file {self.file} "
                      f"(available: {list(self.maggies.keys())}).")
            log.critical(errmsg)
            raise KeyError(errmsg)

        maggies = self.maggies[photsys] # [ntemplate, nz, nband]

        zvar = np.log10(1. + redshift) if self.logspace else redshift
        zvar = np.clip(zvar, self.zgrid_interp[0], self.zgrid_interp[-1])

        idx = int(np.clip(np.searchsorted(self.zgrid_interp, zvar), 1, len(self.zgrid_interp) - 1))
        z0, z1 = self.zgrid_interp[idx - 1], self.zgrid_interp[idx]
        frac = 0. if z1 == z0 else (zvar - z0) / (z1 - z0)

        return maggies[:, idx - 1, :] + frac * (maggies[:, idx, :] - maggies[:, idx - 1, :])


# global structure with single-copy grid data, initially empty
bg_data = BayesianGrid()


def _initialize_fastbayes_workers(fphotofile=None, gridfile=None):
    """MPPool initializer: populate ``sc_data.photometry`` and ``bg_data`` in each worker.

    ``sc_data.photometry`` is populated directly (rather than via
    ``sc_data.initialize()``) so this mode never loads the stellar template
    basis or emission-line tables, neither of which it needs.

    """
    sc_data.photometry = Photometry(fphotofile=fphotofile)
    bg_data.load(gridfile)

    if list(bg_data.bands) != list(sc_data.photometry.bands):
        errmsg = ('Band mismatch between the photometry configuration and the Bayesian grid file; '
                 'they must be built from the same (or a compatible) fphoto configuration file.')
        log.critical(errmsg)
        raise ValueError(errmsg)


def _weighted_percentile(values, weights, percentiles):
    """Weighted percentiles of ``values`` via linear interpolation of the weighted CDF."""
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w)
    cw /= cw[-1]
    return np.interp(np.asarray(percentiles) / 100., cw, v)


def _weighted_mode(values, weights):
    """Weighted mode: the (exact) grid value with the largest total weight."""
    uniq, inv = np.unique(values, return_inverse=True)
    binweights = np.bincount(inv, weights=weights)
    return uniq[np.argmax(binweights)]


def get_fastbayes_dtype(topk=0):
    """Build the output dtype for the per-object Bayesian-fitting results.

    Parameters
    ----------
    topk : :class:`int`, optional
        Number of top-weight grid templates to store per object (sparse
        joint posterior). Disabled when ``0`` (default).

    Returns
    -------
    :class:`numpy.dtype`

    """
    cols = []
    for pname in PARAM_NAMES:
        cols += [(f'{pname}_BEST', 'f4'), (f'{pname}_MEAN', 'f4'), (f'{pname}_MODE', 'f4'),
                (f'{pname}_P16', 'f4'), (f'{pname}_P50', 'f4'), (f'{pname}_P84', 'f4')]
    cols += [('CHI2MIN', 'f4'), ('NBANDS', 'i2')]
    if topk > 0:
        cols += [('TOPK_INDEX', 'i4', (topk,)), ('TOPK_WEIGHT', 'f4', (topk,))]
    return np.dtype(cols)


def fastbayes_one(iobj, data, meta, fastbayes_dtype, topk=0, uncertainty_floor=0.01):
    """Fit one object's broadband photometry against the Bayesian grid.

    Parameters
    ----------
    iobj : :class:`int`
        Index of the object in the input list, used for log messages.
    data : :class:`dict`
        Per-object data dictionary from
        :meth:`fastspecfit.io.DESISpectra.read` (``fastphot=True``).
    meta : :class:`astropy.table.Row`
        Metadata row for this object.
    fastbayes_dtype : :class:`numpy.dtype`
        Output dtype, from :func:`get_fastbayes_dtype`.
    topk : :class:`int`, optional
        Number of top-weight grid templates to store (0 disables).
    uncertainty_floor : :class:`float`, optional
        Minimum fractional photometric uncertainty added in quadrature.

    Returns
    -------
    meta : :class:`astropy.table.Row`
        Updated metadata row with observed photometry filled in.
    result : :class:`numpy.ndarray`
        Bayesian-fitting output row.

    """
    from fastspecfit.io import one_spectrum

    phot = sc_data.photometry

    log.info(f'Bayesian fitting object {iobj} [{phot.uniqueid_col.lower()} '
            f'{data["uniqueid"]}, z={data["redshift"]:.6f}].')

    one_spectrum(data, meta, uncertainty_floor=uncertainty_floor, fastphot=True, synthphot=True)

    # Copy parsed photometry from the 'data' dictionary to the 'meta' table
    # (mirroring fastspec_one's convention).
    flux = data['photometry']['nanomaggies']
    fluxivar = data['photometry']['nanomaggies_ivar']
    for iband, band in enumerate(phot.bands):
        meta[f'FLUX_{band.upper()}'] = flux[iband]
        meta[f'FLUX_IVAR_{band.upper()}'] = fluxivar[iband]

    flam = np.asarray(data['photometry']['flam'], dtype='f8')
    flam_ivar = np.asarray(data['photometry']['flam_ivar'], dtype='f8') * phot.bands_to_fit
    lambda_eff = np.asarray(data['photometry']['lambda_eff'], dtype='f8')

    model_maggies = bg_data.interpolate_at_z(data['photsys'], data['redshift']) # [ntemplate, nband]
    model_flam = Photometry.get_photflam(model_maggies, lambda_eff) # [ntemplate, nband]

    # Closed-form, non-negative, chi2-minimizing mass amplitude per template.
    numer = (model_flam * (flam_ivar * flam)[np.newaxis, :]).sum(axis=1)
    denom = (model_flam**2 * flam_ivar[np.newaxis, :]).sum(axis=1)
    amplitude = np.divide(numer, denom, out=np.zeros_like(numer), where=denom > 0.)
    amplitude = np.clip(amplitude, 0., None)

    resid = flam[np.newaxis, :] - amplitude[:, np.newaxis] * model_flam
    chi2 = (flam_ivar[np.newaxis, :] * resid**2).sum(axis=1)

    chi2min = np.min(chi2)
    weight = np.exp(-0.5 * (chi2 - chi2min))
    weight /= weight.sum()

    ibest = np.argmin(chi2)
    logmass = np.log10(np.clip(amplitude, 1e-30, None))

    result = np.zeros(1, dtype=fastbayes_dtype)[0]
    for pname in PARAM_NAMES:
        vals = logmass if pname == 'LOGMASS' else np.asarray(bg_data.meta[_GRID_COLUMN[pname]])
        result[f'{pname}_BEST'] = vals[ibest]
        result[f'{pname}_MEAN'] = np.sum(weight * vals)
        result[f'{pname}_MODE'] = _weighted_mode(vals, weight)
        p16, p50, p84 = _weighted_percentile(vals, weight, (16., 50., 84.))
        result[f'{pname}_P16'] = p16
        result[f'{pname}_P50'] = p50
        result[f'{pname}_P84'] = p84

    result['CHI2MIN'] = chi2min
    result['NBANDS'] = int(np.sum(flam_ivar > 0.))

    if topk > 0:
        idx = np.argsort(weight)[::-1][:topk]
        result['TOPK_INDEX'][:len(idx)] = idx.astype('i4')
        result['TOPK_WEIGHT'][:len(idx)] = weight[idx].astype('f4')

    return meta, result


def write_fastbayes(meta, results, outfile, gridfile, fphotofile, topk=0):
    """Write the Bayesian-fitting output to a multi-extension FITS file.

    Parameters
    ----------
    meta : :class:`astropy.table.Table`
        Output metadata table (see :func:`fastspecfit.io.create_output_meta`).
    results : :class:`astropy.table.Table`
        Output fitting-results table.
    outfile : :class:`str`
        Full path of the output FITS file.
    gridfile : :class:`str`
        Full path of the Bayesian grid file used for fitting.
    fphotofile : :class:`str`
        Full path of the photometry configuration file used.
    topk : :class:`int`, optional
        Number of top-weight grid templates stored per object, if any.

    """
    from astropy.io import fits

    hduprim = fits.PrimaryHDU()
    hduprim.header['GRIDFILE'] = os.path.abspath(str(gridfile))
    hduprim.header['FPHOTO'] = os.path.abspath(str(fphotofile)) if fphotofile else ''
    hduprim.header['TOPK'] = topk

    hdumeta = fits.convenience.table_to_hdu(meta)
    hdumeta.header['EXTNAME'] = 'METADATA'

    hduresults = fits.convenience.table_to_hdu(results)
    hduresults.header['EXTNAME'] = 'FASTBAYES'

    hx = fits.HDUList([hduprim, hdumeta, hduresults])

    outdir = os.path.dirname(outfile)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    log.info(f'Writing results for {len(meta)} objects to {outfile}')
    hx.writeto(outfile, overwrite=True)


def parse(options=None):
    """Parse input arguments to the ``fastbayes`` script.

    Parameters
    ----------
    options : list of str or None, optional
        Command-line argument strings. If ``None``, reads from ``sys.argv``.

    Returns
    -------
    args : :class:`argparse.Namespace`

    """
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('redrockfiles', nargs='+', help='Full path to input redrock file(s).')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename (required).')
    parser.add_argument('--gridfile', type=str, required=True,
                        help='Full path to the Bayesian photometry grid file '
                        '(output of bin/build-bayesian-photometry).')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to process in each file, zero-indexed.')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--input-redshifts', type=str, default=None,
                        help='Comma-separated list of input redshifts corresponding to the (required) --targetids input.')
    parser.add_argument('--zmin', type=float, default=None, help='Override the default minimum redshift required for modeling.')
    parser.add_argument('--topk', type=int, default=0,
                        help='Number of top-weight grid templates to save per object (sparse joint posterior); 0 disables.')
    parser.add_argument('--use-quasarnet', dest='use_quasarnet', default=False, action='store_true',
                        help='Use QuasarNet to improve QSO redshifts.')
    parser.add_argument('--fphotodir', type=str, default=None, help='Top-level location of the source photometry.')
    parser.add_argument('--fphotofile', type=str, default=None,
                        help='Photometric configuration file (default: taken from the grid file header).')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--redux_dir', type=str, default=None, help='Optional full path $DESI_SPECTRO_REDUX.')
    parser.add_argument('--specproddir', type=str, default=None, help='Optional directory name for the spectroscopic production.')
    parser.add_argument('--uncertainty-floor', type=float, default=0.01,
                        help='Minimum fractional uncertainty to add in quadrature to the formal inverse variance.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    log.info('fastbayes {}'.format(' '.join(options)))

    return parser.parse_args(options)


def fastbayes(args=None, mp_pool=None):
    """Main fastbayes engine: read, fit, and write Bayesian grid photometric results.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or list of str or None, optional
        Pre-parsed arguments or raw argument list. If ``None``, reads from
        ``sys.argv``.
    mp_pool : :class:`fastspecfit.util.MPPool` or None, optional
        Pre-built worker pool; a new one is created (and closed) when
        ``None``.

    Returns
    -------
    :class:`int`
        Exit code (0 on success).

    """
    from astropy.table import vstack
    from fastspecfit.io import DESISpectra, create_output_meta, create_output_table

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    envlist = []
    if args.redux_dir is None:
        envlist.append('DESI_SPECTRO_REDUX')
    if args.mapdir is None:
        envlist.append('DUST_DIR')
    if args.fphotodir is None:
        envlist.append('FPHOTO_DIR')
    for env in envlist:
        if env not in os.environ:
            errmsg = f'Mandatory environment variable {env} missing.'
            log.critical(errmsg)
            raise KeyError(errmsg)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    targetids = None
    input_redshifts = None
    if args.targetids is not None:
        targetids = [int(x) for x in args.targetids.split(',')]
        if args.input_redshifts is not None:
            input_redshifts = [float(x) for x in args.input_redshifts.split(',')]
            if len(input_redshifts) != len(targetids):
                errmsg = 'targetids and input_redshifts must have the same number of elements.'
                log.critical(errmsg)
                raise ValueError(errmsg)

    gridfile = os.path.expandvars(args.gridfile)
    if not os.path.isfile(gridfile):
        errmsg = f'Bayesian grid file {gridfile} not found.'
        log.critical(errmsg)
        raise IOError(errmsg)

    fphotofile = args.fphotofile
    if fphotofile is None:
        with fitsio.FITS(gridfile) as F:
            fphotofile = F[0].read_header().get('FPHOTO')

    init_argdict = {'fphotofile': fphotofile, 'gridfile': gridfile}

    t0 = time.time()
    _initialize_fastbayes_workers(**init_argdict)
    log.info(fsftime('bg_data_init', time.time() - t0))

    _own_pool = False
    if mp_pool is None:
        mp_pool = MPPool(args.mp, initializer=_initialize_fastbayes_workers, init_argdict=init_argdict)
        _own_pool = True

    phot = sc_data.photometry
    cosmo = TabulatedDESI()

    Spec = DESISpectra(phot=phot, cosmo=cosmo, fphotodir=args.fphotodir,
                       mapdir=args.mapdir, redux_dir=args.redux_dir)

    Spec.gather_metadata(args.redrockfiles, firsttarget=args.firsttarget,
                         targetids=targetids, input_redshifts=input_redshifts,
                         ntargets=args.ntargets, zmin=args.zmin,
                         redrockfile_prefix=args.redrockfile_prefix,
                         use_quasarnet=args.use_quasarnet, specprod_dir=args.specproddir)
    if len(Spec.specfiles) == 0:
        return 0

    data, meta = Spec.read(phot, fastphot=True)

    nobj = len(meta)
    fastbayes_dtype = get_fastbayes_dtype(topk=args.topk)

    fitargs = [{
        'iobj': iobj,
        'data': data[iobj],
        'meta': meta[iobj],
        'fastbayes_dtype': fastbayes_dtype,
        'topk': args.topk,
        'uncertainty_floor': args.uncertainty_floor,
    } for iobj in range(nobj)]

    t0 = time.time()
    out = mp_pool.starmap(fastbayes_one, fitargs)
    out = list(zip(*out))

    outmeta = create_output_meta(vstack(out[0]), phot=phot, fastphot=True)

    units = {f'AGE_{stat}': 'Gyr' for stat in ('BEST', 'MEAN', 'MODE', 'P16', 'P50', 'P84')}
    results = create_output_table(out[1], outmeta, units)

    if _own_pool:
        mp_pool.close()

    log.info(fsftime('fastbayes_all', time.time() - t0, context=f'nobj={nobj}'))

    write_fastbayes(outmeta, results, outfile=args.outfile, gridfile=gridfile,
                  fphotofile=phot.fphotofile, topk=args.topk)

    return 0
