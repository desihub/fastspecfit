"""
fastspecfit.fastspecfit
=======================

See sandbox/running-fastspecfit for examples.

"""
import os, sys, time
import numpy as np

from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.singlecopy import sc_data, _initialize_sc_data
from fastspecfit.util import BoxedScalar, MPPool, NMONTE_DEFAULT, fsftime
from fastspecfit.templates import VDISP_NOMINAL, VDISP_BOUNDS

def make_init_sc_args(args, fastphot=False, fitstack=False):
    """Build the sc_data.initialize() kwargs from a parsed args Namespace.

    Parameters
    ----------
    args : :class:`argparse.Namespace`
        Parsed command-line arguments (from :func:`parse` or from
        ``mpi-fastspecfit``).
    fastphot : bool, optional
        Passed through to :meth:`~fastspecfit.singlecopy.Singletons.initialize`.
    fitstack : bool, optional
        Passed through to :meth:`~fastspecfit.singlecopy.Singletons.initialize`.

    Returns
    -------
    dict
        Keyword arguments suitable for ``sc_data.initialize(**...)``.

    """
    return {
        'emlines_file':      getattr(args, 'emlinesfile', None),
        'fphotofile':        getattr(args, 'fphotofile', None),
        'fastphot':          fastphot,
        'fitstack':          fitstack,
        'ignore_photometry': getattr(args, 'ignore_photometry', False),
        'template_file':     getattr(args, 'templates', None),
        'template_version':  getattr(args, 'templateversion', None),
        'template_imf':      getattr(args, 'imf', None),
        'log_verbose':       getattr(args, 'verbose', False),
        'vdisp_nominal':     getattr(args, 'vdisp_nominal', VDISP_NOMINAL),
        'vdisp_bounds':      getattr(args, 'vdisp_bounds', VDISP_BOUNDS),
    }


def parse(options=None, rank=0):
    """Parse input arguments to fastspec and fastphot scripts.

    Parameters
    ----------
    options : list of str or None, optional
        Command-line argument strings. If ``None``, reads from ``sys.argv``.
    rank : int, optional
        MPI rank of the calling process, used to suppress log output on
        non-zero ranks. Defaults to 0.

    Returns
    -------
    args : :class:`argparse.Namespace`
        Parsed command-line arguments.

    """
    import argparse, sys
    from fastspecfit.templates import Templates

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('redrockfiles', nargs='+', help='Full path to input redrock file(s).')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename (required).')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads per MPI rank.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file, zero-indexed.')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--input-redshifts', type=str, default=None, help='Comma-separated list of input redshifts corresponding to the (required) --targetids input.')
    parser.add_argument('--input-seeds', type=str, default=None, help='Comma-separated list of input random-number seeds corresponding to the (required) --targetids input.')
    parser.add_argument('--seed', type=int, default=1, help='Random seed for Monte Carlo reproducibility; ignored if --input-seeds is passed.')
    parser.add_argument('--nmonte', type=int, default=NMONTE_DEFAULT, help='Number of Monte Carlo realizations.')
    parser.add_argument('--vdisp-nominal', type=float, default=VDISP_NOMINAL, help='Nominal (default) velocity dispersion in km/s.')
    parser.add_argument('--vdisp-bounds', type=float, default=VDISP_BOUNDS, nargs=2, help='Nominal (default) velocity dispersion in km/s.')
    parser.add_argument('--zmin', type=float, default=None, help='Override the default minimum redshift required for modeling.')
    parser.add_argument('--no-broadlinefit', default=True, action='store_false', dest='broadlinefit',
                        help='Do not model broad Balmer and helium line-emission.')
    parser.add_argument('--ignore-photometry', default=False, action='store_true', help='Ignore the broadband photometry during model fitting.')
    parser.add_argument('--ignore-quasarnet', dest='use_quasarnet', default=True, action='store_false', help='Do not use QuasarNet to improve QSO redshifts.')
    parser.add_argument('--constrain-age', action='store_true', help='Constrain the age of the SED.')
    parser.add_argument('--no-smooth-continuum', action='store_true', help='Do not fit the smooth continuum.')
    parser.add_argument('--imf', type=str, default=Templates.DEFAULT_IMF, help='Initial mass function.')
    parser.add_argument('--templateversion', type=str, default=Templates.DEFAULT_TEMPLATEVERSION, help='Template version number.')
    parser.add_argument('--templates', type=str, default=None, help='Optional full path and filename to the templates.')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--specfile-prefix', type=str, default='coadd-', help='Prefix of the spectral file(s).')
    parser.add_argument('--qnfile-prefix', type=str, default='qso_qn-', help='Prefix of the QuasarNet afterburner file(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--fphotodir', type=str, default=None, help='Top-level location of the source photometry.')
    parser.add_argument('--fphotofile', type=str, default=None, help='Photometric information file.')
    parser.add_argument('--emlinesfile', type=str, default=None, help='Emission line parameter file.')
    parser.add_argument('--redux_dir', type=str, default=None, help='Optional full path $DESI_SPECTRO_REDUX.')
    parser.add_argument('--specproddir', type=str, default=None, help='Optional directory name for the spectroscopic production.')
    parser.add_argument('--uncertainty-floor', type=float, default=0.01, help='Minimum fractional uncertainty to add in quadrature to the formal inverse variance spectrum.')
    parser.add_argument('--minsnr-balmer-broad', type=float, default=2.5, help='Minimum broad Balmer S/N to force broad+narrow-line model.')
    parser.add_argument('--debug-plots', action='store_true', help='Generate a variety of debugging plots (written to $PWD).')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    if rank == 0:
        log.info(f'fastspec {" ".join(options)}')

    return parser.parse_args(options)


def fastspec_one(iobj, data, meta, fastfit_dtype, specphot_dtype, broadlinefit=True,
                 fastphot=False, fitstack=False, constrain_age=False,
                 no_smooth_continuum=False, debug_plots=False, uncertainty_floor=0.01,
                 minsnr_balmer_broad=2.5, nmonte=NMONTE_DEFAULT, seed=1):
    """Fit the continuum and emission lines for a single DESI object.

    Parameters
    ----------
    iobj : int
        Index of the object in the input list, used for log messages.
    data : dict
        Per-object data dictionary from :class:`fastspecfit.io.DESISpectra`.
    meta : :class:`astropy.table.Row`
        Metadata row for this object.
    fastfit_dtype : :class:`numpy.dtype`
        NumPy dtype for the emission-line fitting output table.
    specphot_dtype : :class:`numpy.dtype`
        NumPy dtype for the spectrophotometric output table.
    broadlinefit : bool, optional
        Attempt to fit broad Balmer components. Defaults to ``True``.
    fastphot : bool, optional
        Fit broadband photometry only (no spectra). Defaults to ``False``.
    fitstack : bool, optional
        Treat input as stacked spectra. Defaults to ``False``.
    constrain_age : bool, optional
        Constrain the stellar population age during fitting. Defaults to
        ``False``.
    no_smooth_continuum : bool, optional
        Skip smooth continuum fitting. Defaults to ``False``.
    debug_plots : bool, optional
        Write diagnostic QA plots to the working directory. Defaults to
        ``False``.
    uncertainty_floor : float, optional
        Minimum fractional uncertainty added in quadrature to the formal
        inverse variance spectrum. Defaults to 0.01.
    minsnr_balmer_broad : float, optional
        Minimum S/N required to adopt the broad Balmer component. Defaults
        to 2.5.
    nmonte : int, optional
        Number of Monte Carlo realizations for uncertainty estimation.
    seed : int, optional
        Random seed for Monte Carlo reproducibility. Defaults to 1.

    Returns
    -------
    meta : :class:`astropy.table.Row`
        Updated metadata row with observed photometry filled in.
    specphot : :class:`numpy.ndarray`
        Spectrophotometric output row.
    fastfit : :class:`numpy.ndarray`
        Emission-line fitting output row.
    emmodel : :class:`astropy.table.Table` or None
        Coadded model spectra table, or ``None`` if ``fastphot=True``.

    """
    from fastspecfit.io import one_spectrum, one_stacked_spectrum
    from fastspecfit.emlines import emline_specfit
    from fastspecfit.continuum import continuum_specfit

    igm = sc_data.igm
    phot = sc_data.photometry
    emline_table = sc_data.emlines.table
    templates = sc_data.templates

    t0 = time.time()

    if fastphot:
        log.info(f'Continuum fitting object {iobj} [{phot.uniqueid_col.lower()} ' + \
                 f'{data["uniqueid"]}, seed {seed}, z={data["redshift"]:.6f}].')
    else:
        log.info(f'Continuum- and emission-line fitting object {iobj} [{phot.uniqueid_col.lower()} ' + \
                 f'{data["uniqueid"]}, seed {seed}, z={data["redshift"]:.6f}].')

    if fitstack:
        one_stacked_spectrum(data, meta, synthphot=True, debug_plots=debug_plots)
    else:
        one_spectrum(data, meta, uncertainty_floor=uncertainty_floor, fastphot=fastphot,
                     synthphot=True, debug_plots=debug_plots)

    # Copy parsed photometry from the 'data' dictionary to the 'meta' table.
    if not fitstack:
        flux = data['photometry']['nanomaggies']
        fluxivar = data['photometry']['nanomaggies_ivar']
        for iband, band in enumerate(phot.bands):
            meta[f'FLUX_{band.upper()}'] = flux[iband]
            meta[f'FLUX_IVAR_{band.upper()}'] = fluxivar[iband]

        if hasattr(phot, 'fiber_bands'):
            fiberflux = data['fiberphot']['nanomaggies']
            fibertotflux = data['fibertotphot']['nanomaggies']
            for iband, band in enumerate(phot.fiber_bands):
                meta[f'FIBERFLUX_{band.upper()}'] = fiberflux[iband]
                meta[f'FIBERTOTFLUX_{band.upper()}'] = fibertotflux[iband]

    # output structures
    fastfit = BoxedScalar(fastfit_dtype)
    specphot = BoxedScalar(specphot_dtype)

    continuummodel, smooth_continuum, continuummodel_monte, specflux_monte = \
        continuum_specfit(data, fastfit, specphot, templates, igm, phot, constrain_age=constrain_age,
                          no_smooth_continuum=no_smooth_continuum, fastphot=fastphot,
                          fitstack=fitstack, debug_plots=debug_plots, nmonte=nmonte,
                          seed=seed)

    # Optionally fit the emission-line spectrum.
    if fastphot:
        emmodel = None
    else:
        emmodel = emline_specfit(data, fastfit, specphot, continuummodel, smooth_continuum,
                                 phot, emline_table, broadlinefit=broadlinefit,
                                 minsnr_balmer_broad=minsnr_balmer_broad,
                                 debug_plots=debug_plots, specflux_monte=specflux_monte,
                                 continuummodel_monte=continuummodel_monte)

    log.info(fsftime('fastspec_one', time.time()-t0,
                     context=f'{phot.uniqueid_col.lower()}={data["uniqueid"]}'))

    return meta, specphot.value, fastfit.value, emmodel


def fastspec(fastphot=False, fitstack=False, args=None, comm=None, verbose=False, mp_pool=None):
    """Main fastspec engine: read, fit, and write results for one or more DESI spectra.

    Parameters
    ----------
    fastphot : bool, optional
        Fit broadband photometry only (no spectra). Defaults to ``False``.
    fitstack : bool, optional
        Treat input as stacked spectra. Defaults to ``False``.
    args : :class:`argparse.Namespace` or list of str or None, optional
        Pre-parsed arguments or raw argument list. If ``None``, reads from
        ``sys.argv``.
    comm : :class:`mpi4py.MPI.Comm` or None, optional
        MPI intracommunicator for parallel execution, or ``None`` for
        single-process mode.
    verbose : bool, optional
        Enable verbose (debug-level) logging. Defaults to ``False``.

    Returns
    -------
    int
        Exit code (0 on success).

    """
    from astropy.table import vstack
    from fastspecfit.io import (DESISpectra, get_output_dtype, create_output_meta,
                                create_output_table, write_fastspecfit)

    if comm:
        rank, size = comm.rank, comm.size
    else:
        rank, size = 0, 1

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args, rank=rank)

    if fitstack:
        args.ignore_photometry = True

    # check for mandatory environment variables
    envlist = []
    if not fitstack and args.redux_dir is None:
        envlist += ['DESI_SPECTRO_REDUX']
    if not fitstack and args.mapdir is None:
        envlist += ['DUST_DIR']
    if not fitstack and args.fphotodir is None:
        envlist += ['FPHOTO_DIR']
    if args.templates is None:
        envlist += ['FTEMPLATES_DIR']
    for env in envlist:
        if not env in os.environ:
            errmsg = f'Mandatory environment variable {env} missing.'
            log.critical(errmsg)
            raise KeyError(errmsg)

    if verbose:
        args.verbose = True

    targetids = None
    input_redshifts = None
    input_seeds = None

    if args.targetids is not None:
        targetids = [int(x) for x in args.targetids.split(',')]
        if args.input_redshifts is not None:
            input_redshifts = [float(x) for x in args.input_redshifts.split(',')]
            if len(input_redshifts) != len(targetids):
                errmsg = 'targetids and input_redshifts must have the same number of elements.'
                log.critical(errmsg)
                raise ValueError(errmsg)
        if args.input_seeds is not None:
            input_seeds = [np.int64(x) for x in args.input_seeds.split(',')]
            if len(input_seeds) != len(targetids):
                errmsg = 'targetids and input_seeds must have the same number of elements.'
                log.critical(errmsg)
                raise ValueError(errmsg)

    init_sc_args = make_init_sc_args(args, fastphot=fastphot, fitstack=fitstack)
    t0 = time.time()
    sc_data.initialize(**init_sc_args)
    log.info(fsftime('sc_data_init', time.time()-t0))

    _own_pool = False
    if rank == 0:
        t0 = time.time()
        if comm is None and mp_pool is None:
            mp_pool = MPPool(args.mp, initializer=_initialize_sc_data,
                             init_argdict=init_sc_args)
            _own_pool = True

        log.debug(fsftime('init_workers', time.time()-t0,
                          context=f'nworkers={args.mp}'))

        log.info(f'Cached stellar templates {sc_data.templates.file}')
        log.info(f'Cached emission-line table {sc_data.emlines.file}')
        log.info(f'Cached photometric filters and parameters {sc_data.photometry.fphotofile}')
        log.info(f'Cached cosmology table {sc_data.cosmology.file}')
        log.info(f'Cached {sc_data.igm.reference} IGM attenuation parameters.')

        # Read the data.
        Spec = DESISpectra(phot=sc_data.photometry, cosmo=sc_data.cosmology,
                           fphotodir=args.fphotodir, mapdir=args.mapdir,
                           redux_dir=args.redux_dir)

        if fitstack:
            data, meta = Spec.read_stacked(args.redrockfiles, firsttarget=args.firsttarget,
                                           stackids=targetids, ntargets=args.ntargets,
                                           constrain_age=args.constrain_age)
        else:
            Spec.gather_metadata(args.redrockfiles, firsttarget=args.firsttarget,
                                 targetids=targetids, input_redshifts=input_redshifts,
                                 ntargets=args.ntargets, zmin=args.zmin,
                                 redrockfile_prefix=args.redrockfile_prefix,
                                 specfile_prefix=args.specfile_prefix, qnfile_prefix=args.qnfile_prefix,
                                 use_quasarnet=args.use_quasarnet, specprod_dir=args.specproddir)
            if len(Spec.specfiles) == 0:
                return 0

            data, meta = Spec.read(sc_data.photometry, fastphot=fastphot, constrain_age=args.constrain_age)

        ntargets = len(meta)
        ncoeff = sc_data.templates.ntemplates
        if fastphot:
            cameras = None
        else:
            cameras = data[0]['cameras']

        fastfit_dtype, fastfit_units = get_output_dtype(
            Spec.specprod, phot=sc_data.photometry,
            linetable=sc_data.emlines.table, ncoeff=ncoeff,
            cameras=cameras, fastphot=fastphot,
            fitstack=fitstack)

        specphot_dtype, specphot_units = get_output_dtype(
            Spec.specprod, phot=sc_data.photometry,
            linetable=sc_data.emlines.table, ncoeff=ncoeff,
            cameras=cameras, fastphot=fastphot,
            fitstack=fitstack, specphot=True)

        # If using Monte Carlo, generate the random seed(s).
        if args.nmonte > 0:
            if input_seeds is not None:
                seeds = input_seeds
            else:
                rng = np.random.default_rng(seed=args.seed)
                seeds = rng.integers(2**32, size=ntargets, dtype=np.int64)
        else:
            seeds = [1] * ntargets

        nobj = len(meta)
        if args.outfile:
            _outbase = os.path.basename(args.outfile).replace('.fits.gz', '').replace('.fits', '')
            for d in data:
                d['outfile_base'] = _outbase

        fitargs = [{
            'iobj':                iobj,
            'data':                data[iobj],
            'meta':                meta[iobj],
            'fastfit_dtype':       fastfit_dtype,
            'specphot_dtype':      specphot_dtype,
            'broadlinefit':        args.broadlinefit,
            'fastphot':            fastphot,
            'fitstack':            fitstack,
            'constrain_age':       args.constrain_age,
            'no_smooth_continuum': args.no_smooth_continuum,
            'debug_plots':         args.debug_plots,
            'uncertainty_floor':   args.uncertainty_floor,
            'minsnr_balmer_broad': args.minsnr_balmer_broad,
            'nmonte':              args.nmonte,
            'seed':                seeds[iobj],
        } for iobj in range(nobj)]


    # Fit in parallel
    t0 = time.time()
    if comm is not None:
        # Rank=0 sends work to all the other ranks (and also does work itself).
        if rank == 0:
            log.info(f'Rank {rank}: distributing {nobj:,d} objects to {comm.size:,d} ranks.')
            fitargs_byrank = np.array_split(fitargs, size)
            for onerank in range(1, size):
                #log.debug(f'Rank 0 sending data on {len(fitargs_byrank[onerank])}/{len(meta)} objects to rank {onerank}')
                comm.send(fitargs_byrank[onerank], dest=onerank)
            fitargs_onerank = fitargs_byrank[rank]
        else:
            fitargs_onerank = comm.recv(source=0)
        #log.debug(f'Rank {rank} received data on {len(fitargs_onerank)} objects from rank 0')

        # Each rank, including rank 0, iterates over each object and then sends
        # the results to rank 0.
        log.info(f'Rank {rank}: fitting {len(fitargs_onerank):,d} objects.')
        t1 = time.time()
        out = []
        for fitarg_onerank in fitargs_onerank:
            out.append(fastspec_one(**fitarg_onerank))
        log.info(fsftime('fit_rank', time.time()-t1,
                         context=f'rank={rank}, nobj={len(out):,d}'))

        if rank > 0:
            #log.debug(f'Rank {rank} sending data on {len(out)} objects to rank 0.')
            comm.send(out, dest=0)
        else:
            for onerank in range(1, size):
                out.extend(comm.recv(source=onerank))
            #log.debug(f'Rank 0 received data on {len(out)} objects.')
            log.info(f'Rank {rank}: collected fitting results for {len(out):,d} ' + \
                     f'objects from {comm.size:,d} ranks.')
    else:
        out = mp_pool.starmap(fastspec_one, fitargs)

    out = list(zip(*out))

    if rank == 0:
        meta = create_output_meta(vstack(out[0]), phot=sc_data.photometry,
                                  fastphot=fastphot, fitstack=fitstack)
        specphot = create_output_table(out[1], meta, specphot_units, fitstack=fitstack)

        if fastphot:
            fastfit = None
            modelspectra = None
        else:
            fastfit = create_output_table(out[2], meta, fastfit_units, fitstack=fitstack)
            modelspectra = vstack(out[3], join_type='exact', metadata_conflicts='error')

        if comm is None and _own_pool:
            mp_pool.close()

        _elapsed = time.time() - t0
        _ncore = max(args.mp, 1)
        _per_obj = _elapsed / _ncore / ntargets
        log.info(fsftime('fit_all', _elapsed,
                         context=f'nobj={ntargets},{_per_obj:.2f}s/obj/core'))

        write_fastspecfit(
            meta, specphot, fastfit, modelspectra=modelspectra,
            outfile=args.outfile, specprod=Spec.specprod, coadd_type=Spec.coadd_type,
            fphotofile=sc_data.photometry.fphotofile,
            template_file=sc_data.templates.file,
            emlinesfile=sc_data.emlines.file, fastphot=fastphot,
            inputz=input_redshifts is not None,
            nmonte=args.nmonte, vdisp_nominal=args.vdisp_nominal,
            vdisp_bounds=args.vdisp_bounds,
            seed=args.seed, inputseeds=input_seeds is not None,
            uncertainty_floor=args.uncertainty_floor,
            minsnr_balmer_broad=args.minsnr_balmer_broad,
            ignore_photometry=args.ignore_photometry,
            broadlinefit=args.broadlinefit, constrain_age=args.constrain_age,
            use_quasarnet=args.use_quasarnet,
            no_smooth_continuum=args.no_smooth_continuum)

        return 0


def fastphot(args=None, comm=None, verbose=False, mp_pool=None):
    """Main fastphot entry point: fit broadband photometry for DESI objects.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or list of str or None, optional
        Pre-parsed arguments or raw argument list. If ``None``, reads from
        ``sys.argv``.
    comm : :class:`mpi4py.MPI.Comm` or None, optional
        MPI intracommunicator for parallel execution, or ``None`` for
        single-process mode.
    verbose : bool, optional
        Enable verbose (debug-level) logging. Defaults to ``False``.
    mp_pool : :class:`fastspecfit.util.MPPool` or None, optional
        Pre-created worker pool to reuse; a new pool is created and closed
        when ``None``.

    Returns
    -------
    int
        Exit code (0 on success).

    """
    return fastspec(fastphot=True, args=args, comm=comm, verbose=verbose, mp_pool=mp_pool)


def stackfit(args=None, comm=None, verbose=False):
    """Entry point to fit generic stacked DESI spectra.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or list of str or None, optional
        Pre-parsed arguments or raw argument list. If ``None``, reads from
        ``sys.argv``.
    comm : :class:`mpi4py.MPI.Comm` or None, optional
        MPI intracommunicator for parallel execution, or ``None`` for
        single-process mode.
    verbose : bool, optional
        Enable verbose (debug-level) logging. Defaults to ``False``.

    Returns
    -------
    int
        Exit code (0 on success).

    """
    return fastspec(fitstack=True, args=args, comm=comm, verbose=verbose, mp_pool=None)
