"""
fastspecfit.fastspecfit
=======================

See sandbox/running-fastspecfit for examples.

"""
import os, time
import numpy as np

from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.singlecopy import sc_data
from fastspecfit.util import BoxedScalar, MPPool


def fastspec_one(iobj, data, meta, out_dtype, broadlinefit=True, fastphot=False,
                 fitstack=False, constrain_age=False, no_smooth_continuum=False,
                 debug_plots=False, uncertainty_floor=0.01, minsnr_balmer_broad=2.5,
                 nmonte=50, seed=1):
    """Run :func:`fastspec` on a single object.

    """
    from fastspecfit.io import one_spectrum, one_stacked_spectrum
    from fastspecfit.emlines import emline_specfit
    from fastspecfit.continuum import continuum_specfit

    igm = sc_data.igm
    phot = sc_data.photometry
    emline_table = sc_data.emlines.table
    templates = sc_data.templates

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
    if not fastphot:
        if not fitstack:
            flux = data['photometry']['nanomaggies']
            fluxivar = data['photometry']['nanomaggies_ivar']
            for iband, band in enumerate(phot.bands):
                meta[f'FLUX_{band.upper()}'] = flux[iband]
                meta[f'FLUX_IVAR_{band.upper()}'] = fluxivar[iband]

            if hasattr(phot, 'fiber_bands'):
                fibertotflux = data['fiberphot']['nanomaggies']
                for iband, band in enumerate(phot.fiber_bands):
                    meta[f'FIBERTOTFLUX_{band.upper()}'] = fibertotflux[iband]

    # output structure
    out = BoxedScalar(out_dtype)

    continuummodel, smooth_continuum, continuummodel_monte, specflux_monte = \
        continuum_specfit(data, out, templates, igm, phot, constrain_age=constrain_age,
                          no_smooth_continuum=no_smooth_continuum, fastphot=fastphot,
                          fitstack=fitstack, debug_plots=debug_plots, nmonte=nmonte,
                          seed=seed)

    # Optionally fit the emission-line spectrum.
    if fastphot:
        emmodel = None
    else:
        emmodel = emline_specfit(data, out, continuummodel, smooth_continuum,
                                 phot, emline_table, broadlinefit=broadlinefit,
                                 minsnr_balmer_broad=minsnr_balmer_broad,
                                 debug_plots=debug_plots, specflux_monte=specflux_monte,
                                 continuummodel_monte=continuummodel_monte)

    return out.value, meta, emmodel


def fastspec(fastphot=False, fitstack=False, args=None, comm=None, verbose=False):
    """Main fastspec script.

    This script is the engine to model one or more DESI spectra. It initializes
    the :class:`FastFit` class, reads the data, fits each spectrum (with the
    option of fitting in parallel), and writes out the results as a
    multi-extension binary FITS table.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line.
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    from astropy.table import vstack
    from fastspecfit.io import (DESISpectra, get_output_dtype, create_output_meta,
                                create_output_table, write_fastspecfit)

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    # check for mandatory environment variables
    envlist = []
    if args.specproddir is None:
        envlist += ['DESI_SPECTRO_REDUX']
    if args.mapdir is None:
        envlist += ['DUST_DIR']
    if args.fphotodir is None:
        envlist += ['FPHOTO_DIR']
    if args.templates is None:
        envlist += ['FTEMPLATES_DIR']
    for env in envlist:
        if not env in os.environ:
            errmsg = f'Mandatory environment variable {env} missing.'
            log.critical(errmsg)
            raise KeyError(errmsg)

    if fitstack:
        args.ignore_photometry = True

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

    # initialize single-copy objects im main process
    init_sc_args = {
        'emlines_file':      args.emlinesfile,
        'fphotofile':        args.fphotofile,
        'fastphot':          fastphot,
        'fitstack':          fitstack,
        'ignore_photometry': args.ignore_photometry,
        'template_file':     args.templates,
        'template_version':  args.templateversion,
        'template_imf':      args.imf,
        'log_verbose':       args.verbose
    }

    sc_data.initialize(**init_sc_args)

    # If multiprocessing, create a pool of worker processes and initialize
    # single-copy objects in each worker.
    if args.mp > 1 and not 'NERSC_HOST' in os.environ:
        import multiprocessing
        multiprocessing.set_start_method('fork')

    t0 = time.time()
    mp_pool = MPPool(args.mp,
                     initializer=sc_data.initialize,
                     init_argdict=init_sc_args)
    log.debug(f'Caching took {time.time()-t0:.5f} seconds.')

    log.info(f'Cached stellar templates {sc_data.templates.file}')
    log.info(f'Cached emission-line table {sc_data.emlines.file}')
    log.info(f'Cached photometric filters and parameters {sc_data.photometry.fphotofile}')
    log.info(f'Cached cosmology table {sc_data.cosmology.file}')
    log.info(f'Cached {sc_data.igm.reference} IGM attenuation parameters.')

    # Read the data.
    Spec = DESISpectra(phot=sc_data.photometry, cosmo=sc_data.cosmology,
                       fphotodir=args.fphotodir, mapdir=args.mapdir)

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
            return

        data, meta = Spec.read(sc_data.photometry, fastphot=fastphot, constrain_age=args.constrain_age)

    ntargets = len(meta)
    ncoeff = sc_data.templates.ntemplates
    if fastphot:
        cameras = None
    else:
        cameras = data[0]['cameras']

    out_dtype, out_units = get_output_dtype(
        Spec.specprod, phot=sc_data.photometry,
        linetable=sc_data.emlines.table, ncoeff=ncoeff,
        cameras=cameras, fastphot=fastphot,
        fitstack=fitstack)

    # If using Monte Carlo, generate the random seed(s).
    if args.nmonte > 0:
        if input_seeds is not None:
            seeds = input_seeds
        else:
            rng = np.random.default_rng(seed=args.seed)
            seeds = rng.integers(2**32, size=ntargets, dtype=np.int64)
    else:
        seeds = [1] * ntargets

    # Fit in parallel
    t0 = time.time()
    fitargs = [{
        'iobj':                iobj,
        'data':                data[iobj],
        'meta':                meta[iobj],
        'out_dtype':           out_dtype,
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
    } for iobj in range(len(meta))]

    _out = mp_pool.starmap(fastspec_one, fitargs)
    out = list(zip(*_out))

    meta = create_output_meta(vstack(out[1]), phot=sc_data.photometry,
                              fastphot=fastphot, fitstack=fitstack)

    results = create_output_table(out[0], meta, out_units,
                                  fitstack=fitstack)

    if fastphot:
        modelspectra = None
    else:
        modelspectra = vstack(out[2], join_type='exact', metadata_conflicts='error')

    # if multiprocessing, clean up workers
    mp_pool.close()

    log.info(f'Fitting {ntargets} object(s) took {time.time()-t0:.2f} seconds.')

    write_fastspecfit(results, meta, modelspectra=modelspectra, outfile=args.outfile,
                      specprod=Spec.specprod, coadd_type=Spec.coadd_type,
                      fphotofile=sc_data.photometry.fphotofile,
                      template_file=sc_data.templates.file,
                      emlinesfile=sc_data.emlines.file, fastphot=fastphot,
                      inputz=input_redshifts is not None,
                      nmonte=args.nmonte, seed=args.seed,
                      inputseeds=input_seeds is not None,
                      uncertainty_floor=args.uncertainty_floor,
                      minsnr_balmer_broad=args.minsnr_balmer_broad,
                      ignore_photometry=args.ignore_photometry,
                      broadlinefit=args.broadlinefit, constrain_age=args.constrain_age,
                      use_quasarnet=args.use_quasarnet,
                      no_smooth_continuum=args.no_smooth_continuum)


def fastphot(args=None, comm=None):
    """Main fastphot script.

    This script is the engine to model the broadband photometry of one or more
    DESI objects. It initializes the :class:`ContinuumFit` class, reads the
    data, fits each object (with the option of fitting in parallel), and writes
    out the results as a multi-extension binary FITS table.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line.
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    fastspec(fastphot=True, args=args, comm=comm)


def stackfit(args=None, comm=None):
    """Wrapper script to fit (model) generic stacked spectra.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line.
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    fastspec(fitstack=True, args=args, comm=comm)


def parse(options=None):
    """Parse input arguments to fastspec and fastphot scripts.

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
    parser.add_argument('--nmonte', type=int, default=50, help='Number of Monte Carlo realizations.')
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
    parser.add_argument('--specproddir', type=str, default=None, help='Optional directory name for the spectroscopic production.')
    parser.add_argument('--uncertainty-floor', type=float, default=0.01, help='Minimum fractional uncertainty to add in quadrature to the formal inverse variance spectrum.')
    parser.add_argument('--minsnr-balmer-broad', type=float, default=2.5, help='Minimum broad Balmer S/N to force broad+narrow-line model.')
    parser.add_argument('--debug-plots', action='store_true', help='Generate a variety of debugging plots (written to $PWD).')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    args = parser.parse_args(options)

    log.info(f'fastspec {" ".join(options)}')

    return args


