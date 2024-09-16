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


def fastspec_one(iobj, data, out_dtype, broadlinefit=True, fastphot=False,
                 constrain_age=False, no_smooth_continuum=False,
                 percamera_models=False, debug_plots=False,
                 minsnr_balmer_broad=2.5):
    """Run :func:`fastspec` on a single object.

    """
    from fastspecfit.emlines import emline_specfit
    from fastspecfit.continuum import continuum_specfit

    igm = sc_data.igm
    phot = sc_data.photometry
    emline_table = sc_data.emlines.table
    templates = sc_data.templates

    if fastphot:
        log.info(f'Continuum fitting object {iobj} [{phot.uniqueid_col.lower()} ' + \
                 f'{data["uniqueid"]}, z={data["redshift"]:.6f}].')
    else:
        log.info(f'Continuum- and emission-line fitting object {iobj} [{phot.uniqueid_col.lower()} ' + \
                 f'{data["uniqueid"]}, z={data["redshift"]:.6f}].')

    # output structure
    out = BoxedScalar(out_dtype)

    continuummodel, smooth_continuum = continuum_specfit(
        data, out, templates, igm, phot, constrain_age=constrain_age,
        no_smooth_continuum=no_smooth_continuum,
        fastphot=fastphot, debug_plots=debug_plots)

    # Optionally fit the emission-line spectrum.
    if fastphot:
        emmodel = None
    else:
        emmodel = emline_specfit(data, out, continuummodel, smooth_continuum,
                                 phot, emline_table,
                                 broadlinefit=broadlinefit,
                                 minsnr_balmer_broad=minsnr_balmer_broad,
                                 percamera_models=percamera_models)

    return out.value, emmodel


def fastspec(fastphot=False, stackfit=False, args=None, comm=None, verbose=False):
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
    from fastspecfit.io import DESISpectra, write_fastspecfit

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

    if stackfit:
        args.ignore_photometry = True

    if verbose:
        args.verbose = True

    if args.verbose:
        from fastspecfit.logger import DEBUG
        log.setLevel(DEBUG)

    targetids = None
    input_redshifts = None

    if args.targetids is not None:
        targetids = [ int(x) for x in args.targetids.split(',') ]
        if args.input_redshifts is not None:
            input_redshifts = [ float(x) for x in args.input_redshifts.split(',') ]
            if len(input_redshifts) != len(targetids):
                errmsg = 'targetids and input_redshifts must have the same number of elements.'
                log.critical(errmsg)
                raise ValueError(errmsg)

    # initialize single-copy objects im main process
    init_sc_args = {
        'emlines_file':      args.emlinesfile,
        'fphotofile':        args.fphotofile,
        'fastphot':          fastphot,
        'stackfit':          stackfit,
        'ignore_photometry': args.ignore_photometry,
        'template_file':     args.templates,
        'template_version':  args.templateversion,
        'template_imf':      args.imf,
    }

    sc_data.initialize(**init_sc_args)

    # if multiprocessing, create a pool of worker processes
    # and initialize single-copy objects in each worker
    #t0 = time.time()
    mp_pool = MPPool(args.mp,
                     initializer=sc_data.initialize,
                     init_argdict=init_sc_args)
    #log.info(f'Caching took {time.time()-t0:.5f} seconds.')

    if sc_data.templates.use_legacy_fitting:
        log.warning(f'Fitting with deprecated spectrophotometric templates (version={sc_data.templates.version})!')

    log.info(f'Cached stellar templates {sc_data.templates.file}')
    log.info(f'Cached emission-line table {sc_data.emlines.file}')
    log.info(f'Cached photometric filters and parameters {sc_data.photometry.fphotofile}')
    log.info(f'Cached cosmology table {sc_data.cosmology.file}')
    log.info(f'Cached {sc_data.igm.reference} IGM attenuation parameters.')

    # Read the data.
    Spec = DESISpectra(phot=sc_data.photometry, cosmo=sc_data.cosmology,
                       fphotodir=args.fphotodir, mapdir=args.mapdir)

    if stackfit:
        data = Spec.read_stacked(mp_pool, args.redrockfiles, firsttarget=args.firsttarget,
                                 stackids=targetids, ntargets=args.ntargets)
    else:
        Spec.gather_metadata(args.redrockfiles, firsttarget=args.firsttarget,
                             targetids=targetids, input_redshifts=input_redshifts,
                             ntargets=args.ntargets, zmin=args.zmin,
                             redrockfile_prefix=args.redrockfile_prefix,
                             specfile_prefix=args.specfile_prefix, qnfile_prefix=args.qnfile_prefix,
                             use_quasarnet=args.use_quasarnet, specprod_dir=args.specproddir)
        if len(Spec.specfiles) == 0:
            return

        data = Spec.read(mp_pool, fastphot=fastphot, debug_plots=args.debug_plots,
                         constrain_age=args.constrain_age)

    ncoeff = sc_data.templates.ntemplates
    out_dtype, out_units = get_output_dtype(Spec.specprod,
                                            phot=sc_data.photometry,
                                            linetable=sc_data.emlines.table,
                                            ncoeff=ncoeff,
                                            fastphot=fastphot, stackfit=stackfit)

    # Fit in parallel
    t0 = time.time()
    fitargs = [{
        'iobj':                iobj,
        'data':                data[iobj],
        'out_dtype':           out_dtype,
        'broadlinefit':        args.broadlinefit,
        'fastphot':            fastphot,
        'constrain_age':       args.constrain_age,
        'no_smooth_continuum': args.no_smooth_continuum,
        'percamera_models':    args.percamera_models,
        'debug_plots':         args.debug_plots,
        'minsnr_balmer_broad': args.minsnr_balmer_broad,
    } for iobj in range(Spec.ntargets)]

    _out = mp_pool.starmap(fastspec_one, fitargs)
    out = list(zip(*_out))

    meta = create_output_meta(Spec.meta, data,
                              phot=sc_data.photometry,
                              fastphot=fastphot, stackfit=stackfit)

    results = create_output_table(out[0], meta, out_units,
                                  stackfit=stackfit)

    if fastphot:
        modelspectra = None
    else:
        from astropy.table import vstack # preserves metadata
        modelspectra = vstack(out[1], join_type='exact', metadata_conflicts='error')

    # if multiprocessing, clean up workers
    mp_pool.close()

    log.info(f'Fitting {Spec.ntargets} object(s) took {time.time()-t0:.2f} seconds.')

    write_fastspecfit(results, meta, modelspectra=modelspectra, outfile=args.outfile,
                      specprod=Spec.specprod, coadd_type=Spec.coadd_type,
                      fphotofile=sc_data.photometry.fphotofile,
                      template_file=sc_data.templates.file,
                      emlinesfile=sc_data.emlines.file, fastphot=fastphot,
                      inputz=input_redshifts is not None,
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
    fastspec(stackfit=True, args=args, comm=comm)


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
    parser.add_argument('--zmin', type=float, default=None, help='Override the default minimum redshift required for modeling.')
    parser.add_argument('--no-broadlinefit', default=True, action='store_false', dest='broadlinefit',
                        help='Do not model broad Balmer and helium line-emission.')
    parser.add_argument('--ignore-photometry', default=False, action='store_true', help='Ignore the broadband photometry during model fitting.')
    parser.add_argument('--ignore-quasarnet', dest='use_quasarnet', default=True, action='store_false', help='Do not use QuasarNet to improve QSO redshifts.')
    parser.add_argument('--constrain-age', action='store_true', help='Constrain the age of the SED.')
    parser.add_argument('--no-smooth-continuum', action='store_true', help='Do not fit the smooth continuum.')
    parser.add_argument('--percamera-models', action='store_true', help='Return the per-camera (not coadded) model spectra.')
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
    parser.add_argument('--minsnr-balmer-broad', type=float, default=3., help='Minimum broad Balmer S/N to force broad+narrow-line model.')
    parser.add_argument('--debug-plots', action='store_true', help='Generate a variety of debugging plots (written to $PWD).')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    args = parser.parse_args(options)

    log.info(f'fastspec {" ".join(options)}')

    return args


def get_output_dtype(specprod, phot, linetable, ncoeff,
                     fastphot=False, stackfit=False):
    """
    Get type of one fastspecfit output data record, along
    with dictionary of units for any fields that have them.

    """
    import astropy.units as u

    out_dtype = []
    out_units = {}

    def add_field(name, dtype, shape=None, unit=None):
        if shape is not None:
            t = (name, dtype, shape) # subarray
        else:
            t = (name, dtype)
        out_dtype.append(t)

        if unit is not None:
            out_units[name] = unit

    add_field('Z', dtype='f8') # redshift
    add_field('COEFF', shape=(ncoeff,), dtype='f4')

    if not fastphot:
        add_field('RCHI2', dtype='f4')      # full-spectrum reduced chi2
        add_field('RCHI2_CONT', dtype='f4') # rchi2 fitting just to the continuum (spec+phot)
    add_field('RCHI2_PHOT', dtype='f4') # rchi2 fitting just to the photometry

    if stackfit:
        for cam in ('BRZ'):
            add_field(f'SNR_{cam}', dtype='f4') # median S/N in each camera
        for cam in ('BRZ'):
            add_field(f'SMOOTHCORR_{cam}', dtype='f4')
    else:
        if not fastphot:
            # if the zeroth object has a fully masked camera, this data model will fail
            #if data is not None:
            #    for cam in data[0]['cameras']:
            #        add_field(f'SNR_{cam.upper()}'), dtype='f4') # median S/N in each camera
            #    for cam in data[0]['cameras']:
            #        add_field(f'SMOOTHCORR_{cam.upper()}'), dtype='f4')
            for cam in ('B', 'R', 'Z'):
                add_field(f'SNR_{cam.upper()}', dtype='f4') # median S/N in each camera
            for cam in ('B', 'R', 'Z'):
                add_field(f'SMOOTHCORR_{cam.upper()}', dtype='f4')

    add_field('VDISP', dtype='f4', unit=u.kilometer/u.second)
    if not fastphot:
        add_field('VDISP_IVAR', dtype='f4', unit=u.second**2/u.kilometer**2)
    add_field('AV', dtype='f4', unit=u.mag)
    add_field('AGE', dtype='f4', unit=u.Gyr)
    add_field('ZZSUN', dtype='f4')
    add_field('LOGMSTAR', dtype='f4', unit=u.solMass)
    add_field('SFR', dtype='f4', unit=u.solMass/u.year)
    #add_field('FAGN', dtype='f4')

    if not fastphot:
        add_field('DN4000', dtype='f4')
        add_field('DN4000_OBS', dtype='f4')
        add_field('DN4000_IVAR', dtype='f4')
    add_field('DN4000_MODEL', dtype='f4')

    if not fastphot:
        # observed-frame photometry synthesized from the spectra
        for band in phot.synth_bands:
            add_field(f'FLUX_SYNTH_{band.upper()}', dtype='f4', unit='nanomaggies')
            #add_field(f'FLUX_SYNTH_IVAR_{band.upper()}'), dtype='f4', unit='1/nanomaggies**2')
        # observed-frame photometry synthesized the best-fitting spectroscopic model
        for band in phot.synth_bands:
            add_field(f'FLUX_SYNTH_SPECMODEL_{band.upper()}', dtype='f4', unit='nanomaggies')
    # observed-frame photometry synthesized the best-fitting continuum model
    for band in phot.bands:
        add_field(f'FLUX_SYNTH_PHOTMODEL_{band.upper()}', dtype='f4', unit='nanomaggies')

    for band, shift in zip(phot.absmag_bands, phot.band_shift):
        band = band.upper()
        shift = int(10*shift)
        add_field(f'ABSMAG{shift:02d}_{band}', dtype='f4', unit=u.mag) # absolute magnitudes
        add_field(f'ABSMAG{shift:02d}_IVAR_{band}', dtype='f4', unit=1/u.mag**2)
        add_field(f'KCORR{shift:02d}_{band}', dtype='f4', unit=u.mag)

    for cflux in ('LOGLNU_1500', 'LOGLNU_2800'):
        add_field(cflux,  dtype='f4', unit=10**(-28)*u.erg/u.second/u.Hz)
    add_field('LOGL_1450', dtype='f4', unit=10**(10)*u.solLum)
    add_field('LOGL_1700', dtype='f4', unit=10**(10)*u.solLum)
    add_field('LOGL_3000', dtype='f4', unit=10**(10)*u.solLum)
    add_field('LOGL_5100', dtype='f4', unit=10**(10)*u.solLum)

    for cflux in ('FLYA_1215_CONT', 'FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT'):
        add_field(cflux, dtype='f4', unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))

    if not fastphot:
        # Add chi2 metrics
        #add_field('DOF',  dtype='i8') # full-spectrum dof
        add_field('RCHI2_LINE', dtype='f4') # reduced chi2 with broad line-emission
        #add_field('NDOF_LINE', dtype='i8') # number of degrees of freedom corresponding to rchi2_line
        #add_field('DOF_BROAD', dtype='i8')
        add_field('DELTA_LINECHI2', dtype='f4') # delta-reduced chi2 with and without broad line-emission
        add_field('DELTA_LINENDOF', dtype=np.int32)

        # aperture corrections
        add_field('APERCORR', dtype='f4') # median aperture correction
        for band in phot.synth_bands:
            add_field(f'APERCORR_{band.upper()}', dtype='f4')

        add_field('NARROW_Z', dtype='f8')
        add_field('NARROW_ZRMS', dtype='f8')
        add_field('BROAD_Z', dtype='f8')
        add_field('BROAD_ZRMS', dtype='f8')
        add_field('UV_Z', dtype='f8')
        add_field('UV_ZRMS', dtype='f8')

        add_field('NARROW_SIGMA', dtype='f4', unit=u.kilometer / u.second)
        add_field('NARROW_SIGMARMS', dtype='f4', unit=u.kilometer / u.second)
        add_field('BROAD_SIGMA', dtype='f4', unit=u.kilometer / u.second)
        add_field('BROAD_SIGMARMS', dtype='f4', unit=u.kilometer / u.second)
        add_field('UV_SIGMA', dtype='f4', unit=u.kilometer / u.second)
        add_field('UV_SIGMARMS', dtype='f4', unit=u.kilometer / u.second)

        # special columns for the fitted doublets
        add_field('MGII_DOUBLET_RATIO', dtype='f4')
        add_field('OII_DOUBLET_RATIO', dtype='f4')
        add_field('SII_DOUBLET_RATIO', dtype='f4')

        for line in linetable['name']:
            line = line.upper()
            add_field(f'{line}_MODELAMP', dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
            add_field(f'{line}_AMP', dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
            add_field(f'{line}_AMP_IVAR', dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2)
            add_field(f'{line}_FLUX', dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2))
            add_field(f'{line}_FLUX_IVAR', dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4/u.erg**2)
            add_field(f'{line}_BOXFLUX', dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2))
            add_field(f'{line}_BOXFLUX_IVAR', dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4/u.erg**2)

            add_field(f'{line}_VSHIFT', dtype='f4',
                                  unit=u.kilometer/u.second)
            add_field(f'{line}_SIGMA', dtype='f4',
                                  unit=u.kilometer / u.second)

            add_field(f'{line}_CONT', dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
            add_field(f'{line}_CONT_IVAR', dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2)
            add_field(f'{line}_EW', dtype='f4',
                                  unit=u.Angstrom)
            add_field(f'{line}_EW_IVAR', dtype='f4',
                                  unit=1/u.Angstrom**2)
            add_field(f'{line}_FLUX_LIMIT', dtype='f4',
                                  unit=u.erg/(u.second*u.cm**2))
            add_field(f'{line}_EW_LIMIT', dtype='f4',
                                  unit=u.Angstrom)
            add_field(f'{line}_CHI2', dtype='f4')
            add_field(f'{line}_NPIX', dtype=np.int32)

    return np.dtype(out_dtype, align=True), out_units


def create_output_meta(input_meta, data, phot,
                       fastphot=False, stackfit=False):
    """Create the fastspecfit output metadata table.

    """
    from fastspecfit.io import TARGETINGBITS
    from astropy.table import Table
    import astropy.units as u

    nobj = len(input_meta)

    # The information stored in the metadata table depends on which spectra
    # were fitted (exposures, nightly coadds, deep coadds).
    if stackfit:
        fluxcols = ['PHOTSYS']
    else:
        fluxcols = []
        if hasattr(phot, 'outcols'):
            fluxcols.extend(phot.outcols)
        if hasattr(phot, 'fiber_bands'):
            fluxcols.extend([f'FIBERFLUX_{band.upper()}' for band in phot.fiber_bands])
            fluxcols.extend([f'FIBERTOTFLUX_{band.upper()}' for band in phot.fiber_bands])
        fluxcols.extend(phot.fluxcols)
        fluxcols.extend(phot.fluxivarcols)
        fluxcols.append('EBV')
        fluxcols.extend([f'MW_TRANSMISSION_{band.upper()}' for band in phot.bands])

    colunit = {'RA': u.deg, 'DEC': u.deg, 'EBV': u.mag}
    for fcol, icol in zip(phot.fluxcols, phot.fluxivarcols):
        colunit[fcol.upper()] = phot.photounits
        colunit[icol.upper()] = f'{phot.photounits}-2'
    if hasattr(phot, 'fiber_bands'):
        for band in phot.fiber_bands:
            band = band.upper()
            colunit[f'FIBERFLUX_{band}'] = phot.photounits
            colunit[f'FIBERTOTFLUX_{band}'] = phot.photounits

    skipcols = fluxcols + ['OBJTYPE', 'TARGET_RA', 'TARGET_DEC', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID', 'RELEASE']

    if stackfit:
        redrockcols = ('Z')
    else:
        redrockcols = ('Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE', 'Z_RR', 'TSNR2_BGS',
                       'TSNR2_LRG', 'TSNR2_ELG', 'TSNR2_QSO', 'TSNR2_LYA')

    meta = Table()
    metacols = set(input_meta.colnames)

    # All of this business is so we can get the columns in the order we want
    # (i.e., the order that matches the data model).
    if stackfit:
        for metacol in ('STACKID', 'SURVEY', 'PROGRAM'):
            if metacol in metacols:
                meta[metacol] = input_meta[metacol]
    else:
        for metacol in ('TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'NIGHT', 'FIBER',
                        'EXPID', 'TILEID_LIST', 'RA', 'DEC', 'COADD_FIBERSTATUS'):
            if metacol in metacols:
                meta[metacol] = input_meta[metacol]

        if 'SURVEY' in meta.colnames:
            if np.any(np.isin(meta['SURVEY'], 'main')) or np.any(np.isin(meta['SURVEY'], 'special')):
                TARGETINGCOLS = TARGETINGBITS['default']
            else:
                TARGETINGCOLS = TARGETINGBITS['all']
        else:
            TARGETINGCOLS = TARGETINGBITS['all']

        for metacol in metacols:
            if metacol in skipcols or metacol in TARGETINGCOLS or metacol in meta.colnames or metacol in redrockcols:
                continue
            else:
                meta[metacol] = input_meta[metacol]

        for bitcol in TARGETINGCOLS:
            if bitcol in metacols:
                meta[bitcol] = input_meta[bitcol]
            else:
                meta[bitcol] = np.zeros(shape=(1,), dtype=np.int64)

    for redrockcol in redrockcols:
        if redrockcol in metacols: # the Z_RR from quasarnet may not be present
            meta[redrockcol] = input_meta[redrockcol]

    for fluxcol in fluxcols:
        meta[fluxcol] = input_meta[fluxcol]

    # assign units to any columns that should have them
    for col in meta.colnames:
        if col in colunit:
            meta[col].unit = colunit[col]

    # copy some values from the input data to the output metadata
    for iobj, _data in enumerate(data):
        if not fastphot:
            if not stackfit:
                if hasattr(phot, 'fiber_bands'):
                    fibertotflux = _data['fiberphot']['nanomaggies']
                    #fibertotfluxivar = _data['fiberphot']['nanomaggies_ivar']
                    for iband, band in enumerate(phot.fiber_bands):
                        meta[f'FIBERTOTFLUX_{band.upper()}'][iobj] = fibertotflux[iband]
                        #result['FIBERTOTFLUX_IVAR_{band.upper()}'] = fibertotfluxivar[iband]

            flux = _data['photometry']['nanomaggies']
            fluxivar = _data['photometry']['nanomaggies_ivar']
            for iband, band in enumerate(phot.bands):
                meta[f'FLUX_{band.upper()}'][iobj] = flux[iband]
                meta[f'FLUX_IVAR_{band.upper()}'][iobj] = fluxivar[iband]

    return meta


def create_output_table(result_records, meta, units, stackfit=False):

    from astropy.table import hstack

    # Initialize the output table from the metadata table.
    metacols = set(meta.colnames)

    if stackfit:
        initcols = ('STACKID', 'SURVEY', 'PROGRAM')
    else:
        initcols = ('TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'NIGHT', 'FIBER', 'EXPID')
    initcols = [col for col in initcols if col in metacols]

    cdata = [meta[col] for col in initcols]
    results = Table()
    results.add_columns(cdata)

    # Now add the measurements. Columns and their dtypes are inferred from the
    # array's dtype.
    results = hstack((results, Table(np.array(result_records), units=units)))

    return results
