#!/usr/bin/env python
"""
fastspecfit.fastspecfit
=======================

See sandbox/running-fastspecfit for examples.

"""
import pdb # for debugging

import os, time
import numpy as np

def _fastspec_one(args):
    """Multiprocessing wrapper."""
    return fastspec_one(*args)

def _assign_units_to_columns(fastfit, metadata, Spec, templates, fastphot, stackfit, log):
    """Assign astropy units to output tables.

    """
    from fastspecfit.io import init_fastspec_output
    
    fastcols = fastfit.colnames
    metacols = metadata.colnames

    T, M = init_fastspec_output(Spec.meta, Spec.specprod, fphoto=Spec.fphoto, 
                                templates=templates, log=log, fastphot=fastphot, 
                                stackfit=stackfit)
    
    for col in T.colnames:
        if col in fastcols:
            fastfit[col].unit = T[col].unit
    for col in M.colnames:
        if col in metacols:
            metadata[col].unit = M[col].unit

def fastspec_one(iobj, data, out, meta, fphoto, templates, log=None,
                 emlinesfile=None, broadlinefit=True, fastphot=False,
                 constrain_age=False, no_smooth_continuum=False,
                 percamera_models=False, debug_plots=False):
    """Multiprocessing wrapper to run :func:`fastspec` on a single object.

    """
    from fastspecfit.io import cache_templates
    from fastspecfit.emlines import emline_specfit
    from fastspecfit.continuum import continuum_specfit

    log.info('Continuum- and emission-line fitting object {} [{} {}, z={:.6f}].'.format(
        iobj, fphoto['uniqueid'].lower(), data['uniqueid'], meta['Z']))

    # Read the templates and then fit the continuum spectrum. Note that 450 A as
    # the minimum wavelength will allow us to synthesize u-band photometry only
    # up to z=5.53, even though some targets are at higher redshift. Handle this
    # case in continuum.ContinuumTools.
    templatecache = cache_templates(templates, log=log, mintemplatewave=450.0,
                                    maxtemplatewave=40e4, fastphot=fastphot)

    continuummodel, smooth_continuum = continuum_specfit(data, out, templatecache, fphoto=fphoto,
                                                         emlinesfile=emlinesfile, constrain_age=constrain_age,
                                                         no_smooth_continuum=no_smooth_continuum,
                                                         fastphot=fastphot, debug_plots=debug_plots,
                                                         log=log)

    # Optionally fit the emission-line spectrum.
    if fastphot:
        emmodel = None
    else:
        emmodel = emline_specfit(data, templatecache, out, continuummodel, smooth_continuum,
                                 fphoto=fphoto, emlinesfile=emlinesfile, broadlinefit=broadlinefit,
                                 percamera_models=percamera_models, log=log)
        
    return out, meta, emmodel

def parse(options=None, log=None):
    """Parse input arguments to fastspec and fastphot scripts.

    """
    import argparse, sys

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('redrockfiles', nargs='+', help='Full path to input redrock file(s).')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename (required).')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads per MPI rank.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file, zero-indexed.') 
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--input-redshifts', type=str, default=None, help='Comma-separated list of input redshifts corresponding to the (required) --targetids input.')
    parser.add_argument('--no-broadlinefit', default=True, action='store_false', dest='broadlinefit',
                        help='Do not allow for broad Balmer and Helium line-fitting.')
    parser.add_argument('--nophoto', action='store_true', help='Do not include the photometry in the model fitting.')
    parser.add_argument('--constrain-age', action='store_true', help='Constrain the age of the SED.')
    parser.add_argument('--no-smooth-continuum', action='store_true', help='Do not fit the smooth continuum.')
    parser.add_argument('--ignore-quasarnet', dest='use_quasarnet', default=True, action='store_false', help='Do not use QuasarNet to improve QSO redshifts.')    
    parser.add_argument('--percamera-models', action='store_true', help='Return the per-camera (not coadded) model spectra.')
    parser.add_argument('--imf', type=str, default='chabrier', help='Initial mass function.')
    parser.add_argument('--templateversion', type=str, default='1.1.0', help='Template version number.')
    parser.add_argument('--templates', type=str, default=None, help='Optional full path and filename to the templates.')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--specfile-prefix', type=str, default='coadd-', help='Prefix of the spectral file(s).')
    parser.add_argument('--qnfile-prefix', type=str, default='qso_qn-', help='Prefix of the QuasarNet afterburner file(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--fphotodir', type=str, default=None, help='Top-level location of the source photometry.')
    parser.add_argument('--fphotofile', type=str, default=None, help='Photometric information file.')
    parser.add_argument('--emlinesfile', type=str, default=None, help='Emission line parameter file.')
    parser.add_argument('--specproddir', type=str, default=None, help='Optional directory name for the spectroscopic production.')
    parser.add_argument('--debug-plots', action='store_true', help='Generate a variety of debugging plots (written to $PWD).')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if log is None:
        from desiutil.log import get_logger
        log = get_logger()
            
    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspec {}'.format(' '.join(options)))

    return args

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
    from astropy.table import Table, vstack
    from desiutil.log import get_logger, DEBUG
    from fastspecfit.io import DESISpectra, write_fastspecfit, init_fastspec_output

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.verbose or verbose:
        log = get_logger(DEBUG)
    else:
        log = get_logger()

    if args.emlinesfile is None:
        from importlib import resources
        emlinesfile = resources.files('fastspecfit').joinpath('data/emlines.ecsv')
    else:
        emlinesfile = args.emlinesfile
        
    if not os.path.isfile(emlinesfile):
        errmsg = f'Emission lines parameter file {emlinesfile} does not exist.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    input_redshifts = None
    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
        if args.input_redshifts:
            input_redshifts = [float(x) for x in args.input_redshifts.split(',')]
            if len(input_redshifts) != len(targetids):
                errmsg = 'targetids and input_redshifts must have the same number of elements.'
                log.critical(errmsg)
                raise ValueError(errmsg)
    else:
        targetids = args.targetids

    # Read the data.
    t0 = time.time()
    Spec = DESISpectra(stackfit=stackfit, fphotodir=args.fphotodir,
                       fphotofile=args.fphotofile, mapdir=args.mapdir)

    if stackfit:
        data = Spec.read_stacked(args.redrockfiles, firsttarget=args.firsttarget,
                                 stackids=targetids, ntargets=args.ntargets,
                                 mp=args.mp)
        args.nophoto = True # force nophoto=True
    else:
        Spec.select(args.redrockfiles, firsttarget=args.firsttarget, targetids=targetids,
                    input_redshifts=input_redshifts, ntargets=args.ntargets,
                    redrockfile_prefix=args.redrockfile_prefix,
                    specfile_prefix=args.specfile_prefix, qnfile_prefix=args.qnfile_prefix,
                    use_quasarnet=args.use_quasarnet, specprod_dir=args.specproddir)
        if len(Spec.specfiles) == 0:
            return
    
        data = Spec.read_and_unpack(fastphot=fastphot, mp=args.mp, verbose=args.verbose)
        
    log.info('Reading and unpacking {} spectra to be fitted took {:.2f} seconds.'.format(
        Spec.ntargets, time.time()-t0))

    # Initialize the output tables.
    if args.templates is None:
        from fastspecfit.io import get_templates_filename
        templates = get_templates_filename(templateversion=args.templateversion, imf=args.imf)
    else:
        templates = args.templates

    out, meta = init_fastspec_output(Spec.meta, Spec.specprod, fphoto=Spec.fphoto, 
                                     templates=templates, data=data, log=log,
                                     emlinesfile=emlinesfile, fastphot=fastphot,
                                     stackfit=stackfit)

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], Spec.fphoto, templates, log,
                emlinesfile, args.broadlinefit, fastphot, args.constrain_age,
                args.no_smooth_continuum, args.percamera_models, args.debug_plots)
                for iobj in np.arange(Spec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastspec_one, fitargs)
    else:
        _out = [fastspec_one(*_fitargs) for _fitargs in fitargs]
        
    _out = list(zip(*_out))
    out = Table(np.hstack(_out[0]))
    meta = Table(np.hstack(_out[1]))

    if fastphot:
        modelspectra = None
    else:
        from astropy.table import vstack
        try:
            # need to vstack to preserve the wavelength metadata 
            modelspectra = vstack(_out[2], metadata_conflicts='error')
        except:
            errmsg = 'Metadata conflict when stacking model spectra.'
            log.critical(errmsg)
            raise ValueError(errmsg)
       
    log.info('Fitting {} object(s) took {:.2f} seconds.'.format(Spec.ntargets, time.time()-t0))

    # Assign units and write out.
    _assign_units_to_columns(out, meta, Spec, templates, fastphot, stackfit, log)

    write_fastspecfit(out, meta, modelspectra=modelspectra, outfile=args.outfile,
                      specprod=Spec.specprod, coadd_type=Spec.coadd_type,
                      fphotofile=Spec.fphotofile, templates=templates,
                      emlinesfile=emlinesfile, fastphot=fastphot,
                      inputz=input_redshifts is not None, nophoto=args.nophoto,
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
