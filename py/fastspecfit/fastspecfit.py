#!/usr/bin/env python
"""
fastspecfit.fastspecfit
=======================

See sandbox/running-fastspecfit for examples.


"""
import pdb # for debugging

import os, time
import numpy as np

from desiutil.log import get_logger
log = get_logger()

## ridiculousness! - this seems to come from healpy, blarg
#import tempfile
#os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

def _fastspec_one(args):
    """Multiprocessing wrapper."""
    return _fastspec_one(*args)

def _fastphot_one(args):
    """Multiprocessing wrapper."""
    return _fastphot_one(*args)

def _desiqa_one(args):
    """Multiprocessing wrapper."""
    return _desiqa_one(*args)

def fastspec_one(iobj, data, out, meta, CFit, EMFit, solve_vdisp=False, verbose=False):
    """Multiprocessing wrapper to run :func:`fastspec` on a single object."""
    
    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()

    cfit, continuummodel, smooth_continuum = CFit.continuum_specfit(data, solve_vdisp=solve_vdisp)
    for col in cfit.colnames:
        out[col] = cfit[col]

    # Copy over the reddening-corrected fluxes -- messy!
    for iband, band in enumerate(CFit.fiber_bands):
        meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    for iband, band in enumerate(CFit.bands):
        meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} [targetid {}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    # Fit the emission-line spectrum.
    t0 = time.time()
    emfit = EMFit.fit(data, continuummodel, smooth_continuum, verbose=verbose)
    for col in emfit.colnames:
        out[col] = emfit[col]
    log.info('Line-fitting object {} (targetid={}) took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    return out, meta

def fastphot_one(iobj, data, out, meta, CFit):
    """Multiprocessing wrapper to run :func:`fastphot` on a single object."""

    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()
    cfit, _ = CFit.continuum_fastphot(data)
    for col in cfit.colnames:
        out[col] = cfit[col]

    # Copy over the reddening-corrected fluxes -- messy!
    for iband, band in enumerate(CFit.fiber_bands):
        meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    for iband, band in enumerate(CFit.bands):
        meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} [targetid {}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))
    
    return out, meta

def desiqa_one(CFit, EMFit, data, fastfit, metadata, coadd_type,
               fastphot=False, outdir=None, outprefix=None):
    """Multiprocessing wrapper to generate QA for a single object."""

    #t0 = time.time()
    if fastphot:
        CFit.qa_fastphot(data, fastfit, metadata, coadd_type=coadd_type,
                         outprefix=outprefix, outdir=outdir)
    else:
        EMFit.qa_fastspec(data, fastfit, metadata, coadd_type=coadd_type,
                          outprefix=outprefix, outdir=outdir)
    #log.info('Building took {:.2f} sec'.format(time.time()-t0))

def parse(options=None):
    """Parse input arguments to fastspec and fastphot scripts.

    """
    import argparse, sys

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('redrockfiles', nargs='*', help='Full path to input redrock file(s).')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename (required).')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads per MPI rank.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file, zero-indexed.') 
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--solve-vdisp', action='store_true', help='Solve for the velocity dispersion (only when using fastspec).')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    # make specprod required until this ticket is addressed--
    # https://github.com/desihub/desispec/issues/1077
    #parser.add_argument('--specprod', type=str, default='everest', #choices=['everest', 'denali', 'daily'],
    #                    help='Spectroscopic production.')
    #parser.add_argument('--coadd-type', type=str, default='healpix', choices=['healpix', 'cumulative', 'pernight', 'perexp'],
    #                    help='Type of spectral coadds corresponding to the input redrockfiles.')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspec {}'.format(' '.join(options)))

    return args

def fastspec(args=None, comm=None):
    """Main fastspec script.

    This script is the engine to model one or more DESI spectra. It initializes
    the :class:`ContinuumFit` and :class:`EMLineFit` classes, reads the data, fits
    each spectrum (with the option of fitting in parallel), and writes out the
    results as a multi-extension binary FITS table.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line. 
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum- and emission-line fitting classes.
    t0 = time.time()
    CFit = ContinuumFit()
    EMFit = EMLineFit()
    Spec = DESISpectra(specprod=args.specprod)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()

    Spec.find_specfiles(args.redrockfiles, firsttarget=args.firsttarget,
                        targetids=targetids, ntargets=args.ntargets)
    if len(Spec.specfiles) == 0:
        return
    data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, remember_coadd=True)

    out, meta = Spec.init_output(CFit=CFit, EMFit=EMFit, fastphot=False)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, EMFit, args.solve_vdisp, args.verbose)
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
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # Write out.
    write_fastspecfit(out, meta, outfile=args.outfile, specprod=Spec.specprod,
                      coadd_type=Spec.coadd_type, fastphot=False)

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
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum-fitting classes.
    t0 = time.time()
    CFit = ContinuumFit()
    Spec = DESISpectra(specprod=args.specprod)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    Spec.find_specfiles(args.redrockfiles, firsttarget=args.firsttarget,
                        targetids=targetids, ntargets=args.ntargets)
    if len(Spec.specfiles) == 0:
        return
    data = Spec.read_and_unpack(CFit, fastphot=True, synthphot=False)
    
    out, meta = Spec.init_output(CFit=CFit, fastphot=True)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit)
               for iobj in np.arange(Spec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastphot_one, fitargs)
    else:
        _out = [fastphot_one(*_fitargs) for _fitargs in fitargs]
    _out = list(zip(*_out))
    out = Table(np.hstack(_out[0]))
    meta = Table(np.hstack(_out[1]))
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # Write out.
    write_fastspecfit(out, meta, outfile=args.outfile, specprod=Spec.specprod,
                      coadd_type=Spec.coadd_type, fastphot=True)
