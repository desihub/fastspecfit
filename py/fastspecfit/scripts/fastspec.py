#!/usr/bin/env python
"""
fastspecfit.scripts.fastspec
============================

FastSpec wrapper. Call with, e.g.,

  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/20210324/redrock-4-80613-thru20210324.fits -o fastspec.fits --targetids 39633345008634465

# nice BGS example
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/20210324/redrock-4-80613-thru20210324.fits -o fastspec.fits --targetids 39633345008634465

  # redrock is wrong!
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80605/redrock-0-80605-deep.fits -o fastspec.fits --targetids 39627652595714901

  # good test of needing smoothing continuum residuals before line-fitting
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80605/redrock-9-80605-deep.fits -o fastspec.fits --targetids 39627658622930703

  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/redrock-0-80613-deep.fits -o fastspec.fits --targetids 39633314155332057
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/redrock-0-80606-deep.fits -o fastspec.fits --ntargets 2

  # nice QSO with broad lines
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80607/redrock-9-80607-deep.fits -o fastspec2.fits --targetids 39633331528141827

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
    return fastspec_one(*args)

def fastspec_one(iobj, data, out, meta, CFit, EMFit, solve_vdisp=False):
    """Fit one spectrum."""
    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()

    cfit, continuummodel, smooth_continuum = CFit.continuum_specfit(data, solve_vdisp=solve_vdisp)
    for col in cfit.colnames:
        out[col] = cfit[col]

    ## Copy over the reddening-corrected fluxes -- messy!
    #for iband, band in enumerate(CFit.fiber_bands):
    #    meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
    #    #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    #for iband, band in enumerate(CFit.bands):
    #    meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
    #    meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} [targetid {}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))
    
    # Fit the emission-line spectrum.
    t0 = time.time()
    emfit = EMFit.fit(data, continuummodel, smooth_continuum)
    for col in emfit.colnames:
        out[col] = emfit[col]
    log.info('Line-fitting object {} (targetid={}) took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    return out, meta

def parse(options=None):
    """Parse input arguments.

    """
    import argparse, sys

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file (0-indexed).')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of target IDs to process.')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing processes per MPI rank or node.')
    #parser.add_argument('--suffix', type=str, default=None, help='Optional suffix for output filename.')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename.')
    parser.add_argument('--solve-vdisp', action='store_true', help='Solve for the velocity disperion.')

    parser.add_argument('--specprod', type=str, default='denali', choices=['everest', 'denali', 'daily'],
                        help='Spectroscopic production to process.')
    parser.add_argument('--coadd-type', type=str, default='cumulative', choices=['cumulative', 'pernight', 'perexp'],
                        help='Type of spectral coadds corresponding to the input redrockfiles.')

    parser.add_argument('redrockfiles', nargs='*', help='Full path to input redrock file(s).')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspec {}'.format(' '.join(options)))

    return args

def main(args=None, comm=None):
    """Main module.

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
    Spec = DESISpectra(specprod=args.specprod, coadd_type=args.coadd_type)
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
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, EMFit, args.solve_vdisp)
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
