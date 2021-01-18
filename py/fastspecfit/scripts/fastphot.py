#!/usr/bin/env python
"""
fastspecfit.scripts.fastphot
============================

Fastphot wrapper. Call with, e.g.,
  fastphot /global/cfs/cdirs/desi/spectro/redux/blanc/tiles/80607/deep/zbest-0-80607-deep.fits -o fastphot.fits --ntargets 2

"""
import pdb # for debugging

import os, time
import numpy as np

from desiutil.log import get_logger
log = get_logger()

# ridiculousness! - this seems to come from healpy, blarg
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

def _fastphot_one(args):
    """Multiprocessing wrapper."""
    return fastphot_one(*args)

def fastphot_one(iobj, data, out, meta, CFit, solve_vdisp=False):
    """Fit one spectrum."""
    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()
    cfit, _ = CFit.continuum_fastphot(data)
    for col in cfit.colnames:
        out[col] = cfit[col]

    # Copy over the reddening-corrected fluxes.
    for iband, band in enumerate(CFit.fiber_bands):
        meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    for iband, band in enumerate(CFit.bands):
        meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} took {:.2f} sec'.format(iobj, time.time()-t0))
    
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

    parser.add_argument('zbestfiles', nargs='*', help='Full path to input zbest file(s).')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastphot {}'.format(' '.join(options)))

    return args

def main(args=None, comm=None):
    """Main module.

    """
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.io import DESISpectra, write_fastspec

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum- and emission-line fitting classes.
    t0 = time.time()
    CFit = ContinuumFit()
    Spec = DESISpectra()
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    Spec.find_specfiles(args.zbestfiles, firsttarget=args.firsttarget,
                        targetids=targetids, ntargets=args.ntargets)
    if len(Spec.specfiles) == 0:
        return
    data = Spec.read_and_unpack(CFit, fastphot=True, synthphot=False)
    
    out, meta = Spec.init_output(CFit=CFit, fastphot=True)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, args.solve_vdisp)
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
    write_fastspec(out, meta, outfile=args.outfile, specprod=Spec.specprod, fastphot=True)
