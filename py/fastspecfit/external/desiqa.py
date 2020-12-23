#!/usr/bin/env python
"""
fastspecfit.external.desiqa
===========================


"""
import pdb # for debugging
import os, sys, time
import numpy as np

from desiutil.log import get_logger
log = get_logger()

# ridiculousness!
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import matplotlib
matplotlib.use('Agg')

import multiprocessing

def _desiqa_one(args):
    """Multiprocessing wrapper."""
    return desiqa_one(*args)

def desiqa_one(indx, data, specfit, photfit, CFit, EMFit, qadir):
    """QA on one spectrum."""
    t0 = time.time()
    continuum = CFit.qa_continuum(data, specfit, photfit, qadir)
    EMFit.qa_emlines(data, specfit, continuum, qadir)
    log.info('Building QA for object {} took {:.2f} sec'.format(indx, time.time()-t0))

def parse(options=None):
    """Parse input arguments.

    """
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required but with sensible defaults
    parser.add_argument('--night', default='20200225', type=str, help='Night to process.')
    parser.add_argument('--tile', default='70502', type=str, help='Tile number to process.')

    # optional inputs
    parser.add_argument('--first', type=int, help='Index of first spectrum to process (0-indexed).')
    parser.add_argument('--last', type=int, help='Index of last spectrum to process (max of nobj-1).')
    parser.add_argument('--targetid', type=np.int64, nargs='*', default=None, help='List of one or more targetids to consider.')

    parser.add_argument('--nproc', default=1, type=int, help='Number of cores.')
    parser.add_argument('--specprod', type=str, default='variance-model', choices=['andes', 'daily', 'variance-model'],
                        help='Spectroscopic production to process.')

    parser.add_argument('--use-vi', action='store_true', help='Select spectra with high-quality visual inspections (VI).')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite any existing files.')
    parser.add_argument('--no-write-spectra', dest='write_spectra', default=True, action='store_false',
                        help='Do not write out the selected spectra for the specified tile and night.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose.')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspecfit_qa {}'.format(' '.join(options)))

    return args

def main(args=None, comm=None):
    """Main module.

    """
    import fitsio
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit
    from fastspecfit.external.desi import read_spectra, unpack_all_spectra

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    fastspecfit_dir = os.getenv('FASTSPECFIT_DATA')
    resultsdir = os.path.join(fastspecfit_dir, 'results', args.specprod)
    qadir = os.path.join(fastspecfit_dir, 'qa', args.specprod)

    specfitfile = os.path.join(resultsdir, 'specfit-{}-{}.fits'.format(
        args.tile, args.night))
    photfitfile = os.path.join(resultsdir, 'photfit-{}-{}.fits'.format(
        args.tile, args.night))
    if not os.path.isfile(specfitfile):
        log.info('Spectroscopic fit file {} not found!'.format(specfitfile))
        return
    if not os.path.isfile(photfitfile):
        log.info('Photometric fit file {} not found!'.format(photfitfile))
        return
    specfit = Table(fitsio.read(specfitfile))
    log.info('Read {} objects from {}'.format(len(specfit), specfitfile))
    photfit = Table(fitsio.read(photfitfile))
    log.info('Read {} objects from {}'.format(len(photfit), photfitfile))
    assert(len(specfit) == len(photfit))

    # Read the data 
    zbest, specobj = read_spectra(tile=args.tile, night=args.night,
                                  specprod=args.specprod,
                                  use_vi=args.use_vi, 
                                  write_spectra=args.write_spectra,
                                  verbose=args.verbose)

    if args.targetid is not None:
        fitindx = np.where(np.isin(zbest['TARGETID'], args.targetid))[0]
        miss = np.where(np.logical_not(np.isin(args.targetid, zbest['TARGETID'])))[0]
        for this in miss:
            log.warning('TARGETID {} not found!'.format(args.targetid[this]))
        if len(fitindx) == 0:
            return
    else:
        if args.first is None:
            args.first = 0
        if args.last is None:
            args.last = len(zbest) - 1
        fitindx = np.arange(args.last - args.first + 1) + args.first
        
    # Initialize the continuum- and emission-line fitting classes.
    CFit = ContinuumFit(nproc=args.nproc, verbose=args.verbose)
    EMFit = EMLineFit()

    # Unpacking with multiprocessing takes a lot longer (maybe pickling takes a
    # long time?) so suppress the `nproc` argument here for now.
    data = unpack_all_spectra(specobj, zbest, CFit, fitindx)#, nproc=args.nproc)
    del specobj, zbest # free memory

    # Build the QA in parallel
    qaargs = [(indx, data[iobj], specfit[indx], photfit[indx], CFit, EMFit, qadir)
               for iobj, indx in enumerate(fitindx)]
    if args.nproc > 1:
        with multiprocessing.Pool(args.nproc) as P:
            P.map(_desiqa_one, qaargs)
    else:
        [desiqa_one(*_qaargs) for _qaargs in qaargs]
