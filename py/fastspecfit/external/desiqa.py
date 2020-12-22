#!/usr/bin/env python
"""
fastspecfit.external.desiqa
===========================


"""
import pdb # for debugging

import os, sys, time
import numpy as np
from desiutil.log import get_logger

# ridiculousness!
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import matplotlib
matplotlib.use('Agg')

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
    parser.add_argument('--nproc', default=1, type=int, help='Number of cores.')
    parser.add_argument('--use-vi', action='store_true', help='Select spectra with high-quality visual inspections (VI).')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite any existing files.')
    parser.add_argument('--no-write-spectra', dest='write_spectra', default=True, action='store_false',
                        help='Do not write out the selected spectra for the specified tile and night.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose.')

    log = get_logger()
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
    from astropy.table import Table
    from desigal.fastspecfit import read_spectra, unpack_all_spectra, ContinuumFit, EMLineFit

    log = get_logger()
    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    for key in ['FASTSPECFIT_DATA', 'FASTSPECFIT_TEMPLATES']:
        if key not in os.environ:
            log.fatal('Required ${} environment variable not set'.format(key))
            raise EnvironmentError('Required ${} environment variable not set'.format(key))

    fastspecfit_dir = os.getenv('FASTSPECFIT_DATA')
    resultsdir = os.path.join(fastspecfit_dir, 'results')
    qadir = os.path.join(fastspecfit_dir, 'qa')
    if not os.path.isdir(qadir):
        os.makedirs(qadir)

    fastspecfitfile = os.path.join(resultsdir, 'fastspecfit-{}-{}.fits'.format(
        args.tile, args.night))
    if not os.path.isfile(fastspecfitfile):
        log.info('Output file {} not found!'.format(fastspecfitfile))
        return
    fastspecfit = Table.read(fastspecfitfile)
    log.info('Read {} objects from {}'.format(len(fastspecfit), fastspecfitfile))

    # Read the data 
    zbest, specobj = read_spectra(tile=args.tile, night=args.night,
                                  use_vi=args.use_vi, 
                                  write_spectra=args.write_spectra,
                                  verbose=args.verbose)

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

    for iobj, indx in enumerate(fitindx):
        continuum = CFit.fnnls_continuum_plot(data[iobj], fastspecfit[indx], qadir=qadir)
        EMFit.emlineplot(data[iobj], fastspecfit[indx], continuum, qadir=qadir)

