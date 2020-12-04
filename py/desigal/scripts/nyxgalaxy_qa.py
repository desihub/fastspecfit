#!/usr/bin/env python
"""Main QA module for nyxgalaxy.

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
    parser.add_argument('--no-write-spectra', dest='write_spectra', default=True, action='store_false',
                        help='Do not write out the selected spectra for the specified tile and night.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose.')

    log = get_logger()
    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('nyxgalaxy_qa {}'.format(' '.join(options)))

    return args

def main(args=None):
    """Main module.

    """
    from astropy.table import Table
    from desigal.nyxgalaxy import read_spectra, unpack_all_spectra, ContinuumFit, EMLineFit

    log = get_logger()
    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    for key in ['NYXGALAXY_DATA', 'NYXGALAXY_TEMPLATES']:
        if key not in os.environ:
            log.fatal('Required ${} environment variable not set'.format(key))
            raise EnvironmentError('Required ${} environment variable not set'.format(key))

    nyxgalaxy_dir = os.getenv('NYXGALAXY_DATA')
    qadir = os.path.join(nyxgalaxy_dir, 'qa')
    if not os.path.isdir(qadir):
        os.makedirs(qadir)

    nyxgalaxyfile = os.path.join(nyxgalaxy_dir, 'nyxgalaxy-{}-{}.fits'.format(args.tile, args.night))
    if not os.path.isfile(nyxgalaxyfile):
        log.info('Output file {} not found!'.format(nyxgalaxyfile))
        return
    nyxgalaxy = Table.read(nyxgalaxyfile)
    log.info('Read {} objects from {}'.format(len(nyxgalaxy), nyxgalaxyfile))

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
        CFit.fnnls_continuum_plot(data[iobj], nyxgalaxy[indx], qadir=qadir)

        #south = True
        #targetid = nyxgalaxy['TARGETID'][indx]
        #continuum = CFit.fnnls_continuum_bestfit(nyxgalaxy['CONTINUUM_COEFF'][indx], specwave=specwave,
        #                                         specres=specres, redshift=zredrock)
        #continuum_fullwave, fullwave = CFit.fnnls_continuum_bestfit(nyxgalaxy['CONTINUUM_PHOT_COEFF'][indx],
        #                                                            redshift=zredrock)
        #emlinemodel = EMFit.emlinemodel_bestfit(specwave, specres, nyxgalaxy[indx])
        #
        ## continuum fit
        #pngfile = os.path.join(qadir, 'continuum-{}-{}-{}.png'.format(args.tile, args.night, targetid))
        #
        #if south:
        #    filters = CFit.decamwise
        #else:
        #    filters = CFit.bassmzlswise
        #filtwave = filters.effective_wavelengths.value
        #
        #CFit.fnnls_continuum_plot(specwave, specflux, specivar, galphot, continuum, 
        #                          continuum_fullwave, fullwave, objinfo, png=pngfile)
        #
        #pdb.set_trace()
        #
        ## emission-line fit
        #pngfile = os.path.join(qadir, 'emlinefit-{}-{}-{}.png'.format(args.tile, args.night, targetid))
        #objinfo = {
        #    'targetid': '{} {}'.format(zbest['TARGETID'][indx], -999),
        #    'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(nyxgalaxy['Z'][indx]),
        #    'linevshift_forbidden': '$\Delta\,v_{{\\rm forbidden}}$={:.1f} km/s'.format(nyxgalaxy['LINEVSHIFT_FORBIDDEN'][indx]),
        #    'linevshift_balmer': '$\Delta\,v_{{\\rm Balmer}}$={:.1f} km/s'.format(nyxgalaxy['LINEVSHIFT_BALMER'][indx]),
        #    'linesigma_forbidden': '$\sigma_{{\\rm forbidden}}$={:.1f} km/s'.format(nyxgalaxy['LINESIGMA_FORBIDDEN'][indx]),
        #    'linesigma_balmer': '$\sigma_{{\\rm Balmer}}$={:.1f} km/s'.format(nyxgalaxy['LINESIGMA_BALMER'][indx]),
        #    }
        #EMFit.emlineplot(specwave, specflux, specivar, continuum,
        #                 emlinemodel, zredrock, objinfo, png=pngfile)
