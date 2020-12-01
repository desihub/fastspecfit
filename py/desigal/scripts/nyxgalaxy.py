#!/usr/bin/env python
"""Main module for nyxgalaxy.

ToDo:
* Generalize to work on coadded spectra.
* Fit to the photometry.
* Solve for vdisp, zcontinuum, and E(B-V).
* Add polynomial correction templates
* Capture ivar=0 problems.
* Correct for MW extinction.
* Fit Mg II.

ELG - 20200228 70005
BGS+MWS - 20200303 70500
BGS+MWS - 20200315 66003

new truth table for 70500 is here--
/global/cfs/cdirs/desi/sv/vi/TruthTables/Andes_reinspection/BGS

https://desi.lbl.gov/trac/wiki/SurveyValidation/TruthTables
https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=5720;filename=DESI_data_042820.pdf

"""
import pdb # for debugging

import os, sys
import numpy as np
from desiutil.log import get_logger

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
    parser.add_argument('--makeqa', action='store_true', help='Build QA output.')

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
        log.info('nyxgalaxy {}'.format(' '.join(options)))

    return args

def main(args=None):
    """Main module.

    """
    from desigal.nyxgalaxy import read_spectra, unpack_all_spectra, init_nyxgalaxy
    from desigal.nyxgalaxy import ContinuumFit, EMLineFit

    log = get_logger()
    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    for key in ['NYXGALAXY_DATA', 'NYXGALAXY_TEMPLATES']:
        if key not in os.environ:
            log.fatal('Required ${} environment variable not set'.format(key))
            raise EnvironmentError('Required ${} environment variable not set'.format(key))

    desigal_dir = os.getenv('NYXGALAXY_DATA')

    # If the output file exists, we're done!
    nyxgalaxyfile = os.path.join(desigal_dir, 'nyxgalaxy-{}-{}.fits'.format(args.tile, args.night))
    if os.path.isfile(nyxgalaxyfile) and not args.overwrite:
        log.info('Output file {} exists; all done!'.format(nyxgalaxyfile))
        return
        
    # Read the data 
    zbest, specobj = read_spectra(tile=args.tile, night=args.night,
                                  use_vi=args.use_vi, 
                                  write_spectra=args.write_spectra,
                                  verbose=args.verbose)

    # Initialize the continuum- and emission-line fitting classes and the output
    # data table.
    CFit = ContinuumFit(nproc=args.nproc, verbose=args.verbose)
    EMFit = EMLineFit()
    nyxgalaxy = init_nyxgalaxy(args.tile, args.night, zbest, specobj.fibermap, CFit)

    # Choose the set of spectra to fit and unpack them.
    if args.first is None:
        args.first = 0
    if args.last is None:
        args.last = len(zbest) - 1
    fitindx = np.arange(args.last - args.first + 1) + args.first

    data = unpack_all_spectra(specobj, zbest, CFit, fitindx, nproc=args.nproc)
    del specobj, zbest # free memory

    pdb.set_trace()
    
    # Fit each object in sequence.
    for iobj in fitindx:

        south = True
        galwave, galflux, galivar, galres, galphot, zredrock = unpack_one_spectrum(
            specobj, zbest, iobj, CFit, south=south)


        # fit the stellar continuum
        contfit, continuum = CFit.fnnls_continuum(galwave, galflux, galivar, galres,
                                                  galphot, zredrock, CFit.linetable)
        for col in ['coeff', 'chi2', 'dof', 'age', 'ebv', 'vdisp', 'z', 'phot_coeff']:
            nyxgalaxy['CONTINUUM_{}'.format(col).upper()][iobj] = contfit[col]

        # fit the emission-line spectrum and populate the output table
        emfit, emlinemodel = EMFit.fit(galwave, galflux, galivar, galres, continuum,
                                       zredrock, verbose=args.verbose)

        for col in ['d4000', 'd4000_model']:
            nyxgalaxy['{}'.format(col).upper()][iobj] = emfit['{}'.format(col)]

        for col in ['linevshift', 'linesigma']:
            for suffix in ['forbidden', 'balmer']:
                nyxgalaxy['{}_{}'.format(col, suffix).upper()][iobj] = emfit['{}_{}'.format(col, suffix)]
                nyxgalaxy['{}_{}_IVAR'.format(col, suffix).upper()][iobj] = emfit['{}_{}_ivar'.format(col, suffix)]

        for line in emfit['linenames']:
            for suffix in ['chi2', 'npix', 'flux_limit', 'ew_limit']:
                nyxgalaxy['{}_{}'.format(line, suffix).upper()][iobj] = emfit['{}_{}'.format(line, suffix)]
            for suffix in ['amp', 'flux', 'boxflux', 'ew', 'cont']:
                nyxgalaxy['{}_{}'.format(line, suffix).upper()][iobj] = emfit['{}_{}'.format(line, suffix)]
                nyxgalaxy['{}_{}_IVAR'.format(line, suffix).upper()][iobj] = emfit['{}_{}_ivar'.format(line, suffix)]

        #for col in nyxgalaxy.colnames[-14:]:
        #    log.info('{:.4f} {}'.format(nyxgalaxy[col][iobj], col))
        #pdb.set_trace()

    # write out
    log.info('Writing {} spectra to {}'.format(len(nyxgalaxy), nyxgalaxyfile))
    nyxgalaxy.write(nyxgalaxyfile, overwrite=True)

    if args.makeqa:
        from astropy.table import Table
        nyxgalaxy = Table.read(nyxgalaxyfile)
        log.info('Read {} objects from {}'.format(len(nyxgalaxy), nyxgalaxyfile))

        qadir = os.path.join(desigal_dir, 'qa')
        if not os.path.isdir(qadir):
            os.makedirs(qadir)

        for iobj in fitindx:
            south = True

            targetid = nyxgalaxy['TARGETID'][iobj]
            galwave, galflux, galivar, galres, galphot, zredrock = _unpack_spectrum(
                specobj, zbest, iobj, CFit, south=south)

            continuum = CFit.fnnls_continuum_bestfit(nyxgalaxy['CONTINUUM_COEFF'][iobj], galwave=galwave,
                                                     galres=galres, redshift=zredrock)
            continuum_fullwave, fullwave = CFit.fnnls_continuum_bestfit(nyxgalaxy['CONTINUUM_PHOT_COEFF'][iobj],
                                                                        redshift=zredrock)

            emlinemodel = EMFit.emlinemodel_bestfit(galwave, galres, nyxgalaxy[iobj])

            # continuum fit
            pngfile = os.path.join(qadir, 'continuum-{}-{}-{}.png'.format(args.tile, args.night, targetid))
            objinfo = {
                'targetid': '{} {}'.format(zbest['TARGETID'][iobj], -999),
                #'targetid': 'TARGETID={} fiber={}'.format(zbest['TARGETID'][iobj], -999),
                'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(nyxgalaxy['CONTINUUM_CHI2'][iobj]),
                'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(nyxgalaxy['Z'][iobj]),
                'znyxgalaxy': '$z_{{\\rm nyxgalaxy}}$={:.6f}'.format(nyxgalaxy['CONTINUUM_Z'][iobj]),
                'age': '<Age>={:.3f} Gyr'.format(nyxgalaxy['CONTINUUM_AGE'][iobj]),
                'vdisp': '$\sigma$={:.1f} km/s'.format(nyxgalaxy['CONTINUUM_VDISP'][iobj]),
                'ebv': 'E(B-V)={:.4f} km/s'.format(nyxgalaxy['CONTINUUM_EBV'][iobj]),
                }

            if south:
                filters = CFit.decamwise
            else:
                filters = CFit.bassmzlswise
            filtwave = filters.effective_wavelengths.value

            CFit.fnnls_continuum_plot(galwave, galflux, galivar, galphot, continuum, 
                                      continuum_fullwave, fullwave, objinfo, png=pngfile)

            pdb.set_trace()

            # emission-line fit
            pngfile = os.path.join(qadir, 'emlinefit-{}-{}-{}.png'.format(args.tile, args.night, targetid))
            objinfo = {
                'targetid': '{} {}'.format(zbest['TARGETID'][iobj], -999),
                'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(nyxgalaxy['Z'][iobj]),
                'linevshift_forbidden': '$\Delta\,v_{{\\rm forbidden}}$={:.1f} km/s'.format(nyxgalaxy['LINEVSHIFT_FORBIDDEN'][iobj]),
                'linevshift_balmer': '$\Delta\,v_{{\\rm Balmer}}$={:.1f} km/s'.format(nyxgalaxy['LINEVSHIFT_BALMER'][iobj]),
                'linesigma_forbidden': '$\sigma_{{\\rm forbidden}}$={:.1f} km/s'.format(nyxgalaxy['LINESIGMA_FORBIDDEN'][iobj]),
                'linesigma_balmer': '$\sigma_{{\\rm Balmer}}$={:.1f} km/s'.format(nyxgalaxy['LINESIGMA_BALMER'][iobj]),
                }
            EMFit.emlineplot(galwave, galflux, galivar, continuum,
                             emlinemodel, zredrock, objinfo, png=pngfile)
