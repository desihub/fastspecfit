#!/usr/bin/env python
"""
fastspecfit.external.desi
=========================

fastspecfit wrapper for DESI

fastspecfit-desi /global/cfs/cdirs/desi/spectro/redux/daily/tiles/80605/20201215/zbest-9-80605-20201215.fits -o photfit.fits --ntargets 2 --mp 1 --photfit

"""
import pdb # for debugging

import os, sys, time
import numpy as np
import multiprocessing

from desiutil.log import get_logger
log = get_logger()

## ridiculousness!
#if False:
#    import tempfile
#    os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
#    import matplotlib
#    matplotlib.use('Agg')

def _fastspecfit_one(args):
    """Multiprocessing wrapper."""
    return fastspecfit_one(*args)

def fastspecfit_one(iobj, data, CFit, EMFit, out, photfit=False, solve_vdisp=False):
    """Fit one spectrum."""
    log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()
    if photfit:
        cfit, _ = CFit.continuum_photfit(data)
    else:
        cfit, continuum = CFit.continuum_specfit(data, solve_vdisp=solve_vdisp)
    for col in cfit.colnames:
        out[col] = cfit[col]
    log.info('Continuum-fitting object {} took {:.2f} sec'.format(iobj, time.time()-t0))
    
    if photfit:
        return out

    # fit the emission-line spectrum
    t0 = time.time()
    emfit = EMFit.fit(data, continuum)
    for col in emfit.colnames:
        out[col] = emfit[col]
    log.info('Line-fitting object {} took {:.2f} sec'.format(iobj, time.time()-t0))
        
    return out

class DESISpectra(object):
    def __init__(self, zbestfile, exposures=False):
        """Class to read in the DESI data needed by fastspecfit.

        """
        import fitsio
        from astropy.table import Table
        #from desispec.spectra import Spectra

        self.zbestfile = zbestfile

        if exposures:
            self.specfile = zbestfile.replace('zbest-', 'spectra-')
        else:
            self.specfile = zbestfile.replace('zbest-', 'coadd-')            

        # Hack! Gather tile and targets info.
        if False:
            tilefile = '/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/sv1-tiles.fits'
            log.info('Reading {}'.format(tilefile))
            tiles = fitsio.read(tilefile)
            thistile = tiles[tiles['TILEID'] == np.int(tile)]

            hpdirname = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/targets/sv1/resolve/dark'
            log.info('Assuming hpdirname {}'.format(hpdirname))

            targets = read_targets_in_tiles(hpdirname, tiles=thistile, quick=True)
            #rr = read_targets_in_cap('/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/targets/sv1/resolve/dark', (36.448, -4.601, 1.5), quick=True)

        # Read the zbest file and fibermap and select the objects to fit.
        zb = fitsio.read(self.zbestfile, 'ZBEST')
        #fm = fitsio.read(zbestfile, 'FIBERMAP')
        fm = fitsio.read(self.specfile, 'FIBERMAP')
        
        if exposures:
            pdb.set_trace()
            _, I, _ = np.intersect1d(fm['TARGETID'], zb['TARGETID'], return_indices=True)            
            fm = fm[I]
            assert(np.all(zb['TARGETID'] == fm['TARGETID']))
        else:
            J = ((zb['Z'] > 0) * (zb['ZWARN'] == 0) * (zb['SPECTYPE'] == 'GALAXY') * 
                 (fm['OBJTYPE'] == 'TGT') * (fm['FIBERSTATUS'] == 0))
            self.fitindx = np.where(J)[0]
            #zb = zb[J]
            #fm = fm[J]

            self.zbest = Table(zb[self.fitindx])

    @staticmethod
    def desitarget_resolve_dec():
        """Default Dec cut to separate targets in BASS/MzLS from DECaLS."""
        dec = 32.375
        return dec

    def read_and_unpack(self, CFit, firsttarget=0, targetids=None, ntargets=None,
                        photfit=False, exposures=False):
        """Unpack and pre-process a single DESI spectrum.
        
        Parameters
        ----------
        specobj : :class:`desispec.spectra.Spectra`
            DESI spectra (in the standard format).
        zbest : :class:`astropy.table.Table`
            Redrock redshift table (row-aligned to `specobj`).
        CFit : :class:`fastspecfit.continuum.ContinuumFit`
            Continuum-fitting class which contains filter curves and some additional
            photometric convenience functions.
        indx : :class:`int`
            Index number (0-indexed) of the spectrum to unpack and pre-process.

        Returns
        -------
        :class:`dict` with the following keys:
            zredrock : :class:`numpy.float64`
                Redrock redshift.
            wave : :class:`list`
                Three-element list of `numpy.ndarray` wavelength vectors, one for
                each camera.    
            flux : :class:`list`    
                Three-element list of `numpy.ndarray` flux spectra, one for each
                camera and corrected for Milky Way extinction.
            ivar : :class:`list`    
                Three-element list of `numpy.ndarray` inverse variance spectra, one
                for each camera.    
            res : :class:`list`
                Three-element list of `desispec.resolution.Resolution` objects, one
                for each camera.
            snr : :class:`numpy.ndarray`
                Median per-pixel signal-to-noise ratio in the grz cameras.
            linemask : :class:`list`
                Three-element list of `numpy.ndarray` boolean emission-line masks,
                one for each camera. This mask is used during continuum-fitting.
            coadd_wave : :class:`numpy.ndarray`
                Coadded wavelength vector with all three cameras combined.
            coadd_flux : :class:`numpy.ndarray`
                Flux corresponding to `coadd_wave`.
            coadd_ivar : :class:`numpy.ndarray`
                Inverse variance corresponding to `coadd_flux`.
            photsys_south : :class:`bool`
                Boolean indicating whether this object is on the south (True) or
                north (False) photometric system based on the declination cut coded
                in `desitarget.io.desispec.resolution.Resolution`.
            phot : :class:`astropy.table.Table`
                Imaging photometry in `grzW1W2`, corrected for Milky Way extinction.
            synthphot : :class:`astropy.table.Table`
                Photometry in `grz` synthesized from the extinction-corrected
                coadded spectra (with a mild extrapolation of the data blueward and
                redward to accommodate the g-band and z-band filter curves,
                respectively.

        Notes
        -----
        Hard-coded to assume that all three cameras (grz) have spectra.

        """
        from desispec.resolution import Resolution
        from desiutil.dust import ext_odonnell
        #from desitarget.io import desitarget_resolve_dec
        from fastspecfit.util import C_LIGHT

        if photfit:
            import fitsio
            from astropy.table import Table
            fibermap = Table(fitsio.read(self.specfile, ext='FIBERMAP'))
        else:
            from desispec.io import read_spectra
            spec = read_spectra(self.specfile)#.select(targets=zb['TARGETID'])
            fibermap = spec.fibermap
            
        if targetids is None:
            targetids = self.zbest['TARGETID'].data

        if ntargets is None:
            ntargets = len(self.zbest)
        if ntargets > len(self.zbest):
            log.warning('ntargets exceeds the number of spectra.')
            raise ValueError

        targetids = targetids[firsttarget:firsttarget+ntargets]
        these = np.where([tid in targetids for tid in self.zbest['TARGETID']])[0]

        self.ntargets = len(these)
        self.zbest = self.zbest[these]
        self.fitindx = self.fitindx[these]
        self.fibermap = fibermap[self.fitindx]
        assert(np.all(self.fibermap['TARGETID'] == self.zbest['TARGETID']))

        # Unpack the subset of objects we care about into a simple dictionary.
        alldata = []
        for igal, indx in enumerate(self.fitindx):
            ra = self.fibermap['TARGET_RA'][igal]
            dec = self.fibermap['TARGET_DEC'][igal]
            ebv = CFit.SFDMap.ebv(ra, dec)

            # Unpack the data and correct for Galactic extinction. Also flag pixels that
            # may be affected by emission lines.
            data = {'zredrock': self.zbest['Z'][igal], 'photsys_south': dec < self.desitarget_resolve_dec()}

            if data['photsys_south']:
                filters = CFit.decam
                allfilters = CFit.decamwise
            else:
                filters = CFit.bassmzls
                allfilters = CFit.bassmzlswise            

            if photfit:
                # Unpack the imaging photometry.

                maggies = np.zeros(CFit.nband)
                ivarmaggies = np.zeros(CFit.nband)
                dust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(allfilters.effective_wavelengths.value, Rv=CFit.RV))

                for iband, band in enumerate(CFit.bands):
                    maggies[iband] = self.fibermap['FLUX_{}'.format(band.upper())][igal] / dust[iband]

                    ivarcol = 'FLUX_IVAR_{}'.format(band.upper())
                    if ivarcol in self.fibermap.colnames:
                        ivarmaggies[iband] = self.fibermap[ivarcol][igal] * dust[iband]**2
                    else:
                        if maggies[iband] > 0:
                            ivarmaggies[iband] = (10.0 / maggies [iband])**2 # constant S/N hack!!

                data['phot'] = CFit.parse_photometry(
                    maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
                    lambda_eff=allfilters.effective_wavelengths.value)
            else:
                data.update({'wave': [], 'flux': [], 'ivar': [], 'res': [],
                             'linemask': [], 'snr': np.zeros(3).astype('f4')})
                for icam, camera in enumerate(spec.bands):
                    dust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(spec.wave[camera], Rv=CFit.RV))       
                    data['wave'].append(spec.wave[camera])
                    data['flux'].append(spec.flux[camera][indx, :] * dust)
                    data['ivar'].append(spec.ivar[camera][indx, :] / dust**2)
                    data['res'].append(Resolution(spec.resolution_data[camera][indx, :, :]))
                    data['snr'][icam] = np.median(spec.flux[camera][indx, :] * np.sqrt(spec.ivar[camera][indx, :]))

                    linemask = np.ones_like(spec.wave[camera], bool)
                    for line in CFit.linetable:
                        zwave = line['restwave'] * (1 + data['zredrock'])
                        I = np.where((spec.wave[camera] >= (zwave - 1.5*CFit.linemask_sigma * zwave / C_LIGHT)) *
                                     (spec.wave[camera] <= (zwave + 1.5*CFit.linemask_sigma * zwave / C_LIGHT)))[0]
                        if len(I) > 0:
                            linemask[I] = False
                    data['linemask'].append(linemask)

                # Synthesize photometry from a quick coadded (inverse-variance
                # weighted) spectrum, being careful to resolve between the north and
                # south.
                coadd_wave = np.unique(np.hstack(data['wave']))
                coadd_flux3d = np.zeros((len(coadd_wave), 3))
                coadd_ivar3d = np.zeros_like(coadd_flux3d)
                for icam in np.arange(len(spec.bands)):
                    I = np.where(np.isin(data['wave'][icam], coadd_wave))[0]
                    J = np.where(np.isin(coadd_wave, data['wave'][icam]))[0]
                    coadd_flux3d[J, icam] = data['flux'][icam][I]
                    coadd_ivar3d[J, icam] = data['ivar'][icam][I]

                coadd_ivar = np.sum(coadd_ivar3d, axis=1)
                coadd_flux = np.zeros_like(coadd_ivar)
                good = np.where(coadd_ivar > 0)[0]
                coadd_flux[good] = np.sum(coadd_ivar3d[good, :] * coadd_flux3d[good, :], axis=1) / coadd_ivar[good]

                #data.update({'coadd_wave': coadd_wave, 'coadd_flux': coadd_flux, 'coadd_ivar': coadd_ivar})
                #del coadd_wave, coadd_ivar, coadd_flux

                padflux, padwave = filters.pad_spectrum(coadd_flux, coadd_wave, method='edge')
                synthmaggies = filters.get_ab_maggies(padflux / CFit.fluxnorm, padwave)
                synthmaggies = synthmaggies.as_array().view('f8')

                # code to synthesize uncertainties from the variance spectrum
                #var, mask = _ivar2var(data['coadd_ivar'])
                #padvar, padwave = filters.pad_spectrum(var[mask], data['coadd_wave'][mask], method='edge')
                #synthvarmaggies = filters.get_ab_maggies(1e-17**2 * padvar, padwave)
                #synthivarmaggies = 1 / synthvarmaggies.as_array().view('f8')[:3] # keep just grz
                #
                #data['synthphot'] = CFit.parse_photometry(
                #    maggies=synthmaggies, lambda_eff=lambda_eff[:3],
                #    ivarmaggies=synthivarmaggies, nanomaggies=False)

                data['synthphot'] = CFit.parse_photometry(
                    maggies=synthmaggies, nanomaggies=False,
                    lambda_eff=filters.effective_wavelengths.value)

            alldata.append(data)

        return alldata

def init_output(tile, night, zbest, fibermap, CFit, EMFit, photfit=False):
    """Initialize the fastspecfit output data table.

    Parameters
    ----------
    tile : :class:`str`
        Tile number.
    night : :class:`str`
        Night on which `tile` was observed.
    zbest : :class:`astropy.table.Table`
        Redrock redshift table (row-aligned to `fibermap`).
    fibermap : :class:`astropy.table.Table`
        Fiber map (row-aligned to `zbest`).
    CFit : :class:`fastspecfit.continuum.ContinuumFit`
        Continuum-fitting class.
    EMFit : :class:`fastspecfit.emlines.EMLineFit`
        Emission-line fitting class.

    Returns
    -------
    

    Notes
    -----

    """
    from astropy.table import Table, Column, hstack

    # Grab info on the emission lines and the continuum.
    nobj = len(zbest)

    out = Table()
    for zbestcol in ['TARGETID', 'Z']:#, 'ZERR']:#, 'SPECTYPE', 'DELTACHI2']
        out[zbestcol] = zbest[zbestcol]
    for fibermapcol in ['FIBER']:
        out[fibermapcol] = fibermap[fibermapcol]
    out.add_column(Column(name='NIGHT', data=np.repeat(night, nobj)), index=0)
    out.add_column(Column(name='TILE', data=np.repeat(tile, nobj)), index=0)

    if photfit:
        out = hstack((out, CFit.init_phot_output(nobj)))
    else:
        out = hstack((out, CFit.init_spec_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

    return out

def parse(options=None):
    """Parse input arguments.

    """
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ## required but with sensible defaults
    #parser.add_argument('--night', default='20200225', type=str, help='Night to process.')
    #parser.add_argument('--tile', default='70502', type=str, help='Tile number to process.')

    # optional inputs
    parser.add_argument('-o', '--outfile', type=str, default='.', help='Full path to output FITS file.')

    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file (0-indexed).')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of target IDs to process.')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing processes per MPI rank or node.')

    #parser.add_argument('--specprod', type=str, default='daily', choices=['andes', 'daily', 'variance-model'],
    #                    help='Spectroscopic production to process.')

    parser.add_argument('--exposures', action='store_true', help='Fit the individual exposures (not the coadds).')
    parser.add_argument('--photfit', action='store_true', help='Fit and write out just the broadband photometry.')
    parser.add_argument('--solve-vdisp', action='store_true', help='Solve for the velocity disperion.')
    parser.add_argument('--verbose', action='store_true', help='Be (more) verbose.')

    #parser.add_argument('--use-vi', action='store_true', help='Select spectra with high-quality visual inspections (VI).')
    #parser.add_argument('--overwrite', action='store_true', help='Overwrite any existing files.')
    #parser.add_argument('--mpi', action='store_true', help='Parallelize using MPI.')
    #parser.add_argument('--no-write-spectra', dest='write_spectra', default=True, action='store_false',
    #                    help='Do not write out the selected spectra for the specified tile and night.')

    parser.add_argument('zbestfiles', nargs='*', help='Input zbest file(s).')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspecfit {}'.format(' '.join(options)))

    return args

def main(args=None, comm=None):
    """Main module.

    """
    from astropy.table import vstack
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    ## Create directories as needed.
    #fastspecfit_dir = os.getenv('FASTSPECFIT_DATA')
    #for subdir in ('spectra', 'results', 'qa'):
    #    outdir = os.path.join(fastspecfit_dir, subdir, args.specprod)
    #    if not os.path.isdir(outdir):
    #        if args.verbose:
    #            log.debug('Creating directory {}'.format(outdir))
    #        os.makedirs(outdir, exist_ok=True)
    
    ## If the output file exists, we're done!
    #resultsdir = os.path.join(fastspecfit_dir, 'results', args.specprod)
    #if args.photfit:
    #    prefix = 'photfit'
    #else:
    #    prefix = 'specfit'
        
    #outfile = os.path.join(resultsdir, '{}-{}-{}.fits'.format(
    #    prefix, args.tile, args.night))
    #if os.path.isfile(outfile) and not args.overwrite:
    #    log.info('Output file {} exists; all done!'.format(outfile))
    #    return

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum- and emission-line fitting classes.
    t0 = time.time()
    CFit = ContinuumFit(verbose=args.verbose)
    EMFit = EMLineFit(verbose=args.verbose)
    #out = init_output(args.tile, args.night, zbest, specobj.fibermap, CFit, EMFit, photfit=args.photfit)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Iteratively read the data.
    for zbestfile in args.zbestfiles:
        t0 = time.time()
        DESISpec = DESISpectra(zbestfile, exposures=args.exposures)

        data = DESISpec.read_and_unpack(CFit, firsttarget=args.firsttarget, targetids=targetids,
                                        ntargets=args.ntargets, exposures=args.exposures,
                                        photfit=args.photfit)
        print('Hacked output!!')
        out = init_output('80605', '20201215', DESISpec.zbest, DESISpec.fibermap, CFit, EMFit, photfit=args.photfit)
        log.info('Reading and unpacking the spectra to be fitted took: {:.2f} sec'.format(time.time()-t0))

    # Fit in parallel
    fitargs = [(iobj, data[iobj], CFit, EMFit, out[iobj], args.photfit, args.solve_vdisp)
               for iobj in np.arange(DESISpec.ntargets)]
    if args.mp > 1:
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastspecfit_one, fitargs)
    else:
        _out = [fastspecfit_one(*_fitargs) for _fitargs in fitargs]
    out = vstack(_out)
    del _out

    # write out
    t0 = time.time()
    log.info('Writing results for {} objects to {}'.format(len(out), args.outfile))
    out.write(args.outfile, overwrite=True)
    log.info('Writing out took {:.2f} sec'.format(time.time()-t0))
    
