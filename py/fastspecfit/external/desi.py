#!/usr/bin/env python
"""
fastspecfit.external.desi
=========================

fastspecfit wrapper for DESI

fastspecfit-desi /global/cfs/cdirs/desi/spectro/redux/daily/tiles/80605/20201215/zbest-9-80605-20201215.fits -o photfit.fits --ntargets 2 --mp 1 --photfit

"""
import pdb # for debugging

import os, time
import numpy as np
import fitsio
from astropy.table import Table

from desiutil.log import get_logger
log = get_logger()

# ridiculousness!
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

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
    def __init__(self, zbestfiles, exposures=False):
        """Class to read in the DESI data needed by fastspecfit.

        """
        #from glob import glob
        from astropy.table import vstack
        from desitarget.io import read_targets_in_tiles

        zbestfiles = np.array(sorted(zbestfiles))
        #zbestfiles = np.array(sorted(glob(os.path.join(zbestdir, 'zbest-?-*-????????.fits'))))
        if len(zbestfiles) == 0:
            log.warning('No zbest-?-*-????????.fits files found in {}'.format(zbestdir))
            raise IOError
        
        self.zbestfiles = zbestfiles
        
        if exposures:
            self.specfiles = [zbestfile.replace('zbest-', 'spectra-') for zbestfile in zbestfiles]
        else:
            self.specfiles = [zbestfile.replace('zbest-', 'coadd-') for zbestfile in zbestfiles]
        
        zbest, fmap = [], []
        self.fitindx = []
        for ifile, (zbestfile, specfile) in enumerate(zip(self.zbestfiles, self.specfiles)):
            # Read the zbest file and fibermap and select the objects to fit.
            zb = fitsio.read(zbestfile, 'ZBEST')
            #fm = fitsio.read(zbestfile, 'FIBERMAP')

            hdr = fitsio.read_header(zbestfile, 'FIBERMAP')
            fm = fitsio.read(specfile, 'FIBERMAP')

            if exposures:
                pdb.set_trace()
                _, I, _ = np.intersect1d(fm['TARGETID'], zb['TARGETID'], return_indices=True)            
                fm = fm[I]
                assert(np.all(zb['TARGETID'] == fm['TARGETID']))
            else:
                J = ((zb['Z'] > 0) * (zb['ZWARN'] == 0) * (zb['SPECTYPE'] == 'GALAXY') * 
                     (fm['OBJTYPE'] == 'TGT') * (fm['FIBERSTATUS'] == 0))
                fitindx = np.where(J)[0]
                self.fitindx.append(fitindx)

                zb = Table(zb[fitindx])
                fm = Table(fm[fitindx])

            # Add additional columns to the fibermap table. Should we pass this
            # forward from mpi-fastspecfit?
            if ifile == 0:
                tileid = fm['TILEID'][0]
                
                tilefile = '/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/sv1-tiles.fits'
                log.info('Reading {}'.format(tilefile))
                tiles = fitsio.read(tilefile)
                thistile = tiles[tiles['TILEID'] == tileid]

                stileid = '{:06d}'.format(tileid)
                fahdr = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/{}/fiberassign-{}.fits.gz'.format(stileid[:3], stileid))
                hpdirname = fahdr['TARG']

                # sometimes this is a KPNO directory!
                if not os.path.isdir(hpdirname):
                    log.info('Targets directory not found {}'.format(hpdirname))
                    hpdirname = os.path.join(os.getenv('DESI_ROOT'), hpdirname.replace('/data/', ''))
                log.info('Reading targets from {}'.format(hpdirname))
                alltargets = read_targets_in_tiles(hpdirname, tiles=thistile, quick=True)

            keep = np.hstack([np.where(tid == alltargets['TARGETID'])[0] for tid in fm['TARGETID']])
            targets = alltargets[keep]
            assert(np.all(targets['TARGETID'] == fm['TARGETID']))
            assert(np.all(targets['FLUX_W1'] == fm['FLUX_W1']))
            
            for col in ['FLUX_IVAR_W1', 'FLUX_IVAR_W2']:
                fm[col] = targets[col]
                
            zbest.append(zb)
            fmap.append(fm)

        self.zbest = zbest
        self.fibermap = fmap
        #self.zbest = vstack(zbest)
        #self.fibermap = vstack(fmap)
        #self.zbest = Table(np.hstack(zbest))

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
        from desispec.io import read_spectra
        #from desitarget.io import desitarget_resolve_dec
        from fastspecfit.util import C_LIGHT

        if targetids is None:
            thesetargetids = None
        else:
            thesetargetids = targetids.copy()

        # Unpack the subset of objects we care about into a simple dictionary.
        alldata, zbests, fibermaps = [], [], []
        for specfile, zbest, fibermap, fitindx in zip(self.specfiles, self.zbest, self.fibermap, self.fitindx):

            # Select the subset of objects to process.
            if thesetargetids is None:
                targetids = zbest['TARGETID'].data
            else:
                targetids = thesetargetids.copy()

            these = np.where([tid in targetids for tid in zbest['TARGETID']])[0]
            if len(these) == 0:
                continue

            if ntargets is None:
                ntargets = len(these)

            if ntargets > len(these):
                log.warning('Number of requested ntargets exceeds the number of spectra on {}.'.format(specfile))
                #raise ValueError
            else:
                # the logic here is not quite right...needs more testing
                these = these[firsttarget:firsttarget+ntargets]
            
            fitindx = fitindx[these]
            zbest = zbest[these]
            fibermap = fibermap[these]
            zbests.append(zbest)
            fibermaps.append(fibermap)

            if not photfit:
                spec = read_spectra(specfile)#.select(targets=zbest['TARGETID'])
                
            for igal, indx in enumerate(fitindx):
                ra = fibermap['TARGET_RA'][igal]
                dec = fibermap['TARGET_DEC'][igal]
                ebv = CFit.SFDMap.ebv(ra, dec)

                # Unpack the data and correct for Galactic extinction. Also flag pixels that
                # may be affected by emission lines.
                data = {'zredrock': zbest['Z'][igal], 'photsys_south': dec < self.desitarget_resolve_dec()}

                if data['photsys_south']:
                    filters = CFit.decam
                    allfilters = CFit.decamwise
                else:
                    filters = CFit.bassmzls
                    allfilters = CFit.bassmzlswise            

                # Unpack the imaging photometry and correct for MW dust.
                if photfit:
                    # all photometry
                    dust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(allfilters.effective_wavelengths.value, Rv=CFit.RV))

                    maggies = np.zeros(len(CFit.bands))
                    ivarmaggies = np.zeros(len(CFit.bands))
                    for iband, band in enumerate(CFit.bands):
                        maggies[iband] = fibermap['FLUX_{}'.format(band.upper())][igal] / dust[iband]
                        ivarmaggies[iband] = fibermap['FLUX_IVAR_{}'.format(band.upper())][igal] * dust[iband]**2

                    data['phot'] = CFit.parse_photometry(CFit.bands,
                        maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
                        lambda_eff=allfilters.effective_wavelengths.value)

                    # fiber fluxes
                    fiberdust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(filters.effective_wavelengths.value, Rv=CFit.RV))

                    fibermaggies = np.zeros(len(CFit.fiber_bands))
                    #ivarfibermaggies = np.zeros(len(CFit.fiber_bands))
                    for iband, band in enumerate(CFit.fiber_bands):
                        fibermaggies[iband] = fibermap['FIBERTOTFLUX_{}'.format(band.upper())][igal] / fiberdust[iband]
                        #ivarfibermaggies[iband] = fibermap['FIBERTOTFLUX_IVAR_{}'.format(band.upper())][igal] * dust[iband]**2

                    data['fiberphot'] = CFit.parse_photometry(CFit.fiber_bands,
                        maggies=fibermaggies, nanomaggies=True,
                        lambda_eff=filters.effective_wavelengths.value)
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
                    #data['synthphot'] = CFit.parse_photometry(CFit.bands,
                    #    maggies=synthmaggies, lambda_eff=lambda_eff[:3],
                    #    ivarmaggies=synthivarmaggies, nanomaggies=False)

                    data['synthphot'] = CFit.parse_photometry(CFit.synth_bands,
                        maggies=synthmaggies, nanomaggies=False,
                        lambda_eff=filters.effective_wavelengths.value)

                alldata.append(data)

        self.zbest = Table(np.hstack(zbests))
        self.fibermap = Table(np.hstack(fibermaps))
        self.ntargets = len(self.zbest)

        return alldata
        
    def init_output(self, CFit, EMFit, photfit=False):
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
        import astropy.units as u
        from astropy.table import hstack

        # Grab info on the emission lines and the continuum.
        nobj = len(self.zbest)

        out = Table()
        for fibermapcol, outcol in zip(['TARGETID', 'TARGET_RA', 'TARGET_DEC'],
                                       ['TARGETID', 'RA', 'DEC']):
            out[outcol] = self.fibermap[fibermapcol]
        for fibermapcol in ['NIGHT', 'TILEID', 'FIBER', 'EXPID']:
            out[fibermapcol] = self.fibermap[fibermapcol]
        for zbestcol in ['Z', 'DELTACHI2']:#, 'ZERR']:#, 'SPECTYPE', ]
            out[zbestcol] = self.zbest[zbestcol]

        # target column names (e.g., DESI_TARGET)
        _targetcols = ['DESI_TARGET' in col or 'BGS_TARGET' in col or 'MWS_TARGET' in col for col in self.fibermap.colnames]
        targetcols = np.array(self.fibermap.colnames)[_targetcols]
        for targetcol in targetcols:
            out[targetcol] = self.fibermap[targetcol]
        
        out['RA'].unit = u.deg
        out['DEC'].unit = u.deg

        if photfit:
            out = hstack((out, CFit.init_phot_output(nobj)))
        else:
            out = hstack((out, CFit.init_spec_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

        return out

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

    parser.add_argument('--exposures', action='store_true', help='Fit the individual exposures (not the coadds).')

    parser.add_argument('--qa', action='store_true', help='Build QA (skips fitting).')
    parser.add_argument('--photfit', action='store_true', help='Fit and write out just the broadband photometry.')
    parser.add_argument('--solve-vdisp', action='store_true', help='Solve for the velocity disperion.')

    parser.add_argument('zbestfiles', nargs='*', help='Full path to input zbest file(s).')

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
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit

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
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    DESISpec = DESISpectra(args.zbestfiles, exposures=args.exposures)

    data = DESISpec.read_and_unpack(CFit, firsttarget=args.firsttarget, targetids=targetids,
                                    ntargets=args.ntargets, exposures=args.exposures,
                                    photfit=args.photfit)
    out = DESISpec.init_output(CFit, EMFit, photfit=args.photfit)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        DESISpec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], CFit, EMFit, out[iobj], args.photfit, args.solve_vdisp)
               for iobj in np.arange(DESISpec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastspecfit_one, fitargs)
    else:
        _out = [fastspecfit_one(*_fitargs) for _fitargs in fitargs]
    out = Table(np.hstack(_out))
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # write out
    t0 = time.time()
    outdir = os.path.dirname(os.path.abspath(args.outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
        
    log.info('Writing results for {} objects to {}'.format(len(out), args.outfile))
    out.write(args.outfile, overwrite=True)
    log.info('Writing out took {:.2f} sec'.format(time.time()-t0))
