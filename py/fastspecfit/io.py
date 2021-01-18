#!/usr/bin/env python
"""
fastspecfit.io
==============

Tools for reading and writing.

"""
import pdb # for debugging

import os, time
import numpy as np
import fitsio
from astropy.table import Table

from desiutil.log import get_logger
log = get_logger()

class DESISpectra(object):
    def __init__(self):
        """Class to read in the DESI data needed by fastspecfit.

        """
        #from desitarget.geomask import pixarea2nside

        self.fiberassign_dir = os.path.join(os.getenv('DESI_ROOT'), 'target', 'fiberassign', 'tiles', 'trunk')

        # read the summary tile file (not sure we need this...)
        tilefile = '/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/sv1-tiles.fits'
        self.tiles = fitsio.read(tilefile)
        log.info('Read {} tiles from {}'.format(len(self.tiles), tilefile))
        
        ## Closest nside to DESI tile area of ~7 deg (from
        ## desitarget.io.read_targets_in_quick).
        #self.tile_nside = pixarea2nside(8.)

    def _get_targetdir(self, tileid):
        """Get the targets catalog used to build a given fiberassign catalog.

        """
        thistile = self.tiles[self.tiles['TILEID'] == tileid]
        stileid = '{:06d}'.format(tileid)
        fahdr = fitsio.read_header(os.path.join(self.fiberassign_dir, stileid[:3], 'fiberassign-{}.fits.gz'.format(stileid)))
        targetdir = fahdr['TARG']

        # sometimes this is a KPNO directory!
        if not os.path.isdir(targetdir):
            log.warning('Targets directory not found {}'.format(targetdir))
            targetdir = os.path.join(os.getenv('DESI_ROOT'), targetdir.replace('/data/', ''))
            if not os.path.isdir(targetdir):
                log.fatal('Targets directory not found {}'.format(targetdir))
                raise IOError

        return targetdir

    def find_specfiles(self, zbestfiles=None, fastfit=None, specprod=None,
                       firsttarget=0, targetids=None, ntargets=None, exposures=False):
        """Initialize the fastspecfit output data table.

        Parameters
        ----------

        Returns
        -------

        Notes
        -----

        fastfit - results table (overrides zbestfiles and then specprod is required

        """
        from glob import glob
        from desimodel.footprint import radec2pix

        if zbestfiles is None and fastfit is None:
            log.warning('At least one of zbestfiles or fastfit (and specprod) are required.')
            raise IOError

        if fastfit is not None:
            if specprod is None:
                log.warning('specprod is required when passing a fastfit results table!')
                raise IOError
            fastfit = Table(fastfit)
            targetids = fastfit['TARGETID'].data
            petals = fastfit['FIBER'].data // 500
            tiles = fastfit['TILEID'].astype(str).data
            nights = fastfit['NIGHT'].astype(str).data

            log.info('keep track of a sorting index which will keep the zbest and fibermap tables, below, row-aligned with the input fastfit table')
            pdb.set_trace()

            zbestfiles = []
            for petal, tile, night in zip(petals, tiles, nights):
                if petal == 2 or petal == 6:
                    log.warning('Temporarily skipping petals 2 & 6!')
                    continue
                zbestfile = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', specprod, 'tiles',
                                         tile, night, 'zbest-{}-{}-{}.fits'.format(petal, tile, night))
                if not os.path.isfile(zbestfile):
                    log.warning('zbestfile {} not found.'.format(zbestfile))
                    raise IOError
                zbestfiles.append(zbestfile)
            zbestfiles = np.array(zbestfiles)

        if len(np.atleast_1d(zbestfiles)) == 0:
            log.warning('You must provide at least one zbestfile!')
            raise IOError

        # Try to glean specprod so we can write it to the output file. This
        # should really be in the file header--- see
        # https://github.com/desihub/desispec/issues/1077
        if specprod is None:
            # stupidly fragile!
            specprod = np.atleast_1d(zbestfiles)[0].replace(os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux'), '').split(os.sep)[1]
            self.specprod = specprod
            log.info('Parsed specprod={}'.format(specprod))
            
        # Should we not sort...?
        #zbestfiles = np.array(set(np.atleast_1d(zbestfiles)))
        zbestfiles = np.array(sorted(set(np.atleast_1d(zbestfiles))))
        log.info('Reading and parsing {} unique zbestfile(s)'.format(len(zbestfiles)))

        #self.fitindx = []
        self.zbest, self.fibermap = [], []
        self.zbestfiles, self.specfiles = [], []
        for zbestfile in np.atleast_1d(zbestfiles):
            if exposures:
                specfile = zbestfile.replace('zbest-', 'spectra-')
            else:
                specfile = zbestfile.replace('zbest-', 'coadd-')

            # If targetids is *not* given we have to choose "good" objects
            # before subselecting (e.g., we don't want sky spectra).
            if targetids is None:
                zb = fitsio.read(zbestfile, 'ZBEST')
                fm = fitsio.read(specfile, 'FIBERMAP')
                #fm = fitsio.read(zbestfile, 'FIBERMAP')
                #hdr = fitsio.read_header(zbestfile, 'FIBERMAP')

                if exposures:
                    log.fatal('Not yet implemented!')
                    _, I, _ = np.intersect1d(fm['TARGETID'], zb['TARGETID'], return_indices=True)            
                    fm = fm[I]
                    assert(np.all(zb['TARGETID'] == fm['TARGETID']))
                else:
                    # The fibermap includes all the spectra that went into the coadd
                    # (based on the unique TARGETID which is in the zbest table).
                    assert(len(zb) == len(fm))
                    fitindx = np.where((zb['Z'] > 0) * (zb['ZWARN'] == 0) * (zb['SPECTYPE'] == 'GALAXY') * 
                         (fm['OBJTYPE'] == 'TGT') * (fm['FIBERSTATUS'] == 0))[0]
                    #self.fitindx.append(fitindx)
            else:
                # We already know we like the input targetids, so no selection
                # needed.
                alltargetids = fitsio.read(zbestfile, 'ZBEST', columns='TARGETID')                
                fitindx = np.where([tid in targetids for tid in alltargetids])[0]
                
            if len(fitindx) == 0:
                log.info('No requested targets found in zbestfile {}'.format(zbestfile))
                continue

            # Do we want just a subset of the available objects?
            if ntargets is None:
                _ntargets = len(fitindx)
            else:
                _ntargets = ntargets
            if _ntargets > len(fitindx):
                log.warning('Number of requested ntargets exceeds the number of targets on {}; reading all of them.'.format(
                    zbestfile))
                #raise ValueError

            __ntargets = len(fitindx)
            fitindx = fitindx[firsttarget:firsttarget+_ntargets]
            if len(fitindx) == 0:
                log.info('All {} targets in zbestfile {} have been dropped (firsttarget={}, ntargets={}).'.format(
                    __ntargets, zbestfile, firsttarget, _ntargets))
                continue
                
            # If firsttarget is a large index then the set can become empty.
            if targetids is None:
                zb = Table(zb[fitindx])
                fm = Table(fm[fitindx])
            else:
                zb = Table(fitsio.read(zbestfile, 'ZBEST', rows=fitindx))
                fm = Table(fitsio.read(specfile, 'FIBERMAP', rows=fitindx))
            assert(np.all(zb['TARGETID'] == fm['TARGETID']))

            self.zbest.append(Table(zb))
            self.fibermap.append(Table(fm))
            self.zbestfiles.append(zbestfile)
            self.specfiles.append(specfile)

        if len(self.zbest) == 0:
            log.warning('No targets read!')
            return

        # The targets catalogs are organized on a much larger spatial scale, so
        # grab additional targeting info outside the main loop. See
        # desitarget.io.read_targets_in_tiles for the algorithm.
        t0 = time.time()
        alltileid = [fm['TILEID'][0] for fm in self.fibermap]
        info = Table(np.hstack([fm['TARGETID', 'TARGET_RA', 'TARGET_DEC'] for fm in self.fibermap]))
        targets = []
        for tileid in set(alltileid):
            targetdir = self._get_targetdir(tileid)
            targetfiles = glob(os.path.join(targetdir, '*-hp-*.fits'))
            filenside = fitsio.read_header(targetfiles[0], ext=1)['FILENSID']
            pixlist = radec2pix(filenside, info['TARGET_RA'], info['TARGET_DEC'])
            targetfiles = [targetfiles[0].split('hp-')[0]+'hp-{}.fits'.format(pix) for pix in set(pixlist)]
            for targetfile in targetfiles:
                alltargets = fitsio.read(targetfile, columns=['TARGETID', 'FLUX_W1', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
                                                              'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
                                                              'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2'])
                alltargets = alltargets[np.isin(alltargets['TARGETID'], info['TARGETID'])]
                targets.append(alltargets)
        targets = Table(np.hstack(targets))

        fmaps = []
        for fm in self.fibermap:
            srt = np.hstack([np.where(tid == targets['TARGETID'])[0] for tid in fm['TARGETID']])
            for col in ['FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R',
                        'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']:
                if col not in fm.colnames:
                    fm[col] = targets[col][srt]
            assert(np.all(fm['TARGETID'] == targets['TARGETID'][srt]))
            fmaps.append(Table(fm))
        log.info('Read and parsed targeting info for {} objects in {:.2f} sec'.format(len(targets), time.time()-t0))

        self.fibermap = fmaps # update

    @staticmethod
    def desitarget_resolve_dec():
        """Default Dec cut to separate targets in BASS/MzLS from DECaLS."""
        dec = 32.375
        return dec

    def read_and_unpack(self, CFit, fastphot=False, exposures=False,
                        synthphot=True, remember_coadd=False):
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

        # Read everything into a simple dictionary.
        t0 = time.time()
        alldata = []
        for specfile, zbest, fibermap in zip(self.specfiles, self.zbest, self.fibermap):
            log.info('Reading {} spectra from {}'.format(len(zbest), specfile))

            # sometimes these are an astropy.table.Row!
            zbest = Table(zbest)
            fibermap = Table(fibermap)

            if not fastphot:
                spec = read_spectra(specfile).select(targets=zbest['TARGETID'])
                assert(np.all(spec.fibermap['TARGETID'] == zbest['TARGETID']))
                assert(np.all(spec.fibermap['TARGETID'] == fibermap['TARGETID']))
            
            for igal in np.arange(len(zbest)):
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
                if fastphot:
                    # all photometry
                    #fibermap['MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']
                    mw_transmission_flux = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(allfilters.effective_wavelengths.value, Rv=CFit.RV))

                    maggies = np.zeros(len(CFit.bands))
                    ivarmaggies = np.zeros(len(CFit.bands))
                    for iband, band in enumerate(CFit.bands):
                        maggies[iband] = fibermap['FLUX_{}'.format(band.upper())][igal] / mw_transmission_flux[iband]
                        ivarmaggies[iband] = fibermap['FLUX_IVAR_{}'.format(band.upper())][igal] * mw_transmission_flux[iband]**2

                    data['phot'] = CFit.parse_photometry(CFit.bands,
                        maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
                        lambda_eff=allfilters.effective_wavelengths.value)

                    # fiber fluxes
                    mw_transmission_fiberflux = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(filters.effective_wavelengths.value, Rv=CFit.RV))

                    fibermaggies = np.zeros(len(CFit.fiber_bands))
                    #ivarfibermaggies = np.zeros(len(CFit.fiber_bands))
                    for iband, band in enumerate(CFit.fiber_bands):
                        fibermaggies[iband] = fibermap['FIBERTOTFLUX_{}'.format(band.upper())][igal] / mw_transmission_fiberflux[iband]
                        #ivarfibermaggies[iband] = fibermap['FIBERTOTFLUX_IVAR_{}'.format(band.upper())][igal] * mw_transmission_fiberflux[iband]**2

                    data['fiberphot'] = CFit.parse_photometry(CFit.fiber_bands,
                        maggies=fibermaggies, nanomaggies=True,
                        lambda_eff=filters.effective_wavelengths.value)
                    
                else:
                    
                    data.update({'wave': [], 'flux': [], 'ivar': [], 'res': [],
                                 'linemask': [], 'snr': np.zeros(3).astype('f4')})
                    for icam, camera in enumerate(spec.bands):
                        mw_transmission_spec = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(spec.wave[camera], Rv=CFit.RV))       
                        data['wave'].append(spec.wave[camera])
                        data['flux'].append(spec.flux[camera][igal, :] / mw_transmission_spec)
                        data['ivar'].append(spec.ivar[camera][igal, :] * mw_transmission_spec**2)
                        data['res'].append(Resolution(spec.resolution_data[camera][igal, :, :]))
                        data['snr'][icam] = np.median(spec.flux[camera][igal, :] * np.sqrt(spec.ivar[camera][igal, :]))

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

                    #import matplotlib.pyplot as plt
                    #plt.clf()
                    #for ii in np.arange(3):
                    #    plt.plot(data['wave'][ii], data['flux'][ii])
                    #plt.plot(coadd_wave, coadd_flux-2, alpha=0.6, color='k')
                    #plt.xlim(5500, 6000)
                    #plt.savefig('test.png')
                    #pdb.set_trace()

                    if remember_coadd:
                        data.update({'coadd_wave': coadd_wave, 'coadd_flux': coadd_flux, 'coadd_ivar': coadd_ivar})

                    if synthphot:
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

        self.zbest = Table(np.hstack(self.zbest))
        self.fibermap = Table(np.hstack(self.fibermap))
        self.ntargets = len(self.zbest)
        log.info('Read data for {} objects in {:.2f} sec'.format(self.ntargets, time.time()-t0))
        
        return alldata
        
    def init_output(self, CFit=None, EMFit=None, fastphot=False):
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
        #for fibermapcol, outcol in zip(['TARGETID', 'TARGET_RA', 'TARGET_DEC'],
        #                               ['TARGETID', 'RA', 'DEC']):
        #    out[outcol] = self.fibermap[fibermapcol]
        for fibermapcol in ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'NIGHT', 'TILEID', 'FIBER', 'EXPID']:
            out[fibermapcol] = self.fibermap[fibermapcol]
        for zbestcol in ['Z', 'DELTACHI2']:#, 'ZERR']:#, 'SPECTYPE', ]
            out[zbestcol] = self.zbest[zbestcol]
        out['PHOTSYS_SOUTH'] = out['TARGET_DEC'] < self.desitarget_resolve_dec()

        # target column names (e.g., DESI_TARGET)
        _targetcols = ['DESI_TARGET' in col or 'BGS_TARGET' in col or 'MWS_TARGET' in col for col in self.fibermap.colnames]
        targetcols = np.array(self.fibermap.colnames)[_targetcols]
        for targetcol in targetcols:
            out[targetcol] = self.fibermap[targetcol]
        
        out['TARGET_RA'].unit = u.deg
        out['TARGET_DEC'].unit = u.deg
        #out['RA'].unit = u.deg
        #out['DEC'].unit = u.deg

        if fastphot:
            out = hstack((out, CFit.init_phot_output(nobj)))
        else:
            out = hstack((out, CFit.init_spec_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

        return out

def write_fastspec(out, outfile=None, specprod=None, fastphot=False):
    """Write out.

    """
    from astropy.io import fits
    
    t0 = time.time()
    outdir = os.path.dirname(os.path.abspath(outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    log.info('Writing results for {} objects to {}'.format(len(out), outfile))
    #out.write(outfile, overwrite=True)
    hduprim = fits.PrimaryHDU()
    hduout = fits.convenience.table_to_hdu(out)
    if fastphot:
        hduout.header['EXTNAME'] = 'FASTPHOT'
    else:
        hduout.header['EXTNAME'] = 'FASTSPEC'
        
    if specprod:
        hduout.header['SPECPROD'] = specprod
    hx = fits.HDUList([hduprim, hduout])
    hx.writeto(outfile, overwrite=True)
    log.info('Writing out took {:.2f} sec'.format(time.time()-t0))
