def read_spectra(tile, night, specprod='andes', use_vi=False, write_spectra=True,
                 overwrite=False, verbose=False):
    """Read the spectra and redshift catalog for a given tile and night.

    Parameters
    ----------
    tile : :class:`str`
        Tile number to analyze.
    night : :class:`str`
        Night on which `tile` was observed.
    specprod : :class:`str`, defaults to `andes`.
        Spectroscopic production to read.
    use_vi : :class:`bool`, optional, defaults to False
        Select the subset of spectra with high-quality visual inspections.
    write_spectra : :class:`bool`, optional, defaults to True
        Write out the selected spectra (useful for testing.)
    overwrite : :class:`bool`, optional, defaults to False
        Overwrite output files if they exist on-disk.
    verbose : :class:`bool`, optional, defaults to False
        Trigger more verbose output.

    Returns
    -------
    :class:`astropy.table.Table`
        Redrock redshift table for the given tile and night.
    :class:`desispec.spectra.Spectra`
        DESI spectra for the given tile and night (in the standard format).

    Notes
    -----
    The spectra from all 10 spectrographs are combined and only the subset of
    galaxy spectra with good redshifts (and, optionally, high-quality visual
    inspections) are returned.

    """
    import fitsio
    from astropy.table import Table, vstack
    from desispec.spectra import Spectra
    import desispec.io
    from desitarget.io import read_targets_in_tiles

    if False:
        hpdirname = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/targets/sv1/resolve/dark'
        print('Assuming hpdirname={}'.format(hpdirname))
        from desimodel.io import load_tiles
        tileinfo = load_tiles()
        tileinfo = tileinfo[tileinfo['TILEID'] == tile]

    data_dir = os.path.join(os.getenv('FASTSPECFIT_DATA'), 'spectra', specprod)
    zbestoutfile = os.path.join(data_dir, 'zbest-{}-{}.fits'.format(tile, night))
    coaddoutfile = os.path.join(data_dir, 'coadd-{}-{}.fits'.format(tile, night))
    if os.path.isfile(coaddoutfile) and not overwrite:
        zbest = Table(fitsio.read(zbestoutfile))
        log.info('Read {} redshifts from {}'.format(len(zbest), zbestoutfile))

        coadd = desispec.io.read_spectra(coaddoutfile)
        log.info('Read {} spectra from {}'.format(len(zbest), coaddoutfile))

        assert(np.all(zbest['TARGETID'] == coadd.fibermap['TARGETID']))
        return zbest, coadd

    log.info('Parsing tile, night {}, {} from specprod {}'.format(tile, night, specprod))

    specprod_dir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', specprod)
    desidatadir = os.path.join(specprod_dir, 'tiles', tile, night)

    # use visual inspection results?
    if use_vi:
        truthdir = os.path.join(os.getenv('DESI_ROOT'), 'sv', 'vi', 'TruthTables')
        tile2file = {
            '66003': os.path.join(truthdir, 'truth_table_BGS_v1.2.csv'),
            '70500': os.path.join(truthdir, 'truth_table_ELG_v1.2_latest.csv'),
            }
        if tile not in tile2file.keys():
            log.warning('No (known) truth file exists for tile {}. Setting use_vi=False.'.format(tile))
            use_vi = False
        else:
            truth = Table.read(truthfile)
            best = np.where(
                (truth['best quality'] >= 2.5) * 
                (truth['Redrock spectype'] == 'GALAXY')# *
                #(truth['Redrock z'] < 0.75) 
            )[0]
            #goodobj = np.where((zb['DELTACHI2'] > 100) * (zb['ZWARN'] == 0) * (zb['SPECTYPE'] == 'GALAXY'))[0]

            log.info('Read {}/{} good redshifts from {}'.format(len(best), len(truth), truthfile))
            truth = truth[best]
    
    zbest = []
    spectra, keepindx = [], []
    if False:
        targets = read_targets_in_tiles(tiles=np.atleast_1d(tileinfo), hpdirname=hpdirname, quick=True)
    
    #for spectro in ('0'):
    for spectro in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
        zbestfile = os.path.join(desidatadir, 'zbest-{}-{}-{}.fits'.format(spectro, tile, night))
        coaddfile = os.path.join(desidatadir, 'coadd-{}-{}-{}.fits'.format(spectro, tile, night))
        if os.path.isfile(zbestfile) and os.path.isfile(coaddfile):
            zb = Table(fitsio.read(zbestfile))
            fmap = fitsio.read(coaddfile, ext='FIBERMAP', columns=['FIBERSTATUS', 'OBJTYPE'])
            if use_vi:
                keep = np.where(np.isin(zb['TARGETID'], truth['TARGETID']))[0]
            else:
                keep = np.where((zb['Z'] > 0) * (zb['ZWARN'] == 0) * (fmap['OBJTYPE'] == 'TGT') *
                                (zb['SPECTYPE'] == 'GALAXY') * (fmap['FIBERSTATUS'] == 0))[0]

            if verbose:
                log.debug('Spectrograph {}: N={}'.format(spectro, len(keep)))
            if len(keep) > 0:
                zbest.append(zb[keep])
                keepindx.append(keep)
                spectra.append(desispec.io.read_spectra(coaddfile))

    if len(zbest) == 0:
        log.fatal('No spectra found!')
        raise ValueError
        
    # combine the spectrographs
    #log.debug('Stacking all the spectra.')
    zbest = vstack(zbest)

    coadd = None
    for camera in ('b', 'r', 'z'):
        wave = spectra[0].wave[camera]
        fm, flux, ivar, mask, res = [], [], [], [], []
        for ii in np.arange(len(spectra)):
            fm.append(spectra[ii].fibermap[keepindx[ii]])
            flux.append(spectra[ii].flux[camera][keepindx[ii], :])
            ivar.append(spectra[ii].ivar[camera][keepindx[ii], :])
            mask.append(spectra[ii].mask[camera][keepindx[ii], :])
            res.append(spectra[ii].resolution_data[camera][keepindx[ii], :, :])

        fm = vstack(fm)
        flux = np.concatenate(flux)
        ivar = np.concatenate(ivar)
        mask = np.concatenate(mask)
        res = np.concatenate(res)

        _coadd = Spectra([camera], {camera: wave}, {camera: flux}, {camera : ivar}, 
                         resolution_data={camera: res}, mask={camera: mask}, 
                         fibermap=fm, single=False)#, meta=meta) 
        if coadd is None:
            coadd = _coadd
        else:
            coadd.update(_coadd)

    assert(np.all(zbest['TARGETID'] == coadd.fibermap['TARGETID']))
    
    log.info('Writing {} redshifts to {}'.format(len(zbest), zbestoutfile))
    zbest.write(zbestoutfile, overwrite=True)

    log.info('Writing {} spectra to {}'.format(len(zbest), coaddoutfile))
    desispec.io.write_spectra(coaddoutfile, coadd)

    return zbest, coadd

def _unpack_one_spectrum(args):
    """Multiprocessing wrapper."""
    return unpack_one_spectrum(*args)

def unpack_one_spectrum(specobj, zbest, CFit, indx):
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
    from desitarget.io import desitarget_resolve_dec
    from fastspecfit.util import C_LIGHT

    cameras = ['b', 'r', 'z']
    ncam = len(cameras)

    ra = specobj.fibermap['TARGET_RA'][indx]
    dec = specobj.fibermap['TARGET_DEC'][indx]
    ebv = CFit.SFDMap.ebv(ra, dec)

    # Unpack the data and correct for Galactic extinction. Also flag pixels that
    # may be affected by emission lines.
    data = {'zredrock': zbest['Z'][indx], 'wave': [], 'flux': [], 'ivar': [],
            'res': [], 'linemask': [], 'snr': np.zeros(3).astype('f4')}
    for icam, camera in enumerate(cameras):
        dust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(specobj.wave[camera], Rv=CFit.RV))       
        data['wave'].append(specobj.wave[camera])
        data['flux'].append(specobj.flux[camera][indx, :] * dust)
        data['ivar'].append(specobj.ivar[camera][indx, :] / dust**2)
        data['res'].append(Resolution(specobj.resolution_data[camera][indx, :, :]))
        data['snr'][icam] = np.median(specobj.flux[camera][indx, :] * np.sqrt(specobj.ivar[camera][indx, :]))

        linemask = np.ones_like(specobj.wave[camera], bool)
        for line in CFit.linetable:
            zwave = line['restwave'] * (1 + data['zredrock'])
            I = np.where((specobj.wave[camera] >= (zwave - 1.5*CFit.linemask_sigma * zwave / C_LIGHT)) *
                         (specobj.wave[camera] <= (zwave + 1.5*CFit.linemask_sigma * zwave / C_LIGHT)))[0]
            if len(I) > 0:
                linemask[I] = False
        data['linemask'].append(linemask)

    # Make a quick coadd using inverse variance weights.
    uspecwave = np.unique(np.hstack(data['wave']))
    uspecflux3d = np.zeros((len(uspecwave), 3))
    uspecivar3d = np.zeros_like(uspecflux3d)
    for icam in np.arange(ncam):
        I = np.where(np.isin(data['wave'][icam], uspecwave))[0]
        J = np.where(np.isin(uspecwave, data['wave'][icam]))[0]
        uspecflux3d[J, icam] = data['flux'][icam][I]
        uspecivar3d[J, icam] = data['ivar'][icam][I]

    uspecivar = np.sum(uspecivar3d, axis=1)
    uspecflux = np.zeros_like(uspecivar)
    good = np.where(uspecivar > 0)[0]
    uspecflux[good] = np.sum(uspecivar3d[good, :] * uspecflux3d[good, :], axis=1) / uspecivar[good]
    data.update({'coadd_wave': uspecwave, 'coadd_flux': uspecflux, 'coadd_ivar': uspecivar})
    del uspecwave, uspecivar, uspecflux

    # Synthesize photometry, resolved into north and south.
    split = desitarget_resolve_dec()
    isouth = specobj.fibermap['TARGET_DEC'][indx] < split
    if isouth:
        filters = CFit.decamwise
    else:
        filters = CFit.bassmzlswise
    lambda_eff = filters.effective_wavelengths.value
    data['photsys_south'] = isouth

    padflux, padwave = filters.pad_spectrum(data['coadd_flux'], data['coadd_wave'], method='edge')
    synthmaggies = filters.get_ab_maggies(1e-17 * padflux, padwave)
    synthmaggies = synthmaggies.as_array().view('f8')[:3] # keep just grz

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
        maggies=synthmaggies, lambda_eff=lambda_eff[:3],
        nanomaggies=False)

    # Unpack the imaging photometry.
    bands = ['g', 'r', 'z', 'W1', 'W2']
    maggies = np.zeros(len(bands))
    ivarmaggies = np.zeros(len(bands))
    dust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(lambda_eff, Rv=CFit.RV))
    
    for iband, band in enumerate(bands):
        maggies[iband] = specobj.fibermap['FLUX_{}'.format(band.upper())][indx] / dust[iband]

        ivarcol = 'FLUX_IVAR_{}'.format(band.upper())
        if ivarcol in specobj.fibermap.colnames:
            ivarmaggies[iband] = specobj.fibermap[ivarcol][indx] * dust[iband]**2
        else:
            if maggies[iband] > 0:
                ivarmaggies[iband] = (10.0 / maggies [iband])**2 # constant S/N hack!!
    
    data['phot'] = CFit.parse_photometry(
        maggies=maggies, lambda_eff=lambda_eff,
        ivarmaggies=ivarmaggies, nanomaggies=True)

    return data

def unpack_all_spectra(specobj, zbest, CFit, fitindx, nproc=1):
    """Wrapper on unpack_one_spectrum to parse all the input spectra.

    Parameters
    ----------
    specobj : :class:`desispec.spectra.Spectra`
        DESI spectra (in the standard format).
    zbest : :class:`astropy.table.Table`
        Redrock redshift table (row-aligned to `specobj`).
    CFit : :class:`fastspecfit.continuum.ContinuumFit`
        Continuum-fitting class.
    fitindx : :class:`int`
        Index numbers of the spectra to unpack and pre-process.
    nproc : :class:`int`, optional, defaults to 1
        Number of cores to use for multiprocessing.
    
    Returns
    -------
    :class:`list`
        List of dictionaries (row-aligned to `fitindx`) populated by
        `unpack_one_spectrum`.
        
    Notes
    -----

    """
    args = [(specobj, zbest, CFit, indx) for indx in fitindx]
    if nproc > 1:
        with multiprocessing.Pool(nproc) as P:
            data = np.vstack(P.map(_unpack_one_spectrum, args))
    else:
        data = [unpack_one_spectrum(*_args) for _args in args]

    return data    

class ContinuumFit(object):
    def __init__(self, metallicity='Z0.0190', minwave=None, maxwave=6e4,
                 nproc=1, verbose=True):
        """Class to model a galaxy stellar continuum.

        Parameters
        ----------
        metallicity : :class:`str`, optional, defaults to `Z0.0190`.
            Stellar metallicity of the SSPs. Currently fixed at solar
            metallicity, Z=0.0190.
        minwave : :class:`float`, optional, defaults to None
            Minimum SSP wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxwave : :class:`float`, optional, defaults to 6e4
            Maximum SSP wavelength to read into memory. 
        nproc : :class:`int`, optional, defaults to 1
            Number of cores to use for multiprocessing.
        verbose : :class:`bool`, optional, defaults to False
            Trigger more verbose output throughout the class.

        Notes
        -----
        Need to document all the attributes.
        
        Plans for improvement (largely in self.fnnls_continuum).
          - Update the continuum redshift using cross-correlation.
          - Don't draw reddening from a flat distribution (try gamma or a custom
            distribution of the form x**2*np.exp(-2*x/scale).

        """
        import fitsio
        from astropy.cosmology import FlatLambdaCDM
        from astropy.table import Table

        from speclite import filters
        from desiutil.dust import SFDMap

        from fastspecfit.emlines import read_emlines

        # pre-compute the luminosity distance on a grid
        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        self.redshift_ref = np.arange(0.0, 5.0, 0.05)
        self.dlum_ref = self.cosmo.luminosity_distance(self.redshift_ref).to(u.pc).value

        self.metallicity = metallicity
        self.Z = float(metallicity[1:])
        self.library = 'CKC14z'
        self.isochrone = 'Padova' # would be nice to get MIST in here
        self.imf = 'Kroupa'

        self.nproc = nproc
        self.verbose = verbose

        self.fluxnorm = 1e17 # normalization factor for the spectra
        self.massnorm = 1e10 # stellar mass normalization factor for the SSPs [Msun]

        # dust and velocity dispersion
        self.SFDMap = SFDMap(scaling=0.86) # SF11 recalibration of the SFD maps
        self.RV = 3.1
        self.dustslope = 0.7

        # Initialize the velocity dispersion and reddening parameters. Make sure
        # the nominal values are in the grid.
        vdispmin, vdispmax, dvdisp, vdisp_nominal = (100.0, 350.0, 20.0, 150.0)
        #vdispmin, vdispmax, dvdisp, vdisp_nominal = (0.0, 0.0, 30.0, 150.0)
        nvdisp = np.ceil((vdispmax - vdispmin) / dvdisp).astype(int)
        if nvdisp == 0:
            nvdisp = 1
        vdisp = np.linspace(vdispmin, vdispmax, nvdisp).astype('f4') # [km/s]

        if not vdisp_nominal in vdisp:
            vdisp = np.sort(np.hstack((vdisp, vdisp_nominal)))
        self.vdisp = vdisp
        self.vdisp_nominal = vdisp_nominal
        self.nvdisp = len(vdisp)

        #AVmin, AVmax, dAV, AV_nominal = (0.0, 0.0, 0.1, 0.0)
        AVmin, AVmax, dAV, AV_nominal = (0.0, 1.0, 0.05, 0.0)
        nAV = np.ceil((AVmax - AVmin) / dAV).astype(int)
        if nAV == 0:
            nAV = 1
        AV = np.linspace(AVmin, AVmax, nAV).astype('f4')
        assert(AV[0] == 0.0) # minimum value has to be zero (assumed in fnnls_continuum)

        if not AV_nominal in AV:
            AV = np.sort(np.hstack((AV, AV_nominal)))        
        self.AV = AV
        self.AV_nominal = AV_nominal
        self.nAV = len(AV)

        # photometry
        self.bands = ['g', 'r', 'z', 'W1', 'W2']
        self.nband = len(self.bands)
        self.decam = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z')
        self.bassmzls = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z')

        self.decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z',
                                              'wise2010-W1', 'wise2010-W2')
        self.bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z',
                                                 'wise2010-W1', 'wise2010-W2')

        # SSPs
        self.sspfile = os.path.join(os.getenv('FASTSPECFIT_TEMPLATES'), 'SSP_{}_{}_{}_{}.fits'.format(
            self.isochrone, self.library, self.imf, self.metallicity))

        if verbose:
            log.info('Reading {}'.format(self.sspfile))
        wave, wavehdr = fitsio.read(self.sspfile, ext='WAVE', header=True)
        flux = fitsio.read(self.sspfile, ext='FLUX')
        sspinfo = Table(fitsio.read(self.sspfile, ext='METADATA'))
        
        # Trim the wavelengths and subselect the number of ages/templates.
        if minwave is None:
            minwave = np.min(wave)
        keep = np.where((wave >= minwave) * (wave <= maxwave))[0]
        sspwave = wave[keep]
        sspflux = flux[keep, ::5]
        sspinfo = sspinfo[::5]
        nage = len(sspinfo)
        npix = len(sspwave)

        self.pixkms = wavehdr['PIXSZBLU'] # pixel size [km/s]

        ## Resample the templates to have constant pixels in velocity /
        ## log-lambda and convolve to the nominal velocity dispersion.
        ## hack here
        #opt_pixkms = 50.0
        #ir_pixkms = 200
        #self.pixkms = opt_pixkms                         # SSP pixel size [km/s]
        #opt_dlogwave = opt_pixkms / C_LIGHT / np.log(10) # SSP pixel size [log-lambda]
        #ir_dlogwave = ir_pixkms / C_LIGHT / np.log(10) 
        #
        #wavesplit = 1e4
        #opt_sspwave = 10**np.arange(np.log10(wave.min()), np.log10(wavesplit), opt_dlogwave)
        #ir_sspwave = 10**np.arange(np.log10(wavesplit), np.log10(wave.max()), ir_dlogwave)
        #sspwave = np.hstack((opt_sspwave, ir_sspwave[1:]))
        #npix = len(sspwave)
        #
        ## None here means no resolution matrix.
        #args = [(flux[:, iage], wave, sspwave, None) for iage in np.arange(nage)]
        #if self.nproc > 1:
        #    with multiprocessing.Pool(self.nproc) as P:
        #        sspflux = np.vstack(P.map(_smooth_and_resample, args)).T # [npix, nage]
        #else:
        #    sspflux = np.vstack([smooth_and_resample(*_args) for _args in args]).T # [npix, nage]

        ## Build and store the nominal attenuation grid based on sspwave and the
        ## grid of AV values.
        #atten = []
        #for AV in self.AV:
        #    atten.append(self.dust_attenuation(sspwave, AV))
        #self.atten = np.stack(atten, axis=-1) # [npix, nAV]

        # Next, precompute a grid of spectra convolved to the nominal velocity
        # dispersion with reddening applied. This isn't quite right redward of
        # ~1 micron where the pixel size changes, but fix that later.
        sspflux_dustvdisp = []
        for AV in self.AV:
            atten = self.dust_attenuation(sspwave, AV)
            _sspflux_dustvdisp = self.convolve_vdisp(sspflux * atten[:, np.newaxis], self.vdisp_nominal)
            sspflux_dustvdisp.append(_sspflux_dustvdisp)
        sspflux_dustvdisp = np.stack(sspflux_dustvdisp, axis=-1) # [npix,nage,nAV]

        self.sspwave = sspwave
        self.sspflux = sspflux                     # no dust, no velocity broadening [npix,nage]
        self.sspflux_dustvdisp = sspflux_dustvdisp # nominal velocity broadening on a grid of A(V) [npix,nage,nAV]
        self.sspinfo = sspinfo
        self.nage = nage
        self.npix = npix

        # table of emission lines to fit
        self.linetable = read_emlines()
        self.linemask_sigma = 150.0 # [km/s]

        # Do a throw-away trapezoidal resampling so we can compile the numba
        # code when instantiating this class.
        #from redrock.rebin import trapz_rebin
        #t0 = time.time()
        #_ = trapz_rebin(np.arange(4), np.ones(4), np.arange(2)+1)
        #print('Initial rebin ', time.time() - t0)

    def init_spec_output(self, nobj=1):
        """Initialize the output data table for this class.

        """
        from astropy.table import Table, Column
        
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        out.add_column(Column(name='CONTINUUM_SNR', length=nobj, shape=(3,), dtype='f4')) # median S/N in each camera

        out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='CONTINUUM_VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        out.add_column(Column(name='CONTINUUM_VDISP_IVAR', length=nobj, dtype='f4', unit=u.second**2/u.kilometer**2))

        out['CONTINUUM_AV'] = self.AV_nominal
        out['CONTINUUM_VDISP'] = self.vdisp_nominal

        if False:
            # continuum fit with *no* dust reddening (to be used as a diagnostic
            # tool to identify potential calibration issues).
            out.add_column(Column(name='CONTINUUM_NODUST_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
            out.add_column(Column(name='CONTINUUM_NODUST_CHI2', length=nobj, dtype='f4')) # reduced chi2
            #out.add_column(Column(name='CONTINUUM_NODUST_AGE', length=nobj, dtype='f4', unit=u.Gyr))

        out.add_column(Column(name='DN4000', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_IVAR', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_NOLINES', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))

        return out

    def init_phot_output(self, nobj=1):
        """Initialize the photometric output data table.

        """
        from astropy.table import Table, Column
        
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        out.add_column(Column(name='PHOT_GRZW1W2', length=nobj, shape=(5,), dtype='f4', unit=u.nanomaggy)) # grzW1W2 photometry
        out.add_column(Column(name='PHOT_GRZW1W2_IVAR', length=nobj, shape=(5,), dtype='f4', unit=u.nanomaggy)) 
        #out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_PHOT_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_PHOT_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_PHOT_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_PHOT_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_PHOT_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_PHOT_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='DN4000_MODEL_PHOT', length=nobj, dtype='f4'))

        return out

    @staticmethod
    def parse_photometry(maggies, lambda_eff, ivarmaggies=None, nanomaggies=True,
                         flam=True, fnu=False, abmag=False):
        """Parse input (nano)maggies to various outputs and pack into a table.

        Parameters
        ----------
        flam - 10-17 erg/s/cm2/A
        fnu - 10-17 erg/s/cm2/Hz
        abmag - AB mag
        nanomaggies - input maggies are actually 1e-9 maggies

        Returns
        -------
        phot - photometric table

        Notes
        -----

        """
        from astropy.table import Table, Column
        
        shp = maggies.shape
        if maggies.ndim == 1:
            nband, ngal = shp[0], 1
        else:
            nband, ngal = shp[0], shp[1]
            
        phot = Table()
        phot.add_column(Column(name='lambda_eff', length=nband, dtype='f4'))
        phot.add_column(Column(name='nanomaggies', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='nanomaggies_ivar', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='flam', length=nband, shape=(ngal, ), dtype='f8')) # note f8!
        phot.add_column(Column(name='flam_ivar', length=nband, shape=(ngal, ), dtype='f8'))
        phot.add_column(Column(name='abmag', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag_ivar', length=nband, shape=(ngal, ), dtype='f4'))

        if ivarmaggies is None:
            ivarmaggies = np.zeros_like(maggies)

        phot['lambda_eff'] = lambda_eff.astype('f4')
        if nanomaggies:
            phot['nanomaggies'] = maggies.astype('f4')
            phot['nanomaggies_ivar'] = ivarmaggies.astype('f4')
        else:
            phot['nanomaggies'] = (maggies * 1e9).astype('f4')
            phot['nanomaggies_ivar'] = (ivarmaggies * 1e-18).astype('f4')

        if nanomaggies:
            nanofactor = 1e-9 # [nanomaggies-->maggies]
        else:
            nanofactor = 1.0

        factor = nanofactor * 10**(-0.4 * 48.6) * C_LIGHT * 1e13 / lambda_eff**2 # [maggies-->erg/s/cm2/A]
        if ngal > 1:
            factor = factor[:, None] # broadcast for the models
        phot['flam'] = (maggies * factor)
        phot['flam_ivar'] = (ivarmaggies / factor**2)

        # approximate the uncertainty as being symmetric in magnitude
        if maggies.ndim > 1:
            igood, jgood = np.unravel_index(np.where(maggies > 0)[0], maggies.shape)
            maggies = maggies[igood, jgood]
            ivarmaggies = ivarmaggies[igood, jgood]
        else:
            igood, jgood = np.where(maggies > 0)[0], [0]
            maggies = maggies[igood]
            ivarmaggies = ivarmaggies[igood]
            
        phot['abmag'][igood, jgood] = (-2.5 * np.log10(nanofactor * maggies)).astype('f4')
        phot['abmag_ivar'][igood, jgood] = (ivarmaggies * (maggies * 0.4 * np.log(10))**2).astype('f4')
        
        return phot

    def convolve_vdisp(self, sspflux, vdisp):
        """Convolve by the velocity dispersion.

        Parameters
        ----------
        sspflux
        vdisp

        Returns
        -------

        Notes
        -----

        """
        from scipy.ndimage import gaussian_filter1d

        if vdisp <= 0.0:
            return sspflux
        sigma = vdisp / self.pixkms # [pixels]
        smoothflux = gaussian_filter1d(sspflux, sigma=sigma, axis=0)
        return smoothflux
    
    def SSP2data(self, _sspflux, _sspwave, redshift=0.0, AV=None, vdisp=None,
                 specwave=None, specres=None, coeff=None, south=True,
                 synthphot=True, nproc=1):
        """Workhorse routine to turn input SSPs into spectra that can be compared to
        real data.

        Redshift, apply the resolution matrix, and resample in wavelength.

        Parameters
        ----------
        redshift
        specwave
        specres
        south
        synthphot - synthesize photometry?

        Returns
        -------
        Vector or 3-element list of [npix, nmodel] spectra.

        Notes
        -----
        This method does none or more of the following:
        - redshifting
        - wavelength resampling
        - apply dust reddening
        - apply velocity dispersion broadening
        - apply the resolution matrix
        - synthesize photometry

        It also naturally handles SSPs which have been precomputed on a grid of
        reddening or velocity dispersion (and therefore have an additional
        dimension). However, if the input grid is 3D, it is reshaped to be 2D
        but then it isn't reshaped back because of the way the photometry table
        is organized (bug or feature?).

        """
        # Are we dealing with a 2D grid [npix,nage] or a 3D grid
        # [npix,nage,nAV] or [npix,nage,nvdisp]?
        sspflux = _sspflux.copy() # why?!?
        sspwave = _sspwave.copy() # why?!?
        ndim = sspflux.ndim
        if ndim == 2:
            npix, nage = sspflux.shape
            nmodel = nage
        elif ndim == 3:
            npix, nage, nprop = sspflux.shape
            nmodel = nage*nprop
            sspflux = sspflux.reshape(npix, nmodel)
        else:
            log.fatal('Input SSPs have an unrecognized number of dimensions, {}'.format(ndim))
            raise ValueError
        
        #t0 = time.time()
        ##sspflux = sspflux.copy().reshape(npix, nmodel)
        #log.info('Copying the data took: {:.2f} sec'.format(time.time()-t0))

        # apply reddening
        if AV:
            atten = self.dust_attenuation(sspwave, AV)
            sspflux *= atten[:, np.newaxis]

        ## broaden for velocity dispersion
        #if vdisp:
        #    sspflux = self.convolve_vdisp(sspflux, vdisp)

        # Apply the redshift factor. The models are normalized to 10 pc, so
        # apply the luminosity distance factor here. Also normalize to a nominal
        # stellar mass.
        #t0 = time.time()
        if redshift:
            zsspwave = sspwave * (1.0 + redshift)
            #dfactor = (10.0 / self.cosmo.luminosity_distance(redshift).to(u.pc).value)**2
            dfactor = (10.0 / np.interp(redshift, self.redshift_ref, self.dlum_ref))**2
            factor = (self.fluxnorm * self.massnorm * dfactor / (1.0 + redshift))[np.newaxis, np.newaxis]
            #factor = np.asarray(self.fluxnorm * self.massnorm * dfactor / (1.0 + redshift))[np.newaxis, np.newaxis]
            #t0 = time.time()
            zsspflux = sspflux * factor
            #zsspflux = self.fluxnorm * self.massnorm * dfactor * sspflux / (1.0 + redshift)
            #print(time.time()-t0)
        else:
            zsspwave = sspwave.copy()
            zsspflux = self.fluxnorm * self.massnorm * sspflux
        #log.info('Cosmology calculations took: {:.2f} sec'.format(time.time()-t0))

        # Optionally synthesize photometry. We assume that velocity broadening,
        # if any, won't impact the measured photometry.
        sspphot = None
        if synthphot:
            if south:
                filters = self.decamwise
            else:
                filters = self.bassmzlswise
            effwave = filters.effective_wavelengths.value

            if ((specwave is None and specres is None and coeff is None) or
               (specwave is not None and specres is not None)):
                #t0 = time.time()
                maggies = filters.get_ab_maggies(zsspflux, zsspwave, axis=0) # speclite.filters wants an [nmodel,npix] array
                maggies = np.vstack(maggies.as_array().tolist()).T
                maggies /= self.fluxnorm * self.massnorm
                sspphot = self.parse_photometry(maggies, effwave, nanomaggies=False)
                #log.info('Synthesizing photometry took: {:.2f} sec'.format(time.time()-t0))
            
        # Are we returning per-camera spectra or a single model? Handle that here.
        #t0 = time.time()
        if specwave is None and specres is None:
            # multiprocess over age
            args = [(zsspflux[:, imodel], zsspwave, specwave, specres, vdisp, self.pixkms)
                    for imodel in np.arange(nmodel)]
            if nproc > 1:
                with multiprocessing.Pool(self.nproc) as P:
                    datasspflux = np.vstack(P.map(_smooth_and_resample, args)).T
            else:
                datasspflux = np.vstack([smooth_and_resample(*_args) for _args in args]).T
                
            if vdisp:
                 datasspflux = self.convolve_vdisp(datasspflux, vdisp)
                 
            # optionally compute the best-fitting model
            if coeff is not None:
                datasspflux = datasspflux.dot(coeff)
                if synthphot:
                    maggies = filters.get_ab_maggies(datasspflux, zsspwave, axis=0)
                    maggies = np.array(maggies.as_array().tolist()[0])
                    maggies /= self.fluxnorm * self.massnorm
                    sspphot = self.parse_photometry(maggies, effwave, nanomaggies=False)
        else:
            # loop over cameras and then multiprocess over age
            datasspflux = []
            for icamera in [0, 1, 2]: # iterate on cameras
                args = [(zsspflux[:, imodel], zsspwave, specwave[icamera], specres[icamera], vdisp, self.pixkms)
                        for imodel in np.arange(nmodel)]
                if nproc > 1:
                    with multiprocessing.Pool(self.nproc) as P:
                        datasspflux.append(np.vstack(P.map(_smooth_and_resample, args)).T)
                else:
                    _datasspflux = np.vstack([smooth_and_resample(*_args) for _args in args]).T

                if vdisp:
                    _datasspflux = self.convolve_vdisp(_datasspflux, vdisp)
                if coeff is not None:
                    _datasspflux = _datasspflux.dot(coeff)
                datasspflux.append(_datasspflux)
                
        #log.info('Resampling took: {:.2f} sec'.format(time.time()-t0))

        return datasspflux, sspphot # vector or 3-element list of [npix,nmodel] spectra

    def get_meanage(self, coeff):
        """Compute the light-weighted age, given a set of coefficients.

        """
        nage = len(coeff)
        age = self.sspinfo['age'][0:nage] # account for age of the universe trimming

        if np.count_nonzero(coeff > 0) == 0:
            log.warning('Coefficients are all zero!')
            meanage = -1.0
            #raise ValueError
        else:
            meanage = np.sum(coeff * age) / np.sum(coeff) / 1e9 # [Gyr]
        
        return meanage

    def younger_than_universe(self, redshift):
        """Return the indices of the SSPs younger than the age of the universe at the
        given redshift.

        """
        return np.where(self.sspinfo['age'] <= self.cosmo.age(redshift).to(u.year).value)[0]

    def _fnnls_parallel(self, modelflux, flux, ivar, xparam=None, debug=False):
        """Wrapper on fnnls to set up the multiprocessing. Works with both spectroscopic
        and photometric input and with both 2D and 3D model spectra.

        To be documented.

        """
        from redrock import fitz
        
        if xparam is not None:
            nn = len(xparam)
        ww = np.sqrt(ivar)
        xx = flux * ww

        # If xparam is None (equivalent to modelflux having just two
        # dimensions, [npix,nage]), assume we are just finding the
        # coefficients at some best-fitting value...
        #if modelflux.ndim == 2:
        if xparam is None:
            ZZ = modelflux * ww[:, np.newaxis]
            warn, coeff, chi2 = fnnls_continuum(ZZ, xx, flux=flux, ivar=ivar,
                                                modelflux=modelflux, get_chi2=True)
            if np.any(warn):
                print('WARNING: fnnls did not converge after 10 iterations.')

            return coeff, chi2

        # ...otherwise multiprocess over the xparam (e.g., AV or vdisp)
        # dimension.
        ZZ = modelflux * ww[:, np.newaxis, np.newaxis] # reshape into [npix/nband,nage,nAV/nvdisp]

        fitargs = [(ZZ[:, :, ii], xx, flux, ivar, modelflux[:, :, ii], None, True) for ii in np.arange(nn)]
        if self.nproc > 1:
            with multiprocessing.Pool(self.nproc) as P:
                rr = P.map(_fnnls_continuum, fitargs)
        else:
            #fnnls_continuum(*fitargs[10])
            #pdb.set_trace()
            rr = [fnnls_continuum(*_fitargs) for _fitargs in fitargs]
        warn, _, chi2grid = list(zip(*rr)) # unpack
        if np.any(warn):
            vals = ','.join(['{:.1f}'.format(xp) for xp in xparam[np.where(warn)[0]]])
            log.warning('fnnls did not converge after 10 iterations for parameter value(s) {}.'.format(vals))
        chi2grid = np.array(chi2grid)

        try:
            imin = fitz.find_minima(chi2grid)[0]
            xbest, xerr, chi2min, warn = fitz.minfit(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2])
        except:
            print('Problem here!', chi2grid)
            imin, xbest, xerr, chi2min, warn = 0, 0.0, 0.0, 0.0, 1

        #if np.all(chi2grid == 0):
        #    imin, xbest, xerr, chi2min, warn = 0, 0.0, 0.0, 0.0, 1
        #else:

        if warn == 0:
            xivar = 1.0 / xerr**2
        else:
            chi2min = 1e6
            xivar = 0.0

        if debug:
            import matplotlib.pyplot as plt
            plt.clf()
            plt.scatter(xparam, chi2grid)
            plt.scatter(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2], color='red')
            #plt.plot(xx, np.polyval([aa, bb, cc], xx), ls='--')
            plt.axvline(x=xbest, color='k')
            if xivar > 0:
                plt.axhline(y=chi2min, color='k')
            plt.yscale('log')
            plt.savefig('qa-chi2min.png')

        return chi2min, xbest, xivar

    def continuum_photfit(self, data):
        """Fit the broad photometry.

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            `unpack_one_spectrum`.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in `init_phot_output`.

        Notes
        -----
        See
          https://github.com/jvendrow/fnnls
          https://github.com/mikeiovine/fast-nnls
        for the fNNLS algorithm(s).

        """
        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_phot_output()

        redshift = data['zredrock']
        #result['CONTINUUM_Z'] = redshift
        result['PHOT_GRZW1W2'] = data['phot']['nanomaggies']
        result['PHOT_GRZW1W2_IVAR'] = data['phot']['nanomaggies_ivar']

        # Prepare the reddened and unreddened SSP templates. Note that we ignore
        # templates which are older than the age of the universe at the galaxy
        # redshift.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustvdisp, zsspphot_dustvdisp = self.SSP2data(
            self.sspflux_dustvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=None, specres=None,
            south=data['photsys_south'], nproc=1)
        log.info('Preparing the models took {:.2f} sec'.format(time.time()-t0))
        
        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = data['phot']['flam_ivar'].data / self.fluxnorm**2
        zsspflam_dustvdisp = zsspphot_dustvdisp['flam'].data * self.fluxnorm * self.massnorm # [nband,nage*nAV]
        assert(np.all(objflamivar >= 0))

        inodust = np.asscalar(np.where(self.AV == 0)[0]) # should always be index 0

        npix, nmodel = zsspflux_dustvdisp.shape
        nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)

        zsspflam_dustvdisp = zsspflam_dustvdisp.reshape(self.nband, nage, self.nAV) # [nband,nage,nAV]

        t0 = time.time()
        AVchi2min, AVbest, AVivar = self._fnnls_parallel(zsspflam_dustvdisp, objflam,
                                                         objflamivar, xparam=self.AV)
        log.info('Fitting the photometry took: {:.2f} sec'.format(time.time()-t0))
        if AVivar > 0:
            log.info('Best-fitting photometric A(V)={:.4f}+/-{:.4f} with chi2={:.3f}'.format(
                AVbest, 1/np.sqrt(AVivar), AVchi2min))
        else:
            AVbest = self.AV_nominal
            log.info('Finding photometric A(V) failed; adopting A(V)={:.4f}'.format(self.AV_nominal))

        # Get the final set of coefficients and chi2 at the best-fitting
        # reddening and nominal velocity dispersion.
        bestsspflux, bestphot = self.SSP2data(self.sspflux_dustvdisp[:, agekeep, inodust], # equivalent to calling with self.sspflux[:, agekeep]
                                              self.sspwave, AV=AVbest, redshift=redshift,
                                              south=data['photsys_south'])
        coeff, chi2min = self._fnnls_parallel(bestphot['flam'].data*self.massnorm*self.fluxnorm,
                                              objflam, objflamivar) # bestphot['flam'] is [nband, nage]

        continuum = bestsspflux.dot(coeff)
        dn4000, _ = get_dn4000(self.sspwave, continuum, rest=True)
        meanage = self.get_meanage(coeff)
        log.info('Photometric DN(4000)={:.3f}, Age={:.2f} Gyr'.format(dn4000, meanage))

        result['CONTINUUM_PHOT_COEFF'][0][:nage] = coeff
        result['CONTINUUM_PHOT_CHI2'][0] = chi2min
        result['CONTINUUM_PHOT_AGE'][0] = meanage
        result['CONTINUUM_PHOT_AV'][0] = AVbest
        result['CONTINUUM_PHOT_AV_IVAR'][0] = AVivar
        result['DN4000_MODEL_PHOT'][0] = dn4000

        return result, continuum
    
    def continuum_specfit(self, data, solve_vdisp=False):
        """Fit the stellar continuum of a single spectrum using fast non-negative
        least-squares fitting (fNNLS).

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            `unpack_one_spectrum`.
        solve_vdisp : :class:`bool`, optional, defaults to False
            Solve for the velocity dispersion.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in `init_fastspecfit`.

        Notes
        -----
        See https://github.com/jvendrow/fnnls for the fNNLS algorithm.

        ToDo:
          - Use cross-correlation to update the redrock redshift.
          - Need to mask more emission lines than we fit (e.g., Mg II).

        """
        from fnnls import fnnls

        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_spec_output()

        redshift = data['zredrock']
        result['CONTINUUM_Z'] = redshift
        result['CONTINUUM_SNR'] = data['snr']

        # Prepare the reddened and unreddened SSP templates. Note that we ignore
        # templates which are older than the age of the universe at the galaxy
        # redshift.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustvdisp, _ = self.SSP2data(
            self.sspflux_dustvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=data['wave'], specres=data['res'],
            nproc=1, synthphot=False)
        log.info('Preparing the models took {:.2f} sec'.format(time.time()-t0))
        
        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        npixpercamera = [len(gw) for gw in data['wave']]
        npixpercam = np.hstack([0, npixpercamera])
        
        specwave = np.hstack(data['wave'])
        specflux = np.hstack(data['flux'])
        specivar = np.hstack(data['ivar']) * np.hstack(data['linemask']) # mask emission lines
        zsspflux_dustvdisp = np.concatenate(zsspflux_dustvdisp, axis=0)  # [npix,nage*nAV]
        assert(np.all(specivar >= 0))

        inodust = np.asscalar(np.where(self.AV == 0)[0]) # should always be index 0

        npix, nmodel = zsspflux_dustvdisp.shape
        nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)

        zsspflux_dustvdisp = zsspflux_dustvdisp.reshape(npix, nage, self.nAV)       # [npix,nage,nAV]

        if False:
            # Fit the spectra with *no* dust reddening so we can identify potential
            # calibration issues (again, at the nominal velocity dispersion).
            t0 = time.time()
            coeff, chi2min = self._fnnls_parallel(zsspflux_dustvdisp[:, :, inodust],
                                                  specflux, specivar)
            log.info('No-dust model fit has chi2={:.3f} and took {:.2f} sec'.format(
                chi2min, time.time()-t0))

            result['CONTINUUM_NODUST_COEFF'][0][0:nage] = coeff
            result['CONTINUUM_NODUST_CHI2'] = chi2min

        # Fit the spectra for reddening using the models convolved to the
        # nominal velocity dispersion and then fit for velocity dispersion.
        t0 = time.time()
        AVchi2min, AVbest, AVivar = self._fnnls_parallel(zsspflux_dustvdisp, specflux, specivar,
                                                         xparam=self.AV, debug=False)
        log.info('Fitting for the reddening took: {:.2f} sec'.format(time.time()-t0))
        if AVivar > 0:
            log.info('Best-fitting spectroscopic A(V)={:.4f}+/-{:.4f} with chi2={:.3f}'.format(
                AVbest, 1/np.sqrt(AVivar), AVchi2min))
        else:
            AVbest = self.AV_nominal
            log.info('Finding spectroscopic A(V) failed; adopting A(V)={:.4f}'.format(
                self.AV_nominal))

        # Optionally build out the model spectra on our grid of velocity
        # dispersion and then solve.
        if solve_vdisp:
            t0 = time.time()
            zsspflux_vdisp = []
            for vdisp in self.vdisp:
                _zsspflux_vdisp, _ = self.SSP2data(self.sspflux[:, agekeep], self.sspwave,
                                                   specwave=data['wave'], specres=data['res'],
                                                   AV=AVbest, vdisp=vdisp, redshift=redshift,
                                                   synthphot=False)
                _zsspflux_vdisp = np.concatenate(_zsspflux_vdisp, axis=0)
                zsspflux_vdisp.append(_zsspflux_vdisp)

            zsspflux_vdisp = np.stack(zsspflux_vdisp, axis=-1) # [npix,nage,nvdisp] at best A(V)
            vdispchi2min, vdispbest, vdispivar = self._fnnls_parallel(zsspflux_vdisp, specflux, specivar,
                                                                      xparam=self.vdisp, debug=False)
            log.info('Fitting for the velocity dispersion took: {:.2f} sec'.format(time.time()-t0))
            if vdispivar > 0:
                log.info('Best-fitting vdisp={:.2f}+/-{:.2f} km/s with chi2={:.3f}'.format(
                    vdispbest, 1/np.sqrt(vdispivar), vdispchi2min))
            else:
                vdispbest = self.vdisp_nominal
                log.info('Finding vdisp failed; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))
        else:
            vdispbest, vdispivar = self.vdisp_nominal, 0.0

        # Get the final set of coefficients and chi2 at the best-fitting
        # reddening and velocity dispersion.
        bestsspflux, bestphot = self.SSP2data(self.sspflux[:, agekeep], self.sspwave,
                                              specwave=data['wave'], specres=data['res'],
                                              AV=AVbest, vdisp=vdispbest, redshift=redshift,
                                              south=data['photsys_south'])
        bestsspflux = np.concatenate(bestsspflux, axis=0)
        coeff, chi2min = self._fnnls_parallel(bestsspflux, specflux, specivar)

        bestfit = bestsspflux.dot(coeff)
        meanage = self.get_meanage(coeff)
        dn4000_model, _ = get_dn4000(specwave, bestfit, redshift=redshift)
        dn4000, dn4000_ivar = get_dn4000(specwave, specflux, specivar, redshift=redshift)
        log.info('Spectroscopic DN(4000)={:.3f}, Age={:.2f} Gyr'.format(dn4000, meanage))

        result['CONTINUUM_COEFF'][0][0:nage] = coeff
        result['CONTINUUM_CHI2'][0] = chi2min
        result['CONTINUUM_AV'][0] = AVbest
        result['CONTINUUM_AV_IVAR'][0] = AVivar
        result['CONTINUUM_VDISP'][0] = vdispbest
        result['CONTINUUM_VDISP_IVAR'][0] = vdispivar
        result['CONTINUUM_AGE'] = meanage
        result['DN4000'][0] = dn4000
        result['DN4000_IVAR'][0] = dn4000_ivar
        result['DN4000_MODEL'][0] = dn4000_model

        # Unpack the continuum into individual cameras.
        continuum = []
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npixpercam[:ii+1])
            jpix = np.sum(npixpercam[:ii+2])
            continuum.append(bestfit[ipix:jpix])

        return result, continuum
    
    def qa_continuum(self, data, specfit, photfit, qadir='.'):
        """QA of the best-fitting continuum.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        from fastspecfit.util import ivar2var
    
        sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]
        ymin, ymax = 1e6, -1e6

        redshift = data['zredrock']

        # rebuild the best-fitting spectroscopic and photometric models
        inodust = np.asscalar(np.where(self.AV == 0)[0]) # should always be index 0            
        agekeep = self.younger_than_universe(redshift)
        continuum, _ = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift, 
                                     specwave=data['wave'], specres=data['res'],
                                     AV=specfit['CONTINUUM_AV'],
                                     vdisp=specfit['CONTINUUM_VDISP'],
                                     coeff=specfit['CONTINUUM_COEFF'],
                                     synthphot=False)
        continuum_phot, bestphot = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift,
                                                 AV=photfit['CONTINUUM_PHOT_AV'],
                                                 coeff=photfit['CONTINUUM_PHOT_COEFF'] * self.massnorm,
                                                 south=data['photsys_south'],
                                                 synthphot=True)
        continuum_wave_phot = self.sspwave * (1 + redshift)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
        for ii in [0, 1, 2]: # iterate over cameras
            sigma, _ = ivar2var(data['ivar'][ii], sigma=True)

            ax1.fill_between(data['wave'][ii]/1e4, data['flux'][ii]-sigma,
                             data['flux'][ii]+sigma, color=col1[ii])
            ax1.plot(data['wave'][ii]/1e4, continuum[ii], color=col2[ii], alpha=1.0)#, color='k')
            #ax1.plot(data['wave'][ii]/1e4, continuum_nodust[ii], alpha=0.5, color='k')

            # get the robust range
            filtflux = median_filter(data['flux'][ii], 5)
            if np.min(filtflux) < ymin:
                ymin = np.min(filtflux) * 0.5
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux) * 1.5

        leg = {
            'targetid': 'targetid={} fiber={}'.format(specfit['TARGETID'], specfit['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(specfit['CONTINUUM_CHI2']),
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(specfit['Z']),
            #'zfastspecfit': '$z_{{\\rm fastspecfit}}$={:.6f}'.format(specfit['CONTINUUM_Z']),
            'z': '$z$={:.6f}'.format(specfit['CONTINUUM_Z']),
            'age': '<Age>={:.3f} Gyr'.format(specfit['CONTINUUM_AGE']),
            }
        if specfit['CONTINUUM_VDISP_IVAR'] == 0:
            leg.update({'vdisp': '$\sigma$={:.1f} km/s'.format(specfit['CONTINUUM_VDISP'])})
        else:
            leg.update({'vdisp': '$\sigma$={:.1f}+/-{:.1f} km/s'.format(
                specfit['CONTINUUM_VDISP'], 1/np.sqrt(specfit['CONTINUUM_VDISP_IVAR']))})
        if specfit['CONTINUUM_AV_IVAR'] == 0:
            leg.update({'AV': '$A(V)$={:.3f} mag'.format(specfit['CONTINUUM_AV'])})
        else:
            leg.update({'AV': '$A(V)$={:.3f}+/-{:.3f} mag'.format(
                specfit['CONTINUUM_AV'], 1/np.sqrt(specfit['CONTINUUM_AV_IVAR']))})

        ax1.text(0.95, 0.92, '{}'.format(leg['targetid']), 
                 ha='right', va='center', transform=ax1.transAxes, fontsize=14)
        #ax1.text(0.95, 0.92, '{} {}'.format(leg['targetid'], leg['zredrock']), 
        #         ha='right', va='center', transform=ax1.transAxes, fontsize=14)
        ax1.text(0.95, 0.86, r'{} {}'.format(leg['z'], leg['chi2']),
                 ha='right', va='center', transform=ax1.transAxes, fontsize=14)
        ax1.text(0.95, 0.80, r'{}'.format(leg['age']),
                 ha='right', va='center', transform=ax1.transAxes, fontsize=14)
        ax1.text(0.95, 0.74, r'{}'.format(leg['AV']),
                 ha='right', va='center', transform=ax1.transAxes, fontsize=14)
        ax1.text(0.95, 0.68, r'{}'.format(leg['vdisp']),
                 ha='right', va='center', transform=ax1.transAxes, fontsize=14)
                    
        ax1.set_xlim(3500/1e4, 10000/1e4)
        ax1.set_ylim(ymin, ymax)
        #ax1.set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
        #ax1.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        ax1.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 

        # add the photometry
        if False:
            for ii in [0, 1, 2]: # iterate over cameras
                #galsigma = 1 / np.sqrt(specivar[ii])
                factor = 1e-17  * specwave[ii]**2 / (C_LIGHT * 1e13) # [10-17 erg/s/cm2/A --> maggies]
                good = np.where(specflux[ii] > 0)[0]
                if len(good) > 0:
                    ax2.plot(specwave[ii][good]/1e4, -2.5*np.log10(specflux[ii][good]*factor[good])-48.6, color=col1[ii])
                    #ax1.fill_between(specwave[ii]/1e4, -2.5*np.log10((specflux[ii]-galsigma) * factor,
                    #                 (specflux[ii]+galsigma) * factor, color=col1[ii])
                #ax2.plot(specwave[ii]/1e4, -2.5*np.log10(continuum[ii]*factor)-48.6, color=col2[ii], alpha=1.0)#, color='k')
                ax2.plot(specwave[ii]/1e4, -2.5*np.log10(continuum[ii]*factor)-48.6, color=col2[ii], alpha=1.0)#, color='k')

        wavemin, wavemax = 0.2, 6.0
        indx = np.where((continuum_wave_phot/1e4 > wavemin) * (continuum_wave_phot/1e4 < wavemax))[0]

        if np.any(continuum_phot <= 0):
            log.warning('Best-fitting photometric continuum is all zeros or negative!')
            continuum_phot_abmag = continuum_phot*0 + np.median(data['phot']['abmag'][:3])
        else:
            factor = 10**(0.4 * 48.6) * continuum_wave_phot**2 / (C_LIGHT * 1e13) / self.fluxnorm / self.massnorm # [erg/s/cm2/A --> maggies]
            continuum_phot_abmag = -2.5*np.log10(continuum_phot * factor)
            ax2.plot(continuum_wave_phot[indx] / 1e4, continuum_phot_abmag[indx], color='gray', zorder=1)

        ax2.scatter(data['synthphot']['lambda_eff']/1e4, data['synthphot']['abmag'], 
                     marker='o', s=130, color='blue', edgecolor='k',
                     label=r'$grz$ (synthesized)', alpha=1.0, zorder=2)
        ax2.scatter(data['phot']['lambda_eff']/1e4, data['phot']['abmag'],
                    marker='s', s=130, facecolor='red', edgecolor='k',
                    label=r'$grzW1W2$ (imaging)', alpha=1.0, zorder=3)
        #abmag_sigma, _ = _ivar2var(data['phot']['abmag_ivar'], sigma=True)
        #ax2.errorbar(data['phot']['lambda_eff']/1e4, data['phot']['abmag'], yerr=abmag_sigma,
        #             fmt='s', markersize=15, markeredgewidth=3, markeredgecolor='k', markerfacecolor='red',
        #             elinewidth=3, ecolor='blue', capsize=3)
        ax2.legend(loc='lower right', fontsize=16)

        dm = 0.75
        good = data['phot']['abmag_ivar'] > 0
        ymin = np.max((np.nanmax(data['phot']['abmag'][good]),
                       np.nanmax(continuum_phot_abmag[indx]))) + dm
        ymax = np.min((np.nanmin(data['phot']['abmag'][good]),
                       np.nanmin(continuum_phot_abmag[indx]))) - dm
        if np.isnan(ymin) or np.isnan(ymax):
            pdb.set_trace()

        ax2.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        ax2.set_ylabel(r'AB Mag') 
        ax2.set_xlim(wavemin, wavemax)
        ax2.set_ylim(ymin, ymax)

        ax2.set_xscale('log')
        ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax2.set_xticks([0.3, 0.4, 0.6, 1.0, 1.5, 2.5, 5.0])

        plt.subplots_adjust(bottom=0.1, right=0.95, top=0.95, wspace=0.12)
        #plt.subplots_adjust(bottom=0.1, right=0.95, top=0.95, wspace=0.17)

        pngfile = os.path.join(qadir, 'continuum-{}-{}-{}.png'.format(
            specfit['TILE'], specfit['NIGHT'], specfit['TARGETID']))
        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()

        return continuum
        
    def _get_uncertainties(self, pcov=None, jac=None, return_covariance=False):
        """Determine the uncertainties from the diagonal terms of the covariance
        matrix. If the covariance matrix is not known, estimate it from the
        Jacobian following
        https://github.com/scipy/scipy/blob/master/scipy/optimize/minpack.py#L805

        """
        if pcov is None:
            from scipy.linalg import svd
            _, s, VT = svd(jac, full_matrices=False)
            threshold = np.finfo(float).eps * max(jac.shape) * s[0]
            s = s[s > threshold]
            VT = VT[:s.size]
            pcov = np.dot(VT.T / s**2, VT)

        unc = np.diag(pcov)**0.5

        if return_covariance:
            return unc, pcov
        else:
            return unc

    def dust_attenuation(self, wave, AV):
        """Compute the dust attenuation curve A(lambda)/A(V) from Charlot & Fall 2000.

        ToDo: add a UV bump and IGM attenuation!
          https://gitlab.lam.fr/cigale/cigale/-/blob/master/pcigale/sed_modules/dustatt_powerlaw.py#L42

        """
        return 10**(-0.4 * AV * (wave / 5500.0)**(-self.dustslope))
        
    def dusty_continuum(self, ebv_and_coeffs, wave, sspflux):
        """Continuum model with dust attenuation."""
        ebv, coeffs = ebv_and_coeffs[0], ebv_and_coeffs[1:]
        atten = self.dust_attenuation(wave, ebv)
        modelflux = (sspflux * atten[:, None]).dot(coeffs) 
        return modelflux
    
    def _dusty_continuum_resid(self, ebv_and_coeffs, wave, flux, isigma, sspflux):
        """Wrapper cost function for scipy.optimize.least_squares."""
        modelflux = self.dusty_continuum(ebv_and_coeffs, wave, sspflux)
        return isigma * (modelflux - flux)
    
    def fit_continuum(self):
        """More detailed continuum model fit, including dust attenuation.

        https://stackoverflow.com/questions/60335524/how-to-use-scipy-optimize-curve-fit-to-use-lists-of-variable/60409667#60409667 
        
        """
        from scipy.optimize import least_squares

        sspflux = np.concatenate(self.redshift_smooth_and_resample(), axis=0)

        wave = np.hstack(self.specwave)
        flux = np.hstack(self.specflux)
        isigma = 1 / np.sqrt(np.hstack(self.specivar))

        _ = self.fnnls_continuum()
        init_params = np.hstack([0.05, self.fnnls_coeffs])
        log.info(init_params)

        params = least_squares(self._dusty_continuum_resid, x0=init_params, #kwargs=sspflux,
                               bounds=(0.0, np.inf), args=(wave, flux, isigma, sspflux),
                               method='trf', tr_options={'regularize': True})
        #params, cov = curve_fit(self._dusty_continuum, wave, flux, sigma=1/np.sqrt(ivar),
        #                        kwargs=sspflux)
        #continuum_fit = fitter(ContinuumModel, wave, flux, 
        unc = self._get_uncertainties(jac=params.jac, return_covariance=False)
        
        log.info(params.x[0], self.ssp.info['age'][params.x[1:] > 1] / 1e9)
        #pdb.set_trace()

