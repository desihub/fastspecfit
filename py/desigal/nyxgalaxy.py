"""
desigal.nyxgalaxy
=================

"""
import pdb # for debugging

import os, time
import numpy as np
import multiprocessing

from scipy.ndimage import gaussian_filter1d
import astropy.units as u
from astropy.table import Table, Column, vstack, join, hstack
from astropy.modeling import Fittable1DModel

from fnnls import fnnls
from desiutil.log import get_logger
#from desispec.interpolation import resample_flux
from redrock.rebin import trapz_rebin

from scipy import constants
C_LIGHT = constants.c / 1000.0 # [km/s]

log = get_logger()

def init_nyxgalaxy(tile, night, zbest, fibermap, CFit, EMFit):
    """Initialize the nyxgalaxy output data table.

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
    CFit : :class:`desigal.nyxgalaxy.ContinuumFit`
        Continuum-fitting class.
    CFit : :class:`desigal.nyxgalaxy.EMLineFit`
        Emission-line fitting class.

    Returns
    -------
    

    Notes
    -----

    """
    # Grab info on the emission lines and the continuum.
    nobj = len(zbest)

    out = Table()
    for zbestcol in ['TARGETID', 'Z']:#, 'ZERR']:#, 'SPECTYPE', 'DELTACHI2']
        out[zbestcol] = zbest[zbestcol]
    for fibermapcol in ['FIBER']:
        out[fibermapcol] = fibermap[fibermapcol]
    out.add_column(Column(name='NIGHT', data=np.repeat(night, nobj)), index=0)
    out.add_column(Column(name='TILE', data=np.repeat(tile, nobj)), index=0)

    out = hstack((out, CFit.init_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

    return out
    
def read_spectra(tile, night, use_vi=False, write_spectra=True, overwrite=False,
                 verbose=False):
    """Read the spectra and redshift catalog for a given tile and night.

    Parameters
    ----------
    tile : :class:`str`
        Tile number to analyze.
    night : :class:`str`
        Night on which `tile` was observed.
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
    from astropy.table import join
    from desispec.spectra import Spectra
    import desispec.io

    data_dir = os.path.join(os.getenv('NYXGALAXY_DATA'), 'spectra')
    if not os.path.isdir(data_dir):
        if verbose:
            log.debug('Creating directory {}'.format(data_dir))
        os.makedirs(data_dir, exist_ok=True)
    
    zbestoutfile = os.path.join(data_dir, 'zbest-{}-{}.fits'.format(tile, night))
    coaddoutfile = os.path.join(data_dir, 'coadd-{}-{}.fits'.format(tile, night))
    if os.path.isfile(coaddoutfile) and not overwrite:
        zbest = Table(fitsio.read(zbestoutfile))
        log.info('Read {} redshifts from {}'.format(len(zbest), zbestoutfile))

        coadd = desispec.io.read_spectra(coaddoutfile)
        log.info('Read {} spectra from {}'.format(len(zbest), coaddoutfile))

        assert(np.all(zbest['TARGETID'] == coadd.fibermap['TARGETID']))

        return zbest, coadd

    log.info('Parsing data from tile, night {}, {}'.format(tile, night))

    specprod_dir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', 'andes')
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
    for spectro in ('0'):
    #for spectro in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
        zbestfile = os.path.join(desidatadir, 'zbest-{}-{}-{}.fits'.format(spectro, tile, night))
        coaddfile = os.path.join(desidatadir, 'coadd-{}-{}-{}.fits'.format(spectro, tile, night))
        if os.path.isfile(zbestfile) and os.path.isfile(coaddfile):
            zb = Table(fitsio.read(zbestfile))
            fiberstatus = fitsio.read(coaddfile, ext='FIBERMAP', columns='FIBERSTATUS')
            if use_vi:
                keep = np.where(np.isin(zb['TARGETID'], truth['TARGETID']))[0]
            else:
                #snr = np.median(specobj.flux['r']*np.sqrt(specobj.ivar['r']), axis=1)
                keep = np.where((zb['Z'] > 0) * (zb['ZWARN'] == 0) *
                                (zb['SPECTYPE'] == 'GALAXY') * (fiberstatus == 0))[0]
                #keep = np.where((zb['ZWARN'] == 0) * (zb['DELTACHI2'] > 50) * (zb['SPECTYPE'] == 'GALAXY'))[0]

            if verbose:
                log.debug('Spectrograph {}: N={}'.format(spectro, len(keep)))
            if len(keep) > 0:
                zbest.append(zb[keep])
                keepindx.append(keep)
                spectra.append(desispec.io.read_spectra(coaddfile))

    if len(zbest) == 0:
        log.fatal('No spectra found for tile {} and night {}!'.format(tile, night))
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

def _ivar2var(ivar, sigma=False):
    """Safely convert an inverse variance to a variance."""
    var = np.zeros_like(ivar)
    goodmask = ivar > 0 # True is good
    if np.count_nonzero(goodmask) == 0:
        log.warning('All values are masked!')
        raise ValueError
    var[goodmask] = 1 / ivar[goodmask]
    if sigma:
        var = np.sqrt(var) # return a sigma
    return var, goodmask

def get_d4000(wave, flam, flam_ivar=None, redshift=None, rest=True):
    """Compute D(4000) and, optionally, the inverse variance.

    Parameters
    ----------
    wave
    flam
    flam_ivar
    redshift
    rest

    Returns
    -------

    Notes
    -----
    If `rest`=``False`` then `redshift` input is required.

    """
    d4000, d4000_ivar = 0.0, 0.0

    if rest:
        flam2fnu =  wave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
    else:
        wave /= (1 + redshift) # [Angstrom]
        flam2fnu = (1 + redshift) * wave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]

    if flam_ivar is None:
        goodmask = np.ones(len(flam), bool) # True is good
    else:
        goodmask = flam_ivar > 0

    indxblu = np.where((wave >= 3850.) * (wave <= 3950.) * goodmask)[0]
    indxred = np.where((wave >= 4000.) * (wave <= 4100.) * goodmask)[0]
    if len(indxblu) < 5 or len(indxred) < 5:
        return d4000, d4000_ivar

    blufactor, redfactor = 3950.0 - 3850.0, 4100.0 - 4000.0
    deltawave = np.gradient(wave) # should be constant...

    fnu = flam * flam2fnu # [erg/s/cm2/Hz]

    numer = blufactor * np.sum(deltawave[indxred] * fnu[indxred])
    denom = redfactor * np.sum(deltawave[indxblu] * fnu[indxblu])
    if denom == 0.0:
        log.warning('D(4000) is ill-defined!')
        return d4000, d4000_ivar
    d4000 =  numer / denom

    if flam_ivar is not None:
        fnu_ivar = flam_ivar / flam2fnu**2
        fnu_var, _ = _ivar2var(fnu_ivar)

        numer_var = blufactor**2 * np.sum(deltawave[indxred] * fnu_var[indxred])
        denom_var = redfactor**2 * np.sum(deltawave[indxblu] * fnu_var[indxblu])
        d4000_var = (numer_var + numer**2 * denom_var) / denom**2
        if d4000_var <= 0:
            log.warning('D(4000) variance is ill-defined!')
            d4000_ivar = 0.0
        else:
            d4000_ivar = 1.0 / d4000_var

    return d4000, d4000_ivar

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
    CFit : :class:`desigal.nyxgalaxy.ContinuumFit`
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

    cameras = ['b', 'r', 'z']
    ncam = len(cameras)

    ra = specobj.fibermap['TARGET_RA'][indx]
    dec = specobj.fibermap['TARGET_DEC'][indx]
    ebv = CFit.SFDMap.ebv(ra, dec, scaling=1.0) # SFD coefficients

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
    for iband, band in enumerate(bands):
        dust = specobj.fibermap['MW_TRANSMISSION_{}'.format(band.upper())][indx]
        maggies[iband] = specobj.fibermap['FLUX_{}'.format(band.upper())][indx] / dust
        ivarmaggies[iband] = specobj.fibermap['FLUX_IVAR_{}'.format(band.upper())][indx] * dust**2
    
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
    CFit : :class:`desigal.nyxgalaxy.ContinuumFit`
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

def _fit_fnnls_continuum(myargs):
    """Multiprocessing wrapper."""
    return fit_fnnls_continuum(*myargs)

def fit_fnnls_continuum(ZZ, xx, flux=None, ivar=None, modelflux=None,
                        support=None, get_chi2=False, printme=None):
    """Fit a continuum using fNNLS. This function is a simple wrapper on fnnls; see
    the ContinuumFit.fnnls_continuum method for documentation.

    """
    if support is None:
        support = np.zeros(0, dtype=int)
    #print(printme)
    warn, coeff, _ = fnnls(ZZ, xx, P_initial=support)
    #if warn:
    #    print('WARNING: fnnls did not converge after 5 iterations.')

    if get_chi2:
        chi2 = np.sum(ivar * (flux - modelflux.dot(coeff))**2)
        chi2 /= np.sum(ivar > 0) # reduced chi2
        return warn, coeff, chi2
    else:
        return warn, coeff

def convolve_vdisp(sspflux, pixkms, vdisp):
    """Convolve by the velocity dispersion.

    Parameters
    ----------
    sspflux
    pixkms - 
    vdisp

    Returns
    -------

    Notes
    -----

    """
    if vdisp <= 0.0:
        return sspflux
    sigma = vdisp / pixkms # [pixels]

    return gaussian_filter1d(sspflux, sigma=sigma, axis=0)

def _smooth_and_resample(args):
    """Multiprocessing wrapper."""
    return smooth_and_resample(*args)

def smooth_and_resample(sspflux, sspwave, specwave=None, specres=None,
                        vdisp=None, pixkms=None):
    """Given a single template, apply the resolution matrix and resample in
    wavelength.
    
    Parameters
    ----------
    sspflux : :class:`numpy.ndarray` [npix]
        Input (model) spectrum.
    sspwave : :class:`numpy.ndarray` [npix]
        Wavelength array corresponding to `sspflux`.
    specwave : :class:`numpy.ndarray` [noutpix], optional, defaults to None
        Desired output wavelength array, usually that of the object being fitted.
    specres : :class:`desispec.resolution.Resolution`, optional, defaults to None 
        Resolution matrix.
    vdisp : :class:`float`, optional, defaults to None
        Velocity dispersion broadening factor [km/s].
    pixkms : :class:`float`, optional, defaults to None
        Pixel size of input spectra [km/s].
    
    Returns
    -------
    :class:`numpy.ndarray` [noutpix]
        Smoothed and resampled flux at the new resolution and wavelength sampling.
        
    Notes
    -----
    This function stands by itself rather than being in a class because we call
    it with multiprocessing, below.

    """
    if specwave is None:
        resampflux = sspflux 
    else:
        #t0 = time.time()
        trim = (sspwave > (specwave.min()-10.0)) * (sspwave < (specwave.max()+10.0))
        #resampflux = trapz_rebin(sspwave, sspflux, specwave)
        resampflux = trapz_rebin(sspwave[trim], sspflux[trim], specwave)
        #print(time.time()-t0)
        #resampflux = resample_flux(specwave, sspwave, sspflux, extrapolate=True)
        #resampflux = np.interp(specwave, sspwave, sspflux)

    ## broaden for velocity dispersion
    #print('Turning off vdisp')
    #if vdisp:
    #    resampflux = convolve_vdisp(resampflux, pixkms, vdisp)

    if specres is None:
        smoothflux = resampflux
    else:
        smoothflux = specres.dot(resampflux)
        
    return smoothflux # [noutpix]

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

        from speclite import filters
        from desiutil.dust import SFDMap

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
        self.SFDMap = SFDMap()
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
        self.nband = 5
        self.decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z',
                                              'wise2010-W1', 'wise2010-W2')
        self.bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z',
                                                 'wise2010-W1', 'wise2010-W2')

        # SSPs
        self.sspfile = os.path.join(os.getenv('NYXGALAXY_TEMPLATES'), 'SSP_{}_{}_{}_{}.fits'.format(
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
            _sspflux_dustvdisp = convolve_vdisp(sspflux * atten[:, np.newaxis], self.pixkms, self.vdisp_nominal)
            sspflux_dustvdisp.append(_sspflux_dustvdisp)
        sspflux_dustvdisp = np.stack(sspflux_dustvdisp, axis=-1) # [npix,nage,nAV]

        self.sspwave = sspwave
        self.sspflux = sspflux                     # no dust, no velocity broadening [npix,nage]
        self.sspflux_dustvdisp = sspflux_dustvdisp # nominal velocity broadening on a grid of A(V) [npix,nage,nAV]
        self.sspinfo = sspinfo
        self.nage = nage
        self.npix = npix

        # table of emission lines to fit
        self.linetable = read_nyxgalaxy_lines()
        self.linemask_sigma = 150.0 # [km/s]

        # Ddo a throw-away trapezoidal resampling so we can compile the numba
        # code when instantiating this class.
        #t0 = time.time()
        _ = trapz_rebin(np.arange(4), np.ones(4), np.arange(2)+1)
        #print('Initial rebin ', time.time() - t0)

    def init_output(self, nobj=1):
        """Initialize the output data table for this class.

        """
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

        # continuum fit with *no* dust reddening (to be used as a diagnostic
        # tool to identify potential calibration issues).
        out.add_column(Column(name='CONTINUUM_NODUST_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_NODUST_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_NODUST_AGE', length=nobj, dtype='f4', unit=u.Gyr))

        out.add_column(Column(name='CONTINUUM_PHOT_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_PHOT_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_PHOT_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_PHOT_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_PHOT_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_PHOT_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))

        out.add_column(Column(name='D4000', length=nobj, dtype='f4'))
        out.add_column(Column(name='D4000_IVAR', length=nobj, dtype='f4'))
        out.add_column(Column(name='D4000_NOLINES', length=nobj, dtype='f4'))
        out.add_column(Column(name='D4000_MODEL', length=nobj, dtype='f4'))
        out.add_column(Column(name='D4000_MODEL_PHOT', length=nobj, dtype='f4'))

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

    def obsolete_convolve_vdisp(self, sspflux, vdisp):
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
                 datasspflux = self.obsolete_convolve_vdisp(datasspflux, vdisp)
                 
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
                    _datasspflux = self.obsolete_convolve_vdisp(_datasspflux, vdisp)
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
            log.fatal('Coefficients are all zero!')
            raise ValueError
        
        meanage = np.sum(coeff * age) / np.sum(coeff) / 1e9 # [Gyr]
        
        return meanage

    def obsolete_get_d4000(self, wave, flam, flam_ivar=None, redshift=None, rest=True):
        """Compute D(4000) and, optionally, the inverse variance.

        Parameters
        ----------
        wave
        flam
        flam_ivar
        redshift
        rest

        Returns
        -------

        Notes
        -----
        If `rest`=``False`` then `redshift` input is required.
        
        """
        d4000, d4000_ivar = 0.0, 0.0

        if rest:
            flam2fnu =  wave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
        else:
            wave /= (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * wave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]

        if flam_ivar is None:
            goodmask = np.ones(len(flam), bool) # True is good
        else:
            goodmask = flam_ivar > 0
            
        indxblu = np.where((wave >= 3850.) * (wave <= 3950.) * goodmask)[0]
        indxred = np.where((wave >= 4000.) * (wave <= 4100.) * goodmask)[0]
        if len(indxblu) < 5 or len(indxred) < 5:
            return d4000, d4000_ivar

        blufactor, redfactor = 3950.0 - 3850.0, 4100.0 - 4000.0
        deltawave = np.gradient(wave) # should be constant...

        fnu = flam * flam2fnu # [erg/s/cm2/Hz]

        numer = blufactor * np.sum(deltawave[indxred] * fnu[indxred])
        denom = redfactor * np.sum(deltawave[indxblu] * fnu[indxblu])
        if denom == 0.0:
            log.warning('D(4000) is ill-defined!')
            return d4000, d4000_ivar
        d4000 =  numer / denom
        
        if flam_ivar is not None:
            fnu_ivar = flam_ivar / flam2fnu**2
            fnu_var, _ = _ivar2var(fnu_ivar)
            
            numer_var = blufactor**2 * np.sum(deltawave[indxred] * fnu_var[indxred])
            denom_var = redfactor**2 * np.sum(deltawave[indxblu] * fnu_var[indxblu])
            d4000_var = (numer_var + numer**2 * denom_var) / denom**2
            if d4000_var <= 0:
                log.warning('D(4000) variance is ill-defined!')
                d4000_ivar = 0.0
            else:
                d4000_ivar = 1.0 / d4000_var

        return d4000, d4000_ivar

    def younger_than_universe(self, redshift):
        """Return the indices of the SSPs younger than the age of the universe at the
        given redshift.

        """
        return np.where(self.sspinfo['age'] <= self.cosmo.age(redshift).to(u.year).value)[0]

    def fnnls_continuum_bestfit(self, coeff, sspflux=None, specwave=None,
                                specres=None, redshift=None, south=True):
        if sspflux is None:
            sspflux, sspphot = self.redshift_smooth_and_resample(redshift, specwave, specres, south=south)
            if specres is None and specwave is None: # ignore per-camera
                bestfit = sspflux.dot(coeff)
            else: # iterate over camera
                bestfit = [_sspflux.dot(coeff) for _sspflux in sspflux]
        else:
            bestfit = sspflux.dot(coeff)

        if specres is None and specwave is None:
            return bestfit, self.sspwave
        else:
            return bestfit

    @staticmethod
    def find_chi2min(xx, chi2):
        """Find the chi2 minimum and fit a parabola centered on it. 

        Parameters
        ----------
        xx : :class:`numpy.ndarray`
            Independent variable of values (e.g., reddening).
        chi2 : :class:`numpy.ndarray`
            Array of chi2 values corresponding to each value of `xx`. 

        Returns
        -------
        :class:`tuple` with three `float` elements
            Chi2 minimum, value at that minimum, and the inverse variance. 

        Notes
        -----
        We minimize chi2 by fitting a parabola to the three points around the
        minimum, a standard approach. This particular algorithm was canibalized
        largely from `redrock.fitz.fitminimum`, written by S. Bailey.

        Poor chi2 minima are flagged with an inverse variance set to zero.

        """
        if len(xx) < 3:
            return (-1.0, -1.0, 0.0)
        
        mindx = np.argmin(chi2)
        if mindx == 0 or mindx == len(chi2)-1:
            return (-1.0, -1.0, 0.0)
        
        # model: y = aa*x**2 + bb*x + cc
        try:
            aa, bb, cc = np.polyfit(xx[mindx-1:mindx+2], chi2[mindx-1:mindx+2], 2)
        except np.linalg.LinAlgError:
            return (-1.0, -1.0, 0.0)

        if aa == 0.0:
            return (-1.0, -1.0, 0.0)
        
        # recast as y = y0 + ((x-x0)/xerr)^2
        xbest = -bb / (2*aa)
        chi2min = -(bb**2) / (4*aa) + cc

        # xivar==0 means a bad fit
        xivar = aa # err = 1/np.sqrt(aa)
        if (xbest <= np.min(xx)) or (np.max(xx) <= xbest):
            xivar = 0.0
        if (chi2min <= 0.):
            xivar = 0.0
        if aa <= 0.0:
            xivar = 0.0

        return (chi2min, xbest, xivar)

    def fnnls_continuum(self, data):
        """Fit the continuum using fast non-negative least-squares fitting (fNNLS).

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            `unpack_one_spectrum`.
        sigma_mask : :class:`float`, optional
            Mask all pixels within +/-1.5`sigma_mask` [km/s] of known emission
            lines during continuum-fitting. Defaults to 300 km/s.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in `init_nyxgalaxy`.

        Notes
        -----
        See https://github.com/jvendrow/fnnls for the fNNLS algorithm.

        ToDo:
          - Use cross-correlation to update the redrock redshift.
          - Need to mask more emission lines than we fit (e.g., Mg II).

        """
        from fnnls import fnnls
        from redrock import fitz

        def _fnnls_parallel(modelflux, flux, ivar, xparam=None, debug=False):
            """Wrapper on fnnls to set up the multiprocessing. Works with both spectroscopic
            and photometric input and with both 2D and 3D model spectra.

            """
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
                warn, coeff, chi2 = fit_fnnls_continuum(ZZ, xx, flux=flux, ivar=ivar,
                                                        modelflux=modelflux, get_chi2=True)
                if np.any(warn):
                    print('WARNING: fnnls did not converge after 5 iterations.')

                return coeff, chi2

            # ...otherwise multiprocess over the xparam (e.g., AV or vdisp)
            # dimension.
            ZZ = modelflux * ww[:, np.newaxis, np.newaxis] # reshape into [npix/nband,nage,nAV/nvdisp]

            fitargs = [(ZZ[:, :, ii], xx, flux, ivar, modelflux[:, :, ii], None, True, xparam[ii]) for ii in np.arange(nn)]
            if self.nproc > 1:
                with multiprocessing.Pool(self.nproc) as P:
                    rr = P.map(_fit_fnnls_continuum, fitargs)
            else:
                #fit_fnnls_continuum(*fitargs[10])
                #pdb.set_trace()
                rr = [fit_fnnls_continuum(*_fitargs) for _fitargs in fitargs]
            warn, _, chi2grid = list(zip(*rr)) # unpack
            if np.any(warn):
                vals = ','.join(['{:.1f}'.format(xp) for xp in xparam[np.where(warn)[0]]])
                log.warning('fnnls did not converge after 5 iterations for parameter value(s) {}.'.format(vals))
            chi2grid = np.array(chi2grid)

            imin = fitz.find_minima(chi2grid)[0]
            xbest, xerr, chi2min, warn = fitz.minfit(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2])
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

        # Initialize the output table; see init_nyxgalaxy for the data model.
        result = self.init_output()

        redshift = data['zredrock']
        result['CONTINUUM_Z'] = redshift
        result['CONTINUUM_SNR'] = data['snr']

        # Prepare the reddened and unreddened SSP templates. Note that we ignore
        # templates which are older than the age of the universe at the galaxy
        # redshift.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustvdisp, zsspphot_dustvdisp = self.SSP2data(
            self.sspflux_dustvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=data['wave'], specres=data['res'],
            south=data['photsys_south'], nproc=1)
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

        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = data['phot']['flam_ivar'].data / self.fluxnorm**2
        zsspflam_dustvdisp = zsspphot_dustvdisp['flam'].data * self.fluxnorm * self.massnorm # [nband,nage*nAV]
        assert(np.all(objflamivar >= 0))

        inodust = np.asscalar(np.where(self.AV == 0)[0]) # should always be index 0

        npix, nmodel = zsspflux_dustvdisp.shape
        nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)

        zsspflux_dustvdisp = zsspflux_dustvdisp.reshape(npix, nage, self.nAV)       # [npix,nage,nAV]
        zsspflam_dustvdisp = zsspflam_dustvdisp.reshape(self.nband, nage, self.nAV) # [nband,nage,nAV]
        
        # [1] First, fit just the photometry, independent of the spectra,
        # including reddening and at our nominal velocity dispersion.
        t0 = time.time()
        AVchi2min, AVbest, AVivar = _fnnls_parallel(zsspflam_dustvdisp, objflam, objflamivar, xparam=self.AV)
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
        coeff, chi2min = _fnnls_parallel(bestphot['flam'].data*self.massnorm*self.fluxnorm,
                                         objflam, objflamivar) # bestphot['flam'] is [nband, nage]

        bestfit = bestsspflux.dot(coeff)
        d4000, _ = get_d4000(self.sspwave, bestfit, rest=True)
        #d4000, _ = self.get_d4000(self.sspwave, bestfit, rest=True)
        meanage = self.get_meanage(coeff)
        log.info('Photometric D(4000)={:.3f}, Age={:.2f} Gyr'.format(d4000, meanage))

        result['CONTINUUM_PHOT_COEFF'][0][:nage] = coeff
        result['CONTINUUM_PHOT_CHI2'][0] = chi2min
        result['CONTINUUM_PHOT_AGE'][0] = meanage
        result['CONTINUUM_PHOT_AV'][0] = AVbest
        result['CONTINUUM_PHOT_AV_IVAR'][0] = AVivar
        result['D4000_MODEL_PHOT'][0] = d4000
            
        # [2] Next, fit the spectra  with *no* dust reddening so we can identify
        # potential calibration issues (again, at the nominal velocity
        # dispersion).
        t0 = time.time()
        coeff, chi2min = _fnnls_parallel(zsspflux_dustvdisp[:, :, inodust], specflux, specivar)
        log.info('No-dust model fit has chi2={:.3f} and took {:.2f} sec'.format(
            chi2min, time.time()-t0))

        result['CONTINUUM_NODUST_COEFF'][0][0:nage] = coeff
        result['CONTINUUM_NODUST_CHI2'] = chi2min

        # [3] Third, fit the spectra for reddening using the models convolved to
        # the nominal velocity dispersion and then fit for velocity dispersion.
        t0 = time.time()
        AVchi2min, AVbest, AVivar = _fnnls_parallel(zsspflux_dustvdisp, specflux, specivar,
                                                    xparam=self.AV, debug=False)
        log.info('Fitting for the reddening took: {:.2f} sec'.format(time.time()-t0))
        if AVivar > 0:
            log.info('Best-fitting spectroscopic A(V)={:.4f}+/-{:.4f} with chi2={:.3f}'.format(
                AVbest, 1/np.sqrt(AVivar), AVchi2min))
        else:
            AVbest = self.AV_nominal
            log.info('Finding spectroscopic A(V) failed; adopting A(V)={:.4f}'.format(
                self.AV_nominal))
        
        # Build out the model spectra on our grid of velocity dispersion
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
        vdispchi2min, vdispbest, vdispivar = _fnnls_parallel(zsspflux_vdisp, specflux, specivar,
                                                             xparam=self.vdisp, debug=False)
        log.info('Fitting for the velocity dispersion took: {:.2f} sec'.format(time.time()-t0))
        if vdispivar > 0:
            log.info('Best-fitting vdisp={:.2f}+/-{:.2f} km/s with chi2={:.3f}'.format(
                vdispbest, 1/np.sqrt(vdispivar), vdispchi2min))
        else:
            vdispbest = self.vdisp_nominal
            log.info('Finding vdisp failed; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))

        # Get the final set of coefficients and chi2 at the best-fitting
        # reddening and velocity dispersion.
        bestsspflux, bestphot = self.SSP2data(self.sspflux[:, agekeep], self.sspwave,
                                              specwave=data['wave'], specres=data['res'],
                                              AV=AVbest, vdisp=vdispbest, redshift=redshift,
                                              south=data['photsys_south'])
        bestsspflux = np.concatenate(bestsspflux, axis=0)
        coeff, chi2min = _fnnls_parallel(bestsspflux, specflux, specivar)

        bestfit = bestsspflux.dot(coeff)
        meanage = self.get_meanage(coeff)
        d4000_model, _ = get_d4000(specwave, bestfit, redshift=redshift)
        d4000, d4000_ivar = get_d4000(specwave, specflux, specivar, redshift=redshift)
        #d4000_model, _ = self.get_d4000(specwave, bestfit, redshift=redshift)
        #d4000, d4000_ivar = self.get_d4000(specwave, specflux, specivar, redshift=redshift)
        log.info('Spectroscopic D(4000)={:.3f}, Age={:.2f} Gyr'.format(d4000, meanage))

        #import matplotlib.pyplot as plt
        #plt.clf()
        #plt.plot(specwave, specflux)
        #plt.plot(specwave, bestfit)
        #plt.savefig('qa-bestfit.png')

        result['CONTINUUM_COEFF'][0][0:nage] = coeff
        result['CONTINUUM_CHI2'][0] = chi2min
        result['CONTINUUM_AV'][0] = AVbest
        result['CONTINUUM_AV_IVAR'][0] = AVivar
        result['CONTINUUM_VDISP'][0] = vdispbest
        result['CONTINUUM_VDISP_IVAR'][0] = vdispivar
        result['CONTINUUM_AGE'] = meanage
        result['D4000'][0] = d4000
        result['D4000_IVAR'][0] = d4000_ivar
        result['D4000_MODEL'][0] = d4000_model

        # Unpack the continuum into individual cameras.
        continuum = []
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npixpercam[:ii+1])
            jpix = np.sum(npixpercam[:ii+2])
            continuum.append(bestfit[ipix:jpix])

        return result, continuum
    
    def fnnls_continuum_plot(self, data, result, qadir='.'):
        """QA of the best-fitting continuum.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

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
                                     AV=result['CONTINUUM_AV'],
                                     vdisp=result['CONTINUUM_VDISP'],
                                     coeff=result['CONTINUUM_COEFF'],
                                     synthphot=False)
        continuum_phot, bestphot = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift,
                                                 AV=result['CONTINUUM_PHOT_AV'],
                                                 coeff=result['CONTINUUM_PHOT_COEFF'] * self.massnorm,
                                                 south=data['photsys_south'],
                                                 synthphot=True)
        continuum_wave_phot = self.sspwave * (1 + redshift)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
        for ii in [0, 1, 2]: # iterate over cameras
            sigma, _ = _ivar2var(data['ivar'][ii], sigma=True)

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
            'targetid': 'targetid={} fiber={}'.format(result['TARGETID'], result['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(result['CONTINUUM_CHI2']),
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(result['Z']),
            #'znyxgalaxy': '$z_{{\\rm nyxgalaxy}}$={:.6f}'.format(result['CONTINUUM_Z']),
            'z': '$z$={:.6f}'.format(result['CONTINUUM_Z']),
            'age': '<Age>={:.3f} Gyr'.format(result['CONTINUUM_AGE']),
            }
        if result['CONTINUUM_VDISP_IVAR'] == 0:
            leg.update({'vdisp': '$\sigma$={:.1f} km/s'.format(result['CONTINUUM_VDISP'])})
        else:
            leg.update({'vdisp': '$\sigma$={:.1f}+/-{:.1f} km/s'.format(
                result['CONTINUUM_VDISP'], 1/np.sqrt(result['CONTINUUM_VDISP_IVAR']))})
        if result['CONTINUUM_AV_IVAR'] == 0:
            leg.update({'AV': '$A(V)$={:.3f} mag'.format(result['CONTINUUM_AV'])})
        else:
            leg.update({'AV': '$A(V)$={:.3f}+/-{:.3f} mag'.format(
                result['CONTINUUM_AV'], 1/np.sqrt(result['CONTINUUM_AV_IVAR']))})

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
            result['TILE'], result['NIGHT'], result['TARGETID']))
        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)

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

def read_nyxgalaxy_lines():
    """Read the set of emission lines of interest.

    ToDo: add lines to mask during continuum-fitting but which we do not want to
    emission-line fit.

    """
    from pkg_resources import resource_filename
    
    linefile = resource_filename('desigal', 'data/nyxgalaxy_lines.ecsv')    
    linetable = Table.read(linefile, format='ascii.ecsv', guess=False)
    
    return linetable    

class EMLineModel(Fittable1DModel):
    """Class to model the emission-line spectra.

    """
    from astropy.modeling import Parameter

    # NB! The order of the parameters here matters! linevshift is the shift of
    # the emission-line velocity (in km/s) with respect to the systemic redshift
    linevshift_forbidden = Parameter(name='linevshift_forbidden', default=0.0, bounds=(-500.0, 500.0)) # [km/s]
    linevshift_balmer = Parameter(name='linevshift_balmer', default=0.0, bounds=(-500.0, 500.0)) # [km/s]

    linesigma_forbidden = Parameter(name='linesigma_forbidden', default=50.0, bounds=(5, 350)) # line-sigma [km/s]
    linesigma_balmer = Parameter(name='linesigma_balmer', default=50.0, bounds=(5, 350)) # line-sigma [km/s]

    # Fragile because the lines are hard-coded--
    oii_3726_amp = Parameter(name='oii_3726_amp', default=1.0)
    oii_3729_amp = Parameter(name='oii_3729_amp', default=1.0)
    oiii_4959_amp = Parameter(name='oiii_4959_amp', default=1.0)
    oiii_5007_amp = Parameter(name='oiii_5007_amp', default=3.0)
    nii_6548_amp = Parameter(name='nii_6548_amp', default=1.0)
    nii_6584_amp = Parameter(name='nii_6584_amp', default=3.0)
    sii_6716_amp = Parameter(name='sii_6716_amp', default=1.0)
    sii_6731_amp = Parameter(name='sii_6731_amp', default=1.0)
    hepsilon_amp = Parameter(name='hepsilon_amp', default=0.5)
    hdelta_amp = Parameter(name='hdelta_amp', default=0.5)
    hgamma_amp = Parameter(name='hgamma_amp', default=0.5)
    hbeta_amp = Parameter(name='hbeta_amp', default=1.0)
    halpha_amp = Parameter(name='halpha_amp', default=3.0)

    # tie the velocity shifts and line-widths together
    def tie_vshift(model):
        return model.linevshift_balmer
    linevshift_forbidden.tied = tie_vshift

    def tie_sigma(model):
        return model.linesigma_balmer
    linesigma_forbidden.tied = tie_sigma

    # tie the [NII] and [OIII] line-strengths together
    def tie_oiii(model):
        return model.oiii_5007_amp / 2.8875
    oiii_4959_amp.tied = tie_oiii

    def tie_nii(model):
        return model.nii_6584_amp / 2.936
    nii_6548_amp.tied = tie_nii
    
    def __init__(self,
                 linevshift_forbidden=linevshift_forbidden.default,
                 linevshift_balmer=linevshift_balmer.default,
                 linesigma_forbidden=linesigma_forbidden.default,
                 linesigma_balmer=linesigma_balmer.default,
                 oii_3726_amp=oii_3726_amp.default, 
                 oii_3729_amp=oii_3729_amp.default, 
                 oiii_4959_amp=oiii_4959_amp.default, 
                 oiii_5007_amp=oiii_5007_amp.default, 
                 nii_6548_amp=nii_6548_amp.default, 
                 nii_6584_amp=nii_6584_amp.default, 
                 sii_6716_amp=sii_6716_amp.default, 
                 sii_6731_amp=sii_6731_amp.default, 
                 hepsilon_amp=hepsilon_amp.default, 
                 hdelta_amp=hdelta_amp.default, 
                 hgamma_amp=hgamma_amp.default, 
                 hbeta_amp=hbeta_amp.default, 
                 halpha_amp=halpha_amp.default,
                 redshift=None,
                 emlineR=None, npixpercamera=None,
                 log10wave=None, **kwargs):
        """Initialize the emission-line model.
        
        emlineR -
        redshift - required keyword
        
        """
        self.linetable = read_nyxgalaxy_lines()
        self.redshift = redshift
        self.emlineR = emlineR
        self.npixpercamera = np.hstack([0, npixpercamera])

        # internal wavelength vector for building the emission-line model
        if log10wave is None:
            pixkms = 10.0
            dlogwave = pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
            log10wave = np.arange(np.log10(3000), np.log10(1e4), dlogwave)
        self.log10wave = log10wave
            
        super(EMLineModel, self).__init__(
            linevshift_forbidden=linevshift_forbidden,
            linevshift_balmer=linevshift_balmer,
            linesigma_forbidden=linesigma_forbidden,
            linesigma_balmer=linesigma_balmer,
            oii_3726_amp=oii_3726_amp,
            oii_3729_amp=oii_3729_amp,
            oiii_4959_amp=oiii_4959_amp,
            oiii_5007_amp=oiii_5007_amp,
            nii_6548_amp=nii_6548_amp,
            nii_6584_amp=nii_6584_amp,
            sii_6716_amp=sii_6716_amp,
            sii_6731_amp=sii_6731_amp,
            hepsilon_amp=hepsilon_amp,
            hdelta_amp=hdelta_amp,
            hgamma_amp=hgamma_amp,
            hbeta_amp=hbeta_amp,
            halpha_amp=halpha_amp, **kwargs)

    def evaluate(self, emlinewave, *args):
        """Evaluate the emission-line model.
        emlineR=None, npixpercamera=None, 

        """ 
        linevshift_forbidden, linevshift_balmer = args[0], args[1]
        linez_forbidden = self.redshift + linevshift_forbidden / C_LIGHT
        linez_balmer = self.redshift + linevshift_balmer / C_LIGHT

        linesigma_forbidden, linesigma_balmer = args[2], args[3]
        log10sigma_forbidden = linesigma_forbidden / C_LIGHT / np.log(10) # line-width [log-10 Angstrom]
        log10sigma_balmer = linesigma_balmer / C_LIGHT / np.log(10)       # line-width [log-10 Angstrom]

        lineamps = args[4:]
        linenames = self.linetable['name'].data
        isbalmers = self.linetable['isbalmer'].data

        # build the emission-line model [erg/s/cm2/A, observed frame]; should we
        # multiprocess this step?
        log10model = np.zeros_like(self.log10wave)
        for linename, lineamp, isbalmer in zip(linenames, lineamps, isbalmers):
            restlinewave = self.linetable[self.linetable['name'] == linename]['restwave'][0]
            if isbalmer:
                log10sigma = log10sigma_balmer
                linezwave = np.log10(restlinewave * (1.0 + linez_balmer)) # redshifted wavelength [log-10 Angstrom]
            else:
                log10sigma = log10sigma_forbidden
                linezwave = np.log10(restlinewave * (1.0 + linez_forbidden)) # redshifted wavelength [log-10 Angstrom]

            ww = np.abs(self.log10wave - linezwave) < 20 * log10sigma
            if np.count_nonzero(ww) > 0:
                #log.info(linename, 10**linezwave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
                log10model[ww] += lineamp * np.exp(-0.5 * (self.log10wave[ww]-linezwave)**2 / log10sigma**2)

        # split into cameras, resample, and convolve with the instrumental
        # resolution
        emlinemodel = []
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(self.npixpercamera[:ii+1])
            jpix = np.sum(self.npixpercamera[:ii+2])
            #_emlinemodel = resample_flux(emlinewave[ipix:jpix], 10**self.log10wave, log10model)
            _emlinemodel = trapz_rebin(10**self.log10wave, log10model, emlinewave[ipix:jpix])
            
            if self.emlineR is not None:
                _emlinemomdel = self.emlineR[ii].dot(_emlinemodel)
            
            #plt.plot(10**_emlinewave, _emlinemodel)
            #plt.plot(10**_emlinewave, self.emlineR[ii].dot(_emlinemodel))
            #plt.xlim(3870, 3920) ; plt.show()
            #pdb.set_trace()
            emlinemodel.append(_emlinemodel)
            
        return np.hstack(emlinemodel)

class EMLineFit(object):
    """Class to fit an emission-line spectrum.

    * https://docs.astropy.org/en/stable/modeling/example-fitting-constraints.html#tied
    * https://docs.astropy.org/en/stable/modeling/new-model.html
    * https://docs.astropy.org/en/stable/modeling/compound-models.html#parameters

    """
    def __init__(self, nball=10, chi2fail=1e6, nproc=1, verbose=False):
        """Write me.
        
        """
        from astropy.modeling import fitting

        self.verbose = verbose
        self.nproc = nproc
        
        self.nball = nball
        self.chi2fail = chi2fail
        self.pixkms = 10.0 # pixel size for internal wavelength array [km/s]

        self.fitter = fitting.LevMarLSQFitter()

    def init_output(self, linetable, nobj=1):
        """Initialize the output data table for this class.

        """
        out = Table()
        out.add_column(Column(name='LINEVSHIFT_FORBIDDEN', length=nobj, dtype='f4'))
        out.add_column(Column(name='LINEVSHIFT_FORBIDDEN_IVAR', length=nobj, dtype='f4'))
        out.add_column(Column(name='LINEVSHIFT_BALMER', length=nobj, dtype='f4'))
        out.add_column(Column(name='LINEVSHIFT_BALMER_IVAR', length=nobj, dtype='f4'))
        out.add_column(Column(name='LINESIGMA_FORBIDDEN', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='LINESIGMA_FORBIDDEN_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        out.add_column(Column(name='LINESIGMA_BALMER', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='LINESIGMA_BALMER_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))

        for line in linetable['name']:
            line = line.upper()
            out.add_column(Column(name='{}_AMP'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))
            out.add_column(Column(name='{}_AMP_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2))
            out.add_column(Column(name='{}_FLUX'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2)))
            out.add_column(Column(name='{}_FLUX_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4/u.erg**2))
            out.add_column(Column(name='{}_BOXFLUX'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2)))
            out.add_column(Column(name='{}_BOXFLUX_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4/u.erg**2))
            out.add_column(Column(name='{}_CONT'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))
            out.add_column(Column(name='{}_CONT_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2))
            out.add_column(Column(name='{}_EW'.format(line), length=nobj, dtype='f4',
                                  unit=u.Angstrom))
            out.add_column(Column(name='{}_EW_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=1/u.Angstrom**2))
            out.add_column(Column(name='{}_FLUX_LIMIT'.format(line), length=nobj, dtype='f4',
                                  unit=u.erg/(u.second*u.cm**2)))
            out.add_column(Column(name='{}_EW_LIMIT'.format(line), length=nobj, dtype='f4',
                                  unit=u.Angstrom))
            out.add_column(Column(name='{}_CHI2'.format(line), length=nobj, dtype='f4'))
            out.add_column(Column(name='{}_NPIX'.format(line), length=nobj, dtype=np.int32))

        return out
        
    def chi2(self, bestfit, emlinewave, emlineflux, emlineivar):
        """Compute the reduced chi^2."""
        dof = len(emlinewave) - len(bestfit.parameters)
        emlinemodel = bestfit(emlinewave)
        chi2 = np.sum(emlineivar * (emlineflux - emlinemodel)**2) / dof
        return chi2

    def emlinemodel_bestfit(self, specwave, specres, nyxgalaxy_table):
        """Wrapper function to get the best-fitting emission-line model
        from an nyxgalaxy table (to be used to build QA).

        """
        npixpercamera = [len(gw) for gw in specwave]

        redshift = nyxgalaxy_table['Z']
        linesigma_forbidden = nyxgalaxy_table['LINESIGMA_FORBIDDEN']
        linesigma_balmer = nyxgalaxy_table['LINESIGMA_BALMER']

        linevshift_forbidden = nyxgalaxy_table['LINEVSHIFT_FORBIDDEN']
        linevshift_balmer = nyxgalaxy_table['LINEVSHIFT_BALMER']
        #linez_forbidden = nyxgalaxy_table['LINEZ_FORBIDDEN']
        #linez_balmer = nyxgalaxy_table['LINEZ_BALMER']
        #linevshift_forbidden = (linez_forbidden - redshift) * C_LIGHT # [km/s]
        #linevshift_balmer = (linez_balmer - redshift) * C_LIGHT # [km/s]

        EMLine = EMLineModel(linevshift_forbidden=linevshift_forbidden,
                             linevshift_balmer=linevshift_balmer,
                             linesigma_forbidden=linesigma_forbidden,
                             linesigma_balmer=linesigma_balmer,
                             redshift=redshift, emlineR=specres,
                             npixpercamera=npixpercamera)
        # skip linevshift_[forbidden,balmer] and linesigma_[forbidden,balmer]
        lineargs = [nyxgalaxy_table[linename.upper()] for linename in EMLine.param_names[4:]] 
        lineargs = [linevshift_forbidden, linevshift_balmer, linesigma_forbidden, linesigma_balmer] + lineargs

        _emlinemodel = EMLine.evaluate(np.hstack(specwave), *lineargs)

        # unpack it
        emlinemodel = []
        npix = np.hstack([0, npixpercamera])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            emlinemodel.append(_emlinemodel[ipix:jpix])

        return emlinemodel
    
    def fit(self, data, continuum):
        """Perform the fit minimization / chi2 minimization.
        
        EMLineModel object
        FC - ContinuumFit object

        need to take into account the instrumental velocity width when computing integrated fluxes
        
        """
        #from scipy import integrate
        from astropy.stats import sigma_clipped_stats
        
        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        npixpercamera = [len(gw) for gw in data['wave']]
        npixpercam = np.hstack([0, npixpercamera])

        redshift = data['zredrock']
        emlinewave = np.hstack(data['wave'])
        emlineivar = np.hstack(data['ivar'])
        specflux = np.hstack(data['flux'])
        emlineflux = specflux - np.hstack(continuum)

        dlogwave = self.pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        log10wave = np.arange(np.log10(3e3), np.log10(1e4), dlogwave)
        #log10wave = np.arange(np.log10(emlinewave.min()), np.log10(emlinewave.max()), dlogwave)
        
        self.EMLineModel = EMLineModel(redshift=redshift,
                                       emlineR=data['res'],
                                       npixpercamera=npixpercamera,
                                       log10wave=log10wave)
        nparam = len(self.EMLineModel.parameters)

        #weights = np.ones_like
        #emlinevar, _ = _ivar2var(emlineivar)
        #weights[np.logical_not(goodmask)] = 1e16
        bestfit = self.fitter(self.EMLineModel, emlinewave, emlineflux, weights=np.sqrt(emlineivar))

        emlinemodel = bestfit(emlinewave)
        chi2 = self.chi2(bestfit, emlinewave, emlineflux, emlineivar).astype('f4')

        # Initialize the output table; see init_nyxgalaxy for the data model.
        result = self.init_output(self.EMLineModel.linetable)

        specflux_nolines = specflux - emlinemodel

        d4000_nolines, _ = get_d4000(emlinewave, specflux_nolines, redshift=redshift)
        result['D4000_NOLINES'] = d4000_nolines

        ## Pack the results in a dictionary and return.
        ## https://gist.github.com/eteq/1f3f0cec9e4f27536d52cd59054c55f2
        #result = {
        #    'converged': False,
        #    'fit_message': self.fitter.fit_info['message'],
        #    'nparam': nparam,
        #    'npix': len(emlinewave),
        #    'dof': len(emlinewave) - len(self.EMLineModel.parameters),
        #    'chi2': chi2,
        #    'linenames': [ll.replace('_amp', '') for ll in self.EMLineModel.param_names[4:]],
        #}
        #for param in bestfit.param_names:
        #    result.update({param: getattr(bestfit, param).value})
        
        # uncertainties
        if self.fitter.fit_info['param_cov'] is not None:
            cov = self.fitter.fit_info['param_cov']
            ivar = 1 / np.diag(cov)
            #result['converged'] = True
        else:
            cov = np.zeros((nparam, nparam))
            ivar = np.zeros(nparam)

        # Need to be careful about uncertainties for tied parameters--
        # https://github.com/astropy/astropy/issues/7202
        
        #err_params = np.sqrt(np.diag(fitter.fit_info['param_cov']))
        #err = model.copy()
        #fitting._fitter_to_model_params(err, err_params)            
            
        count = 0
        for ii, pp in enumerate(bestfit.param_names):
            pinfo = getattr(bestfit, pp)
            iinfo = getattr(self.EMLineModel, pp)

            # need to think about this more deeply
            #if pinfo.value == iinfo.value: # not fitted
            #    result.update({pinfo.name: np.float(0.0)})
            #else:
            #    result.update({pinfo.name: pinfo.value.astype('f4')})
            #result.update({pinfo.name: pinfo.value.astype('f4')})
            result[pinfo.name.upper()] = pinfo.value.astype('f4')
                
            if pinfo.fixed:
                #result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
                result['{}_IVAR'.format(pinfo.name.upper())] = pinfo.value.astype('f4')
            elif pinfo.tied:
                # hack! see https://github.com/astropy/astropy/issues/7202
                #result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
                result['{}_IVAR'.format(pinfo.name.upper())] = np.float32(0.0)
            else:
                result['{}_IVAR'.format(pinfo.name.upper())] = ivar[count].astype('f4')
                #result.update({'{}_ivar'.format(pinfo.name): ivar[count].astype('f4')})
                count += 1

            #if 'forbidden' in pinfo.name:
            #    pdb.set_trace()

        # hack for tied parameters---gotta be a better way to do this
        #result['oiii_4959_amp_ivar'] = result['oiii_5007_amp_ivar'] * 2.8875**2
        result['OIII_4959_AMP_IVAR'] = result['OIII_5007_AMP_IVAR'] * 2.8875**2
        result['NII_6548_AMP_IVAR'] = result['NII_6548_AMP_IVAR'] * 2.936**2
        result['LINEVSHIFT_FORBIDDEN_IVAR'] = result['LINEVSHIFT_BALMER_IVAR']
        result['LINESIGMA_FORBIDDEN_IVAR'] = result['LINESIGMA_BALMER_IVAR']

        ## convert the vshifts to redshifts
        #result['linez_forbidden'] = redshift + result['linevshift_forbidden'] / C_LIGHT
        #result['linez_balmer'] = redshift + result['linevshift_balmer'] / C_LIGHT
        #result['linez_forbidden_ivar'] = result['linevshift_forbidden_ivar'] * C_LIGHT**2
        #result['linez_balmer_ivar'] = result['linevshift_balmer_ivar'] * C_LIGHT**2

        # now loop back through and if ivar==0 then set the parameter value to zero
        if False:
            if self.fitter.fit_info['param_cov'] is not None:
                for pp in bestfit.param_names:
                    if result['{}_IVAR'.format(pp.upper())] == 0.0:
                        result[pp] = np.float(0.0)

        # get continuum fluxes, EWs, and upper limits
        sigma_cont = 150.0
        for oneline in self.EMLineModel.linetable:
            line = oneline['name'].upper()

            # get the emission-line flux
            zwave = oneline['restwave'] * (1 + redshift)
            if oneline['isbalmer']:
                linesigma = result['LINESIGMA_BALMER']
            else:
                linesigma = result['LINESIGMA_FORBIDDEN']

            linesigma_ang = zwave * linesigma / C_LIGHT # [observed-frame Angstrom]
            norm = np.sqrt(2.0 * np.pi) * linesigma_ang

            result['{}_FLUX'.format(line)] = result['{}_AMP'.format(line)] * norm
            result['{}_FLUX_IVAR'.format(line)] = result['{}_AMP_IVAR'.format(line)] / norm**2

            # boxcar integration, chi2, and number of pixels
            lineindx = np.where((emlinewave > (zwave - 3*sigma_cont * zwave / C_LIGHT)) *
                                (emlinewave < (zwave + 3.*sigma_cont * zwave / C_LIGHT)) *
                                (emlineivar > 0))[0]
            npix = len(lineindx)
            if npix > 0:
                dof = npix - 3 # ??? [redshift, sigma, and amplitude]
                chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-emlinemodel[lineindx])**2) / dof
                boxflux = np.sum(emlineflux[lineindx])
                boxflux_ivar = 1 / np.sum(1 / emlineivar[lineindx])
            else:
                npix, chi2, boxflux, boxflux_ivar = 0.0, 1e6, 0.0, 0.0

            result['{}_NPIX'.format(line)] = npix
            result['{}_CHI2'.format(line)] = chi2
            result['{}_BOXFLUX'.format(line)] = boxflux
            result['{}_BOXFLUX_IVAR'.format(line)] = boxflux_ivar
            
            # get the continuum and EWs
            indxlo = np.where((emlinewave > (zwave - 10*sigma_cont * zwave / C_LIGHT)) *
                              (emlinewave < (zwave - 3.*sigma_cont * zwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indxhi = np.where((emlinewave < (zwave + 10*sigma_cont * zwave / C_LIGHT)) *
                              (emlinewave > (zwave + 3.*sigma_cont * zwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indx = np.hstack((indxlo, indxhi))

            if len(indx) > 5: # require at least 5 pixels to get the continuum level
                _, cmed, csig = sigma_clipped_stats(specflux_nolines[indx], sigma=3.0)
                civar = (np.sqrt(len(indx)) / csig)**2
            else:
                cmed, civar = 0.0, 0.0
                
            result['{}_CONT'.format(line)] = cmed
            result['{}_CONT_IVAR'.format(line)] = civar

            if result['{}_CONT_IVAR'.format(line)] != 0.0:
                factor = (1 + redshift) / result['{}_CONT'.format(line)]
                ew = result['{}_FLUX'.format(line)] * factor   # rest frame [A]
                ewivar = result['{}_FLUX_IVAR'.format(line)] / factor**2

                # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar)
                ewlimit = fluxlimit * factor
            else:
                ew, ewivar, fluxlimit, ewlimit = 0.0, 0.0, 0.0, 0.0

            result['{}_EW'.format(line)] = ew
            result['{}_EW_IVAR'.format(line)] = ewivar
            result['{}_FLUX_LIMIT'.format(line)] = fluxlimit
            result['{}_EW_LIMIT'.format(line)] = ewlimit

            # simple QA
            if 'alpha' in line and False:
                import matplotlib.pyplot as plt
                _indx = np.where((emlinewave > (zwave - 15*sigma_cont * zwave / C_LIGHT)) *
                                (emlinewave < (zwave + 15*sigma_cont * zwave / C_LIGHT)))[0]
                plt.plot(emlinewave[_indx], emlineflux[_indx])
                plt.plot(emlinewave[_indx], specflux_nolines[_indx])
                plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
                plt.axhline(y=cmed, color='k')
                plt.axhline(y=cmed+csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.axhline(y=cmed-csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.savefig('junk.png')

        return result
    
    def emlineplot(self, data, result, continuum, qadir='.'):
        """Plot the emission-line spectrum and best-fitting model.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        redshift = result['Z']
        _emlinemodel = self.emlinemodel_bestfit(data['wave'], data['res'], result)

        sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]

        leg = {
            'targetid': '{} {}'.format(result['TARGETID'], -999),
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            'linevshift_forbidden': '$\Delta\,v_{{\\rm forbidden}}$={:.1f} km/s'.format(result['LINEVSHIFT_FORBIDDEN']),
            'linevshift_balmer': '$\Delta\,v_{{\\rm Balmer}}$={:.1f} km/s'.format(result['LINEVSHIFT_BALMER']),
            'linesigma_forbidden': '$\sigma_{{\\rm forbidden}}$={:.1f} km/s'.format(result['LINESIGMA_FORBIDDEN']),
            'linesigma_balmer': '$\sigma_{{\\rm Balmer}}$={:.1f} km/s'.format(result['LINESIGMA_BALMER']),
            }

        #fig, ax = plt.subplots(1, 4, figsize=(16, 10))#, sharey=True)
        fig = plt.figure(figsize=(16, 16))
        gs = fig.add_gridspec(3, 4, height_ratios=[4, 2, 2])
        #gs = fig.add_gridspec(2, 4, gridspec_kw={'width_ratios': 1.0, 'height_ratios': 0.5})

        # full spectrum
        bigax = fig.add_subplot(gs[0, :])

        ymin, ymax = 1e6, -1e6
        for ii in [0, 1, 2]: # iterate over cameras
            emlinewave = data['wave'][ii]
            emlineflux = data['flux'][ii] - continuum[ii]
            emlinemodel = _emlinemodel[ii]

            emlinesigma, good = _ivar2var(data['ivar'][ii], sigma=True)
            emlinewave = emlinewave[good]
            emlineflux = emlineflux[good]
            emlinesigma = emlinesigma[good]
            emlinemodel = emlinemodel[good]

            bigax.fill_between(emlinewave, emlineflux-emlinesigma,
                               emlineflux+emlinesigma, color=col1[ii], alpha=0.7)
            bigax.plot(emlinewave, emlinemodel, color=col2[ii], lw=2)

            # get the robust range
            sigflux, filtflux = np.std(emlineflux), median_filter(emlineflux, 5)

            if -3 * sigflux < ymin:
                ymin = -2 * sigflux
            #if np.min(filtflux) < ymin:
            #    ymin = np.min(filtflux)
            if np.min(emlinemodel) < ymin:
                ymin = 0.8 * np.min(emlinemodel)
                
            if 5 * sigflux > ymax:
                ymax = 4 * sigflux
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux)
            if np.max(emlinemodel) > ymax:
                ymax = np.max(emlinemodel) * 1.2

        bigax.text(0.95, 0.92, '{}'.format(leg['targetid']), 
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
        bigax.text(0.95, 0.86, r'{}'.format(leg['zredrock']),
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
        bigax.text(0.95, 0.80, r'{} {}'.format(leg['linevshift_balmer'], leg['linevshift_forbidden']),
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
        bigax.text(0.95, 0.74, r'{} {}'.format(leg['linesigma_balmer'], leg['linesigma_forbidden']),
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
                
        bigax.set_xlim(3500, 10000)
        bigax.set_ylim(ymin, ymax)
        
        #bigax.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 
        
        # zoom in on individual emission lines - use linetable!
        sig = 500.0 # [km/s]

        meanwaves = [np.mean([3730,3727]), 3971, 4103, 4342, 4863, np.mean([4960, 5008]), 6565, np.mean([6718.294, 6732.673])]
        deltawaves = [0, 0, 0, 0, 0, (5007-4959)/2, (6585-6550)/2, (6733-6718)/2]
        linenames = [r'[OII]', r'H$\epsilon$', r'H$\delta$', r'H$\gamma$', r'H$\beta$', r'[OIII]', r'H$\alpha$+[NII]', r'[SII]']
        nline = len(meanwaves)

        removelabels = np.ones(nline, bool)
        ymin, ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
        
        ax = []
        for iax, (meanwave, deltawave, linename) in enumerate(zip(meanwaves, deltawaves, linenames)):
            if iax < 4:
                xx = fig.add_subplot(gs[1, iax])
            else:
                xx = fig.add_subplot(gs[2, iax-4])
            ax.append(xx)

            wmin = (meanwave-deltawave)*(1+redshift)-2.5*sig*meanwave/3e5
            wmax = (meanwave+deltawave)*(1+redshift)+2.5*sig*meanwave/3e5

            # iterate over cameras
            for ii in [0, 1, 2]: # iterate over cameras
                emlinewave = data['wave'][ii]
                emlineflux = data['flux'][ii] - continuum[ii]
                emlinemodel = _emlinemodel[ii]

                emlinesigma, good = _ivar2var(data['ivar'][ii], sigma=True)
                emlinewave = emlinewave[good]
                emlineflux = emlineflux[good]
                emlinesigma = emlinesigma[good]
                emlinemodel = emlinemodel[good]

                indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
                #log.info(ii, linename, len(indx))
                if len(indx) > 1:
                    removelabels[iax] = False
                    xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                                    emlineflux[indx]+emlinesigma[indx], color=col1[ii], alpha=0.5)
                    xx.plot(emlinewave[indx], emlinemodel[indx], color=col2[ii], lw=3)

                    # get the robust range
                    sigflux, filtflux = np.std(emlineflux[indx]), median_filter(emlineflux[indx], 3)

                    _ymin, _ymax = -1.5 * sigflux, 3 * sigflux
                    if np.max(emlinemodel[indx]) > _ymax:
                        _ymax = np.max(emlinemodel[indx]) * 1.2
                    if np.max(filtflux) > _ymax:
                        _ymax = np.max(filtflux)
                    if np.min(emlinemodel[indx]) < _ymin:
                        _ymin = 0.8 * np.min(emlinemodel[indx])
                        
                    if _ymax > ymax[iax]:
                        ymax[iax] = _ymax
                    if _ymin < ymin[iax]:
                        ymin[iax] = _ymin
                    
                xx.text(0.08, 0.9, linename, ha='left', va='center',
                        transform=xx.transAxes, fontsize=20)
                    
        for iax, xx in enumerate(ax):
            if removelabels[iax]:
                xx.set_ylim(0, 1)
                xx.set_xticklabels([])
                xx.set_yticklabels([])
            else:
                xx.set_ylim(ymin[iax], ymax[iax])
                xlim = xx.get_xlim()
                #log.info(linenames[iax], xlim, np.diff(xlim))
                xx.xaxis.set_major_locator(ticker.MaxNLocator(3))
                #xx.xaxis.set_major_locator(ticker.MultipleLocator(20)) # wavelength spacing of ticks [Angstrom]
                #if iax == 2:
                #    pdb.set_trace()

        # common axis labels
        tp, bt, lf, rt = 0.95, 0.09, 0.12, 0.95
        
        fig.text(lf-0.07, (tp-bt)/2+bt,
                 r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
                 ha='center', va='center', rotation='vertical')
        fig.text((rt-lf)/2+lf, bt-0.06, r'Observed-frame Wavelength ($\AA$)',
                 ha='center', va='center')
            
        plt.subplots_adjust(wspace=0.27, top=tp, bottom=bt, left=lf, right=rt, hspace=0.22)

        pngfile = os.path.join(qadir, 'emlinefit-{}-{}-{}.png'.format(
            result['TILE'], result['NIGHT'], result['TARGETID']))
        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
