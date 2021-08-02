"""
fastspecfit.continuum
=====================

Methods and tools for continuum-fitting.

"""
import pdb # for debugging

import os, time
import numpy as np

import astropy.units as u

from fastspecfit.util import C_LIGHT
from desiutil.log import get_logger
log = get_logger()

def _fnnls_continuum(myargs):
    """Multiprocessing wrapper."""
    return fnnls_continuum(*myargs)

def fnnls_continuum(ZZ, xx, flux=None, ivar=None, modelflux=None,
                    support=None, get_chi2=False, jvendrow=False):
    """Fit a continuum using fNNLS. This function is a simple wrapper on fnnls; see
    the ContinuumFit.fnnls_continuum method for documentation.

        Mapping between mikeiovine fnnls(AtA, Aty) and jvendrow fnnls(Z, x) inputs:
          Z [mxn] --> A [mxn]
          x [mx1] --> y [mx1]

        And mikeiovine wants:
          A^T * A
          A^T * y

          AtA = A.T.dot(A)
          Aty = A.T.dot(y)

    """
    if jvendrow:
        from fnnls import fnnls

        if support is None:
            support = np.zeros(0, dtype=int)
            
        try:
            warn, coeff, _ = fnnls(ZZ, xx)#, P_initial=support)
        except:
            log.warning('fnnls failed to converge.')
            warn, coeff = True, np.zeros(modelflux.shape[1])
    else:
        #from fastnnls import fnnls
        from fastspecfit.fnnls import fnnls

        AtA = ZZ.T.dot(ZZ)
        Aty = ZZ.T.dot(xx)
        coeff = fnnls(AtA, Aty)
        warn = False
        
    #if warn:
    #    print('WARNING: fnnls did not converge after 5 iterations.')

    if get_chi2:
        chi2 = np.sum(ivar * (flux - modelflux.dot(coeff))**2)
        chi2 /= np.sum(ivar > 0) # reduced chi2
        return warn, coeff, chi2
    else:
        return warn, coeff
    
class ContinuumTools(object):
    def __init__(self, metallicity='Z0.0190', minwave=None, maxwave=6e4, seed=1):
        """Tools for dealing with stellar continua..

        """
        import fitsio
        from astropy.cosmology import FlatLambdaCDM
        from astropy.table import Table, Column

        from speclite import filters
        from desiutil.dust import SFDMap
        from fastspecfit.emlines import read_emlines
        from fastspecfit.io import FASTSPECFIT_TEMPLATES_NERSC, DUST_DIR_NERSC

        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        # pre-compute the luminosity distance on a grid
        #self.redshift_ref = np.arange(0.0, 5.0, 0.05)
        #self.dlum_ref = self.cosmo.luminosity_distance(self.redshift_ref).to(u.pc).value

        self.fluxnorm = 1e17 # normalization factor for the spectra
        self.massnorm = 1e10 # stellar mass normalization factor for the SSPs [Msun]

        self.metallicity = metallicity
        self.Z = float(metallicity[1:])
        self.library = 'CKC14z'
        self.isochrone = 'Padova' # would be nice to get MIST in here
        self.imf = 'Kroupa'

        # dust maps
        mapdir = os.path.join(os.environ.get('DUST_DIR', DUST_DIR_NERSC), 'maps')
        self.SFDMap = SFDMap(scaling=1.0, mapdir=mapdir)
        #self.SFDMap = SFDMap(scaling=0.86, mapdir=mapdir) # SF11 recalibration of the SFD maps
        self.RV = 3.1
        self.dustslope = 0.7

        # SSPs
        templates_dir = os.environ.get('FASTSPECFIT_TEMPLATES', FASTSPECFIT_TEMPLATES_NERSC)
        self.sspfile = os.path.join(templates_dir, 'SSP_{}_{}_{}_{}.fits'.format(
            self.isochrone, self.library, self.imf, self.metallicity))
        if not os.path.isfile(self.sspfile):
            log.warning('SSP templates file not found {}'.format(self.sspfile))
            raise IOError

        log.info('Reading {}'.format(self.sspfile))
        wave, wavehdr = fitsio.read(self.sspfile, ext='WAVE', header=True)
        flux = fitsio.read(self.sspfile, ext='FLUX')
        sspinfo = Table(fitsio.read(self.sspfile, ext='METADATA'))
        
        # Trim the wavelengths and select the number/ages of the templates.
        # https://www.sdss.org/dr14/spectro/galaxy_mpajhu
        if minwave is None:
            minwave = np.min(wave)
        keep = np.where((wave >= minwave) * (wave <= maxwave))[0]
        sspwave = wave[keep]

        myages = np.array([0.005, 0.025, 0.1, 0.2, 0.6, 0.9, 1.4, 2.5, 5, 10.0, 13.0])*1e9
        iage = np.array([np.argmin(np.abs(sspinfo['age']-myage)) for myage in myages])
        sspflux = flux[:, iage][keep, :] # flux[keep, ::5]
        sspinfo = sspinfo[iage]

        nage = len(sspinfo)
        npix = len(sspwave)

        self.pixkms = wavehdr['PIXSZBLU'] # pixel size [km/s]

        # add AGN templates here?
        if False:
            # https://www.aanda.org/articles/aa/pdf/2017/08/aa30378-16.pdf
            # F_nu \propto \nu^(-alpha) or F_lambda \propto \lambda^(alpha-2)
            self.agn_lambda0 = 4020.0 # [Angstrom]a
            self.agn_alpha = [0.5, 1.0, 1.5, 2.0]
            nagn = len(self.agn_alpha)
            
            agnflux = np.zeros((npix, nagn), 'f4')
            #import matplotlib.pyplot as plt
            for ii, alpha in enumerate(self.agn_alpha):
                agnflux[:, ii] = sspwave**(alpha-2) / self.agn_lambda0**(alpha-2)
                #plt.plot(sspwave, agnflux[:, ii])
            #plt.xlim(3000, 9000) ; plt.ylim(0.1, 2.2) ; plt.savefig('junk.png')

            sspflux = np.vstack((agnflux.T, sspflux.T)).T

            sspinfo = Table(np.hstack([sspinfo[:nagn], sspinfo]))
            sspinfo.add_column(Column(name='agn_alpha', length=nagn+nage, dtype='f4'))
            sspinfo['age'][:nagn] = 0.0
            sspinfo['mstar'][:nagn] = 0.0
            sspinfo['lbol'][:nagn] = 0.0
            sspinfo['agn_alpha'][:nagn] = self.agn_alpha

            nage = len(sspinfo)

        self.sspwave = sspwave
        self.sspflux = sspflux                     # no dust, no velocity broadening [npix,nage]
        self.sspinfo = sspinfo
        self.nage = nage
        self.npix = npix

        # emission lines
        self.linetable = read_emlines()

        # photometry
        self.bands = ['g', 'r', 'z', 'W1', 'W2']
        self.synth_bands = ['g', 'r', 'z'] # for synthesized photometry
        self.fiber_bands = ['g', 'r', 'z'] # for fiber fluxes

        self.decam = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z')
        self.bassmzls = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z')

        self.decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z',
                                              'wise2010-W1', 'wise2010-W2')
        self.bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z',
                                                 'wise2010-W1', 'wise2010-W2')

        # rest-frame filters
        self.absmag_bands = ['U', 'B', 'V', 'sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z', 'W1']
        self.absmag_filters = filters.load_filters('bessell-U', 'bessell-B', 'bessell-V',
                                                   'sdss2010-u', 'sdss2010-g', 'sdss2010-r',
                                                   'sdss2010-i', 'sdss2010-z', 'wise2010-W1')
        self.absmag_bandshift = np.array([0.0, 0.0, 0.0,
                                          0.1, 0.1, 0.1,
                                          0.1, 0.1, 0.0])

        # used in one place...
        self.rand = np.random.RandomState(seed=seed)

    @staticmethod
    def get_dn4000(wave, flam, flam_ivar=None, redshift=None, rest=True):
        """Compute DN(4000) and, optionally, the inverse variance.

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
        from fastspecfit.util import ivar2var

        dn4000, dn4000_ivar = 0.0, 0.0

        if rest:
            flam2fnu =  wave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
        else:
            wave = np.copy(wave)
            wave /= (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * wave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]

        if flam_ivar is None:
            goodmask = np.ones(len(flam), bool) # True is good
        else:
            goodmask = flam_ivar > 0

        indxblu = np.where((wave >= 3850.) * (wave <= 3950.) * goodmask)[0]
        indxred = np.where((wave >= 4000.) * (wave <= 4100.) * goodmask)[0]
        if len(indxblu) < 5 or len(indxred) < 5:
            return dn4000, dn4000_ivar

        blufactor, redfactor = 3950.0 - 3850.0, 4100.0 - 4000.0
        deltawave = np.gradient(wave) # should be constant...

        fnu = flam * flam2fnu # [erg/s/cm2/Hz]

        numer = blufactor * np.sum(deltawave[indxred] * fnu[indxred])
        denom = redfactor * np.sum(deltawave[indxblu] * fnu[indxblu])
        if denom == 0.0:
            log.warning('DN(4000) is ill-defined!')
            return dn4000, dn4000_ivar
        dn4000 =  numer / denom

        if flam_ivar is not None:
            fnu_ivar = flam_ivar / flam2fnu**2
            fnu_var, _ = ivar2var(fnu_ivar)

            numer_var = blufactor**2 * np.sum(deltawave[indxred] * fnu_var[indxred])
            denom_var = redfactor**2 * np.sum(deltawave[indxblu] * fnu_var[indxblu])
            dn4000_var = (numer_var + numer**2 * denom_var) / denom**2
            if dn4000_var <= 0:
                log.warning('DN(4000) variance is ill-defined!')
                dn4000_ivar = 0.0
            else:
                dn4000_ivar = 1.0 / dn4000_var

        return dn4000, dn4000_ivar

    @staticmethod
    def parse_photometry(bands, maggies, lambda_eff, ivarmaggies=None,
                         nanomaggies=True, nsigma=1.0):
        """Parse input (nano)maggies to various outputs and pack into a table.

        Parameters
        ----------
        flam - 10-17 erg/s/cm2/A
        fnu - 10-17 erg/s/cm2/Hz
        abmag - AB mag
        nanomaggies - input maggies are actually 1e-9 maggies

        nsigma - magnitude limit 

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
        phot.add_column(Column(name='band', data=bands))
        phot.add_column(Column(name='lambda_eff', length=nband, dtype='f4'))
        phot.add_column(Column(name='nanomaggies', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='nanomaggies_ivar', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='flam', length=nband, shape=(ngal, ), dtype='f8')) # note f8!
        phot.add_column(Column(name='flam_ivar', length=nband, shape=(ngal, ), dtype='f8'))
        phot.add_column(Column(name='abmag', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag_ivar', length=nband, shape=(ngal, ), dtype='f4'))
        #phot.add_column(Column(name='abmag_err', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag_brighterr', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag_fainterr', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag_limit', length=nband, shape=(ngal, ), dtype='f4'))

        if ivarmaggies is None:
            ivarmaggies = np.zeros_like(maggies)

        ## Gaia-only targets all have grz=-99 fluxes (we now cut these out in
        ## io.DESISpectra.find_specfiles)
        #if np.all(maggies==-99):
        #    log.warning('Gaia-only targets not supported.')
        #    raise ValueError

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

        # deal with measurements
        good = np.where(maggies > 0)[0]        
        if len(good) > 0:
            if maggies.ndim > 1:
                igood, jgood = np.unravel_index(good, maggies.shape)
                goodmaggies = maggies[igood, jgood]                
            else:
                igood, jgood = good, [0]
                goodmaggies = maggies[igood]
            phot['abmag'][igood, jgood] = (-2.5 * np.log10(nanofactor * goodmaggies)).astype('f4')
        
        # deal with the uncertainties
        snr = maggies * np.sqrt(ivarmaggies)
        good = np.where(snr > nsigma)[0]
        upper = np.where((ivarmaggies > 0) * (snr <= nsigma))[0]
        if maggies.ndim > 1:
            if len(upper) > 0:
                iupper, jupper = np.unravel_index(upper, maggies.shape)
                abmag_limit = +2.5 * np.log10(np.sqrt(ivarmaggies[iupper, jupper]) / nsigma) # note "+" instead of 1/ivarmaggies
                
            igood, jgood = np.unravel_index(good, maggies.shape)
            maggies = maggies[igood, jgood]
            ivarmaggies = ivarmaggies[igood, jgood]
            errmaggies = 1 / np.sqrt(ivarmaggies)
            #fracerr = 1 / snr[igood, jgood]
        else:
            if len(upper) > 0:
                iupper, jupper = upper, [0]
                abmag_limit = +2.5 * np.log10(np.sqrt(ivarmaggies[iupper]) / nsigma)
                
            igood, jgood = good, [0]
            maggies = maggies[igood]
            ivarmaggies = ivarmaggies[igood]
            errmaggies = 1 / np.sqrt(ivarmaggies)
            #fracerr = 1 / snr[igood]

        # significant detections
        if len(good) > 0:
            phot['abmag_brighterr'][igood, jgood] = errmaggies / (0.4 * np.log(10) * (maggies+errmaggies)).astype('f4') # bright end (flux upper limit)
            phot['abmag_fainterr'][igood, jgood] = errmaggies / (0.4 * np.log(10) * (maggies-errmaggies)).astype('f4') # faint end (flux lower limit)
            #phot['abmag_loerr'][igood, jgood] = +2.5 * np.log10(1 + fracerr) # bright end (flux upper limit)
            #phot['abmag_uperr'][igood, jgood] = +2.5 * np.log10(1 - fracerr) # faint end (flux lower limit)
            #test = 2.5 * np.log(np.exp(1)) * fracerr # symmetric in magnitude (approx)

            # approximate the uncertainty as being symmetric in magnitude
            phot['abmag_ivar'][igood, jgood] = (ivarmaggies * (maggies * 0.4 * np.log(10))**2).astype('f4')
            
        if len(upper) > 0:
            phot['abmag_limit'][iupper, jupper] = abmag_limit.astype('f4')
            
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
    
    def dust_attenuation(self, wave, AV):
        """Compute the dust attenuation curve A(lambda)/A(V) from Charlot & Fall 2000.

        ToDo: add a UV bump and IGM attenuation!
          https://gitlab.lam.fr/cigale/cigale/-/blob/master/pcigale/sed_modules/dustatt_powerlaw.py#L42

        """
        return 10**(-0.4 * AV * (wave / 5500.0)**(-self.dustslope))

    def smooth_and_resample(self, sspflux, sspwave, specwave=None, specres=None):
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
        from redrock.rebin import trapz_rebin

        if specwave is None:
            resampflux = sspflux 
        else:
            trim = (sspwave > (specwave.min()-10.0)) * (sspwave < (specwave.max()+10.0))
            resampflux = trapz_rebin(sspwave[trim], sspflux[trim], specwave)
            #try:
            #    resampflux = trapz_rebin(sspwave[trim], sspflux[trim], specwave)
            #except:
            #    pdb.set_trace()

        if specres is None:
            smoothflux = resampflux
        else:
            smoothflux = specres.dot(resampflux)

        return smoothflux # [noutpix]
    
    def SSP2data(self, _sspflux, _sspwave, redshift=0.0, AV=None, vdisp=None,
                 cameras=['b', 'r', 'z'], specwave=None, specres=None, coeff=None,
                 south=True, synthphot=True):
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
            dfactor = (10.0 / self.cosmo.luminosity_distance(redshift).to(u.pc).value)**2
            #dfactor = (10.0 / np.interp(redshift, self.redshift_ref, self.dlum_ref))**2
            factor = (self.fluxnorm * self.massnorm * dfactor / (1.0 + redshift))[np.newaxis, np.newaxis]
            zsspflux = sspflux * factor
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
                sspphot = self.parse_photometry(self.bands, maggies, effwave, nanomaggies=False)
                #log.info('Synthesizing photometry took: {:.2f} sec'.format(time.time()-t0))
            
        # Are we returning per-camera spectra or a single model? Handle that here.
        #t0 = time.time()
        if specwave is None and specres is None:
            datasspflux = []
            for imodel in np.arange(nmodel):
                datasspflux.append(self.smooth_and_resample(zsspflux[:, imodel], zsspwave))
            datasspflux = np.vstack(datasspflux).T
                
            if vdisp:
                 datasspflux = self.convolve_vdisp(datasspflux, vdisp)
                 
            # optionally compute the best-fitting model
            if coeff is not None:
                datasspflux = datasspflux.dot(coeff)
                if synthphot:
                    maggies = filters.get_ab_maggies(datasspflux, zsspwave, axis=0)
                    maggies = np.array(maggies.as_array().tolist()[0])
                    maggies /= self.fluxnorm * self.massnorm
                    sspphot = self.parse_photometry(self.bands, maggies, effwave, nanomaggies=False)
        else:
            # loop over cameras and then multiprocess over age
            datasspflux = []
            for icamera in np.arange(len(cameras)): # iterate on cameras
                _datasspflux = []
                for imodel in np.arange(nmodel):
                    _datasspflux.append(self.smooth_and_resample(
                        zsspflux[:, imodel], zsspwave, specwave=specwave[icamera],
                        specres=specres[icamera]))
                _datasspflux = np.vstack(_datasspflux).T
                if vdisp:
                    _datasspflux = self.convolve_vdisp(_datasspflux, vdisp)
                if coeff is not None:
                    _datasspflux = _datasspflux.dot(coeff)
                datasspflux.append(_datasspflux)
                
        #log.info('Resampling took: {:.2f} sec'.format(time.time()-t0))

        return datasspflux, sspphot # vector or 3-element list of [npix,nmodel] spectra

    def smooth_residuals(self, continuummodel, specwave, specflux, specivar,
                         linemask, percamera=False, binwave=0.8*160):
        """Derive a median-smoothed correction to the continuum fit in order to pick up
        any unmodeled flux, before emission-line fitting.

        Parameters
        ----------

        Returns
        -------

        Notes
        -----
        There are a few different algorithms in here, but the default one is to
        do a very simple median-smoothing on the camera-stacked
        (percamera=False) spectrum, which prevents any crazy discontinuities
        between the cameras.
        
        https://github.com/moustakas/moustakas-projects/blob/master/ages/ppxf/ages_gandalf_specfit.pro#L138-L145

        """
        from scipy.ndimage import median_filter

        def robust_median(wave, flux, ivar, binwave):
            from scipy import ptp
            from scipy.stats import sigmaclip
            #from scipy.signal import savgol_filter
            
            smoothflux = np.zeros_like(flux)
            minwave, maxwave = np.min(wave), np.max(wave)
            nbin = int(ptp(wave) / binwave)
            wavebins = np.linspace(minwave, maxwave, nbin)
            idx  = np.digitize(wave, wavebins)
            for kk in np.arange(nbin):
                I = (idx == kk)
                J = I * (ivar > 0)
                #print(kk, np.sum(I), np.sum(J))
                if np.sum(J) > 10:
                    clipflux, _, _ = sigmaclip(flux[J], low=5.0, high=5.0)
                    smoothflux[I] = np.median(clipflux)
                    #print(kk, np.sum(J), np.sum(I), np.median(wave[I]), smoothflux[I][0])
            return smoothflux
                
        medbin = np.int(binwave * 2)
        #print(binwave, medbin)

        if percamera:
            smooth_continuum = []
            for icam in np.arange(len(specflux)): # iterate over cameras        
                residuals = specflux[icam] - continuummodel[icam]
                if False:
                    smooth1 = robust_median(specwave[icam], residuals, specivar[icam], binwave)
                    #smooth1 = robust_median(specwave[icam], residuals, specivar[icam] * linemask[icam], binwave)
                    smooth2 = median_filter(smooth1, medbin, mode='nearest')
                    smooth_continuum.append(smooth2)
                else:
                    # Fragile...replace the line-affected pixels with noise to make the
                    # smoothing better-behaved.
                    pix_nolines = linemask[icam]      # unaffected by a line = True
                    pix_emlines = np.logical_not(pix_nolines) # affected by line = True
                    residuals[pix_emlines] = (self.rand.normal(np.count_nonzero(pix_emlines)) *
                                              np.std(residuals[pix_nolines]) +
                                              np.median(residuals[pix_nolines]))
                    smooth1 = median_filter(residuals, medbin, mode='nearest')
                    #smooth2 = savgol_filter(smooth1, 151, 2)
                    smooth2 = median_filter(smooth1, medbin//2, mode='nearest')
                    smooth_continuum.append(smooth2)
        else:
            residuals = specflux - continuummodel
            pix_nolines = linemask      # unaffected by a line = True
            pix_emlines = np.logical_not(pix_nolines) # affected by line = True
            
            residuals[pix_emlines] = (self.rand.normal(np.count_nonzero(pix_emlines)) *
                                      np.std(residuals[pix_nolines]) + np.median(residuals[pix_nolines]))
            smooth1 = median_filter(residuals, medbin, mode='nearest')
            smooth_continuum = median_filter(smooth1, medbin//2, mode='nearest')

        return smooth_continuum

class ContinuumFit(ContinuumTools):
    def __init__(self, metallicity='Z0.0190', minwave=None, maxwave=6e4):
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

        Notes
        -----
        Need to document all the attributes.
        
        Plans for improvement (largely in self.fnnls_continuum).
          - Update the continuum redshift using cross-correlation.
          - Don't draw reddening from a flat distribution (try gamma or a custom
            distribution of the form x**2*np.exp(-2*x/scale).

        """
        super(ContinuumFit, self).__init__(metallicity=metallicity, minwave=minwave, maxwave=maxwave)
        
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
        AVmin, AVmax, dAV, AV_nominal = (0.0, 1.5, 0.1, 0.0)
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

        # Next, precompute a grid of spectra convolved to the nominal velocity
        # dispersion with reddening applied. This isn't quite right redward of
        # ~1 micron where the pixel size changes, but fix that later.
        sspflux_dustvdisp = []
        for AV in self.AV:
            atten = self.dust_attenuation(self.sspwave, AV)
            _sspflux_dustvdisp = self.convolve_vdisp(self.sspflux * atten[:, np.newaxis], self.vdisp_nominal)
            sspflux_dustvdisp.append(_sspflux_dustvdisp)

        # nominal velocity broadening on a grid of A(V) [npix,nage,nAV]
        self.sspflux_dustvdisp = np.stack(sspflux_dustvdisp, axis=-1) # [npix,nage,nAV]

        # table of emission lines to fit
        self.linemask_sigma = 200.0 # [km/s]
        self.linemask_sigma_broad = 1000.0 # [km/s]

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
        out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='CONTINUUM_VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        out.add_column(Column(name='CONTINUUM_VDISP_IVAR', length=nobj, dtype='f4', unit=u.second**2/u.kilometer**2))
        for cam in ['B', 'R', 'Z']:
            out.add_column(Column(name='CONTINUUM_SNR_{}'.format(cam), length=nobj, dtype='f4')) # median S/N in each camera
        #out.add_column(Column(name='CONTINUUM_SNR', length=nobj, shape=(3,), dtype='f4')) # median S/N in each camera

        # maximum correction to the median-smoothed continuum
        for cam in ['B', 'R', 'Z']:
            out.add_column(Column(name='CONTINUUM_SMOOTHCORR_{}'.format(cam), length=nobj, dtype='f4')) 
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
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))

        return out

    def init_phot_output(self, nobj=1):
        """Initialize the photometric output data table.

        """
        from astropy.table import Table, Column
        
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        #out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))

        if False:
            for band in self.fiber_bands:
                out.add_column(Column(name='FIBERTOTFLUX_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.nanomaggy)) # observed-frame fiber photometry
                #out.add_column(Column(name='FIBERTOTFLUX_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/u.nanomaggy**2))
            for band in self.bands:
                out.add_column(Column(name='FLUX_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.nanomaggy)) # observed-frame photometry
                out.add_column(Column(name='FLUX_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/u.nanomaggy**2))
                
        for band in self.absmag_bands:
            out.add_column(Column(name='KCORR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.mag))
            out.add_column(Column(name='ABSMAG_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.mag)) # absolute magnitudes
            out.add_column(Column(name='ABSMAG_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/u.mag**2))

        return out

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

    def kcorr_and_absmag(self, data, continuum, coeff, band_shift=0.0, snrmin=2.0):
        """Computer K-corrections and absolute magnitudes.

        # To get the absolute r-band magnitude we would do:
        #   M_r = m_X_obs + 2.5*log10(r_synth_rest/X_synth_obs)
        # where X is the redshifted bandpass

        """
        redshift = data['zredrock']
        
        if data['photsys'] == 'S':
            filters_in = self.decamwise
        else:
            filters_in = self.bassmzlswise
        filters_out = self.absmag_filters
        nout = len(filters_out)

        lambda_in = filters_in.effective_wavelengths.value
        lambda_out = filters_out.effective_wavelengths.value / (1 + band_shift)

        # redshifted wavelength array and distance modulus
        zsspwave = self.sspwave * (1 + redshift)
        dmod = self.cosmo.distmod(redshift).value

        maggies = data['phot']['nanomaggies'].data * 1e-9
        ivarmaggies = data['phot']['nanomaggies_ivar'].data / 1e-9**2

        # input bandpasses, observed frame; maggies and bestmaggies should be
        # very close.
        bestmaggies = filters_in.get_ab_maggies(continuum / self.fluxnorm, zsspwave)
        bestmaggies = np.array(bestmaggies.as_array().tolist()[0])

        # output bandpasses, rest frame -- need to shift the filter curves
        # blueward by a factor of 1+band_shift!
        synth_outmaggies_rest = filters_out.get_ab_maggies(continuum * (1 + redshift) / self.fluxnorm, self.sspwave) 
        synth_outmaggies_rest = np.array(synth_outmaggies_rest.as_array().tolist()[0])

        # output bandpasses, observed frame
        synth_outmaggies_obs = filters_out.get_ab_maggies(continuum / self.fluxnorm, zsspwave)
        synth_outmaggies_obs = np.array(synth_outmaggies_obs.as_array().tolist()[0])

        absmag = np.zeros(nout, dtype='f4')
        ivarabsmag = np.zeros(nout, dtype='f4')
        kcorr = np.zeros(nout, dtype='f4')
        for jj in np.arange(nout):
            lambdadist = np.abs(lambda_in / (1 + redshift) - lambda_out[jj])
            # K-correct from the nearest "good" bandpass (to minimizes the K-correction)
            #oband = np.argmin(lambdadist)
            #oband = np.argmin(lambdadist + (ivarmaggies == 0)*1e10)
            oband = np.argmin(lambdadist + (maggies*np.sqrt(ivarmaggies) < snrmin)*1e10)
            kcorr[jj] = + 2.5 * np.log10(synth_outmaggies_rest[jj] / bestmaggies[oband])

            # m_R = M_Q + DM(z) + K_QR(z) or
            # M_Q = m_R - DM(z) - K_QR(z)
            if maggies[oband] * np.sqrt(ivarmaggies[oband]) > snrmin:
            #if (maggies[oband] > 0) and (ivarmaggies[oband]) > 0:
                absmag[jj] = -2.5 * np.log10(maggies[oband]) - dmod - kcorr[jj]
                ivarabsmag[jj] = maggies[oband]**2 * ivarmaggies[oband] * (0.4 * np.log(10.))**2
            else:
                # if we use synthesized photometry then ivarabsmag is zero
                # (which should never happen?)
                absmag[jj] = -2.5 * np.log10(synth_outmaggies_rest[jj]) - dmod

        # get the stellar mass
        if False:
            nage = len(coeff)
            dfactor = (10.0 / self.cosmo.luminosity_distance(redshift).to(u.pc).value)**2        
            mstar = self.sspinfo['mstar'][:nage].dot(coeff) * self.massnorm * self.fluxnorm * dfactor / (1+redshift)
        
        # From Taylor+11, eq 8
        #https://researchportal.port.ac.uk/ws/files/328938/MNRAS_2011_Taylor_1587_620.pdf
        #mstar = 1.15 + 0.7*(absmag[1]-absmag[3]) - 0.4*absmag[3]

        return kcorr, absmag, ivarabsmag

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

    def continuum_fastphot(self, data):
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

        # Prepare the reddened and unreddened SSP templates. Note that we ignore
        # templates which are older than the age of the universe at the galaxy
        # redshift.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustvdisp, zsspphot_dustvdisp = self.SSP2data(
            self.sspflux_dustvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=None, specres=None,
            south=data['photsys'] == 'S')
        log.info('Preparing the models took {:.2f} sec'.format(time.time()-t0))
        
        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = data['phot']['flam_ivar'].data / self.fluxnorm**2
        zsspflam_dustvdisp = zsspphot_dustvdisp['flam'].data * self.fluxnorm * self.massnorm # [nband,nage*nAV]
        assert(np.all(objflamivar >= 0))

        inodust = np.asscalar(np.where(self.AV == 0)[0]) # should always be index 0

        npix, nmodel = zsspflux_dustvdisp.shape
        nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)

        zsspflam_dustvdisp = zsspflam_dustvdisp.reshape(len(self.bands), nage, self.nAV) # [nband,nage,nAV]

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
                                              south=data['photsys'] == 'S')
        coeff, chi2min = self._fnnls_parallel(bestphot['flam'].data*self.massnorm*self.fluxnorm,
                                              objflam, objflamivar) # bestphot['flam'] is [nband, nage]
        continuummodel = bestsspflux.dot(coeff)

        # Compute DN4000, K-corrections, and rest-frame quantities.
        if np.count_nonzero(coeff > 0) == 0:
            log.warning('Continuum coefficients are all zero!')
            chi2min, dn4000, meanage = 1e6, -1.0, -1.0
            kcorr = np.zeros(len(self.absmag_bands))
            absmag = np.zeros(len(self.absmag_bands))-99.0
            ivarabsmag = np.zeros(len(self.absmag_bands))
        else:
            dn4000, _ = self.get_dn4000(self.sspwave, continuummodel, rest=True)
            meanage = self.get_meanage(coeff)
            kcorr, absmag, ivarabsmag = self.kcorr_and_absmag(data, continuummodel, coeff)

            log.info('Photometric DN(4000)={:.3f}, Age={:.2f} Gyr, Mr={:.2f} mag'.format(
                dn4000, meanage, absmag[1]))

        # Pack it up and return.
        result['CONTINUUM_COEFF'][0][:nage] = coeff
        result['CONTINUUM_CHI2'][0] = chi2min
        result['CONTINUUM_AGE'][0] = meanage
        result['CONTINUUM_AV'][0] = AVbest
        result['CONTINUUM_AV_IVAR'][0] = AVivar
        result['DN4000_MODEL'][0] = dn4000
        if False:
            for iband, band in enumerate(self.fiber_bands):
                result['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
                #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
            for iband, band in enumerate(self.bands):
                result['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
                result['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        for iband, band in enumerate(self.absmag_bands):
            result['KCORR_{}'.format(band.upper())] = kcorr[iband]
            result['ABSMAG_{}'.format(band.upper())] = absmag[iband]
            result['ABSMAG_IVAR_{}'.format(band.upper())] = ivarabsmag[iband]

        return result, continuummodel
    
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
        ToDo:
          - Use cross-correlation to update the redrock redshift.
          - Need to mask more emission lines than we fit (e.g., Mg II).

        """
        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_spec_output()

        redshift = data['zredrock']
        result['CONTINUUM_Z'] = redshift
        for icam, cam in enumerate(data['cameras']):
            result['CONTINUUM_SNR_{}'.format(cam.upper())] = data['snr'][icam]

        # Prepare the reddened and unreddened SSP templates. Note that we ignore
        # templates which are older than the age of the universe at the galaxy
        # redshift.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustvdisp, _ = self.SSP2data(
            self.sspflux_dustvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=data['wave'], specres=data['res'],
            cameras=data['cameras'], synthphot=False)
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
                                                   cameras=data['cameras'], synthphot=False)
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
                                              cameras=data['cameras'], south=data['photsys'] == 'S')
        bestsspflux = np.concatenate(bestsspflux, axis=0)
        coeff, chi2min = self._fnnls_parallel(bestsspflux, specflux, specivar)

        # Get the mean age and DN(4000).
        bestfit = bestsspflux.dot(coeff)
        meanage = self.get_meanage(coeff)
        dn4000_model, _ = self.get_dn4000(specwave, bestfit, redshift=redshift, rest=False)
        dn4000, dn4000_ivar = self.get_dn4000(specwave, specflux, specivar, redshift=redshift, rest=False)
        log.info('Spectroscopic DN(4000)={:.3f}, Age={:.2f} Gyr'.format(dn4000, meanage))

        # Do a quick median-smoothing of the stellar continuum-subtracted
        # residuals, to help with the emission-line fitting.
        _smooth_continuum = self.smooth_residuals(
            bestfit, specwave, specflux, np.hstack(data['ivar']),
            np.hstack(data['linemask']))
        #smooth_continuum = self.smooth_residuals(
        #    bestfit, data['coadd_wave'], data['coadd_flux'],
        #    data['coadd_ivar'], data['coadd_linemask'],
        #    npixpercam)
        
        # Unpack the continuum into individual cameras.
        continuummodel = []
        smooth_continuum = []
        for icam in np.arange(len(data['cameras'])): # iterate over cameras
            ipix = np.sum(npixpercam[:icam+1])
            jpix = np.sum(npixpercam[:icam+2])
            continuummodel.append(bestfit[ipix:jpix])
            smooth_continuum.append(_smooth_continuum[ipix:jpix])

        ## Like above, but with per-camera smoothing.
        #smooth_continuum = self.smooth_residuals(
        #    continuummodel, data['wave'], data['flux'],
        #    data['ivar'], data['linemask'], percamera=False)

        if False:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1)
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                resid = data['flux'][icam]-continuummodel[icam]
                ax[0].plot(data['wave'][icam], resid)
                ax[1].plot(data['wave'][icam], resid-smooth_continuum[icam])
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                resid = data['flux'][icam]-continuummodel[icam]
                pix_emlines = np.logical_not(data['linemask'][icam]) # affected by line = True
                ax[0].scatter(data['wave'][icam][pix_emlines], resid[pix_emlines], s=30, color='red')
                ax[0].plot(data['wave'][icam], smooth_continuum[icam], color='k', alpha=0.7, lw=2)
            plt.savefig('junk.png')
            pdb.set_trace()

        # Pack it in and return.
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

        for icam, cam in enumerate(data['cameras']):
            nonzero = np.abs(continuummodel[icam]) > 1e-5
            if np.sum(nonzero) > 0:
                corr = np.mean(smooth_continuum[icam][nonzero] / continuummodel[icam][nonzero])
                result['CONTINUUM_SMOOTHCORR_{}'.format(cam.upper())] = corr
                
        log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
            100*result['CONTINUUM_SMOOTHCORR_B'][0], 100*result['CONTINUUM_SMOOTHCORR_R'][0],
            100*result['CONTINUUM_SMOOTHCORR_Z'][0]))
        
        return result, continuummodel, smooth_continuum
    
    def qa_fastphot(self, data, fastphot, metadata, coadd_type='healpix',
                    outdir=None, outprefix=None):
        """QA of the best-fitting continuum.

        """
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        from fastspecfit.util import ivar2var
    
        sns.set(context='talk', style='ticks', font_scale=1.2)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]
        ymin, ymax = 1e6, -1e6

        redshift = metadata['Z']

        if metadata['PHOTSYS'] == 'S':
            filters = self.decam
            allfilters = self.decamwise
        else:
            filters = self.bassmzls
            allfilters = self.bassmzlswise

        if outdir is None:
            outdir = '.'
        if outprefix is None:
            outprefix = 'fastphot'

        if coadd_type == 'healpix':
            title = 'Survey/Program/HealPix: {}/{}/{}, TargetID: {}'.format(
                    metadata['SURVEY'], metadata['FAPRGRM'], metadata['HPXPIXEL'], metadata['TARGETID'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['FAPRGRM'], metadata['HPXPIXEL'], metadata['TARGETID']))
        elif coadd_type == 'cumulative':
            title = 'Tile/thruNight: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['THRUNIGHT'], metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-thru{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['THRUNIGHT'], metadata['TARGETID']))
        elif coadd_type == 'pernight':
            title = 'Tile/Night: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID'],
                    metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID']))
        elif coadd_type == 'perexp':
            title = 'Tile/Night/Expid: {}/{}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['EXPID'],
                    metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-exp{:08d}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'],
                    metadata['EXPID'], metadata['TARGETID']))
        else:
            pass

        # rebuild the best-fitting photometric model fit
        continuum_phot, _ = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift,
                                          AV=fastphot['CONTINUUM_AV'],
                                          coeff=fastphot['CONTINUUM_COEFF'] * self.massnorm,
                                          synthphot=False)
        continuum_wave_phot = self.sspwave * (1 + redshift)

        wavemin, wavemax = 0.2, 6.0
        indx = np.where((continuum_wave_phot/1e4 > wavemin) * (continuum_wave_phot/1e4 < wavemax))[0]     

        phot = self.parse_photometry(self.bands,
                                     maggies=np.array([metadata['FLUX_{}'.format(band.upper())] for band in self.bands]),
                                     ivarmaggies=np.array([metadata['FLUX_IVAR_{}'.format(band.upper())] for band in self.bands]),
                                     lambda_eff=allfilters.effective_wavelengths.value)
        fiberphot = self.parse_photometry(self.fiber_bands,
                                          maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in self.fiber_bands]),
                                          lambda_eff=filters.effective_wavelengths.value)
        
        #fibertotphot = self.parse_photometry(self.fiber_bands,
        #                                     maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in self.fiber_bands]),
        #                                     lambda_eff=filters.effective_wavelengths.value)
        #if specfit:
        #    synthphot = self.parse_photometry(self.synth_bands,
        #                                      maggies=np.array([specfit['FLUX_SYNTH_{}'.format(band.upper())] for band in self.synth_bands]),
        #                                      lambda_eff=filters.effective_wavelengths.value)
        #    synthmodelphot = self.parse_photometry(self.synth_bands,
        #                                           maggies=np.array([specfit['FLUX_SYNTH_MODEL_{}'.format(band.upper())] for band in self.synth_bands]),
        #                                           lambda_eff=filters.effective_wavelengths.value)
        #else:
        #    synthphot, synthmodelphot = None, None
        
        fig, ax = plt.subplots(figsize=(12, 8))

        if np.any(continuum_phot <= 0):
            log.warning('Best-fitting photometric continuum is all zeros or negative!')
            continuum_phot_abmag = continuum_phot*0 + np.median(fiberphot['abmag'])
        else:
            factor = 10**(0.4 * 48.6) * continuum_wave_phot**2 / (C_LIGHT * 1e13) / self.fluxnorm / self.massnorm # [erg/s/cm2/A --> maggies]
            continuum_phot_abmag = -2.5*np.log10(continuum_phot * factor)
            ax.plot(continuum_wave_phot[indx] / 1e4, continuum_phot_abmag[indx], color='gray', zorder=1)

        # we have to set the limits *before* we call errorbar, below!
        dm = 0.75
        good = phot['abmag_ivar'] > 0
        ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(continuum_phot_abmag[indx]))) + dm
        ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(continuum_phot_abmag[indx]))) - dm
        if ymin > 31:
            ymin = 31
        if np.isnan(ymin) or np.isnan(ymax):
            raise('Problem here!')

        ax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        #ax.set_ylabel(r'AB mag') 
        ax.set_ylabel(r'Apparent Brightness (AB mag)') 
        ax.set_xlim(wavemin, wavemax)
        ax.set_ylim(ymin, ymax)

        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.set_xticks([0.3, 0.4, 0.6, 1.0, 1.5, 3.0, 5.0])

        ax.set_title(title, fontsize=20)

        # integrated flux / photometry
        #ax.scatter(phot['lambda_eff']/1e4, phot['abmag'],
        #           marker='s', s=130, facecolor='red', edgecolor='k',
        #           label=r'$grzW1W2$ (imaging)', alpha=1.0, zorder=3)
        abmag = np.squeeze(phot['abmag'])
        abmag_limit = np.squeeze(phot['abmag_limit'])
        abmag_fainterr = np.squeeze(phot['abmag_fainterr'])
        abmag_brighterr = np.squeeze(phot['abmag_brighterr'])
        yerr = np.squeeze([abmag_fainterr, abmag_brighterr])

        lolims = abmag_limit > 0
        #lolims[[2, 4]] = True
        if np.count_nonzero(lolims) > 0:
            abmag[lolims] = abmag_limit[lolims]

        ax.errorbar(phot['lambda_eff']/1e4, abmag, lolims=lolims,
                    yerr=yerr,
                    fmt='s', markersize=11, markeredgewidth=3, markeredgecolor='k',
                    markerfacecolor='red', elinewidth=3, ecolor='red', capsize=4,
                    label=r'$grz\,W_{1}W_{2}$', zorder=2)

        good = np.where(fiberphot['abmag'] > 0)[0]
        if len(good) > 0:
            ax.scatter(fiberphot['lambda_eff'][good]/1e4, fiberphot['abmag'][good],
                        marker='o', s=150, facecolor='blue', edgecolor='k',
                        label=r'$grz$ (fiberflux)', alpha=0.9, zorder=5)

        if False:
            good = np.where(fibertotphot['abmag'] > 0)[0]
            if len(good) > 0:
                ax.scatter(fibertotphot['lambda_eff'][good]/1e4, fibertotphot['abmag'][good],
                            marker='^', s=150, facecolor='orange', edgecolor='k',
                            label=r'$grz$ (total fiberflux)', alpha=0.9, zorder=6)
                
        #if synthphot:
        #    ax.scatter(synthmodelphot['lambda_eff']/1e4, synthmodelphot['abmag'], 
        #               marker='s', s=175, color='green', #edgecolor='k',
        #               label=r'$grz$ (spectrum, synthesized)', alpha=0.7, zorder=3)
        #    ax.scatter(synthphot['lambda_eff']/1e4, synthphot['abmag'], 
        #               marker='o', s=130, color='blue', edgecolor='k',
        #               label=r'$grz$ (spectral model, synthesized)', alpha=1.0, zorder=4)

        leg = ax.legend(loc='lower left', fontsize=16)
        #for hndl in leg.legendHandles:
        #    hndl.set_markersize(8)

        if coadd_type == 'healpix':
            targetid_str = str(metadata['TARGETID'])
        else:
            targetid_str = '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),

        leg = {
            'targetid': targetid_str,
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastphot['CONTINUUM_CHI2']),
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            #'zfastfastphot': '$z_{{\\rm fastfastphot}}$={:.6f}'.format(fastphot['CONTINUUM_Z']),
            #'z': '$z$={:.6f}'.format(fastphot['CONTINUUM_Z']),
            'age': '<Age>={:.3f} Gyr'.format(fastphot['CONTINUUM_AGE']),
            'absmag_r': '$M_{{^{{0.0}}r}}={:.2f}$'.format(fastphot['ABSMAG_SDSS_R']),
            'absmag_gr': '$^{{0.0}}(g-r)={:.3f}$'.format(fastphot['ABSMAG_SDSS_G']-fastphot['ABSMAG_SDSS_R']),
            }
        if fastphot['CONTINUUM_AV_IVAR'] == 0:
            leg.update({'AV': '$A(V)$={:.2f} mag'.format(fastphot['CONTINUUM_AV'])})
        else:
            leg.update({'AV': '$A(V)$={:.2f}+/-{:.2f} mag'.format(
                fastphot['CONTINUUM_AV'], 1/np.sqrt(fastphot['CONTINUUM_AV_IVAR']))})

        bbox = dict(boxstyle='round', facecolor='gray', alpha=0.25)
        legfntsz = 20

        txt = '\n'.join((
            r'{}'.format(leg['absmag_r']),
            r'{}'.format(leg['absmag_gr'])
            ))
        
        legxpos, legypos = 0.04, 0.94
        ax.text(legxpos, legypos, txt, ha='left', va='top',
                transform=ax.transAxes, fontsize=legfntsz,
                bbox=bbox)
        
        txt = '\n'.join((
            r'{}'.format(leg['zredrock']),
            r'{} {}'.format(leg['chi2'], leg['age']),
            r'{}'.format(leg['AV']),
            ))
        
        legxpos, legypos = 0.98, 0.06
        ax.text(legxpos, legypos, txt, ha='right', va='bottom',
                transform=ax.transAxes, fontsize=legfntsz,
                bbox=bbox)

        plt.subplots_adjust(bottom=0.14, right=0.95, top=0.93)

        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()
