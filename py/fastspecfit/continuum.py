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
from desiutil.log import get_logger, DEBUG
log = get_logger(DEBUG)

def _fnnls_continuum(myargs):
    """Multiprocessing wrapper for fnnls_continuum."""
    return fnnls_continuum(*myargs)

def fnnls_continuum(ZZ, xx, flux=None, ivar=None, modelflux=None, get_chi2=False):
    """Fit a stellar continuum using fNNLS. 

    Parameters
    ----------
    ZZ : :class:`~numpy.ndarray`
        Array.
    xx : :class:`~numpy.ndarray`
        Array.
    flux : :class:`~numpy.ndarray`, optional, defaults to ``None``
        Input flux spectrum. 
    ivar : :class:`~numpy.ndarray`, optional, defaults to ``None``
        Input inverse variance spectrum corresponding to ``flux``.
    modelflux : :class:`~numpy.ndarray`, optional, defaults to ``None``
        Input model flux spectrum. 

    Returns
    -------
    :class:`bool`
        Boolean flag indicating whether the non-negative fit did not converge.
    :class:`~numpy.ndarray`
        Coefficients of the best-fitting spectrum.
    :class:`float`
        Reduced chi-squared of the fit. Only returned if ``get_chi2=True``.

    Notes
    -----
    - This function is a simple wrapper on fastspecfit.fnnls.fnnls(); see the
      ContinuumFit.fnnls_continuum method for documentation.
    - The arguments ``flux``, ``ivar`` and ``modelflux`` are only used when
      ``get_chi2=True``.

    """
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
    """Tools for dealing with stellar continua.

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

    .. note::
        Need to document all the attributes.

    """
    def __init__(self, metallicity='Z0.0190', minwave=None, maxwave=6e4, 
                 seed=1, mapdir=None):

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
        if mapdir is None:
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
        self.sspflux = sspflux # no dust, no velocity broadening [npix,nage]
        self.sspinfo = sspinfo
        self.nage = nage
        self.npix = npix

        # emission lines
        self.linetable = read_emlines()

        self.linemask_sigma_narrow = 200.0  # [km/s]
        self.linemask_sigma_balmer = 1000.0 # [km/s]
        self.linemask_sigma_broad = 2000.0  # [km/s]

        # photometry
        self.bands = np.array(['g', 'r', 'z', 'W1', 'W2'])
        self.synth_bands = np.array(['g', 'r', 'z']) # for synthesized photometry
        self.fiber_bands = np.array(['g', 'r', 'z']) # for fiber fluxes

        self.decam = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z')
        self.bassmzls = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z')

        self.decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z',
                                              'wise2010-W1', 'wise2010-W2')
        self.bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z',
                                                 'wise2010-W1', 'wise2010-W2')

        self.bands_to_fit = np.ones(len(self.bands), bool)
        self.bands_to_fit[self.bands == 'W2'] = False # drop W2

        # rest-frame filters
        self.absmag_bands = ['U', 'B', 'V', 'sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z', 'W1']
        self.absmag_filters = filters.load_filters('bessell-U', 'bessell-B', 'bessell-V',
                                                   'sdss2010-u', 'sdss2010-g', 'sdss2010-r',
                                                   'sdss2010-i', 'sdss2010-z', 'wise2010-W1')
        self.absmag_bandshift = np.array([0.0, 0.0, 0.0,
                                          0.0, 0.0, 0.0,
                                          0.0, 0.0, 0.0])
        #self.absmag_bandshift = np.array([0.0, 0.0, 0.0,
        #                                  0.1, 0.1, 0.1,
        #                                  0.1, 0.1, 0.0])
        
        self.min_uncertainty = np.array([0.01, 0.01, 0.01, 0.02, 0.02])

        # used in one place...
        #self.rand = np.random.RandomState(seed=seed)

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
                         nanomaggies=True, nsigma=1.0, min_uncertainty=None):
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

        if min_uncertainty is None:
            min_uncertainty = np.zeros_like(maggies)

        # Gaia-only targets can sometimes have grz=-99.
        if np.any(ivarmaggies < 0) or np.any(maggies == -99.0):
            errmsg = 'All ivarmaggies must be zero or positive!'
            log.critical(errmsg)
            raise ValueError(errmsg)

        phot['lambda_eff'] = lambda_eff#.astype('f4')
        if nanomaggies:
            phot['nanomaggies'] = maggies#.astype('f4')
            phot['nanomaggies_ivar'] = ivarmaggies#.astype('f4')
        else:
            phot['nanomaggies'] = (maggies * 1e9)#.astype('f4')
            phot['nanomaggies_ivar'] = (ivarmaggies * 1e-18)#.astype('f4')

        if nanomaggies:
            nanofactor = 1e-9 # [nanomaggies-->maggies]
        else:
            nanofactor = 1.0

        #pdb.set_trace()
        #factor=(2.5/alog(10.))
        #err=factor/sqrt(maggies_ivar[k,igood])/maggies[k,igood]
        #err2=err^2+minerrors[k]^2
        #maggies_ivar[k,igood]=factor^2/(maggies[k,igood]^2*err2)

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
            phot['abmag'][igood, jgood] = (-2.5 * np.log10(nanofactor * goodmaggies))#.astype('f4')
        
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
            phot['abmag_brighterr'][igood, jgood] = errmaggies / (0.4 * np.log(10) * (maggies+errmaggies))#.astype('f4') # bright end (flux upper limit)
            phot['abmag_fainterr'][igood, jgood] = errmaggies / (0.4 * np.log(10) * (maggies-errmaggies))#.astype('f4') # faint end (flux lower limit)
            #phot['abmag_loerr'][igood, jgood] = +2.5 * np.log10(1 + fracerr) # bright end (flux upper limit)
            #phot['abmag_uperr'][igood, jgood] = +2.5 * np.log10(1 - fracerr) # faint end (flux lower limit)
            #test = 2.5 * np.log(np.exp(1)) * fracerr # symmetric in magnitude (approx)

            # approximate the uncertainty as being symmetric in magnitude
            phot['abmag_ivar'][igood, jgood] = (ivarmaggies * (maggies * 0.4 * np.log(10))**2)#.astype('f4')
            
        if len(upper) > 0:
            phot['abmag_limit'][iupper, jupper] = abmag_limit#.astype('f4')
            
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

    def smooth_continuum(self, wave, flux, ivar, redshift, medbin=150, 
                         smooth_window=50, smooth_step=10, maskkms_uv=3000.0, 
                         maskkms_balmer=1000.0, maskkms_narrow=200.0, 
                         linemask=None, png=None):
        """Build a smooth, nonparametric continuum spectrum.

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
            Observed-frame wavelength array.
        flux : :class:`numpy.ndarray` [npix]
            Spectrum corresponding to `wave`.
        ivar : :class:`numpy.ndarray` [npix]
            Inverse variance spectrum corresponding to `flux`.
        redshift : :class:`float`
            Object redshift.
        medbin : :class:`int`, optional, defaults to 150 pixels
            Width of the median-smoothing kernel in pixels; a magic number.
        smooth_window : :class:`int`, optional, defaults to 50 pixels
            Width of the sliding window used to compute the iteratively clipped
            statistics (mean, median, sigma); a magic number. Note: the nominal
            extraction width (0.8 A) and observed-frame wavelength range
            (3600-9800 A) corresponds to pixels that are 66-24 km/s. So
            `smooth_window` of 50 corresponds to 3300-1200 km/s, which is
            conservative for all but the broadest lines. A magic number. 
        smooth_step : :class:`int`, optional, defaults to 10 pixels
            Width of the step size when computing smoothed statistics; a magic
            number.
        maskkms_uv : :class:`float`, optional, defaults to 3000 km/s
            Masking width for UV emission lines. Pixels within +/-3*maskkms_uv
            are masked before median-smoothing.
        maskkms_balmer : :class:`float`, optional, defaults to 3000 km/s
            Like `maskkms_uv` but for Balmer lines.
        maskkms_narrow : :class:`float`, optional, defaults to 300 km/s
            Like `maskkms_uv` but for narrow, forbidden lines.
        linemask : :class:`numpy.ndarray` of type :class:`bool`, optional, defaults to `None`
            Boolean mask with the same number of pixels as `wave` where `True`
            means a pixel is (possibly) affected by an emission line
            (specifically a strong line which likely cannot be median-smoothed).
        png : :class:`str`, optional, defaults to `None`
            Generate a simple QA plot and write it out to this filename.

        Returns
        -------
        smooth :class:`numpy.ndarray` [npix]
            Smooth continuum spectrum which can be subtracted from `flux` in
            order to create a pure emission-line spectrum.
        smoothsigma :class:`numpy.ndarray` [npix]
            Smooth one-sigma uncertainty spectrum.

        """
        from scipy.ndimage import median_filter
        from numpy.lib.stride_tricks import sliding_window_view
        from astropy.stats import sigma_clip

        npix = len(wave)

        # If we're not given a linemask, make a conservative one.
        if linemask is None:
            linemask = np.zeros(npix, bool) # True = (possibly) affected by emission line

            nsig = 3

            # select just strong lines
            zlinewaves = self.linetable['restwave'] * (1 + redshift)
            inrange = (zlinewaves > np.min(wave)) * (zlinewaves < np.max(wave))
            if np.sum(inrange) > 0:
                linetable = self.linetable[inrange]
                linetable = linetable[linetable['amp'] >= 1]
                if len(linetable) > 0:
                    for oneline in linetable:
                        zlinewave = oneline['restwave'] * (1 + redshift)
                        if oneline['isbroad']:
                            if oneline['isbalmer']:
                                sigma = maskkms_balmer
                            else:
                                sigma = maskkms_uv
                        else:
                            sigma = maskkms_narrow
                    
                        sigma *= zlinewave / C_LIGHT # [km/s --> Angstrom]
                        I = (wave >= (zlinewave - nsig*sigma)) * (wave <= (zlinewave + nsig*sigma))
                        if len(I) > 0:
                            linemask[I] = True

            # ToDo: mask Ly-a (1215 A) here.

        if len(linemask) != npix:
            errmsg = 'Linemask must have the same number of pixels as the input spectrum.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        # Build the smooth (line-free) continuum by computing statistics in a
        # sliding window, accounting for masked pixels and trying to be smart
        # about broad lines. See:
        #   https://stackoverflow.com/questions/41851044/python-median-filter-for-1d-numpy-array
        #   https://numpy.org/devdocs/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html
        
        wave_win = sliding_window_view(wave, window_shape=smooth_window)
        flux_win = sliding_window_view(flux, window_shape=smooth_window)
        ivar_win = sliding_window_view(ivar, window_shape=smooth_window)
        noline_win = sliding_window_view(np.logical_not(linemask), window_shape=smooth_window)

        smooth_wave, smooth_flux, smooth_sigma, smooth_mask = [], [], [], []
        for swave, sflux, sivar, noline in zip(wave_win[::smooth_step],
                                               flux_win[::smooth_step],
                                               ivar_win[::smooth_step],
                                               noline_win[::smooth_step]):

            # if there are fewer than 10 good pixels after accounting for the
            # line-mask, skip this window.
            sflux = sflux[noline]
            if len(sflux) < 10:
                smooth_mask.append(True)
                continue
            swave = swave[noline]
            sivar = sivar[noline]

            cflux = sigma_clip(sflux, sigma=2.0, cenfunc='median', stdfunc='std', masked=False, grow=1.5)
            if np.sum(np.isfinite(cflux)) < 10:
                smooth_mask.append(True)
                continue

            I = np.isfinite(cflux) # should never be fully masked!
            smooth_wave.append(np.mean(swave[I]))
            smooth_mask.append(False)

            # simple median and sigma
            sig = np.std(cflux[I])
            mn = np.median(cflux[I])

            ## inverse-variance weighted mean and sigma
            #norm = np.sum(sivar[I])
            #mn = np.sum(sivar[I] * cflux[I]) / norm # weighted mean
            #sig = np.sqrt(np.sum(sivar[I] * (cflux[I] - mn)**2) / norm) # weighted sigma

            smooth_sigma.append(sig)
            smooth_flux.append(mn)

        smooth_wave = np.array(smooth_wave)
        smooth_sigma = np.array(smooth_sigma)
        smooth_flux = np.array(smooth_flux)
        smooth_mask = np.array(smooth_mask)

        smooth_flux = np.interp(wave, smooth_wave, smooth_flux)
        smooth_sigma = np.interp(wave, smooth_wave, smooth_sigma)

        smooth = median_filter(smooth_flux, medbin, mode='nearest')
        smoothsigma = median_filter(smooth_sigma, medbin, mode='nearest')

        # Optional QA.
        if png:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
            ax[0].plot(wave, flux)
            ax[0].scatter(wave[linemask], flux[linemask], s=10, marker='s', color='k', zorder=2)
            ax[0].plot(wave, smooth, color='red')

            ax[1].plot(wave, flux - smooth)
            ax[1].axhline(y=0, color='k')

            #for xx in ax:
                #xx.set_xlim(3800, 4300)
                #xx.set_xlim(5200, 6050)
                #xx.set_xlim(4500, 5500)
                #xx.set_xlim(6300, 6900)
            #for xx in ax:
            #    xx.set_ylim(-5, 100)
            zlinewaves = self.linetable['restwave'] * (1 + redshift)
            linenames = self.linetable['name']
            inrange = np.where((zlinewaves > np.min(wave)) * (zlinewaves < np.max(wave)))[0]
            if len(inrange) > 0:
                for linename, zlinewave in zip(linenames[inrange], zlinewaves[inrange]):
                    #print(linename, zlinewave)
                    for xx in ax:
                        xx.axvline(x=zlinewave, color='gray')
    
            fig.savefig(png)

        return smooth, smoothsigma
    
    def estimate_linesigmas(self, wave, flux, ivar, redshift=0.0, png=None, refit=True):
        """Estimate the velocity width from potentially strong, isolated lines.
    
        """
        def get_linesigma(zlinewaves, init_linesigma, label='Line', ax=None):
    
            from scipy.optimize import curve_fit
            
            linesigma, linesigma_snr = 0.0, 0.0

            if ax:
                _empty = True
            
            inrange = (zlinewaves > np.min(wave)) * (zlinewaves < np.max(wave))
            if np.sum(inrange) > 0:
                stackdvel, stackflux, stackivar, contflux = [], [], [], []
                for zlinewave in zlinewaves[inrange]:
                    I = ((wave >= (zlinewave - 5*init_linesigma * zlinewave / C_LIGHT)) *
                         (wave <= (zlinewave + 5*init_linesigma * zlinewave / C_LIGHT)) *
                         (ivar > 0))
                    J = np.logical_or(
                        (wave > (zlinewave - 8*init_linesigma * zlinewave / C_LIGHT)) *
                        (wave < (zlinewave - 5*init_linesigma * zlinewave / C_LIGHT)) *
                        (ivar > 0),
                        (wave < (zlinewave + 8*init_linesigma * zlinewave / C_LIGHT)) *
                        (wave > (zlinewave + 5*init_linesigma * zlinewave / C_LIGHT)) *
                        (ivar > 0))
    
                    if (np.sum(I) > 3) and np.max(flux[I]*ivar[I]) > 1:
                        stackdvel.append((wave[I] - zlinewave) / zlinewave * C_LIGHT)
                        norm = np.percentile(flux[I], 99)
                        if norm <= 0:
                            norm = 1.0
                        stackflux.append(flux[I] / norm)
                        stackivar.append(ivar[I] * norm**2)
                        if np.sum(J) > 3:
                            contflux.append(flux[J] / norm) # continuum pixels
                            #contflux.append(np.std(flux[J]) / norm) # error in the mean
                            #contflux.append(np.std(flux[J]) / np.sqrt(np.sum(J)) / norm) # error in the mean
                        else:
                            contflux.append(flux[I] / norm) # shouldn't happen...
                            #contflux.append(np.std(flux[I]) / norm) # shouldn't happen...
                            #contflux.append(np.std(flux[I]) / np.sqrt(np.sum(I)) / norm) # shouldn't happen...
    
                if len(stackflux) > 0: 
                    stackdvel = np.hstack(stackdvel)
                    stackflux = np.hstack(stackflux)
                    stackivar = np.hstack(stackivar)
                    contflux = np.hstack(contflux)
    
                    if len(stackflux) > 10: # require at least 10 pixels
                        #onegauss = lambda x, amp, sigma: amp * np.exp(-0.5 * x**2 / sigma**2) # no pedestal
                        onegauss = lambda x, amp, sigma, const: amp * np.exp(-0.5 * x**2 / sigma**2) + const
                        #onegauss = lambda x, amp, sigma, const, slope: amp * np.exp(-0.5 * x**2 / sigma**2) + const + slope*x
        
                        stacksigma = 1 / np.sqrt(stackivar)
                        try:
                            popt, _ = curve_fit(onegauss, xdata=stackdvel, ydata=stackflux,
                                                sigma=stacksigma, p0=[1.0, init_linesigma, 0.0])
                                                #sigma=stacksigma, p0=[1.0, init_linesigma, np.median(stackflux)])
                                                #sigma=stacksigma, p0=[1.0, sigma, np.median(stackflux), 0.0])
                            popt[1] = np.abs(popt[1])
                            if popt[0] > 0 and popt[1] > 0:
                                linesigma = popt[1]
                                robust_std = np.diff(np.percentile(contflux, [25, 75]))[0] / 1.349 # robust sigma
                                #robust_std = np.std(contflux)
                                if robust_std > 0:
                                    linesigma_snr = popt[0] / robust_std
                                else:
                                    linesigma_snr = 0.0
                            else:
                                popt = None
                        except RuntimeError:
                            popt = None

                        if ax:
                            _label = r'{} $\sigma$={:.0f} km/s S/N={:.1f}'.format(label, linesigma, linesigma_snr)
                            ax.scatter(stackdvel, stackflux, s=10, label=_label)
                            if popt is not None:
                                srt = np.argsort(stackdvel)
                                linemodel = onegauss(stackdvel[srt], *popt)
                                ax.plot(stackdvel[srt], linemodel, color='k')#, label='Gaussian Model')
                            else:
                                linemodel = stackflux * 0

                            #_min, _max = np.percentile(stackflux, [5, 95])
                            _max = np.max([np.max(linemodel), 1.05*np.percentile(stackflux, 99)])

                            ax.set_ylim(-2*np.median(contflux), _max)
                            if linesigma > 0:
                                if linesigma < np.max(stackdvel):
                                    ax.set_xlim(-5*linesigma, +5*linesigma)

                            ax.set_xlabel(r'$\Delta v$ (km/s)')
                            ax.set_ylabel('Relative Flux')
                            ax.legend(loc='upper left', fontsize=8, frameon=False)
                            _empty = False
                        
            log.debug('{} masking sigma={:.3f} and S/N={:.3f}'.format(label, linesigma, linesigma_snr))

            if ax and _empty:
                ax.plot([0, 0], [0, 0], label='{}-No Data'.format(label))
                ax.axes.xaxis.set_visible(False)
                ax.axes.yaxis.set_visible(False)
                ax.legend(loc='upper left', fontsize=10)
            
            return linesigma, linesigma_snr
    
        linesigma_snr_min = 1.5
        init_linesigma_balmer = 1000.0
        init_linesigma_narrow = 200.0
        init_linesigma_uv = 2000.0

        if png:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 3, figsize=(12, 4))
        else:
            ax = [None] * 3
            
        # [OII] doublet, [OIII] 4959,5007
        zlinewaves = np.array([3728.483, 4960.295, 5008.239]) * (1 + redshift)
        linesigma_narrow, linesigma_narrow_snr = get_linesigma(zlinewaves, init_linesigma_narrow, 
                                                               label='Forbidden', #label='[OII]+[OIII]',
                                                               ax=ax[0])

        # refit with the new value
        if refit:
            if (linesigma_narrow > init_linesigma_narrow) and (linesigma_narrow < 5*init_linesigma_narrow) and (linesigma_narrow_snr > linesigma_snr_min):
                if ax[0] is not None:
                    ax[0].clear()
                linesigma_narrow, linesigma_narrow_snr = get_linesigma(
                    zlinewaves, linesigma_narrow, label='Forbidden', ax=ax[0])

        if (linesigma_narrow < 50) or (linesigma_narrow > 5*init_linesigma_narrow) or (linesigma_narrow_snr < linesigma_snr_min):
            linesigma_narrow_snr = 0.0
            linesigma_narrow = init_linesigma_narrow
    
        # Hbeta, Halpha
        zlinewaves = np.array([4862.683, 6564.613]) * (1 + redshift)
        linesigma_balmer, linesigma_balmer_snr = get_linesigma(zlinewaves, init_linesigma_balmer, 
                                                               label='Balmer', #label=r'H$\alpha$+H$\beta$',
                                                               ax=ax[1])
        # refit with the new value
        if refit:
            if (linesigma_balmer > init_linesigma_balmer) and (linesigma_balmer < 5*init_linesigma_balmer) and (linesigma_balmer_snr > linesigma_snr_min): 
                if ax[1] is not None:
                    ax[1].clear()
                linesigma_balmer, linesigma_balmer_snr = get_linesigma(zlinewaves, linesigma_balmer, 
                                                                       label='Balmer', #label=r'H$\alpha$+H$\beta$',
                                                                       ax=ax[1])
                
        # if no good fit, should we use narrow or Balmer??
        if (linesigma_balmer < 50) or (linesigma_balmer > 5*init_linesigma_balmer) or (linesigma_balmer_snr < linesigma_snr_min):
            linesigma_balmer_snr = 0.0
            linesigma_balmer = init_linesigma_balmer
            #linesigma_balmer = init_linesigma_narrow 
    
        # Lya, SiIV doublet, CIV doublet, CIII], MgII doublet
        zlinewaves = np.array([1549.4795, 2799.942]) * (1 + redshift)
        #zlinewaves = np.array([1215.670, 1398.2625, 1549.4795, 1908.734, 2799.942]) * (1 + redshift)
        linesigma_uv, linesigma_uv_snr = get_linesigma(zlinewaves, init_linesigma_uv, 
                                                       label='UV/Broad', ax=ax[2])
        
        # refit with the new value
        if refit:
            if (linesigma_uv > init_linesigma_uv) and (linesigma_uv < 5*init_linesigma_uv) and (linesigma_uv_snr > linesigma_snr_min): 
                if ax[2] is not None:
                    ax[2].clear()
                linesigma_uv, linesigma_uv_snr = get_linesigma(zlinewaves, linesigma_uv, 
                                                               label='UV/Broad', ax=ax[2])
                
        if (linesigma_uv < 300) or (linesigma_uv > 5*init_linesigma_uv) or (linesigma_uv_snr < linesigma_snr_min):
            linesigma_uv_snr = 0.0
            linesigma_uv = init_linesigma_uv

        if png:
            fig.subplots_adjust(left=0.1, bottom=0.15, wspace=0.2, right=0.95, top=0.95)
            fig.savefig(png)

        return (linesigma_narrow, linesigma_balmer, linesigma_uv,
                linesigma_narrow_snr, linesigma_balmer_snr, linesigma_uv_snr)

    def build_linemask(self, wave, flux, ivar, redshift=0.0, nsig=5.0):
        """Generate a mask which identifies pixels impacted by emission lines.

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
            Observed-frame wavelength array.
        flux : :class:`numpy.ndarray` [npix]
            Spectrum corresponding to `wave`.
        ivar : :class:`numpy.ndarray` [npix]
            Inverse variance spectrum corresponding to `flux`.
        redshift : :class:`float`, optional, defaults to 0.0
            Object redshift.
        nsig : :class:`float`, optional, defaults to 5.0
            Mask pixels which are within +/-`nsig`-sigma of each emission line,
            where `sigma` is the line-width in km/s.

        Returns
        -------
        smooth :class:`numpy.ndarray` [npix]
            Smooth continuum spectrum which can be subtracted from `flux` in
            order to create a pure emission-line spectrum.

        
        :class:`dict` with the following keys:
            linemask : :class:`list`
            linemask_all : :class:`list`
            linename : :class:`list`
            linepix : :class:`list`
            contpix : :class:`list`

        Notes
        -----
        Code exists to generate some QA but the code is not exposed.

        """
        # Initially, mask aggressively, especially the Balmer lines.
        png = None
        #png = 'smooth.png'
        #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/smooth.png'
        smooth, smoothsigma = self.smooth_continuum(wave, flux, ivar, redshift, maskkms_uv=5000.0,
                                                    maskkms_balmer=5000.0, maskkms_narrow=500.0,
                                                    png=png)

        # Get a better estimate of the Balmer, forbidden, and UV/QSO line-widths.
        png = None
        #png = 'linesigma.png'
        #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/linesigma.png'
        linesigma_narrow, linesigma_balmer, linesigma_uv, linesigma_narrow_snr, linesigma_balmer_snr, linesigma_uv_snr = \
          self.estimate_linesigmas(wave, flux-smooth, ivar, redshift, png=png)

        # Next, build the emission-line mask.
        linemask = np.zeros_like(wave, bool)      # True = affected by possible emission line.
        linemask_strong = np.zeros_like(linemask) # True = affected by strong emission lines.
    
        linenames = np.hstack(('Lya', self.linetable['name'])) # include Lyman-alpha
        zlinewaves = np.hstack((1215.0, self.linetable['restwave'])) * (1 + redshift)
        lineamps = np.hstack((1.0, self.linetable['amp']))
        isbroads = np.hstack((True, self.linetable['isbroad'] * (self.linetable['isbalmer'] == False)))
        isbalmers = np.hstack((False, self.linetable['isbalmer'] * (self.linetable['isbroad'] == False)))
    
        png = None
        #png = 'linemask.png'
        #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/linemask.png'
        snr_strong = 3.0
    
        inrange = (zlinewaves > np.min(wave)) * (zlinewaves < np.max(wave))
        nline = np.sum(inrange)
        if nline > 0:
            # Index I for building the line-mask; J for estimating the local
            # continuum (to be used in self.smooth_continuum).
    
            # initial line-mask
            for _linename, zlinewave, lineamp, isbroad, isbalmer in zip(
                    linenames[inrange], zlinewaves[inrange], lineamps[inrange],
                    isbroads[inrange], isbalmers[inrange]):
                if isbroad:
                    linesigma = linesigma_uv
                elif isbalmer or 'broad' in _linename:
                    linesigma = linesigma_balmer
                else:
                    linesigma = linesigma_narrow
                    
                sigma = linesigma * zlinewave / C_LIGHT # [km/s --> Angstrom]
                I = (wave >= (zlinewave - nsig*sigma)) * (wave <= (zlinewave + nsig*sigma))
                if np.sum(I) > 0:
                    linemask[I] = True

                    #if 'broad' in _linename:
                    #    linemask_strong[I] = True

                    # Now find "strong" lines using a constant sigma so that we
                    # focus on the center of the wavelength range of the
                    # line. For example: if sigma is too big, like 2000 km/s
                    # then we can sometimes be tricked into flagging lines (like
                    # the Helium lines) as strong because we pick up the high
                    # S/N of adjacent truly strong lines.
                    if linesigma > 500.0:
                        linesigma_strong = 500.0
                    else:
                        linesigma_strong = linesigma
                    sigma_strong = linesigma_strong * zlinewave / C_LIGHT # [km/s --> Angstrom]
                    J = (ivar > 0) * (smoothsigma > 0) * (wave >= (zlinewave - sigma_strong)) * (wave <= (zlinewave + sigma_strong))
                    if np.sum(J) > 0:
                        snr = (flux[J] - smooth[J]) / smoothsigma[J]
                        # require peak S/N>3 and at least 5 pixels with S/N>3
                        #print(_linename, zlinewave, np.percentile(snr, 98), np.sum(snr > snr_strong))
                        if np.percentile(snr, 98) > snr_strong and np.sum(snr > snr_strong) > 5:
                            linemask_strong[I] = True
                        #if '4686' in _linename:
                        #    pdb.set_trace()

            # now get the continuum, too
            if png:
                import matplotlib.pyplot as plt
                nrows = np.ceil(nline/4).astype(int)
                fig, ax = plt.subplots(nrows, 4, figsize=(8, 2*nrows))
                ax = ax.flatten()
            else:
                ax = [None] * nline

            linepix, contpix, linename = [], [], []        
            for _linename, zlinewave, lineamp, isbroad, isbalmer, xx in zip(
                    linenames[inrange], zlinewaves[inrange], lineamps[inrange],
                    isbroads[inrange], isbalmers[inrange], ax):
                
                if isbroad:
                    sigma = linesigma_uv
                elif isbalmer or 'broad' in _linename:
                    sigma = linesigma_balmer
                else:
                    sigma = linesigma_narrow

                sigma *= zlinewave / C_LIGHT # [km/s --> Angstrom]
                I = (wave >= (zlinewave - nsig*sigma)) * (wave <= (zlinewave + nsig*sigma))

                # get the pixels of the local continuum
                Jblu = (wave > (zlinewave - 2*nsig*sigma)) * (wave < (zlinewave - nsig*sigma)) * (linemask_strong == False)
                Jred = (wave < (zlinewave + 2*nsig*sigma)) * (wave > (zlinewave + nsig*sigma)) * (linemask_strong == False)
                J = np.logical_or(Jblu, Jred)

                if np.sum(J) < 10: # go further out
                    Jblu = (wave > (zlinewave - 3*nsig*sigma)) * (wave < (zlinewave - nsig*sigma)) * (linemask_strong == False)
                    Jred = (wave < (zlinewave + 3*nsig*sigma)) * (wave > (zlinewave + nsig*sigma)) * (linemask_strong == False)
                    J = np.logical_or(Jblu, Jred)
                
                if np.sum(J) < 10: # drop the linemask_ condition
                    Jblu = (wave > (zlinewave - 2*nsig*sigma)) * (wave < (zlinewave - nsig*sigma))
                    Jred = (wave < (zlinewave + 2*nsig*sigma)) * (wave > (zlinewave + nsig*sigma))
                    J = np.logical_or(Jblu, Jred)

                #print(_linename, np.sum(I), np.sum(J))
                if np.sum(I) > 0 and np.sum(J) > 0:
                    linename.append(_linename)
                    linepix.append(I)
                    contpix.append(J)
    
                    if png:
                        _Jblu = np.where((wave > (zlinewave - 2*nsig*sigma)) * (wave < (zlinewave - nsig*sigma)))[0]
                        _Jred = np.where((wave < (zlinewave + 2*nsig*sigma)) * (wave > (zlinewave + nsig*sigma)))[0]
                        if len(_Jblu) == 0:
                            _Jblu = [0]
                        if len(_Jred) == 0:
                            _Jred = [len(wave)-1]
                        plotwave, plotflux = wave[_Jblu[0]:_Jred[-1]], flux[_Jblu[0]:_Jred[-1]]
                        xx.plot(plotwave, plotflux, label=_linename)
                        #xx.plot(np.hstack((wave[Jblu], wave[I], wave[Jred])), 
                        #        np.hstack((flux[Jblu], flux[I], flux[Jred])), label=_linename)
                        #xx.plot(wave[I], flux[I], label=_linename)
                        xx.scatter(wave[I], flux[I], s=10, color='orange', marker='s')
                        if np.sum(linemask_strong[I]) > 0:
                            xx.scatter(wave[I][linemask_strong[I]], flux[I][linemask_strong[I]], s=10, color='green', marker='s')
                        #if '4471' in _linename:
                        #    pdb.set_trace()
                        
                        xx.scatter(wave[Jblu], flux[Jblu], color='blue', s=10)
                        xx.scatter(wave[Jred], flux[Jred], color='red', s=10)
                        xx.set_ylim(np.min(plotflux), np.max(flux[I]))
                        xx.legend(frameon=False, fontsize=10, loc='upper left')
        
            linemask_dict = {'linemask_all': linemask, 'linemask': linemask_strong,
                             'linename': linename, 'linepix': linepix, 'contpix': contpix,
                             'linesigma_narrow': linesigma_narrow, 'linesigma_narrow_snr': linesigma_narrow_snr, 
                             'linesigma_balmer': linesigma_balmer, 'linesigma_balmer_snr': linesigma_balmer_snr, 
                             'linesigma_uv': linesigma_uv, 'linesigma_uv_snr': linesigma_uv_snr, 
                             }
    
            if png:
                fig.savefig(png)
    
        else:
            linemask_dict = {'linemask_all': [], 'linemask': [],
                             'linename': [], 'linepix': [], 'contpix': []}

        return linemask_dict

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
            errmsg = 'Input SSPs have an unrecognized number of dimensions, {}'.format(ndim)
            log.critical(errmsg)
            raise ValueError(errmsg)
        
        #t0 = time.time()
        ##sspflux = sspflux.copy().reshape(npix, nmodel)
        #log.info('Copying the data took: {:.2f} sec'.format(time.time()-t0))

        # optionally apply reddening
        if AV:
            atten = self.dust_attenuation(sspwave, AV)
            sspflux *= atten[:, np.newaxis]

        # broaden for velocity dispersion
        if vdisp:
            sspflux = self.convolve_vdisp(sspflux, vdisp)

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

            ## This is wrong...
            #if vdisp:
            #     datasspflux = self.convolve_vdisp(datasspflux, vdisp)
                 
            # optionally compute the best-fitting model
            if coeff is not None:
                datasspflux = datasspflux.dot(coeff)
                if synthphot:
                    maggies = filters.get_ab_maggies(datasspflux, zsspwave, axis=0)
                    maggies = np.array(maggies.as_array().tolist()[0])
                    maggies /= self.fluxnorm * self.massnorm
                    sspphot = self.parse_photometry(self.bands, maggies, effwave, nanomaggies=False)
        else:
            # loop over cameras
            datasspflux = []
            for icamera in np.arange(len(cameras)): # iterate on cameras
                _datasspflux = []
                for imodel in np.arange(nmodel):
                    _datasspflux.append(self.smooth_and_resample(
                        zsspflux[:, imodel], zsspwave, specwave=specwave[icamera],
                        specres=specres[icamera]))
                _datasspflux = np.vstack(_datasspflux).T
                ## This is wrong...
                #if vdisp:
                #    _datasspflux = self.convolve_vdisp(_datasspflux, vdisp)
                if coeff is not None:
                    _datasspflux = _datasspflux.dot(coeff)
                datasspflux.append(_datasspflux)
                
        #log.info('Resampling took: {:.2f} sec'.format(time.time()-t0))

        return datasspflux, sspphot # vector or 3-element list of [npix,nmodel] spectra

class ContinuumFit(ContinuumTools):
    def __init__(self, metallicity='Z0.0190', minwave=None, maxwave=6e4, 
                 nolegend=False, solve_vdisp=False, mapdir=None):
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
        super(ContinuumFit, self).__init__(metallicity=metallicity, minwave=minwave, 
                                           maxwave=maxwave, mapdir=mapdir)

        self.nolegend = nolegend

        # Initialize the velocity dispersion and reddening parameters. Make sure
        # the nominal values are in the grid.
        self.solve_vdisp = solve_vdisp
        
        #vdispmin, vdispmax, dvdisp, vdisp_nominal = (20.0, 300.0, 20.0, 75.0)
        vdispmin, vdispmax, dvdisp, vdisp_nominal = (100.0, 350.0, 20.0, 150.0)
        #vdispmin, vdispmax, dvdisp, vdisp_nominal = (0.0, 0.0, 30.0, 150.0)
        nvdisp = int(np.ceil((vdispmax - vdispmin) / dvdisp))
        if nvdisp == 0:
            nvdisp = 1
        vdisp = np.linspace(vdispmin, vdispmax, nvdisp)#.astype('f4') # [km/s]

        if not vdisp_nominal in vdisp:
            vdisp = np.sort(np.hstack((vdisp, vdisp_nominal)))
        self.vdisp = vdisp
        self.vdisp_nominal = vdisp_nominal
        self.nvdisp = len(vdisp)

        #AVmin, AVmax, dAV, AV_nominal = (0.0, 0.0, 0.1, 0.0)
        AVmin, AVmax, dAV, AV_nominal = (0.0, 1.5, 0.1, 0.0)
        #AVmin, AVmax, dAV, AV_nominal = (0.0, 1.5, 0.05, 0.0)
        nAV = int(np.ceil((AVmax - AVmin) / dAV))
        if nAV == 0:
            nAV = 1
        AV = np.linspace(AVmin, AVmax, nAV)#.astype('f4')
        assert(AV[0] == 0.0) # minimum value has to be zero (assumed in fnnls_continuum)

        if not AV_nominal in AV:
            AV = np.sort(np.hstack((AV, AV_nominal)))        
        self.AV = AV
        self.AV_nominal = AV_nominal
        self.nAV = len(AV)

        # Next, precompute a grid of spectra convolved to the nominal velocity
        # dispersion with reddening applied. This isn't quite right redward of
        # ~1 micron where the pixel size changes, but fix that later.
        sspflux_dustnomvdisp = []
        for AV in self.AV:
            atten = self.dust_attenuation(self.sspwave, AV)
            _sspflux_dustnomvdisp = self.convolve_vdisp(self.sspflux * atten[:, np.newaxis], self.vdisp_nominal)
            sspflux_dustnomvdisp.append(_sspflux_dustnomvdisp)

        #import matplotlib.pyplot as plt
        #plt.plot(self.sspwave, self.sspflux[:, 5])
        #plt.plot(self.sspwave, sspflux_dustnomvdisp[3][:, 5])
        #plt.xlim(3500, 4700)
        #plt.savefig('desi-users/ioannis/tmp/test-vdisp.png')

        # nominal velocity broadening on a grid of A(V) [npix,nage,nAV]
        self.sspflux_dustnomvdisp = np.stack(sspflux_dustnomvdisp, axis=-1) # [npix,nage,nAV]

        # Finally, optionally precompute a grid of spectra with nominal
        # reddening on a grid of velocity dispersion. Again, this isn't quite
        # right redward of ~1 micron where the pixel size changes.
        if self.solve_vdisp:
            sspflux_vdispnomdust = []
            for vdisp in self.vdisp:
                _sspflux_vdispnomdust = self.convolve_vdisp(self.sspflux, vdisp)
                sspflux_vdispnomdust.append(_sspflux_vdispnomdust)

            #import matplotlib.pyplot as plt
            #plt.plot(self.sspwave, self.sspflux[:, 5])
            #plt.plot(self.sspwave, sspflux_vdispnomdust[3][:, 5])
            #plt.plot(self.sspwave, sspflux_vdispnomdust[7][:, 5])
            #plt.xlim(3500, 4700)
            #plt.savefig('desi-users/ioannis/tmp/test-vdisp.png')
    
            # nominal dust on a grid of velocity broadening [npix,nage,nvdisp]
            self.sspflux_vdispnomdust = np.stack(sspflux_vdispnomdust, axis=-1) # [npix,nage,nvdisp]

    def init_spec_output(self, nobj=1):
        """Initialize the output data table for this class.

        """
        from astropy.table import Table, Column
        
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_RCHI2', length=nobj, dtype='f4')) # reduced chi2
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
        out.add_column(Column(name='CONTINUUM_RCHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))
        
        # observed-frame photometry synthesized from the best-fitting continuum model fit
        for band in self.bands:
            out.add_column(Column(name='FLUX_SYNTH_MODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.nanomaggy))

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

    def kcorr_and_absmag(self, data, continuum, coeff, snrmin=2.0):
        """Computer K-corrections and absolute magnitudes.

        # To get the absolute r-band magnitude we would do:
        #   M_r = m_X_obs + 2.5*log10(r_synth_rest/X_synth_obs)
        # where X is the redshifted bandpass

        """
        redshift = data['zredrock']
        band_shift = self.absmag_bandshift
        bands_to_fit = self.bands_to_fit
        
        if data['photsys'] == 'S':
            filters_in = self.decamwise
        else:
            filters_in = self.bassmzlswise
        filters_out = self.absmag_filters
        nout = len(filters_out)

        # shift the bandpasses blueward by a factor of 1+band_shift
        lambda_in = filters_in.effective_wavelengths.value
        lambda_out = filters_out.effective_wavelengths.value / (1 + band_shift)

        # redshifted wavelength array and distance modulus
        zsspwave = self.sspwave * (1 + redshift)
        dmod = self.cosmo.distmod(redshift).value

        maggies = data['phot']['nanomaggies'].data * 1e-9
        ivarmaggies = (data['phot']['nanomaggies_ivar'].data / 1e-9**2) * bands_to_fit

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

        return kcorr, absmag, ivarabsmag, bestmaggies

    def _fnnls_parallel(self, modelflux, flux, ivar, xparam=None, debug=False):
        """Wrapper on fnnls to set up the multiprocessing. 

        Works with both spectroscopic and photometric input and with both 2D and
        3D model spectra.

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

        fitargs = [(ZZ[:, :, ii], xx, flux, ivar, modelflux[:, :, ii], True) for ii in np.arange(nn)]
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
            #plt.ylim(chi2min*0.95, np.max(chi2grid[imin-1:imin+2])*1.05)
            #plt.plot(xx, np.polyval([aa, bb, cc], xx), ls='--')
            plt.axvline(x=xbest, color='k')
            if xivar > 0:
                plt.axhline(y=chi2min, color='k')
            plt.yscale('log')
            plt.savefig('desi-users/ioannis/tmp/qa-chi2min.png')

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

        .. note::

            See https://github.com/mikeiovine/fast-nnls for the fNNLS algorithm(s).

        """
        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_phot_output()

        redshift = data['zredrock']
        #result['CONTINUUM_Z'] = redshift

        # Prepare the reddened and unreddened SSP templates by redshifting and
        # normalizing. Note that we ignore templates which are older than the
        # age of the universe at the galaxy redshift.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustnomvdisp, zsspphot_dustnomvdisp = self.SSP2data(
            self.sspflux_dustnomvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=None, specres=None,
            south=data['photsys'] == 'S')
        log.info('Preparing the models took {:.2f} sec'.format(time.time()-t0))
        
        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = (data['phot']['flam_ivar'].data / self.fluxnorm**2) * self.bands_to_fit
        assert(np.all(objflamivar >= 0))

        if np.all(objflamivar == 0): # can happen for secondary targets
            log.info('All photometry is masked or not available!')
            AVbest, AVivar = self.AV_nominal, 0.0
            nage = self.nage
            coeff = np.zeros(self.nage)
            continuummodel = np.zeros(len(self.sspwave))
        else:
            zsspflam_dustnomvdisp = zsspphot_dustnomvdisp['flam'].data * self.fluxnorm * self.massnorm # [nband,nage*nAV]

            inodust = np.ndarray.item(np.where(self.AV == 0)[0]) # should always be index 0

            npix, nmodel = zsspflux_dustnomvdisp.shape
            nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)

            zsspflam_dustnomvdisp = zsspflam_dustnomvdisp.reshape(len(self.bands), nage, self.nAV) # [nband,nage,nAV]

            t0 = time.time()
            AVchi2min, AVbest, AVivar = self._fnnls_parallel(zsspflam_dustnomvdisp, objflam, objflamivar,
                                                             xparam=self.AV, debug=False)
            log.info('Fitting the photometry took: {:.2f} sec'.format(time.time()-t0))
            if AVivar > 0:
                log.info('Best-fitting photometric A(V)={:.4f}+/-{:.4f} with chi2={:.3f}'.format(
                    AVbest, 1/np.sqrt(AVivar), AVchi2min))
            else:
                AVbest = self.AV_nominal
                log.info('Finding photometric A(V) failed; adopting A(V)={:.4f}'.format(self.AV_nominal))

            # Get the final set of coefficients and chi2 at the best-fitting
            # reddening and nominal velocity dispersion.
            bestsspflux, bestphot = self.SSP2data(self.sspflux_dustnomvdisp[:, agekeep, inodust], # equivalent to calling with self.sspflux[:, agekeep]
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
            synth_bestmaggies = np.zeros(len(self.bands))
        else:
            dn4000, _ = self.get_dn4000(self.sspwave, continuummodel, rest=True)
            meanage = self.get_meanage(coeff)
            kcorr, absmag, ivarabsmag, synth_bestmaggies = self.kcorr_and_absmag(data, continuummodel, coeff)

            # convert to nanomaggies
            synth_bestmaggies *= 1e9

            log.info('Photometric DN(4000)={:.3f}, Age={:.2f} Gyr, Mr={:.2f} mag'.format(
                dn4000, meanage, absmag[1]))

        # Pack it up and return.
        result['CONTINUUM_COEFF'][0][:nage] = coeff
        result['CONTINUUM_RCHI2'][0] = chi2min
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
        for iband, band in enumerate(self.bands):
            result['FLUX_SYNTH_MODEL_{}'.format(band.upper())] = synth_bestmaggies[iband]

        return result, continuummodel
    
    def continuum_specfit(self, data):
        """Fit the non-negative stellar continuum of a single spectrum.

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            :func:`fastspecfit.io.DESISpectra.read_and_unpack`.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in :func:`self.init_spec_output`.

        Notes
        -----
          - Consider using cross-correlation to update the redrock redshift.

        """
        result = self.init_spec_output()

        redshift = data['zredrock']
        result['CONTINUUM_Z'] = redshift
        for icam, cam in enumerate(data['cameras']):
            result['CONTINUUM_SNR_{}'.format(cam.upper())] = data['snr'][icam]

        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        npixpercamera = [len(gw) for gw in data['wave']]
        npixpercam = np.hstack([0, npixpercamera])
        
        specwave = np.hstack(data['wave'])
        specflux = np.hstack(data['flux'])
        specivar = np.hstack(data['ivar']) * np.logical_not(np.hstack(data['linemask'])) # mask emission lines
        if np.all(specivar == 0) or np.any(specivar < 0):
            errmsg = 'All pixels are masked or some inverse variances are negative!'
            log.critical(errmsg)
            raise ValueError(errmsg)
        
        # Prepare the reddened and unreddened SSP templates by redshifting and
        # normalizing. Note that we ignore templates which are older than the
        # age of the universe at the redshift of the object.
        agekeep = self.younger_than_universe(redshift)
        t0 = time.time()
        zsspflux_dustnomvdisp, _ = self.SSP2data(
            self.sspflux_dustnomvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=data['wave'], specres=data['res'],
            cameras=data['cameras'], synthphot=False)
        zsspflux_dustnomvdisp = np.concatenate(zsspflux_dustnomvdisp, axis=0)  # [npix,nage*nAV]
        npix, nmodel = zsspflux_dustnomvdisp.shape
        nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)
        zsspflux_dustnomvdisp = zsspflux_dustnomvdisp.reshape(npix, nage, self.nAV)       # [npix,nage,nAV]
        log.info('Preparing the models took {:.2f} sec'.format(time.time()-t0))

        ## Fit the spectra with *no* dust reddening so we can identify potential
        ## calibration issues (again, at the nominal velocity dispersion).
        #t0 = time.time()
        #inodust = np.ndarray.item(np.where(self.AV == 0)[0]) # should always be index 0 here
        #coeff, chi2min = self._fnnls_parallel(zsspflux_dustnomvdisp[:, :, inodust],
        #                                      specflux, specivar)
        #log.info('No-dust model fit has chi2={:.3f} and took {:.2f} sec'.format(
        #    chi2min, time.time()-t0))
        #
        #result['CONTINUUM_NODUST_COEFF'][0][0:nage] = coeff
        #result['CONTINUUM_NODUST_CHI2'] = chi2min

        # Fit the spectra for reddening using the models convolved to the
        # nominal velocity dispersion.
        t0 = time.time()
        AVchi2min, AVbest, AVivar = self._fnnls_parallel(zsspflux_dustnomvdisp, specflux, specivar,
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
        if self.solve_vdisp:
            # Here!
            pdb.set_trace()
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

        # Do median-smooth the continuum-subtracted spectrum, to help with the
        # emission-line fitting, particularly for AGN.
        #log.warning('Skipping smoothing residuals!')
        #_smooth_continuum = np.zeros_like(bestfit)

        png = None
        #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/smooth-continuum.png'
        #linemask = np.hstack(data['linemask_all'])
        linemask = np.hstack(data['linemask'])
        if np.all(coeff == 0):
            _smooth_continuum = np.zeros_like(bestfit)
        else:
            _smooth_continuum, _ = self.smooth_continuum(specwave, specflux - bestfit, specivar,
                                                         redshift, linemask=linemask, png=png)

        #_residuals = specflux - bestfit
        #residuals = [_residuals[np.sum(npixpercam[:icam+1]):np.sum(npixpercam[:icam+2])] for icam in np.arange(len(data['cameras']))]
        #_smooth_continuum, _ = self.smooth_continuum(data['wave'], residuals, data['ivar'], png=png,
        #                                             redshift, linemask_dict=data['linemask_dict'])

        #_smooth_continuum = self.smooth_residuals(
        #    residuals, data['wave'], data['ivar'],
        #    data['linemask'], data['linepix'], data['contpix'])
        
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

        # Pack it in and return.
        result['CONTINUUM_COEFF'][0][0:nage] = coeff
        result['CONTINUUM_RCHI2'][0] = chi2min
        result['CONTINUUM_AV'][0] = AVbest
        result['CONTINUUM_AV_IVAR'][0] = AVivar
        result['CONTINUUM_VDISP'][0] = vdispbest
        result['CONTINUUM_VDISP_IVAR'][0] = vdispivar
        result['CONTINUUM_AGE'] = meanage
        result['DN4000'][0] = dn4000
        result['DN4000_IVAR'][0] = dn4000_ivar
        result['DN4000_MODEL'][0] = dn4000_model

        for icam, cam in enumerate(data['cameras']):
            nonzero = continuummodel[icam] != 0
            #nonzero = np.abs(continuummodel[icam]) > 1e-5
            if np.sum(nonzero) > 0:
                corr = np.median(smooth_continuum[icam][nonzero] / continuummodel[icam][nonzero])
                result['CONTINUUM_SMOOTHCORR_{}'.format(cam.upper())] = corr * 100 # [%]

        log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
            result['CONTINUUM_SMOOTHCORR_B'][0], result['CONTINUUM_SMOOTHCORR_R'][0],
            result['CONTINUUM_SMOOTHCORR_Z'][0]))
        
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
                    metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
        elif coadd_type == 'cumulative':
            title = 'Tile/ThruNight: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], coadd_type, metadata['TARGETID']))
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
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
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
        if np.sum(good) > 0:
            ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(continuum_phot_abmag[indx]))) + dm
            ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(continuum_phot_abmag[indx]))) - dm
        else:
            good = phot['abmag'] > 0
            if np.sum(good) > 0:
                ymin = np.nanmax(phot['abmag'][good]) + dm
                ymax = np.nanmin(phot['abmag'][good]) - dm
            else:
                ymin, ymax = [30, 20]
            
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

        if not self.nolegend:
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

        dofit = np.where((abmag > 0) * self.bands_to_fit)[0]
        if len(dofit) > 0:
            ax.errorbar(phot['lambda_eff'][dofit]/1e4, abmag[dofit], lolims=lolims[dofit],
                        yerr=yerr[:, dofit],
                        fmt='s', markersize=11, markeredgewidth=3, markeredgecolor='k',
                        markerfacecolor='red', elinewidth=3, ecolor='red', capsize=4,
                        label=r'$grz\,W_{1}W_{2}$', zorder=2)
            
        ignorefit = np.where((abmag > 0) * (self.bands_to_fit == False))[0]
        if len(ignorefit) > 0:
            ax.errorbar(phot['lambda_eff'][ignorefit]/1e4, abmag[ignorefit], lolims=lolims[ignorefit],
                        yerr=yerr[:, ignorefit],
                        fmt='s', markersize=11, markeredgewidth=3, markeredgecolor='k',
                        markerfacecolor='none', elinewidth=3, ecolor='k', capsize=4)

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
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastphot['CONTINUUM_RCHI2']),
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

        if not self.nolegend:
            legxpos, legypos = 0.04, 0.94
            txt = '\n'.join((
                r'{}'.format(leg['absmag_r']),
                r'{}'.format(leg['absmag_gr'])
                ))
            ax.text(legxpos, legypos, txt, ha='left', va='top',
                    transform=ax.transAxes, fontsize=legfntsz,
                    bbox=bbox)
            
            legxpos, legypos = 0.98, 0.06
            txt = '\n'.join((
                r'{}'.format(leg['zredrock']),
                r'{} {}'.format(leg['chi2'], leg['age']),
                r'{}'.format(leg['AV']),
                ))
            ax.text(legxpos, legypos, txt, ha='right', va='bottom',
                    transform=ax.transAxes, fontsize=legfntsz,
                    bbox=bbox)
    
        plt.subplots_adjust(bottom=0.14, right=0.95, top=0.93)

        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()
