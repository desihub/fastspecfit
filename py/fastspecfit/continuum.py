"""
fastspecfit.continuum
=====================

Methods and tools for continuum-fitting.

"""
import pdb # for debugging

import os, time
import numpy as np
from astropy.table import Table

from fastspecfit.io import FLUXNORM
from fastspecfit.util import TabulatedDESI, C_LIGHT, quantile, median
from fastspecfit.inoue14 import Inoue14

# SPS template constants (used by build-templates)
PIXKMS_BLU = 25.  # [km/s]
PIXKMS_RED = 100. # [km/s]
PIXKMS_WAVESPLIT = 1e4 # [Angstrom]


def convolve_vdisp(templateflux, vdisp, pixsize_kms=PIXKMS_BLU, limit=None):
    """Convolve an input spectrum to a desired velocity dispersion in km/s.

    Parameters
    ----------
    templateflux : :class:`numpy.ndarray` [npix, nmodel]
        One- or two-dimensional input model spectra.
    vdisp : :class:`float`
        Desired velocity dispersion.
    pixsize_kms : :class:`float`
        Pixel size of `templateflux` in km/s.
    limit : :class:`int`
        Only smooth up to the pixel position (in `templateflux`) specified by
        this parameter.

    Returns
    -------
    :class:`numpy.ndarray` [npix, nmodel]
        Gaussian-smoothed model spectra.

    """
    from scipy.ndimage import gaussian_filter1d
    
    # Convolve by the velocity dispersion.
    if vdisp <= 0.:
        output = templateflux.copy()
    else:
        output = np.empty_like(templateflux)
        sigma = vdisp / pixsize_kms # [pixels]

        if limit is None:
            limit = templateflux.shape[0]

        if templateflux.ndim == 1:
            gaussian_filter1d(templateflux[:limit], sigma=sigma, output=output[:limit])
            output[limit:] = templateflux[limit:]
        else:
            if limit is None:
                limit = templateflux.shape[0]
            gaussian_filter1d(templateflux[:limit, :], sigma=sigma, axis=0,
                              output=output[:limit:, :])
            output[limit:, :] = templateflux[limit:, :]

    return output


def _smooth_continuum(wave, flux, ivar, linemask, redshift=0., png=None, log=None):
    """Build a smooth, nonparametric continuum spectrum.

    Parameters
    ----------
    wave : :class:`numpy.ndarray` [npix]
        Observed-frame wavelength array.
    flux : :class:`numpy.ndarray` [npix]
        Spectrum corresponding to `wave`.
    ivar : :class:`numpy.ndarray` [npix]
        Inverse variance spectrum corresponding to `flux`.
    linemask : :class:`numpy.ndarray` of type :class:`bool`
        Boolean mask with the same number of pixels as `wave` where `True`
        means a pixel is affected by an emission line.
    camerapix : :class:`list` [ncamera]
        List of starting and ending :class:`int` pixel indices indicating
         the range of each camera.
    png : :class:`str`
        Generate a simple QA plot and write it out to this filename.

    Returns
    -------
    smooth :class:`numpy.ndarray` [npix]
        Smooth continuum spectrum which can be subtracted from `flux` in
        order to create a pure emission-line spectrum.

    """
    from scipy.interpolate import splev, splrep
        
    if log is None:
        from desiutil.log import get_logger
        log = get_logger()

    if len(linemask) != len(wave):
        errmsg = 'Linemask must have the same number of pixels as the input spectrum.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    mask = linemask
    mask[wave < (1215. * (1. + redshift))] = True

    srt = np.argsort(wave)

    #import numpy.ma as ma
    #mawave = ma.masked_array(wave)
    #mawave[mask] = ma.masked
    #clumps = ma.clump_unmasked(mawave)
    #edges = ma.notmasked_edges(mawave)

    good = np.logical_not(linemask[srt])
    ww = np.sqrt(good * ivar[srt])

    from scipy.interpolate import UnivariateSpline
    spl = UnivariateSpline(wave[srt], flux[srt], w=ww, k=3, ext=1)

    knots = spl.get_knots()
    smooth = spl(wave[srt])

    #linemask[wave < (1215. * (1. + 2.6916212))] = True
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(wave, flux, alpha=0.75)
    plt.plot(wave, smooth, color='blue')
    plt.scatter(knots, np.interp(knots, wave[srt], flux[srt]), marker='s', s=75, color='red')
    plt.savefig('junk.png')

    plt.scatter(wave[mask], flux[mask], color='k', marker='s', zorder=2)

    #plt.xlim(4000, 5000) ; 
    
    pdb.set_trace()


    mm = np.sum(good)
    ss = mm - np.sqrt(2 * mm)

    pdb.set_trace()

    spl = splrep(wave[srt], flux[srt], w=w, k=3, s=ss)

    if np.any(np.isnan(spl[1])):
        errmsg = 'Spline-fitting has failed!'
        log.critical(errmsg)
        raise ValueError(errmsg)

    smooth = splev(wave, spl)

    ## Treat each camera independently. Assume wavelengths are sorted!
    #if camerapix is not None:
    #    smooth = []
    #    for pix in camerapix:
    #        spl = splrep(wave[pix[0]:pix[1]], flux[pix[0]:pix[1]],
    #                     w=np.sqrt(ivar[pix[0]:pix[1]]), k=3)
    #        smooth.append(splev(wave[pix[0]:pix[1]], spl))
    #    smooth = np.hstack(smooth)
    #else:
    #    spl = splrep(wave, flux, w=np.sqrt(ivar), k=3)
    #    smooth = splev(wave, spl)


    # Optional QA.
    if png:
        import matplotlib.pyplot as plt
        import seaborn as sns

        resid = flux - smooth        
        noise = np.ptp(quantile(resid[~linemask], (0.25, 0.75))) / 1.349 # robust sigma

        sns.set(context='talk', style='ticks', font_scale=0.8)

        fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        ax[0].plot(wave / 1e4, flux, alpha=0.75, label='Data')
        ax[0].scatter(wave[linemask] / 1e4, flux[linemask], s=10, marker='s', 
                      color='k', zorder=2, alpha=0.5, label='Line-masked pixel')
        ax[0].plot(wave / 1e4, smooth, color='red', label='Smooth continuum')
        ax[0].set_ylim(np.min((-5. * noise, quantile(flux, 0.05))), 
                       np.max((5. * noise, 1.5 * quantile(flux, 0.975))))
        ax[0].set_ylabel('Continuum-subtracted Spectrum')
        ax[0].legend(fontsize=12)

        ax[1].plot(wave / 1e4, resid, alpha=0.75)#, label='Residuals')
        ax[1].axhline(y=0, color='k')
        ax[1].set_ylim(np.min((-5. * noise, quantile(resid, 0.05))), 
                       np.max((5. * noise, 1.5 * quantile(resid, 0.975))))
        ax[1].set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
        ax[1].set_ylabel('Residual Spectrum')
        #ax[1].legend(fontsize=10)

        fig.tight_layout()
        fig.savefig(png)#, bbox_inches='tight')
        plt.close()

    return smooth
    

class Tools(TabulatedDESI):
    """Class to load filters and containing filter- and dust-related methods.

    """
    def __init__(self, ignore_photometry=False, fphoto=None, load_filters=True, 
                 log=None, verbose=False):
        """
        Parameters
        ----------
        ignore_photometry : :class:`bool`
            Boolean flag indicating whether or not to ignore the broadband
            photometry.
        fphoto : :class:`dict`
            Photometry dictionary containing information on bandpasses, filter
            response curves, etc.
        load_filters : :class:`bool`
            Load the filters by reading them from disk.
        verbose : :class:`bool`
            Include more verbose (debugging) output.
        log : :class:`desiutil.log.logger`
            Logger object.

        """
        super(Tools, self).__init__()
        
        from speclite import filters

        if log is None:
            from desiutil.log import get_logger, DEBUG
            if verbose:
                log = get_logger(DEBUG)
            else:
                log = get_logger()
        self.log = log
        
        if fphoto is not None:
            keys = fphoto.keys()
            self.uniqueid = fphoto['uniqueid']
            self.photounits = fphoto['photounits']
            if 'readcols' in keys:
                self.readcols = np.array(fphoto['readcols'])
            if 'dropcols' in keys:
                self.dropcols = np.array(fphoto['dropcols'])
            if 'outcols' in keys:
                self.outcols = np.array(fphoto['outcols'])
            self.bands = np.array(fphoto['bands'])
            self.bands_to_fit = np.array(fphoto['bands_to_fit'])
            self.fluxcols = np.array(fphoto['fluxcols'])
            self.fluxivarcols = np.array(fphoto['fluxivarcols'])
            self.min_uncertainty = np.array(fphoto['min_uncertainty'])
            
            if 'legacysurveydr' in keys:
                self.legacysurveydr = fphoto['legacysurveydr']
            if 'viewer_layer' in keys:
                self.viewer_layer = fphoto['viewer_layer']
            if 'viewer_pixscale' in keys:
                self.viewer_pixscale = fphoto['viewer_pixscale']
            if 'synth_bands' in keys:
                self.synth_bands = np.array(fphoto['synth_bands'])
            if 'fiber_bands' in keys:
                self.fiber_bands = np.array(fphoto['fiber_bands'])

            self.absmag_bands = np.array(fphoto['absmag_bands'])
            self.band_shift = np.array(fphoto['band_shift'])

            if load_filters:
                # If fphoto['filters'] is a dictionary, then assume that there
                # are N/S filters (as indicated by photsys).
                self.filters = {}
                for key in fphoto['filters'].keys():
                    self.filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                                for filtname in fphoto['filters'][key]])
                self.synth_filters = {}
                for key in fphoto['synth_filters'].keys():
                    self.synth_filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                                      for filtname in fphoto['synth_filters'][key]])
                if hasattr(self, 'fiber_bands'):
                    self.fiber_filters = {}
                    for key in fphoto['fiber_filters'].keys():
                        self.fiber_filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                                          for filtname in fphoto['fiber_filters'][key]])
                # Simple list of filters.
                self.absmag_filters = filters.FilterSequence([filters.load_filter(filtname) for filtname in fphoto['absmag_filters']])

        else:
            # This should never happen in production because we read the default
            # fphoto file in io.DESISpectra.__init__().
            self.uniqueid = 'TARGETID'
            self.photounits = 'nanomaggies'
            self.readcols = np.array(['TARGETID', 'RA', 'DEC', 'RELEASE', 'LS_ID',
                                      'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',
                                      'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z'])
            self.outcols = np.array(['PHOTSYS', 'LS_ID'])
            self.bands = np.array(['g', 'r', 'z', 'W1', 'W2', 'W3', 'W4'])
            self.bands_to_fit = np.ones(len(self.bands), bool)
            self.fluxcols = np.array(['FLUX_G', 'FLUX_R', 'FLUX_Z',
                                      'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4'])
            self.fluxivarcols = np.array(['FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                                          'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FLUX_IVAR_W3', 'FLUX_IVAR_W4'])
            self.min_uncertainty = np.array([0.02, 0.02, 0.02, 0.05, 0.05, 0.05, 0.05]) # mag
            
            self.legacysurveydr = 'dr9'
            self.viewer_layer = 'ls-dr9'
            self.viewer_pixscale = 0.262
            self.synth_bands = np.array(['g', 'r', 'z']) # for synthesized photometry
            self.fiber_bands = np.array(['g', 'r', 'z']) # for fiber fluxes

            self.absmag_bands = ['decam_g', 'decam_r', 'decam_z',
                                 'U', 'B', 'V',
                                 'sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z',
                                 'W1']#, 'W2']
            self.band_shift = [1.0, 1.0, 1.0,
                               0.0, 0.0, 0.0,
                               0.1, 0.1, 0.1, 0.1, 0.1,
                               0.1]#, 0.1]

            if load_filters:
                self.filters = {'N': filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z', 
                                                          'wise2010-W1', 'wise2010-W2', 'wise2010-W3', 'wise2010-W4'),
                                'S': filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z', 
                                                          'wise2010-W1', 'wise2010-W2', 'wise2010-W3', 'wise2010-W4')}
                self.synth_filters = {'N': filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z'),
                                      'S': filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z')}
                self.fiber_filters = self.synth_filters

                self.absmag_filters = filters.FilterSequence((
                    filters.load_filter('decam2014-g'), filters.load_filter('decam2014-r'), filters.load_filter('decam2014-z'), 
                    filters.load_filter('bessell-U'), filters.load_filter('bessell-B'), filters.load_filter('bessell-V'), 
                    filters.load_filter('sdss2010-u'), filters.load_filter('sdss2010-g'), filters.load_filter('sdss2010-r'),
                    filters.load_filter('sdss2010-i'), filters.load_filter('sdss2010-z'),
                    filters.load_filter('wise2010-W1')))#, filters.load_filter('wise2010-W2')))

        if len(self.absmag_bands) != len(self.band_shift):
            errmsg = 'absmag_bands and band_shift must have the same number of elements.'
            self.log.critical(errmsg)
            raise ValueError(errmsg)
        
        if self.photounits != 'nanomaggies':
            errmsg = 'nanomaggies is the only currently supported photometric unit!'
            self.log.critical(errmsg)
            raise ValueError(errmsg)

        # Do not fit the photometry.
        if ignore_photometry:
            self.bands_to_fit *= [False]


    def restframe_photometry(self, redshift, zmodelflux, zmodelwave, maggies, ivarmaggies,
                             filters_in, absmag_filters, band_shift=None, snrmin=2.,
                             dmod=None, cosmo=None):
        """Compute K-corrections and rest-frame photometry for a single object.
        
        Parameters
        ----------
        redshift : :class:`float`
           Galaxy or QSO redshift.
        zmodelwave : `numpy.ndarray`
           Observed-frame (redshifted) model wavelength array.
        zmodelflux : `numpy.ndarray`
           Observed-frame model spectrum.
        maggies : `numpy.ndarray`
           Input photometric fluxes in the `filters_in` bandpasses.
        ivarmaggies : `numpy.ndarray`
           Inverse variance photometry corresponding to `maggies`.
        filters_in : `speclite.filters.FilterSequence`
           Input filter curves.
        absmag_filters : `speclite.filters.FilterSequence`
           Filter curves corresponding to desired bandpasses.
        band_shift : `numpy.ndarray` or `None`
           Band-shift each bandpass in `absmag_filters` by this amount.
        snrmin : :class:`float`, defaults to 2.
           Minimum signal-to-noise ratio in the input photometry (`maggies`) in
           order for that bandpass to be used to compute a K-correction.
        dmod : :class:`float` or `None`
           Distance modulus corresponding to `redshift`. Not needed if `cosmo` is
           provided.
        cosmo : `fastspecfit.util.TabulatedDESI` or `None`
           Cosmological model class needed to compute the distance modulus.
        log : `desiutil.log`
           Logging object.
    
        Returns
        -------
        kcorr : `numpy.ndarray`
           K-corrections for each bandpass in `absmag_filters`.
        absmag : `numpy.ndarray`
           Absolute magnitudes, band-shifted according to `band_shift` (if
           provided) for each bandpass in `absmag_filters`. 
        ivarabsmag : `numpy.ndarray`
           Inverse variance corresponding to `absmag`.
        synth_absmag : `numpy.ndarray`
           Like `absmag`, but entirely based on synthesized photometry.
        synth_maggies_in : `numpy.ndarray`
           Synthesized input photometry (should closely match `maggies` if the
           model fit is good).
        
        Notes
        -----
        By default, the K-correction is computed by finding the observed-frame
        bandpass closest in wavelength (and with a minimum signal-to-noise ratio) to
        the desired band-shifted absolute magnitude bandpass. In other words, by
        default we endeavor to minimize the K-correction. The inverse variance,
        `ivarabsmag`, is derived from the inverse variance of the K-corrected
        photometry. If no bandpass is available then `ivarabsmag` is set to zero and
        `absmag` is derived from the synthesized rest-frame photometry.
        
        """
        from speclite import filters
        
        nabs = len(absmag_filters)
            
        if redshift <= 0.0:
            errmsg = 'Input redshift not defined, zero, or negative!'
            self.log.warning(errmsg)
            kcorr        = np.zeros(nabs, dtype='f8')
            absmag       = np.zeros(nabs, dtype='f8')
            ivarabsmag   = np.zeros(nabs, dtype='f8')
            synth_absmag = np.zeros(nabs, dtype='f8')
            synth_maggies_in = np.zeros(len(maggies))
            return kcorr, absmag, ivarabsmag, synth_absmag, synth_maggies_in

        if cosmo is None:
            from fastspecfit.util import TabulatedDESI
            cosmo = TabulatedDESI()

        if dmod is None:
            dmod = cosmo.distance_modulus(redshift)
    
        modelwave = zmodelwave / (1. + redshift)
        lambda_in = filters_in.effective_wavelengths.value

        if band_shift is None:
            band_shift = np.zeros_like(lambda_in)

        # input bandpasses, observed frame; maggies and synth_maggies_in should be
        # very close.
        synth_maggies_in = self.get_ab_maggies(filters_in,
                                               zmodelflux / FLUXNORM,
                                               zmodelwave,
                                               log=self.log)
        filters_out = \
            filters.FilterSequence( [ f.create_shifted(band_shift=bs) for f, bs in zip(absmag_filters, band_shift) ])
        lambda_out = filters_out.effective_wavelengths.value
        
        # Multiply by (1+z) to convert the best-fitting model to the "rest
        # frame".
        synth_outmaggies_rest = self.get_ab_maggies(filters_out,
                                                    zmodelflux * (1. + redshift) / FLUXNORM,
                                                    modelwave,
                                                    log=self.log)
        
        synth_absmag = -2.5 * np.log10(synth_outmaggies_rest) - dmod

        # K-correct from the nearest "good" bandpass (to minimizes the K-correction)
        oband = np.empty(nabs, dtype=np.int32)
        for jj in range(nabs):
            lambdadist = np.abs(lambda_in / (1. + redshift) - lambda_out[jj])
            oband[jj] = np.argmin(lambdadist + (maggies * np.sqrt(ivarmaggies) < snrmin) * 1e10)

        kcorr = + 2.5 * np.log10(synth_outmaggies_rest / synth_maggies_in[oband])
        
        # m_R = M_Q + DM(z) + K_QR(z) or
        # M_Q = m_R - DM(z) - K_QR(z)
        absmag = -2.5 * np.log10(maggies[oband]) - dmod - kcorr

        C = 0.8483036976765437 # (0.4 * np.log(10.))**2
        ivarabsmag = maggies[oband]**2 * ivarmaggies[oband] * C
    
        # if we use synthesized photometry then ivarabsmag is zero
        # (which should never happen?)
        I = (maggies[oband] * np.sqrt(ivarmaggies[oband]) <= snrmin)
        absmag[I] = synth_absmag[I]
        ivarabsmag[I] = 0.
    
        return kcorr, absmag, ivarabsmag, synth_absmag, synth_maggies_in


    @staticmethod
    def get_ab_maggies_fast(filters, flux, wave, log=None):
        """Like `self.get_ab_maggies()`, but by-passing the units and error-checking
        used by `speclite`. Specifically, we assume that the filter response
        function always lies within the wavelength range of the input spectrum.

        flux and wave are assumed to be in erg/s/cm2/A and A, respectively. 

        """
        # AB reference spctrum in erg/s/cm2/Hz times the speed of light in A/s
        # and converted to erg/s/cm2/A.
        abflam = 3.631e-20 * C_LIGHT * 1e13 / wave**2

        maggies = np.zeros(len(filters))
        for ifilt, filt in enumerate(filters):
            lo = np.searchsorted(wave, filt.wavelength[0]-2., 'right')
            hi = np.searchsorted(wave, filt.wavelength[-1]+2., 'left')
            resp = np.interp(wave[lo:hi], filt.wavelength, filt.response, left=0., right=0.)
            numer = np.trapz(resp * flux[lo:hi] * wave[lo:hi], x=wave[lo:hi])
            denom = np.trapz(resp * abflam[lo:hi] * wave[lo:hi], x=wave[lo:hi])
            maggies[ifilt] = numer / denom

        return maggies


    @staticmethod
    def get_ab_maggies(filters, flux, wave, log=None):
        
        try:
            maggies0 = filters.get_ab_maggies(flux, wave)
        except:
            # pad in case of an object at very high redshift (z > 5.5)
            if log is not None:
                log.warning('Padding model spectrum due to insufficient wavelength coverage to synthesize photometry.') 
            padflux, padwave = filters.pad_spectrum(flux, wave, axis=0, method='edge')
            maggies0 = filters.get_ab_maggies(padflux, padwave)
        
        if len(maggies0) == 1:
            maggies = np.fromiter(maggies0.values(), np.float64)
        else:
            maggies = np.empty((len(maggies0.colnames), len(maggies0)), dtype=np.float64)
            for i, col in enumerate(maggies0.values()):
                maggies[i, :] = col.value
        
        return maggies

    
    @staticmethod
    def to_nanomaggies(maggies):
        return maggies * 1e9

    
    @staticmethod
    def get_photflam(maggies, lambda_eff, nanomaggies=True):
        
        shp = maggies.shape
        if maggies.ndim == 1:
            ngal = 1
        else:
            ngal = shp[1]
        
        if nanomaggies:
            nanofactor = 1e-9 # [nanomaggies-->maggies]
        else:
            nanofactor = 1.0

        factor = nanofactor * 10**(-0.4 * 48.6) * C_LIGHT * 1e13 / lambda_eff**2 # [maggies-->erg/s/cm2/A]
        if ngal > 1:
            factor = factor[:, None] # broadcast for the models
        return maggies * factor

    
    @staticmethod
    def parse_photometry(bands, maggies, lambda_eff, ivarmaggies=None,
                         nanomaggies=True, nsigma=2.0, min_uncertainty=None,
                         log=None, get_abmag=False):
        """Parse input (nano)maggies to various outputs and pack into a table.

        Parameters
        ----------
        flam - 10-17 erg/s/cm2/A
        fnu - 10-17 erg/s/cm2/Hz
        abmag - AB mag
        nanomaggies - input maggies are actually 1e-9 maggies

        nsigma - magnitude limit 

        qa - true iff table will be used by fastqa (which needs
        columns that fastspec does not)
        
        Returns
        -------
        phot - photometric table

        Notes
        -----

        """
        if log is None:
            from desiutil.log import get_logger
            log = get_logger()
        
        if ivarmaggies is None:
            ivarmaggies = np.zeros_like(maggies)
            
        # Gaia-only targets can sometimes have grz=-99.
        if np.any(ivarmaggies < 0.) or np.any(maggies == -99.0):
            errmsg = 'All ivarmaggies must be zero or positive!'
            log.critical(errmsg)
            raise ValueError(errmsg)

        if nanomaggies:
            nanofactor = 1e-9 # [nanomaggies-->maggies]
            nmg = maggies
            nmg_ivar = ivarmaggies.copy()
        else:
            nanofactor = 1.0
            nmg = maggies * 1e9
            nmg_ivar = ivarmaggies * 1e-18

        if get_abmag:
            # compute columns used only by fastqa
            abmag           = np.zeros_like(maggies)
            abmag_limit     = np.zeros_like(maggies)
            abmag_brighterr = np.zeros_like(maggies)
            abmag_fainterr  = np.zeros_like(maggies)
            abmag_ivar      = np.zeros_like(maggies)
            
            # deal with measurements
            good = (maggies > 0.)
            abmag[good] = -2.5 * np.log10(nanofactor * maggies[good])
            
            # deal with upper limits
            snr = maggies * np.sqrt(ivarmaggies)
            upper = ((ivarmaggies > 0.) & (snr <= nsigma))
            abmag_limit[upper] = - 2.5 * np.log10(nanofactor * nsigma / np.sqrt(ivarmaggies[upper]))
            
            # significant detections
            C = 0.4 * np.log(10)
            good = (snr > nsigma)
            maggies_good = maggies[good]
            ivarmaggies_good = ivarmaggies[good]
            errmaggies = 1. / np.sqrt(ivarmaggies_good)
            abmag_brighterr[good] = errmaggies / (C * (maggies_good + errmaggies)) # bright end (flux upper limit)
            abmag_fainterr[good]  = errmaggies / (C * (maggies_good - errmaggies)) # faint end (flux lower limit)
            abmag_ivar[good]      = ivarmaggies_good * (C * maggies_good)**2
            
        # Add a minimum uncertainty in quadrature **but only for flam**, which
        # is used in the fitting.
        if min_uncertainty is not None:
            log.info('Propagating minimum photometric uncertainties (mag): [{}]'.format(
                ' '.join(min_uncertainty.astype(str))))
            good = ((maggies != 0.) & (ivarmaggies > 0.))
            maggies_good = maggies[good]
            factor = 2.5 / np.log(10.)
            magerr = factor / (np.sqrt(ivarmaggies[good]) * maggies_good)
            magerr2 = magerr**2 + min_uncertainty[good]**2
            ivarmaggies[good] = factor**2 / (maggies_good**2 * magerr2)

        factor = nanofactor * 10**(-0.4 * 48.6) * C_LIGHT * 1e13 / lambda_eff**2 # [maggies-->erg/s/cm2/A]
        
        ngal = 1 if maggies.ndim == 1 else maggies.shape[1]
        if ngal > 1:
            factor = factor[:, None] # broadcast for the models
        flam = maggies * factor
        flam_ivar = ivarmaggies / factor**2

        data = {
            'band':             bands,
            'lambda_eff':       lambda_eff, 
            'nanomaggies':      nmg,
            'nanomaggies_ivar': nmg_ivar,
            'flam':             flam,
            'flam_ivar':        flam_ivar,
        }
        dtypes = [
            bands.dtype,
            'f4',
            'f4',
            'f4',
            'f8', # flam
            'f8', # flam_ivar
        ]
        
        if get_abmag:
            # add columns used only by fastqa
            data_qa = {
                'abmag':            abmag,
                'abmag_ivar':       abmag_ivar,
                'abmag_brighterr':  abmag_brighterr,
                'abmag_fainterr':   abmag_fainterr,
                'abmag_limit':      abmag_limit,
            }
            data |= data_qa

            dtypes_qa = [
                'f4',
                'f4',
                'f4',
                'f4',
                'f4',
            ]
            dtypes.extend(dtypes_qa)            
            
        phot = Table(data=data, dtype=dtypes)
        
        return phot

    
    @staticmethod
    def get_dn4000(wave, flam, flam_ivar=None, redshift=None, rest=True, log=None):
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

        Require full wavelength coverage over the definition of the index.

        See eq. 11 in Bruzual 1983
        (https://articles.adsabs.harvard.edu/pdf/1983ApJ...273..105B) but with
        the "narrow" definition of Balogh et al. 1999.

        """
        from fastspecfit.util import ivar2var

        if log is None:
            from desiutil.log import get_logger
            log = get_logger()
    
        dn4000, dn4000_ivar = 0.0, 0.0

        if rest is False or redshift is not None:
            restwave = wave / (1. + redshift) # [Angstrom]
            flam2fnu = (1. + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
        else:
            restwave = wave
            flam2fnu =  restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]

        # Require a 2-Angstrom pad around the break definition.
        wpad = 2.0
        if np.min(restwave) > (3850-wpad) or np.max(restwave) < (4100+wpad):
            log.warning('Too little wavelength coverage to compute Dn(4000)')
            return dn4000, dn4000_ivar

        fnu = flam * flam2fnu # [erg/s/cm2/Hz]

        if flam_ivar is not None:
            fnu_ivar = flam_ivar / flam2fnu**2
        else:
            fnu_ivar = np.ones_like(flam) # uniform weights

        def _integrate(wave, flux, ivar, w1, w2):
            from scipy import integrate, interpolate
            # trim for speed
            I = (wave > (w1-wpad)) * (wave < (w2+wpad))
            J = np.logical_and(I, ivar > 0)
            # Require no more than 20% of pixels are masked.
            if np.sum(J) / np.sum(I) < 0.8:
                log.warning('More than 20% of pixels in Dn(4000) definition are masked.')
                return 0.0
            wave = wave[J]
            flux = flux[J]
            ivar = ivar[J]
            # should never have to extrapolate
            f = interpolate.interp1d(wave, flux, assume_sorted=False, bounds_error=True)
            f1 = f(w1)
            f2 = f(w2)
            i = interpolate.interp1d(wave, ivar, assume_sorted=False, bounds_error=True)
            i1 = i(w1)
            i2 = i(w2)
            # insert the boundary wavelengths then integrate
            I = np.where((wave > w1) * (wave < w2))[0]
            wave = np.insert(wave[I], [0, len(I)], [w1, w2])
            flux = np.insert(flux[I], [0, len(I)], [f1, f2])
            ivar = np.insert(ivar[I], [0, len(I)], [i1, i2])
            weight = integrate.simps(x=wave, y=ivar)
            index = integrate.simps(x=wave, y=flux*ivar) / weight
            index_var = 1 / weight
            return index, index_var

        blufactor = 3950.0 - 3850.0
        redfactor = 4100.0 - 4000.0
        try:
            # yes, blue wavelength go with red integral bounds
            numer, numer_var = _integrate(restwave, fnu, fnu_ivar, 4000, 4100)
            denom, denom_var = _integrate(restwave, fnu, fnu_ivar, 3850, 3950)
        except:
            log.warning('Integration failed when computing DN(4000).')
            return dn4000, dn4000_ivar

        if denom == 0.0 or numer == 0.0:
            log.warning('DN(4000) is ill-defined or could not be computed.')
            return dn4000, dn4000_ivar
        
        dn4000 =  (blufactor / redfactor) * numer / denom
        if flam_ivar is not None:
            dn4000_ivar = (1.0 / (dn4000**2)) / (denom_var / (denom**2) + numer_var / (numer**2))
    
        return dn4000, dn4000_ivar

    
    @staticmethod
    def _linepix_and_contpix(wave, ivar, linetable, linesigmas, linevshifts=None, 
                             patchMap=None, redshift=0., nsigma=4., minlinesigma=50., 
                             mincontpix=11, get_contpix=True, log=None):
        """Support routine to determine the pixels potentially containing emission lines
        and the corresponding (adjacent) continuum.

        linesigmas in km/s
        minlinesigma in kms - minimum line-sigma for the purposes of masking
        mincontpix - minimum number of continuum pixels per line

        """
        def _get_linepix(zlinewave, sigma):
            # line-emission
            I = (wave > (zlinewave - nsigma * sigma)) * (wave < (zlinewave + nsigma * sigma)) * (ivar > 0.)
            return np.where(I)[0]

        def _get_contpix(zlinewaves, sigmas, nsigma_factor=2., linemask=None, lya=False):
            # never use continuum pixels blueward of Lyman-alpha
            if lya:
                minwave = np.min(zlinewaves)
            else:
                minwave = np.min(zlinewaves - nsigma_factor * nsigma * sigmas)
            maxwave = np.max(zlinewaves + nsigma_factor * nsigma * sigmas)
            if linemask is None:
                J = (wave > minwave) * (wave < maxwave) * (ivar > 0.)
            else:
                J = (wave > minwave) * (wave < maxwave) * (ivar > 0.) * ~linemask
            return np.where(J)[0]


        if log is None:
            from desiutil.log import get_logger
            log = get_logger()

        if linevshifts is None:
            linevshifts = np.zeros_like(linesigmas)

        if len(linevshifts) != len(linesigmas):
            errmsg = 'Mismatch in linevshifts and linesigmas dimensions.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        linemask = np.zeros_like(wave, bool) # True - affected by possible emission line
        linenames = linetable['name'].value

        zlinewaves = linetable['restwave'].value * (1. + redshift + linevshifts / C_LIGHT)
        linesigmas[(linesigmas > 0.) * (linesigmas < minlinesigma)] = minlinesigma
        linesigmas_ang = linesigmas * zlinewaves / C_LIGHT # [km/s --> Angstrom]

        pix = {'linepix': {}, 'contpix': {}}
        if patchMap is not None:
            pix.update({'patch_contpix': {}, 'dropped': [], 'merged': []})
            patchids = list(patchMap.keys())
            npatch = len(patchids)

        # Initial set of pixels that may contain emission lines.
        for linename, zlinewave, sigma in zip(linenames, zlinewaves, linesigmas_ang):
            # skip fixed (e.g., hbeta_broad) lines
            if sigma <= 0.:
                continue 
            I = _get_linepix(zlinewave, sigma)
            if len(I) > 0:
                linemask[I] = True
                pix['linepix'][linename] = I

        # skip contpix
        if not get_contpix:
            return pix

        if patchMap is None:
            for linename, zlinewave, sigma in zip(linenames, zlinewaves, linesigmas_ang):
                # skip fixed (e.g., hbeta_broad) or fully masked lines
                if sigma <= 0. or not linename in pix['linepix'].keys():
                    continue
                lya = linename == 'lyalpha'
                J = _get_contpix(zlinewave, sigma, nsigma_factor=2., linemask=linemask, lya=lya)
                # go further out
                if len(J) < mincontpix:
                    J = _get_contpix(zlinewave, sigma, nsigma_factor=2.5, linemask=linemask, lya=lya)
                # drop the linemask_ condition; dangerous??
                if len(J) < mincontpix: 
                    #log.debug(f'Dropping linemask condition for {linename} with width ' + \
                    #          f'{nsigma*sigma:.3f}-{3*nsigma*sigma:.3f} Angstrom')

                    # Note the smaller nsigma_factor (to stay closer to the
                    # line); remove the pixels already assigned to this
                    # line.
                    J = _get_contpix(zlinewave, sigma, nsigma_factor=2., linemask=None, lya=lya)
                    J = np.delete(J, np.isin(J, pix['linepix'][linename]))

                if len(J) > 0:
                    pix['contpix'][linename] = J
        else:
            # Now get the corresponding continuum pixels in "patches."
            for patchid in patchids:
                # Not all patchlines will be in the 'linepix' dictionary
                # because, e.g., the broad Balmer lines have linesigma=0 when
                # fitting the narrow-only linemodel. An entire patch can also
                # get dropped if a large portion of the spectrum is fully masked
                # (ivar==0).
                patchlines = patchMap[patchid][0]
                keep = np.isin(patchlines, list(pix['linepix'].keys()))
                if np.count_nonzero(keep) == 0:
                    pix['dropped'].append(patchid)
                    log.info(f'Dropping patch {patchid} ({len(patchlines)} lines fully masked).')
                    continue

                patchlines = patchlines[keep]
                I = patchMap[patchid][1][keep]

                zlinewaves_patch = zlinewaves[I]
                sigmas_patch = linesigmas_ang[I]
                lya = 'lyalpha' in patchlines

                J = _get_contpix(zlinewaves_patch, sigmas_patch, nsigma_factor=2.,
                                 linemask=linemask, lya=lya)

                # go further out
                if len(J) < mincontpix: 
                    J = _get_contpix(zlinewaves_patch, sigmas_patch, nsigma_factor=2.5,
                                     linemask=linemask, lya=lya)

                if len(J) > 0:
                    # all lines in this patch get the same continuum indices
                    mn, mx = np.min(J), np.max(J)
                    for patchline in patchlines:
                        pix['contpix'][patchline] = J
                        
                        # Make sure the left/right edges of the patch include all
                        # the emission lines on this patch.
                        mn = np.min((mn, np.min(pix['linepix'][patchline])))
                        mx = np.max((mx, np.max(pix['linepix'][patchline])))
                    J = np.unique(np.hstack((mn, J, mx)))

                    pix['patch_contpix'][patchid] = J # updated vector

            # Loop back through and merge patches that overlap by at least
            # mincontpix.
            patchkeys = pix['patch_contpix'].keys()
            for ipatch, patchid in enumerate(patchids[:npatch-1]):
                if patchid in patchkeys and patchids[ipatch+1] in patchkeys:
                    left = pix['patch_contpix'][patchid]
                    rght = pix['patch_contpix'][patchids[ipatch+1]]
                    incommon = np.intersect1d(left, rght, assume_unique=True)
                    if len(incommon) > mincontpix:
                        log.debug(f'Merging patches {patchid} and {patchids[ipatch+1]}')
                        newcontpix = np.union1d(left, rght)
                        newpatchid = patchids[ipatch] + patchids[ipatch+1]
                        del pix['patch_contpix'][patchids[ipatch]]
                        del pix['patch_contpix'][patchids[ipatch+1]]
                        pix['patch_contpix'][newpatchid] = newcontpix
                        pix['merged'].append(newpatchid)

        if patchMap is not None:
            if len(pix['dropped']) > 0:
                pix['dropped'] = np.hstack(pix['dropped'])
            if len(pix['merged']) > 0:
                pix['merged'] = np.hstack(pix['merged'])

        # make sure we haven't mucked up our indexing.
        if not pix['linepix'].keys() == pix['contpix'].keys():
            errmsg = 'Mismatch in linepix and contpix!'
            log.critical(errmsg)
            raise ValueError(errmsg)

        return pix


class ContinuumTools(Tools):
    """Tools for dealing with spectral continua.

    """
    def __init__(self, emlinesfile=None, **kwargs):
        """
        Parameters
        ----------
        emlinesfile : :class:`str`
            Full path to the emission lines parameter file to read.
        kwargs : :class:`dict`
            Optional keyword arguments passed to `Tools` :class:`class`.
    
        """
        from fastspecfit.io import read_emlines

        super(ContinuumTools, self).__init__(**kwargs)

        self.emlinesfile = emlinesfile
        self.linetable = read_emlines(emlinesfile=emlinesfile)
        self.igm = Inoue14()

        self.massnorm = 1e10       # stellar mass normalization factor [Msun]
        self.dust_power = -0.7     # power-law slope
        self.dust_normwave = 5500. # pivot wavelength

        # The dust attenuation curve will be cached once we have the template
        # wavelength vector.
        self.dust_klambda = None   

        # The redshift-dependent factors will be cached once we have the
        # redshift and the redshifted template wavelength vector.
        self.zfactors = None       

        # The pixel position in templatewave corresponding to PIXKMS_WAVESPLIT
        # will be cached once we have the wavelength vector.
        self.pixpos_wavesplit = None 

        
    @staticmethod
    def smooth_continuum(*args, **kwargs):
        return _smooth_continuum(*args, **kwargs)

    
    def build_linemask_patches(self, wave, flux, ivar, resolution_matrix, redshift=0.0, 
                               uniqueid=0, initsigma_broad=3000., initsigma_narrow=150., 
                               initsigma_balmer_broad=1000., initvshift_broad=0., 
                               initvshift_narrow=0., initvshift_balmer_broad=0.,
                               minsnr_balmer_broad=1.5, niter=2, debug_plots=False):
        """Generate a mask which identifies pixels impacted by emission lines.

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
            Observed-frame wavelength array.
        flux : :class:`numpy.ndarray` [npix]
            Spectrum corresponding to `wave`.
        ivar : :class:`numpy.ndarray` [npix]
            Inverse variance spectrum corresponding to `flux`.
        resolution_matrix : :class:`list`
            List of :class:`fastspecfit.emline_fit.sparse_rep.ResMatrix`
            objects, one per camera.
        redshift : :class:`float`
            Object redshift.
        uniqueid : :class:`int`
            Unique identification number.
        initsigma_broad : :class:`float`
            Initial guess of the broad-line emission-line width in km/s.
        initsigma_narrow : :class:`float`
            Like `initsigma_broad` but for the forbidden lines.
        initsigma_balmer_broad : :class:`float`
            Like `initsigma_broad` but for the broad Balmer lines.
        initvshift_broad : :class:`float`
            Initial guess of the broad-line velocity shift in km/s.
        initvshift_narrow : :class:`float`
            Like `initvshift_broad` but for the forbidden lines.
        initvshift_balmer_broad : :class:`float`
            Like `initvshift_broad` but for the broad Balmer lines.
        minsnr_balmer_broad : :class:`float`
            Minimum signal-to-noise ratio for identifying a significant broad
            Balmer line.
        niter : :class:`int`
            Number of fit-in-patches iterations.
        debug_plots : :class:`bool`
            Optionally generate and write out various debugging plots.

        Returns
        -------
        :class:`dict` with line-masking arrays

        """
        from astropy.table import vstack
        from fastspecfit.util import quantile
        from fastspecfit.emlines import EMFitTools, ParamType

        def _make_patchTable(patchids):
            """Initialize the table containing information on each patch."""
            patchids = np.atleast_1d(patchids)
            npatch = len(patchids)
            continuum_patches = Table()
            continuum_patches['patchid'] = patchids
            continuum_patches['endpts'] = np.zeros((npatch,2), int) # starting index relative to coadd_wave
            continuum_patches['pivotwave'] = np.zeros(npatch)
            continuum_patches['slope'] = np.zeros(npatch)
            continuum_patches['intercept'] = np.zeros(npatch)
            continuum_patches['slope_bounds'] = np.broadcast_to([-1e2, +1e2], (npatch, 2))
            continuum_patches['intercept_bounds'] = np.broadcast_to([-1e5, +1e5], (npatch, 2))
            continuum_patches['balmerbroad'] = np.zeros(npatch, bool) # does this patch have a broad Balmer line?
            return continuum_patches

        
        def _fit_patches(continuum_patches, patchMap, linemodel, debug_plots=False, 
                         testBalmerBroad=False, minsnr=1.5, modelname='', 
                         png=None):
            """Iteratively fit all the lines in patches."""
            
            linesigmas = np.zeros(nline)
            linesigmas[EMFit.isBroad] = initsigma_broad
            linesigmas[EMFit.isNarrow] = initsigma_narrow
            if testBalmerBroad:
                linesigmas[EMFit.isBalmerBroad] = initsigma_balmer_broad

            linevshifts = np.zeros_like(linesigmas)
            linevshifts[EMFit.isBroad] = initvshift_broad
            linevshifts[EMFit.isNarrow] = initvshift_narrow
            if testBalmerBroad:
                linevshifts[EMFit.isBalmerBroad] = initvshift_balmer_broad

            initial_guesses = None

            t0 = time.time()
            for iiter in range(niter):

                # Build the line and continuum masks (only for lines in range).
                pix = self._linepix_and_contpix(wave, ivar, linetable_inrange, 
                                                linesigmas[EMFit.line_in_range], 
                                                linevshifts=linevshifts[EMFit.line_in_range], 
                                                patchMap=patchMap, redshift=redshift,
                                                log=self.log)

                # Check for fully dropped patches, which can happen if large
                # parts of the spectrum are masked.
                if len(pix['dropped']) > 0:
                    for patchid in np.atleast_1d(pix['dropped']):
                        drop = np.where(continuum_patches['patchid'] == patchid)[0][0]
                        continuum_patches.remove_row(drop)
                        patchMap.pop(patchid)

                # In the case that patches have been merged, update the
                # continuum_patches table and patchMap dictionary.
                if len(pix['merged']) > 0:
                    for newpatchid in np.atleast_1d(pix['merged']):
                        oldpatchids = list(newpatchid)
                        O = np.where(np.isin(continuum_patches['patchid'].value, oldpatchids))[0]

                        # update continuum_patches
                        new_continuum_patch = _make_patchTable(newpatchid)
                        new_continuum_patch['pivotwave'] = np.mean(continuum_patches[O]['pivotwave'])
                        new_continuum_patch['balmerbroad'] = np.any(continuum_patches[O]['balmerbroad'])
                        continuum_patches.remove_rows(O)
                        continuum_patches = vstack((continuum_patches, new_continuum_patch))

                        # update patchMap
                        newlines, newI, newJ = [], [], []
                        for oldpatchid in oldpatchids:
                            newlines.append(patchMap[oldpatchid][0])
                            newI.append(patchMap[oldpatchid][1])
                            newJ.append(patchMap[oldpatchid][2])
                            del patchMap[oldpatchid]
                        patchMap[newpatchid] = (np.hstack(newlines), np.hstack(newI), np.hstack(newJ))


                linemask = np.zeros(len(wave), bool) # False=masked
    
                # Determine the edges of each patch based on the continuum
                # (line-free) pixels of all the lines on that patch.
                for ipatch, patchid in enumerate(pix['patch_contpix'].keys()):
                    contindx = pix['patch_contpix'][patchid]
                    continuum_patches['endpts'][ipatch] = (np.min(contindx), np.max(contindx)) # should be sorted...
                    # initial guesses and bounds, but check for pathological distributions
                    lo, med, hi = quantile(flux[contindx], [0.05, 0.5, 0.95])
                    if lo < med and lo < hi:
                        continuum_patches['intercept'][ipatch] = med
                        continuum_patches['intercept_bounds'][ipatch] = [lo, hi]

                    # unmask the continuum patches
                    linemask[contindx] = True # True=unmasked

                # unmask the lines
                for line in pix['linepix'].keys():
                    linemask[pix['linepix'][line]] = True

                # only fit the pixels in the patches
                weights = np.sqrt(ivar * linemask)
    
                # Get initial guesses on the line-emission on the first iteration.
                if initial_guesses is None:
                    initial_guesses, param_bounds = EMFit._initial_guesses_and_bounds(
                        pix['linepix'], flux, contpix=pix['contpix'],
                        subtract_local_continuum=True,
                        initial_linesigma_broad=initsigma_broad, 
                        initial_linesigma_narrow=initsigma_narrow, 
                        initial_linesigma_balmer_broad=initsigma_balmer_broad,
                        initial_linevshift_broad=0.,
                        initial_linevshift_narrow=0.,
                        initial_linevshift_balmer_broad=0.,
                        log=self.log,
                    )

                # fit!
                linefit, contfit = EMFit.optimize(linemodel, initial_guesses,
                                                  param_bounds, wave,
                                                  flux, weights, redshift,
                                                  resolution_matrix, camerapix,
                                                  continuum_patches=continuum_patches, 
                                                  log=self.log)

                # Update the initial guesses as well as linesigmas and
                # linevshifts (for _linepix_and_contpix, at the top of the
                # iteration loop).
                if iiter < niter-1:
                    initial_guesses = linefit['value'].value.copy()
                    linevshifts = initial_guesses[EMFit.line_table['params'][:, ParamType.VSHIFT]]
                    linesigmas = initial_guesses[EMFit.line_table['params'][:, ParamType.SIGMA]]


            # Build the best-fitting model and estimate the S/N of each line.
            parameters = linefit['value'].value.copy()
            parameters[EMFit.doublet_idx] *= parameters[EMFit.doublet_src]

            lineamps = parameters[EMFit.line_table['params'][:, ParamType.AMPLITUDE]]
            linevshifts = parameters[EMFit.line_table['params'][:, ParamType.VSHIFT]]
            linesigmas = parameters[EMFit.line_table['params'][:, ParamType.SIGMA]]

            bestfit = EMFit.bestfit(linefit, redshift, wave, resolution_matrix, 
                                    camerapix, continuum_patches=contfit)
            residuals = flux - bestfit

            linesnrs = np.zeros_like(lineamps)
            noises = np.zeros(len(contfit))
            for ipatch, (patchid, endpts) in enumerate(contfit.iterrows('patchid', 'endpts')):
                s, e = endpts
                noise = np.ptp(quantile(residuals[s:e], (0.25, 0.75))) / 1.349 # robust sigma                
                noises[ipatch] = noise
                lineindx = patchMap[patchid][2] # index back to full line_table
                if noise != 0:
                    linesnrs[lineindx] = lineamps[lineindx] / noise

            # Derive the final linesigmas and linevshifts, and the maximum S/N
            # of each type of line. If the line isn't well-measured, drop the
            # S/N condition.
            strong = linesnrs > minsnr
            Ifree = linefit[EMFit.line_table['params'][:, ParamType.SIGMA]]['free']
            isBroad = EMFit.isBroad * Ifree * strong
            isNarrow = EMFit.isNarrow * Ifree * strong
            isBalmerBroad = EMFit.isBalmerBroad * Ifree * strong

            if np.any(isBroad):
                linesigma_broad = np.atleast_1d(linesigmas[isBroad])[0] # all values should be the same
                linevshift_broad = np.atleast_1d(linevshifts[isBroad])[0]
                maxsnr_broad = np.max(linesnrs[isBroad])
            else:
                isBroad = EMFit.isBroad * Ifree
                if np.any(isBroad):
                    linesigma_broad = np.atleast_1d(linesigmas[isBroad])[0]
                    linevshift_broad = np.atleast_1d(linevshifts[isBroad])[0]
                    maxsnr_broad = np.max(linesnrs[isBroad])
                else:
                    linesigma_broad = initsigma_broad # 0.
                    linevshift_broad = initvshift_broad
                    maxsnr_broad = 0.

            if np.any(isNarrow):
                linesigma_narrow = np.atleast_1d(linesigmas[isNarrow])[0]
                linevshift_narrow = np.atleast_1d(linevshifts[isNarrow])[0]
                maxsnr_narrow = np.max(linesnrs[isNarrow])
            else:
                isNarrow = EMFit.isNarrow * Ifree
                if np.any(isNarrow):
                    linesigma_narrow = np.atleast_1d(linesigmas[isNarrow])[0] 
                    linevshift_narrow = np.atleast_1d(linevshifts[isNarrow])[0]
                    maxsnr_narrow = np.max(linesnrs[isNarrow])
                else:
                    linesigma_narrow = initsigma_narrow
                    linevshift_narrow = initvshift_narrow
                    maxsnr_narrow = 0.

            if np.any(isBalmerBroad):
                linesigma_balmer_broad = np.atleast_1d(linesigmas[isBalmerBroad])[0]
                linevshift_balmer_broad = np.atleast_1d(linevshifts[isBalmerBroad])[0]
                maxsnr_balmer_broad = np.max(linesnrs[isBalmerBroad])
            else:
                isBalmerBroad = EMFit.isBalmerBroad * Ifree
                if np.any(isBalmerBroad):
                    linesigma_balmer_broad = np.atleast_1d(linesigmas[isBalmerBroad])[0]
                    linevshift_balmer_broad = np.atleast_1d(linevshifts[isBalmerBroad])[0]
                    maxsnr_balmer_broad = np.max(linesnrs[isBalmerBroad])
                else:
                    # we do not want to overmask if a broad Balmer line isn't detected 
                    linesigma_balmer_broad = 0. # initsigma_balmer_broad
                    linevshift_balmer_broad = initvshift_balmer_broad
                    maxsnr_balmer_broad = 0.

            final_linesigmas = (linesigma_broad, linesigma_narrow, linesigma_balmer_broad)
            final_linevshifts = (linevshift_broad, linevshift_narrow, linevshift_balmer_broad)
            maxsnrs = (maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad)

            self.log.info(f'Broad: S/N={maxsnr_broad:.1f}, (sigma,dv)=({linesigma_broad:.0f},{linevshift_broad:.0f}) km/s; ' + \
                          f'Narrow: S/N={maxsnr_narrow:.1f}, ({linesigma_narrow:.0f},{linevshift_narrow:.0f}) km/s; '+ \
                          f'Balmer Broad: S/N={maxsnr_balmer_broad:.1f}, ({linesigma_balmer_broad:.0f},{linevshift_balmer_broad:.0f}) km/s.')
            self.log.info(f'Fitting {modelname} with patches ({niter} iterations) took {time.time()-t0:.4f} seconds.')

            # optionally build a QA figure
            if debug_plots:
                import matplotlib.pyplot as plt
                import seaborn as sns
                from fastspecfit.emline_fit import EMLine_MultiLines

                sns.set(context='talk', style='ticks', font_scale=0.8)

                npatch = len(contfit)
                ncols = 3
                nrows = int(np.ceil(npatch / ncols))
    
                lines = EMLine_MultiLines(parameters, wave, redshift,
                                          linetable['restwave'].value,
                                          resolution_matrix, camerapix)
    
                def _niceline(line):
                    match line:
                        case 'lyalpha':
                            return r'S/N(Ly$\alpha$)='
                        case 'nv_1240':
                            return r'S/N(NV$\lambda1240$)='
                        case 'oi_1304':
                            return r'S/N(OI$\lambda1304$)='
                        case 'civ_1549': 
                            return r'S/N(CIV$\lambda1549$)='
                        case 'ciii_1908': 
                            return r'S/N(CIII]$\lambda1908$)='
                        case 'mgii_2796': 
                            return r'S/N(MgII$\lambda2796$)='
                        case 'mgii_2803': 
                            return r'S/N(MgII$\lambda2803$)='
                        case 'oii_3726': 
                            return r'S/N([OII]$\lambda3726$)='
                        case 'oii_3729': 
                            return r'S/N([OII]$\lambda3729$)='
                        case 'hgamma': 
                            return r'S/N(H$\gamma$)='
                        case 'hgamma_broad':
                            return r'S/N(H$\gamma_{b}$)='
                        case 'hbeta': 
                            return r'S/N(H$\beta$)='
                        case 'hbeta_broad': 
                            return r'S/N(H$\beta_{b}$)='
                        case 'oiii_4959': 
                            return r'S/N([OIII]$\lambda4959$)='
                        case 'oiii_5007': 
                            return r'S/N([OIII]$\lambda5007$)='
                        case 'nii_6548': 
                            return r'S/N([NII]$\lambda6548$)='
                        case 'halpha': 
                            return r'S/N(H$\alpha$)='
                        case 'halpha_broad': 
                            return r'S/N(H$\alpha_{b}$)='
                        case 'nii_6584': 
                            return r'S/N([NII]$\lambda6584$)='
                        case 'sii_6716': 
                            return r'S/N([SII]$\lambda6716$)='
                        case 'sii_6731':
                            return r'S/N([SII]$\lambda6731$)='
                        case _:
                            self.log.warning('New "niceline" needed.')
                            return 'Line'


                fig, ax = plt.subplots(nrows, ncols, figsize=(5.5*ncols, 5.5*nrows))
                for ipatch, ((patchid, endpts, slope, intercept, pivotwave), xx) in enumerate(
                        zip(contfit.iterrows('patchid', 'endpts', 'slope', 'intercept', 'pivotwave'), ax.flat)):
                    # get the starting and ending pixels first
                    s, e = endpts

                    for iline in patchMap[patchid][2]:
                        (ls, le), _ = lines.getLine(iline)
                        if ls != le: # fixed / dropped lines
                            s = np.min((s, ls))
                            e = np.max((e, le))
                    
                    xx.plot(wave[s:e] / 1e4, flux[s:e], color='gray')
                    xx.plot(wave[s:e] / 1e4, bestfit[s:e], color='k', ls='-', alpha=0.75)
                    cmodel = slope * (wave[s:e]-pivotwave) + intercept
                    xx.plot(wave[s:e] / 1e4, cmodel+noises[ipatch], color='gray', lw=1, ls='-')
                    xx.plot(wave[s:e] / 1e4, cmodel, color='k', lw=2, ls='--')
                    xx.plot(wave[s:e] / 1e4, cmodel-noises[ipatch], color='gray', lw=1, ls='-')
                    for line, iline in zip(patchMap[patchid][0], patchMap[patchid][2]):
                        (ls, le), profile = lines.getLine(iline)
                        if ls != le: # skip fixed lines
                            label = _niceline(line)+f'{linesnrs[iline]:.1f}; '+r'$\sigma$='+f'{linesigmas[iline]:.0f}'+' km/s'
                            xx.plot(wave[ls:le] / 1e4, profile, alpha=0.75, label=label)
                    if len(patchMap[patchid][0]) > 4:
                        nlegcol = 1
                        yfactor = 1.4
                    else:
                        nlegcol = 1
                        yfactor = 1.3
                    ymin = -1.2 * noises[ipatch] 
                    ymax = yfactor * np.max((quantile(flux[s:e], 0.99), np.max(bestfit[s:e])))
                    xx.set_ylim(ymin, ymax)
                    xx.legend(loc='upper left', fontsize=8, ncols=nlegcol)
                    xx.set_title(f'Patch {patchid}')
                for rem in np.arange(ncols*nrows-ipatch-1)+ipatch+1:
                    ax.flat[rem].axis('off')

                if ax.ndim == 1:
                    ulpos = ax[0].get_position()
                    llpos = ax[0].get_position()
                    lrpos = ax[-1].get_position()
                    dxlabel = 0.08
                    bottom = 0.14
                else:
                    ulpos = ax[0, 0].get_position()
                    llpos = ax[-1, 0].get_position()
                    lrpos = ax[-1, -1].get_position()
                    dxlabel = 0.07
                    bottom = 0.11

                xpos = (lrpos.x1 - llpos.x0) / 2. + llpos.x0
                ypos = llpos.y0 - dxlabel
                fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                         ha='center', va='center')
    
                xpos = ulpos.x0 - 0.09
                ypos = (ulpos.y1 - llpos.y0) / 2. + llpos.y0
                fig.text(xpos, ypos, r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                         ha='center', va='center', rotation=90)

                fig.subplots_adjust(left=0.08, right=0.97, bottom=bottom, top=0.92, wspace=0.23, hspace=0.3)
                
                if png:
                    fig.savefig(png, bbox_inches='tight')
                    plt.close()

            return linefit, contfit, final_linesigmas, final_linevshifts, maxsnrs


        # main function begins here
        camerapix = np.array([[0, len(wave)]]) # one camera

        # Read just the strong lines and determine which lines are in range of the camera.
        EMFit = EMFitTools(emlinesfile=self.emlinesfile, uniqueid=uniqueid, stronglines=True)
        EMFit.compute_inrange_lines(redshift, wavelims=(np.min(wave), np.max(wave)))
    
        # Build the narrow and narrow+broad emission-line models.
        linemodel_broad, linemodel_nobroad = EMFit.build_linemodels(separate_oiii_fit=False)
        
        # ToDo: are there ever *no* "strong" lines in range?
        linetable = EMFit.line_table
        linetable_inrange = linetable[EMFit.line_in_range]
        nline = len(linetable)

        # Initialize the continuum_patches table for all patches in range.
        patchids = np.unique(linetable_inrange['patch'])
        continuum_patches = _make_patchTable(patchids)

        # Get the mapping between patchid and the set of lines belonging to each
        # patch, and pivotwave.
        patchMap = {}
        for ipatch, patchid in enumerate(patchids):
            I = np.where(linetable_inrange['patch'] == patchid)[0] # index into fitted / in-range lines
            J = np.where(linetable['patch'] == patchid)[0]         # index into *all* lines
            patchlines = linetable_inrange['name'][I].value
            patchMap[patchid] = (patchlines, I, J)
            linewaves = linetable_inrange['restwave'][I] * (1. + redshift)
            pivotwave = 0.5 * (np.min(linewaves) + np.max(linewaves)) # midpoint
            continuum_patches['pivotwave'][ipatch] = pivotwave
            # is there a broad Balmer line on this patch?
            continuum_patches['balmerbroad'][ipatch] = np.any(EMFit.isBalmerBroad_noHelium_Strong[EMFit.line_in_range][I])
            
        # Need to pass copies of continuum_patches and patchMap because they can
        # get modified dynamically by _fit_patches.
        linefit_nobroad, contfit_nobroad, linesigmas_nobroad, linevshifts_nobroad, maxsnrs_nobroad = \
            _fit_patches(continuum_patches.copy(), patchMap.copy(), 
                         linemodel_nobroad, testBalmerBroad=False, 
                         debug_plots=debug_plots, modelname='narrow lines only',
                         png=f'qa-patches-nobroad-{uniqueid}.png')

        # Only fit with broad Balmer lines if at least one patch contains a
        # broad line.
        B = contfit_nobroad['balmerbroad']
        if np.any(B):
            linefit_broad, contfit_broad, linesigmas_broad, linevshifts_broad, maxsnrs_broad = \
                _fit_patches(continuum_patches.copy(), patchMap.copy(), 
                             linemodel_broad, testBalmerBroad=True, 
                             debug_plots=debug_plots, modelname='narrow+broad lines',
                             png=f'qa-patches-broad-{uniqueid}.png')

            # if a broad Balmer line is well-detected, take its linewidth
            if maxsnrs_broad[2] > minsnr_balmer_broad: 
                self.log.info(f'Adopting broad Balmer-line masking: S/N(broad Balmer) ' + \
                              f'{maxsnrs_broad[2]:.1f} > {minsnr_balmer_broad:.1f}.')
                finalsigma_broad, finalsigma_narrow, finalsigma_balmer_broad = linesigmas_broad
                finalvshift_broad, finalvshift_narrow, finalvshift_balmer_broad = linevshifts_broad
                maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad = maxsnrs_broad
            else:
                self.log.info(f'Adopting narrow Balmer-line masking: S/N(broad Balmer) ' + \
                              f'{maxsnrs_broad[2]:.1f} < {minsnr_balmer_broad:.1f}.')
                finalsigma_broad, finalsigma_narrow, finalsigma_balmer_broad = linesigmas_nobroad
                finalvshift_broad, finalvshift_narrow, finalvshift_balmer_broad = linevshifts_nobroad
                maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad = maxsnrs_nobroad
        else:
            self.log.info(f'Adopting narrow Balmer-line masking: no Balmer lines in wavelength range.')
            finalsigma_broad, finalsigma_narrow, finalsigma_balmer_broad = linesigmas_nobroad
            finalvshift_broad, finalvshift_narrow, finalvshift_balmer_broad = linevshifts_nobroad
            maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad = maxsnrs_nobroad

        # Build the final pixel mask for *all* lines using our current best
        # knowledge of the broad Balmer lines....(comment continued below)
        EMFit = EMFitTools(emlinesfile=self.emlinesfile, uniqueid=uniqueid, stronglines=False)
        EMFit.compute_inrange_lines(redshift, wavelims=(np.min(wave), np.max(wave)))

        linesigmas = np.zeros(len(EMFit.line_table))
        linesigmas[EMFit.isBroad] = finalsigma_broad
        linesigmas[EMFit.isNarrow] = finalsigma_narrow
        linesigmas[EMFit.isBalmerBroad] = finalsigma_balmer_broad

        linevshifts = np.zeros_like(linesigmas)
        linevshifts[EMFit.isBroad] = finalvshift_broad
        linevshifts[EMFit.isNarrow] = finalvshift_narrow
        linevshifts[EMFit.isBalmerBroad] = finalvshift_balmer_broad

        pix = self._linepix_and_contpix(wave, ivar, 
                                        EMFit.line_table[EMFit.line_in_range], 
                                        linesigmas[EMFit.line_in_range], 
                                        linevshifts=linevshifts[EMFit.line_in_range],
                                        patchMap=None, redshift=redshift)

        # Build another QA figure
        if debug_plots:
            import matplotlib.pyplot as plt
            import seaborn as sns

            sns.set(context='talk', style='ticks', font_scale=0.8)

            png = f'qa-linemask-{uniqueid}.png'

            linenames = list(pix['linepix'].keys())
            zlinewaves = EMFit.line_table[EMFit.line_in_range]['restwave'] * (1. + redshift)
            nline = len(linenames)

            ncols = 5
            nrows = int(np.ceil(nline / ncols))

            fig, ax = plt.subplots(nrows, ncols, figsize=(3*ncols, 2*nrows))

            for iline, (linename, xx) in enumerate(zip(linenames, ax.flat)):
                linepix = pix['linepix'][linename]
                contpix = pix['contpix'][linename]
                s = np.min((np.min(contpix), np.min(linepix)))
                e = np.max((np.max(contpix), np.max(linepix)))

                ylim = quantile(flux[s:e], [0.01, 0.99])

                xx.plot(wave[s:e] / 1e4, flux[s:e], color='gray', alpha=0.5)
                xx.scatter(wave[contpix] / 1e4, flux[contpix], color='blue', s=10, marker='s')
                xx.scatter(wave[linepix] / 1e4, flux[linepix], color='orange', s=10, marker='o', alpha=0.7)
                xx.axvline(x=zlinewaves[iline] / 1e4, color='k', ls='--', lw=2)
                xx.set_ylim(ylim[0], 1.3 * ylim[1])
                xx.text(0.1, 0.9, linename, ha='left', va='center', transform=xx.transAxes, 
                        fontsize=8)

            for rem in np.arange(ncols*nrows-iline-1)+iline+1:
                ax.flat[rem].axis('off')

            if ax.ndim == 1:
                ulpos = ax[0].get_position()
                urpos = ax[-1].get_position()
                llpos = ax[0].get_position()
                lrpos = ax[-1].get_position()
                top = 0.92
                bottom = 0.14
                dytitle = 0.13
                dyxlabel = 0.15
            else:
                ulpos = ax[0, 0].get_position()
                urpos = ax[0, -1].get_position()
                llpos = ax[-1, 0].get_position()
                lrpos = ax[-1, -1].get_position()
                top = 0.95
                bottom = 0.07
                dytitle = 0.09
                dyxlabel = 0.08

            xpos = (lrpos.x1 - llpos.x0) / 2. + llpos.x0
            ypos = llpos.y0 - dyxlabel
            fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                     ha='center', va='center')

            xpos = ulpos.x0 - 0.1
            ypos = (ulpos.y1 - llpos.y0) / 2. + llpos.y0# + 0.03
            fig.text(xpos, ypos, r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                     ha='center', va='center', rotation=90)

            xpos = (urpos.x1 - ulpos.x0) / 2. + ulpos.x0
            ypos = ulpos.y1 + dytitle
            fig.text(xpos, ypos, f'LinePix/ContPix: {uniqueid}', ha='center', va='center')

            fig.subplots_adjust(left=0.06, right=0.97, bottom=bottom, top=top, wspace=0.23, hspace=0.3)
            fig.savefig(png, bbox_inches='tight')
            plt.close()

        # (comment continued from above) ...but reset the broad Balmer
        # line-width to a minimum value and make another linepix mask. We need
        # to do this so that emlines.emline_specfit has a chance to remeasure
        # the broad Balmer lines after continuum-subtraction.
        if finalsigma_balmer_broad < initsigma_balmer_broad:
            finalsigma_balmer_broad = initsigma_balmer_broad
            linesigmas[EMFit.isBalmerBroad] = finalsigma_balmer_broad
            linevshifts[EMFit.isBalmerBroad] = finalvshift_balmer_broad
            newpix = self._linepix_and_contpix(wave, ivar, 
                                               EMFit.line_table[EMFit.line_in_range], 
                                               linesigmas[EMFit.line_in_range], 
                                               linevshifts=linevshifts[EMFit.line_in_range],
                                               patchMap=None, redshift=redshift)
            linepix = newpix['linepix']
        else:
            linepix = pix['linepix']
           
        out = {
            'linesigma_broad': finalsigma_broad, 
            'linesigma_narrow': finalsigma_narrow, 
            'linesigma_balmer_broad': finalsigma_balmer_broad, # updated value
            'linevshift_broad': finalvshift_broad, 
            'linevshift_narrow': finalvshift_narrow, 
            'linevshift_balmer_broad': finalvshift_balmer_broad, # updated value
            'maxsnr_broad': maxsnr_broad,
            'maxsnr_narrow': maxsnr_narrow,
            'maxsnr_balmer_broad': maxsnr_balmer_broad,
            'balmerbroad': np.any(contfit_nobroad['balmerbroad']), # True = one or more broad Balmer line in range
            'coadd_linepix': linepix,
        }

        return out


    @staticmethod
    def resample_and_smooth(templateflux, templatewave, specwave, specres):
        """Match an input template spectrum to the DESI wavelength sampling and
        instrumental resolution.

        Parameters
        ----------
        templateflux : :class:`numpy.ndarray` [npix]
            Input (model) spectrum.
        templatewave : :class:`numpy.ndarray` [npix]
            Wavelength array corresponding to `templateflux`.
        specwave : :class:`numpy.ndarray` [noutpix]
            Desired output wavelength array.
        specres : :class:`desispec.resolution.Resolution`
            Resolution matrix.

        Returns
        -------
        :class:`numpy.ndarray` [noutpix]
            Resampled and smoothed model spectrum.

        """
        from fastspecfit.util import trapz_rebin

        # Find lo, hi s.t. templatewave[lo:hi] is strictly within given
        # bounds. Should these pixel positions be cached??
        lo = np.searchsorted(templatewave, specwave[0]-2., 'right')
        hi = np.searchsorted(templatewave, specwave[-1]+2., 'left')
        resampflux = trapz_rebin(templatewave[lo:hi], templateflux[lo:hi], specwave)
        smoothflux = specres.dot(resampflux)

        return smoothflux

    
    @staticmethod
    def convolve_vdisp(templateflux, vdisp, pixsize_kms=PIXKMS_BLU, limit=None):
        return convolve_vdisp(templateflux, vdisp, pixsize_kms=pixsize_kms, limit=limit)


    def templates2data(self, templateflux, templatewave, redshift=0.0, dluminosity=None,
                       vdisp=None, cameras=['b', 'r', 'z'], specwave=None, specres=None, 
                       specmask=None, coeff=None, photsys=None, synthphot=True, 
                       stack_cameras=False, flamphot=False, debug=False, get_abmag=False):
        """Deprecated. Work-horse method for turning input templates into data-like
        spectra and photometry.

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
        This method does (n)one or more of the following:
        - redshifting
        - wavelength resampling
        - apply velocity dispersion broadening
        - apply the resolution matrix
        - synthesize photometry

        It also naturally handles templates which have been precomputed on a
        grid of velocity dispersion (and therefore have an additional
        dimension). However, if the input grid is 3D, it is reshaped to be 2D
        but then it isn't reshaped back because of the way the photometry table
        is organized (bug or feature?).

        """
        # Are we dealing with a 2D grid [npix,nage] or a 3D grid
        # [npix,nage,nAV] or [npix,nage,nvdisp]?
        ndim = templateflux.ndim
        if ndim == 2:
            npix, nsed = templateflux.shape
            nmodel = nsed
        elif ndim == 3:
            npix, nsed, nprop = templateflux.shape
            nmodel = nsed*nprop
            templateflux = templateflux.reshape(npix, nmodel)
        else:
            errmsg = f'Input templates have an unrecognized number of dimensions, {ndim}'
            self.log.critical(errmsg)
            raise ValueError(errmsg)
        
        # broaden for velocity dispersion but only out to ~1 micron
        if vdisp is not None:
            hi = np.searchsorted(templatewave, PIXKMS_WAVESPLIT, 'left')
            vd_templateflux = self.convolve_vdisp(templateflux, vdisp, 
                                                  pixsize_kms=PIXKMS_BLU, limit=hi)
        else:
            vd_templateflux = templateflux
            
        # Apply the redshift factor. The models are normalized to 10 pc, so
        # apply the luminosity distance factor here. Also normalize to a nominal
        # stellar mass.
        if redshift > 0.:
            # FIXME: try accelerating this whole thing in Numba
            ztemplatewave = templatewave * (1. + redshift)
            T = self.igm.full_IGM(redshift, ztemplatewave) 
            T *= FLUXNORM * self.massnorm * (10. / (1e6 * dluminosity))**2 / (1. + redshift)
            ztemplateflux = vd_templateflux * T[:, np.newaxis]
        else:
            errmsg = 'Input redshift not defined, zero, or negative!'
            self.log.warning(errmsg)
            ztemplatewave = templatewave
            T = FLUXNORM * self.massnorm
            ztemplateflux = vd_templateflux * T
        
        # Optionally synthesize photometry.
        templatephot = None
        if synthphot:
            if photsys is not None:
                filters = self.filters[photsys]
            else:
                filters = self.filters
            effwave = filters.effective_wavelengths.value

            if ((specwave is None and specres is None and coeff is None) or
               (specwave is not None and specres is not None)):
                
                maggies = self.get_ab_maggies(filters,
                                              ztemplateflux.T,
                                              ztemplatewave,
                                              log=self.log) / \
                                              (FLUXNORM * self.massnorm)
                if flamphot:                
                    templatephot = self.get_photflam(maggies, effwave, nanomaggies=False)
                else:
                    templatephot = self.parse_photometry(self.bands, maggies, effwave, 
                                                         nanomaggies=False, get_abmag=get_abmag, 
                                                         log=self.log)
                
        # Are we returning per-camera spectra or a single model? Handle that here.
        if specwave is None and specres is None:
            # cannot smooth/resample
            datatemplateflux = ztemplateflux
            
            # optionally compute the best-fitting model
            if coeff is not None:
                datatemplateflux = datatemplateflux.dot(coeff)
                if synthphot:
                    maggies = self.get_ab_maggies(filters,
                                                  datatemplateflux.T,
                                                  ztemplatewave,
                                                  log=self.log) / \
                                                  (FLUXNORM * self.massnorm)
                    
                    if flamphot:
                        templatephot = self.get_photflam(maggies, effwave, nanomaggies=False)
                    else:
                        templatephot = self.parse_photometry(self.bands, maggies, effwave, 
                                                             nanomaggies=False, get_abmag=get_abmag, 
                                                             log=self.log)

        else:
            # loop over cameras
            datatemplateflux = []
            for icamera in range(len(cameras)): # iterate on cameras
                #from scipy.ndimage import binary_dilation
                ## interpolate over pixels where the resolution matrix is masked
                #if specmask is not None and np.any(specmask[icamera] != 0):
                #    sw_mask = binary_dilation(specmask[icamera] != 0, iterations=2)
                #else:
                #    sw_mask = None
                
                _datatemplateflux = np.empty((len(specwave[icamera]), nmodel),
                                             dtype=ztemplateflux.dtype)
                for imodel in range(nmodel):
                    resampflux = self.resample_and_smooth(ztemplateflux[:, imodel],
                                                          ztemplatewave, 
                                                          specwave=specwave[icamera],
                                                          specres=specres[icamera])

                    ## interpolate over pixels where the resolution matrix is masked
                    ## FixMe -- move this to io.read_and_unpack
                    #if sw_mask is not None:
                    #    resampflux[sw_mask] = np.interp(specwave[icamera][sw_mask],
                    #                                    ztemplatewave,
                    #                                    ztemplateflux[:, imodel])
                        
                    _datatemplateflux[:,imodel] = resampflux
                
                #_datatemplateflux = np.vstack(_datatemplateflux).T
                if coeff is not None:
                    _datatemplateflux = _datatemplateflux.dot(coeff)
                datatemplateflux.append(_datatemplateflux)
                
        # Optionally stack and reshape (used in fitting).
        if stack_cameras:
            datatemplateflux = np.concatenate(datatemplateflux, axis=0)  # [npix,nsed*nprop] or [npix,nsed]
            if ndim == 3:
                nwavepix = np.sum([ len(sw) for sw in specwave ])
                datatemplateflux = datatemplateflux.reshape(nwavepix, nsed, nprop) # [npix,nsed,nprop]
                
        return datatemplateflux, templatephot # vector or 3-element list of [npix,nmodel] spectra


    @staticmethod
    def call_nnls(modelflux, flux, ivar, xparam=None, debug=False,
                  interpolate_coeff=False, xlabel=None, png=None, 
                  log=None):
        """Deprecated. Wrapper on `scipy.optimize.nnls`.
    
        Works with both spectroscopic and photometric input and with both 2D and
        3D model spectra.
    
        interpolate_coeff - return the interpolated coefficients when exploring
          an array or grid of xparam
    
        """
        from scipy.optimize import nnls
        from fastspecfit.util import find_minima, minfit
        
        if log is None:
            from desiutil.log import get_logger
            log = get_logger()
    
        if xparam is not None:
            nn = len(xparam)
    
        inverr = np.sqrt(ivar)
        bvector = flux * inverr
    
        # If xparam is None (equivalent to modelflux having just two
        # dimensions, [npix,nage]), assume we are just finding the
        # coefficients at some best-fitting value...
        if xparam is None:
            Amatrix = modelflux * inverr[:, np.newaxis]
            try:
                coeff, rnorm = nnls(A=Amatrix, b=bvector)
            except RuntimeError:
                coeff, _ = nnls(A=Amatrix, b=bvector, maxiter=Amatrix.shape[1] * 100)                
            chi2 = np.sum(ivar * (flux - modelflux.dot(coeff))**2)
            return coeff, chi2
    
        # ...otherwise iterate over the xparam (e.g., AV or vdisp) dimension.
        Amatrix = modelflux * inverr[:, np.newaxis, np.newaxis] # reshape into [npix/nband,nage,nAV/nvdisp]
        coeff, chi2grid = [], []
        for ii in range(nn):
            _coeff, _ = nnls(A=Amatrix[:, :, ii], b=bvector)
            chi2 = np.sum(ivar * (flux - modelflux[:, :, ii].dot(_coeff))**2)
            coeff.append(_coeff)
            chi2grid.append(chi2)
        coeff = np.array(coeff)
        chi2grid = np.array(chi2grid)
        
        try:
            imin = find_minima(chi2grid)[0]
            if debug:
                xbest, xerr, chi2min, warn, (a, b, c) = minfit(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2], return_coeff=True)
            else:
                xbest, xerr, chi2min, warn = minfit(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2])
        except:
            errmsg = 'A problem was encountered minimizing chi2.'
            log.warning(errmsg)
            imin, xbest, xerr, chi2min, warn, (a, b, c) = 0, 0.0, 0.0, 0.0, 1, (0., 0., 0.)
    
        if warn == 0:
            xivar = 1.0 / xerr**2
        else:
            chi2min = 0.0
            xivar = 0.0
            xbest = xparam[0]
    
        # optionally interpolate the coefficients
        if interpolate_coeff:
            from scipy.interpolate import interp1d
            if xbest == xparam[0]:
                bestcoeff = coeff[0, :]
            else:
                xindx = range(len(xparam))
                f = interp1d(xindx, coeff, axis=0)
                bestcoeff = f(np.interp(xbest, xparam, xindx))
        else:
            bestcoeff = None
    
            # interpolate the coefficients
            #np.interp(xbest, xparam, np.arange(len(xparam)))            
    
        if debug:
            if xivar > 0:
                leg = r'${:.1f}\pm{:.1f}$'.format(xbest, 1 / np.sqrt(xivar))
                #leg = r'${:.3f}\pm{:.3f}\ (\chi^2_{{min}}={:.2f})$'.format(xbest, 1/np.sqrt(xivar), chi2min)
            else:
                leg = r'${:.3f}$'.format(xbest)
                
            if xlabel:
                leg = f'{xlabel} = {leg}'
                
            import matplotlib.pyplot as plt
            import seaborn as sns
            sns.set(context='talk', style='ticks', font_scale=0.8)
            fig, ax = plt.subplots()
            ax.scatter(xparam, chi2grid-chi2min, marker='s', s=50, color='gray', edgecolor='k')
            #ax.scatter(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2]-chi2min, color='red')
            if not np.all(np.array([a, b, c]) == 0):
                ax.plot(xparam, np.polyval([a, b, c], xparam)-chi2min, lw=2, ls='--')
            #ax.axvline(x=xbest, color='k')
            #if xivar > 0:
            #    ax.axhline(y=chi2min, color='k')
            if xlabel:
                ax.set_xlabel(xlabel)
                #ax.text(0.03, 0.9, '{}={}'.format(xlabel, leg), ha='left',
                #        va='center', transform=ax.transAxes)
            ax.text(0.9, 0.9, leg, ha='right', va='center', transform=ax.transAxes)
            ax.set_ylabel(r'$\Delta\chi^2$')
            #ax.set_ylabel(r'$\chi^2 - {:.1f}$'.format(chi2min))
            #fig.savefig('qa-chi2min.png')
            fig.tight_layout()
            if png:
                log.info(f'Writing {png}')
                fig.savefig(png)
                plt.close()

        return chi2min, xbest, xivar, bestcoeff

    
    def klambda(self, wave):
        """Construct the total-to-selective attenuation curve, k(lambda).

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
            Input rest-frame wavelength array in Angstrom.

        Returns
        -------
        :class:`numpy.ndarray` [npix]
            Total-to-selective attenuation curve, k(lambda).

        Notes
        -----
        ToDo: support more sophisticated dust models.

        """
        return (wave / self.dust_normwave)**self.dust_power


    def _cache_klambda(self, *args):
        """Convenience method to cache the dust attenuation curve.

        """
        if self.dust_klambda is None:
            self.dust_klambda = self.klambda(*args)


    def _cache_zfactors(self, ztemplatewave, redshift, dluminosity=None):
        """Convenience method to cache the factors that depend on redshift and the
        redshifted wavelength array, including the attenuation due to the IGM.

        ztemplatewave : :class:`numpy.ndarray` [npix]
            Redshifted wavelength array.
        redshift : :class:`float`
            Object redshift.
        dluminosity : :class:`float`
            Luminosity distance corresponding to `redshift`.

        """
        if self.zfactors is None:
            if dluminosity is None:
                dluminosity = self.luminosity_distance(redshift)

            T = self.igm.full_IGM(redshift, ztemplatewave)
            T *= FLUXNORM * self.massnorm * (10. / (1e6 * dluminosity))**2 / (1. + redshift)

            self.zfactors = T


    def _cache_pixpos_wavesplit(self, templatewave):
        """Cache the pixel position in wave corresponding to `PIXKMS_WAVESPLIT`.

        Parameters
        ----------
        templatewave : :class:`numpy.ndarray`
            Input rest-frame wavelength array in Angstrom.

        """
        if self.pixpos_wavesplit is None:
            self.pixpos_wavesplit = np.searchsorted(templatewave, PIXKMS_WAVESPLIT, 'left')


    def build_stellar_continuum(self, templatewave, templateflux, templatecoeff,
                                dustflux=None, specwave=None, specres=None, 
                                camerapix=None, photsys=None, redshift=0., 
                                ebv=0., vdisp=None, synthphot=True, flamphot=False, 
                                get_abmag=False, filters=None):
        """Build a stellar continuum model.

        Parameters
        ----------
        templatewave : :class:`numpy.ndarray` [npix]
            Rest-frame template (model) wavelength array in Angstrom.
        templateflux : :class:`numpy.ndarray` [npix, ntemplate]
            Rest-frame, native-resolution template spectra corresponding to
            `templatewave`.
        templatecoeff : :class:`numpy.ndarray` [ntemplate]
            Column vector of positive coefficients corresponding to each
            template.
        dustflux : :class:`numpy.ndarray` [npix]
            Infrared dust emission spectrum. If provided, energy-balance is used
            to compute the normalization of this spectrum.
        specwave : :class:`numpy.ndarray` [nwave]
            Output (DESI) observed-frame wavelength array in Angstrom.
        specres : :class:`list` [ncamera]
            List of :class:`desispec.resolution.Resolution` objects, one per
            camera.
        camerapix : :class:`list` [ncamera]
            List of starting and ending :class:`int` pixel indices indicating
            the range of each camera.
        photsys : :class:`str`
            Photometric system.
        redshift : :class:`float`
            Object redshift.
        ebv : :class:`float`
            Dust extinction parameter, E(B-V), in mag.
        vdisp : :class:`float`
            Velocity dispersion in km/s. If `None`, do not convolve to the
            specified velocity dispersion (usually because `templateflux` has
            already been smoothed to some nominal value).
        synthphot : :class:`bool`
            Synthesize photometry from the observed-frame spectrum.
        flamphot : :class:`bool`
            Return synthesized photometry in f_lambda (10^17 erg/s/cm2/A) units,
            which are used in fitting. Otherwise, return an
            :class:`astropy.table.Table` with additional bandpass information.
        get_abmag : :class:`bool`
            Add AB magnitudes to the synthesized photometry table (only
            triggered if `synthphot=True` and `flamphot=False`. Usually only
            used for QA.
        filters : :class:`list` or :class:`speclite.filters.FilterSequence` [nfilt]
            Optionally override the filter curves stored in the `filters`
            attribute of this class.

        Returns
        -------
        fullmodel : :class:`numpy.ndarray` [npix]
            Full-wavelength, native-resolution, observed-frame model spectrum.
        modelflux : :class:`numpy.ndarray` [nwave]
            Observed-frame model spectrum at the instrumental resolution and
            wavelength sampling given by `specwave` and `specres`.
        modelphot : :class:`numpy.ndarray` or :class:`astropy.table.Table` [nfilt]
            Synthesized model photometry as an array or a table, depending on
            `flamphot`.

        Notes
        -----
        This method will compute 1-3 quantities depending on the input
        provided. 

        * It will always return `fullmodel`.
        * If `synthphot=True` it will synthesize observed-frame photometry from
          `fullmodel` and return `modelphot`.
        * If `specwave`, `specres`, and `camerapix` are all provided, it will
          generate a DESI-like spectrum and return `modelflux`.

        ToDo: `scipy.integrate.simpson` may not be the fastest integration
        method, so its use here is a placeholder.

        """
        from scipy.integrate import simpson

        # [1] - Compute the non-negative sum of the templates.
        fullmodel = templateflux.dot(templatecoeff)

        # [2] - Optionally convolve to the desired velocity dispersion.
        if vdisp is not None:
            fullmodel = self.convolve_vdisp(fullmodel, vdisp, 
                                            pixsize_kms=PIXKMS_BLU, 
                                            limit=self.pixpos_wavesplit)

        # [3] - Apply dust attenuation; ToDo: allow age-dependent
        # attenuation. Also compute the bolometric luminosity before and after
        # attenuation but only if we have dustflux.
        if self.dust_klambda is None:
            self._cache_klambda(templatewave)

        if dustflux is not None:
            lbol0 = simpson(fullmodel, x=templatewave)

        if ebv > 0.:
            atten = 10.**(-0.4 * ebv * self.dust_klambda)
            fullmodel *= atten
            if dustflux is not None:
                lbolabs = simpson(fullmodel, x=templatewave)

        # [4] - Optionally add the dust emission spectrum, conserving the total
        # (absorbed + re-emitted) energy. NB: (1) dustflux must be normalized to
        # unity already (see io.cache_templates); and (2) we are ignoring dust
        # self-absorption, which should be mostly negligible.
        if dustflux is not None:
            fullmodel += dustflux * (lbol0 - lbolabs)

        # [5] - Redshift factors.
        ztemplatewave = templatewave * (1. + redshift)
        if self.zfactors is None:
            dluminosity = self.luminosity_distance(redshift)
            self._cache_zfactors(ztemplatewave, redshift=redshift, 
                                 dluminosity=dluminosity)
        fullmodel *= self.zfactors

        # [6] - Optionally synthesize photometry
        if synthphot:
            if filters is None:
                if photsys is None:
                    filters = self.filters
                else:
                    filters = self.filters[photsys]

            effwave = filters.effective_wavelengths.value
            modelmaggies = self.get_ab_maggies_fast(filters, fullmodel,
                                                    ztemplatewave, log=self.log)
            if flamphot:                
                modelphot = self.get_photflam(modelmaggies, effwave, nanomaggies=False)
            else:
                modelmaggies /= FLUXNORM * self.massnorm
                modelphot = self.parse_photometry(self.bands, modelmaggies, effwave, 
                                                  nanomaggies=False, get_abmag=get_abmag, 
                                                  log=self.log)
        else:
            modelphot = None

        # [7] - Optionally handle the multi-camera spectroscopy. For each
        # camera, resample and multiply by the instrumental resolution.
        if specwave is not None and specres is not None and camerapix is not None:
            modelflux = np.empty(len(specwave))
            for icam, pix in enumerate(camerapix):
                modelflux[pix[0]:pix[1]] = self.resample_and_smooth(
                    fullmodel, ztemplatewave, specwave[pix[0]:pix[1]], 
                    specres[icam])
        else:
            modelflux = None

        return fullmodel, modelflux, modelphot

    
    def _stellar_objective(self, params, templatewave, templateflux, dustflux, 
                           objflam, objflamivar, specwave, specflux, specivar, 
                           specres, camerapix, photsys, redshift, fit_vdisp, 
                           synthphot):
        """Objective function for fitting a stellar continuum.
        
        """
        if fit_vdisp:
            ebv, vdisp  = params[:2]
            templatecoeff = params[2:]
        else:
            ebv = params[:1]
            templatecoeff = params[1:]
            vdisp = None

        fullmodel, modelflux, modelflam = self.build_stellar_continuum(
            templatewave, templateflux, templatecoeff, 
            dustflux=dustflux, specwave=specwave, specres=specres, 
            camerapix=camerapix, photsys=photsys, redshift=redshift, 
            ebv=ebv, vdisp=vdisp, synthphot=synthphot, flamphot=True)

        # Compute the residuals based on the information we have.
        if modelflux is None and modelflam is None:
            errmsg = 'Catastrophic problem constructing the model!'
            self.log.critical(errmsg)
            raise ValueError(errmsg)

        if modelflux is not None and modelflam is not None:
            spec_resid = np.sqrt(specivar) * (specflux - modelflux)
            phot_resid = np.sqrt(objflamivar) * (objflam - modelflam)
            resid = np.hstack((spec_resid, phot_resid))
        elif modelflux is not None and modelflam is None:
            resid = np.sqrt(specivar) * (specflux - modelflux)
        elif modelflux is None and modelflam is not None:
            resid = np.sqrt(objflamivar) * (objflam - modelflam)

        return resid
    
    
    def fit_stellar_continuum(self, templatewave, templateflux, dustflux=None,
                              objflam=None, objflamivar=None, specwave=None, 
                              specflux=None, specivar=None, specres=None, 
                              camerapix=None, photsys=None, redshift=0., 
                              ebv_guess=0.05, vdisp_guess=125., 
                              coeff_guess=None, ebv_bounds=(0., 3.), 
                              vdisp_bounds=(75., 500.), fit_vdisp=True):
        """Fit a stellar continuum using bounded non-linear least-squares.

        Parameters
        ----------
        templatewave : :class:`numpy.ndarray` [npix]
            Wavelength array corresponding to `templateflux`.
        templateflux : :class:`numpy.ndarray` [npix, ntemplate]
            Grid of input (model) spectra.
        dustflux : :class:`numpy.ndarray` [npix]
            Input (model) infrared dust spectrum.
        objflam : :class:`numpy.ndarray` [nfilt]
            Observed photometric fluxes in 10**-17 erg/s/cm2/A.
        objflamivar : :class:`numpy.ndarray` [nfilt]
            Inverse variance fluxes corresponding to `objflam`.
        specwave : :class:`numpy.ndarray` [nwave]
            Spectroscopic observed-frame wavelength array in Angstrom. The
            pixels corresponding to multiple cameras should be encoded in the
            `camerapix` array.
        specflux : :class:`numpy.ndarray` [nwave]
            Observed-frame spectrum in 10**-17 erg/s/cm2/A corresponding to
            `specwave`.
        specres : :class:`list` of :class:`desispec.resolution.Resolution` objects
            List of resolution matrices, with one element per camera (as encoded
            in `camerapix`).
        camerapix : :class:`list` [ncamera]
            List of a two-element list of starting and ending pixel numbers for
            each camera. E.g., [[0, 1100], [1050, 2200]] would indicate two
            camera with 50 pixels of overlap between them.
        photsys : :class:`str`
            Photometric system.
        redshift : :class:`float`
            Input object redshift.
        ebv_guess : :class:`float`
            Guess scalar value of the dust attenuation.
        vdisp_guess : :class:`float`
            Guess scalar value of the velocity dispersion.
        coeff_guess : :class:`numpy.ndarray` [nfilt]
            Guess of the template coefficients.
        ebv_bounds : :class:`tuple`
            Two-element array of minimum and maximum allowable values of the
            reddening, E(B-V).
        vdisp_bounds : :class:`tuple`
            Two-element array of minimum and maximum allowable values of the
            velocity dispersion; only used if `fit_vdisp=True`.
        fit_vdisp : :class:`bool`
            Optionally solve for the velocity dispersion.

        Returns
        -------
        ebv : :class:`float`
            Maximum-likelihood dust extinction parameter in mag.
        vdisp : :class:`float`
            Maximum-likelihood velocity dispersion in km/s.
        templatecoeff : :class:`numpy.ndarray` [ntemplate]
            Column vector of maximum-likelihood template coefficients.

        Notes
        -----
        This method supports several different fitting 'modes', depending on the
        inputs that are provided.

        The full set of possible modeling modes are: (1a) spectroscopy only,
        solving for velocity dispersion; (1b) spectroscopy only, fixed velocity
        dispersion; (2) spectroscopy and photometry, fixed velocity dispersion;
        (3) photometry only, fixed velocity dispersion.

        In all cases, we solve for dust attenuation via the E(B-V) parameter and
        we also include IGM attenuation.

        The required inputs for each mode are as follows:

        * Mode (1a): `specwave`, `specflux`, `specres`, `camerapix`, and
          `fit_vdisp=True`.
        * Mode (1b): `specwave`, `specflux`, `specres`, `camerapix`, and
          `fit_vdisp=False`.
        * Mode (2): `specwave`, `specflux`, `specres`, `camerapix`, `objflam`,
          `objflamivar`, and `fit_vdisp=False`.
        * Mode (3): `objflam`, `objflamivar`, and `fit_vdisp=False`.

        """
        from scipy.optimize import least_squares
    
        npix, ntemplate = templateflux.shape

        # Unpack the input data to infer the fitting "mode" and the objective
        # function.
        farg = [templatewave, templateflux, dustflux, objflam, 
                objflamivar, specwave, specflux, specivar, 
                specres, camerapix, photsys, redshift, fit_vdisp, ]

        if coeff_guess is None:
            coeff_guess = np.ones(ntemplate)
        else:
            if len(coeff_guess) != ntemplate:
                errmsg = 'Mismatching dimensions between coeff_guess and ntemplate!'
                self.log.critical(errmsg)
                raise ValueError(errmsg)
        
        if fit_vdisp == True:
            initial_guesses = np.hstack(([ebv_guess, vdisp_guess], coeff_guess))
            bounds = [ebv_bounds, vdisp_bounds, ] + [(0., 1e5)] * ntemplate
            #xscale = np.hstack(([0.1, 50.], np.ones(ntemplate) * 1e-1))
        else:
            initial_guesses = np.hstack((ebv_guess, coeff_guess))
            bounds = [ebv_bounds, ] + [(0., 1e5)] * ntemplate
            #xscale = np.hstack((0.1, np.ones(ntemplate) * 1e-1))

        # Don't synthesize photometry unless we have to.
        if objflam is not None and objflamivar is not None:
            synthphot = True
        else:
            synthphot = False
        farg += [synthphot, ]

        # NB: `x_scale` has been set to `jac` here to help with the numerical
        # convergence. There may be faster ways, of course...
        fit_info = least_squares(self._stellar_objective, initial_guesses, args=farg, 
                                 bounds=tuple(zip(*bounds)), method='trf',
                                 tr_solver='lsmr', tr_options={'regularize': True}, 
                                 x_scale='jac', max_nfev=5000, xtol=1e-10)#, verbose=2)
        bestparams = fit_info.x

        if fit_vdisp:
            ebv, vdisp  = bestparams[:2]
            templatecoeff = bestparams[2:]
        else:
            ebv = bestparams[0]
            templatecoeff = bestparams[1:]
            vdisp = vdisp_guess

        return ebv, vdisp, templatecoeff


    def stellar_continuum_chi2(self, templatecoeff, templatewave, templateflux, 
                               dustflux=None, objflam=None, objflamivar=None, 
                               specwave=None, specflux=None, specivar=None, 
                               specres=None, camerapix=None, photsys=None, 
                               redshift=0., ebv=0., vdisp=125., fit_vdisp=True):
        """Compute the reduced spectroscopic and/or photometric chi2.

        """
        if objflam is not None and objflamivar is not None:
            synthphot = True
        else:
            synthphot = False

        fullmodel, modelflux, modelflam = self.build_stellar_continuum(
            templatewave, templateflux, templatecoeff, dustflux=dustflux,
            specwave=specwave, specres=specres, camerapix=camerapix, 
            photsys=photsys, redshift=redshift, ebv=ebv, vdisp=vdisp, 
            synthphot=synthphot, flamphot=True)

        # ebv is always a free parameter
        nfree = len(templatecoeff) + 1 + 1 * fit_vdisp

        def _get_rchi2(chi2, ndof, nfree):
            """Guard against ndof=nfree."""
            if ndof > nfree:
                return chi2 / (ndof - nfree)
            else:
                return chi2 / ndof # ???
            

        # Compute chi2 based on the information we have.
        rchi2_spec, rchi2_phot, rchi2_tot = None, None, None

        if modelflux is not None:
            chi2_spec = np.sum(specivar * (specflux - modelflux)**2)
            ndof_spec = np.sum(specivar > 0)
            rchi2_spec = _get_rchi2(chi2_spec, ndof_spec, nfree)
            rchi2_tot = rchi2_spec

        if modelflam is not None:
            chi2_phot = np.sum(objflamivar * (objflam - modelflam)**2)
            ndof_phot = np.sum(objflamivar > 0)
            rchi2_phot = _get_rchi2(chi2_phot, ndof_phot, nfree)
            rchi2_tot = rchi2_phot

        if modelflux is not None and modelflam is not None:
            rchi2_tot = _get_rchi2(chi2_spec + chi2_phot, ndof_spec + ndof_phot, nfree)

        return fullmodel, rchi2_spec, rchi2_phot, rchi2_tot


    def continuum_fluxes(self, data, templatewave, continuum):
        """Compute rest-frame luminosities and observed-frame continuum fluxes.

        """
        from scipy.ndimage import median_filter
        from scipy.stats import sigmaclip
        
        def get_cflux(cwave):
            lo = np.searchsorted(templatewave, cwave - 20., 'right')
            hi = np.searchsorted(templatewave, cwave + 20., 'left')
            
            ksize = 200
            lo2 = np.maximum(0,              lo - ksize//2)
            hi2 = np.minimum(len(continuum), hi + ksize//2)
            smooth = median_filter(continuum[lo2:hi2], ksize)
            
            clipflux, _, _ = sigmaclip(smooth[lo-lo2:hi-lo2], low=1.5, high=3)
            return median(clipflux) # [flux in 10**-17 erg/s/cm2/A]

        
        redshift = data['zredrock']
        if redshift <= 0.0:
            errmsg = 'Input redshift not defined, zero, or negative!'
            self.log.warning(errmsg)
            return {}, {}
        
        dlum = data['dluminosity']

        # compute the model continuum flux at 1500 and 2800 A (to facilitate UV
        # luminosity-based SFRs) and at the positions of strong nebular emission
        # lines [OII], Hbeta, [OIII], and Halpha
        
        dfactor = (1. + redshift) * 4. * np.pi * (3.08567758e24 * dlum)**2 / FLUXNORM
        
        lums = {}
        cwaves = (1500.0, 2800.0, 1450., 1700., 3000., 5100.)
        labels = ('LOGLNU_1500', 'LOGLNU_2800', 'LOGL_1450', 'LOGL_1700', 'LOGL_3000', 'LOGL_5100')
        for cwave, label in zip(cwaves, labels):
            cflux = get_cflux(cwave) * dfactor # [monochromatic luminosity in erg/s/A]
            
            if 'LOGL_' in label:
                norm = 1e10
                cflux *= cwave / 3.846e33 / norm # [luminosity in 10**10 L_sun]
            else:
                # Convert the UV fluxes to rest-frame luminosity in erg/s/Hz. This
                # luminosity can be converted into a SFR using, e.g., Kennicutt+98,
                # SFR=1.4e-28 * L_UV
                norm = 1e28
                cflux *= cwave**2 / (C_LIGHT * 1e13) / norm # [monochromatic luminosity in 10**(-28) erg/s/Hz]

            if cflux > 0.:
                lums[label] = np.log10(cflux) # * u.erg/(u.second*u.Hz)
        
        cfluxes = {}
        cwaves = (1215.67, 3728.483, 4862.683, 5008.239, 6564.613)
        labels = ('FLYA_1215_CONT', 'FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT')
        for cwave, label in zip(cwaves, labels):
            cfluxes[label] = get_cflux(cwave) # * u.erg/(u.second*u.cm**2*u.Angstrom)

        return lums, cfluxes


    def kcorr_and_absmag(self, data, templatewave, continuum, snrmin=2.0):
        """Compute K-corrections and absolute magnitudes.

        """
        redshift = data['zredrock']
        
        # distance modulus, luminosity distance, and redshifted wavelength array
        dmod = data['dmodulus']
        ztemplatewave = templatewave * (1. + redshift)
        
        filters_in = self.filters[data['photsys']]

        maggies = data['phot']['nanomaggies'].data * 1e-9
        ivarmaggies = (data['phot']['nanomaggies_ivar'].data / 1e-9**2) * self.bands_to_fit

        kcorr, absmag, ivarabsmag, synth_absmag, synth_maggies_in = self.restframe_photometry(
            redshift, continuum, ztemplatewave, maggies, ivarmaggies,
            filters_in=filters_in, absmag_filters=self.absmag_filters,
            band_shift=self.band_shift, dmod=data['dmodulus'],
            snrmin=snrmin)
        
        return kcorr, absmag, ivarabsmag, synth_absmag, synth_maggies_in


def _younger_than_universe(age, tuniv, agepad=0.5):
    """Return the indices of the templates younger than the age of the universe
    (plus an agepadding amount) at the given redshift. age in yr, agepad and
    tuniv in Gyr

    """
    return np.where(age <= 1e9 * (agepad + tuniv))[0]


def continuum_specfit(data, result, templatecache, fphoto=None, emlinesfile=None,
                      constrain_age=False, no_smooth_continuum=False, 
                      ignore_photometry=False, fastphot=False, log=None, 
                      debug_plots=False, verbose=False):
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
        in :func:`init_output`.

    Notes
    -----
      - Consider using cross-correlation to update the redrock redshift.
      - We solve for velocity dispersion if ...

    """
    if log is None:
        from desiutil.log import get_logger, DEBUG
        if verbose:
            log = get_logger(DEBUG)
        else:
            log = get_logger()
            
    tall = time.time()

    CTools = ContinuumTools(fphoto=fphoto, emlinesfile=emlinesfile,
                            ignore_photometry=ignore_photometry, 
                            log=log)

    redshift = result['Z']
    if redshift <= 0.:
        log.warning('Input redshift not defined, zero, or negative!')

    photometry = data['phot']
    objflam = photometry['flam'].data * FLUXNORM
    objflamivar = (photometry['flam_ivar'].data / FLUXNORM**2) * CTools.bands_to_fit
    assert(np.all(objflamivar >= 0.))

    if np.any(CTools.bands_to_fit):
        # Require at least one photometric optical band; do not just fit the IR
        # because we will not be able to compute the aperture correction, below.
        lambda_eff = photometry['lambda_eff']
        opt = ((lambda_eff > 3e3) * (lambda_eff < 1e4))
        if np.all(objflamivar[opt] == 0.):
            log.warning('All optical bands are masked; masking all photometry.')
            objflamivar[:] = 0.0

    # Optionally ignore templates which are older than the age of the
    # universe at the redshift of the object.
    #nsed = len(templatecache['templateinfo'])
    if constrain_age:
        agekeep = _younger_than_universe(templatecache['templateinfo']['age'], data['tuniv'])
    else:
        agekeep = np.arange(len(templatecache['templateinfo']))
    nage = len(agekeep)

    ztemplatewave = templatecache['templatewave'] * (1. + redshift)
    vdisp_nominal = templatecache['vdisp_nominal']
    ebv_guess = 0.05

    # Cache k(lambda), the redshift-dependent factors (incl. IGM attenuation),
    # and a pixel position used when solving for velocity dispersion.
    CTools._cache_klambda(templatecache['templatewave'])
    CTools._cache_pixpos_wavesplit(templatecache['templatewave'])
    CTools._cache_zfactors(ztemplatewave, redshift=redshift, 
                           dluminosity=data['dluminosity'])

    # Photometry-only fitting.
    if fastphot:
        vdispbest, vdispivar = vdisp_nominal, 0.0
        log.info(f'Adopting nominal vdisp={vdisp_nominal:.0f} km/s.')

        if np.all(objflamivar == 0.):
            log.info('All photometry is masked.')
            coeff = np.zeros(nage, 'f4') # nage not nsed
            rchi2_cont, rchi2_phot = 0., 0.
            dn4000_model = 0.
            sedmodel = np.zeros(len(templatecache['templatewave']))
        else:
            # Get the coefficients and chi2 at the nominal velocity dispersion. 
            t0 = time.time()
            # maintain backwards-compatibility
            if templatecache['oldtemplates']:
                sedtemplates, sedphot_flam = CTools.templates2data(
                    templatecache['templateflux_nomvdisp'][:, agekeep],
                    templatecache['templatewave'], flamphot=True,
                    redshift=redshift, dluminosity=data['dluminosity'],
                    vdisp=None, synthphot=True, photsys=data['photsys'])
                sedflam = sedphot_flam * CTools.massnorm * FLUXNORM
                coeff, rchi2_phot = CTools.call_nnls(sedflam, objflam, objflamivar)
                rchi2_phot /= np.sum(objflamivar > 0) # dof???
            else:
                ebv, _, coeff = CTools.fit_stellar_continuum(
                    templatecache['templatewave'], 
                    templatecache['templateflux_nomvdisp'][:, agekeep], # [npix,nsed]
                    dustflux=templatecache['dustflux'], 
                    objflam=objflam, objflamivar=objflamivar, 
                    photsys=data['photsys'], redshift=redshift, 
                    ebv_guess=ebv_guess, vdisp_guess=vdispbest, 
                    fit_vdisp=False)

                sedmodel, _, rchi2_phot, rchi2_cont = CTools.stellar_continuum_chi2(
                    coeff, templatecache['templatewave'], 
                    templatecache['templateflux_nomvdisp'][:, agekeep],
                    dustflux=templatecache['dustflux'], 
                    objflam=objflam, objflamivar=objflamivar, 
                    photsys=data['photsys'], redshift=redshift, 
                    ebv=ebv, vdisp=None, fit_vdisp=False)

            log.info(f'Fitting {nage} models took {time.time()-t0:.2f} seconds.')
            
            if np.all(coeff == 0.):
                log.warning('Continuum coefficients are all zero.')
                sedmodel = np.zeros(len(templatecache['templatewave']))
                dn4000_model = 0.
            else:
                # Measure Dn(4000) from the line-free model.
                if templatecache['oldtemplates']:
                    sedmodel = sedtemplates.dot(coeff)
                    sedtemplates_nolines, _ = CTools.templates2data(
                        templatecache['templateflux_nolines_nomvdisp'][:, agekeep],
                        templatecache['templatewave'],
                        redshift=redshift, dluminosity=data['dluminosity'],
                        vdisp=None, synthphot=False)
                    sedmodel_nolines = sedtemplates_nolines.dot(coeff)
                else:
                    sedmodel_nolines, _, _ = CTools.build_stellar_continuum(                       
                        templatecache['templatewave'], 
                        templatecache['templateflux_nolines_nomvdisp'][:, agekeep],
                        coeff, dustflux=None, photsys=data['photsys'], 
                        redshift=redshift, ebv=ebv, vdisp=vdispbest, synthphot=False)

                    log.info(f'Best-fitting E(B-V)={ebv:.3f} mag.')
                    
                dn4000_model, _ = CTools.get_dn4000(templatecache['templatewave'],
                                                    sedmodel_nolines, rest=True, 
                                                    log=log)
                log.info(f'Model Dn(4000)={dn4000_model:.3f}.')
    else:
        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        specwave = np.hstack(data['wave'])
        specflux = np.hstack(data['flux'])
        flamivar = np.hstack(data['ivar'])
        specivar = flamivar * np.logical_not(np.hstack(data['linemask'])) # mask emission lines

        if np.all(specivar == 0.) or np.any(specivar < 0.):
            specivar = flamivar # not great...
            if np.all(specivar == 0.) or np.any(specivar < 0.):
                errmsg = 'All pixels are masked or some inverse variances are negative!'
                log.critical(errmsg)
                raise ValueError(errmsg)

        npix = len(specwave)

        # We'll need the filters for the aperture correction, below.
        filters_in = CTools.synth_filters[data['photsys']]

        # Solve for the velocity dispersion if the wavelength coverage is
        # sufficient.
        restwave = specwave / (1. + redshift)
        Ivdisp = np.where((specivar > 0) * (restwave > 3500.) * (restwave < 5500.))[0]
        compute_vdisp = (len(Ivdisp) > 0) and (np.ptp(restwave[Ivdisp]) > 500.)

        # stacked spectra do not have all three cameras
        if 'SNR_B' in result.columns and 'SNR_R' in result.columns and 'SNR_Z' in result.columns:
            log.info('S/N_b={:.2f}, S/N_r={:.2f}, S/N_z={:.2f}, rest wavelength coverage={:.0f}-{:.0f} A.'.format(
                result['SNR_B'], result['SNR_R'], result['SNR_Z'], restwave[0], restwave[-1]))

        # Maintain backwards compatibility. With the old templates, the velocity
        # dispersion and aperture corrections are determined separately, so we
        # separate that code out from the new templates, where they are
        # determined simultatneously.
        if templatecache['oldtemplates']:
            if compute_vdisp:
                t0 = time.time()
                ztemplateflux_vdisp, _ = CTools.templates2data(
                    templatecache['vdispflux'], templatecache['vdispwave'], # [npix,vdispnsed,nvdisp]
                    redshift=redshift, dluminosity=data['dluminosity'],
                    specwave=data['wave'], specres=data['res'],
                    cameras=data['cameras'], synthphot=False, stack_cameras=True)
                vdispchi2min, vdispbest, vdispivar, _ = CTools.call_nnls(
                    ztemplateflux_vdisp[Ivdisp, :, :], 
                    specflux[Ivdisp], specivar[Ivdisp],
                    xparam=templatecache['vdisp'], xlabel=r'$\sigma$ (km/s)', log=log,
                    debug=debug_plots, png='deltachi2-vdisp.png')
                log.info(f'Fitting for the velocity dispersion took {time.time()-t0:.2f} seconds.')

                if vdispivar > 0.:
                    # Require vdisp to be measured with S/N>1, which protects
                    # against tiny ivar becomming infinite in the output table.
                    vdispsnr = vdispbest * np.sqrt(vdispivar)
                    if vdispsnr < 1:
                        log.warning(f'vdisp signal-to-noise {vdispsnr:.2f} is less than ' + \
                                    f'one; adopting vdisp={vdisp_nominal:.0f} km/s.')
                        vdispbest, vdispivar = vdisp_nominal, 0.
                    else:
                        log.info(f'Best-fitting vdisp={vdispbest:.1f}+/-{1./np.sqrt(vdispivar):.1f} km/s.')
                else:
                    vdispbest = vdisp_nominal
                    log.info(f'Finding vdisp failed; adopting vdisp={vdisp_nominal:.0f} km/s.')
            else:
                vdispbest, vdispivar = vdisp_nominal, 0.
                log.info(f'Insufficient wavelength covereage to compute vdisp; adopting nominal vdisp={vdispbest:.2f} km/s')
    
            # Derive the aperture correction. 
            t0 = time.time()
    
            # First, do a quick fit of the DESI spectrum (including
            # line-emission templates) so we can synthesize photometry from a
            # noiseless model.
            if vdispbest == vdisp_nominal:
                # Use the cached templates.
                use_vdisp = None
                input_templateflux = templatecache['templateflux_nomvdisp'][:, agekeep]
                input_templateflux_nolines = templatecache['templateflux_nolines_nomvdisp'][:, agekeep]
            else:
                use_vdisp = vdispbest
                input_templateflux = templatecache['templateflux'][:, agekeep]
                input_templateflux_nolines = templatecache['templateflux_nolines'][:, agekeep]
    
            desitemplates, desitemplatephot_flam = CTools.templates2data(
                input_templateflux, templatecache['templatewave'],
                redshift=redshift, dluminosity=data['dluminosity'],
                specwave=data['wave'], specres=data['res'], specmask=data['mask'], 
                vdisp=use_vdisp, cameras=data['cameras'], stack_cameras=True, 
                synthphot=True, flamphot=True, photsys=data['photsys'])
            desitemplateflam = desitemplatephot_flam * CTools.massnorm * FLUXNORM
    
            apercorrs, apercorr = np.zeros(len(CTools.synth_bands), 'f4'), 0.
            
            sedtemplates, _ = CTools.templates2data(
                input_templateflux, templatecache['templatewave'],
                vdisp=use_vdisp, redshift=redshift,
                dluminosity=data['dluminosity'], synthphot=False)
    
            if not np.any(CTools.bands_to_fit):
                log.info('Skipping aperture correction since no bands were fit.')
                apercorrs, apercorr = np.ones(len(CTools.synth_bands), 'f4'), 1.
            else:
                # Fit using the templates with line-emission.
                quickcoeff, _ = CTools.call_nnls(desitemplates, specflux, specivar)
                if np.all(quickcoeff == 0.):
                    log.warning('Quick continuum coefficients are all zero.')
                else:
                    # Synthesize grz photometry from the full-wavelength SED to make
                    # sure we get the z-band correct.
                    nanomaggies = photometry['nanomaggies'].value
                    numer = np.hstack([nanomaggies[photometry['band'] == band] for band in CTools.synth_bands])
    
                    quicksedflux = sedtemplates.dot(quickcoeff)
                    quickmaggies = CTools.get_ab_maggies(filters_in,
                                                         quicksedflux / FLUXNORM,
                                                         ztemplatewave,
                                                         log=log)
                    
                    denom = CTools.to_nanomaggies(quickmaggies)
                    
                    I = (numer > 0.) * (denom > 0.)
                    if np.any(I):
                        apercorrs[I] = numer[I] / denom[I]
                    I = apercorrs > 0.
                    if np.any(I):
                        apercorr = median(apercorrs[I])
                        
                log.info(f'Median aperture correction = {apercorr:.3f} [{np.min(apercorrs):.3f}-{np.max(apercorrs):.3f}].')
                
                if apercorr <= 0.:
                    log.warning('Aperture correction not well-defined; adopting 1.0.')
                    apercorr = 1.
    
            data['apercorr'] = apercorr # needed for the line-fitting
    
            # Perform the final fit using the line-free templates in the spectrum
            # (since we mask those pixels) but the photometry synthesized from the
            # templates with lines.
            desitemplates_nolines, _ = CTools.templates2data(
                input_templateflux_nolines, templatecache['templatewave'], redshift=redshift,
                dluminosity=data['dluminosity'],
                specwave=data['wave'], specres=data['res'], specmask=data['mask'], 
                vdisp=use_vdisp, cameras=data['cameras'], stack_cameras=True, 
                synthphot=False)
    
            coeff, rchi2_cont = CTools.call_nnls(np.vstack((desitemplateflam, desitemplates_nolines)),
                                                 np.hstack((objflam, specflux * apercorr)),
                                                 np.hstack((objflamivar, specivar / apercorr**2)))
    
            # full-continuum fitting rchi2
            rchi2_cont /= (np.sum(objflamivar > 0) + np.sum(specivar > 0)) # dof???
            log.info(f'Final fitting with {nage} models took {time.time()-t0:.2f} seconds.')
    
            # rchi2 fitting just to the photometry, for analysis purposes
            rchi2_phot = np.sum(objflamivar * (objflam - desitemplateflam.dot(coeff))**2)
            if np.any(objflamivar > 0):
                rchi2_phot /= np.sum(objflamivar > 0)
    
            # Compute the full-wavelength best-fitting model.
            if np.all(coeff == 0):
                log.warning('Continuum coefficients are all zero.')
                sedmodel = np.zeros(len(templatecache['templatewave']), 'f4')
                desimodel = np.zeros_like(specflux)
                desimodel_nolines = np.zeros_like(specflux)
                dn4000_model = 0.0
            else:
                sedmodel = sedtemplates.dot(coeff)
                desimodel = desitemplates.dot(coeff)
                desimodel_nolines = desitemplates_nolines.dot(coeff)
    
                # Measure Dn(4000) from the line-free model.
                sedtemplates_nolines, _ = CTools.templates2data(
                    input_templateflux_nolines, templatecache['templatewave'], 
                    vdisp=use_vdisp, redshift=redshift, dluminosity=data['dluminosity'],
                    synthphot=False)
                sedmodel_nolines = sedtemplates_nolines.dot(coeff)
               
                dn4000_model, _ = CTools.get_dn4000(templatecache['templatewave'], 
                                                    sedmodel_nolines, rest=True, log=log)
            vdisp = use_vdisp
        else: # new templates start here

            # First, estimate the aperture correction from a (noiseless) *model*
            # of the spectrum (using the nominal velocity dispersion).
            apercorrs, apercorr = np.ones(len(CTools.synth_bands), 'f4'), 1.

            if np.any(CTools.bands_to_fit):
                t0 = time.time()
                ebv, _, coeff = CTools.fit_stellar_continuum(
                    templatecache['templatewave'], 
                    templatecache['templateflux_nomvdisp'][:, agekeep], # [npix,nsed]
                    dustflux=None, 
                    specwave=specwave, specflux=specflux, specivar=specivar,
                    specres=data['res'], camerapix=data['camerapix'],
                    redshift=redshift, ebv_guess=ebv_guess, 
                    vdisp_guess=None, fit_vdisp=False)

                if np.all(coeff == 0.):
                    log.warning('Unable to estimate aperture correction because continuum coefficients are all zero; adopting 1.0.')
                    coeff_guess = np.ones(nage)
                else:
                    sedmodel, _, sedflam = CTools.build_stellar_continuum(                       
                        templatecache['templatewave'], 
                        templatecache['templateflux_nomvdisp'][:, agekeep], coeff, 
                        dustflux=None, 
                        redshift=redshift, ebv=ebv, vdisp=vdisp_nominal, 
                        synthphot=True, flamphot=True, filters=filters_in)
    
                    objflam_aper = FLUXNORM * photometry[np.isin(photometry['band'], CTools.synth_bands)]['flam'].value

                    I = (objflam_aper > 0.) * (sedflam > 0.)
                    if np.any(I):
                        apercorrs[I] = objflam_aper[I] / sedflam[I]
    
                    I = apercorrs > 0.
                    if np.any(I):
                        apercorr = median(apercorrs[I])
    
                        log.info(f'Median aperture correction = {apercorr:.3f} ' + \
                                 f'[{np.min(apercorrs):.3f}-{np.max(apercorrs):.3f}].')
                        if apercorr <= 0.:
                            log.warning('Aperture correction not well-defined; adopting 1.0.')
                            apercorr = 1.

                    coeff_guess = coeff

                log.info(f'Deriving the aperture correction took {time.time()-t0:.2f} seconds.')
            else:
                log.info('Skipping aperture correction since no bands were fit.')
                coeff_guess = np.ones(nage)

            data['apercorr'] = apercorr # needed for the line-fitting

            # Now do the full spectrophotometric fit.
            t0 = time.time()
            ebv, vdisp, coeff = CTools.fit_stellar_continuum(
                templatecache['templatewave'], 
                templatecache['templateflux'][:, agekeep], # [npix,nsed]
                dustflux=templatecache['dustflux'], 
                specwave=specwave, specflux=specflux*apercorr, 
                specivar=specivar/apercorr**2,
                specres=data['res'], camerapix=data['camerapix'],
                objflam=objflam, objflamivar=objflamivar, 
                photsys=data['photsys'], redshift=redshift, 
                ebv_guess=ebv, vdisp_guess=vdisp_nominal,
                coeff_guess=coeff_guess, fit_vdisp=compute_vdisp)
            log.info(f'Final spectrophotometric fitting with {nage} models took {time.time()-t0:.2f} seconds.')

            # ToDo: capture case where vdisp (and also maybe ebv) hits its bounds.

            if np.all(coeff == 0.):
                log.warning('Continuum coefficients are all zero.')
                sedmodel = np.zeros(len(templatecache['templatewave']))
                dn4000_model = 0.
                rchi2_phot, rchi2_cont = 0., 0.
            else:
                vdispivar = 0. # ToDo
                ebvivar = 0.   # ToDo

                #log.info(f'Best-fitting E(B-V)={ebv:.3f}+/-{1./np.sqrt(ebvivar):.3f} mag.')
                log.info(f'Best-fitting E(B-V)={ebv:.3f} mag.')
    
                if compute_vdisp:
                    #log.info(f'Best-fitting vdisp={vdisp:.1f}+/-{1./np.sqrt(vdispivar):.1f} km/s.')
                    log.info(f'Best-fitting vdisp={vdisp:.1f} km/s.')
                else:
                    vdisp, vdispivar = vdisp_nominal, 0.
                    log.info(f'Insufficient wavelength covereage to compute vdisp; adopting nominal vdisp={vdisp:.1f} km/s')

                # rebuild the best-fitting model with and without line-emission
                if vdisp == vdisp_nominal:
                    # Use the cached templates.
                    use_vdisp = None
                    input_templateflux = templatecache['templateflux_nomvdisp'][:, agekeep]
                    input_templateflux_nolines = templatecache['templateflux_nolines_nomvdisp'][:, agekeep]
                else:
                    use_vdisp = vdisp
                    input_templateflux = templatecache['templateflux'][:, agekeep]
                    input_templateflux_nolines = templatecache['templateflux_nolines'][:, agekeep]

                sedmodel, _, _ = CTools.build_stellar_continuum(                       
                    templatecache['templatewave'], 
                    input_templateflux, coeff, 
                    dustflux=templatecache['dustflux'], 
                    redshift=redshift, ebv=ebv, 
                    vdisp=use_vdisp, synthphot=False)

                _, rchi2_spec, rchi2_phot, rchi2_cont = CTools.stellar_continuum_chi2(
                    coeff, templatecache['templatewave'], 
                    input_templateflux,
                    dustflux=templatecache['dustflux'], 
                    specwave=specwave, specflux=specflux*apercorr, 
                    specivar=specivar/apercorr**2,
                    specres=data['res'], camerapix=data['camerapix'],
                    objflam=objflam, objflamivar=objflamivar, 
                    photsys=data['photsys'], redshift=redshift, 
                    ebv=ebv, vdisp=use_vdisp, fit_vdisp=compute_vdisp)

                sedmodel_nolines, desimodel_nolines, _ = CTools.build_stellar_continuum(                       
                    templatecache['templatewave'], 
                    input_templateflux_nolines, coeff, 
                    dustflux=None, 
                    specwave=specwave, specres=data['res'], 
                    camerapix=data['camerapix'],
                    redshift=redshift, ebv=ebv, 
                    vdisp=vdisp, synthphot=False)

                dn4000_model, _ = CTools.get_dn4000(templatecache['templatewave'], 
                                                    sedmodel_nolines, rest=True, log=log)

        # Get DN(4000). Specivar is line-masked so we can't use it!
        dn4000, dn4000_ivar = CTools.get_dn4000(specwave, specflux, flam_ivar=flamivar, 
                                                redshift=redshift, rest=False, log=log)
        
        if dn4000_ivar > 0:
            log.info(f'Spectroscopic DN(4000)={dn4000:.3f}+/-{1./np.sqrt(dn4000_ivar):.3f}, ' + \
                     f'Model Dn(4000)={dn4000_model:.3f}')
        else:
            log.info(f'Spectroscopic DN(4000)={dn4000:.3f}, Model Dn(4000)={dn4000_model:.3f}')

        if np.all(coeff == 0.) or no_smooth_continuum:
            _smoothcontinuum = np.zeros_like(specwave)
        else:
            # Need to be careful we don't pass a large negative residual
            # where there are gaps in the data.
            residuals = apercorr * specflux - desimodel_nolines
            I = (specflux == 0.) * (specivar == 0.)
            if np.any(I):
                residuals[I] = 0.

            if debug_plots:
                png = f'qa-smooth-continuum-{data["uniqueid"]}.png'
            else:
                png = None

            linemask = np.hstack(data['linemask'])
            _smoothcontinuum = CTools.smooth_continuum(
                specwave, residuals, flamivar / apercorr**2, linemask=linemask, 
                redshift=redshift, png=png, log=log)

        # Unpack the continuum into individual cameras.
        continuummodel, smoothcontinuum = [], []
        for camerapix in data['camerapix']:
            continuummodel.append(desimodel_nolines[camerapix[0]:camerapix[1]])
            smoothcontinuum.append(_smoothcontinuum[camerapix[0]:camerapix[1]])

        for icam, cam in enumerate(data['cameras']):
            nonzero = continuummodel[icam] != 0
            if np.sum(nonzero) > 0:
                corr = median(smoothcontinuum[icam][nonzero] / continuummodel[icam][nonzero])
                result['SMOOTHCORR_{}'.format(cam.upper())] = corr * 100 # [%]

        if 'SMOOTHCORR_B' in result.columns and 'SMOOTHCORR_R' in result.columns and 'SMOOTHCORR_Z' in result.columns:
            log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
                result['SMOOTHCORR_B'], result['SMOOTHCORR_R'], result['SMOOTHCORR_Z']))
    
    # Compute K-corrections, rest-frame quantities, and physical properties.
    if np.all(coeff == 0):
        kcorr = np.zeros(len(CTools.absmag_bands))
        absmag = np.zeros(len(CTools.absmag_bands))#-99.0
        ivarabsmag = np.zeros(len(CTools.absmag_bands))
        synth_bestmaggies = np.zeros(len(CTools.bands))
        lums, cfluxes = {}, {}

        AV, age, zzsun, logmstar, sfr = 0.0, 0.0, 0.0, 0.0, 0.0
    else:
        kcorr, absmag, ivarabsmag, _, synth_bestmaggies = CTools.kcorr_and_absmag(
            data, templatecache['templatewave'], sedmodel)
        lums, cfluxes = CTools.continuum_fluxes(data, templatecache['templatewave'], sedmodel)

        # get the SPS properties
        mstars = templatecache['templateinfo']['mstar'][agekeep] # [current mass in stars, Msun]
        masstot = coeff.dot(mstars)
        coefftot = np.sum(coeff)
        logmstar = np.log10(CTools.massnorm * masstot)
        zzsun = np.log10(coeff.dot(mstars * 10.**templatecache['templateinfo']['zzsun'][agekeep]) / masstot) # mass-weighted
        age = coeff.dot(templatecache['templateinfo']['age'][agekeep]) / coefftot / 1e9          # luminosity-weighted [Gyr]
        #age = coeff.dot(mstars * templatecache['templateinfo']['age'][agekeep]) / masstot / 1e9 # mass-weighted [Gyr]
        sfr = CTools.massnorm * coeff.dot(templatecache['templateinfo']['sfr'][agekeep])                           # [Msun/yr]
        if templatecache['oldtemplates']:
            AV = coeff.dot(templatecache['templateinfo']['av'][agekeep]) / coefftot # luminosity-weighted [mag]
        else:
            AV = ebv * CTools.klambda(5500.) # [mag]

    rindx = np.argmin(np.abs(CTools.absmag_filters.effective_wavelengths.value / (1.+CTools.band_shift) - 5600))
    log.info(f'log(M/Msun)={logmstar:.2f}, M{CTools.absmag_bands[rindx]}={absmag[rindx]:.2f} mag, ' + \
             f'A(V)={AV:.3f}, Age={age:.3f} Gyr, SFR={sfr:.3f} Msun/yr, Z/Zsun={zzsun:.3f}')

    # Pack it in and return.
    result['COEFF'][agekeep] = coeff
    result['RCHI2_PHOT'] = rchi2_phot
    result['VDISP'] = vdisp # * u.kilometer/u.second
    result['AV'] = AV # * u.mag
    result['AGE'] = age # * u.Gyr
    result['ZZSUN'] = zzsun
    result['LOGMSTAR'] = logmstar
    result['SFR'] = sfr
    result['DN4000_MODEL'] = dn4000_model

    for iband, (band, shift) in enumerate(zip(CTools.absmag_bands, CTools.band_shift)):
        result['KCORR{:02d}_{}'.format(int(10*shift), band.upper())] = kcorr[iband] # * u.mag
        result['ABSMAG{:02d}_{}'.format(int(10*shift), band.upper())] = absmag[iband] # * u.mag
        result['ABSMAG{:02d}_IVAR_{}'.format(int(10*shift), band.upper())] = ivarabsmag[iband] # / (u.mag**2)
    for iband, band in enumerate(CTools.bands):
        result['FLUX_SYNTH_PHOTMODEL_{}'.format(band.upper())] = 1e9 * synth_bestmaggies[iband] # * u.nanomaggy
    if bool(lums):
        for lum in lums.keys():
            result[lum] = lums[lum]
    if bool(cfluxes):
        for cflux in cfluxes.keys():
            result[cflux] = cfluxes[cflux]

    if not fastphot:
        result['RCHI2_CONT'] = rchi2_cont
        result['VDISP_IVAR'] = vdispivar # * (u.second/u.kilometer)**2
        
        result['APERCORR'] = apercorr
        for iband, band in enumerate(CTools.synth_bands):
            result['APERCORR_{}'.format(band.upper())] = apercorrs[iband]
        result['DN4000_OBS'] = dn4000
        result['DN4000_IVAR'] = dn4000_ivar

    log.info('Continuum-fitting took {:.2f} seconds.'.format(time.time()-tall))

    if fastphot:
        return sedmodel, None
    else:
        # divide out the aperture correction
        continuummodel = [_continuummodel / apercorr for _continuummodel in continuummodel]
        smoothcontinuum = [_smoothcontinuum / apercorr for _smoothcontinuum in smoothcontinuum]
        return continuummodel, smoothcontinuum

