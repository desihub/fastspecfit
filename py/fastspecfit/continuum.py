"""
fastspecfit.continuum
=====================

Methods and tools for continuum-fitting.

"""
import pdb # for debugging

import time
import numpy as np
import numba

from fastspecfit.logger import log
from fastspecfit.photometry import Photometry
from fastspecfit.templates import Templates
from fastspecfit.util import (C_LIGHT, FLUXNORM,
    trapz_rebin, trapz_rebin_pre, quantile,
    median, sigmaclip)

class ContinuumTools(object):
    """Tools for dealing with spectral continua.

    """
    def __init__(self, igm, phot, templates, data, fastphot=False):

        self.phot = phot
        self.igm = igm  # only needed by legacy fitting
        self.templates = templates
        self.data = data

        self.massnorm = 1e10       # stellar mass normalization factor [Msun]

        self.pixpos_wavesplit = templates.pixpos_wavesplit
        self.lg_atten         = np.log(10.) * (-0.4 * templates.dust_klambda)

        # Cache the redshift-dependent factors (incl. IGM attenuation),
        redshift = data['zredrock']
        self.ztemplatewave = templates.wave * (1. + redshift)
        self.zfactors = self.get_zfactors(igm,
                                          self.ztemplatewave,
                                          redshift=redshift,
                                          dluminosity=data['dluminosity'])

        # Get preprocessing data to accelerate continuum_to_photometry()
        # but ONLY when it is called with the default filters=None
        photsys = self.data['photsys']
        if photsys is None:
            filters = self.phot.filters
        else:
            filters = self.phot.filters[photsys]

        self.phot_pre = (
            filters,
            filters.effective_wavelengths.value,
            Photometry.get_ab_maggies_pre(filters, self.ztemplatewave)
        )

        if not fastphot:
            self.wavelen = np.sum([len(w) for w in self.data['wave']])

            # get preprocessing data to accelerate resampling
            rspre = []
            for wave in self.data['wave']:
                rspre.append(trapz_rebin_pre(wave))
            self.spec_pre = tuple(rspre)


    def get_zfactors(self, igm, ztemplatewave, redshift, dluminosity):
        """Convenience method to cache the factors that depend on redshift and the
        redshifted wavelength array, including the attenuation due to the IGM.

        ztemplatewave : :class:`numpy.ndarray` [npix]
            Redshifted wavelength array.
        redshift : :class:`float`
            Object redshift.
        dluminosity : :class:`float`
            Luminosity distance corresponding to `redshift`.

        """
        T = igm.full_IGM(redshift, ztemplatewave)
        T *= FLUXNORM * self.massnorm * (10. / (1e6 * dluminosity))**2 / (1. + redshift)

        return T


    @staticmethod
    def smooth_continuum(wave, flux, ivar, linemask, camerapix=None, medbin=175,
                         smooth_window=75, smooth_step=125, png=None):
        """Build a smooth, nonparametric continuum spectrum.

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
           Observed-frame wavelength array.
        flux : :class:`numpy.ndarray` [npix]
           Spectrum corresponding to `wave`.
        ivar : :class:`numpy.ndarray` [npix]
           Inverse variance spectrum corresponding to `flux`.
        linemask : :class:`numpy.ndarray` of type :class:`bool`, optional, defaults to `None`
           Boolean mask with the same number of pixels as `wave` where `True`
           means a pixel is (possibly) affected by an emission line
           (specifically a strong line which likely cannot be median-smoothed).
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
        png : :class:`str`, optional, defaults to `None`
           Generate a simple QA plot and write it out to this filename.

        Returns
        -------
        smooth :class:`numpy.ndarray` [npix]
        Smooth continuum spectrum which can be subtracted from `flux` in
        order to create a pure emission-line spectrum.

        """
        from numpy.lib.stride_tricks import sliding_window_view
        from scipy.ndimage import median_filter
        from scipy.interpolate import make_interp_spline

        npix = len(wave)

        if len(linemask) != npix:
            errmsg = 'Linemask must have the same number of pixels as the input spectrum.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        # If camerapix is given, mask the 10 (presumably noisy) pixels from the edge
        # of each per-camera spectrum.
        if camerapix is not None:
            for campix in camerapix:
                ivar[campix[0]:campix[0]+10] = 0.
                ivar[campix[1]-10:campix[1]] = 0.

        # Build the smooth (line-free) continuum by computing statistics in a
        # sliding window, accounting for masked pixels and trying to be smart
        # about broad lines. See:
        #   https://stackoverflow.com/questions/41851044/python-median-filter-for-1d-numpy-array
        #   https://numpy.org/devdocs/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html
        wave_win = sliding_window_view(wave, window_shape=smooth_window)
        flux_win = sliding_window_view(flux, window_shape=smooth_window)
        ivar_win = sliding_window_view(ivar, window_shape=smooth_window)
        noline_win = sliding_window_view(np.logical_not(linemask), window_shape=smooth_window)

        nminpix = 15

        smooth_wave, smooth_flux, smooth_sigma = [], [], []
        for swave, sflux, sivar, noline in zip(wave_win[::smooth_step],
                                               flux_win[::smooth_step],
                                               ivar_win[::smooth_step],
                                               noline_win[::smooth_step]):

            # if there are fewer than XX good pixels after accounting for the
            # line-mask, skip this window.
            sflux = sflux[noline]
            if len(sflux) < nminpix:
                continue

            cflux, clo, chi = sigmaclip(sflux, low=2.0, high=2.0)
            if len(cflux) < nminpix:
                continue

            sivar = sivar[noline]
            # Toss out regions with too little good data.
            if np.sum(sivar > 0) < nminpix:
                continue

            sig = np.std(cflux) # simple median and sigma
            mn = median(cflux)

            # One more check for crummy spectral regions.
            if mn == 0.0:
                continue

            swave = swave[noline]
            I = ((sflux >= clo) & (sflux <= chi))
            smooth_wave.append(np.mean(swave[I]))

            smooth_sigma.append(sig)
            smooth_flux.append(mn)

        smooth_wave = np.array(smooth_wave)
        smooth_sigma = np.array(smooth_sigma)
        smooth_flux = np.array(smooth_flux)

        # For debugging.
        if png:
            _smooth_wave = smooth_wave.copy()
            _smooth_flux = smooth_flux.copy()

        # corner case for very wacky spectra
        if len(smooth_flux) == 0:
            smooth_flux = flux
            #smooth_sigma = flux * 0 + np.std(flux)
        else:
            #smooth_flux = np.interp(wave, smooth_wave, smooth_flux)
            #smooth_sigma = np.interp(wave, smooth_wave, smooth_sigma)
            _, uindx = np.unique(smooth_wave, return_index=True)
            srt = np.argsort(smooth_wave[uindx])
            bspl_flux = make_interp_spline(smooth_wave[uindx][srt], smooth_flux[uindx][srt], k=1)
            smooth_flux = bspl_flux(wave)

            # check for extrapolation
            blu = np.where(wave < np.min(smooth_wave[srt]))[0]
            red = np.where(wave > np.max(smooth_wave[srt]))[0]
            if len(blu) > 0:
                smooth_flux[:blu[-1]] = smooth_flux[blu[-1]]
            if len(red) > 0:
                smooth_flux[red[0]:] = smooth_flux[red[0]]

        smooth = median_filter(smooth_flux, medbin, mode='nearest')

        # very important!
        Z = (flux == 0.0) * (ivar == 0.0)
        if np.sum(Z) > 0:
            smooth[Z] = 0.0

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
            #ax[0].plot(_smooth_wave / 1e4, _smooth_flux, color='orange')
            #ax[0].plot(wave, median_filter(flux, medbin, mode='nearest'), color='k', lw=2)
            ax[0].set_ylim(np.min((-5. * noise, quantile(flux, 0.05))),
                           np.max((5. * noise, 1.5 * quantile(flux, 0.975))))
            ax[0].set_ylabel('Continuum-subtracted Spectrum')
            #ax[0].scatter(_smooth_wave, _smooth_flux, color='orange', marker='s', ls='-', s=20)
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


    def continuum_fluxes(self, continuum):
        """Compute rest-frame luminosities and observed-frame continuum fluxes.

        """
        from scipy.ndimage import median_filter

        def get_cflux(cwave):
            templatewave = self.templates.wave
            lo = np.searchsorted(templatewave, cwave - 20., 'right')
            hi = np.searchsorted(templatewave, cwave + 20., 'left')

            ksize = 200
            lo2 = np.maximum(0,              lo - ksize//2)
            hi2 = np.minimum(len(continuum), hi + ksize//2)
            smooth = median_filter(continuum[lo2:hi2], ksize)

            clipflux, _, _ = sigmaclip(smooth[lo-lo2:hi-lo2], low=1.5, high=3)
            return median(clipflux) # [flux in 10**-17 erg/s/cm2/A]

        redshift = self.data['zredrock']
        if redshift <= 0.0:
            log.warning('Input redshift not defined, zero, or negative!')
            return {}, {}

        # compute the model continuum flux at 1500 and 2800 A (to facilitate UV
        # luminosity-based SFRs) and at the positions of strong nebular emission
        # lines [OII], Hbeta, [OIII], and Halpha
        dlum = self.data['dluminosity']
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


    # FIXME:
    # * templatewave, redshift, luminosity, cameras, specres, specmask, photsys
    #    all come from stored object data now
    # * specmask ALSO comes from data, but we don't always use it when it exists,
    #    so it needs to remain an argument
    # * replace flamphot with phottable to get it out of this file and into qa
    # * introduce synthspec flag and use it instead of checking for specwave/specres
    # * don't handle redshift <= 0 here -- caller checks it in fastspecfit, and
    #     caller should check it in fastqa if not being done now
    # after we get all of that working, we can consider internal simplifications
    def templates2data(self, templateflux, templatewave, redshift=0., dluminosity=None,
                       vdisp=None, cameras=['b','r','z'], specwave=None, specres=None,
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
            log.critical(errmsg)
            raise ValueError(errmsg)

        # broaden for velocity dispersion but only out to ~1 micron
        if vdisp is not None:
            vd_templateflux = Templates.convolve_vdisp(templateflux, vdisp,
                                                       pixsize_kms=Templates.PIXKMS_BLU,
                                                       limit=self.pixpos_wavesplit)
        else:
            vd_templateflux = templateflux

        # Apply the redshift factor. The models are normalized to 10 pc, so
        # apply the luminosity distance factor here. Also normalize to a nominal
        # stellar mass.
        if redshift > 0.:
            ztemplatewave = templatewave * (1. + redshift)
            T = self.igm.full_IGM(redshift, ztemplatewave)
            T *= FLUXNORM * self.massnorm * (10. / (1e6 * dluminosity))**2 / (1. + redshift)
            ztemplateflux = vd_templateflux * T[:, np.newaxis]
        else:
            log.warning('Input redshift not defined, zero, or negative!')
            ztemplatewave = templatewave
            T = FLUXNORM * self.massnorm
            ztemplateflux = vd_templateflux * T

        # Optionally synthesize photometry.
        templatephot = None
        if synthphot:
            if photsys is not None:
                filters = self.phot.filters[photsys]
            else:
                filters = self.phot.filters
            effwave = filters.effective_wavelengths.value

            if ((specwave is None and specres is None and coeff is None) or
               (specwave is not None and specres is not None)):

                maggies = Photometry.get_ab_maggies(filters,
                                                    ztemplateflux.T,
                                                    ztemplatewave) / \
                                                    (FLUXNORM * self.massnorm)
                if flamphot:
                    templatephot = Photometry.get_photflam(maggies, effwave).T
                else:
                    templatephot = Photometry.parse_photometry(self.phot.bands, maggies, effwave,
                                                               nanomaggies=False, get_abmag=get_abmag)

        # Are we returning per-camera spectra or a single model? Handle that here.
        if specwave is None and specres is None:
            # cannot smooth/resample
            datatemplateflux = ztemplateflux

            # optionally compute the best-fitting model
            if coeff is not None:
                datatemplateflux = datatemplateflux.dot(coeff)
                if synthphot:
                    maggies = Photometry.get_ab_maggies(filters,
                                                        datatemplateflux.T,
                                                        ztemplatewave) / \
                                                        (FLUXNORM * self.massnorm)

                    if flamphot:
                        templatephot = Photometry.get_photflam(maggies, effwave).T
                    else:
                        templatephot = Photometry.parse_photometry(self.phot.bands, maggies, effwave,
                                                                   nanomaggies=False, get_abmag=get_abmag)

        else:
            # loop over cameras
            datatemplateflux = []
            for icam in range(len(cameras)): # iterate on cameras
                _datatemplateflux = np.empty((len(specwave[icam]), nmodel),
                                             dtype=ztemplateflux.dtype)
                for imodel in range(nmodel):
                    resampflux = trapz_rebin(ztemplatewave,
                                             ztemplateflux[:, imodel],
                                             specwave=specwave[icam],
                                             pre=self.spec_pre[icam])
                    _datatemplateflux[:, imodel] = specres[icam].dot(resampflux)

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
                  interpolate_coeff=False, xlabel=None, png=None):
        """Deprecated. Wrapper on `scipy.optimize.nnls`.

        Works with both spectroscopic and photometric input and with both 2D and
        3D model spectra.

        interpolate_coeff - return the interpolated coefficients when exploring
          an array or grid of xparam

        """
        from scipy.optimize import nnls
        from fastspecfit.util import find_minima, minfit

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
            log.warning('A problem was encountered minimizing chi2.')
            imin, xbest, xerr, chi2min, warn, (a, b, c) = 0, 0.0, 0.0, 0.0, 1, (0., 0., 0.)

        if warn == 0:
            xivar = 1.0 / xerr**2
        else:
            chi2min = 0.0
            xivar = 0.0
            xbest = xparam[0]

        # optionally interpolate the coefficients
        if interpolate_coeff:
            if xbest == xparam[0]:
                bestcoeff = coeff[0, :]
            else:
                indxbest = np.interp(xbest, xparam, range(len(xparam)))
                bestcoeff = np.interp(indxbest, xindx, coeff)
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


    # compute attenuated version of a model spectrum,
    # including contribution of dust emission
    @staticmethod
    @numba.jit(nopython=True, fastmath=True, nogil=True)
    def attenuate(M, A, zfactors, wave, dustflux):

        # Concurrently replace M by M * (atten ** ebv) and
        # compute (by trapezoidal integration) integral of
        # difference of bolometric luminosities before and after
        # attenuation at each wavelength.
        #
        # The integration is equivalent to
        # lbol0   = trapz(M, x=wave)
        # lbolabs = trapz(M*(atten**ebv), x=wave)
        # lbol_diff = lbol0 - lbolabs

        ma = M[0] * A[0]
        prev_d = M[0] - ma
        M[0] = ma

        lbol_diff = 0.
        for i in range(len(M) - 1):
            ma = M[i+1] * A[i+1]
            d = M[i+1] - ma
            M[i+1] = ma

            lbol_diff += (wave[i+1] - wave[i]) * (d + prev_d)

            prev_d  = d
        lbol_diff *= 0.5

        # final result is
        # (M * (atten ** ebv) + lbol_diff * dustflux) * zfactors
        for i in range(len(M)):
            M[i] = (M[i] + lbol_diff * dustflux[i]) * zfactors[i]

    # compute attenuated version of a model spectrum M,
    # without dust emission
    @staticmethod
    @numba.jit(nopython=True, fastmath=True, nogil=True)
    def attenuate_nodust(M, A, zfactors):
        # final result is
        # M * (atten ** ebv) * zfactors
        for i in range(len(M)):
            M[i] *= A[i] * zfactors[i]


    def build_stellar_continuum(self, templateflux, templatecoeff,
                                ebv, vdisp=None, dust_emission=True):

        """Build a stellar continuum model.

        Parameters
        ----------
        templateflux : :class:`numpy.ndarray` [npix, ntemplates]
            Rest-frame, native-resolution template spectra corresponding to
            `templatewave`.
        templatecoeff : :class:`numpy.ndarray` [ntemplates]
            Column vector of positive coefficients corresponding to each
            template.
        ebv : :class:`float`
            Dust extinction parameter, E(B-V), in mag.
        vdisp : :class:`float`
            Velocity dispersion in km/s. If `None`, do not convolve to the
            specified velocity dispersion (usually because `templateflux` has
            already been smoothed to some nominal value).
        dust_emission : :class:`bool`
            Model impact of infrared dust emission spectrum. Energy-balance is used
            to compute the normalization of this spectrum.

        Returns
        -------
        contmodel : :class:`numpy.ndarray` [npix]
            Full-wavelength, native-resolution, observed-frame model spectrum.

        """
        # [1] - Compute the weighted sum of the templates.
        contmodel = templateflux.dot(templatecoeff)

        # [2] - Optionally convolve to the desired velocity dispersion.
        if vdisp is not None:
            contmodel = Templates.convolve_vdisp(contmodel, vdisp,
                                                 pixsize_kms=Templates.PIXKMS_BLU,
                                                 limit=self.pixpos_wavesplit)

        # [3] - Apply dust attenuation; ToDo: allow age-dependent
        # attenuation. Also compute the bolometric luminosity before and after
        # attenuation but only if we have dustflux.

        # [4] - Optionally add the dust emission spectrum, conserving the total
        # (absorbed + re-emitted) energy. NB: (1) dustflux must be normalized to
        # unity already (see templates.py); and (2) we are ignoring dust
        # self-absorption, which should be mostly negligible.

        # [5] - Redshift factors.

        # do this part in Numpy because it is very slow
        # in Numba unless accelerated transcendentals are
        # available via, e.g., Intel SVML.
        A = self.lg_atten * ebv
        np.exp(A, out=A)

        if dust_emission:
            self.attenuate(contmodel, A,
                           self.zfactors,
                           self.templates.wave,
                           self.templates.dustflux)
        else:
            self.attenuate_nodust(contmodel, A,
                                  self.zfactors)

        return contmodel


    def continuum_to_spectroscopy(self, contmodel):
        """
        Synthesize spectroscopy from a continuum model.

        Parameters
        ----------
        contmodel : :class:`numpy.ndarray` [npix]
            Full-wavelength, native-resolution, observed-frame model spectrum.

        Returns
        -------
        modelflux : :class:`numpy.ndarray` [nwave]
            Observed-frame model spectrum at the instrumental resolution and
            wavelength sampling given by `specwave` and `specres`.

        """
        camerapix = self.data['camerapix']
        specwave  = self.data['wave']
        specres   = self.data['res']

        modelflux = np.empty(self.wavelen)

        for icam, (s, e) in enumerate(camerapix):
            resampflux = trapz_rebin(self.ztemplatewave,
                                     contmodel,
                                     specwave[icam],
                                     pre=self.spec_pre[icam])
            modelflux[s:e] = specres[icam] @ resampflux

        return modelflux


    def continuum_to_photometry(self, contmodel,
                                filters=None,
                                phottable=False,
                                get_abmag=False):
        """
        Synthesize photometry from a continuum model.

        Parameters
        ----------
        contmodel : :class:`numpy.ndarray` [npix]
            Full-wavelength, native-resolution, observed-frame model spectrum.
        filters : :class:`list` or :class:`speclite.filters.FilterSequence` [nfilt]
            Optionally override the filter curves stored in the `filters`
            attribute of the global Photometry object.
        phottable : :class:`bool`
            Return an :class:`astropy.table.Table` with additional bandpass information.
            Otherwise, return synthesized photometry in f_lambda (10^17 erg/s/cm2/A) units,
            which are used in fitting  Only true for QA.
        get_abmag : :class:`bool`
            Add AB magnitudes to the synthesized photometry table (only
            applies if `phottable=True`. Only used for QA.

        Returns
        -------
        modelphot : :class:`numpy.ndarray` or :class:`astropy.table.Table` [nfilt]
            Synthesized model photometry as an array or a table, depending on
            `phottable`.

        """
        if filters is None:
            filters, effwave, maggies_pre = self.phot_pre
        else:
            effwave = filters.effective_wavelengths.value
            maggies_pre = None

        modelmaggies = Photometry.get_ab_maggies_unchecked(filters,
                                                           contmodel,
                                                           self.ztemplatewave,
                                                           pre=maggies_pre)

        if not phottable:
            modelphot = Photometry.get_photflam(modelmaggies, effwave)
        else:
            modelmaggies /= FLUXNORM * self.massnorm
            modelphot = Photometry.parse_photometry(self.phot.bands, modelmaggies, effwave,
                                                    nanomaggies=False, get_abmag=get_abmag)
        return modelphot


    def _stellar_objective(self, params, templateflux,
                           dust_emission, fit_vdisp,
                           objflam, objflamistd,
                           specflux, specistd,
                           synthphot, synthspec):
        """Objective function for fitting a stellar continuum.

        """
        assert (synthphot or synthspec), "request for empty residuals!"

        if fit_vdisp:
            ebv, vdisp    = params[:2]
            templatecoeff = params[2:]
        else:
            ebv           = params[0]
            vdisp         = None
            templatecoeff = params[1:]

        # FIXME: write results directly into resid array instead of
        # allocating new
        fullmodel = self.build_stellar_continuum(
            templateflux, templatecoeff,
            ebv=ebv, vdisp=vdisp, dust_emission=dust_emission)

        # save the full model each time we compute the objective;
        # after optimization, the final full model will be
        # saved here. (And yes, it works even if we are using
        # finite-differencing; the last computation of the
        # objective occurs after the last computation of the
        # Jacobian in least_squares().)
        self.optimizer_saved_contmodel = fullmodel

        # Compute residuals versus provided spectroscopy
        # and/or photometry. Allocate a residual array
        # big enough to hold whatever we compute.
        spec_reslen = len(specflux) if synthspec else 0
        phot_reslen = len(objflam)  if synthphot else 0
        resid = np.empty(spec_reslen + phot_reslen)
        resid_split = spec_reslen

        if synthspec:
            modelflux = self.continuum_to_spectroscopy(fullmodel)
            spec_resid = resid[:resid_split]
            spec_resid[:]  = specflux
            spec_resid    -= modelflux
            spec_resid    *= specistd

        if synthphot:
            modelflam = self.continuum_to_photometry(fullmodel)
            phot_resid = resid[resid_split:]
            phot_resid[:] = objflam
            phot_resid   -= modelflam
            phot_resid   *= objflamistd

        return resid


    def fit_stellar_continuum(self, templateflux, fit_vdisp,
                              vdisp_guess=125., ebv_guess=0.05,
                              coeff_guess=None,
                              vdisp_bounds=(75., 500.), ebv_bounds=(0., 3.),
                              dust_emission=True,
                              objflam=None, objflamistd=None,
                              specflux=None, specistd=None,
                              synthphot=False, synthspec=False):
        """Fit a stellar continuum using bounded non-linear least-squares.

        Parameters
        ----------
        templateflux : :class:`numpy.ndarray` [npix, ntemplate]
            Grid of input (model) spectra.
        fit_vdisp : :class:`bool`
            If true, solve for the velocity dispersion;
            if false, use a nominal dispersion.
        vdisp_guess : :class:`float`
            Guess for scalar value of the velocity dispersion if fitting.
        ebv_guess : :class:`float`
            Guess scalar value of the dust attenuation.
        coeff_guess : :class:`numpy.ndarray` [ntemplates]
            Guess of the template coefficients.
        vdisp_bounds : :class:`tuple`
            Two-element list of minimum and maximum allowable values of the
            velocity dispersion; only used if `fit_vdisp=True`.
        ebv_bounds : :class:`tuple`
            Two-element list of minimum and maximum allowable values of the
            reddening, E(B-V).
        dust_emission : :class:`bool`
            Model impact of infrared dust emission spectrum. Energy-balance is used
            to compute the normalization of this spectrum.
        objflam: :class: `numpy.ndarray`
            Measured object photometry (used if fitting photometry)
        objflamistd: :class: `numpy.ndarray`
            Sqrt of inverse variance of objflam (used if fitting photometry)
        specflux : :class:`numpy.ndarray` [nwave]
            Observed-frame spectrum in 10**-17 erg/s/cm2/A corresponding to
            `specwave` (used if fitting spectroscopy)
        specfluxistd : :class:`numpy.ndarray` [nwave]
            Sqrt of inverse variance of the observed-frame spectrum
            (used if fitting spectroscopy)
        synthphot: :class:`bool`
            True iff fitting objective includes object photometry
        synthspec: :class:`bool`
            True iff fitting objective includes observed spectrum

        Returns
        -------
        ebv : :class:`float`
            Maximum-likelihood dust extinction parameter in mag.
        vdisp : :class:`float`
            Maximum-likelihood velocity dispersion in km/s, or
            nominal velocity if it was not fitted.
        templatecoeff : :class:`numpy.ndarray` [ntemplate]
            Column vector of maximum-likelihood template coefficients.
        resid: :class:`numpy.ndarray`
            Vector of residuals at final parameter values

        Notes
        -----
        This method supports several different fitting 'modes', depending on the
        whether vdisp is fitted or left nominal, and whether the caller requests
        fitting a model against spectroscopy, photometry, or both.

        In all cases, we solve for dust attenuation via the E(B-V) parameter and
        we also include IGM attenuation.

        """
        from scipy.optimize import least_squares

        ntemplates = templateflux.shape[1]

        # Unpack the input data to infer the fitting "mode" and the objective
        # function.
        farg = {
            'templateflux':  templateflux,
            'dust_emission': dust_emission,
            'fit_vdisp':     fit_vdisp,
            'objflam':       objflam,
            'objflamistd':   objflamistd,
            'specflux':      specflux,
            'specistd':      specistd,
            'synthphot':     synthphot,
            'synthspec':     synthspec,
        }

        if coeff_guess is None:
            coeff_guess = np.ones(ntemplates)
        else:
            if len(coeff_guess) != ntemplates:
                errmsg = 'Mismatching dimensions between coeff_guess and ntemplates!'
                log.critical(errmsg)
                raise ValueError(errmsg)

        coeff_bounds = (0., 1e5)

        if fit_vdisp == True:
            initial_guesses = np.array((ebv_guess, vdisp_guess))
            bounds = [ebv_bounds, vdisp_bounds]
        else:
            initial_guesses = np.array((ebv_guess,))
            bounds = [ebv_bounds]

        initial_guesses = np.concatenate((initial_guesses, coeff_guess))
        bounds = bounds + [coeff_bounds] * ntemplates
        #xscale = np.hstack(([0.1, 50.], np.ones(ntemplates) * 1e-1))

        # NB: `x_scale` has been set to `jac` here to help with the numerical
        # convergence. There may be faster ways, of course...
        fit_info = least_squares(self._stellar_objective, initial_guesses, kwargs=farg,
                                 bounds=tuple(zip(*bounds)), method='trf',
                                 tr_solver='lsmr', tr_options={'regularize': True},
                                 x_scale='jac', max_nfev=5000, ftol=1e-5, xtol=1e-10)#, verbose=2)
        bestparams = fit_info.x
        resid      = fit_info.fun

        if fit_vdisp:
            ebv, vdisp    = bestparams[:2]
            templatecoeff = bestparams[2:]
        else:
            ebv           = bestparams[0]
            templatecoeff = bestparams[1:]
            vdisp = self.templates.vdisp_nominal

        return ebv, vdisp, templatecoeff, resid


    def stellar_continuum_chi2(self, resid,
                               ncoeff, vdisp_fitted,
                               split = 0,
                               ndof_spec = 0, ndof_phot = 0):
        """Compute the reduced spectroscopic and/or photometric chi2.

        resid:
           Vector of residuals from least-squares fitting
        ncoeff:
           Number of template coefficients fitted
        vdisp_fitted:
           True if the velocity dispersion was fitted
        split:
           Boundary between initial spectroscopy elements of
           residual and final photometry elements
        ndof_spec:
           Number of spectroscopy degrees of freedom
        ndof_phot:
           Number of photometry degrees of freedom

        """
        # ebv is always a free parameter
        nfree = ncoeff + 1 + int(vdisp_fitted)

        def _get_rchi2(chi2, ndof, nfree):
            """Guard against ndof=nfree."""
            if ndof > nfree:
                return chi2 / (ndof - nfree)
            else:
                return chi2 / ndof # ???

        if ndof_spec > 0:
            resid_spec = resid[:split]
            chi2_spec = resid_spec.dot(resid_spec)
            rchi2_spec = _get_rchi2(chi2_spec, ndof_spec, nfree)
        else:
            chi2_spec = 0.
            rchi2_spec = 0.

        if ndof_phot > 0:
            resid_phot = resid[split:]
            chi2_phot = resid_phot.dot(resid_phot)
            rchi2_phot = _get_rchi2(chi2_phot, ndof_phot, nfree)
        else:
            chi2_phot = 0.
            rchi2_phot = 0.

        rchi2_tot = _get_rchi2(chi2_spec + chi2_phot, ndof_spec + ndof_phot, nfree)

        return rchi2_spec, rchi2_phot, rchi2_tot


def continuum_specfit(data, result, templates, igm, phot,
                      constrain_age=False, no_smooth_continuum=False,
                      fastphot=False, debug_plots=False):
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
    def younger_than_universe(age, tuniv, agepad=0.5):
        """Return the indices of the templates younger than the age of the universe
        (plus an agepadding amount) at the given redshift. age in yr, agepad and
        tuniv in Gyr
        """
        return np.where(age <= 1e9 * (agepad + tuniv))[0]


    tall = time.time()

    CTools = ContinuumTools(igm, phot, templates, data, fastphot)

    redshift = data['zredrock']
    if redshift <= 0.:
        log.warning('Input redshift not defined, zero, or negative!')

    photometry = data['phot']
    objflam = photometry['flam'].value * FLUXNORM
    objflamivar = (photometry['flam_ivar'].value / FLUXNORM**2) * phot.bands_to_fit
    assert(np.all(objflamivar >= 0.))

    if np.any(phot.bands_to_fit):
        # Require at least one photometric optical band; do not just fit the IR
        # because we will not be able to compute the aperture correction, below.
        lambda_eff = photometry['lambda_eff'].value
        opt = ((lambda_eff > 3e3) & (lambda_eff < 1e4))
        if np.all(objflamivar[opt] == 0.):
            log.warning('All optical bands are masked; masking all photometry.')
            objflamivar[:] = 0.0

    objflamistd = np.sqrt(objflamivar)

    # Optionally ignore templates which are older than the age of the
    # universe at the redshift of the object.
    if constrain_age:
        agekeep = younger_than_universe(templates.info['age'].value, data['tuniv'])
    else:
        agekeep = np.arange(templates.ntemplates)
    nage = len(agekeep)

    vdisp_nominal = templates.vdisp_nominal
    ebv_guess = 0.05

    # Photometry-only fitting.
    if fastphot:
        log.info(f'Adopting nominal vdisp={vdisp_nominal:.0f} km/s.')
        vdisp = vdisp_nominal

        if np.all(objflamivar == 0.):
            log.info('All photometry is masked.')
            coeff = np.zeros(nage) # nage not nsed
            rchi2_cont, rchi2_phot = 0., 0.
            dn4000_model = 0.
            sedmodel = np.zeros(len(templates.wave))
        else:
            # Get the coefficients and chi2 at the nominal velocity dispersion.
            t0 = time.time()

            # maintain backwards-compatibility
            if templates.use_legacy_fitting:
                sedtemplates, sedphot_flam = CTools.templates2data(
                    templates.flux_nomvdisp[:, agekeep],
                    templates.wave, flamphot=True,
                    redshift=redshift, dluminosity=data['dluminosity'],
                    vdisp=None, synthphot=True, photsys=data['photsys'])
                sedflam = sedphot_flam * CTools.massnorm * FLUXNORM
                coeff, rchi2_phot = CTools.call_nnls(sedflam, objflam, objflamivar)
                rchi2_phot /= np.sum(objflamivar > 0.) # dof???
            else:
                ebv, _, coeff, resid = CTools.fit_stellar_continuum(
                    templates.flux_nomvdisp[:, agekeep], # [npix,nsed]
                    fit_vdisp=False,
                    vdisp_guess=vdisp_nominal, ebv_guess=ebv_guess,
                    objflam=objflam, objflamistd=objflamistd,
                    synthphot=True, synthspec=False
                )

                _, rchi2_phot, rchi2_cont = CTools.stellar_continuum_chi2(
                    resid, ncoeff=len(coeff), vdisp_fitted=False,
                    ndof_phot=np.sum(objflamivar > 0.)
                )

                sedmodel = CTools.optimizer_saved_contmodel

            log.info(f'Fitting {nage} models took {time.time()-t0:.2f} seconds.')

            if np.all(coeff == 0.):
                log.warning('Continuum coefficients are all zero.')
                sedmodel = np.zeros(len(templates.wave))
                dn4000_model = 0.
            else:
                # Measure Dn(4000) from the line-free model.
                if templates.use_legacy_fitting:
                    sedmodel = sedtemplates.dot(coeff)
                    sedtemplates_nolines, _ = CTools.templates2data(
                        templates.flux_nolines_nomvdisp[:, agekeep],
                        templates.wave,
                        redshift=redshift, dluminosity=data['dluminosity'],
                        vdisp=None, synthphot=False)
                    sedmodel_nolines = sedtemplates_nolines.dot(coeff)
                else:
                    sedmodel_nolines = CTools.build_stellar_continuum(
                        templates.flux_nolines_nomvdisp[:, agekeep], coeff,
                        ebv=ebv, vdisp=None,
                        dust_emission=False,
                    )

                    log.info(f'Best-fitting E(B-V)={ebv:.3f} mag.')

                dn4000_model, _ = Photometry.get_dn4000(templates.wave,
                                                        sedmodel_nolines, rest=True)
                log.info(f'Model Dn(4000)={dn4000_model:.3f}.')
    else:
        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        specwave = np.hstack(data['wave'])
        specflux = np.hstack(data['flux'])
        flamivar = np.hstack(data['ivar'])
        specivar = flamivar * ~np.hstack(data['linemask']) # mask emission lines

        if np.all(specivar == 0.) or np.any(specivar < 0.):
            specivar = flamivar # not great...
            if np.all(specivar == 0.) or np.any(specivar < 0.):
                errmsg = 'All pixels are masked or some inverse variances are negative!'
                log.critical(errmsg)
                raise ValueError(errmsg)

        specistd = np.sqrt(specivar)

        npix = len(specwave)

        # We'll need the filters for the aperture correction, below.
        filters_in = phot.synth_filters[data['photsys']]

        # Solve for the velocity dispersion if the wavelength coverage is
        # sufficient.
        restwave = specwave / (1. + redshift)
        Ivdisp = np.where((specivar > 0) & (restwave > 3500.) & (restwave < 5500.))[0]
        compute_vdisp = (len(Ivdisp) > 0) and (np.ptp(restwave[Ivdisp]) > 500.)

        if len(data['cameras']) == 3:
            log.info('S/N_{}={:.2f}, S/N_{}={:.2f}, S/N_{}={:.2f}, rest wavelength coverage={:.0f}-{:.0f} A.'.format(
                data['cameras'][0], data['snr'][0],
                data['cameras'][1], data['snr'][1],
                data['cameras'][2], data['snr'][2],
                restwave[0], restwave[-1]))

        # Maintain backwards compatibility. With the old templates, the velocity
        # dispersion and aperture corrections are determined separately, so we
        # separate that code out from the new templates, where they are
        # determined simultatneously.
        if templates.use_legacy_fitting:
            if compute_vdisp:
                t0 = time.time()
                ztemplateflux_vdisp, _ = CTools.templates2data(
                    templates.vdispflux, templates.vdispwave, # [npix,vdispnsed,nvdisp]
                    redshift=redshift, dluminosity=data['dluminosity'],
                    specwave=data['wave'], specres=data['res'],
                    cameras=data['cameras'], synthphot=False, stack_cameras=True)
                vdispchi2min, vdispbest, vdispivar, _ = CTools.call_nnls(
                    ztemplateflux_vdisp[Ivdisp, :, :],
                    specflux[Ivdisp], specivar[Ivdisp],
                    xparam=templates.vdisp, xlabel=r'$\sigma$ (km/s)',
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
                input_templateflux         = templates.flux_nomvdisp[:, agekeep]
                input_templateflux_nolines = templates.flux_nolines_nomvdisp[:, agekeep]
            else:
                use_vdisp = vdispbest
                input_templateflux         = templates.flux[:, agekeep]
                input_templateflux_nolines = templates.flux_nolines[:, agekeep]

            desitemplates, desitemplatephot_flam = CTools.templates2data(
                input_templateflux, templates.wave,
                redshift=redshift, dluminosity=data['dluminosity'],
                specwave=data['wave'], specres=data['res'], specmask=data['mask'],
                vdisp=use_vdisp, cameras=data['cameras'], stack_cameras=True,
                synthphot=True, flamphot=True, photsys=data['photsys'])
            desitemplateflam = desitemplatephot_flam * CTools.massnorm * FLUXNORM

            apercorrs, apercorr = np.zeros(len(phot.synth_bands)), 0.

            sedtemplates, _ = CTools.templates2data(
                input_templateflux, templates.wave,
                vdisp=use_vdisp, redshift=redshift,
                dluminosity=data['dluminosity'], synthphot=False)

            if not np.any(phot.bands_to_fit):
                log.info('Skipping aperture correction since no bands were fit.')
                apercorrs, apercorr = np.ones(len(phot.synth_bands)), 1.
            else:
                # Fit using the templates with line-emission.
                quickcoeff, _ = CTools.call_nnls(desitemplates, specflux, specivar)
                if np.all(quickcoeff == 0.):
                    log.warning('Quick continuum coefficients are all zero.')
                else:
                    ztemplatewave = CTools.ztemplatewave

                    # Synthesize grz photometry from the full-wavelength SED to make
                    # sure we get the z-band correct.
                    nanomaggies = photometry['nanomaggies'].value
                    numer = np.hstack([nanomaggies[photometry['band'] == band] for band in phot.synth_bands])

                    quicksedflux = sedtemplates.dot(quickcoeff)
                    quickmaggies = Photometry.get_ab_maggies(filters_in,
                                                             quicksedflux / FLUXNORM,
                                                             ztemplatewave)

                    denom = Photometry.to_nanomaggies(quickmaggies)

                    I = ((numer > 0.) & (denom > 0.))
                    if np.any(I):
                        apercorrs[I] = numer[I] / denom[I]
                    I = (apercorrs > 0.)
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
                input_templateflux_nolines, templates.wave, redshift=redshift,
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
                sedmodel = np.zeros(len(templates.wave))
                desimodel = np.zeros_like(specflux)
                desimodel_nolines = np.zeros_like(specflux)
                dn4000_model = 0.0
            else:
                sedmodel = sedtemplates.dot(coeff)
                desimodel = desitemplates.dot(coeff)
                desimodel_nolines = desitemplates_nolines.dot(coeff)

                # Measure Dn(4000) from the line-free model.
                sedtemplates_nolines, _ = CTools.templates2data(
                    input_templateflux_nolines, templates.wave,
                    vdisp=use_vdisp, redshift=redshift, dluminosity=data['dluminosity'],
                    synthphot=False)
                sedmodel_nolines = sedtemplates_nolines.dot(coeff)

                dn4000_model, _ = Photometry.get_dn4000(templates.wave,
                                                        sedmodel_nolines, rest=True)

            if use_vdisp is None:
                vdisp = vdisp_nominal
            else:
                vdisp = use_vdisp

        else: # new templates start here

            # First, estimate the aperture correction from a (noiseless) *model*
            # of the spectrum (using the nominal velocity dispersion).
            apercorrs, apercorr = np.ones(len(phot.synth_bands)), 1.

            if np.any(phot.bands_to_fit):
                t0 = time.time()

                ebv, _, coeff, _ = CTools.fit_stellar_continuum(
                    templates.flux_nomvdisp[:, agekeep], # [npix,nsed]
                    dust_emission=False,  fit_vdisp=False,
                    vdisp_guess=None, ebv_guess=ebv_guess,
                    specflux=specflux, specistd=specistd,
                    synthphot=False, synthspec=True
                )

                if np.all(coeff == 0.):
                    log.warning('Unable to estimate aperture correction because continuum coefficients are all zero; adopting 1.0.')
                    coeff_guess = np.ones(nage)
                else:
                    sedflam = CTools.continuum_to_photometry(
                        CTools.optimizer_saved_contmodel,
                        filters=filters_in)

                    objflam_aper = FLUXNORM * photometry[np.isin(photometry['band'], phot.synth_bands)]['flam'].value

                    I = ((objflam_aper > 0.) & (sedflam > 0.))
                    if np.any(I):
                        apercorrs[I] = objflam_aper[I] / sedflam[I]

                    I = (apercorrs > 0.)
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

            if not compute_vdisp:
                # Use the cached templates with nominal velocity dispersion
                input_templateflux         = templates.flux_nomvdisp[:, agekeep]
                input_templateflux_nolines = templates.flux_nolines_nomvdisp[:, agekeep]
            else:
                input_templateflux         = templates.flux[:, agekeep]
                input_templateflux_nolines = templates.flux_nolines[:, agekeep]

            ebv, vdisp, coeff, resid = CTools.fit_stellar_continuum(
                input_templateflux, # [npix,nage]
                fit_vdisp=compute_vdisp,
                ebv_guess=ebv, vdisp_guess=vdisp_nominal,
                coeff_guess=coeff_guess,
                objflam=objflam, objflamistd=objflamistd,
                specflux=specflux*apercorr, specistd=specistd/apercorr,
                synthphot=True, synthspec=True)

            log.info(f'Final spectrophotometric fitting with {nage} models took {time.time()-t0:.2f} seconds.')

            # ToDo: capture case where vdisp (and also maybe ebv) hits its bounds.
            if np.all(coeff == 0.):
                log.warning('Continuum coefficients are all zero.')
                sedmodel = np.zeros(len(templates.wave))
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
                    log.info(f'Insufficient wavelength covereage to compute vdisp; adopting nominal vdisp={vdisp:.1f} km/s')

                _, rchi2_phot, rchi2_cont = CTools.stellar_continuum_chi2(
                    resid, ncoeff=len(coeff), vdisp_fitted=compute_vdisp,
                    split=len(specflux),
                    ndof_spec=np.sum(specivar > 0.),
                    ndof_phot=np.sum(objflamivar > 0.)
                )

                # get the best-fitting model with and without line-emission
                sedmodel = CTools.optimizer_saved_contmodel

                sedmodel_nolines = CTools.build_stellar_continuum(
                    input_templateflux_nolines, coeff,
                    ebv=ebv, vdisp=vdisp, dust_emission=False
                )
                desimodel_nolines = CTools.continuum_to_spectroscopy(sedmodel_nolines)

                dn4000_model, _ = Photometry.get_dn4000(templates.wave,
                                                        sedmodel_nolines, rest=True)

        # Get DN(4000). Specivar is line-masked so we can't use it!
        dn4000, dn4000_ivar = Photometry.get_dn4000(specwave, specflux, flam_ivar=flamivar,
                                                    redshift=redshift, rest=False)

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
            residuals = specflux*apercorr - desimodel_nolines
            I = ((specflux == 0.) & (specivar == 0.))
            residuals[I] = 0.

            if debug_plots:
                png = f'qa-smooth-continuum-{data["uniqueid"]}.png'
            else:
                png = None

            linemask = np.hstack(data['linemask'])
            _smoothcontinuum = CTools.smooth_continuum(
                specwave, residuals, specivar / apercorr**2, linemask,
                camerapix=data['camerapix'], png=png)

        # Unpack the continuum into individual cameras.
        continuummodel, smoothcontinuum = [], []
        for campix in data['camerapix']:
            s, e = campix
            continuummodel.append(desimodel_nolines[s:e])
            smoothcontinuum.append(_smoothcontinuum[s:e])

        for icam, cam in enumerate(data['cameras']):
            nonzero = (continuummodel[icam] != 0)
            if np.any(nonzero):
                corr = median(smoothcontinuum[icam][nonzero] / continuummodel[icam][nonzero])
                result[f'SMOOTHCORR_{cam.upper()}'] = corr * 100 # [%]

        if len(data['cameras']) == 3:
            log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
                result['SMOOTHCORR_B'], result['SMOOTHCORR_R'], result['SMOOTHCORR_Z']))

    result['COEFF'][agekeep] = coeff
    result['RCHI2_PHOT'] = rchi2_phot
    result['VDISP'] = vdisp # * u.kilometer/u.second

    if not fastphot:
        result['RCHI2_CONT'] = rchi2_cont
        result['VDISP_IVAR'] = vdispivar # * (u.second/u.kilometer)**2

        result['APERCORR'] = apercorr
        for iband, band in enumerate(phot.synth_bands):
            result[f'APERCORR_{band.upper()}'] = apercorrs[iband]
        result['DN4000_OBS'] = dn4000
        result['DN4000_IVAR'] = dn4000_ivar
        result['DN4000_MODEL'] = dn4000_model

    # Compute K-corrections, rest-frame quantities, and physical properties.
    if not np.all(coeff == 0):
        kcorr, absmag, ivarabsmag, synth_bestmaggies = phot.kcorr_and_absmag(
            photometry['nanomaggies'].value, photometry['nanomaggies_ivar'].value,
            redshift, data['dmodulus'], data['photsys'],
            CTools.ztemplatewave, sedmodel
        )
        lums, cfluxes = CTools.continuum_fluxes(sedmodel)

        for iband, (band, shift) in enumerate(zip(phot.absmag_bands, phot.band_shift)):
            band = band.upper()
            shift = int(10*shift)
            result[f'KCORR{shift:02d}_{band}'] = kcorr[iband] # * u.mag
            result[f'ABSMAG{shift:02d}_{band}'] = absmag[iband] # * u.mag
            result[f'ABSMAG{shift:02d}_IVAR_{band}'] = ivarabsmag[iband] # / (u.mag**2)
            for iband, band in enumerate(phot.bands):
                result[f'FLUX_SYNTH_PHOTMODEL_{band.upper()}'] = 1e9 * synth_bestmaggies[iband] # * u.nanomaggy
            for lum in lums:
                result[lum] = lums[lum]
            for cflux in cfluxes:
                result[cflux] = cfluxes[cflux]

        # get the SPS properties
        tinfo = templates.info[agekeep]
        mstars = tinfo['mstar'] # [current mass in stars, Msun]
        masstot = coeff.dot(mstars)
        coefftot = np.sum(coeff)
        logmstar = np.log10(CTools.massnorm * masstot)
        zzsun = np.log10(coeff.dot(mstars * 10.**tinfo['zzsun']) / masstot) # mass-weighted
        age = coeff.dot(tinfo['age']) / coefftot / 1e9           # luminosity-weighted [Gyr]
        #age = coeff.dot(mstars * tinfo['age']) / masstot / 1e9        # mass-weighted [Gyr]
        sfr = CTools.massnorm * coeff.dot(tinfo['sfr'])                          # [Msun/yr]
        if templates.use_legacy_fitting:
            AV = coeff.dot(tinfo['av']) / coefftot # luminosity-weighted [mag]
        else:
            AV = ebv * Templates.klambda(5500.) # [mag]

        result['AV'] = AV # * u.mag
        result['AGE'] = age # * u.Gyr
        result['ZZSUN'] = zzsun
        result['LOGMSTAR'] = logmstar
        result['SFR'] = sfr

        rindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1.+phot.band_shift) - 5600))
        log.info(f'log(M/Msun)={logmstar:.2f}, M{phot.absmag_bands[rindx]}={absmag[rindx]:.2f} mag, ' + \
                 f'A(V)={AV:.3f}, Age={age:.3f} Gyr, SFR={sfr:.3f} Msun/yr, Z/Zsun={zzsun:.3f}')

    log.info(f'Continuum-fitting took {time.time()-tall:.2f} seconds.')

    if fastphot:
        return sedmodel, None
    else:
        # divide out the aperture correction
        continuummodel  = [cm / apercorr for cm in continuummodel ]
        smoothcontinuum = [sc / apercorr for sc in smoothcontinuum]
        return continuummodel, smoothcontinuum
