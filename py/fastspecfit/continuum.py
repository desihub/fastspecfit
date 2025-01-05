"""
fastspecfit.continuum
=====================
Methods and tools for continuum-fitting.

"""
import time
import numpy as np
from numba import jit

from fastspecfit.logger import log
from fastspecfit.photometry import Photometry
from fastspecfit.templates import Templates
from fastspecfit.util import (
    C_LIGHT, TINY, F32MAX, quantile, median, var2ivar,
    trapz_rebin, trapz_rebin_pre)


class ContinuumTools(object):
    """Tools for dealing with spectral continua.

    """
    def __init__(self, data, templates, phot, igm, tauv_guess=0.1,
                 vdisp_guess=250., tauv_bounds=(0., 2.),
                 vdisp_bounds=(75., 500.), vdisp_nbin=5,
                 fluxnorm=1e17, massnorm=1e10, fastphot=False,
                 constrain_age=False):

        self.phot = phot
        self.templates = templates
        self.data = data

        self.massnorm = massnorm  # stellar mass normalization factor [Msun]
        self.fluxnorm = fluxnorm  # flux normalization factor [erg/s/cm2/A]

        self.tauv_guess = tauv_guess
        self.vdisp_guess = vdisp_guess
        self.tauv_bounds = tauv_bounds
        self.vdisp_bounds = vdisp_bounds
        self.vdisp_grid = np.linspace(vdisp_bounds[0], vdisp_bounds[1], vdisp_nbin)

        # Cache the redshift-dependent factors (incl. IGM attenuation),
        redshift = data['redshift']
        self.ztemplatewave = templates.wave * (1. + redshift)
        self.zfactors = self.get_zfactors(
            igm, self.ztemplatewave, redshift=redshift,
            dluminosity=data['dluminosity'])

        # Optionally ignore templates which are older than the age of the
        # universe at the redshift of the object.
        if constrain_age:
            self.agekeep = self._younger_than_universe(templates.info['age'].value, data['tuniv'])
            self.nage = len(self.agekeep)
        else:
            # use default slice instead of arange to avoid copying templates
            self.agekeep = slice(None, None)
            self.nage = templates.ntemplates

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
            rspre = [ trapz_rebin_pre(w) for w in self.data['wave'] ]
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
        T *= self.fluxnorm * self.massnorm * (10. / (1e6 * dluminosity))**2 / (1. + redshift)

        return T

    @staticmethod
    def _younger_than_universe(age, tuniv, agepad=0.5):
        """Return the indices of the templates younger than the age of the universe
        (plus an agepadding amount) at the given redshift. age in yr, agepad and
        tuniv in Gyr

        """
        return np.where(age <= 1e9 * (agepad + tuniv))[0]


    @staticmethod
    def smooth_continuum(wave, flux, ivar, linemask, camerapix,
                         uniqueid=0, smooth_window=75, smooth_step=125,
                         clip_sigma=2., nminpix=15, nmaskpix=9,
                         debug_plots=False):
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
        smooth_flux :class:`numpy.ndarray` [npix]
        Smooth continuum spectrum which can be subtracted from `flux` in
        order to create a pure emission-line spectrum.

        """
        from numpy.lib.stride_tricks import sliding_window_view
        from scipy.ndimage import median_filter
        from scipy.interpolate import UnivariateSpline
        from fastspecfit.util import sigmaclip

        npix = len(wave)
        if len(linemask) != npix:
            errmsg = 'Linemask must have the same number of pixels as the input spectrum.'
            log.critical(errmsg)
            raise ValueError(errmsg)


        def _smooth_percamera(camwave, camflux, camivar, camlinemask):

            # Mask nmaskpix (presumably noisy) pixels from the edge
            # of each per-camera spectrum.
            cammask = (camlinemask | (camivar <= 0.))
            cammask[:nmaskpix] = True
            cammask[-nmaskpix:] = True

            # Build the smooth (line-free) continuum by computing statistics in a
            # sliding window, accounting for masked pixels and trying to be smart
            # about broad lines. See:
            #   https://stackoverflow.com/questions/41851044/python-median-filter-for-1d-numpy-array
            #   https://numpy.org/devdocs/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html
            wave_win = sliding_window_view(camwave, window_shape=smooth_window)
            flux_win = sliding_window_view(camflux, window_shape=smooth_window)
            ivar_win = sliding_window_view(camivar, window_shape=smooth_window)
            nomask_win = sliding_window_view(np.logical_not(cammask), window_shape=smooth_window)

            swave, sflux, sisig = [], [], []
            for wwave, wflux, wivar, wnomask in zip(
                    wave_win[::smooth_step], flux_win[::smooth_step],
                    ivar_win[::smooth_step], nomask_win[::smooth_step]):

                # If there are fewer than nminpix good pixels after all
                # masking, discard the window.
                umflux = wflux[wnomask]
                if len(umflux) < nminpix:
                    continue

                cflux, _ = sigmaclip(umflux, low=clip_sigma, high=clip_sigma)
                if len(cflux) < nminpix:
                    continue

                mn, clo, chi = quantile(cflux, (0.5, 0.25, 0.75)) # robust stats
                sig = (chi - clo) / 1.349 # robust sigma

                # One more check for crummy spectral regions.
                if mn == 0. or sig <= 0.:
                    continue

                umwave = wwave[wnomask]
                swave.append(np.mean(umwave))
                sflux.append(mn)
                sisig.append(1. / sig) # inverse sigma

            swave = np.array(swave)
            sflux = np.array(sflux)
            sisig = np.array(sisig)

            # corner case for very wacky spectra
            if len(sflux) == 0:
                smoothflux = np.zeros_like(camflux)
            else:
                ## remove duplicate wavelength values, which should never
                ## happen...
                #_, uindx = np.unique(swave, return_index=True)
                #swave = swave[uindx]
                #sflux = sflux[uindx]
                #sisig = sisig[uindx]

                # We supply estimates local inverse stddev in each window
                # (i.e., how noisy the data is there) so that variation is
                # down-weighted in noisier regions. Note: ext=3 means constant
                # extrapolation.
                if len(swave) > 3:
                    spl_flux = UnivariateSpline(swave, sflux, w=sisig, ext=3, k=2)
                    smoothflux = spl_flux(camwave)
                else:
                    smoothflux = np.zeros_like(camflux)

                # evaluate on the original wavelength vector

            # very important!
            smoothflux[(camflux == 0.) & (camivar == 0.)] = 0.

            return swave, sflux, smoothflux

        smooth_wave, smooth_flux, smoothcontinuum = [], [], []
        for ss, ee in camerapix:
            smooth_wave1, smooth_flux1, smoothcontinuum1 = _smooth_percamera(
                wave[ss:ee], flux[ss:ee], ivar[ss:ee], linemask[ss:ee])
            smooth_wave.append(smooth_wave1)
            smooth_flux.append(smooth_flux1)
            smoothcontinuum.append(smoothcontinuum1)
        smooth_wave = np.hstack(smooth_wave)
        smooth_flux = np.hstack(smooth_flux)
        smoothcontinuum = np.hstack(smoothcontinuum)

        # Optional QA.
        if debug_plots:
            import numpy.ma as ma
            import matplotlib.pyplot as plt
            import seaborn as sns

            pngfile = f'qa-smooth-continuum-{uniqueid}.png'
            sns.set(context='talk', style='ticks', font_scale=0.7)

            srt = np.argsort(wave)

            resid = flux - smoothcontinuum
            noise = np.ptp(quantile(resid[~linemask], (0.25, 0.75))) / 1.349 # robust sigma

            msk = ma.array(linemask)
            msk.mask = linemask
            clumps_masked = ma.clump_masked(msk)
            clumps_unmasked = ma.clump_unmasked(msk)

            fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
            for iclump, clump in enumerate(clumps_unmasked):
                if iclump == 0:
                    label = 'Unmasked Flux'
                else:
                    label = None
                ax[0].plot(wave[srt][clump] / 1e4, flux[srt][clump], color='grey',
                           alpha=0.5, lw=0.5, label=label)
            for iclump, clump in enumerate(clumps_masked):
                if iclump == 0:
                    label = 'Masked Flux'
                else:
                    label = None
                ax[0].plot(wave[srt][clump] / 1e4, flux[srt][clump], alpha=0.3, lw=0.5,
                           color='blue', label=label)
            ax[0].scatter(smooth_wave / 1e4, smooth_flux, edgecolor='k', color='orange',
                          marker='s', alpha=0.8, s=20, zorder=3,
                          label='Smooth Data')
            for icam, (ss, ee) in enumerate(camerapix):
                if icam == 0:
                    label = 'Smooth Model'
                else:
                    label = None
                ax[0].plot(wave[ss:ee] / 1e4, smoothcontinuum[ss:ee], color='red',
                           zorder=4, ls='-', lw=2, alpha=0.6, label=label)

            ax[0].set_ylim(np.min((-5. * noise, quantile(flux, 0.05))),
                           np.max((5. * noise, 1.5 * quantile(flux, 0.975))))
            ax[0].set_ylabel('Continuum-subtracted Flux\n' + \
                             r'$(10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$')
            leg = ax[0].legend(fontsize=10, loc='upper left')
            for line in leg.get_lines():
                line.set_linewidth(2)

            #ax[1].plot(wave[srt] / 1e4, resid[srt], alpha=0.75, lw=0.5)
            for clump in clumps_unmasked:
                ax[1].plot(wave[srt][clump] / 1e4, resid[srt][clump], color='grey',
                           alpha=0.5, lw=0.5)
            for clump in clumps_masked:
                ax[1].plot(wave[srt][clump] / 1e4, resid[srt][clump], alpha=0.3,
                           lw=0.5, color='blue')
            ax[1].axhline(y=0, color='k', ls='--', lw=2)
            ax[1].set_ylim(np.min((-5. * noise, quantile(resid, 0.05))),
                           np.max((5. * noise, 1.5 * quantile(resid, 0.975))))
            ax[1].set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
            ax[1].set_ylabel('Residual Flux\n' + r'$(10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$')
            #ax[1].legend(fontsize=10)

            ax[0].set_title(f'Smooth Continuum: {uniqueid}')
            fig.tight_layout()
            fig.savefig(pngfile)#, bbox_inches='tight')
            plt.close()
            log.info(f'Wrote {pngfile}')

        return smoothcontinuum


    @staticmethod
    def lums_keys():
        """Simple method which defines the rest-frame luminosity keys and
        wavelengths.

        """
        keys = ('LOGL_1450', 'LOGLNU_1500', 'LOGL_1700',
                'LOGLNU_2800', 'LOGL_3000', 'LOGL_5100')
        waves = (1450., 1500., 1700., 2800., 3000., 5100.)
        return keys, waves


    @staticmethod
    def cfluxes_keys():
        """Simple method which defines the observed-frame continuum flux keys
        and wavelengths.

        """
        keys = ('FLYA_1215_CONT', 'FOII_3727_CONT', 'FHBETA_CONT',
                'FOIII_5007_CONT', 'FHALPHA_CONT')
        waves = (1215.67, 3728.48, 4862.71, 5008.24, 6564.6)
        return keys, waves


    def continuum_fluxes(self, continuum, uniqueid=0, width1=50., width2=100.,
                         debug_plots=False):
        """Compute rest-frame luminosities and observed-frame continuum fluxes.

        """
        from fastspecfit.util import sigmaclip

        def _get_cflux(cwave, linear_fit=False, siglo=2., sighi=2.,
                       ignore_core=False, return_slope=False):
            # continuum flux in 10**-17 erg/s/cm2/A

            lo = np.searchsorted(templatewave, cwave - width1, 'right')
            hi = np.searchsorted(templatewave, cwave + width1, 'left')
            lo2 = np.searchsorted(templatewave, cwave - width2, 'right')
            hi2 = np.searchsorted(templatewave, cwave + width2, 'left')

            if ignore_core:
                ylo, lomask = sigmaclip(continuum[lo2:lo], low=siglo, high=sighi)
                yhi, himask = sigmaclip(continuum[hi:hi2], low=siglo, high=sighi)
                xlo = templatewave[lo2:lo][lomask]
                xhi = templatewave[hi:hi2][himask]
                xfit = np.hstack((xlo, xhi))
                yfit = np.hstack((ylo, yhi))
            else:
                yfit, mask = sigmaclip(continuum[lo2:hi2], low=siglo, high=sighi)
                xfit = templatewave[lo2:hi2][mask]

            if linear_fit:
                slope, cflux = np.polyfit(xfit - cwave, yfit, 1)
            else:
                slope = None
                cflux = median(yfit)

            if return_slope:
                return cflux, slope
            else:
                return cflux

        llabels, lcwaves = self.lums_keys()
        flabels, fcwaves = self.cfluxes_keys()

        lums = np.zeros(len(lcwaves))
        cfluxes = np.zeros(len(fcwaves))

        redshift = self.data['redshift']
        if redshift <= 0.0:
            log.warning('Input redshift not defined, zero, or negative!')
            return lums, cfluxes

        templatewave = self.templates.wave

        # compute the model continuum flux at 1500 and 2800 A (to facilitate UV
        # luminosity-based SFRs) and at the positions of strong nebular emission
        # lines [OII], Hbeta, [OIII], and Halpha
        dlum = self.data['dluminosity']
        dfactor = (1. + redshift) * 4. * np.pi * (3.08567758e24 * dlum)**2 / self.fluxnorm

        for ilum, (cwave, label) in enumerate(zip(lcwaves, llabels)):
            cflux = _get_cflux(cwave, linear_fit=True) * dfactor # [monochromatic luminosity in erg/s/A]

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
                lums[ilum] = np.log10(cflux) # * u.erg/(u.second*u.Hz)

        for icflux, (cwave, label) in enumerate(zip(fcwaves, flabels)):
            if 'FLYA' in label or 'FHBETA' in label or 'FHALPHA' in label:
                ignore_core = True
            else:
                ignore_core = False
            cfluxes[icflux] = _get_cflux(cwave, linear_fit=True, ignore_core=ignore_core)

        # simple QA
        if debug_plots:
            import matplotlib.pyplot as plt
            import seaborn as sns

            pngfile = f'qa-continuum-fluxes-{uniqueid}.png'
            sns.set(context='talk', style='ticks', font_scale=0.6)

            templatewave = self.templates.wave

            labels = np.hstack((llabels, flabels))
            cwaves = np.hstack((lcwaves, fcwaves))
            #linear_fits = np.hstack(([True] * len(lcwaves), [False] * len(fcwaves)))
            linear_fits = [True] * len(labels)
            ncwaves = len(cwaves)

            ncols = 3
            nrows = int(np.ceil(ncwaves / ncols))

            fig, ax = plt.subplots(nrows, ncols, figsize=(3*ncols, 2*nrows))
            for iwave, (cwave, label, linear_fit, xx) in enumerate(zip(cwaves, labels, linear_fits, ax.flat)):
                lo = np.searchsorted(templatewave, cwave - width1, 'right')
                hi = np.searchsorted(templatewave, cwave + width1, 'left')
                lo2 = np.searchsorted(templatewave, cwave - width2, 'right')
                hi2 = np.searchsorted(templatewave, cwave + width2, 'left')
                lo3 = np.searchsorted(templatewave, cwave - 300., 'right')
                hi3 = np.searchsorted(templatewave, cwave + 300., 'left')

                if 'FLYA' in label or 'FHBETA' in label or 'FHALPHA' in label:
                    ignore_core = True
                else:
                    ignore_core = False
                cflux, slope = _get_cflux(cwave, linear_fit=linear_fit, return_slope=True,
                                          ignore_core=ignore_core)

                xx.plot(templatewave[lo3:hi3] / 1e4, continuum[lo3:hi3])
                if slope is not None:
                    xx.plot(templatewave[lo2:hi2] / 1e4, cflux + slope * (templatewave[lo2:hi2] - cwave),
                            color='k', lw=2, ls='-')
                xx.axvline(x=(cwave - width1)  / 1e4, lw=1, ls='--', color='k')
                xx.axvline(x=(cwave + width1) / 1e4, lw=1, ls='--', color='k')
                xx.axvline(x=(cwave - width2)  / 1e4, lw=1, ls='-', color='k')
                xx.axvline(x=(cwave + width2) / 1e4, lw=1, ls='-', color='k')
                xx.axhline(y=cflux, lw=2, ls='-', color='k')
                xx.axvline(x=cwave, lw=2, ls='-', color='k')
                xx.scatter([cwave / 1e4] * 2, [cflux] * 2, zorder=10, marker='s',
                           color='red', edgecolor='k', s=70)
                xx.set_xlim(templatewave[lo3] / 1e4, templatewave[hi3] / 1e4)
                xx.set_ylim(min(continuum[lo2:hi2]) * 0.9, max(continuum[lo2:hi2]) * 1.1)
                xx.text(0.05, 0.9, label, ha='left', va='center',
                        transform=xx.transAxes, fontsize=8)

            for rem in range(iwave+1, ncols*nrows):
                ax.flat[rem].axis('off')

            ulpos = ax[0, 0].get_position()
            urpos = ax[0, -1].get_position()
            llpos = ax[-1, 0].get_position()
            lrpos = ax[-1, -1].get_position()
            dxlabel = 0.07
            bottom = 0.11
            top = 0.9
            dytitle = 0.06

            xpos = (lrpos.x1 - llpos.x0) / 2. + llpos.x0
            ypos = llpos.y0 - dxlabel
            fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                     ha='center', va='center')

            xpos = ulpos.x0 - 0.09
            ypos = (ulpos.y1 - llpos.y0) / 2. + llpos.y0
            fig.text(xpos, ypos, r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                     ha='center', va='center', rotation=90)

            xpos = (urpos.x1 - ulpos.x0) / 2. + ulpos.x0
            ypos = ulpos.y1 + dytitle
            fig.text(xpos, ypos, f'Continuum Fluxes: {uniqueid}', ha='center', va='center')

            fig.subplots_adjust(left=0.1, right=0.97, bottom=bottom, top=top, wspace=0.23, hspace=0.3)
            fig.savefig(pngfile)#, bbox_inches='tight')
            plt.close()

        return lums, cfluxes


    @staticmethod
    @jit(nopython=True, nogil=True, fastmath=True, cache=True)
    def attenuate(M, A, zfactors, wave, dustflux):
        """
        Compute attenuated version of a model spectrum,
        including contribution of dust emission.

        """
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


    @staticmethod
    @jit(nopython=True, nogil=True, fastmath=True, cache=True)
    def attenuate_nodust(M, A, zfactors):
        """
        Compute attenuated version of a model spectrum M,
        without dust emission.

        """
        # final result is
        # M * (atten ** ebv) * zfactors
        for i in range(len(M)):
            M[i] *= A[i] * zfactors[i]


    def build_stellar_continuum(self, templateflux, templatecoeff,
                                tauv, vdisp=None, conv_pre=None,
                                dust_emission=True):

        """Build a stellar continuum model.

        Parameters
        ----------
        templateflux : :class:`numpy.ndarray` [ntemplates, npix]
            Rest-frame, native-resolution template spectra corresponding to
            `templatewave`.
        templatecoeff : :class:`numpy.ndarray` [ntemplates]
            Column vector of positive coefficients corresponding to each
            template.
        tauv : :class:`float`
            V-band optical depth, tau(V).
        vdisp : :class:`float`
            Velocity dispersion in km/s. If `None`, do not convolve to the
            specified velocity dispersion (usually because `templateflux` has
            already been smoothed to some nominal value).
        conv_pre: :class:`tuple` or None
            Optional preprocessing data to accelerate template convolution with vdisp
            (may be present only if vdisp is not None).
        dust_emission : :class:`bool`
            Model impact of infrared dust emission spectrum. Energy-balance is used
            to compute the normalization of this spectrum.

        Returns
        -------
        contmodel : :class:`numpy.ndarray` [npix]
            Full-wavelength, native-resolution, observed-frame model spectrum.

        """
        if conv_pre is None or vdisp > Templates.MAX_PRE_VDISP:
            # Compute the weighted sum of the templates.
            contmodel = templatecoeff.dot(templateflux)

            # Optionally convolve to the desired velocity dispersion.
            if vdisp is not None:
                contmodel = self.templates.convolve_vdisp(contmodel, vdisp)
        else:
            # if conv_pre is present, it contains flux values for non-convolved
            # regions of template fluxes, plus FTs of tempaltes for convolved
            # region.  Both must be combined using template coefficients.
            flux_lohi, ft_flux_mid, fft_len = conv_pre

            # Compute the weighted sum of the templates.
            cont_lohi   = templatecoeff.dot(flux_lohi)
            ft_cont_mid = templatecoeff.dot(ft_flux_mid)

            # Convolve to the desired velocity dispersion. Use the vdisp
            # convolution that takes precomputed FT of flux for convolved
            # region.
            flux_len = templateflux.shape[1]
            contmodel = self.templates.convolve_vdisp_from_pre(
                cont_lohi, ft_cont_mid, flux_len, fft_len, vdisp)

            # sanity check for debugging
            #contmodel0 = templateflux.dot(templatecoeff)
            #contmodel0 = self.templates.convolve_vdisp(contmodel0, vdisp)
            #print("DIFF ", np.max(np.abs(contmodel - contmodel0)))

        # [3] - Apply dust attenuation; ToDo: allow age-dependent
        # attenuation. Also compute the bolometric luminosity before and after
        # attenuation but only if we have dustflux.

        # [4] - Optionally add the dust emission spectrum, conserving the total
        # (absorbed + re-emitted) energy. NB: (1) dustflux must be normalized to
        # unity already (see templates.py); and (2) we are ignoring dust
        # self-absorption, which should be mostly negligible.

        # [5] - Redshift factors.

        # Do this part in Numpy because it is very slow in Numba unless
        # accelerated transcendentals are available via, e.g., Intel SVML.
        A = -tauv * self.templates.dust_klambda
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


    def continuum_to_spectroscopy(self, contmodel, interp=False):
        """
        Synthesize spectroscopy from a continuum model.

        Parameters
        ----------
        contmodel : :class:`numpy.ndarray` [npix]
            Full-wavelength, native-resolution, observed-frame model spectrum.
        interp : :class:`bool`
            For cosmetic (plotting) purposes, interpolate over masked pixels.

        Returns
        -------
        modelflux : :class:`numpy.ndarray` [nwave]
            Observed-frame model spectrum at the instrumental resolution and
            wavelength sampling given by `specres` and `specwave`.

        """
        camerapix = self.data['camerapix']
        specwave = self.data['wave']
        specres = self.data['res']
        specmask = self.data['mask']

        modelflux = np.empty(self.wavelen)

        for icam, (s, e) in enumerate(camerapix):
            resampflux = trapz_rebin(
                self.ztemplatewave, contmodel,
                specwave[icam], pre=self.spec_pre[icam])
            specres[icam].dot(resampflux, out=modelflux[s:e])

            # optionally interpolate the model over masked pixels, for cosmetic
            # purposes
            if interp and np.any(specmask[icam]):
                mask = specmask[icam]
                modelflux[s:e][mask] = np.interp(specwave[icam][mask], specwave[icam][~mask], modelflux[s:e][~mask])

        return modelflux


    def continuum_to_photometry(self, contmodel, filters=None,
                                phottable=False, get_abmag=False):
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

        modelmaggies = Photometry.get_ab_maggies_unchecked(
            filters, contmodel, self.ztemplatewave, pre=maggies_pre)

        if not phottable:
            modelphot = Photometry.get_photflam(modelmaggies, effwave)
        else:
            modelmaggies /= self.fluxnorm * self.massnorm
            modelphot = Photometry.parse_photometry(self.phot.bands, modelmaggies, effwave,
                                                    nanomaggies=False, get_abmag=get_abmag)
        return modelphot


    def _stellar_objective(self, params, templateflux, dust_emission,
                           fit_vdisp, conv_pre, objflam, objflamistd,
                           specflux, specistd, synthphot, synthspec):
        """Objective function for fitting a stellar continuum.

        """
        assert (synthphot or synthspec), "request for empty residuals!"

        if fit_vdisp:
            tauv, vdisp = params[:2]
            templatecoeff = params[2:]
        else:
            tauv = params[0]
            vdisp = None
            templatecoeff = params[1:]

        fullmodel = self.build_stellar_continuum(
            templateflux, templatecoeff, tauv=tauv,
            vdisp=vdisp, conv_pre=conv_pre,
            dust_emission=dust_emission)

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
            np.subtract(modelflux, specflux, spec_resid)
            spec_resid *= specistd

        if synthphot:
            modelflam = self.continuum_to_photometry(fullmodel)
            phot_resid = resid[resid_split:]
            np.subtract(modelflam, objflam, phot_resid)
            phot_resid *= objflamistd

        return resid


    def fit_stellar_continuum(self, templateflux, fit_vdisp=False, conv_pre=None,
                              vdisp_guess=None, tauv_guess=None, vdisp_bounds=None,
                              tauv_bounds=None, coeff_guess=None, dust_emission=True,
                              objflam=None, objflamistd=None, specflux=None,
                              specistd=None, synthphot=False, synthspec=False):
        """Fit a stellar continuum using bounded non-linear least-squares.

        Parameters
        ----------
        templateflux : :class:`numpy.ndarray` [ntemplate, npix]
            Grid of input (model) spectra.
        fit_vdisp : :class:`bool`
            If `True`, solve for the velocity dispersion;
            if `False`, use a nominal dispersion.
        conv_pre : :class:`tuple` of None
            If not None, preprocessing data for convolving templateflux
            with vdisp values.  (Occurs only if fit_vdisp is True.)
        vdisp_guess : :class:`float`
            Guess for scalar value of the velocity dispersion if fitting.
        tauv_guess : :class:`float`
            Guess scalar value of the dust optical depth.
        coeff_guess : :class:`numpy.ndarray` [ntemplates]
            Guess of the template coefficients.
        vdisp_bounds : :class:`tuple`
            Two-element list of minimum and maximum allowable values of the
            velocity dispersion; only used if `fit_vdisp=True`.
        tauv_bounds : :class:`tuple`
            Two-element list of minimum and maximum allowable values of the
            V-band optical depth, tau(V).
        dust_emission : :class:`bool`
            Model impact of infrared dust emission spectrum. Energy-balance is used
            to compute the normalization of this spectrum.
        objflam: :class: `numpy.ndarray`
            Measured object photometry (used if fitting photometry).
        objflamistd: :class: `numpy.ndarray`
            Sqrt of inverse variance of objflam (used if fitting photometry).
        specflux : :class:`numpy.ndarray` [nwave]
            Observed-frame spectrum in 10**-17 erg/s/cm2/A corresponding to
            `specwave` (used if fitting spectroscopy).
        specfluxistd : :class:`numpy.ndarray` [nwave]
            Sqrt of inverse variance of the observed-frame spectrum
            (used if fitting spectroscopy).
        synthphot: :class:`bool`
            True iff fitting objective includes object photometry.
        synthspec: :class:`bool`
            True iff fitting objective includes observed spectrum.

        Returns
        -------
        tauv : :class:`float`
            Maximum-likelihood V-band optical depth.
        vdisp : :class:`float`
            Maximum-likelihood velocity dispersion in km/s, or
            nominal velocity if it was not fitted.
        templatecoeff : :class:`numpy.ndarray` [ntemplate]
            Column vector of maximum-likelihood template coefficients.
        resid: :class:`numpy.ndarray`
            Vector of residuals at final parameter values.

        Notes
        -----
        This method supports several different fitting 'modes', depending on the
        whether vdisp is fitted or left nominal, and whether the caller requests
        fitting a model against spectroscopy, photometry, or both.

        In all cases, we solve for dust attenuation via the tau(V) parameter and
        we also include IGM attenuation.

        """
        from scipy.optimize import least_squares

        if tauv_guess is None:
            tauv_guess = self.tauv_guess
        if vdisp_guess is None:
            vdisp_guess = self.vdisp_guess
        if tauv_bounds is None:
            tauv_bounds = self.tauv_bounds
        if vdisp_bounds is None:
            vdisp_bounds = self.vdisp_bounds

        ntemplates = templateflux.shape[0]

        # Unpack the input data to infer the fitting "mode" and the objective
        # function.
        farg = {
            'templateflux':  templateflux,
            'dust_emission': dust_emission,
            'fit_vdisp':     fit_vdisp,
            'conv_pre':      conv_pre,
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

        coeff_bounds = (0., 1e6)

        if fit_vdisp:
            initial_guesses = np.array((tauv_guess, vdisp_guess))
            bounds = [tauv_bounds, vdisp_bounds]
        else:
            initial_guesses = np.array((tauv_guess,))
            bounds = [tauv_bounds]

        initial_guesses = np.concatenate((initial_guesses, coeff_guess))
        bounds = bounds + [coeff_bounds] * ntemplates
        #xscale = np.hstack(([0.1, 50.], np.ones(ntemplates) * 1e-1))

        # NB: `x_scale` has been set to `jac` here to help with the numerical
        # convergence. There may be faster ways, of course...
        fit_info = least_squares(self._stellar_objective, initial_guesses, kwargs=farg,
                                 bounds=tuple(zip(*bounds)), method='trf',
                                 tr_solver='exact', tr_options={'regularize': True},
                                 x_scale='jac', max_nfev=5000, ftol=1e-6, xtol=1e-10)#, verbose=2)
        bestparams = fit_info.x
        resid      = fit_info.fun

        if fit_vdisp:
            tauv, vdisp = bestparams[:2]
            templatecoeff = bestparams[2:]
        else:
            tauv = bestparams[0]
            templatecoeff = bestparams[1:]
            vdisp = self.templates.vdisp_nominal

        return tauv, vdisp, templatecoeff, resid


    def stellar_continuum_chi2(self, resid, ncoeff, vdisp_fitted,
                               split=0, ndof_spec=0, ndof_phot=0):
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
        # tauv is always a free parameter
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


def can_compute_vdisp(redshift, specwave, min_restrange=(3800., 4800.), fit_restrange=(3800., 6000.)):
    """Determine if we can solve for the velocity dispersion.

    min_restrange - minimum required rest wavelength range
    fit_restrange - fitted rest wavelength range

    """
    restwave = specwave / (1. + redshift)
    minwave = np.min(restwave)
    maxwave = np.max(restwave)

    compute_vdisp = (minwave < min_restrange[0]) and (maxwave > min_restrange[1])

    if compute_vdisp:
        s = np.searchsorted(restwave, fit_restrange[0], 'left')
        e = np.searchsorted(restwave, fit_restrange[1], 'left')
        log.debug(f'Solving for vdisp: min(restwave)={minwave:.0f}<{min_restrange[0]:.0f} A, ' + \
                  f'max(restwave)={maxwave:.0f}>{min_restrange[1]:.0f} A')
    else:
        s, e = 0, 0

    return compute_vdisp, (s, e)


def continuum_fastphot(redshift, objflam, objflamivar, CTools, uniqueid=0,
                       nmonte=10, rng=None, debug_plots=False):
    """Model the broadband photometry.

    """
    data = CTools.data
    templates = CTools.templates
    agekeep = CTools.agekeep
    nage = CTools.nage

    vdisp = templates.vdisp_nominal

    ndof_phot = np.sum(objflamivar > 0.)

    if ndof_phot == 0:
        log.info('All photometry is masked.')
        coeff = np.zeros(nage) # nage not nsed
        rchi2_phot = 0.
        sedmodel = np.zeros(len(templates.wave))
        sedmodel_nolines = np.zeros(len(templates.wave))
        dn4000_model = 0.

        coeff_monte = None
        tauv_monte = None
        sedmodel_monte = None
        sedmodel_nolines_monte = None

        tauv_ivar = 0.
        dn4000_model_ivar = 0.

    else:
        objflamistd = np.sqrt(objflamivar)

        t0 = time.time()
        def do_fit(objflam):
            tauv, _, coeff, resid = CTools.fit_stellar_continuum(
                templates.flux_nomvdisp[agekeep, :], fit_vdisp=False,
                objflam=objflam, objflamistd=objflamistd, synthphot=True,
                synthspec=False)

            if np.all(coeff == 0.):
                log.warning('Continuum coefficients are all zero.')
                sedmodel = np.zeros(templates.npix)
                sedmodel_nolines = np.zeros(templates.npix)
                dn4000_model = 0.
            else:
                # Get the best-fitting model with and without line-emission.
                sedmodel = CTools.optimizer_saved_contmodel
                sedmodel_nolines = CTools.build_stellar_continuum(
                    templates.flux_nolines_nomvdisp[agekeep, :], coeff,
                    tauv=tauv, vdisp=None, dust_emission=False)

                # Measure Dn(4000) from the line-free model.
                dn4000_model, _ = Photometry.get_dn4000(
                    templates.wave, sedmodel_nolines, rest=True)

            return (tauv, coeff, sedmodel, sedmodel_nolines, dn4000_model, resid)

        (tauv, coeff, sedmodel, sedmodel_nolines, dn4000_model, resid) = do_fit(objflam)

        if np.all(coeff == 0.):
            rchi2_phot = 0.
        else:
            _, rchi2_phot, _ = CTools.stellar_continuum_chi2(
                resid, ncoeff=len(coeff), vdisp_fitted=False,
                ndof_phot=ndof_phot)

        dt = time.time()-t0
        log.info(f'Fitting {nage} models took {dt:.2f} seconds ' + \
                 f'[rchi2_phot={rchi2_phot:.1f}, ndof={ndof_phot:.0f}].')

        # Monte Carlo to get tauv_ivar and coeff_monte.
        if nmonte > 0:
            objflamstd = np.zeros_like(objflamistd)
            I = objflamistd > 0.
            objflamstd[I] = 1. / objflamistd[I]
            objflam_monte = rng.normal(objflam[np.newaxis, :],
                                       objflamstd[np.newaxis, :],
                                       size=(nmonte, len(objflam)))

            res = [do_fit(*args) for args in zip(objflam_monte)]
            (tauv_monte, coeff_monte, sedmodel_monte, sedmodel_nolines_monte,
             dn4000_model_monte, _) = tuple(zip(*res))

            tauv_ivar = var2ivar(np.var(tauv_monte))
            dn4000_model_ivar = var2ivar(np.var(dn4000_model_monte))

            msg = []
            for label, units, val, val_ivar in zip(
                    ['Model Dn(4000)', 'tau(V)'], ['', ' mag'],
                    [dn4000_model, tauv], [dn4000_model_ivar, tauv_ivar]):
                var_msg = f'+/-{1./np.sqrt(val_ivar):.3f}' if val_ivar > 0. else ''
                msg.append(f'{label}={val:.3f}{var_msg}{units}')
            msg.append(f'vdisp={vdisp:.0f} km/s')
            log.info(' '.join(msg))

    return (coeff, coeff_monte, rchi2_phot, tauv, tauv_monte, tauv_ivar, vdisp,
            dn4000_model, dn4000_model_ivar, sedmodel, sedmodel_monte,
            sedmodel_nolines, sedmodel_nolines_monte)


def vdisp_by_chi2scan(CTools, templates, uniqueid, specflux, specwave,
                      specistd, fitmask, agekeep, deltachi2min=25.,
                      fit_for_min=False, debug_plots=False):
    """Determine the velocity dispersion via a chi2 scan.

    """
    from fastspecfit.util import find_minima, minfit

    ngrid = len(CTools.vdisp_grid)
    chi2grid = np.zeros(ngrid)
    for iv, vdisp1 in enumerate(CTools.vdisp_grid):
        # convolve the templates at the derived vdisp and fit
        input_templateflux_nolines = templates.convolve_vdisp(
            templates.flux_nolines[agekeep, :], vdisp1)
        tauv, _, coeff, resid1 = CTools.fit_stellar_continuum(
            input_templateflux_nolines, fit_vdisp=False, conv_pre=None,
            #tauv_bounds=(0., 2.),
            specflux=specflux, specistd=specistd*fitmask,
            dust_emission=False, synthspec=True)
        chi2grid[iv] = resid1.dot(resid1)

    # Require the peak-to-peak delta-chi2 to be at least deltachi2min and the
    # minimum to not be on either endpoint.
    imin = find_minima(chi2grid)[0]
    deltachi2 = np.ptp(chi2grid)
    if deltachi2 < deltachi2min or imin == 0 or imin == ngrid-1:
        vdisp_init = CTools.vdisp_grid[imin]
        vdisp = templates.vdisp_nominal
        vdisp_ivar = 0.
        if deltachi2 < deltachi2min:
            log.info('Initial velocity dispersion fit failed: delta-chi2=' + \
                     f'{deltachi2:.0f}<{deltachi2min:.0f}')
        else:
            log.info('Initial velocity dispersion fit failed: vdisp_init=' + \
                     f'{vdisp_init:.0f} km/s a boundary value.')
    else:
        vdisp = CTools.vdisp_grid[imin]
        vdisp_ivar = 1. # =! 0.
        log.info('Initial velocity dispersion fit succeeded: delta-chi2=' + \
                 f'{deltachi2:.0f}>{deltachi2min:.0f}, vdisp_init={vdisp:.0f} km/s')

    # Optionally fit for the minimum (best) value (only useful with a dense
    # velocity dispersion grid and so deprecated by default).
    if fit_for_min:
        vdisp, vdisp_sigma, chi2min, warn, (a, b, c) = minfit(
            CTools.vdisp_grid[imin-1:imin+2], chi2grid[imin-1:imin+2],
            return_coeff=True)

        # Did fitting fail?
        if vdisp < 0.:
            vdisp = templates.vdisp_nominal
            vdisp_ivar = 0.
            chi2min = 0.
        else:
            vdisp_ivar = var2ivar(vdisp_sigma, sigma=True)

        if debug_plots:
            import matplotlib.pyplot as plt
            import seaborn as sns
            sns.set(context='talk', style='ticks', font_scale=0.8)

            pngfile = f'qa-vdisp-chi2scan-{uniqueid}.png'

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(CTools.vdisp_grid, chi2grid-chi2min, marker='s', s=50, color='gray', edgecolor='k')
            if not np.all(np.array([a, b, c]) == 0.):
                yquad = np.polyval([a, b, c], CTools.vdisp_grid[imin-3:imin+4])-chi2min
                ax.plot(CTools.vdisp_grid[imin-3:imin+4], yquad, lw=2, ls='--')
                #ax.set_ylim(-0.1*np.max(yquad), np.min((3.*np.max(yquad), np.max(chi2grid-chi2min))))
            ax.set_xlabel(r'$\sigma_{star}$ (km/s)')
            if vdisp_ivar > 0:
                txt = r'$\sigma_{star}$='+f'{vdisp:.0f}'+r'$\pm$'+f'{vdisp_sigma:.0f} km/s'
                ax.set_ylabel(r'$\Delta\chi^2$')
            else:
                txt = r'$\sigma_{star}$='+f'{vdisp:.0f} km/s'
                ax.set_ylabel(r'$\chi^2$')
            ax.text(0.9, 0.9, txt, ha='right', va='center', transform=ax.transAxes)
            ax.set_title(r'Velocity Dispersion $\chi^2$ Scan: '+f'{uniqueid}')
            fig.tight_layout()
            fig.savefig(pngfile)#, bbox_inches='tight')
            plt.close()
            log.info(f'Wrote {pngfile}')

    return vdisp, vdisp_ivar


def _continuum_nominal_vdisp(CTools, templates, specflux, specwave,
                             specistd, agekeep, compute_chi2=False):
    """Support routine to fit a spectrum at the nominal velocity dispersion.

    """
    tauv, vdisp, coeff, resid = CTools.fit_stellar_continuum(
        templates.flux_nolines_nomvdisp[agekeep, :], fit_vdisp=False,
        conv_pre=None, specflux=specflux, specistd=specistd,
        dust_emission=False, synthspec=True)
    contmodel = CTools.optimizer_saved_contmodel.copy() # copy needed??

    if compute_chi2:
        chi2 = resid.dot(resid)
    else:
        chi2 = 1e6

    return tauv, vdisp, coeff, contmodel, chi2


def continuum_fastspec(redshift, objflam, objflamivar, CTools, nmonte=10,
                       rng=None, uniqueid=0, no_smooth_continuum=False,
                       debug_plots=False):
    """Jointly model the spectroscopy and broadband photometry.

    """
    from fastspecfit.util import find_minima, minfit

    data = CTools.data
    phot = CTools.phot
    templates = CTools.templates
    agekeep = CTools.agekeep
    nage = CTools.nage

    # Combine all three cameras; we will unpack them to build the
    # best-fitting model (per-camera) below.
    specwave = np.hstack(data['wave'])
    specflux = np.hstack(data['flux'])
    specivar_nolinemask = np.hstack(data['ivar'])
    slinemask = np.hstack(data['linemask'])
    specivar = specivar_nolinemask * np.logical_not(slinemask) # mask emission lines
    npix = len(specwave)

    if np.all(specivar == 0.) or np.any(specivar < 0.):
        errmsg = 'All pixels are masked or some inverse variances are negative!'
        log.critical(errmsg)
        raise ValueError(errmsg)

    ncam = len(data['cameras'])
    snrmsg = f"Median spectral S/N_{data['cameras'][0]}={data['snr'][0]:.2f}"
    for icam in range(1, ncam):
        snrmsg += f" S/N_{data['cameras'][icam]}={data['snr'][icam]:.2f}"
    log.info(snrmsg)

    specistd = np.sqrt(specivar)
    objflamistd = np.sqrt(objflamivar)
    ndof_cont = np.sum(specivar > 0.)
    ndof_phot = np.sum(objflamivar > 0.)

    if nmonte > 0 and (ndof_cont > 0 or ndof_phot > 0):
        # ndof_cont==0 should never happen...

        # Use specivar_nolinemask here rather than specivar because we are
        # going to use the same set of realizations in line-fitting.
        specstd = np.zeros_like(specivar_nolinemask)
        I = specivar_nolinemask > 0.
        specstd[I] = 1. / np.sqrt(specivar_nolinemask[I])
        specflux_monte = rng.normal(specflux[np.newaxis, :],
                                    specstd[np.newaxis, :],
                                    size=(nmonte, len(specflux)))

        if ndof_phot > 0:
            objflamstd = np.zeros_like(objflamistd)
            I = objflamistd > 0.
            objflamstd[I] = 1. / objflamistd[I]
            objflam_monte = rng.normal(objflam[np.newaxis, :],
                                       objflamstd[np.newaxis, :],
                                       size=(nmonte, len(objflam)))
        else:
            # should be all zeros
            objflam_monte = np.repeat(objflam, nmonte).reshape(nmonte, len(objflam))
    else:
        specflux_monte = None

    # Attempt to solve for the velocity dispersion based on the rest-wavelength coverage.
    compute_vdisp, (vdisp_s, vdisp_e) = can_compute_vdisp(redshift, specwave)

    if not compute_vdisp:
        # Fit to the cached templates at the nominal velocity dispersion.
        tauv, vdisp, coeff, contmodel, _ = _continuum_nominal_vdisp(
            CTools, templates, specflux, specwave,
            specistd, agekeep, compute_chi2=False)

        vdisp_ivar = 0.

        input_templateflux = templates.flux_nomvdisp[agekeep, :]
        input_templateflux_nolines = templates.flux_nolines_nomvdisp[agekeep, :]

        log.debug('Insufficient wavelength coverage to compute velocity ' + \
                  f'dispersion; adopting {vdisp:.0f} km/s')
    else:
        t0 = time.time()

        # Fit for the velocity dispersion over a restricted wavelength range.
        fitmask = np.zeros(len(specflux), bool)
        fitmask[vdisp_s:vdisp_e] = True

        # First, perform a basic chi2 scan over a limited set of vdisp values.
        vdisp, vdisp_ivar = vdisp_by_chi2scan(
            CTools, templates, uniqueid, specflux, specwave,
            specistd, fitmask, agekeep, deltachi2min=25.,
            fit_for_min=False, debug_plots=debug_plots)

        # If the scan is unsuccessful, adopt the nominal velocity dispersion
        # and continue....
        if vdisp_ivar == 0.:
            tauv, vdisp, coeff, contmodel, _ = _continuum_nominal_vdisp(
                CTools, templates, specflux, specwave,
                specistd, agekeep, compute_chi2=False)

            input_templateflux = templates.flux_nomvdisp[agekeep, :]
            input_templateflux_nolines = templates.flux_nolines_nomvdisp[agekeep, :]
        else:
            # ...otherwise fit for the maximum likelihood value.

            def do_fit_vdisp(specflux):
                tauv, vdisp, coeff, resid = CTools.fit_stellar_continuum(
                    templates.flux_nolines[agekeep, :], fit_vdisp=True,
                    conv_pre=input_conv_pre_nolines, specflux=specflux,
                    specistd=specistd*fitmask, dust_emission=False, synthspec=True)
                age = coeff.dot(templates.info['age']) / np.sum(coeff) / 1e9 # luminosity-weighted [Gyr]
                return (tauv, vdisp, coeff, age, resid)

            input_conv_pre_nolines = templates.conv_pre_select(
                templates.conv_pre_nolines, agekeep)

            (tauv, vdisp, coeff, age, resid) = do_fit_vdisp(specflux)

            contmodel = CTools.optimizer_saved_contmodel

            # Get the templates, coefficients, and model at the derived vdisp.
            input_templateflux = templates.convolve_vdisp(
                templates.flux[agekeep, :], vdisp)
            input_templateflux_nolines = templates.convolve_vdisp(
                templates.flux_nolines[agekeep, :], vdisp)

            # Monte Carlo to get vdisp_ivar (and the diagnostic plot, if
            # requested).
            if specflux_monte is not None:
                res = [do_fit_vdisp(sf) for sf in specflux_monte]
                (tauv_monte, vdisp_monte, _, age_monte, _) = tuple(zip(*res))

                vdisp_ivar = var2ivar(np.var(vdisp_monte))

                if debug_plots and vdisp_ivar > 0.:
                    import matplotlib.pyplot as plt
                    import corner as cn
                    import seaborn as sns

                    pngfile = f'qa-vdisp-{uniqueid}.png'

                    sns.set(context='talk', style='ticks', font_scale=0.6)
                    colors = sns.color_palette()

                    dkw = {'color': colors[1], 'ms': 10, 'alpha': 0.75, 'mec': 'k'}
                    hkw = {'fill': True, 'alpha': 0.75, 'color': 'gray',
                           'align': 'left', 'edgecolor': 'k'}
                    ndim = 3

                    tauv_sigma = np.std(tauv_monte)
                    age_sigma = np.std(age_monte)
                    vdisp_sigma = np.std(vdisp_monte)

                    plotdata = np.vstack((vdisp_monte, tauv_monte, age_monte)).T
                    truths = [vdisp, tauv, age]
                    labels = [r'$\sigma_{star}$ (km/s)', r'$\tau_{V}$', 'Age (Gyr)']
                    titles = [r'$\sigma_{star}$='+f'{vdisp:.0f}'+r'$\pm$'+f'{vdisp_sigma:.0f} km/s',
                              r'$\tau_{V}$='+f'{tauv:.2f}'+r'$\pm$'+f'{tauv_sigma:.2f}',
                              f'Age={age:.2f}'+r'$\pm$'+f'{age_sigma:.2f} Gyr']
                    sig = [max(5.*vdisp_sigma, 3.), max(5.*tauv_sigma, 0.005), max(5.*age_sigma, 0.005)]
                    ranges = ((vdisp-sig[0], vdisp+sig[0]), (tauv-sig[1], tauv+sig[1]), (age-sig[2], age+sig[2]))

                    bins = nmonte // 3
                    if bins < 10:
                        bins = 10
                    fig = cn.corner(plotdata, bins=bins, smooth=None, plot_density=False,
                                    plot_contours=False, range=ranges,
                                    data_kwargs=dkw, hist_kwargs=hkw, labels=labels)
                    ax = np.array(fig.axes).reshape((ndim, ndim))
                    for ii, mlval, sig in zip(range(ndim), (vdisp, tauv, age), (vdisp_sigma, tauv_sigma, age_sigma)):
                        ax[ii, ii].axvline(mlval, color=colors[0], lw=2, ls='-')
                        ax[ii, ii].axvline(mlval+sig, color=colors[0], lw=1, ls='--')
                        ax[ii, ii].axvline(mlval-sig, color=colors[0], lw=1, ls='--')
                        ax[ii, ii].set_title(titles[ii])
                        if ii == 0:
                            ax[ii, ii].set_ylabel('Number of\nRealizations')
                        else:
                            xx = ax[ii, ii].twinx()
                            xx.set_yticklabels([])
                            xx.set_ylabel('Number of\nRealizations')
                            xx.tick_params(right=False)
                    for yi in range(ndim):
                        for xi in range(yi):
                            ax[yi, xi].axvline(truths[xi], color=colors[0], lw=1, ls='-', alpha=0.75)
                            ax[yi, xi].axhline(truths[yi], color=colors[0], lw=1, ls='-', alpha=0.75)
                    fig.suptitle(f'Velocity Dispersion: {uniqueid}')

                    fig.subplots_adjust(left=0.13, right=0.9, bottom=0.13, top=0.91, wspace=0.14, hspace=0.14)
                    fig.savefig(pngfile)#, bbox_inches='tight')
                    plt.close()
                    log.info(f'Wrote {pngfile}')

        log.debug(f'Estimating the velocity dispersion took {time.time()-t0:.2f} seconds.')

    # Next, estimate the aperture correction.
    apercorrs = np.ones(len(phot.synth_bands))
    median_apercorr = 1.

    if not np.any(phot.bands_to_fit):
        log.info('Skipping aperture correction since no bands were fit.')
    else:
        if np.all(coeff == 0.):
            log.warning('Unable to estimate aperture correction because ' + \
                        'continuum coefficients are all zero; adopting 1.0.')
        else:
            sedflam = CTools.continuum_to_photometry(
                contmodel, filters=phot.synth_filters[data['photsys']])

            I = np.isin(data['photometry']['band'], phot.synth_bands)
            objflam_aper = CTools.fluxnorm * data['photometry'][I]['flam'].value

            I = ((objflam_aper > 0.) & (sedflam > 0.))
            if np.any(I):
                apercorrs[I] = objflam_aper[I] / sedflam[I]

            I = (apercorrs > 0.)
            if np.any(I):
                median_apercorr = median(apercorrs[I])

                if median_apercorr <= 0.:
                    log.warning('Aperture correction not well-defined; adopting 1.0.')
                    median_apercorr = 1.
                else:
                    log.info(f'Median aperture correction {median_apercorr:.3f} ' + \
                             f'[{np.min(apercorrs):.3f}-{np.max(apercorrs):.3f}].')

    # Now do the full spectrophotometric fit with the velocity dispersion
    # fixed and the bounds on tauv relaxed.
    t0 = time.time()

    def do_fit_full(objflam, specflux):
        tauv, _, coeff, resid = CTools.fit_stellar_continuum(
            input_templateflux, fit_vdisp=False, conv_pre=None,
            objflam=objflam, objflamistd=objflamistd,
            specflux=specflux*median_apercorr,
            specistd=specistd/median_apercorr,
            synthphot=True, synthspec=True)

        if np.all(coeff == 0.):
            log.warning('Continuum coefficients are all zero.')
            sedmodel = np.zeros(templates.npix)
            sedmodel_nolines = np.zeros(templates.npix)
            desimodel_nolines = np.zeros(len(specflux))
            dn4000_model = 0.
        else:
            # Get the best-fitting model with and without line-emission. Set
            # dust_emission=False for sedmodel_nolines since we only use it to get
            # Dn(4000) and the UV/optical continuum fluxes.
            sedmodel = CTools.optimizer_saved_contmodel
            sedmodel_nolines = CTools.build_stellar_continuum(
                input_templateflux_nolines, coeff, tauv=tauv,
                vdisp=None, conv_pre=None, dust_emission=False)
            desimodel_nolines = CTools.continuum_to_spectroscopy(sedmodel_nolines)

            # Get DN(4000).
            dn4000_model, _ = Photometry.get_dn4000(
                templates.wave, sedmodel_nolines, rest=True)

        return (tauv, coeff, sedmodel, sedmodel_nolines,
                desimodel_nolines, dn4000_model, resid)

    (tauv, coeff, sedmodel, sedmodel_nolines, desimodel_nolines, dn4000_model, resid) = \
        do_fit_full(objflam, specflux)

    if np.all(coeff == 0.):
        rchi2_cont = 0.
        rchi2_phot = 0.
    else:
        _, rchi2_phot, rchi2_cont = CTools.stellar_continuum_chi2(
        resid, ncoeff=nage, vdisp_fitted=False, split=len(specflux),
        ndof_spec=ndof_cont, ndof_phot=ndof_phot)

    log.debug(f'Fitting {nage} models took {time.time()-t0:.2f} seconds ' + \
              f'[rchi2_cont={rchi2_cont:.1f}, ndof={ndof_cont:.0f}; ' + \
              f'rchi2_phot={rchi2_phot:.1f}, ndof={ndof_phot:.0f}].')

    if specflux_monte is not None:
        res = [do_fit_full(*args) for args in zip(objflam_monte, specflux_monte)]

        (tauv_monte, coeff_monte, sedmodel_monte, sedmodel_nolines_monte,
         desimodel_nolines_monte, dn4000_model_monte, _) = tuple(zip(*res))

        continuummodel_monte = np.vstack(desimodel_nolines_monte)

        tauv_ivar = var2ivar(np.var(tauv_monte))
        dn4000_model_ivar = var2ivar(np.var(dn4000_model_monte))
    else:
        coeff_monte = None
        tauv_monte = None
        sedmodel_monte = None
        sedmodel_nolines_monte = None
        desimodel_nolines_monte = None
        continuummodel_monte = None

        tauv_ivar = 0.
        dn4000_model_ivar = 0.

    dn4000, dn4000_ivar = Photometry.get_dn4000(
        specwave, specflux, flam_ivar=specivar_nolinemask,
        redshift=redshift, rest=False)

    # Get the smooth continuum.
    t0 = time.time()

    smoothstats = np.zeros(len(data['camerapix']))

    if np.all(coeff == 0.) or no_smooth_continuum:
        smoothcontinuum = np.zeros_like(specwave)
    else:
        # Need to be careful we don't pass a large negative residual
        # where there are gaps in the data.
        residuals = specflux * median_apercorr - desimodel_nolines
        I = ((specflux == 0.) & (specivar == 0.))
        residuals[I] = 0.

        smoothcontinuum = CTools.smooth_continuum(
            specwave, residuals, specivar / median_apercorr**2,
            slinemask, uniqueid=data['uniqueid'],
            camerapix=data['camerapix'], debug_plots=debug_plots)

        for icam, (ss, ee) in enumerate(data['camerapix']):
            I = ((specflux[ss:ee] != 0.) & (specivar[ss:ee] != 0.) & (smoothcontinuum[ss:ee] != 0.))
            if np.count_nonzero(I) > 3: # require three good pixels to compute the mean
                smoothstats[icam] = median(smoothcontinuum[ss:ee][I] / specflux[ss:ee][I])

    log.debug(f'Deriving the smooth continuum took {time.time()-t0:.2f} seconds.')

    return (coeff, coeff_monte, rchi2_cont, rchi2_phot, median_apercorr, apercorrs,
            tauv, tauv_monte, tauv_ivar, vdisp, vdisp_ivar, dn4000, dn4000_ivar,
            dn4000_model, dn4000_model_ivar, sedmodel, sedmodel_nolines,
            desimodel_nolines, smoothcontinuum, smoothstats, specflux_monte,
            sedmodel_monte, sedmodel_nolines_monte, continuummodel_monte)


def continuum_specfit(data, fastfit, specphot, templates, igm, phot,
                      nmonte=10, seed=1, constrain_age=False,
                      no_smooth_continuum=False, fitstack=False,
                      fastphot=False, debug_plots=False):
    """Fit the non-negative stellar continuum of a single spectrum.

    Parameters
    ----------
    data : :class:`dict`
        Dictionary of input spectroscopy (plus ancillary data) populated by
        :func:`fastspecfit.io.DESISpectra.read`.

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
    tall = time.time()

    redshift = data['redshift']
    if redshift <= 0.:
        log.warning('Input redshift not defined, zero, or negative!')

    if fitstack:
        FLUXNORM = 1.
    else:
        from fastspecfit.util import FLUXNORM

    objflam = data['photometry']['flam'].value * FLUXNORM
    objflamivar = (data['photometry']['flam_ivar'].value / FLUXNORM**2) * phot.bands_to_fit

    if np.any(phot.bands_to_fit):
        # Require at least one *optical* photometric band; do not just fit the
        # IR because we will not be able to compute the aperture correction.
        lambda_eff = data['photometry']['lambda_eff'].value
        opt = ((lambda_eff > 3e3) & (lambda_eff < 1e4))
        if np.all(objflamivar[opt] == 0.):
            log.warning('All optical bands are masked; masking all photometry.')
            objflamivar[:] = 0.

    # Instantiate the continuum tools class.
    CTools = ContinuumTools(data, templates, phot, igm, fastphot=fastphot,
                            vdisp_guess=templates.vdisp_nominal,
                            fluxnorm=FLUXNORM, constrain_age=constrain_age)

    # Instantiate the random-number generator.
    if nmonte > 0:
        rng = np.random.default_rng(seed=seed)
    else:
        rng = None

    if fastphot:
        # Photometry-only fitting.
        (coeff, coeff_monte, rchi2_phot, tauv, tauv_monte, tauv_ivar, vdisp, dn4000_model,
         dn4000_model_ivar, sedmodel, sedmodel_monte, sedmodel_nolines, sedmodel_nolines_monte) = \
            continuum_fastphot(redshift, objflam, objflamivar, CTools,
                               uniqueid=data['uniqueid'], debug_plots=debug_plots,
                               nmonte=nmonte, rng=rng)
    else:
        (coeff, coeff_monte, rchi2_cont, rchi2_phot, median_apercorr, apercorrs,
         tauv, tauv_monte, tauv_ivar, vdisp, vdisp_ivar, dn4000, dn4000_ivar,
         dn4000_model, dn4000_model_ivar, sedmodel, sedmodel_nolines, continuummodel,
         smoothcontinuum, smoothstats, specflux_monte, sedmodel_monte,
         sedmodel_nolines_monte, continuummodel_monte) = \
             continuum_fastspec(redshift, objflam, objflamivar, CTools,
                                nmonte=nmonte, rng=rng, uniqueid=data['uniqueid'],
                                debug_plots=debug_plots, no_smooth_continuum=no_smooth_continuum)

        data['apercorr'] = median_apercorr # needed for the line-fitting

        # populate the output table
        for icam, cam in enumerate(np.atleast_1d(data['cameras'])):
            fastfit[f'SNR_{cam.upper()}'] = data['snr'][icam]

        msg = ['Smooth continuum correction:']
        for cam, corr in zip(np.atleast_1d(data['cameras']), smoothstats):
            fastfit[f'SMOOTHCORR_{cam.upper()}'] = corr * 100. # [%]
            msg.append(f'{cam}={100.*corr:.3f}%')
        log.info(' '.join(msg))

    #result['Z'] = redshift
    specphot['SEED'] = seed
    specphot['COEFF'][CTools.agekeep] = coeff
    specphot['RCHI2_PHOT'] = rchi2_phot
    specphot['VDISP'] = vdisp # * u.kilometer/u.second
    specphot['DN4000_MODEL'] = dn4000_model
    specphot['DN4000_MODEL_IVAR'] = dn4000_model_ivar

    if not fastphot:
        specphot['RCHI2_CONT'] = rchi2_cont

        # add the initial line-masking parameters to the output table
        fastfit['INIT_SIGMA_UV'] = data['linesigma_broad']
        fastfit['INIT_SIGMA_NARROW'] = data['linesigma_narrow']
        fastfit['INIT_SIGMA_BALMER'] = data['linesigma_balmer_broad']
        fastfit['INIT_VSHIFT_UV'] = data['linevshift_broad']
        fastfit['INIT_VSHIFT_NARROW'] = data['linevshift_narrow']
        fastfit['INIT_VSHIFT_BALMER'] = data['linevshift_balmer_broad']
        fastfit['INIT_BALMER_BROAD'] = data['balmerbroad']

        fastfit['APERCORR'] = median_apercorr
        for iband, band in enumerate(phot.synth_bands):
            fastfit[f'APERCORR_{band.upper()}'] = apercorrs[iband]
        specphot['DN4000_OBS'] = dn4000
        specphot['DN4000_IVAR'] = dn4000_ivar
        specphot['VDISP_IVAR'] = vdisp_ivar # * (u.second/u.kilometer)**2

    # Compute K-corrections, rest-frame quantities, and physical properties.
    if not np.all(coeff == 0.):

        def do_kcorr(sedmodel, sedmodel_nolines, debug_plots=False):
            synth_absmag, synth_maggies_rest = phot.synth_absmag(
                redshift, data['dmodulus'], CTools.ztemplatewave,
                sedmodel)

            kcorr, absmag, ivarabsmag, synth_bestmaggies = phot.kcorr_and_absmag(
                data['photometry']['nanomaggies'].value,
                data['photometry']['nanomaggies_ivar'].value,
                redshift, data['dmodulus'], data['photsys'],
                CTools.ztemplatewave, sedmodel, synth_absmag,
                synth_maggies_rest)

            lums, cfluxes = CTools.continuum_fluxes(
                sedmodel_nolines, uniqueid=data['uniqueid'],
                debug_plots=debug_plots)

            return (synth_absmag, synth_maggies_rest, kcorr, absmag,
                    ivarabsmag, synth_bestmaggies, lums, cfluxes)

        (synth_absmag, synth_maggies_rest, kcorr, absmag, ivarabsmag, synth_bestmaggies, lums, cfluxes) = \
            do_kcorr(sedmodel, sedmodel_nolines, debug_plots=debug_plots)

        for iband, (band, shift) in enumerate(zip(phot.absmag_bands, phot.band_shift)):
            band = band.upper()
            shift = int(10*shift)
            specphot[f'KCORR{shift:02d}_{band}'] = kcorr[iband] # * u.mag
            specphot[f'ABSMAG{shift:02d}_{band}'] = absmag[iband] # * u.mag
            specphot[f'ABSMAG{shift:02d}_SYNTH_{band}'] = synth_absmag[iband] # * u.mag
            specphot[f'ABSMAG{shift:02d}_IVAR_{band}'] = ivarabsmag[iband] # / (u.mag**2)

        for iband, band in enumerate(phot.bands):
            specphot[f'FLUX_SYNTH_PHOTMODEL_{band.upper()}'] = 1e9 * synth_bestmaggies[iband] # * u.nanomaggy

        lumskeys, _ = CTools.lums_keys()
        for ikey, key in enumerate(lumskeys):
            specphot[key] = lums[ikey]

        cfluxeskeys, _ = CTools.cfluxes_keys()
        for ikey, key in enumerate(cfluxeskeys):
            specphot[key] = cfluxes[ikey]

        # Get the variance via Monte Carlo.
        if sedmodel_monte is not None:
            res = [do_kcorr(sm, snm, False) for sm, snm in zip(sedmodel_monte, sedmodel_nolines_monte)]
            (synth_absmag_monte, _, _, _, _, _, lums_monte, cfluxes_monte) = tuple(zip(*res))

            synth_absmag_var = np.var(synth_absmag_monte, axis=0)
            for band, shift, var in zip(phot.absmag_bands, phot.band_shift, synth_absmag_var):
                if var > TINY:
                    band = band.upper()
                    shift = int(10*shift)
                    specphot[f'ABSMAG{shift:02d}_SYNTH_IVAR_{band}'] = 1. / var

            lums_var = np.var(lums_monte, axis=0)
            for lumkey, var in zip(lumskeys, lums_var):
                if var > TINY:
                    specphot[f'{lumkey}_IVAR'] = 1. / var

            cfluxes_var = np.var(cfluxes_monte, axis=0)
            for cfluxkey, var in zip(cfluxeskeys, cfluxes_var):
                if var > TINY:
                    specphot[f'{cfluxkey}_IVAR'] = 1. / var

        # get the SPS properties
        def _get_sps_properties(coeff):
            tinfo = templates.info[CTools.agekeep]
            mstars = tinfo['mstar'] # [current mass in stars, Msun]
            masstot = coeff.dot(mstars)
            coefftot = np.sum(coeff)
            logmstar = np.log10(CTools.massnorm * masstot)
            zzsun = np.log10(coeff.dot(mstars * 10.**tinfo['zzsun']) / masstot) # mass-weighted
            age = coeff.dot(tinfo['age']) / coefftot / 1e9           # luminosity-weighted [Gyr]
            #age = coeff.dot(mstars * tinfo['age']) / masstot / 1e9  # mass-weighted [Gyr]
            sfr = CTools.massnorm * coeff.dot(tinfo['sfr'])          # [Msun/yr]
            return age, zzsun, logmstar, sfr

        age, zzsun, logmstar, sfr = _get_sps_properties(coeff)
        specphot['TAUV'] = tauv
        specphot['TAUV_IVAR'] = tauv_ivar
        specphot['AGE'] = age
        specphot['ZZSUN'] = zzsun
        specphot['LOGMSTAR'] = logmstar
        specphot['SFR'] = sfr

        if coeff_monte is not None:
            res = [_get_sps_properties(c) for c in coeff_monte]
            age_monte, zzsun_monte, logmstar_monte, sfr_monte = tuple(zip(*res))

            for val_monte, col in zip([age_monte, zzsun_monte, logmstar_monte, sfr_monte],
                                      ['AGE_IVAR', 'ZZSUN_IVAR', 'LOGMSTAR_IVAR', 'SFR_IVAR']):
                val_ivar = var2ivar(np.var(val_monte))
                if val_ivar < F32MAX:
                    specphot[col] = val_ivar

            # optional debugging plot
            if debug_plots:
                import matplotlib.pyplot as plt
                import corner as cn
                import seaborn as sns

                pngfile = f'qa-sps-properties-{data["uniqueid"]}.png'

                sns.set(context='talk', style='ticks', font_scale=0.6)
                colors = sns.color_palette()

                dkw = {'color': colors[1], 'ms': 10, 'alpha': 0.75, 'mec': 'k'}
                hkw = {'fill': True, 'alpha': 0.75, 'color': 'gray',
                       'align': 'left', 'edgecolor': 'k'}
                ndim = 5

                zzsun_sigma = np.std(zzsun_monte)
                tauv_sigma = np.std(tauv_monte)
                sfr_sigma = np.std(sfr_monte)
                logmstar_sigma = np.std(logmstar_monte)
                age_sigma = np.std(age_monte)

                plotdata = np.vstack((zzsun_monte, tauv_monte, sfr_monte, logmstar_monte, age_monte)).T
                truths = [zzsun, tauv, sfr, logmstar, age]
                labels = [r'$Z/Z_{\odot}$', r'$\tau_{V}$', r'SFR ($M_{\odot}/\mathrm{yr}$)',
                          '\n'+r'$\log_{10}(M/M_{\odot})$', 'Age (Gyr)']
                titles = [r'$Z/Z_{\odot}$='+f'{zzsun:.1f}'+r'$\pm$'+f'{zzsun_sigma:.1f}',
                          r'$\tau_{V}$='+f'{tauv:.2f}'+r'$\pm$'+f'{tauv_sigma:.2f}',
                          r'SFR='+f'{sfr:.1f}'+r'$\pm$'+f'{sfr_sigma:.1f}'+r' $M_{\odot}/\mathrm{yr}$',
                          r'$\log_{10}(M/M_{\odot})$='+f'{logmstar:.2f}'+r'$\pm$'+f'{logmstar_sigma:.2f}',
                          f'Age={age:.2f}'+r'$\pm$'+f'{age_sigma:.2f} Gyr']
                sig = [max(5.*zzsun_sigma, 0.1), max(5.*tauv_sigma, 0.005), max(5.*sfr_sigma, 3),
                       max(5.*logmstar_sigma, 0.1), max(5.*age_sigma, 0.005)]
                ranges = [(prop-sig1, prop+sig1) for prop, sig1 in zip([zzsun, tauv, sfr, logmstar, age], sig)]

                bins = nmonte // 3
                if bins < 10:
                    bins = 10
                fig = cn.corner(plotdata, bins=bins, smooth=None, plot_density=False,
                                plot_contours=False, range=ranges,
                                data_kwargs=dkw, hist_kwargs=hkw, labels=labels)
                ax = np.array(fig.axes).reshape((ndim, ndim))
                for ii, mlval, sig in zip(range(ndim), (zzsun, tauv, sfr, logmstar, age),
                                          (zzsun_sigma, tauv_sigma, sfr_sigma, logmstar_sigma, age_sigma)):
                    ax[ii, ii].axvline(mlval, color=colors[0], lw=2, ls='-')
                    ax[ii, ii].axvline(mlval+sig, color=colors[0], lw=1, ls='--')
                    ax[ii, ii].axvline(mlval-sig, color=colors[0], lw=1, ls='--')
                    ax[ii, ii].set_title(titles[ii])
                    if ii == 0:
                        ax[ii, ii].set_ylabel('Number of\nRealizations')
                    else:
                        xx = ax[ii, ii].twinx()
                        xx.set_yticklabels([])
                        xx.set_ylabel('Number of\nRealizations')
                        xx.tick_params(right=False)
                for yi in range(ndim):
                    for xi in range(yi):
                        ax[yi, xi].axvline(truths[xi], color=colors[0], lw=1, ls='-', alpha=0.75)
                        ax[yi, xi].axhline(truths[yi], color=colors[0], lw=1, ls='-', alpha=0.75)
                fig.suptitle(f'SPS Properties: {data["uniqueid"]}')

                fig.subplots_adjust(left=0.1, right=0.92, bottom=0.1, top=0.95, wspace=0.14, hspace=0.14)
                fig.savefig(pngfile)#, bbox_inches='tight')
                plt.close()
                log.info(f'Wrote {pngfile}')


        msg = []
        for label, units, val, col in zip(['vdisp', 'log(M/Msun)', 'tau(V)', 'Age', 'SFR', 'Z/Zsun'],
                                          [' km/s', '', '', ' Gyr', ' Msun/yr', ''],
                                          [vdisp, logmstar, tauv, age, sfr, zzsun],
                                          ['VDISP', 'LOGMSTAR', 'TAUV', 'AGE', 'SFR', 'ZZSUN']):
            ivarcol = f'{col}_IVAR'
            if ivarcol in specphot.value.dtype.names:
                val_ivar = specphot[ivarcol]
                var_msg = f'+/-{1./np.sqrt(val_ivar):.2f}' if val_ivar > 0. else ''
            else:
                var_msg = ''
            msg.append(f'{label}={val:.2f}{var_msg}{units}')
        log.info(' '.join(msg))

    log.debug(f'Continuum-fitting took {time.time()-tall:.2f} seconds.')

    if fastphot:
        return sedmodel, None, None, None
    else:
        # divide out the aperture correction
        continuummodel /= median_apercorr
        smoothcontinuum /= median_apercorr
        if continuummodel_monte is not None:
            continuummodel_monte /= median_apercorr

        return continuummodel, smoothcontinuum, continuummodel_monte, specflux_monte
