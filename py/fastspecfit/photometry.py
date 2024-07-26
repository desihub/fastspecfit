import numpy as np
from astropy.table import Table

from fastspecfit.logger import log

from fastspecfit.cosmo import TabulatedDESI
from fastspecfit.util import C_LIGHT, FLUXNORM

class Photometry(object):
    """Class to load filters and containing filter- and dust-related methods.

    """
    def __init__(self, fphotofile=None, stackfit=False, ignore_photometry=False):
        """
        Parameters
        ----------
        ignore_photometry : :class:`bool`
            Boolean flag indicating whether or not to ignore the broadband
            photometry.
        log : :class:`desiutil.log.logger`
            Logger object.

        """
        from speclite import filters
        import yaml

        if fphotofile is None:
            from importlib import resources
            if stackfit:
                fphotofile = resources.files('fastspecfit').joinpath('data/stacked-phot.yaml')
            else:
                fphotofile = resources.files('fastspecfit').joinpath('data/legacysurvey-dr9.yaml')
        self.fphotofile = fphotofile
        
        try:
            with open(fphotofile, 'r') as F:
                fphoto = yaml.safe_load(F)
        except:
            errmsg = f'Unable to read parameter file {fphotofile}'
            log.critical(errmsg)
            raise ValueError(errmsg)

        self.uniqueid = fphoto['uniqueid']
        self.photounits = fphoto['photounits']

        if 'readcols' in fphoto:
            self.readcols = np.array(fphoto['readcols'])
        if 'dropcols' in fphoto:
            self.dropcols = np.array(fphoto['dropcols'])
        if 'outcols' in fphoto:
            self.outcols = np.array(fphoto['outcols'])
            
        self.bands = np.array(fphoto['bands'])
        self.bands_to_fit = np.array(fphoto['bands_to_fit'])
        self.fluxcols = np.array(fphoto['fluxcols'])
        self.fluxivarcols = np.array(fphoto['fluxivarcols'])
        self.min_uncertainty = np.array(fphoto['min_uncertainty'])
            
        if 'legacysurveydr' in fphoto:
            self.legacysurveydr = fphoto['legacysurveydr']
        if 'viewer_layer' in fphoto:
            self.viewer_layer = fphoto['viewer_layer']
        if 'viewer_pixscale' in fphoto:
            self.viewer_pixscale = fphoto['viewer_pixscale']
        if 'synth_bands' in fphoto:
            self.synth_bands = np.array(fphoto['synth_bands'])
        if 'fiber_bands' in fphoto:
            self.fiber_bands = np.array(fphoto['fiber_bands'])
        
        self.absmag_bands = np.array(fphoto['absmag_bands'])
        self.band_shift = np.array(fphoto['band_shift'])

        # If fphoto['filters'] is a dictionary, then assume that there
        # are N/S filters (as indicated by photsys).
        self.filters = {}
        for key in fphoto['filters']:
            self.filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                        for filtname in fphoto['filters'][key]])
        self.synth_filters = {}
        for key in fphoto['synth_filters']:
            self.synth_filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                              for filtname in fphoto['synth_filters'][key]])
        if 'fiber_bands' in fphoto:
            self.fiber_filters = {}
            for key in fphoto['fiber_filters']:
                self.fiber_filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                                  for filtname in fphoto['fiber_filters'][key]])
        # Simple list of filters.
        self.absmag_filters = filters.FilterSequence([filters.load_filter(filtname) for filtname in fphoto['absmag_filters']])

        if len(self.absmag_bands) != len(self.band_shift):
            errmsg = 'absmag_bands and band_shift must have the same number of elements.'
            log.critical(errmsg)
            raise ValueError(errmsg)
        
        if self.photounits != 'nanomaggies':
            errmsg = 'nanomaggies is the only currently supported photometric unit!'
            log.critical(errmsg)
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
            log.warning(errmsg)
            kcorr        = np.zeros(nabs, dtype='f8')
            absmag       = np.zeros(nabs, dtype='f8')
            ivarabsmag   = np.zeros(nabs, dtype='f8')
            synth_absmag = np.zeros(nabs, dtype='f8')
            synth_maggies_in = np.zeros(len(maggies))
            return kcorr, absmag, ivarabsmag, synth_absmag, synth_maggies_in

        if cosmo is None:
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
                                               zmodelwave)
        filters_out = \
            filters.FilterSequence( [ f.create_shifted(band_shift=bs) for f, bs in zip(absmag_filters, band_shift) ])
        lambda_out = filters_out.effective_wavelengths.value
        
        # Multiply by (1+z) to convert the best-fitting model to the "rest
        # frame".
        synth_outmaggies_rest = self.get_ab_maggies(filters_out,
                                                    zmodelflux * (1. + redshift) / FLUXNORM,
                                                    modelwave)
        
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
    def get_ab_maggies_fast(filters, flux, wave):
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
    def get_ab_maggies(filters, flux, wave):
        
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
                         get_abmag=False):
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

        Require full wavelength coverage over the definition of the index.

        See eq. 11 in Bruzual 1983
        (https://articles.adsabs.harvard.edu/pdf/1983ApJ...273..105B) but with
        the "narrow" definition of Balogh et al. 1999.

        """
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
