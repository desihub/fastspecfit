"""
fastspecfit.photometry
======================

Tools for handling filters and photometry calculations.

"""
import os
import numpy as np
import fitsio
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.util import trapz, C_LIGHT, FLUXNORM


class Photometry(object):
    """Class to load filters and containing filter- and dust-related methods.

    """
    def __init__(self, fphotofile=None, fitstack=False, ignore_photometry=False):

        """
        Parameters
        ----------
        ignore_photometry : :class:`bool`
            Boolean flag indicating whether or not to ignore the broadband
            photometry.

        """
        from speclite import filters
        import yaml

        if fphotofile is None:
            from importlib import resources
            if fitstack:
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
            raise IOError(errmsg)

        self.uniqueid_col = fphoto['uniqueid']
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

        # Deprecated - trim excessive filter wavelengths where the response is
        # effectively ~zero.
        #trim = {
        #    'decam2014-u': [None, 4100.],
        #    'decam2014-g': [3850., 5700.],
        #    'decam2014-r': [5500., 7300.],
        #    'decam2014-i': [6750., 8750.],
        #    'decam2014-z': [8200., 10200.],
        #    'decam2014-Y': [9300., 10800.],
        #    'BASS-g': [None, 5750.],
        #    'BASS-r': [None, None],
        #    'MzLS-z': [8200., 10300.],
        #}
        #
        #def trim_filter(filt):
        #    """Trim a filter."""
        #    lo, hi = 0, filt.wavelength.size
        #    if trim[filtname][0] is not None:
        #        lo = np.searchsorted(filt.wavelength, trim[filtname][0], 'left')
        #    if trim[filtname][1] is not None:
        #        hi = np.searchsorted(filt.wavelength, trim[filtname][1], 'left')
        #    filtwave = filt.wavelength[lo:hi]
        #    filtresp = filt.response[lo:hi]
        #    # response has to go to zero
        #    filtwave = np.hstack((filtwave[0]-0.1, filtwave, filtwave[-1]+0.1))
        #    filtresp = np.hstack((0., filtresp, 0.))
        #    return filters.FilterResponse(filtwave, filtresp, filt.meta)


        # If fphoto['filters'] is a dictionary, then assume that there
        # are N/S filters (as indicated by photsys).
        self.filters = {}
        for key in fphoto['filters']:
            filts = []
            for filtname in fphoto['filters'][key]:
                filt = filters.load_filter(filtname)
                #if filtname in trim.keys():
                #    filt = trim_filter(filt)
                filts.append(filt)
            self.filters[key] = filters.FilterSequence(filts)

        self.synth_filters = {}
        for key in fphoto['synth_filters']:
            self.synth_filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                              for filtname in fphoto['synth_filters'][key]])
        if 'fiber_bands' in fphoto:
            self.fiber_filters = {}
            for key in fphoto['fiber_filters']:
                self.fiber_filters[key] = filters.FilterSequence([filters.load_filter(filtname)
                                                                  for filtname in fphoto['fiber_filters'][key]])

        # absmag filters
        self.absmag_filters = filters.FilterSequence([filters.load_filter(filtname) for filtname in fphoto['absmag_filters']])

        # shifted absmag filters for use in kcorr_and_absmag
        self.filters_out = \
            filters.FilterSequence( [ f.create_shifted(band_shift=bs) for f, bs in zip(self.absmag_filters, self.band_shift) ])


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


    def synth_absmag(self, redshift, dmod, zmodelwave, zmodelflux):
        """Synthesize absolute magnitudes from the best-fitting SED.

        Parameters
        ----------
        redshift : :class:`float`
           Galaxy or QSO redshift.
        dmod : :class:`float`
           Distance modulus corresponding to `redshift`.
        zmodelwave : `numpy.ndarray`
           Observed-frame (redshifted) model wavelength array.
        zmodelflux : `numpy.ndarray`
           Observed-frame (redshifted) model spectrum.

        Returns
        -------
        synth_absmag : `numpy.ndarray`
           Absolute magnitudes based on synthesized photometry.
        synth_maggies_rest : `numpy.ndarray`
           Synthesized rest-frame photometry.

        """
        if redshift <= 0.:
            log.warning('Input redshift not defined, zero, or negative!')
            nabs = len(self.absmag_filters)
            synth_absmag = np.zeros(nabs, dtype='f8')
            synth_maggies_rest = np.zeros(nabs, dtype='f8')
            return synth_absmag, synth_maggies_rest

        # Multiply by (1+z) to convert the best-fitting model to the "rest frame".
        synth_maggies_rest = self.get_ab_maggies_unchecked(
            self.filters_out, zmodelflux * (1. + redshift) / FLUXNORM,
            zmodelwave / (1. + redshift))
        synth_absmag = -2.5 * np.log10(synth_maggies_rest) - dmod

        return synth_absmag, synth_maggies_rest


    def kcorr_and_absmag(self, nanomaggies, ivar_nanomaggies, redshift, dmod,
                         photsys, zmodelwave, zmodelflux, synth_absmag,
                         synth_maggies_rest, snrmin=2.):
        """Compute K-corrected rest-frame photometry.

        Parameters
        ----------
        nanomaggies : `numpy.ndarray`
           Input photometric fluxes in the `filters_obs` bandpasses.
        ivar_nanomaggies : `numpy.ndarray`
           Inverse variance photometry corresponding to `nanomaggies`.
        redshift : :class:`float`
           Galaxy or QSO redshift.
        dmod : :class:`float`
           Distance modulus corresponding to `redshift`.
        zmodelwave : `numpy.ndarray`
           Observed-frame (redshifted) model wavelength array.
        zmodelflux : `numpy.ndarray`
           Observed-frame (redshifted) model spectrum.
        synth_absmag : `numpy.ndarray`
           Absolute magnitudes based on synthesized photometry.
        synth_maggies_rest : `numpy.ndarray`
           Synthesized rest-frame photometry.
        snrmin : :class:`float`, defaults to 2.
           Minimum signal-to-noise ratio in the input photometry (`maggies`) in
           order for that bandpass to be used to compute a K-correction.

        Returns
        -------
        kcorr : `numpy.ndarray`
           K-corrections for each bandpass in `absmag_filters`.
        absmag : `numpy.ndarray`
           Absolute magnitudes, band-shifted according to `band_shift` (if
           provided) for each bandpass in `absmag_filters`.
        ivarabsmag : `numpy.ndarray`
           Inverse variance corresponding to `absmag`.
        synth_maggies_obs : `numpy.ndarray`
           Synthesized observed-frame photometry.

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
        nabs = len(self.absmag_filters)
        if redshift <= 0.:
            log.warning('Input redshift not defined, zero, or negative!')
            kcorr = np.zeros(nabs, dtype='f8')
            absmag = np.zeros(nabs, dtype='f8')
            ivarabsmag = np.zeros(nabs, dtype='f8')
            synth_maggies_obs = np.zeros(len(nanomaggies))
            return kcorr, absmag, ivarabsmag, synth_maggies_obs

        maggies = nanomaggies * 1e-9
        ivarmaggies = (ivar_nanomaggies / 1e-9**2) * self.bands_to_fit

        # Input bandpasses, observed frame; maggies and synth_maggies_obs
        # should be very close.
        filters_obs = self.filters[photsys]
        lambda_obs = filters_obs.effective_wavelengths.value
        lambda_out = self.filters_out.effective_wavelengths.value

        # Synthesize observed-frame photometry (should be close to maggies).
        synth_maggies_obs = self.get_ab_maggies_unchecked(
            filters_obs, zmodelflux / FLUXNORM, zmodelwave)

        # K-correct from the nearest "good" bandpass (to minimizes the K-correction)
        oband = np.empty(nabs, dtype=np.int16)
        for jj in range(nabs):
            lambdadist = np.abs(lambda_obs / (1. + redshift) - lambda_out[jj])
            oband[jj] = np.argmin(lambdadist + (maggies * np.sqrt(ivarmaggies) < snrmin) * 1e10)

        kcorr = + 2.5 * np.log10(synth_maggies_rest / synth_maggies_obs[oband])

        # m_R = M_Q + DM(z) + K_QR(z) or
        # M_Q = m_R - DM(z) - K_QR(z)
        absmag = np.copy(synth_absmag)
        ivarabsmag = np.zeros_like(absmag)

        # if we use synthesized photometry then ivarabsmag is zero
        # (which should never happen?)
        I = (maggies[oband] * np.sqrt(ivarmaggies[oband]) > snrmin)
        if np.any(I):
            C = 0.8483036976765437 # (0.4 * np.log(10.))**2
            absmag[I] = -2.5 * np.log10(maggies[oband[I]]) - dmod - kcorr[I]
            ivarabsmag[I] = maggies[oband[I]]**2 * ivarmaggies[oband[I]] * C

        return kcorr, absmag, ivarabsmag, synth_maggies_obs


    @staticmethod
    def get_ab_maggies_pre(filters, wave):
        """
        Compute preprocessing data for get_ab_maggies_unchecked() for given filter list
        and target wavelength.

        """
        # AB reference spctrum in erg/s/cm2/Hz times the speed of light in A/s
        # and converted to erg/s/cm2/A.
        abflam = 3.631e-20 * C_LIGHT * 1e13 / wave**2

        pre = []
        for ifilt, filt in enumerate(filters):
            if wave[0] > filt.wavelength[0] or filt.wavelength[-1] > wave[-1]:
                #print(filt.name, wave[0], wave[-1], filt.wavelength[0], filt.wavelength[-1])
                raise RuntimeError('Filter boundaries exceed wavelength array')

            # NB: if we assume that no padding is needed, wave extends
            # strictly past filt_wavelength, so this is safe
            lo = np.searchsorted(wave, filt.wavelength[ 0], 'right')
            hi = np.searchsorted(wave, filt.wavelength[-1], 'left') + 1
            resp = np.interp(wave[lo:hi], filt.wavelength, filt.response, left=0., right=0.) * wave[lo:hi]
            idenom = 1. / trapz(resp * abflam[lo:hi], x=wave[lo:hi])

            pre.append((lo, hi, resp, idenom))

        return tuple(pre)


    @staticmethod
    def get_ab_maggies_unchecked(filters, flux, wave, pre=None):
        """Like `get_ab_maggies()`, but by-passing the units and, more
        importantly, the padding and interpolation of wave/flux that
        speclite does.  We assume that the response function for each
        filter lies strictly within the bounds of wave, and that the
        response functions don't change so fast that we would need to
        interpolate wave to get an accurate integral.

        When wave comes from a stellar template, it has a very large
        wavelength range, so these assumptions are reasonable.  When
        wave comes from an actual camera, however, the filter
        responses are known to exceed the cameras' range of observed
        wavelengths.

        flux and wave are assumed to be in erg/s/cm2/A and A, respectively.

        """
        if pre == None:
            pre = Photometry.get_ab_maggies_pre(filters, wave)

        maggies = np.empty(len(pre))

        for ifilt, filtpre in enumerate(pre):
            lo, hi, resp, idenom = filtpre
            numer = trapz(resp * flux[lo:hi], x=wave[lo:hi])
            maggies[ifilt] = numer * idenom

        return maggies


    @staticmethod
    def get_ab_maggies(filters, flux, wave):
        """This version of get_ab_maggies() is robust to wavelength vectors
        that do not entirely cover one the response range of one or more
        filters.

        """
        try:
            if flux.ndim > 1:
                nflux = flux.shape[0]
                maggies = np.empty((nflux, len(filters)))
                for ii in range(nflux):
                    maggies[ii, :] = Photometry.get_ab_maggies_unchecked(filters, flux[ii, :], wave)
            else:
                maggies = Photometry.get_ab_maggies_unchecked(filters, flux, wave)
        except:
            # pad in case of an object at very high redshift (z > 5.5)
            log.debug('Padding input spectrum due to insufficient wavelength coverage to synthesize photometry.')
            padflux, padwave = filters.pad_spectrum(flux, wave, axis=0, method='edge')

            if flux.ndim > 1:
                nflux = padflux.shape[0]
                maggies = np.empty((nflux, len(filters)))

                for ii in range(nflux):
                    maggies[ii, :] = Photometry.get_ab_maggies_unchecked(filters, padflux[ii, :], padwave)
            else:
                maggies = Photometry.get_ab_maggies_unchecked(filters, padflux, padwave)

        return maggies


    @staticmethod
    def to_nanomaggies(maggies):
        return maggies * 1e9


    @staticmethod
    def get_photflam(maggies, lambda_eff):
        factor = 10.**(-0.4 * 48.6) * C_LIGHT * 1e13 / lambda_eff**2 # [maggies-->erg/s/cm2/A]
        return maggies * factor


    @staticmethod
    def parse_photometry(bands, maggies, lambda_eff, ivarmaggies=None,
                         nanomaggies=True, nsigma=2., min_uncertainty=None,
                         get_abmag=False):
        """Parse input (nano)maggies to various outputs and pack into a table.

        Parameters
        ----------
        flam - 10-17 erg/s/cm2/A
        fnu - 10-17 erg/s/cm2/Hz
        abmag - AB mag
        nanomaggies - input maggies are actually 1e-9 maggies

        nsigma - magnitude limit

        get_abmag - true iff table will be used by fastqa (which needs
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
        if np.any(ivarmaggies < 0.) or np.any(maggies == -99.):
            errmsg = 'All ivarmaggies must be zero or positive!'
            log.critical(errmsg)
            raise ValueError(errmsg)

        if nanomaggies:
            nanofactor = 1e-9 # [nanomaggies-->maggies]
            nmg = maggies
            nmg_ivar = ivarmaggies.copy()
        else:
            nanofactor = 1.
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
            C = 0.4 * np.log(10.)
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
            log.debug('Propagating minimum photometric uncertainties (mag): [{}]'.format(
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
        dn4000, dn4000_ivar = 0., 0.

        if rest is False or redshift is not None:
            restwave = wave / (1. + redshift) # [Angstrom]
            flam2fnu = (1. + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
        else:
            restwave = wave
            flam2fnu =  restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]

        # Require a 2-Angstrom pad around the break definition.
        wpad = 2.
        if np.min(restwave) > (3850.-wpad) or np.max(restwave) < (4100.+wpad):
            log.debug('Too little wavelength coverage to compute Dn(4000).')
            return dn4000, dn4000_ivar

        fnu = flam * flam2fnu # [erg/s/cm2/Hz]

        if flam_ivar is not None:
            fnu_ivar = flam_ivar / flam2fnu**2
        else:
            fnu_ivar = np.ones_like(flam) # uniform weights

        def _integrate(wave, flux, ivar, w1, w2):
            from scipy import integrate
            # trim for speed
            I = ((wave > w1-wpad) & (wave < w2+wpad))
            J = np.logical_and(I, ivar > 0.)
            if np.sum(I) == 0:
                return 0., 0.
            if np.sum(J) / np.sum(I) < 0.9:
                log.warning('More than 10% of pixels in Dn(4000) definition are masked.')
                return 0., 0.
            wave = wave[J]
            flux = flux[J]
            ivar = ivar[J]
            srt = np.argsort(wave)
            # should never have to extrapolate
            f1, f2 = np.interp((w1, w2), wave[srt], flux[srt])
            i1, i2 = np.interp((w1, w2), wave[srt], ivar[srt])
            # insert the boundary wavelengths then integrate
            I = ((wave > w1) & (wave < w2))
            wave = np.hstack((w1, wave[I], w2))
            flux = np.hstack((f1, flux[I], f2))
            ivar = np.hstack((i1, ivar[I], i2))
            #index_var = 1. / trapz(ivar, x=wave)
            #index = trapz(flux*ivar, x=wave) * index_var
            index_var = 1. / integrate.simpson(y=ivar, x=wave)
            index = integrate.simpson(y=flux*ivar, x=wave) * index_var
            return index, index_var

        blufactor = 3950. - 3850.
        redfactor = 4100. - 4000.

        try:
            # yes, blue wavelength go with red integral bounds
            numer, numer_var = _integrate(restwave, fnu, fnu_ivar, 4000., 4100.)
            denom, denom_var = _integrate(restwave, fnu, fnu_ivar, 3850., 3950.)
        except:
            log.warning('Integration failed when computing DN(4000).')
            return dn4000, dn4000_ivar

        if denom == 0. or numer == 0.:
            log.warning('DN(4000) is ill-defined or could not be computed.')
            return dn4000, dn4000_ivar

        dn4000 = (blufactor / redfactor) * numer / denom
        if flam_ivar is not None:
            dn4000_ivar = (1. / (dn4000**2)) / (denom_var / (denom**2) + numer_var / (numer**2))

        return dn4000, dn4000_ivar


def tractorphot_datamodel(from_file=False, datarelease='dr9'):
    """Initialize the tractorphot data model for a given Legacy Surveys data
    release.

    Args:
        from_file (bool, optional): read the datamodel from a file on-disk.
        datarelease (str, optional): data release to read; currently only `dr9`
          and `dr10` are supported.

    Returns an `astropy.table.Table` which follows the Tractor catalog
    datamodel for the given data release.

    """
    if from_file:
        from desispec.io.meta import get_desi_root_readonly
        desi_root = get_desi_root_readonly()

        datamodel_file = f'{desi_root}/external/legacysurvey/{datarelease}/south/tractor/000/tractor-0001m002.fits'
        datamodel = Table(fitsio.read(datamodel_file, rows=0, upper=True))
        for col in datamodel.colnames:
            datamodel[col] = np.zeros(datamodel[col].shape, dtype=datamodel[col].dtype)

        #for col in datamodel.colnames:
        #   print("('{}', {}, '{}'),".format(col, datamodel[col].shape, datamodel[col].dtype))
    else:
        if datarelease.lower() == 'dr9':
            COLS = [
                ('RELEASE', (1,), '>i2'),
                ('BRICKID', (1,), '>i4'),
                ('BRICKNAME', (1,), '<U8'),
                ('OBJID', (1,), '>i4'),
                ('BRICK_PRIMARY', (1,), 'bool'),
                ('MASKBITS', (1,), '>i2'),
                ('FITBITS', (1,), '>i2'),
                ('TYPE', (1,), '<U3'),
                ('RA', (1,), '>f8'),
                ('DEC', (1,), '>f8'),
                ('RA_IVAR', (1,), '>f4'),
                ('DEC_IVAR', (1,), '>f4'),
                ('BX', (1,), '>f4'),
                ('BY', (1,), '>f4'),
                ('DCHISQ', (1, 5), '>f4'),
                ('EBV', (1,), '>f4'),
                ('MJD_MIN', (1,), '>f8'),
                ('MJD_MAX', (1,), '>f8'),
                ('REF_CAT', (1,), '<U2'),
                ('REF_ID', (1,), '>i8'),
                ('PMRA', (1,), '>f4'),
                ('PMDEC', (1,), '>f4'),
                ('PARALLAX', (1,), '>f4'),
                ('PMRA_IVAR', (1,), '>f4'),
                ('PMDEC_IVAR', (1,), '>f4'),
                ('PARALLAX_IVAR', (1,), '>f4'),
                ('REF_EPOCH', (1,), '>f4'),
                ('GAIA_PHOT_G_MEAN_MAG', (1,), '>f4'),
                ('GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR', (1,), '>f4'),
                ('GAIA_PHOT_G_N_OBS', (1,), '>i2'),
                ('GAIA_PHOT_BP_MEAN_MAG', (1,), '>f4'),
                ('GAIA_PHOT_BP_MEAN_FLUX_OVER_ERROR', (1,), '>f4'),
                ('GAIA_PHOT_BP_N_OBS', (1,), '>i2'),
                ('GAIA_PHOT_RP_MEAN_MAG', (1,), '>f4'),
                ('GAIA_PHOT_RP_MEAN_FLUX_OVER_ERROR', (1,), '>f4'),
                ('GAIA_PHOT_RP_N_OBS', (1,), '>i2'),
                ('GAIA_PHOT_VARIABLE_FLAG', (1,), 'bool'),
                ('GAIA_ASTROMETRIC_EXCESS_NOISE', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_EXCESS_NOISE_SIG', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_N_OBS_AL', (1,), '>i2'),
                ('GAIA_ASTROMETRIC_N_GOOD_OBS_AL', (1,), '>i2'),
                ('GAIA_ASTROMETRIC_WEIGHT_AL', (1,), '>f4'),
                ('GAIA_DUPLICATED_SOURCE', (1,), 'bool'),
                ('GAIA_A_G_VAL', (1,), '>f4'),
                ('GAIA_E_BP_MIN_RP_VAL', (1,), '>f4'),
                ('GAIA_PHOT_BP_RP_EXCESS_FACTOR', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_SIGMA5D_MAX', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_PARAMS_SOLVED', (1,), 'uint8'),
                ('FLUX_G', (1,), '>f4'),
                ('FLUX_R', (1,), '>f4'),
                ('FLUX_Z', (1,), '>f4'),
                ('FLUX_W1', (1,), '>f4'),
                ('FLUX_W2', (1,), '>f4'),
                ('FLUX_W3', (1,), '>f4'),
                ('FLUX_W4', (1,), '>f4'),
                ('FLUX_IVAR_G', (1,), '>f4'),
                ('FLUX_IVAR_R', (1,), '>f4'),
                ('FLUX_IVAR_Z', (1,), '>f4'),
                ('FLUX_IVAR_W1', (1,), '>f4'),
                ('FLUX_IVAR_W2', (1,), '>f4'),
                ('FLUX_IVAR_W3', (1,), '>f4'),
                ('FLUX_IVAR_W4', (1,), '>f4'),
                ('FIBERFLUX_G', (1,), '>f4'),
                ('FIBERFLUX_R', (1,), '>f4'),
                ('FIBERFLUX_Z', (1,), '>f4'),
                ('FIBERTOTFLUX_G', (1,), '>f4'),
                ('FIBERTOTFLUX_R', (1,), '>f4'),
                ('FIBERTOTFLUX_Z', (1,), '>f4'),
                ('APFLUX_G', (1, 8), '>f4'),
                ('APFLUX_R', (1, 8), '>f4'),
                ('APFLUX_Z', (1, 8), '>f4'),
                ('APFLUX_RESID_G', (1, 8), '>f4'),
                ('APFLUX_RESID_R', (1, 8), '>f4'),
                ('APFLUX_RESID_Z', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_G', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_R', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_Z', (1, 8), '>f4'),
                ('APFLUX_IVAR_G', (1, 8), '>f4'),
                ('APFLUX_IVAR_R', (1, 8), '>f4'),
                ('APFLUX_IVAR_Z', (1, 8), '>f4'),
                ('APFLUX_MASKED_G', (1, 8), '>f4'),
                ('APFLUX_MASKED_R', (1, 8), '>f4'),
                ('APFLUX_MASKED_Z', (1, 8), '>f4'),
                ('APFLUX_W1', (1, 5), '>f4'),
                ('APFLUX_W2', (1, 5), '>f4'),
                ('APFLUX_W3', (1, 5), '>f4'),
                ('APFLUX_W4', (1, 5), '>f4'),
                ('APFLUX_RESID_W1', (1, 5), '>f4'),
                ('APFLUX_RESID_W2', (1, 5), '>f4'),
                ('APFLUX_RESID_W3', (1, 5), '>f4'),
                ('APFLUX_RESID_W4', (1, 5), '>f4'),
                ('APFLUX_IVAR_W1', (1, 5), '>f4'),
                ('APFLUX_IVAR_W2', (1, 5), '>f4'),
                ('APFLUX_IVAR_W3', (1, 5), '>f4'),
                ('APFLUX_IVAR_W4', (1, 5), '>f4'),
                ('MW_TRANSMISSION_G', (1,), '>f4'),
                ('MW_TRANSMISSION_R', (1,), '>f4'),
                ('MW_TRANSMISSION_Z', (1,), '>f4'),
                ('MW_TRANSMISSION_W1', (1,), '>f4'),
                ('MW_TRANSMISSION_W2', (1,), '>f4'),
                ('MW_TRANSMISSION_W3', (1,), '>f4'),
                ('MW_TRANSMISSION_W4', (1,), '>f4'),
                ('NOBS_G', (1,), '>i2'),
                ('NOBS_R', (1,), '>i2'),
                ('NOBS_Z', (1,), '>i2'),
                ('NOBS_W1', (1,), '>i2'),
                ('NOBS_W2', (1,), '>i2'),
                ('NOBS_W3', (1,), '>i2'),
                ('NOBS_W4', (1,), '>i2'),
                ('RCHISQ_G', (1,), '>f4'),
                ('RCHISQ_R', (1,), '>f4'),
                ('RCHISQ_Z', (1,), '>f4'),
                ('RCHISQ_W1', (1,), '>f4'),
                ('RCHISQ_W2', (1,), '>f4'),
                ('RCHISQ_W3', (1,), '>f4'),
                ('RCHISQ_W4', (1,), '>f4'),
                ('FRACFLUX_G', (1,), '>f4'),
                ('FRACFLUX_R', (1,), '>f4'),
                ('FRACFLUX_Z', (1,), '>f4'),
                ('FRACFLUX_W1', (1,), '>f4'),
                ('FRACFLUX_W2', (1,), '>f4'),
                ('FRACFLUX_W3', (1,), '>f4'),
                ('FRACFLUX_W4', (1,), '>f4'),
                ('FRACMASKED_G', (1,), '>f4'),
                ('FRACMASKED_R', (1,), '>f4'),
                ('FRACMASKED_Z', (1,), '>f4'),
                ('FRACIN_G', (1,), '>f4'),
                ('FRACIN_R', (1,), '>f4'),
                ('FRACIN_Z', (1,), '>f4'),
                ('ANYMASK_G', (1,), '>i2'),
                ('ANYMASK_R', (1,), '>i2'),
                ('ANYMASK_Z', (1,), '>i2'),
                ('ALLMASK_G', (1,), '>i2'),
                ('ALLMASK_R', (1,), '>i2'),
                ('ALLMASK_Z', (1,), '>i2'),
                ('WISEMASK_W1', (1,), 'uint8'),
                ('WISEMASK_W2', (1,), 'uint8'),
                ('PSFSIZE_G', (1,), '>f4'),
                ('PSFSIZE_R', (1,), '>f4'),
                ('PSFSIZE_Z', (1,), '>f4'),
                ('PSFDEPTH_G', (1,), '>f4'),
                ('PSFDEPTH_R', (1,), '>f4'),
                ('PSFDEPTH_Z', (1,), '>f4'),
                ('GALDEPTH_G', (1,), '>f4'),
                ('GALDEPTH_R', (1,), '>f4'),
                ('GALDEPTH_Z', (1,), '>f4'),
                ('NEA_G', (1,), '>f4'),
                ('NEA_R', (1,), '>f4'),
                ('NEA_Z', (1,), '>f4'),
                ('BLOB_NEA_G', (1,), '>f4'),
                ('BLOB_NEA_R', (1,), '>f4'),
                ('BLOB_NEA_Z', (1,), '>f4'),
                ('PSFDEPTH_W1', (1,), '>f4'),
                ('PSFDEPTH_W2', (1,), '>f4'),
                ('PSFDEPTH_W3', (1,), '>f4'),
                ('PSFDEPTH_W4', (1,), '>f4'),
                ('WISE_COADD_ID', (1,), '<U8'),
                ('WISE_X', (1,), '>f4'),
                ('WISE_Y', (1,), '>f4'),
                ('LC_FLUX_W1', (1, 15), '>f4'),
                ('LC_FLUX_W2', (1, 15), '>f4'),
                ('LC_FLUX_IVAR_W1', (1, 15), '>f4'),
                ('LC_FLUX_IVAR_W2', (1, 15), '>f4'),
                ('LC_NOBS_W1', (1, 15), '>i2'),
                ('LC_NOBS_W2', (1, 15), '>i2'),
                ('LC_FRACFLUX_W1', (1, 15), '>f4'),
                ('LC_FRACFLUX_W2', (1, 15), '>f4'),
                ('LC_RCHISQ_W1', (1, 15), '>f4'),
                ('LC_RCHISQ_W2', (1, 15), '>f4'),
                ('LC_MJD_W1', (1, 15), '>f8'),
                ('LC_MJD_W2', (1, 15), '>f8'),
                ('LC_EPOCH_INDEX_W1', (1, 15), '>i2'),
                ('LC_EPOCH_INDEX_W2', (1, 15), '>i2'),
                ('SERSIC', (1,), '>f4'),
                ('SERSIC_IVAR', (1,), '>f4'),
                ('SHAPE_R', (1,), '>f4'),
                ('SHAPE_R_IVAR', (1,), '>f4'),
                ('SHAPE_E1', (1,), '>f4'),
                ('SHAPE_E1_IVAR', (1,), '>f4'),
                ('SHAPE_E2', (1,), '>f4'),
                ('SHAPE_E2_IVAR', (1,), '>f4'),
                # added columns
                ('LS_ID', (1,), '>i8'),
                ('TARGETID', (1,), '>i8'),
                ]
        elif datarelease.lower() == 'dr10':
            COLS = [
                ('RELEASE', (1,), '>i2'),
                ('BRICKID', (1,), '>i4'),
                ('BRICKNAME', (1,), '<U8'),
                ('OBJID', (1,), '>i4'),
                ('BRICK_PRIMARY', (1,), 'bool'),
                ('MASKBITS', (1,), '>i2'),
                ('FITBITS', (1,), '>i2'),
                ('TYPE', (1,), '<U3'),
                ('RA', (1,), '>f8'),
                ('DEC', (1,), '>f8'),
                ('RA_IVAR', (1,), '>f4'),
                ('DEC_IVAR', (1,), '>f4'),
                ('BX', (1,), '>f4'),
                ('BY', (1,), '>f4'),
                ('DCHISQ', (1, 5), '>f4'),
                ('EBV', (1,), '>f4'),
                ('MJD_MIN', (1,), '>f8'),
                ('MJD_MAX', (1,), '>f8'),
                ('REF_CAT', (1,), '<U2'),
                ('REF_ID', (1,), '>i8'),
                ('PMRA', (1,), '>f4'),
                ('PMDEC', (1,), '>f4'),
                ('PARALLAX', (1,), '>f4'),
                ('PMRA_IVAR', (1,), '>f4'),
                ('PMDEC_IVAR', (1,), '>f4'),
                ('PARALLAX_IVAR', (1,), '>f4'),
                ('REF_EPOCH', (1,), '>f4'),
                ('GAIA_PHOT_G_MEAN_MAG', (1,), '>f4'),
                ('GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR', (1,), '>f4'),
                ('GAIA_PHOT_G_N_OBS', (1,), '>i2'),
                ('GAIA_PHOT_BP_MEAN_MAG', (1,), '>f4'),
                ('GAIA_PHOT_BP_MEAN_FLUX_OVER_ERROR', (1,), '>f4'),
                ('GAIA_PHOT_BP_N_OBS', (1,), '>i2'),
                ('GAIA_PHOT_RP_MEAN_MAG', (1,), '>f4'),
                ('GAIA_PHOT_RP_MEAN_FLUX_OVER_ERROR', (1,), '>f4'),
                ('GAIA_PHOT_RP_N_OBS', (1,), '>i2'),
                ('GAIA_PHOT_VARIABLE_FLAG', (1,), 'bool'),
                ('GAIA_ASTROMETRIC_EXCESS_NOISE', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_EXCESS_NOISE_SIG', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_N_OBS_AL', (1,), '>i2'),
                ('GAIA_ASTROMETRIC_N_GOOD_OBS_AL', (1,), '>i2'),
                ('GAIA_ASTROMETRIC_WEIGHT_AL', (1,), '>f4'),
                ('GAIA_DUPLICATED_SOURCE', (1,), 'bool'),
                ('GAIA_A_G_VAL', (1,), '>f4'),
                ('GAIA_E_BP_MIN_RP_VAL', (1,), '>f4'),
                ('GAIA_PHOT_BP_RP_EXCESS_FACTOR', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_SIGMA5D_MAX', (1,), '>f4'),
                ('GAIA_ASTROMETRIC_PARAMS_SOLVED', (1,), 'uint8'),
                ('FLUX_G', (1,), '>f4'),
                ('FLUX_R', (1,), '>f4'),
                ('FLUX_I', (1,), '>f4'),
                ('FLUX_Z', (1,), '>f4'),
                ('FLUX_W1', (1,), '>f4'),
                ('FLUX_W2', (1,), '>f4'),
                ('FLUX_W3', (1,), '>f4'),
                ('FLUX_W4', (1,), '>f4'),
                ('FLUX_IVAR_G', (1,), '>f4'),
                ('FLUX_IVAR_R', (1,), '>f4'),
                ('FLUX_IVAR_I', (1,), '>f4'),
                ('FLUX_IVAR_Z', (1,), '>f4'),
                ('FLUX_IVAR_W1', (1,), '>f4'),
                ('FLUX_IVAR_W2', (1,), '>f4'),
                ('FLUX_IVAR_W3', (1,), '>f4'),
                ('FLUX_IVAR_W4', (1,), '>f4'),
                ('FIBERFLUX_G', (1,), '>f4'),
                ('FIBERFLUX_R', (1,), '>f4'),
                ('FIBERFLUX_I', (1,), '>f4'),
                ('FIBERFLUX_Z', (1,), '>f4'),
                ('FIBERTOTFLUX_G', (1,), '>f4'),
                ('FIBERTOTFLUX_R', (1,), '>f4'),
                ('FIBERTOTFLUX_I', (1,), '>f4'),
                ('FIBERTOTFLUX_Z', (1,), '>f4'),
                ('APFLUX_G', (1, 8), '>f4'),
                ('APFLUX_R', (1, 8), '>f4'),
                ('APFLUX_I', (1, 8), '>f4'),
                ('APFLUX_Z', (1, 8), '>f4'),
                ('APFLUX_RESID_G', (1, 8), '>f4'),
                ('APFLUX_RESID_R', (1, 8), '>f4'),
                ('APFLUX_RESID_I', (1, 8), '>f4'),
                ('APFLUX_RESID_Z', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_G', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_R', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_I', (1, 8), '>f4'),
                ('APFLUX_BLOBRESID_Z', (1, 8), '>f4'),
                ('APFLUX_IVAR_G', (1, 8), '>f4'),
                ('APFLUX_IVAR_R', (1, 8), '>f4'),
                ('APFLUX_IVAR_I', (1, 8), '>f4'),
                ('APFLUX_IVAR_Z', (1, 8), '>f4'),
                ('APFLUX_MASKED_G', (1, 8), '>f4'),
                ('APFLUX_MASKED_R', (1, 8), '>f4'),
                ('APFLUX_MASKED_I', (1, 8), '>f4'),
                ('APFLUX_MASKED_Z', (1, 8), '>f4'),
                ('APFLUX_W1', (1, 5), '>f4'),
                ('APFLUX_W2', (1, 5), '>f4'),
                ('APFLUX_W3', (1, 5), '>f4'),
                ('APFLUX_W4', (1, 5), '>f4'),
                ('APFLUX_RESID_W1', (1, 5), '>f4'),
                ('APFLUX_RESID_W2', (1, 5), '>f4'),
                ('APFLUX_RESID_W3', (1, 5), '>f4'),
                ('APFLUX_RESID_W4', (1, 5), '>f4'),
                ('APFLUX_IVAR_W1', (1, 5), '>f4'),
                ('APFLUX_IVAR_W2', (1, 5), '>f4'),
                ('APFLUX_IVAR_W3', (1, 5), '>f4'),
                ('APFLUX_IVAR_W4', (1, 5), '>f4'),
                ('MW_TRANSMISSION_G', (1,), '>f4'),
                ('MW_TRANSMISSION_R', (1,), '>f4'),
                ('MW_TRANSMISSION_I', (1,), '>f4'),
                ('MW_TRANSMISSION_Z', (1,), '>f4'),
                ('MW_TRANSMISSION_W1', (1,), '>f4'),
                ('MW_TRANSMISSION_W2', (1,), '>f4'),
                ('MW_TRANSMISSION_W3', (1,), '>f4'),
                ('MW_TRANSMISSION_W4', (1,), '>f4'),
                ('NOBS_G', (1,), '>i2'),
                ('NOBS_R', (1,), '>i2'),
                ('NOBS_I', (1,), '>i2'),
                ('NOBS_Z', (1,), '>i2'),
                ('NOBS_W1', (1,), '>i2'),
                ('NOBS_W2', (1,), '>i2'),
                ('NOBS_W3', (1,), '>i2'),
                ('NOBS_W4', (1,), '>i2'),
                ('RCHISQ_G', (1,), '>f4'),
                ('RCHISQ_R', (1,), '>f4'),
                ('RCHISQ_I', (1,), '>f4'),
                ('RCHISQ_Z', (1,), '>f4'),
                ('RCHISQ_W1', (1,), '>f4'),
                ('RCHISQ_W2', (1,), '>f4'),
                ('RCHISQ_W3', (1,), '>f4'),
                ('RCHISQ_W4', (1,), '>f4'),
                ('FRACFLUX_G', (1,), '>f4'),
                ('FRACFLUX_R', (1,), '>f4'),
                ('FRACFLUX_I', (1,), '>f4'),
                ('FRACFLUX_Z', (1,), '>f4'),
                ('FRACFLUX_W1', (1,), '>f4'),
                ('FRACFLUX_W2', (1,), '>f4'),
                ('FRACFLUX_W3', (1,), '>f4'),
                ('FRACFLUX_W4', (1,), '>f4'),
                ('FRACMASKED_G', (1,), '>f4'),
                ('FRACMASKED_R', (1,), '>f4'),
                ('FRACMASKED_I', (1,), '>f4'),
                ('FRACMASKED_Z', (1,), '>f4'),
                ('FRACIN_G', (1,), '>f4'),
                ('FRACIN_R', (1,), '>f4'),
                ('FRACIN_I', (1,), '>f4'),
                ('FRACIN_Z', (1,), '>f4'),
                ('NGOOD_G', (1,), '>i2'),
                ('NGOOD_R', (1,), '>i2'),
                ('NGOOD_I', (1,), '>i2'),
                ('NGOOD_Z', (1,), '>i2'),
                ('ANYMASK_G', (1,), '>i2'),
                ('ANYMASK_R', (1,), '>i2'),
                ('ANYMASK_I', (1,), '>i2'),
                ('ANYMASK_Z', (1,), '>i2'),
                ('ALLMASK_G', (1,), '>i2'),
                ('ALLMASK_R', (1,), '>i2'),
                ('ALLMASK_I', (1,), '>i2'),
                ('ALLMASK_Z', (1,), '>i2'),
                ('WISEMASK_W1', (1,), 'uint8'),
                ('WISEMASK_W2', (1,), 'uint8'),
                ('PSFSIZE_G', (1,), '>f4'),
                ('PSFSIZE_R', (1,), '>f4'),
                ('PSFSIZE_I', (1,), '>f4'),
                ('PSFSIZE_Z', (1,), '>f4'),
                ('PSFDEPTH_G', (1,), '>f4'),
                ('PSFDEPTH_R', (1,), '>f4'),
                ('PSFDEPTH_I', (1,), '>f4'),
                ('PSFDEPTH_Z', (1,), '>f4'),
                ('GALDEPTH_G', (1,), '>f4'),
                ('GALDEPTH_R', (1,), '>f4'),
                ('GALDEPTH_I', (1,), '>f4'),
                ('GALDEPTH_Z', (1,), '>f4'),
                ('NEA_G', (1,), '>f4'),
                ('NEA_R', (1,), '>f4'),
                ('NEA_I', (1,), '>f4'),
                ('NEA_Z', (1,), '>f4'),
                ('BLOB_NEA_G', (1,), '>f4'),
                ('BLOB_NEA_R', (1,), '>f4'),
                ('BLOB_NEA_I', (1,), '>f4'),
                ('BLOB_NEA_Z', (1,), '>f4'),
                ('PSFDEPTH_W1', (1,), '>f4'),
                ('PSFDEPTH_W2', (1,), '>f4'),
                ('PSFDEPTH_W3', (1,), '>f4'),
                ('PSFDEPTH_W4', (1,), '>f4'),
                ('WISE_COADD_ID', (1,), '<U8'),
                ('WISE_X', (1,), '>f4'),
                ('WISE_Y', (1,), '>f4'),
                ('LC_FLUX_W1', (1, 17), '>f4'),
                ('LC_FLUX_W2', (1, 17), '>f4'),
                ('LC_FLUX_IVAR_W1', (1, 17), '>f4'),
                ('LC_FLUX_IVAR_W2', (1, 17), '>f4'),
                ('LC_NOBS_W1', (1, 17), '>i2'),
                ('LC_NOBS_W2', (1, 17), '>i2'),
                ('LC_FRACFLUX_W1', (1, 17), '>f4'),
                ('LC_FRACFLUX_W2', (1, 17), '>f4'),
                ('LC_RCHISQ_W1', (1, 17), '>f4'),
                ('LC_RCHISQ_W2', (1, 17), '>f4'),
                ('LC_MJD_W1', (1, 17), '>f8'),
                ('LC_MJD_W2', (1, 17), '>f8'),
                ('LC_EPOCH_INDEX_W1', (1, 17), '>i2'),
                ('LC_EPOCH_INDEX_W2', (1, 17), '>i2'),
                ('SERSIC', (1,), '>f4'),
                ('SERSIC_IVAR', (1,), '>f4'),
                ('SHAPE_R', (1,), '>f4'),
                ('SHAPE_R_IVAR', (1,), '>f4'),
                ('SHAPE_E1', (1,), '>f4'),
                ('SHAPE_E1_IVAR', (1,), '>f4'),
                ('SHAPE_E2', (1,), '>f4'),
                ('SHAPE_E2_IVAR', (1,), '>f4'),
                # added columns
                ('LS_ID', (1,), '>i8'),
                ('TARGETID', (1,), '>i8'),
                ]
        else:
            errmsg = f'Unrecognized data release {datarelease}; only dr9 and dr10 currently supported.'
            log.critical(errmsg)
            raise IOError(errmsg)

        datamodel = Table()
        for col in COLS:
            datamodel[col[0]] = np.zeros(shape=col[1], dtype=col[2])

    return datamodel


def _gather_tractorphot_onebrick(input_cat, legacysurveydir, radius_match, racolumn, deccolumn,
                                 datamodel, restrict_region):
    """Support routine for gather_tractorphot.

    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from desitarget import geomask
    from desitarget.io import desitarget_resolve_dec

    assert(np.all(input_cat['BRICKNAME'] == input_cat['BRICKNAME'][0]))
    brick = input_cat['BRICKNAME'][0]

    idr9 = np.where((input_cat['RELEASE'] > 0) * (input_cat['BRICKID'] > 0) *
                    (input_cat['BRICK_OBJID'] > 0) * (input_cat['PHOTSYS'] != ''))[0]
    ipos = np.delete(np.arange(len(input_cat)), idr9)

    out = Table(np.hstack(np.repeat(datamodel, len(np.atleast_1d(input_cat)))))
    out['TARGETID'] = input_cat['TARGETID']

    # DR9 targeting photometry exists
    if len(idr9) > 0:
        assert(np.all(input_cat['PHOTSYS'][idr9] == input_cat['PHOTSYS'][idr9][0]))

        # find the catalog
        photsys = input_cat['PHOTSYS'][idr9][0]

        if photsys == 'S':
            region = 'south'
        elif photsys == 'N':
            region = 'north'

        #raslice = np.array(['{:06d}'.format(int(ra*1000))[:3] for ra in input_cat['RA']])
        tractorfile = os.path.join(legacysurveydir, region, 'tractor', brick[:3], f'tractor-{brick}.fits')

        if not os.path.isfile(tractorfile):
            errmsg = f'Unable to find Tractor catalog {tractorfile}'
            log.critical(errmsg)
            raise IOError(errmsg)

        # Some commissioning and SV targets can have brick_primary==False, so don't require it here.
        #<Table length=1>
        #     TARGETID     BRICKNAME BRICKID BRICK_OBJID RELEASE CMX_TARGET DESI_TARGET   SV1_DESI_TARGET   SV2_DESI_TARGET SV3_DESI_TARGET SCND_TARGET
        #      int64          str8    int32     int32     int16    int64       int64           int64             int64           int64         int64
        #----------------- --------- ------- ----------- ------- ---------- ----------- ------------------- --------------- --------------- -----------
        #39628509856927757  0352p315  503252        4109    9010          0           0 2305843009213693952               0               0           0
        #<Table length=1>
        #     TARGETID         TARGET_RA          TARGET_DEC     TILEID SURVEY PROGRAM
        #      int64            float64            float64       int32   str7    str6
        #----------------- ------------------ ------------------ ------ ------ -------
        #39628509856927757 35.333944142134406 31.496490061792002  80611    sv1  bright

        _tractor = fitsio.read(tractorfile, columns=['OBJID', 'BRICK_PRIMARY'], upper=True)
        #I = np.where(_tractor['BRICK_PRIMARY'] * np.isin(_tractor['OBJID'], input_cat['BRICK_OBJID']))[0]
        I = np.where(np.isin(_tractor['OBJID'], input_cat['BRICK_OBJID'][idr9]))[0]

        ## Some secondary programs have BRICKNAME!='' and BRICK_OBJID==0 (i.e.,
        ## not populated). However, there should always be a match here because
        ## we "repair" brick_objid in the main function.
        #if len(I) == 0:
        #    return Table()

        tractor_dr9 = Table(fitsio.read(tractorfile, rows=I, upper=True))

        # sort explicitly in order to ensure order
        srt = geomask.match_to(tractor_dr9['OBJID'], input_cat['BRICK_OBJID'][idr9])
        tractor_dr9 = tractor_dr9[srt]
        assert(np.all((tractor_dr9['BRICKID'] == input_cat['BRICKID'][idr9])*(tractor_dr9['OBJID'] == input_cat['BRICK_OBJID'][idr9])))

        tractor_dr9['LS_ID'] = np.int64(0) # will be filled in at the end
        tractor_dr9['TARGETID'] = input_cat['TARGETID'][idr9]

        out[idr9] = tractor_dr9
        del tractor_dr9

    # use positional matching
    if len(ipos) > 0:
        rad = radius_match * u.arcsec

        # resolve north/south unless restrict region is set
        if restrict_region is not None:
            tractorfile = os.path.join(legacysurveydir, restrict_region, 'tractor', brick[:3], f'tractor-{brick}.fits')
            if not os.path.isfile(tractorfile):
                return out
        else:
            tractorfile_north = os.path.join(legacysurveydir, 'north', 'tractor', brick[:3], f'tractor-{brick}.fits')
            tractorfile_south = os.path.join(legacysurveydir, 'south', 'tractor', brick[:3], f'tractor-{brick}.fits')
            if os.path.isfile(tractorfile_north) and not os.path.isfile(tractorfile_south):
                tractorfile = tractorfile_north
            elif not os.path.isfile(tractorfile_north) and os.path.isfile(tractorfile_south):
                tractorfile = tractorfile_south
            elif os.path.isfile(tractorfile_north) and os.path.isfile(tractorfile_south):
                if np.median(input_cat[deccolumn][ipos]) < desitarget_resolve_dec():
                    tractorfile = tractorfile_south
                else:
                    tractorfile = tractorfile_north
            elif not os.path.isfile(tractorfile_north) and not os.path.isfile(tractorfile_south):
                return out

        _tractor = fitsio.read(tractorfile, columns=['RA', 'DEC', 'BRICK_PRIMARY'], upper=True)
        iprimary = np.where(_tractor['BRICK_PRIMARY'])[0] # only primary targets
        if len(iprimary) == 0:
            log.warning(f'No primary photometric targets on brick {brick}.')
        else:
            _tractor = _tractor[iprimary]
            coord_tractor = SkyCoord(ra=_tractor['RA']*u.deg, dec=_tractor['DEC']*u.deg)
            # Some targets can appear twice (with different targetids), so
            # to make sure we do it right, we have to loop. Example:
            #
            #     TARGETID    SURVEY PROGRAM     TARGET_RA          TARGET_DEC    OBJID BRICKID RELEASE  SKY  GAIADR    RA     DEC   GROUP BRICKNAME
            #      int64       str7    str6       float64            float64      int64  int64   int64  int64 int64  float64 float64 int64    str8
            # --------------- ------ ------- ------------------ ----------------- ----- ------- ------- ----- ------ ------- ------- ----- ---------
            # 234545047666699    sv1   other 150.31145983340912 2.587887211205909    11  345369      53     0      0     0.0     0.0     0  1503p025
            # 243341140688909    sv1   other 150.31145983340912 2.587887211205909    13  345369      55     0      0     0.0     0.0     0  1503p025

            for indx_cat, (ra, dec, targetid) in enumerate(zip(input_cat[racolumn][ipos],
                                                               input_cat[deccolumn][ipos],
                                                               input_cat['TARGETID'][ipos])):
                coord_cat = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
                indx_tractor, d2d, _ = coord_cat.match_to_catalog_sky(coord_tractor)
                if d2d < rad:
                    _tractor = Table(fitsio.read(tractorfile, rows=iprimary[indx_tractor], upper=True))
                    _tractor['LS_ID'] = np.int64(0) # will be filled in at the end
                    _tractor['TARGETID'] = targetid
                    out[ipos[indx_cat]] = _tractor[0]

    # Add a unique DR9 identifier.
    out['LS_ID'] = (out['RELEASE'].astype(np.int64) << 40) | (out['BRICKID'].astype(np.int64) << 16) | (out['OBJID'].astype(np.int64))

    assert(np.all(input_cat['TARGETID'] == out['TARGETID']))

    return out


def gather_tractorphot(input_cat, racolumn='TARGET_RA', deccolumn='TARGET_DEC',
                       legacysurveydir=None, dr9dir=None, radius_match=1.0,
                       restrict_region=None, columns=None):
    """Retrieve the Tractor catalog for all the objects in this catalog (one brick).

    Args:
        input_cat (astropy.table.Table): input table with the following
          (required) columns: TARGETID, TARGET_RA, TARGET_DEC. Additional
          optional columns that will ensure proper matching are BRICKNAME,
          RELEASE, PHOTSYS, BRICKID, and BRICK_OBJID.
        legacysurveydir (str): full path to the location of the Tractor catalogs
        dr9dir (str): relegated keyword; please use `legacysurveydir`
        radius_match (float, arcsec): matching radius (default, 1 arcsec)
        restrict_region (str): when positional matching, restrict the region to
          check for photometry, otherwise check both 'north' and 'south'
          (default None)
        columns (str array): return this subset of columns

    Returns a table of Tractor photometry. Matches are identified either using
    BRICKID and BRICK_OBJID or using positional matching (1 arcsec radius).

    """
    from desitarget.targets import decode_targetid
    from desiutil.brick import brickname

    if len(input_cat) == 0:
        log.warning('No objects in input catalog.')
        return Table()

    for col in ['TARGETID', racolumn, deccolumn]:
        if col not in input_cat.colnames:
            errmsg = f'Missing required input column {col}'
            log.critical(errmsg)
            raise ValueError(errmsg)

    # If these columns don't exist, add them with blank entries:
    COLS = [('RELEASE', (1,), '>i2'), ('BRICKID', (1,), '>i4'),
            ('BRICKNAME', (1,), '<U8'), ('BRICK_OBJID', (1,), '>i4'),
            ('PHOTSYS', (1,), '<U1')]
    for col in COLS:
        if col[0] not in input_cat.colnames:
            input_cat[col[0]] = np.zeros(col[1], dtype=col[2])

    if dr9dir is not None:
        log.warning('Keyword dr9dir is relegated; please use legacysurveydir.')
        legacysurveydir = dr9dir

    if legacysurveydir is None:
        from desispec.io.meta import get_desi_root_readonly
        desi_root = get_desi_root_readonly()
        legacysurveydir = os.path.join(desi_root, 'external', 'legacysurvey', 'dr9')

    if not os.path.isdir(legacysurveydir):
        errmsg = f'Legacy Surveys directory {legacysurveydir} not found.'
        log.critical(errmsg)
        raise IOError(errmsg)

    if 'dr9' in legacysurveydir:
        datarelease = 'dr9'
    elif 'dr10' in legacysurveydir:
        datarelease = 'dr10'
    else:
        errmsg = f'Unable to determine data release from {legacysurveydir}; falling back to DR9.'
        log.warning(errmsg)
        datarelease = 'dr9'

    if restrict_region is not None:
        if restrict_region not in ['north', 'south']:
            errmsg = f"Optional input restrict_region must be either 'north' or 'south'."
            log.critical(errmsg)
            raise ValueError(errmsg)

    ## Some secondary programs (e.g., 39632961435338613, 39632966921487347)
    ## have BRICKNAME!='' & BRICKID!=0, but BRICK_OBJID==0. Unpack those here
    ## using decode_targetid.
    #idecode = np.where((input_cat['BRICKNAME'] != '') * (input_cat['BRICK_OBJID'] == 0))[0]
    #if len(idecode) > 0:
    #    log.debug('Inferring BRICK_OBJID for {} objects using decode_targetid'.format(len(idecode)))
    #    new_objid, new_brickid, _, _, _, _ = decode_targetid(input_cat['TARGETID'][idecode])
    #    assert(np.all(new_brickid == input_cat['BRICKID'][idecode]))
    #    input_cat['BRICK_OBJID'][idecode] = new_objid

    # BRICKNAME can sometimes be blank; fix that here. NB: this step has to come
    # *after* the decode step, above!
    inobrickname = np.where(input_cat['BRICKNAME'] == '')[0]
    if len(inobrickname) > 0:
        log.debug(f'Inferring brickname for {len(inobrickname):,d} objects')
        input_cat['BRICKNAME'][inobrickname] = brickname(input_cat[racolumn][inobrickname],
                                                         input_cat[deccolumn][inobrickname])

    # Split into unique brickname(s) and initialize the data model.
    bricknames = input_cat['BRICKNAME']
    datamodel = tractorphot_datamodel(datarelease=datarelease)

    out = Table(np.hstack(np.repeat(datamodel, len(np.atleast_1d(input_cat)))))
    for onebrickname in set(bricknames):
        I = np.where(onebrickname == bricknames)[0]
        out[I] = _gather_tractorphot_onebrick(input_cat[I], legacysurveydir, radius_match, racolumn,
                                              deccolumn, datamodel, restrict_region)

    if 'RELEASE' in input_cat.colnames:
        _, _, check_release, _, _, _ = decode_targetid(input_cat['TARGETID'])
        bug = np.where(out['RELEASE'] != check_release)[0]
        if len(bug) > 0:
            input_cat['BRICKNAME'][bug] = brickname(input_cat[racolumn][bug], input_cat[deccolumn][bug])
            input_cat['RELEASE'][bug] = 0
            input_cat['BRICKID'][bug] = 0
            input_cat['BRICK_OBJID'][bug] = 0
            input_cat['PHOTSYS'][bug] = ''

            bugout = Table(np.hstack(np.repeat(datamodel, len(bug))))
            for onebrickname in set(input_cat['BRICKNAME'][bug]):
                I = np.where(onebrickname == input_cat['BRICKNAME'][bug])[0]
                bugout[I] = _gather_tractorphot_onebrick(input_cat[bug][I], legacysurveydir, radius_match, racolumn, deccolumn,
                                                         datamodel, restrict_region)

    if columns is not None:
        if type(columns) is not list:
            columns = columns.tolist()
        out = out[columns]

    return out
