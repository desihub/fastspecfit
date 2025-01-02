"""
fastspecfit.templates
=====================

Tools for handling templates.

"""
import os
import numpy as np
import numba

import fitsio
from astropy.table import Table

from fastspecfit.logger import log

class Templates(object):

    # SPS template constants (used by build-templates)
    # https://github.com/cconroy20/fsps/tree/master/SPECTRA/C3K#readme
    PIXKMS = 25.  # [km/s]
    PIXKMS_BOUNDS = (2750., 9100.)

    AGN_PIXKMS = 75.  # [km/s]
    AGN_PIXKMS_BOUNDS = (1075., 3090.)

    DEFAULT_TEMPLATEVERSION = '2.0.0'
    DEFAULT_IMF = 'chabrier'

    # highest vdisp for which we attempt to use cached FFTs
    MAX_PRE_VDISP = 500.

    def __init__(self, template_file=None, template_version=None, imf=None,
                 mintemplatewave=None, maxtemplatewave=40e4, vdisp_nominal=250.,
                 fastphot=False, read_linefluxes=False):
        """"
        Read the templates into a dictionary.

        Parameters
        ----------
        template_file : :class:`str`
            Full path to the templates to read. Defaults to the output of
            :class:`get_templates_filename`.
        template_version : :class:`str`
            Version of the templates.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.
        mintemplatewave : :class:`float`, optional, defaults to None
            Minimum template wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxtemplatewave : :class:`float`, optional, defaults to 6e4
            Maximum template wavelength to read into memory.

        """
        self.init_ffts()

        if template_file is None:
            if template_version is None:
                template_version = Templates.DEFAULT_TEMPLATEVERSION
            if imf is None:
                imf = Templates.DEFAULT_IMF
            template_file = self.get_templates_filename(template_version=template_version, imf=imf)

        if not os.path.isfile(template_file):
            errmsg = f'Templates file {template_file} not found.'
            log.critical(errmsg)
            raise IOError(errmsg)

        self.file = template_file

        T = fitsio.FITS(template_file)
        templatewave     = T['WAVE'].read()        # [npix]
        wavehdr          = T['WAVE'].read_header() # [npix]
        templateflux     = T['FLUX'].read()        # [npix x nsed]
        templatelineflux = T['LINEFLUX'].read()    # [npix x nsed]
        templateinfo     = T['METADATA'].read()
        templatehdr      = T['METADATA'].read_header()

        templateflux     = np.transpose(templateflux).copy()
        templatelineflux = np.transpose(templatelineflux).copy()

        self.version = T[0].read_header()['VERSION']

        self.imf = templatehdr['IMF']
        self.ntemplates = len(templateinfo)

        if mintemplatewave is None:
            mintemplatewave = np.min(templatewave)

        keeplo = np.searchsorted(templatewave, mintemplatewave, 'left')
        keephi = np.searchsorted(templatewave, maxtemplatewave, 'right')
        self.wave = templatewave[keeplo:keephi]
        self.flux = templateflux[:, keeplo:keephi]
        self.flux_nolines = self.flux - templatelineflux[:, keeplo:keephi]
        self.npix = len(self.wave)

        # dust attenuation curve
        self.dust_klambda = Templates.klambda(self.wave)
        self.vdisp_nominal = vdisp_nominal # [km/s]

        pixkms_bounds = np.searchsorted(self.wave, Templates.PIXKMS_BOUNDS, 'left')
        self.pixkms_bounds = pixkms_bounds

        self.conv_pre = self.convolve_vdisp_pre(self.flux)
        self.flux_nomvdisp = self.convolve_vdisp(self.flux, vdisp_nominal)

        self.conv_pre_nolines = self.convolve_vdisp_pre(self.flux_nolines)
        self.flux_nolines_nomvdisp = self.convolve_vdisp(self.flux_nolines, vdisp_nominal)

        self.info = Table(templateinfo)

        if 'DUSTFLUX' in T and 'AGNFLUX' in T:
            from fastspecfit.util import trapz

            # make sure fluxes are normalized to unity
            dustflux = T['DUSTFLUX'].read()
            #dustflux /= trapz(dustflux, x=templatewave) # should already be 1.0
            self.dustflux = dustflux[keeplo:keephi]

            #dusthdr = T['DUSTFLUX'].read_header()
            #self.qpah     = dusthdr['QPAH']
            #self.umin     = dusthdr['UMIN']
            #self.gamma    = dusthdr['GAMMA']

            # construct the AGN wavelength vector
            iragnflux = T['AGNFLUX'].read()
            iragnwave = T['AGNWAVE'].read()
            #iragnflux  /= trapz(iragnflux, x=iragnwave) # should already be 1.0

            trim = np.searchsorted(iragnwave, 1e4, 'left') # hack...
            iragnflux = iragnflux[trim:]
            iragnwave = iragnwave[trim:]

            feflux = T['FEFLUX'].read()
            fewave = T['FEWAVE'].read()

            febounds = np.searchsorted(templatewave, Templates.AGN_PIXKMS_BOUNDS, 'left')
            irbounds = np.searchsorted(templatewave, iragnwave[0], 'left')

            agnwave = np.hstack((templatewave[:febounds[0]], fewave,
                                 templatewave[febounds[1]:irbounds],
                                 iragnwave))
            #self.agnwave = agnwave

            #agnhdr = T['AGNFLUX'].read_header()
            #self.agntau   = agnhdr['AGNTAU']
        else:
            errmsg = f'Templates file {template_file} missing mandatory extensions DUSTFLUX and AGNFLUX.'
            log.critical(errmsg)
            raise IOError(errmsg)

        # Read the model emission-line fluxes; only present for
        # template_version>=1.1.1 and generally only useful to a power-user.
        if read_linefluxes:
            self.linewaves  = T['LINEWAVES'].read()
            self.linefluxes = T['LINEFLUXES'].read()


    def init_ffts(self):
        """
        Use MKL's accelerated FFT backend in SciPy if available.
        If other FFT libraries are acceptable, we can also
        check for them here.

        Set the appropriate convolution function to call
        for non-cached convolutions in convolve_vdisp(),
        based on what we believe to be fastest.

        """
        import scipy.fft as sc_fft
        import scipy.signal as sc_sig
        from importlib.util import find_spec

        if find_spec("mkl_fft") is not None:
            import mkl_fft._scipy_fft_backend as be
            sc_fft.set_global_backend(be)

            self.convolve = sc_sig.convolve
            log.debug('Using mkl_fft library for FFTs')
        else:
            self.convolve = sc_sig.oaconvolve


    @staticmethod
    def get_templates_filename(template_version, imf):
        """Get the templates filename.

        """
        template_dir = os.path.expandvars(os.environ.get('FTEMPLATES_DIR'))
        template_file = os.path.join(template_dir, template_version, f'ftemplates-{imf}-{template_version}.fits')
        return template_file


    def convolve_vdisp_pre(self, templateflux):
        """Given a 2D array of of one or more template fluxes, return a
        preprocessing structure to accelerate future vdisp convolutions via the
        FFT. This structure is unpacked in
        ContinuumTools.build_stellar_continuum()

        Parameters
        ----------
        templateflux: :class:`np.ndarray`
            [ntemplates x nwavelengths] array of templates

        Returns
        -------
        preprocessing structure

        """
        import scipy.fft as sc_fft

        # determine largest kernel we will support
        # based on the maximum supported vdisp.
        pixsize_kms = Templates.PIXKMS
        sigma = Templates.MAX_PRE_VDISP / pixsize_kms # [pixels]
        radius = Templates._gaussian_radius(sigma)
        kernel_size = 2*radius + 1

        lo, hi = self.pixkms_bounds

        # extract the un-convolved ranges of templateflux
        flux_lo = templateflux[:, :lo]
        flux_hi = templateflux[:, hi:]
        flux_mid = templateflux[:, lo:hi]

        mid_len = flux_mid.shape[1]

        fft_len = sc_fft.next_fast_len(mid_len + kernel_size - 1,
                                       real=True)
        ft_flux_mid = sc_fft.rfft(flux_mid, n=fft_len)

        return (np.hstack((flux_lo, flux_hi)), ft_flux_mid, fft_len)


    @staticmethod
    def conv_pre_select(conv_pre, rows):
        """
        Filter a set of preprocessing data for vdisp convolution
        to use only the specified rows (templates)

        Parameters
        ----------
        conv_pre: :class:`tuple` or None
            Preprocessing structure for vdisp convolution (may be None)
        rows: :class:`int`
            Which rows (templates) from the preprocessing
            data should we keep?

        Returns
        -------
        Revised preprocessing data with only the selected rows

        """
        if conv_pre is None:
            return None
        else:
            flux_lohi, ft_flux_mid, fft_len = conv_pre
            return (flux_lohi[rows, :], ft_flux_mid[rows, :], fft_len)


    def convolve_vdisp_from_pre(self, flux_lohi, ft_flux_mid, flux_len, fft_len, vdisp):
        """
        Convolve an array of fluxes with a velocity dispersion, using
        precomputed Fourier transform information for the fluxes.  The
        raw array is not passed as an argument; instead, it is decomposed
        into a central region, for which we receive only the Fourier transform
        in ft_flux_mid, and peripheral regions, for which we receive the raw fluxes
        in flux_lohi.

        Parameters
        ----------
        flux_lohi: :class:`np.ndarray` (float64)
            1D array of fluxes for wavelengths outside the wavelength range
            defined by self.pixkms_bounds.  Lower and upper ranges are
            concatenated together.
        ft_flux_mid: :class:`np.ndarray` (complex128)
            FT of fluxes for range definded by self.pixkms_bounds
        flux_len: :class:`int`
            Total length of flux array to be returned
        fft_len: :class:`int`
            Length of the FFT for ft_flux_mid
        vdisp: :class:`float64`
            Velocity dispersion to convolve with fluxes

        Returns
        -------
        1D array of flux_len fluxes, equivalent to what would be
        computed for the raw array of input fluxes by convolve_vdisp().

        """

        import scipy.fft as sc_fft

        assert vdisp <= Templates.MAX_PRE_VDISP

        output = np.empty(flux_len)

        pixsize_kms = Templates.PIXKMS
        sigma = vdisp / pixsize_kms # [pixels]

        radius = Templates._gaussian_radius(sigma)
        kernel = Templates._gaussian_kernel1d(sigma, radius)

        # compute FFT of Gaussian kernel, then complete convolution
        ft_kernel = sc_fft.rfft(kernel, n=fft_len)
        conv      = sc_fft.irfft(ft_flux_mid * ft_kernel, n=fft_len)

        lo, hi = self.pixkms_bounds

        # extract middle of convolution (eqv to mode='same')
        s = len(kernel)//2
        e = s + hi - lo
        output[lo:hi] = conv[s:e]

        output[:lo] = flux_lohi[:lo]
        output[hi:] = flux_lohi[lo:]

        return output


    def convolve_vdisp(self, templateflux, vdisp):
        """
        Convolve one or more arrays of fluxes with a velocity dispersion.

        Parameters
        ----------
        templateflux: :class:`np.ndarray`
           Either a 1D template array of fluxes, or a 2D array of size
          [ntemplates x nfluxes] representing ntemplates fluxes
        vdisp: :class:`float64`
           Velocity dispersion to convolve with fluxes

        Returns
        -------
        Array with convlution of template(s) with vdisp.  Only
        fluxes in the range self.pixkms_bounds are convolved;
        the rest are copied unchanged.

        """
        from scipy.signal import oaconvolve

        # Convolve by the velocity dispersion.
        if vdisp <= 0.:
            output = templateflux.copy()
        else:
            output = np.empty_like(templateflux)
            pixsize_kms = Templates.PIXKMS
            sigma = vdisp / pixsize_kms # [pixels]

            radius = Templates._gaussian_radius(sigma)
            kernel = Templates._gaussian_kernel1d(sigma, radius)

            lo, hi = self.pixkms_bounds

            if templateflux.ndim == 1:
                output[lo:hi] = self.convolve(
                    templateflux[lo:hi], kernel, mode='same')
                output[:lo] = templateflux[:lo]
                output[hi:] = templateflux[hi:]
            else:
                ntemplates = templateflux.shape[0]
                for ii in range(ntemplates):
                    output[ii, lo:hi] = self.convolve(
                        templateflux[ii, lo:hi], kernel, mode='same')
                output[:, :lo] = templateflux[:, :lo]
                output[:, hi:] = templateflux[:, hi:]

        return output


    @staticmethod
    def _gaussian_radius(sigma):
        """
        Compute the radius of a Gaussian filter with sttdev sigma.

        Note: truncation removes very small tail values of Gaussian
        to limit size of filter.

        """
        truncate = 4.
        return int(truncate * sigma + 0.5)


    @staticmethod
    def _gaussian_kernel1d(sigma, radius, order=0):
        """
        Computes a 1-D Gaussian convolution kernel.

        Borrowed from scipy.ndimage.  Width of kernel
        is 2 * radius + 1.

        order == k > 0 --> compute kth derivative of kernel

        """
        sigma2 = sigma * sigma
        x = np.arange(-radius, radius+1, dtype=np.float64)
        phi_x = np.exp(-0.5 / sigma2 * x ** 2)
        phi_x /= phi_x.sum()

        if order == 0:
            return phi_x
        else:
            # f(x) = q(x) * phi(x) = q(x) * exp(p(x))
            # f'(x) = (q'(x) + q(x) * p'(x)) * phi(x)
            # p'(x) = -1 / sigma ** 2
            # Implement q'(x) + q(x) * p'(x) as a matrix operator and apply to the
            # coefficients of q(x)
            exponent_range = np.arange(order + 1)
            q = np.zeros(order + 1)
            q[0] = 1
            D = np.diag(exponent_range[1:], 1)  # D @ q(x) = q'(x)
            P = np.diag(np.ones(order)/-sigma2, -1)  # P @ q(x) = q(x) * p'(x)
            Q_deriv = D + P
            for _ in range(order):
                q = Q_deriv.dot(q)
            q = (x[:, None] ** exponent_range).dot(q)
            return q * phi_x


    @staticmethod
    def klambda(wave):
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
        dust_power = -0.7     # power-law slope
        dust_normwave = 5500. # pivot wavelength

        return (wave / dust_normwave)**dust_power
