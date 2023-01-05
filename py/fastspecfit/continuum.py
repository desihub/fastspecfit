"""
fastspecfit.continuum
=====================

Methods and tools for continuum-fitting.

"""
import pdb # for debugging

import os, time
import numpy as np

import astropy.units as u
from astropy.table import Table, Column

from fastspecfit.util import C_LIGHT

#import numba
#@numba.jit(nopython=True)
def nnls(A, b, maxiter=None, eps=1e-7, v=False):
    # https://en.wikipedia.org/wiki/Non-negative_least_squares#Algorithms
    # not sure what eps should be set to
    m, n = A.shape
    if b.shape[0] != m:
        raise ValueError()#f"Shape mismatch: {A.shape} {b.shape}. Expected: (m, n) (m) ")
    if maxiter is None:
        # this is what scipy does when maxiter is not specified
        maxiter = 3 * n
    # init
    P = np.zeros(n).astype(bool)
    R = np.ones(n).astype(bool)
    x = np.zeros(n)
    s = np.zeros(n)
    # compute these once ahead of time
    ATA = A.T.dot(A)
    ATb = A.T.dot(b)
    # main loop
    w = ATb - ATA.dot(x)
    j = np.argmax(R * w)
    c = 0 # iteration count
    # while R != {} and max(w[R]) > eps
    while np.any(R) and w[j] > eps:
        #if v: print(f"{c=}", f"{P=}\n {R=}\n {w=}\n {j=}")
        # add j to P, remove j from R
        P[j], R[j] = True, False
        s[P] = np.linalg.inv(ATA[P][:, P]).dot(ATb[P])
        s[R] = 0
        d = 0 # inner loop iteration count, for debugging
        #if v: print(f"{c=}", f"{P=}\n {R=}\n {s=}")
        # while P != {} and min(s[P]) < eps
        # make sure P is not empty before checking min s[P]
        while np.any(P) and np.min(s[P]) < eps:
            i = P & (s < eps)
            #if v: print(f" {d=}", f"  {P=}\n  {i=}\n  {s=}\n  {x=}")
            a = np.min(x[i] / (x[i] - s[i]))
            x = x + a * (s - x)
            j = P & (x < eps)
            R[j], P[j] = True, False
            s[P] = np.linalg.inv(ATA[P][:, P]).dot(ATb[P])
            s[R] = 0
            d += 1
        x[:] = s
        # w = A.T.dot(b - A.dot(x))
        w = ATb - ATA.dot(x)
        j = np.argmax(R * w)
        #if v: print(f"{c=}", f"{P=}\n {R=}\n {w=}\n {j=}")
        c += 1
        if c >= maxiter:
            break
    res = np.linalg.norm(A.dot(x) - b)
    return x, res

class ContinuumTools(object):
    """Tools for dealing with stellar continua.

    Parameters
    ----------
    ssptemplates : :class:`str`, optional
        Full path to the SSP templates used for continuum-fitting.
    sspversion : :class:`str`, optional, defaults to `v1.0`
        Version of the SSP templates.
    mapdir : :class:`str`, optional
        Full path to the Milky Way dust maps.
    metallicity : :class:`str`, optional, defaults to `Z0.0190`.
        Stellar metallicity of the SSPs. Currently fixed at solar
        metallicity, Z=0.0190.
    minsspwave : :class:`float`, optional, defaults to None
        Minimum SSP wavelength to read into memory. If ``None``, the minimum
        available wavelength is used (around 100 Angstrom).
    maxsspwave : :class:`float`, optional, defaults to 6e4
        Maximum SSP wavelength to read into memory. 

    .. note::
        Need to document all the attributes.

    """
    def __init__(self, ssptemplates=None, sspversion='v1.0', minsspwave=None,
                 maxsspwave=40e4, mapdir=None, verbose=False):

        import fitsio
        from astropy.cosmology import FlatLambdaCDM

        from speclite import filters
        from desiutil.dust import SFDMap
        from desiutil.log import get_logger, DEBUG

        from fastspecfit.emlines import read_emlines
        from fastspecfit.io import FASTSPECFIT_TEMPLATES_NERSC, DUST_DIR_NERSC, DR9_DIR_NERSC

        if verbose:
            self.log = get_logger(DEBUG)
        else:
            self.log = get_logger()

        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        # pre-compute the luminosity distance on a grid
        #self.redshift_ref = np.arange(0.0, 5.0, 0.05)
        #self.dlum_ref = self.cosmo.luminosity_distance(self.redshift_ref).to(u.pc).value

        self.fluxnorm = 1e17 # normalization factor for the spectra
        self.massnorm = 1e10 # stellar mass normalization factor for the SSPs [Msun]

        # dust maps
        if mapdir is None:
            mapdir = os.path.join(os.environ.get('DUST_DIR', DUST_DIR_NERSC), 'maps')
        self.SFDMap = SFDMap(scaling=1.0, mapdir=mapdir)
        #self.SFDMap = SFDMap(scaling=0.86, mapdir=mapdir) # SF11 recalibration of the SFD maps
        self.RV = 3.1

        # SSPs
        if ssptemplates is not None:
            self.ssptemplates = ssptemplates
        else:
            templates_dir = os.environ.get('FASTSPECFIT_TEMPLATES', FASTSPECFIT_TEMPLATES_NERSC)
            self.ssptemplates = os.path.join(templates_dir, 'fastspecfit-templates-{}.fits'.format(sspversion))
        if not os.path.isfile(self.ssptemplates):
            errmsg = 'SSP templates file not found {}'.format(self.ssptemplates)
            self.log.critical(errmsg)
            raise IOError(errmsg)

        self.log.info('Reading {}'.format(self.ssptemplates))
        wave, wavehdr = fitsio.read(self.ssptemplates, ext='WAVE', header=True) # [npix]
        sspflux = fitsio.read(self.ssptemplates, ext='FLUX')  # [npix,nsed]
        sspflux_vdisp, vdisphdr = fitsio.read(self.ssptemplates, ext='FLUXVDISP', header=True) # [npix,nsed,nvdisp]
        sspinfo = Table(fitsio.read(self.ssptemplates, ext='METADATA'))

        # Trim the wavelengths and select the number/ages of the templates.
        # https://www.sdss.org/dr14/spectro/galaxy_mpajhu
        if minsspwave is None:
            minsspwave = np.min(wave)
        wavekeep = np.where((wave >= minsspwave) * (wave <= maxsspwave))[0]

        self.sspwave = wave[wavekeep]
        self.sspflux = sspflux[wavekeep, :]
        self.sspflux_vdisp = sspflux_vdisp[wavekeep, :, :]
        self.sspinfo = sspinfo
        self.nsed = len(sspinfo)
        self.npix = len(wavekeep)

        self.continuum_pixkms = wavehdr['PIXSZBLU'] # pixel size [km/s]
        self.pixkms_wavesplit = wavehdr['PIXSZSPT'] # wavelength where the pixel size changes [A]

        # see bin/build-fsps-templates
        vdisp_nominal = 150.0
        nvdisp = int(np.ceil((vdisphdr['VDISPMAX'] - vdisphdr['VDISPMIN']) / vdisphdr['VDISPRES'])) + 1
        vdisp = np.linspace(vdisphdr['VDISPMIN'], vdisphdr['VDISPMAX'], nvdisp)

        if not vdisp_nominal in vdisp:
            errmsg = 'Nominal velocity dispersion is not in velocity dispersion vector.'
            self.log.critical(errmsg)
            raise ValueError(errmsg)
        
        self.vdisp = vdisp
        self.vdisp_nominal = vdisp_nominal
        self.vdisp_nominal_indx = np.where(vdisp == vdisp_nominal)[0]
        self.nvdisp = nvdisp

        # emission lines
        self.linetable = read_emlines()

        self.linemask_sigma_narrow = 200.0  # [km/s]
        self.linemask_sigma_balmer = 1000.0 # [km/s]
        self.linemask_sigma_broad = 2000.0  # [km/s]

        # photometry
        self.bands = np.array(['g', 'r', 'z', 'W1', 'W2', 'W3', 'W4'])
        self.synth_bands = np.array(['g', 'r', 'z']) # for synthesized photometry
        self.fiber_bands = np.array(['g', 'r', 'z']) # for fiber fluxes

        self.decam = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z')
        self.bassmzls = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z')

        self.decamwise = filters.load_filters(
            'decam2014-g', 'decam2014-r', 'decam2014-z', 'wise2010-W1', 'wise2010-W2', 'wise2010-W3', 'wise2010-W4')
        self.bassmzlswise = filters.load_filters(
            'BASS-g', 'BASS-r', 'MzLS-z', 'wise2010-W1', 'wise2010-W2', 'wise2010-W3', 'wise2010-W4')

        self.bands_to_fit = np.ones(len(self.bands), bool)
        #for B in ['W2', 'W3', 'W4']:
        #for B in ['W4']:
        #    self.bands_to_fit[self.bands == B] = False # drop W2-W4

        # rest-frame filters
        self.absmag_bands = ['U', 'B', 'V', 'sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z', 'W1']

        self.absmag_bands_00 = ['U', 'B', 'V', 'W1'] # band_shift=0.0
        self.absmag_bands_01 = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z'] # band_shift=0.1

        self.absmag_filters_00 = filters.FilterSequence((
            filters.load_filter('bessell-U'), filters.load_filter('bessell-B'),
            filters.load_filter('bessell-V'), filters.load_filter('wise2010-W1')
            ))
        
        self.absmag_filters_01 = filters.FilterSequence((
            filters.load_filter('sdss2010-u'),
            filters.load_filter('sdss2010-g'),
            filters.load_filter('sdss2010-r'),
            filters.load_filter('sdss2010-i'),
            filters.load_filter('sdss2010-z'),
            ))

        self.min_uncertainty = np.array([0.02, 0.02, 0.02, 0.05, 0.05, 0.05, 0.05]) # mag

    def get_dn4000(self, wave, flam, flam_ivar=None, redshift=None, rest=True):
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

        dn4000, dn4000_ivar = 0.0, 0.0

        if rest:
            restwave = wave
            flam2fnu =  restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
        else:
            restwave = wave / (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]

        # Require a 2-Angstrom pad around the break definition.
        wpad = 2.0
        if np.min(restwave) > (3850-wpad) or np.max(restwave) < (4100+wpad):
            self.log.warning('Too little wavelength coverage to compute Dn(4000)')
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
                self.log.warning('More than 20% of pixels in Dn(4000) definition are masked.')
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
            self.log.warning('Integration failed when computing DN(4000).')
            return dn4000, dn4000_ivar

        if denom == 0.0 or numer == 0.0:
            self.log.warning('DN(4000) is ill-defined or could not be computed.')
            return dn4000, dn4000_ivar
        
        dn4000 =  (blufactor / redfactor) * numer / denom
        if flam_ivar is not None:
            dn4000_ivar = (1.0 / (dn4000**2)) / (denom_var / (denom**2) + numer_var / (numer**2))

        return dn4000, dn4000_ivar

    def parse_photometry(self, bands, maggies, lambda_eff, ivarmaggies=None,
                         nanomaggies=True, nsigma=2.0, min_uncertainty=None,
                         debug=False):
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

        # Gaia-only targets can sometimes have grz=-99.
        if np.any(ivarmaggies < 0) or np.any(maggies == -99.0):
            errmsg = 'All ivarmaggies must be zero or positive!'
            self.log.critical(errmsg)
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

        #print('Hack!!')
        #if debug:
        #    maggies[3:5, 4:7] = 0.0

        dims = maggies.shape
        flatmaggies = maggies.flatten()
        flativarmaggies = ivarmaggies.flatten()

        abmag = np.zeros_like(flatmaggies)
        abmag_limit = np.zeros_like(flatmaggies)
        abmag_brighterr = np.zeros_like(flatmaggies)
        abmag_fainterr = np.zeros_like(flatmaggies)
        abmag_ivar = np.zeros_like(flatmaggies)
        
        # deal with measurements
        good = np.where(flatmaggies > 0)[0]
        if len(good) > 0:
            abmag[good] = -2.5 * np.log10(nanofactor * flatmaggies[good])

        # deal with upper limits
        snr = flatmaggies * np.sqrt(flativarmaggies)
        upper = np.where((flativarmaggies > 0) * (snr <= nsigma))[0]
        if len(upper) > 0:
            abmag_limit[upper] = - 2.5 * np.log10(nanofactor * nsigma / np.sqrt(flativarmaggies[upper]))

        # significant detections
        good = np.where(snr > nsigma)[0]
        if len(good) > 0:
            errmaggies = 1 / np.sqrt(flativarmaggies[good])
            abmag_brighterr[good] = errmaggies / (0.4 * np.log(10) * (flatmaggies[good]+errmaggies))#.astype('f4') # bright end (flux upper limit)
            abmag_fainterr[good] = errmaggies / (0.4 * np.log(10) * (flatmaggies[good]-errmaggies))#.astype('f4') # faint end (flux lower limit)
            abmag_ivar[good] = (flativarmaggies[good] * (flatmaggies[good] * 0.4 * np.log(10))**2)#.astype('f4')

        phot['abmag'] = abmag.reshape(dims)
        phot['abmag_limit'] = abmag_limit.reshape(dims)
        phot['abmag_brighterr'] = abmag_brighterr.reshape(dims)
        phot['abmag_fainterr'] = abmag_fainterr.reshape(dims)
        phot['abmag_ivar'] = abmag_ivar.reshape(dims)
            
        # Add a minimum uncertainty in quadrature **but only for flam**, which
        # is used in the fitting.
        if min_uncertainty is not None:
            self.log.debug('Propagating minimum photometric uncertainties (mag): [{}]'.format(
                ' '.join(min_uncertainty.astype(str))))
            good = np.where((maggies != 0) * (ivarmaggies > 0))[0]
            if len(good) > 0:
                factor = 2.5 / np.log(10.)
                magerr = factor / (np.sqrt(ivarmaggies[good]) * maggies[good])
                magerr2 = magerr**2 + min_uncertainty[good]**2
                ivarmaggies[good] = factor**2 / (maggies[good]**2 * magerr2)

        factor = nanofactor * 10**(-0.4 * 48.6) * C_LIGHT * 1e13 / lambda_eff**2 # [maggies-->erg/s/cm2/A]
        if ngal > 1:
            factor = factor[:, None] # broadcast for the models
        phot['flam'] = (maggies * factor)
        phot['flam_ivar'] = (ivarmaggies / factor**2)

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
        sigma = vdisp / self.continuum_pixkms # [pixels]

        smoothflux = gaussian_filter1d(sspflux, sigma=sigma, axis=0)

        return smoothflux
    
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

            # Special: mask Ly-a (1215 A)
            zlinewave = 1215.0 * (1 + redshift)
            if (zlinewave > np.min(wave)) * (zlinewave < np.max(wave)):
                sigma = maskkms_uv * zlinewave / C_LIGHT # [km/s --> Angstrom]
                I = (wave >= (zlinewave - nsig*sigma)) * (wave <= (zlinewave + nsig*sigma))
                if len(I) > 0:
                    linemask[I] = True

        if len(linemask) != npix:
            errmsg = 'Linemask must have the same number of pixels as the input spectrum.'
            self.log.critical(errmsg)
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

        # corner case for very wacky spectra
        if len(smooth_flux) == 0:
            smooth_flux = flux
            smooth_sigma = flux * 0 + np.std(flux)
        else:
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

            for xx in ax:
                #xx.set_xlim(3800, 4300)
                #xx.set_xlim(5200, 6050)
                xx.set_xlim(7000, 9000)
                #xx.set_xlim(6300, 6900)
            for xx in ax:
                xx.set_ylim(-1, 8)
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
                        
            self.log.info('{} masking sigma={:.3f} km/s and S/N={:.3f}'.format(label, linesigma, linesigma_snr))

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
        if refit and linesigma_narrow_snr > 0:
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
        if refit and linesigma_balmer_snr > 0:
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
        if refit and linesigma_uv_snr > 0:
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

    def build_linemask(self, wave, flux, ivar, redshift=0.0, nsig=7.0):
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
                    elif linesigma < 100.0:
                        linesigma_strong = 100.0
                    else:
                        linesigma_strong = linesigma
                    sigma_strong = linesigma_strong * zlinewave / C_LIGHT # [km/s --> Angstrom]
                    J = (ivar > 0) * (smoothsigma > 0) * (wave >= (zlinewave - sigma_strong)) * (wave <= (zlinewave + sigma_strong))
                    if np.sum(J) > 0:
                        snr = (flux[J] - smooth[J]) / smoothsigma[J]
                        # require peak S/N>3 and at least 5 pixels with S/N>3
                        #print(_linename, zlinewave, np.percentile(snr, 98), np.sum(snr > snr_strong))
                        if len(snr) > 5:
                            if np.percentile(snr, 98) > snr_strong and np.sum(snr > snr_strong) > 5:
                                linemask_strong[I] = True
                        else:
                            # Very narrow, strong lines can have fewer than 5
                            # pixels but if they're all S/N>3 then flag this
                            # line here.
                            if np.all(snr > snr_strong):
                                linemask_strong[I] = True
                        # Always identify Lya as "strong"
                        if _linename == 'Lya':
                            linemask_strong[I] = True


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
                        xx.plot(plotwave, plotflux, label=_linename, color='gray')
                        #xx.plot(np.hstack((wave[Jblu], wave[I], wave[Jred])), 
                        #        np.hstack((flux[Jblu], flux[I], flux[Jred])), label=_linename)
                        #xx.plot(wave[I], flux[I], label=_linename)
                        xx.scatter(wave[I], flux[I], s=10, color='orange', marker='s')
                        if np.sum(linemask_strong[I]) > 0:
                            xx.scatter(wave[I][linemask_strong[I]], flux[I][linemask_strong[I]], s=15, color='k', marker='x')
                            
                        xx.scatter(wave[Jblu], flux[Jblu], color='blue', s=10)
                        xx.scatter(wave[Jred], flux[Jred], color='red', s=10)
                        xx.set_ylim(np.min(plotflux), np.max(flux[I]))
                        xx.legend(frameon=False, fontsize=10, loc='upper left')
        
            linemask_dict = {'linemask_all': linemask, 'linemask': linemask_strong, # note we make linemask_strong the default
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

        Returns
        -------
        :class:`numpy.ndarray` [noutpix]
            Smoothed and resampled flux at the new resolution and wavelength sampling.

        Notes
        -----
        This function stands by itself rather than being in a class because we call
        it with multiprocessing, below.

        """
        from fastspecfit.util import trapz_rebin

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
    
    def SSP2data(self, _sspflux, _sspwave, redshift=0.0, vdisp=None,
                 cameras=['b', 'r', 'z'], specwave=None, specres=None, coeff=None,
                 south=True, synthphot=True, debug=False):
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
        - apply velocity dispersion broadening
        - apply the resolution matrix
        - synthesize photometry

        It also naturally handles SSPs which have been precomputed on a grid of
        velocity dispersion (and therefore have an additional
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
            self.log.critical(errmsg)
            raise ValueError(errmsg)
        
        #t0 = time.time()
        ##sspflux = sspflux.copy().reshape(npix, nmodel)
        #self.log.info('Copying the data took: {:.2f} seconds.'.format(time.time()-t0))

        # broaden for velocity dispersion but only out to ~1 micron
        if vdisp is not None:
            #sspflux = self.convolve_vdisp(sspflux, vdisp)
            I = np.where(sspwave < self.pixkms_wavesplit)[0]
            sspflux[I, :] = self.convolve_vdisp(sspflux[I, :], vdisp)

        # Apply the redshift factor. The models are normalized to 10 pc, so
        # apply the luminosity distance factor here. Also normalize to a nominal
        # stellar mass.
        if redshift:
            zsspwave = sspwave * (1.0 + redshift)
            dfactor = (10.0 / self.cosmo.luminosity_distance(redshift).to(u.pc).value)**2
            #dfactor = (10.0 / np.interp(redshift, self.redshift_ref, self.dlum_ref))**2
            factor = (self.fluxnorm * self.massnorm * dfactor / (1.0 + redshift))[np.newaxis, np.newaxis]
            zsspflux = sspflux * factor
        else:
            zsspwave = sspwave.copy() # ???
            zsspflux = self.fluxnorm * self.massnorm * sspflux

        # Optionally synthesize photometry.
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
                sspphot = self.parse_photometry(self.bands, maggies, effwave, nanomaggies=False, debug=debug)
                #self.log.info('Synthesizing photometry took: {:.2f} seconds.'.format(time.time()-t0))

        # Are we returning per-camera spectra or a single model? Handle that here.
        #t0 = time.time()
        if specwave is None and specres is None:
            datasspflux = []
            for imodel in np.arange(nmodel):
                datasspflux.append(self.smooth_and_resample(zsspflux[:, imodel], zsspwave))
            datasspflux = np.vstack(datasspflux).T

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
                if coeff is not None:
                    _datasspflux = _datasspflux.dot(coeff)
                datasspflux.append(_datasspflux)
                
        #self.log.info('Resampling took: {:.2f} seconds.'.format(time.time()-t0))

        return datasspflux, sspphot # vector or 3-element list of [npix,nmodel] spectra

