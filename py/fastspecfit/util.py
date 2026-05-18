"""
fastspecfit.util
================

General utilities.

"""
import numpy as np
from numba import jit

from fastspecfit.logger import log

try:  # this fails when building the documentation
    from scipy import constants
    C_LIGHT = constants.c / 1000.0  # [km/s]
except:
    C_LIGHT = 299792.458  # [km/s]

FLUXNORM = 1e17  # flux normalization factor for all DESI spectra [erg/s/cm2/A]

TINY = np.nextafter(0, 1, dtype=np.float32)
SQTINY = np.sqrt(TINY)
F32MAX = np.finfo(np.float32).max


class BoxedScalar(object):
    """A zero-initialized NumPy structured scalar that can be passed by reference.

    Parameters
    ----------
    dtype : dtype-like
        NumPy dtype of the scalar value.

    Attributes
    ----------
    value : numpy scalar
        The boxed value; access via ``.value`` to unbox.

    """
    def __init__(self, dtype):
        self.value = np.zeros(1, dtype=dtype)[0]

    def __getitem__(self, key):
        return self.value[key]

    def __setitem__(self, key, v):
        self.value[key] = v


class MPPool(object):
    """Parallel execution pool that falls back to sequential for a single worker.

    Unlike :class:`multiprocessing.Pool`, :meth:`starmap` accepts a list of
    keyword-argument dictionaries rather than positional arguments.

    Parameters
    ----------
    nworkers : :class:`int`
        Number of worker processes. When 1, work runs in the current process.
    initializer : callable or None, optional
        Function to call in each worker subprocess at startup.
    init_argdict : :class:`dict` or None, optional
        Keyword arguments passed to ``initializer`` at startup.

    """
    def __init__(self, nworkers, initializer=None, init_argdict=None):
        initfunc = None if initializer is None else self.apply_to_dict

        # If multiprocessing, create a pool of worker processes and initialize
        # single-copy objects in each worker.
        if nworkers > 1:
            from multiprocessing import Pool
            self.pool = Pool(nworkers,
                             initializer=initfunc,
                             initargs=(initializer, init_argdict,))
        else:
            self.pool = None


    def starmap(self, func, argdicts):
        """Apply ``func`` to each keyword-argument dict in ``argdicts``."""
        # we cannot pickle a local function, so we must pass
        # both func and the argument dictionary to the subprocess
        # worker and have it apply one to the other.
        wrapped_args = [ ( func, a, ) for a in argdicts ]

        if self.pool is not None:
            out = self.pool.starmap(self.apply_to_dict, wrapped_args)
        else:
            from itertools import starmap
            out = starmap(self.apply_to_dict, wrapped_args)

        return out

    def close(self):
        """Close the multiprocessing pool if one was created."""
        if self.pool is not None:
            self.pool.close()

    @staticmethod
    def apply_to_dict(f, argdict):
        return f(**argdict)


class ZWarningMask(object):
    """Redrock ``ZWARN`` bitmask definitions (from Redrock 0.15.4).

    Not all flags are currently used by fastspecfit.

    """
    SKY               = 2**0  #- sky fiber
    LITTLE_COVERAGE   = 2**1  #- too little wavelength coverage
    SMALL_DELTA_CHI2  = 2**2  #- chi-squared of best fit is too close to that of second best
    NEGATIVE_MODEL    = 2**3  #- synthetic spectrum is negative
    MANY_OUTLIERS     = 2**4  #- fraction of points more than 5 sigma away from best model is too large (>0.05)
    Z_FITLIMIT        = 2**5  #- chi-squared minimum at edge of the redshift fitting range
    NEGATIVE_EMISSION = 2**6  #- a QSO line exhibits negative emission, triggered only in QSO spectra, if  C_IV, C_III, Mg_II, H_beta, or H_alpha has LINEAREA + 3 * LINEAREA_ERR < 0
    UNPLUGGED         = 2**7  #- the fiber was unplugged/broken, so no spectrum obtained
    BAD_TARGET        = 2**8  #- catastrophically bad targeting data
    NODATA            = 2**9  #- No data for this fiber, e.g. because spectrograph was broken during this exposure (ivar=0 for all pixels)
    BAD_MINFIT        = 2**10 #- Bad parabola fit to the chi2 minimum
    POORDATA          = 2**11 #- Poor input data quality but try fitting anyway

    @classmethod
    def flags(cls):
        flagmask = list()
        for key, value in cls.__dict__.items():
            if not key.startswith('_') and key.isupper():
                flagmask.append((key, value))

        import numpy as np
        isort = np.argsort([x[1] for x in flagmask])
        flagmask = [flagmask[i] for i in isort]
        return flagmask


def mwdust_transmission(ebv, filtername):
    """Convert SFD E(B-V) to Milky Way dust transmission in a given bandpass.

    Parameters
    ----------
    ebv : :class:`float` or array-like
        SFD E(B-V) reddening value(s).
    filtername : :class:`str`
        Filter name, e.g. ``'decam2014-r'``.

    Returns
    -------
    transmission : :class:`float` or :class:`numpy.ndarray`
        Milky Way dust transmission (0–1), same shape as ``ebv``.

    Notes
    -----
    Uses tabulated k_X = A(X)/E(B-V) values where
    transmission = 10^{-0.4 * k_X * ebv}. Based on
    ``desiutil.dust.mwdust_transmission``.

    """
    k_X = {
        # From https://desi.lbl.gov/trac/wiki/ImagingStandardBandpass
        # DECam u  3881.6   3.994
        # DECam g  4830.8   3.212
        # DECam r  6409.0   2.164
        # DECam i  7787.5   1.591
        # DECam z  9142.7   1.211
        # DECam Y  9854.5   1.063
        # BASS g  4772.1   3.258
        # BASS r  6383.6   2.176
        # MzLS z  9185.1   1.199
        # Consistent with the synthetic magnitudes and function dust_transmission
        'BASS-g': 3.258,
        'BASS-r': 2.176,
        'MzLS-z': 1.199,
        'decam2014-u': 3.994,
        'decam2014-g': 3.212,
        'decam2014-r': 2.164,
        'decam2014-i': 1.591,
        'decam2014-z': 1.211,
        'decam2014-Y': 1.063,
        # Anand Raichoor, private communication
        'cfht_megacam-u': 4.01,
        'cfht_megacam-ustar': 4.12,
        'odin-N419': 4.324,
        'odin-N501': 3.540,
        'odin-N673': 2.438,
        'hsc2017-g': 3.24,
        'hsc2017-r': 2.276,
        'hsc2017-r2': 2.276,
        'hsc2017-i': 1.633,
        'hsc2017-i2': 1.633,
        'hsc2017-z': 1.263,
        'hsc2017-y': 1.075,
        'suprime-IB427': 4.202,
        'suprime-IB464': 3.894,
        'suprime-IB484': 3.694,
        'suprime-IB505': 3.490,
        'suprime-IB527': 3.304,
        # Add WISE from
        # https://github.com/dstndstn/tractor/blob/main/tractor/sfd.py#L23-L35
        'wise2010-W1': 0.184,
        'wise2010-W2': 0.113,
        'wise2010-W3': 0.0241,
        'wise2010-W4': 0.00910,
        }

    if filtername not in k_X:
        errmsg = f'Filtername {filtername} is missing from dictionary of known bandpasses!'
        log.critical(errmsg)
        raise ValueError(errmsg)

    A_X = k_X[filtername] * ebv

    transmission = 10.**(-0.4 * A_X)

    return transmission


def var2ivar(var, sigma=False):
    """Safely convert a scalar variance (or standard deviation) to an inverse variance.

    Parameters
    ----------
    var : :class:`float`
        Variance, or standard deviation when ``sigma=True``.
    sigma : :class:`bool`, optional
        If ``True``, treat ``var`` as a standard deviation. Default is
        ``False``.

    Returns
    -------
    ivar : :class:`float`
        Inverse variance; 0 when ``var`` is too small to invert safely.

    """
    if sigma:
        if var > SQTINY:
            ivar = 1. / var**2
        else:
            ivar = 0.
    else:
        if var > TINY:
            ivar = 1. / var
        else:
            ivar = 0.
    return ivar


def ivar2var(ivar, clip=1e-8, sigma=False, allmasked_ok=False):
    """Safely convert an inverse variance array to variance (or sigma).

    Parameters
    ----------
    ivar : :class:`numpy.ndarray`
        Inverse variance array.
    clip : :class:`float`, optional
        Minimum ``ivar`` value treated as valid. Default is 1e-8.
    sigma : :class:`bool`, optional
        If ``True``, return the square root (i.e. standard deviation).
        Default is ``False``.
    allmasked_ok : :class:`bool`, optional
        If ``True``, return zeros rather than raising when all pixels are
        masked. Default is ``False``.

    Returns
    -------
    var : :class:`numpy.ndarray`
        Variance (or sigma) array; zero where ``ivar <= clip``.
    goodmask : :class:`numpy.ndarray` of bool
        ``True`` where ``ivar > clip``.

    """
    var = np.zeros_like(ivar)
    goodmask = ivar > clip # True is good
    if np.count_nonzero(goodmask) == 0:
        # Try clipping at zero.
        goodmask = ivar > 0. # True is good
        if np.count_nonzero(goodmask) == 0:
            if allmasked_ok:
                return var, goodmask
            errmsg = 'All values are masked!'
            log.critical(errmsg)
            raise ValueError(errmsg)
    var[goodmask] = 1. / ivar[goodmask]
    if sigma:
        var = np.sqrt(var) # return a sigma
    return var, goodmask


# currently unused - JDB
def air2vac(airwave):
    """Convert an air wavelength in Angstroms to vacuum wavelength."""
    if airwave <= 0:
        raise ValueError('Input wavelength is not defined.')
    ss = 1e4 / airwave
    nn = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - ss**2) + 0.0001599740894897 / (38.92568793293 - ss**2)
    return airwave * nn


def find_minima(x):
    """Return indices of local minima of ``x``, including edges.

    Indices are sorted in ascending order of ``x`` value. Conservative with
    repeated values: ``find_minima([1,1,1,2,2,2])`` returns ``[0,1,2,4,5]``.

    Parameters
    ----------
    x : array-like
        Input data array.

    Returns
    -------
    ii : :class:`numpy.ndarray` of int
        Indices of local minima, sorted by ``x`` value (smallest first).

    """
    x = np.asarray(x)
    ii = np.where(np.r_[True, x[1:]<=x[:-1]] & np.r_[x[:-1]<=x[1:], True])[0]

    jj = np.argsort(x[ii])

    return ii[jj]


def minfit(x, y, return_coeff=False):
    """Fit a parabola y = y0 + ((x - x0) / xerr)^2 to find the minimum.

    Parameters
    ----------
    x : array-like
        x values.
    y : array-like
        y values.
    return_coeff : :class:`bool`, optional
        If ``True``, also return the raw polynomial coefficients ``(a, b,
        c)``. Default is ``False``.

    Returns
    -------
    x0 : :class:`float`
        x position of the parabola minimum (-1 on failure).
    xerr : :class:`float`
        Half-width of the parabola at unit height (-1 on failure).
    y0 : :class:`float`
        Minimum y value (-1 on failure).
    zwarn : :class:`int`
        Zero on success; :attr:`ZWarningMask.BAD_MINFIT` on failure.
    coeff : :class:`tuple`, optional
        Raw ``(a, b, c)`` polynomial coefficients; only present when
        ``return_coeff=True``.

    """
    a, b, c = 0., 0., 0.
    if len(x) < 3:
        if return_coeff:
            return (-1,-1,-1,ZWarningMask.BAD_MINFIT,(a,b,c))
        else:
            return (-1,-1,-1,ZWarningMask.BAD_MINFIT)

    try:
        #- y = a x^2 + b x + c
        a,b,c = np.polyfit(x,y,2)
    except np.linalg.LinAlgError:
        if return_coeff:
            return (-1,-1,-1,ZWarningMask.BAD_MINFIT,(a,b,c))
        else:
            return (-1,-1,-1,ZWarningMask.BAD_MINFIT)

    if a == 0.0:
        if return_coeff:
            return (-1,-1,-1,ZWarningMask.BAD_MINFIT,(a,b,c))
        else:
            return (-1,-1,-1,ZWarningMask.BAD_MINFIT)

    #- recast as y = y0 + ((x-x0)/xerr)^2
    x0 = -b / (2*a)
    y0 = -(b**2) / (4*a) + c

    zwarn = 0
    if (x0 <= np.min(x)) or (np.max(x) <= x0):
        zwarn |= ZWarningMask.BAD_MINFIT
    if (y0<=0.):
        zwarn |= ZWarningMask.BAD_MINFIT

    if a > 0.0:
        xerr = 1 / np.sqrt(a)
    else:
        xerr = 1 / np.sqrt(-a)
        zwarn |= ZWarningMask.BAD_MINFIT

    if return_coeff:
        return (x0, xerr, y0, zwarn,(a,b,c))
    else:
        return (x0, xerr, y0, zwarn)


#
# sigma clipping in Numba
# Rewrite basic sigma clipping to avoid
# array copies and redundant summation
# on each iteration
#
@jit(nopython=True, nogil=True, cache=True)
def sigmaclip(c, low=3., high=3.):
    """Iterative sigma-clipping; returns the clipped array and a boolean mask."""
    n  = len(c)
    mask = np.full(n, True, dtype=np.bool_)

    s  = np.sum(c)
    s2 = np.sum(c*c)

    delta = 1
    while delta > 0:
        mean = s/n

        # differences very close to zero can end up negative
        # due to limited floating-point precision
        std = np.sqrt(np.maximum(s2/n - mean*mean, 0.))
        clo = mean - std * low
        chi = mean + std * high

        if std == 0.: # don't mask everything
            break

        n0 = n
        for j, cval in enumerate(c):
            if mask[j] and (cval < clo or cval > chi):
                mask[j] = False
                n  -= 1
                s  -= cval
                s2 -= cval * cval

        delta = n0-n

    return c[mask], mask


# Numba's quantile impl is much faster
# than Numpy's standard version
@jit(nopython=True, nogil=True, cache=True)
def quantile(A, q):
    """Compute quantile(s) ``q`` of array ``A`` (Numba-accelerated)."""
    return np.quantile(A, q)


# Numba's median impl is also faster
@jit(nopython=True, nogil=True, cache=True)
def median(A):
    """Compute the median of array ``A`` (Numba-accelerated)."""
    return np.median(A)


# Open-coded Numba trapz is much faster than np.traz
@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def trapz(y, x):
    """Trapezoidal integration of ``y`` over ``x`` (Numba-accelerated)."""
    res = 0.
    for i in range(len(x) - 1):
        res += (x[i+1] - x[i]) * (y[i+1] + y[i])
    return 0.5 * res


def trapz_rebin(src_x, src_y, bin_centers, out=None, pre=None):
    """Resample ``src_y`` into bins centered at ``bin_centers`` (area-conserving).

    Parameters
    ----------
    src_x : :class:`numpy.ndarray`
        x coordinates of the input signal.
    src_y : :class:`numpy.ndarray`
        Input signal sampled at ``src_x``.
    bin_centers : :class:`numpy.ndarray`
        Center x coordinates of the output bins.
    out : :class:`numpy.ndarray` or None, optional
        Pre-allocated output array; allocated if ``None``.
    pre : :class:`tuple` or None, optional
        Preprocessing data from :func:`trapz_rebin_pre`; computed from
        ``bin_centers`` when ``None``.

    Returns
    -------
    out : :class:`numpy.ndarray`
        Resampled signal at the center of each output bin.

    """
    if pre == None:
        pre = trapz_rebin_pre(bin_centers)

    edges, ibw = pre

    return _trapz_rebin(src_x, src_y, edges, ibw, out)


def trapz_rebin_pre(bin_centers):
    """Precompute bin edges and inverse widths for :func:`trapz_rebin`.

    Parameters
    ----------
    bin_centers : :class:`numpy.ndarray`
        Center x coordinates of the output bins.

    Returns
    -------
    pre : :class:`tuple`
        ``(edges, ibw)`` where ``edges`` are the bin edges and ``ibw`` is
        the array of inverse bin widths. Pass as the ``pre`` argument to
        :func:`trapz_rebin` to avoid recomputation.

    """

    edges = centers2edges(bin_centers)
    ibw = 1. / np.diff(edges)

    return (edges, ibw)


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _trapz_rebin(x, y, edges, ibw, out):
    """Numba-accelerated trapezoidal rebinning core."""
    # interpolated value of y at edge_x, which lies between x[j] and x[j+1]
    def y_at(edge_x, j): # j: largest j s.t. x[j] < edge_x
        return y[j] + (edge_x - x[j]) * (y[j+1] - y[j]) / (x[j+1] - x[j])

    if edges[0] < x[0] or x[-1] < edges[-1]:
        raise ValueError('edges must lie strictly within x range')

    nbins = len(edges) - 1

    results = np.empty(nbins, dtype=y.dtype) if out is None else out

    # greatest j s.t. x[j] < edges[0]
    j = np.searchsorted(x, edges[0], 'right')

    xlo = edges[0]
    ylo = y_at(xlo, j)

    # loop invariant: on entry to iteration i,
    #   x[j] is greatest x < edges[i]
    #   xlo = edges[i]
    #   ylo = y_at(edges[i], j)
    for i in range(nbins):

        a = 0.

        while x[j+1] < edges[i+1]:
            xhi = x[j+1]
            yhi = y[j+1]

            # add area from prev boundary to x[j+1]
            a += (xhi - xlo) * (yhi + ylo)

            xlo = xhi
            ylo = yhi

            j += 1

        # partial area up to edge i+1
        xhi = edges[i+1]
        yhi = y_at(xhi, j)

        a += (xhi - xlo) * (yhi + ylo)

        xlo = xhi
        ylo = yhi

        results[i] = 0.5 * ibw[i] * a

    return results


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _trapz_rebin_batch(x, Y, edges, ibw, out):
    """Apply :func:`_trapz_rebin` to each row of ``Y`` (ntemplates × npix).

    Parameters
    ----------
    x : :class:`numpy.ndarray`, shape (npix,)
    Y : :class:`numpy.ndarray`, shape (ntemplates, npix)
    edges, ibw : precomputed bin edges and inverse widths from :func:`trapz_rebin_pre`
    out : :class:`numpy.ndarray`, shape (ntemplates, nbins)

    """
    for t in range(Y.shape[0]):
        _trapz_rebin(x, Y[t], edges, ibw, out[t])


@jit(nopython=True, nogil=True, cache=True)
def centers2edges(centers):
    """Convert bin centers to bin edges by linear extrapolation at the boundaries.

    Parameters
    ----------
    centers : :class:`numpy.ndarray`
        Bin center coordinates.

    Returns
    -------
    edges : :class:`numpy.ndarray`
        Bin edge coordinates, length ``len(centers) + 1``.

    """

    edges = np.empty(len(centers) + 1, dtype=centers.dtype)

    #- Interior edges are just points half way between bin centers
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])

    #- edge edges are extrapolation of interior bin sizes
    edges[0]  = centers[0]  - (centers[1]  - edges[1])
    edges[-1] = centers[-1] + (centers[-1] - edges[-2])

    return edges


def radec2pix(nside, ra, dec):
    """Convert RA/Dec to HEALPix nested pixel numbers.

    Parameters
    ----------
    nside : :class:`int`
        HEALPix ``nside`` parameter, must be ``2**k`` for 0 < k < 30.
    ra : :class:`float` or array-like
        Right ascension in degrees.
    dec : :class:`float` or array-like
        Declination in degrees.

    Returns
    -------
    pix : :class:`numpy.ndarray` of int
        HEALPix pixel numbers in the nested scheme.

    """
    import healpy as hp
    theta, phi = np.radians(90-dec), np.radians(ra)
    if np.isnan(np.sum(theta)) :
        raise ValueError("some NaN theta values")

    if np.sum((theta < 0)|(theta > np.pi))>0 :
        raise ValueError("some theta values are outside [0,pi]: {}".format(theta[(theta < 0)|(theta > np.pi)]))

    return hp.ang2pix(nside, theta, phi, nest=True)
