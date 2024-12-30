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
    """
    A BoxedScalar is an item of an Numpy
    structured scalar type that is initialized
    to all zeros and can then be passed
    around by reference.  Access the .value
    field to unbox the scalar.

    """

    def __init__(self, dtype):
        self.value = np.zeros(1, dtype=dtype)[0]

    def __getitem__(self, key):
        return self.value[key]

    def __setitem__(self, key, v):
        self.value[key] = v


class MPPool(object):
    """
    A Pool encapsulates paraallel execution with a
    multiprocessing.Pool, falling back to sequential
    execution in the current process if just one worker
    is requested.

    Unlike multiprocessing.Pool, our starmap() function
    takes a list of keyword argument dictionaries,
    rather than a list of positional arguments.

    """
    def __init__(self, nworkers, initializer=None, init_argdict=None):
        """
        create a pool with nworkers workers, using the current
        process if nworkers is 1.  If initiializer is not None,
        apply this function to the arguments in keyword dictionary
        init_argdict on startup in each each worker subprocess.

        """
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
        """
        apply function func to each of a list of inputs, represented
        as a list of keyword argument dictionaries.

        """
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
        """
        close our multiprocess pool if we created one

        """
        if self.pool is not None:
            self.pool.close()

    @staticmethod
    def apply_to_dict(f, argdict):
        return f(**argdict)


class ZWarningMask(object):
    """
    Mask bit definitions for zwarning.
    Taken from Redrock/0.15.4
    WARNING on the warnings: not all of these are implemented yet.

    #- TODO: Consider using something like desispec.maskbits to provide a more
    #- convenient wrapper class (probably copy it here; don't make a dependency)
    #- That class as-is would bring in a yaml dependency.
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
    """Convert SFD E(B-V) value to dust transmission 0-1 given the bandpass.

    Args:
        ebv (float or array-like): SFD E(B-V) value(s)
        filtername (str): Filter name, e.g., 'decam2014-r'.

    Returns:
        Scalar or array (same as ebv input), Milky Way dust transmission 0-1.

    Note:

        This function tabulates the total-to-selective extinction ratio,
        k_X=A(X)/E(B-V) for many different filter bandpasses, X, where
        A(X)=-2.5*log10(transmission in X band). And so given the ebv, it
        returns mwdust_transmission=10**(-0.4*k_X*ebv).

    Returns:
        scalar, total extinction A(band) = -2.5*log10(transmission(band))

    Notes:
        Based on `desiutil.dust.mwdust_transmission`.

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
    """Simple function to safely turn a (scalar, for now) variance into an
    inverse variance.

    if sigma=True then assume that `var` is a standard deviation

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


def ivar2var(ivar, clip=1e-3, sigma=False, allmasked_ok=False):
    """Safely convert an inverse variance to a variance. Note that we clip at 1e-3
    by default, not zero.

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
    """http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion"""
    if airwave <= 0:
        raise ValueError('Input wavelength is not defined.')
    ss = 1e4 / airwave
    nn = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - ss**2) + 0.0001599740894897 / (38.92568793293 - ss**2)
    return airwave * nn


def find_minima(x):
    """Return indices of local minima of x, including edges.

    The indices are sorted small to large.

    Note:
        this is somewhat conservative in the case of repeated values:
        find_minima([1,1,1,2,2,2]) -> [0,1,2,4,5]

    Args:
        x (array-like): The data array.

    Returns:
        (array): The indices.

    """
    x = np.asarray(x)
    ii = np.where(np.r_[True, x[1:]<=x[:-1]] & np.r_[x[:-1]<=x[1:], True])[0]

    jj = np.argsort(x[ii])

    return ii[jj]


def minfit(x, y, return_coeff=False):
    """Fits y = y0 + ((x-x0)/xerr)**2

    See redrock.zwarning.ZWarningMask.BAD_MINFIT for zwarn failure flags

    Args:
        x (array): x values.
        y (array): y values.

    Returns:
        (tuple):  (x0, xerr, y0, zwarn) where zwarn=0 is good fit.

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
    return np.quantile(A, q)


# Numba's median impl is also faster
@jit(nopython=True, nogil=True, cache=True)
def median(A):
    return np.median(A)


# Open-coded Numba trapz is much faster than np.traz
@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def trapz(y, x):
    res = 0.
    for i in range(len(x) - 1):
        res += (x[i+1] - x[i]) * (y[i+1] + y[i])
    return 0.5 * res


def trapz_rebin(src_x, src_y, bin_centers, out=None, pre=None):
    """
    Resample an signal src_y sampled at points src_x into
    into a series of x-axis bins, in an area-conserving way

    Parameters
    ----------
    src_y : :class:`numpy.ndarray` [npix]
        Input signal
    src_x : :class:`numpy.ndarray` [npix]
        array of x coordinates corresponing to src_y
    bin_centers : :class:`numpy.ndarray` [noutpix]
        center x coordinates of output bins
    out:  :class:`numpy.ndarray` or None
        if not None, an array of correct size and type
        that will receive the result of the computation.
        If None, a new array will be allocated.
    pre :
        preprocessing data computed by trapz_rebin_pre(),
        if available.  If not None, it must correspond
        to the input bin_centers array.

    Returns
    -------
    :class:`numpy.ndarray` [noutpix]
        Resampled signal at centers of each bin

    """
    if pre == None:
        pre = trapz_rebin_pre(bin_centers)

    edges, ibw = pre

    return _trapz_rebin(src_x, src_y, edges, ibw, out)


def trapz_rebin_pre(bin_centers):
    """
    Perform preprocessing on the bin centers passed to trapz_rebin()
    to avoid some computations each time it is called with these
    same centers.

    Parameters
    ----------
    bin_centers : :class:`numpy.ndarray` [npix]
        array of bin centers for trapz_rebin()

    Returns
    -------
    A preprocessing data structure that can be passed as the
    'pre' argument of trapz_rebin whenever it is called with
    the given bin_centers (for *any* src_x and src_y)

    """

    edges = centers2edges(bin_centers)
    ibw = 1. / np.diff(edges)

    return (edges, ibw)


@jit(nopython=True, nogil=True, fastmath=True, cache=True)
def _trapz_rebin(x, y, edges, ibw, out):
    """
    Trapezoidal rebinning

    This implementation is optimized for compilation with Numba.  It
    agrees with the earlier fastspecfit implementation to within around
    1e-12 and is somewhat faster.  It also tolerates the first bin being
    arbitrarily greater than the first x without performance loss. ibw is
    the array of inverse bin widths.

    We require, as did the original version of the code, that the values of x
    extend strictly beyond the edges array in both directions.

    """
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


@jit(nopython=True, nogil=True, cache=True)
def centers2edges(centers):
    """
    Convert bin centers to bin edges, guessing at what you probably meant.

    Parameters
    ----------
    centers: :class:`numpy.ndarray`
      bin centers

    Returns
    -------
    array: bin edges, length = len(centers) + 1

    """

    edges = np.empty(len(centers) + 1, dtype=centers.dtype)

    #- Interior edges are just points half way between bin centers
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])

    #- edge edges are extrapolation of interior bin sizes
    edges[0]  = centers[0]  - (centers[1]  - edges[1])
    edges[-1] = centers[-1] + (centers[-1] - edges[-2])

    return edges
