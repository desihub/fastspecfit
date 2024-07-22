"""
fastspecfit.util
================

General utilities.

"""
import numpy as np
import numba

from desiutil.log import get_logger
log = get_logger()

try: # this fails when building the documentation
    from scipy import constants
    C_LIGHT = constants.c / 1000.0 # [km/s]
except:
    C_LIGHT = 299792.458 # [km/s]

FLUXNORM = 1e17 # flux normalization factor for all DESI spectra [erg/s/cm2/A]

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

    if filtername not in k_X.keys():
        errmsg = f'Filtername {filtername} is missing from dictionary of known bandpasses!'
        log.critical(errmsg)
        raise ValueError(errmsg)

    A_X = k_X[filtername] * ebv

    transmission = 10**(-0.4 * A_X)
    
    return transmission

def ivar2var(ivar, clip=1e-3, sigma=False, allmasked_ok=False):
    """Safely convert an inverse variance to a variance. Note that we clip at 1e-3
    by default, not zero.
    
    """
    var = np.zeros_like(ivar)
    goodmask = ivar > clip # True is good
    if np.count_nonzero(goodmask) == 0:
        # Try clipping at zero.
        goodmask = ivar > 0 # True is good
        if np.count_nonzero(goodmask) == 0:
            if allmasked_ok:
                return var, goodmask
            errmsg = 'All values are masked!'
            log.critical(errmsg)
            raise ValueError(errmsg)
    var[goodmask] = 1 / ivar[goodmask]
    if sigma:
        var = np.sqrt(var) # return a sigma
    return var, goodmask

def air2vac(airwave):
    """http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion"""
    if airwave <= 0:
        raise ValueError('Input wavelength is not defined.')
    ss = 1e4 / airwave
    nn = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - ss**2) + 0.0001599740894897 / (38.92568793293 - ss**2)
    return airwave * nn


@numba.jit(nopython=True, nogil=True)
def centers2edges(centers):
    """Convert bin centers to bin edges, guessing at what you probably meant

    Args:
        centers (array): bin centers,

    Returns:
        array: bin edges, lenth = len(centers) + 1

    """
    
    edges = np.empty(len(centers) + 1, dtype=np.float64)
    
    #- Interior edges are just points half way between bin centers
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    
    #- edge edges are extrapolation of interior bin sizes
    edges[0]  = centers[0]  - (centers[1]  - edges[1])
    edges[-1] = centers[-1] + (centers[-1] - edges[-2])

    return edges


# This code is purposely written in a very "C-like" way.  The logic
# being that it may help numba optimization and also makes it easier
# if it ever needs to be ported to Cython.  Actually Cython versions
# of this code have already been tested and shown to perform no better
# than numba on Intel haswell and KNL architectures.

@numba.jit(nopython=True, fastmath=True, nogil=True)
def _trapz_rebin(x, y, edges):
    '''
    Numba-friendly version of trapezoidal rebinning

    See redrock.rebin.trapz_rebin() for input descriptions.
    `results` is pre-allocated array of length len(edges)-1 to keep results
    '''
    nbins = len(edges) - 1
        
    results = np.empty(nbins, dtype=np.float64)
    
    j = 0  #- index counter for inputs
    for i in range(nbins):
        
        results[i] = 0.
        
        # seek least j s.t. edges[i] < x[j]
        while x[j] <= edges[i]:
            j += 1
            
        #- What is the y value where the interpolation crossed the edge?
        yedge = y[j-1] + (edges[i] - x[j-1]) * (y[j] - y[j-1]) / (x[j] - x[j-1])
        
        #- Is this sample inside this bin?
        if x[j] < edges[i+1]:
            results[i] += (y[j] + yedge) * (x[j] - edges[i]) # 2*area
            
            #- Continue with interior bins
            while x[j+1] < edges[i+1]:
                j += 1
                results[i] += (y[j] + y[j-1]) * (x[j] - x[j-1]) # 2*area

            #- Next sample will be outside this bin; handle upper edge
            yedge = y[j] + (edges[i+1] - x[j]) * (y[j+1] - y[j]) / (x[j+1] - x[j])
            results[i] += (yedge + y[j]) * (edges[i+1] - x[j]) # 2*area
            
        #- Otherwise the samples span over this bin
        else:
            ylo = y[j] + (edges[i]   - x[j]) * (y[j] - y[j-1]) / (x[j] - x[j-1])
            yhi = y[j] + (edges[i+1] - x[j]) * (y[j] - y[j-1]) / (x[j] - x[j-1])
            results[i] += (ylo + yhi) * (edges[i+1] - edges[i]) # 2*area
    
        results[i] /= 2 * (edges[i+1] - edges[i])
        
    return results


def trapz_rebin(x, y, xnew, edges=None):
    """Rebin y(x) flux density using trapezoidal integration between bin edges

    Notes:
        y is interpreted as a density, as is the output, e.g.

        >>> x = np.arange(10)
        >>> y = np.ones(10)
        >>> trapz_rebin(x, y, edges=[0,2,4,6,8])  #- density still 1, not 2
        array([ 1.,  1.,  1.,  1.])

    Args:
        x (array): input x values.
        y (array): input y values.
        edges (array): (optional) new bin edges.

    Returns:
        array: integrated results with len(results) = len(edges)-1

    Raises:
        ValueError: if edges are outside the range of x or if len(x) != len(y)

    """
    if edges is None:
        edges = centers2edges(xnew)
#    else:
#        edges = np.asarray(edges)
    
    if edges[0] < x[0] or x[-1] < edges[-1]:
        raise ValueError('edges must be within input x range')

    return _trapz_rebin(x, y, edges)


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

    #- The following bits are reserved for experiment-specific post-redrock
    #- afterburner updates to ZWARN; redrock commits to *not* setting these bits
    RESERVED16        = 2**16
    RESERVED17        = 2**17
    RESERVED18        = 2**18
    RESERVED19        = 2**19
    RESERVED20        = 2**20
    RESERVED21        = 2**21
    RESERVED22        = 2**22
    RESERVED23        = 2**23

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

class TabulatedDESI(object):
    """
    Class to load tabulated z->E(z) and z->comoving_radial_distance(z) relations within DESI fiducial cosmology
    (in LSS/data/desi_fiducial_cosmology.dat) and perform the (linear) interpolations at any z.

    >>> cosmo = TabulatedDESI()
    >>> distance = cosmo.comoving_radial_distance([0.1, 0.2])
    >>> efunc = cosmo.efunc(0.3)

    The cosmology is defined in https://github.com/abacusorg/AbacusSummit/blob/master/Cosmologies/abacus_cosm000/CLASS.ini
    and the tabulated file was obtained using https://github.com/adematti/cosmoprimo/blob/main/cosmoprimo/fiducial.py.

    Note
    ----
    Redshift interpolation range is [0, 100].

    """
    def __init__(self):
        from importlib import resources
        cosmofile = resources.files('fastspecfit').joinpath('data/desi_fiducial_cosmology.dat')

        self._z, self._efunc, self._comoving_radial_distance = np.loadtxt(cosmofile, comments='#', usecols=None, unpack=True)

        self.H0 = 100.0
        self.h = self.H0 / 100
        self.hubble_time = 3.08567758e19 / 3.15576e16 / self.H0 # Hubble time [Gyr]

    def efunc(self, z):
        r"""Return :math:`E(z)`, where the Hubble parameter is defined as :math:`H(z) = H_{0} E(z)`, unitless."""
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any():
            raise ValueError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._efunc, left=None, right=None)

    def comoving_radial_distance(self, z):
        r"""Return comoving radial distance, in :math:`\mathrm{Mpc}/h`."""
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any():
            raise ValueError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._comoving_radial_distance, left=None, right=None)

    def luminosity_distance(self, z):
        r"""Return luminosity distance, in :math:`\mathrm{Mpc}/h`."""
        return self.comoving_radial_distance(z) * (1.0+z)

    def distance_modulus(self, z):
        """Return the distance modulus at the given redshift (Hogg Eq. 24)."""
        return 5. * np.log10(self.luminosity_distance(z)) + 25

    def universe_age(self, z):
        """Return the age of the universe at the given redshift.

        """
        from scipy.integrate import quad

        def _agefunc(z):
            return 1.0 / self.efunc(z) / (1.0 + z)
        
        if np.isscalar(z):
            integ, _ =  quad(_agefunc, z, self._z[-1])
            return integ * self.hubble_time
        else:
            age = []
            for _z in z:
                integ, _ =  quad(_agefunc, _z, self._z[-1])
                age.append(integ * self.hubble_time)
            return np.array(age)

#
# sigma clipping in Numba
# Rewrite basic sigma clipping to avoid
# array copies and redundant summation
# on each iteration
#
@numba.jit(nopython=True, nogil=True)
def sigmaclip(c, low=3., high=3.):

    n  = len(c)
    mask = np.full(n, True, dtype=np.bool_)
    
    s  = np.sum(c)
    s2 = np.sum(c*c)
    
    delta = 1
    while delta > 0:
        mean = s/n
        std  = np.sqrt(s2/n - mean*mean)
        clo = mean - std * low
        chi = mean + std * high
        
        n0 = n
        for j, cval in enumerate(c):
            if mask[j] and (cval < clo or cval > chi):
                mask[j] = False
                n  -= 1
                s  -= cval
                s2 -= cval * cval

        delta = n0-n
        
    return c[mask], clo, chi


def quantile(A, q):
    return _quantile(A, q)


def median(A):
    return _median(A)


# Numba's quantile impl is much faster
# than Numpy's standard version
@numba.jit(nopython=True, nogil=True)
def _quantile(A, q):
    return np.quantile(A, q)


# Numba's median impl is also faster
@numba.jit(nopython=True, nogil=True)
def _median(A):
    return np.median(A)
