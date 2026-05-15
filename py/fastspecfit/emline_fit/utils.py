import numpy as np
from math import erf, erfc

from numba import jit

from fastspecfit.util import C_LIGHT


# Do not bother computing normal PDF/CDF if more than this many
# standard deviations from mean.
MAX_SDEV = 5.


@jit(nopython=True, nogil=True, cache=True)
def norm_pdf(a):
    """PDF of the standard normal distribution.

    Parameters
    ----------
    a : float
        Evaluation point.

    Returns
    -------
    float
        Probability density at ``a``.

    """

    SQRT_2PI = np.sqrt(2 * np.pi)

    return 1/SQRT_2PI * np.exp(-0.5 * a**2)


@jit(nopython=True, nogil=True, cache=True)
def norm_cdf(a):
    """Approximate the CDF of the standard normal distribution.

    Parameters
    ----------
    a : float
        Upper integration limit.

    Returns
    -------
    float
        Integral of the standard normal PDF from :math:`-\\infty` to ``a``.

    """

    SQRT1_2 = 1.0 / np.sqrt(2)

    z = np.abs(a)

    # Optimization (currently disabled because it is not needed): If
    # |a| > MAX_SDEV, treat the value as extreme and return 0 or 1 as
    # appropriate.

    #if z > MAX_SDEV: # short-circuit extreme values
    #    if a > 0:
    #        y = 1.
    #    else:
    #        y = 0.
    if z < 1.:
        y = 0.5 + 0.5 * erf(a * SQRT1_2)
    else:
        y = 0.5 * erfc(z * SQRT1_2)
        if a > 0:
            y = 1.0 - y

    return y


@jit(nopython=True, nogil=True, cache=True)
def max_buffer_width(log_obs_bin_edges, line_sigmas, padding=0):
    """Estimate the maximum number of nonzero bins possible for any spectral line.

    Parameters
    ----------
    log_obs_bin_edges : :class:`np.ndarray`
        Log wavelengths of all observed bin edges.
    line_sigmas : :class:`np.ndarray`
        Gaussian widths of all spectral lines in km/s.
    padding : int, optional
        Extra bins to add on each side for future expansion. Defaults to 0.

    Returns
    -------
    int
        Safe upper bound on the number of bins with nonzero flux for any line.

    """

    # Find the largest width sigma for any line, and
    # allocate enough space for twice that much width
    # in bins, given the smallest observed bin width.
    # Add padding and a little fudge factor to be safe.
    max_width = \
        int(2*MAX_SDEV*np.max(line_sigmas/C_LIGHT) / \
            np.min(np.diff(log_obs_bin_edges))) + \
            2*padding + 4
    return max_width
