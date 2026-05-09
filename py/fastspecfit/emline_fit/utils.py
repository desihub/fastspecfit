import numpy as np
from math import erf, erfc

from numba import jit

from fastspecfit.util import C_LIGHT


# Do not bother computing normal PDF/CDF if more than this many
# standard deviations from mean.
MAX_SDEV = 5.


@jit(nopython=True, nogil=True, cache=True)
def norm_pdf(a):
    """
    PDF of standard normal distribution at a point a
    """

    SQRT_2PI = np.sqrt(2 * np.pi)

    return 1/SQRT_2PI * np.exp(-0.5 * a**2)


@jit(nopython=True, nogil=True, cache=True)
def norm_cdf(a):
    """
    Approximate the integral of a standard normal PDF from -infty to a.
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
def max_buffer_width(log_obs_bin_edges, line_sigmas, sigma0_angstrom, padding=0):
    # sigma_eff >= sigma0; use sigma0 as a floor so narrow/dropped lines
    # still get enough buffer for the fiducial pre-convolved profile.
    min_log_bin = np.min(np.diff(log_obs_bin_edges))
    max_sigma_line = np.max(line_sigmas / C_LIGHT)

    # conservative sigma_eff: quadrature sum at the blue end (smallest lambda),
    # where sigma0 in log-lambda is largest
    sigma0_max = sigma0_angstrom / np.exp(log_obs_bin_edges[0])
    sigma_eff_max = np.sqrt(max_sigma_line**2 + sigma0_max**2)

    max_width = int(2 * MAX_SDEV * sigma_eff_max / min_log_bin) + 2 * padding + 4
    return max_width
