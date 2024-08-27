import numpy as np
from math import erf, erfc

from numba import jit

from fastspecfit.util import C_LIGHT

# Do not bother computing normal PDF/CDF if more than this many 
# standard deviations from mean.
MAX_SDEV = 5.

#
# norm_pdf()
# PDF of standard normal distribution at a point a
#
@jit(nopython=True, fastmath=False, nogil=True)
def norm_pdf(a):

    SQRT_2PI = np.sqrt(2 * np.pi)
    
    return 1/SQRT_2PI * np.exp(-0.5 * a**2)

#
# norm_cdf()
# Approximate the integral of a standard normal PDF from -infty to a.
#
# Optimization (currently disabled because it is not needed): If
# |a| > MAX_SDEV, treat the value as extreme and return 0 or 1 as
# appropriate.
#
@jit(nopython=True, fastmath=False, nogil=True)
def norm_cdf(a):

    SQRT1_2 = 1.0 / np.sqrt(2)
    
    z = np.abs(a)

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


#
# max_buffer_width()
# Compute a safe estimate of the number of nonzero bin fluxes possible
# for a line spanning a subrange of bins with edges log_obs_bin_edges,
# assuming the line's width is one of the values in line_sigmas.
# Optionally add 2*padding to allow future expansion to left and right.
#
@jit(nopython=True, fastmath=False, nogil=True)
def max_buffer_width(log_obs_bin_edges, line_sigmas, padding=0):

    # Find the largest width sigma for any line, and
    # allocate enough space for twice that much width
    # in bins, given the smallest observed bin width.
    # Add padding and a little fudge factor to be safe.
    max_width = \
        int(2*MAX_SDEV*np.max(line_sigmas/C_LIGHT) / \
            np.min(np.diff(log_obs_bin_edges))) + \
            2*padding + 4
    return max_width
