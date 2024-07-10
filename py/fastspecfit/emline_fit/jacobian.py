import numpy as np

from numba import jit

from fastspecfit.util import C_LIGHT

from .utils import (
    MAX_SDEV,
    max_buffer_width,
    norm_cdf,
    norm_pdf
)

#
# emline_model_jacobian() 
#
# Compute the Jacobian of the function computed in emlines_model().
# Inputs are as for emlines_model().
#
# RETURNS:
# Sparse Jacobian as tuple (endpts, dd), where
#  column j has nonzero values in interval [ endpts[j,0] , endpts[j,1] )
#  which are stored in dd[j].
#
@jit(nopython=True, fastmath=False, nogil=True)
def emline_model_jacobian(line_amplitudes, line_vshifts, line_sigmas,
                          log_obs_bin_edges,
                          ibin_widths,
                          redshift,
                          line_wavelengths,
                          padding):

    SQRT_2PI = np.sqrt(2*np.pi)

    nbins = len(log_obs_bin_edges) - 1

    # buffers for per-parameter calculations, sized large
    # enough for line's max possible range [s .. e)
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, padding)
    
    nlines = len(line_wavelengths)
    dd     = np.empty((3 * nlines, max_width), dtype=line_amplitudes.dtype)
    endpts = np.zeros((nlines, 2), dtype=np.int32)

    starts = endpts[:,0]
    ends   = endpts[:,1]
    
    # compute partial derivatives for avg values of all Gaussians
    # inside each bin. For each Gaussian, we only compute
    # contributions for bins where it is non-negligible.
    for j in range(len(line_wavelengths)):
        
        # line width
        sigma = line_sigmas[j] / C_LIGHT
        
        c0 = SQRT_2PI * np.exp(0.5 * sigma**2)
        
        # wavelength shift for spectral lines
        line_shift = 1. + redshift + line_vshifts[j] / C_LIGHT
        shifted_line     = line_wavelengths[j] * line_shift
        log_shifted_line = np.log(shifted_line)
        
        # leftmost edge i that needs a value (> 0) for this line 
        lo = np.searchsorted(log_obs_bin_edges,
                             log_shifted_line - MAX_SDEV * sigma,
                             side="left")
        
        # leftmost edge i that does *not* need a value (== 1) for this line
        hi = np.searchsorted(log_obs_bin_edges,
                             log_shifted_line + MAX_SDEV * sigma,
                             side="right")

        # check if entire Gaussian is outside bounds of log_obs_bin_edges
        if hi == 0 or lo == len(log_obs_bin_edges): 
            continue
        
        nedges = hi - lo + 2 # compute values at edges [lo - 1 ... hi]
        
        # Compute contribs of each line to each partial derivative in place.
        # No sharing of params between peaks means that we never have to
        # add contributions from two peaks to same line.
        dda_vals = dd[           j]
        ddv_vals = dd[nlines   + j]
        dds_vals = dd[2*nlines + j]

        offset = (log_shifted_line / sigma + sigma) if sigma > 0. else 0.
        
        c = c0 * line_wavelengths[j]
        A = c / C_LIGHT * line_amplitudes[j]
        
        # vals[i] --> edge i + lo - 1
        
        dda_vals[0] = 0. # edge lo - 1
        ddv_vals[0] = 0.
        dds_vals[0] = 0.
        
        for i in range(1, nedges - 1):
            
            # x - offset = (log(lambda_j) - mu_i)/sigma - sigma,
            # the argument of the Gaussian integral
            
            x = log_obs_bin_edges[i+lo-1]/sigma - offset
            pdf = norm_pdf(x)
            cdf = norm_cdf(x)
            
            dda_vals[i] = c * line_shift * sigma * cdf
            ddv_vals[i] = A * (sigma * cdf - pdf)
            dds_vals[i] = A * line_shift * \
                ((1 + sigma**2) * cdf - (x + 2*sigma) * pdf)
            
        dda_vals[nedges - 1] = c * line_shift * sigma     # edge hi
        ddv_vals[nedges - 1] = A * sigma
        dds_vals[nedges - 1] = A * line_shift * (1 + sigma**2)

        # convert *_vals[i] to partial derivatives for bin i+lo-1
        # (last value in each array is garbage)
        # we get values for bins lo-1 to hi-1 inclusive
        for i in range(nedges - 1):
            dda_vals[i] = (dda_vals[i+1] - dda_vals[i]) * ibin_widths[i+lo]
            ddv_vals[i] = (ddv_vals[i+1] - ddv_vals[i]) * ibin_widths[i+lo]
            dds_vals[i] = (dds_vals[i+1] - dds_vals[i]) * ibin_widths[i+lo]
        
        # starts[j] is set to first valid bin
        if lo == 0:
            # bin low - 1 is before start of requested bins,
            # and its true left endpt value is unknown
            starts[j] = lo

            # discard bin lo - 1 in dd*_vals
            dda_vals[0:nedges - 1] = dda_vals[1:nedges]
            ddv_vals[0:nedges - 1] = ddv_vals[1:nedges]
            dds_vals[0:nedges - 1] = dds_vals[1:nedges]
        else:
            # bin lo - 1 is valid
            starts[j] = lo - 1
        
        # ends[j] is set one past last valid bin
        if hi == nbins + 1:
            # bin hi - 1 is one past end of requested bins,
            # and its true right endpt value is unknown
            ends[j] = hi - 1
        else:
            # bin hi - 1 is valid
            ends[j] = hi

    all_endpts = tile_2d(endpts, 3)
    
    # for lines with zero amplitude, partial derivatives
    # w/r to their vshifts and sigmas are zero.
    for i, amp in enumerate(line_amplitudes):
        if amp == 0.:
            all_endpts[i + nlines  ] = np.array([0, 0])
            all_endpts[i + 2*nlines] = np.array([0, 0])

    # for lines with zero width, partial derivatives
    # w/r to their amplitudes are zero.
    for i, sig in enumerate(line_sigmas):
        if sig == 0.:
            all_endpts[i] = np.array([0, 0])
    
    return (all_endpts, dd)    


# horizontally tile a 2D array n times
# replaces np.tile, which is not supported by Numba,
@jit(nopython=True, fastmath=False, nogil=True)
def tile_2d(a, n):
    sz = a.shape[0]
    r = np.empty((n * sz, a.shape[1]), dtype=a.dtype)
    for i in range(n):
        r[i*sz:(i+1)*sz,:] = a
    return r


# compute partial Jacobian associated with just the patch parameters
# in the sparse column form used in sparse_rep.py.  This Jacobian
# is independent of both the line model and the particular choices
# of slope/intercept for each patch, so can be computed once when
# we set up the optimization.
#
# FIXME: do we need to do this for each camera and apply resolution?
@staticmethod
@jit(nopython=True, fastmath=False, nogil=True)
def patch_jacobian(obs_bin_centers,
                   obs_weights,
                   patch_endpts,
                   patch_pivotwave):
    
    nPatches = patch_endpts.shape[0]

    #
    # make two copies of the endpts array (one each
    # for slopes and intercepts) and compute the
    # maximum width of any patch.
    #
    
    endpts = np.empty((2*nPatches, 2), dtype=np.int32)
    maxPatchWidth = 0

    for i in range(nPatches):
        s, e = patch_endpts[i]
        
        endpts[i]          = (s, e)
        endpts[i+nPatches] = (s, e)
        
        maxPatchWidth = np.maximum(maxPatchWidth, e - s)

    #
    # Compute weighted partial derivatives of
    # flux w/r to slope and intercept of each
    # patch.  These derivatives are nonzero only
    # within the boundaries of the patch.
    #
    M = np.empty((2*nPatches, maxPatchWidth))
    for i in range(nPatches):
        s, e = endpts[i]

        # dobj/dslope for patch
        M[i,:e-s] = \
            (obs_bin_centers[s:e] - patch_pivotwave[i]) * \
            obs_weights[s:e]
        
        # dobj/dintercept for patch
        M[i + nPatches, :e-s] = obs_weights[s:e] # 1. x obs_weights
        
    return (endpts, M)
