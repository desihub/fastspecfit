#
# Calculation of objective function and
# Jacobian for emission line fitting
#

import numpy as np
from math import erf, erfc

from numba import jit

from .params_mapping import ParamsMapping
from .sparse_rep import EMLineJacobian

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

"""
#
# for debugging -- compute singular values of Jacobian
# and estimate its condition number.
#
def jacEstimateConditionNumber(J):

    from scipi.sparse.linalg import svds
    
    _, nFreeParms = J[0]
    
    try:
        svs = svds(J, return_singular_vectors=False,
                   k=nFreeParms - 1, which="LM")
        sv0 = svds(J, return_singular_vectors=False,
                   k=1, which="SM")[0]
        cond = svs[-1] / sv0
        print(np.hstack((sv0, svs)))       
        print(f"cond(J) = {cond:.3e}")
        
    except:
        print("Failed to compute Jacobian condition number")
"""

###################################################################

#
# Objective function for least-squares optimization
#
# Build the emline model as described above and compute the weighted
# vector of residuals between the modeled fluxes and the observations.
#
def objective(free_parameters,
              bin_data,
              obs_fluxes,
              obs_weights,
              redshift,
              line_wavelengths,
              resolution_matrices,
              camerapix,
              params_mapping):
    
    #
    # expand free paramters into complete
    # parameter array, handling tied params
    # and doublets
    #
    parameters = params_mapping.mapFreeToFull(free_parameters)
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)
    
    log_obs_bin_edges, ibin_widths = bin_data
    
    model_fluxes = np.empty_like(obs_fluxes, dtype=obs_fluxes.dtype)
    
    for icam, campix in enumerate(camerapix):
        
        # start and end for obs fluxes of camera icam
        s, e = campix

        # Actual bin widths are in ibw[1..e-s].
        # Setup guarantees that ibw[0] and
        # ibw[e-s+1] are dummy values
        ibw = ibin_widths[s:e+2]
        
        mf = emline_model(line_wavelengths,
                          lineamps, linevshifts, linesigmas,
                          log_obs_bin_edges[s+icam:e+icam+1],
                          redshift,
                          ibw)
        
        # convolve model with resolution matrix and store in
        # this camera's subrange of model_fluxes
        resolution_matrices[icam].matvec(mf, model_fluxes[s:e])
        
    return obs_weights * (model_fluxes - obs_fluxes) # residuals


#
# emline_model() 
#
# Given a fixed set of spectral lines and known redshift, and estimates
# for the amplitude, width, and velocity shift of each line, compute the
# average flux values observed in each of a series of spectral bins
# whose log wavelengths are bounded by log_obs_bin_edges.
#
# INPUTS:
#   line_wavelengths -- array of nominal wavelengths for all fitted lines
#   line_amplitudes  -- amplitude for each fitted line
#   line_vshifts     -- additional velocity shift for each fitted lines
#   line_sigmas      -- width of Gaussian profile for each fitted lines
#
#   log_obs_bin_edges -- natural logs of observed wavelength bin edges
#   redshift          -- red shift of observed spectrum
#   ibin_widths       -- one over widths of each observed wavelength bin
#
# RETURNS:
#   vector of average fluxes in each observed wavelength bin
#
@jit(nopython=True, fastmath=False, nogil=True)
def emline_model(line_wavelengths,
                 line_amplitudes, line_vshifts, line_sigmas,
                 log_obs_bin_edges,
                 redshift,
                 ibin_widths):
    
    # temporary buffer for per-line calculations, sized large
    # enough for whatever we may need to compute ( [s..e) )
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas)
    bin_vals = np.empty(max_width, dtype=line_amplitudes.dtype)
    
    # output per-bin fluxes
    # entry i corresponds bin i-1
    # entry 0 is a dummy in case s == -1
    # last entry is a dummy in case e == nbins + 1
    nbins = len(log_obs_bin_edges) - 1
    model_fluxes = np.zeros(nbins + 2, dtype=line_amplitudes.dtype)
    
    # compute total area of all Gaussians inside each bin.
    # For each Gaussian, we only compute area contributions
    # for bins where it is non-negligible.
    for j in range(len(line_wavelengths)):
        
        if line_amplitudes[j] == 0.:
            continue
        
        s, e = emline_model_core(line_wavelengths[j],
                                 line_amplitudes[j],
                                 line_vshifts[j],
                                 line_sigmas[j],
                                 log_obs_bin_edges,
                                 redshift,
                                 ibin_widths,
                                 bin_vals)
        
        # add bin avgs for this peak to bins [s, e)
        # (fluxes are offset by 1 to allow for s == -1
        # and extended to allow for e = nbins + 1;
        # these bin values may be ignored)
        model_fluxes[s+1:e+1] += bin_vals[:e-s]
        
    # trim off left and right dummy values before returning
    return model_fluxes[1:-1]


#
# emline_model_core()
# Compute the flux contribution of one spectral line to a model
# spectrum.  We compute the range of bins where the line
# contributes flux, then write those contributions to the 
# output array edge_vals and return half-open range of bins
# [s, e) to which it contributes.
#
# INPUTS:
#   line_wavelength -- nominal wavelength of line
#   line_amplitude  -- amplitude for spectral line
#   line_vshift     -- additional velocity shift for line
#   line_sigma      -- width of Gaussian profile for line
#
#   log_obs_bin_edges -- natural log of observed wavelength bin edges
#   redshift          -- red shift of observed spectrum
#   ibin_widths       -- one over widths of each observed wavelength bin
#
#   vals            -- array in which to record values of
#                      any nonzero bin fluxes for this line
#
# RETURNS:
#   tuple (s,e) such that nonzero fluxes occur in bins [s,e)
#    (if s == e, no nonzero fluxes were computed)
#   vals[0:e-s] contains these fluxes
#
# NB: s can be as low as -1, while e can be as high as nbins + 1.
# In these cases, the first (resp last) returned fluxes in
# vals are invalid and may be discarded. The caller should
# handle these edge cases,
# 
@jit(nopython=True, fastmath=False, nogil=True)
def emline_model_core(line_wavelength,
                      line_amplitude,
                      line_vshift,
                      line_sigma,
                      log_obs_bin_edges,
                      redshift,
                      ibin_widths,
                      vals):
    
    SQRT_2PI = np.sqrt(2 * np.pi)
    
    # line width
    sigma = line_sigma / C_LIGHT
    c = SQRT_2PI * sigma * np.exp(0.5 * sigma**2)
    
    # wavelength shift for spectral lines
    line_shift = 1. + redshift + line_vshift / C_LIGHT
    
    shifted_line     = line_wavelength * line_shift
    log_shifted_line = np.log(shifted_line)
    
    # leftmost edge i that needs a value (> 0) for this line 
    lo = np.searchsorted(log_obs_bin_edges,
                         log_shifted_line - MAX_SDEV * sigma,
                         side="left")
     
    # leftmost edge i that does *not* need a value (== 1) for this line
    hi = np.searchsorted(log_obs_bin_edges,
                         log_shifted_line + MAX_SDEV * sigma,
                         side="right")
    
    if hi == lo:  # entire Gaussian is outside bounds of log_obs_bin_edges
        return (0,0)
    
    nedges = hi - lo + 2  # compute values at edges [lo - 1 ... hi]
    
    A = c * line_amplitude * shifted_line
    offset = log_shifted_line / sigma + sigma
    
    # vals[i] --> edge i + lo - 1
    
    vals[0] = 0. # edge lo - 1
    
    for i in range(1, nedges-1):
        
        # x = (log(lambda_j) - mu_i)/sigma - sigma,
        # the argument of the Gaussian integral
        
        x = log_obs_bin_edges[i+lo-1]/sigma - offset
        vals[i] = A * norm_cdf(x)
        
    vals[nedges-1] = A  # edge hi
    
    # convert vals[i] to avg of bin i+lo-1 (last value is garbage)
    # we get values for bins lo-1 to hi-1 inclusive
    for i in range(nedges-1):
        vals[i] = (vals[i+1] - vals[i]) * ibin_widths[i+lo]
    
    return (lo - 1, hi)


##################################################################################

#
# compatibility entry point to compute modeled fluxes
# (FIXME: we should merge most of this code with the copy in objective())
#
def build_model(redshift,
                line_amplitudes,
                line_vshifts,
                line_sigmas,
                line_wavelengths,
                obs_bin_centers,
                resolution_matrices,
                camerapix):
    
    log_obs_bin_edges, ibin_widths = prepare_bins(obs_bin_centers, camerapix)

    # below here is common code with objective()
    model_fluxes = np.empty_like(obs_bin_centers, dtype=obs_bin_centers.dtype)
    
    for icam, campix in enumerate(camerapix):
        
        # start and end for obs fluxes of camera icam
        s, e = campix
        
        # Actual bin widths are in ibw[1..e-s].
        # Setup guarantees that ibw[0] and
        # ibw[e-s+1] are dummy values
        ibw = ibin_widths[s:e+2]
        
        mf = emline_model(line_wavelengths,
                          line_amplitudes, line_vshifts, line_sigmas,
                          log_obs_bin_edges[s+icam:e+icam+1],
                          redshift,
                          ibw)
        
        # convolve model with resolution matrix and store in
        # this camera's subrange of model_fluxes
        resolution_matrices[icam].matvec(mf, model_fluxes[s:e])
        
    return model_fluxes


    
#
# find_peak_amplitudes()
# Given fitted parameters for all emission lines, report for
# each line the largest flux that it contributes to any
# observed bin.
#
# INPUTS:
#  parameters -- full array of fitted line parameters
#  bin_data   -- preprocessed bin data in the same form
#                taken by the optimizer objective
#  redshift   -- object redshift
#  line_wavelengths -- wavelengths of all fitted lines
#  resolution_matrices -- list of sparse resolution matrices
#                         in same form taken by objective
#  camerapix  -- wavelength ranges for each camera
#
# RETURNS:
#   an array of maximum amplitudes for each line
#
def find_peak_amplitudes(parameters,
                         bin_data,
                         redshift,
                         line_wavelengths,
                         resolution_matrices,
                         camerapix):

    #
    # expand free paramters into complete
    # parameter array, handling tied params
    # and doublets
    #
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)
    
    log_obs_bin_edges, ibin_widths = bin_data
    
    max_amps = np.zeros_like(line_wavelengths, dtype=lineamps.dtype)
    
    for icam, campix in enumerate(camerapix):
        
        # start and end for obs fluxes of camera icam
        s, e = campix
        
        # Actual bin widths are in ibw[1..e-s].
        # Setup guarantees that ibw[0] and
        # ibw[e-s+1] are not out of bounds.
        ibw = ibin_widths[s:e+1]

        # compute model waveform for each spectral line
        line_models = emline_perline_models(line_wavelengths,
                                            lineamps, linevshifts, linesigmas,
                                            log_obs_bin_edges[s+icam:e+icam+1],
                                            redshift,
                                            ibw,
                                            resolution_matrices[icam].ndiag())
        
        # convolve each line's waveform with resolution matrix
        line_models = mulWMJ(np.ones(e - s),
                             resolution_matrices[icam].data,
                             line_models)
        
        # find highest flux for each line; if it's
        # bigger than any seen for that line so far,
        # update the line's global max
        update_line_maxima(max_amps, line_models)
        
    return max_amps


#
# update_line_maxima()
# Given an array of line waveforms and an array of
# maximum amplitudes observed for each line so far,
# if the waveform contains a higher amplitude than
# the previous maximum, update the maximum in place.
#
@jit(nopython=True, fastmath=False, nogil=True)
def update_line_maxima(max_amps, line_models):

    endpts, vals = line_models
    
    # find the highest flux for each peak; if it's
    # bigger than any seen so far, update global max
    for i in range(vals.shape[0]):
        ps, pe = endpts[i]
        if pe > ps:
            max_amps[i] = np.maximum(max_amps[i],
                                    np.max(vals[i,:pe-ps]))
            

#
# emline_perline_models()
# Given parameters for a set of lines, compute for each
# line individually its waveform in the range of bins
# described by log_bin_edges.  The waveforms are returned
# sparsely in the same forma as the Jacobian below.
#
# This function shares the core computation of emline_model()
# but does not collapse the lines into one composite.
#
@jit(nopython=True, fastmath=False, nogil=True)
def emline_perline_models(line_wavelengths,
                          line_amplitudes, line_vshifts, line_sigmas,
                          log_obs_bin_edges,
                          redshift,
                          ibin_widths,
                          padding):
    
    nbins = len(log_obs_bin_edges) - 1
    
    # buffers for per-line calculations, sized large
    # enough for max possible range [s .. e)
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, padding)
    
    nlines = len(line_wavelengths)
    line_profiles = np.empty((nlines, max_width), dtype=line_amplitudes.dtype)
    endpts        = np.zeros((nlines, 2), dtype=np.int32)
    
    # compute total area of all Gaussians inside each bin.
    # For each Gaussian, we only compute area contributions
    # for bins where it is non-negligible.
    for j in range(nlines):

        if line_amplitudes[j] == 0.:
            continue
        
        bin_vals = line_profiles[j]
        
        s, e = emline_model_core(line_wavelengths[j],
                                 line_amplitudes[j],
                                 line_vshifts[j],
                                 line_sigmas[j],
                                 log_obs_bin_edges,
                                 redshift,
                                 ibin_widths,
                                 bin_vals)
        
        if s == -1: 
            # bin s is before start of requested bins,
            # and its true left endpt value is unknown
            s = 0
            
            # discard bin lo - 1 in bin_vals
            bin_vals[0:e-s-1] = bin_vals[1:e-s]
            
        if e == nbins + 1:
            # bin e-1 is one past end of requested bins,
            # and its true right endpt value is unknown
            e -= 1

        endpts[j] = np.array([s,e])

    return (endpts, line_profiles)


##########################################################################

#
# Jacobian of objective function for least-squares EMLine
# optimization. The result of the detailed calculation is converted
# to a sparse matrix, since it is extremely sparse, to speed up
# subsequent matrix-vector multiplies in the optimizer.
#
def jacobian(free_parameters,
             bin_data,
             obs_fluxes, # not used, but must match objective
             obs_weights,
             redshift,
             line_wavelengths,
             resolution_matrices,
             camerapix,
             params_mapping):
    
    #
    # expand free paramters into complete
    # parameter array, handling tied params
    # and doublets
    #
    
    parameters = params_mapping.mapFreeToFull(free_parameters)
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)

    log_obs_bin_edges, ibin_widths = bin_data
    
    J_S = params_mapping.getJacobian(free_parameters)

    jacs = []
    for icam, campix in enumerate(camerapix):
        s, e = campix

        # Actual bin widths are in ibw[1..e-s].
        # Setup guarantees that ibw[0] and
        # ibw[e-s+1] are not out of bounds.
        ibw = ibin_widths[s:e+1]

        idealJac = \
            emline_model_jacobian(lineamps, linevshifts, linesigmas,
                                  log_obs_bin_edges[s+icam:e+icam+1],
                                  ibw,
                                  redshift,
                                  line_wavelengths,
                                  resolution_matrices[icam].ndiag())
        
        # ignore any columns corresponding to fixed parameters
        endpts = idealJac[0]
        endpts[params_mapping.fixedMask(), :] = (0,0)
        
        jacs.append( mulWMJ(obs_weights[s:e],
                            resolution_matrices[icam].data,
                            idealJac) )
    
    nBins = np.sum(np.diff(camerapix))
    nFreeParms = len(free_parameters)
    nParms = len(parameters)
    J =  EMLineJacobian((nBins, nFreeParms), nParms,
                        camerapix, jacs, J_S)

    return J

    
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
                    
        if hi == lo:  # Gaussian is entirely outside bounds of log_obs_bin_edges
            continue

        nedges = hi - lo + 2 # compute values at edges [lo - 1 ... hi]
        
        # Compute contribs of each line to each partial derivative in place.
        # No sharing of params between peaks means that we never have to
        # add contributions from two peaks to same line.
        dda_vals = dd[           j]
        ddv_vals = dd[nlines   + j]
        dds_vals = dd[2*nlines + j]
        
        offset = log_shifted_line / sigma + sigma
        
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

##############################################################################

#
# mulWMJ()
# Compute the sparse matrix product P = WMJ, where
#   W is a diagonal matrix (represented by a vector)
#   M is a resolution matrix (the .data field of a
#     ResMatrix, since Numba can't handle classes)
#   J is a column-sparse matrix giving the nonzero
#     values in one contiguous range per column
#
# We assume that J jac has enough free space in
# each column to store the results of the multiply
# and update it in place, returning the same
# jac as our return value.
#
@jit(nopython=True, fastmath=False, nogil=True)
def mulWMJ(w, M, jac):

    endpts, J = jac
    
    nbins, ndiag = M.shape
    ncol, maxColSize = J.shape
    
    hdiag = ndiag//2

    # temporary buffer for each column of WMJ
    buf = np.empty(maxColSize, dtype=J.dtype)
        
    for j in range(ncol):
        # boundaries of nonzero entries
        # in jth column of J
        s, e = endpts[j]
        
        if s == e: # no nonzero values in column j
            continue
        
        # boundaries of entries in jth column of P
        # impacted by matrix multiply
        imin = np.maximum(s - hdiag, 0)
        imax = np.minimum(e + hdiag, nbins) # one past last impacted entry
        
        for i in range(imin, imax):
            
            # boundaries of interval of k where both
            # M[i, k] and J[k, j] are nonzero.
            kmin = np.maximum(i - hdiag,     s)
            kmax = np.minimum(i + hdiag, e - 1)
            
            acc = 0.
            for k in range(kmin, kmax + 1):
                acc += M[i, k - i + hdiag] * J[j, k - s]
            buf[i - imin] = acc * w[i]

        # write result back to J and update endpts for col
        newS = np.maximum(imin, 0)
        newE = np.minimum(imax, nbins)
        J[j, :newE - newS] = buf[:newE - newS]
        endpts[j] = np.array([newS, newE])
    
    return (endpts, J)


###############################################################################

#
# prepare_bins
# Convert N bin centers to the info needed by the optimizer,
# returned as an opaque tuple.  Given N bins, return values
# include
#
# - N+1 log bin edges.  Edges are placed halfway between centers,
# with extrapolation at the ends.  We return the natural log of
# each edge's wavelength, since that is what model and Jacobian
# building need.
# - N inverse bin widths, in an array padded by one cell on
#   the left and right.
#
@jit(nopython=True, fastmath=False, nogil=True)
def prepare_bins(centers, camerapix):

    ncameras = camerapix.shape[0]
    edges = np.empty(len(centers) + ncameras, dtype=centers.dtype)
    
    for icam, campix in enumerate(camerapix):
        
        s, e = campix
        icenters = centers[s:e]
        
        #- interior edges are just points half way between bin centers
        int_edges = 0.5 * (icenters[:-1] + icenters[1:])
        
        #- exterior edges are extrapolation of interior bin sizes
        edge_l = icenters[ 0] - (icenters[ 1] - int_edges[ 0])
        edge_r = icenters[-1] + (icenters[-1] - int_edges[-1])

        edges[s + icam]              = edge_l
        edges[s + icam + 1:e + icam] = int_edges
        edges[e + icam]              = edge_r

    # dummies before and after widths are needed
    # for corner cases in edge -> bin computation
    ibin_widths = np.empty(len(centers) + 2, dtype=centers.dtype)
    ibin_widths[0]  = 0.
    for i in range(len(centers)):
        ibin_widths[i+1] = 1. / (edges[i+1] - edges[i])
    ibin_widths[-1] = 0.
    
    return (np.log(edges), ibin_widths)
