import numpy as np

from numba import jit

from fastspecfit.util import C_LIGHT

from .utils import (
    MAX_SDEV,
    max_buffer_width,
    norm_cdf
)

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
#   line_parameters  -- parameters of each fitted line
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
                 line_parameters,
                 log_obs_bin_edges,
                 redshift,
                 ibin_widths):

    line_amplitudes, line_vshifts, line_sigmas = \
        np.split(line_parameters, 3)
        
    # output per-bin fluxes
    nbins = len(log_obs_bin_edges) - 1
    model_fluxes = np.zeros(nbins, dtype=line_amplitudes.dtype)

    # temporary buffer for per-line calculations, sized large
    # enough for whatever we may need to compute ( [s..e) )
    bin_vals = np.empty(nbins + 2, dtype=line_amplitudes.dtype)

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
        model_fluxes[s:e] += bin_vals[:e-s]
    
    return model_fluxes
    

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
                          line_parameters,
                          log_obs_bin_edges,
                          redshift,
                          ibin_widths,
                          padding):
    
    line_amplitudes, line_vshifts, line_sigmas = \
        np.split(line_parameters, 3)
    
    nbins = len(log_obs_bin_edges) - 1
    
    # buffers for per-line calculations, sized large
    # enough for max possible range [s .. e), plus extra
    # padding as directed by caller.
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

        endpts[j,0] = s
        endpts[j,1] = e
        
    return (endpts, line_profiles)


#
# emline_model_core()
# Compute the flux contribution of one spectral line to a model
# spectrum.  We compute the range of bins where the line
# contributes flux, then write those contributions to the 
# output array vals and return half-open range of bins
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
# NB: vals MUST have length at least two more than the largest number
# of bins that could be computed.  Entries in vals after the returned 
# range may be set to arbitrary values.
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

    # check if entire Gaussian is outside bounds of log_obs_bin_edges
    if hi == 0 or lo == len(log_obs_bin_edges): 
        return (0,0)
    
    nedges = hi - lo + 2  # compute values at edges [lo - 1 ... hi]
    
    A = c * line_amplitude * shifted_line
    offset = (log_shifted_line / sigma + sigma)  if sigma > 0. else 0.
    
    # vals[i] --> edge i + lo - 1
    
    vals[0] = 0. # edge lo - 1
    
    for i in range(1, nedges-1):
        
        # x = (log(lambda_j) - mu_i)/sigma - sigma,
        # the argument of the Gaussian integral
        
        x = log_obs_bin_edges[i+lo-1]/sigma - offset
        vals[i] = A * norm_cdf(x)
        
    vals[nedges-1] = A  # edge hi
    
    # Compute avg of bin i+lo-1 for 0 <= i < nedges - 1.
    # But:
    #  * if lo == 0, bin lo - 1 is not defined, so skip it
    #  * if hi == |log_obs_bin_edges|, bin hi - 1 is not defined,
    #      so skip it
    # Implement skips by modifying range of computation loop.
    #
    # After loop, vals contains weighted area of bins
    # lo - 1 + dlo .. hi - dhi, plus up to two cells of
    # trailing garbage.
    #
    # We ensure that values for all valid bins are written
    # starting at vals[0], and return the range of indices
    # of the computed bins w/r to the input bin array.
    
    dlo = 1 if lo == 0                      else 0
    dhi = 1 if hi == len(log_obs_bin_edges) else 0
    
    for i in range(dlo, nedges - 1 - dhi):
        vals[i-dlo] = (vals[i+1] - vals[i]) * ibin_widths[i+lo]
    
    return (lo - 1 + dlo, hi - dhi)
