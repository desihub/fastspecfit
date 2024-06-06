"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import pdb # for debugging

import os, time
import numpy as np
from numba import jit
from math import erf, erfc
from astropy.table import Table

from fastspecfit.io import FLUXNORM
from fastspecfit.continuum import Filters
from fastspecfit.util import C_LIGHT

# Do not bother computing normal PDF/CDF if more than this many 
# standard deviations from mean.
MAX_SDEV = 5.


def read_emlines(emlinesfile=None):
    """Read the set of emission lines of interest.

    """
    if emlinesfile is None:
        from importlib import resources
        emlinesfile = resources.files('fastspecfit').joinpath('data/emlines.ecsv')

    try:
        linetable = Table.read(emlinesfile, format='ascii.ecsv', guess=False)
    except: 
        from desiutil.log import get_logger
        log = get_logger()
        errmsg = f'Problem reading emission lines parameter file {emlinesfile}.'
        log.critical(errmsg)
        raise ValueError(errmsg)
    
    return linetable    


@jit(nopython=True, fastmath=False, nogil=True)
def norm_pdf(a):
    """PDF of standard normal distribution at a point a.

    """
    SQRT_2PI = np.sqrt(2 * np.pi)
    
    return 1/SQRT_2PI * np.exp(-0.5 * a**2)


@jit(nopython=True, fastmath=False, nogil=True)
def norm_cdf(a):
    """Approximate the integral of a standard normal PDF from -infty to a.
    
    Optimization (currently disabled because it is not needed): If |a| >
    MAX_SDEV, treat the value as extreme and return 0 or 1 as appropriate.

    """
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


@jit(nopython=True, fastmath=False, nogil=True)
def max_buffer_width(log_obs_bin_edges, line_sigmas, padding=0):
    """Compute a safe estimate of the number of nonzero bin fluxes possible for a
    line spanning a subrange of bins with edges log_obs_bin_edges, assuming the
    line's width is one of the values in line_sigmas.  Optionally add 2*padding
    to allow future expansion to left and right.

    """

    # Find the largest width sigma for any line, and allocate enough space for
    # twice that much width in bins, given the smallest observed bin width.  Add
    # padding and a little fudge factor to be safe.
    max_width = \
        int(2*MAX_SDEV*np.max(line_sigmas/C_LIGHT) / \
            np.min(np.diff(log_obs_bin_edges))) + \
            2*padding + 4
    return max_width


def objective(free_parameters,
              bin_data,
              obs_fluxes,
              obs_weights,
              redshift,
              line_wavelengths,
              resolution_matrices,
              camerapix,
              params_mapping):
    """Objective function for least-squares optimization
    
    Build the emline model as described above and compute the weighted vector of
    residuals between the modeled fluxes and the observations.

    """

    # Expand free paramters into complete parameter array, handling tied params
    # and doublets.
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


@jit(nopython=True, fastmath=False, nogil=True)
def emline_model(line_wavelengths,
                 line_amplitudes, line_vshifts, line_sigmas,
                 log_obs_bin_edges,
                 redshift,
                 ibin_widths):
    """Given a fixed set of spectral lines and known redshift, and estimates for
    the amplitude, width, and velocity shift of each line, compute the average
    flux values observed in each of a series of spectral bins whose log
    wavelengths are bounded by log_obs_bin_edges.
   
    INPUTS:
      line_wavelengths -- array of nominal wavelengths for all fitted lines
      line_amplitudes  -- amplitude for each fitted line
      line_vshifts     -- additional velocity shift for each fitted lines
      line_sigmas      -- width of Gaussian profile for each fitted lines
   
      log_obs_bin_edges -- natural logs of observed wavelength bin edges
      redshift          -- red shift of observed spectrum
      ibin_widths       -- one over widths of each observed wavelength bin
   
    RETURNS:
      vector of average fluxes in each observed wavelength bin

    """
    # Temporary buffer for per-line calculations, sized large enough for
    # whatever we may need to compute ( [s..e) ).
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas)
    bin_vals = np.empty(max_width, dtype=line_amplitudes.dtype)
    
    # output per-bin fluxes
    #   entry i corresponds bin i-1
    #   entry 0 is a dummy in case s == -1
    #   last entry is a dummy in case e == nbins + 1
    nbins = len(log_obs_bin_edges) - 1
    model_fluxes = np.zeros(nbins + 2, dtype=line_amplitudes.dtype)
    
    # Compute total area of all Gaussians inside each bin.  For each Gaussian,
    # we only compute area contributions for bins where it is non-negligible.
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
        
        # Add bin avgs for this peak to bins [s, e). Fluxes are offset by 1 to
        # allow for s == -1 and extended to allow for e = nbins + 1; these bin
        # values may be ignored.
        model_fluxes[s+1:e+1] += bin_vals[:e-s]
        
    # Trim off left and right dummy values before returning.
    return model_fluxes[1:-1]


@jit(nopython=True, fastmath=False, nogil=True)
def emline_model_core(line_wavelength,
                      line_amplitude,
                      line_vshift,
                      line_sigma,
                      log_obs_bin_edges,
                      redshift,
                      ibin_widths,
                      vals):
    """Compute the flux contribution of one spectral line to a model spectrum.  We
    compute the range of bins where the line contributes flux, then write those
    contributions to the output array edge_vals and return half-open range of
    bins [s, e) to which it contributes.
    
    INPUTS:
      line_wavelength -- nominal wavelength of line
      line_amplitude  -- amplitude for spectral line
      line_vshift     -- additional velocity shift for line
      line_sigma      -- width of Gaussian profile for line
    
      log_obs_bin_edges -- natural log of observed wavelength bin edges
      redshift          -- red shift of observed spectrum
      ibin_widths       -- one over widths of each observed wavelength bin
    
      vals            -- array in which to record values of
                         any nonzero bin fluxes for this line
    
    RETURNS:
      tuple (s,e) such that nonzero fluxes occur in bins [s,e)
       (if s == e, no nonzero fluxes were computed)
      vals[0:e-s] contains these fluxes
    
    NB: s can be as low as -1, while e can be as high as nbins + 1.
    In these cases, the first (resp last) returned fluxes in
    vals are invalid and may be discarded. The caller should
    handle these edge cases,

    """    
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


def build_model(redshift,
                line_amplitudes,
                line_vshifts,
                line_sigmas,
                line_wavelengths,
                obs_bin_centers,
                resolution_matrices,
                camerapix):
    """Compatibility entry point to compute modeled fluxes.  (FIXME: we should
    merge most of this code with the copy in objective())

    """
    log_obs_bin_edges, ibin_widths = prepare_bins(obs_bin_centers, camerapix)

    # Below here is common code with objective().
    model_fluxes = np.empty_like(obs_bin_centers, dtype=obs_bin_centers.dtype)
    
    for icam, campix in enumerate(camerapix):
        
        # start and end for obs fluxes of camera icam
        s, e = campix
        
        # Actual bin widths are in ibw[1..e-s].  Setup guarantees that ibw[0]
        # and ibw[e-s+1] are dummy values.
        ibw = ibin_widths[s:e+2]
        
        mf = emline_model(line_wavelengths,
                          line_amplitudes, line_vshifts, line_sigmas,
                          log_obs_bin_edges[s+icam:e+icam+1],
                          redshift,
                          ibw)
        
        # Convolve model with resolution matrix and store in this camera's
        # subrange of model_fluxes.
        resolution_matrices[icam].matvec(mf, model_fluxes[s:e])
        
    return model_fluxes


def find_peak_amplitudes(parameters,
                         bin_data,
                         redshift,
                         line_wavelengths,
                         resolution_matrices,
                         camerapix):
    """Given fitted parameters for all emission lines, report for each line the
    largest flux that it contributes to any observed bin.
    
    INPUTS:
     parameters -- full array of fitted line parameters
     bin_data   -- preprocessed bin data in the same form
                   taken by the optimizer objective
     redshift   -- object redshift
     line_wavelengths -- wavelengths of all fitted lines
     resolution_matrices -- list of sparse resolution matrices
                            in same form taken by objective
     camerapix  -- wavelength ranges for each camera
    
    RETURNS:
      an array of maximum amplitudes for each line

    """    
    # Expand free paramters into complete parameter array, handling tied params
    # and doublets.
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)
    
    log_obs_bin_edges, ibin_widths = bin_data
    
    max_amps = np.zeros_like(line_wavelengths, dtype=lineamps.dtype)
    
    for icam, campix in enumerate(camerapix):
        
        # start and end for obs fluxes of camera icam
        s, e = campix
        
        # Actual bin widths are in ibw[1..e-s]. Setup guarantees that ibw[0] and
        # ibw[e-s+1] are not out of bounds.
        ibw = ibin_widths[s:e+1]

        # Compute model waveform for each spectral line.
        line_models = emline_perline_models(line_wavelengths,
                                            lineamps, linevshifts, linesigmas,
                                            log_obs_bin_edges[s+icam:e+icam+1],
                                            redshift,
                                            ibw,
                                            resolution_matrices[icam].ndiag())
        
        # Convolve each line's waveform with resolution matrix.
        line_models = mulWMJ(np.ones(e - s),
                             resolution_matrices[icam].data,
                             line_models)
        
        # Find highest flux for each line; if it's bigger than any seen for that
        # line so far, update the line's global max.
        update_line_maxima(max_amps, line_models)
        
    return max_amps


@jit(nopython=True, fastmath=False, nogil=True)
def update_line_maxima(max_amps, line_models):
    """Given an array of line waveforms and an array of maximum amplitudes observed
    for each line so far, if the waveform contains a higher amplitude than the
    previous maximum, update the maximum in place.

    """
    endpts, vals = line_models
    
    # Find the highest flux for each peak; if it's bigger than any seen so far,
    # update global max.
    for i in range(vals.shape[0]):
        ps, pe = endpts[i]
        if pe > ps:
            max_amps[i] = np.maximum(max_amps[i], np.max(vals[i,:pe-ps]))
            

@jit(nopython=True, fastmath=False, nogil=True)
def emline_perline_models(line_wavelengths,
                          line_amplitudes, line_vshifts, line_sigmas,
                          log_obs_bin_edges,
                          redshift,
                          ibin_widths,
                          padding):
    """Given parameters for a set of lines, compute for each line individually its
    waveform in the range of bins described by log_bin_edges.  The waveforms are
    returned sparsely in the same forma as the Jacobian below.

    This function shares the core computation of emline_model() but does not
    collapse the lines into one composite.

    """
    nbins = len(log_obs_bin_edges) - 1
    
    # Buffers for per-line calculations, sized large enough for max possible
    # range [s .. e).
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, padding)
    
    nlines = len(line_wavelengths)
    line_profiles = np.empty((nlines, max_width), dtype=line_amplitudes.dtype)
    endpts        = np.zeros((nlines, 2), dtype=np.int32)
    
    # Compute total area of all Gaussians inside each bin.  For each Gaussian,
    # we only compute area contributions for bins where it is non-negligible.
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
            # Bin s is before start of requested bins, and its true left endpt
            # value is unknown.
            s = 0
            
            # discard bin lo - 1 in bin_vals
            bin_vals[0:e-s-1] = bin_vals[1:e-s]
            
        if e == nbins + 1:
            # bin e-1 is one past end of requested bins,
            # and its true right endpt value is unknown
            e -= 1

        endpts[j] = np.array([s,e])

    return (endpts, line_profiles)


def jacobian(free_parameters,
             bin_data,
             obs_fluxes, # not used, but must match objective
             obs_weights,
             redshift,
             line_wavelengths,
             resolution_matrices,
             camerapix,
             params_mapping):
    """
    Jacobian of objective function for least-squares EMLine
    optimization. The result of the detailed calculation is converted
    to a sparse matrix, since it is extremely sparse, to speed up
    subsequent matrix-vector multiplies in the optimizer.

    """
    from fastspecfit.fitting import EMLineJacobian
    
    # Expand free paramters into complete parameter array, handling tied params
    # and doublets.
    parameters = params_mapping.mapFreeToFull(free_parameters)
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)

    log_obs_bin_edges, ibin_widths = bin_data
    
    J_S = params_mapping.getJacobian(free_parameters)

    jacs = []
    for icam, campix in enumerate(camerapix):
        s, e = campix

        # Actual bin widths are in ibw[1..e-s].  Setup guarantees that ibw[0]
        # and ibw[e-s+1] are not out of bounds.
        ibw = ibin_widths[s:e+1]

        idealJac = \
            emline_model_jacobian(lineamps, linevshifts, linesigmas,
                                  log_obs_bin_edges[s+icam:e+icam+1],
                                  ibw,
                                  redshift,
                                  line_wavelengths,
                                  resolution_matrices[icam].ndiag())
        
        # Ignore any columns corresponding to fixed parameters.
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


@jit(nopython=True, fastmath=False, nogil=True)
def emline_model_jacobian(line_amplitudes, line_vshifts, line_sigmas,
                          log_obs_bin_edges,
                          ibin_widths,
                          redshift,
                          line_wavelengths,
                          padding):
    """
    Compute the Jacobian of the function computed in emlines_model().
    Inputs are as for emlines_model().
    
    RETURNS:
    Sparse Jacobian as tuple (endpts, dd), where
     column j has nonzero values in interval [ endpts[j,0] , endpts[j,1] )
     which are stored in dd[j].

    """
    SQRT_2PI = np.sqrt(2*np.pi)

    nbins = len(log_obs_bin_edges) - 1

    # Buffers for per-parameter calculations, sized large enough for line's max
    # possible range [s .. e).
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, padding)
    
    nlines = len(line_wavelengths)
    dd     = np.empty((3 * nlines, max_width), dtype=line_amplitudes.dtype)
    endpts = np.zeros((nlines, 2), dtype=np.int32)

    starts = endpts[:,0]
    ends   = endpts[:,1]
    
    # Compute partial derivatives for avg values of all Gaussians inside each
    # bin. For each Gaussian, we only compute contributions for bins where it is
    # non-negligible.
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


@jit(nopython=True, fastmath=False, nogil=True)
def tile_2d(a, n):
    """Horizontally tile a 2D array n times # replaces np.tile, which is supported
    by Numba.

    """
    sz = a.shape[0]
    r = np.empty((n * sz, a.shape[1]), dtype=a.dtype)
    for i in range(n):
        r[i*sz:(i+1)*sz,:] = a
    return r


@jit(nopython=True, fastmath=False, nogil=True)
def mulWMJ(w, M, jac):
    """Compute the sparse matrix product P = WMJ, where
    
      W is a diagonal matrix (represented by a vector)
      M is a resolution matrix (the .data field of a
        ResMatrix, since Numba can't handle classes)
      J is a column-sparse matrix giving the nonzero
        values in one contiguous range per column
    
    We assume that J jac has enough free space in each column to store the
    results of the multiply and update it in place, returning the same jac as
    our return value.

    """    
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


@jit(nopython=True, fastmath=False, nogil=True)
def prepare_bins(centers, camerapix):
    """Convert N bin centers to the info needed by the optimizer, returned as an
    opaque tuple.  Given N bins, return values include
    
    - N+1 log bin edges.  Edges are placed halfway between centers,
    with extrapolation at the ends.  We return the natural log of
    each edge's wavelength, since that is what model and Jacobian
    building need.
    - N inverse bin widths, in an array padded by one cell on
      the left and right.

    """
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


def _drop_params(parameters, emfit, linemodel, Ifree, log):
    """Drop dubious free parameters after fitting.

    """
    # Conditions for dropping a parameter (all parameters, not just those
    # being fitted):
    # --negative amplitude or sigma
    # --parameter at its default value (fit failed, right??)
    # --parameter within 0.1% of its bounds
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
    notfixed = np.logical_not(linemodel['fixed'])

    # drop any negative amplitude or sigma parameter that is not fixed 
    drop1 = np.hstack((lineamps < 0, np.zeros(len(linevshifts), bool), linesigmas <= 0)) * notfixed
    
    # Require equality, not np.isclose, because the optimization can be very
    # small (<1e-6) but still significant, especially for the doublet
    # ratios. If linesigma is dropped this way, make sure the corresponding
    # line-amplitude is dropped, too (see MgII 2796 on
    # sv1-bright-17680-39627622543528153).
    drop2 = np.zeros(len(parameters), bool)
        
    # if any amplitude is zero, drop the corresponding sigma and vshift
    amp_param_bool = emfit.amp_param_bool[Ifree]
    I = np.where(parameters[Ifree][amp_param_bool] == 0.)[0]
    if len(I) > 0:
        _Ifree = np.zeros(len(parameters), bool)
        _Ifree[Ifree] = True
        for pp in linemodel[Ifree][amp_param_bool][I]['param_name']:
            J = np.where(_Ifree * (linemodel['param_name'] == pp.replace('_amp', '_sigma')))[0]
            drop2[J] = True
            K = np.where(_Ifree * (linemodel['param_name'] == pp.replace('_amp', '_vshift')))[0]
            drop2[K] = True
            #print(pp, J, K, np.sum(drop2))

    # drop amplitudes for any lines tied to a line with a dropped sigma
    sigmadropped = np.where(emfit.sigma_param_bool * drop2)[0]
    for lineindx, dropline in zip(sigmadropped, linemodel[sigmadropped]['linename']):
        # Check whether lines are tied to this line. If so, find the
        # corresponding amplitude and drop that, too.
        T = linemodel['tiedtoparam'] == lineindx
        for tiedline in set(linemodel['linename'][T]):
            drop2[linemodel['param_name'] == f'{tiedline}_amp'] = True
        drop2[linemodel['param_name'] == f'{dropline}_amp'] = True

    # drop amplitudes for any lines tied to a line with a dropped vshift
    vshiftdropped = np.where(emfit.vshift_param_bool * drop2)[0]
    for lineindx, dropline in zip(vshiftdropped, linemodel[vshiftdropped]['linename']):
        # Check whether lines are tied to this line. If so, find the
        # corresponding amplitude and drop that, too.
        T = linemodel['tiedtoparam'] == lineindx
        for tiedline in set(linemodel['linename'][T]):
            drop2[linemodel['param_name'] == f'{tiedline}_amp'] = True
        drop2[linemodel['param_name'] == f'{dropline}_amp'] = True

    # drop any non-fixed parameters outside their bounds
    # It's OK for parameters to be *at* their bounds.
    drop3 = np.zeros(len(parameters), bool)
    drop3[Ifree] = np.logical_or(parameters[Ifree] < linemodel['bounds'][Ifree, 0], 
                                 parameters[Ifree] > linemodel['bounds'][Ifree, 1])
    drop3 *= notfixed
    
    log.debug(f'Dropping {np.sum(drop1)} negative-amplitude lines.') # linewidth can't be negative
    log.debug(f'Dropping {np.sum(drop2)} sigma,vshift parameters of zero-amplitude lines.')
    log.debug(f'Dropping {np.sum(drop3)} parameters which are out-of-bounds.')
    Idrop = np.where(np.logical_or.reduce((drop1, drop2, drop3)))[0]
    
    if len(Idrop) > 0:
        log.debug(f'  Dropping {len(Idrop)} unique parameters.')
        parameters[Idrop] = 0.0


def emfit_bestfit(emfit, linemodel, redshift, emlinewave, resolution_matrix, camerapix, fast=False):
    """Construct the best-fitting emission-line spectrum from a linemodel."""
    
    parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, \
                 doubletpair, linewaves) = emfit._linemodel_to_parameters(linemodel, emfit.fit_linetable)
    
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
    
    # doublets
    lineamps[doubletindx] *= lineamps[doubletpair]

    if fast:
        from FasterSpecFit import EMLine_build_model
        emlinemodel = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas,
                                         linewaves, emlinewave, resolution_matrix, camerapix)

    else:
        from fastspecfit.emlines import build_emline_model
        emlinemodel = build_emline_model(emfit.dlog10wave, redshift, lineamps, linevshifts, linesigmas,
                                         linewaves, emlinewave, resolution_matrix, camerapix)

    return emlinemodel


def optimize(emfit, linemodel,
             obs_bin_centers,
             obs_bin_fluxes,
             obs_weights,
             redshift,
             resolution_matrices,
             camerapix,
             get_finalamp=False,
             log=None,
             verbose=False,
             debug=False):
    """Optimization routine.

    """
    from scipy.optimize import least_squares
    from fastspecfit.fitting import ParamsMapping
    
    if log is None:
        from desiutil.log import get_logger, DEBUG
        if verbose:
            log = get_logger(DEBUG)
        else:
            log = get_logger()
            
    parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, doubletpair, \
                 line_wavelengths) = emfit._linemodel_to_parameters(linemodel, emfit.fit_linetable)
    
    log.debug(f"Optimizing {len(Ifree)} free parameters")
    
    # corner case where all lines are out of the wavelength range, which can
    # happen at high redshift and with the red camera masked, e.g.,
    # iron/main/dark/6642/39633239580608311).
    initial_guesses = parameters[Ifree]

    params_mapping = ParamsMapping(parameters, Ifree,
                                   Itied, tiedtoparam, tiedfactor,
                                   doubletindx, doubletpair)
    
    bin_data = prepare_bins(obs_bin_centers, camerapix)
    
    farg = (bin_data, obs_bin_fluxes, obs_weights, redshift,
            line_wavelengths, tuple(resolution_matrices), camerapix,
            params_mapping)
    
    if len(Ifree) == 0:
        fit_info = {'nfev': 0, 'status': 0}
    else:
        try:
            fit_info = least_squares(objective, initial_guesses, jac=jacobian, args=farg, max_nfev=5000, 
                                     xtol=1e-10, ftol=1e-5, #x_scale='jac' gtol=1e-10,
                                     tr_solver='lsmr', tr_options={'maxiter': 1000, 'regularize': True},
                                     method='trf', bounds=tuple(zip(*bounds)),) # verbose=2)
            parameters[Ifree] = fit_info.x
        except:
            if emfit.uniqueid:
                errmsg = f'Problem in scipy.optimize.least_squares for {emfit.uniqueid}.'
            else:
                errmsg = 'Problem in scipy.optimize.least_squares.'
                log.critical(errmsg)
                raise RuntimeError(errmsg)

    # Drop (zero out) any dubious free parameters.
    _drop_params(parameters, emfit, linemodel, Ifree, log)
    
    # At this point, parameters contains correct *free* and *fixed* values, but
    # we need to update *tied* values to reflect any changes to free params.  We
    # do *not* apply doublet rules, as other code expects us to return a params
    # array with doublet ratios as ratios, not amplitudes.
    for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
        parameters[I] = parameters[indx] * factor
    
    out_linemodel = linemodel.copy()
    out_linemodel['value'] = parameters.copy() # so we don't munge it below
    out_linemodel.meta['nfev'] = fit_info['nfev']
    out_linemodel.meta['status'] = fit_info['status']
    
    if get_finalamp:
        
        lineamps, _, __ = np.array_split(parameters, 3) # 3 parameters per line

        # apply doublet rules
        lineamps[doubletindx] *= lineamps[doubletpair]
        
        # calculate the observed maximum amplitude for each
        # fitted spectral line after convolution with the resolution
        # matrix.
        peaks = find_peak_amplitudes(parameters,
                                     bin_data,
                                     redshift,
                                     line_wavelengths,
                                     resolution_matrices,
                                     camerapix)

        # FIXME: is len(lineamps) == # lines obtainable from line table?
        out_linemodel['obsvalue'][:len(lineamps)] = peaks
        
    return out_linemodel

#import numba
#@numba.jit(nopython=True)
def build_emline_model(dlog10wave, redshift, lineamps, linevshifts, linesigmas, 
                       linewaves, emlinewave, resolution_matrix, camerapix=None,
                       debug=False):

    """Given parameters, build the model emission-line spectrum.

    ToDo: can this be optimized using numba?

    """
    from fastspecfit.util import trapz_rebin

    #log10model = np.zeros_like(log10wave) # [erg/s/cm2/A, observed frame]
    log10wave = [] # initialize

    # Cut to lines with non-zero amplitudes.
    #I = linesigmas > 0
    I = lineamps > 0
    if np.count_nonzero(I) > 0:
        linevshifts = linevshifts[I]
        linesigmas = linesigmas[I]
        lineamps = lineamps[I]
        linewaves = linewaves[I]

        # demand at least 20 km/s for rendering the model
        if np.any(linesigmas) < 20.:
            linesigmas[linesigmas<20.] = 20.

        # line-width [log-10 Angstrom] and redshifted wavelength [log-10 Angstrom]
        log10sigmas = linesigmas / C_LIGHT / np.log(10) 
        linezwaves = np.log10(linewaves * (1.0 + redshift + linevshifts / C_LIGHT))

        # Build the wavelength vector on-the-fly.
        if camerapix is None:
            minwave = emlinewave[0][0]-2.
            maxwave = emlinewave[-1][-1]+2.
        else:
            minwave = emlinewave[0]-2.
            maxwave = emlinewave[-1]+2.
        
        log10wave = []
        for linezwave, log10sigma in zip(linezwaves, log10sigmas):
            log10wave.append(np.arange(linezwave - (5 * log10sigma), linezwave + (5 * log10sigma), dlog10wave))
        log10wave = np.hstack([np.log10(minwave), np.log10(maxwave), ] + log10wave)
        S = np.argsort(log10wave)
        log10wave = log10wave[S]
        log10model = np.zeros_like(log10wave)
        
        for lineamp, linezwave, log10sigma in zip(lineamps, linezwaves, log10sigmas):
            J = np.abs(log10wave - linezwave) < (5 * log10sigma)
            #print(lineamp, 10**linezwave, 10**log10wave[J].min(), 10**log10wave[J].max())
            log10model[J] += lineamp * np.exp(-0.5 * (log10wave[J]-linezwave)**2 / log10sigma**2)

    # Optionally split into cameras, resample, and convolve with the resolution
    # matrix.
    emlinemodel = []
    if camerapix is None:
        for icam, specwave in enumerate(emlinewave):
            if len(log10wave) > 0:
                _emlinemodel = trapz_rebin(10**log10wave, log10model, specwave)
                _emlinemodel = resolution_matrix[icam].dot(_emlinemodel)
                _emlinemodel[_emlinemodel < 0] = 0. # resolution matrix is occasionally negative
            else:
                _emlinemodel = specwave * 0
            emlinemodel.append(_emlinemodel)
        return emlinemodel
    else:
        for icam, campix in enumerate(camerapix):
            # can be zero-length if all lines are dropped
            if len(log10wave) > 0:
                _emlinemodel = trapz_rebin(10**log10wave, log10model, emlinewave[campix[0]:campix[1]])
                _emlinemodel = resolution_matrix[icam].dot(_emlinemodel)
                _emlinemodel[_emlinemodel < 0] = 0. # resolution matrix is occasionally negative
            else:
                _emlinemodel = emlinewave[campix[0]:campix[1]] * 0
            emlinemodel.append(_emlinemodel)
        return np.hstack(emlinemodel)

def _objective_function(free_parameters, emlinewave, emlineflux, weights, redshift, 
                        dlog10wave, resolution_matrix, camerapix, parameters, Ifree, 
                        Itied, tiedtoparam, tiedfactor, doubletindx, doubletpair, 
                        linewaves):
    """The parameters array should only contain free (not tied or fixed) parameters."""

    # Parameters have to be allowed to exceed their bounds in the optimization
    # function, otherwise they get stuck at the boundary.

    # Handle tied parameters and bounds. Only check bounds on the free
    # parameters.
    #for I, (value, bnd) in enumerate(zip(free_parameters, bounds)):
    #    if value < bnd[0]:
    #        free_parameters[I] = bnd[0]
    #    if value > bnd[1]:
    #        free_parameters[I] = bnd[1]

    #print(free_parameters)
    parameters[Ifree] = free_parameters

    if len(Itied) > 0:
        for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
            parameters[I] = parameters[indx] * factor

    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line

    # doublets
    lineamps[doubletindx] *= lineamps[doubletpair]

    # Build the emission-line model.
    emlinemodel = build_emline_model(dlog10wave, redshift, lineamps, linevshifts, 
                                     linesigmas, linewaves, emlinewave, 
                                     resolution_matrix, camerapix)

    if weights is None:
        residuals = emlinemodel - emlineflux
    else:
        residuals = weights * (emlinemodel - emlineflux)

    return residuals

class EMFitTools(Filters):
    def __init__(self, fphoto=None, emlinesfile=None, uniqueid=None,
                 minsnr_balmer_broad=3.):
        """Class to model a galaxy stellar continuum.

        Parameters
        ----------
        templates : :class:`str`, optional
            Full path to the templates used for continuum-fitting.
        mintemplatewave : :class:`float`, optional, defaults to None
            Minimum template wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxtemplatewave : :class:`float`, optional, defaults to 6e4
            Maximum template wavelength to read into memory. 
        chi2_default : :class:`float`, optional, defaults to 0.0.
            Default chi2 value for a emission line not fitted.
        maxiter : :class:`int`, optional, defaults to 5000.
            Maximum number of iterations.
        accuracy : :class:`float`, optional, defaults to 0.01.
            Fitting accuracy.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.

        Notes
        -----
        Need to document all the attributes.
        
        Plans for improvement:
          - Update the continuum redshift using cross-correlation.
          - Don't draw reddening from a flat distribution (try gamma or a custom
            distribution of the form x**2*np.exp(-2*x/scale).

        """
        super(EMFitTools, self).__init__(fphoto=fphoto)

        self.uniqueid = uniqueid

        self.linetable = read_emlines(emlinesfile=emlinesfile)

        self.emwave_pixkms = 5.                                   # pixel size for internal wavelength array [km/s]
        self.dlog10wave = self.emwave_pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]

        # default line-sigma for computing upper limits
        self.limitsigma_narrow = 75.0
        self.limitsigma_broad = 1200.0 
        self.wavepad = 2.5 # Angstrom

        # Establish the names of the parameters and doublets here, at
        # initialization, because we use them when instantiating the best-fit
        # model not just when fitting.
        doublet_names = ['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio']
        doublet_pairs = ['mgii_2803_amp', 'oii_3729_amp', 'sii_6716_amp']

        param_names = []
        for param in ['amp', 'vshift', 'sigma']:
            for linename in self.linetable['name'].data:
                param_name = linename+'_'+param
                # Use doublet-ratio parameters for close or physical
                # doublets. Note that any changes here need to be propagated to
                # the XX method, which needs to "know" about these doublets.
                if param_name == 'mgii_2796_amp':
                    param_name = 'mgii_doublet_ratio' # MgII 2796/2803
                if param_name == 'oii_3726_amp':
                    param_name = 'oii_doublet_ratio'  # [OII] 3726/3729
                if param_name == 'sii_6731_amp':
                    param_name = 'sii_doublet_ratio'  # [SII] 6731/6716
                param_names.append(param_name)
        self.param_names = np.hstack(param_names)
        self.amp_param_bool = np.array(['_amp' in pp for pp in self.param_names])
        self.amp_balmer_bool = np.array(['_amp' in pp and 'hei_' not in pp and 'heii_' not in pp for pp in self.param_names]) # no helium lines
        self.sigma_param_bool = np.array(['_sigma' in pp for pp in self.param_names])
        self.vshift_param_bool = np.array(['_vshift' in pp for pp in self.param_names])

        self.doubletindx = np.hstack([np.where(self.param_names == doublet)[0] for doublet in doublet_names])
        self.doubletpair = np.hstack([np.where(self.param_names == pair)[0] for pair in doublet_pairs])

        self.minsigma_balmer_broad = 250. # minimum broad-line sigma [km/s]
        self.minsnr_balmer_broad = minsnr_balmer_broad # minimum broad-line S/N

    def summarize_linemodel(self, linemodel):
        """Simple function to summarize an input linemodel."""
        def _print(linenames):
            for linename in linenames:
                for param in ['amp', 'sigma', 'vshift']:
                    I = np.where(self.param_names == linename+'_'+param)[0]
                    if len(I) == 1:
                        I = I[0]
                        if linemodel['tiedtoparam'][I] == -1:
                            if linemodel['fixed'][I]:
                                print('{:25s} is NOT FITTED'.format(linename+'_'+param))
                            else:
                                print('{:25s} untied'.format(linename+'_'+param))
                        else:
                            if linemodel['fixed'][I]:
                                print('{:25s} tied to {:25s} with factor {:.4f} and FIXED'.format(
                                    linename+'_'+param, self.param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))
                            else:
                                print('{:25s} tied to {:25s} with factor {:.4f}'.format(
                                    linename+'_'+param, self.param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))
                                
        linenames = self.fit_linetable['name'].data

        print('---------------------')
        print('UV/QSO (broad) lines:')
        print('---------------------')
        _print(linenames[(self.fit_linetable['isbroad'] == True) * (self.fit_linetable['isbalmer'] == False)])
        print()
        print('--------------------------')
        print('Broad Balmer+helium lines:')
        print('--------------------------')
        _print(linenames[(self.fit_linetable['isbroad'] == True) * (self.fit_linetable['isbalmer'] == True)])
        print()
        print('---------------------------')
        print('Narrow Balmer+helium lines:')
        print('---------------------------')
        _print(linenames[(self.fit_linetable['isbroad'] == False) * (self.fit_linetable['isbalmer'] == True)])
        print()
        print('----------------')
        print('Forbidden lines:')
        print('----------------')
        _print(linenames[(self.fit_linetable['isbroad'] == False) * (self.fit_linetable['isbalmer'] == False)])

    def build_linemodels(self, redshift, wavelims=[3000, 10000], verbose=False, strict_broadmodel=True):
        """Build all the multi-parameter emission-line models we will use.
    
        """
        def _fix_parameters(linemodel, verbose=False):
            """Set the "fixed" attribute for all the parameters in a given linemodel."""
            # First loop through all tied parameters and set fixed to the
            # parameter it's tied to.
            I = np.where(linemodel['tiedtoparam'] != -1)[0] # should always have len(I)>0
            alltied = linemodel[I]['tiedtoparam']
            utied = np.unique(alltied)
            for tied in utied:
                J = tied == alltied
                if verbose:
                    print('Tying {} to {}'.format(' '.join(linemodel[I][J]['param_name']), linemodel[tied]['param_name']))
                linemodel[I][J]['fixed'] = linemodel[tied]['fixed']
            if verbose:
                print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

            # Next, fix out-of-range lines but not those that are in the 'utied'
            # array---those out-of-range lines need to be in the optimization
            # list because the in-range lines depend on them.
            outofrange = fit_linetable['inrange'] == False
            if np.sum(outofrange) > 0: # should always be true
                for linename in fit_linetable['name'][outofrange]:
                    for param in ['amp', 'vshift', 'sigma']:
                        param_name = linename+'_'+param
                        I = np.where(linemodel['param_name'] == param_name)[0]
                        if len(I) > 0:
                            if I in utied:
                                if verbose:
                                    print('Not fixing out-of-range parameter {}'.format(param_name))
                            else:
                                linemodel['fixed'][I] |= True
                if verbose:
                    print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                    print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                    #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

                # Finally loop through each 'utied' line and if all the lines
                # tied to it are fixed, then fix that line, too.
                for tied in utied:
                    if linemodel['param_name'][tied] and np.all(linemodel[linemodel['tiedtoparam'] == tied]['fixed']):
                        if outofrange[linemodel['linename'][tied] == fit_linetable['name']]:
                            if verbose:
                                print('Fixing {} because line is out of range and all tied lines are fixed: {}'.format(
                                    linemodel['param_name'][tied], ' '.join(linemodel[linemodel['tiedtoparam'] == tied]['param_name'])))
                            linemodel[tied]['fixed'] = True

                # Also handle the doublets.
                I = np.where(linemodel['doubletpair'] != -1)[0]
                if len(I) > 0:
                    for doublet in linemodel[I]['doubletpair']:
                        J = doublet == linemodel['doubletpair']
                        if linemodel[doublet]['fixed']:
                            linemodel['fixed'][J] = True

                if verbose:
                    print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                    print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                    #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

        initvshift = 1.0
        vmaxshift_narrow = 500.0
        vmaxshift_broad = 2500.0 # 3000.0
        vmaxshift_balmer_broad = 2500.
    
        minsigma_narrow = 1.0
        maxsigma_narrow = 750.0 # 500.0

        minsigma_broad = 1.0 # 100.
        maxsigma_broad = 1e4

        minsigma_balmer_broad = 1.0 # 100.0 # minsigma_narrow
        maxsigma_balmer_broad = maxsigma_broad
    
        # Be very careful about changing the default broad line-sigma. Smaller
        # values like 1500 km/s (which is arguably more sensible) can lead to
        # low-amplitude broad lines in a bunch of normal star-forming galaxy
        # spectra. (They act to "suck up" local continuum variations.) Also
        # recall that if it's well-measured, we use the initial line-sigma in
        # estimate_linesigma, which is a better initial guess.
        initsigma_narrow = 75.0 # 260.0 # 75.0
        initsigma_broad = 3000.0  
    
        initamp = 0.0
        #minamp = 0.0
        minamp = -1e2
        maxamp = +1e5
        minamp_balmer_broad = minamp # 0.0
        maxamp_balmer_broad = maxamp
    
        # Specialized parameters on the MgII, [OII], and [SII] doublet ratios. See
        # https://github.com/desihub/fastspecfit/issues/39. Be sure to set
        # self.doublet_names, below, and also note that any change in the order of
        # these lines has to be handled in _emline_spectrum!
        init_mgii_doublet = 0.5 # MgII 2796/2803
        init_oii_doublet = 0.74 # [OII] 3726/3729
        init_sii_doublet = 0.74 # [SII] 6731/6716

        bounds_mgii_doublet = [0.0, 10.0] 
        bounds_oii_doublet = [0.0, 2.0] # [0.5, 1.5] # [0.66, 1.4]
        bounds_sii_doublet = [0.0, 2.0] # [0.5, 1.5] # [0.67, 1.2]
    
        # Create a new line-fitting table which contains the redshift-dependent
        # quantities for this object.
        fit_linetable = Table()
        fit_linetable['name'] = self.linetable['name']
        fit_linetable['isbalmer'] = self.linetable['isbalmer']
        fit_linetable['ishelium'] = self.linetable['ishelium']
        fit_linetable['isbroad'] = self.linetable['isbroad']
        fit_linetable['restwave'] = self.linetable['restwave']
        fit_linetable['zwave'] = self.linetable['restwave'].data * (1 + redshift)
        fit_linetable['inrange'] = ((fit_linetable['zwave'] > (wavelims[0]+self.wavepad)) * 
                                    (fit_linetable['zwave'] < (wavelims[1]-self.wavepad)))
        self.fit_linetable = fit_linetable
        
        linenames = fit_linetable['name'].data
        param_names = self.param_names
        nparam = len(param_names)

        # Model 1 -- here, parameters are minimally tied together for the final
        # fit and only lines outside the wavelength range are fixed. Includes
        # broad lines.
        linemodel_broad = Table()
        linemodel_broad['param_name'] = param_names
        linemodel_broad['index'] = np.arange(nparam).astype(np.int32)
        linemodel_broad['linename'] = np.tile(linenames, 3) # 3 parameters per line
        linemodel_broad['isbalmer'] = np.zeros(nparam, bool)
        linemodel_broad['ishelium'] = np.zeros(nparam, bool)
        linemodel_broad['isbroad'] = np.zeros(nparam, bool)
        linemodel_broad['tiedfactor'] = np.zeros(nparam, 'f8')
        linemodel_broad['tiedtoparam'] = np.zeros(nparam, np.int16)-1
        linemodel_broad['doubletpair'] = np.zeros(nparam, np.int16)-1
        linemodel_broad['fixed'] = np.zeros(nparam, bool)
        linemodel_broad['bounds'] = np.zeros((nparam, 2), 'f8')
        linemodel_broad['initial'] = np.zeros(nparam, 'f8')
        linemodel_broad['value'] = np.zeros(nparam, 'f8')
        linemodel_broad['obsvalue'] = np.zeros(nparam, 'f8')
        linemodel_broad['civar'] = np.zeros(nparam, 'f8') # continuum inverse variance

        linemodel_broad['doubletpair'][self.doubletindx] = self.doubletpair

        # Build the relationship of "tied" parameters. In the 'tied' array, the
        # non-zero value is the multiplicative factor by which the parameter
        # represented in the 'tiedtoparam' index should be multiplied.
    
        # Physical doublets and lines in the same ionization species should have
        # their velocity shifts and line-widths always tied. In addition, set fixed
        # doublet-ratios here. Note that these constraints must be set on *all*
        # lines, not just those in range.
    
        for iline, linename in enumerate(linenames):
            linemodel_broad['isbalmer'][linemodel_broad['linename'] == linename] = fit_linetable[fit_linetable['name'] == linename]['isbalmer']
            linemodel_broad['ishelium'][linemodel_broad['linename'] == linename] = fit_linetable[fit_linetable['name'] == linename]['ishelium']
            linemodel_broad['isbroad'][linemodel_broad['linename'] == linename] = fit_linetable[fit_linetable['name'] == linename]['isbroad']
            
            # initial values and bounds - broad He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp_balmer_broad, maxamp_balmer_broad], 
                                                   [minsigma_balmer_broad, maxsigma_balmer_broad],
                                                   [-vmaxshift_balmer_broad, +vmaxshift_balmer_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - narrow He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - broad UV/QSO lines (non-Balmer)
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_broad, maxsigma_broad],
                                                   [-vmaxshift_broad, +vmaxshift_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - forbidden lines
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # tie parameters

            # broad He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]
            #print('Releasing the narrow Balmer lines!')
            # narrow He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False and linename != 'halpha':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
            # other lines
            if linename == 'mgii_2796':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'mgii_2803_'+param)[0]
            if linename == 'nev_3346' or linename == 'nev_3426': # should [NeIII] 3869 be tied to [NeV]???
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'neiii_3869_'+param)[0]
            if linename == 'oii_3726':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_3729_'+param)[0]
            # Tentative! Tie auroral lines to [OIII] 4363 but maybe we shouldn't tie [OI] 6300 here...
            if linename == 'nii_5755' or linename == 'oi_6300' or linename == 'siii_6312':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_4363_'+param)[0]
            if linename == 'oiii_4959':
                """
                [O3] (4-->2): airwave: 4958.9097 vacwave: 4960.2937 emissivity: 1.172e-21
                [O3] (4-->3): airwave: 5006.8417 vacwave: 5008.2383 emissivity: 3.497e-21
                """
                linemodel_broad['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9839 # 2.8875
                linemodel_broad['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oiii_5007_amp')[0]
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
            if linename == 'nii_6548':
                """
                [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.02198e-21
                [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.94901e-21
                """
                linemodel_broad['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9421 # 2.936
                linemodel_broad['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'nii_6584_amp')[0]
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nii_6584_'+param)[0]
            if linename == 'sii_6731':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'sii_6716_'+param)[0]
            if linename == 'oii_7330':
                """
                [O2] (5-->2): airwave: 7318.9185 vacwave: 7320.9350 emissivity: 8.18137e-24
                [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 2.40519e-23
                [O2] (5-->3): airwave: 7329.6613 vacwave: 7331.6807 emissivity: 1.35614e-23
                [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.27488e-23
                """
                linemodel_broad['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 1.2251
                linemodel_broad['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oii_7320_amp')[0]
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_7320_'+param)[0]
            if linename == 'siii_9069':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'siii_9532_'+param)[0]
            # Tentative! Tie SiIII] 1892 to CIII] 1908 because they're so close in wavelength.
            if linename == 'siliii_1892':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'ciii_1908_'+param)[0]

            # Tie all the forbidden and narrow Balmer+helium lines *except
            # [OIII] 4959,5007* to [NII] 6584 when we have broad lines. The
            # [OIII] doublet frequently has an outflow component, so fit it
            # separately. See the discussion at
            # https://github.com/desihub/fastspecfit/issues/160
            if strict_broadmodel:
                if fit_linetable['isbroad'][iline] == False and linename != 'nii_6584' and linename != 'oiii_4959' and linename != 'oiii_5007':
                    for param in ['sigma', 'vshift']:
                        linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                        linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nii_6584_'+param)[0]
                        
                #if fit_linetable['isbroad'][iline] == False and linename != 'oiii_5007':
                #    for param in ['sigma', 'vshift']:
                #        linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                #        linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
                
                ## Tie all forbidden lines to [OIII] 5007; the narrow Balmer and
                ## helium lines are separately tied together.
                #if fit_linetable['isbroad'][iline] == False and fit_linetable['isbalmer'][iline] == False and linename != 'oiii_5007':
                #    for param in ['sigma']:
                #        linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                #        linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
                    
        # Finally set the initial values and bounds on the doublet ratio parameters.
        for param, bounds, default in zip(['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio'],
                                          [bounds_mgii_doublet, bounds_oii_doublet, bounds_sii_doublet],
                                          [init_mgii_doublet, init_oii_doublet, init_sii_doublet]):
            linemodel_broad['initial'][linemodel_broad['param_name'] == param] = default
            linemodel_broad['bounds'][linemodel_broad['param_name'] == param] = bounds
                    
        # Assign fixed=True to parameters which are outside the wavelength range
        # except those that are tied to other lines.
        _fix_parameters(linemodel_broad, verbose=False)

        assert(np.all(linemodel_broad['tiedtoparam'][linemodel_broad['tiedfactor'] != 0] != -1))
        # It's OK for the doublet ratios to be bounded at zero.
        #assert(len(linemodel_broad[np.sum(linemodel_broad['bounds'] == [0.0, 0.0], axis=1) > 0]) == 0)
    
        #_print_linemodel(linemodel_broad)
        #linemodel_broad[np.logical_and(linemodel_broad['fixed'] == False, linemodel_broad['tiedtoparam'] == -1)]

        # Model 2 - like linemodel, but broad lines have been fixed at zero.
        linemodel_nobroad = linemodel_broad.copy()
        linemodel_nobroad['fixed'] = False # reset

        for iline, linename in enumerate(linenames):
            if linename == 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    linemodel_nobroad['fixed'][param_names == linename+'_'+param] = True

            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

            if strict_broadmodel:
                # Tie the forbidden lines to [OIII] 5007.
                if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline] == False and linename != 'oiii_5007':
                    for param in ['sigma', 'vshift']:
                        linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                        linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
                        
                # Tie narrow Balmer and helium lines together.
                if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False:
                    if linename == 'halpha':
                        for param in ['sigma', 'vshift']:
                            linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 0.0
                            linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = -1
                    else:
                        for param in ['sigma', 'vshift']:
                            linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                            linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
                
        #linemodel_nobroad[np.logical_and(linemodel_nobroad['fixed'] == False, linemodel_nobroad['tiedtoparam'] == -1)]

        _fix_parameters(linemodel_nobroad, verbose=False)
        assert(np.all(linemodel_nobroad['tiedtoparam'][linemodel_nobroad['tiedfactor'] != 0] != -1))

        return linemodel_broad, linemodel_nobroad

    def initial_guesses_and_bounds(self, data, emlinewave, emlineflux, log):
        """For all lines in the wavelength range of the data, get a good initial guess
        on the amplitudes and line-widths. This step is critical for cases like,
        e.g., 39633354915582193 (tile 80613, petal 05), which has strong narrow
        lines.

        """
        initial_guesses, param_bounds = {}, {}
        #init_amplitudes, init_sigmas = {}, {}

        coadd_sigma = data['smoothsigma'] # robust estimate of the variance in the spectrum
        coadd_emlineflux = np.interp(data['coadd_wave'], emlinewave, emlineflux)
    
        for linename, linepix, contpix in zip(data['coadd_linename'], data['coadd_linepix'], data['coadd_contpix']):
            ## skip the physical doublets
            #if not hasattr(self.EMLineModel, '{}_amp'.format(linename)):
            #    continue

            npix = len(linepix)
            if npix > 5:
                mnpx, mxpx = linepix[npix//2]-3, linepix[npix//2]+3
                if mnpx < 0:
                    mnpx = 0
                if mxpx > linepix[-1]:
                    mxpx = linepix[-1]
                amp = np.max(coadd_emlineflux[mnpx:mxpx])
            else:
                amp = np.percentile(coadd_emlineflux[linepix], 97.5)
            if amp < 0:
                amp = np.abs(amp)

            noise = np.mean(coadd_sigma[linepix])
            if noise == 0.:
                civar = 0.
                errmsg = 'Noise estimate for line {} is zero!'.format(linename)
                log.warning(errmsg)
                #raise ValueError(errmsg)
            else:
                civar = 1. / noise**2

            # update the bounds on the line-amplitude
            #bounds = [-np.min(np.abs(coadd_emlineflux[linepix])), 3*np.max(coadd_emlineflux[linepix])]
            mx = 5*np.max(coadd_emlineflux[linepix])
            if mx < 0: # ???
                mx = 5*np.max(np.abs(coadd_emlineflux[linepix]))
            
            iline = self.linetable[self.linetable['name'] == linename]
            bounds = [-1.5*np.min(np.abs(coadd_emlineflux[linepix])), mx]

            # In extremely rare cases many of the pixels are zero, in which case
            # bounds[0] becomes zero, which is bad (e.g.,
            # iron/main/dark/27054/39627811564029314). Fix that here.
            if np.abs(bounds[0]) == 0.0:
                N = coadd_emlineflux[linepix] != 0
                if np.sum(N) > 0:
                    bounds[0] = -1.5*np.min(np.abs(coadd_emlineflux[linepix][N]))
                if np.abs(bounds[0]) == 0.0:
                    bounds[0] = -1e-3 # ??
            
            if (bounds[0] > bounds[1]) or (amp < bounds[0]) or (amp > bounds[1]):
                log.warning('Initial amplitude is outside its bound for line {}.'.format(linename))
                amp = np.diff(bounds)/2 + bounds[0]
                # Should never happen.
                if (bounds[0] > bounds[1]) or (amp < bounds[0]) or (amp > bounds[1]):
                    errmsg = 'Initial amplitude is outside its bound for line {}.'.format(linename)
                    self.log.critical(errmsg)
                    raise ValueError(errmsg)

            initial_guesses[linename+'_amp'] = amp
            param_bounds[linename+'_amp'] = bounds
            initial_guesses[linename+'_civar'] = civar

        # Now update the linewidth but here we need to loop over *all* lines
        # (not just those in range). E.g., if H-alpha is out of range we need to
        # set its initial value correctly since other lines are tied to it
        # (e.g., main-bright-32406-39628257196245904).
        for iline in self.linetable:
            linename = iline['name']
            if iline['isbroad']:
                if iline['isbalmer']: # broad Balmer lines
                    if data['linesigma_balmer'] > data['linesigma_narrow']:
                        initial_guesses[linename+'_sigma'] = data['linesigma_balmer']
                    else:
                        initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
                else: # broad UV/QSO lines
                    initial_guesses[linename+'_sigma'] = data['linesigma_uv']
            else:
                # prefer narrow over Balmer
                initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
    
        return initial_guesses, param_bounds

    @staticmethod
    def _linemodel_to_parameters(linemodel, fit_linetable):
        """Convert a linemodel model to a list of emission-line parameters."""

        linesplit = np.array_split(linemodel['index'], 3) # 3 parameters per line
        #linesplit = (np.arange(3) + 1) * len(linemodel) // 3 # 3 parameters per line
        lineamps = linemodel['value'][linesplit[0]].data
        linevshifts = linemodel['value'][linesplit[1]].data
        linesigmas = linemodel['value'][linesplit[2]].data

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        #doublet = np.where(linemodel['doubletpair'] != -1)[0]
        #lineamps[doublet] *= linemodel['value'][linemodel['doubletpair'][doublet]]
        parameters = np.hstack((lineamps, linevshifts, linesigmas))

        linewaves = fit_linetable['restwave'].data
        #lineinrange = fit_linetable['inrange'].data

        #Itied = np.where((linemodel['tiedtoparam'] != -1))[0]
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        Ifree = np.where((linemodel['tiedtoparam'] == -1) * (linemodel['fixed'] == False))[0]

        tiedtoparam = linemodel['tiedtoparam'][Itied].data
        tiedfactor = linemodel['tiedfactor'][Itied].data
        bounds = linemodel['bounds'][Ifree].data

        doubletindx = np.where(linemodel['doubletpair'] != -1)[0]
        doubletpair = linemodel['doubletpair'][doubletindx].data

        parameter_extras = (Ifree, Itied, tiedtoparam, tiedfactor, bounds,
                            doubletindx, doubletpair, linewaves)

        return parameters, parameter_extras

    @staticmethod
    def populate_linemodel(linemodel, initial_guesses, param_bounds, log):
        """Populate an input linemodel with initial guesses and parameter bounds, taking
        into account fixed parameters.

        """
        # Set initial values and bounds.
        for iparam, param in enumerate(linemodel['param_name']):
            if param in initial_guesses.keys():
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0 # always set fixed parameter to zero
                else:
                    # Make sure the initial guess for the narrow Balmer+helium
                    # line is smaller than the guess for the broad-line model.
                    if linemodel[iparam]['isbalmer'] and 'sigma' in param:
                        if linemodel[iparam]['isbroad']:
                            linemodel['initial'][iparam] = 1.1 * initial_guesses[param]
                        else:
                            linemodel['initial'][iparam] = initial_guesses[param]
                    else:
                        linemodel['initial'][iparam] = initial_guesses[param]
                    if param in param_bounds.keys():
                        linemodel['bounds'][iparam] = param_bounds[param]
                        # set the lower boundary on broad lines to be XX times the local noise
                        if linemodel['isbalmer'][iparam] and linemodel['isbroad'][iparam]:
                            civarkey = linemodel['linename'][iparam]+'_civar'
                            linemodel['civar'][iparam] = initial_guesses[civarkey]

                            #broadbalmer_snrmin = 3.
                            #linemodel['bounds'][iparam][0] = broadbalmer_snrmin * initial_guesses[noisekey]
                            #if linemodel['initial'][iparam] < linemodel['bounds'][iparam][0]:
                            #    linemodel['initial'][iparam] = linemodel['bounds'][iparam][0]
                            #if linemodel['initial'][iparam] > linemodel['bounds'][iparam][1]:
                            #    linemodel['initial'][iparam] = linemodel['bounds'][iparam][1]
            else:
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0
                else:
                    linemodel['initial'][iparam] = 1.0

            # Check bounds for free parameters but do not crash.
            if linemodel['fixed'][iparam] == False and linemodel['tiedtoparam'][iparam] == -1:
                toosml = linemodel['initial'][iparam] < linemodel['bounds'][iparam, 0]
                toobig = linemodel['initial'][iparam] > linemodel['bounds'][iparam, 1]
                if toosml:
                    errmsg = 'Initial parameter {} is outside its bound, {:.2f} < {:.2f}.'.format(
                        param, linemodel['initial'][iparam], linemodel['bounds'][iparam, 0])
                    log.warning(errmsg)
                    #raise ValueError(errmsg)
                    linemodel['initial'][iparam] = linemodel['bounds'][iparam, 0]
                if toobig:
                    errmsg = 'Initial parameter {} is outside its bound, {:.2f} > {:.2f}.'.format(
                        param, linemodel['initial'][iparam], linemodel['bounds'][iparam, 1])
                    log.warning(errmsg)
                    #raise ValueError(errmsg)
                    linemodel['initial'][iparam] = linemodel['bounds'][iparam, 1]

        # Now loop back through and ensure that tied relationships are enforced.
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        if len(Itied) > 0:
            for iparam, param in enumerate(linemodel['param_name'][Itied]):
                tieindx = linemodel['tiedtoparam'][Itied[iparam]]
                tiefactor = linemodel['tiedfactor'][Itied[iparam]]
                #log.info('{} tied to {} with factor {:.4f}'.format(param, linemodel[tieindx]['param_name'], tiefactor))
                linemodel['initial'][Itied[iparam]] = linemodel[tieindx]['initial'] * tiefactor

        linemodel['value'] = linemodel['initial'] # copy

    def optimize(self, linemodel, emlinewave, emlineflux, weights, redshift,
                 resolution_matrix, camerapix, log=None, get_finalamp=False,
                 verbose=False, debug=False):
        """Wrapper to call the least-squares minimization given a linemodel.

        """
        from scipy.optimize import least_squares
        from fastspecfit.util import trapz_rebin        

        if log is None:
            from desiutil.log import get_logger, DEBUG
            if verbose:
                log = get_logger(DEBUG)
            else:
                log = get_logger()

        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, doubletpair, \
                     linewaves) = self._linemodel_to_parameters(linemodel, self.fit_linetable)
        log.debug('Optimizing {} free parameters'.format(len(Ifree)))

        farg = (emlinewave, emlineflux, weights, redshift, self.dlog10wave, 
                resolution_matrix, camerapix, parameters, ) + \
                (Ifree, Itied, tiedtoparam, tiedfactor, doubletindx, 
                 doubletpair, linewaves)

        # corner case where all lines are out of the wavelength range, which can
        # happen at high redshift and with the red camera masked, e.g.,
        # iron/main/dark/6642/39633239580608311).
        initial_guesses = parameters[Ifree]

        if len(Ifree) == 0:
            fit_info = {'nfev': 0, 'status': 0}
        else:
            try:
                fit_info = least_squares(_objective_function, initial_guesses, args=farg, max_nfev=5000, 
                                         xtol=1e-10, #x_scale='jac', #ftol=1e-10, gtol=1e-10,
                                         tr_solver='lsmr', tr_options={'regularize': True},
                                         method='trf', bounds=tuple(zip(*bounds)))#, verbose=2)
                parameters[Ifree] = fit_info.x
            except:
                if self.uniqueid:
                    errmsg = 'Problem in scipy.optimize.least_squares for {}.'.format(self.uniqueid)
                else:
                    errmsg = 'Problem in scipy.optimize.least_squares.'
                log.critical(errmsg)
                raise RuntimeError(errmsg)

            # If the narrow-line sigma didn't change by more than ~one km/s from
            # its initial guess, then something has gone awry, so perturb the
            # initial guess by 10% and try again. Examples where this matters:
            #   fuji-sv3-bright-28119-39628390055022140
            #   fuji-sv3-dark-25960-1092092734472204
            S = np.where(self.sigma_param_bool[Ifree] * (linemodel['isbroad'][Ifree] == False))[0]
            if len(S) > 0:
                sig_init = initial_guesses[S]
                sig_final = parameters[Ifree][S]
                G = np.abs(sig_init - sig_final) < 1.
                if np.any(G):
                    log.warning(f'Poor convergence on line-sigma for {self.uniqueid}; perturbing initial guess and refitting.')
                    initial_guesses[S[G]] *= 0.9
                    try:
                        fit_info = least_squares(_objective_function, initial_guesses, args=farg, max_nfev=5000, 
                                                 xtol=1e-10, #x_scale='jac', #ftol=1e-10, gtol=1e-10,
                                                 tr_solver='lsmr', tr_options={'regularize': True},
                                                 method='trf', bounds=tuple(zip(*bounds)))#, verbose=2)
                        parameters[Ifree] = fit_info.x
                    except:
                        if self.uniqueid:
                            errmsg = f'Problem in scipy.optimize.least_squares for {self.uniqueid}.'
                        else:
                            errmsg = 'Problem in scipy.optimize.least_squares.'
                        log.critical(errmsg)
                        raise RuntimeError(errmsg)

        # Conditions for dropping a parameter (all parameters, not just those
        # being fitted):
        # --negative amplitude or sigma
        # --parameter at its default value (fit failed, right??)
        # --parameter within 0.1% of its bounds
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
        notfixed = np.logical_not(linemodel['fixed'])

        drop1 = np.hstack((lineamps < 0, np.zeros(len(linevshifts), bool), linesigmas <= 0)) * notfixed

        # Require equality, not np.isclose, because the optimization can be very
        # small (<1e-6) but still significant, especially for the doublet
        # ratios. If linesigma is dropped this way, make sure the corresponding
        # line-amplitude is dropped, too (see MgII 2796 on
        # sv1-bright-17680-39627622543528153).
        drop2 = np.zeros(len(parameters), bool)

        amp_param_bool = self.amp_param_bool[Ifree]
        I = np.where(parameters[Ifree][amp_param_bool] == 0)[0]
        if len(I) > 0:
            _Ifree = np.zeros(len(parameters), bool)
            _Ifree[Ifree] = True
            for pp in linemodel[Ifree][amp_param_bool][I]['param_name']:
                J = np.where(_Ifree * (linemodel['param_name'] == pp.replace('_amp', '_sigma')))[0]
                if len(J) > 0:
                    drop2[J] = True
                K = np.where(_Ifree * (linemodel['param_name'] == pp.replace('_amp', '_vshift')))[0]
                if len(K) > 0:
                    drop2[K] = True
                #print(pp, J, K, np.sum(drop2))
        
        #drop2[self.amp_param_bool] = parameters[self.amp_param_bool] == 0.0
        #drop2[Ifree] = parameters[Ifree] == linemodel['value'][Ifree]
        #drop2 *= notfixed
        
        sigmadropped = np.where(self.sigma_param_bool * drop2)[0]
        if len(sigmadropped) > 0:
            for lineindx, dropline in zip(sigmadropped, linemodel[sigmadropped]['linename']):
                # Check whether lines are tied to this line. If so, find the
                # corresponding amplitude and drop that, too.
                T = linemodel['tiedtoparam'] == lineindx
                if np.any(T):
                    for tiedline in set(linemodel['linename'][T]):
                        drop2[linemodel['param_name'] == '{}_amp'.format(tiedline)] = True
                drop2[linemodel['param_name'] == '{}_amp'.format(dropline)] = True

        vshiftdropped = np.where(self.vshift_param_bool * drop2)[0]
        if len(vshiftdropped) > 0:
            for lineindx, dropline in zip(vshiftdropped, linemodel[vshiftdropped]['linename']):
                # Check whether lines are tied to this line. If so, find the
                # corresponding amplitude and drop that, too.
                T = linemodel['tiedtoparam'] == lineindx
                if np.any(T):
                    for tiedline in set(linemodel['linename'][T]):
                        drop2[linemodel['param_name'] == '{}_amp'.format(tiedline)] = True
                drop2[linemodel['param_name'] == '{}_amp'.format(dropline)] = True

        # It's OK for parameters to be *at* their bounds.
        drop3 = np.zeros(len(parameters), bool)
        drop3[Ifree] = np.logical_or(parameters[Ifree] < linemodel['bounds'][Ifree, 0], 
                                     parameters[Ifree] > linemodel['bounds'][Ifree, 1])
        drop3 *= notfixed
        
        log.debug('Dropping {} negative-amplitude lines.'.format(np.sum(drop1))) # linewidth can't be negative
        log.debug('Dropping {} sigma,vshift parameters of zero-amplitude lines.'.format(np.sum(drop2)))
        log.debug('Dropping {} parameters which are out-of-bounds.'.format(np.sum(drop3)))
        Idrop = np.where(np.logical_or.reduce((drop1, drop2, drop3)))[0]
        
        if len(Idrop) > 0:
            log.debug('  Dropping {} unique parameters.'.format(len(Idrop)))
            parameters[Idrop] = 0.0
            
        # apply tied constraints
        if len(Itied) > 0:
            for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
                parameters[I] = parameters[indx] * factor

        # Now loop back through and drop Broad balmer lines that:
        #   (1) are narrower than their narrow-line counterparts;
        #   (2) have a narrow line whose amplitude is smaller than that of the broad line
        #      --> Deprecated! main-dark-32303-39628176678192981 is an example
        #          of an object where there's a broad H-alpha line but no other
        #          forbidden lines!

        out_linemodel = linemodel.copy()
        out_linemodel['value'] = parameters
        out_linemodel.meta['nfev'] = fit_info['nfev']
        out_linemodel.meta['status'] = fit_info['status']

        # Get the final line-amplitudes, after resampling and convolution (see
        # https://github.com/desihub/fastspecfit/issues/139). Some repeated code
        # from build_emline_model...
        if get_finalamp:
            lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
            lineindxs = np.arange(len(lineamps))
            
            I = lineamps > 0
            if np.count_nonzero(I) > 0:
                linevshifts = linevshifts[I]
                linesigmas = linesigmas[I]
                lineamps = lineamps[I]
                linewaves = linewaves[I]
                lineindxs = lineindxs[I]
    
                # demand at least 20 km/s for rendering the model
                if np.any(linesigmas) < 20.:
                    linesigmas[linesigmas<20.] = 20.

                if camerapix is None:
                    minwave = emlinewave[0][0]-2.
                    maxwave = emlinewave[-1][-1]+2.
                else:
                    minwave = emlinewave[0]-2.
                    maxwave = emlinewave[-1]+2.
                    
                _emlinewave = []
                for icam, campix in enumerate(camerapix):
                    _emlinewave.append(emlinewave[campix[0]:campix[1]])
        
                # line-width [log-10 Angstrom] and redshifted wavelength [log-10 Angstrom]
                log10sigmas = linesigmas / C_LIGHT / np.log(10)                
                linezwaves = np.log10(linewaves * (1.0 + redshift + linevshifts / C_LIGHT))
        
                for lineindx, lineamp, linezwave, log10sigma in zip(lineindxs, lineamps, linezwaves, log10sigmas):
                    log10wave = np.arange(linezwave - (5 * log10sigma), linezwave + (5 * log10sigma), self.dlog10wave)
                    log10wave = np.hstack((np.log10(minwave), log10wave, np.log10(maxwave)))
                    log10model = lineamp * np.exp(-0.5 * (log10wave-linezwave)**2 / log10sigma**2)
                    # Determine which camera we're on and then resample and
                    # convolve with the resolution matrix.
                    icam = np.argmin([np.abs((np.max(emwave)-np.min(emwave))/2+np.min(emwave)-10**linezwave) for emwave in _emlinewave])
                    model_resamp = trapz_rebin(10**log10wave, log10model, _emlinewave[icam])
                    model_convol = resolution_matrix[icam].dot(model_resamp)
                    out_linemodel['obsvalue'][lineindx] = np.max(model_convol)
                    #if out_linemodel[lineindx]['param_name'] == 'oiii_5007_amp':
                    #    import matplotlib.pyplot as plt
                    #    bestfit = self.bestfit(out_linemodel, redshift, emlinewave, resolution_matrix, camerapix)
                    #    plt.clf()
                    #    plt.plot(emlinewave, emlineflux, label='data', color='gray', lw=4)
                    #    plt.plot(emlinewave, bestfit, label='bestfit', ls='--', lw=3, alpha=0.7, color='k')
                    #    plt.plot(10**self.log10wave, log10model, label='hires model')
                    #    plt.plot(_emlinewave[icam], model_resamp, label='resamp')
                    #    plt.plot(_emlinewave[icam], model_convol, label='convol', lw=2)
                    #    plt.xlim(5386, 5394)
                    #    plt.legend()
                    #    plt.savefig('desi-users/ioannis/tmp/junk.png')
        
        return out_linemodel

    @staticmethod
    def chi2(linemodel, emlinewave, emlineflux, emlineivar, emlinemodel,
             continuum_model=None, return_dof=False):
        """Compute the reduced chi^2."""

        nfree = np.count_nonzero((linemodel['fixed'] == False) * (linemodel['tiedtoparam'] == -1))
        dof = np.count_nonzero(emlineivar > 0) - nfree

        if dof > 0:
            if continuum_model is None:
                model = emlinemodel
            else:
                model = emlinemodel + continuum_model
            chi2 = np.sum(emlineivar * (emlineflux - model)**2) / dof
        else:
            chi2 = 0.0
            
        if return_dof:
            return chi2, dof, nfree
        else:
            return chi2

    def bestfit(self, linemodel, redshift, emlinewave, resolution_matrix, camerapix):
        """Construct the best-fitting emission-line spectrum from a linemodel."""

        from fastspecfit.emlines import build_emline_model

        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, \
                     doubletpair, linewaves) = self._linemodel_to_parameters(linemodel, self.fit_linetable)

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line

        # doublets
        lineamps[doubletindx] *= lineamps[doubletpair]

        emlinemodel = build_emline_model(self.dlog10wave, redshift, lineamps, linevshifts, linesigmas,
                                         linewaves, emlinewave, resolution_matrix, camerapix)

        return emlinemodel

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table, redshift=None, 
                            snrcut=None):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (used for QA and elsewhere).

        """
        from fastspecfit.emlines import build_emline_model

        if redshift is None:
            redshift = fastspecfit_table['Z']
        
        linewaves = self.linetable['restwave'].data

        parameters = []
        for param in self.param_names:
            if '_amp' in param:
                param = param.replace('_amp', '_modelamp')
            parameters.append(fastspecfit_table[param.upper()])
        #parameters = [fastspecfit_table[param.upper()] for param in self.param_names]

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line    

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        lineamps[self.doubletindx] *= lineamps[self.doubletpair]

        if snrcut is not None:
            lineamps_ivar = [fastspecfit_table[param.upper()+'_AMP_IVAR'] for param in self.linetable['name']]
            lineamps[lineamps * np.sqrt(lineamps_ivar) < snrcut] = 0.
            
        emlinemodel = build_emline_model(self.dlog10wave, redshift, lineamps, 
                                         linevshifts, linesigmas, linewaves, 
                                         specwave, specres, None)

        return emlinemodel

    def populate_emtable(self, result, finalfit, finalmodel, emlinewave, emlineflux,
                         emlineivar, oemlineivar, specflux_nolines, redshift,
                         resolution_matrix, camerapix, log, nminpix=7, nsigma=3.):
        """Populate the output table with the emission-line measurements.

        """
        from scipy.stats import sigmaclip
        from fastspecfit.util import centers2edges
        from scipy.special import erf

        for param in finalfit:
            val = param['value']
            obsval = param['obsvalue']

            # special case the tied doublets
            if param['param_name'] == 'oii_doublet_ratio':
                result['OII_DOUBLET_RATIO'] = val
                result['OII_3726_MODELAMP'] = val * finalfit[param['doubletpair']]['value']
                result['OII_3726_AMP'] = val * finalfit[param['doubletpair']]['obsvalue']
            elif param['param_name'] == 'sii_doublet_ratio':
                result['SII_DOUBLET_RATIO'] = val
                result['SII_6731_MODELAMP'] = val * finalfit[param['doubletpair']]['value']
                result['SII_6731_AMP'] = val * finalfit[param['doubletpair']]['obsvalue']
            elif param['param_name'] == 'mgii_doublet_ratio':
                result['MGII_DOUBLET_RATIO'] = val
                result['MGII_2796_MODELAMP'] = val * finalfit[param['doubletpair']]['value']
                result['MGII_2796_AMP'] = val * finalfit[param['doubletpair']]['obsvalue']
            else:
                if '_amp' in param['param_name']:
                    result[param['param_name'].upper().replace('AMP', 'MODELAMP')] = val
                    result[param['param_name'].upper()] = obsval
                else:
                    result[param['param_name'].upper()] = val                    

        gausscorr = erf(nsigma / np.sqrt(2))      # correct for the flux outside of +/-nsigma
        dpixwave = np.median(np.diff(emlinewave)) # median pixel size [Angstrom]

        # Where the cameras overlap, we have to account for the variable pixel
        # size by sorting in wavelength.
        Wsrt = np.argsort(emlinewave)
        dwaves = np.diff(centers2edges(emlinewave[Wsrt]))

        # zero out all out-of-range lines
        for oneline in self.fit_linetable[~self.fit_linetable['inrange']]:
            linename = oneline['name'].upper()
            #print(linename, result['{}_AMP'.format(linename)], result['{}_MODELAMP'.format(linename)],
            #      result['{}_SIGMA'.format(linename)], result['{}_VSHIFT'.format(linename)])
            result['{}_AMP'.format(linename)] = 0.0
            result['{}_MODELAMP'.format(linename)] = 0.0
            result['{}_VSHIFT'.format(linename)] = 0.0
            result['{}_SIGMA'.format(linename)] = 0.0

        # get continuum fluxes, EWs, and upper limits
        narrow_sigmas, broad_sigmas, uv_sigmas = [], [], []
        narrow_redshifts, broad_redshifts, uv_redshifts = [], [], []
        for oneline in self.fit_linetable[self.fit_linetable['inrange']]:

            linename = oneline['name'].upper()
            linez = redshift + result['{}_VSHIFT'.format(linename)] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)
            linesigma = result['{}_SIGMA'.format(linename)] # [km/s]

            # if the line was dropped, use a default sigma value
            if linesigma == 0:
                if oneline['isbroad']:
                    if oneline['isbalmer']:
                        linesigma = self.limitsigma_narrow
                    else:
                        linesigma = self.limitsigma_broad
                else:
                    linesigma = self.limitsigma_broad

            linesigma_ang = linesigma * linezwave / C_LIGHT    # [observed-frame Angstrom]

            # require at least 2 pixels
            if linesigma_ang < 2 * dpixwave:
                linesigma_ang_window = 2 * dpixwave
                _gausscorr = 1.
            else:
                linesigma_ang_window = linesigma_ang
                _gausscorr = gausscorr

            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next line.
            lineindx = np.where((emlinewave >= (linezwave - nsigma*linesigma_ang_window)) *
                                (emlinewave <= (linezwave + nsigma*linesigma_ang_window)))[0]
            
            if len(lineindx) > 0 and np.sum(oemlineivar[lineindx] == 0) / len(lineindx) > 0.3: # use original ivar
                result['{}_AMP'.format(linename)] = 0.0
                result['{}_MODELAMP'.format(linename)] = 0.0
                result['{}_VSHIFT'.format(linename)] = 0.0
                result['{}_SIGMA'.format(linename)] = 0.0
            else:
                # number of pixels, chi2, and boxcar integration
                lineindx = np.where((emlinewave >= (linezwave - nsigma*linesigma_ang_window)) *
                                    (emlinewave <= (linezwave + nsigma*linesigma_ang_window)) *
                                    (emlineivar > 0))[0]

                npix = len(lineindx)
                result['{}_NPIX'.format(linename)] = npix
    
                if npix >= nminpix: # magic number: required at least XX unmasked pixels centered on the line
                    
                    if np.any(emlineivar[lineindx] == 0):
                        errmsg = 'Ivar should never be zero within an emission line!'
                        log.critical(errmsg)
                        raise ValueError(errmsg)

                    lineindx_Wsrt = np.where((emlinewave[Wsrt] >= (linezwave - nsigma*linesigma_ang_window)) *
                                             (emlinewave[Wsrt] <= (linezwave + nsigma*linesigma_ang_window)) *
                                             (emlineivar[Wsrt] > 0))[0]

                    # boxcar integration of the flux
                    #dwave = np.median(np.abs(np.diff(emlinewave_edges[lineindx])))
                    boxflux = np.sum(emlineflux[Wsrt][lineindx_Wsrt] * dwaves[lineindx_Wsrt])
                    boxflux_ivar = 1 / np.sum((1 / emlineivar[Wsrt][lineindx_Wsrt]) * dwaves[lineindx_Wsrt]**2)

                    result['{}_BOXFLUX'.format(linename)] = boxflux # * u.erg/(u.second*u.cm**2)
                    result['{}_BOXFLUX_IVAR'.format(linename)] = boxflux_ivar # * u.second**2*u.cm**4/u.erg**2
                    
                    # Get the uncertainty in the line-amplitude based on the scatter
                    # in the pixel values from the emission-line subtracted
                    # spectrum.
                    amp_sigma = np.diff(np.percentile(specflux_nolines[lineindx], [25, 75]))[0] / 1.349 # robust sigma
                    #clipflux, _, _ = sigmaclip(specflux_nolines[lineindx], low=3, high=3)
                    #amp_sigma = np.std(clipflux)
                    if amp_sigma > 0:
                        result['{}_AMP_IVAR'.format(linename)] = 1 / amp_sigma**2 # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
    
                    # require amp > 0 (line not dropped) to compute the flux and chi2
                    if result['{}_MODELAMP'.format(linename)] > 0:
    
                        result['{}_CHI2'.format(linename)] = np.sum(emlineivar[lineindx] * (emlineflux[lineindx] - finalmodel[lineindx])**2)
                        
                        lineprofile = build_emline_model(self.dlog10wave, redshift,
                                                         np.array([result['{}_MODELAMP'.format(linename)]]),
                                                         np.array([result['{}_VSHIFT'.format(linename)]]),
                                                         np.array([result['{}_SIGMA'.format(linename)]]),
                                                         np.array([oneline['restwave']]), emlinewave,
                                                         resolution_matrix, camerapix, debug=False)
                        
                        if np.sum(lineprofile) == 0. or np.any(lineprofile) < 0.:
                            errmsg = 'Line-profile should never be zero or negative!'
                            log.critical(errmsg)
                            raise ValueError(errmsg)
                        
                        # theoretical Gaussian line-flux and the (wrong) weighted flux_ivar
                        #flux = result['{}_MODELAMP'.format(linename)] * np.sqrt(2.0 * np.pi) * linesigma_ang # * u.Angstrom
                        #flux_ivar = np.sum(lineprofile[lineindx])**2 / np.sum(lineprofile[lineindx]**2 / emlineivar[lineindx])

                        # matched-filter (maximum-likelihood) Gaussian flux
                        pro_j = lineprofile[Wsrt][lineindx_Wsrt] / np.sum(lineprofile[Wsrt][lineindx_Wsrt])
                        I = pro_j > 0. # very narrow lines can have a profile that goes to zero
                        weight_j = (pro_j[I]**2 / dwaves[lineindx_Wsrt][I]**2) * emlineivar[Wsrt][lineindx_Wsrt][I]
                        flux = np.sum(weight_j * dwaves[lineindx_Wsrt][I] * lineprofile[Wsrt][lineindx_Wsrt][I] / pro_j[I]) / np.sum(weight_j)
                        flux_ivar = np.sum(weight_j)

                        # correction for missing flux
                        flux /= _gausscorr
                        flux_ivar *= _gausscorr**2

                        result['{}_FLUX'.format(linename)] = flux
                        result['{}_FLUX_IVAR'.format(linename)] = flux_ivar # * u.second**2*u.cm**4/u.erg**2
    
                        # keep track of sigma and z but only using XX-sigma lines
                        linesnr = result['{}_AMP'.format(linename)] * np.sqrt(result['{}_AMP_IVAR'.format(linename)])
                        #print(linename, result['{}_AMP'.format(linename)], amp_sigma, linesnr)
                        if linesnr > 1.5:
                            if oneline['isbroad']: # includes UV and broad Balmer lines
                                if oneline['isbalmer']:
                                    broad_sigmas.append(linesigma)
                                    broad_redshifts.append(linez)
                                else:
                                    uv_sigmas.append(linesigma)
                                    uv_redshifts.append(linez)
                            else:
                                narrow_sigmas.append(linesigma)
                                narrow_redshifts.append(linez)
    
                # next, get the continuum, the inverse variance in the line-amplitude, and the EW
                indxlo = np.where((emlinewave > (linezwave - 10*linesigma_ang_window)) *
                                  (emlinewave < (linezwave - 3.*linesigma_ang_window)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indxhi = np.where((emlinewave < (linezwave + 10*linesigma_ang_window)) *
                                  (emlinewave > (linezwave + 3.*linesigma_ang_window)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indx = np.hstack((indxlo, indxhi))

                if len(indx) >= nminpix: # require at least XX pixels to get the continuum level
                    #_, cmed, csig = sigma_clipped_stats(specflux_nolines[indx], sigma=3.0)
                    clipflux, _, _ = sigmaclip(specflux_nolines[indx], low=3, high=3)
                    # corner case: if a portion of a camera is masked
                    if len(clipflux) > 0:
                        #cmed, csig = np.mean(clipflux), np.std(clipflux)
                        cmed = np.median(clipflux)
                        csig = np.diff(np.percentile(clipflux, [25, 75])) / 1.349 # robust sigma
                        if csig > 0:
                            civar = (np.sqrt(len(indx)) / csig)**2
                        else:
                            civar = 0.0
                    else:
                        cmed, civar = 0.0, 0.0
    
                    result['{}_CONT'.format(linename)] = cmed # * u.erg/(u.second*u.cm**2*u.Angstrom)
                    result['{}_CONT_IVAR'.format(linename)] = civar # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
    
                if result['{}_CONT'.format(linename)] != 0.0 and result['{}_CONT_IVAR'.format(linename)] != 0.0:
                    lineflux = result['{}_FLUX'.format(linename)]
                    #linefluxivar = result['{}_BOXFLUX_IVAR'.format(linename)]
                    linefluxivar = result['{}_FLUX_IVAR'.format(linename)]
                    if lineflux > 0 and linefluxivar > 0:
                        # add the uncertainties in flux and the continuum in quadrature
                        ew = lineflux / cmed / (1 + redshift) # rest frame [A]
                        ewivar = (1+redshift)**2 / (1 / (cmed**2 * linefluxivar) + lineflux**2 / (cmed**4 * civar))
                    else:
                        ew, ewivar = 0.0, 0.0
                        
                    # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                    fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar) # * u.erg/(u.second*u.cm**2)
                    ewlimit = fluxlimit * cmed / (1+redshift)
    
                    result['{}_EW'.format(linename)] = ew
                    result['{}_EW_IVAR'.format(linename)] = ewivar
                    result['{}_FLUX_LIMIT'.format(linename)] = fluxlimit 
                    result['{}_EW_LIMIT'.format(linename)] = ewlimit

            if 'debug' in log.name:
                for col in ('VSHIFT', 'SIGMA', 'MODELAMP', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)]))
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'BOXFLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT', 'EW_LIMIT'):
                    log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)]))
                print()
                #log.debug(' ')
    
            ## simple QA
            #if linename == 'OIII_5007':
            #    import matplotlib.pyplot as plt
            #    _indx = np.arange(indx[-1]-indx[0])+indx[0]
            #    # continuum bandpasses and statistics
            #    plt.clf()
            #    plt.plot(emlinewave[_indx], specflux_nolines[_indx], color='gray')
            #    plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
            #    plt.axhline(y=cmed, color='k')
            #    plt.axhline(y=cmed+1/np.sqrt(civar), color='k', ls='--')
            #    plt.axhline(y=cmed-1/np.sqrt(civar), color='k', ls='--')
            #    plt.savefig('desi-users/ioannis/tmp/junk.png')
            #
            #    # emission-line integration
            #    plt.clf()
            #    plt.plot(emlinewave[_indx], emlineflux[_indx], color='gray')
            #    plt.plot(emlinewave[_indx], finalmodel[_indx], color='red')
            #    #plt.plot(emlinewave[_indx], specflux_nolines[_indx], color='orange', alpha=0.5)
            #    plt.axvline(x=emlinewave[lineindx[0]], color='blue')
            #    plt.axvline(x=emlinewave[lineindx[-1]], color='blue')
            #    plt.axhline(y=0, color='k', ls='--')
            #    plt.axhline(y=amp_sigma, color='k', ls='--')
            #    plt.axhline(y=2*amp_sigma, color='k', ls='--')
            #    plt.axhline(y=3*amp_sigma, color='k', ls='--')
            #    plt.axhline(y=result['{}_AMP'.format(linename)], color='k', ls='-')
            #    plt.savefig('desi-users/ioannis/tmp/junk2.png')

        # Clean up the doublets whose amplitudes were tied in the fitting since
        # they may have been zeroed out in the clean-up, above. This should be
        # smarter.
        if result['OIII_5007_MODELAMP'] == 0.0 and result['OIII_5007_NPIX'] > 0:
            result['OIII_4959_MODELAMP'] = 0.0
            result['OIII_4959_AMP'] = 0.0
            result['OIII_4959_FLUX'] = 0.0
            result['OIII_4959_EW'] = 0.0
        if result['NII_6584_MODELAMP'] == 0.0 and result['NII_6584_NPIX'] > 0:
            result['NII_6548_MODELAMP'] = 0.0
            result['NII_6548_AMP'] = 0.0
            result['NII_6548_FLUX'] = 0.0
            result['NII_6548_EW'] = 0.0
        if result['OII_7320_MODELAMP'] == 0.0 and result['OII_7320_NPIX'] > 0:
            result['OII_7330_MODELAMP'] = 0.0
            result['OII_7330_AMP'] = 0.0
            result['OII_7330_FLUX'] = 0.0
            result['OII_7330_EW'] = 0.0
        if result['MGII_2796_MODELAMP'] == 0.0 and result['MGII_2803_MODELAMP'] == 0.0:
            result['MGII_DOUBLET_RATIO'] = 0.0
        if result['OII_3726_MODELAMP'] == 0.0 and result['OII_3729_MODELAMP'] == 0.0:
            result['OII_DOUBLET_RATIO'] = 0.0
        if result['SII_6716_MODELAMP'] == 0.0 and result['SII_6731_MODELAMP'] == 0.0:
            result['SII_DOUBLET_RATIO'] = 0.0

        if 'debug' in log.name:
            for col in ('MGII_DOUBLET_RATIO', 'OII_DOUBLET_RATIO', 'SII_DOUBLET_RATIO'):
                log.debug('{}: {:.4f}'.format(col, result[col]))
            #log.debug(' ')
            print()

        # get the average emission-line redshifts and velocity widths
        if len(narrow_redshifts) > 0:
            result['NARROW_Z'] = np.median(narrow_redshifts)
            result['NARROW_SIGMA'] = np.median(narrow_sigmas) # * u.kilometer / u.second
            result['NARROW_ZRMS'] = np.std(narrow_redshifts)
            result['NARROW_SIGMARMS'] = np.std(narrow_sigmas)
        else:
            result['NARROW_Z'] = redshift
            
        if len(broad_redshifts) > 0:
            result['BROAD_Z'] = np.median(broad_redshifts)
            result['BROAD_SIGMA'] = np.median(broad_sigmas) # * u.kilometer / u.second
            result['BROAD_ZRMS'] = np.std(broad_redshifts)
            result['BROAD_SIGMARMS'] = np.std(broad_sigmas)
        else:
            result['BROAD_Z'] = redshift

        if len(uv_redshifts) > 0:
            result['UV_Z'] = np.median(uv_redshifts)
            result['UV_SIGMA'] = np.median(uv_sigmas) # * u.kilometer / u.second
            result['UV_ZRMS'] = np.std(uv_redshifts)
            result['UV_SIGMARMS'] = np.std(uv_sigmas)
        else:
            result['UV_Z'] = redshift

        # fragile
        if 'debug' in log.name:
            for line in ('NARROW', 'BROAD', 'UV'):
                log.debug('{}_Z: {:.9f}+/-{:.9f}'.format(line, result['{}_Z'.format(line)], result['{}_ZRMS'.format(line)]))
                log.debug('{}_SIGMA: {:.3f}+/-{:.3f}'.format(line, result['{}_SIGMA'.format(line)], result['{}_SIGMARMS'.format(line)]))

    def synthphot_spectrum(self, data, result, modelwave, modelflux):
        """Synthesize photometry from the best-fitting model (continuum+emission lines).

        """
        filters = self.synth_filters[data['photsys']]

        # Pad (simply) in wavelength...
        padflux, padwave = filters.pad_spectrum(modelflux, modelwave, method='edge')
        synthmaggies = filters.get_ab_maggies(padflux / FLUXNORM, padwave)
        synthmaggies = synthmaggies.as_array().view('f8')
        model_synthphot = self.parse_photometry(self.synth_bands, maggies=synthmaggies,
                                                nanomaggies=False,
                                                lambda_eff=filters.effective_wavelengths.value)

        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_{}'.format(band.upper())] = data['synthphot']['nanomaggies'][iband] # * 'nanomaggies'
            #result['FLUX_SYNTH_IVAR_{}'.format(band.upper())] = data['synthphot']['nanomaggies_ivar'][iband]
        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_SPECMODEL_{}'.format(band.upper())] = model_synthphot['nanomaggies'][iband] # * 'nanomaggies'

    def qa_fastspec(self, data, fastspec, metadata, coadd_type='healpix',
                    spec_wavelims=(3550, 9900), phot_wavelims=(0.1, 35),
                    fastphot=False, outprefix=None, outdir=None):
        """QA plot the emission-line spectrum and best-fitting model.

        """
        import subprocess
        from scipy.ndimage import median_filter

        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        from matplotlib import colors
        from matplotlib.patches import Circle, Rectangle, ConnectionPatch
        from matplotlib.lines import Line2D
        import matplotlib.gridspec as gridspec

        import astropy.units as u
        from astropy.io import fits
        from astropy.wcs import WCS
        import seaborn as sns
        from PIL import Image, ImageDraw

        from fastspecfit.util import ivar2var
        from fastspecfit.emlines import build_emline_model
        
        Image.MAX_IMAGE_PIXELS = None

        sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['dodgerblue', 'darkseagreen', 'orangered']]
        col2 = [colors.to_hex(col) for col in ['darkblue', 'darkgreen', 'darkred']]
        col3 = [colors.to_hex(col) for col in ['blue', 'green', 'red']]

        photcol1 = colors.to_hex('darkorange')
        #photcol1 = colors.to_hex('darkblue') # 'darkgreen', 'darkred', 'dodgerblue', 'darkseagreen', 'orangered']]
        
        if outdir is None:
            outdir = '.'
        if outprefix is None:
            outprefix = 'fastspec'

        if metadata['PHOTSYS'] == 'S':
            filters = self.decam
            allfilters = self.decamwise
        else:
            filters = self.bassmzls
            allfilters = self.bassmzlswise

        #survey = None
        #if 'SURVEY' not in metadata.colnames:
        #    if np.any([metadata['DESI_TARGET'] > 0, metadata['BGS_TARGET'] > 0, metadata['MWS_TARGET'] > 0, metadata['SCND_TARGET'] > 0]):
        #        survey = 'Main'
        #    else:
        #        for checksurvey in ['SV1', 'SV2', 'SV3']:
        #            #bit = []
        #            for targ in ['DESI', 'BGS', 'MWS', 'SCND']:
        #                if metadata['{}_{}_TARGET'.format(checksurvey, targ)] > 0:
        #                    survey = checksurvey.lower()
        #                    #bit += ['{}_{}_TARGET: {}'.format(checksurvey, targ, metadata['{}_{}_TARGET'.format(checksurvey, targ)])]
        #            if survey is not None:
        #                break
        #        if survey is None:
        #            if metadata['CMX_TARGET'] == 0:
        #                errmsg = 'Unable to determine survey!'
        #                log.critical(errmsg)
        #                raise ValueError(errmsg)
        #            else:
        #                survey = 'CMX'
        #else:

        if coadd_type == 'healpix':
            target = [
                'Survey/Program/Healpix: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX']),
                'TargetID: {}'.format(metadata['TARGETID']),
                ]
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
        elif coadd_type == 'cumulative':
            target = [
                'Tile/Night/Fiber: {}/{}/{}'.format(metadata['TILEID'], metadata['NIGHT'], metadata['FIBER']),
                #'Survey/Program/Tile: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['TILEID']),
                #'Night/Fiber: {}/{}'.format(metadata['NIGHT'], metadata['FIBER']),
                'TargetID: {}'.format(metadata['TARGETID']),
            ]
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], coadd_type, metadata['TARGETID']))
        elif coadd_type == 'pernight':
            target = [
                #'Survey/Program/Tile: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['TILEID']),
                'Tile/Night/Fiber: {}/{}/{}'.format(metadata['TILEID'], metadata['NIGHT'], metadata['FIBER']),
                'TargetID: {}'.format(metadata['TARGETID']),
            ]
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID']))
        elif coadd_type == 'perexp':
            target = [
                #'Survey/Program/Tile: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['TILEID']),
                'Tile/Night/Fiber: {}/{}/{}'.format(metadata['TILEID'], metadata['NIGHT'], metadata['FIBER']),
                'Night/Fiber: {}/{}'.format(metadata['NIGHT'], metadata['FIBER']),
                'TargetID: {}'.format(metadata['TARGETID']),
            ]
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'],
                    metadata['EXPID'], metadata['TARGETID']))
        elif coadd_type == 'custom':
            target = [
                'Survey/Program/Healpix: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX']),
                'TargetID: {}'.format(metadata['TARGETID']),
                ]
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
        else:
            errmsg = 'Unrecognized coadd_type {}!'.format(coadd_type)
            log.critical(errmsg)
            raise ValueError(errmsg)

        #print('Hack to skip existing files!')
        #if os.path.isfile(pngfile):
        #    return

        #target += bit

        apercorr = fastspec['APERCORR']
        redshift = fastspec['Z']

        leg = {
            'radec': '$(\\alpha,\\delta)=({:.7f},{:.6f})$'.format(metadata['RA'], metadata['DEC']),
            #'targetid': '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),

            'z': '$z={:.7f}$'.format(redshift),
            'zwarn': '$z_{{\\rm warn}}={}$'.format(metadata['ZWARN']),
            'dn4000_model': '$D_{{n}}(4000)_{{\\rm model}}={:.3f}$'.format(fastspec['DN4000_MODEL']),

            'rchi2': '$\\chi^{{2}}_{{\\nu, \\rm specphot}}$={:.2f}'.format(fastspec['RCHI2']), 
            'rchi2_cont': '$\\chi^{{2}}_{{\\nu, \\rm cont}}$={:.2f}'.format(fastspec['RCHI2_CONT']),
            'rchi2_phot': '$\\chi^{{2}}_{{\\nu, \\rm phot}}$={:.2f}'.format(fastspec['RCHI2_PHOT']),

            'age': 'Age$={:.3f}$ Gyr'.format(fastspec['AGE']),
            'AV': '$A_{{V}}={:.3f}$ mag'.format(fastspec['AV']),
            'mstar': '$\\log_{{10}}(M/M_{{\odot}})={:.3f}$'.format(fastspec['LOGMSTAR']),
            'sfr': '${{\\rm SFR}}={:.1f}\ M_{{\odot}}/{{\\rm yr}}$'.format(fastspec['SFR']),
            #'fagn': '$f_{{\\rm AGN}}={:.3f}$'.format(fastspec['FAGN']),
            'zzsun': '$Z/Z_{{\\odot}}={:.3f}$'.format(fastspec['ZZSUN']),

            'absmag_r': '$M_{{0.1r}}={:.2f}$'.format(fastspec['ABSMAG_SDSS_R']),
            'absmag_gr': '$^{{0.1}}(g-r)={:.3f}$'.format(fastspec['ABSMAG_SDSS_G']-fastspec['ABSMAG_SDSS_R']),
            'absmag_rz': '$^{{0.1}}(r-z)={:.3f}$'.format(fastspec['ABSMAG_SDSS_R']-fastspec['ABSMAG_SDSS_Z']),       
            }

        if redshift != metadata['Z_RR']:
            leg['zredrock'] = '$z_{{\\rm Redrock}}={:.7f}$'.format(metadata['Z_RR'])

        if fastspec['VDISP_IVAR'] > 0:
            leg['vdisp'] = '$\\sigma_{{star}}={:.0f}\pm{:.0f}$ km/s'.format(fastspec['VDISP'], 1/np.sqrt(fastspec['VDISP_IVAR']))
        else:
            leg['vdisp'] = '$\\sigma_{{star}}={:g}$ km/s'.format(fastspec['VDISP'])

        if fastspec['DN4000_IVAR'] > 0:
            leg['dn4000_spec'] = '$D_{{n}}(4000)_{{\\rm data}}={:.3f}$'.format(fastspec['DN4000'])
            #leg.update({'dn4000_spec': '$D_{{n}}(4000)_{{\\rm spec}}={:.3f}\pm{:.3f}$'.format(fastspec['DN4000'], 1/np.sqrt(fastspec['DN4000_IVAR']))})

        # kinematics
        if fastspec['NARROW_Z'] != redshift:
            if fastspec['NARROW_ZRMS'] > 0:
                leg['dv_narrow'] = '$\\Delta v_{{\\rm narrow}}={:.0f}\pm{:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['NARROW_Z']-redshift), C_LIGHT*fastspec['NARROW_ZRMS'])
            else:
                leg['dv_narrow'] = '$\\Delta v_{{\\rm narrow}}={:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['NARROW_Z']-redshift))
        if fastspec['NARROW_SIGMA'] != 0.0:
            if fastspec['NARROW_SIGMARMS'] > 0:
                leg['sigma_narrow'] = '$\\sigma_{{\\rm narrow}}={:.0f}\pm{:.0f}$ km/s'.format(
                    fastspec['NARROW_SIGMA'], fastspec['NARROW_SIGMARMS'])
            else:
                leg['sigma_narrow'] = '$\\sigma_{{\\rm narrow}}={:.0f}$ km/s'.format(fastspec['NARROW_SIGMA'])

        snrcut = 1.5
        leg_broad, leg_narrow, leg_uv = {}, {}, {}

        if fastspec['UV_Z'] != redshift:
            if fastspec['UV_ZRMS'] > 0:
                leg_uv['dv_uv'] = '$\\Delta v_{{\\rm UV}}={:.0f}\pm{:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['UV_Z']-redshift), C_LIGHT*fastspec['UV_ZRMS'])
            else:
                leg_uv['dv_uv'] = '$\\Delta v_{{\\rm UV}}={:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['UV_Z']-redshift))
        if fastspec['UV_SIGMA'] != 0.0:
            if fastspec['UV_SIGMARMS'] > 0:
                leg_uv['sigma_uv'] = '$\\sigma_{{\\rm UV}}={:.0f}\pm{:.0f}$ km/s'.format(
                    fastspec['UV_SIGMA'], fastspec['UV_SIGMARMS'])
            else:
                leg_uv['sigma_uv'] = '$\\sigma_{{\\rm UV}}={:.0f}$ km/s'.format(fastspec['UV_SIGMA'])
        if fastspec['BROAD_Z'] != redshift:
            if fastspec['BROAD_ZRMS'] > 0:
                leg_broad['dv_broad'] = '$\\Delta v_{{\\rm broad}}={:.0f}\pm{:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['BROAD_Z']-redshift), C_LIGHT*fastspec['BROAD_ZRMS'])
            else:
                leg_broad['dv_broad'] = '$\\Delta v_{{\\rm broad}}={:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['BROAD_Z']-redshift))
        if fastspec['BROAD_SIGMA'] != 0.0:
            if fastspec['BROAD_SIGMARMS'] > 0:
                leg_broad['sigma_broad'] = '$\\sigma_{{\\rm broad}}={:.0f}\pm{:.0f}$ km/s'.format(
                    fastspec['BROAD_SIGMA'], fastspec['BROAD_SIGMARMS'])
            else:
                leg_broad['sigma_broad'] = '$\\sigma_{{\\rm broad}}={:.0f}$ km/s'.format(fastspec['BROAD_SIGMA'])

        # emission lines

        # UV
        if fastspec['LYALPHA_AMP']*np.sqrt(fastspec['LYALPHA_AMP_IVAR']) > snrcut:
            leg_uv['ewlya'] = 'EW(Ly$\\alpha)={:.1f}\ \\AA$'.format(fastspec['LYALPHA_EW'])
        if fastspec['CIV_1549_AMP']*np.sqrt(fastspec['CIV_1549_AMP_IVAR']) > snrcut:
            leg_uv['ewciv'] = 'EW(CIV)$={:.1f}\ \\AA$'.format(fastspec['CIV_1549_EW'])
        if fastspec['CIII_1908_AMP']*np.sqrt(fastspec['CIII_1908_AMP_IVAR']) > snrcut:
            leg_uv['ewciii'] = 'EW(CIII])$={:.1f}\ \\AA$'.format(fastspec['CIII_1908_EW'])
        if (fastspec['MGII_2796_AMP']*np.sqrt(fastspec['MGII_2796_AMP_IVAR']) > snrcut or
            fastspec['MGII_2803_AMP']*np.sqrt(fastspec['MGII_2803_AMP_IVAR']) > snrcut):
            leg_uv['ewmgii'] = 'EW(MgII)$={:.1f}\ \\AA$'.format(fastspec['MGII_2796_EW']+fastspec['MGII_2803_EW'])
            leg_uv['mgii_doublet'] = 'MgII $\lambda2796/\lambda2803={:.3f}$'.format(fastspec['MGII_DOUBLET_RATIO'])

        leg_broad['linerchi2'] = '$\\chi^{{2}}_{{\\nu,\\rm line}}$={:.2f}'.format(fastspec['RCHI2_LINE'])
        leg_broad['deltarchi2'] = '$\\Delta\\chi^{{2}}_{{\\nu,\\rm nobroad}}={:.2f}$'.format(fastspec['DELTA_LINERCHI2'])

        # choose one broad Balmer line
        if fastspec['HALPHA_BROAD_AMP']*np.sqrt(fastspec['HALPHA_BROAD_AMP_IVAR']) > snrcut:
            leg_broad['ewbalmer_broad'] = 'EW(H$\\alpha)_{{\\rm broad}}={:.1f}\ \\AA$'.format(fastspec['HALPHA_BROAD_EW'])
        elif fastspec['HBETA_BROAD_AMP']*np.sqrt(fastspec['HBETA_BROAD_AMP_IVAR']) > snrcut:
            leg_broad['ewbalmer_broad'] = 'EW(H$\\beta)_{{\\rm broad}}={:.1f}\ \\AA$'.format(fastspec['HBETA_BROAD_EW'])
        elif fastspec['HGAMMA_BROAD_AMP']*np.sqrt(fastspec['HGAMMA_BROAD_AMP_IVAR']) > snrcut:
            leg_broad['ewbalmer_broad'] = 'EW(H$\\gamma)_{{\\rm broad}}={:.1f}\ \\AA$'.format(fastspec['HGAMMA_BROAD_EW'])

        if (fastspec['OII_3726_AMP']*np.sqrt(fastspec['OII_3726_AMP_IVAR']) > snrcut or 
            fastspec['OII_3729_AMP']*np.sqrt(fastspec['OII_3729_AMP_IVAR']) > snrcut):
            #leg_narrow['ewoii'] = 'EW([OII] $\lambda\lambda3726,29)={:.1f}\,\\AA$'.format(fastspec['OII_3726_EW']+fastspec['OII_3729_EW'])
            leg_narrow['ewoii'] = 'EW([OII])$={:.1f}\ \\AA$'.format(fastspec['OII_3726_EW']+fastspec['OII_3729_EW'])

        if fastspec['OIII_5007_AMP']*np.sqrt(fastspec['OIII_5007_AMP_IVAR']) > snrcut:
            leg_narrow['ewoiii'] = 'EW([OIII])$={:.1f}\,\\AA$'.format(fastspec['OIII_5007_EW'])
            #leg_narrow['ewoiii'] = 'EW([OIII] $\lambda5007={:.1f}\,\\AA$'.format(fastspec['OIII_5007_EW'])

        # choose one Balmer line
        if fastspec['HALPHA_AMP']*np.sqrt(fastspec['HALPHA_AMP_IVAR']) > snrcut:
            leg_narrow['ewbalmer_narrow'] = 'EW(H$\\alpha)={:.1f}\ \\AA$'.format(fastspec['HALPHA_EW'])
        elif fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut:
            leg_narrow['ewbalmer_narrow'] = 'EW(H$\\beta)={:.1f}\ \\AA$'.format(fastspec['HBETA_EW'])
        elif fastspec['HGAMMA_AMP']*np.sqrt(fastspec['HGAMMA_AMP_IVAR']) > snrcut:
            leg_narrow['ewbalmer_narrow'] = 'EW(H$\\gamma)={:.1f}\ \\AA$'.format(fastspec['HGAMMA_EW'])

        if (fastspec['HALPHA_AMP']*np.sqrt(fastspec['HALPHA_AMP_IVAR']) > snrcut and 
            fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut):
            leg_narrow['hahb'] = '${{\\rm H}}\\alpha/{{\\rm H}}\\beta={:.3f}$'.format(fastspec['HALPHA_FLUX']/fastspec['HBETA_FLUX'])
        if 'hahb' not in leg_narrow.keys() and (fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut and 
            fastspec['HGAMMA_AMP']*np.sqrt(fastspec['HGAMMA_AMP_IVAR']) > snrcut):
            leg_narrow['hbhg'] = '${{\\rm H}}\\beta/{{\\rm H}}\\gamma={:.3f}$'.format(fastspec['HBETA_FLUX']/fastspec['HGAMMA_FLUX'])
        if (fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut and 
            fastspec['OIII_5007_AMP']*np.sqrt(fastspec['OIII_5007_AMP_IVAR']) > snrcut and 
            fastspec['HBETA_FLUX'] > 0 and fastspec['OIII_5007_FLUX'] > 0):
            leg_narrow['oiiihb'] = '$\\log_{{10}}({{\\rm [OIII]/H}}\\beta)={:.3f}$'.format(np.log10(fastspec['OIII_5007_FLUX']/fastspec['HBETA_FLUX']))
        if (fastspec['HALPHA_AMP']*np.sqrt(fastspec['HALPHA_AMP_IVAR']) > snrcut and 
            fastspec['NII_6584_AMP']*np.sqrt(fastspec['NII_6584_AMP_IVAR']) > snrcut and 
            fastspec['HALPHA_FLUX'] > 0 and fastspec['NII_6584_FLUX'] > 0):
            leg_narrow['niiha'] = '$\\log_{{10}}({{\\rm [NII]/H}}\\alpha)={:.3f}$'.format(np.log10(fastspec['NII_6584_FLUX']/fastspec['HALPHA_FLUX']))

        if (fastspec['OII_3726_AMP']*np.sqrt(fastspec['OII_3726_AMP_IVAR']) > snrcut or 
            fastspec['OII_3729_AMP']*np.sqrt(fastspec['OII_3729_AMP_IVAR']) > snrcut):
            #if fastspec['OII_DOUBLET_RATIO'] != 0:
            leg_narrow['oii_doublet'] = '[OII] $\lambda3726/\lambda3729={:.3f}$'.format(fastspec['OII_DOUBLET_RATIO'])

        if fastspec['SII_6716_AMP']*np.sqrt(fastspec['SII_6716_AMP_IVAR']) > snrcut or fastspec['SII_6731_AMP']*np.sqrt(fastspec['SII_6731_AMP_IVAR']) > snrcut:
            #if fastspec['SII_DOUBLET_RATIO'] != 0:
            leg_narrow['sii_doublet'] = '[SII] $\lambda6731/\lambda6716={:.3f}$'.format(fastspec['SII_DOUBLET_RATIO'])

        # rebuild the best-fitting broadband photometric model
        sedmodel, sedphot = self.templates2data(
            self.templateflux, self.templatewave, 
            redshift=redshift, synthphot=True, 
            coeff=fastspec['COEFF'] * self.massnorm)
        sedwave = self.templatewave * (1 + redshift)

        phot = self.parse_photometry(self.bands,
                                     maggies=np.array([metadata['FLUX_{}'.format(band.upper())] for band in self.bands]),
                                     ivarmaggies=np.array([metadata['FLUX_IVAR_{}'.format(band.upper())] for band in self.bands]),
                                     lambda_eff=allfilters.effective_wavelengths.value,
                                     min_uncertainty=self.min_uncertainty)
        fiberphot = self.parse_photometry(self.fiber_bands,
                                          maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in self.fiber_bands]),
                                          lambda_eff=filters.effective_wavelengths.value)

        indx_phot = np.where((sedmodel > 0) * (sedwave/1e4 > phot_wavelims[0]) * 
                             (sedwave/1e4 < phot_wavelims[1]))[0]
        sedwave = sedwave[indx_phot]
        sedmodel = sedmodel[indx_phot]

        # Rebuild the best-fitting spectroscopic model; prefix "desi" means
        # "per-camera" and prefix "full" has the cameras h-stacked.
        fullwave = np.hstack(data['wave'])

        desicontinuum, _ = self.templates2data(self.templateflux_nolines, self.templatewave, 
                                               redshift=redshift, synthphot=False,
                                               specwave=data['wave'], specres=data['res'],
                                               specmask=data['mask'], cameras=data['cameras'],
                                               vdisp=fastspec['VDISP'],
                                               coeff=fastspec['COEFF'])

        # remove the aperture correction
        desicontinuum = [_desicontinuum / apercorr for _desicontinuum in desicontinuum]
        fullcontinuum = np.hstack(desicontinuum)

         # Need to be careful we don't pass a large negative residual where
         # there are gaps in the data.
        desiresiduals = []
        for icam in np.arange(len(data['cameras'])):
            resid = data['flux'][icam] - desicontinuum[icam]
            I = (data['flux'][icam] == 0.0) * (data['flux'][icam] == 0.0)
            if np.any(I):
                resid[I] = 0.0
            desiresiduals.append(resid)
        
        if np.all(fastspec['COEFF'] == 0):
            fullsmoothcontinuum = np.zeros_like(fullwave)
        else:
            fullsmoothcontinuum, _ = self.smooth_continuum(
                fullwave, np.hstack(desiresiduals), np.hstack(data['ivar']), 
                redshift=redshift, linemask=np.hstack(data['linemask']))

        desismoothcontinuum = []
        for campix in data['camerapix']:
            desismoothcontinuum.append(fullsmoothcontinuum[campix[0]:campix[1]])

        # full model spectrum + individual line-spectra
        desiemlines = self.emlinemodel_bestfit(data['wave'], data['res'], fastspec)

        desiemlines_oneline = []
        inrange = ( (self.linetable['restwave'] * (1+redshift) > np.min(fullwave)) *
                    (self.linetable['restwave'] * (1+redshift) < np.max(fullwave)) )
         
        for oneline in self.linetable[inrange]: # for all lines in range
            linename = oneline['name'].upper()
            amp = fastspec['{}_AMP'.format(linename)]
            if amp != 0:
                desiemlines_oneline1 = build_emline_model(
                    self.dlog10wave, redshift, np.array([amp]),
                    np.array([fastspec['{}_VSHIFT'.format(linename)]]),
                    np.array([fastspec['{}_SIGMA'.format(linename)]]),
                    np.array([oneline['restwave']]), data['wave'], data['res'])
                desiemlines_oneline.append(desiemlines_oneline1)

        # Grab the viewer cutout.
        pixscale = 0.262
        width = int(30 / pixscale)   # =1 arcmin
        height = int(width / 1.3) # 3:2 aspect ratio

        hdr = fits.Header()
        hdr['NAXIS'] = 2
        hdr['NAXIS1'] = width
        hdr['NAXIS2'] = height
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['CRVAL1'] = metadata['RA']
        hdr['CRVAL2'] = metadata['DEC']
        hdr['CRPIX1'] = width/2+0.5
        hdr['CRPIX2'] = height/2+0.5
        hdr['CD1_1'] = -pixscale/3600
        hdr['CD1_2'] = 0.0
        hdr['CD2_1'] = 0.0
        hdr['CD2_2'] = +pixscale/3600
        wcs = WCS(hdr)

        cutoutpng = os.path.join('/tmp', 'tmp.'+os.path.basename(pngfile))
        if not os.path.isfile(cutoutpng):
            cmd = 'wget -O {outfile} https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&width={width}&height={height}&layer=ls-dr9'
            cmd = cmd.format(outfile=cutoutpng, ra=metadata['RA'], dec=metadata['DEC'],
                             width=width, height=height)
            print(cmd)
            cuterr = subprocess.call(cmd.split())
            if cuterr != 0:
                errmsg = 'Something went wrong retrieving the png cutout'
                log.warning(errmsg)
                #log.critical(errmsg)
                #raise ValueError(errmsg)
        else:
            cuterr = 0

        # QA choices

        # 8 columns: 3 for the SED, 5 for the spectra, and 8 for the lines
        # 8 rows: 4 for the SED, 2 each for the spectra, 1 gap, and 3 for the lines
        ngaprows = 1
        nlinerows = 6
        nlinecols = 3
        nrows = 9 + ngaprows
        ncols = 8

        fullheight = 18 # inches
        fullwidth = 24

        height_ratios = np.hstack(([1.0]*3, 0.25, [1.0]*6)) # small gap
        width_ratios = np.hstack(([1.0]*5, [1.0]*3))
    
        fig = plt.figure(figsize=(fullwidth, fullheight))
        gs = fig.add_gridspec(nrows, ncols, height_ratios=height_ratios, width_ratios=width_ratios)

        cutax = fig.add_subplot(gs[0:3, 5:8], projection=wcs) # rows x cols
        sedax = fig.add_subplot(gs[0:3, 0:5])
        specax = fig.add_subplot(gs[4:8, 0:5])
        
        legxpos, legypos, legypos2, legfntsz1, legfntsz = 0.98, 0.94, 0.05, 16, 18
        bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.15)
        bbox2 = dict(boxstyle='round', facecolor='lightgray', alpha=0.7)

        # viewer cutout
        cuterr2 = 0
        if cuterr == 0:
            try:
                with Image.open(cutoutpng) as im:
                    sz = im.size
                    cutax.imshow(im, origin='lower')#, interpolation='nearest')
                cuterr2 = 0
            except:
                cuterr2 = 1
                errmsg = 'Something went wrong reading the png cutout'
                log.warning(errmsg)
                
            cutax.set_xlabel('RA [J2000]')
            cutax.set_ylabel('Dec [J2000]')
    
            cutax.coords[1].set_ticks_position('r')
            cutax.coords[1].set_ticklabel_position('r')
            cutax.coords[1].set_axislabel_position('r')
    
            if metadata['DEC'] > 0:
                sgn = '+'
            else:
                sgn = ''
                
            cutax.text(0.04, 0.95, '$(\\alpha,\\delta)$=({:.7f}, {}{:.6f})'.format(metadata['RA'], sgn, metadata['DEC']),
                       ha='left', va='top', color='k', fontsize=18, bbox=bbox2,
                       transform=cutax.transAxes)

            if cuterr2 == 0:
                cutax.add_artist(Circle((sz[0] / 2, sz[1] / 2), radius=1.5/2/pixscale, facecolor='none', # DESI fiber=1.5 arcsec diameter
                                        edgecolor='red', ls='-', alpha=0.8))#, label='3" diameter'))
                cutax.add_artist(Circle((sz[0] / 2, sz[1] / 2), radius=10/2/pixscale, facecolor='none',
                                        edgecolor='red', ls='--', alpha=0.8))#, label='15" diameter'))
                handles = [Line2D([0], [0], color='red', lw=2, ls='-', label='1.5 arcsec'),
                           Line2D([0], [0], color='red', lw=2, ls='--', label='10 arcsec')]
    
                cutax.legend(handles=handles, loc='lower left', fontsize=18, facecolor='lightgray')
    
        # plot the full spectrum + best-fitting (total) model
        spec_ymin, spec_ymax = 1e6, -1e6

        specax.plot(fullwave/1e4, fullsmoothcontinuum, color='gray', alpha=0.4)
        specax.plot(fullwave/1e4, fullcontinuum, color='k', alpha=0.6)

        desimodelspec = []
        for icam in np.arange(len(data['cameras'])): # iterate over cameras
            wave = data['wave'][icam]
            flux = data['flux'][icam]
            modelflux = desiemlines[icam] + desicontinuum[icam] + desismoothcontinuum[icam]

            sigma, camgood = ivar2var(data['ivar'][icam], sigma=True, allmasked_ok=True, clip=0)

            wave = wave[camgood]
            flux = flux[camgood]
            sigma = sigma[camgood]
            modelflux = modelflux[camgood]

            desimodelspec.append(apercorr * (desicontinuum[icam] + desiemlines[icam]))
            
            # get the robust range
            filtflux = median_filter(flux, 51, mode='nearest')
            if np.sum(camgood) > 0:
                sigflux = np.diff(np.percentile(flux - modelflux, [25, 75]))[0] / 1.349 # robust
                if -2 * sigflux < spec_ymin:
                    spec_ymin = -2 * sigflux
                if 6 * sigflux > spec_ymax:
                    spec_ymax = 6 * sigflux
                if np.max(filtflux) > spec_ymax:
                    #print(icam, spec_ymax, np.max(filtflux), np.max(filtflux) * 1.2)
                    spec_ymax = np.max(filtflux) * 1.25
                if np.max(modelflux) > spec_ymax:
                    spec_ymax = np.max(modelflux) * 1.25
                #print(spec_ymin, spec_ymax)
        
            #specax.fill_between(wave, flux-sigma, flux+sigma, color=col1[icam], alpha=0.2)
            specax.plot(wave/1e4, flux, color=col1[icam], alpha=0.8)
            specax.plot(wave/1e4, modelflux, color=col2[icam], lw=2, alpha=0.8)

        fullmodelspec = np.hstack(desimodelspec)

        specax.spines[['top']].set_visible(False)        
        specax.set_xlim(spec_wavelims[0]/1e4, spec_wavelims[1]/1e4)
        specax.set_ylim(spec_ymin, spec_ymax)
        specax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        #specax.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        specax.set_ylabel(r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$')

        # photometric SED   
        abmag_good = phot['abmag_ivar'] > 0
        abmag_goodlim = phot['abmag_limit'] > 0
        
        if len(sedmodel) == 0:
            self.log.warning('Best-fitting photometric continuum is all zeros or negative!')
            if np.sum(abmag_good) > 0:
                medmag = np.median(phot['abmag'][abmag_good])
            elif np.sum(abmag_goodlim) > 0:
                medmag = np.median(phot['abmag_limit'][abmag_goodlim])
            else:
                medmag = 0.0
            sedmodel_abmag = np.zeros_like(self.templatewave) + medmag
        else:
            factor = 10**(0.4 * 48.6) * sedwave**2 / (C_LIGHT * 1e13) / FLUXNORM / self.massnorm # [erg/s/cm2/A --> maggies]
            sedmodel_abmag = -2.5*np.log10(sedmodel * factor)
            sedax.plot(sedwave / 1e4, sedmodel_abmag, color='slategrey', alpha=0.9, zorder=1)
    
        sedax.scatter(sedphot['lambda_eff']/1e4, sedphot['abmag'], marker='s', 
                      s=400, color='k', facecolor='none', linewidth=2, alpha=1.0, zorder=2)

        #factor = 10**(0.4 * 48.6) * fullwave**2 / (C_LIGHT * 1e13) / FLUXNORM # [erg/s/cm2/A --> maggies]
        #good = fullmodelspec > 0
        #sedax.plot(fullwave[good]/1e4, -2.5*np.log10(fullmodelspec[good]*factor[good]), color='purple', alpha=0.8)
        for icam in np.arange(len(data['cameras'])):
            factor = 10**(0.4 * 48.6) * data['wave'][icam]**2 / (C_LIGHT * 1e13) / FLUXNORM # [erg/s/cm2/A --> maggies]
            good = desimodelspec[icam] > 0
            sedax.plot(data['wave'][icam][good]/1e4, -2.5*np.log10(desimodelspec[icam][good]*factor[good]), color=col2[icam], alpha=0.8)

        # we have to set the limits *before* we call errorbar, below!
        dm = 1.5
        sed_ymin = np.nanmax(sedmodel_abmag) + dm
        sed_ymax = np.nanmin(sedmodel_abmag) - dm
        if np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) > 0:
            sed_ymin = np.max((np.nanmax(phot['abmag'][abmag_good]), np.nanmax(phot['abmag_limit'][abmag_goodlim]), np.nanmax(sedmodel_abmag))) + dm
            sed_ymax = np.min((np.nanmin(phot['abmag'][abmag_good]), np.nanmin(phot['abmag_limit'][abmag_goodlim]), np.nanmin(sedmodel_abmag))) - dm
        elif np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) == 0:
            sed_ymin = np.max((np.nanmax(phot['abmag'][abmag_good]), np.nanmax(sedmodel_abmag))) + dm
            sed_ymax = np.min((np.nanmin(phot['abmag'][abmag_good]), np.nanmin(sedmodel_abmag))) - dm
        elif np.sum(abmag_good) == 0 and np.sum(abmag_goodlim) > 0:
            sed_ymin = np.max((np.nanmax(phot['abmag_limit'][abmag_goodlim]), np.nanmax(sedmodel_abmag))) + dm
            sed_ymax = np.min((np.nanmin(phot['abmag_limit'][abmag_goodlim]), np.nanmin(sedmodel_abmag))) - dm
        else:
            abmag_good = phot['abmag'] > 0
            abmag_goodlim = phot['abmag_limit'] > 0
            if np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) > 0:
                sed_ymin = np.max((np.nanmax(phot['abmag'][abmag_good]), np.nanmax(phot['abmag_limit'][abmag_goodlim]))) + dm
                sed_ymax = np.min((np.nanmin(phot['abmag'][abmag_good]), np.nanmin(phot['abmag_limit'][abmag_goodlim]))) - dm
            elif np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) == 0:                
                sed_ymin = np.nanmax(phot['abmag'][abmag_good]) + dm
                sed_ymax = np.nanmin(phot['abmag'][abmag_good]) - dm
            elif np.sum(abmag_good) == 0 and np.sum(abmag_goodlim) > 0:
                sed_ymin = np.nanmax(phot['abmag_limit'][abmag_goodlim]) + dm
                sed_ymax = np.nanmin(phot['abmag_limit'][abmag_goodlim]) - dm
            #else:
            #    sed_ymin, sed_ymax = [30, 20]
            
        if sed_ymin > 30:
            sed_ymin = 30
        if np.isnan(sed_ymin) or np.isnan(sed_ymax):
            raise('Problem here!')
        print(sed_ymin, sed_ymax)

        #sedax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        sedax.set_xlim(phot_wavelims[0], phot_wavelims[1])
        sedax.set_xscale('log')
        sedax.set_ylabel('AB mag') 
        #sedax.set_ylabel(r'Apparent Brightness (AB mag)') 
        sedax.set_ylim(sed_ymin, sed_ymax)

        @ticker.FuncFormatter
        def major_formatter(x, pos):
            if (x >= 0.01) and (x < 0.1):
                return f'{x:.2f}'
            elif (x >= 0.1) and (x < 1):
                return f'{x:.1f}'
            else:
                return f'{x:.0f}'
        
        sedax.xaxis.set_major_formatter(major_formatter)
        obsticks = np.array([0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 20.0])
        sedax.set_xticks(obsticks)

        sedax_twin = sedax.twiny()
        sedax_twin.set_xlim(phot_wavelims[0]/(1+redshift), phot_wavelims[1]/(1+redshift))
        sedax_twin.set_xscale('log')
        #sedax_twin.set_xlabel(r'Rest-frame Wavelength ($\mu$m)') 

        sedax_twin.xaxis.set_major_formatter(major_formatter)
        restticks = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 15.0, 20.0])
        restticks = restticks[(restticks >= phot_wavelims[0]/(1+redshift)) * (restticks <= phot_wavelims[1]/(1+redshift))]
        sedax_twin.set_xticks(restticks)

        # integrated flux / photometry
        abmag = np.squeeze(phot['abmag'])
        abmag_limit = np.squeeze(phot['abmag_limit'])
        abmag_fainterr = np.squeeze(phot['abmag_fainterr'])
        abmag_brighterr = np.squeeze(phot['abmag_brighterr'])
        yerr = np.squeeze([abmag_brighterr, abmag_fainterr])

        markersize = 14

        dofit = np.where(self.bands_to_fit)[0]
        if len(dofit) > 0:
            good = np.where((abmag[dofit] > 0) * (abmag_limit[dofit] == 0))[0]
            upper = np.where(abmag_limit[dofit] > 0)[0]
            if len(good) > 0:
                sedax.errorbar(phot['lambda_eff'][dofit][good]/1e4, abmag[dofit][good],
                               yerr=yerr[:, dofit[good]],
                               fmt='o', markersize=markersize, markeredgewidth=1, markeredgecolor='k',
                               markerfacecolor=photcol1, elinewidth=3, ecolor=photcol1, capsize=4,
                               label=r'$grz\,W_{1}W_{2}W_{3}W_{4}$', zorder=2, alpha=1.0)
            if len(upper) > 0:
                sedax.errorbar(phot['lambda_eff'][dofit][upper]/1e4, abmag_limit[dofit][upper],
                               lolims=True, yerr=0.75,
                               fmt='o', markersize=markersize, markeredgewidth=3, markeredgecolor='k',
                               markerfacecolor=photcol1, elinewidth=3, ecolor=photcol1, capsize=4, alpha=0.7)
    
        ignorefit = np.where(self.bands_to_fit == False)[0]
        if len(ignorefit) > 0:
            good = np.where((abmag[ignorefit] > 0) * (abmag_limit[ignorefit] == 0))[0]
            upper = np.where(abmag_limit[ignorefit] > 0)[0]
            if len(good) > 0:
                sedax.errorbar(phot['lambda_eff'][ignorefit][good]/1e4, abmag[ignorefit][good],
                               yerr=yerr[:, ignorefit[good]],
                               fmt='o', markersize=markersize, markeredgewidth=3, markeredgecolor='k',
                               markerfacecolor='none', elinewidth=3, ecolor=photcol1, capsize=4, alpha=0.7)
            if len(upper) > 0:
                sedax.errorbar(phot['lambda_eff'][ignorefit][upper]/1e4, abmag_limit[ignorefit][upper],
                               lolims=True, yerr=0.75, fmt='o', markersize=markersize, markeredgewidth=3,
                               markeredgecolor='k', markerfacecolor='none', elinewidth=3,
                               ecolor=photcol1, capsize=5, alpha=0.7)

        # Label the DESI wavelength range and the aperture correction.
        sedax.plot([np.min(fullwave)/1e4, np.max(fullwave)/1e4], [sed_ymin-1, sed_ymin-1],
                   lw=2, ls='-', color='gray', marker='s')#, alpha=0.5)
        sedax.text(((np.max(fullwave)-np.min(fullwave))/2+np.min(fullwave)*0.8)/1e4, sed_ymin-1.7,
                   'DESI x {:.2f}'.format(apercorr), ha='center', va='center', fontsize=16,
                   color='k')

        txt = '\n'.join((leg['rchi2_cont'], leg['rchi2_phot'], leg['rchi2']))
        sedax.text(0.02, 0.94, txt, ha='left', va='top',
                    transform=sedax.transAxes, fontsize=legfntsz)#, bbox=bbox)

        txt = '\n'.join((
            #r'{}'.format(leg['fagn']),
            r'{}'.format(leg['zzsun']),
            r'{}'.format(leg['AV']),
            r'{}'.format(leg['sfr']),
            r'{}'.format(leg['age']),
            r'{}'.format(leg['mstar']),
            ))
        sedax.text(legxpos, legypos2, txt, ha='right', va='bottom',
                    transform=sedax.transAxes, fontsize=legfntsz1, bbox=bbox)

        # draw lines connecting the SED and spectral plots
        sedax.add_artist(ConnectionPatch(xyA=(spec_wavelims[0]/1e4, sed_ymin), 
                                         xyB=(spec_wavelims[0]/1e4, spec_ymax), 
                                         coordsA='data', coordsB='data',
                                         axesA=sedax, axesB=specax, color='k'))
        sedax.add_artist(ConnectionPatch(xyA=(spec_wavelims[1]/1e4, sed_ymin), 
                                         xyB=(spec_wavelims[1]/1e4, spec_ymax), 
                                         coordsA='data', coordsB='data',
                                         axesA=sedax, axesB=specax, color='k'))

        # zoom in on individual emission lines - use linetable!
        linetable = self.linetable
        inrange = (linetable['restwave'] * (1+redshift) > np.min(fullwave)) * (linetable['restwave'] * (1+redshift) < np.max(fullwave))
        linetable = linetable[inrange]

        nline = len(set(linetable['plotgroup']))

        plotsig_default = 200.0 # [km/s]
        plotsig_default_balmer = 500.0 # [km/s]
        plotsig_default_broad = 2000.0 # [km/s]

        meanwaves, deltawaves, sigmas, linenames = [], [], [], []
        for plotgroup in set(linetable['plotgroup']):
            I = np.where(plotgroup == linetable['plotgroup'])[0]
            linename = linetable['nicename'][I[0]].replace('-', ' ')
            linenames.append(linename)
            meanwaves.append(np.mean(linetable['restwave'][I]))
            deltawaves.append((np.max(linetable['restwave'][I]) - np.min(linetable['restwave'][I])) / 2)
        
            sigmas1 = np.array([fastspec['{}_SIGMA'.format(line.upper())] for line in linetable[I]['name']])
            sigmas1 = sigmas1[sigmas1 > 0]
            if len(sigmas1) > 0:
                plotsig = 1.5*np.mean(sigmas1)
                if plotsig < 50:
                    plotsig = 50.0
            else:
                if np.any(linetable['isbroad'][I]):
                    if np.any(linetable['isbalmer'][I]):
                        plotsig = fastspec['BROAD_SIGMA']
                        if plotsig < 50:
                            plotsig = fastspec['NARROW_SIGMA']
                            if plotsig < 50:
                                plotsig = plotsig_default
                                #plotsig = plotsig_default_broad
                    else:
                        plotsig = fastspec['UV_SIGMA']                    
                        if plotsig < 50:
                            plotsig = plotsig_default_broad
                else:
                    plotsig = fastspec['NARROW_SIGMA']
                    if plotsig < 50:
                        plotsig = plotsig_default
            sigmas.append(plotsig)
        
        srt = np.argsort(meanwaves)
        meanwaves = np.hstack(meanwaves)[srt]
        deltawaves = np.hstack(deltawaves)[srt]
        sigmas = np.hstack(sigmas)[srt]
        linenames = np.hstack(linenames)[srt]

        # Add the linenames to the spectrum plot.
        for meanwave, linename in zip(meanwaves*(1+redshift), linenames):
            #print(meanwave, ymax_spec)
            if meanwave > spec_wavelims[0] and meanwave < spec_wavelims[1]:
                if 'SiIII' in linename:
                    thislinename = '\n'+linename.replace('+', '+\n  ')
                elif '4363' in linename:
                    thislinename = linename+'\n'
                else:
                    thislinename = linename
                specax.text(meanwave/1e4, spec_ymax, thislinename, ha='center', va='top',
                            rotation=270, fontsize=12, alpha=0.5)
    
        removelabels = np.ones(nline, bool)
        line_ymin, line_ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
        
        ax, irow, colshift = [], 4, 5 # skip the gap row
        for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):
            icol = iax % nlinecols
            icol += colshift
            if iax > 0 and iax % nlinecols == 0:
                irow += 1
            #print(iax, irow, icol)
        
            xx = fig.add_subplot(gs[irow, icol])
            ax.append(xx)
        
            wmin = (meanwave - deltawave) * (1+redshift) - 6 * sig * meanwave * (1+redshift) / C_LIGHT
            wmax = (meanwave + deltawave) * (1+redshift) + 6 * sig * meanwave * (1+redshift) / C_LIGHT
            #print(linename, wmin, wmax)
        
            # iterate over cameras
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                emlinewave = data['wave'][icam]
                emlineflux = data['flux'][icam] - desicontinuum[icam] - desismoothcontinuum[icam]
                emlinemodel = desiemlines[icam]
        
                emlinesigma, good = ivar2var(data['ivar'][icam], sigma=True, allmasked_ok=True, clip=0)
                emlinewave = emlinewave[good]
                emlineflux = emlineflux[good]
                emlinesigma = emlinesigma[good]
                emlinemodel = emlinemodel[good]
        
                emlinemodel_oneline = []
                for desiemlines_oneline1 in desiemlines_oneline:
                    emlinemodel_oneline.append(desiemlines_oneline1[icam][good])
                    
                indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
                if len(indx) > 1:
                    removelabels[iax] = False
                    xx.plot(emlinewave[indx]/1e4, emlineflux[indx], color=col1[icam], alpha=0.5)
                    #xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                    #                emlineflux[indx]+emlinesigma[indx], color=col1[icam], alpha=0.5)
                    # plot the individual lines first then the total model
                    for emlinemodel_oneline1 in emlinemodel_oneline:
                        if np.sum(emlinemodel_oneline1[indx]) > 0:
                            #P = emlinemodel_oneline1[indx] > 0
                            xx.plot(emlinewave[indx]/1e4, emlinemodel_oneline1[indx], lw=1, alpha=0.8, color=col2[icam])
                    xx.plot(emlinewave[indx]/1e4, emlinemodel[indx], color=col2[icam], lw=2, alpha=0.8)

                    #xx.plot(emlinewave[indx], emlineflux[indx]-emlinemodel[indx], color='gray', alpha=0.3)
                    #xx.axhline(y=0, color='gray', ls='--')
        
                    # get the robust range
                    sigflux = np.std(emlineflux[indx])
                    filtflux = median_filter(emlineflux[indx], 3, mode='nearest')
        
                    #_line_ymin, _line_ymax = -1.5 * sigflux, 4 * sigflux
                    #if np.max(emlinemodel[indx]) > _line_ymax:
                    #    _line_ymax = np.max(emlinemodel[indx]) * 1.3
                    _line_ymin, _line_ymax = -1.5 * sigflux, np.max(emlinemodel[indx]) * 1.4
                    if 4 * sigflux > _line_ymax:
                        _line_ymax = 4 * sigflux
                    if np.max(filtflux) > _line_ymax:
                        _line_ymax = np.max(filtflux)
                    if np.min(emlinemodel[indx]) < _line_ymin:
                        _line_ymin = 0.8 * np.min(emlinemodel[indx])
                    if _line_ymax > line_ymax[iax]:
                        line_ymax[iax] = _line_ymax
                    if _line_ymin < line_ymin[iax]:
                        line_ymin[iax] = _line_ymin
                    #print(linename, line_ymin[iax], line_ymax[iax])
        
                    xx.set_xlim(wmin/1e4, wmax/1e4)
                    
                xx.text(0.03, 0.89, linename, ha='left', va='center',
                        transform=xx.transAxes, fontsize=12)
                xx.tick_params(axis='x', labelsize=16)
                xx.tick_params(axis='y', labelsize=16)
                
        for iax, xx in enumerate(ax):
            if removelabels[iax]:
                xx.set_ylim(0, 1)
                xx.set_xticklabels([])
                xx.set_yticklabels([])
            else:
                xx.set_yticklabels([])
                xx.set_ylim(line_ymin[iax], line_ymax[iax])
                xx_twin = xx.twinx()
                xx_twin.set_ylim(line_ymin[iax], line_ymax[iax])
                xlim = xx.get_xlim()
                xx.xaxis.set_major_locator(ticker.MaxNLocator(2))
        
        plt.subplots_adjust(wspace=0.4, top=0.9, bottom=0.1, left=0.07, right=0.92, hspace=0.33)

        # common axis labels
        if len(ax) > 0:
            if len(ax) == 2:
                xx = fig.add_subplot(gs[irow, icol+1])
                xx.axis('off')
                ax.append(xx)
            elif len(ax) == 1:
                xx = fig.add_subplot(gs[irow, icol+1])
                xx.axis('off')
                ax.append(xx)
                xx = fig.add_subplot(gs[irow, icol+2])
                xx.axis('off')
                ax.append(xx)
            
            ulpos = ax[0].get_position()
            lpos = ax[nline-1].get_position()
            urpos = ax[2].get_position()
            xpos = (urpos.x1 - ulpos.x0) / 2 + ulpos.x0# + 0.03
                
            ypos = lpos.y0 - 0.04
            fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                     ha='center', va='center', fontsize=24)
    
            xpos = urpos.x1 + 0.05
            ypos = (urpos.y1 - lpos.y0) / 2 + lpos.y0# + 0.03
            fig.text(xpos, ypos, 
                     r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                     ha='center', va='center', rotation=270, fontsize=24)

        ppos = sedax.get_position()
        xpos = (ppos.x1 - ppos.x0) / 2 + ppos.x0
        ypos = ppos.y1 + 0.03
        fig.text(xpos, ypos, r'Rest-frame Wavelength ($\mu$m)',
                 ha='center', va='bottom', fontsize=24)

        fig.text(0.647, 0.925, '\n'.join(target), ha='left', va='bottom',
                 fontsize=23, linespacing=1.4)

        #fig.suptitle(title, fontsize=22)

        # add some key results about the object at the bottom of the figure
        legkeys = leg.keys()

        legfntsz = 17
        #toppos, startpos, deltapos = 0.21, 0.04, 0.13
        toppos, leftpos, rightpos, adjust = 0.21, 0.03, 0.62, 0.01
        nbox = 2 + 1*bool(leg_narrow) + bool(leg_broad)
        boxpos = np.arange(nbox) * (rightpos - leftpos)/nbox + leftpos

        txt = [
            r'{}'.format(leg['z']),
        ]
        if 'zredrock' in legkeys:
            txt += [r'{}'.format(leg['zredrock'])]

        txt += [            
            #r'{}'.format(leg['zwarn']),
            r'{}'.format(leg['vdisp']),
        ]
        if 'dv_narrow' in legkeys or 'dv_uv' in leg_uv.keys() or 'dv_broad' in leg_broad.keys():
            txt += ['']

        if 'dv_narrow' in legkeys:
            txt += [
                r'{}'.format(leg['sigma_narrow']),
                r'{}'.format(leg['dv_narrow']),
            ]
        if 'dv_uv' in leg_uv.keys():
            txt += [
                r'{}'.format(leg_uv['sigma_uv']),
                r'{}'.format(leg_uv['dv_uv']),
            ]
            _ = leg_uv.pop('sigma_uv')
            _ = leg_uv.pop('dv_uv')
        if 'dv_broad' in leg_broad.keys():
            txt += [
                r'{}'.format(leg_broad['sigma_broad']),
                r'{}'.format(leg_broad['dv_broad']),
            ]
            _ = leg_broad.pop('sigma_broad')
            _ = leg_broad.pop('dv_broad')

        ibox = 0
        #fig.text(startpos, toppos, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz, 
        fig.text(boxpos[ibox], toppos, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz, 
                 bbox=bbox, linespacing=1.4)
        ibox += 1

        txt = [
            r'{}'.format(leg['absmag_r']),
            r'{}'.format(leg['absmag_gr']),
            r'{}'.format(leg['absmag_rz']), 
            '',
            r'{}'.format(leg['dn4000_model'])
        ]
        if 'dn4000_spec' in legkeys:
            txt += [r'{}'.format(leg['dn4000_spec'])]
            
        #fig.text(startpos+deltapos, toppos, '\n'.join(txt), ha='left', va='top', 
        fig.text(boxpos[ibox], toppos, '\n'.join(txt), ha='left', va='top', 
                 fontsize=legfntsz, bbox=bbox, linespacing=1.4)
        ibox += 1        

        #factor = 2
        if bool(leg_narrow):
            txt = []
            for key in leg_narrow.keys():
                txt += [r'{}'.format(leg_narrow[key])]
            #fig.text(startpos+deltapos*factor, toppos, '\n'.join(txt), ha='left', va='top',
            fig.text(boxpos[ibox]-adjust*2, toppos, '\n'.join(txt), ha='left', va='top',
                     fontsize=legfntsz, bbox=bbox, linespacing=1.4)
            ibox += 1        
            #factor += 1.25

        if bool(leg_broad):
            txt = []
            for key in leg_broad.keys():
                txt += [r'{}'.format(leg_broad[key])]

        if bool(leg_uv):
            if bool(leg_broad):
                txt += ['']
            else:
                txt = []
            for key in leg_uv.keys():
                txt += [r'{}'.format(leg_uv[key])]

        if bool(leg_uv) or bool(leg_broad):
            #fig.text(startpos+deltapos*factor, toppos, '\n'.join(txt), ha='left', va='top',
            fig.text(boxpos[ibox]-adjust*1, toppos, '\n'.join(txt), ha='left', va='top',
                     fontsize=legfntsz, bbox=bbox, linespacing=1.4)

        self.log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)#, dpi=150)
        plt.close()

def emline_specfit(data, result, continuummodel, smooth_continuum,
                   minsnr_balmer_broad=3., fphoto=None, emlinesfile=None,
                   synthphot=True, broadlinefit=True,
                   percamera_models=False, log=None, verbose=False):
    """Perform the fit minimization / chi2 minimization.

    Parameters
    ----------
    data
    continuummodel
    smooth_continuum
    synthphot
    verbose
    broadlinefit

    Returns
    -------
    results
    modelflux
 
    """
    from astropy.table import Column, vstack
    from fastspecfit.util import ivar2var

    tall = time.time()

    if log is None:
        from desiutil.log import get_logger, DEBUG
        if verbose:
            log = get_logger(DEBUG)
        else:
            log = get_logger()

    EMFit = EMFitTools(emlinesfile=emlinesfile, fphoto=fphoto, uniqueid=data['uniqueid'],
                       minsnr_balmer_broad=minsnr_balmer_broad)

    # Combine all three cameras; we will unpack them to build the
    # best-fitting model (per-camera) below.
    redshift = data['zredrock']
    emlinewave = np.hstack(data['wave'])
    oemlineivar = np.hstack(data['ivar'])
    specflux = np.hstack(data['flux'])
    #resolution_matrix = data['res']
    resolution_matrix = data['res_fast']
    camerapix = data['camerapix']

    continuummodelflux = np.hstack(continuummodel)
    smoothcontinuummodelflux = np.hstack(smooth_continuum)
    emlineflux = specflux - continuummodelflux - smoothcontinuummodelflux

    emlineivar = np.copy(oemlineivar)
    emlinevar, emlinegood = ivar2var(emlineivar, clip=1e-3)
    emlinebad = np.logical_not(emlinegood)

    # This is a (dangerous???) hack.
    if np.sum(emlinebad) > 0:
        emlineivar[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineivar[emlinegood])
        emlineflux[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineflux[emlinegood]) # ???

    weights = np.sqrt(emlineivar)

    # Build all the emission-line models for this object.
    linemodel_broad, linemodel_nobroad = EMFit.build_linemodels(
        redshift, wavelims=(np.min(emlinewave)+5, np.max(emlinewave)-5),
        verbose=False, strict_broadmodel=True)
    #EMFit.summarize_linemodel(linemodel_broad)
    #EMFit.summarize_linemodel(linemodel_nobroad)

    # Get initial guesses on the parameters and populate the two "initial"
    # linemodels; the "final" linemodels will be initialized with the
    # best-fitting parameters from the initial round of fitting.
    initial_guesses, param_bounds = EMFit.initial_guesses_and_bounds(
        data, emlinewave, emlineflux, log)

    EMFit.populate_linemodel(linemodel_nobroad, initial_guesses, param_bounds, log)
    EMFit.populate_linemodel(linemodel_broad, initial_guesses, param_bounds, log)

    # Initial fit - initial_linemodel_nobroad
    t0 = time.time()
    fit_nobroad = optimize(EMFit, linemodel_nobroad, emlinewave, emlineflux, weights, redshift,
                           resolution_matrix, camerapix, log=log, debug=False, get_finalamp=True)
    model_nobroad = emfit_bestfit(EMFit, fit_nobroad, redshift, emlinewave, resolution_matrix, camerapix)
    chi2_nobroad, ndof_nobroad, nfree_nobroad = EMFit.chi2(fit_nobroad, emlinewave, emlineflux, emlineivar, model_nobroad, return_dof=True)
    log.info('Line-fitting with no broad lines and {} free parameters took {:.4f} seconds [niter={}, rchi2={:.4f}].'.format(
        nfree_nobroad, time.time()-t0, fit_nobroad.meta['nfev'], chi2_nobroad))
    
    # Now try adding broad Balmer and helium lines and see if we improve the
    # chi2.
    if broadlinefit:
        # Gather the pixels around the broad Balmer lines and the corresponding
        # linemodel table.
        balmer_pix, balmer_linemodel_broad, balmer_linemodel_nobroad = [], [], []
        for icam in np.arange(len(data['cameras'])):
            pixoffset = int(np.sum(data['npixpercamera'][:icam]))
            if len(data['linename'][icam]) > 0:
                I = (linemodel_nobroad['isbalmer'] * (linemodel_nobroad['ishelium'] == False) *
                     linemodel_nobroad['isbroad'] * np.isin(linemodel_nobroad['linename'], data['linename'][icam]))
                _balmer_linemodel_broad = linemodel_broad[I]
                _balmer_linemodel_nobroad = linemodel_nobroad[I]
                balmer_linemodel_broad.append(_balmer_linemodel_broad)
                balmer_linemodel_nobroad.append(_balmer_linemodel_nobroad)
                if len(_balmer_linemodel_broad) > 0: # use balmer_linemodel_broad not balmer_linemodel_nobroad
                    I = np.where(np.isin(data['linename'][icam], _balmer_linemodel_broad['linename']))[0]
                    for ii in I:
                        #print(data['linename'][icam][ii])
                        balmer_pix.append(data['linepix'][icam][ii] + pixoffset)
                        
        if len(balmer_pix) > 0:
            t0 = time.time()
            fit_broad = optimize(EMFit, linemodel_broad, emlinewave, emlineflux, weights, 
                                 redshift, resolution_matrix, camerapix, log=log,
                                 debug=False, get_finalamp=True)
            model_broad = emfit_bestfit(EMFit, fit_broad, redshift, emlinewave, resolution_matrix, camerapix)
            #fit_broad = EMFit.optimize(linemodel_broad, emlinewave, emlineflux, weights, 
            #                           redshift, resolution_matrix, camerapix, log=log,
            #                           debug=False, get_finalamp=True)
            #model_broad = EMFit.bestfit(fit_broad, redshift, emlinewave, resolution_matrix, camerapix)
            chi2_broad, ndof_broad, nfree_broad = EMFit.chi2(fit_broad, emlinewave, emlineflux, emlineivar, model_broad, return_dof=True)
            log.info('Line-fitting with broad lines and {} free parameters took {:.4f} seconds [niter={}, rchi2={:.4f}].'.format(
                nfree_broad, time.time()-t0, fit_broad.meta['nfev'], chi2_broad))

            # compute delta-chi2 around just the Balmer lines
            balmer_pix = np.hstack(balmer_pix)
            balmer_linemodel_broad = vstack(balmer_linemodel_broad)

            balmer_nfree_broad = (np.count_nonzero((balmer_linemodel_broad['fixed'] == False) *
                                                   (balmer_linemodel_broad['tiedtoparam'] == -1)))
            balmer_ndof_broad = np.count_nonzero(emlineivar[balmer_pix] > 0) - balmer_nfree_broad

            balmer_linemodel_nobroad = vstack(balmer_linemodel_nobroad)
            balmer_nfree_nobroad = (np.count_nonzero((balmer_linemodel_nobroad['fixed'] == False) *
                                                     (balmer_linemodel_nobroad['tiedtoparam'] == -1)))
            balmer_ndof_nobroad = np.count_nonzero(emlineivar[balmer_pix] > 0) - balmer_nfree_nobroad

            linechi2_balmer_broad = np.sum(emlineivar[balmer_pix] * (emlineflux[balmer_pix] - model_broad[balmer_pix])**2)
            linechi2_balmer_nobroad = np.sum(emlineivar[balmer_pix] * (emlineflux[balmer_pix] - model_nobroad[balmer_pix])**2)
            delta_linechi2_balmer = linechi2_balmer_nobroad - linechi2_balmer_broad
            delta_linendof_balmer = balmer_ndof_nobroad - balmer_ndof_broad

            # Choose broad-line model only if:
            # --delta-chi2 > delta-ndof
            # --broad_sigma < narrow_sigma
            # --broad_sigma < 250

            dchi2test = delta_linechi2_balmer > delta_linendof_balmer
            Hanarrow = fit_broad['param_name'] == 'halpha_sigma' # Balmer lines are tied to H-alpha even if out of range
            Habroad = fit_broad['param_name'] == 'halpha_broad_sigma'
            Bbroad = fit_broad['isbalmer'] * fit_broad['isbroad'] * (fit_broad['fixed'] == False) * EMFit.amp_balmer_bool
            broadsnr = fit_broad[Bbroad]['obsvalue'].data * np.sqrt(fit_broad[Bbroad]['civar'].data)

            sigtest1 = fit_broad[Habroad]['value'][0] > EMFit.minsigma_balmer_broad
            sigtest2 = (fit_broad[Habroad]['value'] > fit_broad[Hanarrow]['value'])[0]
            if len(broadsnr) == 0:
                broadsnrtest = False
                _broadsnr = 0.
            elif len(broadsnr) == 1:
                broadsnrtest =  broadsnr[-1] > EMFit.minsnr_balmer_broad
                _broadsnr = 'S/N ({}) = {:.1f}'.format(fit_broad[Bbroad]['linename'][-1], broadsnr[-1])
            else:
                broadsnrtest =  np.any(broadsnr[-2:] > EMFit.minsnr_balmer_broad)
                _broadsnr = 'S/N ({}) = {:.1f}, S/N ({}) = {:.1f}'.format(
                    fit_broad[Bbroad]['linename'][-2], broadsnr[-2], fit_broad[Bbroad]['linename'][-1], broadsnr[-1])

            if dchi2test and sigtest1 and sigtest2 and broadsnrtest:
                log.info('Adopting broad-line model:')
                log.info('  delta-chi2={:.1f} > delta-ndof={:.0f}'.format(delta_linechi2_balmer, delta_linendof_balmer))
                log.info('  sigma_broad={:.1f} km/s, sigma_narrow={:.1f} km/s'.format(fit_broad[Habroad]['value'][0], fit_broad[Hanarrow]['value'][0]))
                if _broadsnr:
                    log.info('  {} > {:.0f}'.format(_broadsnr, EMFit.minsnr_balmer_broad))
                finalfit, finalmodel, finalchi2 = fit_broad, model_broad, chi2_broad
            else:
                if dchi2test == False:
                    log.info('Dropping broad-line model: delta-chi2={:.1f} < delta-ndof={:.0f}'.format(
                        delta_linechi2_balmer, delta_linendof_balmer))
                elif sigtest1 == False:
                    log.info('Dropping broad-line model: Halpha_broad_sigma {:.1f} km/s < {:.0f} km/s (delta-chi2={:.1f}, delta-ndof={:.0f}).'.format(
                        fit_broad[Habroad]['value'][0], EMFit.minsigma_balmer_broad, delta_linechi2_balmer, delta_linendof_balmer))
                elif sigtest2 == False:
                    log.info('Dropping broad-line model: Halpha_broad_sigma {:.1f} km/s < Halpha_narrow_sigma {:.1f} km/s (delta-chi2={:.1f}, delta-ndof={:.0f}).'.format(
                        fit_broad[Habroad]['value'][0], fit_broad[Hanarrow]['value'][0], delta_linechi2_balmer, delta_linendof_balmer))
                elif broadsnrtest == False:
                    log.info('Dropping broad-line model: {} < {:.0f}'.format(_broadsnr, EMFit.minsnr_balmer_broad))
                finalfit, finalmodel, finalchi2 = fit_nobroad, model_nobroad, chi2_nobroad
        else:
            log.info('Insufficient Balmer lines to test the broad-line model.')
            finalfit, finalmodel, finalchi2 = fit_nobroad, model_nobroad, chi2_nobroad
            delta_linechi2_balmer, delta_linendof_balmer = 0, np.int32(0)
    else:
        log.info('Skipping broad-line fitting (broadlinefit=False).')
        finalfit, finalmodel, finalchi2 = fit_nobroad, model_nobroad, chi2_nobroad
        delta_linechi2_balmer, delta_linendof_balmer = 0, np.int32(0)

    # Residual spectrum with no emission lines.
    specflux_nolines = specflux - finalmodel

    # Now fill the output table.
    EMFit.populate_emtable(result, finalfit, finalmodel, emlinewave, emlineflux,
                           emlineivar, oemlineivar, specflux_nolines, redshift,
                           resolution_matrix, camerapix, log)

    # Build the model spectra.
    emmodel = np.hstack(EMFit.emlinemodel_bestfit(data['wave'], data['res'], result, redshift=redshift))

    result['RCHI2_LINE'] = finalchi2
    #result['NDOF_LINE'] = finalndof
    result['DELTA_LINECHI2'] = delta_linechi2_balmer # chi2_nobroad - chi2_broad
    result['DELTA_LINENDOF'] = delta_linendof_balmer # ndof_nobroad - ndof_broad

    # full-fit reduced chi2
    rchi2 = np.sum(oemlineivar * (specflux - (continuummodelflux + smoothcontinuummodelflux + emmodel))**2)
    rchi2 /= np.sum(oemlineivar > 0) # dof??
    result['RCHI2'] = rchi2

    # I believe that all the elements of the coadd_wave vector are contained
    # within one or more of the per-camera wavelength vectors, and so we
    # should be able to simply map our model spectra with no
    # interpolation. However, because of round-off, etc., it's probably
    # easiest to use np.interp.

    # package together the final output models for writing; assume constant
    # dispersion in wavelength!
    minwave, maxwave, dwave = np.min(data['coadd_wave']), np.max(data['coadd_wave']), np.diff(data['coadd_wave'][:2])[0]
    minwave = float(int(minwave * 1000) / 1000)
    maxwave = float(int(maxwave * 1000) / 1000)
    dwave = float(int(np.round(dwave * 1000)) / 1000)
    npix = int(np.round((maxwave-minwave)/dwave)) + 1
    modelwave = minwave + dwave * np.arange(npix)

    modelspectra = Table()
    # all these header cards need to be 2-element tuples (value, comment),
    # otherwise io.write_fastspecfit will crash
    modelspectra.meta['NAXIS1'] = (npix, 'number of pixels')
    modelspectra.meta['NAXIS2'] = (npix, 'number of models')
    modelspectra.meta['NAXIS3'] = (npix, 'number of objects')
    modelspectra.meta['BUNIT'] = ('10**-17 erg/(s cm2 Angstrom)', 'flux unit')
    modelspectra.meta['CUNIT1'] = ('Angstrom', 'wavelength unit')
    modelspectra.meta['CTYPE1'] = ('WAVE', 'type of axis')
    modelspectra.meta['CRVAL1'] = (minwave, 'wavelength of pixel CRPIX1 (Angstrom)')
    modelspectra.meta['CRPIX1'] = (0, '0-indexed pixel number corresponding to CRVAL1')
    modelspectra.meta['CDELT1'] = (dwave, 'pixel size (Angstrom)')
    modelspectra.meta['DC-FLAG'] = (0, '0 = linear wavelength vector')
    modelspectra.meta['AIRORVAC'] = ('vac', 'wavelengths in vacuum (vac)')

    wavesrt = np.argsort(emlinewave)
    modelcontinuum = np.interp(modelwave, emlinewave[wavesrt], continuummodelflux[wavesrt]).reshape(1, npix)
    modelsmoothcontinuum = np.interp(modelwave, emlinewave[wavesrt], smoothcontinuummodelflux[wavesrt]).reshape(1, npix)
    modelemspectrum = np.interp(modelwave, emlinewave[wavesrt], emmodel[wavesrt]).reshape(1, npix)
    
    modelspectra.add_column(Column(name='CONTINUUM', dtype='f4', data=modelcontinuum))
    modelspectra.add_column(Column(name='SMOOTHCONTINUUM', dtype='f4', data=modelsmoothcontinuum))
    modelspectra.add_column(Column(name='EMLINEMODEL', dtype='f4', data=modelemspectrum))

    # Finally, optionally synthesize photometry (excluding the
    # smoothcontinuum!) and measure Dn(4000) from the line-free spectrum.
    if synthphot:
        modelflux = modelcontinuum[0, :] + modelemspectrum[0, :]
        EMFit.synthphot_spectrum(data, result, modelwave, modelflux)

    # measure DN(4000) without the emission lines
    if result['DN4000_IVAR'] > 0:
        fluxnolines = data['coadd_flux'] - modelemspectrum[0, :]
        dn4000_nolines, _ = EMFit.get_dn4000(modelwave, fluxnolines, redshift=redshift, log=log, rest=False)
        log.info('Dn(4000)={:.3f} in the emission-line subtracted spectrum.'.format(dn4000_nolines))
        result['DN4000'] = dn4000_nolines

        # Simple QA of the Dn(4000) estimate.
        if False:
            import matplotlib.pyplot as plt

            dn4000, dn4000_obs, dn4000_model, dn4000_ivar = result['DN4000'], result['DN4000_OBS'], result['DN4000_MODEL'], result['DN4000_IVAR']
            print(dn4000, dn4000_obs, dn4000_model, 1/np.sqrt(dn4000_ivar))
    
            restwave = modelwave / (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
            fnu_obs = data['coadd_flux'] * flam2fnu # [erg/s/cm2/Hz]
            fnu = fluxnolines * flam2fnu # [erg/s/cm2/Hz]

            fnu_model = modelcontinuum[0, :] * flam2fnu
            fnu_fullmodel = modelflux * flam2fnu
            
            fnu_ivar = data['coadd_ivar'] / flam2fnu**2            
            fnu_sigma, fnu_mask = ivar2var(fnu_ivar, sigma=True)
    
            I = (restwave > 3835) * (restwave < 4115)
            J = (restwave > 3835) * (restwave < 4115) * fnu_mask
    
            fig, ax = plt.subplots()
            ax.fill_between(restwave[I], fnu_obs[I]-fnu_sigma[I], fnu_obs[I]+fnu_sigma[I],
                            label='Observed Dn(4000)={:.3f}+/-{:.3f}'.format(dn4000_obs, 1/np.sqrt(dn4000_ivar)))
            ax.plot(restwave[I], fnu[I], color='blue', label='Line-free Dn(4000)={:.3f}+/-{:.3f}'.format(
                dn4000, 1/np.sqrt(dn4000_ivar)))
            ax.plot(restwave[I], fnu_fullmodel[I], color='k', label='Model Dn(4000)={:.3f}'.format(dn4000_model))
            ax.plot(restwave[I], fnu_model[I], color='red', label='Model Dn(4000)={:.3f}'.format(dn4000_model))
            ylim = ax.get_ylim()
            ax.fill_between([3850, 3950], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.fill_between([4000, 4100], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.set_xlabel(r'Rest Wavelength ($\AA$)')
            ax.set_ylabel(r'$F_{\nu}$ (erg/s/cm2/Hz)')
            ax.legend()
            fig.savefig('desi-users/ioannis/tmp/qa-dn4000.png')

    log.info('Emission-line fitting took {:.2f} seconds.'.format(time.time()-tall))

    if percamera_models:
        errmsg = 'percamera-models option not yet implemented.'
        log.critical(errmsg)
        raise NotImplementedError(errmsg)

    return modelspectra

