#
# Calculation of objective function and
# Jacobian for emission line fitting
#

import numpy as np
from math import erf, erfc

from numba import jit

from .params_mapping import ParamsMapping
from .sparse_rep import EMLineJacobian

from .utils import (
    MAX_SDEV,
    norm_pdf,
    norm_cdf,
    max_buffer_width
)

from .model import (
    emline_model,
    emline_perline_models
)

from .jacobian import emline_model_jacobian

###################################################################

class EMLine_Objective(object):

    def __init__(self,
                 obs_bin_centers,
                 obs_fluxes,
                 obs_weights,
                 redshift,
                 line_wavelengths,
                 resolution_matrices,
                 camerapix,
                 params_mapping):

        self.dtype = obs_fluxes.dtype
        
        self.obs_fluxes = obs_fluxes
        self.obs_weights = obs_weights
        self.redshift = redshift
        self.line_wavelengths = line_wavelengths
        self.resolution_matrices = tuple(resolution_matrices)
        self.camerapix = camerapix
        self.params_mapping = params_mapping
        
        self.log_obs_bin_edges, self.ibin_widths = \
            prepare_bins(obs_bin_centers, camerapix)
                
    #
    # Objective function for least-squares optimization
    #
    # Build the emline model as described above and compute the weighted
    # vector of residuals between the modeled fluxes and the observations.
    #
    def objective(self, free_parameters):
            
        #
        # expand free paramters into complete
        # parameter array, handling tied params
        # and doublets
        #
        parameters = self.params_mapping.mapFreeToFull(free_parameters)
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)
        
        model_fluxes = np.empty_like(self.obs_fluxes, dtype=self.dtype)
        
        for icam, campix in enumerate(self.camerapix):
            
            # start and end for obs fluxes of camera icam
            s, e = campix
            
            # Actual bin widths are in ibw[1..e-s].
            # Setup guarantees that ibw[0] and
            # ibw[e-s+1] are dummy values
            ibw = self.ibin_widths[s:e+2]

            mf = emline_model(self.line_wavelengths,
                              lineamps, linevshifts, linesigmas,
                              self.log_obs_bin_edges[s+icam:e+icam+1],
                              self.redshift,
                              ibw)
            
            # convolve model with resolution matrix and store in
            # this camera's subrange of model_fluxes
            self.resolution_matrices[icam].matvec(mf, model_fluxes[s:e])
        
        return self.obs_weights * (model_fluxes - self.obs_fluxes) # residuals


    #
    # Jacobian of objective function for least-squares EMLine
    # optimization. The result of the detailed calculation is converted
    # to a sparse matrix, since it is extremely sparse, to speed up
    # subsequent matrix-vector multiplies in the optimizer.
    #
    def jacobian(self, free_parameters):
        
        #
        # expand free paramters into complete
        # parameter array, handling tied params
        # and doublets
        #
        
        parameters = self.params_mapping.mapFreeToFull(free_parameters)
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3)

        J_S = self.params_mapping.getJacobian(free_parameters)

        jacs = []
        for icam, campix in enumerate(self.camerapix):
            s, e = campix

            # Actual bin widths are in ibw[1..e-s].
            # Setup guarantees that ibw[0] and
            # ibw[e-s+1] are not out of bounds.
            ibw = self.ibin_widths[s:e+1]
            
            idealJac = \
                emline_model_jacobian(lineamps, linevshifts, linesigmas,
                                      self.log_obs_bin_edges[s+icam:e+icam+1],
                                      ibw,
                                      self.redshift,
                                      self.line_wavelengths,
                                      self.resolution_matrices[icam].ndiag())
            
            # ignore any columns corresponding to fixed parameters
            endpts = idealJac[0]
            endpts[self.params_mapping.fixedMask(), :] = (0,0)
        
            jacs.append( mulWMJ(self.obs_weights[s:e],
                                self.resolution_matrices[icam].data,
                                idealJac) )
            
        nBins = np.sum(np.diff(self.camerapix))
        nFreeParms = len(free_parameters)
        nParms = len(parameters)
        J =  EMLineJacobian((nBins, nFreeParms), nParms,
                            self.camerapix, jacs, J_S)

        return J
    

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
                         obs_bin_centers,
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
    
    log_obs_bin_edges, ibin_widths = prepare_bins(obs_bin_centers,
                                                  camerapix)
    
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

            
##########################################################################

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
