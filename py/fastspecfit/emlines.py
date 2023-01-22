"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import pdb # for debugging

import os, time
import numpy as np
import numba

from astropy.table import Table, Column

def read_emlines():
    """Read the set of emission lines of interest.

    """
    from pkg_resources import resource_filename
    
    linefile = resource_filename('fastspecfit', 'data/emlines.ecsv')    
    linetable = Table.read(linefile, format='ascii.ecsv', guess=False)
    
    return linetable    

#@numba.jit(nopython=True)
def build_emline_model(log10wave, redshift, lineamps, linevshifts, linesigmas, 
                       linewaves, emlinewave, resolution_matrix, camerapix=None):
    """Given parameters, build the model emission-line spectrum.

    ToDo: can this be optimized using numba?

    """
    from fastspecfit.util import trapz_rebin, C_LIGHT

    log10model = np.zeros_like(log10wave) # [erg/s/cm2/A, observed frame]

    # Cut to lines with non-zero amplitudes.
    #I = linesigmas > 0
    I = lineamps > 0
    if np.count_nonzero(I) > 0:
        linevshifts = linevshifts[I]
        linesigmas = linesigmas[I]
        lineamps = lineamps[I]
        linewaves = linewaves[I]

        # line-width [log-10 Angstrom] and redshifted wavelength [log-10 Angstrom]
        log10sigmas = linesigmas / C_LIGHT / np.log(10) 
        linezwaves = np.log10(linewaves * (1.0 + redshift + linevshifts / C_LIGHT))

        for lineamp, linezwave, log10sigma in zip(lineamps, linezwaves, log10sigmas):
            J = np.abs(log10wave - linezwave) < (8 * log10sigma) # cut to pixels within +/-N-sigma
            if np.count_nonzero(J) > 0:
                #print(lineamp, 10**linezwave, 10**log10wave[J].min(), 10**log10wave[J].max())
                log10model[J] = log10model[J] + lineamp * np.exp(-0.5 * (log10wave[J]-linezwave)**2 / log10sigma**2)

    # Optionally split into cameras, resample, and convolve with the resolution
    # matrix.
    emlinemodel = []
    if camerapix is None:
        for icam, specwave in enumerate(emlinewave):
            _emlinemodel = trapz_rebin(10**log10wave, log10model, specwave)
            _emlinemomdel = resolution_matrix[icam].dot(_emlinemodel)
            emlinemodel.append(_emlinemodel)
        return emlinemodel
    else:
        for icam, campix in enumerate(camerapix):
            _emlinemodel = trapz_rebin(10**log10wave, log10model, emlinewave[campix[0]:campix[1]])
            _emlinemomdel = resolution_matrix[icam].dot(_emlinemodel)
            emlinemodel.append(_emlinemodel)
        return np.hstack(emlinemodel)

def _objective_function(free_parameters, emlinewave, emlineflux, weights, redshift, 
                        log10wave, resolution_matrix, camerapix, parameters, Ifree, 
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
    emlinemodel = build_emline_model(log10wave, redshift, lineamps, linevshifts, 
                                     linesigmas, linewaves, emlinewave, 
                                     resolution_matrix, camerapix)

    if weights is None:
        residuals = emlinemodel - emlineflux
    else:
        residuals = weights * (emlinemodel - emlineflux)

    return residuals
