import numpy as np

from numba import jit

from fastspecfit.util import C_LIGHT

from .utils import (
    MAX_SDEV,
    max_buffer_width,
)


@jit(nopython=True, nogil=True, cache=True)
def emline_model(line_wavelengths,
                 line_parameters,
                 log_obs_bin_centers,
                 redshift,
                 sigma_G):
    """Compute the Gaussian-convolved emission-line model at each observed pixel.

    Parameters
    ----------
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    log_obs_bin_centers : :class:`numpy.ndarray`
        Natural logs of observed wavelength bin centers.
    redshift : float
        Redshift of the observed spectrum.
    sigma_G : float
        Gaussian pre-convolution width in Angstroms (must match the value
        used to build the deconvolved resolution matrix W).

    Returns
    -------
    :class:`numpy.ndarray`
        Model flux at each observed pixel center.

    """

    line_amplitudes, line_vshifts, line_sigmas = \
        np.split(line_parameters, 3)

    nbins = len(log_obs_bin_centers)
    model_fluxes = np.zeros(nbins, dtype=line_amplitudes.dtype)

    # temporary buffer for per-line calculations
    bin_vals = np.empty(nbins, dtype=line_amplitudes.dtype)

    for j in range(len(line_wavelengths)):

        if line_amplitudes[j] == 0.:
            continue

        s, e = emline_model_core(line_wavelengths[j],
                                 line_amplitudes[j],
                                 line_vshifts[j],
                                 line_sigmas[j],
                                 log_obs_bin_centers,
                                 redshift,
                                 sigma_G,
                                 bin_vals)

        model_fluxes[s:e] += bin_vals[:e-s]

    return model_fluxes


@jit(nopython=True, nogil=True, cache=True)
def emline_perline_models(line_wavelengths,
                          line_parameters,
                          log_obs_bin_centers,
                          redshift,
                          sigma_G,
                          padding):
    """Compute the per-line emission-line flux profiles sparsely.

    Shares the core computation of :func:`emline_model` but returns
    individual line profiles rather than a collapsed composite.

    Parameters
    ----------
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    log_obs_bin_centers : :class:`numpy.ndarray`
        Natural logs of observed wavelength bin centers.
    redshift : float
        Redshift of the observed spectrum.
    sigma_G : float
        Gaussian pre-convolution width in Angstroms.
    padding : int
        Extra entries to pad each sparse row for later use.

    Returns
    -------
    endpts : :class:`numpy.ndarray` of int, shape (nlines, 2)
        Start and end bin indices of the nonzero range for each line.
    profiles : :class:`numpy.ndarray`, shape (nlines, max_width)
        Flux values within each line's nonzero bin range.

    """

    line_amplitudes, line_vshifts, line_sigmas = \
        np.split(line_parameters, 3)

    max_width = max_buffer_width(log_obs_bin_centers, line_sigmas, padding, sigma_G)

    nlines = len(line_wavelengths)
    line_profiles = np.empty((nlines, max_width), dtype=line_amplitudes.dtype)
    endpts        = np.zeros((nlines, 2), dtype=np.int32)

    for j in range(nlines):

        if line_amplitudes[j] == 0.:
            continue

        bin_vals = line_profiles[j]

        s, e = emline_model_core(line_wavelengths[j],
                                 line_amplitudes[j],
                                 line_vshifts[j],
                                 line_sigmas[j],
                                 log_obs_bin_centers,
                                 redshift,
                                 sigma_G,
                                 bin_vals)

        endpts[j,0] = s
        endpts[j,1] = e

    return (endpts, line_profiles)


@jit(nopython=True, nogil=True, cache=True)
def emline_model_core(line_wavelength,
                      line_amplitude,
                      line_vshift,
                      line_sigma,
                      log_obs_bin_centers,
                      redshift,
                      sigma_G,
                      vals):
    """Compute the Gaussian-convolved flux of one spectral line at each pixel.

    Evaluates ``F_sigma(lambda_j) = A * (sigma/sigma_eff) * exp(-x_j^2 /
    (2*sigma_eff^2))`` at each pixel center ``lambda_j`` in the nonzero
    range, where ``sigma_eff = sqrt(sigma^2 + (sigma_G/lambda*)^2)`` accounts
    for the Gaussian pre-convolution that the deconvolved resolution matrix W
    expects as input.

    Parameters
    ----------
    line_wavelength : float
        Nominal rest-frame wavelength of the line in Angstroms.
    line_amplitude : float
        Amplitude of the line.
    line_vshift : float
        Velocity shift of the line in km/s.
    line_sigma : float
        Intrinsic Gaussian width of the line in km/s.
    log_obs_bin_centers : :class:`numpy.ndarray`
        Natural logs of observed wavelength pixel centers.
    redshift : float
        Redshift of the observed spectrum.
    sigma_G : float
        Gaussian pre-convolution width in Angstroms (must match the value
        used to build the deconvolved resolution matrix W).
    vals : :class:`numpy.ndarray`
        Output array in which nonzero pixel fluxes are written starting at
        ``vals[0]``.

    Returns
    -------
    s : int
        Index of the first pixel with nonzero flux.
    e : int
        One past the index of the last pixel with nonzero flux.
        ``vals[0:e-s]`` contains the fluxes. If ``s == e``, no nonzero
        fluxes were found.

    """

    # intrinsic line width in log-lambda units
    sigma = line_sigma / C_LIGHT

    # observed line center
    line_shift       = 1. + redshift + line_vshift / C_LIGHT
    shifted_line     = line_wavelength * line_shift
    log_shifted_line = np.log(shifted_line)

    # effective sigma in log-lambda after Gaussian pre-convolution;
    # sigma_G / lambda* converts the linear-space sigma to log-space
    sigma_G_log = sigma_G / shifted_line
    sigma_eff   = np.sqrt(sigma**2 + sigma_G_log**2)

    # amplitude factor A * (sigma / sigma_eff)
    amp_eff = line_amplitude * sigma / sigma_eff

    # find range of pixels where line is non-negligible
    lo = np.searchsorted(log_obs_bin_centers,
                         log_shifted_line - MAX_SDEV * sigma_eff,
                         side="left")
    hi = np.searchsorted(log_obs_bin_centers,
                         log_shifted_line + MAX_SDEV * sigma_eff,
                         side="right")

    if hi == 0 or lo == len(log_obs_bin_centers):
        return (0, 0)

    inv_2_sigma_eff_sq = 0.5 / sigma_eff**2

    for i in range(lo, hi):
        x = log_obs_bin_centers[i] - log_shifted_line
        vals[i - lo] = amp_eff * np.exp(-x**2 * inv_2_sigma_eff_sq)

    return (lo, hi)
