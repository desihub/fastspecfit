import numpy as np
from numba import jit

from fastspecfit.util import C_LIGHT
from fastspecfit.resolution import SIGMA0_ANGSTROM  # single source of truth
from .utils import MAX_SDEV, max_buffer_width


@jit(nopython=True, nogil=True, cache=True)
def emline_model(line_wavelengths,
                 line_parameters,
                 log_obs_bin_edges,
                 redshift):
    """
    Given a fixed set of spectral lines and known redshift, and estimates
    for the amplitude, width, and velocity shift of each line, compute the
    flux of the pre-convolved model F_sigma evaluated at the center of each
    observed spectral bin.  The result is intended to be multiplied by the
    deconvolved resolution matrix W (not the raw R).

    Parameters
    ----------
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    line_parameters : :class:`np.ndarray`
      Parameters of each fitted line.
    log_obs_bin_edges : :class:`np.ndarray` [# obs wavelength bins + 1]
      Natural logs of observed wavelength bin edges.
    redshift : :class:`np.float64`
      Redshift of observed spectrum.

    Returns
    -------
    :class:`np.ndarray` of length [# obs wavelength bins] with
    F_sigma evaluated at each bin center.

    """
    line_amplitudes, line_vshifts, line_sigmas = np.split(line_parameters, 3)

    nbins = len(log_obs_bin_edges) - 1
    model_fluxes = np.zeros(nbins, dtype=line_amplitudes.dtype)
    bin_vals = np.empty(nbins + 2, dtype=line_amplitudes.dtype)

    for j in range(len(line_wavelengths)):
        if line_amplitudes[j] == 0.:
            continue
        s, e = emline_model_core(line_wavelengths[j],
                                 line_amplitudes[j],
                                 line_vshifts[j],
                                 line_sigmas[j],
                                 log_obs_bin_edges,
                                 redshift,
                                 bin_vals)
        model_fluxes[s:e] += bin_vals[:e-s]

    return model_fluxes


@jit(nopython=True, nogil=True, cache=True)
def emline_perline_models(line_wavelengths,
                          line_parameters,
                          log_obs_bin_edges,
                          redshift,
                          padding):
    """
    Given parameters for a set of lines, compute for each line
    individually its waveform (F_sigma at bin centers) in the range of
    bins described by log_obs_bin_edges.

    Parameters
    ----------
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    line_parameters : :class:`np.ndarray`
      Parameters of each fitted line.
    log_obs_bin_edges : :class:`np.ndarray` [# obs wavelength bins + 1]
      Natural logs of observed wavelength bin edges.
    redshift : :class:`np.float64`
      Redshift of observed spectrum.
    padding : :class:`int`
      Number of entries by which to pad each sparse row for later use.

    Returns
    -------
    Row-sparse matrix of size [# lines x # obs wavelength bins]
    as tuple (endpts, profiles).

    """
    line_amplitudes, line_vshifts, line_sigmas = np.split(line_parameters, 3)

    nbins = len(log_obs_bin_edges) - 1
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, SIGMA0_ANGSTROM, padding)

    nlines = len(line_wavelengths)
    line_profiles = np.empty((nlines, max_width), dtype=line_amplitudes.dtype)
    endpts = np.zeros((nlines, 2), dtype=np.int32)

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
                                 bin_vals)
        endpts[j, 0] = s
        endpts[j, 1] = e

    return (endpts, line_profiles)


@jit(nopython=True, nogil=True, cache=True)
def emline_model_core(line_wavelength,
                      line_amplitude,
                      line_vshift,
                      line_sigma,
                      log_obs_bin_edges,
                      redshift,
                      vals):
    """
    Compute F_sigma(lambda_j) for one spectral line at each bin center j
    in its support.

    The model evaluates the intrinsic Gaussian pre-convolved with the
    fiducial Gaussian of width SIGMA0_ANGSTROM (matching the deconvolved
    resolution matrix W).  For an intrinsic line with amplitude A and
    log-lambda width sigma_line, the pre-convolved profile is:

        F_sigma(lambda_j) = A * (sigma_line / sigma_eff)
                              * exp(-0.5 * (log(lambda_j) - log(lambda_0))^2
                                    / sigma_eff^2)

    where sigma_eff^2 = sigma_line^2 + sigma0^2 and sigma0 is the fiducial
    width in log-lambda space.

    The output is intended for multiplication by W (not R); see resolution.py.

    Parameters
    ----------
    line_wavelength : :class:`np.float64`
      Nominal rest-frame wavelength of line.
    line_amplitude : :class:`np.float64`
      Amplitude of the intrinsic Gaussian in log-lambda space.
    line_vshift : :class:`np.float64`
      Velocity shift of line [km/s].
    line_sigma : :class:`np.float64`
      Intrinsic Gaussian width of line [km/s].
    log_obs_bin_edges : :class:`np.ndarray` [# obs wavelength bins + 1]
      Natural logs of observed wavelength bin edges.
    redshift : :class:`np.float64`
      Redshift of observed spectrum.
    vals : :class:`np.ndarray` (output)
      Array in which to write nonzero F_sigma values for this line.

    Returns
    -------
    Tuple (s, e) such that all nonzero values are in bins [s, e).
    vals[0:e-s] contains these values.

    """
    # intrinsic line width in log-lambda
    sigma_line = line_sigma / C_LIGHT

    # wavelength of line center in observed frame
    line_shift   = 1. + redshift + line_vshift / C_LIGHT
    shifted_line = line_wavelength * line_shift
    log_shifted_line = np.log(shifted_line)

    # fiducial pre-convolution width in log-lambda (matches W construction)
    sigma0 = SIGMA0_ANGSTROM / shifted_line

    # effective width: quadrature sum of intrinsic and fiducial widths
    sigma_eff = np.sqrt(sigma_line * sigma_line + sigma0 * sigma0)

    # support bounds using sigma_eff (wider than sigma_line for narrow lines)
    lo = np.searchsorted(log_obs_bin_edges,
                         log_shifted_line - MAX_SDEV * sigma_eff,
                         side="left")
    hi = np.searchsorted(log_obs_bin_edges,
                         log_shifted_line + MAX_SDEV * sigma_eff,
                         side="right")

    if hi == 0 or lo == len(log_obs_bin_edges):
        return (0, 0)

    # amplitude of the pre-convolved (effective) Gaussian;
    # flux-conserving: total flux = A * sqrt(2pi) * sigma_line * shifted_line
    A_eff = line_amplitude * sigma_line / sigma_eff

    dlo = 1 if lo == 0                      else 0
    dhi = 1 if hi == len(log_obs_bin_edges) else 0

    s = lo - 1 + dlo
    e = hi - dhi

    inv_2sigma_eff_sq = 0.5 / (sigma_eff * sigma_eff)

    for i in range(s, e):
        # log-wavelength at bin center (midpoint of log edges)
        log_lambda_center = 0.5 * (log_obs_bin_edges[i] + log_obs_bin_edges[i + 1])
        dx = log_lambda_center - log_shifted_line
        vals[i - s] = A_eff * np.exp(-dx * dx * inv_2sigma_eff_sq)

    return (s, e)
