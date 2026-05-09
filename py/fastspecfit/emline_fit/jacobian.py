import numpy as np
from numba import jit

from fastspecfit.util import C_LIGHT
from fastspecfit.resolution import SIGMA0_ANGSTROM
from .utils import MAX_SDEV, max_buffer_width


@jit(nopython=True, nogil=True, cache=True)
def emline_model_jacobian(line_parameters,
                          log_obs_bin_edges,
                          redshift,
                          line_wavelengths,
                          padding):
    """
    Compute Jacobian of the function computed in emline_model().

    Parameters
    ----------
    line_parameters : :class:`np.ndarray`
      Parameters of all fitted lines.
    log_obs_bin_edges : :class:`np.ndarray` [# wavelength bins + 1]
      Natural logs of observed wavelength bin edges.
    redshift : :class:`np.float64`
      Redshift of observed spectrum.
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    padding : :class:`int`
      Number of entries to be added to size of each sparse row
      allocated for output Jacobian, for later use.

    Returns
    -------
    :class:`tuple` (endpts, dd)
      Sparse Jacobian, where column j has nonzero values in interval
      [ endpts[j,0] , endpts[j,1] ], which are stored in dd[j].

    """
    line_amplitudes, line_vshifts, line_sigmas = np.split(line_parameters, 3)

    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, padding)

    nlines = len(line_wavelengths)
    dd     = np.empty((3 * nlines, max_width), dtype=line_amplitudes.dtype)
    endpts = np.zeros((3 * nlines,         2), dtype=np.int32)

    starts = endpts[:, 0]
    ends   = endpts[:, 1]

    for j in range(nlines):

        if line_amplitudes[j] == 0.:
            continue

        sigma_line = line_sigmas[j] / C_LIGHT

        line_shift       = 1. + redshift + line_vshifts[j] / C_LIGHT
        shifted_line     = line_wavelengths[j] * line_shift
        log_shifted_line = np.log(shifted_line)

        sigma0   = SIGMA0_ANGSTROM / shifted_line
        sigma_eff = np.sqrt(sigma_line * sigma_line + sigma0 * sigma0)

        inv_sigma_eff_sq = 1.0 / (sigma_eff * sigma_eff)

        # support bounds using sigma_eff (wider than sigma_line for narrow lines)
        lo = np.searchsorted(log_obs_bin_edges,
                             log_shifted_line - MAX_SDEV * sigma_eff,
                             side="left")
        hi = np.searchsorted(log_obs_bin_edges,
                             log_shifted_line + MAX_SDEV * sigma_eff,
                             side="right")

        if hi == 0 or lo == len(log_obs_bin_edges):
            continue

        dlo = 1 if lo == 0                      else 0
        dhi = 1 if hi == len(log_obs_bin_edges) else 0

        s = lo - 1 + dlo
        e = hi - dhi

        dda_vals = dd[           j]
        ddv_vals = dd[nlines   + j]
        dds_vals = dd[2*nlines + j]

        # factors shared across bins
        # ∂f/∂amplitude = (sigma_line/sigma_eff) * G_j
        r = sigma_line / sigma_eff

        # ∂f/∂vshift = A_eff * G_j * x_j * inv_sigma_eff_sq / (C_LIGHT * line_shift)
        A_eff     = line_amplitudes[j] * r
        ddv_const = A_eff * inv_sigma_eff_sq / (C_LIGHT * line_shift)

        # ∂f/∂sigma = amp * G_j / (C_LIGHT * sigma_eff) * inv_sigma_eff_sq
        #             * (sigma0² + sigma_line² * x_j² * inv_sigma_eff_sq)
        sigma0_sq      = sigma0 * sigma0
        sigma_line_sq  = sigma_line * sigma_line
        dds_const      = line_amplitudes[j] / (C_LIGHT * sigma_eff) * inv_sigma_eff_sq

        for i in range(s, e):
            log_lambda_center = 0.5 * (log_obs_bin_edges[i] + log_obs_bin_edges[i + 1])
            x = log_lambda_center - log_shifted_line
            G = np.exp(-0.5 * x * x * inv_sigma_eff_sq)

            dda_vals[i - s] = r * G
            ddv_vals[i - s] = ddv_const * x * G
            dds_vals[i - s] = dds_const * (sigma0_sq + sigma_line_sq * x * x * inv_sigma_eff_sq) * G

        starts[j] = s
        ends[j]   = e

    # replicate first third of endpts twice more, since same support
    # applies to all three parameters of each line
    for i in range(1, 3):
        endpts[i*nlines:(i+1)*nlines, :] = endpts[:nlines, :]

    # for lines with zero amplitude, vshift and sigma derivatives are zero
    for i, amp in enumerate(line_amplitudes):
        if amp == 0.:
            endpts[i + nlines,   :] = 0
            endpts[i + 2*nlines, :] = 0

    return (endpts, dd)


@staticmethod
@jit(nopython=True, nogil=True, cache=True)
def patch_jacobian(obs_bin_centers,
                   obs_weights,
                   patch_endpts,
                   patch_pivotwave):
    """
    Compute partial Jacobian associated with just the patch parameters
    in the sparse column form used in sparse_rep.py.  This Jacobian
    is independent of both the line model and the particular choices
    of slope/intercept for each patch, so can be computed once when
    we set up the optimization.

    Parameters
    ----------
    obs_bin_centers : :class:`np.ndarray` [# wavelength bins]
      Center of each observed wavelength bin.
    obs_weights : :class:`np.ndarray` [# wavelength bins]
      Weights for each observed wavelength bin.
    patch_endpts : :class:`np.ndarray` [# patches x 2]
      Endpoints of each patch in wavelength array.
    patch_pivotwave : :class:`np.ndarray` [# patches]
      Wavelength offset for fitted affine params of each patch .

    Returns
    -------
    :class:`tuple` (endpts, M)
      Sparse Jacobian, where column j has nonzero values in interval
      [ endpts[j,0] , endpts[j,1] ], which are stored in M[j].

    """

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
