import numpy as np

from numba import jit

from fastspecfit.util import C_LIGHT

from .utils import (
    MAX_SDEV,
    max_buffer_width,
)


@jit(nopython=True, nogil=True, cache=True)
def emline_model_jacobian(line_parameters,
                          log_obs_bin_centers,
                          redshift,
                          line_wavelengths,
                          sigma_G,
                          padding):
    """Compute the sparse Jacobian of :func:`emline_model`.

    Parameters
    ----------
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    log_obs_bin_centers : :class:`numpy.ndarray`
        Natural logs of observed wavelength pixel centers.
    redshift : float
        Redshift of the observed spectrum.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    sigma_G : float
        Gaussian pre-convolution width in Angstroms.
    padding : int
        Extra entries to pad each sparse row for later use.

    Returns
    -------
    endpts : :class:`numpy.ndarray` of int, shape (3*nlines, 2)
        Start and end pixel indices of the nonzero range for each parameter.
    dd : :class:`numpy.ndarray`, shape (3*nlines, max_width)
        Partial derivative values within each parameter's nonzero pixel range.

    """

    line_amplitudes, line_vshifts, line_sigmas = \
        np.split(line_parameters, 3)

    max_width = max_buffer_width(log_obs_bin_centers, line_sigmas, padding, sigma_G)

    nlines = len(line_wavelengths)
    dd     = np.empty((3 * nlines, max_width), dtype=line_amplitudes.dtype)
    endpts = np.zeros((3 * nlines,         2), dtype=np.int32)

    starts = endpts[:, 0]
    ends   = endpts[:, 1]

    for j in range(nlines):

        amp = line_amplitudes[j]
        if amp == 0.:
            continue

        sigma = line_sigmas[j] / C_LIGHT

        line_shift       = 1. + redshift + line_vshifts[j] / C_LIGHT
        shifted_line     = line_wavelengths[j] * line_shift
        log_shifted_line = np.log(shifted_line)

        sigma_G_log        = sigma_G / shifted_line
        sigma_eff          = np.sqrt(sigma**2 + sigma_G_log**2)
        amp_eff            = amp * sigma / sigma_eff
        amp_over_sigma_eff = amp / sigma_eff  # finite even when sigma=0

        lo = np.searchsorted(log_obs_bin_centers,
                             log_shifted_line - MAX_SDEV * sigma_eff,
                             side="left")
        hi = np.searchsorted(log_obs_bin_centers,
                             log_shifted_line + MAX_SDEV * sigma_eff,
                             side="right")

        if hi == 0 or lo == len(log_obs_bin_centers):
            continue

        dda_vals = dd[           j]
        ddv_vals = dd[nlines   + j]
        dds_vals = dd[2*nlines + j]

        inv_2_sigma_eff_sq = 0.5 / sigma_eff**2
        inv_sigma_eff_sq   = 1.  / sigma_eff**2

        for i in range(lo, hi):
            x     = log_obs_bin_centers[i] - log_shifted_line
            gauss = np.exp(-x**2 * inv_2_sigma_eff_sq)
            g     = amp_eff * gauss

            # dF/d(amplitude)
            dda_vals[i - lo] = g / amp

            # dF/d(vshift): d(log_shifted_line)/d(vshift) = 1/C_LIGHT,
            # so dx/d(vshift) = -1/C_LIGHT
            ddv_vals[i - lo] = g * x * inv_sigma_eff_sq / C_LIGHT

            # dF/d(line_sigma): written as (amp/sigma_eff)*exp*(...) rather than
            # g*(1/sigma - ...) to avoid 1/sigma -> inf when line_sigmas[j]=0;
            # the two forms are algebraically identical for sigma > 0
            dds_vals[i - lo] = amp_over_sigma_eff * gauss / C_LIGHT * \
                (1. - sigma**2 * (1. - x**2 * inv_sigma_eff_sq) * inv_sigma_eff_sq)

        starts[j] = lo
        ends[j]   = hi

    # replicate first third of endpts twice more, since the same pixel
    # range applies to all three parameters of each line
    for i in range(1, 3):
        endpts[i*nlines:(i+1)*nlines, :] = endpts[:nlines, :]

    # for lines with zero amplitude, derivatives w.r.t. vshift and sigma are zero
    for i, amp in enumerate(line_amplitudes):
        if amp == 0.:
            endpts[i + nlines,   :] = 0
            endpts[i + 2*nlines, :] = 0

    # for lines with zero width, derivatives w.r.t. amplitude are zero
    for i, sig in enumerate(line_sigmas):
        if sig == 0.:
            endpts[i, :] = 0

    return (endpts, dd)


@staticmethod
@jit(nopython=True, nogil=True, cache=True)
def patch_jacobian(obs_bin_centers,
                   obs_weights,
                   patch_endpts,
                   patch_pivotwave):
    """Compute the patch pedestal Jacobian in sparse column form.

    This Jacobian depends only on the patch geometry (not the line model
    or slope/intercept values), so it can be computed once at setup time.

    Parameters
    ----------
    obs_bin_centers : :class:`numpy.ndarray`
        Center wavelength of each observed bin in Angstroms.
    obs_weights : :class:`numpy.ndarray`
        Weights for each observed wavelength bin.
    patch_endpts : :class:`numpy.ndarray` of int, shape (npatches, 2)
        Start and end bin indices for each patch.
    patch_pivotwave : :class:`numpy.ndarray`
        Pivot wavelength for the affine pedestal of each patch.

    Returns
    -------
    endpts : :class:`numpy.ndarray` of int, shape (2*npatches, 2)
        Start and end bin indices for slope and intercept columns.
    M : :class:`numpy.ndarray`, shape (2*npatches, max_width)
        Weighted partial derivatives with respect to slope and intercept
        for each patch.

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
