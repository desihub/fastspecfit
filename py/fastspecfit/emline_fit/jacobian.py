import numpy as np

from numba import jit

from fastspecfit.util import C_LIGHT

from .utils import (
    MAX_SDEV,
    max_buffer_width,
    norm_cdf,
    norm_pdf
)


@jit(nopython=True, nogil=True, cache=True)
def emline_model_jacobian(line_parameters,
                          log_obs_bin_edges,
                          ibin_widths,
                          redshift,
                          line_wavelengths,
                          padding):
    """Compute the sparse Jacobian of :func:`emline_model`.

    Parameters
    ----------
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    log_obs_bin_edges : :class:`numpy.ndarray`
        Natural logs of observed wavelength bin edges.
    ibin_widths : :class:`numpy.ndarray`
        Inverse widths of each observed wavelength bin.
    redshift : float
        Redshift of the observed spectrum.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    padding : int
        Extra entries to pad each sparse row for later use.

    Returns
    -------
    endpts : :class:`numpy.ndarray` of int, shape (3*nlines, 2)
        Start and end bin indices of the nonzero range for each parameter.
    dd : :class:`numpy.ndarray`, shape (3*nlines, max_width)
        Partial derivative values within each parameter's nonzero bin range.

    """

    SQRT_2PI = np.sqrt(2*np.pi)

    nbins = len(log_obs_bin_edges) - 1

    line_amplitudes, line_vshifts, line_sigmas = \
        np.split(line_parameters, 3)

    # buffers for per-parameter calculations, sized large
    # enough for max possible range [s .. e), plus extra
    # padding as directed by caller.
    max_width = max_buffer_width(log_obs_bin_edges, line_sigmas, padding)

    nlines = len(line_wavelengths)
    dd     = np.empty((3 * nlines, max_width), dtype=line_amplitudes.dtype)
    endpts = np.zeros((3 * nlines,         2), dtype=np.int32)

    starts = endpts[:, 0]
    ends   = endpts[:, 1]

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

        nedges = hi - lo + 2  # compute values at edges [lo - 1 ... hi]

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

        # Compute partial derivs for bin i+lo-1 for 0 <= i < nedges - 1.
        # But:
        #  * if lo == 0, bin lo - 1 is not defined, so skip it
        #  * if hi == |log_obs_bin_edges|, bin hi - 1 is not defined,
        #      so skip it
        # Implement skips by modifying range of computation loop.
        #
        # After loop, vals contains weighted partial derives of bins
        # lo - 1 + dlo .. hi - dhi, plus up to two cells of
        # trailing garbage.

        dlo = 1 if lo == 0                      else 0
        dhi = 1 if hi == len(log_obs_bin_edges) else 0

        for i in range(dlo, nedges - 1 - dhi):
            dda_vals[i-dlo] = (dda_vals[i+1] - dda_vals[i]) * ibin_widths[i+lo]
            ddv_vals[i-dlo] = (ddv_vals[i+1] - ddv_vals[i]) * ibin_widths[i+lo]
            dds_vals[i-dlo] = (dds_vals[i+1] - dds_vals[i]) * ibin_widths[i+lo]

        # starts[j] is set to first valid bin
        starts[j] = lo - 1 + dlo

        # ends[j] is set one past last valid bin
        ends[j]   = hi - dhi

    # replicate first third of endpts (which is what we
    # set above) twice more, since same endpts apply to
    # all three params of each line
    for i in range(1,3):
        endpts[i*nlines:(i+1)*nlines,:] = endpts[:nlines,:]

    # for lines with zero amplitude,
    # partial derivatives w/r to their vshifts
    # and sigmas are zero
    for i, amp in enumerate(line_amplitudes):
        if amp == 0.:
            endpts[i + nlines,  :] = 0
            endpts[i + 2*nlines,:] = 0

    # for lines with zero width, partial derivatives
    # w/r to their amplitudes are zero.
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
