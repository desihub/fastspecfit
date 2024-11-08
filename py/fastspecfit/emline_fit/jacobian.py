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
    """
    Compute Jacobian of the function computed in emline_model().

    Parameters
    ----------
    line_parameters : :class:`np.ndarray`
      Parameters of all fitted lines.
    log_obs_bin_edges : :class:`np.ndarray` [# wavelength bins + 1]
      Natural logs of observed wavelength bin edges.
    ibin_widths : :class:np.ndarray` [# wavelength bins]
      Inverse of width of each observed wavelength bin.
    redshift : :class:`np.float64`
      Red shift of observed spectrum.
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    padding : :class:`int`:
      Number of entries to be added to size of each sparse row
      allocated for output Jacobian, for later use.

    Returns
    -------
    :class:`tuple` (endpts, dd)
      Sparse Jacobian, where column j has nonzero values in interval
      [ endpts[j,0] , endpts[j,1] ], which are stored in dd[j].

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
