import numpy as np

from numba import jit

from fastspecfit.util import C_LIGHT


# Do not bother computing the Gaussian outside this many standard deviations.
MAX_SDEV = 5.


@jit(nopython=True, nogil=True, cache=True)
def max_buffer_width(log_obs_bin_centers, line_sigmas, padding=0, sigma_G_angstrom=0.):
    """Estimate the maximum number of nonzero bins possible for any spectral line.

    Parameters
    ----------
    log_obs_bin_centers : :class:`np.ndarray`
        Log wavelengths of all observed bin centers.
    line_sigmas : :class:`np.ndarray`
        Gaussian widths of all spectral lines in km/s.
    padding : int, optional
        Extra bins to add on each side for future expansion. Defaults to 0.
    sigma_G_angstrom : float, optional
        Gaussian pre-convolution width in Angstroms. Defaults to 0.

    Returns
    -------
    int
        Safe upper bound on the number of bins with nonzero flux for any line.

    """

    # Effective max sigma in log-lambda, including Gaussian pre-convolution floor.
    sigma_max = np.max(line_sigmas) / C_LIGHT
    lambda_min = np.exp(log_obs_bin_centers[0])
    sigma_G_log = sigma_G_angstrom / lambda_min
    sigma_eff_max = np.sqrt(sigma_max**2 + sigma_G_log**2)

    max_width = \
        int(2*MAX_SDEV*sigma_eff_max / \
            np.min(np.diff(log_obs_bin_centers))) + \
            2*padding + 4
    return max_width
