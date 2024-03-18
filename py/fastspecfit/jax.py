"""
fastspecfit.jax
===============

Jax and JaxOpt fitting utilities.

Here's a sketch of my idea for computing the results of
build_emline_model() without having to sample values at a bunch of
wavelengths and then interpolate them to the desired set of bins.

* The gauss_area() function takes advantage of the inherent
vectorizability of the primitive ops to compute the contributions to a
particular bin of each Gaussian in parallel, and then sum them.

* The gauss_area_vectorized() function uses Jax's vmap mechanism to
further vectorize  gauss_area() so that you can give it arrays of bin
endpoints (rather than just two single endpoints) and get an array back.

There is some tweaking of the parameters that the real
build_emline_model() does that I didn't do, but that would be pretty
trivial to add.

"""
import jax
import jax.numpy as jnp
import numpy as np
from jax.scipy.stats.norm import cdf as normalCDF

from desiutil.log import get_logger
log = get_logger()

try: # this fails when building the documentation
    from scipy import constants
    C_LIGHT = constants.c / 1000.0 # [km/s]
except:
    C_LIGHT = 299792.458 # [km/s]

# Input: endpoints (xl, xr), arrays of Gaussian parameters A, mu, sigma
# Returns: 
def gauss_area(xl, xr, A, mu, sigma):
    """Compute the total area of all Gaussians between input endpoints `xl` and `xr`.

    Parameters
    ----------
    x1 : :class:`jax.numpy.ndarray`
        Pixel position of left endpoint(s).
    xr : :class:`jax.numpy.ndarray`
        Pixel position of right endpoint(s).
    A : :class:`jax.numpy.ndarray`
        Amplitude(s).
    mu : :class:`
        Gaussian central pixel(s).
    sigma : :class:`jax.numpy.ndarray`
        Gaussian dispersion.

    Returns
    -------

    """
    scale = jnp.sqrt(2 * jnp.pi) * sigma * A
    loCDF = normalCDF(xl, mu, sigma)
    hiCDF = normalCDF(xr, mu, sigma)
    return jnp.sum(scale * (hiCDF - loCDF))

# Input: two arrays of endpoints xl and xr, arrays of Gaussian
# parameters A, mu, sigma
# Returns: result[i] = total area of all Gaussians between xl[i], xr[i]
# gauss_area_vectorized = jax.vmap(gauss_area, in_axes=[0, 0, None, None, None])

def build_emline_model(E, A, mu, sigma):
    """Given array of bin endpoints, compute total areas between successive pairs of
    endpoints.

    """
    return gauss_area_vectorized(E[:-1], E[1:], A, mu, sigma)
