"""
fastspecfit.util
================

General utilities.

"""
import numpy as np

try: # this fails when building the documentation
    from scipy import constants
    C_LIGHT = constants.c / 1000.0 # [km/s]
except:
    C_LIGHT = 299792.458

def ivar2var(ivar, sigma=False):
    """Safely convert an inverse variance to a variance."""
    var = np.zeros_like(ivar)
    goodmask = ivar > 0 # True is good
    if np.count_nonzero(goodmask) == 0:
        log.warning('All values are masked!')
        raise ValueError
    var[goodmask] = 1 / ivar[goodmask]
    if sigma:
        var = np.sqrt(var) # return a sigma
    return var, goodmask

