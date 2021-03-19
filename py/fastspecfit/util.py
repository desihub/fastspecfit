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

def air2vac(airwave):
    """http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion"""
    if airwave <= 0:
        raise ValueError('Input wavelength is not defined.')
    ss = 1e4 / airwave
    nn = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - ss**2) + 0.0001599740894897 / (38.92568793293 - ss**2)
    return airwave * nn
