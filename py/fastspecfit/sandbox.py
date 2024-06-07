"""
fastspecfit.sandbox
===================

Sandbox code.

"""
import numpy as np
import numba

from desiutil.log import get_logger
log = get_logger()

def fit_continuum(templateflux, dluminosity, redshift):
    """Fit the stellar continuum using bounded non-linear least-squares.

    templateflux - 
    
    """
    
