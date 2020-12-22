"""
desigal.fastspecfit
===================

"""
import pdb # for debugging

import os, time
import numpy as np
import multiprocessing

from scipy.ndimage import gaussian_filter1d
import astropy.units as u
from astropy.table import Table, Column, vstack, join, hstack
from astropy.modeling import Fittable1DModel

from fnnls import fnnls
#from desispec.interpolation import resample_flux
from redrock.rebin import trapz_rebin

from scipy import constants
C_LIGHT = constants.c / 1000.0 # [km/s]

from desiutil.log import get_logger
log = get_logger()

