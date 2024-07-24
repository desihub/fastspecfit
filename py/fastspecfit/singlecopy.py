#
# Single-copy (per process) data structures read from file
#

import os

from astropy.table import Table

from fastspecfit.cosmo      import TabulatedDESI
from fastspecfit.inoue14    import Inoue14
from fastspecfit.photometry import Photometry

class Singletons(object):

    def __init__(self):
        pass
    
    def setup(self,
              emlines_file=None,
              fphotofile=None,
              stackfit=False,
              ignore_photometry=False):
        
        # IGM model
        self.igm = Inoue14()

        # fiducial cosmology
        self.cosmology = TabulatedDESI()
        
        # emission line table
        self.emline_table =  read_emlines(emlines_file=emlines_file)
        
        # photometry
        self.photometry = Photometry(fphotofile,
                                     stackfit=False,
                                     ignore_photometry=False)

        
# global structure with single-copy data, initially empty
sc_data = Singletons()

# initialization entry point
def initialize_sc_data(*args):
    sc_data.setup(*args)
    
        
def read_emlines(emlines_file=None):
    """Read the set of emission lines of interest.

    """
    if emlines_file is None:
        # use default emlines location
        from importlib import resources
        emlines_file = resources.files('fastspecfit').joinpath('data/emlines.ecsv')
        
    if not os.path.isfile(emlines_file):
        errmsg = f'Emission lines parameter file {emlines_file} does not exist.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    try:
        linetable = Table.read(emlines_file, format='ascii.ecsv', guess=False)
    except: 
        errmsg = f'Problem reading emission lines parameter file {emlines_file}.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    linetable.sort('restwave')
            
    return linetable
