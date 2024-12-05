"""
fastspecfit.linetable
=====================

Emission line table, read from a file.

"""
import os
from astropy.table import Table

from fastspecfit.logger import log

class LineTable(object):

    def __init__(self, emlines_file=None):

        """Read the set of emission lines of interest.

        """
        if emlines_file is None:
            # use default emlines location
            from importlib import resources
            emlines_file = resources.files('fastspecfit').joinpath('data/emlines.ecsv')

        if not os.path.isfile(emlines_file):
            errmsg = f'Emission lines parameter file {emlines_file} does not exist.'
            log.critical(errmsg)
            raise FileNotFoundError(errmsg)

        self.file = emlines_file

        try:
            linetable = Table.read(emlines_file, format='ascii.ecsv', guess=False)
        except:
            errmsg = f'Problem reading emission lines parameter file {emlines_file}.'
            log.critical(errmsg)
            raise IOError(errmsg)

        linetable.sort('restwave')

        self.table = linetable
