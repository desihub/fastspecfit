"""
fastspecfit.linetable
=====================

Emission line table, read from a file.

"""
import os
from astropy.table import Table

from fastspecfit.logger import log

class LineTable(object):
    """Emission line table used for spectral fitting.

    Parameters
    ----------
    emlines_file : str, optional
        Path to an ECSV emission-line parameter file. Defaults to the
        bundled ``data/emlines.ecsv``.

    Attributes
    ----------
    file : str
        Path to the emission-line file that was read.
    table : :class:`astropy.table.Table`
        Emission-line parameters.

    """
    def __init__(self, emlines_file=None):
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
