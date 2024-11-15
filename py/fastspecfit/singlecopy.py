"""
fastspecfit.singlecopy
======================

Single-copy (per process) data structures read from files.

"""
from fastspecfit.cosmo import TabulatedDESI
from fastspecfit.igm import Inoue14
from fastspecfit.photometry import Photometry
from fastspecfit.linetable import LineTable
from fastspecfit.templates import Templates
from fastspecfit.logger import log, DEBUG

class Singletons(object):

    def __init__(self):
        pass

    def initialize(self,
                   emlines_file=None,
                   fphotofile=None,
                   fastphot=False,
                   fitstack=False,
                   ignore_photometry=False,
                   template_file=None,
                   template_version=None,
                   template_imf=None,
                   log_verbose=False,
    ):

        # adjust logging level if requested
        if log_verbose:
            log.setLevel(DEBUG)

        # templates for continnuum fitting
        self.templates = Templates(template_file=template_file,
                                   template_version=template_version,
                                   imf=template_imf,
                                   mintemplatewave=None,
                                   maxtemplatewave=40e4,
                                   fastphot=fastphot)
        log.debug(f'Cached stellar templates {self.templates.file}')

        # emission line table
        self.emlines = LineTable(emlines_file)
        log.debug(f'Cached emission-line table {self.emlines.file}')

        # photometry
        self.photometry = Photometry(fphotofile, fitstack,
                                     ignore_photometry)
        log.debug(f'Cached photometric filters and parameters {self.photometry.fphotofile}')

        # fiducial cosmology
        self.cosmology = TabulatedDESI()
        log.debug(f'Cached cosmology table {self.cosmology.file}')

        # IGM model
        self.igm = Inoue14()
        log.debug(f'Cached {self.igm.reference} IGM attenuation parameters.')


# global structure with single-copy data, initially empty
sc_data = Singletons()
