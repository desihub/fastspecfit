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
from fastspecfit.logger import log

class Singletons(object):

    def __init__(self):
        pass

    def initialize(self,
                   emlines_file=None,
                   fphotofile=None,
                   fastphot=False,
                   stackfit=False,
                   ignore_photometry=False,
                   template_file=None,
                   template_version=None,
                   template_imf=None):

        # templates for continnuum fitting
        self.templates = Templates(template_file=template_file,
                                   template_version=template_version,
                                   imf=template_imf,
                                   mintemplatewave=None,
                                   maxtemplatewave=40e4,
                                   fastphot=fastphot)
        log.info(f'Cached stellar templates {self.templates.file}')

        # emission line table
        self.emlines = LineTable(emlines_file)
        log.info(f'Cached emission-line table {self.emlines.file}')

        # photometry
        self.photometry = Photometry(fphotofile,
                                     stackfit,
                                     ignore_photometry)
        log.info(f'Cached photometric filters and parameters {self.photometry.fphotofile}')

        # fiducial cosmology
        self.cosmology = TabulatedDESI()
        log.info(f'Cached cosmology table {self.cosmology.file}')

        # IGM model
        self.igm = Inoue14()
        log.info('Cached IGM attenuation parameters.')


# global structure with single-copy data, initially empty
sc_data = Singletons()
