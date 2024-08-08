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

        # IGM model
        self.igm = Inoue14()

        # fiducial cosmology
        self.cosmology = TabulatedDESI()

        # emission line table
        self.emlines = LineTable(emlines_file)

        # photometry
        self.photometry = Photometry(fphotofile,
                                     stackfit,
                                     ignore_photometry)

        # templates for continnuum fitting
        self.templates = Templates(template_file=template_file,
                                   template_version=template_version,
                                   imf=template_imf,
                                   mintemplatewave=None,
                                   maxtemplatewave=40e4,
                                   fastphot=fastphot)


# global structure with single-copy data, initially empty
sc_data = Singletons()
