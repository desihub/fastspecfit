"""
fastspecfit.singlecopy
======================

Single-copy (per process) data structures read from files.

"""
from fastspecfit.cosmo import TabulatedDESI
from fastspecfit.igm import Inoue14
from fastspecfit.photometry import Photometry
from fastspecfit.linetable import LineTable
from fastspecfit.templates import Templates, VDISP_NOMINAL, VDISP_BOUNDS
from fastspecfit.logger import log, DEBUG

class Singletons(object):
    """Container for per-process singleton data structures.

    Holds global shared objects (templates, emission lines, photometry,
    cosmology, IGM model) that are read from disk once at startup and
    shared across all worker threads/processes within a single MPI rank.

    """
    def __init__(self):
        pass

    def initialize(self,
                   emlines_file=None,
                   fphotofile=None,
                   fastphot=False,
                   vdisp_nominal=VDISP_NOMINAL,
                   vdisp_bounds=VDISP_BOUNDS,
                   fitstack=False,
                   ignore_photometry=False,
                   template_file=None,
                   template_version=None,
                   template_imf=None,
                   log_verbose=False,
    ):
        """Load all singleton data structures from disk.

        Parameters
        ----------
        emlines_file : :class:`str` or None, optional
            Path to the emission-line parameter file; uses the bundled
            default when ``None``.
        fphotofile : :class:`str` or None, optional
            Path to the photometric configuration YAML file; uses the
            bundled DR9 default when ``None``.
        fastphot : :class:`bool`, optional
            If ``True``, load templates in photometry-only mode. Default
            is ``False``.
        vdisp_nominal : :class:`float`, optional
            Nominal velocity dispersion in km/s used to pre-cache FFTs.
        vdisp_bounds : tuple of float, optional
            ``(min, max)`` velocity dispersion bounds in km/s.
        fitstack : :class:`bool`, optional
            If ``True``, use the stacked-spectra photometry configuration.
            Default is ``False``.
        ignore_photometry : :class:`bool`, optional
            If ``True``, disable photometric fitting. Default is ``False``.
        template_file : :class:`str` or None, optional
            Full path to the SPS template FITS file; auto-detected when
            ``None``.
        template_version : :class:`str` or None, optional
            Template version string; used when ``template_file`` is
            ``None``.
        template_imf : :class:`str` or None, optional
            Initial mass function name for template selection.
        log_verbose : :class:`bool`, optional
            If ``True``, set the logger level to ``DEBUG``. Default is
            ``False``.

        """
        if log_verbose:
            log.setLevel(DEBUG)

        key = (emlines_file, fphotofile, fastphot, fitstack, ignore_photometry,
               template_file, template_version, template_imf,
               vdisp_nominal, tuple(vdisp_bounds) if vdisp_bounds is not None else None)
        if getattr(self, '_init_key', None) == key:
            return
        self._init_key = key

        # templates for continuum fitting
        self.templates = Templates(template_file=template_file,
                                   template_version=template_version,
                                   imf=template_imf,
                                   mintemplatewave=None,
                                   maxtemplatewave=40e4,
                                   vdisp_nominal=vdisp_nominal,
                                   vdisp_bounds=vdisp_bounds,
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
