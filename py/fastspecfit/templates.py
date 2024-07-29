import os
import numpy as np

import fitsio
from astropy.table import Table

from fastspecfit.logger import log

class Templates(object):

    # SPS template constants (used by build-templates)
    PIXKMS_BLU = 25.  # [km/s]
    PIXKMS_RED = 100. # [km/s]
    PIXKMS_WAVESPLIT = 1e4 # [Angstrom]

    FTEMPLATES_DIR_NERSC = '/global/cfs/cdirs/desi/science/gqp/templates/fastspecfit'

    DEFAULT_TEMPLATEVERSION = '2.0.0'
    DEFAULT_IMF = 'chabrier'

    def __init__(self, template_file=None, template_version=None, imf=None,
                 mintemplatewave=None, maxtemplatewave=40e4, vdisp_nominal=125.,
                 fastphot=False, read_linefluxes=False):
        """"
        Read the templates into a dictionary.
    
        Parameters
        ----------
        template_file : :class:`str`
            Full path to the templates to read. Defaults to the output of
            :class:`get_templates_filename`.
        template_version : :class:`str`
            Version of the templates.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.
        mintemplatewave : :class:`float`, optional, defaults to None
            Minimum template wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxtemplatewave : :class:`float`, optional, defaults to 6e4
            Maximum template wavelength to read into memory.

        """
        if template_file is None:
            if template_version is None:
                template_version = Templates.DEFAULT_TEMPLATEVERSION
            if imf is None:
                imf = Templates.DEFAULT_IMF
            template_file = self.get_templates_filename(template_version=template_version, imf=imf)
        
        if not os.path.isfile(template_file):
            errmsg = f'Templates file {template_file} not found.'
            log.critical(errmsg)
            raise IOError(errmsg)

        self.file = template_file
        
        T = fitsio.FITS(template_file)
        templatewave     = T['WAVE'].read()        # [npix]
        wavehdr          = T['WAVE'].read_header() # [npix]
        templateflux     = T['FLUX'].read()        # [npix x nsed]
        templatelineflux = T['LINEFLUX'].read()    # [npix x nsed]
        templateinfo     = T['METADATA'].read()
        templatehdr      = T['METADATA'].read_header()
        vdisphdr         = T['VDISPFLUX'].read_header()

        self.imf = templatehdr['IMF']
        self.ntemplates = len(templateinfo)
        
        if mintemplatewave is None:
            mintemplatewave = np.min(templatewave)
        
        keeplo = np.searchsorted(templatewave, mintemplatewave, 'left')
        keephi = np.searchsorted(templatewave, maxtemplatewave, 'right')
        self.wave = templatewave[keeplo:keephi]
        self.flux = templateflux[keeplo:keephi, :]
        self.flux_nolines = self.flux - templatelineflux[keeplo:keephi, :]

        # Cache a copy of the line-free templates at the nominal velocity
        # dispersion (needed for fastphot as well).
        if 'VDISPNOM' in vdisphdr: # older templates do not have this header card
            vdisp_nominal = vdisphdr['VDISPNOM'] # [km/s]
        self.vdisp_nominal = vdisp_nominal

        hi = np.searchsorted(self.wave, Templates.PIXKMS_WAVESPLIT, 'left')
        self.pixpos_wavesplit = hi
        
        # dust attenuation curve
        self.dust_klambda = Templates.klambda(self.wave)
        
        self.flux_nomvdisp = self.convolve_vdisp(
            self.flux, vdisp_nominal, limit=hi,
            pixsize_kms=Templates.PIXKMS_BLU)
        
        self.flux_nolines_nomvdisp = self.convolve_vdisp(
            self.flux_nolines, vdisp_nominal, limit=hi,
            pixsize_kms=Templates.PIXKMS_BLU)
        
        self.info = Table(templateinfo)

        # maintain backwards compatibility with older templates (versions <2.0.0)
        self.use_legacy_fitting = ('VDISPMIN' in vdisphdr)
        
        if self.use_legacy_fitting:
            if not fastphot:
                vdispwave = T['VDISPWAVE'].read()
                vdispflux = T['VDISPFLUX'].read() # [nvdisppix,nvdispsed,nvdisp]

                # see bin/build-fsps-templates
                if vdisphdr['VDISPRES'] > 0.:
                    nvdisp = int(np.ceil((vdisphdr['VDISPMAX'] - vdisphdr['VDISPMIN']) / vdisphdr['VDISPRES'])) + 1
                    vdisp = np.linspace(vdisphdr['VDISPMIN'], vdisphdr['VDISPMAX'], nvdisp)
                else:
                    vdisp = np.atleast_1d(vdisphdr['VDISPMIN'])
        
                if not vdisp_nominal in vdisp:
                    errmsg = 'Nominal velocity dispersion is not in velocity dispersion vector.'
                    log.critical(errmsg)
                    raise ValueError(errmsg)

                self.vdisp = vdisp
                self.vdisp_nominal_index = np.where(vdisp == vdisp_nominal)[0] 
                self.vdispflux = vdispflux
                self.vdispwave = vdispwave

            # Read the model emission-line fluxes; only present for
            # template_version>=1.1.1 and generally only useful to a power-user.
            if read_linefluxes:
                self.linewaves  = T['LINEWAVES'].read()
                self.linefluxes = T['LINEFLUXES'].read()

        else:
            if 'DUSTFLUX' in T and 'AGNFLUX' in T:
                dustflux = T['DUSTFLUX'].read()
                
                # make sure fluxes are normalized to unity
                dustflux /= np.trapz(dustflux, x=templatewave) # should already be 1.0
                self.dustflux = dustflux[keeplo:keephi]
                
                #dusthdr = T['DUSTFLUX'].read_header()
                #self.qpah     = dusthdr['QPAH']
                #self.umin     = dusthdr['umin']
                #self.gamma    = dusthdr['GAMMA']

                #agnflux = T['AGNFLUX'].read()
                #agnflux  /= np.trapz(agnflux, x=templatewave) # should already be 1.0
                #self.agnflux  = agnflux[keeplo:keephi]
                
                #agnhdr = T['AGNFLUX'].read_header()
                #self.agntau   = agnhdr['AGNTAU']
            else:
                errmsg = f'Templates file {template_file} missing mandatory extensions DUSTFLUX and AGNFLUX.'
                log.critical(errmsg)
                raise IOError(errmsg)
    
    
    @staticmethod
    def get_templates_filename(template_version, imf):
        """Get the templates filename. """
        template_dir = os.path.expandvars(os.environ.get('FTEMPLATES_DIR', Templates.FTEMPLATES_DIR_NERSC))
        template_file = os.path.join(template_dir, template_version, f'ftemplates-{imf}-{template_version}.fits')
        return template_file

    
    @staticmethod
    def convolve_vdisp(templateflux, vdisp, pixsize_kms=None, limit=None):
        """Convolve an input spectrum to a desired velocity dispersion in km/s.
        
        Parameters
        ----------
        templateflux : :class:`numpy.ndarray` [npix, nmodel]
            One- or two-dimensional input model spectra.
        vdisp : :class:`float`
            Desired velocity dispersion.
        pixsize_kms : :class:`float`
            Pixel size of `templateflux` in km/s.
        limit : :class:`int`
            Only smooth up to the pixel position (in `templateflux`) specified by
            this parameter.

        Returns
        -------
        :class:`numpy.ndarray` [npix, nmodel]
            Gaussian-smoothed model spectra
        """
        
        from scipy.ndimage import gaussian_filter1d

        if pixsize_kms is None:
            pixsize_kms = Templates.PIXKMS_BLUE
            
        # Convolve by the velocity dispersion.
        if vdisp <= 0.:
            output = templateflux.copy()
        else:
            output = np.empty_like(templateflux)
            sigma = vdisp / pixsize_kms # [pixels]
            
            if limit is None:
                limit = templateflux.shape[0]
                
            if templateflux.ndim == 1:
                gaussian_filter1d(templateflux[:limit], sigma=sigma, output=output[:limit])
                output[limit:] = templateflux[limit:]
            else:
                gaussian_filter1d(templateflux[:limit, :], sigma=sigma, axis=0,
                                  output=output[:limit:, :])
                output[limit:, :] = templateflux[limit:, :]

        return output

    @staticmethod
    def klambda(wave):
        """Construct the total-to-selective attenuation curve, k(lambda).

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
            Input rest-frame wavelength array in Angstrom.

        Returns
        -------
        :class:`numpy.ndarray` [npix]
            Total-to-selective attenuation curve, k(lambda).

        Notes
        -----
        ToDo: support more sophisticated dust models.

        """
        dust_power = -0.7     # power-law slope
        dust_normwave = 5500. # pivot wavelength

        return (wave / dust_normwave)**dust_power 
