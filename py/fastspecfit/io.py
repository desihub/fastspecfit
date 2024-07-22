#!/usr/bin/env python
"""
fastspecfit.io
==============

Tools for reading DESI spectra and reading and writing fastspecfit files.

"""
import pdb # for debugging

import os, time
import numpy as np
import fitsio
from astropy.table import Table

from fastspecfit.cosmo import TabulatedDESI
from fastspecfit.photometry import Photometry
from fastspecfit.util import FLUXNORM, ZWarningMask

from desiutil.log import get_logger
log = get_logger()

# Default environment variables.
DESI_ROOT_NERSC = '/global/cfs/cdirs/desi'
DUST_DIR_NERSC = '/global/cfs/cdirs/cosmo/data/dust/v0_1'
FPHOTO_DIR_NERSC = '/global/cfs/cdirs/desi/external/legacysurvey/dr9'
FTEMPLATES_DIR_NERSC = '/global/cfs/cdirs/desi/science/gqp/templates/fastspecfit'

# list of all possible targeting bit columns
TARGETINGBITS = {
    'all': ['CMX_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',
             'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET',
             'SV2_DESI_TARGET', 'SV2_BGS_TARGET', 'SV2_MWS_TARGET',
             'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET',
             'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET'],
    'fuji': ['CMX_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',
             'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET',
             'SV2_DESI_TARGET', 'SV2_BGS_TARGET', 'SV2_MWS_TARGET',
             'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET',
             'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET'],
    'default': ['DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET'],
    }

# fibermap and exp_fibermap columns to read
FMCOLS = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'COADD_FIBERSTATUS', 'OBJTYPE',
          'PHOTSYS', 'RELEASE', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID',
          #'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
          #'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
          #'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
          #'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2'
          ]
#FMCOLS = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'COADD_FIBERSTATUS', 'OBJTYPE']

EXPFMCOLS = {
    'perexp': ['TARGETID', 'TILEID', 'FIBER', 'EXPID'],
    'pernight': ['TARGETID', 'TILEID', 'FIBER'],
    'cumulative': ['TARGETID', 'TILEID', 'FIBER'],
    'healpix': ['TARGETID', 'TILEID'], # tileid will be an array
    'custom': ['TARGETID', 'TILEID'], # tileid will be an array
    }

# redshift columns to read
REDSHIFTCOLS = ['TARGETID', 'Z', 'ZWARN', 'SPECTYPE', 'DELTACHI2']

# tsnr columns to read
TSNR2COLS = ['TSNR2_BGS', 'TSNR2_LRG', 'TSNR2_ELG', 'TSNR2_QSO', 'TSNR2_LYA']

# quasarnet afterburner columns to read
QNCOLS = ['TARGETID', 'Z_NEW', 'IS_QSO_QN_NEW_RR', 'C_LYA', 'C_CIV',
          'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']
QNLINES = ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']

DEFAULT_TEMPLATEVERSION = '2.0.0'
DEFAULT_IMF = 'chabrier'


def _unpack_one_spectrum(args):
    """Multiprocessing wrapper."""
    return unpack_one_spectrum(*args)


def unpack_one_spectrum(iobj, specdata, meta, ebv, fphoto, fastphot,
                        synthphot, ignore_photometry, debug_plots, log):
    """Unpack the data for a single object and correct for Galactic extinction. Also
    flag pixels which may be affected by emission lines.

    """
    from scipy.ndimage import binary_dilation
    from desispec.resolution import Resolution
    from fastspecfit.emline_fit import EMLine_Resolution
    from fastspecfit.util import mwdust_transmission
    from fastspecfit.continuum import ContinuumTools

    phot = Photometry(fphoto=fphoto, ignore_photometry=ignore_photometry)
    CTools = ContinuumTools(phot=phot)
    
    log.info(f'Pre-processing object {iobj} [targetid {meta[phot.uniqueid]} z={meta["Z"]:.6f}].')
    
    RV = 3.1
    meta['EBV'] = ebv

    filters = phot.filters[specdata['photsys']]
    synth_filters = phot.synth_filters[specdata['photsys']]
    if hasattr(phot, 'fiber_filters'):
        fiber_filters = phot.fiber_filters[specdata['photsys']]
    else:
        fiber_filters = None
        
    # Unpack the imaging photometry and correct for MW dust.
    if fiber_filters is not None:
        # fiber fluxes
        mw_transmission_fiberflux = np.array([mwdust_transmission(ebv, filtername) for filtername in
                                                                  phot.fiber_filters[specdata['photsys']].names])
        fibermaggies = np.zeros(len(phot.fiber_bands))
        fibertotmaggies = np.zeros(len(phot.fiber_bands))
        #ivarfibermaggies = np.zeros(len(phot.fiber_bands))
        for iband, band in enumerate(phot.fiber_bands):
            fibermaggies[iband] = meta[f'FIBERFLUX_{band.upper()}'] / mw_transmission_fiberflux[iband]
            fibertotmaggies[iband] = meta[f'FIBERTOTFLUX_{band.upper()}'] / mw_transmission_fiberflux[iband]
        
        specdata['fiberphot'] = Photometry.parse_photometry(phot.fiber_bands,
            maggies=fibermaggies, nanomaggies=True,
            lambda_eff=fiber_filters.effective_wavelengths.value, log=log)
        specdata['fibertotphot'] = Photometry.parse_photometry(phot.fiber_bands,
            maggies=fibertotmaggies, nanomaggies=True,
            lambda_eff=fiber_filters.effective_wavelengths.value, log=log)

    # total fluxes
    mw_transmission_flux = np.array([mwdust_transmission(ebv, filtername) for filtername in 
                                     phot.filters[specdata['photsys']].names])
    for band, mwdust in zip(phot.bands, mw_transmission_flux):
        meta[f'MW_TRANSMISSION_{band.upper()}'] = mwdust
            
    maggies = np.zeros(len(phot.bands))
    ivarmaggies = np.zeros(len(phot.bands))
    for iband, (fluxcol, ivarcol) in enumerate(zip(phot.fluxcols, phot.fluxivarcols)):
        maggies[iband] = meta[fluxcol.upper()] / mw_transmission_flux[iband]
        ivarmaggies[iband] = meta[ivarcol.upper()] * mw_transmission_flux[iband]**2
        
    if not np.all(ivarmaggies >= 0.):
        errmsg = 'Some ivarmaggies are negative!'
        log.critical(errmsg)
        raise ValueError(errmsg)
    
    specdata['phot'] = Photometry.parse_photometry(
        phot.bands, maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
        lambda_eff=filters.effective_wavelengths.value,
        min_uncertainty=phot.min_uncertainty, log=log)
        
    if not fastphot:
        from desiutil.dust import dust_transmission
    
        specdata.update({'wave': [], 
                         'flux': [], 
                         'ivar': [], 
                         'mask': [],
                         'res': [], 
                         'res_fast': [], 
                         'snr': np.zeros(3, 'f4'),
                         'linemask': [],
                         'linepix': [],
                         #'linename_balmer_broad': [],
                         #'linepix_balmer_broad': [],
                         })
    
        cameras, npixpercamera = [], []
        for icam, camera in enumerate(specdata['cameras']):
            # Check whether the camera is fully masked.
            if np.sum(specdata['ivar0'][icam]) == 0:
                log.warning(f'Dropping fully masked camera {camera}.')
            else:
                ivar = specdata['ivar0'][icam]
                mask = specdata['mask0'][icam]

                # always mask the first and last pixels
                mask[0] = 1
                mask[-1] = 1

                # In the pipeline, if mask!=0 that does not mean ivar==0, but we
                # want to be more aggressive about masking here.
                ivar[mask != 0] = 0.

                if np.all(ivar == 0.):
                    log.warning(f'Dropping fully masked camera {camera}.')
                else:
                    # interpolate over pixels where the resolution matrix is masked
                    I = mask != 0
                    #I = binary_dilation(mask != 0, iterations=2)
                    if np.any(I):
                        J = np.where(np.logical_not(I))[0]
                        I = np.where(I)[0]
                        res = specdata['res0'][icam]
                        for irow in range(res.shape[0]):
                            res[irow, I] = np.interp(I, J, res[irow, J])

                    # should we also interpolate over the coadded resolution matrix??

                    cameras.append(camera)
                    npixpercamera.append(len(specdata['wave0'][icam])) # number of pixels in this camera

                    # Compute the SNR before we correct for dust.
                    specdata['snr'][icam] = np.median(specdata['flux0'][icam] * np.sqrt(ivar))
                    
                    #mw_transmission_spec = 10**(-0.4 * ebv * RV * ext_odonnell(wave[camera], Rv=RV))
                    mw_transmission_spec = dust_transmission(specdata['wave0'][icam], ebv, Rv=RV)
                    specdata['flux'].append(specdata['flux0'][icam] / mw_transmission_spec)
                    specdata['ivar'].append(ivar * mw_transmission_spec**2)
                    specdata['wave'].append(specdata['wave0'][icam])
                    specdata['mask'].append(specdata['mask0'][icam])
                    # do we need 'res' anymore??
                    specdata['res'].append(Resolution(res))
                    specdata['res_fast'].append(EMLine_Resolution(res))
                    
        if len(cameras) == 0:
            errmsg = 'No good data, which should never happen.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        # clean up the data dictionary
        for key in ['wave0', 'flux0', 'ivar0', 'mask0', 'res0']:
            del specdata[key]

        # Pre-compute some convenience variables for "un-hstacking"
        # an "hstacked" spectrum.
        specdata['cameras'] = cameras
        specdata['npixpercamera'] = npixpercamera
        
        ncam = len(specdata['cameras'])
        npixpercam = np.hstack([0, npixpercamera])
        specdata['camerapix'] = np.zeros((ncam, 2), np.int16)
        for icam in np.arange(ncam):
            specdata['camerapix'][icam, :] = [np.sum(npixpercam[:icam+1]), np.sum(npixpercam[:icam+2])]
                                
        # use the coadded spectrum to build a robust emission-line mask
        pix = CTools.build_linemask_patches(
            specdata['coadd_wave'], specdata['coadd_flux'],
            specdata['coadd_ivar'], specdata['coadd_res_fast'], 
            uniqueid=specdata['uniqueid'], redshift=specdata['zredrock'],
            debug_plots=debug_plots)

        # Map the pixels belonging to individual emission lines onto the
        # original per-camera spectra. This works, but maybe there's a better
        # way to do it?
        for icam in np.arange(ncam):
            camlinepix = {}
            camlinemask = np.zeros(npixpercamera[icam], bool)
            for linename in pix['coadd_linepix'].keys():
                linepix = pix['coadd_linepix'][linename]
                # if the line is entirely off this camera, skip it
                oncam = np.where((specdata["coadd_wave"][linepix] >= np.min(specdata['wave'][icam])) *
                                 (specdata["coadd_wave"][linepix] <= np.max(specdata['wave'][icam])))[0]
                if len(oncam) == 0:
                    continue
                I = np.searchsorted(specdata['wave'][icam], specdata['coadd_wave'][linepix[oncam]])
                #print(f'Line {linename:20}: adding {len(I):02d} pixels to camera {icam}')
                camlinemask[I] = True
                camlinepix[linename] = I
            #print()
            specdata['linemask'].append(camlinemask)
            specdata['linepix'].append(camlinepix)

        specdata.update(pix)
        del pix
    
        # Optionally synthesize photometry from the coadded spectrum.
        if synthphot and synth_filters is not None:
            synthmaggies = Photometry.get_ab_maggies(synth_filters,
                                                     specdata['coadd_flux'] / FLUXNORM,
                                                     specdata['coadd_wave'])
            
            specdata['synthphot'] = Photometry.parse_photometry(
                phot.synth_bands, maggies=synthmaggies, nanomaggies=False,
                lambda_eff=synth_filters.effective_wavelengths.value, log=log)

    return specdata, meta


def _unpack_one_stacked_spectrum(args):
    """Multiprocessing wrapper."""
    return unpack_one_stacked_spectrum(*args)


def unpack_one_stacked_spectrum(iobj, specdata, meta, fphoto, synthphot,
                                ignore_photometry, log):
    """Unpack the data for a single stacked spectrum. Also flag pixels which may be
    affected by emission lines.

    """
    from fastspecfit.continuum import ContinuumTools

    phot = Photometry(fphoto=fphoto, ignore_photometry=ignore_photometry)
    CTools = ContinuumTools(phot)
    
    log.info(f'Pre-processing object {iobj} [stackid {meta[phot.uniqueid]} z={meta["Z"]:.6f}].')

    filters = phot.filters[specdata['photsys']]
    synth_filters = phot.synth_filters[specdata['photsys']]

    # Dummy imaging photometry.
    maggies = np.zeros(len(phot.bands))
    ivarmaggies = np.zeros(len(phot.bands))
    
    specdata['phot'] = Photometry.parse_photometry(
        phot.bands, maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
        lambda_eff=filters.effective_wavelengths.value,
        min_uncertainty=phot.min_uncertainty, log=log)
    
    specdata.update({'linemask': [], 'linemask_all': [], 'linename': [],
                     'linepix': [], 'contpix': [],
                     'wave': [], 'flux': [], 'ivar': [], 'mask': [], 'res': [], 
                     'snr': np.zeros(1, 'f4'),
                     })

    cameras, npixpercamera = [], []
    for icam, camera in enumerate(specdata['cameras']):
        # Check whether the camera is fully masked.
        if np.sum(specdata['ivar0'][icam]) == 0:
            log.warning(f'Dropping fully masked camera {camera}.')
        else:
            ivar = specdata['ivar0'][icam]
            mask = specdata['mask0'][icam]

            # always mask the first and last pixels
            mask[0] = 1
            mask[-1] = 1

            # In the pipeline, if mask!=0 that does not mean ivar==0, but we
            # want to be more aggressive about masking here.
            ivar[mask != 0] = 0

            if np.all(ivar == 0):
                log.warning(f'Dropping fully masked camera {camera}.')
            else:
                cameras.append(camera)
                npixpercamera.append(len(specdata['wave0'][icam])) # number of pixels in this camera

                # Compute the SNR before we correct for dust.
                specdata['snr'][icam] = np.median(specdata['flux0'][icam] * np.sqrt(ivar))
                specdata['flux'].append(specdata['flux0'][icam])
                specdata['ivar'].append(ivar)
                specdata['wave'].append(specdata['wave0'][icam])
                specdata['mask'].append(specdata['mask0'][icam])
                specdata['res'].append(specdata['res0'][icam])

    if len(cameras) == 0:
        errmsg = 'No good data, which should never happen.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    # Pre-compute some convenience variables for "un-hstacking"
    # an "hstacked" spectrum.
    specdata['cameras'] = cameras
    specdata['npixpercamera'] = npixpercamera
    
    ncam = len(specdata['cameras'])
    npixpercam = np.hstack([0, npixpercamera])
    specdata['camerapix'] = npixpercam.reshape(ncam, 2)

    # clean up the data dictionary
    for key in ['wave0', 'flux0', 'ivar0', 'mask0', 'res0']:
        del specdata[key]

    coadd_linemask_dict = CTools.build_linemask_patches(specdata['coadd_wave'], specdata['coadd_flux'],
                                                        specdata['coadd_ivar'], redshift=specdata['zredrock'],
                                                        linetable=CTools.linetable)
    specdata['coadd_linename'] = coadd_linemask_dict['linename']
    specdata['coadd_linepix'] = [np.where(lpix)[0] for lpix in coadd_linemask_dict['linepix']]
    specdata['coadd_contpix'] = [np.where(cpix)[0] for cpix in coadd_linemask_dict['contpix']]

    specdata['linesigma_narrow'] = coadd_linemask_dict['linesigma_narrow']
    specdata['linesigma_balmer'] = coadd_linemask_dict['linesigma_balmer']
    specdata['linesigma_uv'] = coadd_linemask_dict['linesigma_uv']

    specdata['linesigma_narrow_snr'] = coadd_linemask_dict['linesigma_narrow_snr']
    specdata['linesigma_balmer_snr'] = coadd_linemask_dict['linesigma_balmer_snr']
    specdata['linesigma_uv_snr'] = coadd_linemask_dict['linesigma_uv_snr']

    specdata['smoothsigma'] = coadd_linemask_dict['smoothsigma']
    
    # Map the pixels belonging to individual emission lines and
    # their local continuum back onto the original per-camera
    # spectra. These lists of arrays are used in
    # continuum.ContinnuumTools.smooth_continuum.
    for icam in np.arange(len(specdata['cameras'])):
        #specdata['smoothflux'].append(np.interp(specdata['wave'][icam], specdata['coadd_wave'], coadd_linemask_dict['smoothflux']))
        specdata['linemask'].append(np.interp(specdata['wave'][icam], 
                                              specdata['coadd_wave'], 
                                              coadd_linemask_dict['linemask']*1) > 0)
        specdata['linemask_all'].append(np.interp(specdata['wave'][icam], 
                                                  specdata['coadd_wave'], 
                                                  coadd_linemask_dict['linemask_all']*1) > 0)
        _linename, _linenpix, _contpix = [], [], []
        for ipix in np.arange(len(coadd_linemask_dict['linepix'])):
            I = np.interp(specdata['wave'][icam], specdata['coadd_wave'], 
                          coadd_linemask_dict['linepix'][ipix]*1) > 0
            J = np.interp(specdata['wave'][icam], specdata['coadd_wave'], 
                          coadd_linemask_dict['contpix'][ipix]*1) > 0
            if np.sum(I) > 3 and np.sum(J) > 3:
                _linename.append(coadd_linemask_dict['linename'][ipix])
                _linenpix.append(np.where(I)[0])
                _contpix.append(np.where(J)[0])
        specdata['linename'].append(_linename)
        specdata['linepix'].append(_linenpix)
        specdata['contpix'].append(_contpix)

    specdata.update({'coadd_linemask': coadd_linemask_dict['linemask'],
                     'coadd_linemask_all': coadd_linemask_dict['linemask_all']})

    # Optionally synthesize photometry from the coadded spectrum.
    if synthphot:
        synthmaggies = Photometry.get_ab_maggies(synth_filters, specdata['coadd_flux'] / FLUXNORM, specdata['coadd_wave'])
        
        specdata['synthphot'] = Photometry.parse_photometry(phot.synth_bands,
            maggies=synthmaggies, nanomaggies=False,
            lambda_eff=synth_filters.effective_wavelengths.value, log=log)

    return specdata, meta


class DESISpectra(object):
    def __init__(self, stackfit=False, redux_dir=None, fiberassign_dir=None,
                 fphotodir=None, fphotofile=None, mapdir=None):
        """Class to read in DESI spectra and associated metadata.

        Parameters
        ----------
        redux_dir : str
            Full path to the location of the reduced DESI data. Optional and
            defaults to `$DESI_SPECTRO_REDUX`.
        
        fiberassign_dir : str
            Full path to the location of the fiberassign files. Optional and
            defaults to `$DESI_ROOT/target/fiberassign/tiles/trunk`.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.

        """
        import yaml

        self.cosmo = TabulatedDESI()
        
        desi_root = os.environ.get('DESI_ROOT', DESI_ROOT_NERSC)

        if redux_dir is None:
            self.redux_dir = os.path.join(desi_root, 'spectro', 'redux')
        else:
            self.redux_dir = redux_dir
            
        if fiberassign_dir is None:
            self.fiberassign_dir = os.path.join(desi_root, 'target', 'fiberassign', 'tiles', 'trunk')
        else:
            self.fiberassign_dir = fiberassign_dir

        if fphotodir is None:
            self.fphotoext = None
            self.fphotodir = os.environ.get('FPHOTO_DIR', FPHOTO_DIR_NERSC)
        else:
            # parse the extension name, if any
            fphotoext = None
            photodir = os.path.dirname(fphotodir)
            photobase = os.path.basename(fphotodir)
            if '[' in photobase and ']' in photobase:
                try:
                    fphotoext = photobase[photobase.find('[')+1:photobase.find(']')]
                    fphotodir = os.path.join(photodir, photobase[:photobase.find('[')])
                except:
                    pass
            self.fphotoext = fphotoext
            self.fphotodir = fphotodir

        if fphotofile is None:
            from importlib import resources
            if stackfit:
                fphotofile = resources.files('fastspecfit').joinpath('data/stacked-phot.yaml')
            else:
                fphotofile = resources.files('fastspecfit').joinpath('data/legacysurvey-dr9.yaml')

        try:
            with open(fphotofile, 'r') as F:
                fphoto = yaml.safe_load(F)
            self.fphotofile = fphotofile
        except:
            errmsg = f'Unable to read parameter file {fphotofile}'
            log.critical(errmsg)
            raise ValueError(errmsg)

        self.fphoto = fphoto

        if mapdir is None:
            self.mapdir = os.path.join(os.environ.get('DUST_DIR', DUST_DIR_NERSC), 'maps')
        else:
            self.mapdir = mapdir

            
    @staticmethod
    def resolve(targets):
        """Resolve which targets are primary in imaging overlap regions.
    
        Parameters
        ----------
        targets : :class:`~numpy.ndarray`
            Rec array of targets. Must have columns "RA" and "DEC" and
            either "RELEASE" or "PHOTSYS" or "TARGETID".
    
        Returns
        -------
        :class:`~numpy.ndarray`
            The original target list trimmed to only objects from the "northern"
            photometry in the northern imaging area and objects from "southern"
            photometry in the southern imaging area.
        
        """
        import healpy as hp
        
        def _isonnorthphotsys(photsys):
            """ If the object is from the northen photometric system """
            # ADM explicitly checking for NoneType. In the past we have had bugs
            # ADM where we forgot to populate variables before passing them.
            if photsys is None:
                msg = "NoneType submitted to _isonnorthphotsys function"
                log.critical(msg)
                raise ValueError(msg)
        
            psftype = np.asarray(photsys)
            # ADM in Python3 these string literals become byte-like
            # ADM so to retain Python2 compatibility we need to check
            # ADM against both bytes and unicode.
            northern = ((photsys == 'N') | (photsys == b'N'))
        
            return northern
                
        # ADM retrieve the photometric system from the RELEASE.
        from desitarget.io import release_to_photsys, desitarget_resolve_dec
        if 'PHOTSYS' in targets.dtype.names:
            photsys = targets["PHOTSYS"]
        else:
            if 'RELEASE' in targets.dtype.names:
                photsys = release_to_photsys(targets["RELEASE"])
            else:
                _, _, release, _, _, _ = decode_targetid(targets["TARGETID"])
                photsys = release_to_photsys(release)
    
        # ADM a flag of which targets are from the 'N' photometry.
        photn = _isonnorthphotsys(photsys)
    
        # ADM grab the declination used to resolve targets.
        split = desitarget_resolve_dec()
    
        # ADM determine which targets are north of the Galactic plane. As
        # ADM a speed-up, bin in ~1 sq.deg. HEALPixels and determine
        # ADM which of those pixels are north of the Galactic plane.
        # ADM We should never be as close as ~1o to the plane.
        from desitarget.geomask import is_in_gal_box, pixarea2nside
        nside = pixarea2nside(1)
        theta, phi = np.radians(90-targets["DEC"]), np.radians(targets["RA"])
        pixnum = hp.ang2pix(nside, theta, phi, nest=True)
        # ADM find the pixels north of the Galactic plane...
        allpix = np.arange(hp.nside2npix(nside))
        theta, phi = hp.pix2ang(nside, allpix, nest=True)
        ra, dec = np.degrees(phi), 90-np.degrees(theta)
        pixn = is_in_gal_box([ra, dec], [0., 360., 0., 90.], radec=True)
        # ADM which targets are in pixels north of the Galactic plane.
        galn = pixn[pixnum]
    
        # ADM which targets are in the northern imaging area.
        arean = (targets["DEC"] >= split) & galn
    
        # ADM retain 'N' targets in 'N' area and 'S' in 'S' area.
        #keep = (photn & arean) | (~photn & ~arean)
        #return targets[keep]

        inorth = photn & arean
        newphotsys = np.array(['S'] * len(targets))
        newphotsys[inorth] = 'N'

        return newphotsys

    
    def select(self, redrockfiles, zmin=None, zmax=None, zwarnmax=None,
               targetids=None, firsttarget=0, ntargets=None,
               input_redshifts=None, specprod_dir=None, use_quasarnet=True,
               redrockfile_prefix='redrock-', specfile_prefix='coadd-',
               qnfile_prefix='qso_qn-'):
        """Select targets for fitting and gather the necessary spectroscopic metadata.

        Parameters
        ----------
        redrockfiles : str or array
            Full path to one or more input Redrock file(s).
        zmin : float
            Minimum redshift of observed targets to select. Defaults to
            0.001. Note that any value less than or equal to zero will raise an
            exception because a positive redshift is needed to compute the
            distance modulus when modeling the broadband photometry.
        zmax : float or `None`
            Maximum redshift of observed targets to select. `None` is equivalent
            to not having an upper redshift limit.
        zwarnmax : int or `None`
            Maximum Redrock zwarn value for selected targets. `None` is
            equivalent to not having any cut on zwarn.
        targetids : int or array or `None`
            Restrict the sample to the set of targetids in this list. If `None`,
            select all targets which satisfy the selection criteria.
        firsttarget : int
            Integer offset of the first object to consider in each file. Useful
            for debugging and testing. Defaults to 0.
        ntargets : int or `None`
            Number of objects to analyze in each file. Useful for debugging and
            testing. If `None`, select all targets which satisfy the selection
            criteria.
        input_redshifts : float or array or `None`
            Input redshifts to use for each input `targetids` input If `None`,
            use the nominal Redrock (or QuasarNet) redshifts.
        use_quasarnet : `bool`
            Use QuasarNet to improve QSO redshifts, if the afterburner file is
            present. Defaults to `True`.
        redrockfile_prefix : str
            Prefix of the `redrockfiles`. Defaults to `redrock-`.
        specfile_prefix : str
            Prefix of the spectroscopic coadds corresponding to the input
            Redrock file(s). Defaults to `coadd-`.
        qnfile_prefix : str
            Prefix of the QuasarNet afterburner file. Defaults to `qso_qn-`.

        Attributes
        ----------
        coadd_type : str
            Type of coadded spectra (healpix, cumulative, pernight, or perexp).
        meta : list of :class:`astropy.table.Table`
            Array of tables (one per input `redrockfile`) with the metadata
            needed to fit the data and to write to the output file(s).
        redrockfiles : str array
            Input Redrock file names.
        specfiles : str array
            Spectroscopic file names corresponding to each each input Redrock file.
        specprod : str
            Spectroscopic production name for the input Redrock file.

        Notes
        -----
        We assume that `specprod` is the same for all input Redrock files,
        although we don't explicitly do this check. Specifically, we only read
        the header of the first file.

        """
        from astropy.table import vstack, hstack
        from desiutil.depend import getdep
        from desitarget import geomask
        from desitarget.targets import main_cmx_or_sv

        if zmin is None:
            zmin = 1e-3

        if zmin <= 0.0:
            errmsg = 'zmin should generally be >= 0; proceed with caution!'
            log.warning(errmsg)
        
        if zmax is None:
            zmax = 99.0

        if zwarnmax is None:
            zwarnmax = 99999
            
        if zmin >= zmax:
            errmsg = 'zmin must be <= zmax.'
            log.critical(errmsg)
            raise ValueError(errmsg)
        
        if redrockfiles is None:
            errmsg = 'At least one redrockfiles file is required.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        if len(np.atleast_1d(redrockfiles)) == 0:
            errmsg = 'No redrockfiles found!'
            log.warning(errmsg)
            raise ValueError(errmsg)

        # Should we not sort...?
        #redrockfiles = np.array(set(np.atleast_1d(redrockfiles)))
        redrockfiles = np.array(sorted(set(np.atleast_1d(redrockfiles))))
        log.info(f'Reading and parsing {len(redrockfiles)} unique redrockfile(s).')

        alltiles = []
        self.redrockfiles, self.specfiles, self.meta, self.surveys = [], [], [], []

        for ired, redrockfile in enumerate(np.atleast_1d(redrockfiles)):
            if not os.path.isfile(redrockfile):
                log.warning(f'File {redrockfile} not found!')
                continue

            if not redrockfile_prefix in redrockfile:
                errmsg = f'Redrockfile {redrockfile} missing standard prefix {redrockfile_prefix}; ' + \
                    'please specify redrockfile_prefix argument.'
                log.critical(errmsg)
                raise ValueError(errmsg)
            
            specfile = os.path.join(os.path.dirname(redrockfile), os.path.basename(redrockfile).replace(redrockfile_prefix, specfile_prefix))
            if not os.path.isfile(specfile):
                log.warning(f'File {specfile} not found!')
                continue
            
            # Can we use the quasarnet afterburner file to improve QSO redshifts?
            qnfile = redrockfile.replace(redrockfile_prefix, qnfile_prefix)
            if os.path.isfile(qnfile) and use_quasarnet and input_redshifts is None:
                use_qn = True
            else:
                use_qn = False

            # Gather some coadd information from the header. Note: this code is
            # only compatible with Fuji & Guadalupe headers and later.
            hdr = fitsio.read_header(specfile, ext=0)

            specprod = getdep(hdr, 'SPECPROD')
            if hasattr(self, 'specprod'):
                if self.specprod != specprod:
                    errmsg = f'specprod must be the same for all input redrock files! {specprod}!={self.specprod}'
                    log.critical(errmsg)
                    raise ValueError(errmsg)
            
            self.specprod = specprod
            
            if 'SPGRP' in hdr:
                self.coadd_type = hdr['SPGRP']
            else:
                errmsg = f'SPGRP header card missing from spectral file {specfile}'
                log.warning(errmsg)
                self.coadd_type = 'custom'

            #log.info('specprod={}, coadd_type={}'.format(self.specprod, self.coadd_type))

            if self.coadd_type == 'healpix':
                survey = hdr['SURVEY']
                program = hdr['PROGRAM']
                healpix = np.int32(hdr['SPGRPVAL'])
                thrunight = None
                log.info('specprod={}, coadd_type={}, survey={}, program={}, healpix={}'.format(
                    self.specprod, self.coadd_type, survey, program, healpix))

                # I'm not sure we need these attributes but if we end up
                # using them then be sure to document them as attributes of
                # the class!
                #self.hpxnside = hdr['HPXNSIDE']
                #self.hpxnest = hdr['HPXNEST']
            elif self.coadd_type == 'custom':
                survey = 'custom'
                program = 'custom'
                healpix = np.int32(0)
                thrunight = None
                log.info('specprod={}, coadd_type={}, survey={}, program={}, healpix={}'.format(
                    self.specprod, self.coadd_type, survey, program, healpix))
            else:
                tileid = np.int32(hdr['TILEID'])
                petal = np.int16(hdr['PETAL'])
                night = np.int32(hdr['NIGHT']) # thrunight for coadd_type==cumulative
                if self.coadd_type == 'perexp':
                    expid = np.int32(hdr['EXPID'])
                    log.info('specprod={}, coadd_type={}, tileid={}, petal={}, night={}, expid={}'.format(
                        self.specprod, self.coadd_type, tileid, petal, night, expid))
                else:
                    expid = None
                    log.info('specprod={}, coadd_type={}, tileid={}, petal={}, night={}'.format(
                        self.specprod, self.coadd_type, tileid, petal, night))

                # cache the tiles file so we can grab the survey and program name appropriate for this tile
                if not hasattr(self, 'tileinfo'):
                    if specprod_dir is None:
                        specprod_dir = os.path.join(self.redux_dir, self.specprod)
                    infofile = os.path.join(specprod_dir, 'tiles-{}.csv'.format(self.specprod))
                    if os.path.isfile(infofile):
                        self.tileinfo = Table.read(infofile)
                        
                if hasattr(self, 'tileinfo'):
                    tileinfo = self.tileinfo[self.tileinfo['TILEID'] == tileid]
                    survey = tileinfo['SURVEY'][0]
                    program = tileinfo['PROGRAM'][0]
                else:
                    survey, program = '', ''

            if survey == 'main' or survey == 'special':
                TARGETINGCOLS = TARGETINGBITS['default']
            else:
                TARGETINGCOLS = TARGETINGBITS['all']
            
            # add targeting columns
            allfmcols = np.array(fitsio.FITS(specfile)['FIBERMAP'].get_colnames())
            READFMCOLS = FMCOLS + [col for col in TARGETINGCOLS if col in allfmcols]
                    
            # If targetids is *not* given we have to choose "good" objects
            # before subselecting (e.g., we don't want sky spectra).
            if targetids is None:
                zb = fitsio.read(redrockfile, 'REDSHIFTS', columns=REDSHIFTCOLS)
                # Are we reading individual exposures or coadds?
                meta = fitsio.read(specfile, 'FIBERMAP', columns=READFMCOLS)
                # Check for uniqueness.
                uu, cc = np.unique(meta['TARGETID'], return_counts=True)
                if np.any(cc > 1):
                    errmsg = 'Found {} duplicate TARGETIDs in {}: {}'.format(
                        np.sum(cc>1), specfile, ' '.join(uu[cc > 1].astype(str)))
                    log.critical(errmsg)
                    raise ValueError(errmsg)
                assert(np.all(zb['TARGETID'] == meta['TARGETID']))
                # need to also update mpi.get_ntargets_one
                if use_qn:
                    # If using QuasarNet, it can happen that zb['Z']<zmin and
                    # therefore the object falls out of the sample before we
                    # have a chance to even read it. So apply the minimum
                    # redshift cut below after we correct the redshift.
                    fitindx = np.where((zb['Z'] < zmax) *
                                       (meta['OBJTYPE'] == 'TGT') * (zb['ZWARN'] <= zwarnmax) *
                                       (zb['ZWARN'] & ZWarningMask.NODATA == 0))[0]
                else:
                    fitindx = np.where((zb['Z'] > zmin) * (zb['Z'] < zmax) *
                                       (meta['OBJTYPE'] == 'TGT') * (zb['ZWARN'] <= zwarnmax) *
                                       (zb['ZWARN'] & ZWarningMask.NODATA == 0))[0]
            else:
                # We already know we like the input targetids, so no selection
                # needed. But make sure there are no duplicates.
                uu, cc = np.unique(targetids, return_counts=True)
                if np.any(cc > 1):
                    errmsg = 'Found {} duplicate TARGETIDs in {}: {}'.format(
                        np.sum(cc>1), specfile, ' '.join(uu[cc > 1].astype(str)))
                    log.critical(errmsg)
                    raise ValueError(errmsg)
                
                alltargetids = fitsio.read(redrockfile, 'REDSHIFTS', columns='TARGETID')
                fitindx = np.where(np.isin(alltargetids, targetids))[0]

            if len(fitindx) == 0:
                log.info(f'No requested targets found in redrockfile {redrockfile}')
                continue

            # Do we want just a subset of the available objects?
            if ntargets is None:
                _ntargets = len(fitindx)
            else:
                _ntargets = ntargets
            if _ntargets > len(fitindx):
                log.warning(f'Number of requested ntargets exceeds the number of targets on {redrockfile}; reading all of them.')

            __ntargets = len(fitindx)
            fitindx = fitindx[firsttarget:firsttarget+_ntargets]
            if len(fitindx) == 0:
                log.info('All {} targets in redrockfile {} have been dropped (firsttarget={}, ntargets={}).'.format(
                    __ntargets, redrockfile, firsttarget, _ntargets))
                continue
                
            # If firsttarget is a large index then the set can become empty.
            if targetids is None:
                zb = Table(zb[fitindx])
                meta = Table(meta[fitindx])
            else:
                zb = Table(fitsio.read(redrockfile, 'REDSHIFTS', rows=fitindx, columns=REDSHIFTCOLS))
                meta = Table(fitsio.read(specfile, 'FIBERMAP', rows=fitindx, columns=READFMCOLS))
                if input_redshifts is not None:
                    log.info(f'Applying {len(input_redshifts)} input_redshifts.')
                    # fitsio doesn't preserve order, so make sure targetids and
                    # input_redshifts are matched
                    srt = np.hstack([np.where(targetids == tid)[0] for tid in zb['TARGETID']])
                    targetids = np.array(targetids)[srt]
                    input_redshifts = np.array(input_redshifts)[srt]
                    assert(np.all(zb['TARGETID'] == targetids))
                    zb['Z'] = input_redshifts
                
            # Update the redrock redshift when quasarnet disagrees **but only
            # for QSO targets**. From Edmond: the QN afterburner is run with a
            # threshold 0.5. With VI, we choose 0.95 as final threshold. Note,
            # the IS_QSO_QN_NEW_RR column contains only QSO for QN which are not
            # QSO for RR.
            zb['Z_RR'] = zb['Z'] # add it at the end
            if use_qn:
                surv_target, surv_mask, surv = main_cmx_or_sv(meta)
                if surv == 'cmx':
                    desi_target = surv_target[0]
                    desi_mask = surv_mask[0]
                    # need to check multiple QSO masks
                    IQSO = []
                    for bitname in desi_mask.names():
                        if 'QSO' in bitname:
                            IQSO.append(np.where(meta[desi_target] & desi_mask[bitname] != 0)[0])
                    if len(IQSO) > 0:
                        IQSO = np.sort(np.unique(np.hstack(IQSO)))
                else:
                    desi_target, bgs_target, mws_target = surv_target
                    desi_mask, bgs_mask, mws_mask = surv_mask
                    IQSO = np.where(meta[desi_target] & desi_mask['QSO'] != 0)[0]
                if len(IQSO) > 0:
                    qn = Table(fitsio.read(qnfile, 'QN_RR', rows=fitindx[IQSO], columns=QNCOLS))
                    assert(np.all(qn['TARGETID'] == meta['TARGETID'][IQSO]))
                    log.info('Updating QSO redshifts using a QN threshold of 0.95.')
                    qn['IS_QSO_QN'] = np.max(np.array([qn[name] for name in QNLINES]), axis=0) > 0.95
                    qn['IS_QSO_QN_NEW_RR'] &= qn['IS_QSO_QN']
                    if np.count_nonzero(qn['IS_QSO_QN_NEW_RR']) > 0:
                        zb['Z'][IQSO[qn['IS_QSO_QN_NEW_RR']]] = qn['Z_NEW'][qn['IS_QSO_QN_NEW_RR']]
                    del qn
                # now apply zmin
                keep = np.where(zb['Z'] > zmin)[0]
                if len(keep) == 0:
                    log.info('No requested targets found in redrockfile {}'.format(redrockfile))
                    continue
                zb = zb[keep]
                meta = meta[keep]
                fitindx = fitindx[keep]

            tsnr2 = Table(fitsio.read(redrockfile, 'TSNR2', rows=fitindx, columns=TSNR2COLS))
            assert(np.all(zb['TARGETID'] == meta['TARGETID']))

            # astropy 5.0 "feature" -- join no longer preserves order, ugh.
            zb.remove_column('TARGETID')
            meta = hstack((zb, meta, tsnr2))
            #meta = join(zb, meta, keys='TARGETID')
            del zb, tsnr2

            # make sure we're sorted
            if targetids is not None:
                srt = geomask.match_to(meta['TARGETID'], targetids)
                meta = meta[srt]
                assert(np.all(meta['TARGETID'] == targetids))

            # Get the unique set of tiles contributing to the coadded spectra
            # from EXP_FIBERMAP.
            expmeta = fitsio.read(specfile, 'EXP_FIBERMAP', columns=EXPFMCOLS[self.coadd_type])
            I = np.isin(expmeta['TARGETID'], meta['TARGETID'])
            if np.count_nonzero(I) == 0:
                errmsg = 'No matching targets in exposure table.'
                log.critical(errmsg)
                raise ValueError(errmsg)
            expmeta = Table(expmeta[I])

            # build the list of tiles that went into each unique target / coadd
            tileid_list = [] # variable length, so need to build the array first
            for tid in meta['TARGETID']:
                I = tid == expmeta['TARGETID']
                tileid_list.append(' '.join(np.unique(expmeta['TILEID'][I]).astype(str)))
                #meta['TILEID_LIST'][M] = ' '.join(np.unique(expmeta['TILEID'][I]).astype(str))
                # store just the zeroth tile for gather_targetphot, below
                if self.coadd_type == 'healpix':
                    alltiles.append(expmeta['TILEID'][I][0]) 
                elif self.coadd_type == 'custom':
                    alltiles.append(expmeta['TILEID'][I][0])
                else:
                    alltiles.append(tileid)
                    
            if self.coadd_type == 'healpix':                    
                meta['TILEID_LIST'] = tileid_list
            elif self.coadd_type == 'custom':
                meta['TILEID_LIST'] = tileid_list

            # Gather additional info about this pixel.
            if self.coadd_type == 'healpix':
                meta['SURVEY'] = survey
                meta['PROGRAM'] = program
                meta['HEALPIX'] = healpix
            elif self.coadd_type == 'custom':
                meta['SURVEY'] = survey
                meta['PROGRAM'] = program
                meta['HEALPIX'] = healpix
            else:
                if hasattr(self, 'tileinfo'):
                    meta['SURVEY'] = survey
                    meta['PROGRAM'] = program
                meta['NIGHT'] = night
                meta['TILEID'] = tileid
                if expid:
                    meta['EXPID'] = expid

                # get the correct fiber number
                if 'FIBER' in expmeta.colnames:
                    meta['FIBER'] = np.zeros(len(meta), dtype=expmeta['FIBER'].dtype)
                    for iobj, tid in enumerate(meta['TARGETID']):
                        iexp = np.where(expmeta['TARGETID'] == tid)[0][0] # zeroth
                        meta['FIBER'][iobj] = expmeta['FIBER'][iexp]

            self.meta.append(Table(meta))
            self.redrockfiles.append(redrockfile)
            self.specfiles.append(specfile)
            self.surveys.append(survey)

        if len(self.meta) == 0:
            log.warning('No targets read!')
            return

        # Use the metadata in the fibermap to retrieve the LS-DR9 source
        # photometry. Note that we have to make a copy of the input_meta table
        # because otherwise BRICKNAME gets "repaired!"
        t0 = time.time()  
        metas = self._gather_photometry(specprod=specprod, alltiles=alltiles)
        self.meta = metas # update
        log.info(f'Gathered photometric metadata in {time.time()-t0:.2f} sec')


    def read_and_unpack(self, fastphot=False, synthphot=True, ignore_photometry=False,
                        constrain_age=False, debug_plots=False, verbose=False, mp=1):
        """Read and unpack selected spectra or broadband photometry.
        
        Parameters
        ----------
        fastphot : bool
            Read and unpack the broadband photometry; otherwise, handle the DESI
            three-camera spectroscopy. Optional; defaults to `False`.
        synthphot : bool
            Synthesize photometry from the coadded optical spectrum. Optional;
            defaults to `True`.
        remember_coadd : bool
            Add the coadded spectrum to the returned dictionary. Optional;
            defaults to `False` (in order to reduce memory usage).

        Returns
        -------
        List of dictionaries (:class:`dict`, one per object) the following keys:
            targetid : numpy.int64
                DESI target ID.
            zredrock : numpy.float64
                Redrock redshift.
            cameras : :class:`list`
                List of camera names present for this spectrum.
            wave : :class:`list`
                Three-element list of `numpy.ndarray` wavelength vectors, one for
                each camera.    
            flux : :class:`list`
                Three-element list of `numpy.ndarray` flux spectra, one for each
                camera and corrected for Milky Way extinction.
            ivar : :class:`list`
                Three-element list of `numpy.ndarray` inverse variance spectra, one
                for each camera.    
            res : :class:`list`
                Three-element list of :class:`desispec.resolution.Resolution`
                objects, one for each camera.
            snr : `numpy.ndarray`
                Median per-pixel signal-to-noise ratio in the grz cameras.
            linemask : :class:`list`
                Three-element list of `numpy.ndarray` boolean emission-line masks,
                one for each camera. This mask is used during continuum-fitting.
            linename : :class:`list`
                Three-element list of emission line names which might be present
                in each of the three DESI cameras.
            linepix : :class:`list`
                Three-element list of pixel indices, one per camera, which were
                identified in :class:`FFit.build_linemask` to belong to emission
                lines.
            contpix : :class:`list`
                Three-element list of pixel indices, one per camera, which were
                identified in :class:`FFit.build_linemask` to not be
                "contaminated" by emission lines.
            coadd_wave : `numpy.ndarray`
                Coadded wavelength vector with all three cameras combined.
            coadd_flux : `numpy.ndarray`
                Flux corresponding to `coadd_wave`.
            coadd_ivar : `numpy.ndarray`
                Inverse variance corresponding to `coadd_flux`.
            photsys : str
                Photometric system.
            phot : `astropy.table.Table`
                Total photometry in `grzW1W2`, corrected for Milky Way extinction.
            fiberphot : `astropy.table.Table`
                Fiber photometry in `grzW1W2`, corrected for Milky Way extinction.
            fibertotphot : `astropy.table.Table`
                Fibertot photometry in `grzW1W2`, corrected for Milky Way extinction.
            synthphot : :class:`astropy.table.Table`
                Photometry in `grz` synthesized from the Galactic
                extinction-corrected coadded spectra (with a mild extrapolation
                of the data blueward and redward to accommodate the g-band and
                z-band filter curves, respectively.

        """
        from astropy.table import vstack
        from desitarget import geomask
        from desispec.coaddition import coadd_cameras
        from desispec.io import read_spectra
        from desiutil.dust import SFDMap
        from desiutil.log import get_logger, DEBUG
        from desispec.resolution import Resolution
        from fastspecfit.emline_fit import EMLine_Resolution
        
        if verbose:
            log = get_logger(DEBUG)
        else:
            log = get_logger()

        phot = Photometry(fphoto=self.fphoto, ignore_photometry=ignore_photometry)
        
        SFD = SFDMap(scaling=1.0, mapdir=self.mapdir)

        alldata = []
        for ispec, (specfile, meta) in enumerate(zip(self.specfiles, self.meta)):
            nobj = len(meta)
            if nobj == 1:
                log.info(f'Reading {nobj} spectrum from {specfile}')
            else:
                log.info(f'Reading {nobj} spectra from {specfile}')

            ebv = SFD.ebv(meta['RA'], meta['DEC'])

            # Pre-compute the luminosity distance, distance modulus, and age of
            # the universe.
            zobj = meta['Z']
            neg = zobj <= 0.
            if np.any(neg):
                errmsg = f'{np.sum(neg)}/{len(zobj)} input redshifts are zero or negative!'
                log.warning(errmsg)
                zobj[neg] = 1e-8

            dlum = self.cosmo.luminosity_distance(zobj)
            dmod = self.cosmo.distance_modulus(zobj)
            if constrain_age:
                tuniv = self.universe_age(zobj)
            else:
                tuniv = np.zeros_like(zobj) + 100.

            if 'PHOTSYS' in meta.colnames:
                photsys = meta['PHOTSYS']
            else:
                photsys = [''] * len(meta)
            
            if fastphot:
                unpackargs = []
                for iobj in np.arange(len(meta)):
                    specdata = {
                        'uniqueid': meta[phot.uniqueid][iobj], 'zredrock': meta['Z'][iobj],
                        'photsys': photsys[iobj], 'dluminosity': dlum[iobj], 
                        'dmodulus': dmod[iobj], 'tuniv': tuniv[iobj],
                        }
                    unpackargs.append((iobj, specdata, meta[iobj], ebv[iobj], self.fphoto,
                                       True, False, ignore_photometry, debug_plots, log))
            else:
                # Don't use .select since meta and spec can be sorted
                # differently if a non-sorted targetids was passed. Do the
                # selection and sort ourselves.
                spec = read_spectra(specfile)#.select(targets=meta[phot.uniqueid])
                srt = geomask.match_to(spec.fibermap[phot.uniqueid], meta['TARGETID'])
                spec = spec[srt]
                assert(np.all(spec.fibermap[phot.uniqueid] == meta[phot.uniqueid]))

                # Coadd across cameras.
                t0 = time.time()                
                coadd_spec = coadd_cameras(spec)
                log.info(f'Coadding across cameras took {time.time()-t0:.2f} seconds.')

                # unpack the desispec.spectra.Spectra objects into simple arrays
                cameras = spec.bands
                coadd_cameras = coadd_spec.bands[0]
                unpackargs = []
                for iobj in np.arange(len(meta)):
                    specdata = {
                        'uniqueid': meta[phot.uniqueid][iobj], 'zredrock': meta['Z'][iobj],
                        'photsys': photsys[iobj], 'cameras': cameras,
                        'dluminosity': dlum[iobj], 'dmodulus': dmod[iobj], 'tuniv': tuniv[iobj],
                        'wave0': [spec.wave[cam] for cam in cameras],
                        'flux0': [spec.flux[cam][iobj, :] for cam in cameras],
                        'ivar0': [spec.ivar[cam][iobj, :] for cam in cameras],
                        # Also track the mask---see https://github.com/desihub/desispec/issues/1389 
                        'mask0': [spec.mask[cam][iobj, :] for cam in cameras],
                        'res0': [spec.resolution_data[cam][iobj, :, :] for cam in cameras],
                        #'res0': [Resolution(spec.resolution_data[cam][iobj, :, :]) for cam in cameras],
                        #'res_fast0': [EMLine_Resolution(spec.resolution_data[cam][iobj, :, :]) for cam in cameras],
                        'coadd_wave': coadd_spec.wave[coadd_cameras],
                        'coadd_flux': coadd_spec.flux[coadd_cameras][iobj, :],
                        'coadd_ivar': coadd_spec.ivar[coadd_cameras][iobj, :],
                        'coadd_res': [Resolution(coadd_spec.resolution_data[coadd_cameras][iobj, :])],
                        'coadd_res_fast': [EMLine_Resolution(coadd_spec.resolution_data[coadd_cameras][iobj, :])],
                        }
                    unpackargs.append((iobj, specdata, meta[iobj], ebv[iobj], self.fphoto, fastphot,
                                       synthphot, ignore_photometry, debug_plots, log))
                    
            if mp > 1:
                import multiprocessing
                with multiprocessing.Pool(mp) as P:
                    out = P.map(_unpack_one_spectrum, unpackargs)
            else:
                out = [unpack_one_spectrum(*_unpackargs) for _unpackargs in unpackargs]
                
            out = list(zip(*out))
            self.meta[ispec] = Table(np.hstack(out[1]))
            alldata.append(out[0])
            del out
    
        alldata = np.concatenate(alldata)
        self.meta = vstack(self.meta)
        self.ntargets = len(self.meta)

        return alldata

    
    def read_stacked(self, stackfiles, firsttarget=0, ntargets=None,
                     stackids=None, synthphot=True, ignore_photometry=False,
                     mp=1):
        """Read one or more stacked spectra.
        
        Parameters
        ----------
        stackfiles : str or array
            Full path to one or more input stacked-spectra file(s).
        stackids : int or array or `None`
            Restrict the sample to the set of stackids in this list. If `None`,
            fit everything.
        firsttarget : int
            Integer offset of the first object to consider in each file. Useful
            for debugging and testing. Defaults to 0.
        ntargets : int or `None`
            Number of objects to analyze in each file. Useful for debugging and
            testing. If `None`, select all targets which satisfy the selection
            criteria.
        synthphot : bool
            Synthesize photometry from the coadded optical spectrum. Optional;
            defaults to `True`.

        Returns
        -------
        List of dictionaries (:class:`dict`, one per object) the following keys:
            targetid : numpy.int64
                DESI target ID.
            zredrock : numpy.float64
                Redrock redshift.
            cameras : :class:`list`
                List of camera names present for this spectrum.
            wave : :class:`list`
                Three-element list of `numpy.ndarray` wavelength vectors, one for
                each camera.    
            flux : :class:`list`
                Three-element list of `numpy.ndarray` flux spectra, one for each
                camera and corrected for Milky Way extinction.
            ivar : :class:`list`
                Three-element list of `numpy.ndarray` inverse variance spectra, one
                for each camera.    
            res : :class:`list`
                Three-element list of :class:`desispec.resolution.Resolution`
                objects, one for each camera.
            snr : `numpy.ndarray`
                Median per-pixel signal-to-noise ratio in the grz cameras.
            linemask : :class:`list`
                Three-element list of `numpy.ndarray` boolean emission-line masks,
                one for each camera. This mask is used during continuum-fitting.
            linename : :class:`list`
                Three-element list of emission line names which might be present
                in each of the three DESI cameras.
            linepix : :class:`list`
                Three-element list of pixel indices, one per camera, which were
                identified in :class:`FFit.build_linemask` to belong to emission
                lines.
            contpix : :class:`list`
                Three-element list of pixel indices, one per camera, which were
                identified in :class:`FFit.build_linemask` to not be
                "contaminated" by emission lines.
            coadd_wave : `numpy.ndarray`
                Coadded wavelength vector with all three cameras combined.
            coadd_flux : `numpy.ndarray`
                Flux corresponding to `coadd_wave`.
            coadd_ivar : `numpy.ndarray`
                Inverse variance corresponding to `coadd_flux`.
            photsys : str
                Photometric system.
            phot : `astropy.table.Table`
                Total photometry in `grzW1W2`, corrected for Milky Way extinction.
            fiberphot : `astropy.table.Table`
                Fiber photometry in `grzW1W2`, corrected for Milky Way extinction.
            fibertotphot : `astropy.table.Table`
                Fibertot photometry in `grzW1W2`, corrected for Milky Way extinction.
            synthphot : :class:`astropy.table.Table`
                Photometry in `grz` synthesized from the Galactic
                extinction-corrected coadded spectra (with a mild extrapolation
                of the data blueward and redward to accommodate the g-band and
                z-band filter curves, respectively.

        """
        from astropy.table import vstack
        from scipy.sparse import identity
        from desispec.resolution import Resolution

        phot = Photometry(fphoto=self.fphoto, ignore_photometry=ignore_photometry)
        
        if stackfiles is None:
            errmsg = 'At least one stackfiles file is required.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        if len(np.atleast_1d(stackfiles)) == 0:
            errmsg = 'No stackfiles found!'
            log.warning(errmsg)
            raise ValueError(errmsg)

        stackfiles = np.array(sorted(set(np.atleast_1d(stackfiles))))
        log.info('Reading and parsing {} unique stackfile(s).'.format(len(stackfiles)))

        self.specprod = 'stacked'
        self.coadd_type = 'stacked'
        survey = 'stacked'
        program = 'stacked'
        healpix = np.int32(0)
        
        READCOLS = ['STACKID', 'Z']
        
        self.stackfiles, self.meta = [], []
        
        for istack, stackfile in enumerate(np.atleast_1d(stackfiles)):
            if not os.path.isfile(stackfile):
                log.warning('File {} not found!'.format(stackfile))
                continue

            # Gather some coadd information from the header.
            #hdr = fitsio.read_header(stackfile, ext=0)

            log.info('specprod={}, coadd_type={}, survey={}, program={}, healpix={}'.format(
                self.specprod, self.coadd_type, survey, program, healpix))
                    
            # If stackids is *not* given, read everything.
            if stackids is None:
                fitindx = np.arange(fitsio.FITS(stackfile)['STACKINFO'].get_nrows())
                #meta = fitsio.read(stackfile, 'STACKINFO', columns=READCOLS)
                #fitindx = np.arange(len(meta))
            else:
                # We already know we like the input stackids, so no selection
                # needed.
                allstackids = fitsio.read(stackfile, 'STACKINFO', columns='STACKID')
                fitindx = np.where([tid in stackids for tid in allstackids])[0]                

            if len(fitindx) == 0:
                log.info('No requested targets found in stackfile {}'.format(stackfile))
                continue

            # Do we want just a subset of the available objects?
            if ntargets is None:
                _ntargets = len(fitindx)
            else:
                _ntargets = ntargets
            if _ntargets > len(fitindx):
                log.warning('Number of requested ntargets exceeds the number of targets on {}; reading all of them.'.format(
                    stackfile))

            __ntargets = len(fitindx)
            fitindx = fitindx[firsttarget:firsttarget+_ntargets]
            if len(fitindx) == 0:
                log.info('All {} targets in stackfile {} have been dropped (firsttarget={}, ntargets={}).'.format(
                    __ntargets, stackfile, firsttarget, _ntargets))
                continue

            # If firsttarget is a large index then the set can become empty.
            meta = Table(fitsio.read(stackfile, 'STACKINFO', rows=fitindx, columns=READCOLS))

            # Check for uniqueness.
            uu, cc = np.unique(meta['STACKID'], return_counts=True)
            if np.any(cc > 1):
                errmsg = 'Found {} duplicate STACKIDs in {}: {}'.format(
                    np.sum(cc>1), stackfile, ' '.join(uu[cc > 1].astype(str)))
                log.critical(errmsg)
                raise ValueError(errmsg)
            
            # Add some columns and append.
            meta['PHOTSYS'] = ''
            meta['SURVEY'] = survey
            meta['PROGRAM'] = program
            meta['HEALPIX'] = healpix

            self.meta.append(Table(meta))
            self.stackfiles.append(stackfile)

        if len(self.meta) == 0:
            log.warning('No targets read!')
            return

        # Now read the data as in self.read_and_unpack (for unstacked spectra).
        alldata = []
        for ispec, (stackfile, meta) in enumerate(zip(self.stackfiles, self.meta)):
            nobj = len(meta)
            if nobj == 1:
                log.info('Reading {} spectrum from {}'.format(nobj, stackfile))
            else:
                log.info('Reading {} spectra from {}'.format(nobj, stackfile))

            # Age of the universe.
            dlum = np.zeros(len(meta['Z']))
            dmod = np.zeros(len(meta['Z']))
            dlum[meta['Z'] > 0.] = self.cosmo.luminosity_distance(meta['Z'][meta['Z'] > 0.])
            dmod[meta['Z'] > 0.] = self.cosmo.distance_modulus(meta['Z'][meta['Z'] > 0.])
            tuniv = self.universe_age(meta['Z'])

            wave = fitsio.read(stackfile, 'WAVE')
            npix = len(wave)

            flux = fitsio.read(stackfile, 'FLUX')
            flux = flux[fitindx, :]
            
            ivar = fitsio.read(stackfile, 'IVAR')
            ivar = ivar[fitindx, :]
            
            # Check if the file contains a resolution matrix, if it does not
            # then use an identity matrix
            if 'RES' in fitsio.FITS(stackfile):
                res = fitsio.read(stackfile, 'RES')
                res = res[fitindx, :, :]
            else:
                res = np.ones((nobj, 1, npix)) # Hack!

            # unpack the desispec.spectra.Spectra objects into simple arrays
            unpackargs = []
            for iobj in np.arange(len(meta)):
                specdata = {
                    'uniqueid': meta[phot.uniqueid][iobj],
                    'zredrock': meta['Z'][iobj],
                    'photsys': meta['PHOTSYS'][iobj],
                    'cameras': ['brz'],
                    'dluminosity': dlum[iobj], 'dmodulus': dmod[iobj],
                    'tuniv': tuniv[iobj],
                    'wave0': [wave],
                    'flux0': [flux[iobj, :]],
                    'ivar0': [ivar[iobj, :]],
                    'mask0': [np.zeros(npix, np.int16)],
                    'res0': [Resolution(res[iobj, :, :])]
                    } 
                specdata.update({
                    'coadd_wave': specdata['wave0'][0],
                    'coadd_flux': specdata['flux0'][0],
                    'coadd_ivar': specdata['ivar0'][0],
                    'coadd_res': specdata['res0'][0],
                    })
                unpackargs.append((iobj, specdata, meta[iobj], self.fphoto, synthphot,
                                   ignore_photometry, log))
                    
            if mp > 1:
                import multiprocessing
                with multiprocessing.Pool(mp) as P:
                    out = P.map(_unpack_one_stacked_spectrum, unpackargs)
            else:
                out = [unpack_one_stacked_spectrum(*_unpackargs) for _unpackargs in unpackargs]
                
            out = list(zip(*out))
            self.meta[ispec] = Table(np.hstack(out[1]))
            alldata.append(out[0])
            del out
    
        alldata = np.concatenate(alldata)
        self.meta = vstack(self.meta)
        self.ntargets = len(self.meta)

        return alldata

    def _gather_photometry(self, specprod=None, alltiles=None):
        """Gather photometry. Unfortunately some of the bandpass information here will
        be repeated (and has to be consistent with) continuum.Fiters.

        """
        from astropy.table import vstack
        from desitarget import geomask
        from desispec.io.photo import gather_tractorphot, gather_targetphot
    
        input_meta = vstack(self.meta).copy()

        F = Photometry(fphoto=self.fphoto, load_filters=False)
        PHOTCOLS = np.unique(np.hstack((F.readcols, F.fluxcols, F.fluxivarcols)))

        # DR9 or DR10
        if hasattr(F, 'legacysurveydr'):
            from desitarget.io import releasedict
            
            legacysurveydr = F.legacysurveydr

            # targeting and Tractor columns to read from disk
            tractor = gather_tractorphot(input_meta, columns=PHOTCOLS, legacysurveydir=self.fphotodir)

            # DR9-specific stuff
            if legacysurveydr.lower() == 'dr9' or legacysurveydr.lower() == 'dr10':
                metas = []
                for meta in self.meta:
                    srt = geomask.match_to(tractor[F.uniqueid], meta[F.uniqueid])
                    assert(np.all(meta[F.uniqueid] == tractor[F.uniqueid][srt]))
                    
                    # The fibermaps in fuji and guadalupe (plus earlier productions) had a
                    # variety of errors. Fix those here using
                    # desispec.io.photo.gather_targetphot.
                    if specprod == 'fuji' or specprod == 'guadalupe': # fragile...
                        input_meta = meta[F.uniqueid, 'TARGET_RA', 'TARGET_DEC']
                        input_meta['TILEID'] = alltiles
                        targets = gather_targetphot(input_meta, fiberassign_dir=self.fiberassign_dir)
                        assert(np.all(input_meta[F.uniqueid] == targets[F.uniqueid]))
                        for col in meta.colnames:
                            if col in targets.colnames:
                                diffcol = meta[col] != targets[col]
                                if np.any(diffcol):
                                    log.warning('Updating column {} in metadata table: {}-->{}.'.format(
                                        col, meta[col][0], targets[col][0]))
                                    meta[col][diffcol] = targets[col][diffcol]
                    srt = geomask.match_to(tractor[F.uniqueid], meta[F.uniqueid])
                    assert(np.all(meta[F.uniqueid] == tractor[F.uniqueid][srt]))
                    
                    # Add the tractor catalog quantities (overwriting columns if necessary).
                    for col in tractor.colnames:
                        meta[col] = tractor[col][srt]
                        
                    # special case for some secondary and ToOs
                    I = (meta['RA'] == 0) * (meta['DEC'] == 0) * (meta['TARGET_RA'] != 0) * (meta['TARGET_DEC'] != 0)
                    if np.sum(I) > 0:
                        meta['RA'][I] = meta['TARGET_RA'][I]
                        meta['DEC'][I] = meta['TARGET_DEC'][I]
                    assert(np.all((meta['RA'] != 0) * (meta['DEC'] != 0)))
                        
                    # try to repair PHOTSYS
                    # https://github.com/desihub/fastspecfit/issues/75
                    I = np.logical_and(meta['PHOTSYS'] != 'N', meta['PHOTSYS'] != 'S') * (meta['RELEASE'] >= 9000)
                    if np.sum(I) > 0:
                        meta['PHOTSYS'][I] = [releasedict[release] if release >= 9000 else '' for release in meta['RELEASE'][I]]
                    I = np.logical_and(meta['PHOTSYS'] != 'N', meta['PHOTSYS'] != 'S')
                    if np.sum(I) > 0:
                        meta['PHOTSYS'][I] = self.resolve(meta[I])
                    I = np.logical_and(meta['PHOTSYS'] != 'N', meta['PHOTSYS'] != 'S')
                    if np.sum(I) > 0:
                        errmsg = 'Unsupported value of PHOTSYS.'
                        log.critical(errmsg)
                        raise ValueError(errmsg)

                    # placeholders (to be added in DESISpectra.read_and_unpack)
                    meta['EBV'] = np.zeros(shape=(1,), dtype='f4')
                    for band in F.bands:
                        meta['MW_TRANSMISSION_{}'.format(band.upper())] = np.ones(shape=(1,), dtype='f4')
                        
                    metas.append(meta)
        else:
            phot = Table(fitsio.read(self.fphotodir, ext=self.fphotoext, columns=PHOTCOLS))
            log.info('Read {:,d} objects from {}'.format(len(phot), self.fphotodir))

            metas = []
            for meta in self.meta:
                srt = geomask.match_to(phot[F.uniqueid], meta[F.uniqueid])
                assert(np.all(meta[F.uniqueid] == phot[F.uniqueid][srt]))
                if hasattr(F, 'dropcols'):
                    meta.remove_columns(F.dropcols)
                for col in phot.colnames:
                    meta[col] = phot[col][srt]
                # placeholders (to be added in DESISpectra.read_and_unpack)
                meta['EBV'] = np.zeros(shape=(1,), dtype='f4')
                for band in F.bands:
                    meta['MW_TRANSMISSION_{}'.format(band.upper())] = np.ones(shape=(1,), dtype='f4')
                metas.append(meta)

        return metas


_cached_linetable = None

def read_emlines(emlinesfile=None):
    """Read the set of emission lines of interest.

    """
    global _cached_linetable
    
    if _cached_linetable is not None:
        linetable = _cached_linetable
    else:
        if emlinesfile is None:
            from importlib import resources
            emlinesfile = resources.files('fastspecfit').joinpath('data/emlines.ecsv')

        try:
            linetable = Table.read(emlinesfile, format='ascii.ecsv', guess=False)
        except: 
            from desiutil.log import get_logger
            log = get_logger()
            errmsg = f'Problem reading emission lines parameter file {emlinesfile}.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        linetable.sort('restwave')
        _cached_linetable = linetable
        
    return linetable    


def init_fastspec_output(input_meta, specprod, fphoto=None, templates=None, 
                         ncoeff=None, data=None, log=None, fastphot=False, 
                         emlinesfile=None, stackfit=False):
    """Initialize the fastspecfit output data and metadata table.

    Parameters
    ----------
    tile : :class:`str`
        Tile number.
    night : :class:`str`
        Night on which `tile` was observed.
    redrock : :class:`astropy.table.Table`
        Redrock redshift table (row-aligned to `fibermap`).
    fibermap : :class:`astropy.table.Table`
        Fiber map (row-aligned to `redrock`).

    Returns
    -------


    Notes
    -----

    Must provide templates or ncoeff.

    """
    import astropy.units as u
    from astropy.table import hstack, Column
    from fastspecfit.emlines import read_emlines        

    if log is None:
        from desiutil.log import get_logger
        log = get_logger()

    linetable = read_emlines(emlinesfile=emlinesfile)

    nobj = len(input_meta)

    # get the number of templates
    if ncoeff is None:
        if not os.path.isfile(templates):
            errmsg = 'Templates file not found {}'.format(templates)
            log.critical(errmsg)
            raise IOError(errmsg)
        
        templatehdr = fitsio.read_header(templates, ext='METADATA')
        ncoeff = templatehdr['NAXIS2']

    # The information stored in the metadata table depends on which spectra
    # were fitted (exposures, nightly coadds, deep coadds).
    if stackfit:
        fluxcols = ['PHOTSYS']
    else:
        fluxcols = []
        if 'outcols' in fphoto.keys():
            fluxcols = fphoto['outcols']
        if 'fiber_bands' in fphoto.keys():
            fluxcols = np.hstack((fluxcols, ['FIBERFLUX_{}'.format(band.upper()) for band in fphoto['fiber_bands']]))
            fluxcols = np.hstack((fluxcols, ['FIBERTOTFLUX_{}'.format(band.upper()) for band in fphoto['fiber_bands']]))
        fluxcols = np.hstack((fluxcols, fphoto['fluxcols'], fphoto['fluxivarcols'], ['EBV']))
        fluxcols = np.hstack((fluxcols, ['MW_TRANSMISSION_{}'.format(band.upper()) for band in fphoto['bands']]))

    colunit = {'RA': u.deg, 'DEC': u.deg, 'EBV': u.mag}
    for fcol, icol in zip(fphoto['fluxcols'], fphoto['fluxivarcols']):
        colunit[fcol.upper()] = fphoto['photounits']
        colunit[icol.upper()] = '{}-2'.format(fphoto['photounits'])
    if 'fiber_bands' in fphoto.keys():
        for band in fphoto['fiber_bands']:
            colunit['FIBERFLUX_{}'.format(band.upper())] = fphoto['photounits']
            colunit['FIBERTOTFLUX_{}'.format(band.upper())] = fphoto['photounits']

    skipcols = np.hstack((['OBJTYPE', 'TARGET_RA', 'TARGET_DEC', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID', 'RELEASE'], fluxcols))

    if stackfit:
        redrockcols = ['Z']
    else:
        redrockcols = ['Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE', 'Z_RR', 'TSNR2_BGS',
                       'TSNR2_LRG', 'TSNR2_ELG', 'TSNR2_QSO', 'TSNR2_LYA']
    
    meta = Table()
    metacols = input_meta.colnames

    # All of this business is so we can get the columns in the order we want
    # (i.e., the order that matches the data model).
    if stackfit:
        for metacol in ['STACKID', 'SURVEY', 'PROGRAM']:
            if metacol in metacols:
                meta[metacol] = input_meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]
    else:
        for metacol in ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'NIGHT', 'FIBER',
                        'EXPID', 'TILEID_LIST', 'RA', 'DEC', 'COADD_FIBERSTATUS']:
            if metacol in metacols:
                meta[metacol] = input_meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]

        if 'SURVEY' in meta.colnames:
            if np.any(np.isin(meta['SURVEY'], 'main')) or np.any(np.isin(meta['SURVEY'], 'special')):
                TARGETINGCOLS = TARGETINGBITS['default']
            else:
                TARGETINGCOLS = TARGETINGBITS['all']
        else:
            TARGETINGCOLS = TARGETINGBITS['all']            
                
        for metacol in metacols:
            if metacol in skipcols or metacol in TARGETINGCOLS or metacol in meta.colnames or metacol in redrockcols:
                continue
            else:
                meta[metacol] = input_meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]
    
        for bitcol in TARGETINGCOLS:
            if bitcol in metacols:
                meta[bitcol] = input_meta[bitcol]
            else:
                meta[bitcol] = np.zeros(shape=(1,), dtype=np.int64)

    for redrockcol in redrockcols:
        if redrockcol in metacols: # the Z_RR from quasarnet may not be present
            meta[redrockcol] = input_meta[redrockcol]
        if redrockcol in colunit.keys():
            meta[redrockcol].unit = colunit[redrockcol]

    for fluxcol in fluxcols:
        meta[fluxcol] = input_meta[fluxcol]
        if fluxcol in colunit.keys():
            meta[fluxcol].unit = colunit[fluxcol]

    # fastspec table
    out = Table()
    if stackfit:
        for col in ['STACKID', 'SURVEY', 'PROGRAM']:
            if col in metacols:
                out[col] = input_meta[col]
    else:
        for col in ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'NIGHT', 'FIBER', 'EXPID']:
            if col in metacols:
                out[col] = input_meta[col]

    out.add_column(Column(name='Z', length=nobj, dtype='f8')) # redshift
    out.add_column(Column(name='COEFF', length=nobj, shape=(ncoeff,), dtype='f4'))

    if not fastphot:
        out.add_column(Column(name='RCHI2', length=nobj, dtype='f4'))      # full-spectrum reduced chi2
        out.add_column(Column(name='RCHI2_CONT', length=nobj, dtype='f4')) # rchi2 fitting just to the continuum (spec+phot)
    out.add_column(Column(name='RCHI2_PHOT', length=nobj, dtype='f4')) # rchi2 fitting just to the photometry

    if stackfit:
        for cam in ['BRZ']:
            out.add_column(Column(name='SNR_{}'.format(cam), length=nobj, dtype='f4')) # median S/N in each camera
        for cam in ['BRZ']:
            out.add_column(Column(name='SMOOTHCORR_{}'.format(cam), length=nobj, dtype='f4'))
    else:
        if not fastphot:
            # if the zeroth object has a fully masked camera, this data model will fail
            #if data is not None:
            #    for cam in data[0]['cameras']:
            #        out.add_column(Column(name='SNR_{}'.format(cam.upper()), length=nobj, dtype='f4')) # median S/N in each camera
            #    for cam in data[0]['cameras']:
            #        out.add_column(Column(name='SMOOTHCORR_{}'.format(cam.upper()), length=nobj, dtype='f4'))
            for cam in ['B', 'R', 'Z']:
                out.add_column(Column(name='SNR_{}'.format(cam.upper()), length=nobj, dtype='f4')) # median S/N in each camera
            for cam in ['B', 'R', 'Z']:
                out.add_column(Column(name='SMOOTHCORR_{}'.format(cam.upper()), length=nobj, dtype='f4'))
                    
    out.add_column(Column(name='VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
    if not fastphot:
        out.add_column(Column(name='VDISP_IVAR', length=nobj, dtype='f4', unit=u.second**2/u.kilometer**2))
    out.add_column(Column(name='AV', length=nobj, dtype='f4', unit=u.mag))
    out.add_column(Column(name='AGE', length=nobj, dtype='f4', unit=u.Gyr))
    out.add_column(Column(name='ZZSUN', length=nobj, dtype='f4'))
    out.add_column(Column(name='LOGMSTAR', length=nobj, dtype='f4', unit=u.solMass))
    out.add_column(Column(name='SFR', length=nobj, dtype='f4', unit=u.solMass/u.year))
    #out.add_column(Column(name='FAGN', length=nobj, dtype='f4'))
    
    if not fastphot:
        out.add_column(Column(name='DN4000', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_OBS', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_IVAR', length=nobj, dtype='f4'))
    out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))

    if not fastphot:
        # observed-frame photometry synthesized from the spectra
        for band in fphoto['synth_bands']:
            out.add_column(Column(name='FLUX_SYNTH_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies')) 
            #out.add_column(Column(name='FLUX_SYNTH_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies-2'))
        # observed-frame photometry synthesized the best-fitting spectroscopic model
        for band in fphoto['synth_bands']:
            out.add_column(Column(name='FLUX_SYNTH_SPECMODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies'))
    # observed-frame photometry synthesized the best-fitting continuum model
    for band in fphoto['bands']:
        out.add_column(Column(name='FLUX_SYNTH_PHOTMODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies'))

    for band, shift in zip(fphoto['absmag_bands'], fphoto['band_shift']):
        out.add_column(Column(name='ABSMAG{:02d}_{}'.format(int(10*shift), band.upper()), length=nobj, dtype='f4', unit=u.mag)) # absolute magnitudes
        out.add_column(Column(name='ABSMAG{:02d}_IVAR_{}'.format(int(10*shift), band.upper()), length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='KCORR{:02d}_{}'.format(int(10*shift), band.upper()), length=nobj, dtype='f4', unit=u.mag))

    for cflux in ['LOGLNU_1500', 'LOGLNU_2800']:
        out.add_column(Column(name=cflux, length=nobj, dtype='f4', unit=10**(-28)*u.erg/u.second/u.Hz))
    out.add_column(Column(name='LOGL_1450', length=nobj, dtype='f4', unit=10**(10)*u.solLum))
    out.add_column(Column(name='LOGL_1700', length=nobj, dtype='f4', unit=10**(10)*u.solLum))
    out.add_column(Column(name='LOGL_3000', length=nobj, dtype='f4', unit=10**(10)*u.solLum))
    out.add_column(Column(name='LOGL_5100', length=nobj, dtype='f4', unit=10**(10)*u.solLum))

    for cflux in ['FLYA_1215_CONT', 'FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT']:
        out.add_column(Column(name=cflux, length=nobj, dtype='f4', unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))

    if not fastphot:
        # Add chi2 metrics
        #out.add_column(Column(name='DOF', length=nobj, dtype='i8')) # full-spectrum dof
        out.add_column(Column(name='RCHI2_LINE', length=nobj, dtype='f4')) # reduced chi2 with broad line-emission
        #out.add_column(Column(name='NDOF_LINE', length=nobj, dtype='i8')) # number of degrees of freedom corresponding to rchi2_line
        #out.add_column(Column(name='DOF_BROAD', length=nobj, dtype='i8'))
        out.add_column(Column(name='DELTA_LINECHI2', length=nobj, dtype='f4')) # delta-reduced chi2 with and without broad line-emission
        out.add_column(Column(name='DELTA_LINENDOF', length=nobj, dtype=np.int32))

        # aperture corrections
        out.add_column(Column(name='APERCORR', length=nobj, dtype='f4')) # median aperture correction
        for band in fphoto['synth_bands']:        
            out.add_column(Column(name='APERCORR_{}'.format(band.upper()), length=nobj, dtype='f4'))

        out.add_column(Column(name='NARROW_Z', length=nobj, dtype='f8'))
        out.add_column(Column(name='NARROW_ZRMS', length=nobj, dtype='f8'))
        out.add_column(Column(name='BROAD_Z', length=nobj, dtype='f8'))
        out.add_column(Column(name='BROAD_ZRMS', length=nobj, dtype='f8'))
        out.add_column(Column(name='UV_Z', length=nobj, dtype='f8'))
        out.add_column(Column(name='UV_ZRMS', length=nobj, dtype='f8'))

        out.add_column(Column(name='NARROW_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='NARROW_SIGMARMS', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='BROAD_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='BROAD_SIGMARMS', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='UV_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='UV_SIGMARMS', length=nobj, dtype='f4', unit=u.kilometer / u.second))

        # special columns for the fitted doublets
        out.add_column(Column(name='MGII_DOUBLET_RATIO', length=nobj, dtype='f4'))
        out.add_column(Column(name='OII_DOUBLET_RATIO', length=nobj, dtype='f4'))
        out.add_column(Column(name='SII_DOUBLET_RATIO', length=nobj, dtype='f4'))

        for line in linetable['name']:
            line = line.upper()
            out.add_column(Column(name='{}_MODELAMP'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))
            out.add_column(Column(name='{}_AMP'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))
            out.add_column(Column(name='{}_AMP_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2))
            out.add_column(Column(name='{}_FLUX'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2)))
            out.add_column(Column(name='{}_FLUX_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4/u.erg**2))
            out.add_column(Column(name='{}_BOXFLUX'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2)))
            out.add_column(Column(name='{}_BOXFLUX_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4/u.erg**2))
            
            out.add_column(Column(name='{}_VSHIFT'.format(line), length=nobj, dtype='f4',
                                  unit=u.kilometer/u.second))
            out.add_column(Column(name='{}_SIGMA'.format(line), length=nobj, dtype='f4',
                                  unit=u.kilometer / u.second))
            
            out.add_column(Column(name='{}_CONT'.format(line), length=nobj, dtype='f4',
                                  unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))
            out.add_column(Column(name='{}_CONT_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2))
            out.add_column(Column(name='{}_EW'.format(line), length=nobj, dtype='f4',
                                  unit=u.Angstrom))
            out.add_column(Column(name='{}_EW_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=1/u.Angstrom**2))
            out.add_column(Column(name='{}_FLUX_LIMIT'.format(line), length=nobj, dtype='f4',
                                  unit=u.erg/(u.second*u.cm**2)))
            out.add_column(Column(name='{}_EW_LIMIT'.format(line), length=nobj, dtype='f4',
                                  unit=u.Angstrom))
            out.add_column(Column(name='{}_CHI2'.format(line), length=nobj, dtype='f4'))
            out.add_column(Column(name='{}_NPIX'.format(line), length=nobj, dtype=np.int32))

    # Optionally copy over some quantities of interest from the data
    # dictionary. (This step is not needed when assigning units to the
    # output tables.)
    if data is not None:
        for iobj, _data in enumerate(data):
            out['Z'][iobj] = _data['zredrock']
            if not fastphot:
                for icam, cam in enumerate(_data['cameras']):
                    out['SNR_{}'.format(cam.upper())][iobj] = _data['snr'][icam]
            if not stackfit:
                if 'fiber_bands' in fphoto.keys():
                    for iband, band in enumerate(fphoto['fiber_bands']):
                        meta['FIBERTOTFLUX_{}'.format(band.upper())][iobj] = _data['fiberphot']['nanomaggies'][iband]
                        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
                for iband, band in enumerate(fphoto['bands']):
                    meta['FLUX_{}'.format(band.upper())][iobj] = _data['phot']['nanomaggies'][iband]
                    meta['FLUX_IVAR_{}'.format(band.upper())][iobj] = _data['phot']['nanomaggies_ivar'][iband]

    return out, meta

def read_fastspecfit(fastfitfile, rows=None, columns=None, read_models=False):
    """Read the fitting results.

    """
    if os.path.isfile(fastfitfile):
        if 'FASTSPEC' in fitsio.FITS(fastfitfile):
            fastphot = False
            ext = 'FASTSPEC'
        else:
            fastphot = True
            ext = 'FASTPHOT'
            
        fastfit = Table(fitsio.read(fastfitfile, ext=ext, rows=rows, columns=columns))
        meta = Table(fitsio.read(fastfitfile, ext='METADATA', rows=rows, columns=columns))
        if read_models and ext == 'FASTSPEC':
            models = fitsio.read(fastfitfile, ext='MODELS')
            if rows is not None:
                models = models[rows, :, :]
        else:
            models = None
        log.info('Read {:,d} object(s) from {}'.format(len(fastfit), fastfitfile))

        # Add specprod to the metadata table so that we can stack across
        # productions (e.g., Fuji+Guadalupe).
        hdr = fitsio.read_header(fastfitfile, ext=0)#, ext='PRIMARY')

        if 'SPECPROD' in hdr:
            specprod = hdr['SPECPROD']
            meta['SPECPROD'] = specprod
            
        if 'COADDTYP' in hdr:
            coadd_type = hdr['COADDTYP']
        else:
            coadd_type = None

        if read_models:
            return fastfit, meta, coadd_type, fastphot, models
        else:
            return fastfit, meta, coadd_type, fastphot
    
    else:
        log.warning('File {} not found.'.format(fastfitfile))
        if read_models:
            return [None]*5
        else:
            return [None]*4

def write_fastspecfit(out, meta, modelspectra=None, outfile=None, specprod=None,
                      coadd_type=None, fphotofile=None, templates=None,
                      emlinesfile=None, fastphot=False, inputz=False,
                      no_smooth_continuum=False, ignore_photometry=False, broadlinefit=True,
                      use_quasarnet=True, constrain_age=False, verbose=True):
    """Write out.

    """
    import gzip, shutil
    from astropy.io import fits
    from desispec.io.util import fitsheader
    from desiutil.depend import add_dependencies, possible_dependencies, setdep

    t0 = time.time()
    outdir = os.path.dirname(os.path.abspath(outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    nobj = len(out)
    if nobj == 1:
        log.info('Writing {} object to {}'.format(nobj, outfile))
    else:
        log.info('Writing {:,d} objects to {}'.format(nobj, outfile))
    
    if outfile.endswith('.gz'):
        tmpfile = outfile[:-3]+'.tmp'
    else:
        tmpfile = outfile+'.tmp'

    if fastphot:
        extname = 'FASTPHOT'
    else:
        extname = 'FASTSPEC'

    out.meta['EXTNAME'] = extname
    meta.meta['EXTNAME'] = 'METADATA'

    primhdr = []
    if specprod:
        primhdr.append(('EXTNAME', 'PRIMARY'))
        primhdr.append(('SPECPROD', (specprod, 'spectroscopic production name')))
    if coadd_type:
        primhdr.append(('COADDTYP', (coadd_type, 'spectral coadd type')))
    primhdr.append(('INPUTZ', (inputz is True, 'input redshifts provided')))
    primhdr.append(('NOSCORR', (no_smooth_continuum is True, 'no smooth continuum correction')))
    primhdr.append(('NOPHOTO', (ignore_photometry is True, 'no fitting to photometry')))
    primhdr.append(('BRDLFIT', (broadlinefit is True, 'carry out broad-line fitting')))
    primhdr.append(('CONSAGE', (constrain_age is True, 'constrain SPS ages')))
    primhdr.append(('USEQNET', (use_quasarnet is True, 'use QuasarNet redshifts')))

    primhdr = fitsheader(primhdr)
    add_dependencies(primhdr, module_names=possible_dependencies+['fastspecfit'],
                     envvar_names=['DESI_ROOT', 'DUST_DIR', 'FTEMPLATES_DIR', 'FPHOTO_DIR'])
    if fphotofile:
        setdep(primhdr, 'FPHOTO_FILE', str(fphotofile))
    if templates:
        setdep(primhdr, 'FTEMPLATES_FILE', os.path.basename(templates))
    if emlinesfile:
        setdep(primhdr, 'EMLINES_FILE', str(emlinesfile))

    hdus = fits.HDUList()
    hdus.append(fits.PrimaryHDU(None, primhdr))
    hdus.append(fits.convenience.table_to_hdu(out))
    hdus.append(fits.convenience.table_to_hdu(meta))

    if modelspectra is not None:
        hdu = fits.ImageHDU(name='MODELS')
        # [nobj, 3, nwave]
        hdu.data = np.swapaxes(np.array([modelspectra['CONTINUUM'].data,
                                         modelspectra['SMOOTHCONTINUUM'].data,
                                         modelspectra['EMLINEMODEL'].data]), 0, 1)
        for key in modelspectra.meta.keys():
            hdu.header[key] = (modelspectra.meta[key][0], modelspectra.meta[key][1]) # all the spectra are identical, right??
                
        hdus.append(hdu)
        
    hdus.writeto(tmpfile, overwrite=True, checksum=True)

    # compress if needed (via another tempfile), otherwise just rename
    if outfile.endswith('.gz'):
        tmpfilegz = outfile[:-3]+'.tmp.gz'
        with open(tmpfile, 'rb') as f_in:
            with gzip.open(tmpfilegz, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.rename(tmpfilegz, outfile)
        os.remove(tmpfile)
    else:
        os.rename(tmpfile, outfile)

    if verbose:
        log.info('Writing out took {:.2f} seconds.'.format(time.time()-t0))

def select(fastfit, metadata, coadd_type, healpixels=None, tiles=None,
           nights=None, return_index=False):
    """Optionally trim to a particular healpix or tile and/or night."""
    keep = np.ones(len(fastfit), bool)
    if coadd_type == 'healpix':
        if healpixels:
            pixelkeep = np.zeros(len(fastfit), bool)
            for healpixel in healpixels:
                pixelkeep = np.logical_or(pixelkeep, metadata['HEALPIX'].astype(str) == healpixel)
            keep = np.logical_and(keep, pixelkeep)
            log.info('Keeping {:,d} objects from healpixels(s) {}'.format(len(fastfit), ','.join(healpixels)))
    else:
        if tiles:
            tilekeep = np.zeros(len(fastfit), bool)
            for tile in tiles:
                tilekeep = np.logical_or(tilekeep, metadata['TILEID'].astype(str) == tile)
            keep = np.logical_and(keep, tilekeep)
            log.info('Keeping {:,d} objects from tile(s) {}'.format(len(fastfit), ','.join(tiles)))
        if nights and 'NIGHT' in metadata:
            nightkeep = np.zeros(len(fastfit), bool)
            for night in nights:
                nightkeep = np.logical_or(nightkeep, metadata['NIGHT'].astype(str) == night)
            keep = np.logical_and(keep, nightkeep)
            log.info('Keeping {:,d} objects from night(s) {}'.format(len(fastfit), ','.join(nights)))
            
    if return_index:
        return np.where(keep)[0]
    else:
        return fastfit[keep], metadata[keep]

def get_templates_filename(templateversion=DEFAULT_TEMPLATEVERSION, imf=DEFAULT_IMF):
    """Get the templates filename. """
    from fastspecfit.io import FTEMPLATES_DIR_NERSC
    templates_dir = os.path.expandvars(os.environ.get('FTEMPLATES_DIR', FTEMPLATES_DIR_NERSC))
    templates = os.path.join(templates_dir, templateversion, f'ftemplates-{imf}-{templateversion}.fits')
    return templates


def get_qa_filename(metadata, coadd_type, outprefix=None, outdir=None,
                    fastphot=False, log=None):
    """Build the QA filename.

    """
    import astropy

    if log is None:
        from desiutil.log import get_logger
        log = get_logger()

    if outdir is None:
        outdir = '.'
        
    if outprefix is None:
        if fastphot:
            outprefix = 'fastphot'
        else:
            outprefix = 'fastspec'

    def _one_filename(_metadata):
        if coadd_type == 'healpix':
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                outprefix, _metadata['SURVEY'], _metadata['PROGRAM'],
                _metadata['HEALPIX'], _metadata['TARGETID']))
        elif coadd_type == 'cumulative':
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                outprefix, _metadata['TILEID'], coadd_type, _metadata['TARGETID']))
        elif coadd_type == 'pernight':
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                outprefix, _metadata['TILEID'], _metadata['NIGHT'], _metadata['TARGETID']))
        elif coadd_type == 'perexp':
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                outprefix, _metadata['TILEID'], _metadata['NIGHT'],
                _metadata['EXPID'], _metadata['TARGETID']))
        elif coadd_type == 'custom':
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                outprefix, _metadata['SURVEY'], _metadata['PROGRAM'],
                _metadata['HEALPIX'], _metadata['TARGETID']))
        elif coadd_type == 'stacked':
            pngfile = os.path.join(outdir, '{}-{}-{}.png'.format(
                outprefix, coadd_type, _metadata['STACKID']))
        else:
            errmsg = 'Unrecognized coadd_type {}!'.format(coadd_type)
            log.critical(errmsg)
            raise ValueError(errmsg)
        return pngfile

    if type(metadata) is astropy.table.row.Row or type(metadata) is np.void:
        pngfile = _one_filename(metadata)
    else:
        pngfile = [_one_filename(_metadata) for _metadata in metadata]
    
    return pngfile


def cache_templates(templates=None, templateversion=DEFAULT_TEMPLATEVERSION, imf=DEFAULT_IMF,
                    mintemplatewave=None, maxtemplatewave=40e4, vdisp_nominal=125.,
                    read_linefluxes=False, fastphot=False, log=None):
    """"Read the templates into a dictionary.
    
    Parameters
    ----------
    templates : :class:`str`
        Full path to the templates to read. Defaults to the output of
        :class:`get_templates_filename`.
    templateversion : :class:`str`
        Version of the templates. Defaults to :class:`DEFAULT_TEMPLATEVERSION`.
    mapdir : :class:`str`, optional
        Full path to the Milky Way dust maps.
    mintemplatewave : :class:`float`, optional, defaults to None
        Minimum template wavelength to read into memory. If ``None``, the minimum
        available wavelength is used (around 100 Angstrom).
    maxtemplatewave : :class:`float`, optional, defaults to 6e4
        Maximum template wavelength to read into memory.

    """
    import fitsio
    from fastspecfit.continuum import convolve_vdisp, PIXKMS_WAVESPLIT, PIXKMS_BLU
    
    if log is None:
        from desiutil.log import get_logger
        log = get_logger()

    if templates is None:
        templates = get_templates_filename(templateversion=templateversion, imf=imf)
        
    if not os.path.isfile(templates):
        errmsg = f'Templates file {templates} not found.'
        log.critical(errmsg)
        raise IOError(errmsg)

    T = fitsio.FITS(templates)

    #log.info('Reading {}'.format(templates))
    wave = T['WAVE'].read() # [npix]
    wavehdr = T['WAVE'].read_header() # [npix]
    templateflux = T['FLUX'].read()  # [npix,nsed]
    templatelineflux = T['LINEFLUX'].read()  # [npix,nsed]
    templateinfo = T['METADATA'].read()
    templatehdr = T['METADATA'].read_header()
    vdisphdr = T['VDISPFLUX'].read_header()

    if mintemplatewave is None:
        mintemplatewave = np.min(wave)
    
    keeplo = np.searchsorted(wave, mintemplatewave, 'left')
    keephi = np.searchsorted(wave, maxtemplatewave, 'right')
    
    templatewave = wave[keeplo:keephi]
    templateflux = templateflux[keeplo:keephi, :]
    templateflux_nolines = templateflux - templatelineflux[keeplo:keephi, :]
    
    # Cache a copy of the line-free templates at the nominal velocity
    # dispersion (needed for fastphot as well).
    if 'VDISPNOM' in vdisphdr: # older templates do not have this header card
        vdisp_nominal = vdisphdr['VDISPNOM'] # [km/s]

    hi = np.searchsorted(templatewave, PIXKMS_WAVESPLIT, 'left')

    templateflux_nolines_nomvdisp = convolve_vdisp(
        templateflux_nolines, vdisp_nominal, limit=hi,
        pixsize_kms=PIXKMS_BLU)
    
    templateflux_nomvdisp = convolve_vdisp(
        templateflux, vdisp_nominal, limit=hi,
        pixsize_kms=PIXKMS_BLU)
    
    # pack into a dictionary
    templatecache = {'imf': templatehdr['IMF'],
                     #'nsed': len(templateinfo), 'npix': len(wavekeep),
                     'vdisp_nominal': vdisp_nominal,
                     'templateinfo': Table(templateinfo),
                     'templatewave': templatewave,
                     'templateflux': templateflux,
                     'templateflux_nomvdisp': templateflux_nomvdisp,
                     'templateflux_nolines': templateflux_nolines,
                     'templateflux_nolines_nomvdisp': templateflux_nolines_nomvdisp,
                     }

    # maintain backwards compatibility with older templates (versions <2.0.0)
    if 'VDISPMIN' in vdisphdr:
        templatecache['oldtemplates'] = True
    else:
        templatecache['oldtemplates'] = False

    if not fastphot:
        
        if templatecache['oldtemplates']:
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

            templatecache.update({
                'vdispflux': vdispflux, 
                'vdispwave': vdispwave,
                'vdisp': vdisp, 
                'vdisp_nominal_indx': np.where(vdisp == vdisp_nominal)[0]
            })

    if not templatecache['oldtemplates']:
        if 'DUSTFLUX' in T and 'AGNFLUX' in T:
            dustflux = T['DUSTFLUX'].read()
            dusthdr = T['DUSTFLUX'].read_header()

            agnflux = T['AGNFLUX'].read()
            agnhdr = T['AGNFLUX'].read_header()

            # make sure both dustflux and agnflux are normalized to unity
            from scipy.integrate import simpson
            dustflux /= simpson(dustflux, x=wave) # should already be 1.0
            agnflux /= simpson(agnflux, x=wave) # should already be 1.0

            templatecache.update({'dustflux': dustflux[keeplo:keephi], 
                                  'qpah': dusthdr['QPAH'],
                                  'umin': dusthdr['umin'],
                                  'gamma': dusthdr['GAMMA'],
                                  'agnflux': agnflux[keeplo:keephi], 
                                  'agntau': agnhdr['AGNTAU'],
                                  })
        else:
            errmsg = f'Templates file {templates} missing mandatory extensions DUSTFLUX and AGNFLUX.'
            log.critical(errmsg)
            raise IOError(errmsg)

    # Read the model emission-line fluxes; only present for
    # templateversion>=1.1.1 and generally only useful to a power-user.
    if read_linefluxes:
        templatecache.update({
            'linefluxes': T['LINEFLUXES'].read(),
            'linewaves':  T['LINEWAVES'].read()
        })

    return templatecache


def one_desi_spectrum(survey, program, healpix, targetid, specprod='fuji',
                      outdir='.', overwrite=False):
    """Utility function to write a single DESI spectrum (e.g., for paper figures or
    unit tests).

    """
    from redrock.external.desi import write_zbest
    from desispec.io import write_spectra, read_spectra
    from fastspecfit.qa import fastqa
    from fastspecfit.fastspecfit import fastspec

    os.environ['SPECPROD'] = specprod # needed to get write_spectra have the correct dependency

    specdir = os.path.join(os.environ.get('DESI_ROOT'), 'spectro', 'redux', specprod, 'healpix',
                           survey, program, str(healpix//100), str(healpix))
    coaddfile = os.path.join(specdir, f'coadd-{survey}-{program}-{healpix}.fits')
    redrockfile = os.path.join(specdir, f'redrock-{survey}-{program}-{healpix}.fits')

    out_coaddfile = os.path.join(outdir, f'coadd-{survey}-{program}-{healpix}-{targetid}.fits')
    out_redrockfile = os.path.join(outdir, f'redrock-{survey}-{program}-{healpix}-{targetid}.fits')
    out_fastfile = os.path.join(outdir, f'fastspec-{survey}-{program}-{healpix}-{targetid}.fits')
    if (os.path.isfile(out_coaddfile) or os.path.isfile(out_redrockfile) or os.path.isfile(out_fastfile)) and not overwrite:
        if os.path.isfile(out_coaddfile):
            print(f'Coadd file {out_coaddfile} exists and overwrite is False')
        if os.path.isfile(out_redrockfile):
            print(f'Redrock file {out_redrockfile} exists and overwrite is False')
        if os.path.isfile(out_fastfile):
            print(f'fastspec file {out_fastfile} exists and overwrite is False')
        return
    
    redhdr = fitsio.read_header(redrockfile)
    zbest = Table.read(redrockfile, 'REDSHIFTS')
    fibermap = Table.read(redrockfile, 'FIBERMAP')
    expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
    tsnr2 = Table.read(redrockfile, 'TSNR2')
    
    spechdr = fitsio.read_header(coaddfile)
    
    zbest = zbest[np.isin(zbest['TARGETID'], targetid)]
    fibermap = fibermap[np.isin(fibermap['TARGETID'], targetid)]
    expfibermap = expfibermap[np.isin(expfibermap['TARGETID'], targetid)]
    tsnr2 = tsnr2[np.isin(tsnr2['TARGETID'], targetid)]
    
    archetype_version = None
    template_version = {redhdr['TEMNAM{:02d}'.format(nn)]: redhdr['TEMVER{:02d}'.format(nn)] for nn in np.arange(10)}

    print(f'Writing {out_redrockfile}')
    write_zbest(out_redrockfile, zbest, fibermap, expfibermap, tsnr2,
                template_version, archetype_version, spec_header=spechdr)
    
    spec = read_spectra(coaddfile).select(targets=targetid)
    print(f'Writing {out_coaddfile}')
    write_spectra(out_coaddfile, spec)

    fastspec(args=f'{out_redrockfile} -o {out_fastfile}'.split())
    cmdargs = f'{out_fastfile} --redrockfiles {out_redrockfile} -o {outdir}'
    if overwrite:
        cmdargs += '--overwrite'
    fastqa(args=cmdargs.split())
