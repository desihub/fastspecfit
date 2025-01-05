"""
fastspecfit.io
==============

Tools for reading DESI spectra and reading and writing fastspecfit files.

"""
import os, time
import numpy as np

import fitsio
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.singlecopy import sc_data
from fastspecfit.photometry import Photometry
from fastspecfit.util import FLUXNORM, ZWarningMask


# list of all possible targeting bit columns
TARGETINGBITS = {
    'all': ('CMX_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',
             'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET',
             'SV2_DESI_TARGET', 'SV2_BGS_TARGET', 'SV2_MWS_TARGET',
             'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET',
             'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET'),
    'fuji': ('CMX_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',
             'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET',
             'SV2_DESI_TARGET', 'SV2_BGS_TARGET', 'SV2_MWS_TARGET',
             'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET',
             'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET'),
    'default': ('DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET'),
    }

# fibermap and exp_fibermap columns to read
FMCOLS = ('TARGETID', 'TARGET_RA', 'TARGET_DEC', 'COADD_FIBERSTATUS', 'OBJTYPE',
          'PHOTSYS', 'RELEASE', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID')

EXPFMCOLS = {
    'perexp':     ('TARGETID', 'TILEID', 'FIBER', 'EXPID'),
    'pernight':   ('TARGETID', 'TILEID', 'FIBER'),
    'cumulative': ('TARGETID', 'TILEID', 'FIBER'),
    'healpix':    ('TARGETID', 'TILEID'), # tileid will be an array
    'custom':     ('TARGETID', 'TILEID'), # tileid will be an array
    }

# redshift columns to read
REDSHIFTCOLS = ('TARGETID', 'Z', 'ZWARN', 'SPECTYPE', 'SUBTYPE', 'DELTACHI2')

# tsnr columns to read
TSNR2COLS = ('TSNR2_BGS', 'TSNR2_LRG', 'TSNR2_ELG', 'TSNR2_QSO', 'TSNR2_LYA')

# quasarnet and MgII afterburner columns to read
QNLINES = ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha', ]
QNCOLS = ['TARGETID', 'Z_NEW', 'ZWARN_NEW', 'IS_QSO_QN_NEW_RR', ] + QNLINES
MGIICOLS = ['TARGETID', 'IS_QSO_MGII']

def one_spectrum(specdata, meta, uncertainty_floor=0.01, RV=3.1,
                 init_sigma_uv=None, init_sigma_narrow=None,
                 init_sigma_balmer=None, init_vshift_uv=None,
                 init_vshift_narrow=None, init_vshift_balmer=None,
                 fastphot=False, synthphot=True, debug_plots=False):
    """Pre-process the data for a single object.

    """
    from fastspecfit.util import mwdust_transmission

    phot = sc_data.photometry
    filters = phot.filters[specdata['photsys']]

    # Process the total fluxes, correcting for MW dust.
    maggies = np.zeros(len(phot.bands))
    ivarmaggies = np.zeros(len(phot.bands))
    for iband, (band, fluxcol, ivarcol) in enumerate(zip(phot.bands, phot.fluxcols, phot.fluxivarcols)):
        transmission = meta[f'MW_TRANSMISSION_{band.upper()}']
        maggies[iband] = meta[fluxcol.upper()] / transmission
        ivarmaggies[iband] = meta[ivarcol.upper()] * transmission**2

    if not np.all(ivarmaggies >= 0.):
        errmsg = 'Some ivarmaggies are negative!'
        log.critical(errmsg)
        raise ValueError(errmsg)

    specdata['photometry'] = Photometry.parse_photometry(
        phot.bands, maggies=maggies, ivarmaggies=ivarmaggies,
        nanomaggies=True, lambda_eff=filters.effective_wavelengths.value,
        min_uncertainty=phot.min_uncertainty)

    # Optionally add the fiber photometry; note that the transmission
    # factors were computed in DESISpectra.read.
    if hasattr(phot, 'fiber_filters'):
        fiber_filters = phot.fiber_filters[specdata['photsys']]

        mw_transmission_fiberflux = specdata['mw_transmission_fiberflux']

        fibermaggies = np.zeros(len(phot.fiber_bands))
        fibertotmaggies = np.zeros(len(phot.fiber_bands))
        #ivarfibermaggies = np.zeros(len(phot.fiber_bands))

        for iband, band in enumerate(phot.fiber_bands):
            band = band.upper()
            fibermaggies[iband] = meta[f'FIBERFLUX_{band}'] / mw_transmission_fiberflux[iband]
            fibertotmaggies[iband] = meta[f'FIBERTOTFLUX_{band}'] / mw_transmission_fiberflux[iband]

        lambda_eff=fiber_filters.effective_wavelengths.value
        specdata['fiberphot'] = Photometry.parse_photometry(phot.fiber_bands,
                                                            maggies=fibermaggies,
                                                            nanomaggies=True,
                                                            lambda_eff=lambda_eff)
        specdata['fibertotphot'] = Photometry.parse_photometry(phot.fiber_bands,
                                                               maggies=fibertotmaggies,
                                                               nanomaggies=True,
                                                               lambda_eff=lambda_eff)

    if not fastphot:
        from desiutil.dust import dust_transmission
        from fastspecfit.util import median, ivar2var
        from fastspecfit.resolution import Resolution
        from fastspecfit.linemasker import LineMasker

        # now process the spectroscopy
        specdata.update({
            'wave': [],
            'flux': [],
            'ivar': [],
            'mask': [],
            'res': [],
            'snr': np.zeros(len(np.atleast_1d(specdata['cameras'])), 'f4'),
            'linemask': [],
            'linepix': [],
        })

        cameras, npixpercamera = [], []
        for icam, camera in enumerate(specdata['cameras']):
            # Check whether the camera is fully masked.
            if np.sum(specdata['ivar0'][icam]) == 0:
                log.warning(f'Dropping fully masked camera {camera} [{specdata["uniqueid"]}].')
            else:
                ivar = specdata['ivar0'][icam]
                mask = specdata['mask0'][icam].astype(bool)

                # always mask the first and last XX pixels
                mask[:3] = True
                mask[-3:] = True

                # In the pipeline, if mask!=0 that does not mean ivar==0, but we
                # want to be more aggressive about masking here.
                ivar[mask] = 0.

                if np.all(ivar == 0.):
                    log.warning(f'Dropping fully masked camera {camera} [{specdata["uniqueid"]}].')
                else:
                    res = specdata['res0'][icam]
                    ## interpolate over pixels where the resolution matrix is masked
                    #if np.any(mask):
                    #    J = np.where(np.logical_not(mask))[0]
                    #    I = np.where(mask)[0]
                    #    for irow in range(res.shape[0]):
                    #        res[irow, I] = np.interp(I, J, res[irow, J])

                    # should we also interpolate over the coadded resolution matrix??

                    # include the minimum uncertainty in quadrature with the input ivar
                    minvar = (uncertainty_floor * specdata['flux0'][icam])**2
                    var, I = ivar2var(ivar)
                    newivar = np.zeros_like(ivar)
                    newivar[I] = 1. / (minvar[I] + var[I])
                    ivar = newivar

                    cameras.append(camera)
                    npixpercamera.append(len(specdata['wave0'][icam])) # number of pixels in this camera

                    # Compute the SNR before we correct for dust.
                    specdata['snr'][icam] = median(specdata['flux0'][icam] * np.sqrt(ivar))

                    #mw_transmission_spec = 10**(-0.4 * ebv * RV * ext_odonnell(wave[camera], Rv=RV))
                    mw_transmission_spec = dust_transmission(specdata['wave0'][icam], meta['EBV'], Rv=RV)
                    specdata['flux'].append(specdata['flux0'][icam] / mw_transmission_spec)
                    specdata['ivar'].append(ivar * mw_transmission_spec**2)
                    specdata['wave'].append(specdata['wave0'][icam])
                    specdata['mask'].append(mask)
                    specdata['res'].append(Resolution(res))

        if len(cameras) == 0:
            errmsg = 'No good data, which should never happen.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        # clean up unused items in data dictionary and
        # freeze lists that will not be further modified
        for key in ('wave', 'flux', 'ivar', 'mask', 'res'):
            del specdata[key + '0']
            specdata[key] = tuple(specdata[key])

        # Pre-compute some convenience variables for "un-hstacking"
        # an "hstacked" spectrum.
        specdata['cameras'] = np.array(cameras)
        specdata['npixpercamera'] = np.array(npixpercamera)

        ncam = len(specdata['cameras'])
        c_ends   = np.cumsum(specdata['npixpercamera'])
        c_starts = c_ends - specdata['npixpercamera']
        specdata['camerapix'] = np.zeros((ncam, 2), np.int32)
        specdata['camerapix'][:, 0] = c_starts
        specdata['camerapix'][:, 1] = c_ends

        # use the coadded spectrum to build a robust emission-line mask
        LM = LineMasker(sc_data.emlines.table)
        pix = LM.build_linemask(
            specdata['coadd_wave'], specdata['coadd_flux'],
            specdata['coadd_ivar'], specdata['coadd_res'],
            uniqueid=specdata['uniqueid'], redshift=specdata['redshift'],
            initsigma_broad=init_sigma_uv,
            initsigma_narrow=init_sigma_narrow,
            initsigma_balmer_broad=init_sigma_balmer,
            initvshift_broad=init_vshift_uv,
            initvshift_narrow=init_vshift_narrow,
            initvshift_balmer_broad=init_vshift_balmer,
            debug_plots=debug_plots)

        # Map the pixels belonging to individual emission lines onto the
        # original per-camera spectra. This works, but maybe there's a better
        # way to do it?
        for icam in range(ncam):
            camlinepix = {}
            camlinemask = np.zeros(specdata['npixpercamera'][icam], bool)
            for linename in pix['coadd_linepix']:
                linepix = pix['coadd_linepix'][linename]
                # if the line is entirely off this camera, skip it
                oncam = ((specdata["coadd_wave"][linepix] >= np.min(specdata['wave'][icam])) &
                         (specdata["coadd_wave"][linepix] <= np.max(specdata['wave'][icam])))
                if not np.any(oncam):
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
        if synthphot:
            synth_filters = phot.synth_filters[specdata['photsys']]
            synthmaggies = Photometry.get_ab_maggies(
                synth_filters, specdata['coadd_flux'] / FLUXNORM, specdata['coadd_wave'])

            specdata['synthphot'] = Photometry.parse_photometry(
                phot.synth_bands, maggies=synthmaggies, nanomaggies=False,
                lambda_eff=synth_filters.effective_wavelengths.value)

    return specdata


def one_stacked_spectrum(specdata, meta, synthphot=True, debug_plots=False):
    """Unpack the data for a single stacked spectrum.

    """
    from fastspecfit.linemasker import LineMasker
    from fastspecfit.util import median

    phot = sc_data.photometry

    filters = phot.filters[specdata['photsys']]
    synth_filters = phot.synth_filters[specdata['photsys']]

    # Dummy imaging photometry.
    maggies = np.zeros(len(phot.bands))
    ivarmaggies = np.zeros(len(phot.bands))

    specdata['photometry'] = Photometry.parse_photometry(
        phot.bands, maggies=maggies, ivarmaggies=ivarmaggies,
        nanomaggies=True, lambda_eff=filters.effective_wavelengths.value,
        min_uncertainty=phot.min_uncertainty)

    specdata.update({
        'wave': [],
        'flux': [],
        'ivar': [],
        'mask': [],
        'res': [],
        'snr': np.zeros(1, 'f4'),
        'linemask': [],
        'linepix': [],
    })

    cameras, npixpercamera = [], []
    for icam, camera in enumerate(specdata['cameras']):
        # Check whether the camera is fully masked.
        if np.sum(specdata['ivar0'][icam]) == 0:
            log.warning(f'Dropping fully masked camera {camera} [{specdata["uniqueid"]}].')
        else:
            ivar = specdata['ivar0'][icam]
            mask = specdata['mask0'][icam]

            # always mask the first and last XX pixels
            mask[:3] = True
            mask[-3:] = True

            # In the pipeline, if mask!=0 that does not mean ivar==0, but we
            # want to be more aggressive about masking here.
            ivar[mask] = 0.

            if np.all(ivar == 0.):
                log.warning(f'Dropping fully masked camera {camera} [{specdata["uniqueid"]}].')
            else:
                cameras.append(camera)
                npixpercamera.append(len(specdata['wave0'][icam])) # number of pixels in this camera

                specdata['snr'][icam] = median(specdata['flux0'][icam] * np.sqrt(ivar))
                specdata['flux'].append(specdata['flux0'][icam])
                specdata['ivar'].append(ivar)
                specdata['wave'].append(specdata['wave0'][icam])
                specdata['mask'].append(mask)
                specdata['res'].append(specdata['res0'][icam])

    if len(cameras) == 0:
        errmsg = 'No good data, which should never happen.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    # clean up unused items in data dictionary and
    # freeze lists that will not be further modified
    for key in ('wave', 'flux', 'ivar', 'mask', 'res'):
        del specdata[key + '0']
        specdata[key] = tuple(specdata[key])

    # Pre-compute some convenience variables for "un-hstacking"
    # an "hstacked" spectrum.
    specdata['cameras'] = np.array(cameras)
    specdata['npixpercamera'] = np.array(npixpercamera)

    ncam = len(specdata['cameras'])
    c_ends   = np.cumsum(specdata['npixpercamera'])
    c_starts = c_ends - specdata['npixpercamera']
    specdata['camerapix'] = np.zeros((ncam, 2), np.int32)
    specdata['camerapix'][:, 0] = c_starts
    specdata['camerapix'][:, 1] = c_ends

    LM = LineMasker(sc_data.emlines.table)
    pix = LM.build_linemask(
        specdata['coadd_wave'], specdata['coadd_flux'],
        specdata['coadd_ivar'], specdata['coadd_res'],
        uniqueid=specdata['uniqueid'], redshift=specdata['redshift'],
        debug_plots=debug_plots)

    # Map the pixels belonging to individual emission lines onto the
    # original per-camera spectra. This works, but maybe there's a better
    # way to do it?
    for icam in range(ncam):
        camlinepix = {}
        camlinemask = np.zeros(specdata['npixpercamera'][icam], bool)
        for linename in pix['coadd_linepix']:
            linepix = pix['coadd_linepix'][linename]
            # if the line is entirely off this camera, skip it
            oncam = ((specdata["coadd_wave"][linepix] >= np.min(specdata['wave'][icam])) &
                     (specdata["coadd_wave"][linepix] <= np.max(specdata['wave'][icam])))
            if not np.any(oncam):
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
    if synthphot:
        synthmaggies = Photometry.get_ab_maggies(
            synth_filters, specdata['coadd_flux'] / FLUXNORM, specdata['coadd_wave'])

        specdata['synthphot'] = Photometry.parse_photometry(
            phot.synth_bands, maggies=synthmaggies, nanomaggies=False,
            lambda_eff=synth_filters.effective_wavelengths.value)

    return specdata


class DESISpectra(object):
    def __init__(self, phot, cosmo, redux_dir=None, fphotodir=None, mapdir=None):
        """Class to read in DESI spectra and associated metadata.

        Parameters
        ----------
        redux_dir : str
            Full path to the location of the reduced DESI data. Optional and
            defaults to `$DESI_SPECTRO_REDUX`.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.

        """
        if redux_dir is None:
            if not 'DESI_SPECTRO_REDUX' in os.environ:
                errmsg = "'DESI_SPECTRO_REDUX' environment variable or redux_dir must be set"
                log.critical(errmsg)
                raise KeyError(errmsg)
            self.redux_dir = os.path.expandvars(os.environ.get('DESI_SPECTRO_REDUX'))
        else:
            self.redux_dir = os.path.expandvars(redux_dir)

        if fphotodir is None:
            self.fphotoext = None
            self.fphotodir = os.path.expandvars(os.environ.get('FPHOTO_DIR'))
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

        if mapdir is None:
            self.mapdir = os.path.join(os.path.expandvars(os.environ.get('DUST_DIR')), 'maps')
        else:
            self.mapdir = mapdir

        self.phot = phot
        self.cosmo = cosmo


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
        arean = ((targets["DEC"] >= split) & galn)

        # ADM retain 'N' targets in 'N' area and 'S' in 'S' area.
        #keep = (photn & arean) | (~photn & ~arean)
        #return targets[keep]

        inorth = (photn & arean)
        newphotsys = np.array(['S'] * len(targets))
        newphotsys[inorth] = 'N'

        return newphotsys


    def gather_metadata(self, redrockfiles, zmin=None, zmax=None, zwarnmax=None,
                        targetids=None, firsttarget=0, ntargets=None,
                        input_redshifts=None, specprod_dir=None, use_quasarnet=True,
                        redrockfile_prefix='redrock-', specfile_prefix='coadd-',
                        qnfile_prefix='qso_qn-', mgiifile_prefix='qso_mgii-'):
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
        mgiifile_prefix : str
            Prefix of the MgII afterburner file. Defaults to `qso_mgii-`.

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

        redrockfiles = set(redrockfiles)
        #log.info(f'Reading and parsing {len(redrockfiles)} unique redrockfile(s).')

        alltiles = []
        self.redrockfiles, self.specfiles, self.meta, self.surveys = [], [], [], []

        t0 = time.time()
        for ired, redrockfile in enumerate(redrockfiles):
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
            mgiifile = redrockfile.replace(redrockfile_prefix, mgiifile_prefix)
            if os.path.isfile(qnfile) and os.path.isfile(mgiifile) and use_quasarnet and input_redshifts is None:
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
                    infofile = os.path.join(specprod_dir, f'tiles-{self.specprod}.csv')
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
            allfmcols = set(fitsio.FITS(specfile)['FIBERMAP'].get_colnames())
            READFMCOLS = list(FMCOLS) + [col for col in TARGETINGCOLS if col in allfmcols]

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
                    # redshift cut below, after we correct the redshift.
                    fitindx = np.where((zb['Z'] < zmax) &
                                       (meta['OBJTYPE'] == 'TGT') & (zb['ZWARN'] <= zwarnmax) &
                                       (zb['ZWARN'] & ZWarningMask.NODATA == 0))[0]
                else:
                    fitindx = np.where((zb['Z'] > zmin) & (zb['Z'] < zmax) &
                                       (meta['OBJTYPE'] == 'TGT') & (zb['ZWARN'] <= zwarnmax) &
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
            # for QSO targets**.
            zb['Z_RR'] = zb['Z'] # add it at the end
            #zb['ZERR_RR'] = zb['ZERR']
            zb['ZWARN_RR'] = zb['ZWARN']
            if use_qn:
                self.update_qso_redshifts(zb, meta, qnfile, mgiifile, fitindx, self.specprod)
                # now apply zmin
                keep = (zb['Z'] > zmin)
                if not np.any(keep):
                    log.info(f'No requested targets found in redrockfile {redrockfile}')
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
            if not np.any(I):
                errmsg = 'No matching targets in exposure table.'
                log.critical(errmsg)
                raise ValueError(errmsg)
            expmeta = Table(expmeta[I])

            # build the list of tiles that went into each unique target / coadd
            tileid_list = [] # variable length, so need to build the array first
            for tid in meta['TARGETID']:
                I = (tid == expmeta['TARGETID'])
                tileid_list.append(' '.join(np.unique(expmeta['TILEID'][I]).astype(str)))
                #meta['TILEID_LIST'][M] = ' '.join(np.unique(expmeta['TILEID'][I]).astype(str))
                # store just the zeroth tile for gather_targetphot, below
                if self.coadd_type == 'healpix' or self.coadd_type == 'custom':
                    alltiles.append(expmeta['TILEID'][I][0])
                else:
                    alltiles.append(tileid)

            if self.coadd_type == 'healpix' or self.coadd_type == 'custom':
                meta['TILEID_LIST'] = tileid_list

            # Gather additional info about this pixel.
            if self.coadd_type == 'healpix' or self.coadd_type == 'custom':
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

                # Get the correct fiber number.
                if 'FIBER' in expmeta.colnames:
                    meta['FIBER'] = np.zeros(len(meta), dtype=expmeta['FIBER'].dtype)
                    _, uindx = np.unique(expmeta['TARGETID'], return_index=True)
                    I = geomask.match_to(expmeta[uindx]['TARGETID'], meta['TARGETID'])
                    assert(np.all(expmeta[uindx][I]['TARGETID'] == meta['TARGETID']))
                    meta['FIBER'] = expmeta[uindx[I]]['FIBER']

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
        #t1 = time.time()
        metas = self._gather_photometry(specprod=specprod, alltiles=alltiles)
        self.meta = metas # update
        #log.info(f'Gathered photometric metadata in {time.time()-t1:.2f} seconds.')
        if len(redrockfiles) > 1:
            log.debug(f'Gathered spectrophotometric metadata for {len(redrockfiles)} unique ' + \
                      f'redrockfiles in {time.time()-t0:.2f} seconds.')
        else:
            log.debug(f'Gathered spectrophotometric metadata for {len(redrockfiles)} unique ' + \
                      f'redrockfile in {time.time()-t0:.2f} seconds.')


    @staticmethod
    def update_qso_redshifts(zb, meta, qnfile, mgiifile, fitindx, specprod):
        """Update QSO redshifts using the afterburners.

        """
        from desitarget.targets import main_cmx_or_sv

        if specprod in ['fuji', 'guadalupe', 'himalayas', 'iron']:
            QNthresh = 0.95
        else:
            # updated for Jura, Kibo, Loa, ...
            QNthresh = 0.99

        surv_target, surv_mask, surv = main_cmx_or_sv(meta, scnd=True)
        if surv == 'cmx':
            desi_target = surv_target[0]
            scnd_target = surv_target[-1]
            desi_mask = surv_mask[0]
            scnd_mask = surv_mask[-1]
            IQSO = ((meta[desi_target] & desi_mask['SV0_QSO'] != 0) |
                    (meta[desi_target] & desi_mask['MINI_SV_QSO'] != 0))
            IWISE_VAR_QSO = np.zeros(len(fitindx), bool)
        else:
            desi_target, bgs_target, mws_target, scnd_target = surv_target
            desi_mask, bgs_mask, mws_mask, scnd_mask = surv_mask
            IQSO = meta[desi_target] & desi_mask['QSO'] != 0
            if 'WISE_VAR_QSO' in scnd_mask.names():
                IWISE_VAR_QSO = meta[scnd_target] & scnd_mask['WISE_VAR_QSO'] != 0
            else:
                IWISE_VAR_QSO = np.zeros(len(meta), bool)
        if np.sum(IQSO) > 0 or np.sum(IWISE_VAR_QSO) > 0:
            qn = Table(fitsio.read(qnfile, 'QN_RR', rows=fitindx, columns=QNCOLS))
            assert(np.all(qn['TARGETID'] == meta['TARGETID']))
            log.debug('Updating QSO redshifts using a QN threshold of 0.99.')
            qn['IS_QSO_QN_099'] = np.max(np.array([qn[name] for name in QNLINES]), axis=0) > QNthresh
            iqso = IQSO * qn['IS_QSO_QN_NEW_RR'] * qn['IS_QSO_QN_099']
            if np.sum(iqso) > 0:
                zb['Z'][iqso] = qn['Z_NEW'][iqso]
                zb['ZWARN'][iqso] = qn['ZWARN_NEW'][iqso]
            if np.sum(IWISE_VAR_QSO) > 0:
                mgii = Table(fitsio.read(mgiifile, 'MGII', rows=fitindx, columns=MGIICOLS))
                assert(np.all(mgii['TARGETID'] == meta['TARGETID']))
                iwise_var_qso = (((zb['SPECTYPE'] == 'QSO') | mgii['IS_QSO_MGII'] | qn['IS_QSO_QN_099']) & (IWISE_VAR_QSO & qn['IS_QSO_QN_NEW_RR']))
                if np.sum(iwise_var_qso) > 0:
                    zb['Z'][iwise_var_qso] = qn['Z_NEW'][iwise_var_qso]
                    #zb['Z_ERR'][iwise_var_qso] = qn['ZERR_NEW'][iwise_var_qso]
                    zb['ZWARN'][iwise_var_qso] = qn['ZWARN_NEW'][iwise_var_qso]
                del mgii
            del qn


    def read(self, photometry, fastphot=False, constrain_age=False):
        """Read selected spectra and/or broadband photometry.

        Parameters
        ----------
        fastphot : bool
            Read the broadband photometry; otherwise, handle the DESI
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
            redshift : numpy.float64
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
                Three-element list of :class:`fastspecfit.resolution.Resolution`
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
        from fastspecfit.resolution import Resolution
        from fastspecfit.util import mwdust_transmission

        t0 = time.time()

        SFD = SFDMap(scaling=1.0, mapdir=self.mapdir)

        uniqueid_col = self.phot.uniqueid_col

        alldata, allmeta = [], []
        for ispecfile, (specfile, meta) in enumerate(zip(self.specfiles, self.meta)):
            nobj = len(meta)
            if nobj == 1:
                log.info(f'Reading {nobj} spectrum from {specfile}')
            else:
                log.info(f'Reading {nobj} spectra from {specfile}')

            # Pre-compute the luminosity distance, distance modulus, and age of
            # the universe.
            redshift = meta['Z'].value
            neg = (redshift <= 0.)
            if np.any(neg):
                errmsg = f'{np.sum(neg)}/{len(redshift)} input redshifts are zero or negative; setting to 1e-8!'
                log.warning(errmsg)
                redshift[neg] = 1e-8

            dlum = self.cosmo.luminosity_distance(redshift)
            dmod = self.cosmo.distance_modulus(redshift)
            if constrain_age:
                tuniv = self.cosmo.universe_age(redshift)
            else:
                tuniv = np.full_like(redshift, 100.)

            # Populate 'meta' with dust and filter-related quantities.
            ebv = SFD.ebv(meta['RA'], meta['DEC'])
            meta['EBV'] = ebv

            if 'PHOTSYS' in meta.colnames:
                photsys = meta['PHOTSYS'].value
            else:
                photsys = [''] * nobj

            if hasattr(photometry, 'fiber_filters'):
                mw_transmission_fiberflux = np.ones((nobj, len(photometry.fiber_bands)))

            for onephotsys in set(photsys):
                I = np.where(onephotsys == photsys)[0]

                filters = photometry.filters[onephotsys]
                for band, filt in zip(photometry.bands, filters.names):
                    meta[f'MW_TRANSMISSION_{band.upper()}'][I] = mwdust_transmission(ebv[I], filt)

                if hasattr(photometry, 'fiber_filters'):
                    for iband, filt in enumerate(photometry.fiber_filters[onephotsys].names):
                        mw_transmission_fiberflux[I, iband] = mwdust_transmission(ebv[I], filt)
                else:
                    mw_transmission_fiberflux = None

            if fastphot:
                for iobj in range(nobj):
                    specdata = {
                        'uniqueid': meta[uniqueid_col][iobj],
                        'redshift': redshift[iobj],
                        'photsys': photsys[iobj],
                        'dluminosity': dlum[iobj],
                        'dmodulus': dmod[iobj],
                        'tuniv': tuniv[iobj],
                        }
                    if mw_transmission_fiberflux is not None:
                        specdata.update({'mw_transmission_fiberflux': mw_transmission_fiberflux[iobj, :]})
                    alldata.append(specdata)
            else:
                # Don't use .select since meta and spec can be sorted
                # differently if a non-sorted targetids was passed. Do the
                # selection and sort ourselves.
                os.environ['DESI_LOGLEVEL'] = 'warning'
                spec = read_spectra(specfile)#.select(targets=meta[uniqueid])

                srt = geomask.match_to(spec.fibermap[uniqueid_col], meta['TARGETID'])
                spec = spec[srt]
                assert(np.all(spec.fibermap[uniqueid_col] == meta[uniqueid_col]))

                # Coadd across cameras.
                t0 = time.time()
                coadd_spec = coadd_cameras(spec)
                os.environ['DESI_LOGLEVEL'] = 'info'
                log.debug(f'Coadding across cameras took {time.time()-t0:.2f} seconds.')

                # unpack the desispec.spectra.Spectra objects into simple arrays
                cams = spec.bands
                coadd_cams = coadd_spec.bands[0]
                for iobj in range(nobj):
                    specdata = {
                        'uniqueid': meta[uniqueid_col][iobj],
                        'redshift': redshift[iobj],
                        'photsys': photsys[iobj],
                        'cameras': cams,
                        'dluminosity': dlum[iobj],
                        'dmodulus': dmod[iobj],
                        'tuniv': tuniv[iobj],
                        'ebv': ebv[iobj],
                        'wave0': [spec.wave[cam] for cam in cams],
                        'flux0': [spec.flux[cam][iobj, :] for cam in cams],
                        'ivar0': [spec.ivar[cam][iobj, :] for cam in cams],
                        # Also track the mask---see https://github.com/desihub/desispec/issues/1389
                        'mask0': [spec.mask[cam][iobj, :] for cam in cams],
                        'res0': [spec.resolution_data[cam][iobj, :, :] for cam in cams],
                        'coadd_wave': coadd_spec.wave[coadd_cams],
                        'coadd_flux': coadd_spec.flux[coadd_cams][iobj, :],
                        'coadd_ivar': coadd_spec.ivar[coadd_cams][iobj, :],
                        'coadd_res': [Resolution(coadd_spec.resolution_data[coadd_cams][iobj, :])],
                    }
                    if mw_transmission_fiberflux is not None:
                        specdata.update({'mw_transmission_fiberflux': mw_transmission_fiberflux[iobj, :]})

                    alldata.append(specdata)

            allmeta.append(meta)

        allmeta = vstack(allmeta)

        return alldata, allmeta


    def read_stacked(self, stackfiles, firsttarget=0, ntargets=None,
                     stackids=None, synthphot=True, constrain_age=False):
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
            redshift : numpy.float64
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
                Three-element list of :class:`fastspecfit.resolution.Resolution`
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
        from fastspecfit.resolution import Resolution

        if stackfiles is None:
            errmsg = 'At least one stackfiles file is required.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        if len(stackfiles) == 0:
            errmsg = 'No stackfiles found!'
            log.warning(errmsg)
            raise ValueError(errmsg)

        t0 = time.time()

        stackfiles = sorted(set(stackfiles))
        log.debug(f'Reading and parsing {len(stackfiles)} unique stackfile(s).')

        self.specprod = 'stacked'
        self.coadd_type = 'stacked'
        survey = 'stacked'
        program = 'stacked'
        healpix = np.int32(0)

        READCOLS = ('STACKID', 'Z')
        uniqueid_col = self.phot.uniqueid_col

        alldata, allmeta = [], []
        for stackfile in stackfiles:
            if not os.path.isfile(stackfile):
                log.warning(f'File {stackfile} not found!')
                continue

            # Gather some coadd information from the header.
            log.info('specprod={}, coadd_type={}, survey={}, program={}, healpix={}'.format(
                self.specprod, self.coadd_type, survey, program, healpix))

            # If stackids is *not* given, read everything.
            if stackids is None:
                fitindx = np.arange(fitsio.FITS(stackfile)['STACKINFO'].get_nrows())
            else:
                # We already know we like the input stackids, so no selection
                # needed.
                allstackids = fitsio.read(stackfile, 'STACKINFO', columns='STACKID')
                fitindx = np.where([tid in stackids for tid in allstackids])[0]

            if len(fitindx) == 0:
                log.info(f'No requested targets found in stackfile {stackfile}')
                continue

            # Do we want just a subset of the available objects?
            if ntargets is None:
                _ntargets = len(fitindx)
            else:
                _ntargets = ntargets
            if _ntargets > len(fitindx):
                log.warning('Number of requested ntargets exceeds the number ' + \
                            f'of targets on {stackfile}; reading all of them.')

            __ntargets = len(fitindx)
            fitindx = fitindx[firsttarget:firsttarget+_ntargets]
            if len(fitindx) == 0:
                log.info(f'All {__ntargets} targets in stackfile {stackfile} have been ' + \
                         f'dropped (firsttarget={firsttarget}, ntargets={_ntargets}).')
                continue

            # If firsttarget is a large index then the set can become empty.
            meta = Table(fitsio.read(stackfile, 'STACKINFO', rows=fitindx, columns=READCOLS))
            nobj = len(meta)
            if nobj == 0:
                log.warning('No targets read!')
                return [], []

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

            # Now read the data as in self.read (for unstacked spectra).
            if nobj == 1:
                log.info(f'Reading 1 spectrum from {stackfile}')
            else:
                log.info(f'Reading {nobj} spectra from {stackfile}')

            # Age of the universe.
            redshift = meta['Z'].value
            neg = (redshift <= 0.)
            if np.any(neg):
                errmsg = f'{np.sum(neg)}/{len(redshift)} input redshifts are zero or negative; setting to 1e-8!'
                log.warning(errmsg)
                redshift[neg] = 1e-8

            dlum = self.cosmo.luminosity_distance(redshift)
            dmod = self.cosmo.distance_modulus(redshift)
            if constrain_age:
                tuniv = self.cosmo.universe_age(redshift)
            else:
                tuniv = np.full_like(redshift, 100.)

            # read the data
            wave = fitsio.read(stackfile, 'WAVE')
            npix = len(wave)

            flux = fitsio.read(stackfile, 'FLUX')
            flux = flux[fitindx, :].astype('f8')

            ivar = fitsio.read(stackfile, 'IVAR')
            ivar = ivar[fitindx, :].astype('f8')

            # Check if the file contains a resolution matrix, if it does not
            # then use an identity matrix.
            if 'RES' in fitsio.FITS(stackfile):
                res = fitsio.read(stackfile, 'RES')
                res = res[fitindx, :, :]
            else:
                log.warning('No resolution matrix found; using identity matrix.')
                res = np.ones((nobj, 1, npix)) # Hack!

            # Ppack the data into a simple dictionary.
            for iobj in range(nobj):
                specdata = {
                    'uniqueid': meta[uniqueid_col][iobj],
                    'redshift': redshift[iobj],
                    'photsys': meta['PHOTSYS'][iobj],
                    'cameras': ['brz'],
                    'dluminosity': dlum[iobj],
                    'dmodulus': dmod[iobj],
                    'tuniv': tuniv[iobj],
                    'wave0': [wave],
                    'flux0': [flux[iobj, :]],
                    'ivar0': [ivar[iobj, :]],
                    'mask0': [np.zeros(npix, bool)],
                    'res0': [Resolution(res[iobj, :, :])],
                }

                specdata.update({
                    'coadd_wave': specdata['wave0'][0],
                    'coadd_flux': specdata['flux0'][0],
                    'coadd_ivar': specdata['ivar0'][0],
                    'coadd_res': specdata['res0'],
                    })
                alldata.append(specdata)

            allmeta.append(meta)

        allmeta = vstack(allmeta)

        return alldata, allmeta


    def _gather_photometry(self, specprod=None, alltiles=None):
        """Gather photometry. Unfortunately some of the bandpass information here will
        be repeated (and has to be consistent with) continuum.Fiters.

        """
        from astropy.table import vstack
        from desitarget import geomask
        from fastspecfit.photometry import gather_tractorphot

        input_meta = vstack(self.meta).copy()

        uniqueid_col = self.phot.uniqueid_col
        PHOTCOLS = np.unique(np.hstack((self.phot.readcols, self.phot.fluxcols, self.phot.fluxivarcols)))

        # DR9 or DR10
        if hasattr(self.phot, 'legacysurveydr'):
            from desitarget.io import releasedict

            legacysurveydr = self.phot.legacysurveydr

            # targeting and Tractor columns to read from disk
            tractor = gather_tractorphot(input_meta, columns=PHOTCOLS, legacysurveydir=self.fphotodir)

            # DR9-specific stuff
            if legacysurveydr.lower() == 'dr9' or legacysurveydr.lower() == 'dr10':
                metas = []
                for meta in self.meta:
                    srt = geomask.match_to(tractor[uniqueid_col], meta[uniqueid_col])
                    assert(np.all(meta[uniqueid_col] == tractor[uniqueid_col][srt]))

                    # The fibermaps in fuji and guadalupe (plus earlier productions) had a
                    # variety of errors. Fix those here using
                    # desispec.io.photo.gather_targetphot.
                    if specprod == 'fuji' or specprod == 'guadalupe':
                        from desispec.io.photo import gather_targetphot

                        for env in ['DESI_ROOT', 'DESI_TARGET', 'DESI_SURVEYOPS', 'FIBER_ASSIGN_DIR']:
                            if not env in os.environ:
                                errmsg = f'For fuji and guadalupe productions, missing mandatory environment variable {env}'
                                log.critical(errmsg)
                                raise KeyError(errmsg)

                        input_meta = meta[uniqueid_col, 'TARGET_RA', 'TARGET_DEC']
                        input_meta['TILEID'] = alltiles
                        targets = gather_targetphot(input_meta)
                        assert(np.all(input_meta[uniqueid_col] == targets[uniqueid_col]))
                        for col in meta.colnames:
                            if col in targets.colnames:
                                diffcol = (meta[col] != targets[col])
                                if np.any(diffcol):
                                    log.warning('Updating column {} in metadata table: {}-->{}.'.format(
                                        col, meta[col][0], targets[col][0]))
                                    meta[col][diffcol] = targets[col][diffcol]
                    srt = geomask.match_to(tractor[uniqueid_col], meta[uniqueid_col])
                    assert(np.all(meta[uniqueid_col] == tractor[uniqueid_col][srt]))

                    # Add the tractor catalog quantities (overwriting columns if necessary).
                    for col in tractor.colnames:
                        meta[col] = tractor[col][srt]

                    # special case for some secondary and ToOs
                    I = ((meta['RA'] == 0) & (meta['DEC'] == 0) &
                         (meta['TARGET_RA'] != 0) & (meta['TARGET_DEC'] != 0))
                    meta['RA'][I]  = meta['TARGET_RA'][I]
                    meta['DEC'][I] = meta['TARGET_DEC'][I]
                    assert(np.all((meta['RA'] != 0) * (meta['DEC'] != 0)))

                    # try to repair PHOTSYS
                    # https://github.com/desihub/fastspecfit/issues/75
                    I = ((meta['PHOTSYS'] != 'N') & (meta['PHOTSYS'] != 'S') & (meta['RELEASE'] >= 9000))
                    meta['PHOTSYS'][I] = [releasedict[release] if release >= 9000 else '' for release in meta['RELEASE'][I]]
                    I = ((meta['PHOTSYS'] != 'N') & (meta['PHOTSYS'] != 'S'))
                    if np.any(I):
                        meta['PHOTSYS'][I] = self.resolve(meta[I])
                    I = ((meta['PHOTSYS'] != 'N') & (meta['PHOTSYS'] != 'S'))
                    if np.any(I):
                        errmsg = 'Unsupported value of PHOTSYS.'
                        log.critical(errmsg)
                        raise ValueError(errmsg)

                    # placeholders (to be added in DESISpectra.read)
                    meta['EBV'] = np.zeros(shape=(1,), dtype='f4')
                    for band in self.phot.bands:
                        meta[f'MW_TRANSMISSION_{band.upper()}'] = np.ones(shape=(1,), dtype='f4')

                    metas.append(meta)
        else:
            phot_tbl = Table(fitsio.read(self.fphotodir, ext=self.fphotoext, columns=PHOTCOLS))
            log.info(f'Read {len(phot_tbl):,d} objects from {self.fphotodir}')

            metas = []
            for meta in self.meta:
                srt = geomask.match_to(phot_tbl[uniqueid_col], meta[uniqueid_col])
                assert(np.all(meta[uniqueid_col] == phot_tbl[uniqueid_col][srt]))
                if hasattr(self.phot, 'dropcols'):
                    meta.remove_columns(self.phot.dropcols)
                for col in phot_tbl.colnames:
                    meta[col] = phot_tbl[col][srt]
                # placeholders (to be added in DESISpectra.read)
                meta['EBV'] = np.zeros(shape=(1,), dtype='f4')
                for band in self.phot.bands:
                    meta[f'MW_TRANSMISSION_{band.upper()}'] = np.ones(shape=(1,), dtype='f4')
                metas.append(meta)

        return metas


def read_fastspecfit(fastfitfile, rows=None, metadata_columns=None, specphot_columns=None,
                     fastspec_columns=None, read_models=False):
    """Read the fitting results.

    """
    if os.path.isfile(fastfitfile):
        F = fitsio.FITS(fastfitfile)
        meta = Table(F['METADATA'].read(rows=rows, columns=metadata_columns))
        specphot = Table(F['SPECPHOT'].read(rows=rows, columns=specphot_columns))

        if 'FASTSPEC' in F:
            fastphot = False
            fastfit = Table(F['FASTSPEC'].read(rows=rows, columns=fastspec_columns))
            if read_models:
                models = F['MODELS'].read()
                if rows is not None:
                    models = models[rows, :, :]
            else:
                models = None
        else:
            fastphot = True
            fastfit = None
            models = None

        log.info(f'Read {len(specphot):,d} object(s) from {fastfitfile}')

        # Add specprod to the metadata table so that we can stack across
        # productions (e.g., Fuji+Guadalupe).
        hdr = F[0].read_header()

        if 'SPECPROD' in hdr:
            specprod = hdr['SPECPROD']
            meta['SPECPROD'] = specprod

        if 'COADDTYP' in hdr:
            coadd_type = hdr['COADDTYP']
        else:
            coadd_type = None

        if read_models:
            return meta, specphot, fastfit, coadd_type, fastphot, models
        else:
            return meta, specphot, fastfit, coadd_type, fastphot

    else:
        log.warning(f'File {fastfitfile} not found.')
        if read_models:
            return [None]*6
        else:
            return [None]*5


def write_fastspecfit(meta, specphot, fastfit, modelspectra=None, outfile=None,
                      specprod=None, coadd_type=None, fphotofile=None,
                      template_file=None, emlinesfile=None, fastphot=False,
                      inputz=False, inputseeds=None, nmonte=10, seed=1,
                      uncertainty_floor=0.01, minsnr_balmer_broad=2.5,
                      no_smooth_continuum=False, ignore_photometry=False,
                      broadlinefit=True, use_quasarnet=True, constrain_age=False,
                      verbose=True):
    """Write out.

    """
    import gzip, shutil
    from astropy.io import fits
    from desispec.io.util import fitsheader
    from desiutil.depend import add_dependencies, possible_dependencies, setdep

    t0 = time.time()
    outdir = os.path.dirname(os.path.abspath(os.path.expanduser(os.path.expandvars(outfile))))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    nobj = len(meta)
    if nobj == 1:
        log.info(f'Writing 1 object to {outfile}')
    else:
        log.info(f'Writing {nobj:,d} objects to {outfile}')

    if outfile.endswith('.gz'):
        tmpfile = outfile[:-3]+'.tmp'
    else:
        tmpfile = outfile+'.tmp'

    # Also update mpi._domerge
    primhdr = []
    if specprod:
        primhdr.append(('EXTNAME', 'PRIMARY'))
        primhdr.append(('SPECPROD', (specprod, 'spectroscopic production name')))
    if coadd_type is not None:
        primhdr.append(('COADDTYP', (coadd_type, 'spectral coadd type')))
    primhdr.append(('INPUTZ', (inputz is True, 'input redshifts provided')))
    primhdr.append(('INPUTS', (inputseeds is True, 'input seeds provided')))
    primhdr.append(('CONSAGE', (constrain_age is True, 'constrain SPS ages')))
    primhdr.append(('USEQNET', (use_quasarnet is True, 'use QuasarNet redshifts')))
    primhdr.append(('NMONTE', (nmonte, 'number of Monte Carlo realizations')))
    primhdr.append(('SEED', (seed, 'random seed for Monte Carlo reproducibility')))
    if not fastphot:
        primhdr.append(('NOSCORR', (no_smooth_continuum is True, 'no smooth continuum correction')))
        primhdr.append(('NOPHOTO', (ignore_photometry is True, 'no fitting to photometry')))
        primhdr.append(('BRDLFIT', (broadlinefit is True, 'carry out broad-line fitting')))
        primhdr.append(('UFLOOR', (uncertainty_floor, 'fractional uncertainty floor')))
        primhdr.append(('SNRBBALM', (minsnr_balmer_broad, 'minimum broad Balmer S/N')))

    primhdr = fitsheader(primhdr)
    add_dependencies(primhdr, module_names=possible_dependencies+['fastspecfit'],
                     envvar_names=('DESI_SPECTRO_REDUX', 'DUST_DIR', 'FTEMPLATES_DIR', 'FPHOTO_DIR'))
    if fphotofile:
        setdep(primhdr, 'FPHOTO_FILE', str(fphotofile))
    if template_file:
        setdep(primhdr, 'FTEMPLATES_FILE', os.path.basename(template_file))
    if emlinesfile:
        setdep(primhdr, 'EMLINES_FILE', str(emlinesfile))

    meta.meta['EXTNAME'] = 'METADATA'
    specphot.meta['EXTNAME'] = 'SPECPHOT'

    hdus = fits.HDUList()
    hdus.append(fits.PrimaryHDU(None, primhdr))
    hdus.append(fits.convenience.table_to_hdu(meta))
    hdus.append(fits.convenience.table_to_hdu(specphot))
    if fastfit is not None:
        fastfit.meta['EXTNAME'] = 'FASTSPEC'
        hdus.append(fits.convenience.table_to_hdu(fastfit))

    if modelspectra is not None:
        hdu = fits.ImageHDU(name='MODELS')
        # [nobj, 3, nwave]
        hdu.data = np.swapaxes(np.array([modelspectra['CONTINUUM'].data,
                                         modelspectra['SMOOTHCONTINUUM'].data,
                                         modelspectra['EMLINEMODEL'].data]), 0, 1)
        for key in modelspectra.meta:
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
        log.debug(f'Writing out took {time.time()-t0:.2f} seconds.')


def get_qa_filename(metadata, coadd_type, outprefix=None, outdir=None,
                    fastphot=False):
    """Build the QA filename.

    """
    import astropy.table.row as row

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

    if type(metadata) is row.Row or type(metadata) is np.void:
        pngfile = _one_filename(metadata)
    else:
        pngfile = [_one_filename(_metadata) for _metadata in metadata]

    return pngfile


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

    specdir = os.path.join(os.environ.get('DESI_SPECTRO_REDUX'), specprod, 'healpix',
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
    template_version = {redhdr[f'TEMNAM{nn:02d}']: redhdr[f'TEMVER{nn:02d}'] for nn in range(10)}

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


def select(metadata, specphot, fastfit=None, coadd_type='healpix',
           healpixels=None, tiles=None, nights=None, return_index=False):
    """Optionally trim to a particular healpix or tile and/or night."""
    nobj = len(metadata)
    if coadd_type == 'healpix':
        if healpixels is not None:
            strpixels = ','.join(healpixels)
            keep = np.isin(metadata['HEALPIX'].astype(str), healpixels)
            log.info(f'Keeping {np.sum(keep):,d}/{nobj:,d} objects from healpixels(s) {strpixels}')
        else:
            keep = np.ones(nobj, bool)
    else:
        if tiles is not None and nights is not None:
            strtiles = ','.join(tiles)
            strnights = ','.join(nights)
            keep = np.isin(metadata['TILEID'].astype(str), tiles) * np.isin(metadata['NIGHT'].astype(str), nights)
            log.info(f'Keeping {np.sum(keep):,d}/{nobj:,d} objects from tile(s) {strtiles} and night(s) {strnights}')
        elif tiles is not None and nights is None:
            strtiles = ','.join(tiles)
            keep = np.isin(metadata['TILEID'].astype(str), tiles)
            log.info(f'Keeping {np.sum(keep):,d}/{nobj:,d} objects from tile(s) {strtiles}')
        elif tiles is None and nights is not None:
            strnights = ','.join(nights)
            keep = np.isin(metadata['NIGHT'].astype(str), nights)
            log.info(f'Keeping {np.sum(keep):,d}/{nobj:,d} objects from night(s) {strnights}')
        else:
            keep = np.ones(nobj, bool) # keep everything

    if return_index:
        return np.where(keep)[0]
    else:
        if fastfit is not None:
            fastfit = fastfit[keep]
        return metadata[keep], specphot[keep], fastfit


def get_output_dtype(specprod, phot, linetable, ncoeff, cameras=['B', 'R', 'Z'],
                     specphot=False, fastphot=False, fitstack=False):
    """
    Get type of one fastspecfit output data record, along
    with dictionary of units for any fields that have them.

    """
    import astropy.units as u

    out_dtype = []
    out_units = {}

    def add_field(name, dtype, shape=None, unit=None):
        if shape is not None:
            t = (name, dtype, shape) # subarray
        else:
            t = (name, dtype)
        out_dtype.append(t)

        if unit is not None:
            out_units[name] = unit

    if specphot:
        #add_field('Z', dtype='f8') # redshift
        add_field('SEED', dtype=np.int64)
        add_field('COEFF', shape=(ncoeff,), dtype='f4')

        if not fastphot:
            add_field('RCHI2', dtype='f4')      # full-spectrum reduced chi2
            add_field('RCHI2_LINE', dtype='f4') # reduced chi2 with broad line-emission
            add_field('RCHI2_CONT', dtype='f4') # rchi2 fitting just to the continuum (spec+phot)
        add_field('RCHI2_PHOT', dtype='f4') # rchi2 fitting just to the photometry

        add_field('VDISP', dtype='f4', unit=u.kilometer/u.second)
        if not fastphot:
            add_field('VDISP_IVAR', dtype='f4', unit=u.second**2/u.kilometer**2)
        add_field('TAUV', dtype='f4')
        add_field('TAUV_IVAR', dtype='f4')
        add_field('AGE', dtype='f4', unit=u.Gyr)
        add_field('AGE_IVAR', dtype='f4', unit=1/u.Gyr**2)
        add_field('ZZSUN', dtype='f4')
        add_field('ZZSUN_IVAR', dtype='f4')
        add_field('LOGMSTAR', dtype='f4', unit=u.solMass)
        add_field('LOGMSTAR_IVAR', dtype='f4', unit=1/u.solMass**2)
        add_field('SFR', dtype='f4', unit=u.solMass/u.year)
        add_field('SFR_IVAR', dtype='f4', unit=u.year**2/u.solMass**2)
        #add_field('FAGN', dtype='f4')

        if not fastphot:
            add_field('DN4000', dtype='f4')
            add_field('DN4000_OBS', dtype='f4')
            add_field('DN4000_IVAR', dtype='f4')
        add_field('DN4000_MODEL', dtype='f4')
        add_field('DN4000_MODEL_IVAR', dtype='f4')

        if not fastphot:
            # observed-frame photometry synthesized from the spectra
            for band in phot.synth_bands:
                add_field(f'FLUX_SYNTH_{band.upper()}', dtype='f4', unit='nanomaggies')
                #add_field(f'FLUX_SYNTH_IVAR_{band.upper()}'), dtype='f4', unit='1/nanomaggies**2')
            # observed-frame photometry synthesized the best-fitting spectroscopic model
            for band in phot.synth_bands:
                add_field(f'FLUX_SYNTH_SPECMODEL_{band.upper()}', dtype='f4', unit='nanomaggies')
        # observed-frame photometry synthesized the best-fitting continuum model
        for band in phot.bands:
            add_field(f'FLUX_SYNTH_PHOTMODEL_{band.upper()}', dtype='f4', unit='nanomaggies')

        for band, shift in zip(phot.absmag_bands, phot.band_shift):
            band = band.upper()
            shift = int(10*shift)
            add_field(f'ABSMAG{shift:02d}_{band}', dtype='f4', unit=u.mag) # absolute magnitudes
            add_field(f'ABSMAG{shift:02d}_IVAR_{band}', dtype='f4', unit=1/u.mag**2)
            add_field(f'ABSMAG{shift:02d}_SYNTH_{band}', dtype='f4', unit=u.mag) # absolute magnitudes
            add_field(f'ABSMAG{shift:02d}_SYNTH_IVAR_{band}', dtype='f4', unit=1/u.mag**2)
        for band, shift in zip(phot.absmag_bands, phot.band_shift):
            band = band.upper()
            shift = int(10*shift)
            add_field(f'KCORR{shift:02d}_{band}', dtype='f4', unit=u.mag)

        for wave in ['1500', '2800']:
            add_field(f'LOGLNU_{wave}',  dtype='f4', unit=10**(-28)*u.erg/u.second/u.Hz)
            add_field(f'LOGLNU_{wave}_IVAR',  dtype='f4', unit=10**(-28)*u.erg/u.second/u.Hz)
        for wave in ['1450', '1700', '3000', '5100']:
            add_field(f'LOGL_{wave}', dtype='f4', unit=10**(10)*u.solLum)
            add_field(f'LOGL_{wave}_IVAR', dtype='f4', unit=10**(10)*u.solLum)
        for line in ['FLYA_1215', 'FOII_3727', 'FHBETA', 'FOIII_5007', 'FHALPHA']:
            add_field(f'{line}_CONT', dtype='f4', unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
            add_field(f'{line}_CONT_IVAR', dtype='f4', unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))

    else:

        if not fastphot:
            for cam in cameras:
                add_field(f'SNR_{cam.upper()}', dtype='f4') # median S/N in each camera
            for cam in cameras:
                add_field(f'SMOOTHCORR_{cam.upper()}', dtype='f4')

            # aperture corrections
            add_field('APERCORR', dtype='f4') # median aperture correction
            for band in phot.synth_bands:
                add_field(f'APERCORR_{band.upper()}', dtype='f4')

        if not fastphot:
            add_field('INIT_SIGMA_UV', dtype='f4', unit=u.kilometer / u.second)
            add_field('INIT_SIGMA_NARROW', dtype='f4', unit=u.kilometer / u.second)
            add_field('INIT_SIGMA_BALMER', dtype='f4', unit=u.kilometer / u.second)
            add_field('INIT_VSHIFT_UV', dtype='f4', unit=u.kilometer / u.second)
            add_field('INIT_VSHIFT_NARROW', dtype='f4', unit=u.kilometer / u.second)
            add_field('INIT_VSHIFT_BALMER', dtype='f4', unit=u.kilometer / u.second)
            add_field('INIT_BALMER_BROAD', dtype=bool)

        if not fastphot:
            # Add chi2 metrics
            #add_field('DOF',  dtype='i8') # full-spectrum dof
            #add_field('NDOF_LINE', dtype='i8') # number of degrees of freedom corresponding to rchi2_line
            #add_field('DOF_BROAD', dtype='i8')
            add_field('DELTA_LINECHI2', dtype='f4') # delta-reduced chi2 with and without broad line-emission
            add_field('DELTA_LINENDOF', dtype=np.int32)

            # special columns for the fitted doublets
            add_field('MGII_DOUBLET_RATIO', dtype='f4')
            add_field('MGII_DOUBLET_RATIO_IVAR', dtype='f4')
            add_field('OII_DOUBLET_RATIO', dtype='f4')
            add_field('OII_DOUBLET_RATIO_IVAR', dtype='f4')
            add_field('OIII_DOUBLET_RATIO', dtype='f4')
            add_field('OIII_DOUBLET_RATIO_IVAR', dtype='f4')
            add_field('NII_DOUBLET_RATIO', dtype='f4')
            add_field('NII_DOUBLET_RATIO_IVAR', dtype='f4')
            add_field('SII_DOUBLET_RATIO', dtype='f4')
            add_field('SII_DOUBLET_RATIO_IVAR', dtype='f4')
            add_field('OIIRED_DOUBLET_RATIO', dtype='f4')
            add_field('OIIRED_DOUBLET_RATIO_IVAR', dtype='f4')

            for line in linetable['name']:
                line = line.upper()
                add_field(f'{line}_MODELAMP', dtype='f4',
                                      unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
                #add_field(f'{line}_MODELAMP_IVAR', dtype='f4',
                #                      unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2)
                add_field(f'{line}_AMP', dtype='f4',
                                      unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
                add_field(f'{line}_AMP_IVAR', dtype='f4',
                                      unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2)
                add_field(f'{line}_FLUX', dtype='f4',
                                      unit=10**(-17)*u.erg/(u.second*u.cm**2))
                #add_field(f'{line}_FLUX_GAUSS_IVAR', dtype='f4',
                #                      unit=10**34*u.second**2*u.cm**4/u.erg**2)
                add_field(f'{line}_FLUX_IVAR', dtype='f4',
                                      unit=10**34*u.second**2*u.cm**4/u.erg**2)
                add_field(f'{line}_BOXFLUX', dtype='f4',
                                      unit=10**(-17)*u.erg/(u.second*u.cm**2))
                add_field(f'{line}_BOXFLUX_IVAR', dtype='f4',
                                      unit=10**34*u.second**2*u.cm**4/u.erg**2)

                add_field(f'{line}_VSHIFT', dtype='f4',
                                      unit=u.kilometer/u.second)
                add_field(f'{line}_VSHIFT_IVAR', dtype='f4',
                                      unit=u.second**2/u.kilometer**2)
                add_field(f'{line}_SIGMA', dtype='f4',
                                      unit=u.kilometer / u.second)
                add_field(f'{line}_SIGMA_IVAR', dtype='f4',
                                      unit=u.second**2/u.kilometer**2)

                add_field(f'{line}_CONT', dtype='f4',
                                      unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom))
                add_field(f'{line}_CONT_IVAR', dtype='f4',
                                      unit=10**34*u.second**2*u.cm**4*u.Angstrom**2/u.erg**2)
                add_field(f'{line}_EW', dtype='f4',
                                      unit=u.Angstrom)
                add_field(f'{line}_EW_IVAR', dtype='f4',
                                      unit=1/u.Angstrom**2)
                add_field(f'{line}_FLUX_LIMIT', dtype='f4',
                                      unit=u.erg/(u.second*u.cm**2))
                #add_field(f'{line}_EW_LIMIT', dtype='f4',
                #                      unit=u.Angstrom)
                add_field(f'{line}_CHI2', dtype='f4')
                add_field(f'{line}_NPIX', dtype=np.int32)

            for line in ['CIV_1549', 'MGII_2800', 'HBETA', 'OIII_5007']:
                for n in range(1, 4):
                    add_field(f'{line}_MOMENT{n}', dtype='f4', unit=u.Angstrom**n)
                    add_field(f'{line}_MOMENT{n}_IVAR', dtype='f4', unit=1/(u.Angstrom**n)**2)

    return np.dtype(out_dtype, align=True), out_units


def create_output_meta(input_meta, phot, fastphot=False, fitstack=False):
    """Create the fastspecfit output metadata table.

    """
    from fastspecfit.io import TARGETINGBITS
    from astropy.table import Table
    import astropy.units as u

    nobj = len(input_meta)

    # The information stored in the metadata table depends on which spectra
    # were fitted (exposures, nightly coadds, deep coadds).
    if fitstack:
        fluxcols = ['PHOTSYS']
    else:
        fluxcols = []
        if hasattr(phot, 'outcols'):
            fluxcols.extend(phot.outcols)
        if hasattr(phot, 'fiber_bands'):
            fluxcols.extend([f'FIBERFLUX_{band.upper()}' for band in phot.fiber_bands])
            fluxcols.extend([f'FIBERTOTFLUX_{band.upper()}' for band in phot.fiber_bands])
        fluxcols.extend(phot.fluxcols)
        fluxcols.extend(phot.fluxivarcols)
        fluxcols.append('EBV')
        fluxcols.extend([f'MW_TRANSMISSION_{band.upper()}' for band in phot.bands])

    colunit = {'RA': u.deg, 'DEC': u.deg, 'EBV': u.mag}
    for fcol, icol in zip(phot.fluxcols, phot.fluxivarcols):
        colunit[fcol.upper()] = phot.photounits
        colunit[icol.upper()] = f'{phot.photounits}-2'
    if hasattr(phot, 'fiber_bands'):
        for band in phot.fiber_bands:
            band = band.upper()
            colunit[f'FIBERFLUX_{band}'] = phot.photounits
            colunit[f'FIBERTOTFLUX_{band}'] = phot.photounits

    skipcols = fluxcols + ['OBJTYPE', 'TARGET_RA', 'TARGET_DEC', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID', 'RELEASE']

    if fitstack:
        redrockcols = ('Z')
    else:
        redrockcols = ('Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE', 'SUBTYPE', 'Z_RR', 'ZWARN_RR',
                       'TSNR2_BGS', 'TSNR2_LRG', 'TSNR2_ELG', 'TSNR2_QSO', 'TSNR2_LYA')

    meta = Table()
    metacols = set(input_meta.colnames)

    # All of this business is so we can get the columns in the order we want
    # (i.e., the order that matches the data model).
    if fitstack:
        for metacol in ('STACKID', 'SURVEY', 'PROGRAM'):
            if metacol in metacols:
                meta[metacol] = input_meta[metacol]
    else:
        for metacol in ('TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'NIGHT', 'FIBER',
                        'EXPID', 'TILEID_LIST', 'RA', 'DEC', 'COADD_FIBERSTATUS'):
            if metacol in metacols:
                meta[metacol] = input_meta[metacol]

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

        for bitcol in TARGETINGCOLS:
            if bitcol in metacols:
                meta[bitcol] = input_meta[bitcol]
            else:
                meta[bitcol] = np.zeros(shape=(1,), dtype=np.int64)

    for redrockcol in redrockcols:
        if redrockcol in metacols: # the Z_RR from quasarnet may not be present
            meta[redrockcol] = input_meta[redrockcol]

    for fluxcol in fluxcols:
        meta[fluxcol] = input_meta[fluxcol]

    # assign units to any columns that should have them
    for col in meta.colnames:
        if col in colunit:
            meta[col].unit = colunit[col]

    return meta


def create_output_table(records, meta, units, fitstack=False):
    """Generate the output FASTSPEC/FASTPHOT or SPECPHOT table.

    """
    from astropy.table import hstack

    # Initialize the output table from the metadata table.
    metacols = set(meta.colnames)

    if fitstack:
        initcols = ('STACKID', 'SURVEY', 'PROGRAM')
    else:
        initcols = ('TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'NIGHT', 'FIBER', 'EXPID')
    initcols = [col for col in initcols if col in metacols]

    cdata = [meta[col] for col in initcols]

    output_table = Table()
    output_table.add_columns(cdata)

    # Now add the measurements. Columns and their dtypes are inferred from the
    # array's dtype.
    output_table = hstack((output_table, Table(np.array(records), units=units)))

    return output_table

