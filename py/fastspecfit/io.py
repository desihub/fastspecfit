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
from astropy.table import Table, vstack, hstack

from desiutil.log import get_logger
log = get_logger()

# Default environment variables.
DESI_ROOT_NERSC = '/global/cfs/cdirs/desi'
DUST_DIR_NERSC = '/global/cfs/cdirs/cosmo/data/dust/v0_1'
DR9_DIR_NERSC = '/global/cfs/cdirs/desi/external/legacysurvey/dr9'
FASTSPECFIT_TEMPLATES_NERSC = '/global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z'

# list of all possible targeting bit columns
TARGETINGBITS = {
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
    'healpix': ['TARGETID', 'TILEID'] # tileid will be an array
    }

# redshift columns to read
REDSHIFTCOLS = ['TARGETID', 'Z', 'ZWARN', 'SPECTYPE', 'DELTACHI2']

# quasarnet afterburner columns to read
QNCOLS = ['TARGETID', 'Z_NEW', 'IS_QSO_QN_NEW_RR', 'C_LYA', 'C_CIV',
          'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']

# targeting and Tractor columns to read from disk
#TARGETCOLS = ['TARGETID', 'RA', 'DEC', 'FLUX_W3', 'FLUX_W4', 'FLUX_IVAR_W3', 'FLUX_IVAR_W4']
TARGETCOLS = ['TARGETID', 'RA', 'DEC',
              'RELEASE', 'LS_ID',
              #'PHOTSYS',
              'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
              'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
              'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4',
              'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
              'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FLUX_IVAR_W3', 'FLUX_IVAR_W4']#,

class DESISpectra(object):
    def __init__(self, redux_dir=None, fiberassign_dir=None, dr9dir=None):
        """Class to read in DESI spectra and associated metadata.

        Parameters
        ----------
        redux_dir : str
            Full path to the location of the reduced DESI data. Optional and
            defaults to `$DESI_SPECTRO_REDUX`.
        
        fiberassign_dir : str
            Full path to the location of the fiberassign files. Optional and
            defaults to `$DESI_ROOT/target/fiberassign/tiles/trunk`.

        """
        desi_root = os.environ.get('DESI_ROOT', DESI_ROOT_NERSC)

        if redux_dir is None:
            self.redux_dir = os.path.join(desi_root, 'spectro', 'redux')
        else:
            self.redux_dir = redux_dir
            
        if fiberassign_dir is None:
            self.fiberassign_dir = os.path.join(desi_root, 'target', 'fiberassign', 'tiles', 'trunk')
        else:
            self.fiberassign_dir = fiberassign_dir

        if dr9dir is None:
            self.dr9dir = DR9_DIR_NERSC
        else:
            self.dr9dir = dr9dir

    def select(self, redrockfiles, zmin=0.001, zmax=None, zwarnmax=None,
               targetids=None, firsttarget=0, ntargets=None,
               use_quasarnet=True, redrockfile_prefix='redrock-',
               specfile_prefix='coadd-', qnfile_prefix='qso_qn-'):
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
        use_quasarnet : `bool`
            Use QuasarNet to improve QSO redshifts, if the afterburner file is
            present. Defaults to `True`.
        redrockfile_prefix : str
            Prefix of the `redrockfiles`. Defaults to `redrock-`.
        specfile_prefix : str
            Prefix of the spectroscopic coadds corresponding to the input
            Redrock file(s). Defaults to `coadd-`.
        qnfile_prefix : str
            Prefix of the QuasarNet afterburner file. Defaults to `coadd-`.

        Attributes
        ----------
        coadd_type : str
            Type of coadded spectra (healpix, cumulative, pernight, or perexp).
        meta : list of :class:`astropy.table.Table`s
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
        from desiutil.depend import getdep
        from desitarget.io import releasedict        
        from desispec.io.photo import gather_tractorphot

        if zmin <= 0.0:
            errmsg = 'zmin must be >= 0.'
            log.critical(errmsg)
            raise ValueError(zmin)
        
        if zmax is None:
            zmax = 99.0

        if zwarnmax is None:
            zwarnmax = 99999
            
        if zmin >= zmax:
            errmsg = 'zmin must be <= zmax.'
            log.critical(errmsg)
            raise ValueError(zmin)
        
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
        log.info('Reading and parsing {} unique redrockfile(s)'.format(len(redrockfiles)))

        alltiles = []
        self.redrockfiles, self.specfiles, self.meta = [], [], []
        
        for ired, redrockfile in enumerate(np.atleast_1d(redrockfiles)):
            if not os.path.isfile(redrockfile):
                log.warning('File {} not found!'.format(redrockfile))
                continue
            
            specfile = redrockfile.replace(redrockfile_prefix, specfile_prefix)
            if not os.path.isfile(specfile):
                log.warning('File {} not found!'.format(specfile))
                continue
            
            # Can we use the quasarnet afterburner file to improve QSO redshifts?
            qnfile = redrockfile.replace(redrockfile_prefix, qnfile_prefix)
            if os.path.isfile(qnfile) and use_quasarnet:
                use_qn = True
            else:
                use_qn = False

            # Gather some coadd information from the header. Note: this code is
            # only compatible with Fuji & Guadalupe headers and later.
            hdr = fitsio.read_header(specfile, ext=0)
            if not 'SPGRP' in hdr:
                errmsg = 'SPGRP header card missing from spectral file {}'.format(specfile)
                log.critical(errmsg)
                raise ValueError(errmsg)

            specprod = getdep(hdr, 'SPECPROD')
            self.specprod = specprod
            if specprod == 'fuji': # EDR
                TARGETINGCOLS = TARGETINGBITS[specprod]
            else:
                TARGETINGCOLS = TARGETINGBITS['default']

            self.coadd_type = hdr['SPGRP']
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

            # add targeting columns
            allfmcols = np.array(fitsio.FITS(specfile)['FIBERMAP'].get_colnames())
            READFMCOLS = FMCOLS + [col for col in TARGETINGCOLS if col in allfmcols]
                    
            # If targetids is *not* given we have to choose "good" objects
            # before subselecting (e.g., we don't want sky spectra).
            if targetids is None:
                zb = fitsio.read(redrockfile, 'REDSHIFTS', columns=REDSHIFTCOLS)
                # Are we reading individual exposures or coadds?
                meta = fitsio.read(specfile, 'FIBERMAP', columns=READFMCOLS)
                assert(np.all(zb['TARGETID'] == meta['TARGETID']))
                # also update mpi.get_ntargets_one
                fitindx = np.where((zb['Z'] > zmin) * (zb['Z'] < zmax) *
                                   (meta['OBJTYPE'] == 'TGT') *
                                   (meta['TARGETID'] > 0) *
                                   (zb['ZWARN'] <= zwarnmax))[0]
            else:
                # We already know we like the input targetids, so no selection
                # needed.
                alltargetids = fitsio.read(redrockfile, 'REDSHIFTS', columns='TARGETID')
                fitindx = np.where([tid in targetids for tid in alltargetids])[0]                
                
            if len(fitindx) == 0:
                log.info('No requested targets found in redrockfile {}'.format(redrockfile))
                continue

            # Do we want just a subset of the available objects?
            if ntargets is None:
                _ntargets = len(fitindx)
            else:
                _ntargets = ntargets
            if _ntargets > len(fitindx):
                log.warning('Number of requested ntargets exceeds the number of targets on {}; reading all of them.'.format(
                    redrockfile))
                #raise ValueError

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
            assert(np.all(zb['TARGETID'] == meta['TARGETID']))

            # Update the redrock redshift when quasarnet disagrees. From Edmond:
            # the QN afterburner is run with a threshold 0.5. With VI, we choose
            # 0.95 as final threshold. Note, the IS_QSO_QN_NEW_RR column
            # contains only QSO for QN which are not QSO for RR.
            if use_qn:
                qn = Table(fitsio.read(qnfile, 'QN_RR', rows=fitindx, columns=QNCOLS))
                assert(np.all(qn['TARGETID'] == meta['TARGETID']))
                log.info('Updating QSO redshifts using a QN threshold of 0.95.')
                qn['IS_QSO_QN'] = np.max(np.array([qn[name] for name in ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']]), axis=0) > 0.95
                qn['IS_QSO_QN_NEW_RR'] &= qn['IS_QSO_QN']
                #zb.add_column(zb['Z'], name='Z_RR', index=2) # add it after 'Z'
                zb['Z_RR'] = zb['Z'] # add it at the end
                if np.any(qn['IS_QSO_QN_NEW_RR']):
                    zb['Z'][qn['IS_QSO_QN_NEW_RR']] = qn['Z_NEW'][qn['IS_QSO_QN_NEW_RR']]
                del qn
            # add an empty Z_RR column?
            #else:
            #    zb['Z_RR'] = zb['Z'] # add it at the end
                
            # astropy 5.0 "feature" -- join no longer preserves order, ugh.
            zb.remove_column('TARGETID')
            meta = hstack((zb, meta))
            #meta = join(zb, meta, keys='TARGETID')
            del zb

            # Get the unique set of tiles contributing to the coadded spectra
            # from EXP_FIBERMAP.
            expmeta = fitsio.read(specfile, 'EXP_FIBERMAP', columns=EXPFMCOLS[self.coadd_type])
            I = np.isin(expmeta['TARGETID'], meta['TARGETID'])
            if np.count_nonzero(I) == 0:
                errmsg = 'No matching targets in exposure table.'
                log.critical(errmsg)
                raise ValueError(errmsg)
            expmeta = Table(expmeta[I])

            #tiles = np.unique(np.atleast_1d(expmeta['TILEID']).data)
            #alltiles.append(tiles)

            # build the list of tiles that went into each unique target / coadd
            tileid_list = [] # variable length, so need to build the array first
            for tid in meta['TARGETID']:
                I = tid == expmeta['TARGETID']
                tileid_list.append(' '.join(np.unique(expmeta['TILEID'][I]).astype(str)))
                #meta['TILEID_LIST'][M] = ' '.join(np.unique(expmeta['TILEID'][I]).astype(str))
                if self.coadd_type == 'healpix':
                    alltiles.append(expmeta['TILEID'][I][0]) # store just the zeroth tile for gather_targetphot, below
                else:
                    alltiles.append(tileid)
            if self.coadd_type == 'healpix':                    
                meta['TILEID_LIST'] = tileid_list

            # Gather additional info about this pixel.
            if self.coadd_type == 'healpix':
                meta['SURVEY'] = survey
                meta['PROGRAM'] = program
                meta['HEALPIX'] = healpix
            else:
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

        if len(self.meta) == 0:
            log.warning('No targets read!')
            return

        # Use the metadata in the fibermap to retrieve the LS-DR9 source
        # photometry.
        t0 = time.time()
        targets = gather_tractorphot(vstack(self.meta), columns=TARGETCOLS, dr9dir=self.dr9dir)

        # bug! https://github.com/desihub/fastspecfit/issues/75
        #from desitarget.io import releasedict
        #for imeta, meta in enumerate(self.meta):
        #    ibug = np.where((meta['RELEASE'] > 0) * (meta['BRICKID'] > 0) * (meta['BRICK_OBJID'] > 0) * (meta['PHOTSYS'] == ''))[0]
        #    if len(ibug) > 0:
        #        meta['PHOTSYS', 'RELEASE', 'BRICKID', 'BRICK_OBJID', 'TARGETID', 'TILEID', 'NIGHT', 'FIBER'][ibug]
        #        from desitarget.targets import decode_targetid
        #        objid, brickid, release, mock, sky, gaia = decode_targetid(meta['TARGETID'][ibug])
        #        meta['PHOTSYS'][ibug] = [releasedict[release] if release >= 9000 else '' for release in meta['RELEASE'][ibug]]
        #        self.meta[imeta] = meta                    
        #targets = gather_tractorphot(vstack(self.meta), columns=TARGETCOLS, dr9dir=self.dr9dir)

        metas = []
        for meta in self.meta:
            srt = np.hstack([np.where(tid == targets['TARGETID'])[0] for tid in meta['TARGETID']])
            assert(np.all(meta['TARGETID'] == targets['TARGETID'][srt]))
            # Prefer the target catalog quantities over those in the fiberassign
            # table, unless the target catalog is zero.
            for col in targets.colnames:
                meta[col] = targets[col][srt]
            # try to repair PHOTSYS
            # https://github.com/desihub/fastspecfit/issues/75
            I = (meta['PHOTSYS'] == '') * (meta['RELEASE'] > 0)
            if np.sum(I) > 0:
                meta['PHOTSYS'][I] = [releasedict[release] if release >= 9000 else '' for release in meta['RELEASE'][I]]
            # special case for some secondary and ToOs
            I = (meta['RA'] == 0) * (meta['DEC'] == 0) * (meta['TARGET_RA'] != 0) * (meta['TARGET_DEC'] != 0)
            if np.sum(I) > 0:
                meta['RA'][I] = meta['TARGET_RA'][I]
                meta['DEC'][I] = meta['TARGET_DEC'][I]
            assert(np.all((meta['RA'] != 0) * (meta['DEC'] != 0)))
            nobj = len(meta)
            # placeholder (to be added in DESISpectra.read_and_unpack)
            for band in ['G', 'R', 'Z', 'W1', 'W2', 'W3', 'W4']:
                meta['MW_TRANSMISSION_{}'.format(band)] = np.ones(shape=(1,), dtype='f4')
            metas.append(meta)
            
        log.info('Gathered photometric metadata for {} objects in {:.2f} sec'.format(len(targets), time.time()-t0))
        self.meta = metas # update

    def read_and_unpack(self, CFit, fastphot=False, synthphot=True,
                        remember_coadd=False):
        """Read and unpack selected spectra or broadband photometry.
        
        Parameters
        ----------
        CFit : :class:`fastspecfit.continuum.ContinuumFit` class
            Continuum-fitting class which contains filter curves and some additional
            photometric convenience functions.
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
                identified in :class:`CFit.build_linemask` to belong to emission
                lines.
            contpix : :class:`list`
                Three-element list of pixel indices, one per camera, which were
                identified in :class:`CFit.build_linemask` to not be
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
        from desispec.resolution import Resolution
        from desispec.coaddition import coadd_cameras
        from desispec.io import read_spectra
        from desiutil.dust import mwdust_transmission, dust_transmission

        # Read everything into a simple dictionary.
        t0 = time.time()
        alldata = []
        for ispec, (specfile, meta) in enumerate(zip(self.specfiles, self.meta)):
            log.info('Reading {} spectra from {}'.format(len(meta), specfile))

            # sometimes this is an astropy.table.Row!
            meta = Table(meta)

            ebv = CFit.SFDMap.ebv(meta['RA'], meta['DEC'])

            if not fastphot:
                spec = read_spectra(specfile).select(targets=meta['TARGETID'])
                ## Sometimes the order is shuffled, which is annoying.
                #if not np.all(spec.fibermap['TARGETID'] == meta['TARGETID']):
                #srt = np.hstack([np.where(tid == targets['TARGETID'])[0] for tid in meta['TARGETID']])
                assert(np.all(spec.fibermap['TARGETID'] == meta['TARGETID']))

                # Coadd across cameras.
                #t1 = time.time()                
                coadd_spec = coadd_cameras(spec)
                coadd_bands = coadd_spec.bands[0]
                #log.info('Coadded the cameras in {:.2f} sec'.format(time.time()-t1))
                
            for igal in np.arange(len(meta)):
                # Unpack the data and correct for Galactic extinction. Also flag pixels that
                # may be affected by emission lines.
                data = {'targetid': meta['TARGETID'][igal], 'zredrock': meta['Z'][igal],
                        'photsys': meta['PHOTSYS'][igal]}#, 'photsys_south': dec < self.desitarget_resolve_dec()}

                if data['photsys'] == 'S':
                    filters = CFit.decam
                    allfilters = CFit.decamwise
                else:
                    filters = CFit.bassmzls
                    allfilters = CFit.bassmzlswise

                # Unpack the imaging photometry and correct for MW dust.

                # all photometry; do not match the Legacy Surveys here because
                # we want the MW dust extinction correction we apply to the
                # spectra to be self-consistent with how we correct the
                # photometry for dust.
                
                #meta['MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']
                #mw_transmission_flux = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(allfilters.effective_wavelengths.value, Rv=CFit.RV))
                if data['photsys'] != '':
                    mw_transmission_flux = np.array([mwdust_transmission(ebv[igal], band, data['photsys'], match_legacy_surveys=False) for band in CFit.bands])
                    for band, mwdust in zip(CFit.bands, mw_transmission_flux):
                        #print(band, mwdust)
                        self.meta[ispec]['MW_TRANSMISSION_{}'.format(band.upper())][igal] = mwdust
                else:
                    mw_transmission_flux = np.array([self.meta[ispec]['MW_TRANSMISSION_{}'.format(band.upper())][igal] for band in CFit.bands])

                maggies = np.zeros(len(CFit.bands))
                ivarmaggies = np.zeros(len(CFit.bands))
                for iband, band in enumerate(CFit.bands):
                    maggies[iband] = meta['FLUX_{}'.format(band.upper())][igal] / mw_transmission_flux[iband]
                    ivarmaggies[iband] = meta['FLUX_IVAR_{}'.format(band.upper())][igal] * mw_transmission_flux[iband]**2
                    
                if not np.all(ivarmaggies >= 0):
                    errmsg = 'Some ivarmaggies are negative!'
                    log.critical(errmsg)
                    raise ValueError(errmsg)

                data['phot'] = CFit.parse_photometry(
                    CFit.bands, maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
                    lambda_eff=allfilters.effective_wavelengths.value,
                    min_uncertainty=CFit.min_uncertainty)
                
                # fiber fluxes
                #mw_transmission_fiberflux = 10**(-0.4 * ebv[igal] * CFit.RV * ext_odonnell(filters.effective_wavelengths.value, Rv=CFit.RV))
                if data['photsys'] != '':                
                    mw_transmission_fiberflux = np.array([mwdust_transmission(ebv[igal], band, data['photsys']) for band in CFit.fiber_bands])
                else:
                    mw_transmission_fiberflux = np.ones(len(CFit.fiber_bands))

                fibermaggies = np.zeros(len(CFit.fiber_bands))
                fibertotmaggies = np.zeros(len(CFit.fiber_bands))
                #ivarfibermaggies = np.zeros(len(CFit.fiber_bands))
                for iband, band in enumerate(CFit.fiber_bands):
                    fibermaggies[iband] = meta['FIBERFLUX_{}'.format(band.upper())][igal] / mw_transmission_fiberflux[iband]
                    fibertotmaggies[iband] = meta['FIBERTOTFLUX_{}'.format(band.upper())][igal] / mw_transmission_fiberflux[iband]
                    #ivarfibermaggies[iband] = meta['FIBERTOTFLUX_IVAR_{}'.format(band.upper())][igal] * mw_transmission_fiberflux[iband]**2

                data['fiberphot'] = CFit.parse_photometry(CFit.fiber_bands,
                    maggies=fibermaggies, nanomaggies=True,
                    lambda_eff=filters.effective_wavelengths.value)
                data['fibertotphot'] = CFit.parse_photometry(CFit.fiber_bands,
                    maggies=fibertotmaggies, nanomaggies=True,
                    lambda_eff=filters.effective_wavelengths.value)

                if not fastphot:
                    data.update({'wave': [], 'flux': [], 'ivar': [], 'res': [],
                                 'linemask': [], 'linemask_all': [],
                                 'linename': [], 'linepix': [], 'contpix': [],
                                 'snr': np.zeros(3, 'f4'),
                                 #'std': np.zeros(3, 'f4'), # emission-line free standard deviation, per-camera
                                 'cameras': spec.bands})
                    for icam, camera in enumerate(data['cameras']):
                        #mw_transmission_spec = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(spec.wave[camera], Rv=CFit.RV))
                        mw_transmission_spec = dust_transmission(spec.wave[camera], ebv[igal], Rv=CFit.RV)
                        data['wave'].append(spec.wave[camera])
                        data['flux'].append(spec.flux[camera][igal, :] / mw_transmission_spec)
                        data['ivar'].append(spec.ivar[camera][igal, :] * mw_transmission_spec**2)
                        data['res'].append(Resolution(spec.resolution_data[camera][igal, :, :]))
                        data['snr'][icam] = np.median(spec.flux[camera][igal, :] * np.sqrt(spec.ivar[camera][igal, :]))

                        #linemask = CFit.build_linemask(spec.wave[camera], redshift=data['zredrock'])
                        #data['linemask'].append(linemask)
                        #data['std'][icam] = np.std(spec.flux[camera][igal, :][linemask])

                    coadd_wave = coadd_spec.wave[coadd_bands]
                    coadd_flux = coadd_spec.flux[coadd_bands][igal, :]
                    coadd_ivar = coadd_spec.ivar[coadd_bands][igal, :]
                    coadd_res = Resolution(coadd_spec.resolution_data[coadd_bands][igal, :])

                    coadd_linemask_dict = CFit.build_linemask(coadd_wave, coadd_flux, coadd_ivar, redshift=data['zredrock'])
                    data['coadd_linename'] = coadd_linemask_dict['linename']
                    data['coadd_linepix'] = [np.where(lpix)[0] for lpix in coadd_linemask_dict['linepix']]
                    data['coadd_contpix'] = [np.where(cpix)[0] for cpix in coadd_linemask_dict['contpix']]
                    
                    data['linesigma_narrow'] = coadd_linemask_dict['linesigma_narrow']
                    data['linesigma_balmer'] = coadd_linemask_dict['linesigma_balmer']
                    data['linesigma_uv'] = coadd_linemask_dict['linesigma_uv']
                    
                    data['linesigma_narrow_snr'] = coadd_linemask_dict['linesigma_narrow_snr']
                    data['linesigma_balmer_snr'] = coadd_linemask_dict['linesigma_balmer_snr']
                    data['linesigma_uv_snr'] = coadd_linemask_dict['linesigma_uv_snr']

                    # Map the pixels belonging to individual emission lines and
                    # their local continuum back onto the original per-camera
                    # spectra. These lists of arrays are used in
                    # continuum.ContinnuumTools.smooth_continuum.
                    for icam in np.arange(len(data['cameras'])):
                        data['linemask'].append(np.interp(data['wave'][icam], coadd_wave, coadd_linemask_dict['linemask']*1) > 0)
                        data['linemask_all'].append(np.interp(data['wave'][icam], coadd_wave, coadd_linemask_dict['linemask_all']*1) > 0)
                        _linename, _linenpix, _contpix = [], [], []
                        for ipix in np.arange(len(coadd_linemask_dict['linepix'])):
                            I = np.interp(data['wave'][icam], coadd_wave, coadd_linemask_dict['linepix'][ipix]*1) > 0
                            J = np.interp(data['wave'][icam], coadd_wave, coadd_linemask_dict['contpix'][ipix]*1) > 0
                            #if '4686' in coadd_linemask_dict['linename'][ipix]:
                            #    pdb.set_trace()
                            if np.sum(I) > 3 and np.sum(J) > 3:
                                _linename.append(coadd_linemask_dict['linename'][ipix])
                                _linenpix.append(np.where(I)[0])
                                _contpix.append(np.where(J)[0])
                        data['linename'].append(_linename)
                        data['linepix'].append(_linenpix)
                        data['contpix'].append(_contpix)
                        #for ipix in np.arange(len(coadd_linemask_dict['contpix'])):
                        #    if icam == 1:
                        #        pdb.set_trace()
                        #    J = np.interp(data['wave'][icam], coadd_wave, coadd_linemask_dict['contpix'][ipix]*1) > 0
                        #    if np.sum(J) > 0:
                        #        _contpix.append(np.where(J)[0])
                        #data['contpix'].append(_contpix)
                        
                    #import matplotlib.pyplot as plt
                    #plt.clf()
                    #for ii in np.arange(3):
                    #    plt.plot(data['wave'][ii], data['flux'][ii])
                    #plt.plot(coadd_wave, coadd_flux-2, alpha=0.6, color='k')
                    #plt.xlim(5500, 6000)
                    #plt.savefig('test.png')
                    #pdb.set_trace()

                    if remember_coadd:
                        data.update({'coadd_wave': coadd_wave, 'coadd_flux': coadd_flux,
                                     'coadd_ivar': coadd_ivar, 'coadd_res': coadd_res,
                                     'coadd_linemask': coadd_linemask_dict['linemask'],
                                     'coadd_linemask_all': coadd_linemask_dict['linemask_all']})
                                     #'linepix': np.where(coadd_linemask_dict['linepix'])[0],
                                     #'contpix': np.where(coadd_linemask_dict['contpix'])[0]})

                    # Optionally synthesize photometry from the coadded spectrum.
                    if synthphot:
                        padflux, padwave = filters.pad_spectrum(coadd_flux, coadd_wave, method='edge')
                        synthmaggies = filters.get_ab_maggies(padflux / CFit.fluxnorm, padwave)
                        synthmaggies = synthmaggies.as_array().view('f8')

                        # code to synthesize uncertainties from the variance spectrum
                        #var, mask = _ivar2var(data['coadd_ivar'])
                        #padvar, padwave = filters.pad_spectrum(var[mask], data['coadd_wave'][mask], method='edge')
                        #synthvarmaggies = filters.get_ab_maggies(1e-17**2 * padvar, padwave)
                        #synthivarmaggies = 1 / synthvarmaggies.as_array().view('f8')[:3] # keep just grz

                        #data['synthphot'] = CFit.parse_photometry(CFit.bands,
                        #    maggies=synthmaggies, lambda_eff=lambda_eff[:3],
                        #    ivarmaggies=synthivarmaggies, nanomaggies=False)

                        data['synthphot'] = CFit.parse_photometry(CFit.synth_bands,
                            maggies=synthmaggies, nanomaggies=False,
                            lambda_eff=filters.effective_wavelengths.value)

                alldata.append(data)

        self.meta = vstack(self.meta)
        #self.meta = Table(np.hstack(self.meta))
        self.ntargets = len(self.meta)
        log.info('Read data for {} objects in {:.2f} sec'.format(self.ntargets, time.time()-t0))

        return alldata
        
    def init_output(self, CFit=None, EMFit=None, fastphot=False):
        """Initialize the fastspecfit output data table.

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
        CFit : :class:`fastspecfit.continuum.ContinuumFit`
            Continuum-fitting class.
        EMFit : :class:`fastspecfit.emlines.EMLineFit`
            Emission-line fitting class.

        Returns
        -------


        Notes
        -----

        """
        import astropy.units as u
        from astropy.table import hstack, Column

        nobj = len(self.meta)

        # The information stored in the metadata table depends on which spectra
        # were fitted (exposures, nightly coadds, deep coadds).
        fluxcols = ['PHOTSYS', 'LS_ID',
                    #'RELEASE',
                    'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',
                    'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
                    'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4',
                    'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                    'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FLUX_IVAR_W3', 'FLUX_IVAR_W4',
                    'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
                    'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 'MW_TRANSMISSION_W3', 'MW_TRANSMISSION_W4']
            
        colunit = {'RA': u.deg, 'DEC': u.deg,
                   'FIBERFLUX_G': 'nanomaggies', 'FIBERFLUX_R': 'nanomaggies', 'FIBERFLUX_Z': 'nanomaggies',
                   'FIBERTOTFLUX_G': 'nanomaggies', 'FIBERTOTFLUX_R': 'nanomaggies', 'FIBERTOTFLUX_Z': 'nanomaggies',
                   'FLUX_G': 'nanomaggies', 'FLUX_R': 'nanomaggies', 'FLUX_Z': 'nanomaggies',
                   'FLUX_W1': 'nanomaggies', 'FLUX_W2': 'nanomaggies', 'FLUX_W3': 'nanomaggies', 'FLUX_W4': 'nanomaggies', 
                   'FLUX_IVAR_G': 'nanomaggies-2', 'FLUX_IVAR_R': 'nanomaggies-2',
                   'FLUX_IVAR_Z': 'nanomaggies-2', 'FLUX_IVAR_W1': 'nanomaggies-2',
                   'FLUX_IVAR_W2': 'nanomaggies-2', 'FLUX_IVAR_W3': 'nanomaggies-2',
                   'FLUX_IVAR_W4': 'nanomaggies-2',
                   }

        skipcols = ['OBJTYPE', 'TARGET_RA', 'TARGET_DEC', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID', 'RELEASE'] + fluxcols
        redrockcols = ['Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE', 'Z_RR']
        
        meta = Table()
        metacols = self.meta.colnames

        # All of this business is so we can get the columns in the order we want
        # (i.e., the order that matches the data model).
        for metacol in ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'FIBER',
                        'NIGHT', 'TILEID_LIST', 'RA', 'DEC', 'COADD_FIBERSTATUS']:
            if metacol in metacols:
                meta[metacol] = self.meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]

        if self.specprod == 'fuji': # EDR
            TARGETINGCOLS = TARGETINGBITS[self.specprod]
        else:
            TARGETINGCOLS = TARGETINGBITS['default']

        for metacol in metacols:
            if metacol in skipcols or metacol in TARGETINGCOLS or metacol in meta.colnames or metacol in redrockcols:
                continue
            else:
                meta[metacol] = self.meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]

        for bitcol in TARGETINGCOLS:
            if bitcol in metacols:
                meta[bitcol] = self.meta[bitcol]
            else:
                meta[bitcol] = np.zeros(shape=(1,), dtype=np.int64)

        for redrockcol in redrockcols:
            if redrockcol in metacols: # the Z_RR from quasarnet may not be present
                meta[redrockcol] = self.meta[redrockcol]
            if redrockcol in colunit.keys():
                meta[redrockcol].unit = colunit[redrockcol]

        for fluxcol in fluxcols:
            meta[fluxcol] = self.meta[fluxcol]
            if fluxcol in colunit.keys():
                meta[fluxcol].unit = colunit[fluxcol]

        out = Table()
        for col in ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'TILEID', 'FIBER', 'NIGHT']:
            if col in metacols:
                out[col] = self.meta[col]
        if fastphot:
            out = hstack((out, CFit.init_phot_output(nobj)))
        else:
            out = hstack((out, CFit.init_spec_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

        return out, meta

def read_fastspecfit(fastfitfile, rows=None, columns=None):
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
        log.info('Read {} object(s) from {}'.format(len(fastfit), fastfitfile))

        # Add specprod to the metadata table so that we can stack across
        # productions (e.g., Fuji+Guadalupe).
        hdr = fitsio.read_header(fastfitfile, ext='PRIMARY')
        if 'SPECPROD' in hdr:
            specprod = hdr['SPECPROD']
            meta['SPECPROD'] = specprod
            
        if 'COADDTYP' in hdr:
            coadd_type = hdr['COADDTYP']
        else:
            coadd_type = None

        return fastfit, meta, coadd_type, fastphot
    
    else:
        log.warning('File {} not found.'.format(fastfitfile))
        return None, None, None

def write_fastspecfit(out, meta, modelspectra=None, outfile=None, specprod=None,
                      coadd_type=None, fastphot=False):
    """Write out.

    """
    import gzip, shutil
    from astropy.io import fits
    from desispec.io.util import fitsheader
    from desiutil.depend import add_dependencies

    t0 = time.time()
    outdir = os.path.dirname(os.path.abspath(outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    log.info('Writing results for {} objects to {}'.format(len(out), outfile))
    
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
        primhdr.append(('COADDTYP', (coadd_type, 'spectral coadd fitted')))

    primhdr = fitsheader(primhdr)
    add_dependencies(primhdr)
    add_dependencies(primhdr, module_names=['fastspecfit'])
    
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

    log.info('Writing out took {:.2f} sec'.format(time.time()-t0))

def select(fastfit, metadata, coadd_type, healpixels=None, tiles=None, nights=None):
    """Optionally trim to a particular healpix or tile and/or night."""
    keep = np.ones(len(fastfit), bool)
    if coadd_type == 'healpix':
        if healpixels:
            pixelkeep = np.zeros(len(fastfit), bool)
            for healpixel in healpixels:
                pixelkeep = np.logical_or(pixelkeep, metadata['HEALPIX'].astype(str) == healpixel)
            keep = np.logical_and(keep, pixelkeep)
            log.info('Keeping {} objects from healpixels(s) {}'.format(len(fastfit), ','.join(healpixels)))
    else:
        if tiles:
            tilekeep = np.zeros(len(fastfit), bool)
            for tile in tiles:
                tilekeep = np.logical_or(tilekeep, metadata['TILEID'].astype(str) == tile)
            keep = np.logical_and(keep, tilekeep)
            log.info('Keeping {} objects from tile(s) {}'.format(len(fastfit), ','.join(tiles)))
        if nights and 'NIGHT' in metadata:
            nightkeep = np.zeros(len(fastfit), bool)
            for night in nights:
                nightkeep = np.logical_or(nightkeep, metadata['NIGHT'].astype(str) == night)
            keep = np.logical_and(keep, nightkeep)
            log.info('Keeping {} objects from night(s) {}'.format(len(fastfit), ','.join(nights)))
    return fastfit[keep], metadata[keep]
