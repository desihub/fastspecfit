#!/usr/bin/env python
"""
fastspecfit.io
==============

Tools for reading and writing.

"""
import pdb # for debugging

import os, time
import numpy as np
import fitsio
from astropy.table import Table, vstack

from desiutil.log import get_logger
log = get_logger()

TARGETINGBITCOLS = [
    #'CMX_TARGET',
    'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET',
    'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET',
    'SV2_DESI_TARGET', 'SV2_BGS_TARGET', 'SV2_MWS_TARGET',
    'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET',
    'SCND_TARGET',
    'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET',
    ]

TARGETCOLS = np.array([
    'PHOTSYS',
    'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
    'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
    'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
    'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
    'FLUX_IVAR_W1', 'FLUX_IVAR_W2'])#,
    #'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
    #'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']

DESI_ROOT_NERSC = '/global/cfs/cdirs/desi'
DUST_DIR_NERSC = '/global/cfs/cdirs/cosmo/data/dust/v0_1'
FASTSPECFIT_TEMPLATES_NERSC = '/global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z'

#SCDNDIR = {'81097': '/global/cfs/cdirs/desi/target/catalogs/dr9/0.53.0/targets/sv2/resolve/dark'}

#VETO_TILES = np.array([81097, 81098, 81099])

class DESISpectra(object):
    def __init__(self, specprod=None):
        """Class to read in the DESI data needed by fastspecfit.

        """
        desi_root = os.environ.get('DESI_ROOT', DESI_ROOT_NERSC)

        if specprod is None:
            log.warning('specprod input is required.')
            raise IOError

        self.specprod = specprod
        
        self.redux_dir = os.path.join(desi_root, 'spectro', 'redux')
        self.fiberassign_dir = os.path.join(desi_root, 'target', 'fiberassign', 'tiles', 'trunk')

        #self.healpix_dir = os.path.join(self.redux_dir, self.specprod, 'healpix')
        #self.tiles_dir = os.path.join(self.redux_dir, self.specprod, 'tiles', self.coadd_type)

    def _get_targetdirs(self, tileid):
        """Get the targets catalog used to build a given fiberassign catalog.

        """
        from astropy.io import fits
        #thistile = self.tiles[self.tiles['TILEID'] == tileid]
        stileid = '{:06d}'.format(tileid)
        fiberfile = os.path.join(self.fiberassign_dir, stileid[:3], 'fiberassign-{}.fits.gz'.format(stileid))
        if not os.path.isfile(fiberfile):
            fiberfile = fiberfile.replace('.gz', '')
            if not os.path.isfile(fiberfile):
                log.warning('Fiber assignment file {} not found!'.format(fiberfile))
        log.info('Reading {} header.'.format(fiberfile))
        # fitsio can't handle CONTINUE header cards!
        #fahdr = fitsio.read_header(fiberfile, ext=0)
        # fastspec /global/cfs/cdirs/desi/spectro/redux/blanc/tiles/80605/20201222/redrock-6-80605-20201222.fits -o /global/cfs/cdirs/desi/spectro/fastspecfit/blanc/tiles/80605/20201222/fastspec-6-80605-20201222.fits --mp 32
        fahdr = fits.getheader(fiberfile, ext=0)
        targetdirs = [fahdr['TARG']]
        for moretarg in ['TARG2', 'TARG3', 'TARG4']:
            if moretarg in fahdr:
                if 'gaia' not in fahdr[moretarg]: # skip
                    targetdirs += [fahdr[moretarg]]
        if 'SCND' in fahdr:
            if fahdr['SCND'].strip() != '-':
                #if tileid in SCDNDIR.keys():
                #    targetdirs += [SCDNDIR[tileid]]
                #else:
                targetdirs += [fahdr['SCND']]

        desi_root = os.environ.get('DESI_ROOT', DESI_ROOT_NERSC)
        for ii, targetdir in enumerate(targetdirs):
            # for secondary targets, targetdir can be a filename
            #if os.path.isfile(targetdir):
            if targetdir[-4:] == 'fits': # fragile...
                targetdir = os.path.dirname(targetdir)
            if not os.path.isdir(targetdir):
                #log.warning('Targets directory not found {}'.format(targetdir))
                # can be a KPNO directory!
                if 'DESIROOT' in targetdir:
                    targetdir = os.path.join(desi_root, targetdir.replace('DESIROOT/', ''))
                if targetdir[:6] == '/data/':
                    targetdir = os.path.join(desi_root, targetdir.replace('/data/', ''))
                
            if os.path.isdir(targetdir):
                log.info('Found targets directory {}'.format(targetdir))
                targetdirs[ii] = targetdir
            else:
                log.warning('Targets directory {} not found.'.format(targetdir))
                #pdb.set_trace()
                #raise IOError
                continue

        # any ToOs?
        if 'TOO' in fahdr:
            # special case, fragile!
            if 'sv3' in fahdr['TOO']:
                ddir = 'sv3'
            else:
                ddir = 'main'
            TOOfile = os.path.join(desi_root, 'survey', 'ops', 'surveyops', 'trunk', 'mtl', ddir, 'ToO', 'ToO.ecsv')
            if os.path.isfile(TOOfile):
                log.info('Found ToO file {}'.format(TOOfile))
            else:
                log.warning('ToO file {} not found.'.format(TOOfile))
                TOOfile = None                
        else:
            TOOfile = None
            
        return targetdirs, TOOfile

    def find_specfiles(self, redrockfiles=None, firsttarget=0, zmin=0.001, zmax=None,
                       zwarnmax=0, coadd_type=None, targetids=None, ntargets=None):
        """Initialize the fastspecfit output data table.

        Parameters
        ----------

        Returns
        -------

        Notes
        -----
        fastfit - results table (overrides redrockfiles and then specprod is required

        """
        from glob import glob
        from astropy.table import Column
        import numpy.ma as ma
        from desimodel.footprint import radec2pix

        if zmax is None:
            zmax = 99.0
        
        if redrockfiles is None:
            log.warning('At least one redrockfiles file is required.')
            raise ValueError
        #if redrockfiles is None and fastfit is None and metadata is None:
        #    log.warning('At least one of redrockfiles or fastfit (and metadata, specprod, and coadd_type) are required.')
        #    raise ValueError

        if len(np.atleast_1d(redrockfiles)) == 0:
            log.warning('No redrockfiles found!')
            raise IOError

        ## Try to glean specprod so we can write it to the output file. This
        ## should really be in the file header--- see
        ## https://github.com/desihub/desispec/issues/1077
        #if specprod is None:            
        #    #import desiutil.depend
        #    #hdr = fitsio.read_header(np.atleast_1d(redrockfiles)[0].replace('redrock-', 'coadd-'))
        #    # stupidly fragile!
        #    specprod = np.atleast_1d(redrockfiles)[0].replace(self.redux_dir, '').split(os.sep)[1]
        #    self.specprod = specprod
        #    log.info('Parsed specprod={}'.format(specprod))

        # Should we not sort...?
        #redrockfiles = np.array(set(np.atleast_1d(redrockfiles)))
        redrockfiles = np.array(sorted(set(np.atleast_1d(redrockfiles))))
        log.info('Reading and parsing {} unique redrockfile(s)'.format(len(redrockfiles)))

        self.redrock, self.meta, self.tiles = [], [], []
        self.redrockfiles, self.specfiles = [], []
        for ired, redrockfile in enumerate(np.atleast_1d(redrockfiles)):
            specfile = redrockfile.replace('redrock-', 'coadd-')
            if not os.path.isfile(redrockfile):
                log.warning('File {} not found!'.format(redrockfile))
                continue
            if not os.path.isfile(specfile):
                log.warning('File {} not found!'.format(specfile))
                continue

            # Try to figure out coadd_type from the first filename. Fragile!
            # Assumes that coadd_type is a scalar... And, really, this should be
            # in a header.
            if coadd_type is not None:
                self.coadd_type = coadd_type
            else:
                if ired == 0:
                    hdr = fitsio.read_header(specfile, ext=0)
                    if 'HPXPIXEL' in hdr:
                        coadd_type = 'healpix'
                        hpxpixel = np.int32(hdr['HPXPIXEL'])
                    else:
                        import re
                        if re.search('-thru20[0-9]+[0-9]+[0-9]+\.fits', redrockfile) is not None:
                            coadd_type = 'cumulative'
                        elif re.search('-20[0-9]+[0-9]+[0-9]+\.fits', redrockfile) is not None:
                            coadd_type = 'pernight'
                        else:
                            coadd_type = 'perexp'
                    self.coadd_type = coadd_type
            log.info('Parsed coadd_type={}'.format(coadd_type))
    
            # Figure out which fibermap columns to put into the metadata
            # table. Note that the fibermap includes all the spectra that went
            # into the coadd (based on the unique TARGETID which is in the redrock
            # table).  See https://github.com/desihub/desispec/issues/1104
            allfmcols = np.array(fitsio.FITS(specfile)['FIBERMAP'].get_colnames())
            fmcols = ['TARGETID', 'TARGET_RA', 'TARGET_DEC',
                      'COADD_FIBERSTATUS', 'OBJTYPE']#,'PHOTSYS']
                      #'PHOTSYS', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
                      #'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
                      #'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
                      #'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                      #'FLUX_IVAR_W1', 'FLUX_IVAR_W2']
            if coadd_type == 'custom':
                expfmcols = ['TARGETID', 'FIBER'] # hack for M31 project
            else:
                expfmcols = ['TARGETID', 'TILEID', 'FIBER']
            
            # add targeting columns
            HAVE_TARGETINGBITCOLS = [col for col in TARGETINGBITCOLS if col in allfmcols]
            fmcols = np.array(fmcols + HAVE_TARGETINGBITCOLS).tolist()
            morecols = np.isin(TARGETCOLS, allfmcols)
            if np.sum(morecols) > 0:
                fmcols = np.hstack((fmcols, TARGETCOLS[morecols]))

            # For the cumulative coadds, NIGHT is defined to be the last night
            # contributing to the coadd.
            if self.coadd_type == 'custom':
                thrunight = None
            elif self.coadd_type == 'healpix':
                thrunight = None
            elif self.coadd_type == 'cumulative':
                thrunight = np.int32(os.path.basename(os.path.dirname(specfile)))
            else:
                thrunight = None
                expfmcols = expfmcols + ['NIGHT']
                if self.coadd_type == 'perexp':
                    expfmcols = expfmcols + ['EXPID']

            ## older fibermaps are missing the WISE inverse variance
            #if 'FLUX_IVAR_W1' in allfmcols:
            #    fmcols = fmcols + ['FLUX_IVAR_W1', 'FLUX_IVAR_W2']

            zbcols = ['TARGETID', 'Z', 'ZWARN', 'SPECTYPE', 'DELTACHI2']

            # If targetids is *not* given we have to choose "good" objects
            # before subselecting (e.g., we don't want sky spectra).
            if targetids is None:
                zb = fitsio.read(redrockfile, 'REDSHIFTS', columns=zbcols)
                # Are we reading individual exposures or coadds?
                meta = fitsio.read(specfile, 'FIBERMAP', columns=fmcols)
                assert(np.all(zb['TARGETID'] == meta['TARGETID']))
                fitindx = np.where(
                    (zb['Z'] > zmin) *
                    (zb['Z'] < zmax) *
                    (zb['ZWARN'] <= zwarnmax) * 
                    (meta['OBJTYPE'] == 'TGT') *
                    #(meta['PHOTSYS'] != 'G') * #(zb['SPECTYPE'] == 'GALAXY') *
                    (meta['COADD_FIBERSTATUS'] == 0)
                    )[0]

                # bug in ToOs -- they don't get targeting info propagated 
                # fastphot /global/cfs/cdirs/desi/spectro/redux/everest/healpix/sv3/bright/259/25969/redrock-sv3-bright-25969.fits -o fastphot.fits --targetids
                check = np.sum(meta[HAVE_TARGETINGBITCOLS][fitindx].tolist(), axis=1)
                good = check != 0
                bad = check == 0
                if np.sum(bad) > 0:
                    tbad = [str(tid) for tid in meta['TARGETID'][fitindx][bad]]
                    log.warning('Missing _TARGET info for the following objects: {}'.format(' '.join(tbad)))
                #pdb.set_trace()
                if np.sum(good) > 0:
                    fitindx = fitindx[good]
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
                zb = Table(fitsio.read(redrockfile, 'REDSHIFTS', rows=fitindx, columns=zbcols))
                meta = Table(fitsio.read(specfile, 'FIBERMAP', rows=fitindx, columns=fmcols))

            assert(np.all(zb['TARGETID'] == meta['TARGETID']))

            # Get the unique set of tiles contributing to the coadded spectra from EXP_FIBERMAP
            expmeta = fitsio.read(specfile, 'EXP_FIBERMAP', columns=expfmcols)
            I = np.isin(expmeta['TARGETID'], meta['TARGETID'])
            if np.count_nonzero(I) == 0:
                log.warning('No matching targets in exposure table.')
                raise ValueError
            expmeta = Table(expmeta[I])
            if coadd_type == 'custom':
                expmeta['TILEID'] = 80715 # M31 hack!

            tiles = np.unique(np.atleast_1d(expmeta['TILEID']).data)

            ## this code is good for vetoing specific tiles
            #if False:
            #    if np.any(np.isin(tiles, VETO_TILES)):
            #        remtargets = np.array(expmeta[np.isin(expmeta['TILEID'], VETO_TILES)]['TARGETID'].tolist())
            #        remindx = np.isin(meta['TARGETID'], remtargets)
            #        if np.sum(remindx) > 0:
            #            log.info('Removing {} targets from VETO_TILES.'.format(np.sum(remindx)))
            #            keepindx = np.where(np.logical_not(remindx))[0]
            #            meta = meta[keepindx]
            #            if len(meta) == 0:
            #                log.info('All targets have been removed from redrockfile {}'.format(redrockfile))
            #                continue
            #            zb = zb[keepindx]
            #            fitindx = fitindx[keepindx]
            #            tiles = tiles[np.logical_not(np.isin(tiles, VETO_TILES))]

            if thrunight:
                meta['THRUNIGHT'] = thrunight

            # Gather additional info about this pixel.
            if coadd_type == 'healpix':
                pixinfo = os.path.basename(redrockfile).split('-') 
                meta['SURVEY'] = pixinfo[1] # a little fragile...
                meta['FAPRGRM'] = pixinfo[2]
                meta['HPXPIXEL'] = hpxpixel
            else:
                # Initialize the columns to get the data type right and then
                # populate them by targetid.
                meta['TILEID'] = expmeta['TILEID'][0]
                meta['FIBER'] = expmeta['FIBER'][0]
                if coadd_type == 'custom':
                    pass
                else:
                    if coadd_type != 'cumulative':
                        meta['NIGHT'] = expmeta['NIGHT'][0]
                    if coadd_type == 'perexp':
                        meta['EXPID'] = expmeta['EXPID'][0]
                        
                for iobj, tid in enumerate(meta['TARGETID']):
                    iexp = np.where(expmeta['TARGETID'] == tid)[0][0] # zeroth
                    meta['TILEID'][iobj] = expmeta['TILEID'][iexp]
                    meta['FIBER'][iobj] = expmeta['FIBER'][iexp]
                    if coadd_type == 'custom':
                        pass
                    else:
                        if coadd_type != 'cumulative':
                            meta['NIGHT'][iobj] = expmeta['NIGHT'][iexp]
                        if coadd_type == 'perexp':
                            meta['EXPID'][iobj] = expmeta['EXPID'][iexp]
    
            self.redrock.append(Table(zb))
            self.meta.append(Table(meta))
            self.tiles.append(tiles)
            self.redrockfiles.append(redrockfile)
            self.specfiles.append(specfile)

        if len(self.redrock) == 0:
            log.warning('No targets read!')
            return

        # The targets catalogs are organized on a much larger spatial scale, so
        # grab additional targeting info outside the main loop. See
        # desitarget.io.read_targets_in_tiles for the algorithm.
        t0 = time.time()

        targetcols = ['TARGETID', 'RA', 'DEC']
        #if not 'FLUX_IVAR_W1' in fmcols:
        #    targetcols = targetcols + ['FLUX_IVAR_W1', 'FLUX_IVAR_W2']
        #targetcols = targetcols + [
        #    'PHOTSYS',
        #    'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
        #    'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
        #    'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
        #    'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
        #    'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
        #    'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
        #    'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']

        info = Table(np.hstack([meta['TARGETID', 'TARGET_RA', 'TARGET_DEC'] for meta in self.meta]))
        
        alltileid = np.unique(self.tiles)
        targetdirs, TOOfiles = [], []
        for tileid in alltileid:
            _targetdirs, TOOfile = self._get_targetdirs(tileid)
            targetdirs.append(_targetdirs)
            if TOOfile is not None:
                TOOfiles.append(TOOfile)

        targets, scndtargets = [], []
        for targetdir in set(np.hstack(targetdirs)):
            # Handle secondary targets, which have a different data model;
            # update on 2021 July 31: these catalogs are missing DR9
            # photometry, so we have to skip them for now.
            if 'secondary' in targetdir:
                #continue                    
                if 'sv1' in targetdir: # special case
                    if 'dedicated' in targetdir:
                        targetfiles = glob(os.path.join(targetdir, 'DC3R2_GAMA_priorities.fits'))
                    else:
                        #targetfiles = ['/global/cfs/cdirs/desi/users/raga19/data/0.51.0/sv1targets_dark_secondary_dr9.fits']
                        targetfiles = glob(os.path.join(targetdir, '*-secondary-dr9photometry.fits'))
                else:
                    targetfiles = glob(os.path.join(targetdir, '*-secondary.fits'))
            else:
                targetfiles = glob(os.path.join(targetdir, '*-hp-*.fits'))
                filenside = fitsio.read_header(targetfiles[0], ext=1)['FILENSID']
                pixlist = radec2pix(filenside, info['TARGET_RA'], info['TARGET_DEC'])
                targetfiles = [targetfiles[0].split('hp-')[0]+'hp-{}.fits'.format(pix) for pix in set(pixlist)]

            for ifile, targetfile in enumerate(targetfiles):
                tinfo = fitsio.FITS(targetfile)
                for _tinfo in tinfo:
                    extname = _tinfo.get_extname()
                    if 'TARGETS' in extname:
                        break
                if extname == '':
                    extname = 1
                alltargetcols = np.array(tinfo[extname].get_colnames())

                # secondary target catalogs are missing some or all of the
                # DR9 photometry columns we need
                # fastspec /global/cfs/cdirs/desi/spectro/redux/everest/healpix/sv3/bright/153/15343/redrock-sv3-bright-15343.fits -o fastspec
                thesecols = np.isin(alltargetcols, TARGETCOLS)
                if np.sum(thesecols) == 0:
                    log.warning('Problem reading target catalog {}'.format(targetfile))

                readtargetcols = np.hstack((targetcols, alltargetcols[thesecols]))
                #readtargetcols = np.hstack((targetcols, TARGETCOLS[thesecols]))
                alltargets = tinfo[extname].read(columns=readtargetcols)
                match = np.isin(alltargets['TARGETID'], info['TARGETID'])
                log.info('Matched {} targets in {}'.format(np.sum(match), targetfile))

                if np.sum(match) > 0:
                    # need to make sure the column order is the same                        
                    if len(targets) > 0 and not np.all(targets[len(targets)-1].dtype.names == alltargets.dtype.names):
                    #if not np.all(readtargetcols == np.array(alltargets.dtype.names)):                            
                        #_alltargets = Table(alltargets[match])[readtargetcols].as_array() # this doesn't work
                        _alltargets = Table(alltargets[match])
                        _newtargets = Table()
                        I = len(targets) - 1
                        for col in targets[I].dtype.names:
                            if col in _alltargets.colnames:
                                _newtargets[col] = _alltargets[col]
                            else:
                                _newtargets.add_column(Column(name=col, dtype=targets[I][col].dtype, length=len(_alltargets)))
                        alltargets = _newtargets.as_array()
                    else:
                        alltargets = alltargets[match]

                    # handle secondary target catalogs in SV3 and main
                    if np.sum(thesecols) < len(TARGETCOLS):
                        scndtargets.append(alltargets)
                    else:
                        targets.append(alltargets)

        # handle ToO targets
        TOOtargets = []
        if len(TOOfiles) > 0:
            for TOOfile in set(np.hstack(TOOfiles)):
                log.info('Reading {}'.format(TOOfile))
                TOO = Table.read(TOOfile) # ecsv file
                
                match = np.isin(TOO['TARGETID'], info['TARGETID'])
                log.info('Matched {} ToOs in {}'.format(np.sum(match), TOOfile))
                if np.sum(match) > 0:
                    # build a Table that looks like the target catalog
                    TOOmatch = TOO[match]
                    nt = len(targets)
                    if nt > 0:
                        _TOOtargets = Table()
                        for col in targets[nt-1].dtype.names:
                            if col in TOOmatch.colnames:
                                _TOOtargets[col] = TOOmatch[col]
                            else:
                                _TOOtargets.add_column(Column(name=col, dtype=targets[nt-1][col].dtype, length=len(TOOmatch)))
                        #TOOtargets.append(_TOOtargets.as_array())
                        TOOtargets.append(_TOOtargets)
                    else:
                        TOOtargets.append(TOOmatch)

        if len(targets) > 0:
            targets = Table(np.hstack(targets))
                
        if len(scndtargets) > 0:
            scndtargets = Table(np.hstack(scndtargets))
            if len(targets) == 0:
                targets = scndtargets
            
        if len(TOOtargets) > 0:
            TOOtargets = vstack(TOOtargets)
            #TOOtargets = Table(np.hstack(TOOtargets))
            if len(targets) > 0:
                targets = vstack((targets, TOOtargets))
            else:
                targets = TOOtargets
                
        #from desitarget.io import releasedict
        #from desitarget.targets import decode_targetid  
        #_, _, releases, _, _, _ = decode_targetid(targets['TARGETID'])  
        #photsys = [releasedict[release] if release >= 9000 else None for release in releases]

        if len(scndtargets) > 0:
            _, uindx = np.unique(scndtargets['TARGETID'], return_index=True) 
            scndtargets = scndtargets[uindx]

            # add the missing columns, if any
            _scndtargets = Table()
            for col in targets.colnames:
                if col in scndtargets.colnames:
                    _scndtargets[col] = scndtargets[col]
                else:
                    _scndtargets.add_column(Column(name=col, dtype=targets[col].dtype, length=len(scndtargets)))
            #scndtargets = _scndtargets.as_array()
            targets = vstack((targets, scndtargets))

        # targets table can include duplicates...
        _, uindx = np.unique(targets['TARGETID'], return_index=True) 
        targets = targets[uindx]

        if len(targets) != len(info):
            miss = np.delete(np.arange(len(info)), np.where(np.isin(info['TARGETID'], targets['TARGETID']))[0])
            log.warning('Missing targeting info for {} objects!'.format(len(info) - len(targets)))
            #_, _, release, _, _, _ = decode_targetid(info[miss])
            raise ValueError

        # convert masked and negative fluxes and ivarfluxes to zero
        for band in ['G', 'R', 'Z', 'W1', 'W2']:
            for prefix in ['', 'FIBER', 'FIBERTOT']:
                for suffix in ['', '_IVAR']:
                    col = '{}FLUX{}_{}'.format(prefix, suffix, band)
                    if col in targets.colnames:
                        if ma.is_masked(targets[col]):
                            targets[col] = ma.filled(targets[col], 0.0)
                        bad = np.logical_or(targets[col] < 0, np.isnan(targets[col]))
                        if np.sum(bad) > 0:
                            targets[col][bad] = 0.0
        if ma.is_masked(targets['PHOTSYS']):
            targets['PHOTSYS'] = ma.filled(targets['PHOTSYS'], '')
                        
        metas = []
        for meta in self.meta:
            srt = np.hstack([np.where(tid == targets['TARGETID'])[0] for tid in meta['TARGETID']])
            assert(np.all(meta['TARGETID'] == targets['TARGETID'][srt]))
            # prefer the target catalog quantities over those in the fiberassign table
            for col in targets.dtype.names:
                meta[col] = targets[col][srt]
                #if col not in meta.colnames:
                #    meta[col] = targets[col][srt]
            meta = Table(meta)
            #pdb.set_trace()
            nobj = len(meta)
            # placeholder (to be added in DESISpectra.read_and_unpack)
            for band in ['G', 'R', 'Z', 'W1', 'W2']:
                meta.add_column(Column(name='MW_TRANSMISSION_{}'.format(band), data=np.ones(nobj, 'f4')))
            metas.append(meta)
        log.info('Read and parsed targeting info for {} objects in {:.2f} sec'.format(len(targets), time.time()-t0))

        self.meta = metas # update

    #@staticmethod
    #def desitarget_resolve_dec():
    #    """Default Dec cut to separate targets in BASS/MzLS from DECaLS."""
    #    dec = 32.375
    #    return dec

    def read_and_unpack(self, CFit, fastphot=False, synthphot=True,
                        remember_coadd=False):
        """Unpack and pre-process a single DESI spectrum.
        
        Parameters
        ----------
        specobj : :class:`desispec.spectra.Spectra`
            DESI spectra (in the standard format).
        redrock : :class:`astropy.table.Table`
            Redrock redshift table (row-aligned to `specobj`).
        CFit : :class:`fastspecfit.continuum.ContinuumFit`
            Continuum-fitting class which contains filter curves and some additional
            photometric convenience functions.
        indx : :class:`int`
            Index number (0-indexed) of the spectrum to unpack and pre-process.

        Returns
        -------
        :class:`dict` with the following keys:
            zredrock : :class:`numpy.float64`
                Redrock redshift.
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
                Three-element list of `desispec.resolution.Resolution` objects, one
                for each camera.
            snr : :class:`numpy.ndarray`
                Median per-pixel signal-to-noise ratio in the grz cameras.
            linemask : :class:`list`
                Three-element list of `numpy.ndarray` boolean emission-line masks,
                one for each camera. This mask is used during continuum-fitting.
            coadd_wave : :class:`numpy.ndarray`
                Coadded wavelength vector with all three cameras combined.
            coadd_flux : :class:`numpy.ndarray`
                Flux corresponding to `coadd_wave`.
            coadd_ivar : :class:`numpy.ndarray`
                Inverse variance corresponding to `coadd_flux`.
            photsys_south : :class:`bool`
                Boolean indicating whether this object is on the south (True) or
                north (False) photometric system based on the declination cut coded
                in `desitarget.io.desispec.resolution.Resolution`.
            phot : :class:`astropy.table.Table`
                Imaging photometry in `grzW1W2`, corrected for Milky Way extinction.
            synthphot : :class:`astropy.table.Table`
                Photometry in `grz` synthesized from the extinction-corrected
                coadded spectra (with a mild extrapolation of the data blueward and
                redward to accommodate the g-band and z-band filter curves,
                respectively.

        Notes
        -----
        Hard-coded to assume that all three cameras (grz) have spectra.

        """
        #from scipy.interpolate import interp1d        
        from desispec.resolution import Resolution
        from desispec.coaddition import coadd_cameras
        from desispec.io import read_spectra # read_tile_spectra # 
        from desiutil.dust import mwdust_transmission, dust_transmission#, ext_odonnell
        from fastspecfit.util import C_LIGHT

        # Read everything into a simple dictionary.
        t0 = time.time()
        alldata = []
        for ispec, (specfile, redrock, meta) in enumerate(zip(self.specfiles, self.redrock, self.meta)):
            log.info('Reading {} spectra from {}'.format(len(redrock), specfile))

            # sometimes these are an astropy.table.Row!
            redrock = Table(redrock)
            meta = Table(meta)

            ebv = CFit.SFDMap.ebv(meta['RA'], meta['DEC'])

            if not fastphot:
                #spec = read_tile_spectra(meta['TILEID'][0], self.coadd_type, specprod=self.specprod,
                #                         targets=redrock['TARGETID'])
                spec = read_spectra(specfile).select(targets=redrock['TARGETID'])
                assert(np.all(spec.fibermap['TARGETID'] == redrock['TARGETID']))
                assert(np.all(spec.fibermap['TARGETID'] == meta['TARGETID']))

                # Coadd across cameras.
                t1 = time.time()                
                coadd_spec = coadd_cameras(spec)
                coadd_bands = coadd_spec.bands[0]
                log.info('Coadded the cameras in {:.2f} sec'.format(time.time()-t1))
                
            for igal in np.arange(len(redrock)):
                # Unpack the data and correct for Galactic extinction. Also flag pixels that
                # may be affected by emission lines.
                data = {'targetid': meta['TARGETID'][igal], 'zredrock': redrock['Z'][igal],
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
                    log.warning('Some ivarmaggies are negative!')
                    raise ValueError

                data['phot'] = CFit.parse_photometry(CFit.bands,
                    maggies=maggies, ivarmaggies=ivarmaggies, nanomaggies=True,
                    lambda_eff=allfilters.effective_wavelengths.value)

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
                                 'linemask': [], 'linename': [], 'linepix': [], 'contpix': [],
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

                    # Coadd across cameras.
                    #coadd_wave = np.unique(np.hstack(data['wave']))
                    #coadd_flux3d = np.zeros((len(coadd_wave), 3))
                    #coadd_ivar3d = np.zeros_like(coadd_flux3d)
                    #coadd_linemask3d = np.ones((len(coadd_wave), 3), bool)
                    #for icam in np.arange(len(data['cameras'])):
                    #    I = np.where(np.isin(data['wave'][icam], coadd_wave))[0]
                    #    J = np.where(np.isin(coadd_wave, data['wave'][icam]))[0]
                    #    coadd_flux3d[J, icam] = data['flux'][icam][I]
                    #    coadd_ivar3d[J, icam] = data['ivar'][icam][I]
                    #    coadd_linemask3d[J, icam] = data['linemask'][icam][I]
                    #
                    #coadd_ivar = np.sum(coadd_ivar3d, axis=1)
                    #coadd_flux = np.zeros_like(coadd_ivar)
                    #good = np.where(coadd_ivar > 0)[0]
                    #coadd_flux[good] = np.sum(coadd_ivar3d[good, :] * coadd_flux3d[good, :], axis=1) / coadd_ivar[good]
                    #coadd_linemask = np.all(coadd_linemask3d, axis=1)
                    coadd_wave = coadd_spec.wave[coadd_bands]
                    coadd_flux = coadd_spec.flux[coadd_bands][igal, :]
                    coadd_ivar = coadd_spec.ivar[coadd_bands][igal, :]
                    coadd_res = Resolution(coadd_spec.resolution_data[coadd_bands][igal, :])
                    coadd_linemask_dict = CFit.build_linemask(coadd_wave, coadd_flux, coadd_ivar, redshift=data['zredrock'])

                    # Map the pixels belonging to individual emission lines and
                    # their local continuum back onto the original per-camera
                    # spectra. These lists of arrays are used in
                    # continuum.ContinnuumTools.smooth_residuals.
                    for icam in np.arange(len(data['cameras'])):
                        data['linemask'].append(np.interp(data['wave'][icam], coadd_wave, coadd_linemask_dict['linemask']*1) > 0)
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
                                     'coadd_linemask': coadd_linemask_dict['linemask']})
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

        self.redrock = Table(np.hstack(self.redrock))
        self.meta = Table(np.hstack(self.meta))
        self.ntargets = len(self.redrock)
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

        nobj = len(self.redrock)

        # The information stored in the metadata table depends on which spectra
        # were fitted (exposures, nightly coadds, deep coadds).
        fluxcols = ['PHOTSYS',
                    'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 
                    'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',
                    'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
                    'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
                    'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2']
        colunit = {'RA': u.deg, 'DEC': u.deg,
                   'FIBERFLUX_G': u.nanomaggy, 'FIBERFLUX_R': u.nanomaggy, 'FIBERFLUX_Z': u.nanomaggy,
                   'FIBERTOTFLUX_G': u.nanomaggy, 'FIBERTOTFLUX_R': u.nanomaggy, 'FIBERTOTFLUX_Z': u.nanomaggy,
                   'FLUX_G': u.nanomaggy, 'FLUX_R': u.nanomaggy,
                   'FLUX_Z': u.nanomaggy, 'FLUX_W1': u.nanomaggy, 'FLUX_W2': u.nanomaggy, 
                   'FLUX_IVAR_G': 1/u.nanomaggy**2, 'FLUX_IVAR_R': 1/u.nanomaggy**2,
                   'FLUX_IVAR_Z': 1/u.nanomaggy**2, 'FLUX_IVAR_W1': 1/u.nanomaggy**2,
                   'FLUX_IVAR_W2': 1/u.nanomaggy**2,
                   }

        skipcols = ['COADD_FIBERSTATUS', 'OBJTYPE', 'TARGET_RA', 'TARGET_DEC'] + fluxcols
        zcols = ['Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE']
        
        meta = Table()
        metacols = self.meta.colnames
        redrockcols = self.redrock.colnames

        # All of this business is so we can get the columns in the order we want
        # (i.e., the order that matches the data model).
        for metacol in ['TARGETID', 'RA', 'DEC', 'TILEID', 'FIBER', 'NIGHT', 'THRUNIGHT']:
            if metacol in metacols:
                meta[metacol] = self.meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]

        for metacol in metacols:
            if metacol in skipcols or metacol in TARGETINGBITCOLS or metacol in meta.colnames:
                continue
            else:
                meta[metacol] = self.meta[metacol]
                if metacol in colunit.keys():
                    meta[metacol].unit = colunit[metacol]

        for bitcol in TARGETINGBITCOLS:
            if bitcol in metacols:
                meta[bitcol] = self.meta[bitcol]
            else:
                meta.add_column(Column(name=bitcol, dtype=np.int64, length=nobj))

        for zcol in zcols:
            meta[zcol] = self.redrock[zcol]
            if zcol in colunit.keys():
                meta[zcol].unit = colunit[zcol]

        for fluxcol in fluxcols:
            meta[fluxcol] = self.meta[fluxcol]
            if fluxcol in colunit.keys():
                meta[fluxcol].unit = colunit[fluxcol]

        out = Table()
        out['TARGETID'] = self.meta['TARGETID']
        if fastphot:
            out = hstack((out, CFit.init_phot_output(nobj)))
        else:
            out = hstack((out, CFit.init_spec_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

        return out, meta

def read_fastspecfit(fastfitfile, fastphot=False, rows=None, columns=None):
    """Read the fitting results.

    """
    if os.path.isfile(fastfitfile):
        if fastphot:
            ext = 'FASTPHOT'
        else:
            ext = 'FASTSPEC'
            
        hdr = fitsio.read_header(fastfitfile, ext='METADATA')
        specprod, coadd_type = hdr['SPECPROD'], hdr['COADDTYP']

        fastfit = Table(fitsio.read(fastfitfile, ext=ext, rows=rows, columns=columns))
        meta = Table(fitsio.read(fastfitfile, ext='METADATA', rows=rows, columns=columns))
        log.info('Read {} object(s) from {} and specprod={}'.format(len(fastfit), fastfitfile, specprod))

        return fastfit, meta, specprod, coadd_type
    
    else:
        log.warning('File {} not found.'.format(fastfitfile))
        return None, None, None

def write_fastspecfit(out, meta, outfile=None, specprod=None,
                      coadd_type=None, survey=None, program=None,
                      fastphot=False):
    """Write out.

    """
    t0 = time.time()
    outdir = os.path.dirname(os.path.abspath(outfile))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    log.info('Writing results for {} objects to {}'.format(len(out), outfile))

    # astropy is really slow!
    #from astropy.io import fits
    #hduprim = fits.PrimaryHDU()
    #
    #hduout = fits.convenience.table_to_hdu(out)
    #if fastphot:
    #    hduout.header['EXTNAME'] = 'FASTPHOT'
    #else:
    #    hduout.header['EXTNAME'] = 'FASTSPEC'
    #
    #hdumeta = fits.convenience.table_to_hdu(meta)
    #hdumeta.header['EXTNAME'] = 'METADATA'
    #
    #if specprod:
    #    hdumeta.header['SPECPROD'] = (specprod, 'spectroscopic production name')
    #if coadd_type:
    #    hdumeta.header['COADDTYP'] = (coadd_type, 'spectral coadd fitted')
    #
    #hx = fits.HDUList([hduprim, hduout, hdumeta])
    #hx.writeto(outfile, overwrite=True, checksum=True)

    hdr = []
    if specprod:
        hdr.append({'name': 'SPECPROD', 'value': specprod, 'comment': 'spectroscopic production name'})
    if coadd_type:
        hdr.append({'name': 'COADDTYP', 'value': coadd_type, 'comment': 'spectral coadd fitted'})
    if survey:
        hdr.append({'name': 'SURVEY', 'value': survey, 'comment': 'survey name'})
    if program:
        hdr.append({'name': 'FAPRGRM', 'value': program, 'comment': 'program name'})

    fitsio.write(outfile, out.as_array(), header=hdr, clobber=True)
    fitsio.write(outfile, meta.as_array(), header=hdr)

    # update the extension name
    if fastphot:
        extname = 'FASTPHOT'
    else:
        extname = 'FASTSPEC'

    with fitsio.FITS(outfile, 'rw') as fits:
        fits[1].write_key('EXTNAME', extname)
        fits[2].write_key('EXTNAME', 'METADATA')
        
    log.info('Writing out took {:.2f} sec'.format(time.time()-t0))

def select(fastfit, metadata, coadd_type, hpxpixels=None, tiles=None, nights=None):
    """Optionally trim to a particular healpix or tile and/or night."""
    keep = np.ones(len(fastfit), bool)
    if coadd_type == 'healpix':
        if hpxpixels:
            pixelkeep = np.zeros(len(fastfit), bool)
            for hpxpixel in hpxpixels:
                pixelkeep = np.logical_or(pixelkeep, metadata['HPXPIXEL'].astype(str) == hpxpixel)
            keep = np.logical_and(keep, pixelkeep)
            log.info('Keeping {} objects from healpixels(s) {}'.format(len(fastfit), ','.join(hpxpixels)))
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
