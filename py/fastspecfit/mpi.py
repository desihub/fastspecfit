"""
fastspecfit.mpi
===============

MPI tools.

"""
import pdb # for debuggin

import os, time
import numpy as np
from glob import glob
import multiprocessing
import fitsio
from astropy.table import Table

from fastspecfit.io import get_qa_filename

from desiutil.log import get_logger
log = get_logger()

def _get_ntargets_one(args):
    return get_ntargets_one(*args)

def get_ntargets_one(specfile, htmldir_root, outdir_root, coadd_type='healpix',
                     makeqa=False, overwrite=False, fastphot=False):
    if makeqa:
        if overwrite:
            ntargets = fitsio.FITS(specfile)[1].get_nrows()
        else:
            outdir = os.path.dirname(specfile).replace(outdir_root, htmldir_root)
            meta = fitsio.read(specfile, 'METADATA', columns=['SURVEY', 'PROGRAM', 'TARGETID', 'HEALPIX'])
            ntargets = 0
            for meta1 in meta:
                pngfile = get_qa_filename(meta1, coadd_type, outdir=outdir, fastphot=fastphot)
                #print(pngfile)
                if not os.path.isfile(pngfile):
                    ntargets += 1
    else:
        from fastspecfit.io import ZWarningMask
        zb = fitsio.read(specfile, 'REDSHIFTS', columns=['Z', 'ZWARN'])
        fm = fitsio.read(specfile, 'FIBERMAP', columns=['TARGETID', 'OBJTYPE'])
        J = ((zb['Z'] > 0.001) * (fm['OBJTYPE'] == 'TGT') * (zb['ZWARN'] & ZWarningMask.NODATA == 0))
        ntargets = np.sum(J)
    return ntargets

def plan(comm=None, specprod=None, specprod_dir=None, coadd_type='healpix',
         survey=None, program=None, healpix=None, tile=None, night=None,
         sample=None, outdir_data='.', mp=1, merge=False, makeqa=False,
         fastphot=False, overwrite=False):

    import fitsio
    from astropy.table import Table, vstack
    from desispec.parallel import weighted_partition
    from fastspecfit.io import DESI_ROOT_NERSC

    t0 = time.time()
    if comm is None:
        rank, size = 0, 1
    else:
        rank, size = comm.rank, comm.size

    if fastphot:
        outprefix = 'fastphot'
        gzip = False
    else:
        outprefix = 'fastspec'
        gzip = True

    desi_root = os.environ.get('DESI_ROOT', DESI_ROOT_NERSC)
    # look for data in the standard location

    if coadd_type == 'healpix':
        subdir = 'healpix'
        if healpix is not None:
            healpixels = healpix.split(',')
    else:
        subdir = 'tiles'

        # figure out which tiles belong to the SV programs
        if tile is None:
            tilefile = os.path.join(desi_root, 'spectro', 'redux', specprod, f'tiles-{specprod}.csv')
            alltileinfo = Table.read(tilefile)
            tileinfo = alltileinfo[['sv' in survey for survey in alltileinfo['SURVEY']]]

            log.info('Retrieved a list of {} {} tiles from {}'.format(
                len(tileinfo), ','.join(sorted(set(tileinfo['SURVEY']))), tilefile))

            tile = np.array(list(set(tileinfo['TILEID'])))

            #if True:
            #    tileinfo = tileinfo[['lrg' in program or 'elg' in program for program in tileinfo['PROGRAM']]]
            #    tile = np.array(list(set(tileinfo['TILEID'])))
            #print(tileinfo)

    if specprod_dir is None:
        specprod_dir = os.path.join(desi_root, 'spectro', 'redux', specprod, subdir)

    outdir = os.path.join(outdir_data, specprod, subdir)
    htmldir = os.path.join(outdir_data, specprod, 'html', subdir)

    def _findfiles(filedir, prefix='redrock', survey=None, program=None, healpix=None,
                   tile=None, night=None, gzip=False, sample=None):
        if gzip:
            fitssuffix = 'fits.gz'
        else:
            fitssuffix = 'fits'
            
        if sample is not None: # special case of an input catalog
            thesefiles, ntargets = [], []
            for onesurvey in sorted(set(sample['SURVEY'].data)):
                S = np.where(onesurvey == sample['SURVEY'])[0]
                for oneprogram in sorted(set(sample['PROGRAM'][S].data)):
                    log.info(f'Building file list for survey={onesurvey} and program={oneprogram}')
                    P = np.where(oneprogram == sample['PROGRAM'][S])[0]
                    uhealpix, _ntargets = np.unique(sample['HEALPIX'][S][P].astype(str), return_counts=True)
                    ntargets.append(_ntargets)
                    for onepix in uhealpix:
                        #ntargets.append(np.sum(onepix == sample['HEALPIX'][S][P].astype(str)))
                        thesefiles.append(os.path.join(filedir, onesurvey, oneprogram, str(int(onepix)//100), onepix,
                                                       f'{prefix}-{onesurvey}-{oneprogram}-{onepix}.{fitssuffix}'))
            if len(thesefiles) > 0:
                #thesefiles = np.array(sorted(np.unique(np.hstack(thesefiles))))
                thesefiles = np.hstack(thesefiles)
                ntargets = np.hstack(ntargets)
            return thesefiles, ntargets
        elif coadd_type == 'healpix':
            thesefiles = []
            for onesurvey in np.atleast_1d(survey):
                for oneprogram in np.atleast_1d(program):
                    log.info(f'Building file list for survey={onesurvey} and program={oneprogram}')
                    if healpix is not None:
                        for onepix in healpixels:
                            _thesefiles = glob(os.path.join(filedir, onesurvey, oneprogram, str(int(onepix)//100), onepix,
                                                            f'{prefix}-{onesurvey}-{oneprogram}-{onepix}.{fitssuffix}'))
                            thesefiles.append(_thesefiles)
                    else:
                        allpix = glob(os.path.join(filedir, onesurvey, oneprogram, '*'))
                        for onepix in allpix:
                            _thesefiles = glob(os.path.join(onepix, '*', f'{prefix}-{onesurvey}-{oneprogram}-*.{fitssuffix}'))
                            thesefiles.append(_thesefiles)
            if len(thesefiles) > 0:
                thesefiles = np.array(sorted(np.unique(np.hstack(thesefiles))))
        elif coadd_type == 'cumulative':
            # Scrape the disk to get the tiles, but since we read the csv file I don't think this ever happens.
            if tile is None:
                tiledirs = np.array(sorted(set(glob(os.path.join(filedir, 'cumulative', '?????')))))
                if len(tiledirs) > 0:
                    tile = [int(os.path.basename(tiledir)) for tiledir in tiledirs]
            if tile is not None:
                thesefiles = []
                for onetile in tile:
                    nightdirs = np.array(sorted(set(glob(os.path.join(filedir, 'cumulative', str(onetile), '????????')))))
                    if len(nightdirs) > 0:
                        # for a given tile, take just the most recent night
                        thisnightdir = nightdirs[-1]
                        thesefiles.append(glob(os.path.join(thisnightdir, f'{prefix}-[0-9]-{onetile}-thru????????.{fitssuffix}')))
                if len(thesefiles) > 0:
                    thesefiles = np.array(sorted(set(np.hstack(thesefiles))))
        elif coadd_type == 'pernight':
            if tile is not None and night is not None:
                thesefiles = []
                for onetile in tile:
                    for onenight in night:
                        thesefiles.append(glob(os.path.join(
                            filedir, 'pernight', str(onetile), str(onenight), f'{prefix}-[0-9]-{onetile}-{onenight}.{fitssuffix}')))
                if len(thesefiles) > 0:
                    thesefiles = np.array(sorted(set(np.hstack(thesefiles))))
            elif tile is not None and night is None:
                thesefiles = np.array(sorted(set(np.hstack([glob(os.path.join(
                    filedir, 'pernight', str(onetile), '????????', f'{prefix}-[0-9]-{onetile}-????????.{fitssuffix}')) for onetile in tile]))))
            elif tile is None and night is not None:
                thesefiles = np.array(sorted(set(np.hstack([glob(os.path.join(
                    filedir, 'pernight', '?????', str(onenight), f'{prefix}-[0-9]-?????-{onenight}.{fitssuffix}')) for onenight in night]))))
            else:
                thesefiles = np.array(sorted(set(glob(os.path.join(
                    filedir, '?????', '????????', f'{prefix}-[0-9]-?????-????????.{fitssuffix}')))))
        elif coadd_type == 'perexp':
            if tile is not None:
                thesefiles = np.array(sorted(set(np.hstack([glob(os.path.join(
                    filedir, 'perexp', str(onetile), '????????', f'{prefix}-[0-9]-{onetile}-exp????????.{fitssuffix}')) for onetile in tile]))))
            else:
                thesefiles = np.array(sorted(set(glob(os.path.join(
                    filedir, 'perexp', '?????', '????????', f'{prefix}-[0-9]-?????-exp????????.{fitssuffix}')))))
        else:
            pass
        return thesefiles

    if merge:
        redrockfiles = None
        outfiles = _findfiles(outdir, prefix=outprefix, survey=survey, program=program, healpix=healpix, tile=tile, night=night, gzip=gzip)
        log.info(f'Found {len(outfiles)} {outprefix} files to be merged.')
    elif makeqa:
        redrockfiles = None
        outfiles = _findfiles(outdir, prefix=outprefix, survey=survey, program=program, healpix=healpix, tile=tile, night=night, gzip=gzip)
        log.info(f'Found {len(outfiles)} {outprefix} files for QA.')
        ntargs = [(outfile, htmldir, outdir, coadd_type, True, overwrite, fastphot) for outfile in outfiles]
    else:
        if sample is not None: # special case of an input catalog
            redrockfiles, ntargets = _findfiles(specprod_dir, prefix='redrock', sample=sample)
        else:
            redrockfiles = _findfiles(specprod_dir, prefix='redrock', survey=survey, program=program,
                                      healpix=healpix, tile=tile, night=night)
        nfile = len(redrockfiles)
        outfiles = []
        for redrockfile in redrockfiles:
            outfile = redrockfile.replace(specprod_dir, outdir).replace('redrock-', f'{outprefix}-')
            if gzip:
                outfile = outfile.replace('.fits', '.fits.gz')
            outfiles.append(outfile) 
        outfiles = np.array(outfiles)
        
        todo = np.ones(len(redrockfiles), bool)
        for ii, outfile in enumerate(outfiles):
            if os.path.isfile(outfile) and not overwrite:
                todo[ii] = False
        if np.count_nonzero(todo) > 0:
            redrockfiles = redrockfiles[todo]
            outfiles = outfiles[todo]
            if sample is not None:
                ntargets = ntargets[todo]
        else: # all done!
            redrockfiles = []
            outfiles = []
            if sample is not None:
                ntargets = np.array([])
                
        log.info(f'Found {len(redrockfiles)}/{nfile} redrockfiles (left) to do.')

        # we already counted ntargets
        if sample is None:
            ntargs = [(redrockfile, None, None, None, False, False, False) for redrockfile in redrockfiles]

    # create groups weighted by the number of targets
    if merge:
        groups = [np.arange(len(outfiles))]
        ntargets = None
    else:
        # we already counted ntargets
        if sample is None:
            if mp > 1:
                with multiprocessing.Pool(mp) as P:
                    ntargets = P.map(_get_ntargets_one, ntargs)
            else:
                ntargets = [get_ntargets_one(*_ntargs) for _ntargs in ntargs]
            ntargets = np.array(ntargets)

        iempty = np.where(ntargets == 0)[0]
        if len(iempty) > 0:
            log.info(f'Skipping {len(iempty)} files with no targets.')

        itodo = np.where(ntargets > 0)[0]
        if len(itodo) > 0:
            ntargets = ntargets[itodo]
            log.info(f'Number of targets left to do: {np.sum(ntargets):,d}.')
            
            if redrockfiles is not None:
                redrockfiles = redrockfiles[itodo]
            if outfiles is not None:
                outfiles = outfiles[itodo]
                
            groups = weighted_partition(ntargets, size)

            ## Assign the sample to ranks to make the ntargets distribution per rank ~flat.
            ## https://stackoverflow.com/questions/33555496/split-array-into-equally-weighted-chunks-based-on-order
            #indices = np.arange(len(ntargets))
            #cumuweight = ntargets.cumsum() / ntargets.sum()
            #idx = np.searchsorted(cumuweight, np.linspace(0, 1, size, endpoint=False)[1:])
            #if len(idx) < size: # can happen in corner cases or with 1 rank
            #    groups = np.array_split(indices, size) # unweighted
            #else:
            #    groups = np.array_split(indices, idx) # weighted
            #for ii in range(size): # sort by weight
            #    srt = np.argsort(ntargets[groups[ii]])
            #    groups[ii] = groups[ii][srt]
        else:
            if redrockfiles is not None:
                redrockfiles = []
            if outfiles is not None:
                outfiles = []
            groups = [np.array([])]

    if merge:
        if len(outfiles) == 0:
            if rank == 0:
                log.debug(f'No {outprefix} files in {outdir} found!')
            return '', list(), list(), list(), None
        return outdir, redrockfiles, outfiles, None, None
    elif makeqa:
        if len(outfiles) == 0:
            if rank == 0:
                log.debug(f'No {outprefix} files in {outdir} left to do!')
            return '', list(), list(), list(), None
        #  hack--build the output directories and pass them in the 'redrockfiles'
        #  position! for coadd_type==cumulative, strip out the 'lastnight' argument
        if coadd_type == 'cumulative':
            redrockfiles = np.array([os.path.dirname(os.path.dirname(outfile)).replace(outdir, htmldir) for outfile in outfiles])
        else:
            redrockfiles = np.array([os.path.dirname(outfile).replace(outdir, htmldir) for outfile in outfiles])
    else:
        if len(redrockfiles) == 0:
            if rank == 0:
                log.info('All files have been processed!')
            return '', list(), list(), list(), None

    return outdir, redrockfiles, outfiles, groups, ntargets

def _read_to_merge_one(args):
    return read_to_merge_one(*args)

def read_to_merge_one(filename, extname):
    info = fitsio.FITS(filename)
    ext = [_info.get_extname() for _info in info]
    if extname not in ext:
        log.warning(f'Missing extension {extname} in file {filename}')
        return [], []
    if 'METADATA' not in ext:
        log.warning(f'Missing extension METADATA in file {filename}')
        return [], []
    out = Table(info[extname].read())
    meta = Table(info['METADATA'].read())
    return out, meta

def _domerge(outfiles, extname='FASTSPEC', survey=None, program=None,
             outprefix=None, specprod=None, coadd_type=None, mergefile=None,
             fastphot=False, mp=1):
    """Support routine to merge a set of input files.

    """
    from astropy.table import vstack
    from desiutil.depend import getdep, hasdep
    from fastspecfit.io import write_fastspecfit
    
    t0 = time.time()
    out, meta = [], []

    mpargs = [[outfile, extname] for outfile in outfiles]
    if mp > 1:
        with multiprocessing.Pool(mp) as P:
            _out = P.map(_read_to_merge_one, mpargs)
    else:
        _out = [read_to_merge_one(*mparg) for mparg in mpargs]
    _out = list(zip(*_out))
    out = vstack(_out[0])
    meta = vstack(_out[1])
    del _out
        
    # sort?
    srt = np.argsort(meta['TARGETID'])
    out = out[srt]
    meta = meta[srt]

    if outprefix:
        log.info(f'Merging {len(out):,d} objects from {len(outfiles)} {outprefix} files took {(time.time()-t0)/60.0:.2f} min.')
    else:
        log.info(f'Merging {len(out):,d} objects from {len(outfiles)} files took {(time.time()-t0)/60.0:.2f} min.')

    hdr = fitsio.read_header(outfiles[0])

    # sensible defaults
    deps = {}
    deps['INPUTZ'] = False
    deps['NOSCORR'] = False
    deps['NOPHOTO'] = False
    deps['BRDLFIT'] = True
    deps['CONSAGE'] = False
    deps['USEQNET'] = True
    for key in deps.keys():
        if key in hdr:
            deps[key] = hdr[key]
            
    deps2 = {}
    deps2['FPHOTO_FILE'] = None
    deps2['FTEMPLATES_FILE'] = None
    deps2['EMLINES_FILE'] = None
    for key in deps2.keys():
        if hasdep(hdr, key):
            deps2[key] = getdep(hdr, key)

    write_fastspecfit(out, meta, outfile=mergefile, specprod=specprod,
                      coadd_type=coadd_type, fastphot=fastphot,
                      fphotofile=deps2['FPHOTO_FILE'], templates=deps2['FTEMPLATES_FILE'],
                      emlinesfile=deps2['EMLINES_FILE'], inputz=deps['INPUTZ'],
                      ignore_photometry=deps['NOPHOTO'], broadlinefit=deps['BRDLFIT'],
                      constrain_age=deps['CONSAGE'], use_quasarnet=deps['USEQNET'],
                      no_smooth_continuum=deps['NOSCORR'])

def merge_fastspecfit(specprod=None, coadd_type=None, survey=None, program=None,
                      healpix=None, tile=None, night=None, outsuffix=None,
                      fastphot=False, specprod_dir=None, outdir_data='.',
                      fastfiles_to_merge=None, 
                      mergedir=None, supermerge=False, overwrite=False, mp=1):
    """Merge all the individual catalogs into a single large catalog. Runs only on
    rank 0.

    supermerge - merge previously merged catalogs

    """
    import fitsio
    from copy import copy
    from astropy.io import fits
    from astropy.table import Table, vstack
    from fastspecfit.mpi import plan
    
    if fastphot:
        outprefix = 'fastphot'
        extname = 'FASTPHOT'
    else:
        outprefix = 'fastspec'
        extname = 'FASTSPEC'

    if mergedir is None:
        mergedir = os.path.join(outdir_data, specprod, 'catalogs')
        if not os.path.isdir(mergedir):
            os.makedirs(mergedir, exist_ok=True)

    if outsuffix is None:
        outsuffix = specprod

    # merge previously merged catalogs into one big catalog (and then return)
    if supermerge:
        if fastfiles_to_merge is None:
            _outfiles = os.path.join(mergedir, f'{outprefix}-{outsuffix}-*.fits*')
            outfiles = glob(_outfiles)
        else:
            _outfiles = fastfiles_to_merge
            outfiles = _outfiles
        if len(outfiles) > 0:
            log.info(f'Merging {len(outfiles):,d} catalogs')
            mergefile = os.path.join(mergedir, f'{outprefix}-{outsuffix}.fits')
            _domerge(outfiles, extname=extname, mergefile=mergefile, outprefix=outprefix,
                     specprod=specprod, coadd_type=coadd_type, fastphot=fastphot, mp=mp)
        else:
            log.info(f'No catalogs found: {_outfiles}')
        return

    if coadd_type == 'healpix':
        if survey is None or program is None:
            log.warning(f'coadd_type={coadd_type} requires survey and program inputs.')
            return
        surveys = copy(survey)
        programs = copy(program)
        for survey in surveys:
            for program in programs:
                mergefile = os.path.join(mergedir, f'{outprefix}-{specprod}-{survey}-{program}.fits')
                if os.path.isfile(mergefile) and not overwrite:
                    log.info(f'Merged output file {mergefile} exists!')
                    continue
                #survey = np.atleast_1d(survey)
                #program = np.atleast_1d(program)
                _, _, outfiles, _, _ = plan(specprod=specprod, survey=survey, program=program, healpix=healpix,
                                            merge=True, fastphot=fastphot, specprod_dir=specprod_dir,
                                            outdir_data=outdir_data, overwrite=overwrite)
                if len(outfiles) > 0:
                    _domerge(outfiles, extname=extname, survey=survey[0], program=program[0],
                             mergefile=mergefile, outprefix=outprefix, specprod=specprod,
                             coadd_type=coadd_type, fastphot=fastphot, mp=mp)
    else:
        mergefile = os.path.join(mergedir, f'{outprefix}-{specprod}-{coadd_type}.fits')
        if os.path.isfile(mergefile) and not overwrite:
            log.info(f'Merged output file {mergefile} exists!')
            return
        _, _, outfiles, _, _ = plan(specprod=specprod, coadd_type=coadd_type, tile=tile, night=night,
                                    merge=True, fastphot=fastphot, specprod_dir=specprod_dir,
                                    outdir_data=outdir_data, overwrite=overwrite)
        if len(outfiles) > 0:
            _domerge(outfiles, extname=extname, mergefile=mergefile, outprefix=outprefix,
                     specprod=specprod, coadd_type=coadd_type, fastphot=fastphot, mp=mp)
