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

from desiutil.log import get_logger
log = get_logger()

def _get_ntargets_one(args):
    return get_ntargets_one(*args)

def get_ntargets_one(specfile, makeqa=False):
    if makeqa:
        ntargets = fitsio.FITS(specfile)[1].get_nrows() # fragile?
    else:
        zb = fitsio.read(specfile, 'REDSHIFTS', columns=['Z'])#, 'ZWARN'])
        fm = fitsio.read(specfile, 'FIBERMAP', columns=['TARGETID', 'OBJTYPE'])#, 'COADD_FIBERSTATUS', 'PHOTSYS'])
        J = ((zb['Z'] > 0.001) *
             #(zb['ZWARN'] <= 4) *
             #(zb['SPECTYPE'] == 'GALAXY') *
             #(fm['PHOTSYS'] != 'G') *
             #(fm['COADD_FIBERSTATUS'] == 0) *
             (fm['OBJTYPE'] == 'TGT') *
             (fm['TARGETID'] > 0)
             )
        ntargets = np.sum(J)
    return ntargets

def weighted_partition(weights, n):
    """
    Partition ``weights`` into ``n`` groups with approximately same sum(weights).

    Args:
        weights: array-like weights
        n: number of groups

    Returns (groups, groupweights):
        groups: list of lists of indices of weights for each group
        groupweights: sum of weights assigned to each group

    
    """
    sumweights = np.zeros(n, dtype=float)
    groups = list()
    for i in range(n):
        groups.append(list())
    weights = np.asarray(weights)
    for i in np.argsort(-weights):
        j = np.argmin(sumweights)
        groups[j].append(i)
        sumweights[j] += weights[i]

    return groups, np.array([np.sum(x) for x in sumweights])

def group_redrockfiles(specfiles, maxnodes=256, comm=None, makeqa=False):
    '''
    Group redrockfiles to balance runtimes

    Args:
        specfiles: list of spectra filepaths

    Options:
        maxnodes: split the spectra into this number of nodes
        comm: MPI communicator

    Returns (groups, ntargets, grouptimes):
      * groups: list of lists of indices to specfiles
      * list of number of targets per group
      * grouptimes: list of expected runtimes for each group

    '''
    import fitsio
    
    if comm is None:
        rank, size = 0, 1
    else:
        rank, size = comm.rank, comm.size

    npix = len(specfiles)
    pixgroups = np.array_split(np.arange(npix), size)
    ntargets = np.zeros(len(pixgroups[rank]), dtype=int)
    for i, j in enumerate(pixgroups[rank]):
        if makeqa:
            ntargets[i] = fitsio.FITS(specfiles[j])[1].get_nrows()
        else:
            #ntargets[i] = fitsio.FITS(specfiles[j])[1].get_nrows()
            zb = fitsio.read(specfile, 'REDSHIFTS', columns=['Z'])#, 'ZWARN'])
            fm = fitsio.read(specfile, 'FIBERMAP', columns=['TARGETID', 'OBJTYPE'])#, 'COADD_FIBERSTATUS', 'PHOTSYS'])
            J = ((zb['Z'] > 0.001) *
                 #(zb['ZWARN'] <= 4) *
                 #(zb['SPECTYPE'] == 'GALAXY') *
                 #(fm['PHOTSYS'] != 'G') *
                 #(fm['COADD_FIBERSTATUS'] == 0) *
                 (fm['OBJTYPE'] == 'TGT') *
                 (fm['TARGETID'] > 0)
                 )
            nt = np.sum(J)
            ntargets[i] = nt
        log.debug(i, j, specfiles[j], ntargets[i])

    if False:
        if comm is not None:
            ntargets = comm.gather(ntargets)
            if rank == 0:
                ntargets = np.concatenate(ntargets)
            ntargets = comm.bcast(ntargets, root=0)

        sumntargets = np.sum(sumntargets)
        runtimes = 30 + 0.4*sumntargets

        # Aim for 25 minutes, but don't exceed maxnodes number of nodes.
        ntime = 25
        if comm is not None:
            numnodes = comm.size
        else:
            numnodes = min(maxnodes, int(np.ceil(np.sum(runtimes)/(ntime*60))))

        groups, grouptimes = weighted_partition(runtimes, numnodes)
        ntargets = np.array([np.sum(ntargets[ii]) for ii in groups])
    else:
        groups = pixgroups
        grouptimes = None

    return groups, ntargets, grouptimes

def backup_logs(logfile):
    '''
    Move logfile -> logfile.0 or logfile.1 or logfile.n as needed

    TODO: make robust against logfile.abc also existing
    '''
    logfiles = glob(logfile+'.*')
    newlog = logfile+'.'+str(len(logfiles))
    assert not os.path.exists(newlog)
    os.rename(logfile, newlog)
    return newlog

def plan(comm=None, specprod=None, specprod_dir=None, coadd_type='healpix',
         survey=None, program=None, healpix=None, tile=None, night=None, 
         outdir_data='.', outdir_html='.', mp=1, merge=False, makeqa=False,
         fastphot=False, overwrite=False):

    import fitsio
    from astropy.table import Table, vstack
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
            tilefile = os.path.join(desi_root, 'spectro', 'redux', specprod, 'tiles-{}.csv'.format(specprod))
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
    #htmldir = os.path.join(outdir_html, specprod, subdir)

    def _findfiles(filedir, prefix='redrock', survey=None, program=None, healpix=None, tile=None, night=None):
        if coadd_type == 'healpix':
            thesefiles = []
            for onesurvey in np.atleast_1d(survey):
                for oneprogram in np.atleast_1d(program):
                    log.info('Building file list for survey={} and program={}'.format(onesurvey, oneprogram))
                    if healpix is not None:
                        for onepix in healpixels:
                            _thesefiles = glob(os.path.join(filedir, onesurvey, oneprogram, str(int(onepix)//100), onepix,
                                                            '{}-{}-{}-*.fits'.format(prefix, onesurvey, oneprogram)))
                            thesefiles.append(_thesefiles)
                    else:
                        allpix = glob(os.path.join(filedir, onesurvey, oneprogram, '*'))
                        for onepix in allpix:
                            _thesefiles = glob(os.path.join(onepix, '*', '{}-{}-{}-*.fits'.format(prefix, onesurvey, oneprogram)))
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
                        thesefiles.append(glob(os.path.join(thisnightdir, '{}-[0-9]-{}-thru????????.fits*'.format(prefix, onetile))))
                thesefiles = np.array(sorted(set(np.hstack(thesefiles))))
        elif coadd_type == 'pernight':
            if tile is not None and night is not None:
                thesefiles = []
                for onetile in tile:
                    for onenight in night:
                        thesefiles.append(glob(os.path.join(
                            filedir, 'pernight', str(onetile), str(onenight), '{}-[0-9]-{}-{}.fits'.format(prefix, onetile, onenight))))
                if len(thesefiles) > 0:
                    thesefiles = np.array(sorted(set(np.hstack(thesefiles))))
            elif tile is not None and night is None:
                thesefiles = np.array(sorted(set(np.hstack([glob(os.path.join(
                    filedir, 'pernight', str(onetile), '????????', '{}-[0-9]-{}-????????.fits'.format(
                    prefix, onetile))) for onetile in tile]))))
            elif tile is None and night is not None:
                thesefiles = np.array(sorted(set(np.hstack([glob(os.path.join(
                    filedir, 'pernight', '?????', str(onenight), '{}-[0-9]-?????-{}.fits'.format(
                    prefix, onenight))) for onenight in night]))))
            else:
                thesefiles = np.array(sorted(set(glob(os.path.join(
                    filedir, '?????', '????????', '{}-[0-9]-?????-????????.fits'.format(prefix))))))
        elif coadd_type == 'perexp':
            if tile is not None:
                thesefiles = np.array(sorted(set(np.hstack([glob(os.path.join(
                    filedir, 'perexp', str(onetile), '????????', '{}-[0-9]-{}-exp????????.fits'.format(
                    prefix, onetile))) for onetile in tile]))))
            else:
                thesefiles = np.array(sorted(set(glob(os.path.join(
                    filedir, 'perexp', '?????', '????????', '{}-[0-9]-?????-exp????????.fits'.format(prefix))))))
        else:
            pass
        return thesefiles

    if merge:
        redrockfiles = None
        outfiles = _findfiles(outdir, prefix=outprefix, survey=survey, program=program, healpix=healpix, tile=tile, night=night)
        log.info('Found {} {} files to be merged.'.format(len(outfiles), outprefix))
    elif makeqa:
        redrockfiles = None
        outfiles = _findfiles(outdir, prefix=outprefix, survey=survey, program=program, healpix=healpix, tile=tile, night=night)
        log.info('Found {} {} files for QA.'.format(len(outfiles), outprefix))
        ntargs = [(outfile, True) for outfile in outfiles]
    else:
        redrockfiles = _findfiles(specprod_dir, prefix='redrock', survey=survey, program=program, healpix=healpix, tile=tile, night=night)
        nfile = len(redrockfiles)
        outfiles = []
        for redrockfile in redrockfiles:
            outfile = redrockfile.replace(specprod_dir, outdir).replace('redrock-', '{}-'.format(outprefix))
            if gzip:
                outfile = outfile.replace('.fits', '.fits.gz')
            outfiles.append(outfile) 
        outfiles = np.array(outfiles)
        
        todo = np.ones(len(redrockfiles), bool)
        for ii, outfile in enumerate(outfiles):
            if os.path.isfile(outfile) and not overwrite:
                todo[ii] = False
        redrockfiles = redrockfiles[todo]
        outfiles = outfiles[todo]
        log.info('Found {}/{} redrockfiles (left) to do.'.format(len(redrockfiles), nfile))
        ntargs = [(redrockfile, False) for redrockfile in redrockfiles]

    # create groups weighted by the number of targets
    if merge:
        groups = [np.arange(len(outfiles))]
        ntargets = None
    else:
        if mp > 1:
            with multiprocessing.Pool(mp) as P:
                ntargets = P.map(_get_ntargets_one, ntargs)
        else:
            ntargets = [get_ntargets_one(*_ntargs) for _ntargs in ntargs]
        ntargets = np.array(ntargets)

        iempty = np.where(ntargets==0)[0]
        if len(iempty) > 0:
            log.info('Skipping {} files with no targets.'.format(len(iempty)))

        itodo = np.where(ntargets > 0)[0]
        if len(itodo) > 0:
            ntargets = ntargets[itodo]
            indices = np.arange(len(ntargets))
            if redrockfiles is not None:
                redrockfiles = redrockfiles[itodo]
            if outfiles is not None:
                outfiles = outfiles[itodo]

            # Assign the sample to ranks to make the ntargets distribution per rank ~flat.
            # https://stackoverflow.com/questions/33555496/split-array-into-equally-weighted-chunks-based-on-order
            cumuweight = ntargets.cumsum() / ntargets.sum()
            idx = np.searchsorted(cumuweight, np.linspace(0, 1, size, endpoint=False)[1:])
            if len(idx) < size: # can happen in corner cases or with 1 rank
                groups = np.array_split(indices, size) # unweighted
            else:
                groups = np.array_split(indices, idx) # weighted
            for ii in range(size): # sort by weight
                srt = np.argsort(ntargets[groups[ii]])
                groups[ii] = groups[ii][srt]
        else:
            groups = [np.array([])]

    #if comm:
    #    outdir = comm.bcast(outdir, root=0)
    #    redrockfiles = comm.bcast(redrockfiles, root=0)
    #    outfiles = comm.bcast(outfiles, root=0)
    #    groups = comm.bcast(groups, root=0)
    #    ntargets = comm.bcast(ntargets, root=0)

    if merge:
        if len(outfiles) == 0:
            if rank == 0:
                log.debug('No {} files in {} found!'.format(outprefix, outdir))
            return '', list(), list(), list(), None
        return outdir, redrockfiles, outfiles, None, None
    elif makeqa:
        if len(outfiles) == 0:
            if rank == 0:
                log.debug('No {} files in {} found!'.format(outprefix, outdir))
            return '', list(), list(), list(), None
        #  hack--build the output directories and pass them in the 'redrockfiles'
        #  position! for coadd_type==cumulative, strip out the 'lastnight' argument
        if coadd_type == 'cumulative':
            #redrockfiles = []
            #for outfile in outfiles:
            #    dd = os.path.split(outfile)
            #    redrockfiles.append(os.path.dirname(dd[0]).replace(outdir, htmldir))
            #    os.path.dirname(dd[0])
            redrockfiles = np.array([os.path.dirname(os.path.dirname(outfile)).replace(outdir, htmldir) for outfile in outfiles])
        else:
            redrockfiles = np.array([os.path.dirname(outfile).replace(outdir, htmldir) for outfile in outfiles])
    else:
        if len(redrockfiles) == 0:
            if rank == 0:
                log.info('All files have been processed!')
            return '', list(), list(), list(), None

    return outdir, redrockfiles, outfiles, groups, ntargets

def merge_fastspecfit(specprod=None, coadd_type=None, survey=None, program=None,
                      healpix=None, tile=None, night=None, outsuffix=None,
                      fastphot=False, specprod_dir=None, outdir_data='.',
                      mergedir=None, supermerge=False, overwrite=False):
    """Merge all the individual catalogs into a single large catalog. Runs only on
    rank 0.

    supermerge - merge previously merged catalogs

    """
    import fitsio
    from copy import copy
    from astropy.io import fits
    from astropy.table import Table, vstack
    from fastspecfit.mpi import plan
    from fastspecfit.io import write_fastspecfit
    
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

    def _domerge(outfiles, extname='FASTSPEC', survey=None, program=None, mergefile=None):
        t0 = time.time()
        out, meta = [], []
        for outfile in outfiles:
            info = fitsio.FITS(outfile)
            ext = [_info.get_extname() for _info in info]
            if extname not in ext:
                log.warning('Missing extension {} in file {}'.format(extname, outfile))
                continue
            if 'METADATA' not in ext:
                log.warning('Missing extension METADATA in file {}'.format(outfile))
                continue
            out.append(Table(info[extname].read()))
            meta.append(Table(info['METADATA'].read()))
        out = vstack(out)
        meta = vstack(meta)

        # sort?
        srt = np.argsort(meta['TARGETID'])
        out = out[srt]
        meta = meta[srt]
        
        log.info('Merging {:,d} objects from {} {} files took {:.2f} min.'.format(
            len(out), len(outfiles), outprefix, (time.time()-t0)/60.0))
        
        write_fastspecfit(out, meta, outfile=mergefile, specprod=specprod,
                          coadd_type=coadd_type, fastphot=fastphot)

    # merge previously merged catalogs into one big catalog (and then return)
    if supermerge:
        _outfiles = os.path.join(mergedir, '{}-{}-*.fits*'.format(outprefix, outsuffix))
        outfiles = glob(_outfiles)
        #print(_outfiles, outfiles)
        if len(outfiles) > 0:
            log.info('Merging {:,d} catalogs'.format(len(outfiles)))
            mergefile = os.path.join(mergedir, '{}-{}.fits'.format(outprefix, outsuffix))
            _domerge(outfiles, extname=extname, mergefile=mergefile)
        else:
            log.info('No catalogs found: {}'.format(_outfiles))
        return

    if coadd_type == 'healpix':
        if survey is None or program is None:
            log.warning('coadd_type={} requires survey and program inputs.'.format(coadd_type))
            return
        surveys = copy(survey)
        programs = copy(program)
        for survey in surveys:
            for program in programs:
                mergefile = os.path.join(mergedir, '{}-{}-{}-{}.fits'.format(outprefix, specprod, survey, program))
                if os.path.isfile(mergefile) and not overwrite:
                    log.info('Merged output file {} exists!'.format(mergefile))
                    continue
                #survey = np.atleast_1d(survey)
                #program = np.atleast_1d(program)
                _, _, outfiles, _, _ = plan(specprod=specprod, survey=survey, program=program, healpix=healpix,
                                            merge=True, fastphot=fastphot, specprod_dir=specprod_dir,
                                            outdir_data=outdir_data, overwrite=overwrite)
                if len(outfiles) > 0:
                    _domerge(outfiles, extname=extname, survey=survey[0], program=program[0], mergefile=mergefile)
    else:
        mergefile = os.path.join(mergedir, '{}-{}-{}.fits'.format(outprefix, specprod, coadd_type))
        if os.path.isfile(mergefile) and not overwrite:
            log.info('Merged output file {} exists!'.format(mergefile))
            return
        _, _, outfiles, _, _ = plan(specprod=specprod, coadd_type=coadd_type, tile=tile, night=night,
                                    merge=True, fastphot=fastphot, specprod_dir=specprod_dir,
                                    outdir_data=outdir_data, overwrite=overwrite)
        if len(outfiles) > 0:
            _domerge(outfiles, extname=extname, mergefile=mergefile)
