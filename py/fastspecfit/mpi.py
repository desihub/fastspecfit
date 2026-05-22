"""
fastspecfit.mpi
===============

MPI tools.

"""
import os, time
import numpy as np
from glob import glob
import multiprocessing
import fitsio
from astropy.table import Table, vstack

from fastspecfit.io import get_qa_filename
from fastspecfit.logger import log
from fastspecfit.util import fsftime


def _get_ntargets_one(args):
    return get_ntargets_one(*args)


def get_ntargets_one(specfile, htmldir_root, outdir_root, coadd_type='healpix',
                     makeqa=False, overwrite=False, fastphot=False):
    """Count the number of targets to process in a single output file."""
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
        from fastspecfit.util import ZWarningMask
        zb = fitsio.read(specfile, 'REDSHIFTS', columns=['Z', 'ZWARN'])
        fm = fitsio.read(specfile, 'FIBERMAP', columns=['TARGETID', 'OBJTYPE'])
        J = ((zb['Z'] > 0.001) * (fm['OBJTYPE'] == 'TGT') * (zb['ZWARN'] & ZWarningMask.NODATA == 0))
        ntargets = np.sum(J)
    return ntargets


def findfiles(filedir, prefix='redrock', coadd_type=None, survey=None,
              program=None, healpix=None, tile=None, night=None,
              gzip=False, sample=None):
    """Find all DESI spectral files matching the given selection criteria.

    Parameters
    ----------
    filedir : :class:`str`
        Root directory to search for files.
    prefix : :class:`str`, optional
        Filename prefix. Default is ``'redrock'``.
    coadd_type : :class:`str` or None, optional
        Coadd type: ``'healpix'``, ``'cumulative'``, ``'pernight'``, or
        ``'perexp'``.
    survey : :class:`str` or array-like of str, optional
        DESI survey name(s), e.g. ``'main'``.
    program : :class:`str` or array-like of str, optional
        DESI program name(s), e.g. ``'dark'``.
    healpix : array-like or None, optional
        Specific HEALPix pixel(s) to include.
    tile : array-like or None, optional
        Specific tile ID(s) to include.
    night : array-like or None, optional
        Specific observing night(s) to include.
    gzip : :class:`bool`, optional
        If ``True``, look for ``.fits.gz`` files. Default is ``False``.
    sample : :class:`astropy.table.Table` or None, optional
        Input sample catalog; when provided, file paths are built directly
        from ``SURVEY``, ``PROGRAM``, and ``HEALPIX`` columns.

    Returns
    -------
    thesefiles : :class:`numpy.ndarray` of str
        Sorted array of matching file paths.

    """
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
                    for onepix in healpix:
                        _thesefiles = os.path.join(filedir, onesurvey, oneprogram, str(int(onepix)//100), onepix,
                                                   f'{prefix}-{onesurvey}-{oneprogram}-{onepix}.{fitssuffix}')
                        thesefiles.append(glob(_thesefiles))
                else:
                    allpix = os.path.join(filedir, onesurvey, oneprogram, '*')
                    for onepix in glob(allpix):
                        _thesefiles = os.path.join(onepix, '*', f'{prefix}-{onesurvey}-{oneprogram}-*.{fitssuffix}')
                        thesefiles.append(glob(_thesefiles))
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


def plan_merge(outdir, outprefix, coadd_type, survey, program, healpix,
               tile, night, sample=None, gzip=False):
    """Build the list of output files to be merged."""
    redrockfiles = None
    if sample is not None: # special case of an input catalog
        outfiles, _ = findfiles(outdir, prefix=outprefix, sample=sample)
    else:
        outfiles = findfiles(outdir, prefix=outprefix, coadd_type=coadd_type,
                             survey=survey, program=program, healpix=healpix,
                             tile=tile, night=night, gzip=gzip)
    log.info(f'Found {len(outfiles)} {outprefix} files to be merged.')
    return redrockfiles, outfiles


def plan_makeqa(outdir, htmldir, outprefix, coadd_type, survey, program,
                healpix, tile, night, sample=None, gzip=False):
    """Build the list of output files and HTML directories for QA generation."""
    outfiles = findfiles(outdir, prefix=outprefix, coadd_type=coadd_type,
                         survey=survey, program=program, healpix=healpix,
                         tile=tile, night=night, gzip=gzip)
    log.info(f'Found {len(outfiles)} {outprefix} files for QA.')

    #  Hack!--build the output directories and pass them in the 'redrockfiles'
    #  position! for coadd_type=='cumulative', strip out the 'lastnight'
    #  argument.
    if len(outfiles) == 0:
        return _, outfiles

    if coadd_type == 'cumulative':
        redrockfiles = np.array([os.path.dirname(os.path.dirname(outfile)).replace(outdir, htmldir)
                                 for outfile in outfiles])
    else:
        redrockfiles = np.array([os.path.dirname(outfile).replace(outdir, htmldir)
                                 for outfile in outfiles])

    return redrockfiles, outfiles


def plan(comm=None, specprod=None, specprod_dir=None, coadd_type='healpix',
         survey=None, program=None, healpix=None, tile=None, night=None,
         sample=None, outdir_data='.', mp=1, merge=False, makeqa=False,
         fastphot=False, overwrite=False):
    """Determine which files still need to be processed.

    Parameters
    ----------
    comm : MPI communicator or None, optional
        MPI communicator for distributed execution. When ``None``, runs
        single-process.
    specprod : :class:`str` or None, optional
        DESI spectroscopic production name (e.g. ``'iron'``).
    specprod_dir : :class:`str` or None, optional
        Override the standard specprod directory path.
    coadd_type : :class:`str`, optional
        Coadd type: ``'healpix'``, ``'cumulative'``, ``'pernight'``, or
        ``'perexp'``. Default is ``'healpix'``.
    survey : :class:`str` or array-like of str, optional
        DESI survey name(s).
    program : :class:`str` or array-like of str, optional
        DESI program name(s).
    healpix : array-like or None, optional
        Specific HEALPix pixel(s) to process.
    tile : array-like or None, optional
        Specific tile ID(s) to process.
    night : array-like or None, optional
        Specific observing night(s) to process.
    sample : :class:`astropy.table.Table` or None, optional
        Input sample catalog; overrides file-system discovery.
    outdir_data : :class:`str`, optional
        Root output directory. Default is ``'.'``.
    mp : :class:`int`, optional
        Number of multiprocessing workers. Default is 1.
    merge : :class:`bool`, optional
        If ``True``, plan for catalog merging rather than fitting. Default is
        ``False``.
    makeqa : :class:`bool`, optional
        If ``True``, plan for QA generation rather than fitting. Default is
        ``False``.
    fastphot : :class:`bool`, optional
        If ``True``, use ``fastphot`` output filenames. Default is ``False``.
    overwrite : :class:`bool`, optional
        If ``True``, include files that already have output. Default is
        ``False``.

    Returns
    -------
    outdir : :class:`str`
        Root output directory.
    redrockfiles : :class:`list` or :class:`numpy.ndarray`
        Input redrock files to process (or ``None``/``[]`` when merging/QA).
    outfiles : :class:`list` or :class:`numpy.ndarray`
        Corresponding output file paths.
    ntargets : :class:`numpy.ndarray` or None
        Number of targets per file; ``None`` when merging.

    """
    if comm:
        rank, size = comm.rank, comm.size
    else:
        rank, size = 0, 1

    if rank > 0:
        outdir = None
        redrockfiles = None
        outfiles = None
        ntargets = None
    else:
        t0 = time.time()
        if fastphot:
            outprefix = 'fastphot'
            gzip = False
        else:
            outprefix = 'fastspec'
            gzip = True

        redux_dir = os.path.expandvars(os.environ.get('DESI_SPECTRO_REDUX'))

        # look for data in the standard location
        if coadd_type == 'healpix':
            subdir = 'healpix'
        else:
            subdir = 'tiles'

            # figure out which tiles belong to the SV programs
            if tile is None:
                tilefile = os.path.join(redux_dir, specprod, f'tiles-{specprod}.csv')
                alltileinfo = Table.read(tilefile)
                tileinfo = alltileinfo[['sv' in survey for survey in alltileinfo['SURVEY']]]

                log.info('Retrieved a list of {} {} tiles from {}'.format(
                    len(tileinfo), ','.join(sorted(set(tileinfo['SURVEY']))), tilefile))

                tile = np.array(list(set(tileinfo['TILEID'])))

        if specprod_dir is None:
            specprod_dir = os.path.join(redux_dir, specprod, subdir)

        outdir = os.path.join(outdir_data, specprod, subdir)
        htmldir = os.path.join(outdir_data, specprod, 'html', subdir)

        if merge:
            redrockfiles, outfiles = plan_merge(
                outdir, outprefix, coadd_type, survey, program,
                healpix, tile, night, sample=sample, gzip=gzip)
            if len(outfiles) == 0:
                log.debug(f'No {outprefix} files in {outdir} found!')
                return '', list(), list(), None
            return outdir, redrockfiles, outfiles, None
        elif makeqa:
            redrockfiles, outfiles = plan_makeqa(
                outdir, htmldir, outprefix, coadd_type, survey, program,
                healpix, tile, night, sample=sample, gzip=gzip)
            if len(outfiles) == 0:
                log.debug(f'No {outprefix} files in {outdir} left to do!')
                return '', list(), list(), None
            ntargs = [(outfile, htmldir, outdir, coadd_type, True, overwrite, fastphot)
                      for outfile in outfiles]
        else:
            if sample is not None: # special case of an input catalog
                redrockfiles, ntargets = findfiles(specprod_dir, prefix='redrock', sample=sample)
            else:
                redrockfiles = findfiles(specprod_dir, prefix='redrock', coadd_type=coadd_type,
                                         survey=survey, program=program,
                                         healpix=healpix, tile=tile, night=night)

            # In principle, we could parallelize this piece of code...
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

            log.info(f'Found {len(redrockfiles):,d}/{nfile:,d} redrockfiles left.')

            # we already counted ntargets
            if sample is None:
                ntargs = [(redrockfile, None, None, None, False, False, False)
                          for redrockfile in redrockfiles]

    # Count the number of targets in each file. Distribute work to the other
    # ranks, either via MPI or multiprocessing. Note that if `sample` exists,
    # we already counted targets.
    if sample is None:
        if comm:
            if rank == 0:
                ntargs_byrank = np.array_split(ntargs, size)
                for onerank in range(1, size):
                    comm.send(ntargs_byrank[onerank], dest=onerank)
                ntargs_onerank = ntargs_byrank[rank]
            else:
                ntargs_onerank = comm.recv(source=0)

            ntargets = []
            for ntarg_onerank in ntargs_onerank:
                ntargets.append(get_ntargets_one(*ntarg_onerank))

            if rank > 0:
                comm.send(ntargets, dest=0)
            else:
                for onerank in range(1, size):
                    ntargets.extend(comm.recv(source=onerank))
                if len(ntargets) > 0:
                    ntargets = np.hstack(ntargets)
        else:
            if mp > 1:
                with multiprocessing.Pool(mp) as P:
                    ntargets = P.map(_get_ntargets_one, ntargs)
            else:
                ntargets = [get_ntargets_one(*_ntargs) for _ntargs in ntargs]
            ntargets = np.array(ntargets)

        if rank == 0:
            if len(ntargets) == 0:
                if redrockfiles is not None:
                    redrockfiles = []
                if outfiles is not None:
                    outfiles = []
            else:
                iempty = np.where(ntargets == 0)[0]
                if len(iempty) > 0:
                    log.info(f'Skipping {len(iempty):,d} redrockfiles with no targets.')

                itodo = np.where(ntargets > 0)[0]
                if len(itodo) > 0:
                    ntargets = ntargets[itodo]
                    log.info(f'Number of targets left: {np.sum(ntargets):,d}.')
                    if redrockfiles is not None:
                        redrockfiles = redrockfiles[itodo]
                        if coadd_type == 'healpix':
                            maxlen = str(len(max(np.atleast_1d(survey), key=len)) +
                                         len(max(np.atleast_1d(program), key=len)))
                            for onesurvey in np.atleast_1d(survey):
                                for oneprogram in np.atleast_1d(program):
                                    I = np.logical_and(np.char.find(redrockfiles, onesurvey) != -1,
                                                       np.char.find(redrockfiles, oneprogram) != -1)
                                    if np.any(I):
                                        one = f'{onesurvey}:{oneprogram}'
                                        log.info(f' {one:>{maxlen}}: {np.sum(I):,d} ' + \
                                                 f'redrockfiles, {np.sum(ntargets[I]):,d} targets')
                    if outfiles is not None:
                        outfiles = outfiles[itodo]
                else:
                    if redrockfiles is not None:
                        redrockfiles = []
                    if outfiles is not None:
                        outfiles = []

            if len(redrockfiles) == 0:
                log.info('All files have been processed!')
                return '', list(), list(), None

    return outdir, redrockfiles, outfiles, ntargets


def _read_to_merge_one(args):
    return read_to_merge_one(*args)


def read_to_merge_one(filename, fastphot):
    """Read metadata, specphot, and fastfit tables from one output file."""
    info = fitsio.FITS(filename)
    meta = Table(info['METADATA'].read())
    specphot = Table(info['SPECPHOT'].read())
    if fastphot:
        fastfit = None
    else:
        fastfit = Table(info['FASTSPEC'].read())
    return meta, specphot, fastfit


def _domerge(outfiles, outprefix=None, specprod=None, coadd_type=None,
             mergefile=None, fastphot=False, split_hdu=False,
             nside_main=None, mp=1):
    """Read and concatenate a set of per-healpix/tile output files into one catalog."""
    from astropy.table import vstack
    from desiutil.depend import getdep, hasdep
    from fastspecfit.io import write_fastspecfit

    t0 = time.time()
    meta, specphot, fastfit = [], [], []

    mpargs = [[outfile, fastphot] for outfile in outfiles]
    if mp > 1:
        with multiprocessing.Pool(mp) as P:
            out = P.map(_read_to_merge_one, mpargs)
    else:
        out = [read_to_merge_one(*mparg) for mparg in mpargs]
    out = list(zip(*out))

    meta = vstack(out[0])
    specphot = vstack(out[1])
    if not fastphot:
        fastfit = vstack(out[2])
    del out

    _ctx = f'nobj={len(meta):,d}, nfiles={len(outfiles)}'
    if outprefix:
        _ctx += f', prefix={outprefix}'
    log.info(fsftime('merge_catalogs', time.time()-t0, context=_ctx))

    hdr = fitsio.read_header(outfiles[0])

    # sensible defaults
    deps = {}
    deps['INPUTZ'] = False
    deps['INPUTS'] = False
    deps['CONSAGE'] = False
    deps['USEQNET'] = True
    deps['NMONTE'] = 50
    deps['SEED'] = 1
    deps['NOSCORR'] = False
    deps['NOPHOTO'] = False
    deps['BRDLFIT'] = True
    deps['UFLOOR'] = 0.01
    deps['SNRBBALM'] = 2.5
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

    write_fastspecfit(meta, specphot, fastfit, modelspectra=None, outfile=mergefile,
                      specprod=specprod, coadd_type=coadd_type, fastphot=fastphot,
                      fphotofile=deps2['FPHOTO_FILE'], template_file=deps2['FTEMPLATES_FILE'],
                      emlinesfile=deps2['EMLINES_FILE'], inputz=deps['INPUTZ'],
                      ignore_photometry=deps['NOPHOTO'], broadlinefit=deps['BRDLFIT'],
                      constrain_age=deps['CONSAGE'], use_quasarnet=deps['USEQNET'],
                      no_smooth_continuum=deps['NOSCORR'], split_hdu=split_hdu,
                      nside=nside_main)


def merge_fastspecfit(specprod=None, coadd_type=None, survey=None, program=None,
                      healpix=None, tile=None, night=None, sample=None, outsuffix=None,
                      fastphot=False, specprod_dir=None, outdir_data='.',
                      split_hdu=False, fastfiles_to_merge=None, merge_suffix=None,
                      mergedir=None, supermerge=False, overwrite=False, nside_main=None,
                      mp=1):
    """Merge per-healpix or per-tile output catalogs into a single catalog.

    Runs only on MPI rank 0. When ``supermerge=True``, merges previously
    merged per-survey/program catalogs into one master file instead.

    Parameters
    ----------
    specprod : :class:`str` or None, optional
        DESI spectroscopic production name.
    coadd_type : :class:`str` or None, optional
        Coadd type: ``'healpix'``, ``'cumulative'``, ``'pernight'``, etc.
    survey : :class:`str` or array-like of str, optional
        DESI survey name(s).
    program : :class:`str` or array-like of str, optional
        DESI program name(s).
    healpix : array-like or None, optional
        Specific HEALPix pixel(s) to merge.
    tile : array-like or None, optional
        Specific tile ID(s) to merge.
    night : array-like or None, optional
        Specific observing night(s) to merge.
    sample : :class:`astropy.table.Table` or None, optional
        Input sample catalog; overrides file-system discovery.
    outsuffix : :class:`str` or None, optional
        Suffix for the merged output filename; defaults to ``specprod``.
    fastphot : :class:`bool`, optional
        If ``True``, merge ``fastphot`` output files. Default is ``False``.
    specprod_dir : :class:`str` or None, optional
        Override the standard specprod directory path.
    outdir_data : :class:`str`, optional
        Root output directory. Default is ``'.'``.
    split_hdu : :class:`bool`, optional
        If ``True``, write main-survey outputs as separate HDUs by
        HEALPix nside. Default is ``False``.
    fastfiles_to_merge : :class:`list` or None, optional
        Explicit list of files to merge when ``supermerge=True``.
    merge_suffix : :class:`str` or None, optional
        Filename suffix for sample-based merges. Default is ``'sample'``.
    mergedir : :class:`str` or None, optional
        Directory for merged output files; defaults to
        ``<outdir_data>/<specprod>/catalogs``.
    supermerge : :class:`bool`, optional
        If ``True``, merge previously merged per-survey/program catalogs
        into one master file. Default is ``False``.
    overwrite : :class:`bool`, optional
        If ``True``, overwrite existing merged output files. Default is
        ``False``.
    nside_main : :class:`int` or None, optional
        HEALPix nside for splitting main-survey output by healpixel.
    mp : :class:`int`, optional
        Number of multiprocessing workers. Default is 1.

    """
    import fitsio
    from copy import copy
    from astropy.io import fits
    from astropy.table import Table, vstack
    from fastspecfit.mpi import plan

    if split_hdu and nside_main:
        errmsg = 'Using split_hdu=True and nside_main!=None is not supported.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    if fastphot:
        outprefix = 'fastphot'
    else:
        outprefix = 'fastspec'

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
            _domerge(outfiles, mergefile=mergefile, outprefix=outprefix,
                     specprod=specprod, coadd_type=coadd_type, fastphot=fastphot,
                     split_hdu=split_hdu, nside_main=nside_main, mp=mp)
        else:
            log.info(f'No catalogs found: {_outfiles}')
        return

    if sample is not None:
        if merge_suffix is None:
            merge_suffix = 'sample'
        if nside_main:
            print('WRITE ME!')
        else:
            mergefile = os.path.join(mergedir, f'{outprefix}-{specprod}-{merge_suffix}.fits')
            if os.path.isfile(mergefile) and not overwrite:
                log.info(f'Merged output file {mergefile} exists!')
                return

        _, _, outfiles, _ = plan(specprod=specprod, sample=sample, merge=True,
                                 fastphot=fastphot, specprod_dir=specprod_dir,
                                 outdir_data=outdir_data, overwrite=overwrite)
        if len(outfiles) > 0:
            _domerge(outfiles, mergefile=mergefile, outprefix=outprefix,
                     specprod=specprod, coadd_type=coadd_type, fastphot=fastphot,
                     split_hdu=split_hdu, nside_main=nside_main, mp=mp)

    elif coadd_type == 'healpix' and sample is None:
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
                _, _, outfiles, _ = plan(specprod=specprod, survey=survey, program=program, healpix=healpix,
                                         merge=True, fastphot=fastphot, specprod_dir=specprod_dir,
                                         outdir_data=outdir_data, overwrite=overwrite)
                if len(outfiles) > 0:
                    # split main/dark and main/bright merge files
                    # either by HDU or nside_main healpixels
                    do_split_hdu = False
                    use_nside_main = None
                    if survey == 'main' and ((program == 'dark') or (program == 'bright')):
                        if nside_main:
                            use_nside_main = nside_main
                        elif split_hdu:
                            do_split_hdu = True
                        
                    _domerge(outfiles, mergefile=mergefile, outprefix=outprefix, specprod=specprod,
                             coadd_type=coadd_type, fastphot=fastphot, split_hdu=do_split_hdu, 
                             nside_main=use_nside_main, mp=mp)
    else:
        mergefile = os.path.join(mergedir, f'{outprefix}-{specprod}-{coadd_type}.fits')
        if os.path.isfile(mergefile) and not overwrite:
            log.info(f'Merged output file {mergefile} exists!')
            return
        _, _, outfiles, _ = plan(specprod=specprod, coadd_type=coadd_type, tile=tile, night=night,
                                 merge=True, fastphot=fastphot, specprod_dir=specprod_dir,
                                 outdir_data=outdir_data, overwrite=overwrite)
        if len(outfiles) > 0:
            _domerge(outfiles, mergefile=mergefile, outprefix=outprefix,
                     specprod=specprod, coadd_type=coadd_type, fastphot=fastphot,
                     split_hdu=split_hdu, nside_main=nside_main, mp=mp)


# ---------------------------------------------------------------------------
# build_cmdargs: originally in bin/mpi-fastspecfit, moved here so that the
# per-file worker (_perfile_run) can call it from the fastspecfit.mpi module.
# ---------------------------------------------------------------------------

def build_cmdargs(args, redrockfile, outfile, sample=None, fastphot=False,
                  input_redshifts=False):
    """Build the command-line argument string for fastspec, fastphot, or fastqa.

    Parameters
    ----------
    args : :class:`argparse.Namespace`
        Parsed command-line arguments.
    redrockfile : :class:`str`
        Path to the input redrock file (or QA output directory for ``--makeqa``).
    outfile : :class:`str`
        Path to the output FITS file (or input file for ``--makeqa``).
    sample : :class:`astropy.table.Table` or None, optional
        Optional target sample table with columns ``{SURVEY, PROGRAM, HEALPIX, TARGETID}``.
    fastphot : :class:`bool`, optional
        If ``True``, build arguments for ``fastphot`` instead of ``fastspec``.
    input_redshifts : :class:`bool`, optional
        If ``True``, include input redshifts from ``sample['Z']``.

    Returns
    -------
    cmd : :class:`str`
        Executable name (``'fastspec'``, ``'fastphot'``, or ``'fastqa'``).
    cmdargs : :class:`str`
        Command-line argument string.
    logfile : :class:`str`
        Path to the log file.

    """
    if args.makeqa:
        cmd = 'fastqa'
        cmdargs = f'{outfile} -o={redrockfile} --mp={args.mp}'
    else:
        cmd = 'fastphot' if fastphot else 'fastspec'
        cmdargs = f'{redrockfile} -o={outfile} --mp={args.mp}'

        if args.ignore_quasarnet:
            cmdargs += ' --ignore-quasarnet'
        if args.ignore_photometry:
            cmdargs += ' --ignore-photometry'
        if args.no_smooth_continuum:
            cmdargs += ' --no-smooth-continuum'
        if args.templates:
            cmdargs += f' --templates={args.templates}'
        if args.templateversion:
            cmdargs += f' --templateversion={args.templateversion}'
        if args.fphotodir:
            cmdargs += f' --fphotodir={args.fphotodir}'
        if args.fphotofile:
            cmdargs += f' --fphotofile={args.fphotofile}'
        if args.emlinesfile:
            cmdargs += f' --emlinesfile={args.emlinesfile}'
        if args.nmonte:
            cmdargs += f' --nmonte={args.nmonte}'
        if args.vdisp_nominal:
            cmdargs += f' --vdisp-nominal={args.vdisp_nominal}'
        if args.vdisp_bounds:
            cmdargs += f' --vdisp-bounds {float(args.vdisp_bounds[0])} {float(args.vdisp_bounds[1])}'
        if args.seed:
            cmdargs += f' --seed={args.seed}'

        if sample is not None:
            _, survey, program, healpix = os.path.basename(redrockfile).split('-')
            healpix = int(healpix.split('.')[0])
            I = ((sample['SURVEY'] == survey) * (sample['PROGRAM'] == program) *
                 (sample['HEALPIX'] == healpix))
            targetids = ','.join(sample[I]['TARGETID'].astype(str))
            cmdargs += f' --targetids={targetids}'
            if input_redshifts:
                inputz = ','.join(sample[I]['Z'].astype(str))
                cmdargs += f' --input-redshifts={inputz}'
        else:
            if args.targetids is not None:
                cmdargs += f' --targetids={args.targetids}'

    if args.verbose:
        cmdargs += ' --verbose'
    if args.ntargets is not None:
        cmdargs += f' --ntargets={args.ntargets}'
    if args.firsttarget is not None:
        cmdargs += f' --firsttarget={args.firsttarget}'

    if args.makeqa:
        logfile = os.path.join(redrockfile, os.path.basename(outfile).replace('.gz', '').replace('.fits', '.log'))
    else:
        logfile = outfile.replace('.gz', '').replace('.fits', '.log')

    return cmd, cmdargs, logfile


# ---------------------------------------------------------------------------
# Per-file parallel worker state and functions (used by MPPool when
# --mp-per-file is set).  Placing these in the fastspecfit.mpi module (rather
# than in bin/mpi-fastspecfit) makes them picklable by fully-qualified name,
# which is required when multiprocessing uses the 'spawn' start method.
#
# _perfile_init() is the pool initializer: called once per worker at startup.
# _perfile_run()  is the per-task function: called once per healpix file.
# ---------------------------------------------------------------------------

_perfile_sample = None
_perfile_rank = None
_perfile_fastphot = False
_perfile_input_redshifts = False


def _perfile_init(init_sc_args, sample, rank, fastphot, input_redshifts):
    """Pool initializer for --mp-per-file mode.

    Loads sc_data once per worker process and stores per-run state in
    module-level globals so that _perfile_run() can access them without
    per-task pickling overhead.
    """
    global _perfile_sample, _perfile_rank, _perfile_fastphot, _perfile_input_redshifts
    from fastspecfit.singlecopy import sc_data
    sc_data.initialize(**init_sc_args)
    _perfile_sample = sample
    _perfile_rank = rank
    _perfile_fastphot = fastphot
    _perfile_input_redshifts = input_redshifts


def _perfile_run(redrockfile, outfile, ntarget, args):
    """Per-file worker task for --mp-per-file mode.

    Reads worker state from module-level globals set by _perfile_init(),
    builds per-file command args (forcing --mp=1 so no inner pool is
    created), and calls fastspec/fastphot for this file.
    """
    import copy
    import fastspecfit.mpi as _mpi
    from desispec.parallel import stdouterr_redirected

    if _mpi._perfile_fastphot:
        from fastspecfit.fastspecfit import fastphot as fast
    else:
        from fastspecfit.fastspecfit import fastspec as fast

    # Force mp=1: all parallelism is at the file level; no inner pool needed.
    args_inner = copy.copy(args)
    args_inner.mp = 1

    cmd, cmdargs, logfile = build_cmdargs(
        args_inner, redrockfile, outfile,
        sample=_mpi._perfile_sample,
        fastphot=_mpi._perfile_fastphot,
        input_redshifts=_mpi._perfile_input_redshifts,
    )

    if _mpi._perfile_rank == 0:
        log.info(f'Rank {_mpi._perfile_rank}: ntargets={ntarget}: {cmd} {cmdargs}')

    if getattr(args, 'dry_run', False):
        return

    try:
        t1 = time.time()
        outdir = os.path.dirname(logfile)
        if not os.path.isdir(outdir):
            os.makedirs(outdir, exist_ok=True)

        if getattr(args, 'nolog', False):
            err = fast(args=cmdargs.split(), comm=None)
        else:
            with stdouterr_redirected(to=logfile, overwrite=args.overwrite, comm=None):
                err = fast(args=cmdargs.split(), comm=None)

        log.info(fsftime('process_file', time.time()-t1,
                         context=f'rank={_mpi._perfile_rank}'))
        if err != 0:
            if not os.path.exists(outfile):
                log.warning(f'Rank {_mpi._perfile_rank} missing {outfile}')
                raise IOError
    except:
        log.warning(f'Rank {_mpi._perfile_rank} raised an exception')
        import traceback
        traceback.print_exc()
