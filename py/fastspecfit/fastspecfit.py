#!/usr/bin/env python
"""
fastspecfit.fastspecfit
=======================

See sandbox/running-fastspecfit for examples.

"""
import pdb # for debugging

import os, time
import numpy as np
import astropy.units as u

from fastspecfit.continuum import ContinuumTools

from desiutil.log import get_logger
log = get_logger()

## ridiculousness! - this seems to come from healpy, blarg
#import tempfile
#os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

def _fastspec_one(args):
    """Multiprocessing wrapper."""
    return fastspec_one(*args)

def _fastphot_one(args):
    """Multiprocessing wrapper."""
    return fastphot_one(*args)

def _desiqa_one(args):
    """Multiprocessing wrapper."""
    return desiqa_one(*args)

def _assign_units_to_columns(fastfit, metadata, Spec, CFit, EMFit=None, fastphot=False):
    """Assign astropy units to output tables."""
    fastcols = fastfit.colnames
    metacols = metadata.colnames

    T, M = Spec.init_output(CFit=CFit, EMFit=EMFit, fastphot=fastphot)
    for col in T.colnames:
        if col in fastcols:
            fastfit[col].unit = T[col].unit
    for col in M.colnames:
        if col in metacols:
            metadata[col].unit = M[col].unit

    if EMFit is not None:
        E = EMFit.init_output(nobj=1)
        for col in E.colnames:
            if col in fastcols:
                fastfit[col].unit = E[col].unit

def fastspec_one(iobj, data, out, meta, CFit, EMFit, broadlinefit=True):
    """Multiprocessing wrapper to run :func:`fastspec` on a single object."""
    
    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()

    #cfit, continuummodel, smooth_continuum = CFit.continuum_specfit(data)
    cfit, continuummodel, smooth_continuum = CFit.continuum_joint_specfit(data)
    for col in cfit.colnames:
        out[col] = cfit[col]

    # Copy over the reddening-corrected fluxes -- messy!
    for iband, band in enumerate(CFit.fiber_bands):
        meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    for iband, band in enumerate(CFit.bands):
        meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} [targetid {}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    # Fit the emission-line spectrum.
    t0 = time.time()
    emfit, emmodel = EMFit.fit(data, continuummodel, smooth_continuum,
                               broadlinefit=broadlinefit)
    for col in emfit.colnames:
        out[col] = emfit[col]
    log.info('Line-fitting object {} [targetid={}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    return out, meta, emmodel

def fastphot_one(iobj, data, out, meta, CFit):
    """Multiprocessing wrapper to run :func:`fastphot` on a single object."""

    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()
    cfit, _ = CFit.continuum_fastphot(data)
    for col in cfit.colnames:
        out[col] = cfit[col]

    # Copy over the reddening-corrected fluxes -- messy!
    for iband, band in enumerate(CFit.fiber_bands):
        meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    for iband, band in enumerate(CFit.bands):
        meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} [targetid {}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))
    
    return out, meta

def desiqa_one(CFit, EMFit, data, fastfit, metadata, coadd_type,
               fastphot=False, outdir=None, outprefix=None, webqa=False):
    """Multiprocessing wrapper to generate QA for a single object."""

    #t0 = time.time()
    if webqa:
        build_webqa(CFit, EMFit, data, fastfit, metadata, coadd_type=coadd_type,
                    outprefix=outprefix, outdir=outdir)
    elif fastphot:
        CFit.qa_fastphot(data, fastfit, metadata, coadd_type=coadd_type,
                         outprefix=outprefix, outdir=outdir)
    else:
        EMFit.qa_joint_fastspec(data, fastfit, metadata, coadd_type=coadd_type,
                                outprefix=outprefix, outdir=outdir)
    #log.info('Building took {:.2f} sec'.format(time.time()-t0))

def parse(options=None):
    """Parse input arguments to fastspec and fastphot scripts.

    """
    import argparse, sys

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('redrockfiles', nargs='*', help='Full path to input redrock file(s).')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename (required).')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads per MPI rank.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file, zero-indexed.') 
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--solve-vdisp', action='store_true', help='Solve for the velocity dispersion (only when using fastspec).')
    parser.add_argument('--no-broadlinefit', default=True, action='store_false', dest='broadlinefit',
                        help='Do not allow for broad Balmer and Helium line-fitting.')
    parser.add_argument('--ssptemplates', type=str, default=None, help='Optional name of the SSP templates.')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--specfile-prefix', type=str, default='coadd-', help='Prefix of the spectral file(s).')
    parser.add_argument('--qnfile-prefix', type=str, default='qso_qn-', help='Prefix of the QuasarNet afterburner file(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--dr9dir', type=str, default=None, help='Optional directory name for the DR9 photometry.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspec {}'.format(' '.join(options)))

    return args

def fastspec(args=None, comm=None):
    """Main fastspec script.

    This script is the engine to model one or more DESI spectra. It initializes
    the :class:`ContinuumFit` and :class:`EMLineFit` classes, reads the data, fits
    each spectrum (with the option of fitting in parallel), and writes out the
    results as a multi-extension binary FITS table.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line. 
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    from astropy.table import Table, vstack
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum- and emission-line fitting classes. Note: trim
    # the wavelengths of the SSPs to optimize compute time.
    t0 = time.time()
    CFit = ContinuumFit(ssptemplates=args.ssptemplates, mapdir=args.mapdir, 
                        verbose=args.verbose, solve_vdisp=args.solve_vdisp, 
                        minwave=500.0,# maxwave=1e4)
                        # dev - increase maxwave
                        maxwave=40e4)
    EMFit = EMLineFit(mapdir=args.mapdir, ssptemplates=args.ssptemplates,
                      verbose=args.verbose)
    Spec = DESISpectra(dr9dir=args.dr9dir)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    Spec.select(args.redrockfiles, firsttarget=args.firsttarget,
                targetids=targetids, ntargets=args.ntargets,
                redrockfile_prefix=args.redrockfile_prefix,
                specfile_prefix=args.specfile_prefix,
                qnfile_prefix=args.qnfile_prefix)
    if len(Spec.specfiles) == 0:
        return

    data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, mp=args.mp)
    log.info('Read data for {} objects in {:.2f} sec'.format(Spec.ntargets, time.time()-t0))

    #np.savetxt('linemask3.txt', np.array([np.hstack(data[0]['wave']), np.hstack(data[0]['flux']), np.hstack(data[0]['ivar'])]).T)

    out, meta = Spec.init_output(CFit=CFit, EMFit=EMFit, fastphot=False)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, EMFit, args.broadlinefit)
               for iobj in np.arange(Spec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastspec_one, fitargs)
    else:
        _out = [fastspec_one(*_fitargs) for _fitargs in fitargs]
    _out = list(zip(*_out))
    out = Table(np.hstack(_out[0]))
    meta = Table(np.hstack(_out[1]))
    try:
        # need to vstack to preserve the wavelength metadata 
        modelspectra = vstack(_out[2], metadata_conflicts='error')
    except:
        errmsg = 'Metadata conflict when stacking model spectra.'
        log.critical(errmsg)
        raise ValueError(errmsg)
       
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # Assign units and write out.
    _assign_units_to_columns(out, meta, Spec, CFit, EMFit=EMFit, fastphot=False)

    write_fastspecfit(out, meta, modelspectra=modelspectra, outfile=args.outfile,
                      specprod=Spec.specprod, coadd_type=Spec.coadd_type,
                      fastphot=False)

def fastphot(args=None, comm=None):
    """Main fastphot script.

    This script is the engine to model the broadband photometry of one or more
    DESI objects. It initializes the :class:`ContinuumFit` class, reads the
    data, fits each object (with the option of fitting in parallel), and writes
    out the results as a multi-extension binary FITS table.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line. 
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum-fitting classes.
    t0 = time.time()
    CFit = ContinuumFit(ssptemplates=args.ssptemplates, mapdir=args.mapdir, 
                        minwave=None, maxwave=40e4, solve_vdisp=False, 
                        cache_vdisp=False, verbose=args.verbose)

    Spec = DESISpectra(dr9dir=args.dr9dir)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    Spec.select(args.redrockfiles, firsttarget=args.firsttarget,
                targetids=targetids, ntargets=args.ntargets,
                redrockfile_prefix=args.redrockfile_prefix,
                specfile_prefix=args.specfile_prefix,
                qnfile_prefix=args.qnfile_prefix)
    if len(Spec.specfiles) == 0:
        return
    data = Spec.read_and_unpack(CFit, fastphot=True, synthphot=False, mp=args.mp)
    
    out, meta = Spec.init_output(CFit=CFit, fastphot=True)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit)
               for iobj in np.arange(Spec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastphot_one, fitargs)
    else:
        _out = [fastphot_one(*_fitargs) for _fitargs in fitargs]
    _out = list(zip(*_out))
    out = Table(np.hstack(_out[0]))
    meta = Table(np.hstack(_out[1]))
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # Assign units and write out.
    _assign_units_to_columns(out, meta, Spec, CFit, fastphot=True)

    write_fastspecfit(out, meta, outfile=args.outfile, specprod=Spec.specprod,
                      coadd_type=Spec.coadd_type, fastphot=True)

def build_webqa(CFit, EMFit, data, fastfit, metadata, coadd_type='healpix',
                spec_wavelims=(3550, 9900), outprefix=None, outdir=None):
    """QA for the web page.

    """
    import subprocess
    from scipy.ndimage import median_filter

    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.ticker as ticker
    from matplotlib.patches import Circle, Rectangle
    from matplotlib.lines import Line2D

    import seaborn as sns
    from astropy.table import Table, Column, vstack
    from PIL import Image, ImageDraw
    
    from fastspecfit.util import ivar2var, C_LIGHT

    sns.set(context='talk', style='ticks', font_scale=1.1)#, rc=rc)

    col1 = [colors.to_hex(col) for col in ['dodgerblue', 'darkseagreen', 'orangered']]
    col2 = [colors.to_hex(col) for col in ['darkblue', 'darkgreen', 'darkred']]
    col3 = [colors.to_hex(col) for col in ['blue', 'green', 'red']]

    photcol1 = colors.to_hex('darkblue') # 'darkgreen', 'darkred', 'dodgerblue', 'darkseagreen', 'orangered']]

    if outdir is None:
        outdir = '.'
    if outprefix is None:
        outprefix = 'fastfit'

    if metadata['PHOTSYS'] == 'S':
        filters = CFit.decam
        allfilters = CFit.decamwise
    else:
        filters = CFit.bassmzls
        allfilters = CFit.bassmzlswise

    if coadd_type == 'healpix':
        title = 'Survey/Program/Healpix: {}/{}/{}, TARGETID: {}'.format(
            metadata['SURVEY'], metadata['PROGRAM'],
            metadata['HEALPIX'], metadata['TARGETID'])
        pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
    elif coadd_type == 'cumulative':
        title = 'Survey/Program/Tile/Fiber: {}/{}/{}/{}, TARGETID: {}'.format(
            metadata['SURVEY'].upper(), metadata['PROGRAM'], metadata['TILEID'], metadata['FIBER'], metadata['TARGETID'])
        #title = '{}/{}: {} on Tile/Fiber: {}/{}'.format(
        #    metadata['SURVEY'].upper(), metadata['PROGRAM'], metadata['TARGETID'], metadata['TILEID'], metadata['FIBER'])
        pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                outprefix, metadata['TILEID'], coadd_type, metadata['TARGETID']))
    elif coadd_type == 'pernight':
        title = 'Tile/Night: {}/{}, TARGETID/Fiber: {}/{}'.format(
                metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID'],
                metadata['FIBER'])
        pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                outprefix, metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID']))
    elif coadd_type == 'perexp':
        title = 'Tile/Night/Expid: {}/{}/{}, TARGETID/Fiber: {}/{}'.format(
                metadata['TILEID'], metadata['NIGHT'], metadata['EXPID'],
                metadata['TARGETID'], metadata['FIBER'])
        pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                outprefix, metadata['TILEID'], metadata['NIGHT'],
                metadata['EXPID'], metadata['TARGETID']))

    redshift = fastfit['CONTINUUM_Z']

    # rebuild the best-fitting broadband photometric model
    inodust = np.ndarray.item(np.where(CFit.AV == 0)[0]) # should always be index 0
    continuum_phot, synthmodelphot = CFit.SSP2data(
        CFit.sspflux_dustnomvdisp[:, :, inodust], CFit.sspwave, redshift=redshift,
        synthphot=True, AV=fastfit['CONTINUUM_AV_PHOT'], #test=True,
        coeff=fastfit['CONTINUUM_COEFF_PHOT'] * CFit.massnorm)

    continuum_wave_phot = CFit.sspwave * (1 + redshift)

    wavemin, wavemax = 0.1, 35.0 # 6.0
    indx_phot = np.where((continuum_wave_phot/1e4 > wavemin) * (continuum_wave_phot/1e4 < wavemax))[0]     

    phot = CFit.parse_photometry(CFit.bands,
                                 maggies=np.array([metadata['FLUX_{}'.format(band.upper())] for band in CFit.bands]),
                                 ivarmaggies=np.array([metadata['FLUX_IVAR_{}'.format(band.upper())] for band in CFit.bands]),
                                 lambda_eff=allfilters.effective_wavelengths.value,
                                 min_uncertainty=CFit.min_uncertainty)
    fiberphot = CFit.parse_photometry(CFit.fiber_bands,
                                      maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in CFit.fiber_bands]),
                                      lambda_eff=filters.effective_wavelengths.value)

    # rebuild the best-fitting spectroscopic model
    stackwave = np.hstack(data['wave'])

    continuum, _ = EMFit.SSP2data(EMFit.sspflux, EMFit.sspwave, redshift=redshift, 
                                 specwave=data['wave'], specres=data['res'],
                                 cameras=data['cameras'],
                                 AV=fastfit['CONTINUUM_AV_SPEC'],
                                 vdisp=fastfit['CONTINUUM_VDISP'],
                                 coeff=fastfit['CONTINUUM_COEFF_SPEC'],
                                 synthphot=False)

    residuals = [data['flux'][icam] - continuum[icam] for icam in np.arange(len(data['cameras']))]
    if np.all(fastfit['CONTINUUM_COEFF_SPEC'] == 0):
        _smooth_continuum = np.zeros_like(stackwave)
    else:
        _smooth_continuum, _ = EMFit.smooth_continuum(np.hstack(data['wave']), np.hstack(residuals),
                                                      np.hstack(data['ivar']), redshift=redshift,
                                                      linemask=np.hstack(data['linemask']))
    smooth_continuum = []
    for campix in data['camerapix']:
        smooth_continuum.append(_smooth_continuum[campix[0]:campix[1]])

    _emlinemodel = EMFit.emlinemodel_bestfit(data['wave'], data['res'], Table(fastfit))

    # individual-line spectra
    _emlinemodel_oneline = []
    for refline in EMFit.linetable: # [EMFit.inrange]: # for all lines in range
        T = Table(fastfit['CONTINUUM_Z', 'MGII_DOUBLET_RATIO', 'OII_DOUBLET_RATIO', 'SII_DOUBLET_RATIO'])
        for oneline in EMFit.linetable: # need all lines for the model
            linename = oneline['name']
            for linecol in ['AMP', 'VSHIFT', 'SIGMA']:
                col = linename.upper()+'_'+linecol
                if linename == refline['name']:
                    T.add_column(Column(name=col, data=fastfit[col], dtype=fastfit[col].dtype))
                else:
                    T.add_column(Column(name=col, data=0.0, dtype=fastfit[col].dtype))
        # special case the parameter doublets
        if refline['name'] == 'mgii_2796':
            T['MGII_2803_AMP'] = fastfit['MGII_2803_AMP']
        if refline['name'] == 'oii_3726':
            T['OII_3729_AMP'] = fastfit['OII_3729_AMP']
        if refline['name'] == 'sii_6731':
            T['SII_6716_AMP'] = fastfit['SII_6716_AMP']
        _emlinemodel_oneline1 = EMFit.emlinemodel_bestfit(data['wave'], data['res'], T)
        if np.sum(np.hstack(_emlinemodel_oneline1)) > 0:
            _emlinemodel_oneline.append(_emlinemodel_oneline1)

    # Grab the viewer cutout.
    width = int(40 / 0.262)   # =1 arcmin
    height = int(width / 1.3) # 3:2 aspect ratio
    
    cutoutpng = os.path.join('/tmp', 'tmp.'+os.path.basename(pngfile))
    if not os.path.isfile(cutoutpng):
        cmd = 'wget -O {outfile} https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&width={width}&height={height}&layer=ls-dr9'
        cmd = cmd.format(outfile=cutoutpng, ra=metadata['RA'], dec=metadata['DEC'],
                         width=width, height=height)
        print(cmd)
        err = subprocess.call(cmd.split())
        if err != 0:
            errmsg = 'Something went wrong retrieving the png cutout'
            log.critical(errmsg)
            raise ValueError(errmsg)

    # QA choices - everything in inches
    npanelrows = 7
    height_onepanel = 1.6
    fullheight = npanelrows*height_onepanel
    
    npanelcols = 13 # 15
    width_onepanel = 1.3
    fullwidth = npanelcols*width_onepanel
    
    nlinerows = 4 
    nlinecols = 6 
    nlinepanels = nlinecols * nlinerows
    
    nspecrows = 3
    nspeccols = 8 
    
    nsedrows = npanelrows - nspecrows
    nsedcols = npanelcols - nlinecols
    
    ncutrows = nspecrows
    ncutcols = npanelcols - nspeccols

    try:
        assert((ncutrows + nsedrows) == npanelrows)
        assert((nspecrows + nlinerows) == npanelrows)
        assert((ncutcols + nspeccols) == npanelcols)
        assert((nsedcols + nlinecols) == npanelcols)
    except:
        errmsg = 'QA setup has gone awry.'
        log.critical(errmsg)
        raise ValueError(errmsg)

    #print(fullwidth, fullheight)
    #print(ncutcols, nspeccols, npanelcols)
    #print(ncutrows, nspecrows, nsedrows, npanelrows)
    
    height_ratios = np.hstack(([1.2]*nspecrows, [1.0]*nlinerows))

    plt.clf()
    fig = plt.figure(figsize=(fullwidth, fullheight))
    gs = fig.add_gridspec(npanelrows, npanelcols)#, height_ratios=height_ratios)#, width_ratios=width_ratios)

    cutax = fig.add_subplot(gs[0:ncutrows, 0:ncutcols])
    sedax = fig.add_subplot(gs[ncutrows:(ncutrows+nsedrows), 0:nsedcols])
    specax = fig.add_subplot(gs[:nspecrows, ncutcols:(ncutcols+nspeccols)])

    # viewer cutout
    Image.MAX_IMAGE_PIXELS = None
    with Image.open(cutoutpng) as im:
        cutax.imshow(im, interpolation='nearest')
        sz = im.size

    if metadata['DEC'] > 0:
        sgn = '+'
    else:
        sgn = ''
        
    bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.7)
    cutax.text(0.04, 0.95, '$(\\alpha,\\delta)$=({:.7f}, {}{:.6f})'.format(metadata['RA'], sgn, metadata['DEC']),
               ha='left', va='top', color='k', fontsize=11, bbox=bbox,
               transform=cutax.transAxes)

    pixscale = 0.262
    cutax.add_artist(Circle((sz[0] / 2, sz[1] / 2), radius=1.5/2/pixscale, facecolor='none', # DESI fiber=1.5 arcsec diameter
                            edgecolor='red', ls='-', alpha=0.5))#, label='3" diameter'))
    cutax.add_artist(Circle((sz[0] / 2, sz[1] / 2), radius=10/2/pixscale, facecolor='none',
                            edgecolor='red', ls='--', alpha=0.5))#, label='15" diameter'))
    handles = [Line2D([0], [0], color='red', lw=2, ls='-', label='1.5 arcsec'),
               Line2D([0], [0], color='red', lw=2, ls='--', label='10 arcsec')]
    
    cutax.get_xaxis().set_visible(False)
    cutax.get_yaxis().set_visible(False)
    cutax.axis('off')
    cutax.autoscale(False)
    cutax.legend(handles=handles, loc='lower left', fontsize=10, facecolor='lightgray')

    # full spectrum + best-fitting continuum model
    ymin_spec, ymax_spec = 1e6, -1e6
    allfullspec, allfullwave = [], []
    for ii in np.arange(len(data['cameras'])): # iterate over cameras
        sigma, good = ivar2var(data['ivar'][ii], sigma=True, allmasked_ok=True)

        specax.fill_between(data['wave'][ii], data['flux'][ii]-sigma,
                            data['flux'][ii]+sigma, color=col1[ii])
        _fullspec = continuum[ii] + _emlinemodel[ii]
        fullspec = continuum[ii] + smooth_continuum[ii] + _emlinemodel[ii]
        specax.plot(data['wave'][ii], fullspec, color=col2[ii])

        allfullspec.append(_fullspec)
        allfullwave.append(data['wave'][ii])
        
        # get the robust range
        filtflux = median_filter(data['flux'][ii], 51, mode='nearest')
        I = data['ivar'][ii] > 0
        if np.sum(I) > 0:
            #ymin_spec, ymax_spec = np.percentile(data['flux'][ii], [5, 98]) * [-1, +1]
            #ymin_spec, ymax_spec = np.percentile(data['flux'][ii], [5, 95])
            sigflux = np.diff(np.percentile(data['flux'][ii], [25, 75]))[0] / 1.349 # robust
            if -3 * sigflux < ymin_spec:
                ymin_spec = -3 * sigflux
            if 5 * sigflux > ymax_spec:
                ymax_spec = 5 * sigflux
            if np.max(filtflux) > ymax_spec:
                ymax_spec = np.max(filtflux)
            if np.max(_emlinemodel[ii]) > ymax_spec:
                ymax_spec = np.max(_emlinemodel[ii]) * 1.2
            #print(ymin_spec, ymax_spec)

    allfullwave = np.hstack(allfullwave)
    allfullspec = np.hstack(allfullspec)
    
    #specax.plot(stackwave, _smooth_continuum, color='gray')
    specax.plot(stackwave, np.hstack(continuum), color='k', alpha=0.3)

    specax.set_xlim(spec_wavelims)
    specax.set_ylim(ymin_spec, ymax_spec)
    specax.set_xticklabels([])
    specax.set_yticklabels([])
    specax.set_xticks([])
    specax.set_yticks([])
    specax.spines[['left', 'bottom', 'top']].set_visible(False)

    specax_twin = specax.twinx()
    specax_twin.set_ylim(ymin_spec, ymax_spec)
    specax_twin.spines[['left', 'bottom', 'top']].set_visible(False)
    
    # photometric SED
    if fastfit['FLUX_SYNTH_R'] > 0:
        apfactor = synthmodelphot[synthmodelphot['band'] == 'r']['nanomaggies'][0]/fastfit['FLUX_SYNTH_R']
    else:
        apfactor = 1.0
    factor = apfactor * 10**(0.4 * 48.6) * allfullwave**2 / (C_LIGHT * 1e13) / CFit.fluxnorm # [erg/s/cm2/A --> maggies]
    good = allfullspec > 0
    sedax.plot(allfullwave[good]/1e4, -2.5*np.log10(allfullspec[good]*factor[good]), color='gray', alpha=0.5)

    if np.all(continuum_phot[indx_phot] <= 0):
        CFit.log.warning('Best-fitting photometric continuum is all zeros or negative!')
        continuum_phot_abmag = continuum_phot*0 + np.median(fiberphot['abmag'])
    else:
        indx_phot = indx_phot[continuum_phot[indx_phot] > 0] # trim zeros
        factor = 10**(0.4 * 48.6) * continuum_wave_phot[indx_phot]**2 / (C_LIGHT * 1e13) / CFit.fluxnorm / CFit.massnorm # [erg/s/cm2/A --> maggies]
        continuum_phot_abmag = -2.5*np.log10(continuum_phot[indx_phot] * factor)
        sedax.plot(continuum_wave_phot[indx_phot] / 1e4, continuum_phot_abmag, color='tan', zorder=1)

    sedax.scatter(synthmodelphot['lambda_eff']/1e4, synthmodelphot['abmag'], 
               marker='s', s=200, color='k', facecolor='none',
               #label=r'$grz$ (spectrum, synthesized)',
               alpha=0.8, zorder=2)
    
    # we have to set the limits *before* we call errorbar, below!
    dm = 1.0
    good = phot['abmag_ivar'] > 0
    goodlim = phot['abmag_limit'] > 0
    if np.sum(good) > 0 and np.sum(goodlim) > 0:
        ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(phot['abmag_limit'][goodlim]), np.nanmax(continuum_phot_abmag))) + dm
        ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(phot['abmag_limit'][goodlim]), np.nanmin(continuum_phot_abmag))) - dm
    elif np.sum(good) > 0 and np.sum(goodlim) == 0:
        ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(continuum_phot_abmag))) + dm
        ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(continuum_phot_abmag))) - dm
    elif np.sum(good) == 0 and np.sum(goodlim) > 0:
        ymin = np.max((np.nanmax(phot['abmag_limit'][goodlim]), np.nanmax(continuum_phot_abmag))) + dm
        ymax = np.min((np.nanmin(phot['abmag_limit'][goodlim]), np.nanmin(continuum_phot_abmag))) - dm
    else:
        good = phot['abmag'] > 0
        goodlim = phot['abmag_limit'] > 0
        if np.sum(good) > 0 and np.sum(goodlim) > 0:
            ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(phot['abmag_limit'][goodlim]))) + dm
            ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(phot['abmag_limit'][goodlim]))) - dm
        elif np.sum(good) > 0 and np.sum(goodlim) == 0:                
            ymin = np.nanmax(phot['abmag'][good]) + dm
            ymax = np.nanmin(phot['abmag'][good]) - dm
        elif np.sum(good) == 0 and np.sum(goodlim) > 0:
            ymin = np.nanmax(phot['abmag_limit'][goodlim]) + dm
            ymax = np.nanmin(phot['abmag_limit'][goodlim]) - dm
        else:
            ymin, ymax = [30, 20]
        
    if ymin > 30:
        ymin = 30
    if np.isnan(ymin) or np.isnan(ymax):
        raise('Problem here!')

    sedax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
    sedax.set_ylabel('AB mag') 
    #sedax.set_ylabel(r'Apparent Brightness (AB mag)') 
    sedax.set_xlim(wavemin, wavemax)
    sedax.set_ylim(ymin, ymax)
    sedax.set_xscale('log')

    @ticker.FuncFormatter
    def major_formatter(x, pos):
        if x > 1:
            return f'{x:.0f}'
        else:
            return f'{x:.1f}'
    
    sedax.xaxis.set_major_formatter(major_formatter)
    sedax.set_xticks([0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 20.0])

    # integrated flux / photometry
    abmag = np.squeeze(phot['abmag'])
    abmag_limit = np.squeeze(phot['abmag_limit'])
    abmag_fainterr = np.squeeze(phot['abmag_fainterr'])
    abmag_brighterr = np.squeeze(phot['abmag_brighterr'])
    yerr = np.squeeze([abmag_brighterr, abmag_fainterr])

    dofit = np.where(CFit.bands_to_fit)[0]
    if len(dofit) > 0:
        good = np.where((abmag[dofit] > 0) * (abmag_limit[dofit] == 0))[0]
        upper = np.where(abmag_limit[dofit] > 0)[0]
        if len(good) > 0:
            sedax.errorbar(phot['lambda_eff'][dofit][good]/1e4, abmag[dofit][good],
                        yerr=yerr[:, dofit[good]],
                        fmt='o', markersize=12, markeredgewidth=3, markeredgecolor=photcol1,
                        markerfacecolor=photcol1, elinewidth=3, ecolor=photcol1, capsize=4,
                        label=r'$grz\,W_{1}W_{2}W_{3}W_{4}$', zorder=2)
        if len(upper) > 0:
            sedax.errorbar(phot['lambda_eff'][dofit][upper]/1e4, abmag_limit[dofit][upper],
                        lolims=True, yerr=0.75,
                        fmt='o', markersize=12, markeredgewidth=3, markeredgecolor=photcol1,
                        markerfacecolor=photcol1, elinewidth=3, ecolor=photcol1, capsize=4)

    ignorefit = np.where(CFit.bands_to_fit == False)[0]
    if len(ignorefit) > 0:
        good = np.where((abmag[ignorefit] > 0) * (abmag_limit[ignorefit] == 0))[0]
        upper = np.where(abmag_limit[ignorefit] > 0)[0]
        if len(good) > 0:
            sedax.errorbar(phot['lambda_eff'][ignorefit][good]/1e4, abmag[ignorefit][good],
                        yerr=yerr[:, ignorefit[good]],
                        fmt='o', markersize=12, markeredgewidth=3, markeredgecolor=photcol1,
                        markerfacecolor='none', elinewidth=3, ecolor=photcol1, capsize=4)
        if len(upper) > 0:
            sedax.errorbar(phot['lambda_eff'][ignorefit][upper]/1e4, abmag_limit[ignorefit][upper],
                        lolims=True, yerr=0.75, fmt='o', markersize=12, markeredgewidth=3,
                        markeredgecolor=photcol1, markerfacecolor='none', elinewidth=3,
                        ecolor=photcol1, capsize=5)

    sedax.plot([spec_wavelims[0]/1e4, spec_wavelims[1]/1e4], [ymin-1, ymin-1],
               lw=2, ls='-', color='gray', marker='s')#, alpha=0.5)
    sedax.text(((spec_wavelims[1]-spec_wavelims[0])/2+spec_wavelims[0]*0.8)/1e4, ymin-1.4,
               'DESI x {:.2f}'.format(apfactor), ha='center', va='center', fontsize=10,
               color='gray')
    
    # zoom in on individual emission lines - use linetable!
    plotsig_default = 200.0 # [km/s]
    plotsig_default_balmer = 500.0 # [km/s]
    plotsig_default_broad = 2000.0 # [km/s]
    
    linetable = EMFit.linetable # custom plotgroups and labels
    #linetable.remove_columns(['isbalmer', 'isbroad', 'amp'])

    if True:
        lya = Table()
        lya['name'] = ['lya']
        lya['restwave'] = [1215.67]
        lya['nicename'] = ['Ly$\\alpha$-$\lambda$1215']
        lya['isbroad'] = [True]
        lya['isbalmer'] = [False]
        lya['plotgroup'] = ['lya_1215']
        linetable = vstack((lya, linetable))

        fastfit = Table(fastfit)
        fastfit['LYA_SIGMA'] = [plotsig_default_broad]
        fastfit = fastfit[0] # stupid astropy Table/Row hack

    linetable = linetable[np.logical_and.reduce((
        linetable['name'] != 'hei_4471',
        #linetable['name'] != 'hei_5876',
        linetable['name'] != 'heii_4686',
        linetable['name'] != 'hei_broad_4471',
        #linetable['name'] != 'hei_broad_5876',
        linetable['name'] != 'heii_broad_4686',
        ))]

    # custom groups
    linetable['plotgroup'][linetable['name'] == 'hgamma'] = 'hgamma'
    linetable['plotgroup'][linetable['name'] == 'hgamma_broad'] = 'hgamma'
    linetable['nicename'][linetable['name'] == 'hgamma'] = 'H$\gamma$-$\lambda$4340'
    linetable['nicename'][linetable['name'] == 'hgamma_broad'] = 'H$\gamma$-$\lambda$4340'

    linetable['plotgroup'][linetable['name'] == 'oiii_4363'] = 'oiii_4363'
    linetable['nicename'][linetable['name'] == 'oiii_4363'] = '[OIII]-$\lambda4363$'

    #linetable['plotgroup'][linetable['name'] == 'oiii_4959'] = 'oiii_4959'
    #linetable['nicename'][linetable['name'] == 'oiii_4959'] = '[OIII]-$\lambda4959$'
    #linetable['plotgroup'][linetable['name'] == 'oiii_5007'] = 'oiii_5007'
    #linetable['nicename'][linetable['name'] == 'oiii_5007'] = '[OIII]-$\lambda5007$'

    #linetable['linetype'] = 

    uplotgroups = np.unique(linetable['plotgroup'].data)
    nline = len(uplotgroups)
    
    meanwaves, deltawaves, sigmas, linenames = [], [], [], []
    for plotgroup in uplotgroups:
        I = np.where(plotgroup == linetable['plotgroup'])[0]
        linenames.append(linetable['nicename'][I[0]].replace('-', ' '))
        meanwaves.append(np.mean(linetable['restwave'][I]))
        deltawaves.append((np.max(linetable['restwave'][I]) - np.min(linetable['restwave'][I])) / 2)
    
        sigmas1 = np.array([fastfit['{}_SIGMA'.format(line.upper())] for line in linetable[I]['name']])
        sigmas1 = sigmas1[sigmas1 > 0]
        if len(sigmas1) > 0:
            plotsig = 1.5*np.mean(sigmas1)
        else:
            if np.any(linetable['isbroad'][I]):
                if np.any(linetable['isbalmer'][I]):
                    plotsig = fastfit['BROAD_SIGMA']
                    if plotsig < 50:
                        plotsig = fastfit['NARROW_SIGMA']
                        if plotsig < 50:
                            plotsig = plotsig_default
                            #plotsig = plotsig_default_broad
                else:
                    plotsig = fastfit['UV_SIGMA']                    
                    if plotsig < 50:
                        plotsig = plotsig_default_broad
            else:
                plotsig = fastfit['NARROW_SIGMA']
                if plotsig < 50:
                    plotsig = plotsig_default
        sigmas.append(plotsig)
    
    #srt = np.argsort(meanwaves)
    order = [
        # UV/QSO - 6
        'lya_1215', 'oi_1304', 'siliv_1396', 'civ_1549', 'siliii_1892_ciii_1908', 'mgii_2796_2803',
        # Balmer + nebular - 12
        'nev_3346', 'nev_3426', 'oii_3726_29', 'neiii_3869_h6', 'hepsilon', 'hdelta',
        'hgamma', 'hbeta',
        #'oiii_4959', 'oiii_5007',
        'oiii_doublet', 'hei_5876',
        'halpha_nii_6548_48', 'sii_6716_31',
        # auroral - 6
        'oiii_4363', 'nii_5755', 'oi_6300_siii_6312', 'oii_7320_30', 'siii_9069', 'siii_9532'
        ]
    srt = np.hstack([np.where(uplotgroups == grp)[0] for grp in order])

    meanwaves = np.hstack(meanwaves)[srt]
    deltawaves = np.hstack(deltawaves)[srt]
    sigmas = np.hstack(sigmas)[srt]
    linenames = np.hstack(linenames)[srt]

    # add linelabels to the specx plot
    #for line in linetable:
    #    meanwave = line['restwave']*(1+redshift)
    #    if meanwave > spec_wavelims[0] and meanwave < spec_wavelims[1]:
    #        specax.text(meanwave, ymax_spec, line['nicename'], ha='center', va='top',
    #                    rotation=270, fontsize=7, alpha=0.5)
    for meanwave, linename in zip(meanwaves*(1+redshift), linenames):
        #print(meanwave, ymax_spec)
        if meanwave > spec_wavelims[0] and meanwave < spec_wavelims[1]:
            if 'SiIII' in linename:
                thislinename = '\n'+linename.replace('+', '+\n  ')
            elif '4363' in linename:
                thislinename = linename+'\n'
            else:
                thislinename = linename
            specax.text(meanwave, ymax_spec, thislinename, ha='center', va='top',
                        rotation=270, fontsize=7, alpha=0.5)
    
    removelabels = np.ones(nline, bool)
    ymin, ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6

    rowoffset = ncutrows
    coloffset = nsedcols

    ax, irow, icol = [], 0, 0
    for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):    
    #for iax in np.arange(nlinepanels):
        icol = iax % nlinecols
        if iax > 0 and iax % nlinecols == 0:
            irow += 1
        #print(iax, icol+coloffset, irow+rowoffset)
    
        xx = fig.add_subplot(gs[irow+rowoffset, icol+coloffset])
        ax.append(xx)
    
        wmin = (meanwave - deltawave) * (1+redshift) - 6 * sig * meanwave * (1+redshift) / C_LIGHT
        wmax = (meanwave + deltawave) * (1+redshift) + 6 * sig * meanwave * (1+redshift) / C_LIGHT
        #print(linename, wmin, wmax)
    
        # iterate over cameras
        for ii in np.arange(len(data['cameras'])): # iterate over cameras
            emlinewave = data['wave'][ii]
            emlineflux = data['flux'][ii] - continuum[ii] - smooth_continuum[ii]
            emlinemodel = _emlinemodel[ii]
    
            emlinesigma, good = ivar2var(data['ivar'][ii], sigma=True, allmasked_ok=True)
            emlinewave = emlinewave[good]
            emlineflux = emlineflux[good]
            emlinesigma = emlinesigma[good]
            emlinemodel = emlinemodel[good]
    
            emlinemodel_oneline = []
            for _emlinemodel_oneline1 in _emlinemodel_oneline:
                emlinemodel_oneline.append(_emlinemodel_oneline1[ii][good])
    
            indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
            if len(indx) > 1:
                removelabels[iax] = False
                xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                                emlineflux[indx]+emlinesigma[indx], color=col1[ii],
                                alpha=0.5)
                # plot the individual lines first then the total model
                for emlinemodel_oneline1 in emlinemodel_oneline:
                    if np.sum(emlinemodel_oneline1[indx]) > 0:
                        #P = emlinemodel_oneline1[indx] > 0
                        xx.plot(emlinewave[indx], emlinemodel_oneline1[indx], lw=1, alpha=0.8, color=col2[ii])
                xx.plot(emlinewave[indx], emlinemodel[indx], color=col2[ii], lw=2)
                    
                # get the robust range
                sigflux = np.std(emlineflux[indx])
                filtflux = median_filter(emlineflux[indx], 3, mode='nearest')
    
                _ymin, _ymax = -1.5 * sigflux, 5 * sigflux
                if np.max(emlinemodel[indx]) > _ymax:
                    _ymax = np.max(emlinemodel[indx]) * 1.3
                if np.max(filtflux) > _ymax:
                    _ymax = np.max(filtflux)
                if np.min(emlinemodel[indx]) < _ymin:
                    _ymin = 0.8 * np.min(emlinemodel[indx])
                    
                if _ymax > ymax[iax]:
                    ymax[iax] = _ymax
                if _ymin < ymin[iax]:
                    ymin[iax] = _ymin
    
                xx.set_xlim(wmin, wmax)

            if 'SiIII' in linename or '6312' in linename:
                thislinename = '\n'+linename.replace('+', '+\n  ')
            elif '2803' in linename:
                thislinename = '\n'+linename.replace(',2803', ',\n           2803')
            else:
                thislinename = linename
            xx.text(0.04, 0.89, thislinename, ha='left', va='center',
                    transform=xx.transAxes, fontsize=7)
            
    for iax, xx in enumerate(ax):
        if removelabels[iax]:
            xx.set_ylim(0, 1)
            xx.set_xticklabels([])
            xx.set_yticklabels([])
        else:
            xx.set_ylim(ymin[iax], ymax[iax])
            #xlim = xx.get_xlim()
            #xx.xaxis.set_major_locator(ticker.MaxNLocator(2))
            xx.set_xticklabels([])
            xx.set_yticklabels([])
        xx.axis('off')

    #plt.subplots_adjust(wspace=0.27, top=tp, bottom=bt, left=lf, right=rt, hspace=0.22)
    plt.subplots_adjust(bottom=0.22, top=0.94, left=0.08, right=0.91)

    # add some key results about the object at the bottom of the figure
    leg = {
        'radec': '$(\\alpha,\\delta)=({:.7f},{:.6f})$'.format(metadata['RA'], metadata['DEC']),
        'z': '$z={:.7f}$'.format(redshift),
        'zredrock': '$z_{{\\rm Redrock}}={:.7f}$'.format(metadata['Z_RR']),
        #'zredrock': '$z_{{\\rm redrock}}={:.7f}$'.format(redshift),
        'dn4000_spec': '$D_{{n}}(4000)_{{\\rm spec}}={:.3f}$'.format(fastfit['DN4000']),
        'dn4000_modelspec': '$D_{{n}}(4000)_{{\\rm spec,model}}={:.3f}$'.format(fastfit['DN4000_MODEL_SPEC']),
        'dn4000_modelphot': '$D_{{n}}(4000)_{{\\rm phot,model}}={:.3f}$'.format(fastfit['DN4000_MODEL_PHOT']),
        'chi2': '$\\chi^{{2}}_{{\\nu}}={:.3f}$'.format(fastfit['CONTINUUM_RCHI2_SPEC']),
        'rchi2': '$\\chi^{{2}}_{{\\nu}}={:.3f}$'.format(fastfit['RCHI2']),
        'deltarchi2': '$\\Delta\\chi^{{2}}_{{\\nu,\\rm broad,narrow}}={:.3f}$'.format(fastfit['DELTA_LINERCHI2']),
        'age_spec': '<Age>$_{{\\rm spec}}={:.3f}$ Gyr'.format(fastfit['CONTINUUM_AGE_SPEC']),
        'age_phot': '<Age>$_{{\\rm phot}}={:.3f}$ Gyr'.format(fastfit['CONTINUUM_AGE_PHOT']),
        'absmag_r': '$M_{{^{{0.1}}r}}={:.2f}$'.format(fastfit['ABSMAG_SDSS_R']),
        'absmag_gr': '$^{{0.1}}(g-r)={:.3f}$'.format(fastfit['ABSMAG_SDSS_G']-fastfit['ABSMAG_SDSS_R']),
        'absmag_rz': '$^{{0.1}}(r-z)={:.3f}$'.format(fastfit['ABSMAG_SDSS_R']-fastfit['ABSMAG_SDSS_Z']),       
        }

    if fastfit['NARROW_Z'] == redshift:
        leg.update({'dv_narrow': '$\\Delta v_{{\\rm narrow}}=$...'})
    else:
        leg.update({'dv_narrow': '$\\Delta v_{{\\rm narrow}}={:.2f}$ km/s'.format(C_LIGHT*(fastfit['NARROW_Z']-redshift))})
    if fastfit['BROAD_Z'] == redshift:
        leg.update({'dv_broad': '$\\Delta v_{{\\rm broad}}=$...'})
    else:
        leg.update({'dv_broad': '$\\Delta v_{{\\rm broad}}={:.2f}$ km/s'.format(C_LIGHT*(fastfit['BROAD_Z']-redshift))})        
    if fastfit['UV_Z'] == redshift:
        leg.update({'dv_uv': '$\\Delta v_{{\\rm UV}}=$...'})
    else:
        leg.update({'dv_uv': '$\\Delta v_{{\\rm UV}}={:.2f}$ km/s'.format(C_LIGHT*(fastfit['UV_Z']-redshift))})

    if fastfit['NARROW_SIGMA'] == 0.0:
        leg.update({'sigma_narrow': '$\\sigma_{{\\rm narrow}}=$...'})
    else:
        leg.update({'sigma_narrow': '$\\sigma_{{\\rm narrow}}={:.1f}$ km/s'.format(fastfit['NARROW_SIGMA'])})
    if fastfit['BROAD_SIGMA'] == 0.0:
        leg.update({'sigma_broad': '$\\sigma_{{\\rm broad}}=$...'})
    else:
        leg.update({'sigma_broad': '$\\sigma_{{\\rm broad}}={:.1f}$ km/s'.format(fastfit['BROAD_SIGMA'])})
    if fastfit['UV_SIGMA'] == 0.0:
        leg.update({'sigma_uv': '$\\sigma_{{\\rm UV}}=$...'})
    else:
        leg.update({'sigma_uv': '$\\sigma_{{\\rm UV}}={:.1f}$ km/s'.format(fastfit['UV_SIGMA'])})
        
    if fastfit['CONTINUUM_VDISP_IVAR'] == 0:
        leg.update({'vdisp': '$\\sigma_{{\\rm star}}={:.1f}$ km/s'.format(fastfit['CONTINUUM_VDISP'])})
    else:
        leg.update({'vdisp': '$\\sigma_{{\\rm star}}={:.1f}\\pm{:.1f}$ km/s'.format(
            fastfit['CONTINUUM_VDISP'], 1/np.sqrt(fastfit['CONTINUUM_VDISP_IVAR']))})

    if fastfit['LOGMSTAR'] > 0:
        leg.update({'mstar': '$\\log_{{10}}\,(M_{{*}}/M_{{\odot}})={:.3f}$'.format(fastfit['LOGMSTAR'])})
    else:
        leg.update({'mstar': '$\\log_{{10}}\,(M_{{*}}/M_{{\odot}})=$...'})

    if fastfit['CONTINUUM_AV_IVAR_SPEC'] == 0:
        leg.update({'AV': '$A(V)={:.3f}$ mag'.format(fastfit['CONTINUUM_AV_SPEC'])})
    else:
        leg.update({'AV': '$A(V)={:.3f}\\pm{:.3f}$ mag'.format(
            fastfit['CONTINUUM_AV_SPEC'], 1/np.sqrt(fastfit['CONTINUUM_AV_IVAR_SPEC']))})

    # emission lines
    if fastfit['CIV_1549_EW'] == 0:
        leg.update({'ewciv': 'EW(CIV)$=$...'})
    else:
        leg.update({'ewciv': 'EW(CIV)$={:.3f}\ \\AA$'.format(fastfit['CIV_1549_EW'])})
        
    if fastfit['MGII_2796_EW'] == 0 and fastfit['MGII_2796_EW'] == 0:
        leg.update({'ewmgii': 'EW(MgII)$=$...'})
    else:
        leg.update({'ewmgii': 'EW(MgII)$={:.3f}\ \\AA$'.format(fastfit['MGII_2796_EW']+fastfit['MGII_2796_EW'])})
        
    if fastfit['HALPHA_EW'] == 0:
        leg.update({'ewha_narrow': 'EW(H$\\alpha)_{{\\rm narrow}}=$...'})
    else:
        leg.update({'ewha_narrow': 'EW(H$\\alpha)_{{\\rm narrow}}={:.2f}\ \\AA$'.format(fastfit['HALPHA_EW'])})
        
    if fastfit['HBETA_EW'] == 0:
        leg.update({'ewhb_narrow': 'EW(H$\\beta)_{{\\rm narrow}}=$...'})
    else:
        leg.update({'ewhb_narrow': 'EW(H$\\beta)_{{\\rm narrow}}={:.2f}\ \\AA$'.format(fastfit['HBETA_EW'])})
        
    if fastfit['HGAMMA_EW'] == 0:
        leg.update({'ewhg_narrow': 'EW(H$\\gamma)_{{\\rm narrow}}=$...'})
    else:
        leg.update({'ewhg_narrow': 'EW(H$\\gamma)_{{\\rm narrow}}={:.2f}\ \\AA$'.format(fastfit['HGAMMA_EW'])})
        
    if fastfit['HALPHA_BROAD_EW'] == 0:
        leg.update({'ewha_broad': 'EW(H$\\alpha)_{{\\rm broad}}=$...'})
    else:
        leg.update({'ewha_broad': 'EW(H$\\alpha)_{{\\rm broad}}={:.2f}\ \\AA$'.format(fastfit['HALPHA_BROAD_EW'])})
        
    if fastfit['HBETA_BROAD_EW'] == 0:
        leg.update({'ewhb_broad': 'EW(H$\\beta)_{{\\rm broad}}=$...'})
    else:
        leg.update({'ewhb_broad': 'EW(H$\\beta)_{{\\rm broad}}={:.2f}\ \\AA$'.format(fastfit['HBETA_BROAD_EW'])})
        
    if fastfit['HGAMMA_BROAD_EW'] == 0:
        leg.update({'ewhg_broad': 'EW(H$\\gamma)_{{\\rm broad}}=$...'})
    else:
        leg.update({'ewhg_broad': 'EW(H$\\gamma)_{{\\rm broad}}={:.2f}\ \\AA$'.format(fastfit['HGAMMA_BROAD_EW'])})
        
    if fastfit['OII_3726_EW'] == 0 and fastfit['OII_3729_EW'] == 0:
        leg.update({'ewoii': 'EW([OII])$=$...'})
    else:
        leg.update({'ewoii': 'EW([OII])$={:.2f}\ \\AA$'.format(fastfit['OII_3726_EW']+fastfit['OII_3729_EW'])})
        
    if fastfit['OIII_5007_EW'] == 0:
        leg.update({'ewoiii': 'EW([OIII])$=$...'})
    else:
        leg.update({'ewoiii': 'EW([OIII])$={:.2f}\ \\AA$'.format(fastfit['OIII_5007_EW'])})
        
    if fastfit['NII_6584_EW'] == 0:
        leg.update({'ewnii': 'EW([NII])$=$...'})
    else:
        leg.update({'ewnii': 'EW([NII])$={:.2f}\ \\AA$'.format(fastfit['NII_6584_EW'])})
        
    if fastfit['SII_6716_EW'] == 0 and fastfit['SII_6731_EW'] == 0:
        leg.update({'ewsii': 'EW([SII])$=$...'})
    else:
        leg.update({'ewsii': 'EW([SII])$={:.2f}\ \\AA$'.format(fastfit['SII_6716_EW']+fastfit['SII_6731_EW'])})
        
    if fastfit['OII_DOUBLET_RATIO'] == 0:
        #leg.update({'oii_doublet': '[OII] doublet ratio$=$...'})
        leg.update({'oii_doublet': '[OII] $\lambda3726/\lambda3729=$...'})
    else:
        #leg.update({'oii_doublet': '[OII] doublet ratio$={:.3f}$'.format(fastfit['OII_DOUBLET_RATIO'])})
        leg.update({'oii_doublet': '[OII] $\lambda3726/\lambda3729={:.3f}$'.format(fastfit['OII_DOUBLET_RATIO'])})

    if fastfit['SII_DOUBLET_RATIO'] == 0:
        #leg.update({'sii_doublet': '[SII] doublet ratio$=$...'})
        leg.update({'sii_doublet': '[SII] $\lambda6731/\lambda6716=$...'})
    else:
        #leg.update({'sii_doublet': '[SII] doublet ratio$={:.3f}$'.format(fastfit['SII_DOUBLET_RATIO'])})
        leg.update({'sii_doublet': '[SII] $\lambda6731/\lambda6716={:.3f}$'.format(fastfit['SII_DOUBLET_RATIO'])})

    # parse the targeting bits
    #from desitarget.targets import main_cmx_or_sv
    #(desi_target, bgs_target, mws_target), mask, survey = main_cmx_or_sv(metadata)
    
    legfntsz, toppos, startpos, deltapos = 12, 0.125, 0.07, 0.13
    txt = '\n'.join((
        r'{}'.format(leg['mstar']),
        r'{}'.format(leg['absmag_r']),
        r'{}'.format(leg['absmag_gr']),
        r'{}'.format(leg['absmag_rz']),
        #r'{}'.format(leg['AV']),
        ))
    fig.text(startpos, toppos, txt, ha='left', va='top', fontsize=legfntsz)

    txt = '\n'.join((
        r'{}'.format(leg['z']),
        r'{}'.format(leg['zredrock']),
        r'{}'.format(leg['dv_narrow']),
        r'{}'.format(leg['dv_broad']),
        r'{}'.format(leg['dv_uv']),
        ))
    fig.text(startpos+deltapos*1, toppos, txt, ha='left', va='top', fontsize=legfntsz)

    txt = '\n'.join((
        r'{}'.format(leg['age_spec']),
        r'{}'.format(leg['age_phot']),
        #r'{}'.format(leg['dn4000_spec']),
        r'{}'.format(leg['dn4000_modelspec']),
        r'{}'.format(leg['dn4000_modelphot']),
        ))
    fig.text(startpos+deltapos*2, toppos, txt, ha='left', va='top', fontsize=legfntsz)

    txt = '\n'.join((
        r'{}'.format(leg['vdisp']),
        r'{}'.format(leg['sigma_narrow']),
        r'{}'.format(leg['sigma_broad']),
        r'{}'.format(leg['sigma_uv']),
        '',
        r'{}'.format(leg['oii_doublet']),
        r'{}'.format(leg['sii_doublet']),
        ))
    fig.text(startpos+deltapos*3.2, toppos+0.05, txt, ha='left', va='top', fontsize=legfntsz)

    txt = '\n'.join((
        r'{}'.format(leg['ewciv']),
        r'{}'.format(leg['ewmgii']),
        '',
        r'{}'.format(leg['ewoii']),
        r'{}'.format(leg['ewoiii']),
        r'{}'.format(leg['ewnii']),
        r'{}'.format(leg['ewsii']),
        ))
    fig.text(startpos+deltapos*3.2+1.2*deltapos, toppos+0.06, txt, ha='left', va='top', fontsize=legfntsz)

    txt = '\n'.join((
        r'{}'.format(leg['ewhg_narrow']),
        r'{}'.format(leg['ewhb_narrow']),
        r'{}'.format(leg['ewha_narrow']),
        '',
        r'{}'.format(leg['ewhg_broad']),
        r'{}'.format(leg['ewhb_broad']),
        r'{}'.format(leg['ewha_broad']),
        ))
    fig.text(startpos+deltapos*3.2+1.2*deltapos+1.1*deltapos, toppos+0.07, txt, ha='left', va='top', fontsize=legfntsz)

    # common axis labels
    ppos = specax.get_position()
    x0, y0, x1, y1 = ppos.x0, ppos.y0, ppos.x1, ppos.y1
    plt.text(x1+0.05, (y1-y0)/2+y0, 
             r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
             ha='center', va='center', 
             transform=fig.transFigure, rotation=270)#, fontsize=)

    # draw a box around the three groups of lines -- note that this code has to
    # come after the subplots_adjust!
    ax = np.array(ax).reshape(nlinerows, nlinecols)

    # UV/QSO
    llpos = ax[0, 0].get_position()
    urpos = ax[0, nlinecols-1].get_position()
    x0, y0, x1, y1 = llpos.x0, llpos.y0, urpos.x1, urpos.y1
    width, height = urpos.x1-llpos.x0, urpos.y1-llpos.y0
    fig.patches.extend([Rectangle((x0, y0), width, height, # (lower-left corner), width, height
                                  fill=False, color='purple', lw=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    plt.text(x1+0.01, (y1-y0)/2+y0, 'UV/QSO', ha='center', va='center',
             transform=fig.transFigure, rotation=270, fontsize=10)

    llpos = ax[2, 0].get_position()
    urpos = ax[1, nlinecols-1].get_position()
    x0, y0, x1, y1 = llpos.x0, llpos.y0, urpos.x1, urpos.y1
    width, height = urpos.x1-llpos.x0, urpos.y1-llpos.y0
    fig.patches.extend([Rectangle((x0, y0), width, height, # (lower-left corner), width, height
                                  fill=False, color='orange', lw=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    plt.text(x1+0.01, (y1-y0)/2+y0, 'Balmer+He+Narrow-Line', ha='center', va='center',
             transform=fig.transFigure, rotation=270, fontsize=10)

    llpos = ax[3, 0].get_position()
    urpos = ax[3, nlinecols-1].get_position()
    x0, y0, x1, y1 = llpos.x0, llpos.y0, urpos.x1, urpos.y1
    width, height = urpos.x1-llpos.x0, urpos.y1-llpos.y0
    fig.patches.extend([Rectangle((x0, y0), width, height, # (lower-left corner), width, height
                                  fill=False, color='firebrick', lw=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    plt.text(x1+0.01, (y1-y0)/2+y0, 'Auroral', ha='center', va='center',
             transform=fig.transFigure, rotation=270, fontsize=10)

    fig.suptitle(title, fontsize=22)

    EMFit.log.info('Writing {}'.format(pngfile))
    fig.savefig(pngfile)
    plt.close()

class FastFit(ContinuumTools):
    def __init__(self, ssptemplates=None, minwave=None, maxwave=40e4, nolegend=False,
                 solve_vdisp=True, constrain_age=False, mapdir=None, verbose=False):
        """Class to model a galaxy stellar continuum.

        Parameters
        ----------
        metallicity : :class:`str`, optional, defaults to `Z0.0190`.
            Stellar metallicity of the SSPs. Currently fixed at solar
            metallicity, Z=0.0190.
        minwave : :class:`float`, optional, defaults to None
            Minimum SSP wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxwave : :class:`float`, optional, defaults to 6e4
            Maximum SSP wavelength to read into memory. 

        Notes
        -----
        Need to document all the attributes.
        
        Plans for improvement:
          - Update the continuum redshift using cross-correlation.
          - Don't draw reddening from a flat distribution (try gamma or a custom
            distribution of the form x**2*np.exp(-2*x/scale).

        """
        super(ContinuumFit, self).__init__(ssptemplates=ssptemplates, minwave=minwave,
                                           maxwave=maxwave, mapdir=mapdir, verbose=verbose)

        self.nolegend = nolegend
        self.constrain_age = constrain_age
        self.solve_vdisp = solve_vdisp

    def init_spec_output(self, nobj=1):
        """Initialize the output data table for this class.

        """
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        out.add_column(Column(name='APCORR', length=nobj, dtype='f4')) # aperture correction
        out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_RCHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        out.add_column(Column(name='VDISP_IVAR', length=nobj, dtype='f4', unit=u.second**2/u.kilometer**2))
        out.add_column(Column(name='AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='ZZSUN', length=nobj, dtype='f4'))
        out.add_column(Column(name='LOGMSTAR', length=nobj, dtype='f4', unit=u.solMass))
        out.add_column(Column(name='SFR', length=nobj, dtype='f4', unit=u.solMass/u.year))
        out.add_column(Column(name='FAGN', length=nobj, dtype='f4'))
        #out.add_column(Column(name='SFR50', length=nobj, dtype='f4', unit=u.solMass/u.year))
        #out.add_column(Column(name='SFR300', length=nobj, dtype='f4', unit=u.solMass/u.year))
        #out.add_column(Column(name='SFR1000', length=nobj, dtype='f4', unit=u.solMass/u.year))
        
        for cam in ['B', 'R', 'Z']:
            out.add_column(Column(name='CONTINUUM_SNR_{}'.format(cam), length=nobj, dtype='f4')) # median S/N in each camera
        #out.add_column(Column(name='CONTINUUM_SNR', length=nobj, shape=(3,), dtype='f4')) # median S/N in each camera

        # maximum correction to the median-smoothed continuum
        for cam in ['B', 'R', 'Z']:
            out.add_column(Column(name='CONTINUUM_SMOOTHCORR_{}'.format(cam), length=nobj, dtype='f4')) 
        out['CONTINUUM_AV'] = 0.0 #* u.mag
        out['CONTINUUM_VDISP'] = self.vdisp_nominal # * u.kilometer/u.second

        #if False:
        #    # continuum fit with *no* dust reddening (to be used as a diagnostic
        #    # tool to identify potential calibration issues).
        #    out.add_column(Column(name='CONTINUUM_NODUST_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        #    out.add_column(Column(name='CONTINUUM_NODUST_CHI2', length=nobj, dtype='f4')) # reduced chi2
        #    #out.add_column(Column(name='CONTINUUM_NODUST_AGE', length=nobj, dtype='f4', unit=u.Gyr))

        out.add_column(Column(name='DN4000', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_IVAR', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))

        return out

    def init_phot_output(self, nobj=1):
        """Initialize the photometric output data table.

        """
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        #out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_RCHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
        out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='CONTINUUM_AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='CONTINUUM_AV_IVAR', length=nobj, dtype='f4', unit=1/u.mag**2))
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))
        
        # observed-frame photometry synthesized from the best-fitting continuum model fit
        for band in self.bands:
            out.add_column(Column(name='FLUX_SYNTH_MODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies'))

        #if False:
        #    for band in self.fiber_bands:
        #        out.add_column(Column(name='FIBERTOTFLUX_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies')) # observed-frame fiber photometry
        #        #out.add_column(Column(name='FIBERTOTFLUX_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/'nanomaggies'**2))
        #    for band in self.bands:
        #        out.add_column(Column(name='FLUX_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies')) # observed-frame photometry
        #        out.add_column(Column(name='FLUX_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/'nanomaggies'**2))
                
        for band in self.absmag_bands:
            out.add_column(Column(name='KCORR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.mag))
            out.add_column(Column(name='ABSMAG_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.mag)) # absolute magnitudes
            out.add_column(Column(name='ABSMAG_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/u.mag**2))

        out.add_column(Column(name='LOGMSTAR', length=nobj, dtype='f4', unit=u.solMass))

        for cflux in ['LOGLNU_1500', 'LOGLNU_2800']:
            out.add_column(Column(name=cflux, length=nobj, dtype='f4', unit=10**(-28)*u.erg/u.second/u.Hz))
        out.add_column(Column(name='LOGL_5100', length=nobj, dtype='f4', unit=10**(10)*u.solLum))
        for cflux in ['FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT']:
            out.add_column(Column(name=cflux, length=nobj, dtype='f4', unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))

        return out

    def get_mean_property(self, physical_property, coeff, agekeep,
                          normalization=None, log10=False):
        """Compute the mean physical properties, given a set of coefficients.

        """
        values = self.sspinfo[physical_property][agekeep] # account for age of the universe trimming

        if np.count_nonzero(coeff > 0) == 0:
            self.log.warning('Coefficients are all zero!')
            meanvalue = 0.0
            #raise ValueError
        else:
            meanvalue = values.dot(coeff)
            # the coefficients include the stellar mass normalization
            if physical_property != 'mstar' and physical_property != 'sfr':
                meanvalue /= np.sum(coeff) 
            if normalization:
                meanvalue /= normalization
            if log10 and meanvalue > 0:
                meanvalue = np.log10(meanvalue)
        
        return meanvalue

    def younger_than_universe(self, redshift):
        """Return the indices of the SSPs younger than the age of the universe at the
        given redshift.

        """
        return np.where(self.sspinfo['age'] <= self.cosmo.age(redshift).to(u.year).value)[0]

    def kcorr_and_absmag(self, data, continuum, coeff, snrmin=2.0):
        """Computer K-corrections, absolute magnitudes, and a simple stellar mass.

        """
        from scipy.stats import sigmaclip
        from scipy.ndimage import median_filter
        
        redshift = data['zredrock']
        
        if data['photsys'] == 'S':
            filters_in = self.decamwise
        else:
            filters_in = self.bassmzlswise
        lambda_in = filters_in.effective_wavelengths.value

        # redshifted wavelength array and distance modulus
        zsspwave = self.sspwave * (1 + redshift)
        dmod = self.cosmo.distmod(redshift).value

        maggies = data['phot']['nanomaggies'].data * 1e-9
        ivarmaggies = (data['phot']['nanomaggies_ivar'].data / 1e-9**2) * self.bands_to_fit # mask W2-W4

        # input bandpasses, observed frame; maggies and bestmaggies should be
        # very close.
        bestmaggies = filters_in.get_ab_maggies(continuum / self.fluxnorm, zsspwave)
        bestmaggies = np.array(bestmaggies.as_array().tolist()[0])

        # need to handle filters with band_shift!=0 separately from those with band_shift==0
        def _kcorr_and_absmag(filters_out, band_shift):
            nout = len(filters_out)

            # note the factor of 1+band_shift
            lambda_out = filters_out.effective_wavelengths.value / (1 + band_shift)

            # Multiply by (1+z) to convert the best-fitting model to the "rest
            # frame" and then divide by 1+band_shift to shift it and the
            # wavelength vector to the band-shifted redshift. Also need one more
            # factor of 1+band_shift in order maintain the AB mag normalization.
            synth_outmaggies_rest = filters_out.get_ab_maggies(continuum * (1 + redshift) / (1 + band_shift) /
                                                               self.fluxnorm, self.sspwave * (1 + band_shift))
            synth_outmaggies_rest = np.array(synth_outmaggies_rest.as_array().tolist()[0]) / (1 + band_shift)
    
            # output bandpasses, observed frame
            synth_outmaggies_obs = filters_out.get_ab_maggies(continuum / self.fluxnorm, zsspwave)
            synth_outmaggies_obs = np.array(synth_outmaggies_obs.as_array().tolist()[0])
    
            absmag = np.zeros(nout, dtype='f4')
            ivarabsmag = np.zeros(nout, dtype='f4')
            kcorr = np.zeros(nout, dtype='f4')
            for jj in np.arange(nout):
                lambdadist = np.abs(lambda_in / (1 + redshift) - lambda_out[jj])
                # K-correct from the nearest "good" bandpass (to minimizes the K-correction)
                #oband = np.argmin(lambdadist)
                #oband = np.argmin(lambdadist + (ivarmaggies == 0)*1e10)
                oband = np.argmin(lambdadist + (maggies*np.sqrt(ivarmaggies) < snrmin)*1e10)
                kcorr[jj] = + 2.5 * np.log10(synth_outmaggies_rest[jj] / bestmaggies[oband])

                # m_R = M_Q + DM(z) + K_QR(z) or
                # M_Q = m_R - DM(z) - K_QR(z)
                if maggies[oband] * np.sqrt(ivarmaggies[oband]) > snrmin:
                #if (maggies[oband] > 0) and (ivarmaggies[oband]) > 0:
                    absmag[jj] = -2.5 * np.log10(maggies[oband]) - dmod - kcorr[jj]
                    ivarabsmag[jj] = maggies[oband]**2 * ivarmaggies[oband] * (0.4 * np.log(10.))**2
                else:
                    # if we use synthesized photometry then ivarabsmag is zero
                    # (which should never happen?)
                    absmag[jj] = -2.5 * np.log10(synth_outmaggies_rest[jj]) - dmod

                check = absmag[jj], -2.5*np.log10(synth_outmaggies_rest[jj]) - dmod
                self.log.debug(check)

            return kcorr, absmag, ivarabsmag

        kcorr_01, absmag_01, ivarabsmag_01 = _kcorr_and_absmag(self.absmag_filters_01, band_shift=0.1)
        kcorr_00, absmag_00, ivarabsmag_00 = _kcorr_and_absmag(self.absmag_filters_00, band_shift=0.0)

        nout = len(self.absmag_bands)
        kcorr = np.zeros(nout, dtype='f4')
        absmag = np.zeros(nout, dtype='f4')
        ivarabsmag = np.zeros(nout, dtype='f4')

        I00 = np.isin(self.absmag_bands, self.absmag_bands_00)
        I01 = np.isin(self.absmag_bands, self.absmag_bands_01)

        kcorr[I00] = kcorr_00
        absmag[I00] = absmag_00
        ivarabsmag[I00] = ivarabsmag_00

        kcorr[I01] = kcorr_01
        absmag[I01] = absmag_01
        ivarabsmag[I01] = ivarabsmag_01

        #print(kcorr_01, absmag_01)
        #pdb.set_trace()
        
        # get the stellar mass
        nage = len(coeff)
        mstar = self.sspinfo['mstar'][:nage].dot(coeff) * self.massnorm
        
        # From Taylor+11, eq 8
        #https://researchportal.port.ac.uk/ws/files/328938/MNRAS_2011_Taylor_1587_620.pdf
        #mstar = 1.15 + 0.7*(absmag[1]-absmag[3]) - 0.4*absmag[3]

        # compute the model continuum flux at 1500 and 2800 A (to facilitate UV
        # luminosity-based SFRs) and at the positions of strong nebular emission
        # lines [OII], Hbeta, [OIII], and Halpha
        dfactor = (1 + redshift) * 4.0 * np.pi * self.cosmo.luminosity_distance(redshift).to(u.cm).value**2 / self.fluxnorm

        lums, cfluxes = {}, {}
        cwaves = [1500.0, 2800.0, 5100.0]
        labels = ['LOGLNU_1500', 'LOGLNU_2800', 'LOGL_5100']
        norms = [1e28, 1e28, 1e10]
        for cwave, norm, label in zip(cwaves, norms, labels):
            J = (self.sspwave > cwave-500) * (self.sspwave < cwave+500)
            I = (self.sspwave[J] > cwave-20) * (self.sspwave[J] < cwave+20)
            smooth = median_filter(continuum[J], 200)
            clipflux, _, _ = sigmaclip(smooth[I], low=1.5, high=3)
            cflux = np.median(clipflux) # [flux in 10**-17 erg/s/cm2/A]
            cflux *= dfactor # [monochromatic luminosity in erg/s/A]
            if label == 'LOGL_5100':
                cflux *= cwave / 3.846e33 / norm # [luminosity in 10**10 L_sun]
            else:
                # Convert the UV fluxes to rest-frame luminosity in erg/s/Hz. This
                # luminosity can be converted into a SFR using, e.g., Kennicutt+98,
                # SFR=1.4e-28 * L_UV
                cflux *= cwave**2 / (C_LIGHT * 1e13) / norm # [monochromatic luminosity in 10**(-28) erg/s/Hz]
            if cflux > 0:
                lums[label] = np.log10(cflux) # * u.erg/(u.second*u.Hz)

        cwaves = [3728.483, 4862.683, 5008.239, 6564.613]
        labels = ['FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT']
        for cwave, label in zip(cwaves, labels):
            J = (self.sspwave > cwave-500) * (self.sspwave < cwave+500)
            I = (self.sspwave[J] > cwave-20) * (self.sspwave[J] < cwave+20)
            smooth = median_filter(continuum[J], 200)
            clipflux, _, _ = sigmaclip(smooth[I], low=1.5, high=3)
            cflux = np.median(clipflux) # [flux in 10**-17 erg/s/cm2/A]
            cfluxes[label] = cflux # * u.erg/(u.second*u.cm**2*u.Angstrom)
            
            #import matplotlib.pyplot as plt
            #print(cwave, cflux)
            #plt.clf()
            #plt.plot(self.sspwave[J], continuum[J])
            #plt.plot(self.sspwave[J], smooth, color='k')
            #plt.axhline(y=cflux, color='red')
            #plt.axvline(x=cwave, color='red')
            #plt.xlim(cwave - 50, cwave + 50)
            #plt.savefig('desi-users/ioannis/tmp/junk.png')
            #pdb.set_trace()

        return kcorr, absmag, ivarabsmag, bestmaggies, mstar, lums, cfluxes

    def _call_nnls(self, modelflux, flux, ivar, xparam=None, debug=False,
                   interpolate_coeff=False, xlabel=None):
        """Wrapper on nnls.

        Works with both spectroscopic and photometric input and with both 2D and
        3D model spectra.

        To be documented.

        interpolate_coeff - return the interpolated coefficients when exploring
          an array or grid of xparam

        """
        from scipy.optimize import nnls
        from fastspecfit.util import find_minima, minfit
        
        if xparam is not None:
            nn = len(xparam)

        inverr = np.sqrt(ivar)
        bvector = flux * inverr

        # If xparam is None (equivalent to modelflux having just two
        # dimensions, [npix,nage]), assume we are just finding the
        # coefficients at some best-fitting value...
        if xparam is None:
            Amatrix = modelflux * inverr[:, np.newaxis]
            try:
                coeff, rnorm = nnls(A=Amatrix, b=bvector)
            except RuntimeError:
                coeff, _ = nnls(A=Amatrix, b=bvector, maxiter=Amatrix.shape[1] * 100)                
            chi2 = np.sum(ivar * (flux - modelflux.dot(coeff))**2)
            return coeff, chi2

        # ...otherwise iterate over the xparam (e.g., AV or vdisp) dimension.
        Amatrix = modelflux * inverr[:, np.newaxis, np.newaxis] # reshape into [npix/nband,nage,nAV/nvdisp]
        coeff, chi2grid = [], []
        for ii in np.arange(nn):
            _coeff, _ = nnls(A=Amatrix[:, :, ii], b=bvector)
            chi2 = np.sum(ivar * (flux - modelflux[:, :, ii].dot(_coeff))**2)
            coeff.append(_coeff)
            chi2grid.append(chi2)
        coeff = np.array(coeff)
        chi2grid = np.array(chi2grid)
        
        try:
            imin = find_minima(chi2grid)[0]
            xbest, xerr, chi2min, warn = minfit(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2])
        except:
            errmsg = 'A problem was encountered minimizing chi2.'
            self.log.warning(errmsg)
            imin, xbest, xerr, chi2min, warn = 0, 0.0, 0.0, 0.0, 1

        if warn == 0:
            xivar = 1.0 / xerr**2
        else:
            chi2min = 0.0
            xivar = 0.0
            xbest = xparam[0]

        # optionally interpolate the coefficients
        if interpolate_coeff:
            from scipy.interpolate import interp1d
            if xbest == xparam[0]:
                bestcoeff = coeff[0, :]
            else:
                xindx = np.arange(len(xparam))
                f = interp1d(xindx, coeff, axis=0)
                bestcoeff = f(np.interp(xbest, xparam, xindx))
        else:
            bestcoeff = None

            # interpolate the coefficients
            #np.interp(xbest, xparam, np.arange(len(xparam)))            

        if debug:
            if xivar > 0:
                leg = r'${:.3f}\pm{:.3f}\ (\chi^2_{{min}}={:.2f})$'.format(xbest, 1/np.sqrt(xivar), chi2min)
            else:
                leg = r'${:.3f}$'.format(xbest)
                
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.scatter(xparam, chi2grid)
            ax.scatter(xparam[imin-1:imin+2], chi2grid[imin-1:imin+2], color='red')
            #ax.set_ylim(chi2min*0.95, np.max(chi2grid[imin-1:imin+2])*1.05)
            #ax.plot(xx, np.polyval([aa, bb, cc], xx), ls='--')
            ax.axvline(x=xbest, color='k')
            if xivar > 0:
                ax.axhline(y=chi2min, color='k')
            #ax.set_yscale('log')
            #ax.set_ylim(chi2min, 63.3)
            if xlabel:
                ax.set_xlabel(xlabel)
                #ax.text(0.03, 0.9, '{}={}'.format(xlabel, leg), ha='left',
                #        va='center', transform=ax.transAxes)
            ax.text(0.03, 0.9, leg, ha='left', va='center', transform=ax.transAxes)
            ax.set_ylabel(r'$\chi^2$')
            fig.savefig('qa-chi2min.png')
            #fig.savefig('desi-users/ioannis/tmp/qa-chi2min.png')

        return chi2min, xbest, xivar, bestcoeff

    def continuum_fastphot(self, data, minbands=3):
        """Fit the broad photometry.

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            `unpack_one_spectrum`.
        minbands : :class:`int`, defaults to 3
            Minimum number of photometric bands required to attempt a fit.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in `init_phot_output`.

        """
        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_phot_output()

        redshift = data['zredrock']

        # Prepare the reddened and unreddened SSP templates by redshifting and
        # normalizing. Note that we ignore templates which are older than the
        # age of the universe at the galaxy redshift.
        if self.constrain_age:
            agekeep = self.younger_than_universe(redshift)
        else:
            agekeep = np.arange(self.nage)
            
        t0 = time.time()
        zsspflux_dustnomvdisp, zsspphot_dustnomvdisp = self.SSP2data(
            self.sspflux_dustnomvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=None, specres=None,
            south=data['photsys'] == 'S')
        self.log.info('Preparing the models took {:.2f} seconds.'.format(time.time()-t0))
        
        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = (data['phot']['flam_ivar'].data / self.fluxnorm**2) * self.bands_to_fit
        assert(np.all(objflamivar >= 0))

        # some secondary targets have no good photometry
        if np.all(objflamivar == 0) or (np.sum(objflamivar > 0) < minbands):
            self.log.warning('All photometry is masked or number of good photometric bands {}<{}'.format(
                np.sum(objflamivar > 0), minbands))
            AVbest, AVivar = self.AV_nominal, 0.0
            nage = self.nage
            coeff = np.zeros(self.nage)
            continuummodel = np.zeros(len(self.sspwave))
        else:
            zsspflam_dustnomvdisp = zsspphot_dustnomvdisp['flam'].data * self.fluxnorm * self.massnorm # [nband,nage*nAV]

            inodust = np.ndarray.item(np.where(self.AV == 0)[0]) # should always be index 0

            npix, nmodel = zsspflux_dustnomvdisp.shape
            nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)

            zsspflam_dustnomvdisp = zsspflam_dustnomvdisp.reshape(len(self.bands), nage, self.nAV) # [nband,nage,nAV]

            t0 = time.time()
            AVchi2min, AVbest, AVivar, _ = self._call_nnls(
                zsspflam_dustnomvdisp, objflam, objflamivar, xparam=self.AV,
                debug=False)
            self.log.info('Fitting the photometry took: {:.2f} seconds.'.format(time.time()-t0))
            if AVivar > 0:
                # protect against tiny ivar from becomming infinite in the output table
                if AVivar < 1e-3:
                    self.log.warning('AV inverse variance is tiny; capping at 1e-3')
                    AVivar = 1e-3
                self.log.info('Best-fitting photometric A(V)={:.4f}+/-{:.4f}'.format(
                    AVbest, 1/np.sqrt(AVivar)))
            else:
                AVbest = self.AV_nominal
                self.log.info('Finding photometric A(V) failed; adopting A(V)={:.4f}'.format(self.AV_nominal))

            # Get the final set of coefficients and chi2 at the best-fitting
            # reddening and nominal velocity dispersion.
            bestsspflux, bestphot = self.SSP2data(self.sspflux_dustnomvdisp[:, agekeep, inodust],
                                                  self.sspwave, AV=AVbest, redshift=redshift,
                                                  south=data['photsys'] == 'S', debug=False)

            coeff, _chi2min = self._call_nnls(bestphot['flam'].data*self.massnorm*self.fluxnorm,
                                              objflam, objflamivar) # bestphot['flam'] is [nband, nage]

            continuummodel = bestsspflux.dot(coeff)
            
            # For some reason I don't understand, the chi2 returned by
            # _call_nnls, above, can be unreasonably small (like 10**-25 or
            # smaller), so we synthesize the photometry one last time with the
            # final coefficients to get the correct chi2 estimate.
            _, synthmodelphot = self.SSP2data(self.sspflux_dustnomvdisp[:, agekeep, inodust],
                                              self.sspwave, redshift=redshift, synthphot=True,
                                              AV=AVbest, coeff=coeff * self.massnorm)
            
            chi2min = np.sum(objflamivar * (objflam - self.fluxnorm * synthmodelphot['flam'])**2)

            dof = np.sum(objflamivar > 0) - 1 # 1 free parameter??
            chi2min /= dof

        # Compute DN4000, K-corrections, and rest-frame quantities.
        if np.count_nonzero(coeff > 0) == 0:
            self.log.warning('Continuum coefficients are all zero or data not fit.')
            chi2min, dn4000, meanage, mstar = 0.0, 0.0, 0.0, 0.0 # -1.0, -1.0, -1.0
            kcorr = np.zeros(len(self.absmag_bands))
            absmag = np.zeros(len(self.absmag_bands))#-99.0
            ivarabsmag = np.zeros(len(self.absmag_bands))
            synth_bestmaggies = np.zeros(len(self.bands))
            lums, cfluxes = {}, {}
        else:
            dn4000, _ = self.get_dn4000(self.sspwave, continuummodel, rest=True)
            meanage = self.get_meanage(coeff)
            kcorr, absmag, ivarabsmag, synth_bestmaggies, mstar, lums, cfluxes = self.kcorr_and_absmag(data, continuummodel, coeff)

            # convert the synthesized photometry to nanomaggies
            synth_bestmaggies *= 1e9

            self.log.info('Photometric DN(4000)={:.3f}, Age={:.2f} Gyr, Mr={:.2f} mag, Mstar={:.4g}'.format(
                dn4000, meanage, absmag[1], mstar))

        # Pack it up and return.
        result['CONTINUUM_COEFF'][0][:nage] = coeff
        result['CONTINUUM_RCHI2'][0] = chi2min
        result['CONTINUUM_AGE'][0] = meanage # * u.Gyr
        result['CONTINUUM_AV'][0] = AVbest # * u.mag
        result['CONTINUUM_AV_IVAR'][0] = AVivar # / (u.mag**2)
        result['DN4000_MODEL'][0] = dn4000
        #if False:
        #    for iband, band in enumerate(self.fiber_bands):
        #        result['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
        #    for iband, band in enumerate(self.bands):
        #        result['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        #        result['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        for iband, band in enumerate(self.absmag_bands):
            result['KCORR_{}'.format(band.upper())] = kcorr[iband] # * u.mag
            result['ABSMAG_{}'.format(band.upper())] = absmag[iband] # * u.mag
            result['ABSMAG_IVAR_{}'.format(band.upper())] = ivarabsmag[iband] # / (u.mag**2)
        for iband, band in enumerate(self.bands):
            result['FLUX_SYNTH_MODEL_{}'.format(band.upper())] = synth_bestmaggies[iband] # * u.nanomaggy

        if mstar > 0:
            result['LOGMSTAR'] = np.log10(mstar) # * u.solMass
        if bool(lums):
            for lum in lums.keys():
                result[lum] = lums[lum]
        if bool(cfluxes):
            for cflux in cfluxes.keys():
                result[cflux] = cfluxes[cflux]
                
        return result, continuummodel
    
    def continuum_specfit(self, data):
        """Fit the non-negative stellar continuum of a single spectrum.

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            :func:`fastspecfit.io.DESISpectra.read_and_unpack`.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in :func:`self.init_spec_output`.

        Notes
        -----
          - Consider using cross-correlation to update the redrock redshift.
          - We solve for velocity dispersion if solve_vdisp=True or ((SNR_B>3 or
            SNR_R>3) and REDSHIFT<1).

        """
        result = self.init_spec_output()

        redshift = data['zredrock']
        result['CONTINUUM_Z'] = redshift
        for icam, cam in enumerate(data['cameras']):
            result['CONTINUUM_SNR_{}'.format(cam.upper())] = data['snr'][icam]

        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        specwave = np.hstack(data['wave'])
        specflux = np.hstack(data['flux'])
        specivar = np.hstack(data['ivar']) * np.logical_not(np.hstack(data['linemask'])) # mask emission lines
        if np.all(specivar == 0) or np.any(specivar < 0):
            specivar = np.hstack(data['ivar']) # not great...
            if np.all(specivar == 0) or np.any(specivar < 0):
                errmsg = 'All pixels are masked or some inverse variances are negative!'
                self.log.critical(errmsg)
                raise ValueError(errmsg)

        # dev - Add the photometry.
        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = (data['phot']['flam_ivar'].data / self.fluxnorm**2) * self.bands_to_fit
        assert(np.all(objflamivar >= 0))

        # Prepare the reddened and unreddened SSP templates by redshifting and
        # normalizing. Note that we ignore templates which are older than the
        # age of the universe at the redshift of the object.
        if self.constrain_age:
            agekeep = self.younger_than_universe(redshift)
        else:
            agekeep = np.arange(self.nage)

        t0 = time.time()
        zsspflux_dustnomvdisp, _ = self.SSP2data(
            self.sspflux_dustnomvdisp[:, agekeep, :], self.sspwave, # [npix,nage,nAV]
            redshift=redshift, specwave=data['wave'], specres=data['res'],
            cameras=data['cameras'], synthphot=False)

        zsspflux_dustnomvdisp = np.concatenate(zsspflux_dustnomvdisp, axis=0)  # [npix,nage*nAV]
        npix, nmodel = zsspflux_dustnomvdisp.shape
        nage = nmodel // self.nAV # accounts for age-of-the-universe constraint (!=self.nage)
        zsspflux_dustnomvdisp = zsspflux_dustnomvdisp.reshape(npix, nage, self.nAV)       # [npix,nage,nAV]
        self.log.info('Preparing the models took {:.2f} seconds.'.format(time.time()-t0))

        # Fit the spectra for reddening using the models convolved to the
        # nominal velocity dispersion.
        t0 = time.time()
        AVchi2min, AVbest, AVivar, AVcoeff = self._call_nnls(
            zsspflux_dustnomvdisp, specflux, specivar, xparam=self.AV,
            debug=False, interpolate_coeff=self.solve_vdisp,
            xlabel=r'$A_V$ (mag)')
        self.log.info('Fitting for the reddening took: {:.2f} seconds.'.format(time.time()-t0))
        if AVivar > 0:
            # protect against tiny ivar from becomming infinite in the output table
            if AVivar < 1e-3:
                self.log.warning('AV inverse variance is tiny; capping at 1e-3')
                AVivar = 1e-3
            self.log.info('Best-fitting spectroscopic A(V)={:.4f}+/-{:.4f}'.format(
                AVbest, 1/np.sqrt(AVivar)))
        else:
            AVbest = self.AV_nominal
            self.log.info('Finding spectroscopic A(V) failed; adopting A(V)={:.4f}'.format(self.AV_nominal))

        # If the S/N is good enough, solve for the velocity dispersion.
        compute_vdisp = ((result['CONTINUUM_SNR_B'] > 3) and (result['CONTINUUM_SNR_R'] > 3)) and (redshift < 1.0)
        if compute_vdisp:
            self.log.info('Solving for velocity dispersion: S/N_B={:.2f}, S/N_R={:.2f}, redshift={:.3f}'.format(
                result['CONTINUUM_SNR_B'][0], result['CONTINUUM_SNR_R'][0], redshift))
            
        if self.solve_vdisp or compute_vdisp:
            t0 = time.time()
            if True:
                zsspflux_vdisp = []
                for vdisp in self.vdisp:
                    _zsspflux_vdisp, _ = self.SSP2data(self.sspflux[:, agekeep], self.sspwave,
                                                       specwave=data['wave'], specres=data['res'],
                                                       AV=AVbest, vdisp=vdisp, redshift=redshift,
                                                       cameras=data['cameras'], synthphot=False)
                    _zsspflux_vdisp = np.concatenate(_zsspflux_vdisp, axis=0)
                    zsspflux_vdisp.append(_zsspflux_vdisp)
                zsspflux_vdisp = np.stack(zsspflux_vdisp, axis=-1) # [npix,nage,nvdisp] at best A(V)
                
                vdispchi2min, vdispbest, vdispivar, _ = self._call_nnls(
                    zsspflux_vdisp, specflux, specivar, xparam=self.vdisp,
                    xlabel=r'$\sigma$ (km/s)', debug=False)
            else:
                # The code below does a refinement of the velocity dispersion around
                # the "best" vdisp based on a very quick-and-dirty chi2 estimation
                # (based on the AVcoeff coefficients found above; no refitting of
                # the coefficients). However, the refined value is no different than
                # the one found using the coarse grid and the uncertainty in the
                # velocity dispersion is ridiculously large, which I don't fully 
                # understand.
                # /global/u2/i/ioannis/code/desihub/fastspecfit-projects/pv-vdisp/fastspecfit-pv-vdisp --targetids 39627665157658710
                zsspflux_vdispnomdust, _ = self.SSP2data(
                    self.sspflux_vdispnomdust[:, agekeep, :], self.sspwave, # [npix,nage,nvdisp]
                    redshift=redshift, specwave=data['wave'], specres=data['res'],
                    cameras=data['cameras'], synthphot=False)
                zsspflux_vdispnomdust = np.concatenate(zsspflux_vdispnomdust, axis=0)  # [npix,nmodel=nage*nvdisp]
                npix, nmodel = zsspflux_vdispnomdust.shape
                nage = nmodel // self.nvdisp # accounts for age-of-the-universe constraint (!=self.nage)
                zsspflux_vdispnomdust = zsspflux_vdispnomdust.reshape(npix, nage, self.nvdisp) # [npix,nage,nvdisp]
    
                # This refits for the coefficients, so it's slower than the "quick" chi2 minimization.
                #vdispchi2min, vdispbest, vdispivar = self._call_nnls(
                #    zsspflux_vdispnomdust, specflux, specivar, xparam=self.vdisp, debug=False)
    
                # Do a quick chi2 minimization over velocity dispersion using the
                # coefficients from the reddening modeling (see equations 7-9 in
                # Benitez+2000). Should really be using broadcasting...
                vdispchi2 = np.zeros(self.nvdisp)
                for iv in np.arange(self.nvdisp):
                    vdispchi2[iv] = np.sum(specivar * (specflux - zsspflux_vdispnomdust[:, :, iv].dot(AVcoeff))**2)
    
                vmindx = np.argmin(vdispchi2)
                if vmindx == 0 or vmindx == self.nvdisp-1: # on the edge; no minimum
                    self.log.info('Finding vdisp failed; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))
                    vdispbest, vdispivar = self.vdisp_nominal, 0.0
                else:
                    # Do a more refined search with +/-XX km/s around the initial minimum.
                    #vdispinit = self.vdisp[vmindx]
                    vdispfine = np.linspace(self.vdisp[vmindx]-10, self.vdisp[vmindx]+10, 15)
                    #vdispfine = np.linspace(self.vdisp[vmindx-1], self.vdisp[vmindx+1], 15)
                    #vdispmin, vdispmax, dvdisp = vdispinit - 10.0, vdispinit + 10.0, 0.01
                    #if vdispmin < 50:
                    #    vdispmin = 50.0
                    #if vdispmax > 400:
                    #    vdispmax = 400
                    #vdispfine = np.arange(vdispmin, vdispmax, dvdisp)
                    nvdispfine = len(vdispfine)
    
                    #atten = self.dust_attenuation(self.sspwave, AVbest)
                    sspflux_vdispfine = []
                    for vdisp in vdispfine:
                        sspflux_vdispfine.append(self.convolve_vdisp(self.sspflux[:, agekeep], vdisp))
                    sspflux_vdispfine = np.stack(sspflux_vdispfine, axis=-1) # [npix,nage,nvdisp]
                    
                    zsspflux_vdispfine, _ = self.SSP2data(sspflux_vdispfine, self.sspwave,
                        redshift=redshift, specwave=data['wave'], specres=data['res'],
                        cameras=data['cameras'], AV=AVbest, synthphot=False)
                    zsspflux_vdispfine = np.concatenate(zsspflux_vdispfine, axis=0)  # [npix,nmodel=nage*nvdisp]
                    npix, nmodel = zsspflux_vdispfine.shape
                    nage = nmodel // nvdispfine # accounts for age-of-the-universe constraint (!=self.nage)
                    zsspflux_vdispfine = zsspflux_vdispfine.reshape(npix, nage, nvdispfine) # [npix,nage,nvdisp]
                    
                    vdispchi2min, vdispbest, vdispivar, _ = self._call_nnls(
                        zsspflux_vdispfine, specflux, specivar, xparam=vdispfine,
                        interpolate_coeff=False, debug=False)
                
            self.log.info('Fitting for the velocity dispersion took: {:.2f} seconds.'.format(time.time()-t0))
            if vdispivar > 0:
                # protect against tiny ivar from becomming infinite in the output table
                if vdispivar < 1e-5:
                    self.log.warning('vdisp inverse variance is tiny; capping at 1e-5')
                    vdispivar = 1e-5
                self.log.info('Best-fitting vdisp={:.2f}+/-{:.2f} km/s'.format(
                    vdispbest, 1/np.sqrt(vdispivar)))
            else:
                vdispbest = self.vdisp_nominal
                self.log.info('Finding vdisp failed; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))
        else:
            self.log.info('Insufficient S/N to compute vdisp; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))
            vdispbest, vdispivar = self.vdisp_nominal, 0.0

        # Get the final set of coefficients and chi2 at the best-fitting
        # reddening and velocity dispersion.
        bestsspflux, _ = self.SSP2data(self.sspflux[:, agekeep], self.sspwave, redshift=redshift,
                                       specwave=data['wave'], specres=data['res'],
                                       AV=AVbest, vdisp=vdispbest, cameras=data['cameras'],
                                       south=data['photsys'] == 'S', synthphot=False)
        bestsspflux = np.concatenate(bestsspflux, axis=0)
        coeff, chi2min = self._call_nnls(bestsspflux, specflux, specivar)
        chi2min /= np.sum(specivar > 0) # dof???


        # Get the light-weighted age and DN(4000).
        bestfit = bestsspflux.dot(coeff)
        meanage = self.get_meanage(coeff, agekeep)

        flam_ivar = np.hstack(data['ivar']) # specivar is line-masked!
        dn4000, dn4000_ivar = self.get_dn4000(specwave, specflux, flam_ivar=flam_ivar, 
                                              redshift=redshift, rest=False)
        dn4000_model, _ = self.get_dn4000(specwave, bestfit, redshift=redshift, rest=False)
        
        if False:
            print(dn4000, dn4000_model, 1/np.sqrt(dn4000_ivar))
            from fastspecfit.util import ivar2var            
            import matplotlib.pyplot as plt

            restwave = specwave / (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
            fnu = specflux * flam2fnu # [erg/s/cm2/Hz]
            fnu_model = bestfit * flam2fnu
            fnu_ivar = flam_ivar / flam2fnu**2            
            fnu_sigma, fnu_mask = ivar2var(fnu_ivar, sigma=True)

            #from astropy.table import Table
            #out = Table()
            #out['WAVE'] = restwave
            #out['FLUX'] = specflux
            #out['IVAR'] = flam_ivar
            #out.write('desi-39633089965589981.fits', overwrite=True)
            
            I = (restwave > 3835) * (restwave < 4115)
            J = (restwave > 3835) * (restwave < 4115) * fnu_mask

            fig, ax = plt.subplots()
            ax.fill_between(restwave[I], fnu[I]-fnu_sigma[I], fnu[I]+fnu_sigma[I],
                            label='Data Dn(4000)={:.3f}+/-{:.3f}'.format(dn4000, 1/np.sqrt(dn4000_ivar)))
            ax.plot(restwave[I], fnu_model[I], color='red', label='Model Dn(4000)={:.3f}'.format(dn4000_model))
            ylim = ax.get_ylim()
            ax.fill_between([3850, 3950], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.fill_between([4000, 4100], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.set_xlabel(r'Rest Wavelength ($\AA$)')
            ax.set_ylabel(r'$F_{\nu}$ (erg/s/cm2/Hz)')
            ax.legend()
            fig.savefig('desi-users/ioannis/tmp/qa-dn4000.png')

        if dn4000_ivar > 0:
            self.log.info('Spectroscopic DN(4000)={:.3f}+/-{:.3f}, Age={:.2f} Gyr'.format(dn4000, 1/np.sqrt(dn4000_ivar), meanage))
        else:
            self.log.info('Spectroscopic DN(4000)={:.3f}, Age={:.2f} Gyr'.format(dn4000, meanage))

        png = None
        #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/smooth-continuum.png'
        #linemask = np.hstack(data['linemask_all'])
        linemask = np.hstack(data['linemask'])
        if np.all(coeff == 0):
            _smooth_continuum = np.zeros_like(bestfit)
        else:
            _smooth_continuum, _ = self.smooth_continuum(specwave, specflux - bestfit, specivar,
                                                         redshift, linemask=linemask, png=png)

        # Unpack the continuum into individual cameras.
        continuummodel, smooth_continuum = [], []
        for camerapix in data['camerapix']:
            continuummodel.append(bestfit[camerapix[0]:camerapix[1]])
            smooth_continuum.append(_smooth_continuum[camerapix[0]:camerapix[1]])

        ## Like above, but with per-camera smoothing.
        #smooth_continuum = self.smooth_residuals(
        #    continuummodel, data['wave'], data['flux'],
        #    data['ivar'], data['linemask'], percamera=False)

        if False:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1)
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                resid = data['flux'][icam]-continuummodel[icam]
                ax[0].plot(data['wave'][icam], resid)
                ax[1].plot(data['wave'][icam], resid-smooth_continuum[icam])
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                resid = data['flux'][icam]-continuummodel[icam]
                pix_emlines = np.logical_not(data['linemask'][icam]) # affected by line = True
                ax[0].scatter(data['wave'][icam][pix_emlines], resid[pix_emlines], s=30, color='red')
                ax[0].plot(data['wave'][icam], smooth_continuum[icam], color='k', alpha=0.7, lw=2)
            plt.savefig('junk.png')

        # Pack it in and return.
        result['CONTINUUM_COEFF'][0][0:nage] = coeff
        result['CONTINUUM_RCHI2'][0] = chi2min
        result['CONTINUUM_AV'][0] = AVbest # * u.mag
        result['CONTINUUM_AV_IVAR'][0] = AVivar # / (u.mag**2)
        result['CONTINUUM_VDISP'][0] = vdispbest # * u.kilometer/u.second
        result['CONTINUUM_VDISP_IVAR'][0] = vdispivar # * (u.second/u.kilometer)**2
        result['CONTINUUM_AGE'] = meanage # * u.Gyr
        result['DN4000'][0] = dn4000
        result['DN4000_IVAR'][0] = dn4000_ivar
        result['DN4000_MODEL'][0] = dn4000_model

        for icam, cam in enumerate(data['cameras']):
            nonzero = continuummodel[icam] != 0
            #nonzero = np.abs(continuummodel[icam]) > 1e-5
            if np.sum(nonzero) > 0:
                corr = np.median(smooth_continuum[icam][nonzero] / continuummodel[icam][nonzero])
                result['CONTINUUM_SMOOTHCORR_{}'.format(cam.upper())] = corr * 100 # [%]

        self.log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
            result['CONTINUUM_SMOOTHCORR_B'][0], result['CONTINUUM_SMOOTHCORR_R'][0],
            result['CONTINUUM_SMOOTHCORR_Z'][0]))
        
        return result, continuummodel, smooth_continuum
    
    def continuum_joint_specfit(self, data):
        """Fit the non-negative stellar continuum of a single spectrum.

        Parameters
        ----------
        data : :class:`dict`
            Dictionary of input spectroscopy (plus ancillary data) populated by
            :func:`fastspecfit.io.DESISpectra.read_and_unpack`.

        Returns
        -------
        :class:`astropy.table.Table`
            Table with all the continuum-fitting results with columns documented
            in :func:`self.init_spec_output`.

        Notes
        -----
          - Consider using cross-correlation to update the redrock redshift.
          - We solve for velocity dispersion if solve_vdisp=True or ((SNR_B>3 or
            SNR_R>3) and REDSHIFT<1).

        """
        result = self.init_spec_output()

        redshift = data['zredrock']
        result['CONTINUUM_Z'] = redshift
        for icam, cam in enumerate(data['cameras']):
            result['CONTINUUM_SNR_{}'.format(cam.upper())] = data['snr'][icam]

        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        specwave = np.hstack(data['wave'])
        specflux = np.hstack(data['flux'])
        specivar = np.hstack(data['ivar']) * np.logical_not(np.hstack(data['linemask'])) # mask emission lines
        if np.all(specivar == 0) or np.any(specivar < 0):
            specivar = np.hstack(data['ivar']) # not great...
            if np.all(specivar == 0) or np.any(specivar < 0):
                errmsg = 'All pixels are masked or some inverse variances are negative!'
                self.log.critical(errmsg)
                raise ValueError(errmsg)

        # dev - Add the photometry.
        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = (data['phot']['flam_ivar'].data / self.fluxnorm**2) * self.bands_to_fit
        assert(np.all(objflamivar >= 0))

        # dev - aperture correction based on the r-band synthesized and observed
        # photometry. Note that the right thing to do is to do a quick-and-dirty
        # fit of the spectrum and used the model-synthesized flux (ignoring any
        # smooth continuum correction), which will guarantee that the aperture
        # correction is always positive.
        apcorr = data['phot']['nanomaggies'][1] / data['synthphot']['nanomaggies'][1]
        data['apcorr'] = apcorr

        # Prepare the reddened and unreddened SSP templates by redshifting and
        # normalizing. Note that we ignore templates which are older than the
        # age of the universe at the redshift of the object.
        if self.constrain_age:
            agekeep = self.younger_than_universe(redshift)
        else:
            agekeep = np.arange(self.nsed)

        nage = len(agekeep)

        # Prepare the spectral and photometric models at the galaxy
        # redshift. And if the wavelength coverage is sufficient, also solve for
        # the velocity dispersion.

        #compute_vdisp = ((result['CONTINUUM_SNR_B'] > 1) and (result['CONTINUUM_SNR_R'] > 1)) and (redshift < 1.0)
        restwave = specwave / (1+redshift)
        Ivdisp = np.where((specivar > 0) * (restwave > 3500.0) * (restwave < 5500.0))[0]
        compute_vdisp = (len(Ivdisp) > 0) * (np.ptp(restwave[Ivdisp]) > 500.0)

        self.log.info('S/N_B={:.2f}, S/N_R={:.2f}, rest wavelength coverage={:.0f}-{:.0f} A.'.format(
            result['CONTINUUM_SNR_B'][0], result['CONTINUUM_SNR_R'][0], restwave[0], restwave[-1]))

        if self.solve_vdisp or compute_vdisp:
            
            # Only use models which likely affect the velocity dispersion.
            Mvdisp = np.where((self.sspinfo[agekeep]['fagn'] == 0) * (self.sspinfo[agekeep]['av'] == 0.0))[0]
            nMvdisp = len(Mvdisp)

            t0 = time.time()
            zsspflux_vdisp, _ = self.SSP2data(
                self.sspflux_vdisp[:, agekeep[Mvdisp], :], self.sspwave, # [npix,nsed,nvdisp]
                redshift=redshift, specwave=data['wave'], specres=data['res'],
                cameras=data['cameras'], synthphot=False)#, south=data['photsys'] == 'S')
            self.log.info('Preparing {}/{} models took {:.2f} seconds.'.format(
                nMvdisp, nage, time.time()-t0))

            zsspflux_vdisp = np.concatenate(zsspflux_vdisp, axis=0)  # [npix,nMvdisp*nvdisp]
            npix, nmodel = zsspflux_vdisp.shape
            zsspflux_vdisp = zsspflux_vdisp.reshape(npix, nMvdisp, self.nvdisp) # [npix,nMvdisp,nvdisp]

            t0 = time.time()
            vdispchi2min, vdispbest, vdispivar, _ = self._call_nnls(
                zsspflux_vdisp[Ivdisp, :, :], specflux[Ivdisp], specivar[Ivdisp],
                xparam=self.vdisp, xlabel=r'$\sigma$ (km/s)', debug=False)#True)
            self.log.info('Fitting for the velocity dispersion took: {:.2f} seconds.'.format(time.time()-t0))

            if vdispivar > 0:
                # protect against tiny ivar from becomming infinite in the output table
                if vdispivar < 1e-5:
                    self.log.warning('vdisp inverse variance is tiny; capping at 1e-5')
                    vdispivar = 1e-5
                self.log.info('Best-fitting vdisp={:.2f}+/-{:.2f} km/s'.format(
                    vdispbest, 1/np.sqrt(vdispivar)))
            else:
                vdispbest = self.vdisp_nominal
                self.log.info('Finding vdisp failed; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))
        else:
            # Fit at the nominal velocity dispersion.
            t0 = time.time()
            zsspflux_novdisp, zsspphot_novdisp = self.SSP2data(
                self.sspflux_vdisp[:, agekeep, self.vdisp_nominal_indx], self.sspwave, # [npix,nsed,nvdisp]
                redshift=redshift, specwave=data['wave'], specres=data['res'],
                cameras=data['cameras'], synthphot=True, south=data['photsys'] == 'S')
            self.log.info('Preparing {} models took {:.2f} seconds.'.format(nage, time.time()-t0))

            zsspflux_vdisp = np.concatenate(zsspflux_vdisp, axis=0)  # [npix,nage*nvdisp]
            npix, nmodel = zsspflux_vdisp.shape
            zsspflux_vdisp = zsspflux_vdisp.reshape(npix, nage, self.nvdisp) # [npix,nage,nvdisp]

            zsspflam_vdisp = zsspphot_vdisp['flam'].data * self.fluxnorm * self.massnorm # [nband,nage*nvdisp]
            zsspflam_vdisp = zsspflam_vdisp.reshape(len(self.bands), nage, self.nvdisp)  # [nband,nage,nvdisp]

            if compute_vdisp:
                self.log.info('Sufficient wavelength covereage to compute vdisp but solve_vdisp=False; adopting nominal vdisp={:.2f} km/s'.format(
                    self.vdisp_nominal))
            else:
                self.log.info('Insufficient wavelength covereage to compute vdisp; adopting nominal vdisp={:.2f} km/s'.format(
                    self.vdisp_nominal))
                
            vdispbest, vdispivar = self.vdisp_nominal, 0.0

        # Get the final set of coefficients and chi2 at the best-fitting velocity dispersion. 
        bestsspflux, bestphot = self.SSP2data(self.sspflux[:, agekeep], self.sspwave, redshift=redshift,
                                              specwave=data['wave'], specres=data['res'],
                                              vdisp=vdispbest, cameras=data['cameras'],
                                              south=data['photsys'] == 'S', synthphot=True)
        bestsspflux = np.concatenate(bestsspflux, axis=0)
        bestflam = bestphot['flam'].data*self.massnorm*self.fluxnorm

        coeff, chi2min = self._call_nnls(np.vstack((bestflam, bestsspflux)),
                                         np.hstack((objflam, specflux*apcorr)),
                                         np.hstack((objflamivar, specivar/apcorr**2)))
        chi2min /= (np.sum(objflamivar > 0) + np.sum(specivar > 0)) # dof???

        # Get the light-weighted physical properties and DN(4000).
        bestfit = bestsspflux.dot(coeff)

        AV = self.get_mean_property('av', coeff, agekeep)                        # [mag]
        age = self.get_mean_property('age', coeff, agekeep, normalization=1e9)   # [Gyr]
        zzsun = self.get_mean_property('zzsun', coeff, agekeep, log10=False)     # [log Zsun]
        fagn = self.get_mean_property('fagn', coeff, agekeep)
        logmstar = self.get_mean_property('mstar', coeff, agekeep, normalization=1/self.massnorm, log10=True) # [Msun]
        sfr = self.get_mean_property('sfr', coeff, agekeep, normalization=1/self.massnorm, log10=False)       # [Msun/yr]

        flam_ivar = np.hstack(data['ivar']) # specivar is line-masked!
        dn4000, dn4000_ivar = self.get_dn4000(specwave, specflux, flam_ivar=flam_ivar, 
                                              redshift=redshift, rest=False)
        dn4000_model, _ = self.get_dn4000(specwave, bestfit, redshift=redshift, rest=False)
        
        if False:
            print(dn4000, dn4000_model, 1/np.sqrt(dn4000_ivar))
            from fastspecfit.util import ivar2var            
            import matplotlib.pyplot as plt

            restwave = specwave / (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
            fnu = specflux * flam2fnu # [erg/s/cm2/Hz]
            fnu_model = bestfit * flam2fnu
            fnu_ivar = flam_ivar / flam2fnu**2            
            fnu_sigma, fnu_mask = ivar2var(fnu_ivar, sigma=True)

            #from astropy.table import Table
            #out = Table()
            #out['WAVE'] = restwave
            #out['FLUX'] = specflux
            #out['IVAR'] = flam_ivar
            #out.write('desi-39633089965589981.fits', overwrite=True)
            
            I = (restwave > 3835) * (restwave < 4115)
            J = (restwave > 3835) * (restwave < 4115) * fnu_mask

            fig, ax = plt.subplots()
            ax.fill_between(restwave[I], fnu[I]-fnu_sigma[I], fnu[I]+fnu_sigma[I],
                            label='Data Dn(4000)={:.3f}+/-{:.3f}'.format(dn4000, 1/np.sqrt(dn4000_ivar)))
            ax.plot(restwave[I], fnu_model[I], color='red', label='Model Dn(4000)={:.3f}'.format(dn4000_model))
            ylim = ax.get_ylim()
            ax.fill_between([3850, 3950], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.fill_between([4000, 4100], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.set_xlabel(r'Rest Wavelength ($\AA$)')
            ax.set_ylabel(r'$F_{\nu}$ (erg/s/cm2/Hz)')
            ax.legend()
            fig.savefig('desi-users/ioannis/tmp/qa-dn4000.png')

        if dn4000_ivar > 0:
            self.log.info('Spectroscopic DN(4000)={:.3f}+/-{:.3f}, Age={:.2f} Gyr, Mstar={:.4g}'.format(
                dn4000, 1/np.sqrt(dn4000_ivar), age, logmstar))
        else:
            self.log.info('Spectroscopic DN(4000)={:.3f}, Age={:.2f} Gyr, Mstar={:.4g}'.format(
                dn4000, age, logmstar))

        png = None
        #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/smooth-continuum.png'
        #linemask = np.hstack(data['linemask_all'])
        linemask = np.hstack(data['linemask'])
        if np.all(coeff == 0):
            _smooth_continuum = np.zeros_like(bestfit)
        else:
            _smooth_continuum, _ = self.smooth_continuum(specwave, apcorr*specflux - bestfit, specivar/apcorr**2,
                                                         redshift, linemask=linemask, png=png)

        # Unpack the continuum into individual cameras.
        continuummodel, smooth_continuum = [], []
        for camerapix in data['camerapix']:
            continuummodel.append(bestfit[camerapix[0]:camerapix[1]])
            smooth_continuum.append(_smooth_continuum[camerapix[0]:camerapix[1]])

        if False:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1)
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                resid = apcorr*data['flux'][icam]-continuummodel[icam]
                ax[0].plot(data['wave'][icam], resid)
                ax[1].plot(data['wave'][icam], resid-smooth_continuum[icam])
            for icam in np.arange(len(data['cameras'])): # iterate over cameras
                resid = apcorr*data['flux'][icam]-continuummodel[icam]
                pix_emlines = np.logical_not(data['linemask'][icam]) # affected by line = True
                ax[0].scatter(data['wave'][icam][pix_emlines], resid[pix_emlines], s=30, color='red')
                ax[0].plot(data['wave'][icam], smooth_continuum[icam], color='k', alpha=0.7, lw=2)
            plt.savefig('junk.png')

        # Pack it in and return.
        result['APCORR'] = apcorr
        result['CONTINUUM_COEFF'][0][0:nage] = coeff
        result['CONTINUUM_RCHI2'][0] = chi2min
        result['VDISP'][0] = vdispbest # * u.kilometer/u.second
        result['VDISP_IVAR'][0] = vdispivar # * (u.second/u.kilometer)**2
        result['AV'][0] = AV # * u.mag
        result['AGE'] = age # * u.Gyr
        result['ZZSUN'] = zzsun
        result['LOGMSTAR'] = logmstar
        result['SFR'] = sfr
        result['FAGN'] = fagn
        result['DN4000'] = dn4000
        result['DN4000_IVAR'] = dn4000_ivar
        result['DN4000_MODEL'] = dn4000_model

        for icam, cam in enumerate(data['cameras']):
            nonzero = continuummodel[icam] != 0
            #nonzero = np.abs(continuummodel[icam]) > 1e-5
            if np.sum(nonzero) > 0:
                corr = np.median(smooth_continuum[icam][nonzero] / continuummodel[icam][nonzero])
                result['CONTINUUM_SMOOTHCORR_{}'.format(cam.upper())] = corr * 100 # [%]

        self.log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
            result['CONTINUUM_SMOOTHCORR_B'][0], result['CONTINUUM_SMOOTHCORR_R'][0],
            result['CONTINUUM_SMOOTHCORR_Z'][0]))
        
        return result, continuummodel, smooth_continuum
    
    def qa_fastphot(self, data, fastphot, metadata, coadd_type='healpix',
                    outdir=None, outprefix=None):
        """QA of the best-fitting continuum.

        """
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        sns.set(context='talk', style='ticks', font_scale=1.1)#, rc=rc)

        col1 = colors.to_hex('darkblue') # 'darkgreen', 'darkred', 'dodgerblue', 'darkseagreen', 'orangered']]
        
        ymin, ymax = 1e6, -1e6

        redshift = metadata['Z']

        if metadata['PHOTSYS'] == 'S':
            filters = self.decam
            allfilters = self.decamwise
        else:
            filters = self.bassmzls
            allfilters = self.bassmzlswise

        if outdir is None:
            outdir = '.'
        if outprefix is None:
            outprefix = 'fastphot'

        if coadd_type == 'healpix':
            title = 'Survey/Program/HealPix: {}/{}/{}, TargetID: {}'.format(
                    metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
        elif coadd_type == 'cumulative':
            title = 'Tile/ThruNight: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], coadd_type, metadata['TARGETID']))
        elif coadd_type == 'pernight':
            title = 'Tile/Night: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID'],
                    metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID']))
        elif coadd_type == 'perexp':
            title = 'Tile/Night/Expid: {}/{}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['EXPID'],
                    metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'],
                    metadata['EXPID'], metadata['TARGETID']))
        else:
            pass

        # rebuild the best-fitting photometric model fit
        inodust = np.ndarray.item(np.where(self.AV == 0)[0]) # should always be index 0        
        continuum_phot, synthmodelphot = self.SSP2data(
            self.sspflux_dustnomvdisp[:, :, inodust], self.sspwave, redshift=redshift,
            synthphot=True, AV=fastphot['CONTINUUM_AV'], #test=True,
            coeff=fastphot['CONTINUUM_COEFF'] * self.massnorm)

        continuum_wave_phot = self.sspwave * (1 + redshift)

        wavemin, wavemax = 0.1, 30.0 # 6.0
        indx = np.where((continuum_wave_phot/1e4 > wavemin) * (continuum_wave_phot/1e4 < wavemax))[0]     

        phot = self.parse_photometry(self.bands,
                                     maggies=np.array([metadata['FLUX_{}'.format(band.upper())] for band in self.bands]),
                                     ivarmaggies=np.array([metadata['FLUX_IVAR_{}'.format(band.upper())] for band in self.bands]),
                                     lambda_eff=allfilters.effective_wavelengths.value,
                                     min_uncertainty=self.min_uncertainty)
        fiberphot = self.parse_photometry(self.fiber_bands,
                                          maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in self.fiber_bands]),
                                          lambda_eff=filters.effective_wavelengths.value)
        
        fig, ax = plt.subplots(figsize=(12, 8))

        if np.all(continuum_phot[indx] <= 0):
            self.log.warning('Best-fitting photometric continuum is all zeros or negative!')
            continuum_phot_abmag = continuum_phot*0 + np.median(fiberphot['abmag'])
        else:
            indx = indx[continuum_phot[indx] > 0] # trim zeros
            factor = 10**(0.4 * 48.6) * continuum_wave_phot[indx]**2 / (C_LIGHT * 1e13) / self.fluxnorm / self.massnorm # [erg/s/cm2/A --> maggies]
            continuum_phot_abmag = -2.5*np.log10(continuum_phot[indx] * factor)
            ax.plot(continuum_wave_phot[indx] / 1e4, continuum_phot_abmag, color='tan', zorder=1)

        ax.scatter(synthmodelphot['lambda_eff']/1e4, synthmodelphot['abmag'], 
                   marker='s', s=200, color='k', facecolor='none',
                   #label=r'$grz$ (spectrum, synthesized)',
                   alpha=0.8, zorder=2)
        
        # we have to set the limits *before* we call errorbar, below!
        dm = 1.0
        good = phot['abmag_ivar'] > 0
        goodlim = phot['abmag_limit'] > 0
        if np.sum(good) > 0 and np.sum(goodlim) > 0:
            ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(phot['abmag_limit'][goodlim]), np.nanmax(continuum_phot_abmag))) + dm
            ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(phot['abmag_limit'][goodlim]), np.nanmin(continuum_phot_abmag))) - dm
        elif np.sum(good) > 0 and np.sum(goodlim) == 0:
            ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(continuum_phot_abmag))) + dm
            ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(continuum_phot_abmag))) - dm
        elif np.sum(good) == 0 and np.sum(goodlim) > 0:
            ymin = np.max((np.nanmax(phot['abmag_limit'][goodlim]), np.nanmax(continuum_phot_abmag))) + dm
            ymax = np.min((np.nanmin(phot['abmag_limit'][goodlim]), np.nanmin(continuum_phot_abmag))) - dm
        else:
            good = phot['abmag'] > 0
            goodlim = phot['abmag_limit'] > 0
            if np.sum(good) > 0 and np.sum(goodlim) > 0:
                ymin = np.max((np.nanmax(phot['abmag'][good]), np.nanmax(phot['abmag_limit'][goodlim]))) + dm
                ymax = np.min((np.nanmin(phot['abmag'][good]), np.nanmin(phot['abmag_limit'][goodlim]))) - dm
            elif np.sum(good) > 0 and np.sum(goodlim) == 0:                
                ymin = np.nanmax(phot['abmag'][good]) + dm
                ymax = np.nanmin(phot['abmag'][good]) - dm
            elif np.sum(good) == 0 and np.sum(goodlim) > 0:
                ymin = np.nanmax(phot['abmag_limit'][goodlim]) + dm
                ymax = np.nanmin(phot['abmag_limit'][goodlim]) - dm
            else:
                ymin, ymax = [30, 20]
            
        if ymin > 30:
            ymin = 30
        if np.isnan(ymin) or np.isnan(ymax):
            raise('Problem here!')

        ax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        #ax.set_ylabel(r'AB mag') 
        ax.set_ylabel(r'Apparent Brightness (AB mag)') 
        ax.set_xlim(wavemin, wavemax)
        ax.set_ylim(ymin, ymax)

        ax.set_xscale('log')

        @ticker.FuncFormatter
        def major_formatter(x, pos):
            if x > 1:
                return f'{x:.0f}'
            else:
                return f'{x:.1f}'
        
        ax.xaxis.set_major_formatter(major_formatter)
        #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
        ax.set_xticks([0.1, 0.2, 0.4, 0.6, 1.0, 1.5, 3.0, 5.0, 10.0, 20.0])

        if not self.nolegend:
            ax.set_title(title, fontsize=20)

        # integrated flux / photometry
        #ax.scatter(phot['lambda_eff']/1e4, phot['abmag'],
        #           marker='s', s=130, facecolor='red', edgecolor='k',
        #           label=r'$grzW1W2$ (imaging)', alpha=1.0, zorder=3)
        abmag = np.squeeze(phot['abmag'])
        abmag_limit = np.squeeze(phot['abmag_limit'])
        abmag_fainterr = np.squeeze(phot['abmag_fainterr'])
        abmag_brighterr = np.squeeze(phot['abmag_brighterr'])
        yerr = np.squeeze([abmag_brighterr, abmag_fainterr])

        dofit = np.where(self.bands_to_fit)[0]
        if len(dofit) > 0:
            good = np.where((abmag[dofit] > 0) * (abmag_limit[dofit] == 0))[0]
            upper = np.where(abmag_limit[dofit] > 0)[0]
            if len(good) > 0:
                ax.errorbar(phot['lambda_eff'][dofit][good]/1e4, abmag[dofit][good],
                            yerr=yerr[:, dofit[good]],
                            fmt='o', markersize=12, markeredgewidth=3, markeredgecolor=col1,
                            markerfacecolor=col1, elinewidth=3, ecolor=col1, capsize=4,
                            label=r'$grz\,W_{1}W_{2}W_{3}W_{4}$', zorder=2)
            if len(upper) > 0:
                ax.errorbar(phot['lambda_eff'][dofit][upper]/1e4, abmag_limit[dofit][upper],
                            lolims=True, yerr=0.75,
                            fmt='o', markersize=12, markeredgewidth=3, markeredgecolor=col1,
                            markerfacecolor=col1, elinewidth=3, ecolor=col1, capsize=4)

        ignorefit = np.where(self.bands_to_fit == False)[0]
        if len(ignorefit) > 0:
            good = np.where((abmag[ignorefit] > 0) * (abmag_limit[ignorefit] == 0))[0]
            upper = np.where(abmag_limit[ignorefit] > 0)[0]
            if len(good) > 0:
                ax.errorbar(phot['lambda_eff'][ignorefit][good]/1e4, abmag[ignorefit][good],
                            yerr=yerr[:, ignorefit[good]],
                            fmt='o', markersize=12, markeredgewidth=3, markeredgecolor=col1,
                            markerfacecolor='none', elinewidth=3, ecolor=col1, capsize=4)
            if len(upper) > 0:
                ax.errorbar(phot['lambda_eff'][ignorefit][upper]/1e4, abmag_limit[ignorefit][upper],
                            lolims=True, yerr=0.75, fmt='o', markersize=12, markeredgewidth=3,
                            markeredgecolor=col1, markerfacecolor='none', elinewidth=3,
                            ecolor=col1, capsize=5)

        if coadd_type == 'healpix':
            targetid_str = str(metadata['TARGETID'])
        else:
            targetid_str = '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),

        leg = {
            'targetid': targetid_str,
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastphot['CONTINUUM_RCHI2']),
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            #'zfastfastphot': r'$z_{{\\rm fastfastphot}}$={:.6f}'.format(fastphot['CONTINUUM_Z']),
            #'z': '$z$={:.6f}'.format(fastphot['CONTINUUM_Z']),
            'age': '<Age>={:.3f} Gyr'.format(fastphot['CONTINUUM_AGE']),
            'absmag_r': '$M_{{^{{0.1}}r}}={:.2f}$'.format(fastphot['ABSMAG_SDSS_R']),
            'absmag_gr': '$^{{0.1}}(g-r)={:.3f}$'.format(fastphot['ABSMAG_SDSS_G']-fastphot['ABSMAG_SDSS_R'])
            }
            
        if fastphot['LOGMSTAR'] > 0:
            leg.update({'mstar': '$\\log_{{10}}\,(M_{{*}}/M_{{\odot}})={:.3f}$'.format(fastphot['LOGMSTAR'])})
        else:
            leg.update({'mstar': '$\\log_{{10}}\,(M_{{*}}/M_{{\odot}})$=...'})

        if fastphot['CONTINUUM_AV_IVAR'] == 0:
            leg.update({'AV': '$A(V)$={:.2f} mag'.format(fastphot['CONTINUUM_AV'])})
        else:
            leg.update({'AV': '$A(V)$={:.3f}+/-{:.3f} mag'.format(
                fastphot['CONTINUUM_AV'], 1/np.sqrt(fastphot['CONTINUUM_AV_IVAR']))})

        bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.25)
        legfntsz = 16

        if not self.nolegend and fastphot['CONTINUUM_RCHI2'] > 0:
            legxpos, legypos = 0.04, 0.94
            txt = '\n'.join((
                r'{}'.format(leg['mstar']),
                r'{}'.format(leg['absmag_r']),
                r'{}'.format(leg['absmag_gr'])
                ))
            ax.text(legxpos, legypos, txt, ha='left', va='top',
                    transform=ax.transAxes, fontsize=legfntsz,
                    bbox=bbox)
            
        if not self.nolegend:
            legxpos, legypos = 0.98, 0.06
            txt = '\n'.join((
                r'{}'.format(leg['zredrock']),
                r'{} {}'.format(leg['chi2'], leg['age']),
                r'{}'.format(leg['AV']),
                ))
            ax.text(legxpos, legypos, txt, ha='right', va='bottom',
                    transform=ax.transAxes, fontsize=legfntsz,
                    bbox=bbox)
    
        plt.subplots_adjust(bottom=0.14, right=0.95, top=0.93)

        self.log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()
