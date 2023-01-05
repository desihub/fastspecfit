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
from astropy.table import Table, Column, vstack

from fastspecfit.util import C_LIGHT
from fastspecfit.continuum import ContinuumTools

from desiutil.log import get_logger
log = get_logger()

## ridiculousness! - this seems to come from healpy, blarg
#import tempfile
#os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

def _fastspec_one(args):
    """Multiprocessing wrapper."""
    return fastspec_one(*args)

def _desiqa_one(args):
    """Multiprocessing wrapper."""
    return desiqa_one(*args)

def _assign_units_to_columns(fastfit, metadata, Spec, FFit, fastphot=False):
    """Assign astropy units to output tables."""
    fastcols = fastfit.colnames
    metacols = metadata.colnames

    T, M = Spec.init_output(FFit=FFit, fastphot=fastphot)
    for col in T.colnames:
        if col in fastcols:
            fastfit[col].unit = T[col].unit
    for col in M.colnames:
        if col in metacols:
            metadata[col].unit = M[col].unit

def fastspec_one(iobj, data, out, meta, FFit, broadlinefit=True, fastphot=False):
    """Multiprocessing wrapper to run :func:`fastspec` on a single object."""
    
    log.info('Working on object {} [targetid={}, z={:.6f}].'.format(
        iobj, meta['TARGETID'], meta['Z']))

    continuummodel, smooth_continuum = FFit.continuum_specfit(data, out, fastphot=fastphot)

    # Fit the emission-line spectrum.
    if fastphot:
        emmodel = None
    else:
        emmodel = FFit.emline_specfit(data, out, continuummodel, smooth_continuum,
                                      broadlinefit=broadlinefit)

    return out, meta, emmodel

def desiqa_one(FFit, data, fastfit, metadata, coadd_type,
               fastphot=False, outdir=None, outprefix=None, webqa=False):
    """Multiprocessing wrapper to generate QA for a single object."""

    #t0 = time.time()
    if webqa:
        build_webqa(FFit, data, fastfit, metadata, coadd_type=coadd_type,
                    outprefix=outprefix, outdir=outdir)
    elif fastphot:
        FFit.qa_fastphot(data, fastfit, metadata, coadd_type=coadd_type,
                         outprefix=outprefix, outdir=outdir)
    else:
        FFit.qa_fastspec(data, fastfit, metadata, coadd_type=coadd_type,
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

def fastspec(fastphot=False, args=None, comm=None):
    """Main fastspec script.

    This script is the engine to model one or more DESI spectra. It initializes
    the :class:`FastFit` class, reads the data, fits each spectrum (with the
    option of fitting in parallel), and writes out the results as a
    multi-extension binary FITS table.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line. 
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    from astropy.table import Table, vstack
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the fitting class. Note: trim the wavelengths of the SSPs to
    # optimize compute time.
    t0 = time.time()
    FFit = FastFit(ssptemplates=args.ssptemplates, mapdir=args.mapdir, 
                   verbose=args.verbose, solve_vdisp=args.solve_vdisp, 
                   minsspwave=500.0, maxsspwave=40e4)
    Spec = DESISpectra(dr9dir=args.dr9dir)
    log.info('Initializing the classes took {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    Spec.select(args.redrockfiles, firsttarget=args.firsttarget,
                targetids=targetids, ntargets=args.ntargets,
                redrockfile_prefix=args.redrockfile_prefix,
                specfile_prefix=args.specfile_prefix,
                qnfile_prefix=args.qnfile_prefix)
    if len(Spec.specfiles) == 0:
        return

    data = Spec.read_and_unpack(FFit, fastphot=fastphot, synthphot=True, mp=args.mp)
    log.info('Reading and unpacking {} spectra to be fitted took {:.2f} seconds.'.format(
        Spec.ntargets, time.time()-t0))

    t0 = time.time()
    out, meta = Spec.init_output(data, FFit=FFit, fastphot=fastphot)
    log.info('Initializing the output tables took {:.2f} seconds.'.format(time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], FFit, args.broadlinefit, fastphot)
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
    if fastphot:
        modelspectra = None
    else:
        try:
            # need to vstack to preserve the wavelength metadata 
            modelspectra = vstack(_out[2], metadata_conflicts='error')
        except:
            errmsg = 'Metadata conflict when stacking model spectra.'
            log.critical(errmsg)
            raise ValueError(errmsg)
       
    log.info('Fitting {} object(s) took {:.2f} seconds.'.format(Spec.ntargets, time.time()-t0))

    # Assign units and write out.
    _assign_units_to_columns(out, meta, Spec, FFit, fastphot=fastphot)

    write_fastspecfit(out, meta, modelspectra=modelspectra, outfile=args.outfile,
                      specprod=Spec.specprod, coadd_type=Spec.coadd_type,
                      fastphot=fastphot)

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
    fastspec(fastphot=True, args=args, comm=comm)

def build_webqa(FFit, data, fastfit, metadata, coadd_type='healpix',
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
    
    from fastspecfit.util import ivar2var

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

    _emlinemodel = EMFit.emlinemodel_bestfit(data['wave'], data['res'], fastfit)

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
    def __init__(self, ssptemplates=None, minsspwave=None, maxsspwave=40e4, 
                 minspecwave=3500.0, maxspecwave=9900.0, chi2_default=0.0, 
                 maxiter=5000, accuracy=1e-2, nolegend=False, solve_vdisp=True, 
                 constrain_age=False, mapdir=None, fastphot=False, verbose=False):
        """Class to model a galaxy stellar continuum.

        Parameters
        ----------
        ssptemplates : :class:`str`, optional
            Full path to the SSP templates used for continuum-fitting.
        minsspwave : :class:`float`, optional, defaults to None
            Minimum SSP wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxsspwave : :class:`float`, optional, defaults to 6e4
            Maximum SSP wavelength to read into memory. 
        minspecwave : :class:`float`, optional, defaults to 3000 A.
            Minimum observed-frame wavelength, which is used internally for the
            forward modeling of the emission-line spectrum.
        maxspecwave : :class:`float`, optional, defaults to 3000 A.
            Like `minspecwave` but the maximum observed-frame wavelength.
        chi2_default : :class:`float`, optional, defaults to 0.0.
            Default chi2 value for a emission line not fitted.
        maxiter : :class:`int`, optional, defaults to 5000.
            Maximum number of iterations.
        accuracy : :class:`float`, optional, defaults to 0.01.
            Fitting accuracy.
        nolegend : :class:`bool`, optional, defaults to `False`.
            Do not render the legend on the QA output.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.

        Notes
        -----
        Need to document all the attributes.
        
        Plans for improvement:
          - Update the continuum redshift using cross-correlation.
          - Don't draw reddening from a flat distribution (try gamma or a custom
            distribution of the form x**2*np.exp(-2*x/scale).

        """
        super(FastFit, self).__init__(ssptemplates=ssptemplates, minsspwave=minsspwave,
                                      maxsspwave=maxsspwave, mapdir=mapdir, fastphot=fastphot,
                                      verbose=verbose)

        self.nolegend = nolegend

        # continuum stuff
        self.constrain_age = constrain_age
        self.solve_vdisp = solve_vdisp

        # emission line stuff
        if not fastphot:
            self.chi2_default = chi2_default
            self.maxiter = maxiter
            self.accuracy = accuracy
            self.nolegend = nolegend
    
            self.emwave_pixkms = 5.0                                  # pixel size for internal wavelength array [km/s]
            self.dlogwave = self.emwave_pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
            self.log10wave = np.arange(np.log10(minspecwave), np.log10(maxspecwave), self.dlogwave)
    
            # default line-sigma for computing upper limits
            self.limitsigma_narrow = 75.0
            self.limitsigma_broad = 1200.0 
            self.wavepad = 5.0 # Angstrom
    
            # Establish the names of the parameters and doublets here, at
            # initialization, because we use them when instantiating the best-fit
            # model not just when fitting.
            doublet_names = ['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio']
            doublet_pairs = ['mgii_2803_amp', 'oii_3729_amp', 'sii_6716_amp']
    
            param_names = []
            for param in ['amp', 'vshift', 'sigma']:
                for linename in self.linetable['name'].data:
                    param_name = linename+'_'+param
                    # Use doublet-ratio parameters for close or physical
                    # doublets. Note that any changes here need to be propagated to
                    # the XX method, which needs to "know" about these doublets.
                    if param_name == 'mgii_2796_amp':
                        param_name = 'mgii_doublet_ratio' # MgII 2796/2803
                    if param_name == 'oii_3726_amp':
                        param_name = 'oii_doublet_ratio'  # [OII] 3726/3729
                    if param_name == 'sii_6731_amp':
                        param_name = 'sii_doublet_ratio'  # [SII] 6731/6716
                    param_names.append(param_name)
            self.param_names = np.hstack(param_names)
    
            self.doubletindx = np.hstack([np.where(self.param_names == doublet)[0] for doublet in doublet_names])
            self.doubletpair = np.hstack([np.where(self.param_names == pair)[0] for pair in doublet_pairs])

    def init_output(self, nobj=1, fastphot=False):
        """Initialize the output data table for this class.

        """
        nssp_coeff = len(self.sspinfo)
        
        out = Table()
        out.add_column(Column(name='APCORR', length=nobj, dtype='f4')) # aperture correction
        out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8')) # redshift
        out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
        out.add_column(Column(name='CONTINUUM_RCHI2', length=nobj, dtype='f4')) # reduced chi2
        #out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))

        if not fastphot:
            for cam in ['B', 'R', 'Z']:
                out.add_column(Column(name='CONTINUUM_SNR_{}'.format(cam), length=nobj, dtype='f4')) # median S/N in each camera
            for cam in ['B', 'R', 'Z']:
                out.add_column(Column(name='CONTINUUM_SMOOTHCORR_{}'.format(cam), length=nobj, dtype='f4')) 

        out.add_column(Column(name='VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        out.add_column(Column(name='VDISP_IVAR', length=nobj, dtype='f4', unit=u.second**2/u.kilometer**2))
        out.add_column(Column(name='AV', length=nobj, dtype='f4', unit=u.mag))
        out.add_column(Column(name='AGE', length=nobj, dtype='f4', unit=u.Gyr))
        out.add_column(Column(name='ZZSUN', length=nobj, dtype='f4'))
        out.add_column(Column(name='LOGMSTAR', length=nobj, dtype='f4', unit=u.solMass))
        out.add_column(Column(name='SFR', length=nobj, dtype='f4', unit=u.solMass/u.year))
        out.add_column(Column(name='FAGN', length=nobj, dtype='f4'))
        
        out.add_column(Column(name='DN4000', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_IVAR', length=nobj, dtype='f4'))
        out.add_column(Column(name='DN4000_MODEL', length=nobj, dtype='f4'))

        # observed-frame photometry synthesized from the spectra
        for band in self.synth_bands:
            out.add_column(Column(name='FLUX_SYNTH_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies')) 
            #out.add_column(Column(name='FLUX_SYNTH_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies-2'))
        # observed-frame photometry synthesized from the best-fitting continuum model fit
        for band in self.synth_bands:
            out.add_column(Column(name='FLUX_SYNTH_MODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies'))

        for band in self.absmag_bands:
            out.add_column(Column(name='KCORR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.mag))
            out.add_column(Column(name='ABSMAG_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.mag)) # absolute magnitudes
            out.add_column(Column(name='ABSMAG_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/u.mag**2))

        for cflux in ['LOGLNU_1500', 'LOGLNU_2800']:
            out.add_column(Column(name=cflux, length=nobj, dtype='f4', unit=10**(-28)*u.erg/u.second/u.Hz))
        out.add_column(Column(name='LOGL_5100', length=nobj, dtype='f4', unit=10**(10)*u.solLum))

        #for cflux in ['FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT']:
        #    out.add_column(Column(name=cflux, length=nobj, dtype='f4', unit=10**(-17)*u.erg/(u.second*u.cm**2*u.Angstrom)))

        if not fastphot:
            # Add chi2 metrics
            out.add_column(Column(name='RCHI2', length=nobj, dtype='f4')) # full-spectrum reduced chi2
            #out.add_column(Column(name='DOF', length=nobj, dtype='i8')) # full-spectrum dof
            out.add_column(Column(name='LINERCHI2_BROAD', length=nobj, dtype='f4')) # reduced chi2 with broad line-emission
            #out.add_column(Column(name='DOF_BROAD', length=nobj, dtype='i8'))
            out.add_column(Column(name='DELTA_LINERCHI2', length=nobj, dtype='f4')) # delta-reduced chi2 with and without broad line-emission
    
            out.add_column(Column(name='NARROW_Z', length=nobj, dtype='f8'))
            #out.add_column(Column(name='NARROW_Z_ERR', length=nobj, dtype='f8'))
            out.add_column(Column(name='BROAD_Z', length=nobj, dtype='f8'))
            #out.add_column(Column(name='BROAD_Z_ERR', length=nobj, dtype='f8'))
            out.add_column(Column(name='UV_Z', length=nobj, dtype='f8'))
            #out.add_column(Column(name='UV_Z_ERR', length=nobj, dtype='f8'))
    
            out.add_column(Column(name='NARROW_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
            #out.add_column(Column(name='NARROW_SIGMA_ERR', length=nobj, dtype='f4', unit=u.kilometer / u.second))
            out.add_column(Column(name='BROAD_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
            #out.add_column(Column(name='BROAD_SIGMA_ERR', length=nobj, dtype='f4', unit=u.kilometer / u.second))
            out.add_column(Column(name='UV_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
            #out.add_column(Column(name='UV_SIGMA_ERR', length=nobj, dtype='f4', unit=u.kilometer / u.second))
    
            # special columns for the fitted doublets
            out.add_column(Column(name='MGII_DOUBLET_RATIO', length=nobj, dtype='f4'))
            out.add_column(Column(name='OII_DOUBLET_RATIO', length=nobj, dtype='f4'))
            out.add_column(Column(name='SII_DOUBLET_RATIO', length=nobj, dtype='f4'))
    
            for line in self.linetable['name']:
                line = line.upper()
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
                out.add_column(Column(name='{}_CHI2'.format(line), data=np.repeat(self.chi2_default, nobj), dtype='f4'))
                out.add_column(Column(name='{}_NPIX'.format(line), length=nobj, dtype=np.int32))

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

        nage = len(coeff)
        
        # From Taylor+11, eq 8
        #mstar = self.sspinfo['mstar'][:nage].dot(coeff) * self.massnorm
        #https://researchportal.port.ac.uk/ws/files/328938/MNRAS_2011_Taylor_1587_620.pdf
        #mstar = 1.15 + 0.7*(absmag[1]-absmag[3]) - 0.4*absmag[3]

        # compute the model continuum flux at 1500 and 2800 A (to facilitate UV
        # luminosity-based SFRs) and at the positions of strong nebular emission
        # lines [OII], Hbeta, [OIII], and Halpha
        dfactor = (1 + redshift) * 4.0 * np.pi * self.cosmo.luminosity_distance(redshift).to(u.cm).value**2 / self.fluxnorm

        lums = {}
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

        #cfluxes = {}
        #cwaves = [3728.483, 4862.683, 5008.239, 6564.613]
        #labels = ['FOII_3727_CONT', 'FHBETA_CONT', 'FOIII_5007_CONT', 'FHALPHA_CONT']
        #for cwave, label in zip(cwaves, labels):
        #    J = (self.sspwave > cwave-500) * (self.sspwave < cwave+500)
        #    I = (self.sspwave[J] > cwave-20) * (self.sspwave[J] < cwave+20)
        #    smooth = median_filter(continuum[J], 200)
        #    clipflux, _, _ = sigmaclip(smooth[I], low=1.5, high=3)
        #    cflux = np.median(clipflux) # [flux in 10**-17 erg/s/cm2/A]
        #    cfluxes[label] = cflux # * u.erg/(u.second*u.cm**2*u.Angstrom)
        #    
        #    #import matplotlib.pyplot as plt
        #    #print(cwave, cflux)
        #    #plt.clf()
        #    #plt.plot(self.sspwave[J], continuum[J])
        #    #plt.plot(self.sspwave[J], smooth, color='k')
        #    #plt.axhline(y=cflux, color='red')
        #    #plt.axvline(x=cwave, color='red')
        #    #plt.xlim(cwave - 50, cwave + 50)
        #    #plt.savefig('desi-users/ioannis/tmp/junk.png')
        #    #pdb.set_trace()

        return kcorr, absmag, ivarabsmag, bestmaggies, lums

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

    def _qa_dn4000(self):
        """Simple QA of the Dn(4000) estimate."""

        from fastspecfit.util import ivar2var            
        import matplotlib.pyplot as plt

        print(dn4000, dn4000_model, 1/np.sqrt(dn4000_ivar))

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

    def continuum_specfit(self, data, result, fastphot=False):
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
            in :func:`self.init_output`.

        Notes
        -----
          - Consider using cross-correlation to update the redrock redshift.
          - We solve for velocity dispersion if solve_vdisp=True or ((SNR_B>3 or
            SNR_R>3) and REDSHIFT<1).

        """
        tall = time.time()

        redshift = result['CONTINUUM_Z']

        objflam = data['phot']['flam'].data * self.fluxnorm
        objflamivar = (data['phot']['flam_ivar'].data / self.fluxnorm**2) * self.bands_to_fit
        assert(np.all(objflamivar >= 0))

        # Optionally ignore templates which are older than the age of the
        # universe at the redshift of the object.
        if self.constrain_age:
            agekeep = self.younger_than_universe(redshift)
        else:
            agekeep = np.arange(self.nsed)
        nage = len(agekeep)

        # Photometry-only fitting.
        if fastphot:
            vdispbest, vdispivar = self.vdisp_nominal, 0.0
            self.log.info('Adopting nominal vdisp={:.2f} km/s.'.format(self.vdisp_nominal))

            # Get the coefficients and chi2 at the nominal velocity dispersion. 
            t0 = time.time()
            bestsspflux, bestphot = self.SSP2data(self.sspflux_vdisp[:, agekeep, self.vdisp_nominal_indx], 
                                                  self.sspwave, redshift=redshift,
                                                  south=data['photsys'] == 'S', synthphot=True)
            bestflam = bestphot['flam'].data*self.massnorm*self.fluxnorm
    
            coeff, chi2min = self._call_nnls(bestflam, objflam, objflamivar)
            chi2min /= np.sum(objflamivar > 0) # dof???
            self.log.info('Fitting {} models took {:.2f} seconds.'.format(
                nage, time.time()-t0))

            smooth_continuum = None
            continuummodel = bestsspflux.dot(coeff)

            dn4000_model, _ = self.get_dn4000(self.sspwave, continuummodel, rest=True)

            self.log.info('Model Dn(4000)={:.3f}.'.format(dn4000_model))

        else:
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
    
            # Aperture correction based on the r-band synthesized and observed
            # photometry. Note that the right thing to do is to do a quick-and-dirty
            # fit of the spectrum and used the model-synthesized flux (ignoring any
            # smooth continuum correction), which will guarantee that the aperture
            # correction is always positive.
            apcorr = data['phot']['nanomaggies'][1] / data['synthphot']['nanomaggies'][1]
            data['apcorr'] = apcorr
    
            # Prepare the spectral and photometric models at the galaxy
            # redshift. And if the wavelength coverage is sufficient, also solve for
            # the velocity dispersion.
    
            #compute_vdisp = ((result['CONTINUUM_SNR_B'] > 1) and (result['CONTINUUM_SNR_R'] > 1)) and (redshift < 1.0)
            restwave = specwave / (1+redshift)
            Ivdisp = np.where((specivar > 0) * (restwave > 3500.0) * (restwave < 5500.0))[0]
            compute_vdisp = (len(Ivdisp) > 0) * (np.ptp(restwave[Ivdisp]) > 500.0)
    
            self.log.info('S/N_B={:.2f}, S/N_R={:.2f}, rest wavelength coverage={:.0f}-{:.0f} A.'.format(
                result['CONTINUUM_SNR_B'], result['CONTINUUM_SNR_R'], restwave[0], restwave[-1]))
    
            if self.solve_vdisp or compute_vdisp:
                
                # Only use models which likely affect the velocity dispersion.
                Mvdisp = np.where((self.sspinfo[agekeep]['fagn'] == 0) * (self.sspinfo[agekeep]['av'] == 0.0))[0]
                nMvdisp = len(Mvdisp)
    
                t0 = time.time()
                zsspflux_vdisp, _ = self.SSP2data(
                    self.sspflux_vdisp[:, agekeep[Mvdisp], :], self.sspwave, # [npix,nsed,nvdisp]
                    redshift=redshift, specwave=data['wave'], specres=data['res'],
                    cameras=data['cameras'], synthphot=False)#, south=data['photsys'] == 'S')
    
                zsspflux_vdisp = np.concatenate(zsspflux_vdisp, axis=0)  # [npix,nMvdisp*nvdisp]
                npix, nmodel = zsspflux_vdisp.shape
                zsspflux_vdisp = zsspflux_vdisp.reshape(npix, nMvdisp, self.nvdisp) # [npix,nMvdisp,nvdisp]
    
                vdispchi2min, vdispbest, vdispivar, _ = self._call_nnls(
                    zsspflux_vdisp[Ivdisp, :, :], specflux[Ivdisp], specivar[Ivdisp],
                    xparam=self.vdisp, xlabel=r'$\sigma$ (km/s)', debug=False)#True)
                self.log.info('Fitting for the velocity dispersion with {}/{} models took {:.2f} seconds.'.format(
                    nMvdisp, nage, time.time()-t0))
    
                if vdispivar > 0:
                    # protect against tiny ivar from becomming infinite in the output table
                    if vdispivar < 1e-5:
                        self.log.warning('vdisp inverse variance is tiny; capping at 1e-5')
                        vdispivar = 1e-5
                    self.log.info('Best-fitting vdisp={:.2f}+/-{:.2f} km/s.'.format(
                        vdispbest, 1/np.sqrt(vdispivar)))
                else:
                    vdispbest = self.vdisp_nominal
                    self.log.info('Finding vdisp failed; adopting vdisp={:.2f} km/s'.format(self.vdisp_nominal))

                # Get the final set of coefficients and chi2 at the best-fitting velocity dispersion. 
                t0 = time.time()
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
                self.log.info('Final fitting with {} models took {:.2f} seconds.'.format(
                    nage, time.time()-t0))
            else:
                if compute_vdisp:
                    self.log.info('Sufficient wavelength covereage to compute vdisp but solve_vdisp=False; adopting nominal vdisp={:.2f} km/s.'.format(
                        self.vdisp_nominal))
                else:
                    self.log.info('Insufficient wavelength covereage to compute vdisp; adopting nominal vdisp={:.2f} km/s'.format(
                        self.vdisp_nominal))
                    
                vdispbest, vdispivar = self.vdisp_nominal, 0.0

                # Get the final set of coefficients and chi2 at the nominal velocity dispersion. 
                t0 = time.time()
                bestsspflux, bestphot = self.SSP2data(self.sspflux_vdisp[:, agekeep, self.vdisp_nominal_indx], 
                                                      self.sspwave, redshift=redshift,
                                                      specwave=data['wave'], specres=data['res'],
                                                      vdisp=None, cameras=data['cameras'],
                                                      south=data['photsys'] == 'S', synthphot=True)
                bestsspflux = np.concatenate(bestsspflux, axis=0)
                bestflam = bestphot['flam'].data*self.massnorm*self.fluxnorm
        
                coeff, chi2min = self._call_nnls(np.vstack((bestflam, bestsspflux)),
                                                 np.hstack((objflam, specflux*apcorr)),
                                                 np.hstack((objflamivar, specivar/apcorr**2)))
                chi2min /= (np.sum(objflamivar > 0) + np.sum(specivar > 0)) # dof???
                self.log.info('Final fitting with {} models took {:.2f} seconds.'.format(
                    nage, time.time()-t0))

            bestfit = bestsspflux.dot(coeff)

            # Get DN(4000).
            flam_ivar = np.hstack(data['ivar']) # specivar is line-masked so we can't use it!
            dn4000, dn4000_ivar = self.get_dn4000(specwave, specflux, flam_ivar=flam_ivar, 
                                                  redshift=redshift, rest=False)

            print('NEED TO GET D4000_MODEL FROM THE BEST-FITTING FULL SED!')

            dn4000_model, _ = self.get_dn4000(specwave, bestfit, redshift=redshift, rest=False)
            #self._qa_dn4000()

            if dn4000_ivar > 0:
                self.log.info('Spectroscopic DN(4000)={:.3f}+/-{:.3f}, Model Dn(4000)={:.3f}'.format(
                    dn4000, 1/np.sqrt(dn4000_ivar), dn4000_model))
            else:
                self.log.info('Spectroscopic DN(4000)={:.3f}, Model Dn(4000)={:.3f}'.format(
                    dn4000, dn4000_model))

            png = None
            #png = '/global/homes/i/ioannis/desi-users/ioannis/tmp/smooth-continuum.png'
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
    
            for icam, cam in enumerate(data['cameras']):
                nonzero = continuummodel[icam] != 0
                #nonzero = np.abs(continuummodel[icam]) > 1e-5
                if np.sum(nonzero) > 0:
                    corr = np.median(smooth_continuum[icam][nonzero] / continuummodel[icam][nonzero])
                    result['CONTINUUM_SMOOTHCORR_{}'.format(cam.upper())] = corr * 100 # [%]
    
            self.log.info('Smooth continuum correction: b={:.3f}%, r={:.3f}%, z={:.3f}%'.format(
                result['CONTINUUM_SMOOTHCORR_B'], result['CONTINUUM_SMOOTHCORR_R'],
                result['CONTINUUM_SMOOTHCORR_Z']))

            #if False:
            #    import matplotlib.pyplot as plt
            #    fig, ax = plt.subplots(2, 1)
            #    for icam in np.arange(len(data['cameras'])): # iterate over cameras
            #        resid = apcorr*data['flux'][icam]-continuummodel[icam]
            #        ax[0].plot(data['wave'][icam], resid)
            #        ax[1].plot(data['wave'][icam], resid-smooth_continuum[icam])
            #    for icam in np.arange(len(data['cameras'])): # iterate over cameras
            #        resid = apcorr*data['flux'][icam]-continuummodel[icam]
            #        pix_emlines = np.logical_not(data['linemask'][icam]) # affected by line = True
            #        ax[0].scatter(data['wave'][icam][pix_emlines], resid[pix_emlines], s=30, color='red')
            #        ax[0].plot(data['wave'][icam], smooth_continuum[icam], color='k', alpha=0.7, lw=2)
            #    plt.savefig('junk.png')

        # # Compute K-corrections, rest-frame quantities, and physical properties.
        kcorr, absmag, ivarabsmag, synth_bestmaggies, lums = self.kcorr_and_absmag(data, continuummodel, coeff)

        AV = self.get_mean_property('av', coeff, agekeep)                        # [mag]
        age = self.get_mean_property('age', coeff, agekeep, normalization=1e9)   # [Gyr]
        zzsun = self.get_mean_property('zzsun', coeff, agekeep, log10=False)     # [log Zsun]
        fagn = self.get_mean_property('fagn', coeff, agekeep)
        logmstar = self.get_mean_property('mstar', coeff, agekeep, normalization=1/self.massnorm, log10=True) # [Msun]
        sfr = self.get_mean_property('sfr', coeff, agekeep, normalization=1/self.massnorm, log10=False)       # [Msun/yr]

        self.log.info('A(V)={:.3f}, Age={:.3f} Gyr, Mstar={:.4g} Msun, SFR={:.3f} Msun/yr, Z/Zsun={:.3f}, fagn={:.3f}'.format(
            AV, age, logmstar, sfr, zzsun, fagn))

        # Pack it in and return.
        result['CONTINUUM_COEFF'][agekeep] = coeff
        result['CONTINUUM_RCHI2'] = chi2min
        result['VDISP'] = vdispbest # * u.kilometer/u.second
        result['VDISP_IVAR'] = vdispivar # * (u.second/u.kilometer)**2
        result['AV'] = AV # * u.mag
        result['AGE'] = age # * u.Gyr
        result['ZZSUN'] = zzsun
        result['LOGMSTAR'] = logmstar
        result['SFR'] = sfr
        result['FAGN'] = fagn
        result['DN4000_MODEL'] = dn4000_model

        if not fastphot:
            result['APCORR'] = apcorr
            result['DN4000'] = dn4000
            result['DN4000_IVAR'] = dn4000_ivar

        log.info('Continuum-fitting took {:.2f} seconds.'.format(time.time()-tall))

        return continuummodel, smooth_continuum

    def build_linemodels(self, redshift, wavelims=[3000, 10000], verbose=False):
        """Build all the multi-parameter emission-line models we will use.
    
        """
        def _print_linemodel(linemodel):
            for linename in linenames:
                for param in ['amp', 'sigma', 'vshift']:
                    I = np.where(param_names == linename+'_'+param)[0]
                    if len(I) == 1:
                        I = I[0]
                        if linemodel['tiedtoparam'][I] == -1:
                            if linemodel['fixed'][I]:
                                print('{:25s} is FIXED'.format(linename+'_'+param))
                        else:
                            if linemodel['fixed'][I]:
                                print('{:25s} tied to {:25s} with factor {:.4f} and FIXED'.format(
                                    linename+'_'+param, param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))
                            else:
                                print('{:25s} tied to {:25s} with factor {:.4f}'.format(
                                    linename+'_'+param, param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))

        def _fix_parameters(linemodel, verbose=False):
            """Set the "fixed" attribute for all the parameters in a given linemodel."""
            # First loop through all tied parameters and set fixed to the
            # parameter it's tied to.
            I = np.where(linemodel['tiedtoparam'] != -1)[0] # should always have len(I)>0
            alltied = linemodel[I]['tiedtoparam']
            utied = np.unique(alltied)
            for tied in utied:
                J = tied == alltied
                if verbose:
                    print('Tying {} to {}'.format(' '.join(linemodel[I][J]['param_name']), linemodel[tied]['param_name']))
                linemodel[I][J]['fixed'] = linemodel[tied]['fixed']
            if verbose:
                print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

            # Next, fix out-of-range lines but not those that are in the 'utied'
            # array---those out-of-range lines need to be in the optimization
            # list because the in-range lines depend on them.
            outofrange = fit_linetable['inrange'] == False
            if np.sum(outofrange) > 0: # should always be true
                for linename in fit_linetable['name'][outofrange]:
                    for param in ['amp', 'vshift', 'sigma']:
                        param_name = linename+'_'+param
                        I = np.where(linemodel['param_name'] == param_name)[0]
                        if I in utied:
                            if verbose:
                                print('Not fixing out-of-range parameter {}'.format(param_name))
                        else:
                            linemodel['fixed'][I] |= True
                if verbose:
                    print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                    print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                    #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

                # Finally loop through each 'utied' line and if all the lines
                # tied to it are fixed, then fix that line, too.
                for tied in utied:
                    if linemodel['param_name'][tied] and np.all(linemodel[linemodel['tiedtoparam'] == tied]['fixed']):
                        if outofrange[linemodel['linename'][tied] == fit_linetable['name']]:
                            if verbose:
                                print('Fixing {} because line is out of range and all tied lines are fixed: {}'.format(
                                    linemodel['param_name'][tied], ' '.join(linemodel[linemodel['tiedtoparam'] == tied]['param_name'])))
                            linemodel[tied]['fixed'] = True

                # Also handle the doublets.
                I = np.where(linemodel['doubletpair'] != -1)[0]
                if len(I) > 0:
                    for doublet in linemodel[I]['doubletpair']:
                        J = doublet == linemodel['doubletpair']
                        if linemodel[doublet]['fixed']:
                            linemodel['fixed'][J] = True

                if verbose:
                    print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                    print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                    #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

        initvshift = 1.0
        vmaxshift_narrow = 500.0
        vmaxshift_broad = 2500.0 # 3000.0
    
        minsigma_narrow = 1.0
        maxsigma_narrow = 750.0 # 500.0

        minsigma_broad = 1.0
        maxsigma_broad = 1e4

        minsigma_balmer_broad = minsigma_narrow
        maxsigma_balmer_broad = maxsigma_broad
    
        # Be very careful about changing the default broad line-sigma. Smaller
        # values like 1500 km/s (which is arguably more sensible) can lead to
        # low-amplitude broad lines in a bunch of normal star-forming galaxy
        # spectra. (They act to "suck up" local continuum variations.) Also
        # recall that if it's well-measured, we use the initial line-sigma in
        # estimate_linesigma, which is a better initial guess.
        initsigma_narrow = 75.0 # 260.0 # 75.0
        initsigma_broad = 3000.0  
    
        initamp = 0.0
        #minamp = 0.0
        minamp = -1e2
        maxamp = +1e5
        minamp_balmer_broad = minamp # 0.0
        maxamp_balmer_broad = maxamp
    
        # Specialized parameters on the MgII, [OII], and [SII] doublet ratios. See
        # https://github.com/desihub/fastspecfit/issues/39. Be sure to set
        # self.doublet_names, below, and also note that any change in the order of
        # these lines has to be handled in _emline_spectrum!
        init_mgii_doublet = 0.5 # MgII 2796/2803
        init_oii_doublet = 0.74 # [OII] 3726/3729
        init_sii_doublet = 0.74 # [SII] 6731/6716

        bounds_mgii_doublet = [0.01, 10.0] 
        bounds_oii_doublet = [0.1, 2.0] # [0.5, 1.5] # [0.66, 1.4]
        bounds_sii_doublet = [0.1, 2.0] # [0.5, 1.5] # [0.67, 1.2]
    
        # Create a new line-fitting table which contains the redshift-dependent
        # quantities for this object.
        fit_linetable = Table()
        fit_linetable['name'] = self.linetable['name']
        fit_linetable['isbalmer'] = self.linetable['isbalmer']
        fit_linetable['isbroad'] = self.linetable['isbroad']
        fit_linetable['restwave'] = self.linetable['restwave']
        fit_linetable['zwave'] = self.linetable['restwave'].data * (1 + redshift)
        fit_linetable['inrange'] = ((fit_linetable['zwave'] > (wavelims[0]+self.wavepad)) * 
                                    (fit_linetable['zwave'] < (wavelims[1]-self.wavepad)))
        self.fit_linetable = fit_linetable
        
        linenames = fit_linetable['name'].data
        param_names = self.param_names
        nparam = len(param_names)

        # Model 1 -- here, parameters are minimally tied together for the final
        # fit and only lines outside the wavelength range are fixed. Includes
        # broad lines.
        final_linemodel = Table()
        final_linemodel['param_name'] = param_names
        final_linemodel['index'] = np.arange(nparam).astype(np.int32)
        final_linemodel['linename'] = np.tile(linenames, 3) # 3 parameters per line
        final_linemodel['tiedfactor'] = np.zeros(nparam, 'f8')
        final_linemodel['tiedtoparam'] = np.zeros(nparam, np.int16)-1
        final_linemodel['doubletpair'] = np.zeros(nparam, np.int16)-1
        final_linemodel['fixed'] = np.zeros(nparam, bool)
        final_linemodel['bounds'] = np.zeros((nparam, 2), 'f8')
        final_linemodel['initial'] = np.zeros(nparam, 'f8')
        final_linemodel['value'] = np.zeros(nparam, 'f8')

        final_linemodel['doubletpair'][self.doubletindx] = self.doubletpair

        # Build the relationship of "tied" parameters. In the 'tied' array, the
        # non-zero value is the multiplicative factor by which the parameter
        # represented in the 'tiedtoparam' index should be multiplied.
    
        # Physical doublets and lines in the same ionization species should have
        # their velocity shifts and line-widths always tied. In addition, set fixed
        # doublet-ratios here. Note that these constraints must be set on *all*
        # lines, not just those in range.
    
        for iline, linename in enumerate(linenames):
            # initial values and bounds - broad He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp_balmer_broad, maxamp_balmer_broad], 
                                                   [minsigma_balmer_broad, maxsigma_balmer_broad],
                                                   [-vmaxshift_broad, +vmaxshift_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - narrow He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - broad UV/QSO lines (non-Balmer)
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_broad, maxsigma_broad],
                                                   [-vmaxshift_broad, +vmaxshift_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - forbidden lines
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # tie parameters

            # broad He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]
            #print('Releasing the narrow Balmer lines!')
            # narrow He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False and linename != 'halpha':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
            # other lines
            if linename == 'mgii_2796':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'mgii_2803_'+param)[0]
            if linename == 'nev_3346' or linename == 'nev_3426': # should [NeIII] 3869 be tied to [NeV]???
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'neiii_3869_'+param)[0]
            if linename == 'oii_3726':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_3729_'+param)[0]
            # Tentative! Tie auroral lines to [OIII] 4363 but maybe we shouldn't tie [OI] 6300 here...
            if linename == 'nii_5755' or linename == 'oi_6300' or linename == 'siii_6312':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_4363_'+param)[0]
            if linename == 'oiii_4959':
                """
                [O3] (4-->2): airwave: 4958.9097 vacwave: 4960.2937 emissivity: 1.172e-21
                [O3] (4-->3): airwave: 5006.8417 vacwave: 5008.2383 emissivity: 3.497e-21
                """
                final_linemodel['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9839 # 2.8875
                final_linemodel['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oiii_5007_amp')[0]
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
            if linename == 'nii_6548':
                """
                [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.02198e-21
                [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.94901e-21
                """
                final_linemodel['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9421 # 2.936
                final_linemodel['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'nii_6584_amp')[0]
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nii_6584_'+param)[0]
            if linename == 'sii_6731':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'sii_6716_'+param)[0]
            if linename == 'oii_7320':
                """
                [O2] (5-->2): airwave: 7318.9185 vacwave: 7320.9350 emissivity: 8.18137e-24
                [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 2.40519e-23
                [O2] (5-->3): airwave: 7329.6613 vacwave: 7331.6807 emissivity: 1.35614e-23
                [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.27488e-23
                """
                final_linemodel['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 1.2251
                final_linemodel['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oii_7330_amp')[0]
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_7330_'+param)[0]
            if linename == 'siii_9069':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'siii_9532_'+param)[0]
            # Tentative! Tie SiIII] 1892 to CIII] 1908 because they're so close in wavelength.
            if linename == 'siliii_1892':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'ciii_1908_'+param)[0]
            
        # Finally set the initial values and bounds on the doublet ratio parameters.
        for param, bounds, default in zip(['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio'],
                                          [bounds_mgii_doublet, bounds_oii_doublet, bounds_sii_doublet],
                                          [init_mgii_doublet, init_oii_doublet, init_sii_doublet]):
            final_linemodel['initial'][final_linemodel['param_name'] == param] = default
            final_linemodel['bounds'][final_linemodel['param_name'] == param] = bounds
                    
        # Assign fixed=True to parameters which are outside the wavelength range
        # except those that are tied to other lines.
        _fix_parameters(final_linemodel, verbose=False)

        assert(np.all(final_linemodel['tiedtoparam'][final_linemodel['tiedfactor'] != 0] != -1))
        assert(len(final_linemodel[np.sum(final_linemodel['bounds'] == [0.0, 0.0], axis=1) > 0]) == 0)
    
        #_print_linemodel(final_linemodel)
        #final_linemodel[np.logical_and(final_linemodel['fixed'] == False, final_linemodel['tiedtoparam'] == -1)]

        # Model 2 - like final_linemodel, but broad lines have been fixed at
        # zero.
        final_linemodel_nobroad = final_linemodel.copy()
        final_linemodel_nobroad['fixed'] = False # reset

        for iline, linename in enumerate(linenames):
            if linename == 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    final_linemodel_nobroad['fixed'][param_names == linename+'_'+param] = True

            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    final_linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

        #final_linemodel_nobroad[np.logical_and(final_linemodel_nobroad['fixed'] == False, final_linemodel_nobroad['tiedtoparam'] == -1)]

        _fix_parameters(final_linemodel_nobroad, verbose=False)

        assert(np.all(final_linemodel_nobroad['tiedtoparam'][final_linemodel_nobroad['tiedfactor'] != 0] != -1))

        # Model 3 - like final_linemodel, but with all the narrow and forbidden
        # lines tied together and all the broad lines tied together.
        initial_linemodel = final_linemodel.copy()
        initial_linemodel['fixed'] = False # reset

        for iline, linename in enumerate(linenames):
            # Tie all forbidden lines and narrow Balmer & He lines to [OIII] 5007.
            if fit_linetable['isbroad'][iline] == False and linename != 'oiii_5007':
                for param in ['sigma', 'vshift']:
                    initial_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    initial_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]

            ## Tie all narrow Balmer+He lines to narrow Halpha.
            #if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False and linename != 'halpha':
            #    for param in ['sigma', 'vshift']:
            #        initial_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
            #        initial_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
    
            # Tie all broad Balmer+He lines to broad Halpha.
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['sigma', 'vshift']:
                    initial_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    initial_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

        _fix_parameters(initial_linemodel, verbose=False)#True)

        assert(np.all(initial_linemodel['tiedtoparam'][initial_linemodel['tiedfactor'] != 0] != -1))

        # Model 4 - like initial_linemodel, but broad lines have been fixed at
        # zero.
        initial_linemodel_nobroad = initial_linemodel.copy()
        initial_linemodel_nobroad['fixed'] = False # reset        

        for iline, linename in enumerate(linenames):
            if linename == 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    initial_linemodel_nobroad['fixed'][param_names == linename+'_'+param] = True

            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    initial_linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    initial_linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

        _fix_parameters(initial_linemodel_nobroad, verbose=False)

        assert(np.all(initial_linemodel_nobroad['tiedtoparam'][initial_linemodel_nobroad['tiedfactor'] != 0] != -1))

        if verbose:
            _print_linemodel(initial_linemodel_nobroad)

        return final_linemodel, final_linemodel_nobroad, initial_linemodel, initial_linemodel_nobroad

    def _initial_guesses_and_bounds(self, data, emlinewave, emlineflux):
        """For all lines in the wavelength range of the data, get a good initial guess
        on the amplitudes and line-widths. This step is critical for cases like,
        e.g., 39633354915582193 (tile 80613, petal 05), which has strong narrow
        lines.

        """
        initial_guesses, param_bounds = {}, {}
        #init_amplitudes, init_sigmas = {}, {}
    
        coadd_emlineflux = np.interp(data['coadd_wave'], emlinewave, emlineflux)
    
        for linename, linepix in zip(data['coadd_linename'], data['coadd_linepix']):
            ## skip the physical doublets
            #if not hasattr(self.EMLineModel, '{}_amp'.format(linename)):
            #    continue

            npix = len(linepix)
            if npix > 5:
                mnpx, mxpx = linepix[npix//2]-3, linepix[npix//2]+3
                if mnpx < 0:
                    mnpx = 0
                if mxpx > linepix[-1]:
                    mxpx = linepix[-1]
                amp = np.max(coadd_emlineflux[mnpx:mxpx])
            else:
                amp = np.percentile(coadd_emlineflux[linepix], 97.5)
                
            if amp < 0:
                amp = np.abs(amp)
    
            # update the bounds on the line-amplitude
            #bounds = [-np.min(np.abs(coadd_emlineflux[linepix])), 3*np.max(coadd_emlineflux[linepix])]
            mx = 5*np.max(coadd_emlineflux[linepix])
            if mx < 0: # ???
                mx = 5*np.max(np.abs(coadd_emlineflux[linepix]))
            
            # force broad Balmer lines to be positive
            iline = self.linetable[self.linetable['name'] == linename]
            if iline['isbroad']:
                if iline['isbalmer']:
                    bounds = [0.0, mx]
                else:
                    # MgII and other UV lines are dropped relatively frequently
                    # due to the lower bound on the amplitude.
                    #bounds = [None, mx]
                    #bounds = [-1e2, mx]
                    #bounds = [0.0, mx]
                    bounds = [-1.5*np.min(np.abs(coadd_emlineflux[linepix])), mx]
            else:
                #bounds = [0.0, mx]
                bounds = [-1.5*np.min(np.abs(coadd_emlineflux[linepix])), mx]

            if (bounds[0] > bounds[1]) or (amp < bounds[0]) or (amp > bounds[1]):
                errmsg = 'Initial amplitude is outside its bound for line {}!'.format(linename)
                self.log.critical(errmsg)
                raise ValueError(errmsg)

            initial_guesses[linename+'_amp'] = amp
            param_bounds[linename+'_amp'] = bounds
    
        # Now update the linewidth but here we need to loop over *all* lines
        # (not just those in range). E.g., if H-alpha is out of range we need to
        # set its initial value correctly since other lines are tied to it
        # (e.g., main-bright-32406-39628257196245904).
        for iline in self.linetable:
            linename = iline['name']
            if iline['isbroad']:
                if iline['isbalmer']: # broad Balmer lines
                    if data['linesigma_balmer'] > data['linesigma_narrow']:
                        initial_guesses[linename+'_sigma'] = data['linesigma_balmer']
                    else:
                        initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
                else: # broad UV/QSO lines
                    initial_guesses[linename+'_sigma'] = data['linesigma_uv']
            else:
                # prefer narrow over Balmer
                initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
    
        return initial_guesses, param_bounds

    def _linemodel_to_parameters(self, linemodel):
        """Convert a linemodel model to a list of emission-line parameters."""

        linesplit = np.array_split(linemodel['index'], 3) # 3 parameters per line
        #linesplit = (np.arange(3) + 1) * len(linemodel) // 3 # 3 parameters per line
        lineamps = linemodel['value'][linesplit[0]].data
        linevshifts = linemodel['value'][linesplit[1]].data
        linesigmas = linemodel['value'][linesplit[2]].data

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        #doublet = np.where(linemodel['doubletpair'] != -1)[0]
        #lineamps[doublet] *= linemodel['value'][linemodel['doubletpair'][doublet]]
        parameters = np.hstack((lineamps, linevshifts, linesigmas))

        linewaves = self.fit_linetable['restwave'].data
        #lineinrange = self.fit_linetable['inrange'].data

        #Itied = np.where((linemodel['tiedtoparam'] != -1))[0]
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        Ifree = np.where((linemodel['tiedtoparam'] == -1) * (linemodel['fixed'] == False))[0]

        tiedtoparam = linemodel['tiedtoparam'][Itied].data
        tiedfactor = linemodel['tiedfactor'][Itied].data
        bounds = linemodel['bounds'][Ifree].data

        doubletindx = np.where(linemodel['doubletpair'] != -1)[0]
        doubletpair = linemodel['doubletpair'][doubletindx].data

        parameter_extras = (Ifree, Itied, tiedtoparam, tiedfactor, bounds,
                            doubletindx, doubletpair, linewaves)

        return parameters, parameter_extras

    def _populate_linemodel(self, linemodel, initial_guesses, param_bounds):
        """Population an input linemodel with initial guesses and parameter bounds,
        taking into account fixed parameters.

        """
        # Set initial values and bounds.
        for iparam, param in enumerate(linemodel['param_name']):
            if param in initial_guesses.keys():
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0 # always set fixed parameter to zero
                else:
                    linemodel['initial'][iparam] = initial_guesses[param]
                    if param in param_bounds.keys():
                        linemodel['bounds'][iparam] = param_bounds[param]
            else:
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0
                else:
                    linemodel['initial'][iparam] = 1.0

            # Check bounds for free parameters but do not crash.
            if linemodel['fixed'][iparam] == False and linemodel['tiedtoparam'][iparam] == -1:
                toosml = linemodel['initial'][iparam] < linemodel['bounds'][iparam, 0]
                toobig = linemodel['initial'][iparam] > linemodel['bounds'][iparam, 1]
                if toosml:
                    errmsg = 'Initial parameter {} is outside its bound, {:.2f} < {:.2f}.'.format(
                        param, linemodel['initial'][iparam], linemodel['bounds'][iparam, 0])
                    self.log.warning(errmsg)
                    #raise ValueError(errmsg)
                    linemodel['initial'][iparam] = linemodel['bounds'][iparam, 0]
                if toobig:
                    errmsg = 'Initial parameter {} is outside its bound, {:.2f} > {:.2f}.'.format(
                        param, linemodel['initial'][iparam], linemodel['bounds'][iparam, 1])
                    self.log.warning(errmsg)
                    #raise ValueError(errmsg)
                    linemodel['initial'][iparam] = linemodel['bounds'][iparam, 1]
                    
        # Now loop back through and ensure that tied relationships are enforced.
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        if len(Itied) > 0:
            for iparam, param in enumerate(linemodel['param_name'][Itied]):
                tieindx = linemodel['tiedtoparam'][Itied[iparam]]
                tiefactor = linemodel['tiedfactor'][Itied[iparam]]
                #self.log.info('{} tied to {} with factor {:.4f}'.format(param, linemodel[tieindx]['param_name'], tiefactor))
                linemodel['initial'][Itied[iparam]] = linemodel[tieindx]['initial'] * tiefactor

        linemodel['value'] = linemodel['initial'] # copy

    def _optimize(self, linemodel, emlinewave, emlineflux, weights, 
                  redshift, resolution_matrix, camerapix, debug=False):
        """Wrapper to call the least-squares minimization given a linemodel.

        """
        from scipy.optimize import least_squares
        from fastspecfit.emlines import _objective_function

        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, doubletpair, \
                     linewaves) = self._linemodel_to_parameters(linemodel)
        self.log.debug('Optimizing {} free parameters'.format(len(Ifree)))

        farg = (emlinewave, emlineflux, weights, redshift, self.log10wave, 
                resolution_matrix, camerapix, parameters, ) + \
                (Ifree, Itied, tiedtoparam, tiedfactor, doubletindx, 
                 doubletpair, linewaves)

        fit_info = least_squares(_objective_function, parameters[Ifree],
                                 args=farg, max_nfev=self.maxiter, 
                                 xtol=self.accuracy, 
                                 #method='lm')
                                 #verbose=2,
                                 tr_solver='lsmr', tr_options={'regularize': True},
                                 method='trf', bounds=tuple(zip(*bounds)))
        parameters[Ifree] = fit_info.x

        # Conditions for dropping a parameter (all parameters, not just those
        # being fitted):
        # --negative amplitude or sigma
        # --parameter at its default value (fit failed, right??)
        # --parameter outside its bounds [should never be needed if method='trf']
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
        notfixed = np.logical_not(linemodel['fixed'])

        drop1 = np.hstack((lineamps < 0, np.zeros(len(linevshifts), bool), linesigmas <= 0)) * notfixed
        
        drop2 = np.zeros(len(parameters), bool)
        drop2[Ifree] = parameters[Ifree] == linemodel['value'][Ifree] # want 'value' here not 'initial'
        drop2 *= notfixed
        
        drop3 = np.zeros(len(parameters), bool)
        drop3[Ifree] = np.logical_or(parameters[Ifree] < linemodel['bounds'][Ifree, 0], 
                                     parameters[Ifree] > linemodel['bounds'][Ifree, 1])
        drop3 *= notfixed
        
        self.log.debug('Dropping {} negative amplitudes or line-widths.'.format(np.sum(drop1)))
        self.log.debug('Dropping {} parameters which were not optimized.'.format(np.sum(drop2)))
        self.log.debug('Dropping {} parameters which are out-of-bounds.'.format(np.sum(drop3)))
        Idrop = np.where(np.logical_or.reduce((drop1, drop2, drop3)))[0]

        if debug:
            pass

        if len(Idrop) > 0:
            self.log.debug('  Dropping {} unique parameters.'.format(len(Idrop)))
            parameters[Idrop] = 0.0

        # apply tied constraints
        if len(Itied) > 0:
            for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
                parameters[I] = parameters[indx] * factor

        # Now loop back through and drop Broad balmer lines that:
        #   (1) are narrower than their narrow-line counterparts;
        #   (2) have a narrow line whose amplitude is smaller than that of the broad line
        #      --> Deprecated! main-dark-32303-39628176678192981 is an example
        #          of an object where there's a broad H-alpha line but no other
        #          forbidden lines!
        
        out_linemodel = linemodel.copy()
        out_linemodel['value'] = parameters
        out_linemodel.meta['nfev'] = fit_info['nfev']

        if False:
            bestfit = self.bestfit(out_linemodel, redshift, emlinewave, resolution_matrix, camerapix)
            import matplotlib.pyplot as plt
            plt.clf()
            plt.plot(emlinewave, emlineflux)
            plt.plot(emlinewave, bestfit)
            #plt.xlim(5800, 6200)
            #plt.xlim(6600, 6950)
            plt.xlim(5050, 5120)
            #plt.xlim(8850, 9050)
            plt.savefig('junk.png')
        
        return out_linemodel

    def chi2(self, linemodel, emlinewave, emlineflux, emlineivar, emlinemodel,
             continuum_model=None, return_dof=False):
        """Compute the reduced chi^2."""

        nfree = np.count_nonzero((linemodel['fixed'] == False) * (linemodel['tiedtoparam'] == -1))
        dof = np.count_nonzero(emlineivar > 0) - nfree

        if dof > 0:
            if continuum_model is None:
                model = emlinemodel
            else:
                model = emlinemodel + continuum_model
            chi2 = np.sum(emlineivar * (emlineflux - model)**2) / dof
        else:
            chi2 = self.chi2_default
            
        if return_dof:
            return chi2, dof
        else:
            return chi2

    def bestfit(self, linemodel, redshift, emlinewave, resolution_matrix, camerapix):
        """Construct the best-fitting emission-line spectrum from a linemodel."""

        from fastspecfit.emlines import build_emline_model

        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, \
                     doubletpair, linewaves) = self._linemodel_to_parameters(linemodel)

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line

        # doublets
        lineamps[doubletindx] *= lineamps[doubletpair]

        emlinemodel = build_emline_model(self.log10wave, redshift, lineamps, 
                                         linevshifts, linesigmas, linewaves, 
                                         emlinewave, resolution_matrix,
                                         camerapix)

        return emlinemodel

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table, redshift=None):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (used for QA and elsewhere).

        """
        from fastspecfit.emlines import build_emline_model

        if redshift is None:
            redshift = fastspecfit_table['CONTINUUM_Z']
        
        linewaves = self.linetable['restwave'].data

        parameters = [fastspecfit_table[param.upper()] for param in self.param_names]

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line    

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        lineamps[self.doubletindx] *= lineamps[self.doubletpair]

        emlinemodel = build_emline_model(self.log10wave, redshift, lineamps, 
                                         linevshifts, linesigmas, linewaves, 
                                         specwave, specres, None)

        return emlinemodel

    def emline_specfit(self, data, result, continuummodel, smooth_continuum, 
                       synthphot=True, broadlinefit=True, verbose=False):
        """Perform the fit minimization / chi2 minimization.

        Parameters
        ----------
        data
        continuummodel
        smooth_continuum
        synthphot
        verbose
        broadlinefit

        Returns
        -------
        results
        modelflux
     
        """
        from fastspecfit.util import ivar2var

        tall = time.time()

        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        redshift = data['zredrock']
        emlinewave = np.hstack(data['wave'])
        oemlineivar = np.hstack(data['ivar'])/data['apcorr']**2
        specflux = np.hstack(data['flux'])*data['apcorr']
        resolution_matrix = data['res']
        camerapix = data['camerapix']

        continuummodelflux = np.hstack(continuummodel)
        smoothcontinuummodelflux = np.hstack(smooth_continuum)
        emlineflux = specflux - continuummodelflux - smoothcontinuummodelflux

        emlineivar = np.copy(oemlineivar)
        emlinevar, emlinegood = ivar2var(emlineivar, clip=1e-3)
        emlinebad = np.logical_not(emlinegood)

        # This is a (dangerous???) hack.
        if np.sum(emlinebad) > 0:
            emlineivar[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineivar[emlinegood])
            emlineflux[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineflux[emlinegood]) # ???

        weights = np.sqrt(emlineivar)

        # Build all the emission-line models for this object.
        final_linemodel, final_linemodel_nobroad, initial_linemodel, initial_linemodel_nobroad = \
            self.build_linemodels(redshift, wavelims=(np.min(emlinewave)+5, np.max(emlinewave)-5),
                                  verbose=False)

        # Get initial guesses on the parameters and populate the two "initial"
        # linemodels; the "final" linemodels will be initialized with the
        # best-fitting parameters from the initial round of fitting.
        initial_guesses, param_bounds = self._initial_guesses_and_bounds(data, emlinewave, emlineflux)

        for linemodel in [initial_linemodel, initial_linemodel_nobroad]:
            self._populate_linemodel(linemodel, initial_guesses, param_bounds)

        # Initial fit - initial_linemodel_nobroad
        t0 = time.time()
        initfit = self._optimize(initial_linemodel_nobroad, emlinewave, emlineflux, 
                                 weights, redshift, resolution_matrix, camerapix, 
                                 debug=False)
        initmodel = self.bestfit(initfit, redshift, emlinewave, resolution_matrix, camerapix)
        initchi2 = self.chi2(initfit, emlinewave, emlineflux, emlineivar, initmodel)
        nfree = np.sum((initfit['fixed'] == False) * (initfit['tiedtoparam'] == -1))
        self.log.info('Initial line-fitting with {} free parameters took {:.2f} seconds [niter={}, rchi2={:.4f}].'.format(
            nfree, time.time()-t0, initfit.meta['nfev'], initchi2))

        ## Now try adding bround Balmer and helium lines and see if we improve
        ## the chi2. First, do we have enough pixels around Halpha and Hbeta to
        ## do this test?
        #broadlinepix = []
        #for icam in np.arange(len(data['cameras'])):
        #    pixoffset = int(np.sum(data['npixpercamera'][:icam]))
        #    for linename, linepix in zip(data['linename'][icam], data['linepix'][icam]):
        #        if linename == 'halpha_broad' or linename == 'hbeta_broad' or linename == 'hgamma_broad':
        #            broadlinepix.append(linepix + pixoffset)
        #            #print(data['wave'][icam][linepix])

        # Require minimum XX pixels.
        if broadlinefit:# and len(broadlinepix) > 0 and len(np.hstack(broadlinepix)) > 10:
            t0 = time.time()
            broadfit = self._optimize(initial_linemodel, emlinewave, emlineflux, weights, 
                                      redshift, resolution_matrix, camerapix, debug=False)
            broadmodel = self.bestfit(broadfit, redshift, emlinewave, resolution_matrix, camerapix)
            broadchi2 = self.chi2(broadfit, emlinewave, emlineflux, emlineivar, broadmodel)
            nfree = np.sum((broadfit['fixed'] == False) * (broadfit['tiedtoparam'] == -1))
            self.log.info('Second (broad) line-fitting with {} free parameters took {:.2f} seconds [niter={}, rchi2={:.4f}].'.format(
                nfree, time.time()-t0, broadfit.meta['nfev'], broadchi2))

            ## Compare chi2 just in and around the broad lines.
            #broadlinepix = np.hstack(broadlinepix)
            #I = ((self.fit_linetable[self.fit_linetable['inrange']]['zwave'] > np.min(emlinewave[broadlinepix])) *
            #     (self.fit_linetable[self.fit_linetable['inrange']]['zwave'] < np.max(emlinewave[broadlinepix])))
            #Iline = initfit[np.isin(broadfit['linename'], self.fit_linetable[self.fit_linetable['inrange']][I]['name'])]
            #Bline = broadfit[np.isin(broadfit['linename'], self.fit_linetable[self.fit_linetable['inrange']][I]['name'])]
            #dof_init = np.count_nonzero(emlineivar[broadlinepix] > 0) - np.count_nonzero((Iline['fixed'] == False) * (Iline['tiedtoparam'] == -1))
            #dof_broad = np.count_nonzero(emlineivar[broadlinepix] > 0) - np.count_nonzero((Bline['fixed'] == False) * (Bline['tiedtoparam'] == -1))
            #if dof_init == 0 or dof_broad == 0:
            #    errmsg = 'Number of degrees of freedom should never be zero: dof_init={}, dof_free={}'.format(
            #        dof_init, dof_broad)
            #    self.log.critical(errmsg)
            #    raise ValueError(errmsg)
            #
            #linechi2_init = np.sum(emlineivar[broadlinepix] * (emlineflux[broadlinepix] - initmodel[broadlinepix])**2) / dof_init
            #linechi2_broad = np.sum(emlineivar[broadlinepix] * (emlineflux[broadlinepix] - broadmodel[broadlinepix])**2) / dof_broad
            #self.log.info('Chi2 with broad lines = {:.5f} and without broad lines = {:.5f} [delta-chi2={:.5f}]'.format(
            #    linechi2_broad, linechi2_init, linechi2_init - linechi2_broad))

            linechi2_broad, linechi2_init = broadchi2, initchi2

            self.log.info('Chi2 with broad lines = {:.5f} and without broad lines = {:.5f} [delta-chi2={:.5f}]'.format(
                linechi2_broad, linechi2_init, linechi2_init - linechi2_broad))

            if linechi2_broad > linechi2_init:
                bestfit = initfit
            else:
                bestfit = broadfit
        else:
            self.log.info('Skipping broad-line fitting.')

            bestfit = initfit
            linechi2_broad, linechi2_init = 1e6, initchi2
            
            #if broadlinefit:
            #    self.log.info('Too few pixels centered on candidate broad emission lines.')
            #else:
            #    self.log.info('Skipping broad-line fitting.')
            #linechi2_init, linechi2_broad = 0.0, 0.0

        # Finally, one more fitting loop with all the line-constraints relaxed
        # but starting from the previous best-fitting values.
        if linechi2_init < (linechi2_broad + 1):
            linemodel = final_linemodel_nobroad
        else:
            linemodel = final_linemodel

        # Populate the new linemodel being careful to handle the fact that the
        # "tied" relationships are very different between the initial and final
        # linemodels.
        linemodel['bounds'] = bestfit['bounds']

        Ifree = np.where(linemodel['fixed'] == False)[0]
        for I in Ifree:
            # copy initial values
            if bestfit['initial'][I] != 0:
                linemodel['initial'][I] = bestfit['initial'][I]
            # copy best-fit values
            if bestfit['value'][I] != 0:
                linemodel['value'][I] = bestfit['value'][I]
            else:
                if bestfit['tiedtoparam'][I] != -1:
                    linemodel['value'][I] = bestfit['value'][bestfit['tiedtoparam'][I]]
                else:
                    linemodel['value'][I] = linemodel['initial'][I]

        # Are the broad and narrow lines swapped? If so, swap them here.
        if not linemodel[linemodel['param_name'] == 'halpha_broad_sigma']['fixed'] and \
          (linemodel[linemodel['param_name'] == 'halpha_broad_sigma']['value'] > 0) and \
          (linemodel[linemodel['param_name'] == 'halpha_broad_sigma']['value'] < linemodel[linemodel['param_name'] == 'halpha_sigma']['value']):
            for linename in self.fit_linetable[self.fit_linetable['isbalmer'] * self.fit_linetable['isbroad']]['name']:
                if not linemodel[linemodel['param_name'] == '{}_sigma'.format(linename)]['fixed']:
                    sigma_broad = linemodel[linemodel['param_name'] == '{}_sigma'.format(linename)]['value']
                    sigma_narrow = linemodel[linemodel['param_name'] == '{}_sigma'.format(linename).replace('_broad', '')]['value']
                    #print(linename, sigma_broad[0], sigma_narrow[0])
                    linemodel['value'][linemodel['param_name'] == '{}_sigma'.format(linename)] = sigma_narrow
                    linemodel['value'][linemodel['param_name'] == '{}_sigma'.format(linename).replace('_broad', '')] = sigma_broad

        # Tied parameters can have initial values of zero if they are fixed
        # (e.g., broad emission lines) but nothing else.
        I = (linemodel['value'][Ifree] == 0) * (linemodel['tiedtoparam'][Ifree] == -1)
        if np.any(I):
            errmsg = 'Initial values should never be zero [targetid={}]!'.format(data['targetid'])
            self.log.critical(errmsg)
            raise ValueError(errmsg)

        #linemodel[linemodel['linename'] == 'halpha']
        #B = np.where(['ne' in param for param in self.param_names])[0]
        #B = np.where(['broad' in param for param in self.param_names])[0]

        t0 = time.time()
        finalfit = self._optimize(linemodel, emlinewave, emlineflux, weights, 
                                  redshift, resolution_matrix, camerapix, 
                                  debug=False)
        finalmodel = self.bestfit(finalfit, redshift, emlinewave, resolution_matrix, camerapix)
        finalchi2 = self.chi2(finalfit, emlinewave, emlineflux, emlineivar, finalmodel)
        nfree = np.sum((finalfit['fixed'] == False) * (finalfit['tiedtoparam'] == -1))
        self.log.info('Final line-fitting with {} free parameters took {:.2f} seconds [niter={}, rchi2={:.4f}].'.format(
            nfree, time.time()-t0, finalfit.meta['nfev'], finalchi2))

        # Residual spectrum with no emission lines.
        specflux_nolines = specflux - finalmodel

        #import matplotlib.pyplot as plt
        #W = (emlinewave>6560)*(emlinewave<6660)
        #plt.clf()
        #plt.plot(emlinewave[W], emlineflux[W], color='gray')
        #plt.plot(emlinewave[W], finalmodel[W], color='orange', alpha=0.7)
        ##plt.plot(emlinewave[W], specflux[W], color='gray')
        ##plt.plot(emlinewave[W], specflux_nolines[W], color='orange', alpha=0.7)
        #plt.savefig('desi-users/ioannis/tmp/junk2.png')

        # Initialize the output table; see init_fastspecfit for the data model.
        result['RCHI2'] = finalchi2
        result['LINERCHI2_BROAD'] = linechi2_broad
        result['DELTA_LINERCHI2'] = linechi2_init - linechi2_broad

        # Now fill the output table.
        self._populate_emtable(result, finalfit, finalmodel, emlinewave, emlineflux,
                               emlineivar, oemlineivar, specflux_nolines, redshift)

        # Build the model spectra.
        emmodel = np.hstack(self.emlinemodel_bestfit(data['wave'], data['res'], result, redshift=redshift))

        # As a consistency check, make sure that the emission-line spectrum
        # rebuilt from the final table is not (very) different from the one
        # based on the best-fitting model parameters.
        #assert(np.all(np.isclose(emmodel, emlinemodel, rtol=1e-4)))
            
        #import matplotlib.pyplot as plt
        #plt.clf()
        #plt.plot(emlinewave, np.hstack(emmodel)-emlinemodel)
        ##plt.plot(emlinewave, emlinemodel)
        #plt.savefig('desi-users/ioannis/tmp/junk.png')

        # I believe that all the elements of the coadd_wave vector are contained
        # within one or more of the per-camera wavelength vectors, and so we
        # should be able to simply map our model spectra with no
        # interpolation. However, because of round-off, etc., it's probably
        # easiest to use np.interp.

        # package together the final output models for writing; assume constant
        # dispersion in wavelength!
        minwave, maxwave, dwave = np.min(data['coadd_wave']), np.max(data['coadd_wave']), np.diff(data['coadd_wave'][:2])[0]
        minwave = float(int(minwave * 1000) / 1000)
        maxwave = float(int(maxwave * 1000) / 1000)
        dwave = float(int(dwave * 1000) / 1000)
        npix = int((maxwave-minwave)/dwave)+1
        modelwave = minwave + dwave * np.arange(npix)

        modelspectra = Table()
        # all these header cards need to be 2-element tuples (value, comment),
        # otherwise io.write_fastspecfit will crash
        modelspectra.meta['NAXIS1'] = (npix, 'number of pixels')
        modelspectra.meta['NAXIS2'] = (npix, 'number of models')
        modelspectra.meta['NAXIS3'] = (npix, 'number of objects')
        modelspectra.meta['BUNIT'] = ('10**-17 erg/(s cm2 Angstrom)', 'flux unit')
        modelspectra.meta['CUNIT1'] = ('Angstrom', 'wavelength unit')
        modelspectra.meta['CTYPE1'] = ('WAVE', 'type of axis')
        modelspectra.meta['CRVAL1'] = (minwave, 'wavelength of pixel CRPIX1 (Angstrom)')
        modelspectra.meta['CRPIX1'] = (0, '0-indexed pixel number corresponding to CRVAL1')
        modelspectra.meta['CDELT1'] = (dwave, 'pixel size (Angstrom)')
        modelspectra.meta['DC-FLAG'] = (0, '0 = linear wavelength vector')
        modelspectra.meta['AIRORVAC'] = ('vac', 'wavelengths in vacuum (vac)')

        modelcontinuum = np.interp(modelwave, emlinewave, continuummodelflux).reshape(1, npix)
        modelsmoothcontinuum = np.interp(modelwave, emlinewave, smoothcontinuummodelflux).reshape(1, npix)
        modelemspectrum = np.interp(modelwave, emlinewave, emmodel).reshape(1, npix)
        
        modelspectra.add_column(Column(name='CONTINUUM', dtype='f4', data=modelcontinuum))
        modelspectra.add_column(Column(name='SMOOTHCONTINUUM', dtype='f4', data=modelsmoothcontinuum))
        modelspectra.add_column(Column(name='EMLINEMODEL', dtype='f4', data=modelemspectrum))

        # Finally, optionally synthesize photometry.
        if synthphot:
            modelflux = modelcontinuum[0, :] + modelsmoothcontinuum[0, :] + modelemspectrum[0, :]
            self._synthphot_spectrum(data, result, modelwave, modelflux)

        log.info('Emission-line fitting took {:.2f} seconds.'.format(time.time()-tall))

        return modelspectra

    def _populate_emtable(self, result, finalfit, finalmodel, emlinewave, emlineflux,
                          emlineivar, oemlineivar, specflux_nolines, redshift):
        """Populate the output table with the emission-line measurements.

        """
        from scipy.stats import sigmaclip

        for param in finalfit:
            val = param['value']
            # special case the tied doublets
            if param['param_name'] == 'oii_doublet_ratio':
                result['OII_DOUBLET_RATIO'] = val
                result['OII_3726_AMP'] = val * finalfit[param['doubletpair']]['value']
            elif param['param_name'] == 'sii_doublet_ratio':
                result['SII_DOUBLET_RATIO'] = val
                result['SII_6731_AMP'] = val * finalfit[param['doubletpair']]['value']
            elif param['param_name'] == 'mgii_doublet_ratio':
                result['MGII_DOUBLET_RATIO'] = val
                result['MGII_2796_AMP'] = val * finalfit[param['doubletpair']]['value']
            else:
                result[param['param_name'].upper()] = val

        # get continuum fluxes, EWs, and upper limits
        narrow_sigmas, broad_sigmas, uv_sigmas = [], [], []
        narrow_redshifts, broad_redshifts, uv_redshifts = [], [], []
        for oneline in self.fit_linetable[self.fit_linetable['inrange']]:

            linename = oneline['name'].upper()
            linez = redshift + result['{}_VSHIFT'.format(linename)] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)
            linesigma = result['{}_SIGMA'.format(linename)] # [km/s]

            # if the line was dropped, use a default sigma value
            if linesigma == 0:
                if oneline['isbroad']:
                    if oneline['isbalmer']:
                        linesigma = self.limitsigma_narrow
                    else:
                        linesigma = self.limitsigma_broad
                else:
                    linesigma = self.limitsigma_broad

            linesigma_ang = linesigma * linezwave / C_LIGHT    # [observed-frame Angstrom]
            #log10sigma = linesigma / C_LIGHT / np.log(10)     # line-width [log-10 Angstrom]            

            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next line.
            lineindx = np.where((emlinewave >= (linezwave - 3.0*linesigma_ang)) *
                                (emlinewave <= (linezwave + 3.0*linesigma_ang)))[0]
            
            if len(lineindx) > 0 and np.sum(oemlineivar[lineindx] == 0) / len(lineindx) > 0.3: # use original ivar
                result['{}_AMP'.format(linename)] = 0.0
                result['{}_VSHIFT'.format(linename)] = 0.0
                result['{}_SIGMA'.format(linename)] = 0.0
            else:
                # number of pixels, chi2, and boxcar integration
                lineindx = np.where((emlinewave >= (linezwave - 3.0*linesigma_ang)) *
                                    (emlinewave <= (linezwave + 3.0*linesigma_ang)) *
                                    (emlineivar > 0))[0]
    
                # can happen if sigma is very small (depending on the wavelength)
                #if (linezwave > np.min(emlinewave)) * (linezwave < np.max(emlinewave)) * len(lineindx) > 0 and len(lineindx) <= 3: 
                if (linezwave > np.min(emlinewave)) * (linezwave < np.max(emlinewave)) * (len(lineindx) <= 3):
                    dwave = emlinewave - linezwave
                    lineindx = np.argmin(np.abs(dwave))
                    if dwave[lineindx] > 0:
                        pad = np.array([-2, -1, 0, +1])
                    else:
                        pad = np.array([-1, 0, +1, +2])
    
                    # check to make sure we don't hit the edges
                    if (lineindx-pad[0]) < 0 or (lineindx+pad[-1]) >= len(emlineivar):
                        lineindx = np.array([])
                    else:
                        lineindx += pad
                        # the padded pixels can have ivar==0
                        good = oemlineivar[lineindx] > 0 # use the original ivar
                        lineindx = lineindx[good]
    
                npix = len(lineindx)
                result['{}_NPIX'.format(linename)] = npix
    
                if npix > 3: # magic number: required at least XX unmasked pixels centered on the line
    
                    if np.any(emlineivar[lineindx] == 0):
                        errmsg = 'Ivar should never be zero within an emission line!'
                        self.log.critical(errmsg)
                        raise ValueError(errmsg)
                        
                    # boxcar integration of the flux; should we weight by the line-profile???
                    boxflux = np.sum(emlineflux[lineindx])                
                    boxflux_ivar = 1 / np.sum(1 / emlineivar[lineindx])
    
                    result['{}_BOXFLUX'.format(linename)] = boxflux # * u.erg/(u.second*u.cm**2)
                    result['{}_BOXFLUX_IVAR'.format(linename)] = boxflux_ivar # * u.second**2*u.cm**4/u.erg**2
                    
                    # Get the uncertainty in the line-amplitude based on the scatter
                    # in the pixel values from the emission-line subtracted
                    # spectrum.
                    amp_sigma = np.diff(np.percentile(specflux_nolines[lineindx], [25, 75]))[0] / 1.349 # robust sigma
                    #clipflux, _, _ = sigmaclip(specflux_nolines[lineindx], low=3, high=3)
                    #amp_sigma = np.std(clipflux)
                    if amp_sigma > 0:
                        result['{}_AMP_IVAR'.format(linename)] = 1 / amp_sigma**2 # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
    
                    #if np.isinf(result['{}_AMP_IVAR'.format(linename)]):
                    #    pdb.set_trace()
    
                    # require amp > 0 (line not dropped) to compute the flux and chi2
                    if result['{}_AMP'.format(linename)] > 0:
    
                        # get the emission-line flux
                        linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang # * u.Angstrom
                        result['{}_FLUX'.format(linename)] = result['{}_AMP'.format(linename)] * linenorm
            
                        #result['{}_FLUX_IVAR'.format(linename)] = result['{}_AMP_IVAR'.format(linename)] / linenorm**2
                        #weight = np.exp(-0.5 * np.log10(emlinewave/linezwave)**2 / log10sigma**2)
                        #weight = (weight / np.max(weight)) > 1e-3
                        #result['{}_FLUX_IVAR'.format(linename)] = 1 / np.sum(1 / emlineivar[weight])
                        result['{}_FLUX_IVAR'.format(linename)] = boxflux_ivar # * u.second**2*u.cm**4/u.erg**2
    
                        dof = npix - 3 # ??? [redshift, sigma, and amplitude]
                        chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-finalmodel[lineindx])**2) / dof
    
                        result['{}_CHI2'.format(linename)] = chi2
    
                        # keep track of sigma and z but only using XX-sigma lines
                        linesnr = result['{}_AMP'.format(linename)] * np.sqrt(result['{}_AMP_IVAR'.format(linename)])
                        #print(linename, result['{}_AMP'.format(linename)], amp_sigma, linesnr)
                        if linesnr > 2:
                            if oneline['isbroad']: # includes UV and broad Balmer lines
                                if oneline['isbalmer']:
                                    broad_sigmas.append(linesigma)
                                    broad_redshifts.append(linez)
                                else:
                                    uv_sigmas.append(linesigma)
                                    uv_redshifts.append(linez)
                            else:
                                narrow_sigmas.append(linesigma)
                                narrow_redshifts.append(linez)
    
                # next, get the continuum, the inverse variance in the line-amplitude, and the EW
                indxlo = np.where((emlinewave > (linezwave - 10*linesigma * linezwave / C_LIGHT)) *
                                  (emlinewave < (linezwave - 3.*linesigma * linezwave / C_LIGHT)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indxhi = np.where((emlinewave < (linezwave + 10*linesigma * linezwave / C_LIGHT)) *
                                  (emlinewave > (linezwave + 3.*linesigma * linezwave / C_LIGHT)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indx = np.hstack((indxlo, indxhi))
    
                if len(indx) >= 3: # require at least XX pixels to get the continuum level
                    #_, cmed, csig = sigma_clipped_stats(specflux_nolines[indx], sigma=3.0)
                    clipflux, _, _ = sigmaclip(specflux_nolines[indx], low=3, high=3)
                    # corner case: if a portion of a camera is masked
                    if len(clipflux) > 0:
                        #cmed, csig = np.mean(clipflux), np.std(clipflux)
                        cmed = np.median(clipflux)
                        csig = np.diff(np.percentile(clipflux, [25, 75])) / 1.349 # robust sigma
                        if csig > 0:
                            civar = (np.sqrt(len(indx)) / csig)**2
                        else:
                            civar = 0.0
                    else:
                        cmed, civar = 0.0, 0.0
    
                    result['{}_CONT'.format(linename)] = cmed # * u.erg/(u.second*u.cm**2*u.Angstrom)
                    result['{}_CONT_IVAR'.format(linename)] = civar # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
    
                if result['{}_CONT'.format(linename)] != 0.0 and result['{}_CONT_IVAR'.format(linename)] != 0.0:
                    factor = 1 / ((1 + redshift) * result['{}_CONT'.format(linename)]) # --> rest frame
                    ew = result['{}_FLUX'.format(linename)] * factor # rest frame [A]
                    ewivar = result['{}_FLUX_IVAR'.format(linename)] / factor**2
    
                    # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                    fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar) # * u.erg/(u.second*u.cm**2)
                    ewlimit = fluxlimit * factor
    
                    result['{}_EW'.format(linename)] = ew
                    result['{}_EW_IVAR'.format(linename)] = ewivar
                    result['{}_FLUX_LIMIT'.format(linename)] = fluxlimit 
                    result['{}_EW_LIMIT'.format(linename)] = ewlimit

            if 'debug' in self.log.name:
                for col in ('VSHIFT', 'SIGMA', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    self.log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)]))
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'BOXFLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT', 'EW_LIMIT'):
                    self.log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)]))
                print()
                #self.log.debug(' ')
    
            ## simple QA
            #if linename == 'OIII_5007':
            #    import matplotlib.pyplot as plt
            #    _indx = np.arange(indx[-1]-indx[0])+indx[0]
            #    # continuum bandpasses and statistics
            #    plt.clf()
            #    plt.plot(emlinewave[_indx], specflux_nolines[_indx], color='gray')
            #    plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
            #    plt.axhline(y=cmed, color='k')
            #    plt.axhline(y=cmed+1/np.sqrt(civar), color='k', ls='--')
            #    plt.axhline(y=cmed-1/np.sqrt(civar), color='k', ls='--')
            #    plt.savefig('desi-users/ioannis/tmp/junk.png')
            #
            #    # emission-line integration
            #    plt.clf()
            #    plt.plot(emlinewave[_indx], emlineflux[_indx], color='gray')
            #    plt.plot(emlinewave[_indx], finalmodel[_indx], color='red')
            #    #plt.plot(emlinewave[_indx], specflux_nolines[_indx], color='orange', alpha=0.5)
            #    plt.axvline(x=emlinewave[lineindx[0]], color='blue')
            #    plt.axvline(x=emlinewave[lineindx[-1]], color='blue')
            #    plt.axhline(y=0, color='k', ls='--')
            #    plt.axhline(y=amp_sigma, color='k', ls='--')
            #    plt.axhline(y=2*amp_sigma, color='k', ls='--')
            #    plt.axhline(y=3*amp_sigma, color='k', ls='--')
            #    plt.axhline(y=result['{}_AMP'.format(linename)], color='k', ls='-')
            #    plt.savefig('desi-users/ioannis/tmp/junk2.png')
            #    pdb.set_trace()

        # Clean up the doublets whose amplitudes were tied in the fitting since
        # they may have been zeroed out in the clean-up, above.
        if result['OIII_5007_AMP'] == 0.0 and result['OIII_5007_NPIX'] > 0:
            result['OIII_4959_AMP'] = 0.0
            result['OIII_4959_FLUX'] = 0.0
            result['OIII_4959_EW'] = 0.0
        if result['NII_6584_AMP'] == 0.0 and result['NII_6584_NPIX'] > 0:
            result['NII_6548_AMP'] = 0.0
            result['NII_6548_FLUX'] = 0.0
            result['NII_6548_EW'] = 0.0
        if result['OII_7320_AMP'] == 0.0 and result['OII_7320_NPIX'] > 0:
            result['OII_7330_AMP'] = 0.0
            result['OII_7330_FLUX'] = 0.0
            result['OII_7330_EW'] = 0.0
        if result['MGII_2796_AMP'] == 0.0 and result['MGII_2803_AMP'] == 0.0:
            result['MGII_DOUBLET_RATIO'] = 0.0
        if result['OII_3726_AMP'] == 0.0 and result['OII_3729_AMP'] == 0.0:
            result['OII_DOUBLET_RATIO'] = 0.0
        if result['SII_6716_AMP'] == 0.0 and result['SII_6731_AMP'] == 0.0:
            result['SII_DOUBLET_RATIO'] = 0.0

        if 'debug' in self.log.name:
            for col in ('MGII_DOUBLET_RATIO', 'OII_DOUBLET_RATIO', 'SII_DOUBLET_RATIO'):
                self.log.debug('{}: {:.4f}'.format(col, result[col]))
            #self.log.debug(' ')
            print()

        # get the average emission-line redshifts and velocity widths
        if len(narrow_redshifts) > 0:
            result['NARROW_Z'] = np.median(narrow_redshifts)
            result['NARROW_SIGMA'] = np.median(narrow_sigmas) # * u.kilometer / u.second
            #result['NARROW_Z_ERR'] = np.std(narrow_redshifts)
            #result['NARROW_SIGMA_ERR'] = np.std(narrow_sigmas)
        else:
            result['NARROW_Z'] = redshift
            
        if len(broad_redshifts) > 0:
            result['BROAD_Z'] = np.median(broad_redshifts)
            result['BROAD_SIGMA'] = np.median(broad_sigmas) # * u.kilometer / u.second
            #result['BROAD_Z_ERR'] = np.std(broad_redshifts)
            #result['BROAD_SIGMA_ERR'] = np.std(broad_sigmas)
        else:
            result['BROAD_Z'] = redshift
            
        if len(uv_redshifts) > 0:
            result['UV_Z'] = np.median(uv_redshifts)
            result['UV_SIGMA'] = np.median(uv_sigmas) # * u.kilometer / u.second
            #result['UV_Z_ERR'] = np.std(uv_redshifts)
            #result['UV_SIGMA_ERR'] = np.std(uv_sigmas)
        else:
            result['UV_Z'] = redshift

        # fragile
        if 'debug' in self.log.name:
            for line in ('NARROW', 'BROAD', 'UV'):
                self.log.debug('{}_Z: {:.12f}'.format(line, result['{}_Z'.format(line)]))
                self.log.debug('{}_SIGMA: {:.3f}'.format(line, result['{}_SIGMA'.format(line)]))

    def _synthphot_spectrum(self, data, result, modelwave, modelflux):
        """Synthesize photometry from the best-fitting model (continuum+emission lines).

        """
        if data['photsys'] == 'S':
            filters = self.decam
        else:
            filters = self.bassmzls

        # Pad (simply) in wavelength...
        padflux, padwave = filters.pad_spectrum(modelflux, modelwave, method='edge')
        synthmaggies = filters.get_ab_maggies(padflux / self.fluxnorm, padwave)
        synthmaggies = synthmaggies.as_array().view('f8')
        model_synthphot = self.parse_photometry(self.synth_bands, maggies=synthmaggies,
                                                nanomaggies=False,
                                                lambda_eff=filters.effective_wavelengths.value)

        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_{}'.format(band.upper())] = data['synthphot']['nanomaggies'][iband] # * 'nanomaggies'
            #result['FLUX_SYNTH_IVAR_{}'.format(band.upper())] = data['synthphot']['nanomaggies_ivar'][iband]
        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_MODEL_{}'.format(band.upper())] = model_synthphot['nanomaggies'][iband] # * 'nanomaggies'

        ## measure DN(4000) without the emission lines
        #if False:
        #    dn4000_nolines, _ = self.get_dn4000(emlinewave, specflux_nolines, redshift=redshift)
        #    result['DN4000_NOLINES'] = dn4000_nolines
    
    def qa_fastspec(self, data, fastspec, metadata, coadd_type='healpix',
                    spec_wavelims=(3550, 9900), outprefix=None, outdir=None):
        """QA plot the emission-line spectrum and best-fitting model.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        from fastspecfit.util import ivar2var

        sns.set(context='talk', style='ticks', font_scale=1.1)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['dodgerblue', 'darkseagreen', 'orangered']]
        col2 = [colors.to_hex(col) for col in ['darkblue', 'darkgreen', 'darkred']]
        col3 = [colors.to_hex(col) for col in ['blue', 'green', 'red']]

        photcol1 = colors.to_hex('darkblue') # 'darkgreen', 'darkred', 'dodgerblue', 'darkseagreen', 'orangered']]
        
        if outdir is None:
            outdir = '.'
        if outprefix is None:
            outprefix = 'fastspec'

        if metadata['PHOTSYS'] == 'S':
            filters = self.decam
            allfilters = self.decamwise
        else:
            filters = self.bassmzls
            allfilters = self.bassmzlswise
    
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
        elif coadd_type == 'custom':
            title = 'Survey/Program/HealPix: {}/{}/{}, TargetID: {}'.format(
                    metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
        else:
            pass

        apcorr = fastspec['APCORR']
        redshift = fastspec['CONTINUUM_Z']

        # rebuild the best-fitting broadband photometric model
        continuum_phot, synthmodelphot = self.SSP2data(
            self.sspflux, self.sspwave, redshift=redshift,
            synthphot=True, #AV=fastspec['CONTINUUM_AV'],
            coeff=fastspec['CONTINUUM_COEFF'] * self.massnorm)

        continuum_wave_phot = self.sspwave * (1 + redshift)
    
        wavemin, wavemax = 0.1, 35 # 6.0
        indx_phot = np.where((continuum_wave_phot/1e4 > wavemin) * (continuum_wave_phot/1e4 < wavemax))[0]     
    
        phot = self.parse_photometry(self.bands,
                                     maggies=np.array([metadata['FLUX_{}'.format(band.upper())] for band in self.bands]),
                                     ivarmaggies=np.array([metadata['FLUX_IVAR_{}'.format(band.upper())] for band in self.bands]),
                                     lambda_eff=allfilters.effective_wavelengths.value,
                                     min_uncertainty=self.min_uncertainty)
        fiberphot = self.parse_photometry(self.fiber_bands,
                                          maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in self.fiber_bands]),
                                          lambda_eff=filters.effective_wavelengths.value)

        # rebuild the best-fitting spectroscopic model
        stackwave = np.hstack(data['wave'])

        continuum, _ = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift, 
                                     specwave=data['wave'], specres=data['res'],
                                     cameras=data['cameras'],
                                     vdisp=fastspec['VDISP'],
                                     coeff=fastspec['CONTINUUM_COEFF'],
                                     synthphot=False)
        
        residuals = [apcorr*data['flux'][icam] - continuum[icam] for icam in np.arange(len(data['cameras']))]
        
        if np.all(fastspec['CONTINUUM_COEFF'] == 0):
            _smooth_continuum = np.zeros_like(stackwave)
        else:
            _smooth_continuum, _ = self.smooth_continuum(np.hstack(data['wave']), np.hstack(residuals),
                                                         np.hstack(data['ivar']), redshift=redshift,
                                                         linemask=np.hstack(data['linemask']))
        
        smooth_continuum = []
        for campix in data['camerapix']:
            smooth_continuum.append(_smooth_continuum[campix[0]:campix[1]])

        _emlinemodel = self.emlinemodel_bestfit(data['wave'], data['res'], fastspec)

        # individual-line spectra
        _emlinemodel_oneline = []
        for refline in self.linetable: # [self.inrange]: # for all lines in range
            T = Table(fastspec['CONTINUUM_Z', 'MGII_DOUBLET_RATIO', 'OII_DOUBLET_RATIO', 'SII_DOUBLET_RATIO'])
            for oneline in self.linetable: # need all lines for the model
                linename = oneline['name']
                for linecol in ['AMP', 'VSHIFT', 'SIGMA']:
                    col = linename.upper()+'_'+linecol
                    if linename == refline['name']:
                        T.add_column(Column(name=col, data=fastspec[col], dtype=fastspec[col].dtype))
                    else:
                        T.add_column(Column(name=col, data=0.0, dtype=fastspec[col].dtype))
            # special case the parameter doublets
            if refline['name'] == 'mgii_2796':
                T['MGII_2803_AMP'] = fastspec['MGII_2803_AMP']
            if refline['name'] == 'oii_3726':
                T['OII_3729_AMP'] = fastspec['OII_3729_AMP']
            if refline['name'] == 'sii_6731':
                T['SII_6716_AMP'] = fastspec['SII_6716_AMP']
            _emlinemodel_oneline1 = self.emlinemodel_bestfit(data['wave'], data['res'], T[0])
            if np.sum(np.hstack(_emlinemodel_oneline1)) > 0:
                _emlinemodel_oneline.append(_emlinemodel_oneline1)

        # QA choices

        # 8 columns: 3 for the SED, 5 for the spectra, and 8 for the lines
        # 8 rows: 4 for the SED, 2 each for the spectra, 1 gap, and 3 for the lines
        ngaprows = 1
        nlinerows = 3
        nlinecols = 8
        nrows = 4 + ngaprows + nlinerows
        ncols = 8
        nlinecols = 6

        fullheight = 18 # inches
        fullwidth = 24

        height_ratios = np.hstack(([1.0]*4, 0.25, [1.0]*nlinerows)) # small gap
    
        fig = plt.figure(figsize=(fullwidth, fullheight))
        gs = fig.add_gridspec(nrows, ncols, height_ratios=height_ratios)#, width_ratios=width_ratios)

        sedax = fig.add_subplot(gs[1:4, 5:]) # rows x cols
        specax1 = fig.add_subplot(gs[:2, :5])
        specax2 = fig.add_subplot(gs[2:4, :5])
        
        # full spectrum + best-fitting continuum model
        leg = {
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            'dv_narrow': '$\\Delta v_{{\\rm narrow}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['NARROW_Z']-redshift)),
            'dv_broad': '$\\Delta v_{{\\rm broad}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['BROAD_Z']-redshift)),
            'dv_uv': '$\\Delta v_{{\\rm UV}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['UV_Z']-redshift)),
            'sigma_narrow': '$\\sigma_{{\\rm narrow}}$={:.1f} km/s'.format(fastspec['NARROW_SIGMA']),
            'sigma_broad': '$\\sigma_{{\\rm broad}}$={:.1f} km/s'.format(fastspec['BROAD_SIGMA']),
            'sigma_uv': '$\\sigma_{{\\rm UV}}$={:.1f} km/s'.format(fastspec['UV_SIGMA']),
            #'targetid': '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),
            'cchi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastspec['CONTINUUM_RCHI2']),
            'rchi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastspec['RCHI2']),
            'deltarchi2': '$\\Delta\\chi^{{2}}_{{\\nu,\\rm broad,narrow}}$={:.3f}'.format(fastspec['DELTA_LINERCHI2']),
            #'zfastfastspec': '$z_{{\\rm fastspecfit}}$={:.6f}'.format(fastspec['CONTINUUM_Z']),
            #'z': '$z$={:.6f}'.format(fastspec['CONTINUUM_Z']),
            'age': '<Age>$={:.3f}$ Gyr'.format(fastspec['AGE']),
            'AV': '$A_{{V}}={:.3f}$ mag'.format(fastspec['AV']),
            'mstar': '$\\log_{{10}}(M_{{*}}/M_{{\odot}})={:.3f}$'.format(fastspec['LOGMSTAR']),
            'sfr': '${{\\rm SFR}}={:.2f}\ M_{{\odot}}/{{\\rm yr}}$'.format(fastspec['SFR']),
            'fagn': '$f_{{\\rm AGN}}={:.3f}$'.format(fastspec['FAGN']),
            'zzsun': '$Z/Z_{{\\odot}}={:.3f}$'.format(fastspec['ZZSUN']),
            }

        if fastspec['VDISP_IVAR'] == 0:
            leg.update({'vdisp': '$\\sigma_{{\\rm star}}={:.1f}$ km/s'.format(fastspec['VDISP'])})
        else:
            leg.update({'vdisp': '$\\sigma_{{\\rm star}}={:.1f}\\pm{:.1f}$ km/s'.format(
                fastspec['VDISP'], 1/np.sqrt(fastspec['VDISP_IVAR']))})
            
        ymin, ymax = 1e6, -1e6

        legxpos, legypos, legypos2, legfntsz = 0.98, 0.94, 0.05, 14 # 20
        bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.25)

        for ii in np.arange(len(data['cameras'])): # iterate over cameras
            sigma, good = ivar2var(data['ivar'][ii]/apcorr**2, sigma=True, allmasked_ok=True)

            specax1.plot(data['wave'][ii], apcorr*data['flux'][ii], color=col1[ii])
            #specax1.fill_between(data['wave'][ii], apcorr*data['flux'][ii]-sigma,
            #                    apcorr*data['flux'][ii]+sigma, color=col1[ii])
            specax1.plot(data['wave'][ii], continuum[ii]+smooth_continuum[ii], color=col2[ii])
            
            # get the robust range
            filtflux = median_filter(apcorr*data['flux'][ii], 31, mode='nearest')
            #filtflux = median_filter(data['flux'][ii] - _emlinemodel[ii], 51, mode='nearest')
            #perc = np.percentile(filtflux[data['ivar'][ii] > 0], [5, 95])
            #sigflux = np.std(apcorr*data['flux'][ii][data['ivar'][ii] > 0])
            I = data['ivar'][ii] > 0
            if np.sum(I) > 0:
                sigflux = np.diff(np.percentile(apcorr*data['flux'][ii][I], [25, 75]))[0] / 1.349 # robust
            else:
                sigflux = 0.0
            #sigflux = np.std(filtflux[data['ivar'][ii] > 0])
            #if -2 * perc[0] < ymin:
            #    ymin = -2 * perc[0]
            #ymin[ii] = np.min(data['flux'][ii])
            #ymax[ii] = np.max(data['flux'][ii])
            #if ymin < _ymin:
            #    ymin = _ymin
            #if ymax > _ymax:
            #    ymax = _ymax
            if -1.5 * sigflux < ymin:
                ymin = -1.5 * sigflux
            #if perc[1] > ymax:
            #    ymax = perc[1]
            #if np.min(filtflux) < ymin:
            #    ymin = np.min(filtflux) * 0.5
            if sigflux * 8 > ymax:
                ymax = sigflux * 10 # 5
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux) * 1.5
            #print(ymin, ymax)

        specax1.plot(stackwave, _smooth_continuum, color='gray')#col3[ii])#, alpha=0.3, lw=2)#, color='k')
        specax1.plot(stackwave, np.hstack(continuum), color='k', alpha=0.1)#col3[ii])#, alpha=0.3, lw=2)#, color='k')

        specax1.text(0.03, 0.9, 'Observed Spectrum + Continuum Model',
                     ha='left', va='center', transform=specax1.transAxes, fontsize=18)
        if not self.nolegend:
            txt = '\n'.join((
                r'{}'.format(leg['zredrock']),
                r'{}'.format(leg['vdisp']),
                ))
            specax1.text(legxpos, legypos, txt, ha='right', va='top',
                        transform=specax1.transAxes, fontsize=legfntsz,
                        bbox=bbox)
            #specax1.set_title(title, fontsize=22)
        
        specax1.set_xlim(spec_wavelims)
        specax1.set_ylim(ymin, ymax)
        specax1.set_xticklabels([])
        #specax1.set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
        #specax1.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 

        # full emission-line spectrum + best-fitting lines
        ymin, ymax = 1e6, -1e6
        allfullspec, allfullwave = [], []        
        for ii in np.arange(len(data['cameras'])): # iterate over cameras
            emlinewave = data['wave'][ii]
            emlineflux = apcorr*data['flux'][ii] - continuum[ii] - smooth_continuum[ii]
            emlinemodel = _emlinemodel[ii]
        
            emlinesigma, good = ivar2var(data['ivar'][ii]/apcorr**2, sigma=True, allmasked_ok=True)
            emlinewave = emlinewave[good]
            emlineflux = emlineflux[good]
            emlinesigma = emlinesigma[good]
            emlinemodel = emlinemodel[good]
        
            specax2.plot(emlinewave, emlineflux, color=col1[ii], alpha=0.7)
            #specax2.fill_between(emlinewave, emlineflux-emlinesigma,
            #                     emlineflux+emlinesigma, color=col1[ii], alpha=0.7)
            specax2.plot(emlinewave, emlinemodel, color=col2[ii], lw=3)

            allfullspec.append(continuum[ii] + _emlinemodel[ii])
            allfullwave.append(data['wave'][ii])
            
            # get the robust range
            filtflux = median_filter(emlineflux, 51, mode='nearest')
            #sigflux = np.std(filtflux)
            #sigflux = np.std(emlineflux)
            if np.sum(good) > 0:
                sigflux = np.diff(np.percentile(emlineflux, [25, 75]))[0] / 1.349 # robust
                if -2 * sigflux < ymin:
                    ymin = -2 * sigflux
                #if np.min(filtflux) < ymin:
                #    ymin = np.min(filtflux)
                #if np.min(emlinemodel) < ymin:
                #    ymin = 0.8 * np.min(emlinemodel)
                if 5 * sigflux > ymax:
                    ymax = 5 * sigflux
                if np.max(filtflux) > ymax:
                    ymax = np.max(filtflux)
                if np.max(emlinemodel) > ymax:
                    ymax = np.max(emlinemodel) * 1.2
                #print(ymin, ymax)
        
        allfullwave = np.hstack(allfullwave)
        allfullspec = np.hstack(allfullspec)
    
        if not self.nolegend:
            txt = '\n'.join((
                r'{} {}'.format(leg['rchi2'], leg['deltarchi2']),
                r'{} {}'.format(leg['dv_narrow'], leg['sigma_narrow']),
                r'{} {}'.format(leg['dv_broad'], leg['sigma_broad']),
                r'{} {}'.format(leg['dv_uv'], leg['sigma_uv']),
                ))
            specax2.text(legxpos, legypos, txt, ha='right', va='top',
                        transform=specax2.transAxes, fontsize=legfntsz,
                        bbox=bbox)
            specax2.text(0.03, 0.9, 'Residual Spectrum + Emission-Line Model',
                         ha='left', va='center', transform=specax2.transAxes,
                         fontsize=18)
                
        specax2.set_xlim(spec_wavelims)
        specax2.set_ylim(ymin, ymax)
        specax2.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 

        # photometric SED   
        if np.all(continuum_phot[indx_phot] <= 0):
            self.log.warning('Best-fitting photometric continuum is all zeros or negative!')
            continuum_phot_abmag = continuum_phot*0 + np.median(fiberphot['abmag'])
        else:
            indx_phot = indx_phot[continuum_phot[indx_phot] > 0] # trim zeros
            factor = 10**(0.4 * 48.6) * continuum_wave_phot[indx_phot]**2 / (C_LIGHT * 1e13) / self.fluxnorm / self.massnorm # [erg/s/cm2/A --> maggies]
            continuum_phot_abmag = -2.5*np.log10(continuum_phot[indx_phot] * factor)
            sedax.plot(continuum_wave_phot[indx_phot] / 1e4, continuum_phot_abmag, color='tan', zorder=1)
    
        sedax.scatter(synthmodelphot['lambda_eff']/1e4, synthmodelphot['abmag'], 
                   marker='s', s=200, color='k', facecolor='none',
                   #label=r'$grz$ (spectrum, synthesized)',
                   alpha=0.8, zorder=2)

        factor = 10**(0.4 * 48.6) * allfullwave**2 / (C_LIGHT * 1e13) / self.fluxnorm # [erg/s/cm2/A --> maggies]
        good = allfullspec > 0
        sedax.plot(allfullwave[good]/1e4, -2.5*np.log10(allfullspec[good]*factor[good]), color='gray', alpha=0.5)
        
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
        sedax.set_xlim(wavemin, wavemax)
        sedax.set_xscale('log')

        #sedax.set_ylabel('AB mag') 
        #sedax.set_ylabel(r'Apparent Brightness (AB mag)') 
        sedax.set_ylim(ymin, ymax)
        sedax.set_yticklabels([])
    
        sedax_twin = sedax.twinx()
        sedax_twin.set_ylabel('AB mag') 
        sedax_twin.set_ylim(ymin, ymax)
    
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
    
        dofit = np.where(self.bands_to_fit)[0]
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
    
        ignorefit = np.where(self.bands_to_fit == False)[0]
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
                   'DESI x {:.2f}'.format(apcorr), ha='center', va='center', fontsize=10,
                   color='gray')

        if not self.nolegend:
            txt = '\n'.join((
                r'{}'.format(leg['cchi2']),
                r'{}'.format(leg['fagn']),
                r'{}'.format(leg['zzsun']),
                r'{}'.format(leg['AV']),
                r'{}'.format(leg['sfr']),
                r'{}'.format(leg['age']),
                r'{}'.format(leg['mstar']),
                ))
            sedax.text(legxpos, legypos2, txt, ha='right', va='bottom',
                        transform=sedax.transAxes, fontsize=legfntsz,
                        bbox=bbox)

        # zoom in on individual emission lines - use linetable!
        linetable = self.linetable
        inrange = (linetable['restwave'] * (1+redshift) > spec_wavelims[0]) * (linetable['restwave'] * (1+redshift) < spec_wavelims[1])
        linetable = linetable[inrange]

        nline = len(set(linetable['plotgroup']))

        plotsig_default = 200.0 # [km/s]
        plotsig_default_balmer = 500.0 # [km/s]
        plotsig_default_broad = 2000.0 # [km/s]
        
        meanwaves, deltawaves, sigmas, linenames = [], [], [], []
        for plotgroup in set(linetable['plotgroup']):
            I = np.where(plotgroup == linetable['plotgroup'])[0]
            linenames.append(linetable['nicename'][I[0]].replace('-', ' '))
            meanwaves.append(np.mean(linetable['restwave'][I]))
            deltawaves.append((np.max(linetable['restwave'][I]) - np.min(linetable['restwave'][I])) / 2)
        
            sigmas1 = np.array([fastspec['{}_SIGMA'.format(line.upper())] for line in linetable[I]['name']])
            sigmas1 = sigmas1[sigmas1 > 0]
            if len(sigmas1) > 0:
                plotsig = 1.5*np.mean(sigmas1)
            else:
                if np.any(linetable['isbroad'][I]):
                    if np.any(linetable['isbalmer'][I]):
                        plotsig = fastspec['BROAD_SIGMA']
                        if plotsig < 50:
                            plotsig = fastspec['NARROW_SIGMA']
                            if plotsig < 50:
                                plotsig = plotsig_default
                                #plotsig = plotsig_default_broad
                    else:
                        plotsig = fastspec['UV_SIGMA']                    
                        if plotsig < 50:
                            plotsig = plotsig_default_broad
                else:
                    plotsig = fastspec['NARROW_SIGMA']
                    if plotsig < 50:
                        plotsig = plotsig_default
            sigmas.append(plotsig)
        
        srt = np.argsort(meanwaves)
        meanwaves = np.hstack(meanwaves)[srt]
        deltawaves = np.hstack(deltawaves)[srt]
        sigmas = np.hstack(sigmas)[srt]
        linenames = np.hstack(linenames)[srt]
        
        removelabels = np.ones(nline, bool)
        ymin, ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
        
        ax, irow, icol = [], 5, 0 # skip the gap row
        for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):
            icol = iax % nlinecols
            if iax > 0 and iax % nlinecols == 0:
                irow += 1
            #print(iax, irow, icol)
        
            xx = fig.add_subplot(gs[irow, icol])
            ax.append(xx)
        
            wmin = (meanwave - deltawave) * (1+redshift) - 6 * sig * meanwave * (1+redshift) / C_LIGHT
            wmax = (meanwave + deltawave) * (1+redshift) + 6 * sig * meanwave * (1+redshift) / C_LIGHT
            #print(linename, wmin, wmax)
        
            # iterate over cameras
            for ii in np.arange(len(data['cameras'])): # iterate over cameras
                emlinewave = data['wave'][ii]
                emlineflux = apcorr*data['flux'][ii] - continuum[ii] - smooth_continuum[ii]
                emlinemodel = _emlinemodel[ii]
        
                emlinesigma, good = ivar2var(data['ivar'][ii]/apcorr**2, sigma=True, allmasked_ok=True)
                emlinewave = emlinewave[good]
                emlineflux = emlineflux[good]
                emlinesigma = emlinesigma[good]
                emlinemodel = emlinemodel[good]
        
                #if ii == 0:
                #    import matplotlib.pyplot as plt ; plt.clf() ; plt.plot(emlinewave, emlineflux) ; plt.plot(emlinewave, emlinemodel) ; plt.xlim(4180, 4210) ; plt.ylim(-15, 17) ; plt.savefig('desi-users/ioannis/tmp/junkg.png')
                    
                emlinemodel_oneline = []
                for _emlinemodel_oneline1 in _emlinemodel_oneline:
                    emlinemodel_oneline.append(_emlinemodel_oneline1[ii][good])
        
                indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
                if len(indx) > 1:
                    removelabels[iax] = False
                    xx.plot(emlinewave[indx], emlineflux[indx], color=col1[ii], alpha=0.5)
                    #xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                    #                emlineflux[indx]+emlinesigma[indx], color=col1[ii], alpha=0.5)
                    # plot the individual lines first then the total model
                    for emlinemodel_oneline1 in emlinemodel_oneline:
                        if np.sum(emlinemodel_oneline1[indx]) > 0:
                            #P = emlinemodel_oneline1[indx] > 0
                            xx.plot(emlinewave[indx], emlinemodel_oneline1[indx], lw=1, alpha=0.8, color=col2[ii])
                    xx.plot(emlinewave[indx], emlinemodel[indx], color=col2[ii], lw=3)
                        
                    #xx.plot(emlinewave[indx], emlineflux[indx]-emlinemodel[indx], color='gray', alpha=0.3)
                    #xx.axhline(y=0, color='gray', ls='--')
        
                    # get the robust range
                    sigflux = np.std(emlineflux[indx])
                    filtflux = median_filter(emlineflux[indx], 3, mode='nearest')
        
                    _ymin, _ymax = -1.5 * sigflux, 4 * sigflux
                    if np.max(emlinemodel[indx]) > _ymax:
                        _ymax = np.max(emlinemodel[indx]) * 1.2
                    if np.max(filtflux) > _ymax:
                        _ymax = np.max(filtflux)
                    if np.min(emlinemodel[indx]) < _ymin:
                        _ymin = 0.8 * np.min(emlinemodel[indx])
                        
                    if _ymax > ymax[iax]:
                        ymax[iax] = _ymax
                    if _ymin < ymin[iax]:
                        ymin[iax] = _ymin
        
                    xx.set_xlim(wmin, wmax)
                    
                xx.text(0.03, 0.89, linename, ha='left', va='center',
                        transform=xx.transAxes, fontsize=12)
                
        for iax, xx in enumerate(ax):
            if removelabels[iax]:
                xx.set_ylim(0, 1)
                xx.set_xticklabels([])
                xx.set_yticklabels([])
            else:
                xx.set_ylim(ymin[iax], ymax[iax])
                xlim = xx.get_xlim()
                xx.xaxis.set_major_locator(ticker.MaxNLocator(2))
        
        # common axis labels
        tp, bt, lf, rt = 0.95, 0.08, 0.07, 0.94
        
        fig.text(lf-0.04, (tp-bt)/2+bt,
                 r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                 ha='center', va='center', rotation='vertical', fontsize=24)
        #fig.text(lf-0.07, (tp-bt)/4+bt,
        #         r'Flux Density ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
        #         ha='center', va='center', rotation='vertical', fontsize=30)
        fig.text((rt-lf)/2+lf, bt-0.05, r'Observed-frame Wavelength ($\AA$)',
                 ha='center', va='center', fontsize=24)

        plt.subplots_adjust(wspace=0.3, top=tp, bottom=bt, left=lf, right=rt, hspace=0.3)
        #plt.subplots_adjust(bottom=0.22, top=0.94, left=0.08, right=0.91)
        #plt.subplots_adjust(bottom=0.1, top=0.95, left=0.08, right=0.9)

        fig.suptitle(title, fontsize=22)

        self.log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()

