#!/usr/bin/env python
"""
fastspecfit.fastspecfit
=======================

See sandbox/running-fastspecfit for examples.

"""
import pdb # for debugging

import os, time
import numpy as np

def _fastspec_one(args):
    """Multiprocessing wrapper."""
    return fastspec_one(*args)

def _desiqa_one(args):
    """Multiprocessing wrapper."""
    return desiqa_one(*args)

def _assign_units_to_columns(fastfit, metadata, Spec, templates, fastphot, stackfit, log):
    """Assign astropy units to output tables.

    """
    from fastspecfit.io import init_fastspec_output
    
    fastcols = fastfit.colnames
    metacols = metadata.colnames

    T, M = init_fastspec_output(Spec.meta, Spec.specprod, templates=templates,
                                log=log, fastphot=fastphot, stackfit=stackfit)
    
    for col in T.colnames:
        if col in fastcols:
            fastfit[col].unit = T[col].unit
    for col in M.colnames:
        if col in metacols:
            metadata[col].unit = M[col].unit

def fastspec_one(iobj, data, out, meta, templates, log=None,
                 minspecwave=3500., maxspecwave=9900., broadlinefit=True,
                 fastphot=False, stackfit=False, nophoto=False, percamera_models=False):
    """Multiprocessing wrapper to run :func:`fastspec` on a single object.

    """
    from fastspecfit.io import cache_templates
    from fastspecfit.emlines import emline_specfit
    from fastspecfit.continuum import continuum_specfit
    
    log.info('Continuum- and emission-line fitting object {} [targetid {}, z={:.6f}].'.format(
        iobj, data['targetid'], meta['Z']))

    # Read the templates and then fit the continuum spectrum. Note that 450 A as
    # the minimum wavelength will allow us to synthesize u-band photometry only
    # up to z=5.53, even though some targets are at higher redshift. Handle this
    # case in continuum.ContinuumTools.
    templatecache = cache_templates(templates, log=log, mintemplatewave=450.0,
                                    maxtemplatewave=40e4, fastphot=fastphot)

    continuummodel, smooth_continuum = continuum_specfit(data, out, templatecache,
                                                         fastphot=fastphot, nophoto=nophoto,
                                                         log=log)

    # Optionally fit the emission-line spectrum.
    if fastphot:
        emmodel = None
    else:
        emmodel = emline_specfit(data, templatecache, out, continuummodel, smooth_continuum,
                                 minspecwave=minspecwave, maxspecwave=maxspecwave,
                                 broadlinefit=broadlinefit, percamera_models=percamera_models,
                                 log=log)
        
    return out, meta, emmodel

def desiqa_one(data, fastfit, metadata, templates, coadd_type,
               minspecwave=3500., maxspecwave=9900., fastphot=False, 
               nophoto=False, stackfit=False, outdir=None, outprefix=None,
               log=None):
    """Multiprocessing wrapper to generate QA for a single object.

    """
    from fastspecfit.io import cache_templates    
    templatecache = cache_templates(templates, log=log, mintemplatewave=450.0,
                                    maxtemplatewave=40e4, fastphot=fastphot)

    qa_fastspec(data, templatecache, fastfit, metadata, coadd_type=coadd_type,
                spec_wavelims=(minspecwave, maxspecwave),
                fastphot=fastphot, nophoto=nophoto, stackfit=stackfit,
                outprefix=outprefix, outdir=outdir, log=log)

def parse(options=None, log=None):
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
    parser.add_argument('--no-broadlinefit', default=True, action='store_false', dest='broadlinefit',
                        help='Do not allow for broad Balmer and Helium line-fitting.')
    parser.add_argument('--nophoto', action='store_true', help='Do not include the photometry in the model fitting.')
    parser.add_argument('--percamera-models', action='store_true', help='Return the per-camera (not coadded) model spectra.')
    parser.add_argument('--imf', type=str, default='chabrier', help='Initial mass function.')
    parser.add_argument('--templateversion', type=str, default='1.0.0', help='Template version number.')
    parser.add_argument('--templates', type=str, default=None, help='Optional full path and filename to the templates.')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--specfile-prefix', type=str, default='coadd-', help='Prefix of the spectral file(s).')
    parser.add_argument('--qnfile-prefix', type=str, default='qso_qn-', help='Prefix of the QuasarNet afterburner file(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--dr9dir', type=str, default=None, help='Optional directory name for the DR9 photometry.')
    parser.add_argument('--specproddir', type=str, default=None, help='Optional directory name for the spectroscopic production.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if log is None:
        from desiutil.log import get_logger
        log = get_logger()
            
    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspec {}'.format(' '.join(options)))

    return args

def fastspec(fastphot=False, stackfit=False, args=None, comm=None, verbose=False):
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
    from desiutil.log import get_logger, DEBUG
    from fastspecfit.io import DESISpectra, write_fastspecfit, init_fastspec_output

    if verbose:
        log = get_logger(DEBUG)
    else:
        log = get_logger()

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args, log=log)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Read the data.
    t0 = time.time()
    Spec = DESISpectra(dr9dir=args.dr9dir, mapdir=args.mapdir)

    if stackfit:
        data = Spec.read_stacked(args.redrockfiles, firsttarget=args.firsttarget,
                                 stackids=targetids, ntargets=args.ntargets,
                                 mp=args.mp)
        args.nophoto = True # force nophoto=True

        # stacked spectra can have variable beginning and ending wavelengths
        minspecwave = np.min(data[0]['coadd_wave']) - 20
        maxspecwave = np.max(data[0]['coadd_wave']) + 20
    else:
        Spec.select(args.redrockfiles, firsttarget=args.firsttarget, targetids=targetids,
                    ntargets=args.ntargets, redrockfile_prefix=args.redrockfile_prefix,
                    specfile_prefix=args.specfile_prefix, qnfile_prefix=args.qnfile_prefix,
                    specprod_dir=args.specproddir)
        if len(Spec.specfiles) == 0:
            return
    
        data = Spec.read_and_unpack(fastphot=fastphot, mp=args.mp)

        # default wavelength ranges
        minspecwave = 3500.
        maxspecwave = 9900.
        
    log.info('Reading and unpacking {} spectra to be fitted took {:.2f} seconds.'.format(
        Spec.ntargets, time.time()-t0))

    # Initialize the output tables.
    if args.templates is None:
        from fastspecfit.io import get_templates_filename
        templates = get_templates_filename(templateversion=args.templateversion, imf=args.imf)
    else:
        templates = args.templates

    out, meta = init_fastspec_output(Spec.meta, Spec.specprod, templates=templates,
                                     data=data, log=log, fastphot=fastphot, stackfit=stackfit)

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], templates, log,
                minspecwave, maxspecwave, args.broadlinefit, fastphot, stackfit,
                args.nophoto, args.percamera_models)
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
        from astropy.table import vstack
        try:
            # need to vstack to preserve the wavelength metadata 
            modelspectra = vstack(_out[2], metadata_conflicts='error')
        except:
            errmsg = 'Metadata conflict when stacking model spectra.'
            log.critical(errmsg)
            raise ValueError(errmsg)
       
    log.info('Fitting {} object(s) took {:.2f} seconds.'.format(Spec.ntargets, time.time()-t0))

    # Assign units and write out.
    _assign_units_to_columns(out, meta, Spec, templates, fastphot, stackfit, log)

    write_fastspecfit(out, meta, modelspectra=modelspectra, outfile=args.outfile,
                      specprod=Spec.specprod, coadd_type=Spec.coadd_type, fastphot=fastphot)

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

def stackfit(args=None, comm=None):
    """Wrapper script to fit (model) generic stacked spectra.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or ``None``
        Required and optional arguments parsed via inputs to the command line. 
    comm : :class:`mpi4py.MPI.MPI.COMM_WORLD` or `None`
        Intracommunicator used with MPI parallelism.

    """
    fastspec(stackfit=True, args=args, comm=comm)

def qa_fastspec(data, templatecache, fastspec, metadata, coadd_type='healpix',
                spec_wavelims=(3550, 9900), phot_wavelims=(0.1, 35),
                fastphot=False, nophoto=False, stackfit=False, outprefix=None,
                outdir=None, log=None):
    """QA plot the emission-line spectrum and best-fitting model.

    """
    import subprocess
    from scipy.ndimage import median_filter

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib import colors
    from matplotlib.patches import Circle, Rectangle, ConnectionPatch
    from matplotlib.lines import Line2D
    import matplotlib.gridspec as gridspec
    import matplotlib.image as mpimg

    import astropy.units as u
    from astropy.io import fits
    from astropy.wcs import WCS
    import seaborn as sns
    from PIL import Image, ImageDraw

    from fastspecfit.util import ivar2var, C_LIGHT
    from fastspecfit.io import FLUXNORM, get_qa_filename
    from fastspecfit.continuum import ContinuumTools
    from fastspecfit.emlines import build_emline_model, EMFitTools

    if log is None:
        from desiutil.log import get_logger
        log = get_logger()
    
    Image.MAX_IMAGE_PIXELS = None

    sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

    if stackfit:
        col1 = [colors.to_hex(col) for col in ['violet']]
        col2 = [colors.to_hex(col) for col in ['purple']]
    else:
        col1 = [colors.to_hex(col) for col in ['dodgerblue', 'darkseagreen', 'orangered']]
        col2 = [colors.to_hex(col) for col in ['darkblue', 'darkgreen', 'darkred']]

    photcol1 = colors.to_hex('darkorange')

    @ticker.FuncFormatter
    def major_formatter(x, pos):
        if (x >= 0.01) and (x < 0.1):
            return f'{x:.2f}'
        elif (x >= 0.1) and (x < 1):
            return f'{x:.1f}'
        else:
            return f'{x:.0f}'
    
    CTools = ContinuumTools(nophoto=nophoto,
                            continuum_pixkms=templatecache['continuum_pixkms'],
                            pixkms_wavesplit=templatecache['pixkms_wavesplit'])
    if not fastphot:
        EMFit = EMFitTools(minspecwave=spec_wavelims[0], maxspecwave=spec_wavelims[1])

    if metadata['PHOTSYS'] == 'S':
        filters = CTools.decam
        allfilters = CTools.decamwise
    else:
        filters = CTools.bassmzls
        allfilters = CTools.bassmzlswise

    redshift = fastspec['Z']

    pngfile = get_qa_filename(metadata, coadd_type, outprefix=outprefix,
                              outdir=outdir, fastphot=fastphot, log=log)

    # some arrays to use for the legend
    if coadd_type == 'healpix':
        target = [
            'Survey/Program/Healpix: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX']),
            'TargetID: {}'.format(metadata['TARGETID']),
            ]
    elif coadd_type == 'cumulative':
        target = [
            'Tile/Night/Fiber: {}/{}/{}'.format(metadata['TILEID'], metadata['NIGHT'], metadata['FIBER']),
            'TargetID: {}'.format(metadata['TARGETID']),
        ]
    elif coadd_type == 'pernight':
        target = [
            'Tile/Night/Fiber: {}/{}/{}'.format(metadata['TILEID'], metadata['NIGHT'], metadata['FIBER']),
            'TargetID: {}'.format(metadata['TARGETID']),
        ]
    elif coadd_type == 'perexp':
        target = [
            'Tile/Night/Fiber: {}/{}/{}'.format(metadata['TILEID'], metadata['NIGHT'], metadata['FIBER']),
            'Night/Fiber: {}/{}'.format(metadata['NIGHT'], metadata['FIBER']),
            'TargetID: {}'.format(metadata['TARGETID']),
        ]
    elif coadd_type == 'custom':
        target = [
            'Survey/Program/Healpix: {}/{}/{}'.format(metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX']),
            'TargetID: {}'.format(metadata['TARGETID']),
            ]
    elif coadd_type == 'stacked':
        target = [
            'StackID: {}'.format(metadata['STACKID']),
            '',
            ]
    else:
        errmsg = 'Unrecognized coadd_type {}!'.format(coadd_type)
        log.critical(errmsg)
        raise ValueError(errmsg)

    leg = {
        'z': '$z={:.7f}$'.format(redshift),

        'rchi2_phot': '$\\chi^{{2}}_{{\\nu, \\rm phot}}$={:.2f}'.format(fastspec['RCHI2_PHOT']),

        'dn4000_model': '$D_{{n}}(4000)_{{\\rm model}}={:.3f}$'.format(fastspec['DN4000_MODEL']),
        'age': 'Age$={:.3f}$ Gyr'.format(fastspec['AGE']),
        'AV': '$A_{{V}}={:.3f}$ mag'.format(fastspec['AV']),
        'mstar': '$\\log_{{10}}(M/M_{{\odot}})={:.3f}$'.format(fastspec['LOGMSTAR']),
        'sfr': '${{\\rm SFR}}={:.1f}\ M_{{\odot}}/{{\\rm yr}}$'.format(fastspec['SFR']),
        #'fagn': '$f_{{\\rm AGN}}={:.3f}$'.format(fastspec['FAGN']),
        'zzsun': '$Z/Z_{{\\odot}}={:.3f}$'.format(fastspec['ZZSUN']),

        'absmag_r': '$M_{{0.1r}}={:.2f}$'.format(fastspec['ABSMAG_SDSS_R']),
        'absmag_gr': '$^{{0.1}}(g-r)={:.3f}$'.format(fastspec['ABSMAG_SDSS_G']-fastspec['ABSMAG_SDSS_R']),
        'absmag_rz': '$^{{0.1}}(r-z)={:.3f}$'.format(fastspec['ABSMAG_SDSS_R']-fastspec['ABSMAG_SDSS_Z']),       
        }

    #leg['radec'] = '$(\\alpha,\\delta)=({:.7f},{:.6f})$'.format(metadata['RA'], metadata['DEC'])
    #leg['zwarn'] = '$z_{{\\rm warn}}={}$'.format(metadata['ZWARN'])

    if fastphot:
        leg['vdisp'] = '$\\sigma_{{star}}={:g}$ km/s'.format(fastspec['VDISP'])
    else:
        if fastspec['VDISP_IVAR'] > 0:
            leg['vdisp'] = '$\\sigma_{{star}}={:.0f}\pm{:.0f}$ km/s'.format(fastspec['VDISP'], 1/np.sqrt(fastspec['VDISP_IVAR']))
        else:
            leg['vdisp'] = '$\\sigma_{{star}}={:g}$ km/s'.format(fastspec['VDISP'])
            
        leg['rchi2'] = '$\\chi^{{2}}_{{\\nu, \\rm specphot}}$={:.2f}'.format(fastspec['RCHI2'])
        leg['rchi2_cont'] = '$\\chi^{{2}}_{{\\nu, \\rm cont}}$={:.2f}'.format(fastspec['RCHI2_CONT'])

    if not stackfit:
        if redshift != metadata['Z_RR']:
            leg['zredrock'] = '$z_{{\\rm Redrock}}={:.7f}$'.format(metadata['Z_RR'])

    if fastphot:
        fontsize1 = 16
        fontsize2 = 22
    else:
        if stackfit:
            fontsize1 = 16
            fontsize2 = 22
        else:
            fontsize1 = 18 # 24
            fontsize2 = 24
        
        apercorr = fastspec['APERCORR']
    
        if fastspec['DN4000_IVAR'] > 0:
            leg['dn4000_spec'] = '$D_{{n}}(4000)_{{\\rm data}}={:.3f}$'.format(fastspec['DN4000'])
            #leg.update({'dn4000_spec': '$D_{{n}}(4000)_{{\\rm spec}}={:.3f}\pm{:.3f}$'.format(fastspec['DN4000'], 1/np.sqrt(fastspec['DN4000_IVAR']))})
    
        # kinematics
        if fastspec['NARROW_Z'] != redshift:
            if fastspec['NARROW_ZRMS'] > 0:
                leg['dv_narrow'] = '$\\Delta v_{{\\rm narrow}}={:.0f}\pm{:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['NARROW_Z']-redshift), C_LIGHT*fastspec['NARROW_ZRMS'])
            else:
                leg['dv_narrow'] = '$\\Delta v_{{\\rm narrow}}={:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['NARROW_Z']-redshift))
        if fastspec['NARROW_SIGMA'] != 0.0:
            if fastspec['NARROW_SIGMARMS'] > 0:
                leg['sigma_narrow'] = '$\\sigma_{{\\rm narrow}}={:.0f}\pm{:.0f}$ km/s'.format(
                    fastspec['NARROW_SIGMA'], fastspec['NARROW_SIGMARMS'])
            else:
                leg['sigma_narrow'] = '$\\sigma_{{\\rm narrow}}={:.0f}$ km/s'.format(fastspec['NARROW_SIGMA'])
    
        snrcut = 1.5
        leg_broad, leg_narrow, leg_uv = {}, {}, {}
    
        if fastspec['UV_Z'] != redshift:
            if fastspec['UV_ZRMS'] > 0:
                leg_uv['dv_uv'] = '$\\Delta v_{{\\rm UV}}={:.0f}\pm{:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['UV_Z']-redshift), C_LIGHT*fastspec['UV_ZRMS'])
            else:
                leg_uv['dv_uv'] = '$\\Delta v_{{\\rm UV}}={:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['UV_Z']-redshift))
        if fastspec['UV_SIGMA'] != 0.0:
            if fastspec['UV_SIGMARMS'] > 0:
                leg_uv['sigma_uv'] = '$\\sigma_{{\\rm UV}}={:.0f}\pm{:.0f}$ km/s'.format(
                    fastspec['UV_SIGMA'], fastspec['UV_SIGMARMS'])
            else:
                leg_uv['sigma_uv'] = '$\\sigma_{{\\rm UV}}={:.0f}$ km/s'.format(fastspec['UV_SIGMA'])
        if fastspec['BROAD_Z'] != redshift:
            if fastspec['BROAD_ZRMS'] > 0:
                leg_broad['dv_broad'] = '$\\Delta v_{{\\rm broad}}={:.0f}\pm{:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['BROAD_Z']-redshift), C_LIGHT*fastspec['BROAD_ZRMS'])
            else:
                leg_broad['dv_broad'] = '$\\Delta v_{{\\rm broad}}={:.0f}$ km/s'.format(
                    C_LIGHT*(fastspec['BROAD_Z']-redshift))
        if fastspec['BROAD_SIGMA'] != 0.0:
            if fastspec['BROAD_SIGMARMS'] > 0:
                leg_broad['sigma_broad'] = '$\\sigma_{{\\rm broad}}={:.0f}\pm{:.0f}$ km/s'.format(
                    fastspec['BROAD_SIGMA'], fastspec['BROAD_SIGMARMS'])
            else:
                leg_broad['sigma_broad'] = '$\\sigma_{{\\rm broad}}={:.0f}$ km/s'.format(fastspec['BROAD_SIGMA'])
    
        # emission lines
    
        # UV
        if fastspec['LYALPHA_AMP']*np.sqrt(fastspec['LYALPHA_AMP_IVAR']) > snrcut:
            leg_uv['ewlya'] = 'EW(Ly$\\alpha)={:.1f}\ \\AA$'.format(fastspec['LYALPHA_EW'])
        if fastspec['CIV_1549_AMP']*np.sqrt(fastspec['CIV_1549_AMP_IVAR']) > snrcut:
            leg_uv['ewciv'] = 'EW(CIV)$={:.1f}\ \\AA$'.format(fastspec['CIV_1549_EW'])
        if fastspec['CIII_1908_AMP']*np.sqrt(fastspec['CIII_1908_AMP_IVAR']) > snrcut:
            leg_uv['ewciii'] = 'EW(CIII])$={:.1f}\ \\AA$'.format(fastspec['CIII_1908_EW'])
        if (fastspec['MGII_2796_AMP']*np.sqrt(fastspec['MGII_2796_AMP_IVAR']) > snrcut or
            fastspec['MGII_2803_AMP']*np.sqrt(fastspec['MGII_2803_AMP_IVAR']) > snrcut):
            leg_uv['ewmgii'] = 'EW(MgII)$={:.1f}\ \\AA$'.format(fastspec['MGII_2796_EW']+fastspec['MGII_2803_EW'])
            leg_uv['mgii_doublet'] = 'MgII $\lambda2796/\lambda2803={:.3f}$'.format(fastspec['MGII_DOUBLET_RATIO'])
    
        leg_broad['linerchi2'] = '$\\chi^{{2}}_{{\\nu,\\rm line}}$={:.2f}'.format(fastspec['RCHI2_LINE'])
        #leg_broad['deltarchi2'] = '$\\chi^{{2}}_{{\\nu,\\rm narrow}}-\\chi^{{2}}_{{\\nu,\\rm narrow+broad}}={:.3f}$'.format(fastspec['DELTA_LINERCHI2'])
        #if fastspec['DELTA_LINERCHI2'] != 0:
        leg_broad['deltarchi2'] = '$\\Delta\\chi^{{2}}_{{\\nu,\\rm nobroad}}={:.2f}$'.format(fastspec['DELTA_LINERCHI2'])
    
        # choose one broad Balmer line
        if fastspec['HALPHA_BROAD_AMP']*np.sqrt(fastspec['HALPHA_BROAD_AMP_IVAR']) > snrcut:
            leg_broad['ewbalmer_broad'] = 'EW(H$\\alpha)_{{\\rm broad}}={:.1f}\ \\AA$'.format(fastspec['HALPHA_BROAD_EW'])
        elif fastspec['HBETA_BROAD_AMP']*np.sqrt(fastspec['HBETA_BROAD_AMP_IVAR']) > snrcut:
            leg_broad['ewbalmer_broad'] = 'EW(H$\\beta)_{{\\rm broad}}={:.1f}\ \\AA$'.format(fastspec['HBETA_BROAD_EW'])
        elif fastspec['HGAMMA_BROAD_AMP']*np.sqrt(fastspec['HGAMMA_BROAD_AMP_IVAR']) > snrcut:
            leg_broad['ewbalmer_broad'] = 'EW(H$\\gamma)_{{\\rm broad}}={:.1f}\ \\AA$'.format(fastspec['HGAMMA_BROAD_EW'])
    
        if (fastspec['OII_3726_AMP']*np.sqrt(fastspec['OII_3726_AMP_IVAR']) > snrcut or 
            fastspec['OII_3729_AMP']*np.sqrt(fastspec['OII_3729_AMP_IVAR']) > snrcut):
            #leg_narrow['ewoii'] = 'EW([OII] $\lambda\lambda3726,29)={:.1f}\,\\AA$'.format(fastspec['OII_3726_EW']+fastspec['OII_3729_EW'])
            leg_narrow['ewoii'] = 'EW([OII])$={:.1f}\ \\AA$'.format(fastspec['OII_3726_EW']+fastspec['OII_3729_EW'])
    
        if fastspec['OIII_5007_AMP']*np.sqrt(fastspec['OIII_5007_AMP_IVAR']) > snrcut:
            leg_narrow['ewoiii'] = 'EW([OIII])$={:.1f}\,\\AA$'.format(fastspec['OIII_5007_EW'])
            #leg_narrow['ewoiii'] = 'EW([OIII] $\lambda5007={:.1f}\,\\AA$'.format(fastspec['OIII_5007_EW'])
    
        # choose one Balmer line
        if fastspec['HALPHA_AMP']*np.sqrt(fastspec['HALPHA_AMP_IVAR']) > snrcut:
            leg_narrow['ewbalmer_narrow'] = 'EW(H$\\alpha)={:.1f}\ \\AA$'.format(fastspec['HALPHA_EW'])
        elif fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut:
            leg_narrow['ewbalmer_narrow'] = 'EW(H$\\beta)={:.1f}\ \\AA$'.format(fastspec['HBETA_EW'])
        elif fastspec['HGAMMA_AMP']*np.sqrt(fastspec['HGAMMA_AMP_IVAR']) > snrcut:
            leg_narrow['ewbalmer_narrow'] = 'EW(H$\\gamma)={:.1f}\ \\AA$'.format(fastspec['HGAMMA_EW'])
    
        if (fastspec['HALPHA_AMP']*np.sqrt(fastspec['HALPHA_AMP_IVAR']) > snrcut and 
            fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut):
            leg_narrow['hahb'] = '${{\\rm H}}\\alpha/{{\\rm H}}\\beta={:.3f}$'.format(fastspec['HALPHA_FLUX']/fastspec['HBETA_FLUX'])
        if 'hahb' not in leg_narrow.keys() and (fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut and 
            fastspec['HGAMMA_AMP']*np.sqrt(fastspec['HGAMMA_AMP_IVAR']) > snrcut):
            leg_narrow['hbhg'] = '${{\\rm H}}\\beta/{{\\rm H}}\\gamma={:.3f}$'.format(fastspec['HBETA_FLUX']/fastspec['HGAMMA_FLUX'])
        if (fastspec['HBETA_AMP']*np.sqrt(fastspec['HBETA_AMP_IVAR']) > snrcut and 
            fastspec['OIII_5007_AMP']*np.sqrt(fastspec['OIII_5007_AMP_IVAR']) > snrcut and 
            fastspec['HBETA_FLUX'] > 0 and fastspec['OIII_5007_FLUX'] > 0):
            leg_narrow['oiiihb'] = '$\\log_{{10}}({{\\rm [OIII]/H}}\\beta)={:.3f}$'.format(np.log10(fastspec['OIII_5007_FLUX']/fastspec['HBETA_FLUX']))
        if (fastspec['HALPHA_AMP']*np.sqrt(fastspec['HALPHA_AMP_IVAR']) > snrcut and 
            fastspec['NII_6584_AMP']*np.sqrt(fastspec['NII_6584_AMP_IVAR']) > snrcut and 
            fastspec['HALPHA_FLUX'] > 0 and fastspec['NII_6584_FLUX'] > 0):
            leg_narrow['niiha'] = '$\\log_{{10}}({{\\rm [NII]/H}}\\alpha)={:.3f}$'.format(np.log10(fastspec['NII_6584_FLUX']/fastspec['HALPHA_FLUX']))
    
        if (fastspec['OII_3726_AMP']*np.sqrt(fastspec['OII_3726_AMP_IVAR']) > snrcut or 
            fastspec['OII_3729_AMP']*np.sqrt(fastspec['OII_3729_AMP_IVAR']) > snrcut):
            #if fastspec['OII_DOUBLET_RATIO'] != 0:
            leg_narrow['oii_doublet'] = '[OII] $\lambda3726/\lambda3729={:.3f}$'.format(fastspec['OII_DOUBLET_RATIO'])
    
        if fastspec['SII_6716_AMP']*np.sqrt(fastspec['SII_6716_AMP_IVAR']) > snrcut or fastspec['SII_6731_AMP']*np.sqrt(fastspec['SII_6731_AMP_IVAR']) > snrcut:
            #if fastspec['SII_DOUBLET_RATIO'] != 0:
            leg_narrow['sii_doublet'] = '[SII] $\lambda6731/\lambda6716={:.3f}$'.format(fastspec['SII_DOUBLET_RATIO'])

    # rebuild the best-fitting broadband photometric model
    if not stackfit:
        sedmodel, sedphot = CTools.templates2data(
                templatecache['templateflux'], templatecache['templatewave'],
                redshift=redshift, dluminosity=data['dluminosity'],
                synthphot=True, coeff=fastspec['COEFF'] * CTools.massnorm)
        sedwave = templatecache['templatewave'] * (1 + redshift)
    
        phot = CTools.parse_photometry(CTools.bands,
                                       maggies=np.array([metadata['FLUX_{}'.format(band.upper())] for band in CTools.bands]),
                                       ivarmaggies=np.array([metadata['FLUX_IVAR_{}'.format(band.upper())] for band in CTools.bands]),
                                       lambda_eff=allfilters.effective_wavelengths.value,
                                       min_uncertainty=CTools.min_uncertainty)
        fiberphot = CTools.parse_photometry(CTools.fiber_bands,
                                            maggies=np.array([metadata['FIBERTOTFLUX_{}'.format(band.upper())] for band in CTools.fiber_bands]),
                                            lambda_eff=filters.effective_wavelengths.value)
    
        indx_phot = np.where((sedmodel > 0) * (sedwave/1e4 > phot_wavelims[0]) * 
                             (sedwave/1e4 < phot_wavelims[1]))[0]
        sedwave = sedwave[indx_phot]
        sedmodel = sedmodel[indx_phot]

    if not fastphot:
        # Rebuild the best-fitting spectroscopic model; prefix "desi" means
        # "per-camera" and prefix "full" has the cameras h-stacked.
        fullwave = np.hstack(data['wave'])
    
        desicontinuum, _ = CTools.templates2data(templatecache['templateflux_nolines'], templatecache['templatewave'],
                                                 redshift=redshift, dluminosity=data['dluminosity'], synthphot=False,
                                                 specwave=data['wave'], specres=data['res'],
                                                 specmask=data['mask'], cameras=data['cameras'],
                                                 vdisp=fastspec['VDISP'],
                                                 coeff=fastspec['COEFF'])
    
        # remove the aperture correction
        desicontinuum = [_desicontinuum / apercorr for _desicontinuum in desicontinuum]
        fullcontinuum = np.hstack(desicontinuum)
    
         # Need to be careful we don't pass a large negative residual where
         # there are gaps in the data.
        desiresiduals = []
        for icam in np.arange(len(data['cameras'])):
            resid = data['flux'][icam] - desicontinuum[icam]
            I = (data['flux'][icam] == 0.0) * (data['flux'][icam] == 0.0)
            if np.any(I):
                resid[I] = 0.0
            desiresiduals.append(resid)
        
        if np.all(fastspec['COEFF'] == 0):
            fullsmoothcontinuum = np.zeros_like(fullwave)
        else:
            fullsmoothcontinuum, _ = CTools.smooth_continuum(
                    fullwave, np.hstack(desiresiduals), np.hstack(data['ivar']), 
                    redshift=redshift, linemask=np.hstack(data['linemask']))
    
        desismoothcontinuum = []
        for campix in data['camerapix']:
            desismoothcontinuum.append(fullsmoothcontinuum[campix[0]:campix[1]])
    
        # full model spectrum + individual line-spectra
        desiemlines = EMFit.emlinemodel_bestfit(data['wave'], data['res'], fastspec)
    
        desiemlines_oneline = []
        inrange = ( (EMFit.linetable['restwave'] * (1+redshift) > np.min(fullwave)) *
                    (EMFit.linetable['restwave'] * (1+redshift) < np.max(fullwave)) )
         
        for oneline in EMFit.linetable[inrange]: # for all lines in range
            linename = oneline['name'].upper()
            amp = fastspec['{}_AMP'.format(linename)]
            if amp != 0:
                desiemlines_oneline1 = build_emline_model(
                    EMFit.log10wave, redshift, np.array([amp]),
                    np.array([fastspec['{}_VSHIFT'.format(linename)]]),
                    np.array([fastspec['{}_SIGMA'.format(linename)]]),
                    np.array([oneline['restwave']]), data['wave'], data['res'])
                desiemlines_oneline.append(desiemlines_oneline1)

    # Grab the viewer cutout.
    if not stackfit:
        pixscale = 0.262
        width = int(30 / pixscale)   # =1 arcmin
        height = int(width / 1.3) # 3:2 aspect ratio
    
        hdr = fits.Header()
        hdr['NAXIS'] = 2
        hdr['NAXIS1'] = width
        hdr['NAXIS2'] = height
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['CRVAL1'] = metadata['RA']
        hdr['CRVAL2'] = metadata['DEC']
        hdr['CRPIX1'] = width/2+0.5
        hdr['CRPIX2'] = height/2+0.5
        hdr['CD1_1'] = -pixscale/3600
        hdr['CD1_2'] = 0.0
        hdr['CD2_1'] = 0.0
        hdr['CD2_2'] = +pixscale/3600
        wcs = WCS(hdr)
    
        cutoutjpeg = os.path.join(outdir, 'tmp.'+os.path.basename(pngfile.replace('.png', '.jpeg')))
        if not os.path.isfile(cutoutjpeg):
            wait = 15 # seconds
            cmd = 'timeout {wait} wget -q -o /dev/null -O {outfile} https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&width={width}&height={height}&layer=ls-dr9'
            cmd = cmd.format(wait=wait, outfile=cutoutjpeg, ra=metadata['RA'],
                             dec=metadata['DEC'], width=width, height=height)
            log.info(cmd)
            try:
                subprocess.check_call(cmd.split(), stderr=subprocess.DEVNULL)
            except:
                log.warning('No cutout from viewer after {} seconds; stopping wget'.format(wait))
                
        try:
            img = mpimg.imread(cutoutjpeg)
        except:
            log.warning('Problem reading cutout for targetid {}.'.format(metadata['TARGETID']))
            img = np.zeros((height, width, 3))
    
        if os.path.isfile(cutoutjpeg):
            os.remove(cutoutjpeg)
        
    # QA choices
    legxpos, legypos, legypos2, legfntsz1, legfntsz = 0.98, 0.94, 0.05, 16, 18
    bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.15)
    bbox2 = dict(boxstyle='round', facecolor='lightgray', alpha=0.7)

    if fastphot:
        fullheight = 9 # inches
        fullwidth = 18
        
        nrows = 3
        ncols = 8

        fig = plt.figure(figsize=(fullwidth, fullheight))
        gs = fig.add_gridspec(nrows, ncols)#, width_ratios=width_ratios)

        cutax = fig.add_subplot(gs[0:2, 5:8], projection=wcs) # rows x cols
        sedax = fig.add_subplot(gs[0:3, 0:5])
    elif stackfit:
        fullheight = 14 # inches
        fullwidth = 24

        # 8 columns: 3 for the spectra, and 5 for the lines
        # 8 rows: 4 for the SED, 2 each for the spectra, 1 gap, and 3 for the lines
        nlinerows = 6
        nlinecols = 3
        nrows = nlinerows
        ncols = 8
    
        #height_ratios = np.hstack(([1.0]*3, 0.25, [1.0]*6))
        #width_ratios = np.hstack(([1.0]*5, [1.0]*3))
    
        fig = plt.figure(figsize=(fullwidth, fullheight))
        gs = fig.add_gridspec(nrows, ncols)#, height_ratios=height_ratios)#, width_ratios=width_ratios)
    
        specax = fig.add_subplot(gs[0:4, 0:5])
    else:
        fullheight = 18 # inches
        fullwidth = 24

        # 8 columns: 3 for the SED, 5 for the spectra, and 8 for the lines
        # 8 rows: 4 for the SED, 2 each for the spectra, 1 gap, and 3 for the lines
        ngaprows = 1
        nlinerows = 6
        nlinecols = 3
        nrows = 9 + ngaprows
        ncols = 8
    
        height_ratios = np.hstack(([1.0]*3, 0.25, [1.0]*6)) # small gap
        #width_ratios = np.hstack(([1.0]*5, [1.0]*3))
    
        fig = plt.figure(figsize=(fullwidth, fullheight))
        gs = fig.add_gridspec(nrows, ncols, height_ratios=height_ratios)#, width_ratios=width_ratios)
    
        cutax = fig.add_subplot(gs[0:3, 5:8], projection=wcs) # rows x cols
        sedax = fig.add_subplot(gs[0:3, 0:5])
        specax = fig.add_subplot(gs[4:8, 0:5])
    
    # viewer cutout
    if not stackfit:
        cutax.imshow(img, origin='lower')#, interpolation='nearest')
        cutax.set_xlabel('RA [J2000]')
        cutax.set_ylabel('Dec [J2000]')
    
        cutax.coords[1].set_ticks_position('r')
        cutax.coords[1].set_ticklabel_position('r')
        cutax.coords[1].set_axislabel_position('r')
    
        if metadata['DEC'] > 0:
            sgn = '+'
        else:
            sgn = ''
            
        cutax.text(0.04, 0.95, '$(\\alpha,\\delta)$=({:.7f}, {}{:.6f})'.format(metadata['RA'], sgn, metadata['DEC']),
                   ha='left', va='top', color='k', fontsize=fontsize1, bbox=bbox2,
                   transform=cutax.transAxes)
    
        sz = img.shape
        cutax.add_artist(Circle((sz[1] / 2, sz[0] / 2), radius=1.5/2/pixscale, facecolor='none', # DESI fiber=1.5 arcsec diameter
                                edgecolor='red', ls='-', alpha=0.8))#, label='3" diameter'))
        cutax.add_artist(Circle((sz[1] / 2, sz[0] / 2), radius=10/2/pixscale, facecolor='none',
                                edgecolor='red', ls='--', alpha=0.8))#, label='15" diameter'))
        handles = [Line2D([0], [0], color='red', lw=2, ls='-', label='1.5 arcsec'),
                   Line2D([0], [0], color='red', lw=2, ls='--', label='10 arcsec')]
    
        cutax.legend(handles=handles, loc='lower left', fontsize=fontsize1, facecolor='lightgray')

    if not fastphot:
        # plot the full spectrum + best-fitting (total) model
        specax.plot(fullwave/1e4, fullsmoothcontinuum, color='gray', alpha=0.4)
        specax.plot(fullwave/1e4, fullcontinuum, color='k', alpha=0.6)

        spec_ymin, spec_ymax = 1e6, -1e6
    
        desimodelspec = []
        for icam in np.arange(len(data['cameras'])): # iterate over cameras
            wave = data['wave'][icam]
            flux = data['flux'][icam]
            modelflux = desiemlines[icam] + desicontinuum[icam] + desismoothcontinuum[icam]
    
            sigma, camgood = ivar2var(data['ivar'][icam], sigma=True, allmasked_ok=True, clip=0)
    
            wave = wave[camgood]
            flux = flux[camgood]
            sigma = sigma[camgood]
            modelflux = modelflux[camgood]
    
            desimodelspec.append(apercorr * (desicontinuum[icam] + desiemlines[icam]))
            
            # get the robust range
            filtflux = median_filter(flux, 51, mode='nearest')
            if np.sum(camgood) > 0:
                sigflux = np.diff(np.percentile(flux - modelflux, [25, 75]))[0] / 1.349 # robust
                if -2 * sigflux < spec_ymin:
                    spec_ymin = -2 * sigflux
                if 6 * sigflux > spec_ymax:
                    spec_ymax = 6 * sigflux
                if np.max(filtflux) > spec_ymax:
                    #print(icam, spec_ymax, np.max(filtflux), np.max(filtflux) * 1.2)
                    spec_ymax = np.max(filtflux) * 1.25
                if np.max(modelflux) > spec_ymax:
                    spec_ymax = np.max(modelflux) * 1.25
                #print(spec_ymin, spec_ymax)
                #pdb.set_trace()
        
            #specax.fill_between(wave, flux-sigma, flux+sigma, color=col1[icam], alpha=0.2)
            specax.plot(wave/1e4, flux, color=col1[icam], alpha=0.8)
            specax.plot(wave/1e4, modelflux, color=col2[icam], lw=2, alpha=0.8)
    
        fullmodelspec = np.hstack(desimodelspec)

        if stackfit:
            specax_twin = specax.twiny()
            specax_twin.set_xlim(spec_wavelims[0]/(1+redshift)/1e4, spec_wavelims[1]/(1+redshift)/1e4)
            specax_twin.xaxis.set_major_formatter(major_formatter)
            restticks = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
            restticks = restticks[(restticks >= spec_wavelims[0]/(1+redshift)/1e4) * (restticks <= spec_wavelims[1]/(1+redshift)/1e4)]
            specax_twin.set_xticks(restticks)
        else:
            specax.spines[['top']].set_visible(False)
            
        specax.set_xlim(spec_wavelims[0]/1e4, spec_wavelims[1]/1e4)
        specax.set_ylim(spec_ymin, spec_ymax)
        specax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        #specax.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        specax.set_ylabel(r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$')

    # photometric SED
    if not stackfit:    
        abmag_good = phot['abmag_ivar'] > 0
        abmag_goodlim = phot['abmag_limit'] > 0
        
        if len(sedmodel) == 0:
            log.warning('Best-fitting photometric continuum is all zeros or negative!')
            if np.sum(abmag_good) > 0:
                medmag = np.median(phot['abmag'][abmag_good])
            elif np.sum(abmag_goodlim) > 0:
                medmag = np.median(phot['abmag_limit'][abmag_goodlim])
            else:
                medmag = 0.0
            sedmodel_abmag = np.zeros_like(templatecache['templatewave']) + medmag
        else:
            factor = 10**(0.4 * 48.6) * sedwave**2 / (C_LIGHT * 1e13) / FLUXNORM / CTools.massnorm # [erg/s/cm2/A --> maggies]
            sedmodel_abmag = -2.5*np.log10(sedmodel * factor)
            sedax.plot(sedwave / 1e4, sedmodel_abmag, color='slategrey', alpha=0.9, zorder=1)
    
        sedax.scatter(sedphot['lambda_eff']/1e4, sedphot['abmag'], marker='s', 
                      s=400, color='k', facecolor='none', linewidth=2, alpha=1.0, zorder=2)
    
        if not fastphot:
            #factor = 10**(0.4 * 48.6) * fullwave**2 / (C_LIGHT * 1e13) / FLUXNORM # [erg/s/cm2/A --> maggies]
            #good = fullmodelspec > 0
            #sedax.plot(fullwave[good]/1e4, -2.5*np.log10(fullmodelspec[good]*factor[good]), color='purple', alpha=0.8)
            for icam in np.arange(len(data['cameras'])):
                factor = 10**(0.4 * 48.6) * data['wave'][icam]**2 / (C_LIGHT * 1e13) / FLUXNORM # [erg/s/cm2/A --> maggies]
                good = desimodelspec[icam] > 0
                sedax.plot(data['wave'][icam][good]/1e4, -2.5*np.log10(desimodelspec[icam][good]*factor[good]), color=col2[icam], alpha=0.8)
    
        # we have to set the limits *before* we call errorbar, below!
        dm = 1.5
        sed_ymin = np.nanmax(sedmodel_abmag) + dm
        sed_ymax = np.nanmin(sedmodel_abmag) - dm
        if np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) > 0:
            sed_ymin = np.max((np.nanmax(phot['abmag'][abmag_good]), np.nanmax(phot['abmag_limit'][abmag_goodlim]), np.nanmax(sedmodel_abmag))) + dm
            sed_ymax = np.min((np.nanmin(phot['abmag'][abmag_good]), np.nanmin(phot['abmag_limit'][abmag_goodlim]), np.nanmin(sedmodel_abmag))) - dm
        elif np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) == 0:
            sed_ymin = np.max((np.nanmax(phot['abmag'][abmag_good]), np.nanmax(sedmodel_abmag))) + dm
            sed_ymax = np.min((np.nanmin(phot['abmag'][abmag_good]), np.nanmin(sedmodel_abmag))) - dm
        elif np.sum(abmag_good) == 0 and np.sum(abmag_goodlim) > 0:
            sed_ymin = np.max((np.nanmax(phot['abmag_limit'][abmag_goodlim]), np.nanmax(sedmodel_abmag))) + dm
            sed_ymax = np.min((np.nanmin(phot['abmag_limit'][abmag_goodlim]), np.nanmin(sedmodel_abmag))) - dm
        else:
            abmag_good = phot['abmag'] > 0
            abmag_goodlim = phot['abmag_limit'] > 0
            if np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) > 0:
                sed_ymin = np.max((np.nanmax(phot['abmag'][abmag_good]), np.nanmax(phot['abmag_limit'][abmag_goodlim]))) + dm
                sed_ymax = np.min((np.nanmin(phot['abmag'][abmag_good]), np.nanmin(phot['abmag_limit'][abmag_goodlim]))) - dm
            elif np.sum(abmag_good) > 0 and np.sum(abmag_goodlim) == 0:                
                sed_ymin = np.nanmax(phot['abmag'][abmag_good]) + dm
                sed_ymax = np.nanmin(phot['abmag'][abmag_good]) - dm
            elif np.sum(abmag_good) == 0 and np.sum(abmag_goodlim) > 0:
                sed_ymin = np.nanmax(phot['abmag_limit'][abmag_goodlim]) + dm
                sed_ymax = np.nanmin(phot['abmag_limit'][abmag_goodlim]) - dm
            #else:
            #    sed_ymin, sed_ymax = [30, 20]
            
        if sed_ymin > 30:
            sed_ymin = 30
        if np.isnan(sed_ymin) or np.isnan(sed_ymax):
            raise('Problem here!')
        #print(sed_ymin, sed_ymax)
    
        #sedax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        sedax.set_xlim(phot_wavelims[0], phot_wavelims[1])
        sedax.set_xscale('log')
        sedax.set_ylabel('AB mag') 
        #sedax.set_ylabel(r'Apparent Brightness (AB mag)') 
        sedax.set_ylim(sed_ymin, sed_ymax)
    
        sedax.xaxis.set_major_formatter(major_formatter)
        obsticks = np.array([0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 20.0])
        sedax.set_xticks(obsticks)
    
        if fastphot:
            sedax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
    
        sedax_twin = sedax.twiny()
        sedax_twin.set_xlim(phot_wavelims[0]/(1+redshift), phot_wavelims[1]/(1+redshift))
        sedax_twin.set_xscale('log')
        #sedax_twin.set_xlabel(r'Rest-frame Wavelength ($\mu$m)') 
    
        sedax_twin.xaxis.set_major_formatter(major_formatter)
        restticks = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 15.0, 20.0])
        restticks = restticks[(restticks >= phot_wavelims[0]/(1+redshift)) * (restticks <= phot_wavelims[1]/(1+redshift))]
        sedax_twin.set_xticks(restticks)
    
        # integrated flux / photometry
        abmag = np.squeeze(phot['abmag'])
        abmag_limit = np.squeeze(phot['abmag_limit'])
        abmag_fainterr = np.squeeze(phot['abmag_fainterr'])
        abmag_brighterr = np.squeeze(phot['abmag_brighterr'])
        yerr = np.squeeze([abmag_brighterr, abmag_fainterr])
    
        markersize = 14
    
        dofit = np.where(CTools.bands_to_fit)[0]
        if len(dofit) > 0:
            good = np.where((abmag[dofit] > 0) * (abmag_limit[dofit] == 0))[0]
            upper = np.where(abmag_limit[dofit] > 0)[0]
            if len(good) > 0:
                sedax.errorbar(phot['lambda_eff'][dofit][good]/1e4, abmag[dofit][good],
                               yerr=yerr[:, dofit[good]],
                               fmt='o', markersize=markersize, markeredgewidth=1, markeredgecolor='k',
                               markerfacecolor=photcol1, elinewidth=3, ecolor=photcol1, capsize=4,
                               label=r'$grz\,W_{1}W_{2}W_{3}W_{4}$', zorder=2, alpha=1.0)
            if len(upper) > 0:
                sedax.errorbar(phot['lambda_eff'][dofit][upper]/1e4, abmag_limit[dofit][upper],
                               lolims=True, yerr=0.75,
                               fmt='o', markersize=markersize, markeredgewidth=3, markeredgecolor='k',
                               markerfacecolor=photcol1, elinewidth=3, ecolor=photcol1, capsize=4, alpha=0.7)
    
        ignorefit = np.where(CTools.bands_to_fit == False)[0]
        if len(ignorefit) > 0:
            good = np.where((abmag[ignorefit] > 0) * (abmag_limit[ignorefit] == 0))[0]
            upper = np.where(abmag_limit[ignorefit] > 0)[0]
            if len(good) > 0:
                sedax.errorbar(phot['lambda_eff'][ignorefit][good]/1e4, abmag[ignorefit][good],
                               yerr=yerr[:, ignorefit[good]],
                               fmt='o', markersize=markersize, markeredgewidth=3, markeredgecolor='k',
                               markerfacecolor='none', elinewidth=3, ecolor=photcol1, capsize=4, alpha=0.7)
            if len(upper) > 0:
                sedax.errorbar(phot['lambda_eff'][ignorefit][upper]/1e4, abmag_limit[ignorefit][upper],
                               lolims=True, yerr=0.75, fmt='o', markersize=markersize, markeredgewidth=3,
                               markeredgecolor='k', markerfacecolor='none', elinewidth=3,
                               ecolor=photcol1, capsize=5, alpha=0.7)
    
        if fastphot:
            txt = leg['rchi2_phot']
        else:
            # Label the DESI wavelength range and the aperture correction.
            sedax.plot([np.min(fullwave)/1e4, np.max(fullwave)/1e4], [sed_ymin-1, sed_ymin-1],
                       lw=2, ls='-', color='gray', marker='s')#, alpha=0.5)
            sedax.text(((np.max(fullwave)-np.min(fullwave))/2+np.min(fullwave)*0.8)/1e4, sed_ymin-1.7,
                       'DESI x {:.2f}'.format(apercorr), ha='center', va='center', fontsize=16,
                       color='k')
    
            txt = '\n'.join((leg['rchi2_cont'], leg['rchi2_phot'], leg['rchi2']))
        
        sedax.text(0.02, 0.94, txt, ha='left', va='top',
                   transform=sedax.transAxes, fontsize=legfntsz)#, bbox=bbox)
    
        txt = '\n'.join((
            #r'{}'.format(leg['fagn']),
            r'{}'.format(leg['zzsun']),
            r'{}'.format(leg['AV']),
            r'{}'.format(leg['sfr']),
            r'{}'.format(leg['age']),
            r'{}'.format(leg['mstar']),
            ))
        sedax.text(legxpos, legypos2, txt, ha='right', va='bottom',
                    transform=sedax.transAxes, fontsize=legfntsz1, bbox=bbox)
    
        if not fastphot:
            # draw lines connecting the SED and spectral plots
            sedax.add_artist(ConnectionPatch(xyA=(spec_wavelims[0]/1e4, sed_ymin), 
                                             xyB=(spec_wavelims[0]/1e4, spec_ymax), 
                                             coordsA='data', coordsB='data',
                                             axesA=sedax, axesB=specax, color='k'))
            sedax.add_artist(ConnectionPatch(xyA=(spec_wavelims[1]/1e4, sed_ymin), 
                                             xyB=(spec_wavelims[1]/1e4, spec_ymax), 
                                             coordsA='data', coordsB='data',
                                             axesA=sedax, axesB=specax, color='k'))
    
    # zoom in on individual emission lines - use linetable!
    if not fastphot:
        linetable = EMFit.linetable
        inrange = (linetable['restwave'] * (1+redshift) > np.min(fullwave)) * (linetable['restwave'] * (1+redshift) < np.max(fullwave))
        linetable = linetable[inrange]
    
        nline = len(set(linetable['plotgroup']))
    
        plotsig_default = 200.0 # [km/s]
        plotsig_default_balmer = 500.0 # [km/s]
        plotsig_default_broad = 2000.0 # [km/s]
    
        meanwaves, deltawaves, sigmas, linenames = [], [], [], []
        for plotgroup in set(linetable['plotgroup']):
            I = np.where(plotgroup == linetable['plotgroup'])[0]
            linename = linetable['nicename'][I[0]].replace('-', ' ')
            linenames.append(linename)
            meanwaves.append(np.mean(linetable['restwave'][I]))
            deltawaves.append((np.max(linetable['restwave'][I]) - np.min(linetable['restwave'][I])) / 2)
        
            sigmas1 = np.array([fastspec['{}_SIGMA'.format(line.upper())] for line in linetable[I]['name']])
            sigmas1 = sigmas1[sigmas1 > 0]
            if len(sigmas1) > 0:
                plotsig = 1.5*np.mean(sigmas1)
                if plotsig < 50:
                    plotsig = 50.0
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
    
        if len(linetable) == 0:
            ax = []
        else:
            srt = np.argsort(meanwaves)
            meanwaves = np.hstack(meanwaves)[srt]
            deltawaves = np.hstack(deltawaves)[srt]
            sigmas = np.hstack(sigmas)[srt]
            linenames = np.hstack(linenames)[srt]
        
            # Add the linenames to the spectrum plot.
            for meanwave, linename in zip(meanwaves*(1+redshift), linenames):
                #print(meanwave, ymax_spec)
                if meanwave > spec_wavelims[0] and meanwave < spec_wavelims[1]:
                    if 'SiIII' in linename:
                        thislinename = '\n'+linename.replace('+', '+\n  ')
                    elif '4363' in linename:
                        thislinename = linename+'\n'
                    else:
                        thislinename = linename
                    if stackfit:
                        specax.text(meanwave/1e4, spec_ymax*0.97, thislinename, ha='center', va='top',
                                    rotation=270, fontsize=12, alpha=0.5)
                    else:
                        specax.text(meanwave/1e4, spec_ymax, thislinename, ha='center', va='top',
                                    rotation=270, fontsize=12, alpha=0.5)
        
            removelabels = np.ones(nline, bool)
            line_ymin, line_ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
    
            if stackfit:
                ax, irow, colshift = [], 0, 5
            else:
                ax, irow, colshift = [], 4, 5 # skip the gap row
                
            for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):
                icol = iax % nlinecols
                icol += colshift
                if iax > 0 and iax % nlinecols == 0:
                    irow += 1
                #print(iax, irow, icol)
            
                xx = fig.add_subplot(gs[irow, icol])
                ax.append(xx)
            
                wmin = (meanwave - deltawave) * (1+redshift) - 6 * sig * meanwave * (1+redshift) / C_LIGHT
                wmax = (meanwave + deltawave) * (1+redshift) + 6 * sig * meanwave * (1+redshift) / C_LIGHT
                #print(linename, wmin, wmax)
            
                # iterate over cameras
                for icam in np.arange(len(data['cameras'])): # iterate over cameras
                    emlinewave = data['wave'][icam]
                    emlineflux = data['flux'][icam] - desicontinuum[icam] - desismoothcontinuum[icam]
                    emlinemodel = desiemlines[icam]
            
                    emlinesigma, good = ivar2var(data['ivar'][icam], sigma=True, allmasked_ok=True, clip=0)
                    emlinewave = emlinewave[good]
                    emlineflux = emlineflux[good]
                    emlinesigma = emlinesigma[good]
                    emlinemodel = emlinemodel[good]
            
                    #if icam == 0:
                    #    import matplotlib.pyplot as plt ; plt.clf() ; plt.plot(emlinewave, emlineflux) ; plt.plot(emlinewave, emlinemodel) ; plt.xlim(4180, 4210) ; plt.ylim(-15, 17) ; plt.savefig('desi-users/ioannis/tmp/junkg.png')
                        
                    emlinemodel_oneline = []
                    for desiemlines_oneline1 in desiemlines_oneline:
                        emlinemodel_oneline.append(desiemlines_oneline1[icam][good])
                        
                    indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
                    if len(indx) > 1:
                        removelabels[iax] = False
                        xx.plot(emlinewave[indx]/1e4, emlineflux[indx], color=col1[icam], alpha=0.5)
                        #xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                        #                emlineflux[indx]+emlinesigma[indx], color=col1[icam], alpha=0.5)
                        # plot the individual lines first then the total model
                        for emlinemodel_oneline1 in emlinemodel_oneline:
                            if np.sum(emlinemodel_oneline1[indx]) > 0:
                                #P = emlinemodel_oneline1[indx] > 0
                                xx.plot(emlinewave[indx]/1e4, emlinemodel_oneline1[indx], lw=1, alpha=0.8, color=col2[icam])
                        xx.plot(emlinewave[indx]/1e4, emlinemodel[indx], color=col2[icam], lw=2, alpha=0.8)
        
                        #xx.plot(emlinewave[indx], emlineflux[indx]-emlinemodel[indx], color='gray', alpha=0.3)
                        #xx.axhline(y=0, color='gray', ls='--')
            
                        # get the robust range
                        sigflux = np.std(emlineflux[indx])
                        filtflux = median_filter(emlineflux[indx], 3, mode='nearest')
            
                        #_line_ymin, _line_ymax = -1.5 * sigflux, 4 * sigflux
                        #if np.max(emlinemodel[indx]) > _line_ymax:
                        #    _line_ymax = np.max(emlinemodel[indx]) * 1.3
                        _line_ymin, _line_ymax = -1.5 * sigflux, np.max(emlinemodel[indx]) * 1.4
                        if 4 * sigflux > _line_ymax:
                            _line_ymax = 4 * sigflux
                        if np.max(filtflux) > _line_ymax:
                            _line_ymax = np.max(filtflux)
                        if np.min(emlinemodel[indx]) < _line_ymin:
                            _line_ymin = 0.8 * np.min(emlinemodel[indx])
                        if _line_ymax > line_ymax[iax]:
                            line_ymax[iax] = _line_ymax
                        if _line_ymin < line_ymin[iax]:
                            line_ymin[iax] = _line_ymin
                        #print(linename, line_ymin[iax], line_ymax[iax])
                        #if linename == '[OII] $\lambda\lambda$3726,29':
                        #    pdb.set_trace()
            
                        xx.set_xlim(wmin/1e4, wmax/1e4)
                        
                    xx.text(0.03, 0.89, linename, ha='left', va='center',
                            transform=xx.transAxes, fontsize=12)
                    xx.tick_params(axis='x', labelsize=16)
                    xx.tick_params(axis='y', labelsize=16)
                    
            for iax, xx in enumerate(ax):
                if removelabels[iax]:
                    xx.set_ylim(0, 1)
                    xx.set_xticklabels([])
                    xx.set_yticklabels([])
                else:
                    xx.set_yticklabels([])
                    xx.set_ylim(line_ymin[iax], line_ymax[iax])
                    xx_twin = xx.twinx()
                    xx_twin.set_ylim(line_ymin[iax], line_ymax[iax])
                    xlim = xx.get_xlim()
                    xx.xaxis.set_major_locator(ticker.MaxNLocator(2))

    if fastphot:
        plt.subplots_adjust(wspace=0.6, top=0.85, bottom=0.13, left=0.07, right=0.86)

    else:
        plt.subplots_adjust(wspace=0.4, top=0.9, bottom=0.1, left=0.07, right=0.92, hspace=0.33)

        # common axis labels
        if len(ax) > 0:
            if len(ax) == 2:
                xx = fig.add_subplot(gs[irow, icol+1])
                xx.axis('off')
                ax.append(xx)
            elif len(ax) == 1:
                xx = fig.add_subplot(gs[irow, icol+1])
                xx.axis('off')
                ax.append(xx)
                xx = fig.add_subplot(gs[irow, icol+2])
                xx.axis('off')
                ax.append(xx)
            
            ulpos = ax[0].get_position()
            lpos = ax[nline-1].get_position()
            urpos = ax[2].get_position()
            xpos = (urpos.x1 - ulpos.x0) / 2 + ulpos.x0# + 0.03
                
            ypos = lpos.y0 - 0.04
            fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                     ha='center', va='center', fontsize=fontsize2)
    
            xpos = urpos.x1 + 0.05
            ypos = (urpos.y1 - lpos.y0) / 2 + lpos.y0# + 0.03
            fig.text(xpos, ypos, 
                     r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                     ha='center', va='center', rotation=270, fontsize=fontsize2)

    legkeys = leg.keys()
    legfntsz = 17

    # rest wavelength plus targeting information
    if fastphot:
        ppos = sedax.get_position()
        xpos = (ppos.x1 - ppos.x0) / 2 + ppos.x0
        ypos = ppos.y1 + 0.07

        fig.text(xpos, ypos, r'Rest-frame Wavelength ($\mu$m)',
                 ha='center', va='bottom', fontsize=fontsize2)
        fig.text(0.58, 0.87, '\n'.join(target), ha='left', va='bottom',
                 fontsize=fontsize2, linespacing=1.4)

        toppos, leftpos = 0.25, 0.59

        txt = [
            r'{}'.format(leg['z']),
        ]
        if 'zredrock' in legkeys:
            txt += [r'{}'.format(leg['zredrock'])]
        txt += [            
            #r'{}'.format(leg['zwarn']),
            #'',
            r'{}'.format(leg['vdisp']),
            '',
            r'{}'.format(leg['dn4000_model']),
        ]
    
        fig.text(leftpos, toppos, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz, 
                 bbox=bbox, linespacing=1.4)

        txt = [            
            r'{}'.format(leg['absmag_r']),
            r'{}'.format(leg['absmag_gr']),
            r'{}'.format(leg['absmag_rz']), 
        ]
    
        fig.text(leftpos+0.18, toppos, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz, 
                 bbox=bbox, linespacing=1.4)
        
    else:
        if stackfit:
            ppos = specax.get_position()
            xpos = (ppos.x1 - ppos.x0) / 2 + ppos.x0
            ypos = ppos.y1 + 0.05

            fig.text(xpos, ypos, r'Rest-frame Wavelength ($\mu$m)',
                     ha='center', va='bottom', fontsize=fontsize2)
            fig.text(0.62, 0.91, '\n'.join(target), ha='left', va='bottom',
                     fontsize=fontsize2, linespacing=1.4)
        else:
            ppos = sedax.get_position()
            xpos = (ppos.x1 - ppos.x0) / 2 + ppos.x0
            ypos = ppos.y1 + 0.03

            fig.text(xpos, ypos, r'Rest-frame Wavelength ($\mu$m)',
                     ha='center', va='bottom', fontsize=fontsize2)
            fig.text(0.647, 0.925, '\n'.join(target), ha='left', va='bottom',
                     fontsize=fontsize2, linespacing=1.4)
        
        # add some key results about the object at the bottom of the figure

        if stackfit:
            #toppos, startpos, deltapos = 0.21, 0.04, 0.13
            toppos, leftpos, rightpos, adjust = 0.27, 0.05, 0.62, 0.01
        else:
            toppos, leftpos, rightpos, adjust = 0.21, 0.03, 0.62, 0.01
        
        nbox = 2 + 1*bool(leg_narrow) + bool(leg_broad)
        boxpos = np.arange(nbox) * (rightpos - leftpos)/nbox + leftpos
    
        txt = [
            r'{}'.format(leg['z']),
        ]
        if 'zredrock' in legkeys:
            txt += [r'{}'.format(leg['zredrock'])]
    
        txt += [            
            #r'{}'.format(leg['zwarn']),
            r'{}'.format(leg['vdisp']),
        ]
        if 'dv_narrow' in legkeys or 'dv_uv' in leg_uv.keys() or 'dv_broad' in leg_broad.keys():
            txt += ['']
    
        if 'dv_narrow' in legkeys:
            txt += [
                r'{}'.format(leg['sigma_narrow']),
                r'{}'.format(leg['dv_narrow']),
            ]
        if 'dv_uv' in leg_uv.keys():
            txt += [
                r'{}'.format(leg_uv['sigma_uv']),
                r'{}'.format(leg_uv['dv_uv']),
            ]
            _ = leg_uv.pop('sigma_uv')
            _ = leg_uv.pop('dv_uv')
        if 'dv_broad' in leg_broad.keys():
            txt += [
                r'{}'.format(leg_broad['sigma_broad']),
                r'{}'.format(leg_broad['dv_broad']),
            ]
            _ = leg_broad.pop('sigma_broad')
            _ = leg_broad.pop('dv_broad')
    
        ibox = 0
        #fig.text(startpos, toppos, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz, 
        fig.text(boxpos[ibox], toppos, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz, 
                 bbox=bbox, linespacing=1.4)
        ibox += 1
    
        txt = [
            r'{}'.format(leg['absmag_r']),
            r'{}'.format(leg['absmag_gr']),
            r'{}'.format(leg['absmag_rz']), 
            '',
            r'{}'.format(leg['dn4000_model'])
        ]
        if 'dn4000_spec' in legkeys:
            txt += [r'{}'.format(leg['dn4000_spec'])]
            
        #fig.text(startpos+deltapos, toppos, '\n'.join(txt), ha='left', va='top', 
        fig.text(boxpos[ibox], toppos, '\n'.join(txt), ha='left', va='top', 
                 fontsize=legfntsz, bbox=bbox, linespacing=1.4)
        ibox += 1        
    
        #factor = 2
        if bool(leg_narrow):
            txt = []
            for key in leg_narrow.keys():
                txt += [r'{}'.format(leg_narrow[key])]
            #fig.text(startpos+deltapos*factor, toppos, '\n'.join(txt), ha='left', va='top',
            fig.text(boxpos[ibox]-adjust*2, toppos, '\n'.join(txt), ha='left', va='top',
                     fontsize=legfntsz, bbox=bbox, linespacing=1.4)
            ibox += 1        
            #factor += 1.25
    
        if bool(leg_broad):
            txt = []
            for key in leg_broad.keys():
                txt += [r'{}'.format(leg_broad[key])]
    
        if bool(leg_uv):
            if bool(leg_broad):
                txt += ['']
            else:
                txt = []
            for key in leg_uv.keys():
                txt += [r'{}'.format(leg_uv[key])]
    
        if bool(leg_uv) or bool(leg_broad):
            #fig.text(startpos+deltapos*factor, toppos, '\n'.join(txt), ha='left', va='top',
            fig.text(boxpos[ibox]-adjust*1, toppos, '\n'.join(txt), ha='left', va='top',
                     fontsize=legfntsz, bbox=bbox, linespacing=1.4)

    log.info('Writing {}'.format(pngfile))
    fig.savefig(pngfile)#, dpi=150)
    plt.close()
