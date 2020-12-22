#!/usr/bin/env python
"""
fastspecfit.external.desi
=========================

fastspecfit wrapper for DESI

ToDo:
* Generalize to work on coadded spectra.
* Fit to the photometry.
* Solve for vdisp, zcontinuum, and E(B-V).
* Add polynomial correction templates
* Capture ivar=0 problems.
* Correct for MW extinction.
* Fit Mg II.

ELG - 20200228 70005
BGS+MWS - 20200303 70500
BGS+MWS - 20200315 66003

new truth table for 70500 is here--
/global/cfs/cdirs/desi/sv/vi/TruthTables/Andes_reinspection/BGS

https://desi.lbl.gov/trac/wiki/SurveyValidation/TruthTables
https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=5720;filename=DESI_data_042820.pdf

time fastspecfit --tile 70500 --night 20200303 --nproc 32 --overwrite

time fastspecfit --tile 67230 --night 20200315 --nproc 32 --overwrite # ELG
time fastspecfit --tile 68001 --night 20200315 --nproc 32 --overwrite # QSO+LRG
time fastspecfit --tile 68002 --night 20200315 --nproc 32 --overwrite # QSO+LRG
time fastspecfit --tile 67142 --night 20200315 --nproc 32 --overwrite # ELG
time fastspecfit --tile 66003 --night 20200315 --nproc 32 --overwrite # BGS+MWS

time fastspecfit-qa --tile 67142 --night 20200315 --first 0 --last 10

"""
import pdb # for debugging

import os, sys, time
import numpy as np
import multiprocessing

from desiutil.log import get_logger
log = get_logger()

# ridiculousness!
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')

def _fastspecfit_one(args):
    """Multiprocessing wrapper."""
    return fastspecfit_one(*args)

def fastspecfit_one(indx, data, CFit, EMFit, out):
    """Fit one."""
    log.info('Continuum-fitting object {}'.format(indx))
    t0 = time.time()
    cfit, continuum = CFit.fnnls_continuum(data)
    log.info('Continuum-fitting object {} took {:.2f} sec'.format(indx, time.time()-t0))

    # fit the emission-line spectrum
    t0 = time.time()
    emfit = EMFit.fit(data, continuum)
    log.info('Line-fitting object {} took {:.2f} sec'.format(indx, time.time()-t0))

    for col in emfit.colnames:
        out[col] = emfit[col]
    for col in cfit.colnames:
        out[col] = cfit[col]
        
    return out

def read_spectra(tile, night, specprod='andes', use_vi=False, write_spectra=True,
                 overwrite=False, verbose=False):
    """Read the spectra and redshift catalog for a given tile and night.

    Parameters
    ----------
    tile : :class:`str`
        Tile number to analyze.
    night : :class:`str`
        Night on which `tile` was observed.
    specprod : :class:`str`, defaults to `andes`.
        Spectroscopic production to read.
    use_vi : :class:`bool`, optional, defaults to False
        Select the subset of spectra with high-quality visual inspections.
    write_spectra : :class:`bool`, optional, defaults to True
        Write out the selected spectra (useful for testing.)
    overwrite : :class:`bool`, optional, defaults to False
        Overwrite output files if they exist on-disk.
    verbose : :class:`bool`, optional, defaults to False
        Trigger more verbose output.

    Returns
    -------
    :class:`astropy.table.Table`
        Redrock redshift table for the given tile and night.
    :class:`desispec.spectra.Spectra`
        DESI spectra for the given tile and night (in the standard format).

    Notes
    -----
    The spectra from all 10 spectrographs are combined and only the subset of
    galaxy spectra with good redshifts (and, optionally, high-quality visual
    inspections) are returned.

    """
    import fitsio
    from astropy.table import Table, vstack
    from desispec.spectra import Spectra
    import desispec.io

    data_dir = os.path.join(os.getenv('FASTSPECFIT_DATA'), 'spectra', specprod)
    zbestoutfile = os.path.join(data_dir, 'zbest-{}-{}.fits'.format(tile, night))
    coaddoutfile = os.path.join(data_dir, 'coadd-{}-{}.fits'.format(tile, night))
    if os.path.isfile(coaddoutfile) and not overwrite:
        zbest = Table(fitsio.read(zbestoutfile))
        log.info('Read {} redshifts from {}'.format(len(zbest), zbestoutfile))

        coadd = desispec.io.read_spectra(coaddoutfile)
        log.info('Read {} spectra from {}'.format(len(zbest), coaddoutfile))

        assert(np.all(zbest['TARGETID'] == coadd.fibermap['TARGETID']))

        return zbest, coadd

    log.info('Parsing tile, night {}, {} from specprod {}'.format(tile, night, specprod))

    specprod_dir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', specprod)
    desidatadir = os.path.join(specprod_dir, 'tiles', tile, night)

    # use visual inspection results?
    if use_vi:
        truthdir = os.path.join(os.getenv('DESI_ROOT'), 'sv', 'vi', 'TruthTables')
        tile2file = {
            '66003': os.path.join(truthdir, 'truth_table_BGS_v1.2.csv'),
            '70500': os.path.join(truthdir, 'truth_table_ELG_v1.2_latest.csv'),
            }
        if tile not in tile2file.keys():
            log.warning('No (known) truth file exists for tile {}. Setting use_vi=False.'.format(tile))
            use_vi = False
        else:
            truth = Table.read(truthfile)
            best = np.where(
                (truth['best quality'] >= 2.5) * 
                (truth['Redrock spectype'] == 'GALAXY')# *
                #(truth['Redrock z'] < 0.75) 
            )[0]
            #goodobj = np.where((zb['DELTACHI2'] > 100) * (zb['ZWARN'] == 0) * (zb['SPECTYPE'] == 'GALAXY'))[0]

            log.info('Read {}/{} good redshifts from {}'.format(len(best), len(truth), truthfile))
            truth = truth[best]
    
    zbest = []
    spectra, keepindx = [], []
    #for spectro in ('0'):
    for spectro in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
        zbestfile = os.path.join(desidatadir, 'zbest-{}-{}-{}.fits'.format(spectro, tile, night))
        coaddfile = os.path.join(desidatadir, 'coadd-{}-{}-{}.fits'.format(spectro, tile, night))
        if os.path.isfile(zbestfile) and os.path.isfile(coaddfile):
            zb = Table(fitsio.read(zbestfile))
            fmap = fitsio.read(coaddfile, ext='FIBERMAP', columns=['FIBERSTATUS', 'OBJTYPE'])
            if use_vi:
                keep = np.where(np.isin(zb['TARGETID'], truth['TARGETID']))[0]
            else:
                keep = np.where((zb['Z'] > 0) * (zb['ZWARN'] == 0) * (fmap['OBJTYPE'] == 'TGT') *
                                (zb['SPECTYPE'] == 'GALAXY') * (fmap['FIBERSTATUS'] == 0))[0]

            if verbose:
                log.debug('Spectrograph {}: N={}'.format(spectro, len(keep)))
            if len(keep) > 0:
                zbest.append(zb[keep])
                keepindx.append(keep)
                spectra.append(desispec.io.read_spectra(coaddfile))

    if len(zbest) == 0:
        log.fatal('No spectra found!')
        raise ValueError
        
    # combine the spectrographs
    #log.debug('Stacking all the spectra.')
    zbest = vstack(zbest)

    coadd = None
    for camera in ('b', 'r', 'z'):
        wave = spectra[0].wave[camera]
        fm, flux, ivar, mask, res = [], [], [], [], []
        for ii in np.arange(len(spectra)):
            fm.append(spectra[ii].fibermap[keepindx[ii]])
            flux.append(spectra[ii].flux[camera][keepindx[ii], :])
            ivar.append(spectra[ii].ivar[camera][keepindx[ii], :])
            mask.append(spectra[ii].mask[camera][keepindx[ii], :])
            res.append(spectra[ii].resolution_data[camera][keepindx[ii], :, :])

        fm = vstack(fm)
        flux = np.concatenate(flux)
        ivar = np.concatenate(ivar)
        mask = np.concatenate(mask)
        res = np.concatenate(res)

        _coadd = Spectra([camera], {camera: wave}, {camera: flux}, {camera : ivar}, 
                         resolution_data={camera: res}, mask={camera: mask}, 
                         fibermap=fm, single=False)#, meta=meta) 
        if coadd is None:
            coadd = _coadd
        else:
            coadd.update(_coadd)

    assert(np.all(zbest['TARGETID'] == coadd.fibermap['TARGETID']))
    
    log.info('Writing {} redshifts to {}'.format(len(zbest), zbestoutfile))
    zbest.write(zbestoutfile, overwrite=True)

    log.info('Writing {} spectra to {}'.format(len(zbest), coaddoutfile))
    desispec.io.write_spectra(coaddoutfile, coadd)

    return zbest, coadd

def _unpack_one_spectrum(args):
    """Multiprocessing wrapper."""
    return unpack_one_spectrum(*args)

def unpack_one_spectrum(specobj, zbest, CFit, indx):
    """Unpack and pre-process a single DESI spectrum.

    Parameters
    ----------
    specobj : :class:`desispec.spectra.Spectra`
        DESI spectra (in the standard format).
    zbest : :class:`astropy.table.Table`
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
    from desispec.resolution import Resolution
    from desiutil.dust import ext_odonnell
    from desitarget.io import desitarget_resolve_dec
    from fastspecfit.util import C_LIGHT

    cameras = ['b', 'r', 'z']
    ncam = len(cameras)

    ra = specobj.fibermap['TARGET_RA'][indx]
    dec = specobj.fibermap['TARGET_DEC'][indx]
    ebv = CFit.SFDMap.ebv(ra, dec, scaling=1.0) # SFD coefficients

    # Unpack the data and correct for Galactic extinction. Also flag pixels that
    # may be affected by emission lines.
    data = {'zredrock': zbest['Z'][indx], 'wave': [], 'flux': [], 'ivar': [],
            'res': [], 'linemask': [], 'snr': np.zeros(3).astype('f4')}
    for icam, camera in enumerate(cameras):
        dust = 10**(-0.4 * ebv * CFit.RV * ext_odonnell(specobj.wave[camera], Rv=CFit.RV))       
        data['wave'].append(specobj.wave[camera])
        data['flux'].append(specobj.flux[camera][indx, :] * dust)
        data['ivar'].append(specobj.ivar[camera][indx, :] / dust**2)
        data['res'].append(Resolution(specobj.resolution_data[camera][indx, :, :]))
        data['snr'][icam] = np.median(specobj.flux[camera][indx, :] * np.sqrt(specobj.ivar[camera][indx, :]))

        linemask = np.ones_like(specobj.wave[camera], bool)
        for line in CFit.linetable:
            zwave = line['restwave'] * (1 + data['zredrock'])
            I = np.where((specobj.wave[camera] >= (zwave - 1.5*CFit.linemask_sigma * zwave / C_LIGHT)) *
                         (specobj.wave[camera] <= (zwave + 1.5*CFit.linemask_sigma * zwave / C_LIGHT)))[0]
            if len(I) > 0:
                linemask[I] = False
        data['linemask'].append(linemask)

    # Make a quick coadd using inverse variance weights.
    uspecwave = np.unique(np.hstack(data['wave']))
    uspecflux3d = np.zeros((len(uspecwave), 3))
    uspecivar3d = np.zeros_like(uspecflux3d)
    for icam in np.arange(ncam):
        I = np.where(np.isin(data['wave'][icam], uspecwave))[0]
        J = np.where(np.isin(uspecwave, data['wave'][icam]))[0]
        uspecflux3d[J, icam] = data['flux'][icam][I]
        uspecivar3d[J, icam] = data['ivar'][icam][I]

    uspecivar = np.sum(uspecivar3d, axis=1)
    uspecflux = np.zeros_like(uspecivar)
    good = np.where(uspecivar > 0)[0]
    uspecflux[good] = np.sum(uspecivar3d[good, :] * uspecflux3d[good, :], axis=1) / uspecivar[good]
    data.update({'coadd_wave': uspecwave, 'coadd_flux': uspecflux, 'coadd_ivar': uspecivar})
    del uspecwave, uspecivar, uspecflux

    # Synthesize photometry, resolved into north and south.
    split = desitarget_resolve_dec()
    isouth = specobj.fibermap['TARGET_DEC'][indx] < split
    if isouth:
        filters = CFit.decamwise
    else:
        filters = CFit.bassmzlswise
    lambda_eff = filters.effective_wavelengths.value
    data['photsys_south'] = isouth

    padflux, padwave = filters.pad_spectrum(data['coadd_flux'], data['coadd_wave'], method='edge')
    synthmaggies = filters.get_ab_maggies(1e-17 * padflux, padwave)
    synthmaggies = synthmaggies.as_array().view('f8')[:3] # keep just grz

    # code to synthesize uncertainties from the variance spectrum
    #var, mask = _ivar2var(data['coadd_ivar'])
    #padvar, padwave = filters.pad_spectrum(var[mask], data['coadd_wave'][mask], method='edge')
    #synthvarmaggies = filters.get_ab_maggies(1e-17**2 * padvar, padwave)
    #synthivarmaggies = 1 / synthvarmaggies.as_array().view('f8')[:3] # keep just grz
    #
    #data['synthphot'] = CFit.parse_photometry(
    #    maggies=synthmaggies, lambda_eff=lambda_eff[:3],
    #    ivarmaggies=synthivarmaggies, nanomaggies=False)

    data['synthphot'] = CFit.parse_photometry(
        maggies=synthmaggies, lambda_eff=lambda_eff[:3],
        nanomaggies=False)

    # Unpack the imaging photometry.
    bands = ['g', 'r', 'z', 'W1', 'W2']
    maggies = np.zeros(len(bands))
    ivarmaggies = np.zeros(len(bands))
    for iband, band in enumerate(bands):
        ivar = specobj.fibermap['FLUX_IVAR_{}'.format(band.upper())][indx]
        if ivar > 0:
            dust = specobj.fibermap['MW_TRANSMISSION_{}'.format(band.upper())][indx]
            maggies[iband] = specobj.fibermap['FLUX_{}'.format(band.upper())][indx] / dust
            ivarmaggies[iband] = ivar * dust**2
    
    data['phot'] = CFit.parse_photometry(
        maggies=maggies, lambda_eff=lambda_eff,
        ivarmaggies=ivarmaggies, nanomaggies=True)

    return data

def unpack_all_spectra(specobj, zbest, CFit, fitindx, nproc=1):
    """Wrapper on unpack_one_spectrum to parse all the input spectra.

    Parameters
    ----------
    specobj : :class:`desispec.spectra.Spectra`
        DESI spectra (in the standard format).
    zbest : :class:`astropy.table.Table`
        Redrock redshift table (row-aligned to `specobj`).
    CFit : :class:`fastspecfit.continuum.ContinuumFit`
        Continuum-fitting class.
    fitindx : :class:`int`
        Index numbers of the spectra to unpack and pre-process.
    nproc : :class:`int`, optional, defaults to 1
        Number of cores to use for multiprocessing.
    
    Returns
    -------
    :class:`list`
        List of dictionaries (row-aligned to `fitindx`) populated by
        `unpack_one_spectrum`.
        
    Notes
    -----

    """
    args = [(specobj, zbest, CFit, indx) for indx in fitindx]
    if nproc > 1:
        with multiprocessing.Pool(nproc) as P:
            data = np.vstack(P.map(_unpack_one_spectrum, args))
    else:
        data = [unpack_one_spectrum(*_args) for _args in args]

    return data    

def init_output(tile, night, zbest, fibermap, CFit, EMFit):
    """Initialize the fastspecfit output data table.

    Parameters
    ----------
    tile : :class:`str`
        Tile number.
    night : :class:`str`
        Night on which `tile` was observed.
    zbest : :class:`astropy.table.Table`
        Redrock redshift table (row-aligned to `fibermap`).
    fibermap : :class:`astropy.table.Table`
        Fiber map (row-aligned to `zbest`).
    CFit : :class:`fastspecfit.continuum.ContinuumFit`
        Continuum-fitting class.
    EMFit : :class:`fastspecfit.emlines.EMLineFit`
        Emission-line fitting class.

    Returns
    -------
    

    Notes
    -----

    """
    from astropy.table import Table, Column, hstack

    # Grab info on the emission lines and the continuum.
    nobj = len(zbest)

    out = Table()
    for zbestcol in ['TARGETID', 'Z']:#, 'ZERR']:#, 'SPECTYPE', 'DELTACHI2']
        out[zbestcol] = zbest[zbestcol]
    for fibermapcol in ['FIBER']:
        out[fibermapcol] = fibermap[fibermapcol]
    out.add_column(Column(name='NIGHT', data=np.repeat(night, nobj)), index=0)
    out.add_column(Column(name='TILE', data=np.repeat(tile, nobj)), index=0)

    out = hstack((out, CFit.init_output(nobj), EMFit.init_output(CFit.linetable, nobj)))

    return out

def parse(options=None):
    """Parse input arguments.

    """
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required but with sensible defaults
    parser.add_argument('--night', default='20200225', type=str, help='Night to process.')
    parser.add_argument('--tile', default='70502', type=str, help='Tile number to process.')

    # optional inputs
    parser.add_argument('--first', type=int, help='Index of first spectrum to process (0-indexed).')
    parser.add_argument('--last', type=int, help='Index of last spectrum to process (max of nobj-1).')
    parser.add_argument('--nproc', default=1, type=int, help='Number of cores.')
    parser.add_argument('--specprod', type=str, default='variance-model', choices=['andes', 'daily', 'variance-model'],
                        help='Spectroscopic production to process.')

    parser.add_argument('--use-vi', action='store_true', help='Select spectra with high-quality visual inspections (VI).')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite any existing files.')
    parser.add_argument('--no-write-spectra', dest='write_spectra', default=True, action='store_false',
                        help='Do not write out the selected spectra for the specified tile and night.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose.')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspecfit {}'.format(' '.join(options)))

    return args

def main(args=None, comm=None):
    """Main module.

    """
    from astropy.table import vstack
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    # Create directories as needed.
    fastspecfit_dir = os.getenv('FASTSPECFIT_DATA')
    for subdir in ('spectra', 'results', 'qa'):
        outdir = os.path.join(fastspecfit_dir, subdir, args.specprod)
        if not os.path.isdir(outdir):
            if args.verbose:
                log.debug('Creating directory {}'.format(outdir))
            os.makedirs(outdir, exist_ok=True)
    
    # If the output file exists, we're done!
    resultsdir = os.path.join(fastspecfit_dir, 'results', args.specprod)
    fastspecfitfile = os.path.join(resultsdir, 'fastspecfit-{}-{}.fits'.format(
        args.tile, args.night))
    if os.path.isfile(fastspecfitfile) and not args.overwrite:
        log.info('Output file {} exists; all done!'.format(fastspecfitfile))
        return
        
    # Read the data 
    t0 = time.time()
    zbest, specobj = read_spectra(tile=args.tile, night=args.night,
                                  specprod=args.specprod,
                                  use_vi=args.use_vi, 
                                  write_spectra=args.write_spectra,
                                  verbose=args.verbose)
    log.info('Reading the data took: {:.2f} sec'.format(time.time()-t0))

    if args.first is None:
        args.first = 0
    if args.last is None:
        args.last = len(zbest) - 1
    if args.first > args.last:
        log.warning('Option --first cannot be larger than --last!')
        raise ValueError
    fitindx = np.arange(args.last - args.first + 1) + args.first

    # Initialize the continuum- and emission-line fitting classes and the output
    # data table.
    t0 = time.time()
    CFit = ContinuumFit(verbose=args.verbose)
    EMFit = EMLineFit(verbose=args.verbose)
    fastspecfit = init_output(args.tile, args.night, zbest, specobj.fibermap, CFit, EMFit)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Unpacking with multiprocessing takes a lot longer (maybe pickling takes a
    # long time?) so suppress the `nproc` argument here for now.
    t0 = time.time()
    data = unpack_all_spectra(specobj, zbest, CFit, fitindx)#, nproc=args.nproc)
    del specobj, zbest # free memory
    log.info('Unpacking the spectra to be fitted took: {:.2f} sec'.format(time.time()-t0))

    # Fit in parallel
    fitargs = [(indx, data[iobj], CFit, EMFit, fastspecfit[indx]) for iobj, indx in enumerate(fitindx)]
    if args.nproc > 1:
        with multiprocessing.Pool(args.nproc) as P:
            out = P.map(_fastspecfit_one, fitargs)
    else:
        out = [fastspecfit_one(*_fitargs) for _fitargs in fitargs]

    fastspecfit = vstack(out)

    # write out
    t0 = time.time()
    log.info('Writing results for {} objects to {}'.format(len(fastspecfit), fastspecfitfile))
    fastspecfit.write(fastspecfitfile, overwrite=True)
    log.info('Writing out took {:.2f} sec'.format(time.time()-t0))
