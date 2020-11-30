"""
desigal.nyxgalaxy
=================

"""
import pdb # for debugging

import os
import numpy as np

import astropy.units as u
from astropy.table import Table, Column, vstack, join
from astropy.modeling import Fittable1DModel

from desiutil.log import get_logger
from desispec.interpolation import resample_flux

from scipy import constants
C_LIGHT = constants.c / 1000.0 # [km/s]

log = get_logger()

def init_nyxgalaxy(tile, night, zbest, CFit):
    """Initialize the output data table.

    CFit - ContinuumFit class

    """
    import astropy.units as u
    nobj = len(zbest)

    # Grab info on the emission lines and the continuum.
    linetable = read_nyxgalaxy_lines()
    nssp_coeff = len(CFit.sspinfo)

    out = Table()
    for zbestcol in ['TARGETID', 'Z']:#, 'ZERR']:#, 'SPECTYPE', 'DELTACHI2']
        out[zbestcol] = zbest[zbestcol]
    out.add_column(Column(name='NIGHT', data=np.repeat(night, nobj)), index=0)
    out.add_column(Column(name='TILE', data=np.repeat(tile, nobj)), index=0)
    out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
    out.add_column(Column(name='CONTINUUM_CHI2', length=nobj, dtype='f4')) # reduced chi2
    out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
    out.add_column(Column(name='CONTINUUM_PHOT_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
    out.add_column(Column(name='CONTINUUM_PHOT_CHI2', length=nobj, dtype='f4')) # reduced chi2
    out.add_column(Column(name='CONTINUUM_PHOT_DOF', length=nobj, dtype=np.int32))
    out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8'))
    out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
    out.add_column(Column(name='CONTINUUM_EBV', length=nobj, dtype='f4', unit=u.mag))
    out.add_column(Column(name='CONTINUUM_VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
    out.add_column(Column(name='D4000', length=nobj, dtype='f4'))
    out.add_column(Column(name='D4000_IVAR', length=nobj, dtype='f4'))
    out.add_column(Column(name='D4000_MODEL', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEVSHIFT_FORBIDDEN', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEVSHIFT_FORBIDDEN_IVAR', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEVSHIFT_BALMER', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEVSHIFT_BALMER_IVAR', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINESIGMA_FORBIDDEN', length=nobj, dtype='f4', unit=u.kilometer / u.second))
    out.add_column(Column(name='LINESIGMA_FORBIDDEN_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
    out.add_column(Column(name='LINESIGMA_BALMER', length=nobj, dtype='f4', unit=u.kilometer / u.second))
    out.add_column(Column(name='LINESIGMA_BALMER_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))

    for line in linetable['name']:
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
        out.add_column(Column(name='{}_CHI2'.format(line), length=nobj, dtype='f4'))
        out.add_column(Column(name='{}_NPIX'.format(line), length=nobj, dtype=np.int32))
        
    return out
    
def read_spectra(tile, night, use_vi=False, write_spectra=True, verbose=False):
    """Read the spectra and redshift catalog for a given tile and night.

    Parameters
    ----------
    tile : :class:`str`
        Tile number to analyze.
    night : :class:`str`
        Night on which `tile` was observed.
    use_vi : :class:`bool`, optional, defaults to False
        Select the subset of spectra with high-quality visual inspections.
    write_spectra : :class:`bool`, optional, defaults to True
        Write out the selected spectra. (Useful for testing.)
    verbose : :class:`bool`, optional, defaults to False
        Trigger more verbose output.

    Returns
    -------
    :class:`astropy.table.Table`
        Redrock redshift table.
    :class:`desispec.spectra.Spectra`
        DESI spectra (in the standard format).

    Notes
    -----
    The spectra from all 10 spectrographs are combined and only the subset of
    galaxy spectra with good redshifts (and, optionally, high-quality visual
    inspections) are returned.

    An optional `overwrite` input could be added to overwrite existing spectra.

    """
    import fitsio
    from astropy.table import join
    from desispec.spectra import Spectra
    import desispec.io

    data_dir = os.path.join(os.getenv('NYXGALAXY_DATA'), 'spectra')
    if not os.path.isdir(data_dir):
        if verbose:
            log.debug('Creating directory {}'.format(data_dir))
        os.makedirs(data_dir, exist_ok=True)
    
    zbestoutfile = os.path.join(data_dir, 'zbest-{}-{}.fits'.format(tile, night))
    coaddoutfile = os.path.join(data_dir, 'coadd-{}-{}.fits'.format(tile, night))
    if os.path.isfile(coaddoutfile) and not overwrite:
        zbest = Table(fitsio.read(zbestoutfile))
        log.info('Read {} redshifts from {}'.format(len(zbest), zbestoutfile))

        coadd = desispec.io.read_spectra(coaddoutfile)
        log.info('Read {} spectra from {}'.format(len(zbest), coaddoutfile))

        return zbest, coadd

    log.info('Parsing data from tile, night {}, {}'.format(tile, night))

    specprod_dir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', 'andes')
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
            if use_vi:
                keep = np.where(np.isin(zb['TARGETID'], truth['TARGETID']))[0]
            else:
                #snr = np.median(specobj.flux['r']*np.sqrt(specobj.ivar['r']), axis=1)
                #keep = np.where((zb['ZWARN'] == 0) * (zb['DELTACHI2'] > 50) * (zb['SPECTYPE'] == 'GALAXY'))[0]
                keep = np.where((zb['ZWARN'] == 0) * (zb['DELTACHI2'] > 50) * (zb['SPECTYPE'] == 'GALAXY'))[0]

            if verbose:
                log.debug('Spectrograph {}: N={}'.format(spectro, len(keep)))
            if len(keep) > 0:
                zbest.append(zb[keep])
                keepindx.append(keep)
                spectra.append(desispec.io.read_spectra(coaddfile))

    if len(zbest) == 0:
        log.fatal('No spectra found for tile {} and night {}!'.format(tile, night))
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
                    fibermap=fm, single=True)#, meta=meta)    
        if coadd is None:
            coadd = _coadd
        else:
            coadd.update(_coadd)

    log.info('Writing {} redshifts to {}'.format(len(zbest), zbestoutfile))
    zbest.write(zbestoutfile, overwrite=True)

    log.info('Writing {} spectra to {}'.format(len(zbest), coaddoutfile))
    desispec.io.write_spectra(coaddoutfile, coadd)

    return zbest, coadd

def _unpack_spectrum(specobj, zbest, iobj, CFit, south=True):
    """Unpack a single spectrum into a list.

    """
    from desispec.resolution import Resolution
    
    if south:
        filters = CFit.decamwise
    else:
        filters = CFit.bassmzlswise

    log.info('need to correct for dust and also maybe pack everything into a dictionary!')
    pdb.set_trace()

    galwave, galflux, galivar, galres = [], [], [], []
    for camera in ('b', 'r', 'z'):
        galwave.append(specobj.wave[camera])
        galflux.append(specobj.flux[camera][iobj, :])
        galivar.append(specobj.ivar[camera][iobj, :])
        galres.append(Resolution(specobj.resolution_data[camera][iobj, :, :]))

    # make a quick coadd using inverse variance weights
    ugalwave = np.unique(np.hstack(galwave))
    ugalflux3d = np.zeros((len(ugalwave), 3))
    ugalivar3d = np.zeros_like(ugalflux3d)
    for icam in [0, 1, 2]:
        I = np.where(np.isin(galwave[icam], ugalwave))[0]
        J = np.where(np.isin(ugalwave, galwave[icam]))[0]
        ugalflux3d[J, icam] = galflux[icam][I]
        ugalivar3d[J, icam] = galivar[icam][I]

    ugalivar = np.sum(ugalivar3d, axis=1)
    ugalflux = np.sum(ugalivar3d * ugalflux3d, axis=1) / ugalivar
    padflux, padwave = filters.pad_spectrum(ugalflux, ugalwave, method='edge')
    filters.get_ab_maggies(1e-17*padflux, padwave)

    if False:
        import matplotlib.pyplot as plt
        for icam in [0, 1, 2]:
            plt.plot(galwave[icam], galflux[icam])
        plt.plot(ugalwave, ugalflux, color='k', alpha=0.7)
        plt.savefig('junk.png')
        pdb.set_trace()

    # unpack the photometry
    fluxcols = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2']
    ivarcols = ['FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2']
    #galphot = Table(specobj.fibermap[[fluxcols, ivarcols][iobj])

    maggies = np.array([specobj.fibermap[col][iobj] for col in fluxcols])
    ivarmaggies = np.array([specobj.fibermap[col][iobj] for col in ivarcols])
    lambda_eff = filters.effective_wavelengths.value
    
    galphot = CFit.convert_photometry(
        maggies=maggies, lambda_eff=lambda_eff,
        ivarmaggies=ivarmaggies, nanomaggies=True)

    zredrock = zbest['Z'][iobj]

    return galwave, galflux, galivar, galres, galphot, zredrock

def _smooth_and_resample(args):
    """Wrapper for multiprocessing."""
    return smooth_and_resample(*args)

def smooth_and_resample(sspflux, sspwave, galwave, galR):
    """
    sspflux[npix] - redshifted SSP template
    sspwave[npix] - redshifted SSP wavelength
    
    """
    if galwave is not None:
        resampflux = resample_flux(galwave, sspwave, sspflux, extrapolate=True)
    else:
        resampflux = sspflux
    if galR is not None:
        smoothflux = galR.dot(resampflux)
    else:
        smoothflux = resampflux
    return smoothflux.T

class ContinuumFit():
    def __init__(self, metallicity='Z0.0190', minwave=0.0, maxwave=6e4, nproc=1,
                 verbose=True):
        """Model the stellar continuum.
    
        specobj - Spectra Class data
        zbest - astropy Table with redshift info
        iobj - index of object to fit
        nproc - number of cores to use for multiprocessing

        * Improved (iterative) continuum-fitting:
          - Implement weighted fitting.
          - Mask pixels around emission lines.
          - Update the continuum redshift using cross-correlation. 
          - Resample to constant log-lamba and solve for the velocity dispersion.
          - Need to be careful because several of these steps will be a function of S/N.
          - The last continuum fit should be a non-linear fit which includes dust attenuation.

          decamwise (speclite.filters instance): DECam2014-[g,r,z] and WISE2010-[W1,W2]
            FilterSequence.
          bassmzlswise (speclite.filters instance): BASS-[g,r], MzLS-z and
            WISE2010-[W1,W2] FilterSequence.
        
        """
        import fitsio
        from astropy.cosmology import FlatLambdaCDM
        from speclite import filters

        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)        

        self.metallicity = metallicity
        self.Z = float(metallicity[1:])
        self.library = 'CKC14z'
        self.isochrone = 'Padova' # would be nice to get MIST in here
        self.imf = 'Kroupa'

        self.nproc = nproc
        self.verbose = verbose

        # Don't hard-code the path!
        self.ssppath = os.getenv('NYXGALAXY_TEMPLATES')
        self.sspfile = os.path.join(self.ssppath, 'SSP_{}_{}_{}_{}.fits'.format(
            self.isochrone, self.library, self.imf, self.metallicity))

        if verbose:
            log.info('Reading {}'.format(self.sspfile))
        wave = fitsio.read(self.sspfile, ext='WAVE')
        flux = fitsio.read(self.sspfile, ext='FLUX')
        sspinfo = Table(fitsio.read(self.sspfile, ext='METADATA'))
        
        # Trim the wavelengths and subselect the number of ages/templates.
        keep = np.where((wave >= minwave) * (wave <= maxwave))[0]
        self.sspwave = wave[keep]
        self.sspflux = flux[keep, ::3]
        self.sspinfo = sspinfo[::3]

        self.nage = len(self.sspinfo['age'])
        self.npix = len(self.sspwave)

        # photometry
        self.decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z',
                                              'wise2010-W1', 'wise2010-W2')
        self.bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z',
                                                 'wise2010-W1', 'wise2010-W2')

        # dust parameters and emission lines
        self.dustslope = 0.7
        
        self.linetable = read_nyxgalaxy_lines()

    @staticmethod
    def convert_photometry(maggies, lambda_eff, ivarmaggies=None, nanomaggies=True,
                           flam=True, fnu=False, abmag=False):
        """Convert maggies to various outputs and pack into a table.

        flam - 10-17 erg/s/cm2/A
        fnu - 10-17 erg/s/cm2/Hz
        abmag - AB mag

        nanomaggies - input maggies are actually 1e-9 maggies

        """
        shp = maggies.shape
        if maggies.ndim == 1:
            nband, ngal = shp[0], 1
        else:
            nband, ngal = shp[0], shp[1]
            
        phot = Table()
        phot.add_column(Column(name='lambda_eff', length=nband, dtype='f4'))
        phot.add_column(Column(name='nanomaggies', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='nanomaggies_ivar', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='flam', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='flam_ivar', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag', length=nband, shape=(ngal, ), dtype='f4'))
        phot.add_column(Column(name='abmag_ivar', length=nband, shape=(ngal, ), dtype='f4'))

        if ivarmaggies is None:
            ivarmaggies = np.zeros_like(maggies)

        phot['lambda_eff'] = lambda_eff.astype('f4')
        if nanomaggies:
            phot['nanomaggies'] = maggies.astype('f4')
            phot['nanomaggies_ivar'] = ivarmaggies.astype('f4')
        else:
            phot['nanomaggies'] = (maggies * 1e9).astype('f4')
            phot['nanomaggies_ivar'] = (ivarmaggies * 1e-18).astype('f4')

        if nanomaggies:
            nanofactor = 1e-9 # [nanomaggies-->maggies]
        else:
            nanofactor = 1.0

        factor = nanofactor * 10**(-0.4 * 48.6) * C_LIGHT * 1e13 / lambda_eff**2 # [maggies-->erg/s/cm2/A]
        if ngal > 1:
            factor = factor[:, None] # broadcast for the models
        phot['flam'] = (maggies * factor).astype('f4')
        phot['flam_ivar'] = (ivarmaggies / factor**2).astype('f4')

        phot['abmag'] = (-2.5 * np.log10(nanofactor * maggies)).astype('f4')

        # approximate the uncertainty as being symmetric in magnitude
        phot['abmag_ivar'] = (ivarmaggies * (maggies * 0.4 * np.log(10))**2).astype('f4')

        return phot
        
    def redshift_smooth_and_resample(self, galwave, galres, redshift, south=True):
        """Redshift, convolve with the spectral resolution, and 
        resample in wavelength.

        ToDo: velocity dispersion smoothing

        phot - photometric table
        
        """
        import multiprocessing

        # Synthesize photometry (only need to do this once).
        if south:
            filters = self.decamwise
        else:
            filters = self.bassmzlswise
        zsspwave = self.sspwave * (1 + redshift)
        zsspflux = self.sspflux / (1 + np.array(redshift).repeat(self.npix)[:, None])

        #zsspflux, zsspwave = list(zip(*args))[:2]
        #zsspflux, zsspwave = np.vstack(zsspflux), zsspwave[0]

        # convert to 10-17 flambda
        maggies = filters.get_ab_maggies(zsspflux.T, zsspwave)
        maggies = np.vstack(maggies.as_array().tolist()).T
        effwave = filters.effective_wavelengths.value

        phot = self.convert_photometry(maggies, effwave, nanomaggies=False)

        # ignore the per-camera spectra
        if galres is None and galwave is None:
            args = [(zsspflux[:, iage], zsspwave, None, None)
                    for iage in np.arange(self.nage)]
            if self.nproc > 1:
                with multiprocessing.Pool(self.nproc) as pool:
                    smoothflux = np.vstack(pool.map(_smooth_and_resample, args)).T
            else:
                smoothflux = np.vstack([smooth_and_resample(*_args) for _args in args]).T
        else:
            # loop over cameras then SSP ages
            smoothflux = []
            for icamera in [0, 1, 2]: # iterate on cameras
                args = [(zsspflux[:, iage], zsspwave, galwave[icamera], galres[icamera])
                        for iage in np.arange(self.nage)]

                if self.nproc > 1:
                    with multiprocessing.Pool(self.nproc) as pool:
                        smoothflux.append(np.vstack(pool.map(_smooth_and_resample, args)).T)
                else:
                    smoothflux.append(np.vstack([smooth_and_resample(*_args) for _args in args]).T)
            
        return smoothflux, phot # [npix, nage]

    def fnnls_continuum_bestfit(self, coeff, sspflux=None, galwave=None,
                                galres=None, redshift=None, south=True):
        if sspflux is None:
            sspflux, sspphot = self.redshift_smooth_and_resample(galwave, galres, redshift, south=south)
            if galres is None and galwave is None: # ignore per-camera
                bestfit = sspflux.dot(coeff)
            else: # iterate over camera
                bestfit = [_sspflux.dot(coeff) for _sspflux in sspflux]
        else:
            bestfit = sspflux.dot(coeff)

        if galres is None and galwave is None:
            return bestfit, self.sspwave
        else:
            return bestfit

    def fnnls_continuum(self, galwave, galflux, galivar, galres,
                        galphot, redshift, linetable, sigma_mask=300.0,
                        use_photometry=True, south=True):
        """Fit the continuum using fast NNLS.
        https://github.com/jvendrow/fnnls

        sigma_mask - emission-line masking sigma [km/s]

        ToDo: mask more emission lines than we fit (e.g., Mg II).
        
        """
        from time import time
        from fnnls import fnnls
        
        sspflux, sspphot = self.redshift_smooth_and_resample(galwave, galres, redshift, south=south) 

        # combine the cameras and fit
        _galwave = np.hstack(galwave)
        _galflux = np.hstack(galflux)
        _galivar = np.hstack(galivar)
        _sspflux = np.concatenate(sspflux, axis=0) # [npix, nage]

        # Mask pixels in and around emission lines.
        emlinemask = np.ones_like(_galivar)
        for line in linetable:
            zwave = line['restwave'] * (1+redshift)
            indx = np.where((_galwave >= (zwave - 1.5*sigma_mask * zwave / C_LIGHT)) *
                            (_galwave <= (zwave + 1.5*sigma_mask * zwave / C_LIGHT)))[0]
            if len(indx) > 0:
                emlinemask[indx] = 0

        # Do a fast initial fit of the stellar continuum.
        log.info('ToDo: add vdisp and ebv, and restrict maxage of templates.')

        # fit with and without photometry
        _modelphot = sspphot['flam'] # [nband, nage]
        _galphot = galphot['flam']
        _galphotivar = galphot['flam_ivar']

        ww = np.sqrt(_galivar * emlinemask)
        ZZ = _sspflux * ww[:, None]
        xx = _galflux * ww

        wwphot = np.sqrt(_galphotivar)
        ZZphot = _modelphot * wwphot[:, None]
        xxphot = _galphot * wwphot

        t0 = time()
        coeff = fnnls(ZZ, xx)[0]
        #coeff = fnnls(np.concatenate((ZZ, ZZphot), axis=0), np.concatenate((xx, xxphot), axis=0))[0]
        coeffphot = fnnls(ZZphot, xxphot)[0]
        dt = time() - t0

        pdb.set_trace()

        # ToDo: fit for dust, the redshift, and velocity dispersion
        zcontinuum, vdisp, ebv = redshift, 0.0, 0.0

        # Need to push the calculation of the best-fitting continuum to a
        # function so we can call it when building QA.
        _continuum, _ = self.fnnls_continuum_bestfit(coeff, _sspflux)
        dof = np.sum(_galivar > 0) - self.nage
        chi2 = np.sum(_galivar * (_galflux - _continuum)**2) / dof

        # Compute the light-weighted age.
        weight = coeff[coeff > 0]
        age = np.sum(weight * self.sspinfo['age'][coeff > 0]) / np.sum(weight) / 1e9 # [Gyr]
        if self.verbose:
            log.info('Continuum fit done in {:.3f} seconds:'.format(dt))
            log.info('  Non-zero templates: {}'.format(len(weight)))
            log.info('  Reduced chi2: {:.3f}'.format(chi2))
            log.info('  Velocity dispersion: {:.3f} km/s'.format(vdisp))
            log.info('  Reddening: {:.3f} mag'.format(ebv))
            log.info('  Light-weighted age: {:.3f} Gyr'.format(age))

        # Unpack the continuum into individual cameras.
        continuum = []
        npixpercamera = [len(gw) for gw in galwave]
        npix = np.hstack([0, npixpercamera])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            continuum.append(_continuum[ipix:jpix])

        # Store the results and return.
        results = {'coeff': coeff, 'chi2': np.float32(chi2), 'dof': dof,
                   'age': np.float32(age), 'ebv': np.float32(ebv),
                   'vdisp': np.float32(vdisp), 'z': zcontinuum,
                   'phot_coeff': coeffphot}
        #results = {'coeff': coeff, 'chi2': np.float32(chi2), 'dof': dof,
        #           'age': np.float32(age)*u.Gyr, 'ebv': np.float32(ebv)*u.mag,
        #           'vdisp': np.float32(vdisp)*u.kilometer/u.second}

        return results, continuum
    
    def fnnls_continuum_plot(self, galwave, galflux, galivar, galphot,
                             continuum, continuum_fullwave, fullwave, objinfo,
                             png=None):
        """QA of the best-fitting continuum.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]

        ymin, ymax = 1e6, -1e6
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
        for ii in [0, 1, 2]: # iterate over cameras
            galsigma = 1 / np.sqrt(galivar[ii])
            ax1.fill_between(galwave[ii], galflux[ii]-galsigma, galflux[ii]+galsigma,
                            color=col1[ii])
            ax1.plot(galwave[ii], continuum[ii], color=col2[ii], alpha=1.0)#, color='k')

            #ax2.fill_between(galwave[ii], galflux[ii]-galsigma, galflux[ii]+galsigma,
            #                color=col1[ii])
            #ax2.plot(galwave[ii], continuum[ii], color=col2[ii], alpha=1.0)#, color='k')

            # get the robust range
            filtflux = median_filter(galflux[ii], 5)
            if np.min(filtflux) < ymin:
                ymin = np.min(filtflux)
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux)

        ax1.text(0.95, 0.92, '{}'.format(objinfo['targetid']), 
                 ha='right', va='center', transform=ax1.transAxes, fontsize=18)
        ax1.text(0.95, 0.86, r'{} {}'.format(objinfo['zredrock'], objinfo['znyxgalaxy']),
                 ha='right', va='center', transform=ax1.transAxes, fontsize=18)
        ax1.text(0.95, 0.80, r'{} {}'.format(objinfo['chi2'], objinfo['vdisp']),
                 ha='right', va='center', transform=ax1.transAxes, fontsize=18)
                    
        ax1.set_xlim(3500, 10000)
        ax1.set_ylim(ymin, ymax)
        ax1.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        ax1.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 

        # add the photometry

        if False:
            for ii in [0, 1, 2]: # iterate over cameras
                #galsigma = 1 / np.sqrt(galivar[ii])
                factor = 1e-17  * galwave[ii]**2 / (C_LIGHT * 1e13) # [10-17 erg/s/cm2/A --> maggies]
                good = np.where(galflux[ii] > 0)[0]
                if len(good) > 0:
                    ax2.plot(galwave[ii][good]/1e4, -2.5*np.log10(galflux[ii][good]*factor[good])-48.6, color=col1[ii])
                    #ax1.fill_between(galwave[ii]/1e4, -2.5*np.log10((galflux[ii]-galsigma) * factor,
                    #                 (galflux[ii]+galsigma) * factor, color=col1[ii])
                #ax2.plot(galwave[ii]/1e4, -2.5*np.log10(continuum[ii]*factor)-48.6, color=col2[ii], alpha=1.0)#, color='k')
                ax2.plot(galwave[ii]/1e4, -2.5*np.log10(continuum[ii]*factor)-48.6, color=col2[ii], alpha=1.0)#, color='k')

        factor = 10**(0.4 * 48.6) * fullwave**2 / (C_LIGHT * 1e13) # [erg/s/cm2/A --> maggies]
        ax2.plot(fullwave/1e4, -2.5*np.log10(continuum_fullwave*factor), color='gray', alpha=0.8)

        #ax2.errorbar(filtwave, flam, yerr=flamsigma, fmt='s')
        #ax2.scatter(galphot['lambda_eff']/1e4, galphot['abmag'], marker='s',
        #            color='blue', s=200, edgecolor='k')
        ax2.errorbar(galphot['lambda_eff']/1e4, galphot['abmag'], yerr=1/np.sqrt(galphot['abmag_ivar']),
                     fmt='s', markersize=15, markeredgewidth=3, markeredgecolor='k', markerfacecolor='red',
                     elinewidth=3, ecolor='blue', capsize=3)

        dm = 0.75
        ymin, ymax = galphot['abmag'].max() + dm, galphot['abmag'].min() - dm

        ax2.set_xlabel(r'Observed-frame Wavelength ($\mu$m)') 
        ax2.set_ylabel(r'AB Mag') 
        ax2.set_xlim(0.2, 6)
        ax2.set_ylim(ymin, ymax)

        ax2.set_xscale('log')
        ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax2.set_xticks([0.3, 0.5, 0.7, 1.0, 1.5, 3.0, 5.0])

        plt.subplots_adjust(bottom=0.1, right=0.95, top=0.95, wspace=0.17)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
        
    def _get_uncertainties(self, pcov=None, jac=None, return_covariance=False):
        """Determine the uncertainties from the diagonal terms of the covariance
        matrix. If the covariance matrix is not known, estimate it from the
        Jacobian following
        https://github.com/scipy/scipy/blob/master/scipy/optimize/minpack.py#L805

        """
        if pcov is None:
            from scipy.linalg import svd
            _, s, VT = svd(jac, full_matrices=False)
            threshold = np.finfo(float).eps * max(jac.shape) * s[0]
            s = s[s > threshold]
            VT = VT[:s.size]
            pcov = np.dot(VT.T / s**2, VT)

        unc = np.diag(pcov)**0.5

        if return_covariance:
            return unc, pcov
        else:
            return unc

    def _get_attenuation(self, wave, ebv):
        """Return the dust attenuation A(lambda)=E(B-V)*k(lambda)

        """
        klambda = 5.9 * (wave / 5500)**(-self.dustslope)
        return 10**(-0.4 * ebv * klambda)
        
    def dusty_continuum(self, ebv_and_coeffs, wave, sspflux):
        """Continuum model with dust attenuation."""
        ebv, coeffs = ebv_and_coeffs[0], ebv_and_coeffs[1:]
        atten = self._get_attenuation(wave, ebv)
        modelflux = (sspflux * atten[:, None]).dot(coeffs) 
        return modelflux
    
    def _dusty_continuum_resid(self, ebv_and_coeffs, wave, flux, isigma, sspflux):
        """Wrapper cost function for scipy.optimize.least_squares."""
        modelflux = self.dusty_continuum(ebv_and_coeffs, wave, sspflux)
        return isigma * (modelflux - flux)
    
    def fit_continuum(self):
        """More detailed continuum model fit, including dust attenuation.

        https://stackoverflow.com/questions/60335524/how-to-use-scipy-optimize-curve-fit-to-use-lists-of-variable/60409667#60409667 
        
        """
        from scipy.optimize import least_squares

        sspflux = np.concatenate(self.redshift_smooth_and_resample(), axis=0)

        wave = np.hstack(self.galwave)
        flux = np.hstack(self.galflux)
        isigma = 1 / np.sqrt(np.hstack(self.galivar))

        _ = self.fnnls_continuum()
        init_params = np.hstack([0.05, self.fnnls_coeffs])
        log.info(init_params)

        params = least_squares(self._dusty_continuum_resid, x0=init_params, #kwargs=sspflux,
                               bounds=(0.0, np.inf), args=(wave, flux, isigma, sspflux),
                               method='trf', tr_options={'regularize': True})
        #params, cov = curve_fit(self._dusty_continuum, wave, flux, sigma=1/np.sqrt(ivar),
        #                        kwargs=sspflux)
        #continuum_fit = fitter(ContinuumModel, wave, flux, 
        unc = self._get_uncertainties(jac=params.jac, return_covariance=False)
        
        log.info(params.x[0], self.ssp.info['age'][params.x[1:] > 1] / 1e9)
        pdb.set_trace()

def read_nyxgalaxy_lines():
    """Read the set of emission lines of interest.

    ToDo: add lines to mask during continuum-fitting but which we do not want to
    emission-line fit.

    """
    from pkg_resources import resource_filename
    
    linefile = resource_filename('desigal', 'data/nyxgalaxy_lines.ecsv')    
    linetable = Table.read(linefile, format='ascii.ecsv', guess=False)
    
    return linetable    

class EMLineModel(Fittable1DModel):
    """Class to model the emission-line spectra.

    """
    from astropy.modeling import Parameter

    # NB! The order of the parameters here matters! linevshift is the shift of
    # the emission-line velocity (in km/s) with respect to the systemic redshift
    linevshift_forbidden = Parameter(name='linevshift_forbidden', default=0.0, bounds=(-500.0, 500.0)) # [km/s]
    linevshift_balmer = Parameter(name='linevshift_balmer', default=0.0, bounds=(-500.0, 500.0)) # [km/s]

    linesigma_forbidden = Parameter(name='linesigma_forbidden', default=50.0, bounds=(5, 350)) # line-sigma [km/s]
    linesigma_balmer = Parameter(name='linesigma_balmer', default=50.0, bounds=(5, 350)) # line-sigma [km/s]

    # Fragile because the lines are hard-coded--
    oii_3726_amp = Parameter(name='oii_3726_amp', default=1.0)
    oii_3729_amp = Parameter(name='oii_3729_amp', default=1.0)
    oiii_4959_amp = Parameter(name='oiii_4959_amp', default=1.0)
    oiii_5007_amp = Parameter(name='oiii_5007_amp', default=3.0)
    nii_6548_amp = Parameter(name='nii_6548_amp', default=1.0)
    nii_6584_amp = Parameter(name='nii_6584_amp', default=3.0)
    sii_6716_amp = Parameter(name='sii_6716_amp', default=1.0)
    sii_6731_amp = Parameter(name='sii_6731_amp', default=1.0)
    hepsilon_amp = Parameter(name='hepsilon_amp', default=0.5)
    hdelta_amp = Parameter(name='hdelta_amp', default=0.5)
    hgamma_amp = Parameter(name='hgamma_amp', default=0.5)
    hbeta_amp = Parameter(name='hbeta_amp', default=1.0)
    halpha_amp = Parameter(name='halpha_amp', default=3.0)

    # tie the velocity shifts and line-widths together
    def tie_vshift(model):
        return model.linevshift_balmer
    linevshift_forbidden.tied = tie_vshift

    def tie_sigma(model):
        return model.linesigma_balmer
    linesigma_forbidden.tied = tie_sigma

    # tie the [NII] and [OIII] line-strengths together
    def tie_oiii(model):
        return model.oiii_5007_amp / 2.8875
    oiii_4959_amp.tied = tie_oiii

    def tie_nii(model):
        return model.nii_6584_amp / 2.936
    nii_6548_amp.tied = tie_nii
    
    def __init__(self,
                 linevshift_forbidden=linevshift_forbidden.default,
                 linevshift_balmer=linevshift_balmer.default,
                 linesigma_forbidden=linesigma_forbidden.default,
                 linesigma_balmer=linesigma_balmer.default,
                 oii_3726_amp=oii_3726_amp.default, 
                 oii_3729_amp=oii_3729_amp.default, 
                 oiii_4959_amp=oiii_4959_amp.default, 
                 oiii_5007_amp=oiii_5007_amp.default, 
                 nii_6548_amp=nii_6548_amp.default, 
                 nii_6584_amp=nii_6584_amp.default, 
                 sii_6716_amp=sii_6716_amp.default, 
                 sii_6731_amp=sii_6731_amp.default, 
                 hepsilon_amp=hepsilon_amp.default, 
                 hdelta_amp=hdelta_amp.default, 
                 hgamma_amp=hgamma_amp.default, 
                 hbeta_amp=hbeta_amp.default, 
                 halpha_amp=halpha_amp.default,
                 redshift=None,
                 emlineR=None, npixpercamera=None,
                 log10wave=None, **kwargs):
        """Initialize the emission-line model.
        
        emlineR -
        redshift - required keyword
        
        """
        self.linetable = read_nyxgalaxy_lines()
        self.redshift = redshift
        self.emlineR = emlineR
        self.npixpercamera = np.hstack([0, npixpercamera])

        # internal wavelength vector for building the emission-line model
        if log10wave is None:
            pixkms = 10.0
            dlogwave = pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
            log10wave = np.arange(np.log10(3000), np.log10(1e4), dlogwave)
        self.log10wave = log10wave
            
        super(EMLineModel, self).__init__(
            linevshift_forbidden=linevshift_forbidden,
            linevshift_balmer=linevshift_balmer,
            linesigma_forbidden=linesigma_forbidden,
            linesigma_balmer=linesigma_balmer,
            oii_3726_amp=oii_3726_amp,
            oii_3729_amp=oii_3729_amp,
            oiii_4959_amp=oiii_4959_amp,
            oiii_5007_amp=oiii_5007_amp,
            nii_6548_amp=nii_6548_amp,
            nii_6584_amp=nii_6584_amp,
            sii_6716_amp=sii_6716_amp,
            sii_6731_amp=sii_6731_amp,
            hepsilon_amp=hepsilon_amp,
            hdelta_amp=hdelta_amp,
            hgamma_amp=hgamma_amp,
            hbeta_amp=hbeta_amp,
            halpha_amp=halpha_amp, **kwargs)

    def evaluate(self, emlinewave, *args):
        """Evaluate the emission-line model.
        emlineR=None, npixpercamera=None, 

        """ 
        linevshift_forbidden, linevshift_balmer = args[0], args[1]
        linez_forbidden = self.redshift + linevshift_forbidden / C_LIGHT
        linez_balmer = self.redshift + linevshift_balmer / C_LIGHT

        linesigma_forbidden, linesigma_balmer = args[2], args[3]
        log10sigma_forbidden = linesigma_forbidden / C_LIGHT / np.log(10) # line-width [log-10 Angstrom]
        log10sigma_balmer = linesigma_balmer / C_LIGHT / np.log(10)       # line-width [log-10 Angstrom]

        lineamps = args[4:]
        linenames = self.linetable['name'].data
        isbalmers = self.linetable['isbalmer'].data

        # build the emission-line model [erg/s/cm2/A, observed frame]; should we
        # multiprocess this step?
        log10model = np.zeros_like(self.log10wave)
        for linename, lineamp, isbalmer in zip(linenames, lineamps, isbalmers):
            restlinewave = self.linetable[self.linetable['name'] == linename]['restwave'][0]
            if isbalmer:
                log10sigma = log10sigma_balmer
                linezwave = np.log10(restlinewave * (1.0 + linez_balmer)) # redshifted wavelength [log-10 Angstrom]
            else:
                log10sigma = log10sigma_forbidden
                linezwave = np.log10(restlinewave * (1.0 + linez_forbidden)) # redshifted wavelength [log-10 Angstrom]

            ww = np.abs(self.log10wave - linezwave) < 20 * log10sigma
            if np.count_nonzero(ww) > 0:
                #log.info(linename, 10**linezwave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
                log10model[ww] += lineamp * np.exp(-0.5 * (self.log10wave[ww]-linezwave)**2 / log10sigma**2)

        # split into cameras, resample, and convolve with the instrumental
        # resolution
        emlinemodel = []
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(self.npixpercamera[:ii+1])
            jpix = np.sum(self.npixpercamera[:ii+2])
            _emlinemodel = resample_flux(emlinewave[ipix:jpix], 10**self.log10wave, log10model)
            
            if self.emlineR is not None:
                _emlinemomdel = self.emlineR[ii].dot(_emlinemodel)
            
            #plt.plot(10**_emlinewave, _emlinemodel)
            #plt.plot(10**_emlinewave, self.emlineR[ii].dot(_emlinemodel))
            #plt.xlim(3870, 3920) ; plt.show()
            #pdb.set_trace()
            emlinemodel.append(_emlinemodel)
            
        return np.hstack(emlinemodel)

class EMLineFit(object):
    """Class to fit an emission-line spectrum.

    * https://docs.astropy.org/en/stable/modeling/example-fitting-constraints.html#tied
    * https://docs.astropy.org/en/stable/modeling/new-model.html
    * https://docs.astropy.org/en/stable/modeling/compound-models.html#parameters

    """
    def __init__(self, nball=10, chi2fail=1e6):
        """Write me.
        
        """
        from astropy.modeling import fitting

        self.nball = nball
        self.chi2fail = chi2fail
        self.pixkms = 10.0 # pixel size for internal wavelength array [km/s]

        self.fitter = fitting.LevMarLSQFitter()
                
    def chi2(self, bestfit, emlinewave, emlineflux, emlineivar):
        """Compute the reduced chi^2."""
        dof = len(emlinewave) - len(bestfit.parameters)
        emlinemodel = bestfit(emlinewave)
        chi2 = np.sum(emlineivar * (emlineflux - emlinemodel)**2) / dof
        return chi2

    def emlinemodel_bestfit(self, galwave, galres, nyxgalaxy_table):
        """Wrapper function to get the best-fitting emission-line model
        from an nyxgalaxy table (to be used to build QA).

        """
        npixpercamera = [len(gw) for gw in galwave]

        redshift = nyxgalaxy_table['Z']
        linesigma_forbidden = nyxgalaxy_table['LINESIGMA_FORBIDDEN']
        linesigma_balmer = nyxgalaxy_table['LINESIGMA_BALMER']

        linevshift_forbidden = nyxgalaxy_table['LINEVSHIFT_FORBIDDEN']
        linevshift_balmer = nyxgalaxy_table['LINEVSHIFT_BALMER']
        #linez_forbidden = nyxgalaxy_table['LINEZ_FORBIDDEN']
        #linez_balmer = nyxgalaxy_table['LINEZ_BALMER']
        #linevshift_forbidden = (linez_forbidden - redshift) * C_LIGHT # [km/s]
        #linevshift_balmer = (linez_balmer - redshift) * C_LIGHT # [km/s]

        EMLine = EMLineModel(linevshift_forbidden=linevshift_forbidden,
                             linevshift_balmer=linevshift_balmer,
                             linesigma_forbidden=linesigma_forbidden,
                             linesigma_balmer=linesigma_balmer,
                             redshift=redshift, emlineR=galres,
                             npixpercamera=npixpercamera)
        # skip linevshift_[forbidden,balmer] and linesigma_[forbidden,balmer]
        lineargs = [nyxgalaxy_table[linename.upper()] for linename in EMLine.param_names[4:]] 
        lineargs = [linevshift_forbidden, linevshift_balmer, linesigma_forbidden, linesigma_balmer] + lineargs

        _emlinemodel = EMLine.evaluate(np.hstack(galwave), *lineargs)

        # unpack it
        emlinemodel = []
        npix = np.hstack([0, npixpercamera])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            emlinemodel.append(_emlinemodel[ipix:jpix])

        return emlinemodel
    
    def fit(self, galwave, galflux, galivar, galres, continuum,
            redshift, verbose=False):
        """Perform the fit minimization / chi2 minimization.
        
        EMLineModel object
        FC - ContinuumFit object

        need to take into account the instrumental velocity width when computing integrated fluxes
        
        """
        #from scipy import integrate
        from astropy.stats import sigma_clipped_stats
        
        npixpercamera = [len(gw) for gw in galwave]

        # we have to stack the per-camera spectra for LevMarLSQFitter
        _galflux = np.hstack(galflux)
        emlinewave = np.hstack(galwave)
        emlineivar = np.hstack(galivar)
        emlineflux = _galflux - np.hstack(continuum)

        dlogwave = self.pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        log10wave = np.arange(np.log10(emlinewave.min()), np.log10(emlinewave.max()), dlogwave)
        
        self.EMLineModel = EMLineModel(redshift=redshift, emlineR=galres,
                                       npixpercamera=npixpercamera,
                                       log10wave=log10wave)

        #weights = np.ones_like
        weights = 1 / np.sqrt(emlineivar)
        bestfit = self.fitter(self.EMLineModel, emlinewave, 
                              emlineflux, weights=weights)

        chi2 = self.chi2(bestfit, emlinewave, emlineflux, emlineivar).astype('f4')

        nparam = len(self.EMLineModel.parameters)
        
        # Pack the results in a dictionary and return.
        # https://gist.github.com/eteq/1f3f0cec9e4f27536d52cd59054c55f2
        result = {
            'converged': False,
            'fit_message': self.fitter.fit_info['message'],
            'nparam': nparam,
            'npix': len(emlinewave),
            'dof': len(emlinewave) - len(self.EMLineModel.parameters),
            'chi2': chi2,
            'linenames': [ll.replace('_amp', '') for ll in self.EMLineModel.param_names[4:]],
        }
        #for param in bestfit.param_names:
        #    result.update({param: getattr(bestfit, param).value})
        
        # uncertainties
        if self.fitter.fit_info['param_cov'] is not None:
            cov = self.fitter.fit_info['param_cov']
            ivar = 1 / np.diag(cov)
            result['converged'] = True
        else:
            cov = np.zeros((nparam, nparam))
            ivar = np.zeros(nparam)

        # Need to be careful about uncertainties for tied parameters--
        # https://github.com/astropy/astropy/issues/7202
        
        #err_params = np.sqrt(np.diag(fitter.fit_info['param_cov']))
        #err = model.copy()
        #fitting._fitter_to_model_params(err, err_params)            
            
        count = 0
        for ii, pp in enumerate(bestfit.param_names):
            pinfo = getattr(bestfit, pp)
            iinfo = getattr(self.EMLineModel, pp)

            # need to think about this more deeply
            #if pinfo.value == iinfo.value: # not fitted
            #    result.update({pinfo.name: np.float(0.0)})
            #else:
            #    result.update({pinfo.name: pinfo.value.astype('f4')})
            result.update({pinfo.name: pinfo.value.astype('f4')})
                
            if pinfo.fixed:
                result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
            elif pinfo.tied:
                # hack! see https://github.com/astropy/astropy/issues/7202
                result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
            else:
                result.update({'{}_ivar'.format(pinfo.name): ivar[count].astype('f4')})
                count += 1

            #if 'forbidden' in pinfo.name:
            #    pdb.set_trace()

        # hack for tied parameters---gotta be a better way to do this
        #result['oiii_4959_amp_ivar'] = result['oiii_5007_amp_ivar'] * 2.8875**2
        result['oiii_4959_amp_ivar'] = result['oiii_5007_amp_ivar'] * 2.8875**2
        result['nii_6548_amp_ivar'] = result['nii_6548_amp_ivar'] * 2.936**2
        result['linevshift_forbidden_ivar'] = result['linevshift_balmer_ivar']
        result['linesigma_forbidden_ivar'] = result['linesigma_balmer_ivar']

        ## convert the vshifts to redshifts
        #result['linez_forbidden'] = redshift + result['linevshift_forbidden'] / C_LIGHT
        #result['linez_balmer'] = redshift + result['linevshift_balmer'] / C_LIGHT
        #result['linez_forbidden_ivar'] = result['linevshift_forbidden_ivar'] * C_LIGHT**2
        #result['linez_balmer_ivar'] = result['linevshift_balmer_ivar'] * C_LIGHT**2

        # now loop back through and if ivar==0 then set the parameter value to zero
        if self.fitter.fit_info['param_cov'] is not None:
            for pp in bestfit.param_names:
                if result['{}_ivar'.format(pp)] == 0.0:
                    result[pp] = 0.0

        # get continuum fluxes, EWs, and upper limits
        emlinemodel = bestfit(emlinewave)
        galflux_nolines = _galflux - emlinemodel

        # measure the 4000-Angstrom break from the data and the model
        restwave = emlinewave / (1 + redshift) # [Angstrom]

        restflam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5)
        restflux_nolines_nu = galflux_nolines * restflam2fnu   # rest-frame, [erg/s/cm2/Hz]
        restcontinuum_nu = np.hstack(continuum) * restflam2fnu

        good = emlineivar > 0
        restemlinevar_nu = np.zeros_like(restcontinuum_nu)
        restemlinevar_nu[good] = restflam2fnu[good]**2 / emlineivar[good] 

        indxblu = np.where((restwave >= 3850) * (restwave <= 3950))[0]
        indxred = np.where((restwave >= 4000) * (restwave <= 4100))[0]
        if len(indxblu) > 5 and len(indxred) > 5:
            blufactor, redfactor = 3950.0 - 3850.0, 4100.0 - 4000.0
            deltawave = np.gradient(restwave)

            numer = blufactor * np.sum(deltawave[indxred] * restflux_nolines_nu[indxred])
            denom = redfactor * np.sum(deltawave[indxblu] * restflux_nolines_nu[indxblu])
            numer_var = blufactor**2 * np.sum(deltawave[indxred] * restemlinevar_nu[indxred])
            denom_var = redfactor**2 * np.sum(deltawave[indxblu] * restemlinevar_nu[indxblu])

            d4000 =  numer / denom
            d4000_var = (numer_var + numer**2 * denom_var) / denom**2
            d4000_ivar = 1 / d4000_var
              
            d4000_model = (blufactor / redfactor) * np.sum(deltawave[indxred] * restcontinuum_nu[indxred]) / \
              np.sum(deltawave[indxblu] * restcontinuum_nu[indxblu])
        else:
            d4000, d4000_ivar, d4000_model = 0.0, 0.0, 0.0
 
        result['d4000'] = d4000
        result['d4000_ivar'] = d4000_ivar
        result['d4000_model'] = d4000_model
              
        sigma_cont = 150.0
        for oneline in self.EMLineModel.linetable:
            line = oneline['name']

            # get the emission-line flux
            zwave = oneline['restwave'] * (1 + redshift)
            if oneline['isbalmer']:
                linesigma = result['linesigma_balmer']
            else:
                linesigma = result['linesigma_forbidden']

            linesigma_ang = zwave * linesigma / C_LIGHT # [observed-frame Angstrom]
            norm = np.sqrt(2.0 * np.pi) * linesigma_ang

            result['{}_flux'.format(line)] = result['{}_amp'.format(line)] * norm
            result['{}_flux_ivar'.format(line)] = result['{}_amp_ivar'.format(line)] / norm**2

            # boxcar integration, chi2, and number of pixels
            lineindx = np.where((emlinewave > (zwave - 3*sigma_cont * zwave / C_LIGHT)) *
                                (emlinewave < (zwave + 3.*sigma_cont * zwave / C_LIGHT)) *
                                (emlineivar > 0))[0]
            npix = len(lineindx)
            if npix > 0:
                dof = npix - 3 # ??? [redshift, sigma, and amplitude]
                chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-emlinemodel[lineindx])**2) / dof
                boxflux = np.sum(emlineflux[lineindx])
                boxflux_ivar = 1 / np.sum(1 / emlineivar[lineindx])
            else:
                npix, chi2, boxflux, boxflux_ivar = 0.0, 1e6, 0.0, 0.0

            result['{}_npix'.format(line)] = npix
            result['{}_chi2'.format(line)] = chi2
            result['{}_boxflux'.format(line)] = boxflux
            result['{}_boxflux_ivar'.format(line)] = boxflux_ivar
            
            # get the continuum and EWs
            indxlo = np.where((emlinewave > (zwave - 10*sigma_cont * zwave / C_LIGHT)) *
                              (emlinewave < (zwave - 3.*sigma_cont * zwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indxhi = np.where((emlinewave < (zwave + 10*sigma_cont * zwave / C_LIGHT)) *
                              (emlinewave > (zwave + 3.*sigma_cont * zwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indx = np.hstack((indxlo, indxhi))

            if len(indx) > 5: # require at least 5 pixels to get the continuum level
                _, cmed, csig = sigma_clipped_stats(galflux_nolines[indx], sigma=3.0)
                civar = (np.sqrt(len(indx)) / csig)**2
            else:
                cmed, civar = 0.0, 0.0
                
            result['{}_cont'.format(line)] = cmed
            result['{}_cont_ivar'.format(line)] = civar

            if result['{}_cont_ivar'.format(line)] != 0.0:
                factor = (1 + redshift) / result['{}_cont'.format(line)]
                ew = result['{}_flux'.format(line)] * factor   # rest frame [A]
                ewivar = result['{}_flux_ivar'.format(line)] / factor**2

                # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar)
                ewlimit = fluxlimit * factor
            else:
                ew, ewivar, fluxlimit, ewlimit = 0.0, 0.0, 0.0, 0.0

            result['{}_ew'.format(line)] = ew
            result['{}_ew_ivar'.format(line)] = ewivar
            result['{}_flux_limit'.format(line)] = fluxlimit
            result['{}_ew_limit'.format(line)] = ewlimit

            # simple QA
            if 'alpha' in line and False:
                import matplotlib.pyplot as plt
                _indx = np.where((emlinewave > (zwave - 15*sigma_cont * zwave / C_LIGHT)) *
                                (emlinewave < (zwave + 15*sigma_cont * zwave / C_LIGHT)))[0]
                plt.plot(emlinewave[_indx], emlineflux[_indx])
                plt.plot(emlinewave[_indx], galflux_nolines[_indx])
                plt.scatter(emlinewave[indx], galflux_nolines[indx], color='red')
                plt.axhline(y=cmed, color='k')
                plt.axhline(y=cmed+csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.axhline(y=cmed-csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.savefig('junk.png')
                #pdb.set_trace()
            
        return result, emlinemodel
    
    def emlineplot(self, galwave, galflux, galivar, continuum,
                   _emlinemodel, redshift, objinfo, png=None):
        """Plot the emission-line spectrum and best-fitting model.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]
        
        #fig, ax = plt.subplots(1, 4, figsize=(16, 10))#, sharey=True)
        fig = plt.figure(figsize=(16, 16))
        gs = fig.add_gridspec(3, 4, height_ratios=[4, 2, 2])
        #gs = fig.add_gridspec(2, 4, gridspec_kw={'width_ratios': 1.0, 'height_ratios': 0.5})

        # full spectrum
        bigax = fig.add_subplot(gs[0, :])

        ymin, ymax = 1e6, -1e6
        for ii in [0, 1, 2]: # iterate over cameras
            emlinewave = galwave[ii]
            emlineflux = galflux[ii] - continuum[ii]
            emlinemodel = _emlinemodel[ii]
            emlinesigma = np.zeros_like(emlinewave)
            good = galivar[ii] > 0
            emlinesigma[good] = 1 / np.sqrt(galivar[ii][good])
            
            bigax.fill_between(emlinewave, emlineflux-emlinesigma, emlineflux+emlinesigma,
                               color=col1[ii], alpha=0.7)
            bigax.plot(emlinewave, emlinemodel, color=col2[ii], lw=2)

            # get the robust range
            sigflux, filtflux = np.std(emlineflux), median_filter(emlineflux, 5)

            if -3 * sigflux < ymin:
                ymin = -2 * sigflux
            #if np.min(filtflux) < ymin:
            #    ymin = np.min(filtflux)
            if np.min(emlinemodel) < ymin:
                ymin = 0.8 * np.min(emlinemodel)
                
            if 5 * sigflux > ymax:
                ymax = 4 * sigflux
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux)
            if np.max(emlinemodel) > ymax:
                ymax = np.max(emlinemodel) * 1.2

        bigax.text(0.95, 0.92, '{}'.format(objinfo['targetid']), 
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
        bigax.text(0.95, 0.86, r'{}'.format(objinfo['zredrock']),
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
        bigax.text(0.95, 0.80, r'{} {}'.format(objinfo['linevshift_balmer'], objinfo['linevshift_forbidden']),
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
        bigax.text(0.95, 0.74, r'{} {}'.format(objinfo['linesigma_balmer'], objinfo['linesigma_forbidden']),
                   ha='right', va='center', transform=bigax.transAxes, fontsize=18)
                
        bigax.set_xlim(3500, 10000)
        bigax.set_ylim(ymin, ymax)
        
        #bigax.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 
        
        # zoom in on individual emission lines - use linetable!
        sig = 500.0 # [km/s]

        meanwaves = [np.mean([3730,3727]), 3971, 4103, 4342, 4863, np.mean([4960, 5008]), 6565, np.mean([6718.294, 6732.673])]
        deltawaves = [0, 0, 0, 0, 0, (5007-4959)/2, (6585-6550)/2, (6733-6718)/2]
        linenames = [r'[OII]', r'H$\epsilon$', r'H$\delta$', r'H$\gamma$', r'H$\beta$', r'[OIII]', r'H$\alpha$+[NII]', r'[SII]']
        nline = len(meanwaves)

        removelabels = np.ones(nline, bool)
        ymin, ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
        
        ax = []
        for iax, (meanwave, deltawave, linename) in enumerate(zip(meanwaves, deltawaves, linenames)):
            if iax < 4:
                xx = fig.add_subplot(gs[1, iax])
            else:
                xx = fig.add_subplot(gs[2, iax-4])
            ax.append(xx)

            wmin = (meanwave-deltawave)*(1+redshift)-2.5*sig*meanwave/3e5
            wmax = (meanwave+deltawave)*(1+redshift)+2.5*sig*meanwave/3e5

            # iterate over cameras
            for ii in [0, 1, 2]: # iterate over cameras
                emlinewave = galwave[ii]
                emlineflux = galflux[ii] - continuum[ii]
                emlinemodel = _emlinemodel[ii]
                emlinesigma = np.zeros_like(emlinewave)
                good = galivar[ii] > 0
                emlinesigma[good] = 1 / np.sqrt(galivar[ii][good])
            
                indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
                #log.info(ii, linename, len(indx))
                if len(indx) > 1:
                    removelabels[iax] = False
                    xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                                    emlineflux[indx]+emlinesigma[indx], color=col1[ii], alpha=0.5)
                    xx.plot(emlinewave[indx], emlinemodel[indx], color=col2[ii], lw=3)

                    # get the robust range
                    sigflux, filtflux = np.std(emlineflux[indx]), median_filter(emlineflux[indx], 3)

                    _ymin, _ymax = -1.5 * sigflux, 3 * sigflux
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
                    
                xx.text(0.08, 0.9, linename, ha='left', va='center',
                        transform=xx.transAxes, fontsize=20)
                    
        for iax, xx in enumerate(ax):
            if removelabels[iax]:
                xx.set_ylim(0, 1)
                xx.set_xticklabels([])
                xx.set_yticklabels([])
            else:
                xx.set_ylim(ymin[iax], ymax[iax])
                xlim = xx.get_xlim()
                #log.info(linenames[iax], xlim, np.diff(xlim))
                xx.xaxis.set_major_locator(ticker.MaxNLocator(3))
                #xx.xaxis.set_major_locator(ticker.MultipleLocator(20)) # wavelength spacing of ticks [Angstrom]
                #if iax == 2:
                #    pdb.set_trace()

        # common axis labels
        tp, bt, lf, rt = 0.95, 0.09, 0.12, 0.95
        
        fig.text(lf-0.07, (tp-bt)/2+bt,
                 r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
                 ha='center', va='center', rotation='vertical')
        fig.text((rt-lf)/2+lf, bt-0.06, r'Observed-frame Wavelength ($\AA$)',
                 ha='center', va='center')
            
        plt.subplots_adjust(wspace=0.27, top=tp, bottom=bt, left=lf, right=rt, hspace=0.22)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
