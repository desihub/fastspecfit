"""
desigal.igalfit
===============

"""
import os, pdb
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from astropy.modeling import Fittable1DModel
from desispec.interpolation import resample_flux

from scipy import constants
C_LIGHT = constants.c / 1000.0 # [km/s]

def init_galfit(tile, night, zbest, CFit):
    """Initialize the output data table.

    CFit - ContinuumFit class

    """
    import astropy.units as u
    nobj = len(zbest)

    # Grab info on the emission lines and the continuum.
    linetable = read_igalfit_lines()
    nssp_coeff = len(CFit.sspinfo)

    out = Table()
    for zbestcol in ['TARGETID', 'Z']:#, 'ZERR']:#, 'SPECTYPE', 'DELTACHI2']
        out[zbestcol] = zbest[zbestcol]
    out.add_column(Column(name='NIGHT', data=np.repeat(night, nobj)), index=0)
    out.add_column(Column(name='TILE', data=np.repeat(tile, nobj)), index=0)
    out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f8'))
    out.add_column(Column(name='CONTINUUM_CHI2', length=nobj, dtype='f4')) # reduced chi2
    out.add_column(Column(name='CONTINUUM_DOF', length=nobj, dtype=np.int32))
    out.add_column(Column(name='CONTINUUM_Z', length=nobj, dtype='f8'))
    out.add_column(Column(name='CONTINUUM_AGE', length=nobj, dtype='f4', unit=u.Gyr))
    out.add_column(Column(name='CONTINUUM_EBV', length=nobj, dtype='f4', unit=u.mag))
    out.add_column(Column(name='CONTINUUM_VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
    out.add_column(Column(name='D4000', length=nobj, dtype='f4'))
    out.add_column(Column(name='D4000_IVAR', length=nobj, dtype='f4'))
    out.add_column(Column(name='D4000_MODEL', length=nobj, dtype='f4'))
    #out.add_column(Column(name='LINEZ', length=nobj, dtype='f4'))
    #out.add_column(Column(name='LINEZ_IVAR', length=nobj, dtype='f4'))
    #out.add_column(Column(name='LINESIGMA', length=nobj, dtype='f4', unit=u.kilometer/u.second))
    #out.add_column(Column(name='LINESIGMA_IVAR', length=nobj, dtype='f4', unit=u.second**2/u.kilometer**2))
    out.add_column(Column(name='LINEZ_FORBIDDEN', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEZ_FORBIDDEN_IVAR', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEZ_BALMER', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINEZ_BALMER_IVAR', length=nobj, dtype='f4'))
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
    
def get_data(tile='66003', night='20200315', overwrite=False, verbose=False):
    """Parse the reduced data for a given tile and night. Select the subset of
    objects with good redshifts and visual inspections.

    ELG - tile='70005', night='20200228'       
    20200315, 66003 - BGS+MWS 
    
    https://desi.lbl.gov/trac/wiki/SurveyValidation/TruthTables
    https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=5720;filename=DESI_data_042820.pdf

    """
    import fitsio
    from astropy.table import join, vstack
    from desispec.spectra import Spectra
    import desispec.io

    desigal_dir = os.getenv('DESIGAL_DATA')
    zbestoutfile = os.path.join(desigal_dir, 'zbest-{}-{}.fits'.format(tile, night))
    coaddoutfile = os.path.join(desigal_dir, 'coadd-{}-{}.fits'.format(tile, night))
    if os.path.isfile(coaddoutfile) and not overwrite:
        zbest = Table(fitsio.read(zbestoutfile))
        print('Read {} redshifts from {}'.format(len(zbest), zbestoutfile))

        coadd = desispec.io.read_spectra(coaddoutfile)
        print('Read {} spectra from {}'.format(len(zbest), coaddoutfile))

        return zbest, coadd

    print('Parsing data from tile, night {}, {}'.format(tile, night))

    truthdir = os.path.join(os.getenv('DESI_ROOT'), 'sv', 'vi', 'TruthTables')
    specprod_dir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', 'andes')
    datadir = os.path.join(specprod_dir, 'tiles', tile, night)
    
    if tile == '66003':
        truthfile = os.path.join(truthdir, 'truth_table_BGS_v1.2.csv')
    elif tile == '70500':
        truthfile = os.path.join(truthdir, 'truth_table_ELG_v1.2_latest.csv')
    else:
        raise ValueError('Fix me.')
    truth = Table.read(truthfile)
    best = np.where(
        (truth['best quality'] >= 2.5) * 
        (truth['Redrock spectype'] == 'GALAXY') *
        (truth['Redrock z'] < 0.75) 
    )[0]
    #goodobj = np.where((zb['DELTACHI2'] > 100) * (zb['ZWARN'] == 0) * (zb['SPECTYPE'] == 'GALAXY'))[0]
    
    print('Read {}/{} good redshifts from {}'.format(len(best), len(truth), truthfile))
    truth = truth[best]
    
    zbest = []
    spectra, keepindx = [], []
    for spectro in ('0'):
    #for spectro in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
        zbestfile = os.path.join(datadir, 'zbest-{}-{}-{}.fits'.format(spectro, tile, night))
        coaddfile = os.path.join(datadir, 'coadd-{}-{}-{}.fits'.format(spectro, tile, night))
        if os.path.isfile(zbestfile) and os.path.isfile(coaddfile):
            zb = Table(fitsio.read(zbestfile))
            keep = np.where(np.isin(zb['TARGETID'], truth['TARGETID']))[0]
            print('Spectrograph {}: N={}'.format(spectro, len(keep)))
            if len(keep) > 0:
                zbest.append(zb[keep])
                keepindx.append(keep)
                spectra.append(desispec.io.read_spectra(coaddfile))

    if len(zbest) == 0:
        raise ValueError('No spectra found for tile {} and night {}!'.format(tile, night))
        
    # combine the spectrographs
    print('Stacking all the spectra.')
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

    print('Writing {} redshifts to {}'.format(len(zbest), zbestoutfile))
    zbest.write(zbestoutfile, overwrite=True)

    print('Writing {} spectra to {}'.format(len(zbest), coaddoutfile))
    desispec.io.write_spectra(coaddoutfile, coadd)

    return zbest, coadd

def _unpack_spectrum(specobj, zbest, iobj):
    """Unpack a single spectrum into a list.

    """
    from desispec.resolution import Resolution
    
    galwave, galflux, galivar, galres = [], [], [], []
    for camera in ('b', 'r', 'z'):
        galwave.append(specobj.wave[camera])
        galflux.append(specobj.flux[camera][iobj, :])
        galivar.append(specobj.ivar[camera][iobj, :])
        galres.append(Resolution(specobj.resolution_data[camera][iobj, :, :]))

    zredrock = zbest['Z'][iobj]

    return galwave, galflux, galivar, galres, zredrock

def _smooth_and_resample(args):
    """Wrapper for multiprocessing."""
    return smooth_and_resample(*args)

def smooth_and_resample(sspflux, sspwave, galwave, galR):
    """
    sspflux[npix] - redshifted SSP template
    sspwave[npix] - redshifted SSP wavelength
    
    """
    return galR.dot(resample_flux(galwave, sspwave, sspflux, extrapolate=True))

class ContinuumFit():
    def __init__(self, metallicity='Z0.0190', minwave=0.0, maxwave=1e4, nproc=1,
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
    
        """
        import fitsio
        from astropy.cosmology import FlatLambdaCDM

        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)        

        self.metallicity = metallicity
        self.Z = float(metallicity[1:])
        self.library = 'CKC14z'
        self.isochrone = 'Padova' # would be nice to get MIST in here
        self.imf = 'Kroupa'

        self.nproc = nproc
        self.verbose = verbose

        # Don't hard-code the path!
        self.ssppath = os.getenv('DESIGAL_TEMPLATES')
        self.sspfile = os.path.join(self.ssppath, 'SSP_{}_{}_{}_{}.fits'.format(
            self.isochrone, self.library, self.imf, self.metallicity))

        if verbose:
            print('Reading {}'.format(self.sspfile))
        wave = fitsio.read(self.sspfile, ext='WAVE')
        flux = fitsio.read(self.sspfile, ext='FLUX')
        sspinfo = Table(fitsio.read(self.sspfile, ext='METADATA'))
        
        # Trim the wavelengths and subselect the number of ages/templates.
        keep = np.where((wave >= minwave) * (wave <= maxwave))[0]
        self.sspwave = wave[keep]
        self.sspflux = flux[keep, ::3]
        self.sspinfo = sspinfo[::3]

        self.nage = len(self.sspinfo['age'])
        self.npix = len(wave)

        # dust parameters and emission lines
        self.dustslope = 0.7
        
        self.linetable = read_igalfit_lines()
        
    def redshift_smooth_and_resample(self, galwave, galres, redshift):
        """Redshift, convolve with the spectral resolution, and 
        resample in wavelength.

        ToDo: velocity dispersion smoothing
        
        """
        import multiprocessing

        # loop over cameras then SSP ages
        smoothflux = []
        for icamera in [0, 1, 2]: # iterate on cameras
            args = [(self.sspflux[:, iage] / (1 + redshift), self.sspwave * (1 + redshift),
                     galwave[icamera], galres[icamera]) for iage in np.arange(self.nage)]
            with multiprocessing.Pool(self.nproc) as pool:
                smoothflux.append(np.array(pool.map(_smooth_and_resample, args)).T)

        #if plot:
        #    import matplotlib.pyplot as plt
        #    ww = 110
        #    plt.plot(self.ssp.wave * (1+self.zredrock), self.ssp.flux[:, ww], color='gray', alpha=0.5)
        #    for ii in [0, 1, 2]: # iterate over cameras
        #        plt.plot(self.galwave[ii], sspflux[ii][:, 110], color='k')
        #    plt.xlim(3900, 4000)
            
        return smoothflux # [npix, nage]

    def fnnls_continuum_bestfit(self, coeff, sspflux=None, galwave=None,
                                galres=None, redshift=None):
        if sspflux is None:
            sspflux = self.redshift_smooth_and_resample(galwave, galres, redshift)
            bestfit = [_sspflux.dot(coeff) for _sspflux in sspflux]
        else:
            bestfit = sspflux.dot(coeff)
            
        return bestfit

    def fnnls_continuum(self, galwave, galflux, galivar, galres, redshift, linetable,
                        sigma_mask=300.0):
        """Fit the continuum using fast NNLS.
        https://github.com/jvendrow/fnnls

        sigma_mask - emission-line masking sigma [km/s]

        ToDo: mask more emission lines than we fit (e.g., Mg II).
        
        """
        from time import time
        from fnnls import fnnls
        
        sspflux = self.redshift_smooth_and_resample(galwave, galres, redshift) 

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
        print('ToDo: add vdisp and ebv, and restrict maxage of templates.')
        t0 = time()
        coeff = fnnls(_sspflux * np.sqrt((_galivar * emlinemask)[:, None]), _galflux * np.sqrt(_galivar))[0]
        dt = time() - t0

        # ToDo: fit for dust, the redshift, and velocity dispersion
        zcontinuum, vdisp, ebv = redshift, 0.0, 0.0

        # Need to push the calculation of the best-fitting continuum to a
        # function so we can call it when building QA.
        _continuum = self.fnnls_continuum_bestfit(coeff, _sspflux)
        dof = np.sum(_galivar > 0) - self.nage
        chi2 = np.sum(_galivar * (_galflux - _continuum)**2) / dof

        # Compute the light-weighted age.
        weight = coeff[coeff > 0]
        age = np.sum(weight * self.sspinfo['age'][coeff > 0]) / np.sum(weight) / 1e9 # [Gyr]
        if self.verbose:
            print('Continuum fit done in {:.3f} seconds:'.format(dt))
            print('  Non-zero templates: {}'.format(len(weight)))
            print('  Reduced chi2: {:.3f}'.format(chi2))
            print('  Velocity dispersion: {:.3f} km/s'.format(vdisp))
            print('  Reddening: {:.3f} mag'.format(ebv))
            print('  Light-weighted age: {:.3f} Gyr'.format(age))

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
                   'vdisp': np.float32(vdisp), 'z': zcontinuum}
        #results = {'coeff': coeff, 'chi2': np.float32(chi2), 'dof': dof,
        #           'age': np.float32(age)*u.Gyr, 'ebv': np.float32(ebv)*u.mag,
        #           'vdisp': np.float32(vdisp)*u.kilometer/u.second}

        return results, continuum
    
    def fnnls_continuum_plot(self, galwave, galflux, galivar, continuum,
                             objinfo, png=None):
        """QA of the best-fitting continuum.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import seaborn as sns

        sns.set(context='talk', style='ticks', font_scale=1.3)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]

        ymin, ymax = 1e6, -1e6
        
        fig, ax = plt.subplots(figsize=(12, 8))
        for ii in [0, 1, 2]: # iterate over cameras
            galsigma = 1 / np.sqrt(galivar[ii])
            ax.fill_between(galwave[ii], galflux[ii]-galsigma, galflux[ii]+galsigma,
                            color=col1[ii])
            ax.plot(galwave[ii], continuum[ii], color=col2[ii], alpha=1.0)#, color='k')

            # get the robust range
            filtflux = median_filter(galflux[ii], 5)
            if np.min(filtflux) < ymin:
                ymin = np.min(filtflux)
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux)

        ax.text(0.95, 0.92, '{}'.format(objinfo['targetid']), 
                ha='right', va='center', transform=ax.transAxes, fontsize=18)
        ax.text(0.95, 0.86, r'{} {}'.format(objinfo['zredrock'], objinfo['zigalfit']),
                ha='right', va='center', transform=ax.transAxes, fontsize=18)
        ax.text(0.95, 0.80, r'{} {}'.format(objinfo['chi2'], objinfo['vdisp']),
                ha='right', va='center', transform=ax.transAxes, fontsize=18)
                    
        ax.set_xlim(3500, 10000)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        ax.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 

        plt.subplots_adjust(bottom=0.15, right=0.95, top=0.95)

        if png:
            print('Writing {}'.format(png))
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
        print(init_params)

        params = least_squares(self._dusty_continuum_resid, x0=init_params, #kwargs=sspflux,
                               bounds=(0.0, np.inf), args=(wave, flux, isigma, sspflux),
                               method='trf', tr_options={'regularize': True})
        #params, cov = curve_fit(self._dusty_continuum, wave, flux, sigma=1/np.sqrt(ivar),
        #                        kwargs=sspflux)
        #continuum_fit = fitter(ContinuumModel, wave, flux, 
        unc = self._get_uncertainties(jac=params.jac, return_covariance=False)
        
        print(params.x[0], self.ssp.info['age'][params.x[1:] > 1] / 1e9)
        pdb.set_trace()

def read_igalfit_lines():
    """Read the set of emission lines of interest.

    ToDo: add lines to mask during continuum-fitting but which we do not want to
    emission-line fit.

    """
    from pkg_resources import resource_filename
    
    linefile = resource_filename('desigal', 'data/igalfit_lines.ecsv')    
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
        self.linetable = read_igalfit_lines()
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
                #print(linename, 10**linezwave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
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

        #emlinemodel = []
        #for ii in [0, 1, 2]: # iterate over cameras
        #    ipix = np.sum(self.npixpercamera[:ii+1])
        #    jpix = np.sum(self.npixpercamera[:ii+2])
        #    _emlinewave = log10wave[ipix:jpix]
        #    _emlinemodel = np.zeros_like(_emlinewave)
        #    
        #    for linename, lineamp in zip(linenames, lineamps):
        #        restlinewave = self.linetable[self.linetable['name'] == linename]['restwave'][0]
        #        linezwave = np.log10(restlinewave * (1.0 + linez)) # redshifted wavelength [log-10 Angstrom]
        #        #lineflux = lineamp * (np.sqrt(2.0 * np.pi) * log10sigma)
        #        #lineamp = lineflux / (np.sqrt(2.0 * np.pi) * log10sigma)
        #    
        #        # Construct the spectrum [erg/s/cm2/A, rest]
        #        ww = np.abs(_emlinewave - linezwave) < 20 * log10sigma
        #        if np.count_nonzero(ww) > 0:
        #            #print(linename, 10**linezwave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
        #            _emlinemodel[ww] += lineamp * np.exp(-0.5 * (_emlinewave[ww]-linezwave)**2 / log10sigma**2)
        #
        #    # resample and (optionally) convolve with the spectral resolution
        #    if self.emlineR is not None:
        #        #pdb.set_trace()
        #        #here  resample_flux(galwave, sspwave, sspflux, extrapolate=True)
        #        _emlinemomdel = self.emlineR[ii].dot(_emlinemodel)
        #    
        #    #plt.plot(10**_emlinewave, _emlinemodel)
        #    #plt.plot(10**_emlinewave, self.emlineR[ii].dot(_emlinemodel))
        #    #plt.xlim(3870, 3920) ; plt.show()
        #    #pdb.set_trace()
        #    emlinemodel.append(_emlinemodel)

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

    def emlinemodel_bestfit(self, galwave, galres, galfit_table):
        """Wrapper function to get the best-fitting emission-line model
        from an igalfit table (to be used to build QA).

        """
        redshift = galfit_table['Z']
        linez_forbidden = galfit_table['LINEZ_FORBIDDEN']
        linez_balmer = galfit_table['LINEZ_BALMER']
        linesigma_forbidden = galfit_table['LINESIGMA_FORBIDDEN']
        linesigma_balmer = galfit_table['LINESIGMA_BALMER']
        
        npixpercamera = [len(gw) for gw in galwave]

        linevshift_forbidden = (linez_forbidden - redshift) * C_LIGHT # [km/s]
        linevshift_balmer = (linez_balmer - redshift) * C_LIGHT # [km/s]

        EMLine = EMLineModel(linevshift_forbidden=linevshift_forbidden,
                             linevshift_balmer=linevshift_balmer,
                             linesigma_forbidden=linesigma_forbidden,
                             linesigma_balmer=linesigma_balmer,
                             redshift=redshift, emlineR=galres,
                             npixpercamera=npixpercamera)
        # skip linevshift_[forbidden,balmer] and linesigma_[forbidden,balmer]
        lineargs = [galfit_table[linename.upper()] for linename in EMLine.param_names[4:]] 
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

        # convert the vshifts to redshifts
        result['linez_forbidden'] = redshift + result['linevshift_forbidden'] / C_LIGHT
        result['linez_balmer'] = redshift + result['linevshift_balmer'] / C_LIGHT

        result['linez_forbidden_ivar'] = result['linevshift_forbidden_ivar'] * C_LIGHT**2
        result['linez_balmer_ivar'] = result['linevshift_balmer_ivar'] * C_LIGHT**2

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
        bigax.text(0.95, 0.80, r'{} {}'.format(objinfo['linez_balmer'], objinfo['linez_forbidden']),
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
                #print(ii, linename, len(indx))
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
                #print(linenames[iax], xlim, np.diff(xlim))
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
            print('Writing {}'.format(png))
            fig.savefig(png)
