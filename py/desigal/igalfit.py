"""
desigal.igalfit
===============

"""
import os
import numpy as np

from scipy import constants
C_LIGHT = constants.c / 1000.0 # [km/s]

def _smooth_and_resample(args):
    """Wrapper for multiprocessing."""
    return smooth_and_resample(*args)

def smooth_and_resample(sspflux, sspwave, galwave, galR):
    """
    sspflux[npix] - redshifted SSP template
    sspwave[npix] - redshifted SSP wavelength
    
    """
    from desispec.interpolation import resample_flux
    return galR.dot(resample_flux(galwave, sspwave, sspflux))

class CKCz14(object):
    """Class to handle the CKCz14 SSPs.
    
    The SSPs have been converted to FITS format using the script fsps2fits.py
    
    """
    def __init__(self, metallicity='Z0.0190', minwave=0.0, maxwave=1e4, verbose=True):
        """Write me.

        """
        import fitsio
        from astropy.table import Table
        from astropy.cosmology import FlatLambdaCDM

        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)        

        self.metallicity = metallicity
        self.Z = float(metallicity[1:])
        self.library = 'CKC14z'
        self.isochrone = 'Padova' # would be nice to get MIST in here
        self.imf = 'Kroupa'

        # Don't hard-code the path!
        self.ssppath = os.getenv('DESIGAL_TEMPLATES')
        self.sspfile = os.path.join(self.ssppath, 'SSP_{}_{}_{}_{}.fits'.format(
            self.isochrone, self.library, self.imf, self.metallicity))

        if verbose:
            print('Reading {}'.format(self.sspfile))
        wave = fitsio.read(self.sspfile, ext='WAVE')
        flux = fitsio.read(self.sspfile, ext='FLUX')
        
        keep = np.where((wave >= minwave) * (wave <= maxwave))[0]
        self.wave = wave[keep]
        self.flux = flux[keep, :]
        self.info = Table(fitsio.read(self.sspfile, ext='METADATA'))

        self.nage = len(self.info['age'])
        self.npix = len(wave)
        
    def smooth_and_resample(self, galwave, galres, redshift, nproc=1):
        """Write me.
        
        galwave - list of wavelength array for each camera
        
        
        """
        oneplusz = 1 + redshift
                    
        # loop over cameras then SSP ages
        smoothflux = []
        for icamera in [0, 1, 2]: # iterate on cameras
            args = [(self.flux[:, iage] / oneplusz, self.wave * oneplusz, galwave[icamera], 
                     galres[icamera]) for iage in np.arange(self.nage)]
            with multiprocessing.Pool(nproc) as pool:
                smoothflux.append(np.array(pool.map(_smooth_and_resample, args)).T)

        return smoothflux

class FitContinuum():
    def __init__(self, ssp, specobj, zbest, iobj=0, nproc=1):
        """Model the stellar continuum.
    
        ssp - CKCz14 Class (fixed metallicity)
        specobj - Spectra Class data
        zbest - astropy Table with redshift info
        iobj - index of object to fit
        nproc - number of cores to use for multiprocessing
    
        """
        from astropy.modeling import Parameter
        from desispec.resolution import Resolution

        self.nproc = nproc
        self.ssp = ssp
        
        # Select and repackage the spectrum we want to fit.
        self.npix, self.galwave, self.galflux, self.galivar, self.galres = [], [], [], [], []
        for camera in ('b', 'r', 'z'):
            self.npix.append(len(specobj.wave[camera]))
            self.galwave.append(specobj.wave[camera])
            self.galflux.append(specobj.flux[camera][iobj, :])
            self.galivar.append(specobj.ivar[camera][iobj, :])
            self.galres.append(Resolution(specobj.resolution_data[camera][iobj, :, :]))
            #self.galres.append(specobj.resolution_data[camera][iobj, :, :])
            
        self.zredrock = zbest['Z'][iobj]
            
    def redshift_smooth_and_resample(self, plot=False):
        """Redshift, convolve with the spectral resolution, and 
        resample in wavelength.
        
        """
        sspflux = self.ssp.smooth_and_resample(self.galwave, 
                                               self.galres, 
                                               self.zredrock, 
                                               nproc=self.nproc) # [npix, nage]
        
        if plot:
            ww = 110
            plt.plot(self.ssp.wave * (1+self.zredrock), self.ssp.flux[:, ww], color='gray', alpha=0.5)
            for ii in [0, 1, 2]: # iterate over cameras
                plt.plot(self.galwave[ii], sspflux[ii][:, 110], color='k')
            plt.xlim(3900, 4000)
            
        return sspflux
        
    def fit_fnnls(self, plot=False):
        """Fit using fast NNLS.
        https://github.com/jvendrow/fnnls
        
        """
        from fnnls import fnnls
        
        sspflux = np.concatenate(self.redshift_smooth_and_resample(), axis=0)
        if True:
            ebv, dustslope = 0.3, -0.7
            galwave = np.hstack(self.galwave)
            klambda = 5.9 * (galwave / 5500)**dustslope
            atten = 10**(-0.4 * ebv * klambda)
            sspflux *= atten[:, None]
        
        galflux = np.hstack(self.galflux)
        galivar = np.hstack(self.galivar)

        # Do a fast initial fit of the stellar continuum.
        # ToDo: mask emission lines!
        self.fnnls_coeffs = fnnls(sspflux * np.sqrt(galivar[:, None]), galflux * np.sqrt(galivar))[0]
        continuum = sspflux.dot(self.fnnls_coeffs)
        print('Compute chi2 and store the coeffs.')
        
        self.continuum = []
        npix = np.hstack([0, self.npix])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            self.continuum.append(continuum[ipix:jpix])

        if plot:
            galwave = np.hstack(self.galwave)

            fig, ax = plt.subplots(figsize=(14, 8))
            for ii in [0, 1, 2]: # iterate over cameras
                ax.plot(self.galwave[ii], self.galflux[ii])
            ax.plot(galwave, continuum, alpha=0.5, color='k')
            #ax.set_xlim(3850, 3950)
            
        return continuum    
        
    def fit_continuum(self, ContinuumModel):
        """More detailed continuum model fit, including dust attenuation.
        
        """
        from astropy.modeling import fitting
                
        fitter = fitting.LevMarLSQFitter()
        
        wave = np.hstack(self.galwave)
        flux = np.hstack(self.galflux)
        ivar = np.hstack(self.galivar)
        
        weights = 1 / np.sqrt(ivar)
        continuum_fit = fitter(ContinuumModel, wave, flux, weights=1/np.sqrt(ivar))

        pdb.set_trace()
