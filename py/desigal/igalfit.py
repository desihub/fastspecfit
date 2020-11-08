"""
desigal.igalfit
===============

"""
import os, pdb
import numpy as np
from astropy.table import Table, Column
from astropy.modeling import Fittable1DModel

from scipy import constants
C_LIGHT = constants.c / 1000.0 # [km/s]

def init_output(tile, night, zbest, CFit):
    """Initialize the output data table.

    CFit - ContinuumFit class

    """
    from astropy.table import QTable
    import astropy.units as u
    nobj = len(zbest)

    # Grab info on the emission lines and the continuum.
    linetable = get_linetable()
    nssp_coeff = len(CFit.sspinfo)

    out = QTable()
    for zbestcol in ['TARGETID', 'Z', 'ZERR']:#, 'SPECTYPE', 'DELTACHI2']
        out[zbestcol] = zbest[zbestcol]
    out.add_column(Column(name='NIGHT', data=np.repeat(night, nobj)), index=0)
    out.add_column(Column(name='TILE', data=np.repeat(tile, nobj)), index=0)
    out.add_column(Column(name='CONTINUUM_COEFF', length=nobj, shape=(nssp_coeff,), dtype='f4'))
    out.add_column(Column(name='CONTINUUM_CHI2', length=nobj, dtype='f4')) # reduced chi2
    out.add_column(Column(name='CONTINUUM_EBV', length=nobj, dtype='f4', unit=u.mag))
    out.add_column(Column(name='CONTINUUM_VDISP', length=nobj, dtype='f4', unit=u.kilometer/u.second))
    out.add_column(Column(name='LINEZ', length=nobj, dtype='f4'))
    out.add_column(Column(name='LINESIGMA', length=nobj, dtype='f4', unit=u.kilometer/u.second))
    #out.add_column(Column(name='LINEZ_FORBIDDEN', length=nobj, dtype='f4'))
    #out.add_column(Column(name='LINEZ_BALMER', length=nobj, dtype='f4'))
    #out.add_column(Column(name='LINESIGMA_FORBIDDEN', length=nobj, dtype='f4', unit=u.kilometer / u.second))
    #out.add_column(Column(name='LINESIGMA_BALMER', length=nobj, dtype='f4', unit=u.kilometer / u.second))

    for line in linetable['name']:
        line = line.upper()
        out.add_column(Column(name='{}_FLUX'.format(line), length=nobj, dtype='f4',
                              unit=u.erg/(u.second*u.cm**2)))
        out.add_column(Column(name='{}_FLUX_IVAR'.format(line), length=nobj, dtype='f4',
                              unit=u.second**2*u.cm**4/u.erg**2))
        out.add_column(Column(name='{}_CONT'.format(line), length=nobj, dtype='f4',
                              unit=u.erg/(u.second*u.cm**2*u.Angstrom)))
        out.add_column(Column(name='{}_CONT_IVAR'.format(line), length=nobj, dtype='f4',
                              unit=u.second**2*u.cm**4*u.Angstrom**2/u.erg**2))
        out.add_column(Column(name='{}_EW'.format(line), length=nobj, dtype='f4',
                              unit=u.Angstrom))
        out.add_column(Column(name='{}_EW_IVAR'.format(line), length=nobj, dtype='f4',
                              unit=1/u.Angstrom**2))
    
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
    from desispec.interpolation import resample_flux
    return galR.dot(resample_flux(galwave, sspwave, sspflux))

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

        # dust parameters
        self.dustslope = 0.7
        
    def redshift_smooth_and_resample(self, galwave, galres, redshift, plot=False):
        """Redshift, convolve with the spectral resolution, and 
        resample in wavelength.
        
        """
        import multiprocessing

        oneplusz = 1 + redshift
                    
        # loop over cameras then SSP ages
        smoothflux = []
        for icamera in [0, 1, 2]: # iterate on cameras
            args = [(self.sspflux[:, iage] / oneplusz, self.sspwave * oneplusz, galwave[icamera], 
                     galres[icamera]) for iage in np.arange(self.nage)]
            with multiprocessing.Pool(self.nproc) as pool:
                smoothflux.append(np.array(pool.map(_smooth_and_resample, args)).T)

        if plot:
            import matplotlib.pyplot as plt
            ww = 110
            plt.plot(self.ssp.wave * (1+self.zredrock), self.ssp.flux[:, ww], color='gray', alpha=0.5)
            for ii in [0, 1, 2]: # iterate over cameras
                plt.plot(self.galwave[ii], sspflux[ii][:, 110], color='k')
            plt.xlim(3900, 4000)
            
        return smoothflux # [npix, nage]

    def fnnls_continuum(self, galwave, galflux, galivar, galres, redshift, plot=False):
        """Fit the continuum using fast NNLS.
        https://github.com/jvendrow/fnnls
        
        """
        from fnnls import fnnls
        
        sspflux = self.redshift_smooth_and_resample(galwave, galres, redshift) 

        # combine the cameras and fit
        _galflux = np.hstack(galflux)
        _galivar = np.hstack(galivar)
        _sspflux = np.concatenate(sspflux, axis=0) # [npix, nage]
        
        # Do a fast initial fit of the stellar continuum.
        print('ToDo: compute reduced chi2, mask emission lines, and compute the light-weighted age.')
        coeffs = fnnls(_sspflux * np.sqrt(_galivar[:, None]), _galflux * np.sqrt(_galivar))[0]

        _continuum = _sspflux.dot(coeffs)

        # light-weighted age
        weight = coeffs[coeffs > 0]
        print(np.sum(weight * self.sspinfo['age'][coeffs > 0]) / np.sum(weight) / 1e9) # [Gyr]

        # Repackage into individual cameras
        continuum = []
        npixpercamera = [len(gw) for gw in galwave]
        npix = np.hstack([0, npixpercamera])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            continuum.append(_continuum[ipix:jpix])

        if plot:
            galwave = np.hstack(self.galwave)
            fig, ax = plt.subplots(figsize=(14, 8))
            for ii in [0, 1, 2]: # iterate over cameras
                ax.plot(self.galwave[ii], self.galflux[ii])
            ax.plot(galwave, continuum, alpha=0.5, color='k')
            #ax.set_xlim(3850, 3950)
            
        return coeffs, continuum

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

#    def emlinemodel(self, params, log10wave):
#        """Emission-line model."""
#        zline, linesigma = params[0], params[1]
#
#        linenames = self.linetable['name'].data
#        linefluxes = params[2:]
#        
#        emlinemodel = []
#        for ii in [0, 1, 2]: # iterate over cameras
#            ipix = np.sum(self.npixpercamera[:ii+1])
#            jpix = np.sum(self.npixpercamera[:ii+2])
#            _emlinewave = log10wave[ipix:jpix]
#            _emlinemodel = np.zeros_like(_emlinewave)
#            
#            for linename, lineflux in zip(linenames, linefluxes):
#                restlinewave = self.linetable[self.linetable['name'] == linename]['restwave'][0]
#                zlinewave = np.log10(restlinewave * (1.0 + zline)) # redshifted wavelength [log-10 Angstrom]
#                log10sigma = linesigma / C_LIGHT / np.log(10)      # line-width [log-10 Angstrom]
#                lineamp = lineflux / (np.sqrt(2.0 * np.pi) * log10sigma)
#            
#                # Construct the spectrum [erg/s/cm2/A, rest]
#                ww = np.abs(_emlinewave - zlinewave) < 20 * log10sigma
#                if np.count_nonzero(ww) > 0:
#                    #print(linename, 10**zlinewave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
#                    _emlinemodel[ww] += lineamp * np.exp(-0.5 * (_emlinewave[ww]-zlinewave)**2 / log10sigma**2)
#
#            # optionally convolve with the spectral resolution
#            if self.emlineR is not None:
#                _emlinemomdel = self.emlineR[ii].dot(_emlinemodel)
#            
#            #plt.plot(10**_emlinewave, _emlinemodel)
#            #plt.plot(10**_emlinewave, self.emlineR[ii].dot(_emlinemodel))
#            #plt.xlim(3870, 3920) ; plt.show()
#            #pdb.set_trace()
#            emlinemodel.append(_emlinemodel)
#
#        return np.hstack(emlinemodel)
#    
#    def _emlinemodel_resid(self, params, log10wave, emlineflux, isigma):
#        """Wrapper cost function for scipy.optimize.least_squares."""
#        emlinemodel = self.emlinemodel(params, log10wave)
#        return isigma * (emlinemodel - emlineflux)

def get_linetable():
    # read this from a file!
    linetable = Table()
    linetable['name'] = ['oii_3726', 'oii_3729', 
                         'oiii_4959', 'oiii_5007', 
                         'nii_6548', 'nii_6584', 
                         'hdelta', 'hgamma', 'hbeta', 'halpha']
    linetable['flux'] = [0.73, 1.0, 
                         1.0, 2.875, 
                         1.0, 2.936, 
                         0.259, 0.468, 1.0, 2.863]
    linetable['restwave'] = [3727.092, 3729.874, 
                             4960.295, 5008.239, 
                             6549.852, 6585.277, 
                             4102.892, 4341.684, 4862.683, 6564.613]
    return linetable    

class EMLineModel(Fittable1DModel):
    """Class to model the emission-line spectra.

    """
    from astropy.modeling import Parameter

    # NB! The order of the parameters here matters!
    zline = Parameter(name='zline', default=0.1, bounds=(-0.05, 2.0)) # line-redshift
    linesigma = Parameter(name='linesigma', default=50.0, bounds=(0.1, 350)) # line-sigma [km/s]

    # Fragile because the lines are hard-coded--
    oii_3726_flux = Parameter(name='oii_3726_flux', default=0.73)
    oii_3729_flux = Parameter(name='oii_3729_flux', default=1.0)
    oiii_4959_flux = Parameter(name='oiii_4959_flux', default=1.0)
    oiii_5007_flux = Parameter(name='oiii_5007_flux', default=2.8875)
    nii_6548_flux = Parameter(name='nii_6548_flux', default=1.0)
    nii_6584_flux = Parameter(name='nii_6584_flux', default=2.936)
    hdelta_flux = Parameter(name='hdelta_flux', default=0.259)
    hgamma_flux = Parameter(name='hgamma_flux', default=0.468)
    hbeta_flux = Parameter(name='hbeta_flux', default=1.0)
    halpha_flux = Parameter(name='halpha_flux', default=2.863)

    # tie the [NII] and [OIII] line-strengths together
    def tie_oiii(model):
        return model.oiii_5007_flux / 2.8875
    oiii_4959_flux.tied = tie_oiii

    def tie_nii(model):
        return model.nii_6584_flux / 2.936
    nii_6548_flux.tied = tie_nii
    
    def __init__(self, zline=zline.default, linesigma=linesigma.default,
                 oii_3726_flux=oii_3726_flux.default, 
                 oii_3729_flux=oii_3729_flux.default, 
                 oiii_4959_flux=oiii_4959_flux.default, 
                 oiii_5007_flux=oiii_5007_flux.default, 
                 nii_6548_flux=nii_6548_flux.default, 
                 nii_6584_flux=nii_6584_flux.default, 
                 hdelta_flux=hdelta_flux.default, 
                 hgamma_flux=hgamma_flux.default, 
                 hbeta_flux=hbeta_flux.default, 
                 halpha_flux=halpha_flux.default, 
                 emlineR=None, npixpercamera=None, **kwargs):
        """Initialize the emission-line model.
        
        emlineR - 
        
        """
        self.linetable = get_linetable()
        self.emlineR = emlineR
        self.npixpercamera = np.hstack([0, npixpercamera])
        
        super(EMLineModel, self).__init__(
            zline=zline, linesigma=linesigma,
            oii_3726_flux=oii_3726_flux,
            oii_3729_flux=oii_3729_flux,
            oiii_4959_flux=oiii_4959_flux,
            oiii_5007_flux=oiii_5007_flux,
            nii_6548_flux=nii_6548_flux,
            nii_6584_flux=nii_6584_flux,
            hdelta_flux=hdelta_flux,
            hgamma_flux=hgamma_flux,
            hbeta_flux=hbeta_flux,
            halpha_flux=halpha_flux, **kwargs)
        
    def evaluate(self, log10wave, *args):
        """Evaluate the emission-line model.
        emlineR=None, npixpercamera=None, 
        """ 
        zline, linesigma = args[0], args[1]

        linenames = self.linetable['name'].data
        linefluxes = args[2:]
        
        emlinemodel = []
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(self.npixpercamera[:ii+1])
            jpix = np.sum(self.npixpercamera[:ii+2])
            _emlinewave = log10wave[ipix:jpix]
            _emlinemodel = np.zeros_like(_emlinewave)
            
            for linename, lineflux in zip(linenames, linefluxes):
                restlinewave = self.linetable[self.linetable['name'] == linename]['restwave'][0]
                zlinewave = np.log10(restlinewave * (1.0 + zline)) # redshifted wavelength [log-10 Angstrom]
                log10sigma = linesigma / C_LIGHT / np.log(10)      # line-width [log-10 Angstrom]
                lineamp = lineflux / (np.sqrt(2.0 * np.pi) * log10sigma)
            
                # Construct the spectrum [erg/s/cm2/A, rest]
                ww = np.abs(_emlinewave - zlinewave) < 20 * log10sigma
                if np.count_nonzero(ww) > 0:
                    #print(linename, 10**zlinewave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
                    _emlinemodel[ww] += lineamp * np.exp(-0.5 * (_emlinewave[ww]-zlinewave)**2 / log10sigma**2)

            # optionally convolve with the spectral resolution
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

        self.fitter = fitting.LevMarLSQFitter()
                
    def chi2(self, bestfit, emlinewave, emlineflux, emlineivar):
        """Compute the reduced chi^2."""
        dof = len(emlinewave) - len(bestfit.parameters)
        emlinemodel = bestfit(emlinewave)
        chi2 = np.sum(emlineivar * (emlineflux - emlinemodel)**2) / dof
        return chi2
    
    def fit(self, galwave, galflux, galivar, galres, continuum,
            redshift, verbose=False):
        """Perform the fit minimization / chi2 minimization.
        
        EMLineModel object
        FC - ContinuumFit object
        
        """
        emlinewave = np.hstack(galwave)
        emlineflux = np.hstack(galflux) - np.hstack(continuum)
        emlineivar = np.hstack(galivar)
        npixpercamera = [len(gw) for gw in galwave]
        
        self.EMLineModel = EMLineModel(zline=redshift, emlineR=galres,
                                       npixpercamera=npixpercamera)

        weights = 1 / np.sqrt(emlineivar)
        bestfit = self.fitter(self.EMLineModel, np.log10(emlinewave), 
                              emlineflux, weights=weights)
        chi2 = self.chi2(bestfit, emlinewave, emlineflux, emlineivar).astype('f4')
        
        # Pack the results in a dictionary and return.
        # https://gist.github.com/eteq/1f3f0cec9e4f27536d52cd59054c55f2
        result = {
            'converged': False,
            'fit_message': self.fitter.fit_info['message'],
            'nparam': len(self.EMLineModel.parameters),
            'npix': len(emlinewave),
            'dof': len(emlinewave) - len(self.EMLineModel.parameters),
            'chi2': chi2,
        }
        #for param in bestfit.param_names:
        #    result.update({param: getattr(bestfit, param).value})
        
        # uncertainties
        if self.fitter.fit_info['param_cov'] is not None:
            cov = self.fitter.fit_info['param_cov']
            unc = np.diag(cov)**0.5
            result['converged'] = True
        else:
            cov = np.zeros((nparams, nparams))
            unc = np.zeros(nparams)

        # Need to be careful about uncertainties for tied parameters--
        # https://github.com/astropy/astropy/issues/7202
        
        #err_params = np.sqrt(np.diag(fitter.fit_info['param_cov']))
        #err = model.copy()
        #fitting._fitter_to_model_params(err, err_params)            
            
        count = 0
        for ii, pp in enumerate(bestfit.param_names):
            pinfo = getattr(bestfit, pp)
            result.update({bestfit.param_names[ii]: bestfit.parameters[ii].astype('f4')})
            if pinfo.fixed:
                result.update({bestfit.param_names[ii]+'_err': np.float32(0.0)})
            elif pinfo.tied:
                pass # hack! see https://github.com/astropy/astropy/issues/7202
            else:
                result.update({bestfit.param_names[ii]+'_err': unc[count].astype('f4')})
                count += 1

        emlinemodel = bestfit(np.log10(emlinewave))
        
        return result, emlinemodel
    
    def emlineplot(self, galwave, galflux, galivar, continuum,
                   emlinemodel, redshift, png=None):
        """Plot the emission-line spectrum and best-fitting model.

        """
        import matplotlib.pyplot as plt
        
        emlinewave = np.hstack(galwave)
        emlineflux = np.hstack(galflux) - np.hstack(continuum)
        emlinesigma = 1 / np.sqrt(np.hstack(galivar))

        fig, ax = plt.subplots(1, 4, figsize=(16, 5))#, sharey=True)
        for xx, minwave, maxwave in zip(ax, (3725, 4850, 4950, 6550), (3730, 4870, 5015, 6595)):
            wmin, wmax = np.array([minwave, maxwave]) * (1+redshift) + np.array([-20, +20])
            indx = np.where((emlinewave > wmin) * (emlinewave < wmax))[0]
            if len(indx) > 5:
                xx.errorbar(emlinewave[indx], emlineflux[indx], yerr=emlinesigma[indx])
                xx.plot(emlinewave[indx], emlinemodel[indx])
            #xx.set_ylim(-3, 13)
        plt.subplots_adjust(wspace=0.05)

        if png:
            print('Writing {}'.format(png))
            fig.savefig(png)
