"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import pdb # for debugging

import os, time
import numpy as np
import multiprocessing

import astropy.units as u
from astropy.modeling import Fittable1DModel

from fastspecfit.util import C_LIGHT
from fastspecfit.continuum import ContinuumTools
from desiutil.log import get_logger
log = get_logger()

def read_emlines():
    """Read the set of emission lines of interest.

    ToDo: add lines to mask during continuum-fitting but which we do not want to
    emission-line fit.

    """
    from astropy.table import Table
    from pkg_resources import resource_filename
    
    linefile = resource_filename('fastspecfit', 'data/emlines.ecsv')    
    linetable = Table.read(linefile, format='ascii.ecsv', guess=False)
    
    return linetable    

def _tie_neon_sigma(model):
    return model.nev_3426_sigma
def _tie_oii_3726_sigma(model):
    return model.oii_3729_sigma
def _tie_oiii_4959_sigma(model):
    return model.oiii_5007_sigma
def _tie_nii_6548_sigma(model):
    return model.nii_6584_sigma
def _tie_sii_6731_sigma(model):
    return model.sii_6716_sigma
def _tie_siii_9530_sigma(model):
    return model.siii_9068_sigma

def _tie_hbeta_sigma(model):
    return model.hbeta_sigma
def _tie_hbeta_vshift(model):
    return model.hbeta_vshift

def _tie_oiii_5007_sigma(model):
    return model.oiii_5007_sigma
def _tie_oiii_5007_vshift(model):
    return model.oiii_5007_vshift

def _tie_mgii_2800_sigma(model):
    return model.mgii_2800_sigma
def _tie_mgii_2800_vshift(model):
    return model.mgii_2800_vshift

def _tie_balmer_lines(model):
    for pp in model.param_names:
        if 'sigma' in pp and pp[0] == 'h' and pp != 'hbeta_sigma':
            getattr(model, pp).tied = _tie_hbeta_sigma
        if 'vshift' in pp and pp[0] == 'h' and pp != 'hbeta_vshift':
            getattr(model, pp).tied = _tie_hbeta_vshift
            
    model.hbeta_sigma.tied = False
    model.hbeta_vshift.tied = False
    
    return model

def _tie_forbidden_lines(model):
    for pp in model.param_names:
        if 'sigma' in pp and pp[0] != 'h' and 'mgii' not in pp and pp != 'oiii_5007_sigma':
            getattr(model, pp).tied = _tie_oiii_5007_sigma
        if 'vshift' in pp and pp[0] != 'h' and 'mgii' not in pp and pp != 'oiii_5007_vshift':
            getattr(model, pp).tied = _tie_oiii_5007_vshift
            
    model.oiii_5007_sigma.tied = False
    model.oiii_5007_vshift.tied = False
    model.mgii_2800_sigma.tied = False
    model.mgii_2800_vshift.tied = False
    
    return model

def _tie_qso_lines(model):
    for pp in model.param_names:
        if 'sigma' in pp and 'mgii' in pp:
            getattr(model, pp).tied = _tie_mgii_2800_sigma
        if 'vshift' in pp and 'mgii' in pp:
            getattr(model, pp).tied = _tie_mgii_2800_vshift
    model.mgii_2800_sigma.tied = False
    model.mgii_2800_vshift.tied = False
    return model

def _tie_all_lines(model):
    for pp in model.param_names:
        if 'sigma' in pp and pp != 'hbeta_sigma':
            getattr(model, pp).tied = _tie_hbeta_sigma
        if 'vshift' in pp and pp != 'hbeta_vshift':
            getattr(model, pp).tied = _tie_hbeta_vshift
    return model

class EMLineModel(Fittable1DModel):
    """Class to model the emission-line spectra.

    """
    from astropy.modeling import Parameter

    # NB! The order of the parameters here matters!
    vmaxshift = 300.0

    maxsigma_narrow = 500.0
    maxsigma_broad = 3000.0

    initsigma_narrow = 75.0
    initsigma_broad = 300.0

    # Fragile because the lines are hard-coded--
    mgii_2800_amp = Parameter(name='mgii_2800_amp', default=3.0)
    #mgii_2800b_amp = Parameter(name='mgii_2800b_amp', default=1.0)
    nev_3346_amp = Parameter(name='nev_3346_amp', default=0.1)
    nev_3426_amp = Parameter(name='nev_3426_amp', default=0.1)
    oii_3726_amp = Parameter(name='oii_3726_amp', default=1.0)
    oii_3729_amp = Parameter(name='oii_3729_amp', default=1.0)
    neiii_3869_amp = Parameter(name='neiii_3869_amp', default=0.3)
    oiii_4959_amp = Parameter(name='oiii_4959_amp', default=1.0)
    oiii_5007_amp = Parameter(name='oiii_5007_amp', default=3.0)
    hepsilon_amp = Parameter(name='hepsilon_amp', default=0.5)
    hdelta_amp = Parameter(name='hdelta_amp', default=0.5)
    hgamma_amp = Parameter(name='hgamma_amp', default=0.5)
    hbeta_amp = Parameter(name='hbeta_amp', default=1.0)
    halpha_amp = Parameter(name='halpha_amp', default=3.0)
    nii_6548_amp = Parameter(name='nii_6548_amp', default=1.0)
    nii_6584_amp = Parameter(name='nii_6584_amp', default=3.0)
    sii_6716_amp = Parameter(name='sii_6716_amp', default=1.0)
    sii_6731_amp = Parameter(name='sii_6731_amp', default=1.0)
    siii_9068_amp = Parameter(name='siii_9068_amp', default=0.3)
    siii_9530_amp = Parameter(name='siii_9530_amp', default=0.3)

    mgii_2800_vshift = Parameter(name='mgii_2800_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    #mgii_2800b_vshift = Parameter(name='mgii_2800b_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    nev_3346_vshift = Parameter(name='nev_3346_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    nev_3426_vshift = Parameter(name='nev_3426_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    oii_3726_vshift = Parameter(name='oii_3726_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    oii_3729_vshift = Parameter(name='oii_3729_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    neiii_3869_vshift = Parameter(name='neiii_3869_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    oiii_4959_vshift = Parameter(name='oiii_4959_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    oiii_5007_vshift = Parameter(name='oiii_5007_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    hepsilon_vshift = Parameter(name='hepsilon_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    hdelta_vshift = Parameter(name='hdelta_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    hgamma_vshift = Parameter(name='hgamma_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    hbeta_vshift = Parameter(name='hbeta_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    halpha_vshift = Parameter(name='halpha_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    nii_6548_vshift = Parameter(name='nii_6548_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    nii_6584_vshift = Parameter(name='nii_6584_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    sii_6716_vshift = Parameter(name='sii_6716_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    sii_6731_vshift = Parameter(name='sii_6731_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    siii_9068_vshift = Parameter(name='siii_9068_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])
    siii_9530_vshift = Parameter(name='siii_9530_vshift', default=0.0, bounds=[-vmaxshift, +vmaxshift])

    mgii_2800_sigma = Parameter(name='mgii_2800_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    #mgii_2800b_sigma = Parameter(name='mgii_2800b_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    nev_3346_sigma = Parameter(name='nev_3346_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    nev_3426_sigma = Parameter(name='nev_3426_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    oii_3726_sigma = Parameter(name='oii_3726_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    oii_3729_sigma = Parameter(name='oii_3729_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    neiii_3869_sigma = Parameter(name='neiii_3869_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    oiii_4959_sigma = Parameter(name='oiii_4959_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    oiii_5007_sigma = Parameter(name='oiii_5007_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    hepsilon_sigma = Parameter(name='hepsilon_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    hdelta_sigma = Parameter(name='hdelta_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    hgamma_sigma = Parameter(name='hgamma_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    hbeta_sigma = Parameter(name='hbeta_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    halpha_sigma = Parameter(name='halpha_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_broad])
    nii_6548_sigma = Parameter(name='nii_6548_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    nii_6584_sigma = Parameter(name='nii_6584_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    sii_6716_sigma = Parameter(name='sii_6716_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    sii_6731_sigma = Parameter(name='sii_6731_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    siii_9068_sigma = Parameter(name='siii_9068_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    siii_9530_sigma = Parameter(name='siii_9530_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])

    # doublet ratios are always tied
    def _tie_nii_amp(model):
        return model.nii_6584_amp / 2.936
    def _tie_oiii_amp(model):
        return model.oiii_5007_amp / 2.8875
    
    nii_6548_amp.tied = _tie_nii_amp
    oiii_4959_amp.tied = _tie_oiii_amp
    
    def __init__(self,
                 mgii_2800_amp=mgii_2800_amp.default,
                 #mgii_2800b_amp=mgii_2800b_amp.default,
                 nev_3346_amp=nev_3346_amp.default,
                 nev_3426_amp=nev_3426_amp.default,
                 oii_3726_amp=oii_3726_amp.default, 
                 oii_3729_amp=oii_3729_amp.default, 
                 neiii_3869_amp=neiii_3869_amp.default,
                 oiii_4959_amp=oiii_4959_amp.default, 
                 oiii_5007_amp=oiii_5007_amp.default, 
                 hepsilon_amp=hepsilon_amp.default, 
                 hdelta_amp=hdelta_amp.default, 
                 hgamma_amp=hgamma_amp.default, 
                 hbeta_amp=hbeta_amp.default, 
                 halpha_amp=halpha_amp.default,
                 nii_6548_amp=nii_6548_amp.default, 
                 nii_6584_amp=nii_6584_amp.default, 
                 sii_6716_amp=sii_6716_amp.default, 
                 sii_6731_amp=sii_6731_amp.default,
                 siii_9068_amp=siii_9068_amp.default,
                 siii_9530_amp=siii_9530_amp.default,
                 
                 mgii_2800_vshift=mgii_2800_vshift.default,
                 #mgii_2800b_vshift=mgii_2800b_vshift.default,
                 nev_3346_vshift=nev_3346_vshift.default,
                 nev_3426_vshift=nev_3426_vshift.default,
                 oii_3726_vshift=oii_3726_vshift.default, 
                 oii_3729_vshift=oii_3729_vshift.default, 
                 neiii_3869_vshift=neiii_3869_vshift.default,
                 oiii_4959_vshift=oiii_4959_vshift.default, 
                 oiii_5007_vshift=oiii_5007_vshift.default, 
                 hepsilon_vshift=hepsilon_vshift.default, 
                 hdelta_vshift=hdelta_vshift.default, 
                 hgamma_vshift=hgamma_vshift.default, 
                 hbeta_vshift=hbeta_vshift.default, 
                 halpha_vshift=halpha_vshift.default,
                 nii_6548_vshift=nii_6548_vshift.default, 
                 nii_6584_vshift=nii_6584_vshift.default, 
                 sii_6716_vshift=sii_6716_vshift.default, 
                 sii_6731_vshift=sii_6731_vshift.default,
                 siii_9068_vshift=siii_9068_vshift.default,
                 siii_9530_vshift=siii_9530_vshift.default,
                 
                 mgii_2800_sigma=mgii_2800_sigma.default,
                 #mgii_2800b_sigma=mgii_2800b_sigma.default,
                 nev_3346_sigma=nev_3346_sigma.default,
                 nev_3426_sigma=nev_3426_sigma.default,
                 oii_3726_sigma=oii_3726_sigma.default, 
                 oii_3729_sigma=oii_3729_sigma.default, 
                 neiii_3869_sigma=neiii_3869_sigma.default,
                 oiii_4959_sigma=oiii_4959_sigma.default, 
                 oiii_5007_sigma=oiii_5007_sigma.default, 
                 hepsilon_sigma=hepsilon_sigma.default, 
                 hdelta_sigma=hdelta_sigma.default, 
                 hgamma_sigma=hgamma_sigma.default, 
                 hbeta_sigma=hbeta_sigma.default, 
                 halpha_sigma=halpha_sigma.default,
                 nii_6548_sigma=nii_6548_sigma.default, 
                 nii_6584_sigma=nii_6584_sigma.default, 
                 sii_6716_sigma=sii_6716_sigma.default, 
                 sii_6731_sigma=sii_6731_sigma.default,
                 siii_9068_sigma=siii_9068_sigma.default,
                 siii_9530_sigma=siii_9530_sigma.default,
                 
                 redshift=None,
                 emlineR=None, npixpercamera=None,
                 log10wave=None, **kwargs):
        """Initialize the emission-line model.
        
        emlineR -
        redshift - required keyword
        
        """
        self.linetable = read_emlines()

        self.nparam_perline = 3 # amplitude, velocity shift, and line-width
        self.nline = len(self.param_names) // self.nparam_perline

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
            mgii_2800_amp=mgii_2800_amp,
            #mgii_2800b_amp=mgii_2800b_amp,
            nev_3346_amp=nev_3346_amp,
            nev_3426_amp=nev_3426_amp,
            oii_3726_amp=oii_3726_amp,
            oii_3729_amp=oii_3729_amp,
            neiii_3869_amp=neiii_3869_amp,
            oiii_4959_amp=oiii_4959_amp,
            oiii_5007_amp=oiii_5007_amp,
            hepsilon_amp=hepsilon_amp,
            hdelta_amp=hdelta_amp,
            hgamma_amp=hgamma_amp,
            hbeta_amp=hbeta_amp,
            halpha_amp=halpha_amp,
            nii_6548_amp=nii_6548_amp,
            nii_6584_amp=nii_6584_amp,
            sii_6716_amp=sii_6716_amp,
            sii_6731_amp=sii_6731_amp,
            siii_9068_amp=siii_9068_amp,
            siii_9530_amp=siii_9530_amp,
            
            mgii_2800_vshift=mgii_2800_vshift,
            #mgii_2800b_vshift=mgii_2800b_vshift,
            nev_3346_vshift=nev_3346_vshift,
            nev_3426_vshift=nev_3426_vshift,
            oii_3726_vshift=oii_3726_vshift,
            oii_3729_vshift=oii_3729_vshift,
            neiii_3869_vshift=neiii_3869_vshift,
            oiii_4959_vshift=oiii_4959_vshift,
            oiii_5007_vshift=oiii_5007_vshift,
            hepsilon_vshift=hepsilon_vshift,
            hdelta_vshift=hdelta_vshift,
            hgamma_vshift=hgamma_vshift,
            hbeta_vshift=hbeta_vshift,
            halpha_vshift=halpha_vshift,
            nii_6548_vshift=nii_6548_vshift,
            nii_6584_vshift=nii_6584_vshift,
            sii_6716_vshift=sii_6716_vshift,
            sii_6731_vshift=sii_6731_vshift,
            siii_9068_vshift=siii_9068_vshift,
            siii_9530_vshift=siii_9530_vshift,
            
            mgii_2800_sigma=mgii_2800_sigma,
            #mgii_2800b_sigma=mgii_2800b_sigma,
            nev_3346_sigma=nev_3346_sigma,
            nev_3426_sigma=nev_3426_sigma,
            oii_3726_sigma=oii_3726_sigma,
            oii_3729_sigma=oii_3729_sigma,
            neiii_3869_sigma=neiii_3869_sigma,
            oiii_4959_sigma=oiii_4959_sigma,
            oiii_5007_sigma=oiii_5007_sigma,
            hepsilon_sigma=hepsilon_sigma,
            hdelta_sigma=hdelta_sigma,
            hgamma_sigma=hgamma_sigma,
            hbeta_sigma=hbeta_sigma,
            halpha_sigma=halpha_sigma,
            nii_6548_sigma=nii_6548_sigma,
            nii_6584_sigma=nii_6584_sigma,
            sii_6716_sigma=sii_6716_sigma,
            sii_6731_sigma=sii_6731_sigma,
            siii_9068_sigma=siii_9068_sigma,
            siii_9530_sigma=siii_9530_sigma,
            
            **kwargs)

    def evaluate(self, emlinewave, *args):
        """Evaluate the emission-line model.
        emlineR=None, npixpercamera=None, 

        """ 
        from redrock.rebin import trapz_rebin

        linenames = np.array([linename.replace('_amp', '') for linename in self.param_names[:self.nline]])
        lineamps = np.hstack(args[0:self.nline])
        linevshifts = np.hstack(args[self.nline:2*self.nline])
        linesigmas = np.hstack(args[2*self.nline:])

        # build the emission-line model [erg/s/cm2/A, observed frame]
        log10model = np.zeros_like(self.log10wave)
        for linename, lineamp, linevshift, linesigma in zip(linenames, lineamps, linevshifts, linesigmas):
            
            iline = np.where(self.linetable['name'] == linename)[0]
            if len(iline) != 1:
                log.warning('No matching line found!')
                raise ValueError
                
            restlinewave = self.linetable[iline]['restwave'][0]
            
            linez = self.redshift + linevshift / C_LIGHT
            linezwave = np.log10(restlinewave * (1.0 + linez)) # redshifted wavelength [log-10 Angstrom]

            log10sigma = linesigma / C_LIGHT / np.log(10)      # line-width [log-10 Angstrom]
            
            ww = np.abs(self.log10wave - linezwave) < (100 * log10sigma)
            if np.count_nonzero(ww) > 0:
                #log.info(linename, 10**linezwave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
                log10model[ww] += lineamp * np.exp(-0.5 * (self.log10wave[ww]-linezwave)**2 / log10sigma**2)

            #print(linename, self.linetable['name'][iline][0], restlinewave, lineamp, linezwave, isbalmer)
            #pdb.set_trace()

        # split into cameras, resample, and convolve with the instrumental
        # resolution
        emlinemodel = []
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(self.npixpercamera[:ii+1])
            jpix = np.sum(self.npixpercamera[:ii+2])
            #_emlinemodel = resample_flux(emlinewave[ipix:jpix], 10**self.log10wave, log10model)
            _emlinemodel = trapz_rebin(10**self.log10wave, log10model, emlinewave[ipix:jpix])
            
            if self.emlineR is not None:
                _emlinemomdel = self.emlineR[ii].dot(_emlinemodel)
            
            #plt.plot(10**_emlinewave, _emlinemodel)
            #plt.plot(10**_emlinewave, self.emlineR[ii].dot(_emlinemodel))
            #plt.xlim(3870, 3920) ; plt.show()
            emlinemodel.append(_emlinemodel)
            
        return np.hstack(emlinemodel)

class EMLineFit(ContinuumTools):
    """Class to fit an emission-line spectrum.

    * https://docs.astropy.org/en/stable/modeling/example-fitting-constraints.html#tied
    * https://docs.astropy.org/en/stable/modeling/new-model.html
    * https://docs.astropy.org/en/stable/modeling/compound-models.html#parameters

    """
    def __init__(self, nball=10, chi2fail=1e6):
        """Write me.
        
        """
        from astropy.modeling import fitting

        super(EMLineFit, self).__init__()
        
        ## methods and attributes for synthesizing photometry
        #self.synth_bands = CFit.synth_bands
        #self.decam = CFit.decam
        #self.bassmzls = CFit.bassmzls
        #self.fluxnorm = CFit.fluxnorm
        #self.parse_photometry = CFit.parse_photometry

        self.nball = nball
        self.chi2fail = chi2fail
        self.pixkms = 10.0 # pixel size for internal wavelength array [km/s]

        self.fitter = fitting.LevMarLSQFitter()

    def init_output(self, linetable, nobj=1):
        """Initialize the output data table for this class.

        """
        import astropy.units as u
        from astropy.table import Table, Column
        
        out = Table()
        if False:
            out.add_column(Column(name='D4000_NOLINES', length=nobj, dtype='f4'))
            
        # observed-frame photometry synthesized from the spectra
        for band in self.synth_bands:
            out.add_column(Column(name='FLUX_SYNTH_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.nanomaggy)) 
            #out.add_column(Column(name='FLUX_SYNTH_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit=1/u.nanomaggy**2))
        # observed-frame photometry synthesized from the best-fitting continuum model fit
        for band in self.synth_bands:
            out.add_column(Column(name='FLUX_SYNTH_MODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit=u.nanomaggy))

        #out.add_column(Column(name='LINEVSHIFT_FORBIDDEN', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        #out.add_column(Column(name='LINEVSHIFT_FORBIDDEN_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        #out.add_column(Column(name='LINEVSHIFT_BALMER', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        #out.add_column(Column(name='LINEVSHIFT_BALMER_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        #out.add_column(Column(name='LINESIGMA_FORBIDDEN', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='LINESIGMA_FORBIDDEN_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        #out.add_column(Column(name='LINESIGMA_BALMER', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='LINESIGMA_BALMER_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))

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
            
            out.add_column(Column(name='{}_VSHIFT'.format(line), length=nobj, dtype='f4',
                                  unit=u.kilometer/u.second))
            out.add_column(Column(name='{}_VSHIFT_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=u.second**2 / u.kilometer**2))
            out.add_column(Column(name='{}_SIGMA'.format(line), length=nobj, dtype='f4',
                                  unit=u.kilometer / u.second))
            out.add_column(Column(name='{}_SIGMA_IVAR'.format(line), length=nobj, dtype='f4',
                                  unit=u.second**2 / u.kilometer**2))
            
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
        
    def chi2(self, bestfit, emlinewave, emlineflux, emlineivar):
        """Compute the reduced chi^2."""
        dof = len(emlinewave) - len(bestfit.parameters)
        emlinemodel = bestfit(emlinewave)
        chi2 = np.sum(emlineivar * (emlineflux - emlinemodel)**2) / dof
        return chi2

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (to be used to build QA).

        """
        npixpercamera = [len(gw) for gw in specwave]

        redshift = fastspecfit_table['CONTINUUM_Z']
        
        #linesigma_forbidden = fastspecfit_table['LINESIGMA_FORBIDDEN']
        #linesigma_balmer = fastspecfit_table['LINESIGMA_BALMER']
        #linevshift_forbidden = fastspecfit_table['LINEVSHIFT_FORBIDDEN']
        #linevshift_balmer = fastspecfit_table['LINEVSHIFT_BALMER']
        
        #linez_forbidden = fastspecfit_table['LINEZ_FORBIDDEN']
        #linez_balmer = fastspecfit_table['LINEZ_BALMER']
        #linevshift_forbidden = (linez_forbidden - redshift) * C_LIGHT # [km/s]
        #linevshift_balmer = (linez_balmer - redshift) * C_LIGHT # [km/s]

        EMLine = EMLineModel(
            #linevshift_forbidden=linevshift_forbidden,
            #linevshift_balmer=linevshift_balmer,
            #linesigma_forbidden=linesigma_forbidden,
            #linesigma_balmer=linesigma_balmer,
            redshift=redshift, emlineR=specres,
            npixpercamera=npixpercamera)
        
        lineargs = [fastspecfit_table[linename.upper()] for linename in EMLine.param_names]#[4:]]
        
        # skip linevshift_[forbidden,balmer] and linesigma_[forbidden,balmer]
        #lineargs = [fastspecfit_table[linename.upper()] for linename in EMLine.param_names[4:]] 
        #lineargs = [linevshift_forbidden, linevshift_balmer, linesigma_forbidden, linesigma_balmer] + lineargs

        _emlinemodel = EMLine.evaluate(np.hstack(specwave), *lineargs)

        # unpack it
        emlinemodel = []
        npix = np.hstack([0, npixpercamera])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            emlinemodel.append(_emlinemodel[ipix:jpix])

        return emlinemodel
    
    def fit(self, data, continuummodel):
        """Perform the fit minimization / chi2 minimization.
        
        EMLineModel object
        FC - ContinuumFit object

        need to take into account the instrumental velocity width when computing integrated fluxes
        
        """
        #from scipy import integrate
        from astropy.stats import sigma_clipped_stats
        
        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        npixpercamera = [len(gw) for gw in data['wave']]
        npixpercam = np.hstack([0, npixpercamera])

        redshift = data['zredrock']
        emlinewave = np.hstack(data['wave'])
        emlineivar = np.hstack(data['ivar'])
        specflux = np.hstack(data['flux'])
        continuummodelflux = np.hstack(continuummodel)
        emlineflux = specflux - continuummodelflux

        dlogwave = self.pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        log10wave = np.arange(np.log10(3e3), np.log10(1e4), dlogwave)
        #log10wave = np.arange(np.log10(emlinewave.min()), np.log10(emlinewave.max()), dlogwave)
        
        self.EMLineModel = EMLineModel(redshift=redshift,
                                       emlineR=data['res'],
                                       npixpercamera=npixpercamera,
                                       log10wave=log10wave)
        nparam = len(self.EMLineModel.parameters)
        #params = np.repeat(self.EMLineModel.parameters, self.nball).reshape(nparam, self.nball)

        # do a fast box-car integration to get the initial line-amplitudes
        sigma_cont = 200.0
        for pp in self.EMLineModel.param_names:
            if getattr(self.EMLineModel, pp).tied:
                print('Skipping {}'.format(pp))
                continue

            if 'amp' in pp:
                pinfo = getattr(self.EMLineModel, pp)
                oneline = self.linetable[self.linetable['name'] == pinfo.name.replace('_amp', '')][0]
                zwave = oneline['restwave'] * (1 + redshift)
                lineindx = np.where((emlinewave > (zwave - 3*sigma_cont * zwave / C_LIGHT)) *
                                    (emlinewave < (zwave + 3.*sigma_cont * zwave / C_LIGHT)) *
                                    (emlineivar > 0))[0]
                if len(lineindx) > 10:
                    lineflux = np.sum(emlineflux[lineindx])
                    #lineflux = np.sum(emlineivar[lineindx] * emlineflux[lineindx]) / np.sum(emlineivar[lineindx])
                    linesigma_ang = zwave * sigma_cont / C_LIGHT # [observed-frame Angstrom]
                    linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang
                    lineamp = np.abs(lineflux / linenorm)

                    #if 'beta' in pinfo.name:
                    #    pdb.set_trace()
                    
                    if not pinfo.tied:
                        setattr(self.EMLineModel, pp, lineamp)
                    print(pinfo.name, len(lineindx), lineflux, lineamp)
                    
        ## reset the tied amplitudes
        #self.EMLineModel.oiii_4959_amp = self.EMLineModel.oiii_4959_amp.tied(self.EMLineModel)
        #self.EMLineModel.nii_6548_amp = self.EMLineModel.nii_6548_amp.tied(self.EMLineModel)
        #self.EMLineModel.linevshift_balmer.fixed = True
        #self.EMLineModel.linevshift_forbidden.fixed = True

        # Fit [1]: tie all lines together
        self.EMLineModel = _tie_all_lines(self.EMLineModel)
        bestfit = self.fitter(self.EMLineModel, emlinewave, emlineflux, weights=np.sqrt(emlineivar), maxiter=100)
        #print(bestfit.parameters)

        # Fit [2]: tie Balmer, narrow forbidden, and QSO/broad lines together, separately
        self.EMLineModel = bestfit
        self.EMLineModel = _tie_balmer_lines(self.EMLineModel)
        self.EMLineModel = _tie_forbidden_lines(self.EMLineModel)
        #self.EMLineModel = _tie_qso_lines(self.EMLineModel)
        bestfit = self.fitter(self.EMLineModel, emlinewave, emlineflux, weights=np.sqrt(emlineivar), maxiter=100)
        
        emlinemodel = bestfit(emlinewave)
        chi2 = self.chi2(bestfit, emlinewave, emlineflux, emlineivar).astype('f4')

        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_output(self.EMLineModel.linetable)
        
        # Synthesize photometry from the best-fitting model (continuum+emission lines).
        if data['photsys'] == 'S':
            filters = self.decam
        else:
            filters = self.bassmzls

        # The wavelengths overlap between the cameras a bit...
        srt = np.argsort(emlinewave)
        padflux, padwave = filters.pad_spectrum((continuummodelflux+emlinemodel)[srt], emlinewave[srt], method='edge')
        synthmaggies = filters.get_ab_maggies(padflux / self.fluxnorm, padwave)
        synthmaggies = synthmaggies.as_array().view('f8')
        model_synthphot = self.parse_photometry(self.synth_bands, maggies=synthmaggies,
                                                nanomaggies=False,
                                                lambda_eff=filters.effective_wavelengths.value)

        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_{}'.format(band.upper())] = data['synthphot']['nanomaggies'][iband]
            #result['FLUX_SYNTH_IVAR_{}'.format(band.upper())] = data['synthphot']['nanomaggies_ivar'][iband]
        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_MODEL_{}'.format(band.upper())] = model_synthphot['nanomaggies'][iband]

        specflux_nolines = specflux - emlinemodel
        
        # measure D(4000) without the emission lines
        if False:
            d4000_nolines, _ = self.get_d4000(emlinewave, specflux_nolines, redshift=redshift)
            result['D4000_NOLINES'] = d4000_nolines

        ## Pack the results in a dictionary and return.
        ## https://gist.github.com/eteq/1f3f0cec9e4f27536d52cd59054c55f2
        #result = {
        #    'converged': False,
        #    'fit_message': self.fitter.fit_info['message'],
        #    'nparam': nparam,
        #    'npix': len(emlinewave),
        #    'dof': len(emlinewave) - len(self.EMLineModel.parameters),
        #    'chi2': chi2,
        #    'linenames': [ll.replace('_amp', '') for ll in self.EMLineModel.param_names[4:]],
        #}
        #for param in bestfit.param_names:
        #    result.update({param: getattr(bestfit, param).value})
        
        # uncertainties
        if self.fitter.fit_info['param_cov'] is not None:
            cov = self.fitter.fit_info['param_cov']
            ivar = 1 / np.diag(cov)
            #result['converged'] = True
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
            #result.update({pinfo.name: pinfo.value.astype('f4')})
            result[pinfo.name.upper()] = pinfo.value.astype('f4')
                
            if pinfo.fixed:
                #result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
                result['{}_IVAR'.format(pinfo.name.upper())] = pinfo.value.astype('f4')
            elif pinfo.tied:
                # hack! see https://github.com/astropy/astropy/issues/7202
                #result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
                result['{}_IVAR'.format(pinfo.name.upper())] = np.float32(0.0)
            else:
                result['{}_IVAR'.format(pinfo.name.upper())] = ivar[count].astype('f4')
                #result.update({'{}_ivar'.format(pinfo.name): ivar[count].astype('f4')})
                count += 1

            #if 'forbidden' in pinfo.name:
            #    pdb.set_trace()

        # hack for tied parameters---gotta be a better way to do this
        #result['oiii_4959_amp_ivar'] = result['oiii_5007_amp_ivar'] * 2.8875**2
        result['OIII_4959_AMP_IVAR'] = result['OIII_5007_AMP_IVAR'] * 2.8875**2
        result['NII_6548_AMP_IVAR'] = result['NII_6548_AMP_IVAR'] * 2.936**2

        #if False:
        #    result['LINEVSHIFT_FORBIDDEN_IVAR'] = result['LINEVSHIFT_BALMER_IVAR']
        #    result['LINESIGMA_FORBIDDEN_IVAR'] = result['LINESIGMA_BALMER_IVAR']

        ## convert the vshifts to redshifts
        #result['linez_forbidden'] = redshift + result['linevshift_forbidden'] / C_LIGHT
        #result['linez_balmer'] = redshift + result['linevshift_balmer'] / C_LIGHT
        #result['linez_forbidden_ivar'] = result['linevshift_forbidden_ivar'] * C_LIGHT**2
        #result['linez_balmer_ivar'] = result['linevshift_balmer_ivar'] * C_LIGHT**2

        # now loop back through and if ivar==0 then set the parameter value to zero
        if False:
            if self.fitter.fit_info['param_cov'] is not None:
                for pp in bestfit.param_names:
                    if result['{}_IVAR'.format(pp.upper())] == 0.0:
                        result[pp] = np.float(0.0)

        # get continuum fluxes, EWs, and upper limits
        sigma_cont = 150.0
        for oneline in self.EMLineModel.linetable:
            line = oneline['name'].upper()

            # get the emission-line flux
            
            #zwave = oneline['restwave'] * (1 + redshift)
            #if oneline['isbalmer']:
            #    linesigma = result['LINESIGMA_BALMER']
            #else:
            #    linesigma = result['LINESIGMA_FORBIDDEN']

            linez = redshift + result['{}_VSHIFT'.format(line)] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)
            
            linesigma = result['{}_SIGMA'.format(line)]     # [km/s]
            linesigma_ang = linezwave * linesigma / C_LIGHT # [observed-frame Angstrom]
            norm = np.sqrt(2.0 * np.pi) * linesigma_ang

            if norm == 0.0:
                pdb.set_trace()

            result['{}_FLUX'.format(line)] = result['{}_AMP'.format(line)] * norm
            result['{}_FLUX_IVAR'.format(line)] = result['{}_AMP_IVAR'.format(line)] / norm**2

            # boxcar integration, chi2, and number of pixels
            lineindx = np.where((emlinewave > (linezwave - 3.*linesigma * linezwave / C_LIGHT)) *
                                (emlinewave < (linezwave + 3.*linesigma * linezwave / C_LIGHT)) *
                                (emlineivar > 0))[0]
            npix = len(lineindx)
            if npix > 0:
                dof = npix - 3 # ??? [redshift, sigma, and amplitude]
                chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-emlinemodel[lineindx])**2) / dof
                boxflux = np.sum(emlineflux[lineindx])
                boxflux_ivar = 1 / np.sum(1 / emlineivar[lineindx])
            else:
                npix, chi2, boxflux, boxflux_ivar = 0.0, 1e6, 0.0, 0.0

            result['{}_NPIX'.format(line)] = npix
            result['{}_CHI2'.format(line)] = chi2
            result['{}_BOXFLUX'.format(line)] = boxflux
            result['{}_BOXFLUX_IVAR'.format(line)] = boxflux_ivar

            #print(line, npix, result['{}_FLUX'.format(line)][0], result['{}_BOXFLUX'.format(line)][0], result['{}_AMP'.format(line)][0])
            print(line, npix, result['{}_SIGMA'.format(line)][0], result['{}_VSHIFT'.format(line)][0],
                  result['{}_FLUX'.format(line)][0], result['{}_AMP'.format(line)][0])
            
            # get the continuum and EWs
            indxlo = np.where((emlinewave > (linezwave - 10*linesigma * linezwave / C_LIGHT)) *
                              (emlinewave < (linezwave - 3.*linesigma * linezwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indxhi = np.where((emlinewave < (zwave + 10*linesigma * linezwave / C_LIGHT)) *
                              (emlinewave > (zwave + 3.*linesigma * linezwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indx = np.hstack((indxlo, indxhi))

            if len(indx) > 5: # require at least 5 pixels to get the continuum level
                _, cmed, csig = sigma_clipped_stats(specflux_nolines[indx], sigma=3.0)
                civar = (np.sqrt(len(indx)) / csig)**2
            else:
                cmed, civar = 0.0, 0.0
                
            result['{}_CONT'.format(line)] = cmed
            result['{}_CONT_IVAR'.format(line)] = civar

            if result['{}_CONT_IVAR'.format(line)] != 0.0:
                #if result['{}_CONT'.format(line)] == 0:
                #    pdb.set_trace()
                factor = (1 + redshift) / result['{}_CONT'.format(line)] # --> rest frame
                ew = result['{}_FLUX'.format(line)] * factor   # rest frame [A]
                ewivar = result['{}_FLUX_IVAR'.format(line)] / factor**2

                # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar)
                ewlimit = fluxlimit * factor
            else:
                ew, ewivar, fluxlimit, ewlimit = 0.0, 0.0, 0.0, 0.0

            result['{}_EW'.format(line)] = ew
            result['{}_EW_IVAR'.format(line)] = ewivar
            result['{}_FLUX_LIMIT'.format(line)] = fluxlimit
            result['{}_EW_LIMIT'.format(line)] = ewlimit

            # simple QA
            if 'alpha' in line and False:
                import matplotlib.pyplot as plt
                _indx = np.where((emlinewave > (zwave - 15*sigma_cont * zwave / C_LIGHT)) *
                                (emlinewave < (zwave + 15*sigma_cont * zwave / C_LIGHT)))[0]
                plt.plot(emlinewave[_indx], emlineflux[_indx])
                plt.plot(emlinewave[_indx], specflux_nolines[_indx])
                plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
                plt.axhline(y=cmed, color='k')
                plt.axhline(y=cmed+csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.axhline(y=cmed-csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.savefig('junk.png')

        return result
    
    def qa_fastspec(self, data, fastspec, metadata, coadd_type='deep',
                    outprefix=None, outdir=None):
        """QA plot the emission-line spectrum and best-fitting model.

        """
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        from fastspecfit.util import ivar2var

        sns.set(context='talk', style='ticks', font_scale=1.5)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]

        if outdir is None:
            outdir = '.'
        if outprefix is None:
            outprefix = 'fastspec'

        if coadd_type == 'deep' or coadd_type == 'all':
            title = 'Tile/Coadd: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], coadd_type, metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], coadd_type, metadata['TARGETID']))
        elif coadd_type == 'night':
            title = 'Tile/Night: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID'],
                    metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'], metadata['TARGETID']))
        else:
            title = 'Tile/Night/Expid: {}/{}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['NIGHT'], metadata['EXPID'],
                    metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['NIGHT'],
                    metadata['EXPID'], metadata['TARGETID']))

        redshift = fastspec['CONTINUUM_Z']

        # rebuild the best-fitting spectroscopic and photometric models
        continuum, _ = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift, 
                                     specwave=data['wave'], specres=data['res'],
                                     AV=fastspec['CONTINUUM_AV'],
                                     vdisp=fastspec['CONTINUUM_VDISP'],
                                     coeff=fastspec['CONTINUUM_COEFF'],
                                     synthphot=False)

        _emlinemodel = self.emlinemodel_bestfit(data['wave'], data['res'], fastspec)

        inches_wide = 16
        inches_fullspec = 6
        inches_perline = inches_fullspec / 2.0
        nlinepanels = 4

        nline = len(set(self.linetable['plotgroup']))
        nlinerows = np.ceil(nline / nlinepanels).astype(int)
        nrows = 2 + nlinerows

        height_ratios = np.hstack([1, 1, [0.5]*nlinerows])

        fig = plt.figure(figsize=(inches_wide, 2*inches_fullspec + inches_perline*nlinerows))
        gs = fig.add_gridspec(nrows, nlinepanels, height_ratios=height_ratios)

        # full spectrum + best-fitting continuum model
        bigax1 = fig.add_subplot(gs[0, :])

        leg = {
            'targetid': '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastspec['CONTINUUM_CHI2']),
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            #'zfastfastspec': '$z_{{\\rm fastspecfit}}$={:.6f}'.format(fastspec['CONTINUUM_Z']),
            #'z': '$z$={:.6f}'.format(fastspec['CONTINUUM_Z']),
            'age': '<Age>={:.3f} Gyr'.format(fastspec['CONTINUUM_AGE']),
            }
        if fastspec['CONTINUUM_VDISP_IVAR'] == 0:
            leg.update({'vdisp': '$\sigma_{{\\rm star}}$={:.1f} km/s'.format(fastspec['CONTINUUM_VDISP'])})
        else:
            leg.update({'vdisp': '$\sigma_{{\\rm star}}$={:.1f}+/-{:.1f} km/s'.format(
                fastspec['CONTINUUM_VDISP'], 1/np.sqrt(fastspec['CONTINUUM_VDISP_IVAR']))})
            
        if fastspec['CONTINUUM_AV_IVAR'] == 0:
            leg.update({'AV': '$A(V)$={:.3f} mag'.format(fastspec['CONTINUUM_AV'])})
        else:
            leg.update({'AV': '$A(V)$={:.3f}+/-{:.3f} mag'.format(
                fastspec['CONTINUUM_AV'], 1/np.sqrt(fastspec['CONTINUUM_AV_IVAR']))})

        ymin, ymax = 1e6, -1e6
        wavelims = (3600, 9800)
        
        for ii in [0, 1, 2]: # iterate over cameras
            sigma, _ = ivar2var(data['ivar'][ii], sigma=True)

            bigax1.fill_between(data['wave'][ii], data['flux'][ii]-sigma,
                             data['flux'][ii]+sigma, color=col1[ii])
            bigax1.plot(data['wave'][ii], continuum[ii], color=col2[ii], alpha=1.0)#, color='k')
            #bigax1.plot(data['wave'][ii], continuum_nodust[ii], alpha=0.5, color='k')

            # get the robust range
            filtflux = median_filter(data['flux'][ii], 5)
            if np.min(filtflux) < ymin:
                ymin = np.min(filtflux) * 0.5
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux) * 1.7

        fntsz = 18
        #bigax1.text(0.98, 0.92, '{}'.format(leg['targetid']), 
        #         ha='right', va='center', transform=bigax1.transAxes, fontsize=fntsz)
        #bigax1.text(0.98, 0.92, '{} {}'.format(leg['targetid'], leg['zredrock']), 
        #         ha='right', va='center', transform=bigax1.transAxes, fontsize=fntsz)
        bigax1.text(0.98, 0.92, r'{}'.format(leg['zredrock']),
                 ha='right', va='center', transform=bigax1.transAxes, fontsize=fntsz)
        bigax1.text(0.98, 0.83, r'{} {}'.format(leg['chi2'], leg['age']),
                 ha='right', va='center', transform=bigax1.transAxes, fontsize=fntsz)
        bigax1.text(0.98, 0.74, r'{}'.format(leg['AV']),
                 ha='right', va='center', transform=bigax1.transAxes, fontsize=fntsz)
        bigax1.text(0.98, 0.65, r'{}'.format(leg['vdisp']),
                 ha='right', va='center', transform=bigax1.transAxes, fontsize=fntsz)
        
        bigax1.text(0.03, 0.9, 'Observed Spectrum + Continuum Model',
                    ha='left', va='center', transform=bigax1.transAxes,
                    fontsize=20)
        bigax1.set_title(title, fontsize=20)
            
        bigax1.set_xlim(wavelims)
        bigax1.set_ylim(ymin, ymax)
        bigax1.set_xticklabels([])
        #bigax1.set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
        #bigax1.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax1.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 

        # full emission-line spectrum + best-fitting lines
        bigax2 = fig.add_subplot(gs[1, :])

        if False:
            leg = {
                #'targetid': '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),
                #'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
                'linevshift_forbidden': '$\Delta\,v_{{\\rm forb}}$={:.1f} km/s'.format(fastspec['LINEVSHIFT_FORBIDDEN']),
                'linevshift_balmer': '$\Delta\,v_{{\\rm Balm}}$={:.1f} km/s'.format(fastspec['LINEVSHIFT_BALMER']),
                'linesigma_forbidden': '$\sigma_{{\\rm forb}}$={:.1f} km/s'.format(fastspec['LINESIGMA_FORBIDDEN']),
                'linesigma_balmer': '$\sigma_{{\\rm Balm}}$={:.1f} km/s'.format(fastspec['LINESIGMA_BALMER']),
                'linevshift': '$\Delta_{{\\rm line}}\,v$={:.1f} km/s'.format(fastspec['LINEVSHIFT_BALMER']),
                'linesigma': '$\sigma_{{\\rm line}}$={:.1f} km/s'.format(fastspec['LINESIGMA_BALMER']),
                }

        ymin, ymax = 1e6, -1e6
        for ii in [0, 1, 2]: # iterate over cameras
            emlinewave = data['wave'][ii]
            emlineflux = data['flux'][ii] - continuum[ii]
            emlinemodel = _emlinemodel[ii]

            emlinesigma, good = ivar2var(data['ivar'][ii], sigma=True)
            emlinewave = emlinewave[good]
            emlineflux = emlineflux[good]
            emlinesigma = emlinesigma[good]
            emlinemodel = emlinemodel[good]

            bigax2.fill_between(emlinewave, emlineflux-emlinesigma,
                               emlineflux+emlinesigma, color=col1[ii], alpha=0.7)
            bigax2.plot(emlinewave, emlinemodel, color=col2[ii], lw=2)

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

        fntsz = 18
        #bigax2.text(0.98, 0.92, '{}'.format(leg['targetid']), 
        #           ha='right', va='center', transform=bigax2.transAxes, fontsize=fntsz)
        #bigax2.text(0.98, 0.86, r'{}'.format(leg['zredrock']),
        #           ha='right', va='center', transform=bigax2.transAxes, fontsize=fntsz)
        
        #bigax2.text(0.98, 0.92, r'{}'.format(leg['linesigma']),
        #           ha='right', va='center', transform=bigax2.transAxes, fontsize=fntsz)
        #bigax2.text(0.98, 0.86, r'{}'.format(leg['linevshift']),
        #           ha='right', va='center', transform=bigax2.transAxes, fontsize=fntsz)
        
        #bigax2.text(0.98, 0.92, r'{} {}'.format(leg['linevshift_balmer'], leg['linevshift_forbidden']),
        #           ha='right', va='center', transform=bigax2.transAxes, fontsize=fntsz)
        #bigax2.text(0.98, 0.86, r'{} {}'.format(leg['linesigma_balmer'], leg['linesigma_forbidden']),
        #           ha='right', va='center', transform=bigax2.transAxes, fontsize=fntsz)
        
        bigax2.text(0.03, 0.9, 'Residual Spectrum + Emission-Line Model',
                    ha='left', va='center', transform=bigax2.transAxes,
                    fontsize=20)
                
        bigax2.set_xlim(wavelims)
        bigax2.set_ylim(ymin, ymax)
        
        #bigax2.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax2.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 
        
        # zoom in on individual emission lines - use linetable!
        sig = 500.0 # [km/s]

        meanwaves, deltawaves, sigmas, linenames = [], [], [], []
        for plotgroup in set(self.linetable['plotgroup']):
            I = np.where(plotgroup == self.linetable['plotgroup'])[0]
            linenames.append(self.linetable['nicename'][I[0]])
            meanwaves.append(np.mean(self.linetable['restwave'][I]))
            deltawaves.append((np.max(self.linetable['restwave'][I]) - np.min(self.linetable['restwave'][I])) / 2)

            #if np.any(self.linetable['isbalmer'][I]):
            #    sigmas.append(np.max((fastspec['LINESIGMA_BALMER'], fastspec['LINESIGMA_FORBIDDEN'])))
            #else:
            #    sigmas.append(fastspec['LINESIGMA_FORBIDDEN'])
            sigmas.append(sig)

        srt = np.argsort(meanwaves)
        meanwaves = np.hstack(meanwaves)[srt]
        deltawaves = np.hstack(deltawaves)[srt]
        sigmas = np.hstack(sigmas)[srt]
        linenames = np.hstack(linenames)[srt]

        removelabels = np.ones(nline, bool)
        ymin, ymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
        
        ax, irow, icol = [], 2, 0
        for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):
            icol = iax % nlinepanels
            if iax > 0 and iax % nlinepanels == 0:
                irow += 1
            #print(iax, irow, icol)

            xx = fig.add_subplot(gs[irow, icol])
            ax.append(xx)

            wmin = (meanwave - deltawave) * (1+redshift) - 8 * sig * meanwave / C_LIGHT
            wmax = (meanwave + deltawave) * (1+redshift) + 8 * sig * meanwave / C_LIGHT
            #print(linename, wmin, wmax)

            # iterate over cameras
            for ii in [0, 1, 2]: # iterate over cameras
                emlinewave = data['wave'][ii]
                emlineflux = data['flux'][ii] - continuum[ii]
                emlinemodel = _emlinemodel[ii]

                emlinesigma, good = ivar2var(data['ivar'][ii], sigma=True)
                emlinewave = emlinewave[good]
                emlineflux = emlineflux[good]
                emlinesigma = emlinesigma[good]
                emlinemodel = emlinemodel[good]

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
                xx.xaxis.set_major_locator(ticker.MaxNLocator(2))
                #xx.xaxis.set_major_locator(ticker.MultipleLocator(20)) # wavelength spacing of ticks [Angstrom]
                #if iax == 2:
                #    pdb.set_trace()

        # common axis labels
        tp, bt, lf, rt = 0.95, 0.09, 0.10, 0.94
        
        fig.text(lf-0.06, (tp-bt)/2+bt,
                 r'Flux Density ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
                 ha='center', va='center', rotation='vertical')
        fig.text((rt-lf)/2+lf, bt-0.05, r'Observed-frame Wavelength ($\AA$)',
                 ha='center', va='center')
            
        plt.subplots_adjust(wspace=0.27, top=tp, bottom=bt, left=lf, right=rt, hspace=0.22)

        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()
