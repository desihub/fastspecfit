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
def _tie_siii_9532_sigma(model):
    return model.siii_9069_sigma

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
    
    return model

def _tie_broad_lines(model):
    """For now, the only *broad* line is MgII, so let it float. However, it's pretty
    fragile that the line-name is hard-coded...

    """
    #for pp in model.param_names:
    #    if 'sigma' in pp and 'mgii' in pp:
    #        getattr(model, pp).tied = _tie_mgii_2800_sigma
    #    if 'vshift' in pp and 'mgii' in pp:
    #        getattr(model, pp).tied = _tie_mgii_2800_vshift
    model.mgii_2800_sigma.tied = False
    model.mgii_2800_vshift.tied = False

    # update / improve the initial guess
    model.mgii_2800_sigma.value = model.initsigma_broad
    model.mgii_2800_sigma._default = model.initsigma_broad
    model.mgii_2800_sigma.bounds = [0.0, model.maxsigma_broad]
    
    return model

def _tie_all_lines(model):
    for pp in model.param_names:
        if 'sigma' in pp and pp != 'hbeta_sigma':
            getattr(model, pp).tied = _tie_hbeta_sigma
        if 'vshift' in pp and pp != 'hbeta_vshift':
            getattr(model, pp).tied = _tie_hbeta_vshift
            
    # reset the amplitudes of the tied doublets; fragile...
    model.oiii_4959_amp = model.oiii_4959_amp.tied(model)
    model.nii_6548_amp = model.nii_6548_amp.tied(model)
            
    return model

class EMLineModel(Fittable1DModel):
    """Class to model the emission-line spectra.

    """
    from astropy.modeling import Parameter

    # NB! The order of the parameters here matters!
    vmaxshift = 300.0
    initvshift = 0.0

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
    siii_9069_amp = Parameter(name='siii_9069_amp', default=0.3)
    siii_9532_amp = Parameter(name='siii_9532_amp', default=0.3)

    mgii_2800_vshift = Parameter(name='mgii_2800_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    #mgii_2800b_vshift = Parameter(name='mgii_2800b_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    nev_3346_vshift = Parameter(name='nev_3346_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    nev_3426_vshift = Parameter(name='nev_3426_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    oii_3726_vshift = Parameter(name='oii_3726_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    oii_3729_vshift = Parameter(name='oii_3729_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    neiii_3869_vshift = Parameter(name='neiii_3869_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    oiii_4959_vshift = Parameter(name='oiii_4959_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    oiii_5007_vshift = Parameter(name='oiii_5007_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    hepsilon_vshift = Parameter(name='hepsilon_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    hdelta_vshift = Parameter(name='hdelta_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    hgamma_vshift = Parameter(name='hgamma_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    hbeta_vshift = Parameter(name='hbeta_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    halpha_vshift = Parameter(name='halpha_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    nii_6548_vshift = Parameter(name='nii_6548_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    nii_6584_vshift = Parameter(name='nii_6584_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    sii_6716_vshift = Parameter(name='sii_6716_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    sii_6731_vshift = Parameter(name='sii_6731_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    siii_9069_vshift = Parameter(name='siii_9069_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])
    siii_9532_vshift = Parameter(name='siii_9532_vshift', default=initvshift, bounds=[-vmaxshift, +vmaxshift])

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
    siii_9069_sigma = Parameter(name='siii_9069_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])
    siii_9532_sigma = Parameter(name='siii_9532_sigma', default=initsigma_narrow, bounds=[0.0, maxsigma_narrow])

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
                 siii_9069_amp=siii_9069_amp.default,
                 siii_9532_amp=siii_9532_amp.default,
                 
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
                 siii_9069_vshift=siii_9069_vshift.default,
                 siii_9532_vshift=siii_9532_vshift.default,
                 
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
                 siii_9069_sigma=siii_9069_sigma.default,
                 siii_9532_sigma=siii_9532_sigma.default,
                 
                 redshift=None,
                 emlineR=None,
                 npixpercamera=None,
                 #goodpixpercam=None,
                 log10wave=None,
                 **kwargs):
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
        #self.goodpixpercam = goodpixpercam

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
            siii_9069_amp=siii_9069_amp,
            siii_9532_amp=siii_9532_amp,
            
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
            siii_9069_vshift=siii_9069_vshift,
            siii_9532_vshift=siii_9532_vshift,
            
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
            siii_9069_sigma=siii_9069_sigma,
            siii_9532_sigma=siii_9532_sigma,
            
            **kwargs)

    def evaluate(self, emlinewave, *args):
        """Evaluate the emission-line model.

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

        # split into cameras, resample, and convolve with the instrumental
        # resolution
        emlinemodel = []
        for ii in [0, 1, 2]: # iterate over cameras
            #ipix = np.sum(self.ngoodpixpercamera[:ii+1]) # unmasked pixels!
            #jpix = np.sum(self.ngoodpixpercamera[:ii+2])

            ipix = np.sum(self.npixpercamera[:ii+1]) # all pixels!
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
        
        self.nball = nball
        self.chi2fail = chi2fail
        self.pixkms = 10.0 # pixel size for internal wavelength array [km/s]

        self.fitter = fitting.LevMarLSQFitter()#calc_uncertainties=True)
        #self.fitter_nouncertainties = fitting.LevMarLSQFitter(calc_uncertainties=False)

    def init_output(self, linetable, nobj=1, chi2_default=1e6):
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

        out.add_column(Column(name='BALMER_Z', length=nobj, dtype='f8'))
        #out.add_column(Column(name='BALMER_Z_ERR', length=nobj, dtype='f8'))
        out.add_column(Column(name='FORBIDDEN_Z', length=nobj, dtype='f8'))
        #out.add_column(Column(name='FORBIDDEN_Z_ERR', length=nobj, dtype='f8'))
        out.add_column(Column(name='BROAD_Z', length=nobj, dtype='f8'))
        #out.add_column(Column(name='BROAD_Z_ERR', length=nobj, dtype='f8'))

        out.add_column(Column(name='BALMER_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='BALMER_SIGMA_ERR', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='FORBIDDEN_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='FORBIDDEN_SIGMA_ERR', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        out.add_column(Column(name='BROAD_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='BROAD_SIGMA_ERR', length=nobj, dtype='f4', unit=u.kilometer / u.second))

        #out.add_column(Column(name='BALMER_VSHIFT', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        #out.add_column(Column(name='BALMER_VSHIFT_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        #out.add_column(Column(name='BALMER_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='BALMER_SIGMA_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        #out.add_column(Column(name='FORBIDDEN_VSHIFT', length=nobj, dtype='f4', unit=u.kilometer/u.second))
        #out.add_column(Column(name='FORBIDDEN_VSHIFT_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))
        #out.add_column(Column(name='FORBIDDEN_SIGMA', length=nobj, dtype='f4', unit=u.kilometer / u.second))
        #out.add_column(Column(name='FORBIDDEN_SIGMA_IVAR', length=nobj, dtype='f4', unit=u.second**2 / u.kilometer**2))

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
            #out.add_column(Column(name='{}_BOXFLUX_IVAR'.format(line), length=nobj, dtype='f4',
            #                      unit=10**34*u.second**2*u.cm**4/u.erg**2))
            
            out.add_column(Column(name='{}_VSHIFT'.format(line), length=nobj, dtype='f4',
                                  unit=u.kilometer/u.second))
            #out.add_column(Column(name='{}_VSHIFT_IVAR'.format(line), length=nobj, dtype='f4',
            #                      unit=u.second**2 / u.kilometer**2))
            out.add_column(Column(name='{}_SIGMA'.format(line), length=nobj, dtype='f4',
                                  unit=u.kilometer / u.second))
            #out.add_column(Column(name='{}_SIGMA_IVAR'.format(line), length=nobj, dtype='f4',
            #                      unit=u.second**2 / u.kilometer**2))
            
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
            out.add_column(Column(name='{}_CHI2'.format(line), data=np.repeat(chi2_default, nobj), dtype='f4'))
            out.add_column(Column(name='{}_NPIX'.format(line), length=nobj, dtype=np.int32))

        return out
        
    def chi2(self, bestfit, emlinewave, emlineflux, emlineivar):
        """Compute the reduced chi^2."""
        dof = np.sum(emlineivar > 0) - len(bestfit.parameters)
        if dof > 0:
            emlinemodel = bestfit(emlinewave)
            chi2 = np.sum(emlineivar * (emlineflux - emlinemodel)**2) / dof
        else:
            chi2 = 1e6
        return chi2

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (to be used to build QA).

        """
        npixpercamera = [len(gw) for gw in specwave]
        redshift = fastspecfit_table['CONTINUUM_Z']
        
        EMLine = EMLineModel(
            redshift=redshift, emlineR=specres,
            npixpercamera=npixpercamera)
        
        lineargs = [fastspecfit_table[linename.upper()] for linename in EMLine.param_names]#[4:]]
        
        _emlinemodel = EMLine.evaluate(np.hstack(specwave), *lineargs)

        # unpack the per-camera spectra
        emlinemodel = []
        npix = np.hstack([0, npixpercamera])
        for ii in [0, 1, 2]: # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            emlinemodel.append(_emlinemodel[ipix:jpix])

        return emlinemodel
    
    def fit(self, data, continuummodel, smooth_continuum):
        """Perform the fit minimization / chi2 minimization.
        
        EMLineModel object
        FC - ContinuumFit object

        ToDo: need to take into account the instrumental velocity width when
        computing integrated fluxes...

        """
        from fastspecfit.util import ivar2var
        from scipy.stats import sigmaclip
        #from astropy.stats import sigma_clipped_stats
        from scipy.ndimage.filters import median_filter

        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        redshift = data['zredrock']
        emlinewave = np.hstack(data['wave'])
        emlineivar = np.hstack(data['ivar'])
        specflux = np.hstack(data['flux'])
        continuummodelflux = np.hstack(continuummodel)
        smoothcontinuummodelflux = np.hstack(smooth_continuum)

        emlineflux = specflux - continuummodelflux - smoothcontinuummodelflux

        emlinevar, emlinegood = ivar2var(emlineivar)

        npixpercamera = [len(gw) for gw in data['wave']] # all pixels
        npixpercam = np.hstack([0, npixpercamera])
        #goodpixpercam = [np.where(iv > 0)[0] for iv in data['ivar']] # unmasked pixels used in the fitting

        #ngoodpixpercamera = [np.count_nonzero(iv) for iv in data['ivar']] # unmasked pixels used in the fitting
        #ngoodpixpercam = np.hstack([0, ngoodpixpercamera])

        dlogwave = self.pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        log10wave = np.arange(np.log10(3e3), np.log10(1e4), dlogwave)
        #log10wave = np.arange(np.log10(emlinewave.min()), np.log10(emlinewave.max()), dlogwave)
        
        self.EMLineModel = EMLineModel(redshift=redshift,
                                       emlineR=data['res'],
                                       npixpercamera=npixpercamera,
                                       #goodpixpercam=goodpixpercam,
                                       log10wave=log10wave)
        nparam = len(self.EMLineModel.parameters)
        #params = np.repeat(self.EMLineModel.parameters, self.nball).reshape(nparam, self.nball)

        # Do a fast box-car integration to get the initial line-amplitudes and
        # line-widths...actually these initial guesses are not really working
        # but keep the code here.
        if False:
            sigma_cont = 200.0
            init_linesigmas = []
            for pp in self.EMLineModel.param_names:
                if getattr(self.EMLineModel, pp).tied:
                    #print('Skipping tied parameter {}'.format(pp))
                    continue

                if 'amp' in pp:
                    pinfo = getattr(self.EMLineModel, pp)
                    linename = pinfo.name.replace('_amp', '')

                    iline = np.where(self.linetable['name'] == linename)[0]
                    if len(iline) != 1:
                        log.warning('No matching line found!')
                        raise ValueError

                    oneline = self.linetable[iline][0]
                    linezwave = oneline['restwave'] * (1 + redshift)
                    lineindx = np.where((emlinewave > (linezwave - 5*sigma_cont * linezwave / C_LIGHT)) *
                                        (emlinewave < (linezwave + 5.*sigma_cont * linezwave / C_LIGHT)) *
                                        #(emlineflux * np.sqrt(emlineivar) > 1)
                                        (emlineivar > 0))[0]

                    linesigma = getattr(self.EMLineModel, '{}_sigma'.format(linename)).default # initial guess
                    #print(linename, linesigma, len(lineindx))

                    lineamp = pinfo.default # backup
                    if len(lineindx) > 10:
                        linesigma_ang = linezwave * sigma_cont / C_LIGHT # [observed-frame Angstrom]
                        linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang

                        lineflux = np.sum(emlineflux[lineindx])
                        lineamp = np.abs(lineflux / linenorm)

                        # estimate the velocity width from potentially strong, isolated lines; fragile!
                        if lineflux > 0 and linename in ['mgii_2800', 'oiii_5007', 'hbeta', 'halpha']:
                            linevar = np.sum(emlineflux[lineindx] * (emlinewave[lineindx] - linezwave)**2) / np.sum(emlineflux[lineindx]) / linezwave * C_LIGHT # [km/s]
                            if linevar > 0:
                                linesigma = np.sqrt(linevar)
                                init_linesigmas.append(linesigma)

                    if not pinfo.tied:# and False:
                        setattr(self.EMLineModel, pp, lineamp)
                    #print(pinfo.name, len(lineindx), lineflux, lineamp, linesigma)

            # update the initial velocity widths
            if len(init_linesigmas) >= 3:
                init_linesigma = np.median(init_linesigmas)
                if init_linesigma > 0 and init_linesigma < 300:
                    for pp in self.EMLineModel.param_names:
                        if 'sigma' in pp:
                            setattr(self.EMLineModel, pp, init_linesigma)

        # Fit [1]: tie all lines together
        self.EMLineModel = _tie_all_lines(self.EMLineModel)
        weights = np.sqrt(emlineivar)
        emlinebad = np.logical_not(emlinegood)
        if np.count_nonzero(emlinebad) > 0:
            weights[emlinebad] = 10*np.max(weights[emlinegood]) # 1e16 # ???
        
        t0 = time.time()
        bestfit_init = self.fitter(self.EMLineModel, emlinewave, emlineflux, weights=weights, maxiter=1000)
        #bestfit_init = self.fitter(self.EMLineModel, emlinewave[emlinegood], emlineflux[emlinegood],
        #                           weights=np.sqrt(emlineivar[emlinegood]), maxiter=1000)
        log.info('Initial line-fitting took {:.2f} sec'.format(time.time()-t0))
        #print(bestfit_init.parameters)

        # Fit [2]: tie Balmer, narrow forbidden, and QSO/broad lines together,
        # separately, and refit.
        self.EMLineModel = bestfit_init
        #self.EMLineModel = _tie_all_lines(self.EMLineModel) # this will reset our updates, above, so don't do it!
        self.EMLineModel = _tie_balmer_lines(self.EMLineModel)
        self.EMLineModel = _tie_forbidden_lines(self.EMLineModel)
        self.EMLineModel = _tie_broad_lines(self.EMLineModel)

        t0 = time.time()        
        bestfit = self.fitter(self.EMLineModel, emlinewave, emlineflux, weights=weights, maxiter=1000)
        #bestfit = self.fitter(self.EMLineModel, emlinewave[emlinegood], emlineflux[emlinegood],
        #                      weights=np.sqrt(emlineivar[emlinegood]), maxiter=1000)
        log.info('Final line-fitting took {:.2f} sec'.format(time.time()-t0))

        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_output(self.EMLineModel.linetable)

        # Populate the output table. First do a pass through the sigma
        # parameters. If sigma is zero, restore the default value and if the
        # amplitude is still at its default, it means the line wasn't fit
        # (right??), so set it to zero.
        for pp in bestfit.param_names:
            if 'sigma' in pp:
                pinfo = getattr(bestfit, pp)
                if pinfo.value == 0.0: # sigma=0, drop the line
                    setattr(bestfit, pp, pinfo.default)

        # Now fill the output table.
        for pp in bestfit.param_names:
            pinfo = getattr(bestfit, pp)
            val = pinfo.value
            #if '9532' in pp:
            #    pdb.set_trace()
            if 'amp' in pp and val == pinfo.default: # line not optimized; drop it!
                val = 0.0
                setattr(bestfit, pp, val)
            result[pinfo.name.upper()] = val

        emlinemodel = bestfit(emlinewave)
        chi2 = self.chi2(bestfit, emlinewave, emlineflux, emlineivar).astype('f4')

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

        ## Determine the uncertainties from the diagonal terms of the covariance
        ## matrix. If the covariance matrix is not known, estimate it from the
        ## Jacobian following:
        ##   https://github.com/scipy/scipy/blob/master/scipy/optimize/minpack.py#L805
        #pcov = self.fitter.fit_info['param_cov']
        #if pcov is None:
        #    pcov = np.zeros((nparam, nparam))
        #    #from scipy.linalg import svd
        #    #fjac = self.fitter.fit_info['fjac']
        #    #if fjac is not None:
        #    #    _, s, VT = svd(fjac, full_matrices=False)
        #    #    threshold = np.finfo(float).eps * max(fjac.shape) * s[0]
        #    #    s = s[s > threshold]
        #    #    VT = VT[:s.size]
        #    #    pcov = np.dot(VT.T / s**2, VT)
        #paramvar = np.diag(pcov)
        #ivar = np.zeros(nparam)

        # Need to be careful about uncertainties for tied parameters--
        # https://github.com/astropy/astropy/issues/7202
        
        #err_params = np.sqrt(np.diag(fitter.fit_info['param_cov']))
        #err = model.copy()
        #fitting._fitter_to_model_params(err, err_params)            

        # populate the output table
        #count = 0
        #for ii, pp in enumerate(bestfit.param_names):
        #    pinfo = getattr(bestfit, pp)
        #    result[pinfo.name.upper()] = pinfo.value.astype('f4')
        #    #iinfo = getattr(self.EMLineModel, pp)
        #
        #    # need to think about this more deeply
        #    #if pinfo.value == iinfo.value: # not fitted
        #    #    result.update({pinfo.name: np.float(0.0)})
        #    #else:
        #    #    result.update({pinfo.name: pinfo.value.astype('f4')})
        #    #result.update({pinfo.name: pinfo.value.astype('f4')})
        #    result[pinfo.name.upper()] = pinfo.value.astype('f4')
        #        
        #    #if pinfo.fixed:
        #    #    #result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
        #    #    result['{}_IVAR'.format(pinfo.name.upper())] = pinfo.value.astype('f4')
        #    #elif pinfo.tied:
        #    #    # hack! see https://github.com/astropy/astropy/issues/7202
        #    #    #result.update({'{}_ivar'.format(pinfo.name): np.float32(0.0)})
        #    #    result['{}_IVAR'.format(pinfo.name.upper())] = np.float32(0.0)
        #    #else:
        #    #    result['{}_IVAR'.format(pinfo.name.upper())] = ivar[count].astype('f4')
        #    #    #result.update({'{}_ivar'.format(pinfo.name): ivar[count].astype('f4')})
        #    #    count += 1

        ## hack for tied parameters---gotta be a better way to do this
        ##result['oiii_4959_amp_ivar'] = result['oiii_5007_amp_ivar'] * 2.8875**2
        #result['OIII_4959_AMP_IVAR'] = result['OIII_5007_AMP_IVAR'] * 2.8875**2
        #result['NII_6548_AMP_IVAR'] = result['NII_6548_AMP_IVAR'] * 2.936**2

        ## now loop back through and if ivar==0 then set the parameter value to zero
        #if False:
        #    if self.fitter.fit_info['param_cov'] is not None:
        #        for pp in bestfit.param_names:
        #            if result['{}_IVAR'.format(pp.upper())] == 0.0:
        #                result[pp] = np.float(0.0)

        # get continuum fluxes, EWs, and upper limits

        verbose = False

        balmer_sigmas, forbidden_sigmas, broad_sigmas = [], [], []
        balmer_redshifts, forbidden_redshifts, broad_redshifts = [], [], []
        for oneline in self.EMLineModel.linetable:

            linename = oneline['name'].upper()
            linez = redshift + result['{}_VSHIFT'.format(linename)][0] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)

            linesigma = result['{}_SIGMA'.format(linename)][0] # [km/s]
            log10sigma = linesigma / C_LIGHT / np.log(10)      # line-width [log-10 Angstrom]            

            # number of pixels, chi2, and boxcar integration
            lineindx = np.where((emlinewave > (linezwave - 4.*linesigma * linezwave / C_LIGHT)) *
                                (emlinewave < (linezwave + 4.*linesigma * linezwave / C_LIGHT)) *
                                (emlineivar > 0))[0]
            npix = len(lineindx)
            if npix > 3: # magic number: required at least XX unmasked pixels centered on the line
                dof = npix - 3 # ??? [redshift, sigma, and amplitude]
                chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-emlinemodel[lineindx])**2) / dof
                boxflux = np.sum(emlineflux[lineindx])
                boxflux_ivar = 1 / np.sum(1 / emlineivar[lineindx])

                # Get the uncertainty in the line-amplitude based on the scatter
                # in the pixel values.
                clipflux, _, _ = sigmaclip(specflux_nolines[lineindx], low=3, high=3)
                amp_sigma = np.std(clipflux)
                if amp_sigma > 0:
                    result['{}_AMP_IVAR'.format(linename)] = 1 / amp_sigma**2

                #if 'MGII' in linename or 'BETA' in linename or linename == 'OII_3729':
                #    import matplotlib.pyplot as plt
                #    plt.clf()
                #    plt.plot(emlinewave[lineindx], emlineflux[lineindx])
                #    plt.plot(emlinewave[lineindx], specflux_nolines[lineindx])
                #    plt.savefig('junk.png')
                #    pdb.set_trace()

                # only use 3-sigma lines
                if result['{}_AMP'.format(linename)] * np.sqrt(result['{}_AMP_IVAR'.format(linename)]) > 3:
                    if oneline['isbalmer']:
                        balmer_sigmas.append(linesigma)
                        balmer_redshifts.append(linez)
                    elif oneline['isbroad']:
                        broad_sigmas.append(linesigma)
                        broad_redshifts.append(linez)
                    else:
                        forbidden_sigmas.append(linesigma)
                        forbidden_redshifts.append(linez)
            else:
                npix, chi2, boxflux, boxflux_ivar = 0, 1e6, 0.0, 0.0
                
            result['{}_NPIX'.format(linename)] = npix
            result['{}_CHI2'.format(linename)] = chi2
            result['{}_BOXFLUX'.format(linename)] = boxflux
            #result['{}_BOXFLUX_IVAR'.format(linename)] = boxflux_ivar

            # get the emission-line flux
            linesigma_ang = linezwave * linesigma / C_LIGHT # [observed-frame Angstrom]
            linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang

            #print('{} sigma={:.3f} km/s npix={} chi2={:.3f} boxflux={:.3f} amp={:.3f}'.format(
            #    linename, linesigma, npix, chi2, boxflux, result['{}_AMP'.format(linename)][0]))
            
            result['{}_FLUX'.format(linename)] = result['{}_AMP'.format(linename)][0] * linenorm

            #result['{}_FLUX_IVAR'.format(linename)] = result['{}_AMP_IVAR'.format(linename)] / linenorm**2
            #weight = np.exp(-0.5 * np.log10(emlinewave/linezwave)**2 / log10sigma**2)
            #weight = (weight / np.max(weight)) > 1e-3
            #result['{}_FLUX_IVAR'.format(linename)] = 1 / np.sum(1 / emlineivar[weight])
            result['{}_FLUX_IVAR'.format(linename)] = boxflux_ivar

            # get the continuum, the inverse variance in the line-amplitude, and the EW
            indxlo = np.where((emlinewave > (linezwave - 10*linesigma * linezwave / C_LIGHT)) *
                              (emlinewave < (linezwave - 3.*linesigma * linezwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indxhi = np.where((emlinewave < (linezwave + 10*linesigma * linezwave / C_LIGHT)) *
                              (emlinewave > (linezwave + 3.*linesigma * linezwave / C_LIGHT)))[0]# *
                              #(emlinemodel == 0))[0]
            indx = np.hstack((indxlo, indxhi))

            cmed, civar = 0.0, 0.0
            if len(indx) > 10: # require at least XX pixels to get the continuum level
                #_, cmed, csig = sigma_clipped_stats(specflux_nolines[indx], sigma=3.0)
                clipflux, _, _ = sigmaclip(specflux_nolines[indx], low=3, high=3)
                cmed, csig = np.median(clipflux), np.std(clipflux)
                if csig > 0:
                    civar = (np.sqrt(len(indx)) / csig)**2
                    #result['{}_AMP_IVAR'.format(linename)] = 1 / csig**2
                
            result['{}_CONT'.format(linename)] = cmed
            result['{}_CONT_IVAR'.format(linename)] = civar

            ew, ewivar, fluxlimit, ewlimit = 0.0, 0.0, 0.0, 0.0
            if result['{}_CONT_IVAR'.format(linename)] != 0.0:
                #if result['{}_CONT'.format(linename)] == 0:
                #    pdb.set_trace()
                factor = (1 + redshift) / result['{}_CONT'.format(linename)] # --> rest frame
                ew = result['{}_FLUX'.format(linename)] * factor   # rest frame [A]
                ewivar = result['{}_FLUX_IVAR'.format(linename)] / factor**2

                # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar)
                ewlimit = fluxlimit * factor

            result['{}_EW'.format(linename)] = ew
            result['{}_EW_IVAR'.format(linename)] = ewivar
            result['{}_FLUX_LIMIT'.format(linename)] = fluxlimit
            result['{}_EW_LIMIT'.format(linename)] = ewlimit

            if verbose:
                for col in ('VSHIFT', 'SIGMA', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    print('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)][0]))
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT', 'EW_LIMIT'):
                    print('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)][0]))
                print()

            # simple QA
            if 'alpha' in linename and False:
                sigma_cont = 150.0
                import matplotlib.pyplot as plt
                _indx = np.where((emlinewave > (linezwave - 15*sigma_cont * linezwave / C_LIGHT)) *
                                (emlinewave < (linezwave + 15*sigma_cont * linezwave / C_LIGHT)))[0]
                plt.plot(emlinewave[_indx], emlineflux[_indx])
                plt.plot(emlinewave[_indx], specflux_nolines[_indx])
                plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
                plt.axhline(y=cmed, color='k')
                plt.axhline(y=cmed+csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.axhline(y=cmed-csig/np.sqrt(len(indx)), color='k', ls='--')
                plt.savefig('junk.png')

        # get the average emission-line redshifts and velocity widths
        if len(balmer_redshifts) > 0:
            result['BALMER_Z'] = np.mean(balmer_redshifts)
            result['BALMER_SIGMA'] = np.mean(balmer_sigmas)
            #result['BALMER_Z_ERR'] = np.std(balmer_redshifts)
            #result['BALMER_SIGMA_ERR'] = np.std(balmer_sigmas)
        else:
            result['BALMER_Z'] = redshift
            
        if len(forbidden_redshifts) > 0:
            result['FORBIDDEN_Z'] = np.mean(forbidden_redshifts)
            result['FORBIDDEN_SIGMA'] = np.mean(forbidden_sigmas)
            #result['FORBIDDEN_Z_ERR'] = np.std(forbidden_redshifts)
            #result['FORBIDDEN_SIGMA_ERR'] = np.std(forbidden_sigmas)
        else:
            result['FORBIDDEN_Z'] = redshift
            
        if len(broad_redshifts) > 0:
            result['BROAD_Z'] = np.mean(broad_redshifts)
            result['BROAD_SIGMA'] = np.mean(broad_sigmas)
            #result['BROAD_Z_ERR'] = np.std(broad_redshifts)
            #result['BROAD_SIGMA_ERR'] = np.std(broad_sigmas)
        else:
            result['BROAD_Z'] = redshift
            
        if verbose:
            for line in ('BALMER', 'FORBIDDEN', 'BROAD'):
                for col in ('Z', 'SIGMA'):
                #for col in ('Z', 'Z_ERR', 'SIGMA', 'SIGMA_ERR'):
                    print('{}_{}: {:.12f}'.format(line, col, result['{}_{}'.format(line, col)][0]))

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
        col3 = [colors.to_hex(col) for col in ['blue', 'green', 'red']]

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
        
        #smooth_continuum = self.smooth_residuals(data['flux'], continuum, data['linemask'])
        smooth_continuum = self.smooth_residuals(
            continuum, data['wave'], data['flux'],
            data['ivar'], data['linemask'])
        
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
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            'dv_balmer': '$\Delta v_{{\\rm Balmer}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['BALMER_Z']-redshift)),
            'dv_forbid': '$\Delta v_{{\\rm forbid}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['FORBIDDEN_Z']-redshift)),
            'dv_broad': '$\Delta v_{{\\rm MgII}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['BROAD_Z']-redshift)),
            #'zbalmer': '$z_{{\\rm Balmer}}$={:.6f}'.format(fastspec['BALMER_Z']),
            #'zforbidden': '$z_{{\\rm forbidden}}$={:.6f}'.format(fastspec['FORBIDDEN_Z']),
            #'zbroad': '$z_{{\\rm MgII}}$={:.6f}'.format(fastspec['BROAD_Z']),
            'sigma_balmer': '$\sigma_{{\\rm Balmer}}$={:.1f} km/s'.format(fastspec['BALMER_SIGMA']),
            'sigma_forbid': '$\sigma_{{\\rm forbid}}$={:.1f} km/s'.format(fastspec['FORBIDDEN_SIGMA']),
            'sigma_broad': '$\sigma_{{\\rm MgII}}$={:.1f} km/s'.format(fastspec['BROAD_SIGMA']),
            #'targetid': '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastspec['CONTINUUM_CHI2']),
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
            leg.update({'AV': '$A(V)$={:.2f}+/-{:.2f} mag'.format(
                fastspec['CONTINUUM_AV'], 1/np.sqrt(fastspec['CONTINUUM_AV_IVAR']))})

        ymin, ymax = 1e6, -1e6
        wavelims = (3600, 9800)

        legxpos, legypos, legfntsz = 0.98, 0.94, 20
        bbox = dict(boxstyle='round', facecolor='gray', alpha=0.25)
        
        for ii in [0, 1, 2]: # iterate over cameras
            sigma, _ = ivar2var(data['ivar'][ii], sigma=True)

            bigax1.fill_between(data['wave'][ii], data['flux'][ii]-sigma,
                             data['flux'][ii]+sigma, color=col1[ii])
            #bigax1.plot(data['wave'][ii], continuum[ii], color=col2[ii], alpha=1.0)#, color='k')
            #bigax1.plot(data['wave'][ii], continuum_nodust[ii], alpha=0.5, color='k')
            bigax1.plot(data['wave'][ii], smooth_continuum[ii], color='gray')#col3[ii])#, alpha=0.3, lw=2)#, color='k')
            bigax1.plot(data['wave'][ii], continuum[ii]+smooth_continuum[ii], color=col2[ii])
            
            # get the robust range
            filtflux = median_filter(data['flux'][ii], 51, mode='nearest')
            #filtflux = median_filter(data['flux'][ii] - _emlinemodel[ii], 51, mode='nearest')
            #perc = np.percentile(filtflux[data['ivar'][ii] > 0], [5, 95])
            sigflux = np.std(data['flux'][ii][data['ivar'][ii] > 0])
            #sigflux = np.std(filtflux[data['ivar'][ii] > 0])
            #if -2 * perc[0] < ymin:
            #    ymin = -2 * perc[0]
            if -2 * sigflux < ymin:
                ymin = -2 * sigflux
            #if perc[1] > ymax:
            #    ymax = perc[1]
            #if np.min(filtflux) < ymin:
            #    ymin = np.min(filtflux) * 0.5
            if sigflux * 5 > ymax:
                ymax = sigflux * 5
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux) * 1.4
            print(ymin, ymax)
            #if ii == 2:
            #    pdb.set_trace()

        txt = '\n'.join((
            r'{}'.format(leg['zredrock']),
            r'{} {}'.format(leg['chi2'], leg['age']),
            r'{}'.format(leg['AV']),
            #r'{}'.format(leg['vdisp']),
            ))
        bigax1.text(legxpos, legypos, txt, ha='right', va='top',
                    transform=bigax1.transAxes, fontsize=legfntsz,
                    bbox=bbox)
        bigax1.text(0.03, 0.9, 'Observed Spectrum + Continuum Model',
                    ha='left', va='center', transform=bigax1.transAxes, fontsize=22)
        bigax1.set_title(title, fontsize=22)
        
        bigax1.set_xlim(wavelims)
        bigax1.set_ylim(ymin, ymax)
        bigax1.set_xticklabels([])
        #bigax1.set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
        #bigax1.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax1.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 

        # full emission-line spectrum + best-fitting lines
        bigax2 = fig.add_subplot(gs[1, :])

        ymin, ymax = 1e6, -1e6
        for ii in [0, 1, 2]: # iterate over cameras
            emlinewave = data['wave'][ii]
            emlineflux = data['flux'][ii] - continuum[ii] - smooth_continuum[ii]
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
            filtflux = median_filter(emlineflux, 51, mode='nearest')
            #sigflux = np.std(filtflux)
            sigflux = np.std(emlineflux)

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
            print(ymin, ymax)

        txt = '\n'.join((
            r'{} {}'.format(leg['dv_balmer'], leg['sigma_balmer']),
            r'{} {}'.format(leg['dv_forbid'], leg['sigma_forbid']),
            r'{} {}'.format(leg['dv_broad'], leg['sigma_broad']),
            ))
        bigax2.text(legxpos, legypos, txt, ha='right', va='top',
                    transform=bigax2.transAxes, fontsize=legfntsz,
                    bbox=bbox)
        bigax2.text(0.03, 0.9, 'Residual Spectrum + Emission-Line Model',
                    ha='left', va='center', transform=bigax2.transAxes,
                    fontsize=22)
                
        bigax2.set_xlim(wavelims)
        bigax2.set_ylim(ymin, ymax)
        
        #bigax2.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax2.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 
        
        # zoom in on individual emission lines - use linetable!
        plotsig_default = 300.0 # [km/s]

        meanwaves, deltawaves, sigmas, linenames = [], [], [], []
        for plotgroup in set(self.linetable['plotgroup']):
            I = np.where(plotgroup == self.linetable['plotgroup'])[0]
            linenames.append(self.linetable['nicename'][I[0]])
            meanwaves.append(np.mean(self.linetable['restwave'][I]))
            deltawaves.append((np.max(self.linetable['restwave'][I]) - np.min(self.linetable['restwave'][I])) / 2)

            plotsig = None
            if np.any(self.linetable['isbroad'][I]):
                plotsig = fastspec['BROAD_SIGMA']
            elif np.any(self.linetable['isbalmer'][I]):
                plotsig = np.max((fastspec['BALMER_SIGMA'], fastspec['FORBIDDEN_SIGMA']))
            else:
                plotsig = fastspec['FORBIDDEN_SIGMA']
            if plotsig:
                sigmas.append(plotsig)
            else:
                sigmas.append(plotsig_default)

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
                emlineflux = data['flux'][ii] - continuum[ii] - smooth_continuum[ii]
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
                    xx.axhline(y=0, color='gray', ls='-')

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
                    
                xx.text(0.08, 0.89, linename, ha='left', va='center',
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
        
        fig.text(lf-0.07, (tp-bt)/2+bt,
                 r'Flux Density ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
                 ha='center', va='center', rotation='vertical')
        fig.text((rt-lf)/2+lf, bt-0.05, r'Observed-frame Wavelength ($\AA$)',
                 ha='center', va='center')
            
        plt.subplots_adjust(wspace=0.27, top=tp, bottom=bt, left=lf, right=rt, hspace=0.22)

        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()
        #pdb.set_trace()
