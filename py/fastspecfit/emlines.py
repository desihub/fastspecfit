"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

python -m cProfile -o fastspec.prof /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/20210324/redrock-4-80613-thru20210324.fits -o fastspec.fits --targetids 39633345008634465 --specprod everest

"""
import pdb # for debugging

import os, time
import numpy as np
import multiprocessing

import astropy.units as u
from astropy.modeling import Fittable1DModel
from astropy.modeling.fitting import LevMarLSQFitter

from redrock.rebin import trapz_rebin
from desispec.interpolation import resample_flux

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

# doublet ratios are always tied
def _tie_nii_amp(model):
    """
    [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.022e-21
    [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.949e-21
    """
    return model.nii_6584_amp / 2.9421 # 2.936

def _tie_oiii_amp(model):
    """
    [O3] (4-->2): airwave: 4958.9097 vacwave: 4960.2937 emissivity: 1.172e-21
    [O3] (4-->3): airwave: 5006.8417 vacwave: 5008.2383 emissivity: 3.497e-21
    """
    return model.oiii_5007_amp / 2.9839 # 2.8875

def _tie_oii_amp(model):
    """
    [O2] (2-->1): airwave: 3728.8145 vacwave: 3729.8750 emissivity: 1.948e-21
    [O2] (3-->1): airwave: 3726.0322 vacwave: 3727.0919 emissivity: 1.444e-21
    """
    return model.oii_3729_amp / 1.3490
    
def _tie_hbeta_sigma(model):
    return model.hbeta_sigma
def _tie_hbeta_vshift(model):
    return model.hbeta_vshift

def _tie_oiii_5007_sigma(model):
    return model.oiii_5007_sigma
def _tie_oiii_5007_vshift(model):
    return model.oiii_5007_vshift

def _tie_mgii_2803_sigma(model):
    return model.mgii_2803_sigma
def _tie_mgii_2803_vshift(model):
    return model.mgii_2803_vshift

def _tie_civ_1550_sigma(model):
    return model.civ_1550_sigma
def _tie_civ_1550_vshift(model):
    return model.civ_1550_vshift

def _tie_siliv_1403_sigma(model):
    return model.siliv_1403_sigma
def _tie_siliv_1403_vshift(model):
    return model.siliv_1403_vshift

def _tie_nv_1243_sigma(model):
    return model.nv_1243_sigma
def _tie_nv_1243_vshift(model):
    return model.nv_1243_vshift

def _tie_lines(model):
    """Tie emission lines together using sensible constrains based on galaxy & QSO
    physics.

    """
    for iline in np.arange(len(model.linetable)):
        linename = model.linetable['name'][iline]

        if model.inrange[iline] == False: # not in the wavelength range
            setattr(model, '{}_amp'.format(linename), 0.0)
            getattr(model, '{}_amp'.format(linename)).fixed = True
            getattr(model, '{}_sigma'.format(linename)).fixed = True
            getattr(model, '{}_vshift'.format(linename)).fixed = True

        # Balmer & He lines
        if model.inrange[iline] and model.linetable['isbalmer'][iline] and linename != 'hbeta':
            getattr(model, '{}_sigma'.format(linename)).tied = _tie_hbeta_sigma
            getattr(model, '{}_vshift'.format(linename)).tied = _tie_hbeta_vshift
            
        # broad lines
        if model.inrange[iline] and model.linetable['isbroad'][iline] and model.linetable['isbalmer'][iline] == False:
            if linename == 'mgii_2796':
                getattr(model, '{}_sigma'.format(linename)).tied = _tie_mgii_2803_sigma
                getattr(model, '{}_vshift'.format(linename)).tied = _tie_mgii_2803_vshift
            if linename == 'civ_1548':
                getattr(model, '{}_sigma'.format(linename)).tied = _tie_civ_1550_sigma
                getattr(model, '{}_vshift'.format(linename)).tied = _tie_civ_1550_vshift
            if linename == 'siliv_1394':
                getattr(model, '{}_sigma'.format(linename)).tied = _tie_siliv_1403_sigma
                getattr(model, '{}_vshift'.format(linename)).tied = _tie_siliv_1403_vshift
            if linename == 'nv_1239':
                getattr(model, '{}_sigma'.format(linename)).tied = _tie_nv_1243_sigma
                getattr(model, '{}_vshift'.format(linename)).tied = _tie_nv_1243_vshift
            
        # forbidden lines
        if (model.inrange[iline] and model.linetable['isbroad'][iline] == False and
            model.linetable['isbalmer'][iline] == False and linename != 'oiii_5007'):
            getattr(model, '{}_sigma'.format(linename)).tied = _tie_oiii_5007_sigma
            getattr(model, '{}_vshift'.format(linename)).tied = _tie_oiii_5007_vshift

    #for pp in model.param_names:
    #    print(getattr(model, pp))
    
    return model

#def _tie_all_lines(model):
#    for pp in model.param_names:
#        if 'sigma' in pp and pp != 'hbeta_sigma':
#            getattr(model, pp).tied = _tie_hbeta_sigma
#        if 'vshift' in pp and pp != 'hbeta_vshift':
#            getattr(model, pp).tied = _tie_hbeta_vshift
#            
#    # reset the amplitudes of the tied doublets; fragile...
#    model.oiii_4959_amp = model.oiii_4959_amp.tied(model)
#    model.nii_6548_amp = model.nii_6548_amp.tied(model)
#            
#    return model

class FastLevMarLSQFitter(LevMarLSQFitter):
    """
    Adopted in part from
      https://github.com/astropy/astropy/blob/v4.2.x/astropy/modeling/fitting.py#L1580-L1671
    
    Copyright (c) 2011-2021, Astropy Developers

    See https://github.com/astropy/astropy/issues/12089 for details.  

    """
    def __init__(self, model):

        import operator        
        from functools import reduce

        self.has_tied = any(model.tied.values())
        self.has_fixed = any(model.fixed.values())
        self.has_bound = any(b != (None, None) for b in model.bounds.values())

        _, fit_param_indices = self.model_to_fit_params(model)

        self._fit_param_names = [model.param_names[iparam] for iparam in np.unique(fit_param_indices)]

        self._fit_istart = []
        self._fit_iend = []
        self._fit_slice = []
        self._fit_ibounds = []
        self._tied_slice = []
        self._tied_func = []

        offset = 0
        for name in self._fit_param_names:
            slice_ = model._param_metrics[name]['slice']
            shape = model._param_metrics[name]['shape']
            # This is determining which range of fps (the fitted parameters) maps
            # to parameters of the model
            size = reduce(operator.mul, shape, 1)
            self._fit_slice.append(slice_)
            self._fit_istart.append(offset)
            self._fit_iend.append(offset + size)

            imin, imax = model.bounds[name]
            self._fit_ibounds.append(((imin, imax) != (None, None), imin, imax))
            offset += size            

        if self.has_tied:
            self._tied_param_names = [_param_name for _param_name in model.param_names if model.tied[_param_name]]
            for name in self._tied_param_names:
                self._tied_slice.append(model._param_metrics[name]['slice'])
                self._tied_func.append(model.tied[name])
        
        super().__init__()

    def objective_function(self, fps, *args):

        model = args[0]
        weights = args[1]
        self.fitter_to_model_params(model, fps)

        meas = args[-1]
        if weights is None:
            return np.ravel(model(*args[2: -1]) - meas)
        else:
            return np.ravel(weights * (model(*args[2: -1]) - meas))

    def fitter_to_model_params(self, model, fps):
        """Constructs the full list of model parameters from the fitted and constrained
        parameters.

        """
        parameters = model.parameters

        if not (self.has_tied or self.has_fixed or self.has_bound):
            # We can just assign directly
            model.parameters = fps
            return

        for name, istart, iend, islice, (ibounds, imin, imax) in zip(
                        self._fit_param_names, self._fit_istart, self._fit_iend,
                        self._fit_slice, self._fit_ibounds):
            values = fps[istart:iend]

            # Check bounds constraints
            if ibounds:
                if imin is not None:
                    values = np.fmax(values, imin)
                if imax is not None:
                    values = np.fmin(values, imax)
                    
            #parameters[islice] = values
            setattr(model, name, values)

        # Update model parameters before calling ``tied`` constraints.

        # This has to be done in a separate loop due to how tied parameters are
        # currently evaluated (the fitted parameters need to actually be *set* on
        # the model first, for use in evaluating the "tied" expression--it might be
        # better to change this at some point
        if self.has_tied:
            for name, islice, func in zip(
                    self._tied_param_names, self._tied_slice,
                    self._tied_func):
                #value = model.tied[name](model)

                # To handle multiple tied constraints, model parameters
                # need to be updated after each iteration.
                value = func(model)
                #parameters[islice] = value
                setattr(model, name, value)
                
    def model_to_fit_params(self, model):
        """Convert a model instance's parameter array to an array that can be used with
        a fitter that doesn't natively support fixed or tied parameters.  In
        particular, it removes fixed/tied parameters from the parameter array.
        These may be a subset of the model parameters, if some of them are held
        constant or tied.

        """
        fitparam_indices = list(range(len(model.param_names)))
        if self.has_fixed or self.has_tied:
            params = list(model.parameters)
            param_metrics = model._param_metrics
            for idx, name in list(enumerate(model.param_names))[::-1]:
                if model.fixed[name] or model.tied[name]:
                    slice_ = param_metrics[name]['slice']
                    del params[slice_]
                    del fitparam_indices[idx]
            return (np.array(params), fitparam_indices)
        return (model.parameters, fitparam_indices)

class EMLineModel(Fittable1DModel):
    """Class to model the emission-line spectra.

    """
    from astropy.modeling import Parameter

    # NB! The order of the parameters here matters!
    vmaxshift_narrow = 300.0
    vmaxshift_broad = 3000.0
    initvshift = 0.0

    minsigma = 30.0
    maxsigma_narrow = 500.0
    maxsigma_broad = 5000.0

    initsigma_narrow = 75.0
    initsigma_broad = 300.0

    minamp = -1e3
    maxamp = +1e5

    # Fragile because the lines are hard-coded--
    nv_1239_amp = Parameter(name='nv_1239_amp', default=3.0, bounds=[minamp, maxamp]) # bounds=[0, maxamp])
    nv_1243_amp = Parameter(name='nv_1243_amp', default=3.0, bounds=[minamp, maxamp])
    oi_1304_amp = Parameter(name='oi_1304_amp', default=3.0, bounds=[minamp, maxamp])
    silii_1307_amp = Parameter(name='silii_1307_amp', default=3.0, bounds=[minamp, maxamp]) # bounds=[0, maxamp])
    siliv_1394_amp = Parameter(name='siliv_1394_amp', default=3.0, bounds=[minamp, maxamp])
    siliv_1403_amp = Parameter(name='siliv_1403_amp', default=3.0, bounds=[minamp, maxamp])
    civ_1548_amp = Parameter(name='civ_1548_amp', default=3.0, bounds=[minamp, maxamp]) # bounds=[0, maxamp])
    civ_1550_amp = Parameter(name='civ_1550_amp', default=3.0, bounds=[minamp, maxamp])
    siliii_1892_amp = Parameter(name='siliii_1892_amp', default=3.0, bounds=[minamp, maxamp])
    ciii_1908_amp = Parameter(name='ciii_1908_amp', default=3.0, bounds=[minamp, maxamp])
    mgii_2796_amp = Parameter(name='mgii_2796_amp', default=3.0, bounds=[minamp, maxamp]) # bounds=[0, maxamp])
    mgii_2803_amp = Parameter(name='mgii_2803_amp', default=3.0, bounds=[minamp, maxamp])
    nev_3346_amp = Parameter(name='nev_3346_amp', default=0.1, bounds=[minamp, maxamp])
    nev_3426_amp = Parameter(name='nev_3426_amp', default=0.1, bounds=[minamp, maxamp])
    oii_3726_amp = Parameter(name='oii_3726_amp', default=1.0, bounds=[minamp, maxamp])
    oii_3729_amp = Parameter(name='oii_3729_amp', default=1.0, bounds=[minamp, maxamp])
    neiii_3869_amp = Parameter(name='neiii_3869_amp', default=0.3, bounds=[minamp, maxamp])
    hei_3889_amp = Parameter(name='hei_3889_amp', default=0.3, bounds=[minamp, maxamp])
    h6_amp = Parameter(name='h6_amp', default=0.3, bounds=[minamp, maxamp])
    hepsilon_amp = Parameter(name='hepsilon_amp', default=0.5, bounds=[minamp, maxamp])
    hdelta_amp = Parameter(name='hdelta_amp', default=0.5, bounds=[minamp, maxamp])
    hgamma_amp = Parameter(name='hgamma_amp', default=0.5, bounds=[minamp, maxamp])
    oiii_4363_amp = Parameter(name='oiii_4363_amp', default=0.3, bounds=[minamp, maxamp])
    hei_4471_amp = Parameter(name='hei_4471_amp', default=0.3, bounds=[minamp, maxamp])
    heii_4686_amp = Parameter(name='heii_4686_amp', default=0.3, bounds=[minamp, maxamp])
    hbeta_amp = Parameter(name='hbeta_amp', default=1.0, bounds=[minamp, maxamp])
    oiii_4959_amp = Parameter(name='oiii_4959_amp', default=1.0, bounds=[minamp, maxamp])
    oiii_5007_amp = Parameter(name='oiii_5007_amp', default=3.0, bounds=[minamp, maxamp])
    nii_5755_amp = Parameter(name='nii_5755_amp', default=0.3, bounds=[minamp, maxamp])
    hei_5876_amp = Parameter(name='hei_5876_amp', default=0.3, bounds=[minamp, maxamp])
    oi_6300_amp = Parameter(name='oi_6300_amp', default=0.3, bounds=[minamp, maxamp])
    nii_6548_amp = Parameter(name='nii_6548_amp', default=1.0, bounds=[minamp, maxamp])
    halpha_amp = Parameter(name='halpha_amp', default=3.0, bounds=[minamp, maxamp])
    nii_6584_amp = Parameter(name='nii_6584_amp', default=3.0, bounds=[minamp, maxamp])
    sii_6716_amp = Parameter(name='sii_6716_amp', default=1.0, bounds=[minamp, maxamp])
    sii_6731_amp = Parameter(name='sii_6731_amp', default=1.0, bounds=[minamp, maxamp])
    siii_9069_amp = Parameter(name='siii_9069_amp', default=0.3, bounds=[minamp, maxamp])
    siii_9532_amp = Parameter(name='siii_9532_amp', default=0.3, bounds=[minamp, maxamp])

    nv_1239_vshift = Parameter(name='nv_1239_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    nv_1243_vshift = Parameter(name='nv_1243_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    oi_1304_vshift = Parameter(name='oi_1304_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    silii_1307_vshift = Parameter(name='silii_1307_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    siliv_1394_vshift = Parameter(name='siliv_1395_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    siliv_1403_vshift = Parameter(name='siliv_1403_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    civ_1548_vshift = Parameter(name='civ_1548_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    civ_1550_vshift = Parameter(name='civ_1550_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    siliii_1892_vshift = Parameter(name='siliii_1892_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    ciii_1908_vshift = Parameter(name='ciii_1908_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    mgii_2796_vshift = Parameter(name='mgii_2796_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    mgii_2803_vshift = Parameter(name='mgii_2803_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    nev_3346_vshift = Parameter(name='nev_3346_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nev_3426_vshift = Parameter(name='nev_3426_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oii_3726_vshift = Parameter(name='oii_3726_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oii_3729_vshift = Parameter(name='oii_3729_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    neiii_3869_vshift = Parameter(name='neiii_3869_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_3889_vshift = Parameter(name='hei_3889_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    h6_vshift = Parameter(name='h6_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hepsilon_vshift = Parameter(name='hepsilon_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hdelta_vshift = Parameter(name='hdelta_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hgamma_vshift = Parameter(name='hgamma_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oiii_4363_vshift = Parameter(name='oiii_4363_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_4471_vshift = Parameter(name='hei_4471_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    heii_4686_vshift = Parameter(name='heii_4686_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hbeta_vshift = Parameter(name='hbeta_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oiii_4959_vshift = Parameter(name='oiii_4959_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oiii_5007_vshift = Parameter(name='oiii_5007_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nii_5755_vshift = Parameter(name='nii_5755_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_5876_vshift = Parameter(name='hei_5876_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oi_6300_vshift = Parameter(name='oi_6300_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nii_6548_vshift = Parameter(name='nii_6548_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    halpha_vshift = Parameter(name='halpha_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nii_6584_vshift = Parameter(name='nii_6584_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    sii_6716_vshift = Parameter(name='sii_6716_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    sii_6731_vshift = Parameter(name='sii_6731_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    siii_9069_vshift = Parameter(name='siii_9069_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    siii_9532_vshift = Parameter(name='siii_9532_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])

    nv_1239_sigma = Parameter(name='nv_1239_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    nv_1243_sigma = Parameter(name='nv_1243_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    oi_1304_sigma = Parameter(name='oi_1304_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    silii_1307_sigma = Parameter(name='silii_1307_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    siliv_1394_sigma = Parameter(name='siliv_1395_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    siliv_1403_sigma = Parameter(name='siliv_1403_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    civ_1548_sigma = Parameter(name='civ_1548_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    civ_1550_sigma = Parameter(name='civ_1550_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    siliii_1892_sigma = Parameter(name='siliii_1892_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    ciii_1908_sigma = Parameter(name='ciii_1908_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    mgii_2796_sigma = Parameter(name='mgii_2796_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    mgii_2803_sigma = Parameter(name='mgii_2803_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    nev_3346_sigma = Parameter(name='nev_3346_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    nev_3426_sigma = Parameter(name='nev_3426_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    oii_3726_sigma = Parameter(name='oii_3726_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    oii_3729_sigma = Parameter(name='oii_3729_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    neiii_3869_sigma = Parameter(name='neiii_3869_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_3889_sigma = Parameter(name='hei_3889_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    h6_sigma = Parameter(name='h6_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    hepsilon_sigma = Parameter(name='hepsilon_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    hdelta_sigma = Parameter(name='hdelta_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    hgamma_sigma = Parameter(name='hgamma_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    oiii_4363_sigma = Parameter(name='oiii_4363_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_4471_sigma = Parameter(name='hei_4471_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    heii_4686_sigma = Parameter(name='heii_4686_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    hbeta_sigma = Parameter(name='hbeta_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    oiii_4959_sigma = Parameter(name='oiii_4959_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    oiii_5007_sigma = Parameter(name='oiii_5007_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    nii_5755_sigma = Parameter(name='nii_5755_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_5876_sigma = Parameter(name='hei_5876_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    oi_6300_sigma = Parameter(name='oi_6300_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    nii_6548_sigma = Parameter(name='nii_6548_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    halpha_sigma = Parameter(name='halpha_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_broad])
    nii_6584_sigma = Parameter(name='nii_6584_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    sii_6716_sigma = Parameter(name='sii_6716_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    sii_6731_sigma = Parameter(name='sii_6731_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    siii_9069_sigma = Parameter(name='siii_9069_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    siii_9532_sigma = Parameter(name='siii_9532_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])

    nii_6548_amp.tied = _tie_nii_amp
    oiii_4959_amp.tied = _tie_oiii_amp
    oii_3726_amp.tied = _tie_oii_amp
    
    def __init__(self,
                 nv_1239_amp=nv_1239_amp.default,
                 nv_1243_amp=nv_1243_amp.default,
                 oi_1304_amp=oi_1304_amp.default,
                 silii_1307_amp=silii_1307_amp.default,
                 siliv_1394_amp=siliv_1394_amp.default,
                 siliv_1403_amp=siliv_1403_amp.default,
                 civ_1548_amp=civ_1548_amp.default,
                 civ_1550_amp=civ_1550_amp.default,
                 siliii_1892_amp=siliii_1892_amp.default,
                 ciii_1908_amp=ciii_1908_amp.default,
                 mgii_2796_amp=mgii_2796_amp.default,
                 mgii_2803_amp=mgii_2803_amp.default,
                 nev_3346_amp=nev_3346_amp.default,
                 nev_3426_amp=nev_3426_amp.default,
                 oii_3726_amp=oii_3726_amp.default,
                 oii_3729_amp=oii_3729_amp.default,
                 neiii_3869_amp=neiii_3869_amp.default,
                 hei_3889_amp=hei_3889_amp.default,
                 h6_amp=h6_amp.default,
                 hepsilon_amp=hepsilon_amp.default,
                 hdelta_amp=hdelta_amp.default,
                 hgamma_amp=hgamma_amp.default,
                 oiii_4363_amp=oiii_4363_amp.default,
                 hei_4471_amp=hei_4471_amp.default,
                 heii_4686_amp=heii_4686_amp.default,
                 hbeta_amp=hbeta_amp.default,
                 oiii_4959_amp=oiii_4959_amp.default,
                 oiii_5007_amp=oiii_5007_amp.default,
                 nii_5755_amp=nii_5755_amp.default,
                 hei_5876_amp=hei_5876_amp.default,
                 oi_6300_amp=oi_6300_amp.default,
                 nii_6548_amp=nii_6548_amp.default,
                 halpha_amp=halpha_amp.default,
                 nii_6584_amp=nii_6584_amp.default,
                 sii_6716_amp=sii_6716_amp.default,
                 sii_6731_amp=sii_6731_amp.default,
                 siii_9069_amp=siii_9069_amp.default,
                 siii_9532_amp=siii_9532_amp.default,
                     
                 nv_1239_vshift=nv_1239_vshift.default,
                 nv_1243_vshift=nv_1243_vshift.default,
                 oi_1304_vshift=oi_1304_vshift.default,
                 silii_1307_vshift=silii_1307_vshift.default,
                 siliv_1394_vshift=siliv_1394_vshift.default,
                 siliv_1403_vshift=siliv_1403_vshift.default,
                 civ_1548_vshift=civ_1548_vshift.default,
                 civ_1550_vshift=civ_1550_vshift.default,
                 siliii_1892_vshift=siliii_1892_vshift.default,
                 ciii_1908_vshift=ciii_1908_vshift.default,
                 mgii_2796_vshift=mgii_2796_vshift.default,
                 mgii_2803_vshift=mgii_2803_vshift.default,
                 nev_3346_vshift=nev_3346_vshift.default,
                 nev_3426_vshift=nev_3426_vshift.default,
                 oii_3726_vshift=oii_3726_vshift.default,
                 oii_3729_vshift=oii_3729_vshift.default,
                 neiii_3869_vshift=neiii_3869_vshift.default,
                 hei_3889_vshift=hei_3889_vshift.default,
                 h6_vshift=h6_vshift.default,
                 hepsilon_vshift=hepsilon_vshift.default,
                 hdelta_vshift=hdelta_vshift.default,
                 hgamma_vshift=hgamma_vshift.default,
                 oiii_4363_vshift=oiii_4363_vshift.default,
                 hei_4471_vshift=hei_4471_vshift.default,
                 heii_4686_vshift=heii_4686_vshift.default,
                 hbeta_vshift=hbeta_vshift.default,
                 oiii_4959_vshift=oiii_4959_vshift.default,
                 oiii_5007_vshift=oiii_5007_vshift.default,
                 nii_5755_vshift=nii_5755_vshift.default,
                 hei_5876_vshift=hei_5876_vshift.default,
                 oi_6300_vshift=oi_6300_vshift.default,
                 nii_6548_vshift=nii_6548_vshift.default,
                 halpha_vshift=halpha_vshift.default,
                 nii_6584_vshift=nii_6584_vshift.default,
                 sii_6716_vshift=sii_6716_vshift.default,
                 sii_6731_vshift=sii_6731_vshift.default,
                 siii_9069_vshift=siii_9069_vshift.default,
                 siii_9532_vshift=siii_9532_vshift.default,
                     
                 nv_1239_sigma=nv_1239_sigma.default,
                 nv_1243_sigma=nv_1243_sigma.default,
                 oi_1304_sigma=oi_1304_sigma.default,
                 silii_1307_sigma=silii_1307_sigma.default,
                 siliv_1394_sigma=siliv_1394_sigma.default,
                 siliv_1403_sigma=siliv_1403_sigma.default,
                 civ_1548_sigma=civ_1548_sigma.default,
                 civ_1550_sigma=civ_1550_sigma.default,
                 siliii_1892_sigma=siliii_1892_sigma.default,
                 ciii_1908_sigma=ciii_1908_sigma.default,
                 mgii_2796_sigma=mgii_2796_sigma.default,
                 mgii_2803_sigma=mgii_2803_sigma.default,
                 nev_3346_sigma=nev_3346_sigma.default,
                 nev_3426_sigma=nev_3426_sigma.default,
                 oii_3726_sigma=oii_3726_sigma.default,
                 oii_3729_sigma=oii_3729_sigma.default,
                 neiii_3869_sigma=neiii_3869_sigma.default,
                 hei_3889_sigma=hei_3889_sigma.default,
                 h6_sigma=h6_sigma.default,
                 hepsilon_sigma=hepsilon_sigma.default,
                 hdelta_sigma=hdelta_sigma.default,
                 hgamma_sigma=hgamma_sigma.default,
                 oiii_4363_sigma=oiii_4363_sigma.default,
                 hei_4471_sigma=hei_4471_sigma.default,
                 heii_4686_sigma=heii_4686_sigma.default,
                 hbeta_sigma=hbeta_sigma.default,
                 oiii_4959_sigma=oiii_4959_sigma.default,
                 oiii_5007_sigma=oiii_5007_sigma.default,
                 nii_5755_sigma=nii_5755_sigma.default,
                 hei_5876_sigma=hei_5876_sigma.default,
                 oi_6300_sigma=oi_6300_sigma.default,
                 nii_6548_sigma=nii_6548_sigma.default,
                 halpha_sigma=halpha_sigma.default,
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
        super().__init__()        

        self.linetable = read_emlines()
        self.nline = len(self.linetable)
        self.linenames = np.array([linename.replace('_amp', '') for linename in self.param_names[:self.nline]])
        self.restlinewaves = self.linetable['restwave'].data

        self.redshift = redshift
        self.emlineR = emlineR
        if npixpercamera is not None:
            self.npixpercamera = np.hstack([0, npixpercamera])

        # internal wavelength vector for building the emission-line model
        if log10wave is None:
            pixkms = 10.0
            dlogwave = pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
            log10wave = np.arange(np.log10(3000), np.log10(1e4), dlogwave)
        self.log10wave = log10wave
            
        wavelims = (10**np.min(log10wave), 10**np.max(log10wave))
        zline = self.linetable['restwave'] * (1 + self.redshift)
        self.inrange = (zline > wavelims[0]) * (zline < wavelims[1])

        # tie lines together
        _tie_lines(self)

    def _emline_spectrum(self, *lineargs):
        """Simple wrapper to build an emission-line spectrum.

        # build the emission-line model [erg/s/cm2/A, observed frame]

        """
        #linenames = self.linenames[self.inrange]
        lineamps = np.hstack(lineargs[0:self.nline])[self.inrange]
        linevshifts = np.hstack(lineargs[self.nline:2*self.nline])[self.inrange]
        linesigmas = np.hstack(lineargs[2*self.nline:])[self.inrange]

        # line-width [log-10 Angstrom] and redshifted wavelength [log-10 Angstrom]
        log10sigmas = linesigmas / C_LIGHT / np.log(10) 
        linezwaves = np.log10(self.restlinewaves[self.inrange] * (1.0 + self.redshift + linevshifts / C_LIGHT)) 

        log10model = np.zeros_like(self.log10wave)
        for lineamp, linezwave, log10sigma in zip(lineamps, linezwaves, log10sigmas):
            ww = np.abs(self.log10wave - linezwave) < (10 * log10sigma)
            if np.sum(ww) > 0:
                #log.info(linename, 10**linezwave, 10**_emlinewave[ww].min(), 10**_emlinewave[ww].max())
                log10model[ww] += lineamp * np.exp(-0.5 * (self.log10wave[ww]-linezwave)**2 / log10sigma**2)

        return log10model
    
    def evaluate(self, emlinewave, *lineargs):
        """Evaluate the emission-line model.

        """ 
        # build the emission-line model [erg/s/cm2/A, observed frame]
        log10model = self._emline_spectrum(*lineargs)

        # optionally split into cameras, resample, and convolve with the instrumental
        # resolution
        emlinemodel = []
        for ii in np.arange(len(self.npixpercamera)-1): # iterate over cameras
            ipix = np.sum(self.npixpercamera[:ii+1]) # all pixels!
            jpix = np.sum(self.npixpercamera[:ii+2])
            #_emlinemodel = resample_flux(emlinewave[ipix:jpix], 10**self.log10wave, log10model)
            #_emlinemodel = trapz_rebin(10**self.log10wave, log10model, emlinewave[ipix:jpix])
            try:
                _emlinemodel = trapz_rebin(10**self.log10wave, log10model, emlinewave[ipix:jpix])
            except:
                _emlinemodel = resample_flux(emlinewave[ipix:jpix], 10**self.log10wave, log10model)

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
    def __init__(self, nball=10, chi2fail=1e6, minwave=3000.0,
                 maxwave=10000.0, pixkms=10.0):
        """Write me.
        
        """
        from astropy.modeling import fitting

        super(EMLineFit, self).__init__()
        
        self.nball = nball
        self.chi2fail = chi2fail

        self.pixkms = 10.0 # pixel size for internal wavelength array [km/s]
        self.dlogwave = pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        self.log10wave = np.arange(np.log10(minwave), np.log10(maxwave), self.dlogwave)
        #self.log10wave = np.arange(np.log10(emlinewave.min()), np.log10(emlinewave.max()), dlogwave)

    def init_output(self, linetable, nobj=1, chi2_default=1e6):
        """Initialize the output data table for this class.

        """
        import astropy.units as u
        from astropy.table import Table, Column
        
        out = Table()
        if False:
            out.add_column(Column(name='DN4000_NOLINES', length=nobj, dtype='f4'))
            
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
            npixpercamera=npixpercamera,
            log10wave=self.log10wave)
        
        lineargs = [fastspecfit_table[linename.upper()] for linename in EMLine.param_names]#[4:]]
        
        _emlinemodel = EMLine.evaluate(np.hstack(specwave), *lineargs)

        # unpack the per-camera spectra
        emlinemodel = []
        npix = np.hstack([0, npixpercamera])
        for ii in np.arange(len(specwave)): # iterate over cameras
            ipix = np.sum(npix[:ii+1])
            jpix = np.sum(npix[:ii+2])
            emlinemodel.append(_emlinemodel[ipix:jpix])

        return emlinemodel
    
    def fit(self, data, continuummodel, smooth_continuum, synthphot=True,
            maxiter=1000, accuracy=1e-2):
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

        npixpercamera = [len(gw) for gw in data['wave']] # all pixels
        npixpercam = np.hstack([0, npixpercamera])

        emlinevar, emlinegood = ivar2var(emlineivar)
        weights = np.sqrt(emlineivar)
        emlinebad = np.logical_not(emlinegood)
        if np.count_nonzero(emlinebad) > 0:
            weights[emlinebad] = 10*np.max(weights[emlinegood]) # 1e16 # ???

        self.EMLineModel = EMLineModel(redshift=redshift,
                                       emlineR=data['res'],
                                       npixpercamera=npixpercamera,
                                       log10wave=self.log10wave)
        
        fitter = FastLevMarLSQFitter(self.EMLineModel)

        t0 = time.time()        
        bestfit = fitter(self.EMLineModel, emlinewave, emlineflux, weights=weights,
                         maxiter=maxiter, acc=accuracy)
        log.info('Line-fitting took {:.2f} sec (niter={})'.format(time.time()-t0, fitter.fit_info['nfev']))

        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_output(self.EMLineModel.linetable)

        # Populate the output table. First do a pass through the sigma
        # parameters. If sigma is zero, restore the default value and if the
        # amplitude is still at its default (or its upper bound!), it means the
        # line wasn't fit or the fit failed (right??), so set it to zero.
        for linename in self.linetable['name'].data:
            if not hasattr(bestfit, '{}_amp'.format(linename)): # line not fitted
                continue
            
            amp = getattr(bestfit, '{}_amp'.format(linename))
            sigma = getattr(bestfit, '{}_sigma'.format(linename))
            vshift = getattr(bestfit, '{}_vshift'.format(linename))

            # drop the line if:
            #  sigma = 0, amp = default, or amp = max bound (not optimized!)
            if ((sigma.value <= sigma.bounds[0]) or (sigma.value >= sigma.bounds[1]) or
                (amp.value == amp.default) or (amp.value >= amp.bounds[1]) or
                (amp.value <= amp.bounds[0])):
                setattr(bestfit, '{}_amp'.format(linename), 0.0)
                setattr(bestfit, '{}_sigma'.format(linename), sigma.default)
                setattr(bestfit, '{}_vshift'.format(linename), vshift.default)

        # special case the tied doublets
        if bestfit.nii_6584_amp == 0.0 and bestfit.nii_6548_amp != 0.0:
            bestfit.nii_6548_amp = 0.0
        if bestfit.oiii_5007_amp == 0.0 and bestfit.oiii_4959_amp != 0.0:
            bestfit.oiii_4959_amp = 0.0
        if bestfit.oii_3726_amp == 0.0 and bestfit.oii_3726_amp != 0.0:
            bestfit.oii_3726_amp = 0.0

        # Now fill the output table.
        for pp in bestfit.param_names:
            pinfo = getattr(bestfit, pp)
            val = pinfo.value
            result[pinfo.name.upper()] = val
                
        emlinemodel = bestfit(emlinewave)
        chi2 = self.chi2(bestfit, emlinewave, emlineflux, emlineivar)

        # Synthesize photometry from the best-fitting model (continuum+emission lines).
        if synthphot:
            if data['photsys'] == 'S':
                filters = self.decam
            else:
                filters = self.bassmzls

            # The wavelengths overlap between the cameras a bit...
            srt = np.argsort(emlinewave)
            padflux, padwave = filters.pad_spectrum((continuummodelflux+smoothcontinuummodelflux+emlinemodel)[srt],
                                                    emlinewave[srt], method='edge')
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

        # measure DN(4000) without the emission lines
        if False:
            dn4000_nolines, _ = self.get_dn4000(emlinewave, specflux_nolines, redshift=redshift)
            result['DN4000_NOLINES'] = dn4000_nolines

        # get continuum fluxes, EWs, and upper limits
        verbose = False

        balmer_sigmas, forbidden_sigmas, broad_sigmas = [], [], []
        balmer_redshifts, forbidden_redshifts, broad_redshifts = [], [], []
        for oneline in self.EMLineModel.linetable[self.EMLineModel.inrange]:

            linename = oneline['name'].upper()
            linez = redshift + result['{}_VSHIFT'.format(linename)][0] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)

            linesigma = result['{}_SIGMA'.format(linename)][0] # [km/s]
            log10sigma = linesigma / C_LIGHT / np.log(10)      # line-width [log-10 Angstrom]            

            # number of pixels, chi2, and boxcar integration
            lineindx = np.where((emlinewave > (linezwave - 2.*linesigma * linezwave / C_LIGHT)) *
                                (emlinewave < (linezwave + 2.*linesigma * linezwave / C_LIGHT)) *
                                (emlineivar > 0))[0]
            npix = len(lineindx)
            #if 'MGII' in linename:
            #    pdb.set_trace()
                
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
                result['{}_AMP'.format(linename)] = 0.0 # overwrite
                npix, chi2, boxflux, boxflux_ivar = 0, 1e6, 0.0, 0.0
                
            result['{}_NPIX'.format(linename)] = npix
            result['{}_CHI2'.format(linename)] = chi2
            result['{}_BOXFLUX'.format(linename)] = boxflux
            #result['{}_BOXFLUX_IVAR'.format(linename)] = boxflux_ivar

            #if '5876' in linename:
            #    pdb.set_trace()
            
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
                              (emlinewave < (linezwave - 3.*linesigma * linezwave / C_LIGHT)) *
                              (emlineivar > 0))[0]
                              #(emlinemodel == 0))[0]
            indxhi = np.where((emlinewave < (linezwave + 10*linesigma * linezwave / C_LIGHT)) *
                              (emlinewave > (linezwave + 3.*linesigma * linezwave / C_LIGHT)) *
                              (emlineivar > 0))[0]
                              #(emlinemodel == 0))[0]
            indx = np.hstack((indxlo, indxhi))

            cmed, civar = 0.0, 0.0
            if len(indx) > 10: # require at least XX pixels to get the continuum level
                #_, cmed, csig = sigma_clipped_stats(specflux_nolines[indx], sigma=3.0)
                clipflux, _, _ = sigmaclip(specflux_nolines[indx], low=3, high=3)
                # corner case: if a portion of a camera is masked
                if len(clipflux) > 0:
                    cmed, csig = np.median(clipflux), np.std(clipflux)
                if csig > 0:
                    civar = (np.sqrt(len(indx)) / csig)**2
                    #result['{}_AMP_IVAR'.format(linename)] = 1 / csig**2

            #if '3726' in linename:
            #    pdb.set_trace()
            result['{}_CONT'.format(linename)] = cmed
            result['{}_CONT_IVAR'.format(linename)] = civar

            ew, ewivar, fluxlimit, ewlimit = 0.0, 0.0, 0.0, 0.0
            if result['{}_CONT'.format(linename)] != 0.0 and result['{}_CONT_IVAR'.format(linename)] != 0.0:
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
            #pdb.set_trace()
            
        return result
    
    def qa_fastspec(self, data, fastspec, metadata, coadd_type='healpix',
                    wavelims=(3600, 9800), outprefix=None, outdir=None):
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

        if coadd_type == 'healpix':
            title = 'Survey/Program/HealPix: {}/{}/{}, TargetID: {}'.format(
                    metadata['SURVEY'], metadata['FAPRGRM'], metadata['HPXPIXEL'], metadata['TARGETID'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['FAPRGRM'], metadata['HPXPIXEL'], metadata['TARGETID']))
        elif coadd_type == 'cumulative':
            title = 'Tile/Coadd: {}/{}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], coadd_type, metadata['TARGETID'], metadata['FIBER'])
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

        redshift = fastspec['CONTINUUM_Z']
        npixpercamera = [len(gw) for gw in data['wave']] # all pixels
        npixpercam = np.hstack([0, npixpercamera])
        
        # rebuild the best-fitting spectroscopic and photometric models
        continuum, _ = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift, 
                                     specwave=data['wave'], specres=data['res'],
                                     cameras=data['cameras'],
                                     AV=fastspec['CONTINUUM_AV'],
                                     vdisp=fastspec['CONTINUUM_VDISP'],
                                     coeff=fastspec['CONTINUUM_COEFF'],
                                     synthphot=False)
        
        residuals = [data['flux'][icam] - continuum[icam] for icam in np.arange(len(data['cameras']))]
        _smooth_continuum = self.smooth_residuals(
            residuals, data['wave'], data['ivar'],
            data['linemask'], data['linepix'], data['contpix'])
        #_smooth_continuum = self.smooth_residuals(
        #    np.hstack(continuum), np.hstack(data['wave']), np.hstack(data['flux']),
        #    np.hstack(data['ivar']), np.hstack(data['linemask']))
        smooth_continuum = []
        for icam in np.arange(len(data['cameras'])): # iterate over cameras
            ipix = np.sum(npixpercam[:icam+1])
            jpix = np.sum(npixpercam[:icam+2])
            smooth_continuum.append(_smooth_continuum[ipix:jpix])

        _emlinemodel = self.emlinemodel_bestfit(data['wave'], data['res'], fastspec)

        # QA choices
        inches_wide = 20
        inches_fullspec = 6
        inches_perline = inches_fullspec / 2.0
        nlinepanels = 5

        nline = len(set(self.linetable['plotgroup']))
        nlinerows = np.int(np.ceil(nline / nlinepanels))
        nrows = 2 + nlinerows

        height_ratios = np.hstack([1, 1, [0.5]*nlinerows])

        fig = plt.figure(figsize=(inches_wide, 2*inches_fullspec + inches_perline*nlinerows))
        gs = fig.add_gridspec(nrows, nlinepanels, height_ratios=height_ratios)

        # full spectrum + best-fitting continuum model
        bigax1 = fig.add_subplot(gs[0, :])

        leg = {
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            'dv_balmer': '$\Delta v_{{\\rm H+He}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['BALMER_Z']-redshift)),
            'dv_forbid': '$\Delta v_{{\\rm forbid}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['FORBIDDEN_Z']-redshift)),
            'dv_broad': '$\Delta v_{{\\rm broad}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['BROAD_Z']-redshift)),
            #'zbalmer': '$z_{{\\rm Balmer}}$={:.6f}'.format(fastspec['BALMER_Z']),
            #'zforbidden': '$z_{{\\rm forbidden}}$={:.6f}'.format(fastspec['FORBIDDEN_Z']),
            #'zbroad': '$z_{{\\rm MgII}}$={:.6f}'.format(fastspec['BROAD_Z']),
            'sigma_balmer': '$\sigma_{{\\rm H+He}}$={:.1f} km/s'.format(fastspec['BALMER_SIGMA']),
            'sigma_forbid': '$\sigma_{{\\rm forbid}}$={:.1f} km/s'.format(fastspec['FORBIDDEN_SIGMA']),
            'sigma_broad': '$\sigma_{{\\rm broad}}$={:.1f} km/s'.format(fastspec['BROAD_SIGMA']),
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

        legxpos, legypos, legfntsz = 0.98, 0.94, 20
        bbox = dict(boxstyle='round', facecolor='gray', alpha=0.25)
        
        for ii in np.arange(len(data['cameras'])): # iterate over cameras
            sigma, _ = ivar2var(data['ivar'][ii], sigma=True)

            bigax1.fill_between(data['wave'][ii], data['flux'][ii]-sigma,
                             data['flux'][ii]+sigma, color=col1[ii])
            #bigax1.plot(data['wave'][ii], continuum[ii], color=col2[ii], alpha=1.0)#, color='k')
            #bigax1.plot(data['wave'][ii], continuum_nodust[ii], alpha=0.5, color='k')
            #bigax1.plot(data['wave'][ii], smooth_continuum[ii], color='gray')#col3[ii])#, alpha=0.3, lw=2)#, color='k')
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
            #print(ymin, ymax)
            #if ii == 2:
            #    pdb.set_trace()

        bigax1.plot(np.hstack(data['wave']), _smooth_continuum, color='gray')#col3[ii])#, alpha=0.3, lw=2)#, color='k')

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
        for ii in np.arange(len(data['cameras'])): # iterate over cameras
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
            #print(ymin, ymax)

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
            linenames.append(self.linetable['nicename'][I[0]].replace('-', ' '))
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

            wmin = (meanwave - deltawave) * (1+redshift) - 5 * sig * meanwave * (1+redshift) / C_LIGHT
            wmax = (meanwave + deltawave) * (1+redshift) + 5 * sig * meanwave * (1+redshift) / C_LIGHT
            #print(linename, wmin, wmax)

            # iterate over cameras
            for ii in np.arange(len(data['cameras'])): # iterate over cameras
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
                    
                xx.text(0.04, 0.89, linename, ha='left', va='center',
                        transform=xx.transAxes, fontsize=16)
                
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
