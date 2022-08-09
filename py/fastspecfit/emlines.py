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
from astropy.modeling.fitting import LevMarLSQFitter

from desispec.interpolation import resample_flux
from fastspecfit.util import trapz_rebin, C_LIGHT
from fastspecfit.continuum import ContinuumTools
from desiutil.log import get_logger, DEBUG
log = get_logger()#DEBUG)

def read_emlines():
    """Read the set of emission lines of interest.

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

#def _tie_oii_blue_amp(model):
#    """
#    [O2] (2-->1): airwave: 3728.8145 vacwave: 3729.8750 emissivity: 1.948e-21
#    [O2] (3-->1): airwave: 3726.0322 vacwave: 3727.0919 emissivity: 1.444e-21
#    """
#    return model.oii_3729_amp / 1.3490
    
def _tie_oii_red_amp(model):
    """
    [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 3.229e-22
    [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.692e-22
    """
    return model.oii_7320_amp / 1.9085
    
#def _tie_siii_amp(model):
#    """
#    [S3] (4-->2): airwave: 9068.6140 vacwave: 9071.1034 emissivity: 8.10672e-21
#    [S3] (4-->3): airwave: 9530.6129 vacwave: 9533.2274 emissivity: 2.00141e-20
#    """
#    return model.siii_9532_amp / 2.4688
    
#def _tie_hbeta_sigma(model):
#    return model.hbeta_sigma
#def _tie_hbeta_vshift(model):
#    return model.hbeta_vshift
#
#def _tie_hbeta_broad_sigma(model):
#    return model.hbeta_broad_sigma
#def _tie_hbeta_broad_vshift(model):
#    return model.hbeta_broad_vshift

def _tie_halpha_sigma(model):
    return model.halpha_sigma
def _tie_halpha_vshift(model):
    return model.halpha_vshift

def _tie_halpha_broad_sigma(model):
    return model.halpha_broad_sigma
def _tie_halpha_broad_vshift(model):
    return model.halpha_broad_vshift

def _tie_oiii_5007_sigma(model):
    return model.oiii_5007_sigma
def _tie_oiii_5007_vshift(model):
    return model.oiii_5007_vshift

def _tie_mgii_2803_sigma(model):
    return model.mgii_2803_sigma
def _tie_mgii_2803_vshift(model):
    return model.mgii_2803_vshift

#def _tie_siliii_1892_sigma(model):
#    return model.siliii_1892_sigma
#def _tie_siliii_1892_vshift(model):
#    return model.siliii_1892_vshift

def _tie_ciii_1908_sigma(model):
    return model.ciii_1908_sigma
def _tie_ciii_1908_vshift(model):
    return model.ciii_1908_vshift

#def _tie_civ_1550_sigma(model):
#    return model.civ_1550_sigma
#def _tie_civ_1550_vshift(model):
#    return model.civ_1550_vshift

#def _tie_siliv_1403_sigma(model):
#    return model.siliv_1403_sigma
#def _tie_siliv_1403_vshift(model):
#    return model.siliv_1403_vshift

#def _tie_nv_1243_sigma(model):
#    return model.nv_1243_sigma
#def _tie_nv_1243_vshift(model):
#    return model.nv_1243_vshift

def _count_free_parameters(model):
    nfree = np.count_nonzero([model.tied[tiedkey] is False and model.fixed[fixedkey] is not False
                              for tiedkey, fixedkey in zip(model.tied.keys(), model.fixed.keys())])
    return nfree

def _tie_lines(model):
    """Tie emission lines together using sensible constrains based on galaxy & QSO
    physics.

    """
    # The following parameters are used as "tied" functions **even if the line
    # is out of the wavelength range** so they should never be fixed!
    neverfix = ['halpha_broad_sigma', 'halpha_broad_vshift',
                'mgii_2803_sigma', 'mgii_2803_vshift',
                'oiii_5007_sigma', 'oiii_5007_vshift']

    # Set the doublet ratios to fixed if both lines are out of range.
    for doublet in model.doublet_names:
        if doublet == 'mgii_2796':
            iline1 = model.linetable['name'] == 'mgii_2796'
            iline2 = model.linetable['name'] == 'mgii_2803'
            if not model.inrange[iline1] and not model.inrange[iline2]:
                model.mgii_doublet_ratio.fixed = True
        elif doublet == 'oii_3726':
            iline1 = model.linetable['name'] == 'oii_3726'
            iline2 = model.linetable['name'] == 'oii_3729'
            if not model.inrange[iline1] and not model.inrange[iline2]:
                model.oii_doublet_ratio.fixed = True
        elif doublet == 'sii_6731':
            iline1 = model.linetable['name'] == 'sii_6716'
            iline2 = model.linetable['name'] == 'sii_6731'
            if not model.inrange[iline1] and not model.inrange[iline2]:
                model.sii_doublet_ratio.fixed = True
        else:
            errmsg = 'Unrecognized doublet -- need to add a case in _tie_lines.'
            log.critical(errmsg)
            raise ValueError(errmsg)

    # Now apply additional constraints **for all lines, not just those in
    # range!**
    for iline in np.arange(len(model.linetable)):

        linename = model.linetable['name'][iline]

        # Tie broad Balmer+He lines to Halpha.
        if model.linetable['isbalmer'][iline] and model.linetable['isbroad'][iline] and linename != 'halpha_broad':
            getattr(model, '{}_sigma'.format(linename)).tied = _tie_halpha_broad_sigma
            getattr(model, '{}_vshift'.format(linename)).tied = _tie_halpha_broad_vshift
        #if model.inrange[iline] and model.linetable['isbalmer'][iline] and model.linetable['isbroad'][iline] and linename != 'hbeta_broad':
        #    getattr(model, '{}_sigma'.format(linename)).tied = _tie_hbeta_broad_sigma
        #    getattr(model, '{}_vshift'.format(linename)).tied = _tie_hbeta_broad_vshift
            
        # other broad lines
        if (model.linetable['isbalmer'][iline] == False) and model.linetable['isbroad'][iline]:
            if linename == 'siliii_1892':
                getattr(model, '{}_sigma'.format(linename)).tied = _tie_ciii_1908_sigma
                getattr(model, '{}_vshift'.format(linename)).tied = _tie_ciii_1908_vshift
            if linename == 'mgii_2796':
                getattr(model, '{}_sigma'.format(linename)).tied = _tie_mgii_2803_sigma
                getattr(model, '{}_vshift'.format(linename)).tied = _tie_mgii_2803_vshift
            
        # forbidden lines and narrow Balmer & He lines
        if model.linetable['isbroad'][iline] == False and linename != 'oiii_5007':
            getattr(model, '{}_sigma'.format(linename)).tied = _tie_oiii_5007_sigma
            getattr(model, '{}_vshift'.format(linename)).tied = _tie_oiii_5007_vshift

        # Now, for lines not in the wavelength range, set all their parameters
        # to zero and fixed **except those in the neverfix list**, which are
        # used as .tied functions below. This step should appear last.
        if not model.inrange[iline]:
            if hasattr(model, '{}_amp'.format(linename)): # doublet ratios have no amplitude
               setattr(model, '{}_amp'.format(linename), 0.0)
               getattr(model, '{}_amp'.format(linename)).fixed = True
            if '{}_sigma'.format(linename) not in neverfix:
                setattr(model, '{}_sigma'.format(linename), 0.0)
                getattr(model, '{}_sigma'.format(linename)).fixed = True
            if '{}_vshift'.format(linename) not in neverfix:
                setattr(model, '{}_vshift'.format(linename), 0.0)
                getattr(model, '{}_vshift'.format(linename)).fixed = True

    #nfree = model.count_free_parameters()
    #log.info('Number of free parameters to optimize: {}'.format(nfree))

    #for pp in model.param_names:
    #    print(getattr(model, pp))
    #pdb.set_trace()
    
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
    vmaxshift_broad = 2500.0 # 3000.0
    initvshift = 1.0

    minsigma = 1.0
    minsigma_broad = 1.0
    maxsigma_narrow = 500.0
    maxsigma_broad = 1e4
    minsigma_balmer_broad = minsigma

    # Be very careful about changing the default broad line-sigma. Smaller
    # values like 1500 km/s (which is arguably more sensible) can lead to
    # low-amplitude broad lines in a bunch of normal star-forming galaxy
    # spectra. (They act to "suck up" local continuum variations.) Also recall
    # that if it's well-measured, we use the initial line-sigma in
    # continuum.ContinuumFit.estimate_linesigma, which is a better initial
    # guess.
    initsigma_narrow = 75.0 # 260.0 # 75.0
    initsigma_broad = 3000.0  

    # default line-sigma for computing upper limits
    limitsigma_narrow = 75.0
    limitsigma_broad = 1200.0 

    #minamp = 0.0
    minamp = -1e2
    maxamp = +1e5
    minamp_balmer_broad = minamp # 0.0

    # Specialized parameters on the MgII, [OII], and [SII] doublet ratios. See
    # https://github.com/desihub/fastspecfit/issues/39. Be sure to set
    # self.doublet_names, below, and also note that any change in the order of
    # these lines has to be handled in _emline_spectrum!
    mgii_doublet_ratio = Parameter(name='mgii_doublet_ratio', default=0.5, bounds=[0.01, 10.0]) # MgII 2796/2803
    oii_doublet_ratio = Parameter(name='oii_doublet_ratio', default=0.74, bounds=[0.5, 1.5]) # [0.66, 1.4]) # [OII] 3726/3729
    sii_doublet_ratio = Parameter(name='sii_doublet_ratio', default=0.74, bounds=[0.5, 1.5]) # [0.67, 1.2]) # [SII] 6731/6716
    mgii_2803_amp = Parameter(name='mgii_2803_amp', default=3.0, bounds=[minamp, maxamp])
    oii_3729_amp = Parameter(name='oii_3729_amp', default=1.0, bounds=[minamp, maxamp])
    sii_6716_amp = Parameter(name='sii_6716_amp', default=1.0, bounds=[minamp, maxamp])

    # Fragile because the lines are hard-coded--
    oi_1304_amp = Parameter(name='oi_1304_amp', default=3.0, bounds=[minamp, maxamp])
    siliv_1396_amp = Parameter(name='siliv_1396_amp', default=3.0, bounds=[minamp, maxamp])
    civ_1549_amp = Parameter(name='civ_1549_amp', default=3.0, bounds=[minamp, maxamp]) # bounds=[0, maxamp])
    siliii_1892_amp = Parameter(name='siliii_1892_amp', default=3.0, bounds=[minamp, maxamp])
    ciii_1908_amp = Parameter(name='ciii_1908_amp', default=3.0, bounds=[minamp, maxamp])
    #mgii_2800_amp = Parameter(name='mgii_2800_amp', default=3.0, bounds=[minamp, maxamp])
    #mgii_2796_amp = Parameter(name='mgii_2796_amp', default=3.0, bounds=[minamp, maxamp])
    #mgii_2803_amp = Parameter(name='mgii_2803_amp', default=3.0, bounds=[minamp, maxamp])
    nev_3346_amp = Parameter(name='nev_3346_amp', default=0.1, bounds=[minamp, maxamp])
    nev_3426_amp = Parameter(name='nev_3426_amp', default=0.1, bounds=[minamp, maxamp])
    neiii_3869_amp = Parameter(name='neiii_3869_amp', default=0.3, bounds=[minamp, maxamp])
    #hei_3889_amp = Parameter(name='hei_3889_amp', default=0.3, bounds=[minamp, maxamp])
    #hei_broad_3889_amp = Parameter(name='hei_broad_3889_amp', default=0.3, bounds=[minamp_balmer_broad, maxamp])
    h6_amp = Parameter(name='h6_amp', default=0.3, bounds=[minamp, maxamp])
    h6_broad_amp = Parameter(name='h6_broad_amp', default=0.3, bounds=[minamp_balmer_broad, maxamp])
    hepsilon_amp = Parameter(name='hepsilon_amp', default=0.5, bounds=[minamp, maxamp])
    hepsilon_broad_amp = Parameter(name='hepsilon_broad_amp', default=0.5, bounds=[minamp_balmer_broad, maxamp])
    hdelta_amp = Parameter(name='hdelta_amp', default=0.5, bounds=[minamp, maxamp])
    hdelta_broad_amp = Parameter(name='hdelta_broad_amp', default=0.5, bounds=[minamp_balmer_broad, maxamp])
    hgamma_amp = Parameter(name='hgamma_amp', default=0.5, bounds=[minamp, maxamp])
    hgamma_broad_amp = Parameter(name='hgamma_broad_amp', default=0.5, bounds=[minamp_balmer_broad, maxamp])
    oiii_4363_amp = Parameter(name='oiii_4363_amp', default=0.3, bounds=[minamp, maxamp])
    hei_4471_amp = Parameter(name='hei_4471_amp', default=0.3, bounds=[minamp, maxamp])
    hei_broad_4471_amp = Parameter(name='hei_broad_4471_amp', default=0.3, bounds=[minamp_balmer_broad, maxamp])
    heii_4686_amp = Parameter(name='heii_4686_amp', default=0.3, bounds=[minamp, maxamp])
    heii_broad_4686_amp = Parameter(name='heii_broad_4686_amp', default=0.3, bounds=[minamp_balmer_broad, maxamp])
    hbeta_amp = Parameter(name='hbeta_amp', default=1.0, bounds=[minamp, maxamp])
    hbeta_broad_amp = Parameter(name='hbeta_broad_amp', default=1.0, bounds=[minamp_balmer_broad, maxamp])
    oiii_4959_amp = Parameter(name='oiii_4959_amp', default=1.0, bounds=[minamp, maxamp])
    oiii_5007_amp = Parameter(name='oiii_5007_amp', default=3.0, bounds=[minamp, maxamp])
    nii_5755_amp = Parameter(name='nii_5755_amp', default=0.3, bounds=[minamp, maxamp])
    hei_5876_amp = Parameter(name='hei_5876_amp', default=0.3, bounds=[minamp, maxamp])
    hei_broad_5876_amp = Parameter(name='hei_broad_5876_amp', default=0.3, bounds=[minamp_balmer_broad, maxamp])
    oi_6300_amp = Parameter(name='oi_6300_amp', default=0.3, bounds=[minamp, maxamp])
    siii_6312_amp = Parameter(name='siii_6312_amp', default=0.3, bounds=[minamp, maxamp])
    nii_6548_amp = Parameter(name='nii_6548_amp', default=1.0, bounds=[minamp, maxamp])
    halpha_amp = Parameter(name='halpha_amp', default=3.0, bounds=[minamp, maxamp])
    halpha_broad_amp = Parameter(name='halpha_broad_amp', default=3.0, bounds=[minamp_balmer_broad, maxamp])
    nii_6584_amp = Parameter(name='nii_6584_amp', default=3.0, bounds=[minamp, maxamp])
    oii_7320_amp = Parameter(name='oii_7320_amp', default=1.0, bounds=[minamp, maxamp])
    oii_7330_amp = Parameter(name='oii_7330_amp', default=1.0, bounds=[minamp, maxamp])
    siii_9069_amp = Parameter(name='siii_9069_amp', default=0.3, bounds=[minamp, maxamp])
    siii_9532_amp = Parameter(name='siii_9532_amp', default=1.0, bounds=[minamp, maxamp])

    mgii_2796_vshift = Parameter(name='mgii_2796_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    oii_3726_vshift = Parameter(name='oii_3726_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    sii_6731_vshift = Parameter(name='sii_6731_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    mgii_2803_vshift = Parameter(name='mgii_2803_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    oii_3729_vshift = Parameter(name='oii_3729_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    sii_6716_vshift = Parameter(name='sii_6716_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    
    oi_1304_vshift = Parameter(name='oi_1304_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    siliv_1396_vshift = Parameter(name='siliv_1396_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    civ_1549_vshift = Parameter(name='civ_1549_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    siliii_1892_vshift = Parameter(name='siliii_1892_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    ciii_1908_vshift = Parameter(name='ciii_1908_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    #mgii_2800_vshift = Parameter(name='mgii_2800_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    #mgii_2796_vshift = Parameter(name='mgii_2796_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    #mgii_2803_vshift = Parameter(name='mgii_2803_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])    
    nev_3346_vshift = Parameter(name='nev_3346_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nev_3426_vshift = Parameter(name='nev_3426_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    neiii_3869_vshift = Parameter(name='neiii_3869_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    #hei_3889_vshift = Parameter(name='hei_3889_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    #hei_broad_3889_vshift = Parameter(name='hei_broad_3889_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    h6_vshift = Parameter(name='h6_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    h6_broad_vshift = Parameter(name='h6_broad_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    hepsilon_vshift = Parameter(name='hepsilon_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hepsilon_broad_vshift = Parameter(name='hepsilon_broad_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    hdelta_vshift = Parameter(name='hdelta_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hdelta_broad_vshift = Parameter(name='hdelta_broad_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    hgamma_vshift = Parameter(name='hgamma_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hgamma_broad_vshift = Parameter(name='hgamma_broad_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    oiii_4363_vshift = Parameter(name='oiii_4363_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_4471_vshift = Parameter(name='hei_4471_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_broad_4471_vshift = Parameter(name='hei_broad_4471_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    heii_4686_vshift = Parameter(name='heii_4686_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    heii_broad_4686_vshift = Parameter(name='heii_broad_4686_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    hbeta_vshift = Parameter(name='hbeta_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hbeta_broad_vshift = Parameter(name='hbeta_broad_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    oiii_4959_vshift = Parameter(name='oiii_4959_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oiii_5007_vshift = Parameter(name='oiii_5007_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nii_5755_vshift = Parameter(name='nii_5755_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_5876_vshift = Parameter(name='hei_5876_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    hei_broad_5876_vshift = Parameter(name='hei_broad_5876_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    oi_6300_vshift = Parameter(name='oi_6300_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    siii_6312_vshift = Parameter(name='siii_6312_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    nii_6548_vshift = Parameter(name='nii_6548_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    halpha_vshift = Parameter(name='halpha_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    halpha_broad_vshift = Parameter(name='halpha_broad_vshift', default=initvshift, bounds=[-vmaxshift_broad, +vmaxshift_broad])
    nii_6584_vshift = Parameter(name='nii_6584_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oii_7320_vshift = Parameter(name='oii_7320_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    oii_7330_vshift = Parameter(name='oii_7330_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    siii_9069_vshift = Parameter(name='siii_9069_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])
    siii_9532_vshift = Parameter(name='siii_9532_vshift', default=initvshift, bounds=[-vmaxshift_narrow, +vmaxshift_narrow])

    mgii_2796_sigma = Parameter(name='mgii_2796_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    oii_3726_sigma = Parameter(name='oii_3726_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    sii_6731_sigma = Parameter(name='sii_6731_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    mgii_2803_sigma = Parameter(name='mgii_2803_sigma', default=initsigma_broad, bounds=[minsigma, maxsigma_broad])
    oii_3729_sigma = Parameter(name='oii_3729_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    sii_6716_sigma = Parameter(name='sii_6716_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])

    oi_1304_sigma = Parameter(name='oi_1304_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    siliv_1396_sigma = Parameter(name='siliv_1396_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    civ_1549_sigma = Parameter(name='civ_1549_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    siliii_1892_sigma = Parameter(name='siliii_1892_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    ciii_1908_sigma = Parameter(name='ciii_1908_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    #mgii_2800_sigma = Parameter(name='mgii_2800_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    #mgii_2796_sigma = Parameter(name='mgii_2796_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])
    #mgii_2803_sigma = Parameter(name='mgii_2803_sigma', default=initsigma_broad, bounds=[minsigma_broad, maxsigma_broad])    
    nev_3346_sigma = Parameter(name='nev_3346_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    nev_3426_sigma = Parameter(name='nev_3426_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    neiii_3869_sigma = Parameter(name='neiii_3869_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    #hei_3889_sigma = Parameter(name='hei_3889_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    #hei_broad_3889_sigma = Parameter(name='hei_broad_3889_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    h6_sigma = Parameter(name='h6_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    h6_broad_sigma = Parameter(name='h6_broad_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    hepsilon_sigma = Parameter(name='hepsilon_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hepsilon_broad_sigma = Parameter(name='hepsilon_broad_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    hdelta_sigma = Parameter(name='hdelta_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hdelta_broad_sigma = Parameter(name='hdelta_broad_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    hgamma_sigma = Parameter(name='hgamma_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hgamma_broad_sigma = Parameter(name='hgamma_broad_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    oiii_4363_sigma = Parameter(name='oiii_4363_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_4471_sigma = Parameter(name='hei_4471_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_broad_4471_sigma = Parameter(name='hei_broad_4471_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    heii_4686_sigma = Parameter(name='heii_4686_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    heii_broad_4686_sigma = Parameter(name='heii_broad_4686_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    hbeta_sigma = Parameter(name='hbeta_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hbeta_broad_sigma = Parameter(name='hbeta_broad_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    oiii_4959_sigma = Parameter(name='oiii_4959_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    oiii_5007_sigma = Parameter(name='oiii_5007_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    nii_5755_sigma = Parameter(name='nii_5755_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_5876_sigma = Parameter(name='hei_5876_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    hei_broad_5876_sigma = Parameter(name='hei_broad_5876_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    oi_6300_sigma = Parameter(name='oi_6300_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    siii_6312_sigma = Parameter(name='siii_6312_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    nii_6548_sigma = Parameter(name='nii_6548_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    halpha_sigma = Parameter(name='halpha_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    halpha_broad_sigma = Parameter(name='halpha_broad_sigma', default=initsigma_broad, bounds=[minsigma_balmer_broad, maxsigma_broad])
    nii_6584_sigma = Parameter(name='nii_6584_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    oii_7320_sigma = Parameter(name='oii_7320_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    oii_7330_sigma = Parameter(name='oii_7330_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    siii_9069_sigma = Parameter(name='siii_9069_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])
    siii_9532_sigma = Parameter(name='siii_9532_sigma', default=initsigma_narrow, bounds=[minsigma, maxsigma_narrow])

    nii_6548_amp.tied = _tie_nii_amp
    oiii_4959_amp.tied = _tie_oiii_amp
    oii_7330_amp.tied = _tie_oii_red_amp
    #siii_9069_amp.tied = _tie_siii_amp
    
    def __init__(self,
                 mgii_doublet_ratio=mgii_doublet_ratio.default,
                 oii_doublet_ratio=oii_doublet_ratio.default,
                 sii_doublet_ratio=sii_doublet_ratio.default,
                 mgii_2803_amp=mgii_2803_amp.default,
                 oii_3729_amp=oii_3729_amp.default,
                 sii_6716_amp=sii_6716_amp.default,
                     
                 oi_1304_amp=oi_1304_amp.default,
                 siliv_1396_amp=siliv_1396_amp.default,
                 civ_1549_amp=civ_1549_amp.default,
                 siliii_1892_amp=siliii_1892_amp.default,
                 ciii_1908_amp=ciii_1908_amp.default,
                 #mgii_2800_amp=mgii_2800_amp.default,
                 #mgii_2796_amp=mgii_2796_amp.default,
                 #mgii_2803_amp=mgii_2803_amp.default,
                 nev_3346_amp=nev_3346_amp.default,
                 nev_3426_amp=nev_3426_amp.default,
                 neiii_3869_amp=neiii_3869_amp.default,
                 #hei_3889_amp=hei_3889_amp.default,
                 #hei_broad_3889_amp=hei_broad_3889_amp.default,
                 h6_amp=h6_amp.default,
                 h6_broad_amp=h6_broad_amp.default,
                 hepsilon_amp=hepsilon_amp.default,
                 hepsilon_broad_amp=hepsilon_broad_amp.default,
                 hdelta_amp=hdelta_amp.default,
                 hdelta_broad_amp=hdelta_broad_amp.default,
                 hgamma_amp=hgamma_amp.default,
                 hgamma_broad_amp=hgamma_broad_amp.default,
                 oiii_4363_amp=oiii_4363_amp.default,
                 hei_4471_amp=hei_4471_amp.default,
                 hei_broad_4471_amp=hei_broad_4471_amp.default,
                 heii_4686_amp=heii_4686_amp.default,
                 heii_broad_4686_amp=heii_broad_4686_amp.default,
                 hbeta_amp=hbeta_amp.default,
                 hbeta_broad_amp=hbeta_broad_amp.default,
                 oiii_4959_amp=oiii_4959_amp.default,
                 oiii_5007_amp=oiii_5007_amp.default,
                 nii_5755_amp=nii_5755_amp.default,
                 hei_5876_amp=hei_5876_amp.default,
                 hei_broad_5876_amp=hei_broad_5876_amp.default,
                 oi_6300_amp=oi_6300_amp.default,
                 siii_6312_amp=siii_6312_amp.default,
                 nii_6548_amp=nii_6548_amp.default,
                 halpha_amp=halpha_amp.default,
                 halpha_broad_amp=halpha_broad_amp.default,
                 nii_6584_amp=nii_6584_amp.default,
                 oii_7320_amp=oii_7320_amp.default,
                 oii_7330_amp=oii_7330_amp.default,
                 siii_9069_amp=siii_9069_amp.default,
                 siii_9532_amp=siii_9532_amp.default,

                 #mgii_2796_vshift=mgii_2796_vshift.default,
                 oii_3726_vshift=oii_3726_vshift.default,
                 sii_6731_vshift=sii_6731_vshift.default,
                 #mgii_2803_vshift=mgii_2803_vshift.default,
                 oii_3729_vshift=oii_3729_vshift.default,
                 sii_6716_vshift=sii_6716_vshift.default,

                 oi_1304_vshift=oi_1304_vshift.default,
                 siliv_1396_vshift=siliv_1396_vshift.default,
                 civ_1549_vshift=civ_1549_vshift.default,
                 siliii_1892_vshift=siliii_1892_vshift.default,
                 ciii_1908_vshift=ciii_1908_vshift.default,
                 #mgii_2800_vshift=mgii_2800_vshift.default,
                 #mgii_2796_vshift=mgii_2796_vshift.default,
                 #mgii_2803_vshift=mgii_2803_vshift.default,
                 nev_3346_vshift=nev_3346_vshift.default,
                 nev_3426_vshift=nev_3426_vshift.default,
                 neiii_3869_vshift=neiii_3869_vshift.default,
                 #hei_3889_vshift=hei_3889_vshift.default,
                 #hei_broad_3889_vshift=hei_broad_3889_vshift.default,
                 h6_vshift=h6_vshift.default,
                 h6_broad_vshift=h6_broad_vshift.default,
                 hepsilon_vshift=hepsilon_vshift.default,
                 hepsilon_broad_vshift=hepsilon_broad_vshift.default,
                 hdelta_vshift=hdelta_vshift.default,
                 hdelta_broad_vshift=hdelta_broad_vshift.default,
                 hgamma_vshift=hgamma_vshift.default,
                 hgamma_broad_vshift=hgamma_broad_vshift.default,
                 oiii_4363_vshift=oiii_4363_vshift.default,
                 hei_4471_vshift=hei_4471_vshift.default,
                 hei_broad_4471_vshift=hei_broad_4471_vshift.default,
                 heii_4686_vshift=heii_4686_vshift.default,
                 heii_broad_4686_vshift=heii_broad_4686_vshift.default,
                 hbeta_vshift=hbeta_vshift.default,
                 hbeta_broad_vshift=hbeta_broad_vshift.default,
                 oiii_4959_vshift=oiii_4959_vshift.default,
                 oiii_5007_vshift=oiii_5007_vshift.default,
                 nii_5755_vshift=nii_5755_vshift.default,
                 hei_5876_vshift=hei_5876_vshift.default,
                 hei_broad_5876_vshift=hei_broad_5876_vshift.default,
                 oi_6300_vshift=oi_6300_vshift.default,
                 siii_6312_vshift=siii_6312_vshift.default,
                 nii_6548_vshift=nii_6548_vshift.default,
                 halpha_vshift=halpha_vshift.default,
                 halpha_broad_vshift=halpha_broad_vshift.default,
                 nii_6584_vshift=nii_6584_vshift.default,
                 oii_7320_vshift=oii_7320_vshift.default,
                 oii_7330_vshift=oii_7330_vshift.default,
                 siii_9069_vshift=siii_9069_vshift.default,
                 siii_9532_vshift=siii_9532_vshift.default,

                 #mgii_2796_sigma=mgii_2796_sigma.default,
                 oii_3726_sigma=oii_3726_sigma.default,
                 sii_6731_sigma=sii_6731_sigma.default,
                 #mgii_2803_sigma=mgii_2803_sigma.default,
                 oii_3729_sigma=oii_3729_sigma.default,
                 sii_6716_sigma=sii_6716_sigma.default,

                 oi_1304_sigma=oi_1304_sigma.default,
                 siliv_1396_sigma=siliv_1396_sigma.default,
                 civ_1549_sigma=civ_1549_sigma.default,
                 siliii_1892_sigma=siliii_1892_sigma.default,
                 ciii_1908_sigma=ciii_1908_sigma.default,
                 #mgii_2800_sigma=mgii_2800_sigma.default,
                 #mgii_2796_sigma=mgii_2796_sigma.default,
                 #mgii_2803_sigma=mgii_2803_sigma.default,
                 nev_3346_sigma=nev_3346_sigma.default,
                 nev_3426_sigma=nev_3426_sigma.default,
                 neiii_3869_sigma=neiii_3869_sigma.default,
                 #hei_3889_sigma=hei_3889_sigma.default,
                 #hei_broad_3889_sigma=hei_broad_3889_sigma.default,
                 h6_sigma=h6_sigma.default,
                 h6_broad_sigma=h6_broad_sigma.default,
                 hepsilon_sigma=hepsilon_sigma.default,
                 hepsilon_broad_sigma=hepsilon_broad_sigma.default,
                 hdelta_sigma=hdelta_sigma.default,
                 hdelta_broad_sigma=hdelta_broad_sigma.default,
                 hgamma_sigma=hgamma_sigma.default,
                 hgamma_broad_sigma=hgamma_broad_sigma.default,
                 oiii_4363_sigma=oiii_4363_sigma.default,
                 hei_4471_sigma=hei_4471_sigma.default,
                 hei_broad_4471_sigma=hei_broad_4471_sigma.default,
                 heii_4686_sigma=heii_4686_sigma.default,
                 heii_broad_4686_sigma=heii_broad_4686_sigma.default,
                 hbeta_sigma=hbeta_sigma.default,
                 hbeta_broad_sigma=hbeta_broad_sigma.default,
                 oiii_4959_sigma=oiii_4959_sigma.default,
                 oiii_5007_sigma=oiii_5007_sigma.default,
                 nii_5755_sigma=nii_5755_sigma.default,
                 hei_5876_sigma=hei_5876_sigma.default,
                 hei_broad_5876_sigma=hei_broad_5876_sigma.default,
                 oi_6300_sigma=oi_6300_sigma.default,
                 siii_6312_sigma=siii_6312_sigma.default,
                 nii_6548_sigma=nii_6548_sigma.default,
                 halpha_sigma=halpha_sigma.default,
                 halpha_broad_sigma=halpha_broad_sigma.default,
                 nii_6584_sigma=nii_6584_sigma.default,
                 oii_7320_sigma=oii_7320_sigma.default,
                 oii_7330_sigma=oii_7330_sigma.default,
                 siii_9069_sigma=siii_9069_sigma.default,
                 siii_9532_sigma=siii_9532_sigma.default,
                     
                 redshift=None,
                 wavelims=[None, None],
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

        # re-arrange the default linetable according to our parameters
        linetable = read_emlines()
        self.nline = len(linetable)
        self.linenames = np.array([linename.replace('_sigma', '') for linename in self.param_names[-self.nline:]])
        
        lineindx = np.hstack([np.where(linetable['name'] == linename)[0] for linename in self.linenames])
        self.linetable = linetable[lineindx]
        
        #self.restlinewaves = np.hstack([self.linetable['restwave'][self.linetable['name'] == linename] for linename in self.linenames])
        self.restlinewaves = self.linetable['restwave'].data

        #self.doublet_names = ['oii_3726', 'sii_6731']
        self.doublet_names = ['mgii_2796', 'oii_3726', 'sii_6731']
        self.ndoublet = len(self.doublet_names)

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

        if wavelims[0] is None:
            wavelims[0] = 10**np.min(log10wave)
        if wavelims[1] is None:
            wavelims[1] = 10**np.max(log10wave)
            
        zline = self.restlinewaves * (1 + self.redshift)
        #self.inrange = (zline > 4600) * (zline < 6300)
        self.inrange = (zline > (wavelims[0]+10)) * (zline < (wavelims[1]-10))

        # tie lines together
        _tie_lines(self)

    def _emline_spectrum(self, *lineargs):
        """Simple wrapper to build an emission-line spectrum.

        # build the emission-line model [erg/s/cm2/A, observed frame]

        """
        log10model = np.zeros_like(self.log10wave)

        # carefull parse the parameters (order matters!)
        linesplit = np.array_split(lineargs, 3) # 3 parameters per line
        linevshifts = np.hstack(linesplit[1])
        linesigmas = np.hstack(linesplit[2])

        # stupidly hard-coded
        _lineamps = np.hstack(linesplit[0])
        doublet_ratios = _lineamps[:self.ndoublet]
        doublet_amps = _lineamps[self.ndoublet:2*self.ndoublet] * doublet_ratios
        lineamps = np.hstack((doublet_amps, _lineamps[self.ndoublet:]))
        #print(doublet_ratios)

        # cut to lines in range and non-zero
        I = self.inrange * (lineamps > 0)
        if np.count_nonzero(I) == 0:
            return log10model
        
        linevshifts = linevshifts[I]
        linesigmas = linesigmas[I]
        lineamps = lineamps[I]

        #linenames = self.linenames[self.inrange]
        #lineamps = np.hstack(lineargs[0:self.nline])[self.inrange]
        #linevshifts = np.hstack(lineargs[self.nline:2*self.nline])[self.inrange]
        #linesigmas = np.hstack(lineargs[2*self.nline:])[self.inrange]

        # line-width [log-10 Angstrom] and redshifted wavelength [log-10 Angstrom]
        log10sigmas = linesigmas / C_LIGHT / np.log(10) 
        linezwaves = np.log10(self.restlinewaves[I] * (1.0 + self.redshift + linevshifts / C_LIGHT)) 

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

    def count_free_parameters(self):
        """Count the number of free (not tied, not fixed) parameters in the model."""
        nfree = np.count_nonzero([not self.tied[tiedkey] and not self.fixed[fixedkey]
                                  for tiedkey, fixedkey in zip(self.tied.keys(), self.fixed.keys())])
        return nfree

class EMLineFit(ContinuumTools):
    """Class to fit an emission-line spectrum.

    * https://docs.astropy.org/en/stable/modeling/example-fitting-constraints.html#tied
    * https://docs.astropy.org/en/stable/modeling/new-model.html
    * https://docs.astropy.org/en/stable/modeling/compound-models.html#parameters

    """
    def __init__(self, nball=10, chi2_default=0.0, minwave=3000.0,
                 maxwave=10000.0, pixkms=10.0, nolegend=False, 
                 ssptemplates=None, mapdir=None):
        """Write me.
        
        """
        from astropy.modeling import fitting

        super(EMLineFit, self).__init__(mapdir=mapdir, ssptemplates=ssptemplates)
        
        self.nball = nball
        self.chi2_default = chi2_default
        self.nolegend = nolegend

        self.emwave_pixkms = 10.0 # pixel size for internal wavelength array [km/s]
        self.dlogwave = pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        self.log10wave = np.arange(np.log10(minwave), np.log10(maxwave), self.dlogwave)
        #self.log10wave = np.arange(np.log10(emlinewave.min()), np.log10(emlinewave.max()), dlogwave)

    def init_output(self, linetable, nobj=1):
        """Initialize the output data table for this class.

        """
        import astropy.units as u
        from astropy.table import Table, Column
        
        out = Table()
        #out.add_column(Column(name='DN4000_NOLINES', length=nobj, dtype='f4'))
            
        # observed-frame photometry synthesized from the spectra
        for band in self.synth_bands:
            out.add_column(Column(name='FLUX_SYNTH_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies')) 
            #out.add_column(Column(name='FLUX_SYNTH_IVAR_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies-2'))
        # observed-frame photometry synthesized from the best-fitting continuum model fit
        for band in self.synth_bands:
            out.add_column(Column(name='FLUX_SYNTH_MODEL_{}'.format(band.upper()), length=nobj, dtype='f4', unit='nanomaggies'))

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
        
    def chi2(self, bestfit, emlinewave, emlineflux, emlineivar, continuum_model=None, return_dof=False):
        """Compute the reduced chi^2."""
        #nfree = len(bestfit.parameters)
        nfree = bestfit.count_free_parameters()
        dof = np.sum(emlineivar > 0) - nfree
        if dof > 0:
            emlinemodel = bestfit(emlinewave)
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

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table, redshift=None):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (to be used to build QA).

        """
        npixpercamera = [len(gw) for gw in specwave]
        if redshift is None:
            redshift = fastspecfit_table['CONTINUUM_Z']
        
        EMLine = EMLineModel(
            redshift=redshift, emlineR=specres,
            npixpercamera=npixpercamera,
            log10wave=self.log10wave)

        # old parameter model
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
            maxiter=5000, accuracy=1e-2, verbose=False):
        """Perform the fit minimization / chi2 minimization.
        
        EMLineModel object
        FC - ContinuumFit object

        ToDo: need to take into account the instrumental velocity width when
        computing integrated fluxes...

        """
        from astropy.table import Table, Column
        from scipy.stats import sigmaclip
        from fastspecfit.util import ivar2var
            
        def _clean_linefit(bestfit, init_amplitudes, init_sigmas):
            """Clean up line-fitting results."""
            
            # If the amplitude is still at its default (or its upper or lower
            # bound!), it means the line wasn't fit or the fit failed (right??), so
            # set it to zero. Important: need to loop over *all* lines (not just
            # those in range).
            initkeys = init_amplitudes.keys()

            for linename in self.linetable['name'].data:
                if not hasattr(bestfit, '{}_amp'.format(linename)): # skip the physical doubleta
                    continue
                amp = getattr(bestfit, '{}_amp'.format(linename))
                if amp.fixed: # line is either out of range or a suppressed line--skip it here
                    continue
                sigma = getattr(bestfit, '{}_sigma'.format(linename))
                vshift = getattr(bestfit, '{}_vshift'.format(linename))
    
                # drop the line if:
                #  sigma = 0, amp = default or amp = max bound (not optimized!) or amp =< 0
                if ((amp.value == amp.default) or
                    (amp.bounds[1] is not None and amp.value >= amp.bounds[1]) or
                    #(amp.value <= amp.bounds[0]) or
                    (amp.value <= 0) or
                    (linename in initkeys and init_amplitudes[linename] == amp.value) or 
                    (sigma.value == sigma.default)
                    #(sigma.value <= sigma.bounds[0]) or (sigma.value >= sigma.bounds[1])
                    ):
                    log.debug('Dropping {} (amp={:.3f}, sigma={:.3f}, vshift={:.3f})'.format(
                        linename, amp.value, sigma.value, vshift.value))

                    setattr(bestfit, '{}_amp'.format(linename), 0.0)
                    setattr(bestfit, '{}_sigma'.format(linename), 0.0) # sigma.default)
                    setattr(bestfit, '{}_vshift'.format(linename), 0.0) #vshift.default)

            # Now loop back through and drop Broad balmer lines that:
            #   (1) are narrower than their narrow-line counterparts;
            #   (2) have a narrow line whose amplitude is smaller than that of the broad line
            #      --> Deprecated! main-dark-32303-39628176678192981 is an example
            #          of an object where there's a broad H-alpha line but no other
            #          forbidden lines!
            
            IB = self.linetable['isbalmer'] * self.linetable['isbroad']
            for linename in self.linetable['name'][IB]:
                amp_broad = getattr(bestfit, '{}_amp'.format(linename))
                if amp_broad.fixed: # line is either out of range or a suppressed line--skip it here
                    continue
                amp = getattr(bestfit, '{}_amp'.format(linename.replace('_broad', ''))) # fragile
                sigma = getattr(bestfit, '{}_sigma'.format(linename.replace('_broad', ''))) # fragile
                sigma_broad = getattr(bestfit, '{}_sigma'.format(linename))
                
                if sigma_broad <= sigma:
                    log.debug('Dropping {} (sigma_broad={:.2f}, sigma_narrow={:.2f})'.format(
                        linename, sigma_broad.value, sigma.value))
                    #vshift = getattr(bestfit, '{}_vshift'.format(linename))
                    setattr(bestfit, '{}_amp'.format(linename), 0.0)
                    setattr(bestfit, '{}_sigma'.format(linename), 0.0) # sigma_broad.default)
                    setattr(bestfit, '{}_vshift'.format(linename), 0.0) # vshift.default)
    
                # Deprecated! -- See note above.
                #if (amp <= 0) or (amp_broad <= 0) and ((amp+amp_broad) < amp_broad):
                #    #vshift = getattr(bestfit, '{}_vshift'.format(linename))
                #    setattr(bestfit, '{}_amp'.format(linename), 0.0)
                #    setattr(bestfit, '{}_sigma'.format(linename), 0.0) # sigma_broad.default)
                #    setattr(bestfit, '{}_vshift'.format(linename), 0.0) # vshift.default)

            return bestfit

        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        redshift = data['zredrock']
        emlinewave = np.hstack(data['wave'])
        oemlineivar = np.hstack(data['ivar'])
        specflux = np.hstack(data['flux'])

        continuummodelflux = np.hstack(continuummodel)
        smoothcontinuummodelflux = np.hstack(smooth_continuum)
        emlineflux = specflux - continuummodelflux - smoothcontinuummodelflux

        npixpercamera = [len(gw) for gw in data['wave']] # all pixels
        npixpercam = np.hstack([0, npixpercamera])

        emlineivar = np.copy(oemlineivar)
        emlinevar, emlinegood = ivar2var(emlineivar, clip=1e-3)
        emlinebad = np.logical_not(emlinegood)
        if np.sum(emlinebad) > 0:
            emlineivar[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineivar[emlinegood])
            emlineflux[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineflux[emlinegood]) # ???
        weights = np.sqrt(emlineivar)

        wavelims = (np.min(emlinewave)+5, np.max(emlinewave)-5)

        self.EMLineModel = EMLineModel(redshift=redshift,
                                       wavelims=wavelims,
                                       emlineR=data['res'],
                                       npixpercamera=npixpercamera,
                                       log10wave=self.log10wave)

        # For all lines in range, do a fast update of the initial
        # line-amplitudes which especially helps with cases like
        # 39633354915582193 (tile 80613, petal 05), which has strong narrow
        # lines.
        init_amplitudes = {}
        coadd_emlineflux = np.interp(data['coadd_wave'], emlinewave, emlineflux)
        for linename, linepix in zip(data['coadd_linename'], data['coadd_linepix']):
            # skip the physical doublets
            if not hasattr(self.EMLineModel, '{}_amp'.format(linename)):
                continue
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
                    bounds = [None, mx]
            else:
                bounds = [-np.min(np.abs(coadd_emlineflux[linepix])), mx]

            setattr(self.EMLineModel, '{}_amp'.format(linename), amp)
            getattr(self.EMLineModel, '{}_amp'.format(linename)).bounds = bounds

            init_amplitudes[linename] = amp

        # Now update the linewidth but here we need to loop over *all* lines
        # (not just those in range). E.g., if H-alpha is out of range we need to
        # set its initial value correctly since other lines are tied to it
        # (e.g., main-bright-32406-39628257196245904).
        init_sigmas = {}
        for iline in self.linetable:
            linename = iline['name']
            # If sigma is zero here, it was set in _tie_lines because the line
            # is out of range and not a "neverfix" line.
            if getattr(self.EMLineModel, '{}_sigma'.format(linename)) == 0:
                continue
            init_sigmas['{}_fixed'.format(linename)] = getattr(self.EMLineModel, '{}_sigma'.format(linename)).fixed
            init_sigmas['{}_vshift'.format(linename)] = getattr(self.EMLineModel, '{}_vshift'.format(linename)).value
            if iline['isbroad']:
                if iline['isbalmer']:
                    if data['linesigma_balmer_snr'] > 0: # broad Balmer lines
                        setattr(self.EMLineModel, '{}_sigma'.format(linename), data['linesigma_balmer'])
                        #getattr(self.EMLineModel, '{}_sigma'.format(linename)).default = data['linesigma_balmer']
                        init_sigmas[linename] = data['linesigma_balmer']
                else:
                    if data['linesigma_uv_snr'] > 0: # broad UV/QSO lines
                        setattr(self.EMLineModel, '{}_sigma'.format(linename), data['linesigma_uv'])
                        #getattr(self.EMLineModel, '{}_sigma'.format(linename)).default = data['linesigma_uv']
                        init_sigmas[linename] = data['linesigma_uv']
            else:
                # prefer narrow over Balmer
                if data['linesigma_narrow_snr'] > 0:
                    setattr(self.EMLineModel, '{}_sigma'.format(linename), data['linesigma_narrow'])
                    #getattr(self.EMLineModel, '{}_sigma'.format(linename)).default = data['linesigma_narrow']
                    init_sigmas[linename] = data['linesigma_narrow']
                #elif data['linesigma_narrow_snr'] == 0 and data['linesigma_narrow_balmer'] > 0:
                #    setattr(self.EMLineModel, '{}_sigma'.format(linename), data['linesigma_balmer'])
                #    getattr(self.EMLineModel, '{}_sigma'.format(linename)).default = data['linesigma_balmer']
                #    init_sigmas[linename] = data['linesigma_balmer']

        #print('HACK!!!!!!!!!!!')
        #self.EMLineModel.halpha_sigma.tied = False
        #self.EMLineModel.nii_6584_sigma.tied = False
        #self.EMLineModel.nii_6548_sigma.tied = False
        #self.EMLineModel.mgii_2800_vshift.bounds = [None, None]
        #self.EMLineModel.mgii_2800_sigma = 5000.0
        #self.EMLineModel.hdelta_amp.bounds = [None, None]
        #self.EMLineModel.mgii_2800_amp.bounds = [None, None]
        #self.EMLineModel.mgii_2803_amp.bounds = [None, 5*5.202380598696401]
        #self.EMLineModel.mgii_2803_sigma.bounds = [None, None]
        #self.EMLineModel.mgii_2803_vshift.bounds = [None, None]

        t0 = time.time()        
        # Initial fit: set the broad-line Balmer amplitudes to zero and
        # fixed. Need to loop over all lines, not just those in range.
        IB = self.linetable['isbroad'] * self.linetable['isbalmer']
        for linename in self.linetable['name'][IB]:
            setattr(self.EMLineModel, '{}_amp'.format(linename), 0.0)
            setattr(self.EMLineModel, '{}_vshift'.format(linename), 0.0)
            setattr(self.EMLineModel, '{}_sigma'.format(linename), 0.0)
            getattr(self.EMLineModel, '{}_amp'.format(linename)).fixed = True
            getattr(self.EMLineModel, '{}_vshift'.format(linename)).fixed = True
            getattr(self.EMLineModel, '{}_sigma'.format(linename)).fixed = True
        nfree = self.EMLineModel.count_free_parameters()
        #for pp in self.EMLineModel.param_names:
        #    print(getattr(self.EMLineModel, pp))

        fitter = FastLevMarLSQFitter(self.EMLineModel)
        initfit = fitter(self.EMLineModel, emlinewave, emlineflux, weights=weights,
                         maxiter=maxiter, acc=accuracy)
        initfit = _clean_linefit(initfit, init_amplitudes, init_sigmas)
        initchi2 = self.chi2(initfit, emlinewave, emlineflux, emlineivar)
        log.info('Initial line-fitting with {} free parameters took {:.2f} sec (niter={}) with chi2={:.3f}'.format(
            nfree, time.time()-t0, fitter.fit_info['nfev'], initchi2))

        # Now try adding bround Balmer and helium lines and see if we improve
        # the chi2. First, do we have enough pixels around Halpha and Hbeta to
        # do this test?
        broadlinepix = []
        for icam in np.arange(len(data['cameras'])):
            pixoffset = int(np.sum(npixpercamera[:icam]))
            for linename, linepix in zip(data['linename'][icam], data['linepix'][icam]):
                if linename == 'halpha_broad' or linename == 'hbeta_broad' or linename == 'hgamma_broad':
                    broadlinepix.append(linepix + pixoffset)
                    #print(data['wave'][icam][linepix])

        # Require minimum XX pixels.
        if len(broadlinepix) > 0 and len(np.hstack(broadlinepix)) > 10: 
            broadlinepix = np.hstack(broadlinepix)

            t0 = time.time()
            # Restore the broad lines. Need to loop over all lines, not just
            # those in range. Do not update the narrow lines in case we hit a
            # local minium in initfit.
            IB = self.linetable['isbroad'] * self.linetable['isbalmer']
            for linename in self.linetable['name'][IB]:
                if linename in init_amplitudes.keys():
                    setattr(self.EMLineModel, '{}_amp'.format(linename), init_amplitudes[linename])
                    getattr(self.EMLineModel, '{}_amp'.format(linename)).fixed = False
                # If sigma is not in the dictionary then keep sigma and vshift equal to zero and fixed.
                if linename in init_sigmas.keys():
                    setattr(self.EMLineModel, '{}_sigma'.format(linename), init_sigmas[linename])
                    setattr(self.EMLineModel, '{}_vshift'.format(linename), init_sigmas['{}_vshift'.format(linename)])
                    getattr(self.EMLineModel, '{}_sigma'.format(linename)).fixed = init_sigmas['{}_fixed'.format(linename)]
                    getattr(self.EMLineModel, '{}_vshift'.format(linename)).fixed = init_sigmas['{}_fixed'.format(linename)]

            #for pp in self.EMLineModel.param_names:
            #    print(getattr(self.EMLineModel, pp))
            nfree = self.EMLineModel.count_free_parameters()
            
            fitter = FastLevMarLSQFitter(self.EMLineModel)
            broadfit = fitter(self.EMLineModel, emlinewave, emlineflux, weights=weights,
                              maxiter=maxiter, acc=accuracy)

            broadfit = _clean_linefit(broadfit, init_amplitudes, init_sigmas)
            broadchi2 = self.chi2(broadfit, emlinewave, emlineflux, emlineivar)
            log.info('Second (broad) line-fitting with {} free parameters took {:.2f} sec (niter={}) with chi2={:.3f}'.format(
                nfree, time.time()-t0, fitter.fit_info['nfev'], broadchi2))

            # Compare chi2 just in and around the broad lines. Need to account
            # for the different number of degrees of freedom!
            initmodel = initfit(emlinewave)
            broadmodel = broadfit(emlinewave)
            linechi2_init = np.sum(emlineivar[broadlinepix] * (emlineflux[broadlinepix] - initmodel[broadlinepix])**2) / len(broadlinepix)
            linechi2_broad = np.sum(emlineivar[broadlinepix] * (emlineflux[broadlinepix] - broadmodel[broadlinepix])**2) / len(broadlinepix)
            log.info('Chi2 with broad lines = {:.5f} and without broad lines = {:.5f}'.format(linechi2_broad, linechi2_init))

            # Make a decision! All check whether all the broad lines get
            # dropped. If so, don't use the broad-line model even if the chi2 is
            # formally lower because it generally means that the narrow lines
            # have not been optimized properly (see, e.g., hdelta in
            # sv1-bright-6541-39633076111803113).
            alldropped = np.all([getattr(broadfit, '{}_amp'.format(linename)).value == 0.0
                                 for linename in self.linetable['name'][IB]])
            if alldropped:
                log.info('All broad lines have been dropped, using narrow-line only model.')
            
            if alldropped or linechi2_broad > linechi2_init:
                bestfit = initfit
            else:
                bestfit = broadfit
        else:
            log.info('Too few pixels centered on candidate broad emission lines.')
            bestfit = initfit
            linechi2_init = 0.0
            linechi2_broad = 0.0

        # Finally, one more fitting loop with all the line-constraints relaxed.
        if False:
            for linename in data['coadd_linename']:
                # skip the physical doublets
                if not hasattr(bestfit, '{}_amp'.format(linename)): 
                    continue
                getattr(bestfit, '{}_amp'.format(linename)).tied = None
                getattr(bestfit, '{}_sigma'.format(linename)).tied = None
                getattr(bestfit, '{}_vshift'.format(linename)).tied = None
                    
            t0 = time.time()        
            fitter = FastLevMarLSQFitter(bestfit)
            finalfit = fitter(bestfit, emlinewave, emlineflux, weights=weights,
                              maxiter=maxiter, acc=accuracy)
            finalfit = _clean_linefit(finalfit, init_amplitudes)
            chi2 = self.chi2(finalfit, emlinewave, emlineflux, emlineivar)
            log.info('Final line-fitting took {:.2f} sec (niter={}) with chi2={:.3f}'.format(
                time.time()-t0, fitter.fit_info['nfev'], chi2))
        else:
            finalfit = bestfit

        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_output(self.EMLineModel.linetable)

        emlinemodel = finalfit(emlinewave)

        # get the full-spectrum reduced chi2
        rchi2, dof = self.chi2(finalfit, emlinewave, specflux, emlineivar, return_dof=True,
                               continuum_model=continuummodelflux + smoothcontinuummodelflux)

        result['RCHI2'] = rchi2
        #result['DOF'] = dof
        result['LINERCHI2_BROAD'] = linechi2_broad
        result['DELTA_LINERCHI2'] = linechi2_init - linechi2_broad

        # Now fill the output table.
        for pp in finalfit.param_names:
            pinfo = getattr(finalfit, pp)
            val = pinfo.value
            # If the parameter was not optimized, set its value to zero.
            if pinfo.fixed:
                val = 0.0
            # special case the tied doublets
            if pp == 'oii_doublet_ratio':
                result['OII_DOUBLET_RATIO'] = val
                result['OII_3726_AMP'] = val * finalfit.oii_3729_amp # * u.erg/(u.second*u.cm**2*u.Angstrom)
            elif pp == 'sii_doublet_ratio':
                result['SII_DOUBLET_RATIO'] = val
                result['SII_6731_AMP'] = val * finalfit.sii_6716_amp # * u.erg/(u.second*u.cm**2*u.Angstrom)
            elif pp == 'mgii_doublet_ratio':
                result['MGII_DOUBLET_RATIO'] = val
                result['MGII_2796_AMP'] = val * finalfit.mgii_2803_amp # * u.erg/(u.second*u.cm**2*u.Angstrom)
            else:
                result[pinfo.name.upper()] = val
                #if 'amp' in pinfo.name:
                #    result[pinfo.name.upper()] = val * u.erg/(u.second*u.cm**2*u.Angstrom)
                #elif 'vshift' in pinfo.name or 'sigma' in pinfo.name:
                #    result[pinfo.name.upper()] = val * u.kilometer / u.second
                #else:
                #    result[pinfo.name.upper()] = val

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
                result['FLUX_SYNTH_{}'.format(band.upper())] = data['synthphot']['nanomaggies'][iband] # * 'nanomaggies'
                #result['FLUX_SYNTH_IVAR_{}'.format(band.upper())] = data['synthphot']['nanomaggies_ivar'][iband]
            for iband, band in enumerate(self.synth_bands):
                result['FLUX_SYNTH_MODEL_{}'.format(band.upper())] = model_synthphot['nanomaggies'][iband] # * 'nanomaggies'

        specflux_nolines = specflux - emlinemodel

        ## measure DN(4000) without the emission lines
        #if False:
        #    dn4000_nolines, _ = self.get_dn4000(emlinewave, specflux_nolines, redshift=redshift)
        #    result['DN4000_NOLINES'] = dn4000_nolines

        # get continuum fluxes, EWs, and upper limits
        narrow_sigmas, broad_sigmas, uv_sigmas = [], [], []
        narrow_redshifts, broad_redshifts, uv_redshifts = [], [], []
        for oneline in self.EMLineModel.linetable[self.EMLineModel.inrange]:

            linename = oneline['name'].upper()
            #linez = redshift + result['{}_VSHIFT'.format(linename)][0].value / C_LIGHT
            linez = redshift + result['{}_VSHIFT'.format(linename)][0] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)

            linesigma = result['{}_SIGMA'.format(linename)][0] # [km/s]
            #linesigma = result['{}_SIGMA'.format(linename)][0].value # [km/s]

            # if the line was dropped, use a default sigma value
            if linesigma == 0:
                if oneline['isbroad']:
                    if oneline['isbalmer']:
                        linesigma = self.EMLineModel.limitsigma_narrow
                    else:
                        linesigma = self.EMLineModel.limitsigma_broad
                else:
                    linesigma = self.EMLineModel.limitsigma_broad

            linesigma_ang = linesigma * linezwave / C_LIGHT    # [observed-frame Angstrom]
            #log10sigma = linesigma / C_LIGHT / np.log(10)     # line-width [log-10 Angstrom]            

            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next line.
            lineindx = np.where((emlinewave >= (linezwave - 2.*linesigma_ang)) *
                                (emlinewave <= (linezwave + 2.*linesigma_ang)))[0]
            if np.sum(oemlineivar[lineindx] == 0) / len(lineindx) > 0.3:
                result['{}_AMP'.format(linename)] = 0.0
                result['{}_VSHIFT'.format(linename)] = 0.0
                result['{}_SIGMA'.format(linename)] = 0.0
            else:
                # number of pixels, chi2, and boxcar integration
                lineindx = np.where((emlinewave >= (linezwave - 2.*linesigma_ang)) *
                                    (emlinewave <= (linezwave + 2.*linesigma_ang)) *
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
                        good = emlineivar[lineindx] > 0
                        lineindx = lineindx[good]
    
                npix = len(lineindx)
                result['{}_NPIX'.format(linename)] = npix
    
                if npix > 3: # magic number: required at least XX unmasked pixels centered on the line
    
                    if np.any(emlineivar[lineindx] == 0):
                        errmsg = 'Ivar should never be zero within an emission line!'
                        log.critical(errmsg)
                        raise ValueError(errmsg)
                        
                    # boxcar integration of the flux
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
                        result['{}_FLUX'.format(linename)] = result['{}_AMP'.format(linename)][0] * linenorm
            
                        #result['{}_FLUX_IVAR'.format(linename)] = result['{}_AMP_IVAR'.format(linename)] / linenorm**2
                        #weight = np.exp(-0.5 * np.log10(emlinewave/linezwave)**2 / log10sigma**2)
                        #weight = (weight / np.max(weight)) > 1e-3
                        #result['{}_FLUX_IVAR'.format(linename)] = 1 / np.sum(1 / emlineivar[weight])
                        result['{}_FLUX_IVAR'.format(linename)] = boxflux_ivar # * u.second**2*u.cm**4/u.erg**2
    
                        dof = npix - 3 # ??? [redshift, sigma, and amplitude]
                        chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-emlinemodel[lineindx])**2) / dof
    
                        result['{}_CHI2'.format(linename)] = chi2
    
                        # keep track of sigma and z but only using XX-sigma lines
                        linesnr = result['{}_AMP'.format(linename)][0] * np.sqrt(result['{}_AMP_IVAR'.format(linename)][0])
                        #print(linename, result['{}_AMP'.format(linename)][0], amp_sigma, linesnr)
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
                                  (emlineivar > 0))[0]
                                  #(emlinemodel == 0))[0]
                indxhi = np.where((emlinewave < (linezwave + 10*linesigma * linezwave / C_LIGHT)) *
                                  (emlinewave > (linezwave + 3.*linesigma * linezwave / C_LIGHT)) *
                                  (emlineivar > 0))[0]
                                  #(emlinemodel == 0))[0]
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
    
                    #if np.abs(civar) > 1e10:
                    #    import matplotlib.pyplot as plt
                    #    plt.clf()
                    #    plt.plot(emlinewave, emlineivar)
                    #    plt.savefig('desi-users/ioannis/tmp/junk.png')
                    #    pdb.set_trace()
    
                    result['{}_CONT'.format(linename)] = cmed # * u.erg/(u.second*u.cm**2*u.Angstrom)
                    result['{}_CONT_IVAR'.format(linename)] = civar # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
    
                if result['{}_CONT'.format(linename)] != 0.0 and result['{}_CONT_IVAR'.format(linename)] != 0.0:
                    factor = (1 + redshift) / result['{}_CONT'.format(linename)] # --> rest frame
                    ew = result['{}_FLUX'.format(linename)] * factor # rest frame [A]
                    ewivar = result['{}_FLUX_IVAR'.format(linename)] / factor**2
    
                    # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                    fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar) # * u.erg/(u.second*u.cm**2)
                    ewlimit = fluxlimit * factor
    
                    result['{}_EW'.format(linename)] = ew
                    result['{}_EW_IVAR'.format(linename)] = ewivar
                    result['{}_FLUX_LIMIT'.format(linename)] = fluxlimit 
                    result['{}_EW_LIMIT'.format(linename)] = ewlimit

            if 'debug' in log.name:
                for col in ('VSHIFT', 'SIGMA', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)][0]))
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'BOXFLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT', 'EW_LIMIT'):
                    log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)][0]))
                log.debug(' ')
    
            ## simple QA
            #if 'alpha' in linename and False:
            #    sigma_cont = 150.0
            #    import matplotlib.pyplot as plt
            #    _indx = np.where((emlinewave > (linezwave - 15*sigma_cont * linezwave / C_LIGHT)) *
            #                    (emlinewave < (linezwave + 15*sigma_cont * linezwave / C_LIGHT)))[0]
            #    plt.plot(emlinewave[_indx], emlineflux[_indx])
            #    plt.plot(emlinewave[_indx], specflux_nolines[_indx])
            #    plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
            #    plt.axhline(y=cmed, color='k')
            #    plt.axhline(y=cmed+csig/np.sqrt(len(indx)), color='k', ls='--')
            #    plt.axhline(y=cmed-csig/np.sqrt(len(indx)), color='k', ls='--')
            #    plt.savefig('junk.png')

        # clean up the doublets; 
        if result['OIII_5007_AMP'] == 0 and result['OIII_5007_NPIX'] > 0:
            result['OIII_4959_AMP'] = 0.0
            result['OIII_4959_FLUX'] = 0.0
            result['OIII_4959_EW'] = 0.0

        if result['NII_6584_AMP'] == 0 and result['NII_6584_NPIX'] > 0:
            result['NII_6548_AMP'] = 0.0
            result['NII_6548_FLUX'] = 0.0
            result['NII_6548_EW'] = 0.0

        #if result['OII_3729_AMP'] == 0 and result['OII_3729_NPIX'] > 0:
        #    result['OII_3726_AMP'] = 0.0
        #    result['OII_3726_FLUX'] = 0.0
        #    result['OII_3726_EW'] = 0.0

        if result['OII_7320_AMP'] == 0 and result['OII_7320_NPIX'] > 0:
            result['OII_7330_AMP'] = 0.0
            result['OII_7330_FLUX'] = 0.0
            result['OII_7330_EW'] = 0.0

        if 'debug' in log.name:
            for col in ('MGII_DOUBLET_RATIO', 'OII_DOUBLET_RATIO', 'SII_DOUBLET_RATIO'):
                log.debug('{}: {:.4f}'.format(col, result[col][0]))
            log.debug(' ')

        # get the average emission-line redshifts and velocity widths
        if len(narrow_redshifts) > 0:
            result['NARROW_Z'] = np.mean(narrow_redshifts)
            result['NARROW_SIGMA'] = np.mean(narrow_sigmas) # * u.kilometer / u.second
            #result['NARROW_Z_ERR'] = np.std(narrow_redshifts)
            #result['NARROW_SIGMA_ERR'] = np.std(narrow_sigmas)
        else:
            result['NARROW_Z'] = redshift
            
        if len(broad_redshifts) > 0:
            result['BROAD_Z'] = np.mean(broad_redshifts)
            result['BROAD_SIGMA'] = np.mean(broad_sigmas) # * u.kilometer / u.second
            #result['BROAD_Z_ERR'] = np.std(broad_redshifts)
            #result['BROAD_SIGMA_ERR'] = np.std(broad_sigmas)
        else:
            result['BROAD_Z'] = redshift
            
        if len(uv_redshifts) > 0:
            result['UV_Z'] = np.mean(uv_redshifts)
            result['UV_SIGMA'] = np.mean(uv_sigmas) # * u.kilometer / u.second
            #result['UV_Z_ERR'] = np.std(uv_redshifts)
            #result['UV_SIGMA_ERR'] = np.std(uv_sigmas)
        else:
            result['UV_Z'] = redshift

        # fragile
        if 'debug' in log.name:
            for line in ('NARROW', 'BROAD', 'UV'):
                log.debug('{}_Z: {:.12f}'.format(line, result['{}_Z'.format(line)][0]))
                log.debug('{}_SIGMA: {:.3f}'.format(line, result['{}_SIGMA'.format(line)][0]))

        # As a consistency check, make sure that the emission-line spectrum
        # rebuilt from the final table is not (very) different from the one
        # based on the best-fitting model parameters.
        emmodel = np.hstack(self.emlinemodel_bestfit(data['wave'], data['res'], result, redshift=redshift))
        #try:
        #    assert(np.all(np.isclose(emmodel, emlinemodel, rtol=1e-4)))
        #except:
        #    pdb.set_trace()
            
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
        
        modelspectra.add_column(Column(name='CONTINUUM', dtype='f4', data=np.interp(modelwave, emlinewave, continuummodelflux).reshape(1, npix)))
        modelspectra.add_column(Column(name='SMOOTHCONTINUUM', dtype='f4', data=np.interp(modelwave, emlinewave, smoothcontinuummodelflux).reshape(1, npix)))
        modelspectra.add_column(Column(name='EMLINEMODEL', dtype='f4', data=np.interp(modelwave, emlinewave, emmodel).reshape(1, npix)))
        
        return result, modelspectra
    
    def qa_fastspec(self, data, fastspec, metadata, coadd_type='healpix',
                    wavelims=(3600, 9800), outprefix=None, outdir=None):
        """QA plot the emission-line spectrum and best-fitting model.

        """
        from astropy.table import Table, Column
        from scipy.ndimage import median_filter
        import matplotlib.pyplot as plt
        from matplotlib import colors
        import matplotlib.ticker as ticker
        import seaborn as sns

        from fastspecfit.util import ivar2var

        sns.set(context='talk', style='ticks', font_scale=1.5)#, rc=rc)

        col1 = [colors.to_hex(col) for col in ['dodgerblue', 'darkseagreen', 'orangered']]
        col2 = [colors.to_hex(col) for col in ['darkblue', 'darkgreen', 'darkred']]
        #col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
        #col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]
        col3 = [colors.to_hex(col) for col in ['blue', 'green', 'red']]

        if outdir is None:
            outdir = '.'
        if outprefix is None:
            outprefix = 'fastspec'

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
            title = 'Tile: {}, TargetID/Fiber: {}/{}'.format(
                    metadata['TILEID'], metadata['TARGETID'], metadata['FIBER'])
            pngfile = os.path.join(outdir, '{}-{}-{}.png'.format(
                    outprefix, metadata['TILEID'], metadata['TARGETID']))
        else:
            pass

        redshift = fastspec['CONTINUUM_Z']
        npixpercamera = [len(gw) for gw in data['wave']] # all pixels
        npixpercam = np.hstack([0, npixpercamera])
        
        stackwave = np.hstack(data['wave'])

        # rebuild the best-fitting spectroscopic and photometric models
        continuum, _ = self.SSP2data(self.sspflux, self.sspwave, redshift=redshift, 
                                     specwave=data['wave'], specres=data['res'],
                                     cameras=data['cameras'],
                                     AV=fastspec['CONTINUUM_AV'],
                                     vdisp=fastspec['CONTINUUM_VDISP'],
                                     coeff=fastspec['CONTINUUM_COEFF'],
                                     synthphot=False)

        residuals = [data['flux'][icam] - continuum[icam] for icam in np.arange(len(data['cameras']))]
        if np.all(fastspec['CONTINUUM_COEFF'] == 0):
            _smooth_continuum = np.zeros_like(stackwave)
        else:
            _smooth_continuum, _ = self.smooth_continuum(np.hstack(data['wave']), np.hstack(residuals),
                                                         np.hstack(data['ivar']), redshift=redshift,
                                                         linemask=np.hstack(data['linemask']))
        smooth_continuum = []
        for icam in np.arange(len(data['cameras'])): # iterate over cameras
            ipix = np.sum(npixpercam[:icam+1])
            jpix = np.sum(npixpercam[:icam+2])
            smooth_continuum.append(_smooth_continuum[ipix:jpix])

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
            _emlinemodel_oneline1 = self.emlinemodel_bestfit(data['wave'], data['res'], T)
            if np.sum(np.hstack(_emlinemodel_oneline1)) > 0:
                _emlinemodel_oneline.append(_emlinemodel_oneline1)

        # QA choices
        inches_wide = 20
        inches_fullspec = 6
        inches_perline = inches_fullspec / 2.0
        nlinepanels = 5

        nline = len(set(self.linetable['plotgroup']))
        nlinerows = int(np.ceil(nline / nlinepanels))
        nrows = 2 + nlinerows

        height_ratios = np.hstack([1, 1, [0.5]*nlinerows])

        fig = plt.figure(figsize=(inches_wide, 2*inches_fullspec + inches_perline*nlinerows))
        gs = fig.add_gridspec(nrows, nlinepanels, height_ratios=height_ratios)

        # full spectrum + best-fitting continuum model
        bigax1 = fig.add_subplot(gs[0, :])

        leg = {
            'zredrock': '$z_{{\\rm redrock}}$={:.6f}'.format(redshift),
            'dv_balmer': '$\\Delta v_{{\\rm narrow}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['NARROW_Z']-redshift)),
            'dv_forbid': '$\\Delta v_{{\\rm broad}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['BROAD_Z']-redshift)),
            'dv_broad': '$\\Delta v_{{\\rm UV}}$={:.2f} km/s'.format(C_LIGHT*(fastspec['UV_Z']-redshift)),
            'sigma_balmer': '$\\sigma_{{\\rm narrow}}$={:.1f} km/s'.format(fastspec['NARROW_SIGMA']),
            'sigma_forbid': '$\\sigma_{{\\rm broad}}$={:.1f} km/s'.format(fastspec['BROAD_SIGMA']),
            'sigma_broad': '$\\sigma_{{\\rm UV}}$={:.1f} km/s'.format(fastspec['UV_SIGMA']),
            #'targetid': '{} {}'.format(metadata['TARGETID'], metadata['FIBER']),
            #'targetid': 'targetid={} fiber={}'.format(metadata['TARGETID'], metadata['FIBER']),
            'chi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastspec['CONTINUUM_RCHI2']),
            'rchi2': '$\\chi^{{2}}_{{\\nu}}$={:.3f}'.format(fastspec['RCHI2']),
            'deltarchi2': '$\\Delta\\chi^{{2}}_{{\\nu,\\rm broad,narrow}}$={:.3f}'.format(fastspec['DELTA_LINERCHI2']),
            #'zfastfastspec': '$z_{{\\rm fastspecfit}}$={:.6f}'.format(fastspec['CONTINUUM_Z']),
            #'z': '$z$={:.6f}'.format(fastspec['CONTINUUM_Z']),
            'age': '<Age>={:.3f} Gyr'.format(fastspec['CONTINUUM_AGE']),
            }

        if fastspec['CONTINUUM_VDISP_IVAR'] == 0:
            leg.update({'vdisp': '$\\sigma_{{\\rm star}}$={:.1f} km/s'.format(fastspec['CONTINUUM_VDISP'])})
        else:
            leg.update({'vdisp': '$\\sigma_{{\\rm star}}$={:.1f}+/-{:.1f} km/s'.format(
                fastspec['CONTINUUM_VDISP'], 1/np.sqrt(fastspec['CONTINUUM_VDISP_IVAR']))})
            
        if fastspec['CONTINUUM_AV_IVAR'] == 0:
            leg.update({'AV': '$A(V)$={:.3f} mag'.format(fastspec['CONTINUUM_AV'])})
        else:
            leg.update({'AV': '$A(V)$={:.3f}+/-{:.3f} mag'.format(
                fastspec['CONTINUUM_AV'], 1/np.sqrt(fastspec['CONTINUUM_AV_IVAR']))})

        ymin, ymax = 1e6, -1e6

        legxpos, legypos, legfntsz = 0.98, 0.94, 20
        bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.25)

        #ymin = np.zeros(len(data['cameras']))
        #ymax = np.zeros(len(data['cameras']))
        for ii in np.arange(len(data['cameras'])): # iterate over cameras
            sigma, good = ivar2var(data['ivar'][ii], sigma=True, allmasked_ok=True)

            bigax1.fill_between(data['wave'][ii], data['flux'][ii]-sigma,
                                data['flux'][ii]+sigma, color=col1[ii])
            #bigax1.plot(data['wave'][ii], continuum_nodust[ii], alpha=0.5, color='k')
            #bigax1.plot(data['wave'][ii], smooth_continuum[ii], color='gray')#col3[ii])#, alpha=0.3, lw=2)#, color='k')
            #bigax1.plot(data['wave'][ii], continuum[ii], color=col2[ii])
            bigax1.plot(data['wave'][ii], continuum[ii]+smooth_continuum[ii], color=col2[ii])
            
            # get the robust range
            filtflux = median_filter(data['flux'][ii], 51, mode='nearest')
            #filtflux = median_filter(data['flux'][ii] - _emlinemodel[ii], 51, mode='nearest')
            #perc = np.percentile(filtflux[data['ivar'][ii] > 0], [5, 95])
            #sigflux = np.std(data['flux'][ii][data['ivar'][ii] > 0])
            I = data['ivar'][ii] > 0
            if np.sum(I) > 0:
                sigflux = np.diff(np.percentile(data['flux'][ii][I], [25, 75]))[0] / 1.349 # robust
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
                ymax = sigflux * 8 # 5
            if np.max(filtflux) > ymax:
                ymax = np.max(filtflux) * 1.4
            #print(ymin, ymax)

        #ymin = np.min(ymin)
        #ymax = np.max(ymax)

        bigax1.plot(stackwave, _smooth_continuum, color='gray')#col3[ii])#, alpha=0.3, lw=2)#, color='k')
        bigax1.plot(stackwave, np.hstack(continuum), color='k', alpha=0.1)#col3[ii])#, alpha=0.3, lw=2)#, color='k')

        bigax1.text(0.03, 0.9, 'Observed Spectrum + Continuum Model',
                    ha='left', va='center', transform=bigax1.transAxes, fontsize=30)
        if not self.nolegend:
            txt = '\n'.join((
                r'{}'.format(leg['zredrock']),
                r'{} {}'.format(leg['chi2'], leg['age']),
                r'{}'.format(leg['AV']),
                r'{}'.format(leg['vdisp']),
                ))
            bigax1.text(legxpos, legypos, txt, ha='right', va='top',
                        transform=bigax1.transAxes, fontsize=legfntsz,
                        bbox=bbox)
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

            emlinesigma, good = ivar2var(data['ivar'][ii], sigma=True, allmasked_ok=True)
            emlinewave = emlinewave[good]
            emlineflux = emlineflux[good]
            emlinesigma = emlinesigma[good]
            emlinemodel = emlinemodel[good]

            bigax2.fill_between(emlinewave, emlineflux-emlinesigma,
                                emlineflux+emlinesigma, color=col1[ii], alpha=0.7)
            bigax2.plot(emlinewave, emlinemodel, color=col2[ii], lw=3)

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

        if not self.nolegend:
            txt = '\n'.join((
                r'{} {}'.format(leg['rchi2'], leg['deltarchi2']),
                r'{} {}'.format(leg['dv_balmer'], leg['sigma_balmer']),
                r'{} {}'.format(leg['dv_forbid'], leg['sigma_forbid']),
                r'{} {}'.format(leg['dv_broad'], leg['sigma_broad']),
                ))
            bigax2.text(legxpos, legypos, txt, ha='right', va='top',
                        transform=bigax2.transAxes, fontsize=legfntsz,
                        bbox=bbox)
            bigax2.text(0.03, 0.9, 'Residual Spectrum + Emission-Line Model',
                        ha='left', va='center', transform=bigax2.transAxes,
                        fontsize=30)
                
        bigax2.set_xlim(wavelims)
        bigax2.set_ylim(ymin, ymax)
        
        #bigax2.set_xlabel(r'Observed-frame Wavelength ($\AA$)') 
        #bigax2.set_ylabel(r'Flux ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)') 
        
        # zoom in on individual emission lines - use linetable!
        plotsig_default = 200.0 # [km/s]
        plotsig_default_balmer = 500.0 # [km/s]
        plotsig_default_broad = 2000.0 # [km/s]

        meanwaves, deltawaves, sigmas, linenames = [], [], [], []
        for plotgroup in set(self.linetable['plotgroup']):
            I = np.where(plotgroup == self.linetable['plotgroup'])[0]
            linenames.append(self.linetable['nicename'][I[0]].replace('-', ' '))
            meanwaves.append(np.mean(self.linetable['restwave'][I]))
            deltawaves.append((np.max(self.linetable['restwave'][I]) - np.min(self.linetable['restwave'][I])) / 2)

            sigmas1 = np.array([fastspec['{}_SIGMA'.format(line.upper())] for line in self.linetable[I]['name']])
            sigmas1 = sigmas1[sigmas1 > 0]
            if len(sigmas1) > 0:
                plotsig = 1.5*np.mean(sigmas1)
            else:
                if np.any(self.linetable['isbroad'][I]):
                    if np.any(self.linetable['isbalmer'][I]):
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
        
        ax, irow, icol = [], 2, 0
        for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):
            icol = iax % nlinepanels
            if iax > 0 and iax % nlinepanels == 0:
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
                emlineflux = data['flux'][ii] - continuum[ii] - smooth_continuum[ii]
                emlinemodel = _emlinemodel[ii]

                emlinesigma, good = ivar2var(data['ivar'][ii], sigma=True, allmasked_ok=True)
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
                    xx.fill_between(emlinewave[indx], emlineflux[indx]-emlinesigma[indx],
                                    emlineflux[indx]+emlinesigma[indx], color=col1[ii],
                                    alpha=0.5)
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
                xx.xaxis.set_major_locator(ticker.MaxNLocator(2))

        # common axis labels
        tp, bt, lf, rt = 0.95, 0.08, 0.10, 0.94
        
        fig.text(lf-0.07, (tp-bt)/4+bt+(tp-bt)/1.9,
                 r'Flux Density ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
                 ha='center', va='center', rotation='vertical', fontsize=30)
        fig.text(lf-0.07, (tp-bt)/4+bt,
                 r'Flux Density ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
                 ha='center', va='center', rotation='vertical', fontsize=30)
        fig.text((rt-lf)/2+lf, bt-0.04, r'Observed-frame Wavelength ($\AA$)',
                 ha='center', va='center', fontsize=30)
            
        plt.subplots_adjust(wspace=0.27, top=tp, bottom=bt, left=lf, right=rt, hspace=0.22)

        log.info('Writing {}'.format(pngfile))
        fig.savefig(pngfile)
        plt.close()

