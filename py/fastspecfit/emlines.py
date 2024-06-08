"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import pdb # for debugging

import os, time
import numpy as np
from numba import jit
from math import erf, erfc
from astropy.table import Table

from fastspecfit.io import FLUXNORM
from fastspecfit.continuum import Filters
from fastspecfit.util import C_LIGHT

from fastspecfit.emline_fit import (
    EMLine_Objective,
    EMLine_find_peak_amplitudes,
    EMLine_build_model,
    EMLine_ParamsMapping,
)


def read_emlines(emlinesfile=None):
    """Read the set of emission lines of interest.

    """
    if emlinesfile is None:
        from importlib import resources
        emlinesfile = resources.files('fastspecfit').joinpath('data/emlines.ecsv')

    try:
        linetable = Table.read(emlinesfile, format='ascii.ecsv', guess=False)
    except: 
        from desiutil.log import get_logger
        log = get_logger()
        errmsg = f'Problem reading emission lines parameter file {emlinesfile}.'
        log.critical(errmsg)
        raise ValueError(errmsg)
    
    return linetable    

class EMFitTools(Filters):
    def __init__(self, fphoto=None, emlinesfile=None, uniqueid=None,
                 minsnr_balmer_broad=3.):
        """Class to model a galaxy stellar continuum.

        Parameters
        ----------
        templates : :class:`str`, optional
            Full path to the templates used for continuum-fitting.
        mintemplatewave : :class:`float`, optional, defaults to None
            Minimum template wavelength to read into memory. If ``None``, the minimum
            available wavelength is used (around 100 Angstrom).
        maxtemplatewave : :class:`float`, optional, defaults to 6e4
            Maximum template wavelength to read into memory. 
        chi2_default : :class:`float`, optional, defaults to 0.0.
            Default chi2 value for a emission line not fitted.
        maxiter : :class:`int`, optional, defaults to 5000.
            Maximum number of iterations.
        accuracy : :class:`float`, optional, defaults to 0.01.
            Fitting accuracy.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.

        Notes
        -----
        Need to document all the attributes.
        
        Plans for improvement:
          - Update the continuum redshift using cross-correlation.
          - Don't draw reddening from a flat distribution (try gamma or a custom
            distribution of the form x**2*np.exp(-2*x/scale).

        """
        super(EMFitTools, self).__init__(fphoto=fphoto)

        self.uniqueid = uniqueid

        self.linetable = read_emlines(emlinesfile=emlinesfile)

        #self.emwave_pixkms = 5.                                   # pixel size for internal wavelength array [km/s]
        #self.dlog10wave = self.emwave_pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]

        # default line-sigma for computing upper limits
        self.limitsigma_narrow = 75.0
        self.limitsigma_broad = 1200.0 
        self.wavepad = 2.5 # Angstrom

        # Establish the names of the parameters and doublets here, at
        # initialization, because we use them when instantiating the best-fit
        # model not just when fitting.
        doublet_names = ['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio']
        doublet_pairs = ['mgii_2803_amp', 'oii_3729_amp', 'sii_6716_amp']

        param_names = []
        for param in ['amp', 'vshift', 'sigma']:
            for linename in self.linetable['name'].data:
                param_name = linename+'_'+param
                # Use doublet-ratio parameters for close or physical
                # doublets. Note that any changes here need to be propagated to
                # the XX method, which needs to "know" about these doublets.
                if param_name == 'mgii_2796_amp':
                    param_name = 'mgii_doublet_ratio' # MgII 2796/2803
                if param_name == 'oii_3726_amp':
                    param_name = 'oii_doublet_ratio'  # [OII] 3726/3729
                if param_name == 'sii_6731_amp':
                    param_name = 'sii_doublet_ratio'  # [SII] 6731/6716
                param_names.append(param_name)
        self.param_names = np.hstack(param_names)
        self.amp_param_bool = np.array(['_amp' in pp for pp in self.param_names])
        self.amp_balmer_bool = np.array(['_amp' in pp and 'hei_' not in pp and 'heii_' not in pp for pp in self.param_names]) # no helium lines
        self.sigma_param_bool = np.array(['_sigma' in pp for pp in self.param_names])
        self.vshift_param_bool = np.array(['_vshift' in pp for pp in self.param_names])

        self.doubletindx = np.hstack([np.where(self.param_names == doublet)[0] for doublet in doublet_names])
        self.doubletpair = np.hstack([np.where(self.param_names == pair)[0] for pair in doublet_pairs])

        self.minsigma_balmer_broad = 250. # minimum broad-line sigma [km/s]
        self.minsnr_balmer_broad = minsnr_balmer_broad # minimum broad-line S/N

    def summarize_linemodel(self, linemodel):
        """Simple function to summarize an input linemodel."""
        def _print(linenames):
            for linename in linenames:
                for param in ['amp', 'sigma', 'vshift']:
                    I = np.where(self.param_names == linename+'_'+param)[0]
                    if len(I) == 1:
                        I = I[0]
                        if linemodel['tiedtoparam'][I] == -1:
                            if linemodel['fixed'][I]:
                                print('{:25s} is NOT FITTED'.format(linename+'_'+param))
                            else:
                                print('{:25s} untied'.format(linename+'_'+param))
                        else:
                            if linemodel['fixed'][I]:
                                print('{:25s} tied to {:25s} with factor {:.4f} and FIXED'.format(
                                    linename+'_'+param, self.param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))
                            else:
                                print('{:25s} tied to {:25s} with factor {:.4f}'.format(
                                    linename+'_'+param, self.param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))
                                
        linenames = self.fit_linetable['name'].data

        print('---------------------')
        print('UV/QSO (broad) lines:')
        print('---------------------')
        _print(linenames[(self.fit_linetable['isbroad'] == True) * (self.fit_linetable['isbalmer'] == False)])
        print()
        print('--------------------------')
        print('Broad Balmer+helium lines:')
        print('--------------------------')
        _print(linenames[(self.fit_linetable['isbroad'] == True) * (self.fit_linetable['isbalmer'] == True)])
        print()
        print('---------------------------')
        print('Narrow Balmer+helium lines:')
        print('---------------------------')
        _print(linenames[(self.fit_linetable['isbroad'] == False) * (self.fit_linetable['isbalmer'] == True)])
        print()
        print('----------------')
        print('Forbidden lines:')
        print('----------------')
        _print(linenames[(self.fit_linetable['isbroad'] == False) * (self.fit_linetable['isbalmer'] == False)])

    def build_linemodels(self, redshift, wavelims=[3000, 10000], verbose=False, strict_broadmodel=True):
        """Build all the multi-parameter emission-line models we will use.
    
        """
        def _fix_parameters(linemodel, verbose=False):
            """Set the "fixed" attribute for all the parameters in a given linemodel."""
            # First loop through all tied parameters and set fixed to the
            # parameter it's tied to.
            I = np.where(linemodel['tiedtoparam'] != -1)[0] # should always have len(I)>0
            alltied = linemodel[I]['tiedtoparam']
            utied = np.unique(alltied)
            for tied in utied:
                J = tied == alltied
                if verbose:
                    print('Tying {} to {}'.format(' '.join(linemodel[I][J]['param_name']), linemodel[tied]['param_name']))
                linemodel[I][J]['fixed'] = linemodel[tied]['fixed']
            if verbose:
                print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

            # Next, fix out-of-range lines but not those that are in the 'utied'
            # array---those out-of-range lines need to be in the optimization
            # list because the in-range lines depend on them.
            outofrange = fit_linetable['inrange'] == False
            if np.sum(outofrange) > 0: # should always be true
                for linename in fit_linetable['name'][outofrange]:
                    for param in ['amp', 'vshift', 'sigma']:
                        param_name = linename+'_'+param
                        I = np.where(linemodel['param_name'] == param_name)[0]
                        if len(I) > 0:
                            if I in utied:
                                if verbose:
                                    print('Not fixing out-of-range parameter {}'.format(param_name))
                            else:
                                linemodel['fixed'][I] |= True
                if verbose:
                    print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                    print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                    #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

                # Finally loop through each 'utied' line and if all the lines
                # tied to it are fixed, then fix that line, too.
                for tied in utied:
                    if linemodel['param_name'][tied] and np.all(linemodel[linemodel['tiedtoparam'] == tied]['fixed']):
                        if outofrange[linemodel['linename'][tied] == fit_linetable['name']]:
                            if verbose:
                                print('Fixing {} because line is out of range and all tied lines are fixed: {}'.format(
                                    linemodel['param_name'][tied], ' '.join(linemodel[linemodel['tiedtoparam'] == tied]['param_name'])))
                            linemodel[tied]['fixed'] = True

                # Also handle the doublets.
                I = np.where(linemodel['doubletpair'] != -1)[0]
                if len(I) > 0:
                    for doublet in linemodel[I]['doubletpair']:
                        J = doublet == linemodel['doubletpair']
                        if linemodel[doublet]['fixed']:
                            linemodel['fixed'][J] = True

                if verbose:
                    print('Number of fixed parameters = {}'.format(np.sum(linemodel['fixed'])))
                    print('Number of free parameters = {}'.format(np.sum(np.logical_and(linemodel['fixed'] == False, linemodel['tiedtoparam'] == -1))))
                    #print('Number of fixed or tied parameters = {}'.format(np.sum(np.logical_or(linemodel['fixed'], linemodel['tiedtoparam'] != -1))))

        initvshift = 1.0
        vmaxshift_narrow = 500.0
        vmaxshift_broad = 2500.0 # 3000.0
        vmaxshift_balmer_broad = 2500.
    
        minsigma_narrow = 1.0
        maxsigma_narrow = 750.0 # 500.0

        minsigma_broad = 1.0 # 100.
        maxsigma_broad = 1e4

        minsigma_balmer_broad = 1.0 # 100.0 # minsigma_narrow
        maxsigma_balmer_broad = maxsigma_broad
    
        # Be very careful about changing the default broad line-sigma. Smaller
        # values like 1500 km/s (which is arguably more sensible) can lead to
        # low-amplitude broad lines in a bunch of normal star-forming galaxy
        # spectra. (They act to "suck up" local continuum variations.) Also
        # recall that if it's well-measured, we use the initial line-sigma in
        # estimate_linesigma, which is a better initial guess.
        initsigma_narrow = 75.0 # 260.0 # 75.0
        initsigma_broad = 3000.0  
    
        initamp = 0.0
        #minamp = 0.0
        minamp = -1e2
        maxamp = +1e5
        minamp_balmer_broad = minamp # 0.0
        maxamp_balmer_broad = maxamp
    
        # Specialized parameters on the MgII, [OII], and [SII] doublet ratios. See
        # https://github.com/desihub/fastspecfit/issues/39. Be sure to set
        # self.doublet_names, below, and also note that any change in the order of
        # these lines has to be handled in _emline_spectrum!
        init_mgii_doublet = 0.5 # MgII 2796/2803
        init_oii_doublet = 0.74 # [OII] 3726/3729
        init_sii_doublet = 0.74 # [SII] 6731/6716

        bounds_mgii_doublet = [0.0, 10.0] 
        bounds_oii_doublet = [0.0, 2.0] # [0.5, 1.5] # [0.66, 1.4]
        bounds_sii_doublet = [0.0, 2.0] # [0.5, 1.5] # [0.67, 1.2]
    
        # Create a new line-fitting table which contains the redshift-dependent
        # quantities for this object.
        fit_linetable = Table()
        fit_linetable['name'] = self.linetable['name']
        fit_linetable['isbalmer'] = self.linetable['isbalmer']
        fit_linetable['ishelium'] = self.linetable['ishelium']
        fit_linetable['isbroad'] = self.linetable['isbroad']
        fit_linetable['restwave'] = self.linetable['restwave']
        fit_linetable['zwave'] = self.linetable['restwave'].data * (1 + redshift)
        fit_linetable['inrange'] = ((fit_linetable['zwave'] > (wavelims[0]+self.wavepad)) * 
                                    (fit_linetable['zwave'] < (wavelims[1]-self.wavepad)))
        self.fit_linetable = fit_linetable
        
        linenames = fit_linetable['name'].data
        param_names = self.param_names
        nparam = len(param_names)

        # Model 1 -- here, parameters are minimally tied together for the final
        # fit and only lines outside the wavelength range are fixed. Includes
        # broad lines.
        linemodel_broad = Table()
        linemodel_broad['param_name'] = param_names
        linemodel_broad['index'] = np.arange(nparam).astype(np.int32)
        linemodel_broad['linename'] = np.tile(linenames, 3) # 3 parameters per line
        linemodel_broad['isbalmer'] = np.zeros(nparam, bool)
        linemodel_broad['ishelium'] = np.zeros(nparam, bool)
        linemodel_broad['isbroad'] = np.zeros(nparam, bool)
        linemodel_broad['tiedfactor'] = np.zeros(nparam, 'f8')
        linemodel_broad['tiedtoparam'] = np.zeros(nparam, np.int16)-1
        linemodel_broad['doubletpair'] = np.zeros(nparam, np.int16)-1
        linemodel_broad['fixed'] = np.zeros(nparam, bool)
        linemodel_broad['bounds'] = np.zeros((nparam, 2), 'f8')
        linemodel_broad['initial'] = np.zeros(nparam, 'f8')
        linemodel_broad['value'] = np.zeros(nparam, 'f8')
        linemodel_broad['obsvalue'] = np.zeros(nparam, 'f8')
        linemodel_broad['civar'] = np.zeros(nparam, 'f8') # continuum inverse variance

        linemodel_broad['doubletpair'][self.doubletindx] = self.doubletpair

        # Build the relationship of "tied" parameters. In the 'tied' array, the
        # non-zero value is the multiplicative factor by which the parameter
        # represented in the 'tiedtoparam' index should be multiplied.
    
        # Physical doublets and lines in the same ionization species should have
        # their velocity shifts and line-widths always tied. In addition, set fixed
        # doublet-ratios here. Note that these constraints must be set on *all*
        # lines, not just those in range.
    
        for iline, linename in enumerate(linenames):
            linemodel_broad['isbalmer'][linemodel_broad['linename'] == linename] = fit_linetable[fit_linetable['name'] == linename]['isbalmer']
            linemodel_broad['ishelium'][linemodel_broad['linename'] == linename] = fit_linetable[fit_linetable['name'] == linename]['ishelium']
            linemodel_broad['isbroad'][linemodel_broad['linename'] == linename] = fit_linetable[fit_linetable['name'] == linename]['isbroad']
            
            # initial values and bounds - broad He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp_balmer_broad, maxamp_balmer_broad], 
                                                   [minsigma_balmer_broad, maxsigma_balmer_broad],
                                                   [-vmaxshift_balmer_broad, +vmaxshift_balmer_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - narrow He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - broad UV/QSO lines (non-Balmer)
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_broad, maxsigma_broad],
                                                   [-vmaxshift_broad, +vmaxshift_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - forbidden lines
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    linemodel_broad['initial'][param_names == linename+'_'+param] = default
                    linemodel_broad['bounds'][param_names == linename+'_'+param] = bounds

            # tie parameters

            # broad He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]
            #print('Releasing the narrow Balmer lines!')
            # narrow He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False and linename != 'halpha':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
            # other lines
            if linename == 'mgii_2796':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'mgii_2803_'+param)[0]
            if linename == 'nev_3346' or linename == 'nev_3426': # should [NeIII] 3869 be tied to [NeV]???
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'neiii_3869_'+param)[0]
            if linename == 'oii_3726':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_3729_'+param)[0]
            # Tentative! Tie auroral lines to [OIII] 4363 but maybe we shouldn't tie [OI] 6300 here...
            if linename == 'nii_5755' or linename == 'oi_6300' or linename == 'siii_6312':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_4363_'+param)[0]
            if linename == 'oiii_4959':
                """
                [O3] (4-->2): airwave: 4958.9097 vacwave: 4960.2937 emissivity: 1.172e-21
                [O3] (4-->3): airwave: 5006.8417 vacwave: 5008.2383 emissivity: 3.497e-21
                """
                linemodel_broad['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9839 # 2.8875
                linemodel_broad['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oiii_5007_amp')[0]
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
            if linename == 'nii_6548':
                """
                [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.02198e-21
                [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.94901e-21
                """
                linemodel_broad['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9421 # 2.936
                linemodel_broad['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'nii_6584_amp')[0]
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nii_6584_'+param)[0]
            if linename == 'sii_6731':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'sii_6716_'+param)[0]
            if linename == 'oii_7330':
                """
                [O2] (5-->2): airwave: 7318.9185 vacwave: 7320.9350 emissivity: 8.18137e-24
                [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 2.40519e-23
                [O2] (5-->3): airwave: 7329.6613 vacwave: 7331.6807 emissivity: 1.35614e-23
                [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.27488e-23
                """
                linemodel_broad['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 1.2251
                linemodel_broad['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oii_7320_amp')[0]
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_7320_'+param)[0]
            if linename == 'siii_9069':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'siii_9532_'+param)[0]
            # Tentative! Tie SiIII] 1892 to CIII] 1908 because they're so close in wavelength.
            if linename == 'siliii_1892':
                for param in ['sigma', 'vshift']:
                    linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'ciii_1908_'+param)[0]

            # Tie all the forbidden and narrow Balmer+helium lines *except
            # [OIII] 4959,5007* to [NII] 6584 when we have broad lines. The
            # [OIII] doublet frequently has an outflow component, so fit it
            # separately. See the discussion at
            # https://github.com/desihub/fastspecfit/issues/160
            if strict_broadmodel:
                if fit_linetable['isbroad'][iline] == False and linename != 'nii_6584' and linename != 'oiii_4959' and linename != 'oiii_5007':
                    for param in ['sigma', 'vshift']:
                        linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                        linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nii_6584_'+param)[0]
                        
                #if fit_linetable['isbroad'][iline] == False and linename != 'oiii_5007':
                #    for param in ['sigma', 'vshift']:
                #        linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                #        linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
                
                ## Tie all forbidden lines to [OIII] 5007; the narrow Balmer and
                ## helium lines are separately tied together.
                #if fit_linetable['isbroad'][iline] == False and fit_linetable['isbalmer'][iline] == False and linename != 'oiii_5007':
                #    for param in ['sigma']:
                #        linemodel_broad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                #        linemodel_broad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
                    
        # Finally set the initial values and bounds on the doublet ratio parameters.
        for param, bounds, default in zip(['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio'],
                                          [bounds_mgii_doublet, bounds_oii_doublet, bounds_sii_doublet],
                                          [init_mgii_doublet, init_oii_doublet, init_sii_doublet]):
            linemodel_broad['initial'][linemodel_broad['param_name'] == param] = default
            linemodel_broad['bounds'][linemodel_broad['param_name'] == param] = bounds
                    
        # Assign fixed=True to parameters which are outside the wavelength range
        # except those that are tied to other lines.
        _fix_parameters(linemodel_broad, verbose=False)

        assert(np.all(linemodel_broad['tiedtoparam'][linemodel_broad['tiedfactor'] != 0] != -1))
        # It's OK for the doublet ratios to be bounded at zero.
        #assert(len(linemodel_broad[np.sum(linemodel_broad['bounds'] == [0.0, 0.0], axis=1) > 0]) == 0)
    
        #_print_linemodel(linemodel_broad)
        #linemodel_broad[np.logical_and(linemodel_broad['fixed'] == False, linemodel_broad['tiedtoparam'] == -1)]

        # Model 2 - like linemodel, but broad lines have been fixed at zero.
        linemodel_nobroad = linemodel_broad.copy()
        linemodel_nobroad['fixed'] = False # reset

        for iline, linename in enumerate(linenames):
            if linename == 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    linemodel_nobroad['fixed'][param_names == linename+'_'+param] = True

            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

            if strict_broadmodel:
                # Tie the forbidden lines to [OIII] 5007.
                if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline] == False and linename != 'oiii_5007':
                    for param in ['sigma', 'vshift']:
                        linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                        linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
                        
                # Tie narrow Balmer and helium lines together.
                if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False:
                    if linename == 'halpha':
                        for param in ['sigma', 'vshift']:
                            linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 0.0
                            linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = -1
                    else:
                        for param in ['sigma', 'vshift']:
                            linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                            linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
                
        #linemodel_nobroad[np.logical_and(linemodel_nobroad['fixed'] == False, linemodel_nobroad['tiedtoparam'] == -1)]

        _fix_parameters(linemodel_nobroad, verbose=False)
        assert(np.all(linemodel_nobroad['tiedtoparam'][linemodel_nobroad['tiedfactor'] != 0] != -1))

        return linemodel_broad, linemodel_nobroad

    def initial_guesses_and_bounds(self, data, emlinewave, emlineflux, log):
        """For all lines in the wavelength range of the data, get a good initial guess
        on the amplitudes and line-widths. This step is critical for cases like,
        e.g., 39633354915582193 (tile 80613, petal 05), which has strong narrow
        lines.

        """
        initial_guesses, param_bounds = {}, {}
        #init_amplitudes, init_sigmas = {}, {}

        coadd_sigma = data['smoothsigma'] # robust estimate of the variance in the spectrum
        coadd_emlineflux = np.interp(data['coadd_wave'], emlinewave, emlineflux)
    
        for linename, linepix, contpix in zip(data['coadd_linename'], data['coadd_linepix'], data['coadd_contpix']):
            ## skip the physical doublets
            #if not hasattr(self.EMLineModel, '{}_amp'.format(linename)):
            #    continue

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
            if amp < 0:
                amp = np.abs(amp)

            noise = np.mean(coadd_sigma[linepix])
            if noise == 0.:
                civar = 0.
                errmsg = 'Noise estimate for line {} is zero!'.format(linename)
                log.warning(errmsg)
                #raise ValueError(errmsg)
            else:
                civar = 1. / noise**2

            # update the bounds on the line-amplitude
            #bounds = [-np.min(np.abs(coadd_emlineflux[linepix])), 3*np.max(coadd_emlineflux[linepix])]
            mx = 5*np.max(coadd_emlineflux[linepix])
            if mx < 0: # ???
                mx = 5*np.max(np.abs(coadd_emlineflux[linepix]))
            
            iline = self.linetable[self.linetable['name'] == linename]
            bounds = [-1.5*np.min(np.abs(coadd_emlineflux[linepix])), mx]

            # In extremely rare cases many of the pixels are zero, in which case
            # bounds[0] becomes zero, which is bad (e.g.,
            # iron/main/dark/27054/39627811564029314). Fix that here.
            if np.abs(bounds[0]) == 0.0:
                N = coadd_emlineflux[linepix] != 0
                if np.sum(N) > 0:
                    bounds[0] = -1.5*np.min(np.abs(coadd_emlineflux[linepix][N]))
                if np.abs(bounds[0]) == 0.0:
                    bounds[0] = -1e-3 # ??
            
            if (bounds[0] > bounds[1]) or (amp < bounds[0]) or (amp > bounds[1]):
                log.warning('Initial amplitude is outside its bound for line {}.'.format(linename))
                amp = np.diff(bounds)/2 + bounds[0]
                # Should never happen.
                if (bounds[0] > bounds[1]) or (amp < bounds[0]) or (amp > bounds[1]):
                    errmsg = 'Initial amplitude is outside its bound for line {}.'.format(linename)
                    self.log.critical(errmsg)
                    raise ValueError(errmsg)

            initial_guesses[linename+'_amp'] = amp
            param_bounds[linename+'_amp'] = bounds
            initial_guesses[linename+'_civar'] = civar

        # Now update the linewidth but here we need to loop over *all* lines
        # (not just those in range). E.g., if H-alpha is out of range we need to
        # set its initial value correctly since other lines are tied to it
        # (e.g., main-bright-32406-39628257196245904).
        for iline in self.linetable:
            linename = iline['name']
            if iline['isbroad']:
                if iline['isbalmer']: # broad Balmer lines
                    if data['linesigma_balmer'] > data['linesigma_narrow']:
                        initial_guesses[linename+'_sigma'] = data['linesigma_balmer']
                    else:
                        initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
                else: # broad UV/QSO lines
                    initial_guesses[linename+'_sigma'] = data['linesigma_uv']
            else:
                # prefer narrow over Balmer
                initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
    
        return initial_guesses, param_bounds

    @staticmethod
    def _linemodel_to_parameters(linemodel, fit_linetable):
        """Convert a linemodel model to a list of emission-line parameters."""

        linesplit = np.array_split(linemodel['index'], 3) # 3 parameters per line
        #linesplit = (np.arange(3) + 1) * len(linemodel) // 3 # 3 parameters per line
        lineamps = linemodel['value'][linesplit[0]].data
        linevshifts = linemodel['value'][linesplit[1]].data
        linesigmas = linemodel['value'][linesplit[2]].data

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        #doublet = np.where(linemodel['doubletpair'] != -1)[0]
        #lineamps[doublet] *= linemodel['value'][linemodel['doubletpair'][doublet]]
        parameters = np.hstack((lineamps, linevshifts, linesigmas))

        linewaves = fit_linetable['restwave'].data
        #lineinrange = fit_linetable['inrange'].data

        #Itied = np.where((linemodel['tiedtoparam'] != -1))[0]
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        Ifree = np.where((linemodel['tiedtoparam'] == -1) * (linemodel['fixed'] == False))[0]

        tiedtoparam = linemodel['tiedtoparam'][Itied].data
        tiedfactor = linemodel['tiedfactor'][Itied].data
        bounds = linemodel['bounds'][Ifree].data

        doubletindx = np.where(linemodel['doubletpair'] != -1)[0]
        doubletpair = linemodel['doubletpair'][doubletindx].data

        parameter_extras = (Ifree, Itied, tiedtoparam, tiedfactor, bounds,
                            doubletindx, doubletpair, linewaves)

        return parameters, parameter_extras

    @staticmethod
    def populate_linemodel(linemodel, initial_guesses, param_bounds, log):
        """Populate an input linemodel with initial guesses and parameter bounds, taking
        into account fixed parameters.

        """
        # Set initial values and bounds.
        for iparam, param in enumerate(linemodel['param_name']):
            if param in initial_guesses.keys():
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0 # always set fixed parameter to zero
                else:
                    # Make sure the initial guess for the narrow Balmer+helium
                    # line is smaller than the guess for the broad-line model.
                    if linemodel[iparam]['isbalmer'] and 'sigma' in param:
                        if linemodel[iparam]['isbroad']:
                            linemodel['initial'][iparam] = 1.1 * initial_guesses[param]
                        else:
                            linemodel['initial'][iparam] = initial_guesses[param]
                    else:
                        linemodel['initial'][iparam] = initial_guesses[param]
                    if param in param_bounds.keys():
                        linemodel['bounds'][iparam] = param_bounds[param]
                        # set the lower boundary on broad lines to be XX times the local noise
                        if linemodel['isbalmer'][iparam] and linemodel['isbroad'][iparam]:
                            civarkey = linemodel['linename'][iparam]+'_civar'
                            linemodel['civar'][iparam] = initial_guesses[civarkey]

                            #broadbalmer_snrmin = 3.
                            #linemodel['bounds'][iparam][0] = broadbalmer_snrmin * initial_guesses[noisekey]
                            #if linemodel['initial'][iparam] < linemodel['bounds'][iparam][0]:
                            #    linemodel['initial'][iparam] = linemodel['bounds'][iparam][0]
                            #if linemodel['initial'][iparam] > linemodel['bounds'][iparam][1]:
                            #    linemodel['initial'][iparam] = linemodel['bounds'][iparam][1]
            else:
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0
                else:
                    linemodel['initial'][iparam] = 1.0

            # Check bounds for free parameters but do not crash.
            if linemodel['fixed'][iparam] == False and linemodel['tiedtoparam'][iparam] == -1:
                toosml = linemodel['initial'][iparam] < linemodel['bounds'][iparam, 0]
                toobig = linemodel['initial'][iparam] > linemodel['bounds'][iparam, 1]
                if toosml:
                    errmsg = 'Initial parameter {} is outside its bound, {:.2f} < {:.2f}.'.format(
                        param, linemodel['initial'][iparam], linemodel['bounds'][iparam, 0])
                    log.warning(errmsg)
                    #raise ValueError(errmsg)
                    linemodel['initial'][iparam] = linemodel['bounds'][iparam, 0]
                if toobig:
                    errmsg = 'Initial parameter {} is outside its bound, {:.2f} > {:.2f}.'.format(
                        param, linemodel['initial'][iparam], linemodel['bounds'][iparam, 1])
                    log.warning(errmsg)
                    #raise ValueError(errmsg)
                    linemodel['initial'][iparam] = linemodel['bounds'][iparam, 1]

        # Now loop back through and ensure that tied relationships are enforced.
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        if len(Itied) > 0:
            for iparam, param in enumerate(linemodel['param_name'][Itied]):
                tieindx = linemodel['tiedtoparam'][Itied[iparam]]
                tiefactor = linemodel['tiedfactor'][Itied[iparam]]
                #log.info('{} tied to {} with factor {:.4f}'.format(param, linemodel[tieindx]['param_name'], tiefactor))
                linemodel['initial'][Itied[iparam]] = linemodel[tieindx]['initial'] * tiefactor

        linemodel['value'] = linemodel['initial'] # copy

        
    def _drop_params(self, parameters, linemodel, Ifree, log):
        """Drop dubious free parameters after fitting.
    
        """
        # Conditions for dropping a parameter (all parameters, not just those
        # being fitted):
        # --negative amplitude or sigma
        # --parameter at its default value (fit failed, right??)
        # --parameter within 0.1% of its bounds
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
        notfixed = np.logical_not(linemodel['fixed'])
    
        # drop any negative amplitude or sigma parameter that is not fixed 
        drop1 = np.hstack((lineamps < 0, np.zeros(len(linevshifts), bool), linesigmas <= 0)) * notfixed
        
        # Require equality, not np.isclose, because the optimization can be very
        # small (<1e-6) but still significant, especially for the doublet
        # ratios. If linesigma is dropped this way, make sure the corresponding
        # line-amplitude is dropped, too (see MgII 2796 on
        # sv1-bright-17680-39627622543528153).
        drop2 = np.zeros(len(parameters), bool)
            
        # if any amplitude is zero, drop the corresponding sigma and vshift
        amp_param_bool = self.amp_param_bool[Ifree]
        I = np.where(parameters[Ifree][amp_param_bool] == 0.)[0]
        if len(I) > 0:
            _Ifree = np.zeros(len(parameters), bool)
            _Ifree[Ifree] = True
            for pp in linemodel[Ifree][amp_param_bool][I]['param_name']:
                J = np.where(_Ifree * (linemodel['param_name'] == pp.replace('_amp', '_sigma')))[0]
                drop2[J] = True
                K = np.where(_Ifree * (linemodel['param_name'] == pp.replace('_amp', '_vshift')))[0]
                drop2[K] = True
                #print(pp, J, K, np.sum(drop2))
    
        # drop amplitudes for any lines tied to a line with a dropped sigma
        sigmadropped = np.where(self.sigma_param_bool * drop2)[0]
        for lineindx, dropline in zip(sigmadropped, linemodel[sigmadropped]['linename']):
            # Check whether lines are tied to this line. If so, find the
            # corresponding amplitude and drop that, too.
            T = linemodel['tiedtoparam'] == lineindx
            for tiedline in set(linemodel['linename'][T]):
                drop2[linemodel['param_name'] == f'{tiedline}_amp'] = True
            drop2[linemodel['param_name'] == f'{dropline}_amp'] = True
    
        # drop amplitudes for any lines tied to a line with a dropped vshift
        vshiftdropped = np.where(self.vshift_param_bool * drop2)[0]
        for lineindx, dropline in zip(vshiftdropped, linemodel[vshiftdropped]['linename']):
            # Check whether lines are tied to this line. If so, find the
            # corresponding amplitude and drop that, too.
            T = linemodel['tiedtoparam'] == lineindx
            for tiedline in set(linemodel['linename'][T]):
                drop2[linemodel['param_name'] == f'{tiedline}_amp'] = True
            drop2[linemodel['param_name'] == f'{dropline}_amp'] = True
    
        # drop any non-fixed parameters outside their bounds
        # It's OK for parameters to be *at* their bounds.
        drop3 = np.zeros(len(parameters), bool)
        drop3[Ifree] = np.logical_or(parameters[Ifree] < linemodel['bounds'][Ifree, 0], 
                                     parameters[Ifree] > linemodel['bounds'][Ifree, 1])
        drop3 *= notfixed
        
        log.debug(f'Dropping {np.sum(drop1)} negative-amplitude lines.') # linewidth can't be negative
        log.debug(f'Dropping {np.sum(drop2)} sigma,vshift parameters of zero-amplitude lines.')
        log.debug(f'Dropping {np.sum(drop3)} parameters which are out-of-bounds.')
        Idrop = np.where(np.logical_or.reduce((drop1, drop2, drop3)))[0]
        
        if len(Idrop) > 0:
            log.debug(f'  Dropping {len(Idrop)} unique parameters.')
            parameters[Idrop] = 0.0


    def optimize(self, linemodel,
                 obs_bin_centers,
                 obs_bin_fluxes,
                 obs_weights,
                 redshift,
                 resolution_matrices,
                 camerapix,
                 get_finalamp=False,
                 log=None,
                 verbose=False,
                 debug=False):
        """Optimization routine.
    
        """
        from scipy.optimize import least_squares
        
        if log is None:
            from desiutil.log import get_logger, DEBUG
            if verbose:
                log = get_logger(DEBUG)
            else:
                log = get_logger()
                
        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, doubletpair, \
                     line_wavelengths) = self._linemodel_to_parameters(linemodel, self.fit_linetable)

        params_mapping = EMLine_ParamsMapping(parameters, Ifree,
                                              Itied, tiedtoparam, tiedfactor,
                                              doubletindx, doubletpair)
        
        log.debug(f"Optimizing {len(Ifree)} free parameters")

        if len(Ifree) == 0:
            fit_info = {'nfev': 0, 'status': 0}
        else:
            # corner case where all lines are out of the wavelength range, which can
            # happen at high redshift and with the red camera masked, e.g.,
            # iron/main/dark/6642/39633239580608311).
            initial_guesses = parameters[Ifree]
    
            obj = EMLine_Objective(obs_bin_centers,
                                   obs_bin_fluxes,
                                   obs_weights,
                                   redshift,
                                   line_wavelengths,
                                   resolution_matrices,
                                   tuple(camerapix), # for more efficient iteration
                                   params_mapping)
        
            try:
                fit_info = least_squares(obj.objective, initial_guesses, jac=obj.jacobian, args=(),
                                         max_nfev=5000, xtol=1e-10, ftol=1e-5, #x_scale='jac' gtol=1e-10,
                                         tr_solver='lsmr', tr_options={'maxiter': 1000, 'regularize': True},
                                         method='trf', bounds=tuple(zip(*bounds)),) # verbose=2)
                parameters[Ifree] = fit_info.x
            except:
                if self.uniqueid:
                    errmsg = f'Problem in scipy.optimize.least_squares for {self.uniqueid}.'
                else:
                    errmsg = 'Problem in scipy.optimize.least_squares.'
                log.critical(errmsg)
                raise RuntimeError(errmsg)
    
            # Drop (zero out) any dubious free parameters.
            self._drop_params(parameters, linemodel, Ifree, log)
        
            # At this point, parameters contains correct *free* and *fixed* values, but
            # we need to update *tied* values to reflect any changes to free params.  We
            # do *not* apply doublet rules, as other code expects us to return a params
            # array with doublet ratios as ratios, not amplitudes.
            for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
                parameters[I] = parameters[indx] * factor
        
        out_linemodel = linemodel.copy()
        out_linemodel['value'] = parameters.copy() # so we don't munge it below
        out_linemodel.meta['nfev'] = fit_info['nfev']
        out_linemodel.meta['status'] = fit_info['status']
        
        if get_finalamp:
            
            lineamps, _, __ = np.array_split(parameters, 3) # 3 parameters per line
    
            # apply doublet rules
            lineamps[doubletindx] *= lineamps[doubletpair]
            
            # calculate the observed maximum amplitude for each
            # fitted spectral line after convolution with the resolution
            # matrix.
            peaks = EMLine_find_peak_amplitudes(parameters,
                                                obs_bin_centers,
                                                redshift,
                                                line_wavelengths,
                                                resolution_matrices,
                                                camerapix)
    
            # FIXME: is len(lineamps) == # lines obtainable from line table?
            out_linemodel['obsvalue'][:len(lineamps)] = peaks
            
        return out_linemodel
    
    @staticmethod
    def chi2(linemodel, emlinewave, emlineflux, emlineivar, emlinemodel,
             continuum_model=None, return_dof=False):
        """Compute the reduced chi^2."""

        nfree = np.count_nonzero((linemodel['fixed'] == False) * (linemodel['tiedtoparam'] == -1))
        dof = np.count_nonzero(emlineivar > 0) - nfree

        if dof > 0:
            if continuum_model is None:
                model = emlinemodel
            else:
                model = emlinemodel + continuum_model
            chi2 = np.sum(emlineivar * (emlineflux - model)**2) / dof
        else:
            chi2 = 0.0
            
        if return_dof:
            return chi2, dof, nfree
        else:
            return chi2


    def bestfit(self, linemodel, redshift, emlinewave, resolution_matrix, camerapix):
        """Construct the best-fitting emission-line spectrum from a linemodel."""
        
        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, \
                     doubletpair, linewaves) = self._linemodel_to_parameters(linemodel, self.fit_linetable)
        
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
        
        # doublets
        lineamps[doubletindx] *= lineamps[doubletpair]
    
        emlinemodel = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas,
                                         linewaves, emlinewave, resolution_matrix, camerapix)
    
    
        return emlinemodel
    

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table, redshift=None, 
                            camerapix=None, snrcut=None):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (used for QA and elsewhere).

        """
        if redshift is None:
            redshift = fastspecfit_table['Z']
        
        linewaves = self.linetable['restwave'].data

        parameters = []
        for param in self.param_names:
            if '_amp' in param:
                param = param.replace('_amp', '_modelamp')
            parameters.append(fastspecfit_table[param.upper()])
        #parameters = [fastspecfit_table[param.upper()] for param in self.param_names]

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line    

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        lineamps[self.doubletindx] *= lineamps[self.doubletpair]

        if snrcut is not None:
            lineamps_ivar = [fastspecfit_table[param.upper()+'_AMP_IVAR'] for param in self.linetable['name']]
            lineamps[lineamps * np.sqrt(lineamps_ivar) < snrcut] = 0.

        emlinemodel = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas,
                                         linewaves, np.hstack(specwave), specres, camerapix)
        
        return emlinemodel

    def populate_emtable(self, result, finalfit, finalmodel, emlinewave, emlineflux,
                         emlineivar, oemlineivar, specflux_nolines, redshift,
                         resolution_matrix, camerapix, log, nminpix=7, nsigma=3.):
        """Populate the output table with the emission-line measurements.

        """
        from scipy.stats import sigmaclip
        from fastspecfit.util import centers2edges
        from scipy.special import erf

        for param in finalfit:
            val = param['value']
            obsval = param['obsvalue']

            # special case the tied doublets
            if param['param_name'] == 'oii_doublet_ratio':
                result['OII_DOUBLET_RATIO'] = val
                result['OII_3726_MODELAMP'] = val * finalfit[param['doubletpair']]['value']
                result['OII_3726_AMP'] = val * finalfit[param['doubletpair']]['obsvalue']
            elif param['param_name'] == 'sii_doublet_ratio':
                result['SII_DOUBLET_RATIO'] = val
                result['SII_6731_MODELAMP'] = val * finalfit[param['doubletpair']]['value']
                result['SII_6731_AMP'] = val * finalfit[param['doubletpair']]['obsvalue']
            elif param['param_name'] == 'mgii_doublet_ratio':
                result['MGII_DOUBLET_RATIO'] = val
                result['MGII_2796_MODELAMP'] = val * finalfit[param['doubletpair']]['value']
                result['MGII_2796_AMP'] = val * finalfit[param['doubletpair']]['obsvalue']
            else:
                if '_amp' in param['param_name']:
                    result[param['param_name'].upper().replace('AMP', 'MODELAMP')] = val
                    result[param['param_name'].upper()] = obsval
                else:
                    result[param['param_name'].upper()] = val                    

        gausscorr = erf(nsigma / np.sqrt(2))      # correct for the flux outside of +/-nsigma
        dpixwave = np.median(np.diff(emlinewave)) # median pixel size [Angstrom]

        # Where the cameras overlap, we have to account for the variable pixel
        # size by sorting in wavelength.
        Wsrt = np.argsort(emlinewave)
        dwaves = np.diff(centers2edges(emlinewave[Wsrt]))

        # zero out all out-of-range lines
        for oneline in self.fit_linetable[~self.fit_linetable['inrange']]:
            linename = oneline['name'].upper()
            #print(linename, result['{}_AMP'.format(linename)], result['{}_MODELAMP'.format(linename)],
            #      result['{}_SIGMA'.format(linename)], result['{}_VSHIFT'.format(linename)])
            result['{}_AMP'.format(linename)] = 0.0
            result['{}_MODELAMP'.format(linename)] = 0.0
            result['{}_VSHIFT'.format(linename)] = 0.0
            result['{}_SIGMA'.format(linename)] = 0.0

        # get continuum fluxes, EWs, and upper limits
        narrow_sigmas, broad_sigmas, uv_sigmas = [], [], []
        narrow_redshifts, broad_redshifts, uv_redshifts = [], [], []
        for oneline in self.fit_linetable[self.fit_linetable['inrange']]:

            linename = oneline['name'].upper()
            linez = redshift + result['{}_VSHIFT'.format(linename)] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)
            linesigma = result['{}_SIGMA'.format(linename)] # [km/s]

            # if the line was dropped, use a default sigma value
            if linesigma == 0:
                if oneline['isbroad']:
                    if oneline['isbalmer']:
                        linesigma = self.limitsigma_narrow
                    else:
                        linesigma = self.limitsigma_broad
                else:
                    linesigma = self.limitsigma_broad

            linesigma_ang = linesigma * linezwave / C_LIGHT    # [observed-frame Angstrom]

            # require at least 2 pixels
            if linesigma_ang < 2 * dpixwave:
                linesigma_ang_window = 2 * dpixwave
                _gausscorr = 1.
            else:
                linesigma_ang_window = linesigma_ang
                _gausscorr = gausscorr

            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next line.
            lineindx = np.where((emlinewave >= (linezwave - nsigma*linesigma_ang_window)) *
                                (emlinewave <= (linezwave + nsigma*linesigma_ang_window)))[0]
            
            if len(lineindx) > 0 and np.sum(oemlineivar[lineindx] == 0) / len(lineindx) > 0.3: # use original ivar
                result['{}_AMP'.format(linename)] = 0.0
                result['{}_MODELAMP'.format(linename)] = 0.0
                result['{}_VSHIFT'.format(linename)] = 0.0
                result['{}_SIGMA'.format(linename)] = 0.0
            else:
                # number of pixels, chi2, and boxcar integration
                lineindx = np.where((emlinewave >= (linezwave - nsigma*linesigma_ang_window)) *
                                    (emlinewave <= (linezwave + nsigma*linesigma_ang_window)) *
                                    (emlineivar > 0))[0]

                npix = len(lineindx)
                result['{}_NPIX'.format(linename)] = npix
    
                if npix >= nminpix: # magic number: required at least XX unmasked pixels centered on the line
                    
                    if np.any(emlineivar[lineindx] == 0):
                        errmsg = 'Ivar should never be zero within an emission line!'
                        log.critical(errmsg)
                        raise ValueError(errmsg)

                    lineindx_Wsrt = np.where((emlinewave[Wsrt] >= (linezwave - nsigma*linesigma_ang_window)) *
                                             (emlinewave[Wsrt] <= (linezwave + nsigma*linesigma_ang_window)) *
                                             (emlineivar[Wsrt] > 0))[0]

                    # boxcar integration of the flux
                    #dwave = np.median(np.abs(np.diff(emlinewave_edges[lineindx])))
                    boxflux = np.sum(emlineflux[Wsrt][lineindx_Wsrt] * dwaves[lineindx_Wsrt])
                    boxflux_ivar = 1 / np.sum((1 / emlineivar[Wsrt][lineindx_Wsrt]) * dwaves[lineindx_Wsrt]**2)

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
    
                    # require amp > 0 (line not dropped) to compute the flux and chi2
                    if result['{}_MODELAMP'.format(linename)] > 0:
    
                        result['{}_CHI2'.format(linename)] = np.sum(emlineivar[lineindx] * (emlineflux[lineindx] - finalmodel[lineindx])**2)

                        print('ToDo: need the per-line model here.')
                        lineprofile = np.ones_like(emlinewave)
                        #lineprofile = build_emline_model(self.dlog10wave, redshift,
                        #                                 np.array([result['{}_MODELAMP'.format(linename)]]),
                        #                                 np.array([result['{}_VSHIFT'.format(linename)]]),
                        #                                 np.array([result['{}_SIGMA'.format(linename)]]),
                        #                                 np.array([oneline['restwave']]), emlinewave,
                        #                                 resolution_matrix, camerapix, debug=False)
                        
                        if np.sum(lineprofile) == 0. or np.any(lineprofile) < 0.:
                            errmsg = 'Line-profile should never be zero or negative!'
                            log.critical(errmsg)
                            raise ValueError(errmsg)
                        
                        # theoretical Gaussian line-flux and the (wrong) weighted flux_ivar
                        #flux = result['{}_MODELAMP'.format(linename)] * np.sqrt(2.0 * np.pi) * linesigma_ang # * u.Angstrom
                        #flux_ivar = np.sum(lineprofile[lineindx])**2 / np.sum(lineprofile[lineindx]**2 / emlineivar[lineindx])

                        # matched-filter (maximum-likelihood) Gaussian flux
                        pro_j = lineprofile[Wsrt][lineindx_Wsrt] / np.sum(lineprofile[Wsrt][lineindx_Wsrt])
                        I = pro_j > 0. # very narrow lines can have a profile that goes to zero
                        weight_j = (pro_j[I]**2 / dwaves[lineindx_Wsrt][I]**2) * emlineivar[Wsrt][lineindx_Wsrt][I]
                        flux = np.sum(weight_j * dwaves[lineindx_Wsrt][I] * lineprofile[Wsrt][lineindx_Wsrt][I] / pro_j[I]) / np.sum(weight_j)
                        flux_ivar = np.sum(weight_j)

                        # correction for missing flux
                        flux /= _gausscorr
                        flux_ivar *= _gausscorr**2

                        result['{}_FLUX'.format(linename)] = flux
                        result['{}_FLUX_IVAR'.format(linename)] = flux_ivar # * u.second**2*u.cm**4/u.erg**2
    
                        # keep track of sigma and z but only using XX-sigma lines
                        linesnr = result['{}_AMP'.format(linename)] * np.sqrt(result['{}_AMP_IVAR'.format(linename)])
                        #print(linename, result['{}_AMP'.format(linename)], amp_sigma, linesnr)
                        if linesnr > 1.5:
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
                indxlo = np.where((emlinewave > (linezwave - 10*linesigma_ang_window)) *
                                  (emlinewave < (linezwave - 3.*linesigma_ang_window)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indxhi = np.where((emlinewave < (linezwave + 10*linesigma_ang_window)) *
                                  (emlinewave > (linezwave + 3.*linesigma_ang_window)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indx = np.hstack((indxlo, indxhi))

                if len(indx) >= nminpix: # require at least XX pixels to get the continuum level
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
    
                    result['{}_CONT'.format(linename)] = cmed # * u.erg/(u.second*u.cm**2*u.Angstrom)
                    result['{}_CONT_IVAR'.format(linename)] = civar # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
    
                if result['{}_CONT'.format(linename)] != 0.0 and result['{}_CONT_IVAR'.format(linename)] != 0.0:
                    lineflux = result['{}_FLUX'.format(linename)]
                    #linefluxivar = result['{}_BOXFLUX_IVAR'.format(linename)]
                    linefluxivar = result['{}_FLUX_IVAR'.format(linename)]
                    if lineflux > 0 and linefluxivar > 0:
                        # add the uncertainties in flux and the continuum in quadrature
                        ew = lineflux / cmed / (1 + redshift) # rest frame [A]
                        ewivar = (1+redshift)**2 / (1 / (cmed**2 * linefluxivar) + lineflux**2 / (cmed**4 * civar))
                    else:
                        ew, ewivar = 0.0, 0.0
                        
                    # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                    fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar) # * u.erg/(u.second*u.cm**2)
                    ewlimit = fluxlimit * cmed / (1+redshift)
    
                    result['{}_EW'.format(linename)] = ew
                    result['{}_EW_IVAR'.format(linename)] = ewivar
                    result['{}_FLUX_LIMIT'.format(linename)] = fluxlimit 
                    result['{}_EW_LIMIT'.format(linename)] = ewlimit

            if 'debug' in log.name:
                for col in ('VSHIFT', 'SIGMA', 'MODELAMP', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)]))
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'BOXFLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT', 'EW_LIMIT'):
                    log.debug('{} {}: {:.4f}'.format(linename, col, result['{}_{}'.format(linename, col)]))
                print()
                #log.debug(' ')
    
            ## simple QA
            #if linename == 'OIII_5007':
            #    import matplotlib.pyplot as plt
            #    _indx = np.arange(indx[-1]-indx[0])+indx[0]
            #    # continuum bandpasses and statistics
            #    plt.clf()
            #    plt.plot(emlinewave[_indx], specflux_nolines[_indx], color='gray')
            #    plt.scatter(emlinewave[indx], specflux_nolines[indx], color='red')
            #    plt.axhline(y=cmed, color='k')
            #    plt.axhline(y=cmed+1/np.sqrt(civar), color='k', ls='--')
            #    plt.axhline(y=cmed-1/np.sqrt(civar), color='k', ls='--')
            #    plt.savefig('desi-users/ioannis/tmp/junk.png')
            #
            #    # emission-line integration
            #    plt.clf()
            #    plt.plot(emlinewave[_indx], emlineflux[_indx], color='gray')
            #    plt.plot(emlinewave[_indx], finalmodel[_indx], color='red')
            #    #plt.plot(emlinewave[_indx], specflux_nolines[_indx], color='orange', alpha=0.5)
            #    plt.axvline(x=emlinewave[lineindx[0]], color='blue')
            #    plt.axvline(x=emlinewave[lineindx[-1]], color='blue')
            #    plt.axhline(y=0, color='k', ls='--')
            #    plt.axhline(y=amp_sigma, color='k', ls='--')
            #    plt.axhline(y=2*amp_sigma, color='k', ls='--')
            #    plt.axhline(y=3*amp_sigma, color='k', ls='--')
            #    plt.axhline(y=result['{}_AMP'.format(linename)], color='k', ls='-')
            #    plt.savefig('desi-users/ioannis/tmp/junk2.png')

        # Clean up the doublets whose amplitudes were tied in the fitting since
        # they may have been zeroed out in the clean-up, above. This should be
        # smarter.
        if result['OIII_5007_MODELAMP'] == 0.0 and result['OIII_5007_NPIX'] > 0:
            result['OIII_4959_MODELAMP'] = 0.0
            result['OIII_4959_AMP'] = 0.0
            result['OIII_4959_FLUX'] = 0.0
            result['OIII_4959_EW'] = 0.0
        if result['NII_6584_MODELAMP'] == 0.0 and result['NII_6584_NPIX'] > 0:
            result['NII_6548_MODELAMP'] = 0.0
            result['NII_6548_AMP'] = 0.0
            result['NII_6548_FLUX'] = 0.0
            result['NII_6548_EW'] = 0.0
        if result['OII_7320_MODELAMP'] == 0.0 and result['OII_7320_NPIX'] > 0:
            result['OII_7330_MODELAMP'] = 0.0
            result['OII_7330_AMP'] = 0.0
            result['OII_7330_FLUX'] = 0.0
            result['OII_7330_EW'] = 0.0
        if result['MGII_2796_MODELAMP'] == 0.0 and result['MGII_2803_MODELAMP'] == 0.0:
            result['MGII_DOUBLET_RATIO'] = 0.0
        if result['OII_3726_MODELAMP'] == 0.0 and result['OII_3729_MODELAMP'] == 0.0:
            result['OII_DOUBLET_RATIO'] = 0.0
        if result['SII_6716_MODELAMP'] == 0.0 and result['SII_6731_MODELAMP'] == 0.0:
            result['SII_DOUBLET_RATIO'] = 0.0

        if 'debug' in log.name:
            for col in ('MGII_DOUBLET_RATIO', 'OII_DOUBLET_RATIO', 'SII_DOUBLET_RATIO'):
                log.debug('{}: {:.4f}'.format(col, result[col]))
            #log.debug(' ')
            print()

        # get the average emission-line redshifts and velocity widths
        if len(narrow_redshifts) > 0:
            result['NARROW_Z'] = np.median(narrow_redshifts)
            result['NARROW_SIGMA'] = np.median(narrow_sigmas) # * u.kilometer / u.second
            result['NARROW_ZRMS'] = np.std(narrow_redshifts)
            result['NARROW_SIGMARMS'] = np.std(narrow_sigmas)
        else:
            result['NARROW_Z'] = redshift
            
        if len(broad_redshifts) > 0:
            result['BROAD_Z'] = np.median(broad_redshifts)
            result['BROAD_SIGMA'] = np.median(broad_sigmas) # * u.kilometer / u.second
            result['BROAD_ZRMS'] = np.std(broad_redshifts)
            result['BROAD_SIGMARMS'] = np.std(broad_sigmas)
        else:
            result['BROAD_Z'] = redshift

        if len(uv_redshifts) > 0:
            result['UV_Z'] = np.median(uv_redshifts)
            result['UV_SIGMA'] = np.median(uv_sigmas) # * u.kilometer / u.second
            result['UV_ZRMS'] = np.std(uv_redshifts)
            result['UV_SIGMARMS'] = np.std(uv_sigmas)
        else:
            result['UV_Z'] = redshift

        # fragile
        if 'debug' in log.name:
            for line in ('NARROW', 'BROAD', 'UV'):
                log.debug('{}_Z: {:.9f}+/-{:.9f}'.format(line, result['{}_Z'.format(line)], result['{}_ZRMS'.format(line)]))
                log.debug('{}_SIGMA: {:.3f}+/-{:.3f}'.format(line, result['{}_SIGMA'.format(line)], result['{}_SIGMARMS'.format(line)]))

    def synthphot_spectrum(self, data, result, modelwave, modelflux):
        """Synthesize photometry from the best-fitting model (continuum+emission lines).

        """
        filters = self.synth_filters[data['photsys']]

        # Pad (simply) in wavelength...
        padflux, padwave = filters.pad_spectrum(modelflux, modelwave, method='edge')
        synthmaggies = filters.get_ab_maggies(padflux / FLUXNORM, padwave)
        synthmaggies = synthmaggies.as_array().view('f8')
        model_synthphot = self.parse_photometry(self.synth_bands, maggies=synthmaggies,
                                                nanomaggies=False,
                                                lambda_eff=filters.effective_wavelengths.value)

        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_{}'.format(band.upper())] = data['synthphot']['nanomaggies'][iband] # * 'nanomaggies'
            #result['FLUX_SYNTH_IVAR_{}'.format(band.upper())] = data['synthphot']['nanomaggies_ivar'][iband]
        for iband, band in enumerate(self.synth_bands):
            result['FLUX_SYNTH_SPECMODEL_{}'.format(band.upper())] = model_synthphot['nanomaggies'][iband] # * 'nanomaggies'


def emline_specfit(data, result, continuummodel, smooth_continuum,
                   minsnr_balmer_broad=3., fphoto=None, emlinesfile=None,
                   synthphot=True, broadlinefit=True,
                   percamera_models=False, log=None, verbose=False):
    """Perform the fit minimization / chi2 minimization.

    Parameters
    ----------
    data
    continuummodel
    smooth_continuum
    synthphot
    verbose
    broadlinefit

    Returns
    -------
    results
    modelflux
 
    """
    from astropy.table import Column, vstack
    from fastspecfit.util import ivar2var

    tall = time.time()

    if log is None:
        from desiutil.log import get_logger, DEBUG
        if verbose:
            log = get_logger(DEBUG)
        else:
            log = get_logger()

    EMFit = EMFitTools(emlinesfile=emlinesfile, fphoto=fphoto, uniqueid=data['uniqueid'],
                       minsnr_balmer_broad=minsnr_balmer_broad)

    # Combine all three cameras; we will unpack them to build the
    # best-fitting model (per-camera) below.
    redshift = data['zredrock']
    emlinewave = np.hstack(data['wave'])
    oemlineivar = np.hstack(data['ivar'])
    specflux = np.hstack(data['flux'])
    resolution_matrix = data['res_fast']
    camerapix = data['camerapix']

    continuummodelflux = np.hstack(continuummodel)
    smoothcontinuummodelflux = np.hstack(smooth_continuum)
    emlineflux = specflux - continuummodelflux - smoothcontinuummodelflux

    emlineivar = np.copy(oemlineivar)
    emlinevar, emlinegood = ivar2var(emlineivar, clip=1e-3)
    emlinebad = np.logical_not(emlinegood)

    # This is a (dangerous???) hack.
    if np.sum(emlinebad) > 0:
        emlineivar[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineivar[emlinegood])
        emlineflux[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineflux[emlinegood]) # ???

    weights = np.sqrt(emlineivar)

    # Build all the emission-line models for this object.
    linemodel_broad, linemodel_nobroad = EMFit.build_linemodels(
        redshift, wavelims=(np.min(emlinewave)+5, np.max(emlinewave)-5),
        verbose=False, strict_broadmodel=True)
    #EMFit.summarize_linemodel(linemodel_broad)
    #EMFit.summarize_linemodel(linemodel_nobroad)

    # Get initial guesses on the parameters and populate the two "initial"
    # linemodels; the "final" linemodels will be initialized with the
    # best-fitting parameters from the initial round of fitting.
    initial_guesses, param_bounds = EMFit.initial_guesses_and_bounds(
        data, emlinewave, emlineflux, log)

    EMFit.populate_linemodel(linemodel_nobroad, initial_guesses, param_bounds, log)
    EMFit.populate_linemodel(linemodel_broad, initial_guesses, param_bounds, log)

    # Initial fit - initial_linemodel_nobroad
    t0 = time.time()
    fit_nobroad = EMFit.optimize(linemodel_nobroad, emlinewave, emlineflux, weights, redshift,
                                 resolution_matrix, camerapix, log=log, debug=False, get_finalamp=True)
    model_nobroad = EMFit.bestfit(fit_nobroad, redshift, emlinewave, resolution_matrix, camerapix)
    chi2_nobroad, ndof_nobroad, nfree_nobroad = EMFit.chi2(fit_nobroad, emlinewave, emlineflux, emlineivar, model_nobroad, return_dof=True)
    log.info('Line-fitting with no broad lines and {} free parameters took {:.4f} seconds [niter={}, rchi2={:.4f}].'.format(
        nfree_nobroad, time.time()-t0, fit_nobroad.meta['nfev'], chi2_nobroad))
    
    # Now try adding broad Balmer and helium lines and see if we improve the
    # chi2.
    if broadlinefit:
        # Gather the pixels around the broad Balmer lines and the corresponding
        # linemodel table.
        balmer_pix, balmer_linemodel_broad, balmer_linemodel_nobroad = [], [], []
        for icam in np.arange(len(data['cameras'])):
            pixoffset = int(np.sum(data['npixpercamera'][:icam]))
            if len(data['linename'][icam]) > 0:
                I = (linemodel_nobroad['isbalmer'] * (linemodel_nobroad['ishelium'] == False) *
                     linemodel_nobroad['isbroad'] * np.isin(linemodel_nobroad['linename'], data['linename'][icam]))
                _balmer_linemodel_broad = linemodel_broad[I]
                _balmer_linemodel_nobroad = linemodel_nobroad[I]
                balmer_linemodel_broad.append(_balmer_linemodel_broad)
                balmer_linemodel_nobroad.append(_balmer_linemodel_nobroad)
                if len(_balmer_linemodel_broad) > 0: # use balmer_linemodel_broad not balmer_linemodel_nobroad
                    I = np.where(np.isin(data['linename'][icam], _balmer_linemodel_broad['linename']))[0]
                    for ii in I:
                        #print(data['linename'][icam][ii])
                        balmer_pix.append(data['linepix'][icam][ii] + pixoffset)
                        
        if len(balmer_pix) > 0:
            t0 = time.time()
            fit_broad = EMFit.optimize(linemodel_broad, emlinewave, emlineflux, weights, 
                                       redshift, resolution_matrix, camerapix, log=log,
                                       debug=False, get_finalamp=True)
            model_broad = EMFit.bestfit(fit_broad, redshift, emlinewave, resolution_matrix, camerapix)
            chi2_broad, ndof_broad, nfree_broad = EMFit.chi2(fit_broad, emlinewave, emlineflux, emlineivar, model_broad, return_dof=True)
            log.info('Line-fitting with broad lines and {} free parameters took {:.4f} seconds [niter={}, rchi2={:.4f}].'.format(
                nfree_broad, time.time()-t0, fit_broad.meta['nfev'], chi2_broad))

            # compute delta-chi2 around just the Balmer lines
            balmer_pix = np.hstack(balmer_pix)
            balmer_linemodel_broad = vstack(balmer_linemodel_broad)

            balmer_nfree_broad = (np.count_nonzero((balmer_linemodel_broad['fixed'] == False) *
                                                   (balmer_linemodel_broad['tiedtoparam'] == -1)))
            balmer_ndof_broad = np.count_nonzero(emlineivar[balmer_pix] > 0) - balmer_nfree_broad

            balmer_linemodel_nobroad = vstack(balmer_linemodel_nobroad)
            balmer_nfree_nobroad = (np.count_nonzero((balmer_linemodel_nobroad['fixed'] == False) *
                                                     (balmer_linemodel_nobroad['tiedtoparam'] == -1)))
            balmer_ndof_nobroad = np.count_nonzero(emlineivar[balmer_pix] > 0) - balmer_nfree_nobroad

            linechi2_balmer_broad = np.sum(emlineivar[balmer_pix] * (emlineflux[balmer_pix] - model_broad[balmer_pix])**2)
            linechi2_balmer_nobroad = np.sum(emlineivar[balmer_pix] * (emlineflux[balmer_pix] - model_nobroad[balmer_pix])**2)
            delta_linechi2_balmer = linechi2_balmer_nobroad - linechi2_balmer_broad
            delta_linendof_balmer = balmer_ndof_nobroad - balmer_ndof_broad

            # Choose broad-line model only if:
            # --delta-chi2 > delta-ndof
            # --broad_sigma < narrow_sigma
            # --broad_sigma < 250

            dchi2test = delta_linechi2_balmer > delta_linendof_balmer
            Hanarrow = fit_broad['param_name'] == 'halpha_sigma' # Balmer lines are tied to H-alpha even if out of range
            Habroad = fit_broad['param_name'] == 'halpha_broad_sigma'
            Bbroad = fit_broad['isbalmer'] * fit_broad['isbroad'] * (fit_broad['fixed'] == False) * EMFit.amp_balmer_bool
            broadsnr = fit_broad[Bbroad]['obsvalue'].data * np.sqrt(fit_broad[Bbroad]['civar'].data)

            sigtest1 = fit_broad[Habroad]['value'][0] > EMFit.minsigma_balmer_broad
            sigtest2 = (fit_broad[Habroad]['value'] > fit_broad[Hanarrow]['value'])[0]
            if len(broadsnr) == 0:
                broadsnrtest = False
                _broadsnr = 0.
            elif len(broadsnr) == 1:
                broadsnrtest =  broadsnr[-1] > EMFit.minsnr_balmer_broad
                _broadsnr = 'S/N ({}) = {:.1f}'.format(fit_broad[Bbroad]['linename'][-1], broadsnr[-1])
            else:
                broadsnrtest =  np.any(broadsnr[-2:] > EMFit.minsnr_balmer_broad)
                _broadsnr = 'S/N ({}) = {:.1f}, S/N ({}) = {:.1f}'.format(
                    fit_broad[Bbroad]['linename'][-2], broadsnr[-2], fit_broad[Bbroad]['linename'][-1], broadsnr[-1])

            if dchi2test and sigtest1 and sigtest2 and broadsnrtest:
                log.info('Adopting broad-line model:')
                log.info('  delta-chi2={:.1f} > delta-ndof={:.0f}'.format(delta_linechi2_balmer, delta_linendof_balmer))
                log.info('  sigma_broad={:.1f} km/s, sigma_narrow={:.1f} km/s'.format(fit_broad[Habroad]['value'][0], fit_broad[Hanarrow]['value'][0]))
                if _broadsnr:
                    log.info('  {} > {:.0f}'.format(_broadsnr, EMFit.minsnr_balmer_broad))
                finalfit, finalmodel, finalchi2 = fit_broad, model_broad, chi2_broad
            else:
                if dchi2test == False:
                    log.info('Dropping broad-line model: delta-chi2={:.1f} < delta-ndof={:.0f}'.format(
                        delta_linechi2_balmer, delta_linendof_balmer))
                elif sigtest1 == False:
                    log.info('Dropping broad-line model: Halpha_broad_sigma {:.1f} km/s < {:.0f} km/s (delta-chi2={:.1f}, delta-ndof={:.0f}).'.format(
                        fit_broad[Habroad]['value'][0], EMFit.minsigma_balmer_broad, delta_linechi2_balmer, delta_linendof_balmer))
                elif sigtest2 == False:
                    log.info('Dropping broad-line model: Halpha_broad_sigma {:.1f} km/s < Halpha_narrow_sigma {:.1f} km/s (delta-chi2={:.1f}, delta-ndof={:.0f}).'.format(
                        fit_broad[Habroad]['value'][0], fit_broad[Hanarrow]['value'][0], delta_linechi2_balmer, delta_linendof_balmer))
                elif broadsnrtest == False:
                    log.info('Dropping broad-line model: {} < {:.0f}'.format(_broadsnr, EMFit.minsnr_balmer_broad))
                finalfit, finalmodel, finalchi2 = fit_nobroad, model_nobroad, chi2_nobroad
        else:
            log.info('Insufficient Balmer lines to test the broad-line model.')
            finalfit, finalmodel, finalchi2 = fit_nobroad, model_nobroad, chi2_nobroad
            delta_linechi2_balmer, delta_linendof_balmer = 0, np.int32(0)
    else:
        log.info('Skipping broad-line fitting (broadlinefit=False).')
        finalfit, finalmodel, finalchi2 = fit_nobroad, model_nobroad, chi2_nobroad
        delta_linechi2_balmer, delta_linendof_balmer = 0, np.int32(0)

    # Residual spectrum with no emission lines.
    specflux_nolines = specflux - finalmodel

    # Now fill the output table.
    EMFit.populate_emtable(result, finalfit, finalmodel, emlinewave, emlineflux,
                           emlineivar, oemlineivar, specflux_nolines, redshift,
                           resolution_matrix, camerapix, log)

    # Build the model spectra.
    emmodel = np.hstack(EMFit.emlinemodel_bestfit(data['wave'], data['res_fast'], result, 
                                                  camerapix=camerapix, redshift=redshift))

    result['RCHI2_LINE'] = finalchi2
    #result['NDOF_LINE'] = finalndof
    result['DELTA_LINECHI2'] = delta_linechi2_balmer # chi2_nobroad - chi2_broad
    result['DELTA_LINENDOF'] = delta_linendof_balmer # ndof_nobroad - ndof_broad

    # full-fit reduced chi2
    rchi2 = np.sum(oemlineivar * (specflux - (continuummodelflux + smoothcontinuummodelflux + emmodel))**2)
    rchi2 /= np.sum(oemlineivar > 0) # dof??
    result['RCHI2'] = rchi2

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
    dwave = float(int(np.round(dwave * 1000)) / 1000)
    npix = int(np.round((maxwave-minwave)/dwave)) + 1
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

    wavesrt = np.argsort(emlinewave)
    modelcontinuum = np.interp(modelwave, emlinewave[wavesrt], continuummodelflux[wavesrt]).reshape(1, npix)
    modelsmoothcontinuum = np.interp(modelwave, emlinewave[wavesrt], smoothcontinuummodelflux[wavesrt]).reshape(1, npix)
    modelemspectrum = np.interp(modelwave, emlinewave[wavesrt], emmodel[wavesrt]).reshape(1, npix)
    
    modelspectra.add_column(Column(name='CONTINUUM', dtype='f4', data=modelcontinuum))
    modelspectra.add_column(Column(name='SMOOTHCONTINUUM', dtype='f4', data=modelsmoothcontinuum))
    modelspectra.add_column(Column(name='EMLINEMODEL', dtype='f4', data=modelemspectrum))

    # Finally, optionally synthesize photometry (excluding the
    # smoothcontinuum!) and measure Dn(4000) from the line-free spectrum.
    if synthphot:
        modelflux = modelcontinuum[0, :] + modelemspectrum[0, :]
        EMFit.synthphot_spectrum(data, result, modelwave, modelflux)

    # measure DN(4000) without the emission lines
    if result['DN4000_IVAR'] > 0:
        fluxnolines = data['coadd_flux'] - modelemspectrum[0, :]
        dn4000_nolines, _ = EMFit.get_dn4000(modelwave, fluxnolines, redshift=redshift, log=log, rest=False)
        log.info('Dn(4000)={:.3f} in the emission-line subtracted spectrum.'.format(dn4000_nolines))
        result['DN4000'] = dn4000_nolines

        # Simple QA of the Dn(4000) estimate.
        if False:
            import matplotlib.pyplot as plt

            dn4000, dn4000_obs, dn4000_model, dn4000_ivar = result['DN4000'], result['DN4000_OBS'], result['DN4000_MODEL'], result['DN4000_IVAR']
            print(dn4000, dn4000_obs, dn4000_model, 1/np.sqrt(dn4000_ivar))
    
            restwave = modelwave / (1 + redshift) # [Angstrom]
            flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
            fnu_obs = data['coadd_flux'] * flam2fnu # [erg/s/cm2/Hz]
            fnu = fluxnolines * flam2fnu # [erg/s/cm2/Hz]

            fnu_model = modelcontinuum[0, :] * flam2fnu
            fnu_fullmodel = modelflux * flam2fnu
            
            fnu_ivar = data['coadd_ivar'] / flam2fnu**2            
            fnu_sigma, fnu_mask = ivar2var(fnu_ivar, sigma=True)
    
            I = (restwave > 3835) * (restwave < 4115)
            J = (restwave > 3835) * (restwave < 4115) * fnu_mask
    
            fig, ax = plt.subplots()
            ax.fill_between(restwave[I], fnu_obs[I]-fnu_sigma[I], fnu_obs[I]+fnu_sigma[I],
                            label='Observed Dn(4000)={:.3f}+/-{:.3f}'.format(dn4000_obs, 1/np.sqrt(dn4000_ivar)))
            ax.plot(restwave[I], fnu[I], color='blue', label='Line-free Dn(4000)={:.3f}+/-{:.3f}'.format(
                dn4000, 1/np.sqrt(dn4000_ivar)))
            ax.plot(restwave[I], fnu_fullmodel[I], color='k', label='Model Dn(4000)={:.3f}'.format(dn4000_model))
            ax.plot(restwave[I], fnu_model[I], color='red', label='Model Dn(4000)={:.3f}'.format(dn4000_model))
            ylim = ax.get_ylim()
            ax.fill_between([3850, 3950], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.fill_between([4000, 4100], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                            color='lightgray', alpha=0.5)
            ax.set_xlabel(r'Rest Wavelength ($\AA$)')
            ax.set_ylabel(r'$F_{\nu}$ (erg/s/cm2/Hz)')
            ax.legend()
            fig.savefig('desi-users/ioannis/tmp/qa-dn4000.png')

    log.info('Emission-line fitting took {:.2f} seconds.'.format(time.time()-tall))

    if percamera_models:
        errmsg = 'percamera-models option not yet implemented.'
        log.critical(errmsg)
        raise NotImplementedError(errmsg)

    return modelspectra

