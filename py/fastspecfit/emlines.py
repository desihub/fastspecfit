"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import pdb # for debugging

from enum import IntEnum
from itertools import chain

import os, time
import numpy as np
from numba import jit

from astropy.table import Table

from fastspecfit.io import (
    read_emlines,
    FLUXNORM,
)

from fastspecfit.continuum import Filters
from fastspecfit.util import C_LIGHT

from fastspecfit.emline_fit import (
    EMLine_Objective,
    EMLine_MultiLines,
    EMLine_find_peak_amplitudes,
    EMLine_build_model,
    EMLine_ParamsMapping,
)

class ParamType(IntEnum):
    AMPLITUDE = 0,
    VSHIFT = 1,
    SIGMA = 2


#def emline_patch_model(redshift, lineamps, linevshifts, linesigmas,
#                       linewaves, emlinewave, resolution_matrix, camerapix):
#    """Given parameters, build the model emission-line spectrum + continuum patches.
#
#    """
#
#    emlinemodel = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas,
#                                     linewaves, emlinewave, resolution_matrix, camerapix)
#
#    return emlinemodel
    
    
def emline_patch_objective(free_parameters, 
                           obs_bin_centers,
                           obs_bin_fluxes, 
                           obs_weights,
                           redshift,
                           line_wavelengths,
                           resolution_matrices,
                           camerapix,
                           allparameters,
                           nFree,
                           isFree,
                           isTied,
                           tiedtoparam,
                           tiedfactor,
                           doublet_idx, 
                           doublet_src,
                           nFree_patch,
                           isFree_patch,
                           continuum_patches,
                           nPatch, 
                           return_model):


    # emission lines
    free_line_parameters = free_parameters[:nFree]
    line_parameters = allparameters[:-2*nPatch]

    line_parameters[isFree] = free_line_parameters
    line_parameters[isTied] = tiedfactor * line_parameters[tiedtoparam]

    lineamps, linevshifts, linesigmas = np.array_split(line_parameters, 3) # 3 parameters per line
    lineamps[doublet_idx] *= lineamps[doublet_src] # doublets

    emlinemodel = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas,
                                     line_wavelengths, obs_bin_centers, resolution_matrices, 
                                     camerapix)
    
    # continuum patches
    free_patch_parameters = free_parameters[nFree:].reshape(2, nFree_patch // 2) # [slope,intercept] per (free) patch

    #import matplotlib.pyplot as plt
    #plt.clf()
    #plt.plot(obs_bin_centers, emlinemodel)

    continuummodel = np.zeros_like(emlinemodel)
    for ipatch, patch in enumerate(continuum_patches):
        s, e = patch['s'], patch['e']
        params = free_patch_parameters[:, ipatch] # [slope, intercept]
        patchmodel = params[0] * (obs_bin_centers[s:e] - patch['pivotwave']) + params[1]
        continuummodel[s:e] += patchmodel
        #plt.plot(obs_bin_centers[s:e], patchmodel, alpha=0.5)

    #plt.savefig('junk.png')

    patchmodel = emlinemodel + continuummodel
    if return_model:
        return np.hstack((lineamps, linevshifts, linesigmas)), free_patch_parameters, patchmodel

    #patchmodel = emline_patch_model(redshift, lineamps, linevshifts, linesigmas,
    #                                line_wavelengths, obs_bin_centers, resolution_matrices, 
    #                                camerapix)


    residuals = obs_weights * (patchmodel - obs_bin_fluxes)

    return residuals

    
class EMFitTools(Filters):
    
    # FIXME: all the work in this function except possibly the super
    # init depends only on the contents of emlinesfile.  We should
    # read this file ONCE and turn it into an object that encapsulates
    # line_table, line_map, and param_table, which are fixed once we
    # read the file.
    
    def __init__(self, uniqueid=None, fphoto=None, emlinesfile=None, stronglines=False):
        
        super(EMFitTools, self).__init__(fphoto=fphoto)
        
        self.uniqueid = uniqueid

        # info about tied doublet lines
        doublet_lines = {
            # indx           source         ratio name               
            'mgii_2796' : ( 'mgii_2803' , 'mgii_doublet_ratio' ) , 
            'oii_3726'  : ( 'oii_3729'  , 'oii_doublet_ratio'  ) ,
            'sii_6731'  : ( 'sii_6716'  , 'sii_doublet_ratio'  ) , 
        }

        self.line_table = read_emlines(emlinesfile=emlinesfile)

        # restrict to just strong lines and assign to patches; should probably
        # push this to the emlines.ecsv file
        if stronglines:
            I = np.where(np.isin(self.line_table['name'], 
                                 ['lyalpha', 'civ_1549', 'ciii_1908', 'mgii_2796', 'mgii_2803', 
                                  'oii_3726', 'oii_3729', 'hdelta', 'hdelta_broad', 'hgamma', 'hgamma_broad',
                                  'hbeta', 'hbeta_broad', 'oiii_4959', 'oiii_5007', 'nii_6548', 'halpha', 
                                  'halpha_broad', 'nii_6584', 'sii_6716', 'sii_6731']))[0]
            self.line_table = self.line_table[I]

            # assign continuum patches (should be in emlines.ecsv)
            continuum_patches = {
                'lyalpha': '0', 'civ_1549': '1', 'ciii_1908': '2', 
                'mgii_2796': '3', 'mgii_2803': '3', 
                'oii_3726': '4', 'oii_3729': '4',  
                'hdelta': '5', 'hdelta_broad': '5', 'hgamma': '5', 'hgamma_broad': '5', 
                'hbeta': '6', 'hbeta_broad': '6', 'oiii_4959': '6', 'oiii_5007': '6', 
                'nii_6548': '7', 'halpha': '7', 'halpha_broad': '7', 'nii_6584': '7', 'sii_6716': '7', 'sii_6731': '7'}
            self.line_table['continuum_patch'] = [continuum_patches[linename] for linename in self.line_table['name']]

        line_names = self.line_table['name'].value

        # mapping to enable fast lookup of line number by name
        self.line_map = { line_name : line_idx for line_idx, line_name in enumerate(line_names) }
        
        # mapping from target -> source for all tied doublets
        doublet_src = np.full(len(self.line_table), -1, dtype=np.int32)
        for doublet_tgt in doublet_lines:
            target_line = self.line_map[doublet_tgt]
            src_line    = self.line_map[ doublet_lines[doublet_tgt][0] ]
            doublet_src[target_line] = src_line
        self.line_table['doublet_src'] = doublet_src

        # create sparse doublet target -> src map; this can
        # be used to map target to src amplitudes as well,
        # because ParamType.AMPLITUDE == 0
        self.doublet_idx = np.where(doublet_src != -1)[0]
        self.doublet_src = doublet_src[self.doublet_idx]
        
        # build parameter names for every line in the line table,
        amp_names = [ f"{line_name}_amp" for line_name in line_names ]
        
        # use ratio names instead of target line names for tied doublet amps
        for doublet_target in doublet_lines:
            amp_names[ self.line_map[doublet_target] ] = \
                doublet_lines[ doublet_target ][1] # ratio name
        
        vshift_names = [ f"{line_name}_vshift" for line_name in line_names ]
        sigma_names  = [ f"{line_name}_sigma"  for line_name in line_names ]

        param_names = list(chain(amp_names, vshift_names, sigma_names))

        # compute type of each parameter in the parameter table
        nlines = len(self.line_table)
        param_types = np.empty(3*nlines, dtype=ParamType)
        for t in ParamType:
            param_types[t*nlines:(t+1)*nlines] = t
        
        self.param_table = Table({
            'name' : param_names,
            'type' : param_types,
            'line' : np.tile(np.arange(nlines, dtype=np.int32), 3), # param's line in line_table
        }, copy=False)
        
        # assign each line in line_table the indices of its 3 params in the name list
        param_idx = np.empty((nlines, 3), dtype=np.int32)
        c = np.arange(nlines, dtype=np.int32)
        param_idx[:,ParamType.AMPLITUDE] = c
        param_idx[:,ParamType.VSHIFT]    = c + nlines
        param_idx[:,ParamType.SIGMA]     = c + 2*nlines
        self.line_table['params'] = param_idx
        
        # needed by emlinemodel_bestfit()
        self.param_table['modelname'] = \
            np.array([ s.replace('_amp', '_modelamp').upper() for s in param_names ])
        

    #
    #  Record which lines are within the limits of the cameras
    # 
    def compute_inrange_lines(self, redshift, wavelims=(3000, 10000)):

        restwave = self.line_table['restwave'].value
        wavepad = 2.5 # Angstrom
        
        self.line_in_range = \
            (restwave > (wavelims[0] + wavepad)/(1 + redshift)) & \
            (restwave < (wavelims[1] - wavepad)/(1 + redshift))
        
    
    #
    # Build emission line model tables, with and without suppression of broad lines.
    # Establish fixed (i.e., forced to zero) params, tying relationships between
    # params, and doublet relationships for each model, as well as the relationship
    # between lines and their parameters, which we record in the line table.
    #
    # Parameter fixing needs to know which lines are within the observed wavelength
    # ranges of the cameras, so we first add this information to the line table.
    #
    def build_linemodels(self, strict_broadmodel=True):
        
        #
        # create_model()
        #
        # Given the tying info for the model and the list of in-range
        # lines, determine which parameters of the model are fixed and
        # free, and create the model's table.
        #
        # We fix all parameters for out-of-range lines to zero, unless
        # the parameter is tied to a parameter for an in-range line.
        # We also allow the caller to specify parameters that should
        # be fixed regardless.
        #
        def create_model(tying_info, forceFixed=[]):

            n_params = len(self.param_table)
            isfixed = np.full(n_params, False, dtype=bool)

            tiedtoparam, tiedfactor = tying_info
            tied_mask   = (tiedtoparam != -1)
            tied_source = tiedtoparam[tied_mask]
            
            # fix any parameters explicitly listed by user and 
            # propagate fixed status to their tied params
            if len(forceFixed) > 0:
                isfixed[forceFixed] = True
                isfixed[tied_mask] = isfixed[tied_source]
            
            # identify all params of out-of-range lines
            out_of_range_lines = ~self.line_in_range
            out_of_range_params = out_of_range_lines[self.param_table['line']]
            
            # for each param, count num of other params tied to it
            n_tied = np.bincount(tied_source, weights=np.ones_like(tied_source), minlength=n_params)
            
            # fix any parameters for an out-of-range line that are not the source of another
            # tied parameter
            isfixed[out_of_range_params & (n_tied == 0)] = True
            # for each param, count # of *fixed* params tied to it
            n_tied_fixed = np.bincount(tied_source, weights=isfixed[tied_mask], minlength=n_params)
            
            # Fix any parameter for an out-of-range line for which all its tied params
            # (if any) are fixed.
            isfixed[out_of_range_params & (n_tied == n_tied_fixed)] = True
            
            # finally, fix any doublet ratio whose source is fixed
            isfixed[self.doublet_idx] = isfixed[self.doublet_src]
            
            # delete tying relationships for fixed parameters
            tiedtoparam[isfixed] = -1

            istied = (tiedtoparam != -1)
            tiedfactor[~istied] = 0. # cosmetic cleanup

            # construct the final linemodel
            linemodel = Table({
                'free': ~isfixed & ~istied,
                'fixed': isfixed,
                'tiedtoparam': tiedtoparam.copy(), # we reuse these later
                'tiedfactor':  tiedfactor.copy(),
            }, copy=False)

            return linemodel

        
        # tie parameters of given line to source line.  We don't
        # tie the amplitude unless a tying factor is given for it
        def tie_line(tying_info, line_params, source_linename, amp_factor=None):

            tiedtoparam, tiedfactor = tying_info
            
            amp, vshift, sigma = line_params
            
            source_line = self.line_map[source_linename]
            src_amp, src_vshift, src_sigma = self.line_table['params'][source_line]
                        
            if amp_factor != None:
                tiedfactor[amp] = amp_factor
                tiedtoparam[amp] = src_amp
                    
            tiedfactor[vshift] = 1.0
            tiedtoparam[vshift] = src_vshift
                
            tiedfactor[sigma] = 1.0
            tiedtoparam[sigma] = src_sigma
        

        # Build the relationship of "tied" parameters. In the 'tied' array, the
        # non-zero value is the multiplicative factor by which the parameter
        # represented in the 'tiedtoparam' index should be multiplied.

        n_params = len(self.param_table)
        tying_info = (
            np.full(n_params, -1, np.int32), # source parameter for tying
            np.empty(n_params, np.float64),  # multiplier between source and tied value 
        )

        # Physical doublets and lines in the same ionization species should have
        # their velocity shifts and line-widths always tied. In addition, set fixed
        # doublet-ratios here. Note that these constraints must be set on *all*
        # lines, not just those in range.
        
        for line_name, line_isbalmer, line_isbroad, line_params in \
                self.line_table.iterrows('name', 'isbalmer', 'isbroad', 'params'):
            
            # broad He + Balmer
            if line_isbalmer and line_isbroad and line_name != 'halpha_broad':
                tie_line(tying_info, line_params, 'halpha_broad')
            # narrow He + Balmer
            elif line_isbalmer and not line_isbroad and line_name != 'halpha':
                tie_line(tying_info, line_params, 'halpha')
            else:
                match line_name:
                    case 'mgii_2796':
                        tie_line(tying_info, line_params, 'mgii_2803')
                    case 'oii_3726':
                        tie_line(tying_info, line_params, 'oii_3729')
                    case 'sii_6731':
                        tie_line(tying_info, line_params, 'sii_6716')
                    case 'nev_3346' | 'nev_3426': # should [NeIII] 3869 be tied to [NeV]???
                        tie_line(tying_info, line_params, 'neiii_3869')

                    case 'nii_5755' | 'oi_6300' | 'siii_6312':                        
                        # Tentative! Tie auroral lines to [OIII] 4363 but maybe we shouldn't tie [OI] 6300 here...
                        tie_line(tying_info, line_params, 'oiii_4363')
                    case 'oiii_4959':
                        """
                        [O3] (4-->2): airwave: 4958.9097 vacwave: 4960.2937 emissivity: 1.172e-21
                        [O3] (4-->3): airwave: 5006.8417 vacwave: 5008.2383 emissivity: 3.497e-21
                        """
                        tie_line(tying_info, line_params, 'oiii_5007', amp_factor = 1.0 / 2.9839)
                    case 'nii_6548':
                        """
                        [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.02198e-21
                        [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.94901e-21
                        """
                        tie_line(tying_info, line_params, 'nii_6584', amp_factor = 1.0 / 2.9421)
                    case 'oii_7330':
                        """
                        [O2] (5-->2): airwave: 7318.9185 vacwave: 7320.9350 emissivity: 8.18137e-24
                        [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 2.40519e-23
                        [O2] (5-->3): airwave: 7329.6613 vacwave: 7331.6807 emissivity: 1.35614e-23
                        [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.27488e-23
                        """
                        tie_line(tying_info, line_params, 'oii_7320', amp_factor = 1.0 / 1.2251)
                    case 'siii_9069':
                        tie_line(tying_info, line_params, 'siii_9532')
                    case 'siliii_1892':
                        # Tentative! Tie SiIII] 1892 to CIII] 1908 because they're so close in wavelength.
                        tie_line(tying_info, line_params, 'ciii_1908')
                    
            # Tie all the forbidden and narrow Balmer+helium lines *except
            # [OIII] 4959,5007* to [NII] 6584 when we have broad lines. The
            # [OIII] doublet frequently has an outflow component, so fit it
            # separately. See the discussion at
            # https://github.com/desihub/fastspecfit/issues/160
            if strict_broadmodel:
                if not line_isbroad and not line_name in { 'nii_6584', 'oiii_4959', 'oiii_5007' }:
                    tie_line(tying_info, line_params, 'nii_6584')
                    
                #if not line_isbroad and line_name != 'oiii_5007':
                #    tie_line(tying_info, line_params, 'oiii_5007')
                
                ## Tie all forbidden lines to [OIII] 5007; the narrow Balmer and
                ## helium lines are separately tied together.
                #if not line_isbroad and not line_isbalmer and line_name != 'oiii_5007'):
                #    tie_line(tying_info, line_params, 'oiii_5007')

        
        linemodel_broad = create_model(tying_info)

        
        # Model 2 - like Model 1, but additionally fix params of all
        # broad lines.  we inherit tying info from Model 1, which we
        # will modify below.
        
        forceFixed = []
        for line_name, line_isbalmer, line_isbroad, line_params in \
                self.line_table.iterrows('name', 'isbalmer', 'isbroad', 'params'):

            if line_name == 'halpha_broad':
                for p in line_params:  # all of amp, vshift, sigma
                    forceFixed.append(p) # fix all of these
            
            if line_isbalmer and line_isbroad and line_name != 'halpha_broad':
                tie_line(tying_info, line_params, 'halpha_broad', amp_factor = 1.0)
                
            if strict_broadmodel:
                # Tie the forbidden lines to [OIII] 5007.
                if not line_isbalmer and not line_isbroad and line_name != 'oiii_5007':
                    tie_line(tying_info, line_params, 'oiii_5007')
                
                # Tie narrow Balmer and helium lines together.
                if line_isbalmer and not line_isbroad:
                    if line_name == 'halpha':
                        tiedtoparam, _ = tying_info
                        
                        _, vshift, sigma = line_params
                        for p in (vshift, sigma):
                            # untie the params of this line
                            tiedtoparam[p] = -1    
                    else:
                        tie_line(tying_info, line_params, 'halpha')
        
        linemodel_nobroad = create_model(tying_info, forceFixed)
        
        return linemodel_broad, linemodel_nobroad
    
    
    # debugging function to print basic contents of a line model
    def summarize_linemodel(self, linemodel):

        """Simple function to summarize an input linemodel."""
        def _print(line_mask):
            for line in np.where(line_mask)[0]:
                line_name   = self.line_table['name'][line]
                line_params = self.line_table['params'][line]
                
                for param in line_params:
                    param_name    = linemodel['param_name'][param]
                    param_isfixed = linemodel['fixed'][param]
                    tiedtoparam   = linemodel['tiedtoparam'][param]
                    
                    if tiedtoparam == -1: # not tied
                        print(f'{param_name:25s} ', end='')
                        print('FIXED' if param_isfixed else 'free')
                    else:
                        source_name = linemodel['param_name'][tiedtoparam]
                        tiedfactor  = linemodel['tiedfactor'][param]
                        print(f'{param_name:25s} tied to {source_name:25s} '
                              f'with factor {tiedfactor:.4f}', end='')
                        print(' and FIXED' if param_isfixed else '')
                        
        line_isbroad  = self.line_table['isbroad']
        line_isbalmer = self.line_table['isbalmer']
        
        print('---------------------')
        print('UV/QSO (broad) lines:')
        print('---------------------')
        _print(line_isbroad & ~line_isbalmer)
        print()
        print('--------------------------')
        print('Broad Balmer+helium lines:')
        print('--------------------------')
        _print(line_isbroad & line_isbalmer)
        print()
        print('---------------------------')
        print('Narrow Balmer+helium lines:')
        print('---------------------------')
        _print(~line_isbroad & line_isbalmer)
        print()
        print('----------------')
        print('Forbidden lines:')
        print('----------------')
        _print(~line_isbroad & ~line_isbalmer)
        
    
    def initial_guesses_and_bounds(self, data, emlinewave, emlineflux, log):
        """For all lines in the wavelength range of the data, get a good initial guess
        on the amplitudes and line-widths. This step is critical for cases like,
        e.g., 39633354915582193 (tile 80613, petal 05), which has strong narrow
        lines.

        """
        initials = np.empty(len(self.param_table), dtype=np.float64)
        bounds   = np.empty((len(self.param_table), 2), dtype=np.float64)
        civars   = np.zeros(len(self.line_table), dtype=np.float64)
        
        #
        # a priori initial guesses and bounds
        #
        
        initvshift = 0.0
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
        
        for line_isbalmer, line_isbroad, line_params in \
                self.line_table.iterrows('isbalmer', 'isbroad', 'params'):
            
            amp, vshift, sigma = line_params

            # initial values and bounds for line's parameters
            initials[amp]    = initamp
            initials[vshift] = initvshift
            
            if line_isbroad:
                initials[sigma] = initsigma_broad
                
                if line_isbalmer: # broad He+Balmer lines
                    bounds[amp] = (minamp_balmer_broad, maxamp_balmer_broad) 
                    bounds[vshift] = (-vmaxshift_balmer_broad, +vmaxshift_balmer_broad)
                    bounds[sigma] = (minsigma_balmer_broad, maxsigma_balmer_broad)
                else: # broad UV/QSO lines (non-Balmer)
                    bounds[amp] = (minamp, maxamp)
                    bounds[vshift] = (-vmaxshift_broad, +vmaxshift_broad)
                    bounds[sigma] = (minsigma_broad, maxsigma_broad)
            else: # narrow He+Balmer lines, and forbidden lines
                initials[sigma] = initsigma_narrow
                
                bounds[amp] = (minamp, maxamp)
                bounds[vshift] = (-vmaxshift_narrow, +vmaxshift_narrow)
                bounds[sigma] = (minsigma_narrow, maxsigma_narrow)        
        
        #
        # Replace a priori initial values with those estimated by initial peak removal,
        # if available.
        #
        
        coadd_sigma = data['smoothsigma'] # robust estimate of the variance in the spectrum
        coadd_emlineflux = np.interp(data['coadd_wave'], emlinewave, emlineflux)
    
        for linename, linepix, contpix in zip(data['coadd_linename'],
                                              data['coadd_linepix'],
                                              data['coadd_contpix']):
            ## skip the physical doublets
            #if not hasattr(self.EMLineModel, f'{linename}_amp'):
            #    continue
            
            npix = len(linepix)
            if npix > 5:
                mnpx, mxpx = linepix[npix//2]-3, linepix[npix//2]+3
                mnpx = np.maximum(mnpx, 0)
                mxpx = np.minimum(mxpx, linepix[-1])
                amp = np.max(coadd_emlineflux[mnpx:mxpx])
            else:
                amp = quantile(coadd_emlineflux[linepix], 0.975)
            amp = np.abs(amp)
            
            noise = np.mean(coadd_sigma[linepix])
            if noise == 0.:
                civar = 0.
                errmsg = f'Noise estimate for line {linename} is zero!'
                log.warning(errmsg)
                #raise ValueError(errmsg)
            else:
                civar = 1. / noise**2
            
            # update the bounds on the line-amplitude
            #bounds = [-np.min(np.abs(coadd_emlineflux[linepix])), 3*np.max(coadd_emlineflux[linepix])]
            mx = 5*np.max(coadd_emlineflux[linepix])
            if mx < 0: # ???
                mx = 5*np.max(np.abs(coadd_emlineflux[linepix]))
            
            bds = np.array(( -1.5*np.min(np.abs(coadd_emlineflux[linepix])), mx ))
            
            # In extremely rare cases many of the pixels are zero, in which case
            # bounds[0] becomes zero, which is bad (e.g.,
            # iron/main/dark/27054/39627811564029314). Fix that here.
            if np.abs(bds[0]) == 0.0:
                N = coadd_emlineflux[linepix] != 0
                if np.sum(N) > 0:
                    bds[0] = -1.5*np.min(np.abs(coadd_emlineflux[linepix][N]))
                if np.abs(bds[0]) == 0.0:
                    bds[0] = -1e-3 # ??
            
            if (bds[0] > bds[1]) or (amp < bds[0]) or (amp > bds[1]):
                log.warning(f'Initial amplitude is outside its bound for line {linename}.')
                amp = np.diff(bounds)/2 + bounds[0]
                # Should never happen.
                if (bds[0] > bds[1]) or (amp < bds[0]) or (amp > bds[1]):
                    errmsg = f'Initial amplitude is outside its bound for line {linename}.'
                    self.log.critical(errmsg)
                    raise ValueError(errmsg)
            
            #
            # record our initial gueses and bounds for amplitude
            # and a guess at civar for the line
            
            line = self.line_map[linename]
            amp_idx = self.line_table['params'][line, ParamType.AMPLITUDE]
            initials[amp_idx] = amp
            bounds[amp_idx] = bds
            civars[line] = civar

        # Now update the linewidth but here we need to loop over *all* lines
        # (not just those in range). E.g., if H-alpha is out of range we need to
        # set its initial value correctly since other lines are tied to it
        # (e.g., main-bright-32406-39628257196245904).

        for isbroad, isbalmer, params in \
                self.line_table.iterrows('isbroad', 'isbalmer', 'params'):
            
            if isbroad:
                if isbalmer: # broad Balmer lines
                    if data['linesigma_balmer'] > data['linesigma_narrow']:
                        isigma = data['linesigma_balmer']
                    else:
                        isigma = data['linesigma_narrow']
                else: # broad UV/QSO lines
                    isigma = data['linesigma_uv']
            else:
                # prefer narrow over Balmer
                isigma = data['linesigma_narrow']

            # Make sure the initial guess for the narrow Balmer+helium
            # line is smaller than the guess for the broad-line model.
            if isbroad and isbalmer:
                isigma *= 1.1
            
            sigma_idx = params[ParamType.SIGMA]
            initials[sigma_idx] = isigma

        # Specialized parameters on the MgII, [OII], and [SII] doublet ratios. See
        # https://github.com/desihub/fastspecfit/issues/39. Be sure to set
        # self.doublet_names and also note that any change in the order of
        # these lines has to be handled in _emline_spectrum!
        doublet_bounds = {
            'mgii_doublet_ratio' : (0.0, 10.0), # MgII 2796/2803
            'oii_doublet_ratio'  : (0.0,  2.0), # [OII] 3726/3729 # (0.5, 1.5) # (0.66, 1.4)
            'sii_doublet_ratio'  : (0.0,  2.0), # [SII] 6731/6716 # (0.5, 1.5) # (0.67, 1.2)
        }

        param_names = self.param_table['name'].value

        # assign special bounds to each tied doublet parameter
        for iparam in self.doublet_idx:
            param_name = param_names[iparam]
            bounds[iparam] = doublet_bounds[param_name]
            initials[iparam] = 1.

        for iparam, param_name in enumerate(self.param_table['name'].value):
            
            # make sure each parameter lies within its bounds
            iv = initials[iparam]
            lb, ub = bounds[iparam]
            if iv < lb:
                errmsg = \
                    f'Initial parameter {param_name} is outside its bound, ' + \
                    f'{iv:.2f} < {lb:.2f}.'
                log.warning(errmsg)
                
                initials[iparam] = lb
            if iv > ub:
                errmsg = \
                    f'Initial parameter {param_name} is outside its bound, ' + \
                    f'{iv:.2f} > {ub:.2f}.'
                log.warning(errmsg)
                
                initials[iparam] = ub
                
        return initials, bounds, civars

        
    def optimize(self, linemodel,
                 initials, param_bounds,
                 obs_bin_centers,
                 obs_bin_fluxes,
                 obs_weights,
                 redshift,
                 resolution_matrices,
                 camerapix,
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

        line_wavelengths = self.line_table['restwave'].value

        isFree      = linemodel['free'].value
        tiedtoparam = linemodel['tiedtoparam'].value
        tiedfactor  = linemodel['tiedfactor'].value
                
        params_mapping = EMLine_ParamsMapping(len(linemodel), isFree,
                                              tiedtoparam, tiedfactor,
                                              self.doublet_idx, self.doublet_src)
        nFree = np.sum(isFree)
        
        log.debug(f"Optimizing {nFree} free parameters")
        
        if nFree == 0:
            # corner case where all lines are out of the wavelength range, which can
            # happen at high redshift and with the red camera masked, e.g.,
            # iron/main/dark/6642/39633239580608311).
            fit_info = {'nfev': 0, 'status': 0}
        else:
            initial_guesses = initials[isFree]
            bound_pairs = param_bounds[isFree, :]
            bounds = (bound_pairs[:,0], bound_pairs[:,1])
            
            obj = EMLine_Objective(obs_bin_centers,
                                   obs_bin_fluxes,
                                   obs_weights,
                                   redshift,
                                   line_wavelengths,
                                   tuple(resolution_matrices),
                                   camerapix, # for more efficient iteration
                                   params_mapping)
        
            try:
                fit_info = least_squares(obj.objective, initial_guesses, jac=obj.jacobian, args=(),
                                         max_nfev=5000, xtol=1e-10, ftol=1e-5, #x_scale='jac' gtol=1e-10,
                                         tr_solver='lsmr', tr_options={'maxiter': 1000, 'regularize': True},
                                         method='trf', bounds=bounds,) # verbose=2)
                free_params = fit_info.x
            except:
                if self.uniqueid:
                    errmsg = f'Problem in scipy.optimize.least_squares for {self.uniqueid}.'
                else:
                    errmsg = 'Problem in scipy.optimize.least_squares.'
                log.critical(errmsg)
                raise RuntimeError(errmsg)
            
            # Drop (zero out) any dubious free parameters.
            self._drop_params(linemodel, isFree, free_params, bounds, log)

            # translate free parame to full param array, but do NOT turn doublet
            # ratios into amplitudes yet, as out_linemodel needs them to be ratios
            parameters = params_mapping.mapFreeToFull(free_params, patchDoublets=False)
            
        linemodel['value'] = parameters.copy() # protect from changes below
        linemodel.meta['nfev'] = fit_info['nfev']
        linemodel.meta['status'] = fit_info['status']

        # convert doublet ratios to amplitudes
        parameters[self.doublet_idx] *= parameters[self.doublet_src]
        
        # calculate the observed maximum amplitude for each
        # fitted spectral line after convolution with the resolution
        # matrix.
        obsvalues = EMLine_find_peak_amplitudes(parameters,
                                                obs_bin_centers,
                                                redshift,
                                                line_wavelengths,
                                                resolution_matrices,
                                                camerapix)
        
        # add observed values as metadata, since they are only
        # relevant to amplitudes, not all parameters
        linemodel.meta['obsvalue'] = obsvalues
        
        return linemodel
    
    
    def optimize_continuum_patches(self, linemodel, 
                               initials, param_bounds,
                               continuum_patches, 
                               patch_initials, patch_param_bounds,
                               obs_bin_centers,
                               obs_bin_fluxes,
                               obs_weights,
                               redshift,
                               resolution_matrices,
                               camerapix,
                               log=None,
                               verbose=False,
                               debug=False):
        """Optimization routine which includes the underlying continuum.
    
        """
        from scipy.optimize import least_squares
        
        if log is None:
            from desiutil.log import get_logger, DEBUG
            if verbose:
                log = get_logger(DEBUG)
            else:
                log = get_logger()

        line_wavelengths = self.line_table['restwave'].value

        isFree      = linemodel['free'].value
        tiedtoparam = linemodel['tiedtoparam'].value
        tiedfactor  = linemodel['tiedfactor'].value
        isTied      = ~isFree * (tiedtoparam != -1)
        nFree = np.sum(isFree)

        # handle the continuum patch parameters
        nPatch = len(continuum_patches)
        isFree_patch = continuum_patches['free'].value
        nFree_patch = 2 * np.sum(isFree_patch) # slope + intercept

        #params_mapping = EMLine_ParamsMapping(len(linemodel), isFree,
        #                                      tiedtoparam, tiedfactor,
        #                                      self.doublet_idx, self.doublet_src)
        
        log.debug(f"Optimizing {nFree} free parameters")

        if nFree + nFree_patch == 0:
            # corner case where all lines are out of the wavelength range, which can
            # happen at high redshift and with the red camera masked, e.g.,
            # iron/main/dark/6642/39633239580608311).
            fit_info = {'nfev': 0, 'status': 0}
        else:
            allparameters  = np.hstack((np.zeros(len(linemodel)),
                                        patch_initials[:nPatch],  # slopes
                                        patch_initials[nPatch:])) # intercepts
            initial_guesses = np.hstack((initials[isFree], 
                                         patch_initials[:nPatch][isFree_patch],  # slopes
                                         patch_initials[nPatch:][isFree_patch])) # intercepts
            bound_pairs = np.vstack((param_bounds[isFree, :], 
                                     patch_param_bounds[:nPatch, :][isFree_patch, :],  # slopes
                                     patch_param_bounds[nPatch:, :][isFree_patch, :])) # intercepts
            bounds = (bound_pairs[:, 0], bound_pairs[:, 1])

            farg = [obs_bin_centers,
                    obs_bin_fluxes,
                    obs_weights,
                    redshift,
                    line_wavelengths,
                    [resolution_matrices],
                    camerapix, 
                    allparameters,
                    nFree,
                    isFree,
                    isTied,
                    tiedtoparam[isTied],
                    tiedfactor[isTied],
                    self.doublet_idx, 
                    self.doublet_src,
                    nFree_patch,
                    isFree_patch,
                    continuum_patches[isFree_patch],
                    nPatch,
                    False]

            try:
                fit_info = least_squares(emline_patch_objective, initial_guesses, args=farg,
                                         max_nfev=5000, xtol=1e-10, ftol=1e-5, #x_scale='jac' gtol=1e-10,
                                         tr_solver='lsmr', tr_options={'maxiter': 1000, 'regularize': True},
                                         method='trf', bounds=bounds, verbose=2)
                free_params = fit_info.x
            except:
                if self.uniqueid:
                    errmsg = f'Problem in scipy.optimize.least_squares for {self.uniqueid}.'
                else:
                    errmsg = 'Problem in scipy.optimize.least_squares.'
                log.critical(errmsg)
                raise RuntimeError(errmsg)


            farg[-1] = True
            lineparams, patchparams, patchmodel = emline_patch_objective(fit_info.x, *farg)
            
            import matplotlib.pyplot as plt
            plt.clf()
            plt.plot(obs_bin_centers, obs_bin_fluxes)
            plt.plot(obs_bin_centers, patchmodel, alpha=0.5, color='red')
            plt.savefig('junk.png')

            pdb.set_trace()
            
            # Drop (zero out) any dubious free parameters.
            self._drop_params(linemodel, isFree, free_params, bounds, log)

            # translate free parame to full param array, but do NOT turn doublet
            # ratios into amplitudes yet, as out_linemodel needs them to be ratios
            #parameters = params_mapping.mapFreeToFull(free_params, patchDoublets=False)
            
        linemodel['value'] = parameters.copy() # protect from changes below
        linemodel.meta['nfev'] = fit_info['nfev']
        linemodel.meta['status'] = fit_info['status']

        # convert doublet ratios to amplitudes
        parameters[self.doublet_idx] *= parameters[self.doublet_src]
        
        # calculate the observed maximum amplitude for each
        # fitted spectral line after convolution with the resolution
        # matrix.
        obsvalues = EMLine_find_peak_amplitudes(parameters,
                                                obs_bin_centers,
                                                redshift,
                                                line_wavelengths,
                                                resolution_matrices,
                                                camerapix)
        
        # add observed values as metadata, since they are only
        # relevant to amplitudes, not all parameters
        linemodel.meta['obsvalue'] = obsvalues
        
        return linemodel
    
    
    def _drop_params(self, linemodel, isFree, free_params, free_bounds, log):
        """Drop dubious free parameters after fitting.
    
        """
        
        param_type = self.param_table['type'][isFree]
        
        # Drop any negative amplitude or nonposiive sigma.  These
        # happen because we do not enforce these conditions in our
        # bounds, for allegedly improved convergence (?).
        drop1 = \
            (((param_type == ParamType.AMPLITUDE) & (free_params <  0.)) | \
             ((param_type == ParamType.SIGMA)     & (free_params <= 0.)))
        
        # Require equality, not np.isclose, because the optimization can be very
        # small (<1e-6) but still significant, especially for the doublet
        # ratios. If linesigma is dropped this way, make sure the corresponding
        # line-amplitude is dropped, too (see MgII 2796 on
        # sv1-bright-17680-39627622543528153).

        drop2 = np.full(len(linemodel), False, bool)
        Ifree = np.where(isFree)[0]
        
        param_line = self.param_table['line']
        line_params = self.line_table['params']
        
        # if any amplitude is zero, drop the corresponding sigma and vshift
        isZeroAmp = \
            (param_type == ParamType.AMPLITUDE) & (free_params == 0.)
        for param in Ifree[isZeroAmp]:
            line = param_line[param]
            vshiftdropped = line_params[line, ParamType.VSHIFT]
            drop2[vshiftdropped] = True
            sigmadropped  = line_params[line, ParamType.SIGMA]
            drop2[sigmadropped] = True
        
        tiedtoparam = linemodel['tiedtoparam']
        
        # drop amplitude for any line tied to a line with a dropped sigma/vshift
        isDroppedSV = (((param_type == ParamType.SIGMA) |
                        (param_type == ParamType.VSHIFT)) &
                       drop2[isFree])
        for param in Ifree[isDroppedSV]:
            line = param_line[tiedtoparam == param]
            ampdropped = line_params[line, ParamType.AMPLITUDE]
            drop2[ampdropped] = True

        drop2 = drop2[isFree]
        
        # This should not happen if our optimizer enforces its bounds
        drop3 = ((free_params < free_bounds[0]) | (free_params > free_bounds[1]))
        
        drop = np.logical_or.reduce((drop1, drop2, drop3))
        
        log.debug(f'Dropping {np.sum(drop1)} impossible line amplitudes/widths.')
        log.debug(f'Dropping {np.sum(drop2)} parameters of or tied to zero-amplitude lines.')
        log.debug(f'Dropping {np.sum(drop3)} parameters which are out-of-bounds.')

        ntotal = np.sum(drop)
        if ntotal > 0:
            log.debug(f'  Dropping {ntotal} unique parameters.')
        
        free_params[drop] = 0.
        
        
    @staticmethod
    def chi2(linemodel, emlinewave, emlineflux, emlineivar, emlinemodel,
             continuum_model=None, return_dof=False):
        """Compute the reduced chi^2."""

        nfree = np.sum(linemodel['free'])
        dof   = np.sum(emlineivar > 0) - nfree
        
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

        linewaves = self.line_table['restwave'].value
        
        parameters = linemodel['value'].copy()

        # convert doublet ratios to amplitudes
        parameters[self.doublet_idx] *= parameters[self.doublet_src]
        
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line
        
        emlinemodel = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas,
                                         linewaves, emlinewave, resolution_matrix, camerapix)
        
        return emlinemodel
    

    def emlinemodel_bestfit(self, fastspecfit_table, specwave, specres, camerapix,
                            redshift=None, snrcut=None):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (used for QA and elsewhere).

        """
        
        linewaves = self.line_table['restwave'].value
        
        parameters = np.array([ fastspecfit_table[param] for param in self.param_table['modelname'] ])
        
        # convert doublet ratios to amplitudes
        parameters[self.doublet_idx] *= parameters[self.doublet_src]
        
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line    
        
        if snrcut is not None:
            line_names = self.line_table['name'].value
            lineamps_ivar = [fastspecfit_table[line_name.upper()+'_AMP_IVAR'] for line_name in line_names]
            lineamps[lineamps * np.sqrt(lineamps_ivar) < snrcut] = 0.

        if redshift is None:
            redshift = fastspecfit_table['Z']

        model_fluxes = EMLine_build_model(redshift, lineamps, linevshifts, linesigmas, linewaves,
                                          np.hstack(specwave), specres, camerapix)
        
        return model_fluxes

        
    def populate_emtable(self, result, finalfit, finalmodel, emlinewave, emlineflux,
                         emlineivar, oemlineivar, specflux_nolines, redshift,
                         resolution_matrices, camerapix, log, nminpix=7, nsigma=3.):
        """Populate the output table with the emission-line measurements.

        """
        
        # return range (lo, hi) s.t. all pixels of A in range [v_lo,
        # v_hi] lie in half-open interval [lo, hi)
        def get_boundaries(A, v_lo, v_hi):
            return np.searchsorted(A, (v_lo, v_hi), side='right')
        
        from math import erf
        from fastspecfit.util import centers2edges, sigmaclip, quantile
        
        line_wavelengths = self.line_table['restwave'].value
        
        parameters = finalfit['value'].value.copy()
        
        # convert doublet ratios to amplitudes
        parameters[self.doublet_idx] *= parameters[self.doublet_src]
        
        line_fluxes = EMLine_MultiLines(parameters,
                                        emlinewave,
                                        redshift,
                                        line_wavelengths,
                                        resolution_matrices,
                                        camerapix)

        values = finalfit['value'].value
        obsvalues = finalfit.meta['obsvalue']
        
        gausscorr = erf(nsigma / np.sqrt(2))      # correct for the flux outside of +/-nsigma
        dpixwave = np.median(np.diff(emlinewave)) # median pixel size [Angstrom]
        
        # Where the cameras overlap, we have to account for the
        # variable pixel size by sorting in wavelength.
        Wsrt = np.argsort(emlinewave)
        
        emlinewave_s  = emlinewave[Wsrt]
        emlineflux_s  = emlineflux[Wsrt]
        emlineivar_s  = emlineivar[Wsrt]
        oemlineivar_s = oemlineivar[Wsrt]
        finalmodel_s  = finalmodel[Wsrt]
        specflux_nolines_s = specflux_nolines[Wsrt]
        
        dwaves = np.diff(centers2edges(emlinewave_s))
        
        # get continuum fluxes, EWs, and upper limits
        narrow_stats, broad_stats, uv_stats = [], [], []
        for iline, (name, restwave, isbroad, isbalmer) in \
                enumerate(self.line_table.iterrows('name', 'restwave', 'isbroad', 'isbalmer')):
            
            linename = name.upper()

            line_amp, line_vshift, line_sigma = \
                self.line_table['params'][iline]
            
            # zero out out-of-range lines
            if not self.line_in_range[iline]:
                obsvalues[line_amp] = 0.
                values[line_amp]    = 0.
                values[line_vshift] = 0.
                values[line_sigma]  = 0.
                continue

            linez = redshift + values[line_vshift] / C_LIGHT
            linezwave = restwave * (1 + linez)
            linesigma = values[line_sigma] # [km/s]
            
            # if the line was dropped, use a default sigma value
            if linesigma == 0:
                limitsigma_narrow = 75.0
                limitsigma_broad = 1200.0
                
                if isbroad and isbalmer:
                    linesigma = limitsigma_narrow
                else:
                    linesigma = limitsigma_broad
                    
            linesigma_ang = linesigma * linezwave / C_LIGHT    # [observed-frame Angstrom]
            
            # require at least 2 pixels
            if linesigma_ang < 2 * dpixwave:
                linesigma_ang_window = 2 * dpixwave
                _gausscorr = 1.
            else:
                linesigma_ang_window = linesigma_ang
                _gausscorr = gausscorr
            
            # emlinewave is NOT sorted because it mixes cameras, so we must check all elts
            # to find pixels near line
            line_s, line_e = get_boundaries(emlinewave_s,
                                            linezwave - nsigma*linesigma_ang_window,
                                            linezwave + nsigma*linesigma_ang_window)
            
            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next line.
            if np.sum(oemlineivar_s[line_s:line_e] == 0.) > 0.3 * (line_e - line_s): # use original ivar
                obsvalues[line_amp] = 0.
                values[line_amp]    = 0.
                values[line_vshift] = 0.
                values[line_sigma]  = 0.
                continue
            
            # number of pixels, chi2, and boxcar integration
            patchindx = line_s + np.where( emlineivar_s[line_s:line_e] > 0. )[0]
            npix = len(patchindx)
            result[f'{linename}_NPIX'] = npix
            
            if npix >= nminpix: # magic number: required at least XX unmasked pixels centered on the line
                
                emlineflux_patch = emlineflux_s[patchindx]
                emlineivar_patch = emlineivar_s[patchindx]
                dwaves_patch     = dwaves[patchindx]
                
                if np.any(emlineivar_patch == 0.):
                    errmsg = 'Ivar should never be zero within an emission line!'
                    log.critical(errmsg)
                    raise ValueError(errmsg)
                
                # boxcar integration of the flux
                boxflux = np.sum(emlineflux_patch * dwaves_patch)
                boxflux_ivar = 1 / np.sum(dwaves_patch**2 / emlineivar_patch)
                result[f'{linename}_BOXFLUX'] = boxflux # * u.erg/(u.second*u.cm**2)
                result[f'{linename}_BOXFLUX_IVAR'] = boxflux_ivar # * u.second**2*u.cm**4/u.erg**2
            
                # Get the uncertainty in the line-amplitude based on the scatter
                # in the pixel values from the emission-line subtracted
                # spectrum.
                n_lo, n_hi = quantile(specflux_nolines_s[patchindx], (0.25, 0.75))
                amp_sigma = (n_hi - n_lo) / 1.349 # robust sigma
                
                amp_ivar = 1/amp_sigma**2 if amp_sigma > 0. else 0.
                if amp_ivar > 0.:
                    result[f'{linename}_AMP_IVAR'] = amp_ivar # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
                
                # require amp > 0 (line not dropped) to compute the flux and chi2
                flux, flux_ivar = 0., 0.
                if values[line_amp] > 0.:
                    finalmodel_patch = finalmodel_s[patchindx]
                    chi2 = np.sum(emlineivar_patch * (emlineflux_patch - finalmodel_patch)**2)
                    result[f'{linename}_CHI2'] = chi2
                    
                    (s, e), flux = line_fluxes.getLine(iline)
                    
                    # we could do this sparsely, but it's slower than allocating
                    # and permuting a big, mostly empty array
                    lineprofile = np.zeros_like(emlinewave)
                    lineprofile[s:e] = flux
                    lineprofile_patch = lineprofile[Wsrt][patchindx]

                    patch_sum = np.sum(lineprofile_patch)
                    if patch_sum == 0. or np.any(lineprofile_patch < 0.):
                        errmsg = 'Line-profile should never be zero or negative!'
                        log.critical(errmsg)
                        raise ValueError(errmsg)
                    
                    # matched-filter (maximum-likelihood) Gaussian flux
                    pro_j = lineprofile_patch / patch_sum
                    I = pro_j > 0. # very narrow lines can have a profile that goes to zero

                    r = pro_j[I] / dwaves_patch[I]
                    weight_j = r * emlineivar_patch[I]
                    flux_ivar = np.sum(r * weight_j)                    
                    flux = np.sum(weight_j * lineprofile_patch[I]) / flux_ivar

                    # correction for missing flux
                    flux      /= _gausscorr
                    flux_ivar *= _gausscorr**2
                    result[f'{linename}_FLUX'] = flux
                    result[f'{linename}_FLUX_IVAR'] = flux_ivar # * u.second**2*u.cm**4/u.erg**2
                    
                    # keep track of sigma and z but only using XX-sigma lines
                    linesnr = obsvalues[line_amp] * np.sqrt(amp_ivar)
                    if linesnr > 1.5:
                        if isbroad: # includes UV and broad Balmer lines
                            if isbalmer:
                                stats = broad_stats
                            else:
                                stats = uv_stats
                        else:
                            stats = narrow_stats
                        stats.append((linesigma, linez))

            
            # next, get the continuum, the inverse variance in the line-amplitude, and the EW
            slo, elo = get_boundaries(emlinewave_s,
                                      linezwave - 10 * linesigma_ang_window,
                                      linezwave -  3 * linesigma_ang_window)
            shi, ehi = get_boundaries(emlinewave_s,
                                      linezwave +  3 * linesigma_ang_window,
                                      linezwave + 10 * linesigma_ang_window)
            
            borderindx = np.hstack((slo + np.where(oemlineivar_s[slo:elo] > 0.)[0],
                                    shi + np.where(oemlineivar_s[shi:ehi] > 0.)[0]))
            
            cmed, civar = 0., 0.
            if len(borderindx) >= nminpix: # require at least XX pixels to get the continuum level
                clipflux, _, _ = sigmaclip(specflux_nolines_s[borderindx], low=3, high=3)
                
                if len(clipflux) > 0:
                    clo, cmed, chi = quantile(clipflux, (0.25, 0.50, 0.75))
                    csig = (chi - clo) / 1.349  # robust sigma
                    civar = (np.sqrt(len(borderindx)) / csig)**2 if csig > 0. else 0.
                    
                result[f'{linename}_CONT'] = cmed # * u.erg/(u.second*u.cm**2*u.Angstrom)
                result[f'{linename}_CONT_IVAR'] = civar # * u.second**2*u.cm**4*u.Angstrom**2/u.erg**2
            
            ew, ewivar = 0., 0.
            if cmed != 0. and civar != 0.:
                if flux > 0. and flux_ivar > 0.:
                    # add the uncertainties in flux and the continuum in quadrature
                    ew = flux / cmed / (1 + redshift) # rest frame [A]
                    ewivar = (1+redshift)**2 / (1 / (cmed**2 * flux_ivar) + flux**2 / (cmed**4 * civar))
                
                # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                fluxlimit = np.sqrt(2 * np.pi) * linesigma_ang / np.sqrt(civar) # * u.erg/(u.second*u.cm**2)
                ewlimit = fluxlimit * cmed / (1+redshift)
                
                result[f'{linename}_EW'] = ew
                result[f'{linename}_EW_IVAR'] = ewivar
                result[f'{linename}_FLUX_LIMIT'] = fluxlimit
                result[f'{linename}_EW_LIMIT'] = ewlimit
                
            ########################################################
            
            if 'debug' in log.name:
                for col in ('VSHIFT', 'SIGMA', 'MODELAMP', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    val = result[f'{linename}_{col}']
                    log.debug(f'{linename} {col}: {val:.4f}')
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'BOXFLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT', 'EW_LIMIT'):
                    val = result[f'{linename}_{col}']
                    log.debug(f'{linename} {col}: {val:.4f}')
                print()
        
        # get the per-group average emission-line redshifts and velocity widths
        for stats, groupname in zip((narrow_stats, broad_stats, uv_stats),
                                    ('NARROW', 'BROAD', 'UV')):
            if len(stats) > 0:
                stats = np.array(stats)
                sigmas    = stats[:,0]
                redshifts = stats[:,1]

                stat_sigma    = np.median(sigmas)  # * u.kilometer / u.second
                stat_sigmarms = np.std(sigmas)

                stat_z    = np.median(redshifts)
                stat_zrms = np.std(redshifts)

                if 'debug' in log.name:
                    log.debug(f'{groupname}_SIGMA: {stat_sigma:.3f}+/-{stat_sigmarms:.3f}')
                    log.debug(f'{groupname}_Z:     {stat_z:.9f}+/-{stat_zrms:.9f}')
                    
                result[f'{groupname}_SIGMA']    = stat_sigma
                result[f'{groupname}_SIGMARMS'] = stat_sigmarms

                result[f'{groupname}_Z']        = stat_z 
                result[f'{groupname}_ZRMS']     = stat_zrms                    
            else:
                result[f'{groupname}_Z'] = redshift
    
                    
                
        # write values of final parameters (after any changes above) to result
        param_names      = self.param_table['name'].value
        param_modelnames = self.param_table['modelname'].value
        param_types      = self.param_table['type'].value
        param_lines      = self.param_table['line'].value
        line_doublet_src = self.line_table['doublet_src'].value
        
        # create result entries for every parameter with its fitted value
        # we need both model amplitude and computed amplitude from
        # peak-finding.
        for iparam in range(len(finalfit)):
            pmodelname = param_modelnames[iparam]
            val = values[iparam]
            
            result[pmodelname] = val

            # observed amplitudes
            if param_types[iparam] == ParamType.AMPLITUDE:
                if line_doublet_src[iparam] == -1:
                    # not a doublet ratio
                    result[param_names[iparam].upper()] = obsvalues[iparam]
                else:
                    # line name of doublet target
                    orig_line = self.line_table['name'][param_lines[iparam]].upper()
                    isrc = line_doublet_src[iparam] # valid for amplitude params
                    
                    result[orig_line + '_MODELAMP'] = val * values[isrc]
                    result[orig_line + '_AMP'     ] = val * obsvalues[isrc] 
                
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
                log.debug(f'{col}: {result[col]:.4f}')
            #log.debug(' ')
            print()
        
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

        synthmag = data['synthphot']['nanomaggies'].value
        model_synthmag = model_synthphot['nanomaggies'].value
        for iband, band in enumerate(self.synth_bands):
            bname =  band.upper()
            result[f'FLUX_SYNTH_{bname}'] = synthmag[iband] # * 'nanomaggies'
            result[f'FLUX_SYNTH_SPECMODEL_{bname}'] = model_synthmag[iband] # * 'nanomaggies'

        
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
    from astropy.table import Column
    from fastspecfit.util import ivar2var

    tall = time.time()

    if log is None:
        from desiutil.log import get_logger, DEBUG
        if verbose:
            log = get_logger(DEBUG)
        else:
            log = get_logger()

            
    minsigma_balmer_broad = 250. # minimum broad-line sigma [km/s]
    
    EMFit = EMFitTools(emlinesfile=emlinesfile, fphoto=fphoto, uniqueid=data['uniqueid'])
    
    redshift    = data['zredrock']
    camerapix    = data['camerapix']
    resolution_matrix = data['res_fast']

    # Combine pixels across all cameras
    emlinewave  = np.hstack(data['wave'])
    oemlineivar = np.hstack(data['ivar'])
    specflux    = np.hstack(data['flux'])
    
    continuummodelflux       = np.hstack(continuummodel)
    smoothcontinuummodelflux = np.hstack(smooth_continuum)
    emlineflux = specflux - continuummodelflux - smoothcontinuummodelflux

    emlineivar = np.copy(oemlineivar)
    _, emlinegood = ivar2var(emlineivar, clip=1e-3)
    emlinebad = ~emlinegood

    # This is a (dangerous???) hack.
    if np.any(emlinebad):
        emlineivar[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineivar[emlinegood])
        emlineflux[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineflux[emlinegood]) # ???
    
    weights = np.sqrt(emlineivar)

    # determine which lines are in range of the camera
    EMFit.compute_inrange_lines(redshift,
                                wavelims=(np.min(emlinewave)+5,
                                          np.max(emlinewave)-5))
    
    # Build all the emission-line models for this object.
    linemodel_broad, linemodel_nobroad = \
        EMFit.build_linemodels(strict_broadmodel=True)
    
    # Get initial guesses on the parameters and populate the two "initial"
    # linemodels; the "final" linemodels will be initialized with the
    # best-fitting parameters from the initial round of fitting.
    initial_guesses, param_bounds, civars = \
        EMFit.initial_guesses_and_bounds(data, emlinewave, emlineflux, log)
    
    # fit spectrum *without* any broad lines
    t0 = time.time()
    fit_nobroad = EMFit.optimize(linemodel_nobroad, initial_guesses, param_bounds,
                                 emlinewave, emlineflux, weights, redshift,
                                 resolution_matrix, camerapix, log=log, debug=False)
    model_nobroad = EMFit.bestfit(fit_nobroad, redshift, emlinewave, resolution_matrix, camerapix)
    chi2_nobroad, ndof_nobroad, nfree_nobroad = EMFit.chi2(fit_nobroad, emlinewave, emlineflux, emlineivar, model_nobroad, return_dof=True)
    log.info('Line-fitting {} with no broad lines and {} free parameters took {:.4f} seconds [niter={}, rchi2={:.4f}].'.format(
        data['uniqueid'], nfree_nobroad, time.time()-t0, fit_nobroad.meta['nfev'], chi2_nobroad))
    
    # Now try adding broad Balmer and helium lines and see if we improve the
    # chi2.
    if broadlinefit:

        pixoffset = np.cumsum(data['npixpercamera']) - data['npixpercamera']
        
        # Bbroad == non-helium broad Balmer lines
        line_isBbroad = \
            np.logical_and.reduce((EMFit.line_table['isbalmer'],
                                   EMFit.line_table['isbroad'],
                                   ~EMFit.line_table['ishelium']))
        
        # Gather the pixels around the broad Balmer lines, and count
        # the number of free parameters associated with these lines
        # in each of the two models

        line_params = EMFit.line_table['params'].value
        
        balmer_pix = []
        balmer_nfree_nobroad = 0
        balmer_nfree_broad   = 0
        
        for icam in range(len(data['cameras'])):
            for line_name, pix in zip(data['linename'][icam], data['linepix'][icam]):
                line = EMFit.line_map[line_name]
                if (line_isBbroad[line]):
                    params = line_params[line]
                    balmer_pix.append(pix + pixoffset[icam])
                    balmer_nfree_nobroad += np.sum(linemodel_nobroad['free'][params])
                    balmer_nfree_broad   += np.sum(linemodel_broad['free'][params])
        
        if len(balmer_pix) > 0:

            # fit spectrum *including* any broad lines
            t0 = time.time()
            fit_broad = EMFit.optimize(linemodel_broad, initial_guesses, param_bounds,
                                       emlinewave, emlineflux, weights, 
                                       redshift, resolution_matrix, camerapix, log=log,
                                       debug=False)
            model_broad = EMFit.bestfit(fit_broad, redshift, emlinewave, resolution_matrix, camerapix)
            chi2_broad, ndof_broad, nfree_broad = EMFit.chi2(fit_broad, emlinewave, emlineflux, emlineivar, model_broad, return_dof=True)
            log.info('Line-fitting {} with broad lines and {} free parameters took {:.4f} seconds [niter={}, rchi2={:.4f}].'.format(
                data['uniqueid'], nfree_broad, time.time()-t0, fit_broad.meta['nfev'], chi2_broad))
            
            # compute delta-chi2 around just the broad, non-helium Balmer lines
            balmer_pix = np.hstack(balmer_pix)
            
            bivar = emlineivar[balmer_pix]
            bflux = emlineflux[balmer_pix]
            nbpix = np.sum(bivar > 0)
            balmer_ndof_broad   = nbpix - balmer_nfree_broad
            balmer_ndof_nobroad = nbpix - balmer_nfree_nobroad
            
            linechi2_balmer_broad   = np.sum(bivar * (bflux - model_broad[balmer_pix])**2)
            linechi2_balmer_nobroad = np.sum(bivar * (bflux - model_nobroad[balmer_pix])**2)
            delta_linechi2_balmer   = linechi2_balmer_nobroad - linechi2_balmer_broad
            delta_linendof_balmer   = balmer_ndof_nobroad - balmer_ndof_broad
            
            # Choose broad-line model only if:
            # --delta-chi2 > delta-ndof
            # --broad_sigma < narrow_sigma
            # --broad_sigma < 250
            
            dchi2test = (delta_linechi2_balmer > delta_linendof_balmer)
            
            Hanarrow_idx = line_params[EMFit.line_map['halpha'],       ParamType.SIGMA]
            Hanarrow = fit_broad['value'][Hanarrow_idx]
            
            Habroad_idx  = line_params[EMFit.line_map['halpha_broad'], ParamType.SIGMA]
            Habroad  = fit_broad['value'][Habroad_idx]
            
            # Bbroad_amps enumerates param indices of all non-fixed amplitudes
            # of lines that are broad, Balmer, and not Helium
            Bbroad_all_amps   = line_params[line_isBbroad, ParamType.AMPLITUDE]
            Bbroad_all_stds   = np.sqrt(civars[line_isBbroad])
            Bbroad_all_lnames = EMFit.line_table['name'][line_isBbroad]
            
            Bbroad_nonFixed  = ~fit_broad['fixed'][Bbroad_all_amps]
            Bbroad_amps   = Bbroad_all_amps[Bbroad_nonFixed]
            Bbroad_stds   = Bbroad_all_stds[Bbroad_nonFixed]
            Bbroad_lnames = Bbroad_all_lnames[Bbroad_nonFixed]
            
            broadsnr = fit_broad.meta['obsvalue'][Bbroad_amps] * Bbroad_stds
            
            sigtest1 = Habroad > minsigma_balmer_broad
            sigtest2 = Habroad > Hanarrow
            if len(broadsnr) == 0:
                broadsnrtest = False
                _broadsnr = '0.'
            elif len(broadsnr) == 1:
                broadsnrtest = (broadsnr[-1] > minsnr_balmer_broad)
                _broadsnr = f'S/N ({Bbroad_lnames[-1]}) = {broadsnr[-1]:.1f}'
            else:
                broadsnrtest = np.any(broadsnr[-2:] > minsnr_balmer_broad)
                _broadsnr = \
                    f'S/N ({Bbroad_lnames[-2]}) = {broadsnr[-2]:.1f}, ' \
                    f'S/N ({Bbroad_lnames[-1]}) = {broadsnr[-1]:.1f}'
                
            if dchi2test and sigtest1 and sigtest2 and broadsnrtest:
                log.info('Adopting broad-line model:')
                log.info(f'  delta-chi2={delta_linechi2_balmer:.1f} > delta-ndof={delta_linendof_balmer:.0f}')
                log.info(f'  sigma_broad={Habroad:.1f} km/s, sigma_narrow={Hanarrow:.1f} km/s')
                if _broadsnr:
                    log.info(f'  {_broadsnr} > {minsnr_balmer_broad:.0f}')
                    
                finalfit, finalmodel, finalchi2 = fit_broad, model_broad, chi2_broad
            else:
                if dchi2test == False:
                    log.info(f'Dropping broad-line model: delta-chi2={delta_linechi2_balmer:.1f} < delta-ndof={delta_linendof_balmer:.0f}')
                elif sigtest1 == False:
                    log.info(f'Dropping broad-line model: Halpha_broad_sigma {Habroad:.1f} km/s < {minsigma_balmer_broad:.0f} km/s '
                             f'(delta-chi2={delta_linechi2_balmer:.1f}, delta-ndof={delta_linendof_balmer:.0f}).')
                elif sigtest2 == False:
                    log.info(f'Dropping broad-line model: Halpha_broad_sigma {Habroad:.1f} km/s < Halpha_narrow_sigma {Hanarrow:.1f} km/s '
                             f'(delta-chi2={delta_linechi2_balmer:.1f}, delta-ndof={delta_linendof_balmer:.0f}).')
                elif broadsnrtest == False:
                    log.info(f'Dropping broad-line model: {_broadsnr} < {minsnr_balmer_broad:.0f}')

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
    
    # Build the model spectrum
    emmodel = EMFit.emlinemodel_bestfit(result,
                                        data['wave'],
                                        data['res_fast'],
                                        camerapix,
                                        redshift=redshift)
    
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
    coadd_waves = data['coadd_wave']
    minwave = np.min(coadd_waves)
    maxwave = np.max(coadd_waves)
    dwave = coadd_waves[1] - coadd_waves[0]
    
    minwave = np.floor(minwave * 1000) / 1000
    maxwave = np.floor(maxwave * 1000) / 1000
    dwave   = np.round(dwave, decimals=3)
    
    npix = int(np.round((maxwave-minwave)/dwave)) + 1
    modelwave = minwave + dwave * np.arange(npix)

    wavesrt = np.argsort(emlinewave)
    sorted_waves = emlinewave[wavesrt]
    modelcontinuum       = np.interp(modelwave, sorted_waves, continuummodelflux[wavesrt])
    modelsmoothcontinuum = np.interp(modelwave, sorted_waves, smoothcontinuummodelflux[wavesrt])
    modelemspectrum      = np.interp(modelwave, sorted_waves, emmodel[wavesrt])
    
    modelspectra = Table(
        # ensure that these columns will stack as rows when
        # we vstack the Tables for different spectra, rather
        # than being concatenated into one long row.
        data=(
            Column(name='CONTINUUM', dtype='f4',
                   data=modelcontinuum.reshape(1, npix)),
            Column(name='SMOOTHCONTINUUM', dtype='f4',
                   data=modelsmoothcontinuum.reshape(1, npix)),
            Column(name='EMLINEMODEL', dtype='f4',
                   data=modelemspectrum.reshape(1, npix))
        ),
        # all these header cards need to be 2-element tuples (value, comment),
        # otherwise io.write_fastspecfit will crash
        meta = {
            'NAXIS1':   (npix, 'number of pixels'),
            'NAXIS2':   (npix, 'number of models'),
            'NAXIS3':   (npix, 'number of objects'),
            'BUNIT':    ('10**-17 erg/(s cm2 Angstrom)', 'flux unit'),
            'CUNIT1':   ('Angstrom', 'wavelength unit'),
            'CTYPE1':   ('WAVE', 'type of axis'),
            'CRVAL1':   (minwave, 'wavelength of pixel CRPIX1 (Angstrom)'),
            'CRPIX1':   (0, '0-indexed pixel number corresponding to CRVAL1'),
            'CDELT1':   (dwave, 'pixel size (Angstrom)'),
            'DC-FLAG':  (0, '0 = linear wavelength vector'),
            'AIRORVAC': ('vac', 'wavelengths in vacuum (vac)')
        },
        copy=False
    )
    
    # Finally, optionally synthesize photometry (excluding the
    # smoothcontinuum!) and measure Dn(4000) from the line-free spectrum.
    if synthphot:
        modelflux = modelcontinuum + modelemspectrum
        EMFit.synthphot_spectrum(data, result, modelwave, modelflux)

    # measure DN(4000) without the emission lines
    if result['DN4000_IVAR'] > 0:
        fluxnolines = data['coadd_flux'] - modelemspectrum
        
        dn4000_nolines, _ = EMFit.get_dn4000(modelwave, fluxnolines, redshift=redshift, log=log, rest=False)
        log.info(f'Dn(4000)={dn4000_nolines:.3f} in the emission-line subtracted spectrum.')
        result['DN4000'] = dn4000_nolines

        # Simple QA of the Dn(4000) estimate.
        """
        import matplotlib.pyplot as plt
        
        dn4000, dn4000_obs, dn4000_model, dn4000_ivar = result['DN4000'], result['DN4000_OBS'], result['DN4000_MODEL'], result['DN4000_IVAR']
        print(dn4000, dn4000_obs, dn4000_model, 1/np.sqrt(dn4000_ivar))
        
        restwave = modelwave / (1 + redshift) # [Angstrom]
        flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) # [erg/s/cm2/A-->erg/s/cm2/Hz, rest]
        fnu_obs = data['coadd_flux'] * flam2fnu # [erg/s/cm2/Hz]
        fnu = fluxnolines * flam2fnu # [erg/s/cm2/Hz]

        fnu_model = modelcontinuum * flam2fnu
        fnu_fullmodel = modelflux * flam2fnu
            
        fnu_ivar = data['coadd_ivar'] / flam2fnu**2            
        fnu_sigma, fnu_mask = ivar2var(fnu_ivar, sigma=True)
    
        I = (restwave > 3835) * (restwave < 4115)
        J = (restwave > 3835) * (restwave < 4115) * fnu_mask
    
        fig, ax = plt.subplots()
        dn4000_sdev = 1/np.sqrt(dn4000_ivar)
        ax.fill_between(restwave[I], fnu_obs[I]-fnu_sigma[I], fnu_obs[I]+fnu_sigma[I],
                        label=f'Observed Dn(4000)={dn4000:.3f}+/-{sdev:.3f}')
        ax.plot(restwave[I], fnu[I], color='blue', label=f'Line-free Dn(4000)={dn4000:.3f}+/-{dn4000_dev:.3f}')
        ax.plot(restwave[I], fnu_fullmodel[I], color='k', label=f'Model Dn(4000)={dn4000_model:.3f}')
        ax.plot(restwave[I], fnu_model[I], color='red', label=f'Model Dn(4000)={dn4000_model:.3f}')
        ylim = ax.get_ylim()
        ax.fill_between([3850, 3950], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                        color='lightgray', alpha=0.5)
        ax.fill_between([4000, 4100], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                        color='lightgray', alpha=0.5)
        ax.set_xlabel(r'Rest Wavelength ($\\AA$)')
        ax.set_ylabel(r'$F_{\\nu}$ (erg/s/cm2/Hz)')
        ax.legend()
        fig.savefig('desi-users/ioannis/tmp/qa-dn4000.png')
        """
        
    log.info(f'Emission-line fitting took {time.time()-tall:.2f} seconds.')

    if percamera_models:
        errmsg = 'percamera-models option not yet implemented.'
        log.critical(errmsg)
        raise NotImplementedError(errmsg)
    
    return modelspectra
