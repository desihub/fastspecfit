"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import time
import numpy as np
from enum import IntEnum
from itertools import chain

from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.photometry import Photometry
from fastspecfit.util import (C_LIGHT, TINY, SQTINY, F32MAX,
                              FLUXNORM, var2ivar)
from fastspecfit.emline_fit import (EMLine_Objective,
    EMLine_MultiLines, EMLine_find_peak_amplitudes,
    EMLine_build_model, EMLine_ParamsMapping)


class ParamType(IntEnum):
    AMPLITUDE = 0,
    VSHIFT = 1,
    SIGMA = 2


class EMFitTools(object):
    def __init__(self, emline_table, uniqueid=None, stronglines=False):

        self.line_table = emline_table
        self.uniqueid = uniqueid

        # restrict to just strong lines and assign to patches
        if stronglines:
            isstrong = self.line_table['isstrong'].value
            self.line_table = self.line_table[isstrong]

        # lines for which we want to measure moments
        self.moment_lines = {'CIV_1549': ['civ_1549'], 'MGII_2800': ['mgii_2796', 'mgii_2803'],
                             'HBETA': ['hbeta'], 'OIII_5007': ['oiii_5007']} #, 'HALPHA': ['halpha']}

        # Build some convenience (Boolean) arrays that we'll use in many places:
        line_isbroad  = self.line_table['isbroad'].value
        line_isbalmer = self.line_table['isbalmer'].value
        line_ishelium = self.line_table['ishelium'].value
        line_isstrong = self.line_table['isstrong'].value

        # broad UV/QSO lines but *not* broad Balmer or helium lines
        self.isBroad = (line_isbroad & ~line_isbalmer)

        # narrow lines (forbidden + Balmer + helium)
        self.isNarrow = ~line_isbroad

        # broad Balmer *and* helium lines
        self.isBalmerBroad = (line_isbroad & line_isbalmer)

        # broad Balmer *excluding* helium lines
        self.isBalmerBroad_noHelium = \
            (self.isBalmerBroad & ~line_ishelium)

        # broad Balmer lines used to test the narrow+broad model
        self.isBalmerBroad_noHelium_Strong = \
            (self.isBalmerBroad_noHelium & line_isstrong)

        line_names = self.line_table['name'].value

        # mapping to enable fast lookup of line number by name
        self.line_map = { line_name : line_idx for line_idx, line_name in enumerate(line_names) }

        # info about tied doublet lines
        doublet_lines = {
            # indx           source         ratio name
            'mgii_2796' : ( 'mgii_2803' , 'mgii_doublet_ratio' ) ,
            'oii_3726'  : ( 'oii_3729'  , 'oii_doublet_ratio'  ) ,
            'sii_6731'  : ( 'sii_6716'  , 'sii_doublet_ratio'  ) ,
        }

        # mapping from target -> source for all tied doublets
        doublet_src = np.full(len(self.line_table), -1, dtype=np.int32)
        for doublet_tgt in doublet_lines:
            target_line = self.line_map[ doublet_tgt ]
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


    def compute_inrange_lines(self, redshift, wavelims=(3600., 9900.), wavepad=5*0.8):
        """ Record which lines are within the limits of the cameras """

        zlinewave = self.line_table['restwave'].value * (1. + redshift)

        self.line_in_range = \
            ((zlinewave > (wavelims[0] + wavepad)) & \
             (zlinewave < (wavelims[1] - wavepad)))


    def build_linemodels(self, separate_oiii_fit=True):
        """Build emission line model tables, with and without
        suppression of broad lines.  Establish fixed (i.e., forced to
        zero) params, tying relationships between params, and doublet
        relationships for each model, as well as the relationship
        between lines and their parameters, which we record in the
        line table.

        Parameter fixing needs to know which lines are within the
        observed wavelength ranges of the cameras, so we first add
        this information to the line table.

        """

        def create_model(tying_info, forceFixed=[]):
            """Given the tying info for the model and the list of
            in-range lines, determine which parameters of the model
            are fixed and free, and create the model's table.

            We fix all parameters for out-of-range lines to zero,
            unless the parameter is tied to a parameter for an
            in-range line.  We also allow the caller to specify
            parameters that should be fixed regardless.

            """

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


        def tie_line(tying_info, line_params, source_linename, amp_factor=None):
            """Tie parameters of given line to source line. We don't tie the
            amplitude unless a tying factor is given for it.

            """
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
                        https://ui.adsabs.harvard.edu/abs/2007AIPC..895..313D/abstract

                        Note: The theoretical [OIII] 4959,5007 doublet *flux*
                        ratio is 2.993, so since we fit in velocity space the
                        constrained amplitude ratio has to be
                        2.993*4960.295/5008.240=2.9643.
                        """
                        tie_line(tying_info, line_params, 'oiii_5007', amp_factor=1./2.9643)
                    case 'nii_6548':
                        """
                        [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.02198e-21
                        [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.94901e-21
                        https://ui.adsabs.harvard.edu/abs/2023AdSpR..71.1219D/abstract

                        Note: The theoretical [NII] 6548,84 doublet *flux*
                        ratio is 3.049, so since we fit in velocity space the
                        constrained amplitude ratio has to be
                        3.049*6549.861/6585.273=3.0326
                        """
                        tie_line(tying_info, line_params, 'nii_6584', amp_factor = 1./3.0326)
                    case 'oii_7330':
                        """
                        [O2] (5-->2): airwave: 7318.9185 vacwave: 7320.9350 emissivity: 8.18137e-24
                        [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 2.40519e-23
                        [O2] (5-->3): airwave: 7329.6613 vacwave: 7331.6807 emissivity: 1.35614e-23
                        [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.27488e-23

                        This quadruplet ratio is sufficently poorly determined
                        that we are not going to apply the wavelength
                        correction used for [OIII] and [NII].

                        """
                        tie_line(tying_info, line_params, 'oii_7320', amp_factor = 1./1.225)
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
            if separate_oiii_fit:
                if not line_isbroad and not line_name in { 'nii_6584', 'oiii_4959', 'oiii_5007' }:
                    tie_line(tying_info, line_params, 'nii_6584')
            else:
                if not line_isbroad and line_name != 'oiii_5007':
                    tie_line(tying_info, line_params, 'oiii_5007')

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

            if separate_oiii_fit:
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


    def summarize_linemodel(self, linemodel):
        """Simple function to summarize an input linemodel."""

        def _print(line_mask):
            for line in np.where(line_mask)[0]:
                line_name   = self.line_table['name'][line]
                line_params = self.line_table['params'][line]

                for param in line_params:
                    param_name    = self.param_table['name'][param]
                    param_isfixed = linemodel['fixed'][param]
                    tiedtoparam   = linemodel['tiedtoparam'][param]

                    if tiedtoparam == -1: # not tied
                        print(f'{param_name:25s} ', end='')
                        print('FIXED' if param_isfixed else 'free')
                    else:
                        source_name = self.param_table['name'][tiedtoparam]
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


    def _initial_guesses_and_bounds(self, linepix, coadd_flux, contpix=None,
                                    initial_linesigma_broad=3000.,
                                    initial_linesigma_narrow=150.,
                                    initial_linesigma_balmer_broad=1000.,
                                    initial_linevshift_broad=0.,
                                    initial_linevshift_narrow=0.,
                                    initial_linevshift_balmer_broad=0.,
                                    subtract_local_continuum=False):
        """For all lines in the wavelength range of the data, get a good initial guess
        on the amplitudes and line-widths. This step is critical for cases like,
        e.g., 39633354915582193 (tile 80613, petal 05), which has strong narrow
        lines.

        linepix - dictionary of indices *defined on the coadded spectrum* for
          all lines in range

        If subtract_local_continuum=True, then contpix is mandatory.

        """
        from fastspecfit.util import quantile, median

        initials = np.empty(len(self.param_table), dtype=np.float64)
        bounds = np.empty((len(self.param_table), 2), dtype=np.float64)

        # a priori initial guesses and bounds
        minsigma_broad = 1. # 100.
        minsigma_narrow = 1.
        minsigma_balmer_broad = 1. # 100.0 # minsigma_narrow

        maxsigma_broad = 1e4
        maxsigma_narrow = 750.
        maxsigma_balmer_broad = 1e4

        maxvshift_broad = 2500.
        maxvshift_narrow = 500.
        maxvshift_balmer_broad = 2500.

        minamp = 0.
        maxamp = +1e5

        for iline, (line_isbalmer, line_isbroad, line_params) in \
            enumerate(self.line_table.iterrows('isbalmer', 'isbroad', 'params')):

            amp, vshift, sigma = line_params

            # initial values and bounds for line's parameters
            initials[amp] = 1.
            bounds[amp] = (minamp, maxamp)

            if line_isbroad:
                if line_isbalmer: # broad He+Balmer lines
                    initials[vshift] = initial_linevshift_balmer_broad
                    initials[sigma] = initial_linesigma_balmer_broad
                    bounds[vshift] = (-maxvshift_balmer_broad, +maxvshift_balmer_broad)
                    bounds[sigma] = (minsigma_balmer_broad, maxsigma_balmer_broad)
                else: # broad UV/QSO lines (non-Balmer)
                    initials[vshift] = initial_linevshift_broad
                    initials[sigma] = initial_linesigma_broad
                    bounds[vshift] = (-maxvshift_broad, +maxvshift_broad)
                    bounds[sigma] = (minsigma_broad, maxsigma_broad)
            else: # narrow He+Balmer lines, and forbidden lines
                initials[vshift] = initial_linevshift_narrow
                initials[sigma] = initial_linesigma_narrow
                bounds[vshift] = (-maxvshift_narrow, +maxvshift_narrow)
                bounds[sigma] = (minsigma_narrow, maxsigma_narrow)

        # Replace a priori initial values based on the data, with optional local
        # continuum subtraction.
        for linename in linepix.keys():

            onelinepix = linepix[linename]
            if contpix is not None:
                onecontpix = contpix[linename]

            if subtract_local_continuum:
                local = median(coadd_flux[onecontpix])
            else:
                local = 0.

            npix = len(onelinepix)
            if npix > 5:
                mnpx = np.maximum(onelinepix[npix//2]-3, 0)
                mxpx = np.minimum(onelinepix[npix//2]+3, onelinepix[-1])
                amp = np.max(coadd_flux[mnpx:mxpx] - local)
            else:
                amp = quantile(coadd_flux[onelinepix], 0.975) - local

            # update the bounds on the line-amplitude
            #bounds = [-np.min(np.abs(coadd_flux[linepix])), 3*np.max(coadd_flux[linepix])]
            mx = 5. * np.max(coadd_flux[onelinepix] - local)

            # record our initial gueses and bounds for the amplitude, unless
            # they are nonsensical
            if mx >= 0. and amp >= 0.:
                line = self.line_map[linename]
                amp_idx = self.line_table['params'][line, ParamType.AMPLITUDE]
                initials[amp_idx] = amp
                bounds[amp_idx] = np.array([0., mx])

        # Specialized parameters on the MgII, [OII], and [SII] doublet ratios. See
        # https://github.com/desihub/fastspecfit/issues/39.
        doublet_bounds = {
            'mgii_doublet_ratio' : (0.0, 10.0), # MgII 2796/2803
            'oii_doublet_ratio'  : (0.0,  2.0), # [OII] 3726/3729 # (0.5, 1.5) # (0.66, 1.4)
            'sii_doublet_ratio'  : (0.0,  2.0), # [SII] 6731/6716 # (0.5, 1.5) # (0.67, 1.2)
        }

        param_names = self.param_table['name'].value

        for iparam in self.doublet_idx:
            param_name = param_names[iparam]
            bounds[iparam] = doublet_bounds[param_name]
            initials[iparam] = 1.

        # Make sure all parameters lie within their respective bounds.
        for iparam, param_name in enumerate(self.param_table['name'].value):
            iv = initials[iparam]
            lb, ub = bounds[iparam]
            if iv < lb:
                errmsg = \
                    f'Initial parameter {param_name} is outside its bound, ' + \
                    f'{iv:.2f} < {lb:.2f}.'
                log.critical(errmsg)
                raise ValueError(errmsg)

            if iv > ub:
                errmsg = \
                    f'Initial parameter {param_name} is outside its bound, ' + \
                    f'{iv:.2f} > {ub:.2f}.'
                log.critical(errmsg)
                raise ValueError(errmsg)

        return initials, bounds


    def optimize(self, linemodel,
                 initials, param_bounds,
                 obs_bin_centers,
                 obs_bin_fluxes,
                 obs_weights,
                 redshift,
                 resolution_matrices,
                 camerapix,
                 continuum_patches=None,
                 debug=False):
        """Optimization routine.

        Returns
        -------
        linemodel passed as input (*not* a copy)
        fit continuum patches if requested

        Modifies
        --------
        linemodel:
          - adds/replaces 'value' column with values of fit line parameters
          - sets metadata 'obsamps' to vector of observed amplitudes
          - sets metadate 'nfev' to number of fcn evaluations by least_squares
        """
        from scipy.optimize import least_squares

        line_wavelengths = self.line_table['restwave'].value

        isFree      = linemodel['free'].value
        tiedtoparam = linemodel['tiedtoparam'].value
        tiedfactor  = linemodel['tiedfactor'].value

        params_mapping = EMLine_ParamsMapping(len(linemodel), isFree,
                                              tiedtoparam, tiedfactor,
                                              self.doublet_idx, self.doublet_src)
        nLineFree = np.sum(isFree)
        nPatches = len(continuum_patches) if continuum_patches is not None else 0
        nPatchFree = 2 * nPatches

        log_str = f"Optimizing {nLineFree} emission-line parameters"
        if nPatchFree > 0:
            log_str += f" and {nPatchFree} continuum patch parameters"
            log.debug(log_str)

        if nLineFree == 0:
            # corner case where all lines are out of the wavelength range, which can
            # happen at high redshift and with the red camera masked, e.g.,
            # iron/main/dark/6642/39633239580608311).
            linemodel.meta['nfev'] = 0
            return linemodel
        else:
            obj = EMLine_Objective(obs_bin_centers,
                                   obs_bin_fluxes,
                                   obs_weights,
                                   redshift,
                                   line_wavelengths,
                                   resolution_matrices,
                                   camerapix,
                                   params_mapping,
                                   continuum_patches=continuum_patches)

            objective = obj.objective
            jac       = obj.jacobian

            initial_guesses = np.empty(nLineFree + nPatchFree)
            bounds          = np.empty((nLineFree + nPatchFree, 2))

            # set line initial values and bounds
            initial_guesses[:nLineFree] = initials[isFree]
            bounds[:nLineFree] = param_bounds[isFree, :]

            if continuum_patches is not None:
                # set patch initial values and bounds
                initial_guesses[nLineFree:nLineFree+nPatches] = continuum_patches['slope']
                initial_guesses[nLineFree+nPatches:]          = continuum_patches['intercept']
                bounds[nLineFree:nLineFree+nPatches]          = continuum_patches['slope_bounds']
                bounds[nLineFree+nPatches:]                   = continuum_patches['intercept_bounds']

            # least_squares wants two arrays, not a 2D array
            bounds = ( bounds[:, 0], bounds[:, 1] )

            fit_info = least_squares(objective, initial_guesses, jac=jac, args=(),
                                     max_nfev=5000, xtol=1e-10, ftol=1e-5, #x_scale='jac' gtol=1e-10,
                                     tr_solver='lsmr', tr_options={'maxiter': 1000, 'regularize': True},
                                     method='trf', bounds=bounds,) # verbose=2)
            free_params = fit_info.x

            if not fit_info.success:
                errmsg = 'least_squares optimizer failed' + \
                    (f' for {self.uniqueid}' if self.uniqueid is not None else '')
                log.critical(errmsg)
                raise RuntimeError(errmsg)
            elif fit_info.status == 0:
                log.warning('optimizer failed to converge')

            # This should never happen if our optimizer enforces its bounds
            if np.any((free_params < bounds[0]) | (free_params > bounds[1])):
                errmsg = "ERROR: final parameters are not within requested bounds"
                log.critical(errmsg)
                raise RunTimeError(errmsg)

            if continuum_patches is not None:
                continuum_patches['slope']     = free_params[nLineFree:nLineFree+nPatches]
                continuum_patches['intercept'] = free_params[nLineFree+nPatches:]

            # translate free parame to full param array, but do NOT turn doublet
            # ratios into amplitudes yet, as out_linemodel needs them to be ratios
            parameters = params_mapping.mapFreeToFull(free_params[:nLineFree], patchDoublets=False)

            linemodel['value'] = parameters.copy() # protect from changes below
            linemodel.meta['nfev'] = fit_info['nfev']

            if continuum_patches is None:
                # convert doublet ratios to amplitudes
                parameters[self.doublet_idx] *= parameters[self.doublet_src]

                # calculate the observed maximum amplitude for each
                # fitted spectral line after convolution with the resolution
                # matrix.
                obsamps = EMLine_find_peak_amplitudes(
                    parameters, obs_bin_centers, redshift,
                    line_wavelengths, resolution_matrices,
                    camerapix)

                # add observed values as metadata, since they are only
                # relevant to amplitudes, not all parameters
                linemodel.meta['obsamps'] = obsamps
                return linemodel
            else:
                return linemodel, continuum_patches


    @staticmethod
    def chi2(linemodel, emlinewave, emlineflux, emlineivar, emlineflux_model,
             continuum_model=None, nfree_patches=0, return_dof=False):
        """Compute the reduced chi^2."""

        nfree = np.sum(linemodel['free'])
        nfree += nfree_patches

        dof = np.sum(emlineivar > 0) - nfree

        if dof > 0:
            if continuum_model is None:
                flux_model = emlineflux_model
            else:
                flux_model = emlineflux_model + continuum_model
            chi2 = np.sum(emlineivar * (emlineflux - flux_model)**2) / dof
        else:
            chi2 = 0.

        if return_dof:
            return chi2, dof, nfree
        else:
            return chi2


    def bestfit(self, linemodel, redshift, emlinewave, resolution_matrix,
                camerapix, continuum_patches=None):
        """Construct the best-fitting emission-line spectrum from a linemodel."""

        line_parameters = linemodel['value'].copy()

        # convert doublet ratios to amplitudes
        line_parameters[self.doublet_idx] *= line_parameters[self.doublet_src]

        linewaves = self.line_table['restwave'].value

        emlineflux_best = EMLine_build_model(redshift, line_parameters, linewaves,
                                             emlinewave, resolution_matrix, camerapix,
                                             continuum_patches=continuum_patches)

        return emlineflux_best


    def emlinemodel_bestfit(self, fastfit, redshift, emlinewave, resolution_matrix,
                            camerapix, snrcut=None):
        """Construct the best-fitting emission-line model
        from a fitted result structure (used below and in QA)"""

        line_parameters = np.array([fastfit[param] for param in self.param_table['modelname'] ])

        # convert doublet ratios to amplitudes
        line_parameters[self.doublet_idx] *= line_parameters[self.doublet_src]

        if snrcut is not None:
            lineamps = line_parameters[:len(self.line_table)] # amplitude parameters
            line_names = self.line_table['name'].value
            lineamps_ivar = [fastfit[line_name.upper()+'_AMP_IVAR'] for line_name in line_names]
            lineamps[lineamps * np.sqrt(lineamps_ivar) < snrcut] = 0.

        linewaves = self.line_table['restwave'].value

        emlineflux_best = EMLine_build_model(redshift, line_parameters, linewaves,
                                             emlinewave, resolution_matrix, camerapix)

        return emlineflux_best


    def populate_emtable(self, fastfit, linemodel, emlineflux_model,
                         emlinewave, emlineflux, emlineivar, oemlineivar,
                         specflux_nolines, redshift, resolution_matrices, camerapix,
                         results_monte=None, nminpix=7, nsigma=3., moment_nsigma=5.,
                         limitsigma_narrow_default=75., limitsigma_broad_default=1200.):
        """Populate the output table with the emission-line measurements.

        """
        from math import erf
        from fastspecfit.util import (centers2edges, sigmaclip, quantile,
                                      median, trapz)

        nline = len(self.line_table)
        nwave = len(emlinewave)
        dpixwave = median(np.diff(emlinewave)) # median pixel size [Angstrom]

        param_modelnames = self.param_table['modelname'].value


        def get_boundaries(A, v_lo, v_hi):
            """Find range (lo, hi) such that all pixels of A in range [v_lo,
            v_hi] lie in half-open interval [lo, hi).

            """
            return np.searchsorted(A, (v_lo, v_hi), side='right')


        def preprocess_linesigma(linesigma, linezwave, isbroad, isbalmer):
            """Pre-process linesigma by using a default value when a line is
            dropped and computing the width (in Angstroms) used to compute
            various quantities like boxflux.

            """
            # if the line was dropped, use a default value
            if linesigma == 0.:
                limitsigma_narrow = limitsigma_narrow_default
                limitsigma_broad = limitsigma_broad_default
                if isbroad and isbalmer:
                    linesigma = limitsigma_narrow
                else:
                    linesigma = limitsigma_broad

            linesigma_ang = linesigma * linezwave / C_LIGHT  # [observed-frame Angstrom]

            # require at least 2 pixels
            if linesigma_ang < 2. * dpixwave:
                linesigma_ang_window = 2. * dpixwave
            else:
                linesigma_ang_window = linesigma_ang

            return linesigma, linesigma_ang, linesigma_ang_window


        def get_continuum_pixels(emlinewave_s, linezwave, linesigma_ang_window):
            """Compute the pixels belonging to the continuum.

            """
            slo, elo = get_boundaries(emlinewave_s,
                                      linezwave - 10. * linesigma_ang_window,
                                      linezwave - 3. * linesigma_ang_window)
            shi, ehi = get_boundaries(emlinewave_s,
                                      linezwave + 3. * linesigma_ang_window,
                                      linezwave + 10. * linesigma_ang_window)
            borderindx = np.hstack((slo + np.where(oemlineivar_s[slo:elo] > 0.)[0],
                                    shi + np.where(oemlineivar_s[shi:ehi] > 0.)[0]))
            return borderindx


        # Where the cameras overlap, we have to account for the
        # variable pixel size by sorting in wavelength.
        Wsrt = np.argsort(emlinewave)

        emlinewave_s = emlinewave[Wsrt]
        emlineflux_s = emlineflux[Wsrt]
        emlineivar_s = emlineivar[Wsrt]
        oemlineivar_s = oemlineivar[Wsrt]
        emlineflux_model_s = emlineflux_model[Wsrt]
        specflux_nolines_s = specflux_nolines[Wsrt]

        dwaves = np.diff(centers2edges(emlinewave_s))

        values = linemodel['value'].value
        obsamps = linemodel.meta['obsamps']

        parameters = values.copy()
        parameters[self.doublet_idx] *= parameters[self.doublet_src]

        if results_monte is not None:
            values_monte, obsamps_monte, emlineflux_monte, specflux_nolines_monte = results_monte

            values_var = np.var(values_monte, axis=0)
            obsamps_var = np.var(obsamps_monte, axis=0)

            parameters_monte = values_monte.copy()
            parameters_monte[:, self.doublet_idx] *= parameters_monte[:, self.doublet_src]

            emlineflux_monte_s = emlineflux_monte[:, Wsrt]
            specflux_nolines_monte_s = specflux_nolines_monte[:, Wsrt]


        # initialize the line-stats table
        line_stats = Table()
        for stat in ['Z', 'SIGMA']:
            for groupname in ['NARROW', 'BROAD', 'UV']:
                line_stats[f'{groupname}_{stat}'] = np.zeros(1, 'f4')
                line_stats[f'{groupname}_{stat}RMS'] = np.zeros(1, 'f4')
        narrow_stats, broad_stats, uv_stats = [], [], []

        # iterate on each line
        for iline, (name, restwave, isbroad, isbalmer) in \
            enumerate(self.line_table.iterrows('name', 'restwave', 'isbroad', 'isbalmer')):

            linename = name.upper()
            line_amp, line_vshift, line_sigma = self.line_table['params'][iline]

            def get_fluxes(values, obsamps, emlineflux_s, specflux_nolines_s,
                           return_extras=False):
                """ Get all the computed fluxes associated with the current line.  Return the
                fluxes along with some intermediate quantities if return_extras is True.  (The
                extras are needed only if we are not using this function in Monte Carlo iteration.)

                """
                linez = redshift + values[line_vshift] / C_LIGHT
                linezwave = restwave * (1. + linez)
                linesigma0 = values[line_sigma] # original value [km/s]

                linesigma, linesigma_ang, linesigma_ang_window = \
                    preprocess_linesigma(linesigma0, linezwave, isbroad, isbalmer)

                line_s, line_e = get_boundaries(emlinewave_s,
                                                linezwave - nsigma * linesigma_ang_window,
                                                linezwave + nsigma * linesigma_ang_window)

                patchindx = line_s + np.where(emlineivar_s[line_s:line_e] > 0.)[0]

                # default values to return if not computed below
                emlineflux_patch = []
                boxflux = 0.
                flux = 0.
                cont, clipflux = 0., []

                # Are the pixels based on the original inverse spectrum fully masked?
                if line_s == line_e or np.sum(oemlineivar_s[line_s:line_e] == 0.) > 0.3 * (line_e - line_s):
                    patchindx = [] # return null patch
                else:
                    if len(patchindx) >= nminpix: # magic number: require at least XX unmasked pixels centered on the line
                        # boxcar integration of the flux
                        dwaves_patch = dwaves[patchindx]
                        emlineflux_patch = emlineflux_s[patchindx]
                        boxflux = np.sum(emlineflux_patch * dwaves_patch)

                        # require amp > 0 (line not dropped) to compute the flux
                        if obsamps[line_amp] > TINY:
                            # analytically integrated flux
                            flux = np.sqrt(2. * np.pi) * parameters[line_amp] * linezwave * linesigma0 / C_LIGHT

                        # next, get the continuum level
                        borderindx = get_continuum_pixels(emlinewave_s, linezwave, linesigma_ang_window)
                        if len(borderindx) >= nminpix: # require at least XX pixels to get the continuum level
                            clipflux, _ = sigmaclip(specflux_nolines_s[borderindx], low=3., high=3.)
                            if len(clipflux) == 0:
                                clipflux = specflux_nolines_s[borderindx]
                            cont = np.mean(clipflux)

                # needed by all code, including Monte Carlo
                res = (boxflux, flux, cont)

                if return_extras:
                    # needed by non-Monte Carlo code
                    extras = (linez, linesigma, linesigma_ang, patchindx, clipflux)
                    return (res, extras)
                else:
                    return res


            # zero out out-of-range lines
            if not self.line_in_range[iline]:
                obsamps[line_amp] = 0.
                parameters[line_amp] = 0.
                values[line_amp] = 0.
                values[line_vshift] = 0.
                values[line_sigma] = 0.
                continue

            # Special-case: populate the results table with the 'free' doublet
            # ratio parameters.
            if 'DOUBLET_RATIO' in param_modelnames[line_amp]:
                fastfit[param_modelnames[line_amp]] = values[line_amp]
                if results_monte is not None:
                    doublet_ivar = var2ivar(values_var[line_amp])
                    if doublet_ivar < F32MAX:
                        fastfit[f'{param_modelnames[line_amp]}_IVAR'] = doublet_ivar

            # Also store the 'tied' doublet ratios (fragile...).
            if linemodel['tiedtoparam'][line_amp] != -1:
                ratio = linemodel['tiedfactor'][line_amp]
                match linename:
                    case 'OIII_4959':
                        col = 'OIII_DOUBLET_RATIO'
                    case 'NII_6548':
                        col = 'NII_DOUBLET_RATIO'
                    case 'OII_7330':
                        col = 'OIIRED_DOUBLET_RATIO'
                    case _:
                        errmsg = 'Unrecognized tied doublet {linename}'
                        log.critical(errmsg)
                        raise ValueError(errmsg)
                fastfit[col] = 1. / ratio
                fastfit[f'{col}_IVAR'] = 0. # not optimized


            (boxflux, flux, cont), extras = get_fluxes(
                values, obsamps, emlineflux_s, specflux_nolines_s,
                return_extras=True)

            (linez, linesigma, linesigma_ang, patchindx, clipflux) = extras

            npix = len(patchindx)
            fastfit[f'{linename}_NPIX'] = npix

            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next
            # line.
            if npix == 0:
                obsamps[line_amp] = 0.
                parameters[line_amp] = 0.
                values[line_amp] = 0.
                values[line_vshift] = 0.
                values[line_sigma] = 0.
                continue

            flux_ivar, cont_ivar = 0., 0. # defaults if not computed below

            if npix >= nminpix: # magic number: require at least XX unmasked pixels centered on the line
                fastfit[f'{linename}_AMP'] = obsamps[line_amp]
                fastfit[f'{linename}_VSHIFT'] = values[line_vshift]
                fastfit[f'{linename}_SIGMA'] = values[line_sigma]
                fastfit[f'{linename}_MODELAMP'] = parameters[line_amp]

                fastfit[f'{linename}_BOXFLUX'] = boxflux
                fastfit[f'{linename}_FLUX'] = flux
                fastfit[f'{linename}_CONT'] = cont

                emlineflux_patch = emlineflux_s[patchindx]
                emlineivar_patch = emlineivar_s[patchindx]
                if np.any(emlineivar_patch == 0.):
                    errmsg = 'Ivar should never be zero within an emission line!'
                    log.critical(errmsg)
                    raise ValueError(errmsg)

                if results_monte is not None:
                    res = [get_fluxes(vv, oo, lf, sfnl) for  vv, oo, lf, sfnl in
                           zip(values_monte, obsamps_monte, emlineflux_monte_s,
                               specflux_nolines_monte_s)]
                    boxflux_monte, flux_monte, cont_monte = tuple(zip(*res))
                    flux_monte = np.array(flux_monte)
                    cont_monte = np.array(cont_monte)

                    # Compute the variance on the line-fitting results.
                    obsamps_ivar = var2ivar(obsamps_var[line_amp])
                    vshift_ivar = var2ivar(values_var[line_vshift])
                    sigma_ivar = var2ivar(values_var[line_sigma])
                    boxflux_ivar = var2ivar(np.var(boxflux_monte))
                    if obsamps_ivar < F32MAX:
                        fastfit[f'{linename}_AMP_IVAR'] = obsamps_ivar
                    if vshift_ivar < F32MAX:
                        fastfit[f'{linename}_VSHIFT_IVAR'] = vshift_ivar
                    if sigma_ivar < F32MAX:
                        fastfit[f'{linename}_SIGMA_IVAR'] = sigma_ivar
                    if boxflux_ivar < F32MAX:
                        fastfit[f'{linename}_BOXFLUX_IVAR'] = boxflux_ivar
                    #fastfit[f'{linename}_MODELAMP_IVAR'] = var2ivar(parameters_var[line_amp])
                else:
                    obsamps_ivar = 0.

                # require amp > 0 (line not dropped) to compute the flux and chi2
                if obsamps[line_amp] > TINY:
                    emlineflux_model_patch = emlineflux_model_s[patchindx]
                    chi2 = np.sum(emlineivar_patch * (emlineflux_patch - emlineflux_model_patch)**2)
                    fastfit[f'{linename}_CHI2'] = chi2

                    if results_monte is not None:
                        flux_ivar = var2ivar(np.var(flux_monte))
                        if flux_ivar < F32MAX:
                            fastfit[f'{linename}_FLUX_IVAR'] = flux_ivar

                    # keep track of sigma and z but only using XX-sigma lines
                    linesnr = obsamps[line_amp] * np.sqrt(obsamps_ivar)
                    if linesnr > 1.5:
                        if isbroad: # includes UV and broad Balmer lines
                            if isbalmer:
                                stats = broad_stats
                            else:
                                stats = uv_stats
                        else:
                            stats = narrow_stats
                        stats.append((linesigma, linez))

            if results_monte is not None:
                cont_ivar = var2ivar(np.var(cont_monte))
                if cont_ivar < F32MAX:
                    fastfit[f'{linename}_CONT_IVAR'] = cont_ivar

            if cont != 0. and cont_ivar > 0.:
                # upper limit on the flux is defined by snrcut*cont_err*sqrt(2*pi)*linesigma
                fluxlimit = np.sqrt(2. * np.pi) * linesigma_ang / np.sqrt(cont_ivar) # * u.erg/(u.second*u.cm**2)
                fastfit[f'{linename}_FLUX_LIMIT'] = fluxlimit

                #ewlimit = fluxlimit * cont / (1.+redshift)
                #fastfit[f'{linename}_EW_LIMIT'] = ewlimit

                if flux > 0. and flux_ivar > 0.:
                    # add the uncertainties in the flux and continuum in quadrature
                    ew = flux / cont / (1. + redshift) # rest frame [A]
                    fastfit[f'{linename}_EW'] = ew

                    if results_monte is not None:
                        I = cont_monte != 0.
                        if np.sum(I) > 2:
                            ew_monte = flux_monte[I] / cont_monte[I] / (1. + redshift) # rest frame [A]
                            ew_ivar = var2ivar(np.var(ew_monte))
                            if ew_ivar < F32MAX:
                                fastfit[f'{linename}_EW_IVAR'] = ew_ivar


        # Measure moments for the set of lines in self.moment_lines. We need a
        # separate loop because for one "line" (MgII) we actually want the full
        # doublet.
        for moment_col, moment_line_names in self.moment_lines.items():
            moment_lines = [ self.line_map[name] for name in moment_line_names if name in self.line_map ]
            if len(moment_lines) == 0:
                continue

            restwave = np.mean(self.line_table['restwave'][moment_lines])
            mline = self.line_table[moment_lines[0]] # take the zeroth line in the case of a doublet

            _, line_vshift, line_sigma = mline['params']
            isbroad, isbalmer = mline['isbroad'], mline['isbalmer']

            def get_moments(values, emlineflux_s):
                """Get first three (non-parametric) moments of the flux
                distribution in a given patch centered on a given line.

                """
                linezwave = restwave * (1. + redshift + values[line_vshift] / C_LIGHT)
                linesigma = values[line_sigma] # [km/s]

                linesigma, _, linesigma_ang_window = preprocess_linesigma(
                    linesigma, linezwave, isbroad, isbalmer)

                ss, ee = get_boundaries(emlinewave_s,
                                        linezwave - moment_nsigma * linesigma_ang_window,
                                        linezwave + moment_nsigma * linesigma_ang_window)

                ww = emlinewave_s[ss:ee]
                ff = emlineflux_s[ss:ee]
                patchnorm = np.sum(ff)
                if patchnorm == 0.: # could happen I guess
                    return 0., 0., 0.
                else:
                    # compute the first three moments of the distribution
                    moment1 = np.sum(ww * ff) / patchnorm               # [Angstrom]
                    moment2 = np.sum((ww-moment1)**2 * ff) / patchnorm  # [Angstrom**2]
                    moment3 = np.sum((ww-moment1)**3 * ff) / patchnorm  # [Angstrom**3]
                    return moment1, moment2, moment3


            moment1, moment2, moment3 = get_moments(values, emlineflux_s)

            for n, mom in enumerate((moment1, moment2, moment3)):
                fastfit[f'{moment_col}_MOMENT{n+1}'] = mom

            if results_monte is not None:
                res = [get_moments(v, ef) for v, ef in zip(values_monte, emlineflux_monte_s)]
                moments_monte = tuple(zip(*res))

                for n, mom_monte in enumerate(moments_monte):
                    mom_ivar = var2ivar(np.var(mom_monte))
                    if mom_ivar < F32MAX:
                        fastfit[f'{moment_col}_MOMENT{n+1}_IVAR'] = mom_ivar

        # get the per-group average emission-line redshifts and velocity widths
        for stats, groupname in zip((narrow_stats, broad_stats, uv_stats),
                                    ('NARROW', 'BROAD', 'UV')):
            if len(stats) > 0:
                stats = np.array(stats)
                sigmas = stats[:, 0]
                redshifts = stats[:, 1]

                line_stats[f'{groupname}_SIGMA'] = np.mean(sigmas)
                line_stats[f'{groupname}_SIGMARMS'] = np.std(sigmas)
                line_stats[f'{groupname}_Z'] = np.mean(redshifts)
                line_stats[f'{groupname}_ZRMS'] = np.std(redshifts)
            else:
                line_stats[f'{groupname}_Z'] = redshift


        import logging
        if log.getEffectiveLevel() == logging.DEBUG:
            for ln in self.line_table['name'].value:
                linename = ln.upper()
                for col in ('VSHIFT', 'SIGMA', 'MODELAMP', 'AMP', 'AMP_IVAR', 'CHI2', 'NPIX'):
                    val = fastfit[f'{linename}_{col}']
                    log.debug(f'{linename} {col}: {val:.4f}')
                for col in ('FLUX', 'BOXFLUX', 'FLUX_IVAR', 'BOXFLUX_IVAR', 'CONT', 'CONT_IVAR', 'EW', 'EW_IVAR', 'FLUX_LIMIT'):
                    val = fastfit[f'{linename}_{col}']
                    log.debug(f'{linename} {col}: {val:.4f}')
                print()

            for lname in ['MGII', 'OII', 'SII']:
                col = f'{lname}_DOUBLET_RATIO'
                val = fastfit[col]
                val_ivar = fastfit[f'{col}_IVAR']
                log.debug(f'{col}: {val:.4f}')
                log.debug(f'{col}_IVAR: {val_ivar:.4f}')
            print()

        return line_stats


def synthphot_spectrum(phot, data, specphot, modelwave, modelflux):
    """Synthesize photometry from the best-fitting model (continuum+emission lines).

    """
    filters = phot.synth_filters[data['photsys']]

    synthmaggies = Photometry.get_ab_maggies(filters, modelflux / FLUXNORM, modelwave)
    model_synthmag = Photometry.to_nanomaggies(synthmaggies)  # units of nanomaggies

    model_synthphot = Photometry.parse_photometry(
        phot.synth_bands, maggies=synthmaggies, nanomaggies=False,
        lambda_eff=filters.effective_wavelengths.value)

    synthmag = data['synthphot']['nanomaggies'].value
    model_synthmag = model_synthphot['nanomaggies'].value
    for iband, band in enumerate(phot.synth_bands):
        bname =  band.upper()
        specphot[f'FLUX_SYNTH_{bname}'] = synthmag[iband] # * 'nanomaggies'
        specphot[f'FLUX_SYNTH_SPECMODEL_{bname}'] = model_synthmag[iband] # * 'nanomaggies'


def build_coadded_models(data, emlinewave, emlineflux_model, continuum_flux,
                         smooth_continuum_flux):
    """Package together the output models, coadded across cameras.

    Assume constant dispersion in wavelength!

    """
    from astropy.table import Column

    # I believe that all the elements of the coadd_wave vector are contained
    # within one or more of the per-camera wavelength vectors, and so we
    # should be able to simply map our model spectra with no
    # interpolation. However, because of round-off, etc., it's probably
    # easiest to use np.interp.
    coadd_waves = data['coadd_wave']
    minwave = np.min(coadd_waves)
    maxwave = np.max(coadd_waves)
    dwave = coadd_waves[1] - coadd_waves[0]

    minwave = np.floor(minwave * 1000.) / 1000
    maxwave = np.floor(maxwave * 1000.) / 1000
    dwave = np.round(dwave, decimals=3)

    npix = int(np.round((maxwave-minwave)/dwave)) + 1
    wave_out = minwave + dwave * np.arange(npix, dtype=np.float64)

    wavesrt = np.argsort(emlinewave)
    sorted_waves = emlinewave[wavesrt]
    continuum_out = np.interp(wave_out, sorted_waves, continuum_flux[wavesrt])
    smooth_continuum_out = np.interp(wave_out, sorted_waves, smooth_continuum_flux[wavesrt])
    emlineflux_out = np.interp(wave_out, sorted_waves, emlineflux_model[wavesrt])

    spectra_out = Table(
        # ensure that these columns will stack as rows when
        # we vstack the Tables for different spectra, rather
        # than being concatenated into one long row.
        data=(
            Column(name='CONTINUUM', dtype='f4',
                   data=continuum_out.reshape(1, npix)),
            Column(name='SMOOTHCONTINUUM', dtype='f4',
                   data=smooth_continuum_out.reshape(1, npix)),
            Column(name='EMLINEMODEL', dtype='f4',
                   data=emlineflux_out.reshape(1, npix))
        ),
        # all these header cards need to be 2-element tuples (value, comment),
        # otherwise io.write_fastspecfit will crash
        meta = {
            'NAXIS1': (npix, 'number of pixels'),
            'NAXIS2': (npix, 'number of models'),
            'NAXIS3': (npix, 'number of objects'),
            'BUNIT':  ('10**-17 erg/(s cm2 Angstrom)', 'flux unit'),
            'CUNIT1': ('Angstrom', 'wavelength unit'),
            'CTYPE1': ('WAVE', 'type of axis'),
            'CRVAL1': (minwave, 'wavelength of pixel CRPIX1 (Angstrom)'),
            'CRPIX1': (0, '0-indexed pixel number corresponding to CRVAL1'),
            'CDELT1': (dwave, 'pixel size (Angstrom)'),
            'DC-FLAG': (0, '0 = linear wavelength vector'),
            'AIRORVAC': ('vac', 'wavelengths in vacuum (vac)')
        },
        copy=False
    )

    return wave_out, continuum_out, emlineflux_out, spectra_out


def test_broad_model(EMFit,
                     linemodel_nobroad, emlineflux_model_nobroad,
                     linemodel_broad, emlineflux_model_broad,
                     emlinewave, emlineflux, emlineivar, redshift,
                     minsnr_balmer_broad, minsigma_balmer_broad):
    """Test the broad Balmer-line hypothesis.

    """
    from fastspecfit.util import quantile
    from fastspecfit.linemasker import LineMasker

    residuals = emlineflux - emlineflux_model_broad

    broad_values = linemodel_broad['value'].value

    line_names = EMFit.line_table['name'].value
    line_params = EMFit.line_table['params'].value

    # get the pixels of the broad Balmer lines
    IBalmer = EMFit.isBalmerBroad_noHelium_Strong

    balmer_linesigmas = broad_values[line_params[IBalmer, ParamType.SIGMA ] ]
    balmer_linevshifts = broad_values[line_params[IBalmer, ParamType.VSHIFT] ]

    balmerpix = LineMasker.linepix_and_contpix(
        emlinewave, emlineivar, EMFit.line_table[IBalmer],
        balmer_linesigmas, get_contpix=False, redshift=redshift)

    balmerlines =  [EMFit.line_map[ln] for ln in balmerpix['linepix']]
    balmerpixels = [px for px in balmerpix['linepix'].values()]

    # Determine how many lines (free parameters) are in wavelengths in and
    # around the Balmer lines, with and without broad lines.
    balmer_nfree_broad = 0
    balmer_nfree_nobroad = 0

    zlinewaves = EMFit.line_table['restwave'] * (1. + redshift)

    balmer_linesnrs = np.zeros(len(balmerlines))
    for iline, (bpix, bline) in enumerate(zip(balmerpixels, balmerlines)):
        bpixwave = emlinewave[bpix]
        line_in_balmerpix = ( (zlinewaves > np.min(bpixwave)) &
                              (zlinewaves < np.max(bpixwave)) )

        for xline in np.where(line_in_balmerpix)[0]:
            params = line_params[xline]
            balmer_nfree_nobroad += np.sum(linemodel_nobroad['free'][params])
            balmer_nfree_broad += np.sum(linemodel_broad['free'][params])

        # get the S/N of the broad Balmer line
        lo, hi = quantile(residuals[bpix], (0.25, 0.75))
        bnoise = (hi - lo) / 1.349 # robust sigma
        bindx = line_params[bline, ParamType.AMPLITUDE]
        if bnoise > 0.:
            balmer_linesnrs[iline] = linemodel_broad.meta['obsamps'][bindx] / bnoise

    # compute delta-chi2 around just the broad, non-helium Balmer lines
    balmerpixels = np.unique(np.hstack(balmerpixels))

    bivar = emlineivar[balmerpixels]
    bflux = emlineflux[balmerpixels]
    nbpix = np.sum(bivar > 0)
    balmer_ndof_broad = nbpix - balmer_nfree_broad
    balmer_ndof_nobroad = nbpix - balmer_nfree_nobroad

    linechi2_balmer_broad = np.sum(bivar * (bflux - emlineflux_model_broad[balmerpixels])**2)
    linechi2_balmer_nobroad = np.sum(bivar * (bflux - emlineflux_model_nobroad[balmerpixels])**2)
    delta_linechi2_balmer = linechi2_balmer_nobroad - linechi2_balmer_broad
    delta_linendof_balmer = balmer_ndof_nobroad - balmer_ndof_broad

    # Choose broad-line model if:
    # --delta-chi2 > delta-ndof
    # --broad_sigma < narrow_sigma
    # --broad_sigma < 250

    dchi2test = (delta_linechi2_balmer > delta_linendof_balmer)

    Hanarrow_idx = line_params[EMFit.line_map['halpha'], ParamType.SIGMA]
    Hanarrow = linemodel_broad['value'][Hanarrow_idx]

    Habroad_idx = line_params[EMFit.line_map['halpha_broad'], ParamType.SIGMA]
    Habroad = linemodel_broad['value'][Habroad_idx]

    sigtest1 = Habroad > minsigma_balmer_broad
    sigtest2 = Habroad > Hanarrow

    if len(balmer_linesnrs) == 1:
        broadsnrtest = (balmer_linesnrs[-1] > minsnr_balmer_broad)
        _broadsnr = f'S/N {line_names[balmerlines[-1]]} = {balmer_linesnrs[-1]:.1f}'
    else:
        broadsnrtest = np.any(balmer_linesnrs[-2:] > minsnr_balmer_broad)
        _broadsnr = \
            f'S/N ({line_names[balmerlines[-2]]}) = {balmer_linesnrs[-2]:.1f}, ' \
            f'S/N ({line_names[balmerlines[-1]]}) = {balmer_linesnrs[-1]:.1f}'

    if dchi2test and sigtest1 and sigtest2 and broadsnrtest:
        adopt_broad = True
        log.info('Adopting broad-line model:')
        log.info(f'  delta-chi2={delta_linechi2_balmer:.1f} > delta-ndof={delta_linendof_balmer:.0f}')
        log.info(f'  sigma_broad={Habroad:.1f} km/s, sigma_narrow={Hanarrow:.1f} km/s')
        if _broadsnr:
            log.info(f'  {_broadsnr} > {minsnr_balmer_broad:.0f}')

    else:
        adopt_broad = False
        if dchi2test == False:
            log.debug(f'Dropping broad-line model: delta-chi2={delta_linechi2_balmer:.1f} < delta-ndof={delta_linendof_balmer:.0f}')
        elif sigtest1 == False:
            log.debug(f'Dropping broad-line model: Halpha_broad_sigma {Habroad:.1f} km/s < {minsigma_balmer_broad:.0f} km/s '
                      f'(delta-chi2={delta_linechi2_balmer:.1f}, delta-ndof={delta_linendof_balmer:.0f}).')
        elif sigtest2 == False:
            log.debug(f'Dropping broad-line model: Halpha_broad_sigma {Habroad:.1f} km/s < Halpha_narrow_sigma {Hanarrow:.1f} km/s '
                      f'(delta-chi2={delta_linechi2_balmer:.1f}, delta-ndof={delta_linendof_balmer:.0f}).')
        elif broadsnrtest == False:
            log.debug(f'Dropping broad-line model: {_broadsnr} < {minsnr_balmer_broad:.0f}')

    return adopt_broad, delta_linechi2_balmer, delta_linendof_balmer


def linefit(EMFit, linemodel, initial_guesses, param_bounds,
            emlinewave, emlineflux, emlineivar, weights, redshift,
            resolution_matrix, camerapix, debug=False):
    """Simple wrapper on EMFit.optimize.
    Modifies linemodel in place.
    """

    EMFit.optimize(linemodel, initial_guesses, param_bounds,
                   emlinewave, emlineflux, weights, redshift,
                   resolution_matrix, camerapix, debug=debug)

    emlineflux_model = EMFit.bestfit(linemodel, redshift, emlinewave,
                                     resolution_matrix, camerapix)
    chi2, ndof, nfree = EMFit.chi2(linemodel,
                                   emlinewave, emlineflux, emlineivar,
                                   emlineflux_model, return_dof=True)

    return emlineflux_model, nfree, chi2


def emline_specfit(data, fastfit, specphot, continuummodel, smooth_continuum,
                   phot, emline_table, minsnr_balmer_broad=2.5,
                   minsigma_balmer_broad=250., continuummodel_monte=None,
                   specflux_monte=None, synthphot=True, broadlinefit=True,
                   debug_plots=False):
    """Perform the fit minimization / chi2 minimization.

    Parameters
    ----------
    data
    continuummodel
    smooth_continuum
    synthphot
    broadlinefit
    minsigma_balmer_broad - minimum broad-line sigma [km/s]

    Returns
    -------
    results
    modelflux

    """
    from astropy.table import vstack
    from fastspecfit.util import ivar2var

    tall = time.time()

    EMFit = EMFitTools(emline_table, uniqueid=data['uniqueid'])

    redshift = data['redshift']
    camerapix = data['camerapix']
    resolution_matrix = data['res']

    # Combine pixels across all cameras
    emlinewave = np.hstack(data['wave'])
    oemlineivar = np.hstack(data['ivar'])
    specflux = np.hstack(data['flux'])

    # portion of actual flux predicted to be due to emission lines
    emlineflux = specflux - continuummodel - smooth_continuum

    emlineivar = np.copy(oemlineivar)
    _, emlinegood = ivar2var(emlineivar, clip=1e-3)
    emlinebad = ~emlinegood

    # This is a (dangerous???) hack.
    if np.any(emlinebad):
        emlineivar[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineivar[emlinegood])
        emlineflux[emlinebad] = np.interp(emlinewave[emlinebad], emlinewave[emlinegood], emlineflux[emlinegood]) # ???

    weights = np.sqrt(emlineivar)

    # Monte Carlo spectrum carried over from continuum-fitting. Assume that the
    # smooth continuum model is the same...
    if specflux_monte is not None:
        nmonte = len(specflux_monte)
        if continuummodel_monte is not None:
            emlineflux_monte = (specflux_monte - continuummodel_monte - \
                                smooth_continuum[np.newaxis, :])
        else:
            emlineflux_monte = (specflux_monte - continuummodel[np.newaxis, :] - \
                                smooth_continuum[np.newaxis, :])
    else:
        nmonte = 0

    # determine which lines are in range of the camera
    EMFit.compute_inrange_lines(redshift, wavelims=(np.min(emlinewave),
                                                    np.max(emlinewave)))

    # Create initial line models for broad and nobroad cases, and
    # get initial guesses for their parameters.  We'll fit both models
    # to the data below and pick the one that fits better.

    linemodel_broad, linemodel_nobroad = EMFit.build_linemodels(separate_oiii_fit=True)
    #EMFit.summarize_linemodel(linemodel_nobroad)
    #EMFit.summarize_linemodel(linemodel_broad)

    coadd_flux = np.interp(data['coadd_wave'], emlinewave, emlineflux)

    initial_guesses, param_bounds = EMFit._initial_guesses_and_bounds(
        data['coadd_linepix'], coadd_flux,
        initial_linesigma_broad=data['linesigma_broad'],
        initial_linesigma_narrow=data['linesigma_narrow'],
        initial_linesigma_balmer_broad=data['linesigma_balmer_broad'],
        initial_linevshift_broad=data['linevshift_broad'],
        initial_linevshift_narrow=data['linevshift_narrow'],
        initial_linevshift_balmer_broad=data['linevshift_balmer_broad'])

    # fit spectrum without broad Balmer lines
    t0 = time.time()

    # updates linemodel_nobroad
    emlineflux_model_nobroad, nfree_nobroad, chi2_nobroad = linefit(
        EMFit, linemodel_nobroad, initial_guesses, param_bounds,
        emlinewave, emlineflux, emlineivar, weights, redshift,
        resolution_matrix, camerapix, debug=False)

    log.debug(f'Line-fitting {data["uniqueid"]} with no broad lines and {nfree_nobroad} free parameters took ' + \
              f'{time.time()-t0:.4f} seconds [niter={linemodel_nobroad.meta["nfev"]}, rchi2={chi2_nobroad:.4f}].')

    # Now try to improve the chi2 by adding broad Balmer lines.  Save the preferred line modela and stats
    if broadlinefit and data['balmerbroad']:
        t0 = time.time()

        # updates linemodel_broad
        emlineflux_model_broad, nfree_broad, chi2_broad = linefit(
            EMFit, linemodel_broad, initial_guesses, param_bounds,
            emlinewave, emlineflux, emlineivar, weights, redshift,
            resolution_matrix, camerapix, debug=False)

        log.debug(f'Line-fitting {data["uniqueid"]} with broad lines and {nfree_broad} free parameters took ' + \
                  f'{time.time()-t0:.4f} seconds [niter={linemodel_broad.meta["nfev"]}, rchi2={chi2_broad:.4f}].')

        adopt_broad, delta_linechi2_balmer, delta_linendof_balmer = \
            test_broad_model(EMFit,
                             linemodel_nobroad, emlineflux_model_nobroad,
                             linemodel_broad, emlineflux_model_broad,
                             emlinewave, emlineflux, emlineivar, redshift,
                             minsnr_balmer_broad, minsigma_balmer_broad)

        if adopt_broad:
            linemodel_pref, emlineflux_model_pref, chi2_pref = linemodel_broad, emlineflux_model_broad, chi2_broad
        else:
            linemodel_pref, emlineflux_model_pref, chi2_pref = linemodel_nobroad, emlineflux_model_nobroad, chi2_nobroad
    else:
        adopt_broad = False
        if not broadlinefit:
            log.info('Skipping broad-line fitting (broadlinefit=False).')
        elif not data['balmerbroad']:
            log.info('Skipping broad-line fitting (no broad Balmer lines in the spectral range).')

        linemodel_pref, emlineflux_model_pref, chi2_pref = linemodel_nobroad, emlineflux_model_nobroad, chi2_nobroad
        delta_linechi2_balmer, delta_linendof_balmer = 0, np.int32(0)

    # Residual spectrum with no emission lines
    specflux_nolines = specflux - emlineflux_model_pref

    if nmonte > 0:
        # Monte Carlo to get the uncertainties on the derived parameters.
        def get_results(emlineflux):
            linemodel = linemodel_pref.copy() # avoid stepping on original line model

            # updates linemodel in place
            emlineflux_model, _, _ = linefit(
                EMFit, linemodel, initial_guesses, param_bounds,
                emlinewave, emlineflux, emlineivar,
                weights, redshift, resolution_matrix, camerapix)
            values = linemodel['value'].value
            obsamps = linemodel.meta['obsamps'] # observed amplitudes
            return (values, obsamps, emlineflux_model)

        res = [get_results(emlf) for emlf in emlineflux_monte]
        values_monte, obsamps_monte, emlineflux_model_monte = \
            tuple(zip(*res))
        values_monte  = np.array(values_monte)  # used as array later
        obsamps_monte = np.array(obsamps_monte) # used as array later

        specflux_nolines_monte = specflux_monte - emlineflux_model_monte

        results_monte = (values_monte, obsamps_monte,
                         emlineflux_monte, specflux_nolines_monte)
    else:
        results_monte = None

    # Now fill the output table.
    line_stats = EMFit.populate_emtable(
        fastfit, linemodel_pref, emlineflux_model_pref, emlinewave,
        emlineflux, emlineivar, oemlineivar, specflux_nolines,
        redshift, resolution_matrix, camerapix, results_monte=results_monte)

    msg = []
    dv = C_LIGHT*(np.array([line_stats['UV_Z'][0], line_stats['BROAD_Z'][0], line_stats['NARROW_Z'][0]])-redshift)
    dverr = C_LIGHT*np.array([line_stats['UV_ZRMS'][0], line_stats['BROAD_ZRMS'][0], line_stats['NARROW_ZRMS'][0]])
    for label, units, val, valerr in zip(
            ['delta(v) UV', 'Balmer broad', 'narrow'],
            [' km/s', ' km/s', ' km/s'], dv, dverr):
        err_msg = f'+/-{valerr:.1f}' if valerr > 0. else ''
        msg.append(f'{label}={val:.1f}{err_msg}{units}')
    log.info(' '.join(msg))

    msg = []
    for label, units, val, valerr in zip(
            ['sigma UV', 'Balmer broad', 'narrow'],
            [' km/s', ' km/s', ' km/s'],
            [line_stats['UV_SIGMA'][0], line_stats['BROAD_SIGMA'][0], line_stats['NARROW_SIGMA'][0]],
            [line_stats['UV_SIGMARMS'][0], line_stats['BROAD_SIGMARMS'][0], line_stats['NARROW_SIGMARMS'][0]]):
        err_msg = f'+/-{valerr:.0f}' if valerr > 0. else ''
        msg.append(f'{label}={val:.0f}{err_msg}{units}')
    log.info(' '.join(msg))

    # Build the model spectrum from the reported parameter values
    emlineflux_model_best = EMFit.emlinemodel_bestfit(
        fastfit, redshift, emlinewave, resolution_matrix, camerapix)

    specphot['RCHI2_LINE'] = chi2_pref
    #fastfit['NDOF_LINE'] = ndof_pref
    fastfit['DELTA_LINECHI2'] = delta_linechi2_balmer  # chi2_nobroad - chi2_broad
    fastfit['DELTA_LINENDOF'] = delta_linendof_balmer  # ndof_nobroad - ndof_broad

    # full-fit reduced chi2
    rchi2 = np.sum(oemlineivar * (specflux - (continuummodel + smooth_continuum + emlineflux_model_best))**2)
    rchi2 /= np.sum(oemlineivar > 0)  # dof??
    specphot['RCHI2'] = rchi2

    # Build the output model spectra.
    wave_out, continuum_out, emlineflux_out, spectra_out = build_coadded_models(
        data, emlinewave, emlineflux_model_best, continuummodel, smooth_continuum)

    # Optionally synthesize photometry (excluding the smoothcontinuum!)
    if synthphot:
        specflux_out = continuum_out + emlineflux_out
        synthphot_spectrum(phot, data, specphot, wave_out, specflux_out)

    # measure DN(4000) without the emission lines
    if specphot['DN4000_IVAR'] > 0.:
        flux_out_nolines = data['coadd_flux'] - emlineflux_out

        dn4000_nolines, _ = Photometry.get_dn4000(wave_out, flux_out_nolines,
                                                  redshift=redshift, rest=False)
        log.info(f'Dn(4000)={dn4000_nolines:.3f} in the emission-line subtracted spectrum.')
        specphot['DN4000'] = dn4000_nolines

        # Simple QA of the Dn(4000) estimate.
        if debug_plots:
            dn4000_ivar = specphot['DN4000_IVAR']
            if dn4000_ivar == 0.:
                log.info('Dn(4000) not measured; unable to generate QA figure.')
            else:
                import matplotlib.pyplot as plt
                import seaborn as sns

                pngfile = f'qa-dn4000-{data["uniqueid"]}.png'
                sns.set(context='talk', style='ticks', font_scale=0.7)

                dn4000, dn4000_obs = specphot['DN4000'], specphot['DN4000_OBS']
                dn4000_model, dn4000_model_ivar = specphot['DN4000_MODEL'], specphot['DN4000_MODEL_IVAR']

                dn4000_sigma = 1. / np.sqrt(dn4000_ivar)
                if dn4000_model_ivar > TINY:
                    dn4000_model_sigma = 1. / np.sqrt(dn4000_model_ivar)
                else:
                    dn4000_model_sigma = 0.

                restwave = wave_out / (1. + redshift) # [Angstrom]
                flam2fnu = (1 + redshift) * restwave**2 / (C_LIGHT * 1e5) * 1e-3 * 1e23 / FLUXNORM # [erg/s/cm2/A-->mJy, rest]
                fnu_obs = data['coadd_flux'] * flam2fnu # [mJy]
                fnu = flux_out_nolines * flam2fnu # [mJy]

                fnu_model = continuum_out * flam2fnu
                #fnu_fullmodel = modelflux * flam2fnu

                fnu_ivar = data['coadd_ivar'] / flam2fnu**2
                fnu_sigma, fnu_mask = ivar2var(fnu_ivar, sigma=True)

                I = (restwave > 3835.) * (restwave < 4115.)
                J = (restwave > 3835.) * (restwave < 4115.) * fnu_mask

                fig, ax = plt.subplots(figsize=(7, 6))
                ax.fill_between(restwave[I], fnu_obs[I]-fnu_sigma[I], fnu_obs[I]+fnu_sigma[I], color='red',
                                alpha=0.5, label=f'Observed Dn(4000)={dn4000:.3f}'+r'$\pm$'+f'{dn4000_sigma:.3f}')
                ax.plot(restwave[I], fnu[I], alpha=0.7, color='k', label=f'Line-free Dn(4000)={dn4000:.3f}' + \
                        r'$\pm$'+f'{dn4000_sigma:.3f}')
                #ax.plot(restwave[I], fnu_fullmodel[I], color='k', label=f'Model Dn(4000)={dn4000_model:.3f}')
                if dn4000_model_sigma > 0.:
                    ax.plot(restwave[I], fnu_model[I], alpha=0.7, label=f'Model Dn(4000)={dn4000_model:.3f}' + \
                            r'$\pm$'+f'{dn4000_model_sigma:.3f}')
                else:
                    ax.plot(restwave[I], fnu_model[I], alpha=0.7, label=f'Model Dn(4000)={dn4000_model:.3f}')
                ylim = ax.get_ylim()
                ax.fill_between([3850, 3950], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                                color='lightgray', alpha=0.5)
                ax.fill_between([4000, 4100], [ylim[0], ylim[0]], [ylim[1], ylim[1]],
                                color='lightgray', alpha=0.5)
                ax.set_xlabel(r'Rest Wavelength ($\AA$)')
                ax.set_ylabel(r'$F_{\nu}$ (mJy)')
                #ax.set_ylabel(r'$F_{\nu}\ ({\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~{\rm Hz}^{-1})$')
                ax.legend(loc='upper left', fontsize=10)
                ax.set_title(f'Dn(4000): {data["uniqueid"]}')

                fig.tight_layout()
                fig.savefig(pngfile)#, bbox_inches='tight')
                plt.close()
                log.info(f'Wrote {pngfile}')

    log.debug(f'Emission-line fitting took {time.time()-tall:.2f} seconds.')

    if debug_plots:
        for name in fastfit.value.dtype.names:
            print(name, fastfit[name])
        print()
        for name in specphot.value.dtype.names:
            print(name, specphot[name])

    return spectra_out
