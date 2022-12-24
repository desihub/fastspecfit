"""
fastspecfit.emlines
===================

Methods and tools for fitting emission lines.

"""
import pdb # for debugging

import os, time
import numpy as np
import multiprocessing
import numba

import astropy.units as u
from astropy.table import Table, Column

from desispec.interpolation import resample_flux
from fastspecfit.util import trapz_rebin, C_LIGHT
from fastspecfit.continuum import ContinuumTools
from desiutil.log import get_logger, DEBUG
log = get_logger(DEBUG)

def read_emlines():
    """Read the set of emission lines of interest.

    """
    from pkg_resources import resource_filename
    
    linefile = resource_filename('fastspecfit', 'data/emlines.ecsv')    
    linetable = Table.read(linefile, format='ascii.ecsv', guess=False)
    
    return linetable    

#@numba.jit(nopython=True)
def _build_emline_model(log10wave, redshift, lineamps, linevshifts, linesigmas, 
                        linewaves, emlinewave, resolution_matrix, camerapix):
    """Given parameters, build the model emission-line spectrum.

    ToDo: can this be optimized using numba?

    """
    log10model = np.zeros_like(log10wave) # [erg/s/cm2/A, observed frame]

    # Cut to lines with non-zero amplitudes.
    I = lineamps > 0
    if np.count_nonzero(I) > 0:
        linevshifts = linevshifts[I]
        linesigmas = linesigmas[I]
        lineamps = lineamps[I]
        linewaves = linewaves[I]

        # line-width [log-10 Angstrom] and redshifted wavelength [log-10 Angstrom]
        log10sigmas = linesigmas / C_LIGHT / np.log(10) 
        linezwaves = np.log10(linewaves * (1.0 + redshift + linevshifts / C_LIGHT))

        for lineamp, linezwave, log10sigma in zip(lineamps, linezwaves, log10sigmas):
            J = np.abs(log10wave - linezwave) < (8 * log10sigma) # cut to pixels within +/-N-sigma
            if np.count_nonzero(J) > 0:
                #print(lineamp, 10**linezwave, 10**log10wave[J].min(), 10**log10wave[J].max())
                log10model[J] = log10model[J] + lineamp * np.exp(-0.5 * (log10wave[J]-linezwave)**2 / log10sigma**2)

    # Optionally split into cameras, resample, and convolve with the resolution
    # matrix.
    emlinemodel = []
    if camerapix is None:
        for icam, specwave in enumerate(emlinewave):
            _emlinemodel = trapz_rebin(10**log10wave, log10model, specwave)
            _emlinemomdel = resolution_matrix[icam].dot(_emlinemodel)
            emlinemodel.append(_emlinemodel)
        return emlinemodel
    else:
        for icam, campix in enumerate(camerapix):
            _emlinemodel = trapz_rebin(10**log10wave, log10model, emlinewave[campix[0]:campix[1]])
            _emlinemomdel = resolution_matrix[icam].dot(_emlinemodel)
            emlinemodel.append(_emlinemodel)
        return np.hstack(emlinemodel)

def _objective_function(free_parameters, emlinewave, emlineflux, weights, redshift, 
                        log10wave, resolution_matrix, camerapix, parameters, Ifree, 
                        Itied, tiedtoparam, tiedfactor, doubletindx, doubletpair, 
                        linewaves):
    """The parameters array should only contain free (not tied or fixed) parameters."""

    # Parameters have to be allowed to exceed their bounds in the optimization
    # function, otherwise they get stuck at the boundary.

    # Handle tied parameters and bounds. Only check bounds on the free
    # parameters.
    #for I, (value, bnd) in enumerate(zip(free_parameters, bounds)):
    #    if value < bnd[0]:
    #        free_parameters[I] = bnd[0]
    #    if value > bnd[1]:
    #        free_parameters[I] = bnd[1]

    parameters[Ifree] = free_parameters

    if len(Itied) > 0:
        for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
            parameters[I] = parameters[indx] * factor

    # Take care of the doublets.
    lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line    
    lineamps[doubletindx] = lineamps[doubletindx] * lineamps[doubletpair]

    # Build the emission-line model.
    emlinemodel = _build_emline_model(log10wave, redshift, lineamps, linevshifts, 
                                      linesigmas, linewaves, emlinewave, 
                                      resolution_matrix, camerapix)

    if weights is None:
        residuals = emlinemodel - emlineflux
    else:
        residuals = weights * (emlinemodel - emlineflux)

    return residuals

class EMLineFit(ContinuumTools):
    """Class to fit an emission-line spectrum.

    """
    def __init__(self, chi2_default=0.0, minwave=3000.0, maxwave=10000.0, 
                 maxiter=5000, accuracy=1e-2, nolegend=False, 
                 ssptemplates=None, mapdir=None):
        """Class to model an emission-line spectrum.

        Parameters
        ----------
        chi2_default : :class:`float`, optional, defaults to 0.0.
            Default chi2 value for a emission line not fitted.
        minwave : :class:`float`, optional, defaults to 3000 A.
            Minimum observed-frame wavelength, which is used internally for the
            forward modeling of the emission-line spectrum.
        maxwave : :class:`float`, optional, defaults to 3000 A.
            Like `minwave` but the maximum observed-frame wavelength.
        maxiter : :class:`int`, optional, defaults to 5000.
            Maximum number of iterations.
        accuracy : :class:`float`, optional, defaults to 0.01.
            Fitting accuracy.
        nolegend : :class:`bool`, optional, defaults to `False`.
            Do not render the legend on the QA output.
        ssptemplates : :class:`str`, optional
            Full path to the SSP templates used for continuum-fitting.
        mapdir : :class:`str`, optional
            Full path to the Milky Way dust maps.

        Notes
        -----
        Need to document all the attributes.
        
        """
        super(EMLineFit, self).__init__(mapdir=mapdir, ssptemplates=ssptemplates)
        
        self.chi2_default = chi2_default
        self.maxiter = maxiter
        self.accuracy = accuracy
        self.nolegend = nolegend

        self.emwave_pixkms = 5.0                                  # pixel size for internal wavelength array [km/s]
        self.dlogwave = self.emwave_pixkms / C_LIGHT / np.log(10) # pixel size [log-lambda]
        self.log10wave = np.arange(np.log10(minwave), np.log10(maxwave), self.dlogwave)

        # default line-sigma for computing upper limits
        self.limitsigma_narrow = 75.0
        self.limitsigma_broad = 1200.0 
        self.wavepad = 5.0 # Angstrom

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

        self.doubletindx = np.hstack([np.where(self.param_names == doublet)[0] for doublet in doublet_names])
        self.doubletpair = np.hstack([np.where(self.param_names == pair)[0] for pair in doublet_pairs])

    def init_output(self, nobj=1):
        """Initialize the output data table.

        """
        import astropy.units as u
        
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

        for line in self.linetable['name']:
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

    def build_linemodels(self, redshift, wavelims=[3000, 10000], verbose=False):
        """Build all the multi-parameter emission-line models we will use.
    
        """
        def _print_linemodel(linemodel):
            for linename in linenames:
                for param in ['amp', 'sigma', 'vshift']:
                    I = np.where(param_names == linename+'_'+param)[0]
                    if len(I) == 1:
                        I = I[0]
                        if linemodel['tiedtoparam'][I] == -1:
                            if linemodel['fixed'][I]:
                                print('{:25s} is FIXED'.format(linename+'_'+param))
                        else:
                            if linemodel['fixed'][I]:
                                print('{:25s} tied to {:25s} with factor {:.4f} and FIXED'.format(
                                    linename+'_'+param, param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))
                            else:
                                print('{:25s} tied to {:25s} with factor {:.4f}'.format(
                                    linename+'_'+param, param_names[linemodel['tiedtoparam'][I]], linemodel['tiedfactor'][I]))

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
                    if np.all(linemodel[linemodel['tiedtoparam'] == tied]['fixed']):
                        if verbose:
                            print('Fixing {} because all tied lines are fixed: {}'.format(
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
        vmaxshift_narrow = 300.0
        vmaxshift_broad = 2500.0 # 3000.0
    
        minsigma_narrow = 1.0
        maxsigma_narrow = 500.0

        minsigma_broad = 1.0
        maxsigma_broad = 1e4

        minsigma_balmer_broad = minsigma_narrow
        maxsigma_balmer_broad = maxsigma_broad
    
        # Be very careful about changing the default broad line-sigma. Smaller
        # values like 1500 km/s (which is arguably more sensible) can lead to
        # low-amplitude broad lines in a bunch of normal star-forming galaxy
        # spectra. (They act to "suck up" local continuum variations.) Also
        # recall that if it's well-measured, we use the initial line-sigma in
        # continuum.ContinuumFit.estimate_linesigma, which is a better initial
        # guess.
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

        bounds_mgii_doublet = [0.01, 10.0] 
        bounds_oii_doublet = [0.5, 1.5] # [0.66, 1.4]
        bounds_sii_doublet = [0.5, 1.5] # [0.67, 1.2]
    
        # Create a new line-fitting table which contains the redshift-dependent
        # quantities for this object.
        fit_linetable = Table()
        fit_linetable['name'] = self.linetable['name']
        fit_linetable['isbalmer'] = self.linetable['isbalmer']
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
        final_linemodel = Table()
        final_linemodel['param_name'] = param_names
        final_linemodel['index'] = np.arange(nparam).astype(np.int32)
        final_linemodel['linename'] = np.tile(linenames, 3) # 3 parameters per line
        final_linemodel['tiedfactor'] = np.zeros(nparam, 'f8')
        final_linemodel['tiedtoparam'] = np.zeros(nparam, np.int16)-1
        final_linemodel['doubletpair'] = np.zeros(nparam, np.int16)-1
        final_linemodel['fixed'] = np.zeros(nparam, bool)
        final_linemodel['bounds'] = np.zeros((nparam, 2), 'f8')
        final_linemodel['initial'] = np.zeros(nparam, 'f8')
        final_linemodel['value'] = np.zeros(nparam, 'f8')

        final_linemodel['doubletpair'][self.doubletindx] = self.doubletpair

        # Build the relationship of "tied" parameters. In the 'tied' array, the
        # non-zero value is the multiplicative factor by which the parameter
        # represented in the 'tiedtoparam' index should be multiplied.
    
        # Physical doublets and lines in the same ionization species should have
        # their velocity shifts and line-widths always tied. In addition, set fixed
        # doublet-ratios here. Note that these constraints must be set on *all*
        # lines, not just those in range.
    
        for iline, linename in enumerate(linenames):
            # initial values and bounds - broad He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp_balmer_broad, maxamp_balmer_broad], 
                                                   [minsigma_balmer_broad, maxsigma_balmer_broad],
                                                   [-vmaxshift_broad, +vmaxshift_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - narrow He+Balmer lines
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - broad UV/QSO lines (non-Balmer)
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline]:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_broad, maxsigma_broad],
                                                   [-vmaxshift_broad, +vmaxshift_broad]],
                                                  [initamp, initsigma_broad, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # initial values and bounds - forbidden lines
            if fit_linetable['isbalmer'][iline] == False and fit_linetable['isbroad'][iline] == False:
                for param, bounds, default in zip(['amp', 'sigma', 'vshift'],
                                                  [[minamp, maxamp], [minsigma_narrow, maxsigma_narrow],
                                                   [-vmaxshift_narrow, +vmaxshift_narrow]],
                                                  [initamp, initsigma_narrow, initvshift]):
                    final_linemodel['initial'][param_names == linename+'_'+param] = default
                    final_linemodel['bounds'][param_names == linename+'_'+param] = bounds

            # tie parameters

            # broad He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]
            # narrow He + Balmer
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False and linename != 'halpha':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
            # other lines
            if linename == 'mgii_2796':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'mgii_2803_'+param)[0]
            if linename == 'nev_3346' or linename == 'neiii_3869': # should [NeIII] 3869 be tied to [NeV]???
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nev_3426_'+param)[0]
            if linename == 'oii_3726':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_3729_'+param)[0]
            # Tentative! Tie auroral lines to [OIII] 4363 but maybe we shouldn't tie [OI] 6300 here...
            if linename == 'nii_5755' or linename == 'oi_6300' or linename == 'siii_6312':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_4363_'+param)[0]
            if linename == 'oiii_4959':
                """
                [O3] (4-->2): airwave: 4958.9097 vacwave: 4960.2937 emissivity: 1.172e-21
                [O3] (4-->3): airwave: 5006.8417 vacwave: 5008.2383 emissivity: 3.497e-21
                """
                final_linemodel['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9839 # 2.8875
                final_linemodel['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oiii_5007_amp')[0]
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]
            if linename == 'nii_6548':
                """
                [N2] (4-->2): airwave: 6548.0488 vacwave: 6549.8578 emissivity: 2.02198e-21
                [N2] (4-->3): airwave: 6583.4511 vacwave: 6585.2696 emissivity: 5.94901e-21
                """
                final_linemodel['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 2.9421 # 2.936
                final_linemodel['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'nii_6584_amp')[0]
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'nii_6584_'+param)[0]
            if linename == 'sii_6731':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'sii_6716_'+param)[0]
            if linename == 'oii_7320':
                """
                [O2] (5-->2): airwave: 7318.9185 vacwave: 7320.9350 emissivity: 8.18137e-24
                [O2] (4-->2): airwave: 7319.9849 vacwave: 7322.0018 emissivity: 2.40519e-23
                [O2] (5-->3): airwave: 7329.6613 vacwave: 7331.6807 emissivity: 1.35614e-23
                [O2] (4-->3): airwave: 7330.7308 vacwave: 7332.7506 emissivity: 1.27488e-23
                """
                final_linemodel['tiedfactor'][param_names == linename+'_amp'] = 1.0 / 1.2251
                final_linemodel['tiedtoparam'][param_names == linename+'_amp'] = np.where(param_names == 'oii_7330_amp')[0]
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oii_7330_'+param)[0]
            if linename == 'siii_9069':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'siii_9532_'+param)[0]
            # Tentative! Tie SiIII] 1892 to CIII] 1908 because they're so close in wavelength.
            if linename == 'siliii_1892':
                for param in ['sigma', 'vshift']:
                    final_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'ciii_1908_'+param)[0]
            
        # Finally set the initial values and bounds on the doublet ratio parameters.
        for param, bounds, default in zip(['mgii_doublet_ratio', 'oii_doublet_ratio', 'sii_doublet_ratio'],
                                          [bounds_mgii_doublet, bounds_oii_doublet, bounds_sii_doublet],
                                          [init_mgii_doublet, init_oii_doublet, init_sii_doublet]):
            final_linemodel['initial'][final_linemodel['param_name'] == param] = default
            final_linemodel['bounds'][final_linemodel['param_name'] == param] = bounds
                    
        # Assign fixed=True to parameters which are outside the wavelength range
        # except those that are tied to other lines.
        _fix_parameters(final_linemodel, verbose=False)

        assert(np.all(final_linemodel['tiedtoparam'][final_linemodel['tiedfactor'] != 0] != -1))
        assert(len(final_linemodel[np.sum(final_linemodel['bounds'] == [0.0, 0.0], axis=1) > 0]) == 0)
    
        #_print_linemodel(final_linemodel)
        #final_linemodel[np.logical_and(final_linemodel['fixed'] == False, final_linemodel['tiedtoparam'] == -1)]

        # Model 2 - like final_linemodel, but broad lines have been fixed at
        # zero.
        final_linemodel_nobroad = final_linemodel.copy()
        final_linemodel_nobroad['fixed'] = False # reset

        for iline, linename in enumerate(linenames):
            if linename == 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    final_linemodel_nobroad['fixed'][param_names == linename+'_'+param] = True

            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    final_linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    final_linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

        #final_linemodel_nobroad[np.logical_and(final_linemodel_nobroad['fixed'] == False, final_linemodel_nobroad['tiedtoparam'] == -1)]

        _fix_parameters(final_linemodel_nobroad, verbose=False)

        assert(np.all(final_linemodel_nobroad['tiedtoparam'][final_linemodel_nobroad['tiedfactor'] != 0] != -1))

        # Model 3 - like final_linemodel, but with all the narrow and forbidden
        # lines tied together and all the broad lines tied together.
        initial_linemodel = final_linemodel.copy()
        initial_linemodel['fixed'] = False # reset

        for iline, linename in enumerate(linenames):
            # Tie all forbidden lines and narrow Balmer & He lines to [OIII] 5007.
            if fit_linetable['isbroad'][iline] == False and linename != 'oiii_5007':
                for param in ['sigma', 'vshift']:
                    initial_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    initial_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'oiii_5007_'+param)[0]

            ## Tie all narrow Balmer+He lines to narrow Halpha.
            #if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] == False and linename != 'halpha':
            #    for param in ['sigma', 'vshift']:
            #        initial_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
            #        initial_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_'+param)[0]
    
            # Tie all broad Balmer+He lines to broad Halpha.
            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['sigma', 'vshift']:
                    initial_linemodel['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    initial_linemodel['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

        _fix_parameters(initial_linemodel, verbose=False)

        assert(np.all(initial_linemodel['tiedtoparam'][initial_linemodel['tiedfactor'] != 0] != -1))

        # Model 4 - like initial_linemodel, but broad lines have been fixed at
        # zero.
        initial_linemodel_nobroad = initial_linemodel.copy()
        initial_linemodel_nobroad['fixed'] = False # reset        

        for iline, linename in enumerate(linenames):
            if linename == 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    initial_linemodel_nobroad['fixed'][param_names == linename+'_'+param] = True

            if fit_linetable['isbalmer'][iline] and fit_linetable['isbroad'][iline] and linename != 'halpha_broad':
                for param in ['amp', 'sigma', 'vshift']:
                    initial_linemodel_nobroad['tiedfactor'][param_names == linename+'_'+param] = 1.0
                    initial_linemodel_nobroad['tiedtoparam'][param_names == linename+'_'+param] = np.where(param_names == 'halpha_broad_'+param)[0]

        _fix_parameters(initial_linemodel_nobroad, verbose=False)

        assert(np.all(initial_linemodel_nobroad['tiedtoparam'][initial_linemodel_nobroad['tiedfactor'] != 0] != -1))

        if verbose:
            _print_linemodel(initial_linemodel_nobroad)

        return final_linemodel, final_linemodel_nobroad, initial_linemodel, initial_linemodel_nobroad

    def _initial_guesses_and_bounds(self, data, emlinewave, emlineflux):
        """For all lines in the wavelength range of the data, get a good initial guess
        on the amplitudes and line-widths. This step is critical for cases like,
        e.g., 39633354915582193 (tile 80613, petal 05), which has strong narrow
        lines.

        """
        initial_guesses, param_bounds = {}, {}
        #init_amplitudes, init_sigmas = {}, {}
    
        coadd_emlineflux = np.interp(data['coadd_wave'], emlinewave, emlineflux)
    
        for linename, linepix in zip(data['coadd_linename'], data['coadd_linepix']):
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
                    #bounds = [None, mx]
                    bounds = [-1e2, mx]
            else:
                bounds = [-np.min(np.abs(coadd_emlineflux[linepix])), mx]

            initial_guesses[linename+'_amp'] = amp
            param_bounds[linename+'_amp'] = bounds
    
        # Now update the linewidth but here we need to loop over *all* lines
        # (not just those in range). E.g., if H-alpha is out of range we need to
        # set its initial value correctly since other lines are tied to it
        # (e.g., main-bright-32406-39628257196245904).
        for iline in self.linetable:
            linename = iline['name']
            if iline['isbroad']:
                if iline['isbalmer']:
                    if data['linesigma_balmer_snr'] > 0: # broad Balmer lines
                        initial_guesses[linename+'_sigma'] = data['linesigma_balmer']
                else:
                    if data['linesigma_uv_snr'] > 0: # broad UV/QSO lines
                        initial_guesses[linename+'_sigma'] = data['linesigma_uv']
            else:
                # prefer narrow over Balmer
                if data['linesigma_narrow_snr'] > 0:
                    initial_guesses[linename+'_sigma'] = data['linesigma_narrow']
    
        return initial_guesses, param_bounds

    def _linemodel_to_parameters(self, linemodel):
        """Convert a linemodel model to a list of emission-line parameters."""

        linesplit = np.array_split(linemodel['index'], 3) # 3 parameters per line
        #linesplit = (np.arange(3) + 1) * len(linemodel) // 3 # 3 parameters per line
        lineamps = linemodel['value'][linesplit[0]].data
        linevshifts = linemodel['value'][linesplit[1]].data
        linesigmas = linemodel['value'][linesplit[2]].data

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        doublet = np.where(linemodel['doubletpair'] != -1)[0]
        lineamps[doublet] *= linemodel['value'][linemodel['doubletpair'][doublet]]

        parameters = np.hstack((lineamps, linevshifts, linesigmas))

        linewaves = self.fit_linetable['restwave'].data
        #lineinrange = self.fit_linetable['inrange'].data

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
    def _populate_linemodel(linemodel, initial_guesses, param_bounds):
        """Population an input linemodel with initial guesses and parameter bounds,
        taking into account fixed parameters.

        """
        # Set initial values and bounds.
        for iparam, param in enumerate(linemodel['param_name']):
            if param in initial_guesses.keys():
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0 # always set fixed parameter to zero
                else:
                    linemodel['initial'][iparam] = initial_guesses[param]
                    if param in param_bounds.keys():
                        linemodel['bounds'][iparam] = param_bounds[param]
            else:
                if linemodel['fixed'][iparam]:
                    linemodel['initial'][iparam] = 0.0

        # Now loop back through and ensure that tied relationships are enforced.
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        if len(Itied) > 0:
            for iparam, param in enumerate(linemodel['param_name'][Itied]):
                tieindx = linemodel['tiedtoparam'][Itied[iparam]]
                tiefactor = linemodel['tiedfactor'][Itied[iparam]]
                #log.info('{} tied to {} with factor {:.4f}'.format(param, linemodel[tieindx]['param_name'], tiefactor))
                linemodel['initial'][Itied[iparam]] = linemodel[tieindx]['initial'] * tiefactor

        linemodel['value'] = linemodel['initial'] # copy

    @staticmethod
    def objective_function(args):
        return _objective_function(*args)

    def _optimize(self, linemodel, emlinewave, emlineflux, weights, 
                  redshift, resolution_matrix, camerapix):
        """Wrapper to call the least-squares minimization given a linemodel.

        """
        from scipy.optimize import least_squares

        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, doubletpair, \
                     linewaves) = self._linemodel_to_parameters(linemodel)
        log.debug('Optimizing {} free parameters'.format(len(Ifree)))

        farg = (emlinewave, emlineflux, weights, redshift, self.log10wave, 
                resolution_matrix, camerapix, parameters, ) + \
                (Ifree, Itied, tiedtoparam, tiedfactor, doubletindx, 
                 doubletpair, linewaves)

        fit_info = least_squares(_objective_function, parameters[Ifree],
                                 args=farg, max_nfev=self.maxiter, 
                                 xtol=self.accuracy, method='lm')

        parameters[Ifree] = fit_info.x

        # Conditions for dropping a parameter (all parameters, not just those
        # being fitted):
        # --negative amplitude or sigma
        # --parameter at its default value (fit failed, right??)
        # --parameter outside its bounds
        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line

        drop1 = np.hstack((lineamps < 0, np.zeros(len(linevshifts), bool), linesigmas < 0))

        drop2 = np.zeros(len(parameters), bool)
        drop2[Ifree] = parameters[Ifree] == linemodel['value'][Ifree] # want 'value' here not 'initial'

        drop3 = np.zeros(len(parameters), bool)
        drop3[Ifree] = np.logical_or(parameters[Ifree] < linemodel['bounds'][Ifree, 0], 
                                     parameters[Ifree] > linemodel['bounds'][Ifree, 1])

        log.debug('Dropping {} negative amplitudes or line-widths.'.format(np.sum(drop1)))
        log.debug('Dropping {} parameters which were not optimized.'.format(np.sum(drop2)))
        log.debug('Dropping {} parameters which are out-of-bounds.'.format(np.sum(drop3)))
        Idrop = np.where(np.logical_or.reduce((drop1, drop2, drop3)))[0]

        if len(Idrop) > 0:
            log.debug('  Dropping {} unique parameters.'.format(len(Idrop)))
            parameters[Idrop] = 0.0

        # apply tied constraints
        if len(Itied) > 0:
            for I, indx, factor in zip(Itied, tiedtoparam, tiedfactor):
                parameters[I] = parameters[indx] * factor

        # Now loop back through and drop Broad balmer lines that:
        #   (1) are narrower than their narrow-line counterparts;
        #   (2) have a narrow line whose amplitude is smaller than that of the broad line
        #      --> Deprecated! main-dark-32303-39628176678192981 is an example
        #          of an object where there's a broad H-alpha line but no other
        #          forbidden lines!
        
        #IB = self.linetable['isbalmer'] * self.linetable['isbroad']
        #for linename in self.linetable['name'][IB]:
        #    amp_broad = getattr(bestfit, '{}_amp'.format(linename))
        #    if amp_broad.fixed: # line is either out of range or a suppressed line--skip it here
        #        continue
        #    amp = getattr(bestfit, '{}_amp'.format(linename.replace('_broad', ''))) # fragile
        #    sigma = getattr(bestfit, '{}_sigma'.format(linename.replace('_broad', ''))) # fragile
        #    sigma_broad = getattr(bestfit, '{}_sigma'.format(linename))
        #    
        #    if sigma_broad <= sigma:
        #        log.debug('Dropping {} (sigma_broad={:.2f}, sigma_narrow={:.2f})'.format(
        #            linename, sigma_broad.value, sigma.value))
        #        #vshift = getattr(bestfit, '{}_vshift'.format(linename))
        #        setattr(bestfit, '{}_amp'.format(linename), 0.0)
        #        setattr(bestfit, '{}_sigma'.format(linename), 0.0) # sigma_broad.default)
        #        setattr(bestfit, '{}_vshift'.format(linename), 0.0) # vshift.default)

        out_linemodel = linemodel.copy()
        out_linemodel['value'] = parameters
        out_linemodel.meta['nfev'] = fit_info['nfev']

        if False:
            bestfit = self.bestfit(out_linemodel, redshift, emlinewave, resolution_matrix, camerapix)
            import matplotlib.pyplot as plt
            plt.clf()
            plt.plot(emlinewave, emlineflux)
            plt.plot(emlinewave, bestfit)
            #plt.xlim(5800, 6200)
            #plt.xlim(6600, 6950)
            plt.xlim(5050, 5120)
            #plt.xlim(8850, 9050)
            plt.savefig('junk.png')
        
        return out_linemodel

    def chi2(self, linemodel, emlinewave, emlineflux, emlineivar, emlinemodel,
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
            chi2 = self.chi2_default
            
        if return_dof:
            return chi2, dof
        else:
            return chi2

    def bestfit(self, linemodel, redshift, emlinewave, resolution_matrix, camerapix):
        """Construct the best-fitting emission-line spectrum from a linemodel."""

        parameters, (Ifree, Itied, tiedtoparam, tiedfactor, bounds, doubletindx, \
                     doubletpair, linewaves) = self._linemodel_to_parameters(linemodel)

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line

        emlinemodel = _build_emline_model(self.log10wave, redshift, lineamps, 
                                          linevshifts, linesigmas, linewaves, 
                                          emlinewave, resolution_matrix,
                                          camerapix)

        return emlinemodel

    def emlinemodel_bestfit(self, specwave, specres, fastspecfit_table, redshift=None):
        """Wrapper function to get the best-fitting emission-line model
        from an fastspecfit table (used for QA and elsewhere).

        """
        if redshift is None:
            redshift = fastspecfit_table['CONTINUUM_Z']
        
        linewaves = self.linetable['restwave'].data

        parameters = [fastspecfit_table[param.upper()][0] for param in self.param_names]

        lineamps, linevshifts, linesigmas = np.array_split(parameters, 3) # 3 parameters per line    

        # Handle the doublets. Note we are implicitly assuming that the
        # amplitude parameters are always in the first third of parameters.
        lineamps[self.doubletindx] *= lineamps[self.doubletpair]

        emlinemodel = _build_emline_model(self.log10wave, redshift, lineamps, 
                                          linevshifts, linesigmas, linewaves, 
                                          specwave, specres, None)

        return emlinemodel

    def fit(self, data, continuummodel, smooth_continuum, synthphot=True,
            verbose=False, broadlinefit=True):
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
        from scipy.stats import sigmaclip
        from fastspecfit.util import ivar2var
            
        # Combine all three cameras; we will unpack them to build the
        # best-fitting model (per-camera) below.
        redshift = data['zredrock']
        emlinewave = np.hstack(data['wave'])
        oemlineivar = np.hstack(data['ivar'])
        specflux = np.hstack(data['flux'])
        resolution_matrix = data['res']
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
        wavelims = (np.min(emlinewave)+5, np.max(emlinewave)-5)

        # Build all the emission-line models for this object.
        final_linemodel, final_linemodel_nobroad, initial_linemodel, initial_linemodel_nobroad = \
            self.build_linemodels(redshift, wavelims=wavelims, verbose=False)

        # Get initial guesses on the parameters and populate the two "initial"
        # linemodels; the "final" linemodels will be initialized with the
        # best-fitting parameters from the initial round of fitting.
        initial_guesses, param_bounds = self._initial_guesses_and_bounds(data, emlinewave, emlineflux)

        for linemodel in [initial_linemodel, initial_linemodel_nobroad]:
            self._populate_linemodel(linemodel, initial_guesses, param_bounds)

        # Initial fit - initial_linemodel_nobroad
        t0 = time.time()
        initfit = self._optimize(initial_linemodel_nobroad, emlinewave, emlineflux, 
                                 weights, redshift, resolution_matrix, camerapix)
        initmodel = self.bestfit(initfit, redshift, emlinewave, resolution_matrix, camerapix)
        initchi2 = self.chi2(initfit, emlinewave, emlineflux, emlineivar, initmodel)
        nfree = np.sum((initfit['fixed'] == False) * (initfit['tiedtoparam'] == -1))
        log.info('Initial line-fitting with {} free parameters took {:.2f} sec (niter={}) with rchi2={:.4f}.'.format(
            nfree, time.time()-t0, initfit.meta['nfev'], initchi2))

        # Now try adding bround Balmer and helium lines and see if we improve
        # the chi2. First, do we have enough pixels around Halpha and Hbeta to
        # do this test?
        broadlinepix = []
        for icam in np.arange(len(data['cameras'])):
            pixoffset = int(np.sum(data['npixpercamera'][:icam]))
            for linename, linepix in zip(data['linename'][icam], data['linepix'][icam]):
                if linename == 'halpha_broad' or linename == 'hbeta_broad' or linename == 'hgamma_broad':
                    broadlinepix.append(linepix + pixoffset)
                    #print(data['wave'][icam][linepix])

        # Require minimum XX pixels.
        if broadlinefit or (len(broadlinepix) > 0 and len(np.hstack(broadlinepix)) > 10): 
            broadlinepix = np.hstack(broadlinepix)

            t0 = time.time()
            broadfit = self._optimize(initial_linemodel, emlinewave, emlineflux, weights, 
                                      redshift, resolution_matrix, camerapix)
            broadmodel = self.bestfit(broadfit, redshift, emlinewave, resolution_matrix, camerapix)
            broadchi2 = self.chi2(broadfit, emlinewave, emlineflux, emlineivar, broadmodel)
            nfree = np.sum((broadfit['fixed'] == False) * (broadfit['tiedtoparam'] == -1))
            log.info('Second (broad) line-fitting with {} free parameters took {:.2f} sec (niter={}) with rchi2={:.4f}'.format(
                nfree, time.time()-t0, broadfit.meta['nfev'], broadchi2))

            # Compare chi2 just in and around the broad lines. Need to account
            # for the different number of degrees of freedom!
            linechi2_init = np.sum(emlineivar[broadlinepix] * (emlineflux[broadlinepix] - initmodel[broadlinepix])**2) / len(broadlinepix)
            linechi2_broad = np.sum(emlineivar[broadlinepix] * (emlineflux[broadlinepix] - broadmodel[broadlinepix])**2) / len(broadlinepix)
            log.info('Chi2 with broad lines = {:.5f} and without broad lines = {:.5f} (delta-chi2={:.5f})'.format(
                linechi2_broad, linechi2_init, linechi2_init - linechi2_broad))

            if linechi2_broad > linechi2_init:
                bestfit = initfit
            else:
                bestfit = broadfit
        else:
            if broadlinefit:
                log.info('Too few pixels centered on candidate broad emission lines.')
            else:
                log.info('Skipping broad-line fitting.')
            bestfit = initfit
            linechi2_init = 0.0
            linechi2_broad = 0.0

        # Finally, one more fitting loop with all the line-constraints relaxed
        # but starting from the previous best-fitting values.
        if linechi2_broad > linechi2_init:
            linemodel = final_linemodel_nobroad
        else:
            linemodel = final_linemode

        linemodel['value'] = bestfit['value']
        Itied = np.where((linemodel['tiedtoparam'] != -1) * (linemodel['fixed'] == False))[0]
        if len(Itied) > 0:
            for I, indx, factor in zip(Itied, linemodel['tiedtoparam'][Itied], linemodel['tiedfactor'][Itied]):
                linemodel[I]['value'] = linemodel[indx]['value'] * factor

        t0 = time.time()
        finalfit = self._optimize(linemodel, emlinewave, emlineflux, weights, 
                                  redshift, resolution_matrix, camerapix)
        finalmodel = self.bestfit(finalfit, redshift, emlinewave, resolution_matrix, camerapix)
        finalchi2 = self.chi2(finalfit, emlinewave, emlineflux, emlineivar, finalmodel)
        nfree = np.sum((finalfit['fixed'] == False) * (finalfit['tiedtoparam'] == -1))
        log.info('Final line-fitting with {} free parameters took {:.2f} sec (niter={}) with rchi2={:.4f}.'.format(
            nfree, time.time()-t0, finalfit.meta['nfev'], finalchi2))

        #I = np.where(['oii_' in param for param in self.param_names])[0]

        # Initialize the output table; see init_fastspecfit for the data model.
        result = self.init_output()

        result['RCHI2'] = finalchi2
        result['LINERCHI2_BROAD'] = linechi2_broad
        result['DELTA_LINERCHI2'] = linechi2_init - linechi2_broad

        # Now fill the output table.
        for param in finalfit:
            # If the parameter was not optimized, set its value to zero.
            val = param['value']
            # special case the tied doublets
            if param['param_name'] == 'oii_doublet_ratio':
                result['OII_DOUBLET_RATIO'] = val
                result['OII_3726_AMP'] = val * finalfit[finalfit['param_name'] == 'oii_3729_amp']['value'] # * u.erg/(u.second*u.cm**2*u.Angstrom)
            elif param['param_name'] == 'sii_doublet_ratio':
                result['SII_DOUBLET_RATIO'] = val
                result['SII_6731_AMP'] = val * finalfit[finalfit['param_name'] == 'sii_6716_amp']['value'] # * u.erg/(u.second*u.cm**2*u.Angstrom)
            elif param['param_name'] == 'mgii_doublet_ratio':
                result['MGII_DOUBLET_RATIO'] = val
                result['MGII_2796_AMP'] = val * finalfit[finalfit['param_name'] == 'mgii_2803_amp']['value'] # * u.erg/(u.second*u.cm**2*u.Angstrom)
            else:
                result[param['param_name'].upper()] = val

        # Synthesize photometry from the best-fitting model (continuum+emission lines).
        if synthphot:
            if data['photsys'] == 'S':
                filters = self.decam
            else:
                filters = self.bassmzls

            # The wavelengths overlap between the cameras a bit...
            srt = np.argsort(emlinewave)
            padflux, padwave = filters.pad_spectrum((continuummodelflux+smoothcontinuummodelflux+finalmodel)[srt],
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

        specflux_nolines = specflux - finalmodel

        ## measure DN(4000) without the emission lines
        #if False:
        #    dn4000_nolines, _ = self.get_dn4000(emlinewave, specflux_nolines, redshift=redshift)
        #    result['DN4000_NOLINES'] = dn4000_nolines

        # get continuum fluxes, EWs, and upper limits
        narrow_sigmas, broad_sigmas, uv_sigmas = [], [], []
        narrow_redshifts, broad_redshifts, uv_redshifts = [], [], []
        for oneline in self.fit_linetable[self.fit_linetable['inrange']]:

            linename = oneline['name'].upper()
            linez = redshift + result['{}_VSHIFT'.format(linename)][0] / C_LIGHT
            linezwave = oneline['restwave'] * (1 + linez)
            linesigma = result['{}_SIGMA'.format(linename)][0] # [km/s]

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
            #log10sigma = linesigma / C_LIGHT / np.log(10)     # line-width [log-10 Angstrom]            

            # Are the pixels based on the original inverse spectrum fully
            # masked? If so, set everything to zero and move onto the next line.
            lineindx = np.where((emlinewave >= (linezwave - 2.*linesigma_ang)) *
                                (emlinewave <= (linezwave + 2.*linesigma_ang)))[0]
            
            if len(lineindx) > 0 and np.sum(oemlineivar[lineindx] == 0) / len(lineindx) > 0.3: # use original ivar
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
                        good = oemlineivar[lineindx] > 0 # use the original ivar
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
                        chi2 = np.sum(emlineivar[lineindx]*(emlineflux[lineindx]-finalmodel[lineindx])**2) / dof
    
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
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
                indxhi = np.where((emlinewave < (linezwave + 10*linesigma * linezwave / C_LIGHT)) *
                                  (emlinewave > (linezwave + 3.*linesigma * linezwave / C_LIGHT)) *
                                  (oemlineivar > 0))[0]
                                  #(finalmodel == 0))[0]
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
                    factor = 1 / ((1 + redshift) * result['{}_CONT'.format(linename)]) # --> rest frame
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
            title = 'Survey/Program/HealPix: {}/{}/{}, TargetID: {}'.format(
                    metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID'])
            pngfile = os.path.join(outdir, '{}-{}-{}-{}-{}.png'.format(
                    outprefix, metadata['SURVEY'], metadata['PROGRAM'], metadata['HEALPIX'], metadata['TARGETID']))
            #title = 'Tile: {}, TargetID/Fiber: {}/{}'.format(
            #        metadata['TILEID'], metadata['TARGETID'], metadata['FIBER'])
            #pngfile = os.path.join(outdir, '{}-{}-{}.png'.format(
            #        outprefix, metadata['TILEID'], metadata['TARGETID']))
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
            print('HACK ON LINEMASK!!!')
            _smooth_continuum, _ = self.smooth_continuum(np.hstack(data['wave']), np.hstack(residuals),
                                                         np.hstack(data['ivar']), redshift=redshift,
                                                         linemask=None)#np.hstack(data['linemask']))
        smooth_continuum = []
        for campix in data['camerapix']:
            smooth_continuum.append(_smooth_continuum[campix[0]:campix[1]])

        _emlinemodel = self.emlinemodel_bestfit(data['wave'], data['res'], Table(fastspec))

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

