"""
fastspecfit.linemasker
======================

Tools for pre-fitting and masking emission lines.

"""
import time
import numpy as np
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.util import C_LIGHT, quantile


class LineMasker(object):
    """Compute spectral line mask for continuum estimation,
       and to get initial guesses at emission line parameters

    """
    def __init__(self, emline_table):
        """
        Parameters
        ----------
        emline_table
            Emission line table
        """

        self.emline_table = emline_table


    @staticmethod
    def linepix_and_contpix(wave, ivar, linetable, linesigmas, linevshifts=None,
                            patchMap=None, redshift=0., nsigma=4., minlinesigma=50.,
                            mincontpix=11, get_contpix=True):
        """Support routine to determine the pixels potentially containing emission lines
        and the corresponding (adjacent) continuum.

        linesigmas in km/s
        minlinesigma in kms - minimum line-sigma for the purposes of masking
        mincontpix - minimum number of continuum pixels per line

        """
        def _get_linepix(zlinewave, sigma):
            # line-emission
            I = (wave > (zlinewave - nsigma * sigma)) * (wave < (zlinewave + nsigma * sigma)) * (ivar > 0.)
            return np.where(I)[0]

        def _get_contpix(zlinewaves, sigmas, nsigma_factor=2., linemask=None, lya=False):
            # never use continuum pixels blueward of Lyman-alpha
            if lya:
                minwave = np.min(zlinewaves)
            else:
                minwave = np.min(zlinewaves - nsigma_factor * nsigma * sigmas)
            maxwave = np.max(zlinewaves + nsigma_factor * nsigma * sigmas)
            if linemask is None:
                J = (wave > minwave) * (wave < maxwave) * (ivar > 0.)
            else:
                J = (wave > minwave) * (wave < maxwave) * (ivar > 0.) * ~linemask
            return np.where(J)[0]

        if linevshifts is None:
            linevshifts = np.zeros_like(linesigmas)

        if len(linevshifts) != len(linesigmas):
            errmsg = 'Mismatch in linevshifts and linesigmas dimensions.'
            log.critical(errmsg)
            raise ValueError(errmsg)

        linemask = np.zeros_like(wave, bool) # True - affected by possible emission line
        linenames = linetable['name'].value

        zlinewaves = linetable['restwave'].value * (1. + redshift + linevshifts / C_LIGHT)
        linesigmas[(linesigmas > 0.) * (linesigmas < minlinesigma)] = minlinesigma
        linesigmas_ang = linesigmas * zlinewaves / C_LIGHT # [km/s --> Angstrom]

        pix = {'linepix': {}, 'contpix': {}}
        if patchMap is not None:
            pix.update({'patch_contpix': {}, 'dropped': [], 'merged': []})
            patchids = list(patchMap.keys())
            npatch = len(patchids)

        # Initial set of pixels that may contain emission lines.
        for linename, zlinewave, sigma in zip(linenames, zlinewaves, linesigmas_ang):
            # skip fixed (e.g., hbeta_broad) lines
            if sigma <= 0.:
                continue 
            I = _get_linepix(zlinewave, sigma)
            if len(I) > 0:
                linemask[I] = True
                pix['linepix'][linename] = I

        # skip contpix
        if not get_contpix:
            return pix

        if patchMap is None:
            for linename, zlinewave, sigma in zip(linenames, zlinewaves, linesigmas_ang):
                # skip fixed (e.g., hbeta_broad) or fully masked lines
                if sigma <= 0. or not linename in pix['linepix'].keys():
                    continue
                lya = linename == 'lyalpha'
                J = _get_contpix(zlinewave, sigma, nsigma_factor=2., linemask=linemask, lya=lya)
                # go further out
                if len(J) < mincontpix:
                    J = _get_contpix(zlinewave, sigma, nsigma_factor=2.5, linemask=linemask, lya=lya)
                # drop the linemask_ condition; dangerous??
                if len(J) < mincontpix:
                    #log.debug(f'Dropping linemask condition for {linename} with width ' + \
                        #          f'{nsigma*sigma:.3f}-{3*nsigma*sigma:.3f} Angstrom')

                    # Note the smaller nsigma_factor (to stay closer to the
                    # line); remove the pixels already assigned to this
                    # line.
                    J = _get_contpix(zlinewave, sigma, nsigma_factor=2., linemask=None, lya=lya)
                    J = np.delete(J, np.isin(J, pix['linepix'][linename]))

                if len(J) > 0:
                    pix['contpix'][linename] = J
        else:
            # Now get the corresponding continuum pixels in "patches."
            for patchid in patchids:
                # Not all patchlines will be in the 'linepix' dictionary
                # because, e.g., the broad Balmer lines have linesigma=0 when
                # fitting the narrow-only linemodel. An entire patch can also
                # get dropped if a large portion of the spectrum is fully masked
                # (ivar==0).
                patchlines = patchMap[patchid][0]
                keep = np.isin(patchlines, list(pix['linepix'].keys()))
                if np.count_nonzero(keep) == 0:
                    pix['dropped'].append(patchid)
                    log.debug(f'Dropping patch {patchid} ({len(patchlines)} lines fully masked).')
                    continue

                patchlines = patchlines[keep]
                I = patchMap[patchid][1][keep]

                zlinewaves_patch = zlinewaves[I]
                sigmas_patch = linesigmas_ang[I]
                lya = 'lyalpha' in patchlines

                J = _get_contpix(zlinewaves_patch, sigmas_patch, nsigma_factor=2.,
                                 linemask=linemask, lya=lya)

                # go further out
                if len(J) < mincontpix:
                    J = _get_contpix(zlinewaves_patch, sigmas_patch, nsigma_factor=2.5,
                                     linemask=linemask, lya=lya)

                if len(J) > 0:
                    # all lines in this patch get the same continuum indices
                    mn, mx = np.min(J), np.max(J)
                    for patchline in patchlines:
                        pix['contpix'][patchline] = J

                        # Make sure the left/right edges of the patch include all
                        # the emission lines on this patch.
                        mn = np.min((mn, np.min(pix['linepix'][patchline])))
                        mx = np.max((mx, np.max(pix['linepix'][patchline])))
                    J = np.unique(np.hstack((mn, J, mx)))

                pix['patch_contpix'][patchid] = J # updated vector

            # Loop back through and merge patches that overlap by at least
            # mincontpix.
            patchkeys = pix['patch_contpix'].keys()
            for ipatch, patchid in enumerate(patchids[:npatch-1]):
                if patchid in patchkeys and patchids[ipatch+1] in patchkeys:
                    left = pix['patch_contpix'][patchid]
                    rght = pix['patch_contpix'][patchids[ipatch+1]]
                    incommon = np.intersect1d(left, rght, assume_unique=True)
                    if len(incommon) > mincontpix:
                        log.debug(f'Merging patches {patchid} and {patchids[ipatch+1]}')
                        newcontpix = np.union1d(left, rght)
                        newpatchid = patchids[ipatch] + patchids[ipatch+1]
                        del pix['patch_contpix'][patchids[ipatch]]
                        del pix['patch_contpix'][patchids[ipatch+1]]
                        pix['patch_contpix'][newpatchid] = newcontpix
                        pix['merged'].append(newpatchid)

            if patchMap is not None:
                if len(pix['dropped']) > 0:
                    pix['dropped'] = np.hstack(pix['dropped'])
                if len(pix['merged']) > 0:
                    pix['merged'] = np.hstack(pix['merged'])

        # make sure we haven't mucked up our indexing.
        if not pix['linepix'].keys() == pix['contpix'].keys():
            errmsg = 'Mismatch in linepix and contpix!'
            log.critical(errmsg)
            raise ValueError(errmsg)

        return pix


    def build_linemask(self, wave, flux, ivar, resolution_matrix, redshift=0.,
                       uniqueid=0, initsigma_broad=3000., initsigma_narrow=150.,
                       initsigma_balmer_broad=1000., initvshift_broad=0.,
                       initvshift_narrow=0., initvshift_balmer_broad=0.,
                       minsnr_balmer_broad=1.5, niter=2, debug_plots=False):
        """Generate a mask which identifies pixels impacted by emission lines.

        Parameters
        ----------
        wave : :class:`numpy.ndarray` [npix]
            Observed-frame wavelength array.
        flux : :class:`numpy.ndarray` [npix]
            Spectrum corresponding to `wave`.
        ivar : :class:`numpy.ndarray` [npix]
            Inverse variance spectrum corresponding to `flux`.
        resolution_matrix : :class:`list`
            List of :class:`fastspecfit.emline_fit.sparse_rep.ResMatrix`
            objects, one per camera.
        redshift : :class:`float`
            Object redshift.
        uniqueid : :class:`int`
            Unique identification number.
        initsigma_broad : :class:`float`
            Initial guess of the broad-line emission-line width in km/s.
        initsigma_narrow : :class:`float`
            Like `initsigma_broad` but for the forbidden lines.
        initsigma_balmer_broad : :class:`float`
            Like `initsigma_broad` but for the broad Balmer lines.
        initvshift_broad : :class:`float`
            Initial guess of the broad-line velocity shift in km/s.
        initvshift_narrow : :class:`float`
            Like `initvshift_broad` but for the forbidden lines.
        initvshift_balmer_broad : :class:`float`
            Like `initvshift_broad` but for the broad Balmer lines.
        minsnr_balmer_broad : :class:`float`
            Minimum signal-to-noise ratio for identifying a significant broad
            Balmer line.
        niter : :class:`int`
            Number of fit-in-patches iterations.
        debug_plots : :class:`bool`
            Optionally generate and write out various debugging plots.

        Returns
        -------
        :class:`dict` with line-masking arrays

        """
        from astropy.table import vstack
        from fastspecfit.emlines import EMFitTools, ParamType

        def _make_patchTable(patchids):
            """Initialize the table containing information on each patch."""
            patchids = np.atleast_1d(patchids)
            npatch = len(patchids)
            continuum_patches = Table()
            continuum_patches['patchid'] = patchids
            continuum_patches['endpts'] = np.zeros((npatch,2), int) # starting index relative to coadd_wave
            continuum_patches['pivotwave'] = np.zeros(npatch)
            continuum_patches['slope'] = np.zeros(npatch)
            continuum_patches['intercept'] = np.zeros(npatch)
            continuum_patches['slope_bounds'] = np.broadcast_to([-1e2, +1e2], (npatch, 2))
            continuum_patches['intercept_bounds'] = np.broadcast_to([-1e5, +1e5], (npatch, 2))
            continuum_patches['balmerbroad'] = np.zeros(npatch, bool) # does this patch have a broad Balmer line?
            return continuum_patches


        def _fit_patches(continuum_patches, patchMap, linemodel, debug_plots=False,
                         testBalmerBroad=False, minsnr=1.5, modelname='',
                         png=None):
            """Iteratively fit all the lines in patches."""

            linesigmas = np.zeros(nline)
            linesigmas[EMFit.isBroad] = initsigma_broad
            linesigmas[EMFit.isNarrow] = initsigma_narrow
            if testBalmerBroad:
                linesigmas[EMFit.isBalmerBroad] = initsigma_balmer_broad

            linevshifts = np.zeros_like(linesigmas)
            linevshifts[EMFit.isBroad] = initvshift_broad
            linevshifts[EMFit.isNarrow] = initvshift_narrow
            if testBalmerBroad:
                linevshifts[EMFit.isBalmerBroad] = initvshift_balmer_broad

            initial_guesses = None

            t0 = time.time()
            for iiter in range(niter):

                # Build the line and continuum masks (only for lines in range).
                pix = self.linepix_and_contpix(wave, ivar, linetable_inrange,
                                               linesigmas[EMFit.line_in_range],
                                               linevshifts=linevshifts[EMFit.line_in_range],
                                               patchMap=patchMap, redshift=redshift)

                # Check for fully dropped patches, which can happen if large
                # parts of the spectrum are masked.
                if len(pix['dropped']) > 0:
                    for patchid in np.atleast_1d(pix['dropped']):
                        drop = np.where(continuum_patches['patchid'] == patchid)[0][0]
                        continuum_patches.remove_row(drop)
                        patchMap.pop(patchid)

                # In the case that patches have been merged, update the
                # continuum_patches table and patchMap dictionary.
                if len(pix['merged']) > 0:
                    for newpatchid in np.atleast_1d(pix['merged']):
                        oldpatchids = list(newpatchid)
                        O = np.where(np.isin(continuum_patches['patchid'].value, oldpatchids))[0]

                        # update continuum_patches
                        new_continuum_patch = _make_patchTable(newpatchid)
                        new_continuum_patch['pivotwave'] = np.mean(continuum_patches[O]['pivotwave'])
                        new_continuum_patch['balmerbroad'] = np.any(continuum_patches[O]['balmerbroad'])
                        continuum_patches.remove_rows(O)
                        continuum_patches = vstack((continuum_patches, new_continuum_patch))

                        # update patchMap
                        newlines, newI, newJ = [], [], []
                        for oldpatchid in oldpatchids:
                            newlines.append(patchMap[oldpatchid][0])
                            newI.append(patchMap[oldpatchid][1])
                            newJ.append(patchMap[oldpatchid][2])
                            del patchMap[oldpatchid]
                        patchMap[newpatchid] = (np.hstack(newlines), np.hstack(newI), np.hstack(newJ))


                linemask = np.zeros(len(wave), bool) # False=masked

                # Determine the edges of each patch based on the continuum
                # (line-free) pixels of all the lines on that patch.
                for ipatch, patchid in enumerate(pix['patch_contpix'].keys()):
                    contindx = pix['patch_contpix'][patchid]
                    continuum_patches['endpts'][ipatch] = (np.min(contindx), np.max(contindx)) # should be sorted...
                    # initial guesses and bounds, but check for pathological distributions
                    lo, med, hi = quantile(flux[contindx], (0.05, 0.5, 0.95))
                    if lo < med and lo < hi:
                        continuum_patches['intercept'][ipatch] = med
                        continuum_patches['intercept_bounds'][ipatch] = [lo, hi]

                    # unmask the continuum patches
                    linemask[contindx] = True # True=unmasked

                # unmask the lines
                for line in pix['linepix'].keys():
                    linemask[pix['linepix'][line]] = True

                # only fit the pixels in the patches
                weights = np.sqrt(ivar * linemask)

                # Get initial guesses on the line-emission on the first iteration.
                if initial_guesses is None:
                    initial_guesses, param_bounds = EMFit._initial_guesses_and_bounds(
                        pix['linepix'], flux, contpix=pix['contpix'],
                        subtract_local_continuum=True,
                        initial_linesigma_broad=initsigma_broad,
                        initial_linesigma_narrow=initsigma_narrow,
                        initial_linesigma_balmer_broad=initsigma_balmer_broad,
                        initial_linevshift_broad=0.,
                        initial_linevshift_narrow=0.,
                        initial_linevshift_balmer_broad=0.,
                    )

                # fit!
                linefit, contfit = EMFit.optimize(linemodel, initial_guesses,
                                                  param_bounds, wave,
                                                  flux, weights, redshift,
                                                  resolution_matrix, camerapix,
                                                  continuum_patches=continuum_patches)

                # Update the initial guesses as well as linesigmas and
                # linevshifts (for linepix_and_contpix, at the top of the
                # iteration loop).
                if iiter < niter-1:
                    initial_guesses = linefit['value'].value.copy()
                    linevshifts = initial_guesses[EMFit.line_table['params'][:, ParamType.VSHIFT]]
                    linesigmas = initial_guesses[EMFit.line_table['params'][:, ParamType.SIGMA]]


            # Build the best-fitting model and estimate the S/N of each line.
            parameters = linefit['value'].value.copy()
            parameters[EMFit.doublet_idx] *= parameters[EMFit.doublet_src]

            lineamps = parameters[EMFit.line_table['params'][:, ParamType.AMPLITUDE]]
            linevshifts = parameters[EMFit.line_table['params'][:, ParamType.VSHIFT]]
            linesigmas = parameters[EMFit.line_table['params'][:, ParamType.SIGMA]]

            bestfit = EMFit.bestfit(linefit, redshift, wave, resolution_matrix,
                                    camerapix, continuum_patches=contfit)
            residuals = flux - bestfit

            linesnrs = np.zeros_like(lineamps)
            noises = np.zeros(len(contfit))
            for ipatch, (patchid, endpts) in enumerate(contfit.iterrows('patchid', 'endpts')):
                s, e = endpts
                noise = np.ptp(quantile(residuals[s:e], (0.25, 0.75))) / 1.349 # robust sigma
                noises[ipatch] = noise
                lineindx = patchMap[patchid][2] # index back to full line_table
                if noise != 0:
                    linesnrs[lineindx] = lineamps[lineindx] / noise

            # Derive the final linesigmas and linevshifts, and the maximum S/N
            # of each type of line. If the line isn't well-measured, drop the
            # S/N condition.
            strong = linesnrs > minsnr
            Ifree = linefit[EMFit.line_table['params'][:, ParamType.SIGMA]]['free']
            isBroad = EMFit.isBroad * Ifree * strong
            isNarrow = EMFit.isNarrow * Ifree * strong
            isBalmerBroad = EMFit.isBalmerBroad * Ifree * strong

            if np.any(isBroad):
                linesigma_broad = np.atleast_1d(linesigmas[isBroad])[0] # all values should be the same
                linevshift_broad = np.atleast_1d(linevshifts[isBroad])[0]
                maxsnr_broad = np.max(linesnrs[isBroad])
            else:
                isBroad = EMFit.isBroad * Ifree
                if np.any(isBroad):
                    linesigma_broad = np.atleast_1d(linesigmas[isBroad])[0]
                    linevshift_broad = np.atleast_1d(linevshifts[isBroad])[0]
                    maxsnr_broad = np.max(linesnrs[isBroad])
                else:
                    linesigma_broad = initsigma_broad # 0.
                    linevshift_broad = initvshift_broad
                    maxsnr_broad = 0.

            if np.any(isNarrow):
                linesigma_narrow = np.atleast_1d(linesigmas[isNarrow])[0]
                linevshift_narrow = np.atleast_1d(linevshifts[isNarrow])[0]
                maxsnr_narrow = np.max(linesnrs[isNarrow])
            else:
                isNarrow = EMFit.isNarrow * Ifree
                if np.any(isNarrow):
                    linesigma_narrow = np.atleast_1d(linesigmas[isNarrow])[0] 
                    linevshift_narrow = np.atleast_1d(linevshifts[isNarrow])[0]
                    maxsnr_narrow = np.max(linesnrs[isNarrow])
                else:
                    linesigma_narrow = initsigma_narrow
                    linevshift_narrow = initvshift_narrow
                    maxsnr_narrow = 0.

            if np.any(isBalmerBroad):
                linesigma_balmer_broad = np.atleast_1d(linesigmas[isBalmerBroad])[0]
                linevshift_balmer_broad = np.atleast_1d(linevshifts[isBalmerBroad])[0]
                maxsnr_balmer_broad = np.max(linesnrs[isBalmerBroad])
            else:
                isBalmerBroad = EMFit.isBalmerBroad * Ifree
                if np.any(isBalmerBroad):
                    linesigma_balmer_broad = np.atleast_1d(linesigmas[isBalmerBroad])[0]
                    linevshift_balmer_broad = np.atleast_1d(linevshifts[isBalmerBroad])[0]
                    maxsnr_balmer_broad = np.max(linesnrs[isBalmerBroad])
                else:
                    # we do not want to overmask if a broad Balmer line isn't detected
                    linesigma_balmer_broad = 0. # initsigma_balmer_broad
                    linevshift_balmer_broad = initvshift_balmer_broad
                    maxsnr_balmer_broad = 0.

            final_linesigmas = (linesigma_broad, linesigma_narrow, linesigma_balmer_broad)
            final_linevshifts = (linevshift_broad, linevshift_narrow, linevshift_balmer_broad)
            maxsnrs = (maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad)

            log.debug(f'Broad: S/N={maxsnr_broad:.1f}, (sigma,dv)=({linesigma_broad:.0f},{linevshift_broad:.0f}) km/s; ' + \
                      f'Narrow: S/N={maxsnr_narrow:.1f}, ({linesigma_narrow:.0f},{linevshift_narrow:.0f}) km/s; '+ \
                      f'Balmer Broad: S/N={maxsnr_balmer_broad:.1f}, ({linesigma_balmer_broad:.0f},{linevshift_balmer_broad:.0f}) km/s.')
            log.debug(f'Fitting {modelname} with patches ({niter} iterations) took {time.time()-t0:.4f} seconds.')

            # optionally build a QA figure
            if debug_plots:
                import matplotlib.pyplot as plt
                import seaborn as sns
                from fastspecfit.emline_fit import EMLine_MultiLines

                sns.set(context='talk', style='ticks', font_scale=0.8)

                npatch = len(contfit)
                ncols = 3
                nrows = int(np.ceil(npatch / ncols))

                lines = EMLine_MultiLines(parameters, wave, redshift,
                                          linetable['restwave'].value,
                                          resolution_matrix, camerapix)

                def _niceline(line):
                    match line:
                        case 'lyalpha':
                            return r'S/N(Ly$\alpha$)='
                        case 'civ_1549':
                            return r'S/N(CIV$\lambda1549$)='
                        case 'ciii_1908':
                            return r'S/N(CIII]$\lambda1908$)='
                        case 'mgii_2796':
                            return r'S/N(MgII$\lambda2796$)='
                        case 'mgii_2803':
                            return r'S/N(MgII$\lambda2803$)='
                        case 'oii_3726':
                            return r'S/N([OII]$\lambda3726$)='
                        case 'oii_3729':
                            return r'S/N([OII]$\lambda3729$)='
                        case 'hgamma':
                            return r'S/N(H$\gamma$)='
                        case 'hgamma_broad':
                            return r'S/N(H$\gamma_{b}$)='
                        case 'hbeta':
                            return r'S/N(H$\beta$)='
                        case 'hbeta_broad':
                            return r'S/N(H$\beta_{b}$)='
                        case 'oiii_4959':
                            return r'S/N([OIII]$\lambda4959$)='
                        case 'oiii_5007':
                            return r'S/N([OIII]$\lambda5007$)='
                        case 'nii_6548':
                            return r'S/N([NII]$\lambda6548$)='
                        case 'halpha':
                            return r'S/N(H$\alpha$)='
                        case 'halpha_broad':
                            return r'S/N(H$\alpha_{b}$)='
                        case 'nii_6584':
                            return r'S/N([NII]$\lambda6584$)='
                        case 'sii_6716':
                            return r'S/N([SII]$\lambda6716$)='
                        case 'sii_6731':
                            return r'S/N([SII]$\lambda6731$)='


                fig, ax = plt.subplots(nrows, ncols, figsize=(5.5*ncols, 5.5*nrows))
                for ipatch, ((patchid, endpts, slope, intercept, pivotwave), xx) in enumerate(
                        zip(contfit.iterrows('patchid', 'endpts', 'slope', 'intercept', 'pivotwave'), ax.flat)):
                    # get the starting and ending pixels first
                    s, e = endpts

                    for iline in patchMap[patchid][2]:
                        (ls, le), _ = lines.getLine(iline)
                        if ls != le: # fixed / dropped lines
                            s = np.min((s, ls))
                            e = np.max((e, le))

                    xx.plot(wave[s:e] / 1e4, flux[s:e], color='gray')
                    xx.plot(wave[s:e] / 1e4, bestfit[s:e], color='k', ls='-', alpha=0.75)
                    cmodel = slope * (wave[s:e]-pivotwave) + intercept
                    xx.plot(wave[s:e] / 1e4, cmodel+noises[ipatch], color='gray', lw=1, ls='-')
                    xx.plot(wave[s:e] / 1e4, cmodel, color='k', lw=2, ls='--')
                    xx.plot(wave[s:e] / 1e4, cmodel-noises[ipatch], color='gray', lw=1, ls='-')
                    for line, iline in zip(patchMap[patchid][0], patchMap[patchid][2]):
                        (ls, le), profile = lines.getLine(iline)
                        if ls != le: # skip fixed lines
                            label = _niceline(line)+f'{linesnrs[iline]:.1f}; '+r'$\sigma$='+f'{linesigmas[iline]:.0f}'+' km/s'
                            xx.plot(wave[ls:le] / 1e4, profile, alpha=0.75, label=label)
                    if len(patchMap[patchid][0]) > 4:
                        nlegcol = 1
                        yfactor = 1.4
                    else:
                        nlegcol = 1
                        yfactor = 1.3
                    ymin = -1.2 * noises[ipatch]
                    ymax = yfactor * np.max((quantile(flux[s:e], 0.99), np.max(bestfit[s:e])))
                    xx.set_ylim(ymin, ymax)
                    xx.legend(loc='upper left', fontsize=8, ncols=nlegcol)
                    xx.set_title(f'Patch {patchid}')
                for rem in range(ipatch+1, ncols*nrowspatch):
                    ax.flat[rem].axis('off')

                if ax.ndim == 1:
                    ulpos = ax[0].get_position()
                    llpos = ax[0].get_position()
                    lrpos = ax[-1].get_position()
                    dxlabel = 0.08
                    bottom = 0.14
                else:
                    ulpos = ax[0, 0].get_position()
                    llpos = ax[-1, 0].get_position()
                    lrpos = ax[-1, -1].get_position()
                    dxlabel = 0.07
                    bottom = 0.11

                xpos = (lrpos.x1 - llpos.x0) / 2. + llpos.x0
                ypos = llpos.y0 - dxlabel
                fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                         ha='center', va='center')

                xpos = ulpos.x0 - 0.09
                ypos = (ulpos.y1 - llpos.y0) / 2. + llpos.y0
                fig.text(xpos, ypos, r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                         ha='center', va='center', rotation=90)

                fig.subplots_adjust(left=0.08, right=0.97, bottom=bottom, top=0.92, wspace=0.23, hspace=0.3)

                if png:
                    fig.savefig(png, bbox_inches='tight')
                    plt.close()

            return linefit, contfit, final_linesigmas, final_linevshifts, maxsnrs


        # main function begins here
        camerapix = np.array([[0, len(wave)]]) # one camera

        # Read just the strong lines and determine which lines are in range of the camera.
        EMFit = EMFitTools(emline_table=self.emline_table, uniqueid=uniqueid, stronglines=True)
        EMFit.compute_inrange_lines(redshift, wavelims=(np.min(wave), np.max(wave)))

        # Build the narrow and narrow+broad emission-line models.
        linemodel_broad, linemodel_nobroad = EMFit.build_linemodels(separate_oiii_fit=False)

        # ToDo: are there ever *no* "strong" lines in range?
        linetable = EMFit.line_table
        linetable_inrange = linetable[EMFit.line_in_range]
        nline = len(linetable)

        # Initialize the continuum_patches table for all patches in range.
        patchids = np.unique(linetable_inrange['patch'])
        continuum_patches = _make_patchTable(patchids)

        # Get the mapping between patchid and the set of lines belonging to each
        # patch, and pivotwave.
        patchMap = {}
        for ipatch, patchid in enumerate(patchids):
            I = np.where(linetable_inrange['patch'] == patchid)[0] # index into fitted / in-range lines
            J = np.where(linetable['patch'] == patchid)[0]         # index into *all* lines
            patchlines = linetable_inrange['name'][I].value
            patchMap[patchid] = (patchlines, I, J)
            linewaves = linetable_inrange['restwave'][I] * (1. + redshift)
            pivotwave = 0.5 * (np.min(linewaves) + np.max(linewaves)) # midpoint
            continuum_patches['pivotwave'][ipatch] = pivotwave
            # is there a broad Balmer line on this patch?
            continuum_patches['balmerbroad'][ipatch] = np.any(EMFit.isBalmerBroad_noHelium_Strong[EMFit.line_in_range][I])

        # Need to pass copies of continuum_patches and patchMap because they can
        # get modified dynamically by _fit_patches.
        linefit_nobroad, contfit_nobroad, linesigmas_nobroad, linevshifts_nobroad, maxsnrs_nobroad = \
            _fit_patches(continuum_patches.copy(), patchMap.copy(),
                         linemodel_nobroad, testBalmerBroad=False,
                         debug_plots=debug_plots, modelname='narrow lines only',
                         png=f'qa-patches-nobroad-{uniqueid}.png')

        # Only fit with broad Balmer lines if at least one patch contains a
        # broad line.
        B = contfit_nobroad['balmerbroad']
        if np.any(B):
            linefit_broad, contfit_broad, linesigmas_broad, linevshifts_broad, maxsnrs_broad = \
                _fit_patches(continuum_patches.copy(), patchMap.copy(),
                             linemodel_broad, testBalmerBroad=True,
                             debug_plots=debug_plots, modelname='narrow+broad lines',
                             png=f'qa-patches-broad-{uniqueid}.png')

            # if a broad Balmer line is well-detected, take its linewidth
            if maxsnrs_broad[2] > minsnr_balmer_broad:
                log.info(f'Adopting broad Balmer-line masking: S/N(broad Balmer) ' + \
                         f'{maxsnrs_broad[2]:.1f} > {minsnr_balmer_broad:.1f}')
                finalsigma_broad, finalsigma_narrow, finalsigma_balmer_broad = linesigmas_broad
                finalvshift_broad, finalvshift_narrow, finalvshift_balmer_broad = linevshifts_broad
                maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad = maxsnrs_broad
            else:
                log.info(f'Adopting narrow Balmer-line masking: S/N(broad Balmer) ' + \
                         f'{maxsnrs_broad[2]:.1f} < {minsnr_balmer_broad:.1f}.')
                finalsigma_broad, finalsigma_narrow, finalsigma_balmer_broad = linesigmas_nobroad
                finalvshift_broad, finalvshift_narrow, finalvshift_balmer_broad = linevshifts_nobroad
                maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad = maxsnrs_nobroad
        else:
            log.info(f'Adopting narrow Balmer-line masking: no Balmer lines in wavelength range.')
            finalsigma_broad, finalsigma_narrow, finalsigma_balmer_broad = linesigmas_nobroad
            finalvshift_broad, finalvshift_narrow, finalvshift_balmer_broad = linevshifts_nobroad
            maxsnr_broad, maxsnr_narrow, maxsnr_balmer_broad = maxsnrs_nobroad

        log.info(f'Masking line-widths: broad {finalsigma_broad:.0f} km/s narrow {finalsigma_narrow:.0f} km/s ' + \
                 f'broad Balmer {finalsigma_balmer_broad:.0f} km/s.')

        # Build the final pixel mask for *all* lines using our current best
        # knowledge of the broad Balmer lines....(comment continued below)
        EMFit = EMFitTools(emline_table=self.emline_table, uniqueid=uniqueid, stronglines=False)
        EMFit.compute_inrange_lines(redshift, wavelims=(np.min(wave), np.max(wave)))

        linesigmas = np.zeros(len(EMFit.line_table))
        linesigmas[EMFit.isBroad] = finalsigma_broad
        linesigmas[EMFit.isNarrow] = finalsigma_narrow
        linesigmas[EMFit.isBalmerBroad] = finalsigma_balmer_broad

        linevshifts = np.zeros_like(linesigmas)
        linevshifts[EMFit.isBroad] = finalvshift_broad
        linevshifts[EMFit.isNarrow] = finalvshift_narrow
        linevshifts[EMFit.isBalmerBroad] = finalvshift_balmer_broad

        pix = self.linepix_and_contpix(wave, ivar,
                                       EMFit.line_table[EMFit.line_in_range],
                                       linesigmas[EMFit.line_in_range],
                                       linevshifts=linevshifts[EMFit.line_in_range],
                                       patchMap=None, redshift=redshift)

        # Build another QA figure
        if debug_plots:
            import matplotlib.pyplot as plt
            import seaborn as sns

            sns.set(context='talk', style='ticks', font_scale=0.8)

            png = f'qa-linemask-{uniqueid}.png'

            linenames = list(pix['linepix'].keys())
            zlinewaves = EMFit.line_table[EMFit.line_in_range]['restwave'] * (1. + redshift)
            nline = len(linenames)

            ncols = 5
            nrows = int(np.ceil(nline / ncols))

            fig, ax = plt.subplots(nrows, ncols, figsize=(3*ncols, 2*nrows))

            for iline, (linename, xx) in enumerate(zip(linenames, ax.flat)):
                linepix = pix['linepix'][linename]
                contpix = pix['contpix'][linename]
                s = np.min((np.min(contpix), np.min(linepix)))
                e = np.max((np.max(contpix), np.max(linepix)))

                ylim = quantile(flux[s:e], (0.01, 0.99))

                xx.plot(wave[s:e] / 1e4, flux[s:e], color='gray', alpha=0.5)
                xx.scatter(wave[contpix] / 1e4, flux[contpix], color='blue', s=10, marker='s')
                xx.scatter(wave[linepix] / 1e4, flux[linepix], color='orange', s=10, marker='o', alpha=0.7)
                xx.axvline(x=zlinewaves[iline] / 1e4, color='k', ls='--', lw=2)
                xx.set_ylim(ylim[0], 1.3 * ylim[1])
                xx.text(0.1, 0.9, linename, ha='left', va='center', transform=xx.transAxes,
                        fontsize=8)

            for rem in range(iline+1, ncols*nrows):
                ax.flat[rem].axis('off')

            if ax.ndim == 1:
                ulpos = ax[0].get_position()
                urpos = ax[-1].get_position()
                llpos = ax[0].get_position()
                lrpos = ax[-1].get_position()
                top = 0.92
                bottom = 0.14
                dytitle = 0.13
                dyxlabel = 0.15
            else:
                ulpos = ax[0, 0].get_position()
                urpos = ax[0, -1].get_position()
                llpos = ax[-1, 0].get_position()
                lrpos = ax[-1, -1].get_position()
                top = 0.95
                bottom = 0.07
                dytitle = 0.09
                dyxlabel = 0.08

            xpos = (lrpos.x1 - llpos.x0) / 2. + llpos.x0
            ypos = llpos.y0 - dyxlabel
            fig.text(xpos, ypos, r'Observed-frame Wavelength ($\mu$m)',
                     ha='center', va='center')

            xpos = ulpos.x0 - 0.1
            ypos = (ulpos.y1 - llpos.y0) / 2. + llpos.y0# + 0.03
            fig.text(xpos, ypos, r'$F_{\lambda}\ (10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1})$',
                     ha='center', va='center', rotation=90)

            xpos = (urpos.x1 - ulpos.x0) / 2. + ulpos.x0
            ypos = ulpos.y1 + dytitle
            fig.text(xpos, ypos, f'LinePix/ContPix: {uniqueid}', ha='center', va='center')

            fig.subplots_adjust(left=0.06, right=0.97, bottom=bottom, top=top, wspace=0.23, hspace=0.3)
            fig.savefig(png, bbox_inches='tight')
            plt.close()

        # (comment continued from above) ...but reset the broad Balmer
        # line-width to a minimum value and make another linepix mask. We need
        # to do this so that emlines.emline_specfit has a chance to remeasure
        # the broad Balmer lines after continuum-subtraction.
        if finalsigma_balmer_broad < initsigma_balmer_broad:
            finalsigma_balmer_broad = initsigma_balmer_broad
            linesigmas[EMFit.isBalmerBroad] = finalsigma_balmer_broad
            linevshifts[EMFit.isBalmerBroad] = finalvshift_balmer_broad
            newpix = self.linepix_and_contpix(wave, ivar,
                                              EMFit.line_table[EMFit.line_in_range],
                                              linesigmas[EMFit.line_in_range],
                                              linevshifts=linevshifts[EMFit.line_in_range],
                                              patchMap=None, redshift=redshift)

            linepix = newpix['linepix']
        else:
            linepix = pix['linepix']

        out = {
            'linesigma_broad': finalsigma_broad,
            'linesigma_narrow': finalsigma_narrow,
            'linesigma_balmer_broad': finalsigma_balmer_broad, # updated value
            'linevshift_broad': finalvshift_broad,
            'linevshift_narrow': finalvshift_narrow,
            'linevshift_balmer_broad': finalvshift_balmer_broad, # updated value
            'maxsnr_broad': maxsnr_broad,
            'maxsnr_narrow': maxsnr_narrow,
            'maxsnr_balmer_broad': maxsnr_balmer_broad,
            'balmerbroad': np.any(contfit_nobroad['balmerbroad']), # True = one or more broad Balmer line in range
            'coadd_linepix': linepix,
        }

        return out

