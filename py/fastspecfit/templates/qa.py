"""
fastspecfit.templates.qa
========================

QA for templates

"""
import pdb

import os
import numpy as np
from astropy.table import Table
from scipy.ndimage import median_filter

from fastspecfit.util import ivar2var, C_LIGHT
from fastspecfit.templates.templates import rebuild_fastspec_spectrum

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

from desiutil.log import get_logger
log = get_logger()

def plot_style(font_scale=1.2):
    import seaborn as sns
    sns.set(context='talk', style='ticks', palette='deep', font_scale=font_scale)#, rc=rc)
    colors = sns.color_palette()
    return sns, colors

def qa_template_colors(phot, template_colors, ntspace=25, png=None):

    if ntspace == 1:
        prefix = 'All '
    else:
        prefix = ''
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    ax1.hexbin(phot['RMAG']-phot['ZMAG'], phot['GMAG']-phot['RMAG'], mincnt=1, bins='log',
               #C=cat['weight'], reduce_C_function=np.sum,
               cmap=cmap)
    ax1.set_xlabel(r'$(r - z)_{\rm obs}$')
    ax1.set_ylabel(r'$(g - r)_{\rm obs}$')
    ax1.set_xlim(rzobslim)
    ax1.set_ylim(grobslim)
    ax1.text(0.05, 0.9, 'Data', ha='left', va='bottom',
             transform=ax1.transAxes, fontsize=14)
    ax1.grid(True)

    #cb = fig.colorbar(hb, ax=ax1)
    #cb.set_label(r'log$_{10}$ (Number of Galaxies)')
    
    for tt in np.arange(0, nt, ntspace):
        ax2.plot(template_colors['rz'][tt, :], template_colors['gr'][tt, :], marker='s', 
                 markersize=5, ls='-', alpha=0.5)
        
    for tt in np.arange(0, nt, ntspace):
        ax2.scatter(template_colors['rz'][tt, 0], template_colors['gr'][tt, 0], marker='o', 
                   facecolors='none', s=40, edgecolors='k',
                   linewidth=1, zorder=10)
        
    ax2.text(0.1, 0.05, 'z=0.0', ha='left', va='bottom',
             transform=ax2.transAxes, fontsize=14)
    ax2.text(0.05, 0.9, '{}Models (z=0.0-1.5, dz=0.1)'.format(prefix), 
             ha='left', va='bottom',
             transform=ax2.transAxes, fontsize=14)
    
    ax2.set_xlim(rzobslim)
    ax2.set_ylim(grobslim)
    ax2.set_xlabel(r'$(r - z)_{\rm obs}$')
    ax2.set_ylabel(r'$(g - r)_{\rm obs}$')
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.tick_right()
    ax2.grid(True)
    
    ax3.hexbin(phot['ZMAG']-phot['W1MAG'], phot['RMAG']-phot['ZMAG'], mincnt=1, bins='log',
               #C=cat['weight'], reduce_C_function=np.sum,
               cmap=cmap)
    ax3.set_ylabel(r'$(r - z)_{\rm obs}$')
    ax3.set_xlabel(r'$(z - W1)_{\rm obs}$')
    ax3.set_ylim(rzobslim)
    ax3.set_xlim(zW1obslim)
    ax3.text(0.05, 0.9, 'Data', ha='left', va='bottom',
             transform=ax3.transAxes, fontsize=14)
    ax3.grid(True)
    
    for tt in np.arange(0, nt, ntspace):
        ax4.plot(template_colors['zW1'][tt, :], template_colors['rz'][tt, :], marker='s', 
                 markersize=5, ls='-', alpha=0.5)
        
    for tt in np.arange(0, nt, ntspace):
        ax4.scatter(template_colors['zW1'][tt, 0], template_colors['rz'][tt, 0], marker='o', 
                   facecolors='none', s=40, edgecolors='k',
                   linewidth=1, zorder=10)
        
    ax4.text(0.05, 0.3, 'z=0.0', ha='left', va='bottom',
             transform=ax4.transAxes, fontsize=14)
    ax4.text(0.05, 0.9, '{}Models (z=0.0-1.5, dz=0.1)'.format(prefix), 
             ha='left', va='bottom',
             transform=ax4.transAxes, fontsize=14)
    
    ax4.set_xlim(zW1obslim)
    ax4.set_ylim(rzobslim)
    ax4.set_ylabel(r'$(r - z)_{\rm obs}$')
    ax4.set_xlabel(r'$(z - W1)_{\rm obs}$')
    ax4.yaxis.set_label_position('right')
    ax4.yaxis.tick_right()
    ax4.grid(True)
    
    plt.subplots_adjust(wspace=0.05, hspace=0.28)
    
    if png:
        log.info('Writing {}'.format(png))
        fig.savefig(png)
        plt.close()

def qa_bpt(targetclass, fastspecfile=None, png=None):
    """QA of the fastspec emission-line spectra.

    """
    from fastspecfit.templates.templates import remove_undetected_lines, read_stacked_fastspec

    sns, _ = plot_style()

    fastmeta, _fastspec = read_stacked_fastspec(fastspecfile, read_spectra=False)
    fastspec = remove_undetected_lines(_fastspec)
    nobj = len(fastmeta)

    def oplot_class(ax, kewley=False, **kwargs):
        if kewley:
            niiha = np.linspace(-1.9, 0.4, 1000)
            oiiihb = 0.61 / (niiha-0.47) + 1.19
        else:
            niiha = np.linspace(-1.9, -0.1, 1000)
            oiiihb = 0.61 / (niiha-0.05) + 1.3
        ax.plot(niiha, oiiihb, **kwargs)

    def _bpt(cc, cclabel='Redshift', vmin=None, vmax=None, png=None):
        fig, ax = plt.subplots(figsize=(10, 7))
        cb = ax.scatter(niiha, oiiihb, c=cc, cmap='jet', vmin=vmin, vmax=vmax)
        oplot_class(ax, kewley=True, color='k', ls='--', lw=3, label='Kewley+01')
        oplot_class(ax, kewley=False, color='k', lw=3, label='Kauffmann+03')
        plt.colorbar(cb, label=cclabel)
        ax.set_xlim(-1.9, 0.7)
        ax.set_ylim(-1.2, 1.5)
        ax.set_xlabel(r'$\log_{10}$ ([NII] $\lambda6584$ / H$\alpha$)')
        ax.set_ylabel(r'$\log_{10}$ ([OIII] $\lambda5007$ / H$\beta$)')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax.legend(fontsize=16, loc='lower left')#, ncol=2)
        plt.subplots_adjust(bottom=0.15, left=0.18, top=0.95, right=0.95)
        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    good = np.where(
        (fastspec['HALPHA_FLUX'] > 0) *
        (fastspec['HBETA_FLUX'] > 0) *
        (fastspec['NII_6584_FLUX'] > 0) *
        (fastspec['OIII_5007_FLUX'] > 0)
        #(fastspec['HALPHA_CHI2'] < 1e4)
    )[0]

    niiha = np.log10(fastspec['NII_6584_FLUX'][good] / fastspec['HALPHA_FLUX'][good])
    oiiihb = np.log10(fastspec['OIII_5007_FLUX'][good] / fastspec['HBETA_FLUX'][good])
    ww = np.where((niiha > -0.05) * (niiha < 0.05) * (oiiihb < -0.5))[0]
    #log.info(fastspec[good][ww]['HALPHA_FLUX', 'NII_6584_FLUX'])

    zz = fastspec['CONTINUUM_Z'][good]
    ewhb = fastspec['HBETA_EW'][good]
    #rW1 = fastmeta['RW1'][good]
    #gr = fastmeta['GR'][good]

    _bpt(zz, 'Redshift', vmin=0, vmax=0.5, png=png.replace('.png', '-redshift.png'))
    _bpt(np.log10(ewhb), r'$\log_{10}\,\mathrm{EW}(\mathrm{H}\beta)$', 
         png=png.replace('.png', '-ewhb.png'))            
    #_bpt(rW1, r'$r-W1$', vmin=-0.3, vmax=0.9, png=png.replace('.png', '-rW1.png'))
    #_bpt(gi, r'$g-i$', vmin=0.6, vmax=1.3, png=png.replace('.png', '-gi.png'))

def qa_fastspec_fullspec(targetclass, fastwave=None, fastflux=None, fastivar=None,
                         fastmeta=None, fastspec=None, fastspecfile=None, CFit=None,
                         EMFit=None, ncol=3, nrow=5, pdffile=None):
    """Full-spectrum QA.
    
    """
    from fastspecfit.util import ivar2var, C_LIGHT
    from fastspecfit.templates.templates import rebuild_fastspec_spectrum, read_stacked_fastspec

    sns, _ = plot_style()        
    
    if CFit is None or EMFit is None:
        from fastspecfit.continuum import ContinuumFit
        from fastspecfit.emlines import EMLineFit    
        CFit = ContinuumFit()
        EMFit = EMLineFit()    

    if fastwave is None:
        fastwave, fastflux, fastivar, fastmeta, fastspec = read_stacked_fastspec(fastspecfile)
        #fastspec = remove_undetected_lines(fastspec, EMFit.linetable, devshift=False)

    if targetclass == 'lrg':
        absmagcol = 'MR'
        colorcol = 'RW1'
        absmaglabel = 'M_{{0.0r}}'
        colorlabel = '^{{0.0}}(r-W1)'
    elif targetclass == 'elg':
        absmagcol = 'MG'
        colorcol = 'GR'
        absmaglabel = 'M_{{0.0g}}'
        colorlabel = '^{{0.0}}(g-r)'
    elif targetclass == 'bgs':
        absmagcol = 'MR'
        colorcol = 'GR'
        absmaglabel = 'M_{{0.0r}}'
        colorlabel = '^{{0.0}}(g-r)'
    else:
        raise NotImplemented
        
    nobj = len(fastmeta)
    icam = 0
        
    zobj = np.unique(fastmeta['ZOBJ'])
    npage = len(zobj)

    inches_wide_perpanel = 4.0
    inches_tall_perpanel = 3.0

    if npage == 1:
        png = True
    else:
        png = False
    
    if pdffile:
        if png:
            pdffile = pdffile.replace('.pdf', '.png')
        else:
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages(pdffile)

    for ipage in np.arange(npage):
        log.info('Building page {}/{}'.format(ipage+1, npage))
        pageindx = np.where(zobj[ipage] == fastmeta['ZOBJ'])[0]

        absmag = sorted(set(fastmeta[absmagcol][pageindx])) # subpage
        nsubpage = len(absmag)

        for isubpage in np.arange(nsubpage):

            subpageindx = np.where((absmag[isubpage] == fastmeta[absmagcol][pageindx]))[0]

            fig, allax = plt.subplots(nrow, ncol, figsize=(inches_wide_perpanel*ncol, inches_tall_perpanel*nrow),
                                      sharex=True, sharey=True)
            for iplot, (indx, ax) in enumerate(zip(pageindx[subpageindx], allax.flatten())):
                #log.info(ipage, isubpage, iplot, len(pageindx), len(subpageindx))

                # rebuild the best-fitting spectrum
                modelwave, continuum, smooth_continuum, emlinemodel, data = rebuild_fastspec_spectrum(
                    fastspec[indx], fastwave, fastflux[indx, :], fastivar[indx, :], CFit, EMFit)

                #sigma, _ = ivar2var(data['ivar'][icam], sigma=True)
                #ax.fill_between(data['wave'][icam], data['flux'][icam]-sigma,
                #                            data['flux'][icam]+sigma, color='skyblue')
                ax.plot(data['wave'][icam], data['flux'][icam], color='skyblue')
                ax.plot(modelwave, continuum+emlinemodel, color='firebrick', alpha=0.5)
                ax.plot(modelwave, continuum, color='blue', alpha=0.5)
                #ax.plot(modelwave, continuum+smooth_continuum, color='gray', alpha=0.3)
                ax.plot(modelwave, smooth_continuum, color='gray', alpha=0.7)

                ymin, ymax = 1e6, -1e6

                filtflux = median_filter(data['flux'][icam], 51, mode='nearest')
                sigflux = np.std(data['flux'][icam][data['ivar'][icam] > 0])
                if -2 * sigflux < ymin:
                    ymin = -2 * sigflux
                if sigflux * 5 > ymax:
                    ymax = sigflux * 5
                if np.max(filtflux) > ymax:
                    ymax = np.max(filtflux) * 1.4

                ax.text(0.96, 0.06, r'${:.2f}<{}<{:.2f}$'.format(
                    fastmeta['{}MIN'.format(colorcol)][indx], colorlabel,
                    fastmeta['{}MAX'.format(colorcol)][indx]),
                    ha='right', va='bottom', transform=ax.transAxes, fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='gray', alpha=0.25))
                ax.text(0.04, 0.96, '\n'.join(( 'N={}, S/N={:.1f}'.format(
                    fastmeta['NOBJ'][indx], fastspec['CONTINUUM_SNR_ALL'][indx]), )),
                    ha='left', va='top', transform=ax.transAxes, fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='gray', alpha=0.25))

                ax.set_xlim(modelwave.min(), modelwave.max())
                ax.set_ylim(ymin, ymax)
                ax.set_xticklabels([])
                ax.set_yticklabels([])

                plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.07, right=0.95, top=0.95, bottom=0.1)

                if iplot == ncol*nrow-1:
                    break
            
            fig.text(0.52, 0.968, r'${:.2f}<z<{:.2f}\ {:.1f}<{}<{:.1f}$'.format(
                fastmeta['ZOBJMIN'][indx], fastmeta['ZOBJMAX'][indx],
                fastmeta['{}MIN'.format(absmagcol)][indx], absmaglabel,
                fastmeta['{}MAX'.format(absmagcol)][indx]),
                ha='center', va='center', fontsize=22)

            for rem in np.arange(ncol*nrow-iplot-1)+iplot+1:
                allax.flatten()[rem].axis('off')
                
            if pdffile and png is False:
                pdf.savefig(fig)
                plt.close()

    if pdffile:
        log.info('Writing {}'.format(pdffile))
        if png:
            fig.savefig(pdffile)
            plt.close()
        else:
            pdf.close()

def qa_fastspec_emlinespec(targetclass, fastwave=None, fastflux=None, fastivar=None,
                           fastmeta=None, fastspec=None, fastspecfile=None, CFit=None,
                           EMFit=None, ncol=3, nrow=5, pdffile=None):
    """QA of the fastspec emission-line spectra.

    """
    from matplotlib.colors import Normalize
    from fastspecfit.templates.templates import remove_undetected_lines
    from fastspecfit.util import ivar2var, C_LIGHT
    from fastspecfit.templates.templates import rebuild_fastspec_spectrum, read_stacked_fastspec

    sns, _ = plot_style()        

    if CFit is None or EMFit is None:
        from fastspecfit.continuum import ContinuumFit
        from fastspecfit.emlines import EMLineFit    
        CFit = ContinuumFit()
        EMFit = EMLineFit()    
    
    if fastwave is None:
        fastwave, fastflux, fastivar, fastmeta, fastspec = read_stacked_fastspec(fastspecfile)
        
    fastspec_fix = remove_undetected_lines(fastspec, EMFit.linetable, devshift=False)

    # plotting preferences
    cmap = plt.cm.get_cmap('jet')
    #cmap = sns.color_palette(as_cmap=True)
    cnorm = Normalize(vmin=np.min(fastmeta['ZOBJ']), vmax=np.max(fastmeta['ZOBJ']))

    inches_wide = 16
    inches_fullspec = 6
    inches_perline = inches_fullspec / 2.0
    nlinepanels = 4

    nline = len(set(EMFit.linetable['plotgroup']))
    nlinerows = np.ceil(nline / nlinepanels).astype(int)
    nrows = 1 + nlinerows

    height_ratios = np.hstack([1, [0.5]*nlinerows])

    plotsig_default = 150.0 # 300.0 # [km/s]
    meanwaves, deltawaves, sigmas, linenames = [], [], [], []
    for plotgroup in set(EMFit.linetable['plotgroup']):
        I = np.where(plotgroup == EMFit.linetable['plotgroup'])[0]
        linenames.append(EMFit.linetable['nicename'][I[0]])
        meanwaves.append(np.mean(EMFit.linetable['restwave'][I]))
        deltawaves.append((np.max(EMFit.linetable['restwave'][I]) - 
                           np.min(EMFit.linetable['restwave'][I])) / 2)
        sigmas.append(plotsig_default)
    srt = np.argsort(meanwaves)
    meanwaves = np.hstack(meanwaves)[srt]
    deltawaves = np.hstack(deltawaves)[srt]
    sigmas = np.hstack(sigmas)[srt]
    linenames = np.hstack(linenames)[srt]

    if targetclass == 'lrg':
        absmagcol = 'MR'
        colorcol = 'RW1'
        absmaglabel = 'M_{{0.0r}}'
        colorlabel = '^{{0.0}}(r-W1)'
    elif targetclass == 'elg':
        absmagcol = 'MG'
        colorcol = 'GR'
        absmaglabel = 'M_{{0.0g}}'
        colorlabel = '^{{0.0}}(g-r)'
    elif targetclass == 'bgs':
        absmagcol = 'MR'
        colorcol = 'GR'
        absmaglabel = 'M_{{0.0r}}'
        colorlabel = '^{{0.0}}(g-r)'
    else:
        raise NotImplemented
        
    # how many pages?
    nobj = len(fastmeta)
    icam = 0
    
    restcolor = np.unique(fastmeta[colorcol])
    npage = len(restcolor)

    if npage == 1:
        png = True
    else:
        png = False
    
    if pdffile:
        if png:
            pdffile = pdffile.replace('.pdf', '.png')
        else:
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages(pdffile)

    # make the plot!
    for ipage in np.arange(npage):
        log.info('Building page {}/{}'.format(ipage+1, npage))

        pageindx = np.where(restcolor[ipage] == fastmeta[colorcol])[0]
        absmag = sorted(set(fastmeta[absmagcol][pageindx])) # subpage
        nsubpage = len(absmag)

        for isubpage in np.arange(nsubpage):#[:1]:#[::2]:

            subpageindx = np.where((absmag[isubpage] == fastmeta[absmagcol][pageindx]))[0]

            fig = plt.figure(figsize=(inches_wide, 2*inches_fullspec + inches_perline*nlinerows))
            gs = fig.add_gridspec(nrows, nlinepanels, height_ratios=height_ratios)

            bigax = fig.add_subplot(gs[0, :])
            ax, irow, icol = [], 1, 0
            for iax in np.arange(nline):
                icol = iax % nlinepanels
                if iax > 0 and iax % nlinepanels == 0:
                    irow += 1
                xx = fig.add_subplot(gs[irow, icol])
                ax.append(xx)
            
            bigymin, bigymax = 1e6, -1e6
            lineymin, lineymax = np.zeros(nline)+1e6, np.zeros(nline)-1e6
            removelabels = np.ones(nline, bool)
        
            for iplot, indx in enumerate(pageindx[subpageindx]):
                #log.info(ipage, isubpage, iplot, len(pageindx), len(subpageindx))

                modelwave, continuum, smooth_continuum, emlinemodel, data = rebuild_fastspec_spectrum(
                    fastspec[indx], fastwave, fastflux[indx, :], fastivar[indx, :], CFit, EMFit)

                #if fastmeta['IBIN'][indx] == 1262:
                #    pdb.set_trace()

                redshift = data['zredrock']
                emlineflux = data['flux'][icam] - continuum - smooth_continuum

                modelwave /= (1+redshift) # rest-frame

                label = 'z=[{:.2f}-{:.2f}] (N={})'.format(
                    fastmeta['ZOBJMIN'][indx], fastmeta['ZOBJMAX'][indx],
                    np.sum(fastmeta['ZOBJ'][pageindx[subpageindx]] == fastmeta['ZOBJ'][indx]))
                #bigax.plot(modelwave/(1+redshift), emlineflux, color='gray')
                bigax.plot(modelwave, emlinemodel, label=label, color=cmap(cnorm(fastmeta['ZOBJ'][indx])))

                if -np.max(emlinemodel)*0.05 < bigymin:
                    bigymin = -np.max(emlinemodel)*0.05
                if np.max(emlinemodel)*1.1 > bigymax:
                    bigymax = np.max(emlinemodel)*1.1
                if np.max(emlinemodel) == 0.0:
                    bigymin, bigymax = 0.0, 1.0

                # zoom in on individual emission lines
                for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(
                    meanwaves, deltawaves, sigmas, linenames)):
                    wmin = (meanwave - deltawave) - 8 * sig * meanwave / C_LIGHT
                    wmax = (meanwave + deltawave) + 8 * sig * meanwave / C_LIGHT
                    lineindx = np.where((modelwave > wmin) * (modelwave < wmax))[0]

                    if len(lineindx) > 1:
                        if np.min(emlinemodel[lineindx]) > 0.0: # at least one line kept (snr>3)
                            removelabels[iax] = False
                            ax[iax].plot(modelwave[lineindx], emlinemodel[lineindx],
                                         color=cmap(cnorm(fastmeta['ZOBJ'][indx])))

                            if -np.max(emlinemodel[lineindx])*0.05 < lineymin[iax]:
                                lineymin[iax] = -np.max(emlinemodel[lineindx])*0.05
                            if np.max(emlinemodel[lineindx]) * 1.1 > lineymax[iax]:
                                lineymax[iax] = np.max(emlinemodel[lineindx]) * 1.1
                            if np.abs(lineymax[iax]-lineymin[iax]) < 1e-2:
                                removelabels[iax] = False

            for iax, xx in enumerate(ax):
                xx.text(0.08, 0.89, linenames[iax], ha='left', va='center',
                        transform=xx.transAxes, fontsize=20)
                if removelabels[iax]:
                    xx.set_ylim(0, 1)
                    xx.set_xticklabels([])
                    xx.set_yticklabels([])
                else:
                    if lineymax[iax] == lineymin[iax]:
                        lineymax[iax] = 1.0
                        
                    xx.set_ylim(lineymin[iax], lineymax[iax])
                    xlim = xx.get_xlim()
                    xx.xaxis.set_major_locator(ticker.MaxNLocator(2))

            # don't repeat the legend labels
            hand, lab = bigax.get_legend_handles_labels()
            ulabels = dict(zip(lab, hand))
            bigax.legend(ulabels.values(), ulabels.keys(), fontsize=18, loc='upper left')
            #bigax.legend(fontsize=18, loc='upper left')

            bigax.set_ylim(bigymin, bigymax)
            bigax.set_xlim(2600, 7200) # 3500, 9300)
            bigax.set_title(r'${:.2f}<{}<{:.2f}\ {:.1f}<{}<{:.1f}$'.format(
                fastmeta['{}MIN'.format(colorcol)][indx], colorlabel,
                fastmeta['{}MAX'.format(colorcol)][indx],
                fastmeta['{}MIN'.format(absmagcol)][indx], absmaglabel,
                fastmeta['{}MAX'.format(absmagcol)][indx]))
            #bigax.set_xlabel('Observed-frame Wavelength ($\AA$)')

            plt.subplots_adjust(wspace=0.28, left=0.07, right=0.95, top=0.95, bottom=0.1)

            if pdffile and png is False:
                pdf.savefig(fig)
                plt.close()

    if pdffile:
        log.info('Writing {}'.format(pdffile))
        if png:
            fig.savefig(pdffile)
            plt.close()
        else:
            pdf.close()

def qa_photometry(targetclass, samplefile=None, png_obs=None, png_rest=None, png_rest_bins=None):
    """QA of the observed- and rest-frame photometry.

    """
    from matplotlib.colors import LogNorm
    from fastspecfit.templates.sample import read_parent_sample, stacking_bins

    sns, _ = plot_style()
    cmap = plt.cm.get_cmap('RdYlBu')
    mincnt = 1

    phot, spec, meta = read_parent_sample(samplefile)
    bins, nbins = stacking_bins(targetclass, verbose=True)
    
    def bgs_obs(phot, png=None):
        robslim = (15, 21.0)
        grobslim = (-0.2, 2.5)
        rzobslim = (-0.5, 1.5)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

        ax1.hexbin(phot['RMAG']-phot['ZMAG'], phot['GMAG']-phot['RMAG'],
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum
                   extent=np.hstack((rzobslim, grobslim)))
        ax1.set_xlabel(r'$(r - z)_{\rm obs}$')
        ax1.set_ylabel(r'$(g - r)_{\rm obs}$')
        ax1.set_xlim(rzobslim)
        ax1.set_ylim(grobslim)

        hb = ax2.hexbin(phot['RMAG'], phot['GMAG']-phot['RMAG'],
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((robslim, grobslim)))
        ax2.set_xlabel(r'$r_{\rm obs}$')
        ax2.set_ylim(grobslim)
        ax2.set_xlim(robslim)

        cax = fig.add_axes([0.88, 0.12, 0.02, 0.83])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False)
        fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')

        for aa in (ax1, ax2):
            aa.grid(True)

        plt.subplots_adjust(left=0.12, top=0.95, right=0.85, bottom=0.19, wspace=0.07)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()
            
    def bgs_rest(phot, meta, bins=None, png=None):
        zlim = (0.0, 0.6)
        Mrlim = (-16, -25)
        grlim = (-0.2, 1.2)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

        ax1.hexbin(meta['Z'], phot['ABSMAG_R'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, Mrlim)))
        ax1.set_ylim(Mrlim)
        ax1.set_xlim(zlim)
        ax1.set_xlabel('Redshift')
        ax1.set_ylabel(r'$M_{0.0r}$')
        #ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))        

        if bins:
            dx, dy = bins['zobj']['del'], bins['Mr']['del']
            [ax1.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['Mr']['grid']]

        ax2.hexbin(meta['Z'], phot['ABSMAG_G']-phot['ABSMAG_R'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, grlim)))
        ax2.set_xlim(zlim)
        ax2.set_ylim(grlim)
        ax2.set_xlabel('Redshift')
        ax2.set_ylabel(r'$^{0.0}(g - r)$')#, labelpad=-10)
        #ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.2))        
        #ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.5))        

        if bins:
            dx, dy = bins['zobj']['del'], bins['gr']['del']
            [ax2.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['gr']['grid']]

        hb = ax3.hexbin(phot['ABSMAG_R'], phot['ABSMAG_G']-phot['ABSMAG_R'], 
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((Mrlim, grlim)))
        ax3.set_xlabel(r'$M_{0.0r}$')
        ax3.set_ylabel(r'$^{0.0}(g - r)$')#, labelpad=-10)
        ax3.set_xlim(Mrlim)
        ax3.set_ylim(grlim)
        #ax3.yaxis.set_major_locator(ticker.MultipleLocator(0.5))        

        if bins:
            dx, dy = bins['Mr']['del'], bins['gr']['del']
            [ax3.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mr']['grid'] for yy in bins['gr']['grid']]
            
        ax4.axis('off')

        cax = fig.add_axes([0.49, 0.12, 0.02, 0.36])
        #cax = fig.add_axes([0.54, 0.4, 0.35, 0.03])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False)
        fig.colorbar(hb, format=formatter, label='Number of Galaxies',
                     cax=cax)#, orientation='horizontal')

        for aa in (ax1, ax2, ax3):
            aa.grid(True)

        plt.subplots_adjust(left=0.1, top=0.95, wspace=0.3, hspace=0.3, right=0.88, bottom=0.13)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    def elg_obs(phot, png=None):
        gobslim = (19.5, 24.5)
        grobslim = (-1.2, 1.2)
        rzobslim = (-1.5, 2.2)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

        ax1.hexbin(phot['RMAG']-phot['ZMAG'], phot['GMAG']-phot['RMAG'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((rzobslim, grobslim)))
        ax1.set_xlabel(r'$(r - z)_{\rm obs}$')
        ax1.set_ylabel(r'$(g - r)_{\rm obs}$')
        ax1.set_xlim(rzobslim)
        ax1.set_ylim(grobslim)

        hb = ax2.hexbin(phot['GMAG'], phot['GMAG']-phot['RMAG'], 
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((gobslim, grobslim)))
        ax2.set_xlabel(r'$g_{\rm obs}$')
        ax2.set_ylim(grobslim)
        ax2.set_xlim(gobslim)

        cax = fig.add_axes([0.88, 0.12, 0.02, 0.83])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False)
        fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')

        for aa in (ax1, ax2):
            aa.grid(True)

        plt.subplots_adjust(left=0.12, top=0.95, right=0.85, bottom=0.19, wspace=0.07)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()
            
    def elg_rest(phot, meta, bins=None, png=None):
        zlim = (0.5, 1.6)
        Mglim = (-18, -25)
        grlim = (-0.5, 1.0)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

        ax1.hexbin(meta['Z'], phot['ABSMAG_G'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, Mglim)))
        ax1.set_ylim(Mglim)
        ax1.set_xlim(zlim)
        ax1.set_xlabel('Redshift')
        ax1.set_ylabel(r'$M_{0.0g}$')
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))        

        if bins:
            dx, dy = bins['zobj']['del'], bins['Mg']['del']
            [ax1.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['Mg']['grid']]

        ax2.hexbin(meta['Z'], phot['ABSMAG_G']-phot['ABSMAG_R'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, grlim)))
        ax2.set_xlim(zlim)
        ax2.set_ylim(grlim)
        ax2.set_xlabel('Redshift')
        ax2.set_ylabel(r'$^{0.0}(g - r)$', labelpad=-10)
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.2))        
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.5))        

        if bins:
            dx, dy = bins['zobj']['del'], bins['gr']['del']
            [ax2.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['gr']['grid']]

        hb = ax3.hexbin(phot['ABSMAG_G'], phot['ABSMAG_G']-phot['ABSMAG_R'], 
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((Mglim, grlim)))
        ax3.set_xlabel(r'$M_{0.0g}$')
        ax3.set_ylabel(r'$^{0.0}(g - r)$', labelpad=-10)
        ax3.set_xlim(Mglim)
        ax3.set_ylim(grlim)
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(0.5))        

        if bins:
            dx, dy = bins['Mg']['del'], bins['gr']['del']
            [ax3.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mg']['grid'] for yy in bins['gr']['grid']]
            
        ax4.axis('off')

        cax = fig.add_axes([0.49, 0.12, 0.02, 0.36])
        #cax = fig.add_axes([0.54, 0.4, 0.35, 0.03])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False)
        fig.colorbar(hb, format=formatter, label='Number of Galaxies',
                     cax=cax)#, orientation='horizontal')

        for aa in (ax1, ax2, ax3):
            aa.grid(True)

        plt.subplots_adjust(left=0.1, top=0.95, wspace=0.3, hspace=0.3, right=0.88, bottom=0.13)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    def lrg_obs(phot, png=None):
        zobslim = (16, 22)
        W1obslim = (16, 21)
        grobslim = (0.0, 4)
        rzobslim = (0.0, 3)
        rW1obslim = (0.7, 4.5)
        zW1obslim = (0, 2.7)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

        ax1.hexbin(phot['RMAG']-phot['W1MAG'], phot['GMAG']-phot['RMAG'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   #norm=LogNorm(vmin=1, vmax=100),
                   extent=np.hstack((rW1obslim, grobslim)))
        ax1.set_xlabel(r'$(r - W1)_{\rm obs}$')
        ax1.set_ylabel(r'$(g - r)_{\rm obs}$')
        ax1.set_xlim(rW1obslim)
        ax1.set_ylim(grobslim)

        ax2.hexbin(phot['ZMAG']-phot['W1MAG'], phot['RMAG']-phot['ZMAG'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zW1obslim, rzobslim)))

        ax2.set_ylabel(r'$(r - z)_{\rm obs}$')
        ax2.set_xlabel(r'$(z - W1)_{\rm obs}$')
        ax2.set_xlim(zW1obslim)
        ax2.set_ylim(rzobslim)
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
        
        ax3.hexbin(phot['ZMAG'], phot['RMAG']-phot['ZMAG'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zobslim, rzobslim)))
        ax3.set_ylabel(r'$(r - z)_{\rm obs}$')
        ax3.set_xlabel(r'$z_{\rm obs}$')
        ax3.set_xlim(zobslim)
        ax3.set_ylim(rzobslim)
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(1))

        hb = ax4.hexbin(phot['W1MAG'], phot['ZMAG']-phot['W1MAG'], 
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((W1obslim, zW1obslim)))
        ax4.set_ylabel(r'$(z - W1)_{\rm obs}$')
        ax4.set_xlabel(r'$W1_{\rm obs}$')
        ax4.set_xlim(W1obslim)
        ax4.set_ylim(zW1obslim)
        ax4.yaxis.set_major_locator(ticker.MultipleLocator(1))

        cax = fig.add_axes([0.88, 0.12, 0.02, 0.83])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
        cb = fig.colorbar(hb, cax=cax, label='Number of Galaxies',
                          format=formatter)#, ticks=[1, 10, 50])

        for aa in (ax1, ax2, ax3, ax4):
            aa.grid(True)

        plt.subplots_adjust(left=0.1, top=0.95, wspace=0.25, hspace=0.32, right=0.85, bottom=0.13)

        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    def lrg_rest(phot, meta, bins=None, png=None):
        zlim = (0.0, 1.2)
        Mrlim = (-19, -25)
        rW1lim = (-1.4, 1.7)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

        ax1.hexbin(meta['Z'], phot['ABSMAG_R'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, Mrlim)))
        ax1.set_ylim(Mrlim)
        ax1.set_xlim(zlim)
        ax1.set_xlabel('Redshift')
        ax1.set_ylabel(r'$M_{0.0r}$')

        if bins:
            #ax1.add_patch(Rectangle((bins['z']['min'], bins['absmag']['min']),
            #                         np.ptp(bins['z']['grid'])+bins['z']['del'], 
            #                         np.ptp(bins['absmag']['grid'])+bins['absmag']['del'],
            #                         facecolor='none', edgecolor='k', lw=3))
            dx, dy = bins['zobj']['del'], bins['Mr']['del']
            [ax1.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['Mr']['grid']]            

        ax2.hexbin(meta['Z'], phot['ABSMAG_R']-phot['ABSMAG_W1'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, rW1lim)))
        ax2.set_xlabel('Redshift')
        ax2.set_ylabel(r'$^{0.0}(r - W1)$')
        ax2.set_ylim(rW1lim)
        ax2.set_xlim(zlim)

        if bins:
            dx, dy = bins['zobj']['del'], bins['rW1']['del']
            [ax2.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['rW1']['grid']]

        hb = ax3.hexbin(phot['ABSMAG_R'], phot['ABSMAG_R']-phot['ABSMAG_W1'],
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((Mrlim, rW1lim)))
        ax3.set_xlabel(r'$M_{0.0r}$')
        ax3.set_ylabel(r'$^{0.0}(r - W1)$')
        ax3.set_xlim(Mrlim)
        ax3.set_ylim(rW1lim)

        if bins:
            dx, dy = bins['Mr']['del'], bins['rW1']['del']
            [ax3.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mr']['grid'] for yy in bins['rW1']['grid']]    

        ax4.axis('off')
        
        cax = fig.add_axes([0.49, 0.12, 0.02, 0.36])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
        fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')

        for aa in (ax1, ax2, ax3):
            aa.grid(True)

        plt.subplots_adjust(left=0.1, top=0.95, wspace=0.3, hspace=0.3, right=0.88, bottom=0.13)
        
        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()
            
    def lrg_rest2(phot, meta, bins=None, png=None):
        zlim, Mrlim, gilim, rW1lim = (0.0, 1.2), (-19, -25), (0.2, 1.6), (-1.4, 1.7)

        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(18, 10))

        ax1.hexbin(meta['Z'], phot['ABSMAG_R'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, Mrlim)))
        ax1.set_ylim(Mrlim)
        ax1.set_xlim(zlim)
        ax1.set_xlabel('Redshift')
        ax1.set_ylabel(r'$M_{0.0r}$')

        if bins:
            #ax1.add_patch(Rectangle((bins['z']['min'], bins['absmag']['min']),
            #                         np.ptp(bins['z']['grid'])+bins['z']['del'], 
            #                         np.ptp(bins['absmag']['grid'])+bins['absmag']['del'],
            #                         facecolor='none', edgecolor='k', lw=3))
            dx, dy = bins['zobj']['del'], bins['Mr']['del']
            [ax1.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['Mr']['grid']]            

        ax2.hexbin(meta['Z'], phot['ABSMAG_G']-phot['ABSMAG_I'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, gilim)))
        ax2.set_xlim(zlim)
        ax2.set_ylim(gilim)
        ax2.set_xlabel('Redshift')
        ax2.set_ylabel(r'$^{0.0}(g - i)$')

        if bins:
            dx, dy = bins['zobj']['del'], bins['gi']['del']
            [ax2.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['gi']['grid']]

        ax3.hexbin(meta['Z'], phot['ABSMAG_R']-phot['ABSMAG_W1'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((zlim, rW1lim)))
        ax3.set_xlabel('Redshift')
        ax3.set_ylabel(r'$^{0.0}(r - W1)$')
        ax3.set_ylim(rW1lim)
        ax3.set_xlim(zlim)

        if bins:
            dx, dy = bins['zobj']['del'], bins['rW1']['del']
            [ax3.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['zobj']['grid'] for yy in bins['rW1']['grid']]

        ax4.hexbin(phot['ABSMAG_R'], phot['ABSMAG_G']-phot['ABSMAG_I'], 
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((Mrlim, gilim)))
        ax4.set_xlabel(r'$M_{0.0r}$')
        ax4.set_ylabel(r'$^{0.0}(g - i)$')
        ax4.set_xlim(Mrlim)
        ax4.set_ylim(gilim)

        if bins:
            dx, dy = bins['Mr']['del'], bins['gi']['del']
            [ax4.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mr']['grid'] for yy in bins['gi']['grid']]

        ax5.hexbin(phot['ABSMAG_R'], phot['ABSMAG_R']-phot['ABSMAG_W1'],
                   mincnt=mincnt, bins='log', cmap=cmap,
                   #C=cat['weight'], reduce_C_function=np.sum,
                   extent=np.hstack((Mrlim, rW1lim)))
        ax5.set_xlabel(r'$M_{0.0r}$')
        ax5.set_ylabel(r'$^{0.0}(r - W1)$')
        ax5.set_xlim(Mrlim)
        ax5.set_ylim(rW1lim)

        if bins:
            dx, dy = bins['Mr']['del'], bins['rW1']['del']
            [ax5.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mr']['grid'] for yy in bins['rW1']['grid']]    

        hb = ax6.hexbin(phot['ABSMAG_R']-phot['ABSMAG_W1'], phot['ABSMAG_G']-phot['ABSMAG_I'],
                        mincnt=mincnt, bins='log', cmap=cmap,
                        #C=cat['weight'], reduce_C_function=np.sum,
                        extent=np.hstack((rW1lim, gilim)))
        ax6.set_xlabel(r'$^{0.0}(r - W1)$')
        ax6.set_ylabel(r'$^{0.0}(g - i)$')
        ax6.set_ylim(gilim)
        ax6.set_xlim(rW1lim)

        if bins:
            dx, dy = bins['rW1']['del'], bins['gi']['del']
            [ax6.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['rW1']['grid'] for yy in bins['gi']['grid']]

        cax = fig.add_axes([0.9, 0.12, 0.02, 0.83])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
        fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')

        for aa in (ax1, ax2, ax3, ax4, ax5, ax6):
            aa.grid(True)

        #plt.subplots_adjust(wspace=0.35, hspace=0.3, right=0.85)
        plt.subplots_adjust(left=0.1, top=0.95, wspace=0.37, hspace=0.3, right=0.88, bottom=0.13)
        
        if png:
            log.info('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()
            
    # make the plots!
    if targetclass == 'lrg':
        if png_obs:
            lrg_obs(phot, png=png_obs)            
        if png_rest:
            lrg_rest(phot, meta, png=png_rest)    
        if png_rest_bins:
            lrg_rest(phot, meta, bins=bins, png=png_rest_bins)
    elif targetclass == 'elg':
        if png_obs:
            elg_obs(phot, png=png_obs)            
        if png_rest:
            elg_rest(phot, meta, png=png_rest)    
        if png_rest_bins:
            elg_rest(phot, meta, bins=bins, png=png_rest_bins)
    elif targetclass == 'bgs':
        if png_obs:
            bgs_obs(phot, png=png_obs)            
        if png_rest:
            bgs_rest(phot, meta, png=png_rest)    
        if png_rest_bins:
            bgs_rest(phot, meta, bins=bins, png=png_rest_bins)
    else:
        pass

#def qa_tilefile(targetclass, remove_vi=True, min_efftime=10.0,
#                specprod='denali', png=None):
#    """Read the set of tiles used for the templates and make a simple QA plot
#    showing the distribution of effective exposure times.
#
#    """
#
#    #from fastspecfit.templates.sample import select_tiles
#    #tileinfo = select_tiles(targetclass, remove_vi=remove_vi, specprod=specprod
#    #                        min_efftime=min_efftime)
#    #tileinfo = Table.read(tilefile)
#    
#    sns, _ = plot_style()
#
#    log.info('Read {} tiles from {}'.format(len(tileinfo), tilefile))
#
#    xlim = (efftime.min(), efftime.max())
#    fig, ax = plt.subplots(figsize=(9, 6))
#    _ = ax.hist(tileinfo['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
#                label='All Tiles (N={})'.format(len(tileinfo)))
#    _ = ax.hist(targtiles['EFFTIME_SPEC'] / 60, bins=50, range=xlim, alpha=0.9,
#                label='{} Tiles (N={})'.format(targetclass.upper(), len(targtiles)))
#
#    if vitiles:
#      _ = ax.hist(vitiles['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
#                  label='VI Tiles (N={})'.format(len(vitiles)))
#    if shallowtiles:
#      _ = ax.hist(shallowtiles['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
#                  label='Shallow (<{:.0f} min) Tiles (N={})'.format(
#                      min_efftime, len(shallowtiles)))
#
#    ax.set_xlabel('Effective Time (spec, min)')
#    ax.set_ylabel('Number of Tiles')
#
#    ax.legend(loc='upper right', fontsize=16)
#
#    plt.subplots_adjust(right=0.95, top=0.95, bottom=0.17)
#
#    if png:
#        log.info('Writing {}'.format(png))
#        fig.savefig(png)
#        plt.close()

def qa_parent_sample(samplefile, tilefile, targetclass='lrg',
                     specprod='denali', png=None):
    """Build QA showing how the parent sample was selected.

    """
    from fastspecfit.templates.sample import read_fastspecfit, read_parent_sample

    sns, _ = plot_style()        

    tilestable = Table.read(tilefile)
    log.info('Read {} tiles from {}'.format(len(tilestable), tilefile))
    
    allphot, allspec, allmeta = read_fastspecfit(
        tilestable, targetclass=targetclass,
        specprod=specprod)

    phot, spec, meta = read_parent_sample(samplefile)

    nall = len(allphot)
    nparent = len(phot)
    log.info('Read {} objects in the parent sample from {}'.format(nparent, samplefile))

    if targetclass == 'lrg':
        zlim = (-0.05, 1.5)
    elif targetclass == 'elg':
        zlim = (-0.05, 1.8)
    elif targetclass == 'bgs':
        zlim = (-0.05, 0.65)
    else:
        pass

    dchi2lim = (0.8, 4.5)
    #fastspec_chi2lim = (-2, 1)
    fastspec_chi2lim = (-0.1, 1)
    fastphot_chi2lim = (-2.5, 4)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))#, sharey=True)
    
    ax1.hist(allmeta['Z'], bins=75, range=zlim, label='All (N={})'.format(nall))
    ax1.hist(meta['Z'], bins=75, range=zlim, alpha=0.7, label='Parent (N={})'.format(nparent))
    ax1.set_xlim(zlim)
    ax1.set_xlabel('Redshift')
    ax1.set_ylabel('Number of {} Targets'.format(targetclass.upper()))

    ax2.hist(np.log10(allphot['CONTINUUM_CHI2']), bins=75, range=fastphot_chi2lim, label='All (N={})'.format(nall))
    ax2.hist(np.log10(phot['CONTINUUM_CHI2']), bins=75, range=fastphot_chi2lim, alpha=0.7, label='Parent (N={})'.format(nparent))
    ax2.set_xlim(fastphot_chi2lim)
    ax2.set_xlabel(r'$\log_{10}\,\chi^{2}_{\nu}$ [fastphot, continuum]')
    #ax2.set_xlabel(r'$\log_{10}\,\chi^{2}_{\nu}$ [$grzW1W2$ model fit]')
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.tick_right()
    ax2.set_ylabel('Number of {} Targets'.format(targetclass.upper()))

    ax3.hist(np.log10(allmeta['DELTACHI2']), bins=75, range=dchi2lim, label='All (N={})'.format(nall))
    ax3.hist(np.log10(meta['DELTACHI2']), bins=75, range=dchi2lim, alpha=0.7, label='Parent (N={})'.format(nparent))
    ax3.set_xlim(dchi2lim)
    ax3.set_xlabel(r'$\log_{10}\,\Delta\chi^{2}$ [redrock]')
    ax3.set_ylabel('Number of {} Targets'.format(targetclass.upper()))

    #ax4.hist(np.log10(np.abs(allspec['CONTINUUM_SMOOTHCORR_B'])), bins=75, range=fastspec_chi2lim)
    #ax4.hist(np.log10(np.abs(spec['CONTINUUM_SMOOTHCORR_B'])), bins=75, range=fastspec_chi2lim, alpha=0.7)
    ax4.hist(np.log10(allspec['CONTINUUM_CHI2']), bins=75, range=fastspec_chi2lim, label='All (N={})'.format(nall))
    ax4.hist(np.log10(spec['CONTINUUM_CHI2']), bins=75, range=fastspec_chi2lim, alpha=0.7, label='Parent (N={})'.format(nparent))
    ax4.set_xlim(fastspec_chi2lim)
    ax4.set_xlabel(r'$\log_{10}\,\chi^{2}_{\nu}$ [fastspec, continuum]')
    ax4.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax4.yaxis.set_label_position('right')
    ax4.yaxis.tick_right()
    ax4.set_ylabel('Number of {} Targets'.format(targetclass.upper()))

    ax4.legend(loc='upper right', fontsize=14)

    plt.subplots_adjust(left=0.14, wspace=0.09, hspace=0.3, right=0.85, top=0.95, bottom=0.15)

    if png:
        log.info('Writing {}'.format(png))
        fig.savefig(png)
        plt.close()

def build_all_qa(targetclass, templatedir, tilefile=None, samplefile=None,
                 stackfile=None, fastspecfile=None, specprod='denali'):

    from fastspecfit.templates.sample import select_tiles

    png = os.path.join(templatedir, 'qa', '{}-tiles.png'.format(targetclass))
    #select_tiles(targetclass, png=png)

    png = os.path.join(templatedir, 'qa', '{}-parent.png'.format(targetclass))
    #qa_parent_sample(samplefile, tilefile, targetclass=targetclass, specprod=specprod, png=png)

    png_obs = os.path.join(templatedir, 'qa', '{}-obs.png'.format(targetclass))
    png_rest = os.path.join(templatedir, 'qa', '{}-rest.png'.format(targetclass))
    png_rest_bins = os.path.join(templatedir, 'qa', '{}-rest-bins.png'.format(targetclass))
    #qa_photometry(targetclass, samplefile=samplefile, png_obs=png_obs,
    #              png_rest=png_rest, png_rest_bins=png_rest_bins)

    pdffile = os.path.join(templatedir, 'qa', '{}-fastspec-fullspec.pdf'.format(targetclass))
    #qa_fastspec_fullspec(targetclass, fastspecfile=fastspecfile, pdffile=pdffile)

    pdffile = os.path.join(templatedir, 'qa', '{}-fastspec-emlinespec.pdf'.format(targetclass))
    #qa_fastspec_emlinespec(targetclass, fastspecfile=fastspecfile, pdffile=pdffile)        

    if targetclass != 'elg': # no lines in redshift range
        png = os.path.join(templatedir, 'qa', '{}-bpt.png'.format(targetclass))
        qa_bpt(targetclass, fastspecfile=fastspecfile, png=png)



    pdb.set_trace()

