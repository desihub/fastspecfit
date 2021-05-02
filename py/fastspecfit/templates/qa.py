"""
fastspecfit.templates.qa
========================

QA for templates

"""
import pdb

import os
import numpy as np
import fitsio
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

def remove_undetected_lines(fastspec, linetable, snrmin=3.0):
    """replace weak or undetected emission lines with their upper limits

    """
    linenames = linetable['name']
    
    for linename in linenames:
        amp = fastspec['{}_AMP'.format(linename.upper())].data
        amp_ivar = fastspec['{}_AMP_IVAR'.format(linename.upper())].data
        cont_ivar = fastspec['{}_CONT_IVAR'.format(linename.upper())].data

        fix = np.where(amp_ivar == 0)[0]
        if len(fix) > 0:
            fastspec['{}_AMP'.format(linename.upper())][fix] = 0.0 # fixes a bug with [OIII] 4959

        # fix and then skip tied doublets
        if 'oiii_4959' in linename or 'nii_6548' in linename: 
            #fix = np.where((fastspec['OIII_4959_AMP_IVAR'] > 0) * (fastspec['OIII_5007_AMP_IVAR'] == 0))[0]
            #if len(fix) > 0:
            #    fastspec['OIII_4959_AMP'][fix] = 0.0
            #    fastspec['OIII_4959_AMP_IVAR'][fix] = 0.0
            continue

        snr = amp * np.sqrt(amp_ivar)
        losnr = np.where(np.logical_or((amp_ivar > 0) * (snr < snrmin), (amp_ivar == 0) * (cont_ivar > 0)))[0]

        if len(losnr) > 0:
            csig = 1 / np.sqrt(cont_ivar[losnr])
            if np.any(csig) < 0:
                pdb.set_trace()

            fastspec['{}_AMP'.format(linename.upper())][losnr] = snrmin * csig
            #fastspec['{}_AMP'.format(linename.upper())][losnr] = 0.0
            #fastspec['{}_AMP_IVAR'.format(linename.upper())][losnr] = 0.0
            #fastspec['{}_AMP_IVAR'.format(linename.upper())][losnr] = 0.0

            if linename == 'oiii_5007':
                fastspec['OIII_4959_AMP'][losnr] = fastspec['OIII_5007_AMP'][losnr].data / 2.875
            if linename == 'nii_6584':
                fastspec['NII_6548_AMP'][losnr] = fastspec['NII_6584_AMP'][losnr].data / 2.936

            #if linename == 'nii_6584':
            #    pdb.set_trace()
                
    # corner case of when the stronger doublet is on the edge of the wavelength range
    fix = np.where((fastspec['OIII_5007_AMP'] == 0) * (fastspec['OIII_4959_AMP'] != 0))[0]
    if len(fix) > 0:
        fastspec['OIII_4959_AMP'][fix] = 0.0
        fastspec['OIII_4959_AMP_IVAR'][fix] = 0.0
        
    fix = np.where((fastspec['NII_6584_AMP'] == 0) * (fastspec['NII_6548_AMP'] != 0))[0]
    if len(fix) > 0:
        fastspec['NII_6548_AMP'][fix] = 0.0
        fastspec['NII_6548_AMP_IVAR'][fix] = 0.0

    for linename in linenames:
        amp = fastspec['{}_AMP'.format(linename.upper())].data
        neg = np.where(amp < 0)[0]
        if len(neg) > 0:
            print('Fix {}'.format(linename))
            pdb.set_trace()

    # go back through and update FLUX and EW based on the new line-amplitudes
    redshift = 0.0 # check this...
    for oneline in linetable:
        linename = oneline['name'].upper()

        amp = fastspec['{}_AMP'.format(linename.upper())].data
        amp_ivar = fastspec['{}_AMP_IVAR'.format(linename.upper())].data
        snr = amp * np.sqrt(amp_ivar)
        hisnr = np.where(snr >= snrmin)[0]
        if len(hisnr) == 0:
            continue

        linez = redshift + fastspec['{}_VSHIFT'.format(linename)][hisnr].data / C_LIGHT
        linezwave = oneline['restwave'] * (1 + linez)

        linesigma = fastspec['{}_SIGMA'.format(linename)][hisnr].data # [km/s]
        log10sigma = linesigma / C_LIGHT / np.log(10)     # line-width [log-10 Angstrom]
            
        # get the emission-line flux
        linesigma_ang = linezwave * linesigma / C_LIGHT # [observed-frame Angstrom]
        linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang

        fastspec['{}_FLUX'.format(linename)][hisnr] = fastspec['{}_AMP'.format(linename)][hisnr].data * linenorm

        cpos = np.where(fastspec['{}_CONT'.format(linename)][hisnr] > 0.0)[0]
        if len(cpos) > 0:
            factor = (1 + redshift) / fastspec['{}_CONT'.format(linename)][hisnr][cpos] # --> rest frame
            fastspec['{}_EW'.format(linename)][hisnr][cpos] = fastspec['{}_FLUX'.format(linename)][hisnr][cpos] * factor   # rest frame [A]

    return fastspec

def qa_bpt(fastspecfile, EMFit, png=None):
    """QA of the fastspec emission-line spectra.

    """
    sns, _ = plot_style()

    fastmeta = Table(fitsio.read(fastspecfile, ext='METADATA'))
    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC'))
    nobj = len(fastmeta)

    fastspec = remove_undetected_lines(fastspec, EMFit.linetable)

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
            print('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    good = np.where(
        (fastspec['HALPHA_FLUX'] > 0) * 
        (fastspec['HBETA_FLUX'] > 0) * 
        (fastspec['NII_6584_FLUX'] > 0) * 
        (fastspec['OIII_5007_FLUX'] > 0) 
        #(fastspec['HALPHA_CHI2'] < 1e4)
    )[0]
    
    zz = fastspec['CONTINUUM_Z'][good]
    rW1 = fastmeta['RW1'][good]
    gi = fastmeta['GI'][good]
    ewhb = fastspec['HBETA_EW'][good]
    
    niiha = np.log10(fastspec['NII_6584_FLUX'][good] / fastspec['HALPHA_FLUX'][good])
    oiiihb = np.log10(fastspec['OIII_5007_FLUX'][good] / fastspec['HBETA_FLUX'][good])
    ww = np.where((niiha > -0.05) * (niiha < 0.05) * (oiiihb < -0.5))[0]
    #print(fastspec[good][ww]['HALPHA_FLUX', 'NII_6584_FLUX'])

    _bpt(zz, 'Redshift', vmin=0, vmax=0.5, png=png.replace('.png', '-redshift.png'))
    _bpt(rW1, r'$r-W1$', vmin=-0.3, vmax=0.9, png=png.replace('.png', '-rW1.png'))
    _bpt(gi, r'$g-i$', vmin=0.6, vmax=1.3, png=png.replace('.png', '-gi.png'))
    _bpt(np.log10(ewhb), r'$\log_{10}\,\mathrm{EW}(\mathrm{H}\beta)$', png=png.replace('.png', '-ewhb.png'))

def qa_fastspec_emlinespec(fastspecfile, CFit, EMFit, pdffile=None):
    """QA of the fastspec emission-line spectra.

    """
    from matplotlib.colors import Normalize
    
    sns, _ = plot_style()

    wave = fitsio.read(fastspecfile, ext='WAVE')
    flux = fitsio.read(fastspecfile, ext='FLUX')
    ivar = fitsio.read(fastspecfile, ext='IVAR')

    fastmeta = Table(fitsio.read(fastspecfile, ext='METADATA'))
    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC'))
    nobj = len(fastmeta)

    #if False:
    #    this = (fastmeta['RW1'] == -0.875) * (fastmeta['ZOBJ'] == 0.85) * (fastmeta['MR'] == -22.25)
    #    fastspec['OIII_4959_AMP', 'OIII_4959_AMP_IVAR', 'OIII_5007_AMP', 'OIII_5007_AMP_IVAR'][this]
    #    pdb.set_trace()
    
    fastspec = remove_undetected_lines(fastspec, EMFit.linetable)

    for linename in EMFit.linetable['name']:
        amp = fastspec['{}_AMP'.format(linename.upper())].data
        neg = np.where(amp < 0)[0]
        if len(neg) > 0:
            print('Fix {}'.format(linename))
            pdb.set_trace()

    ncol, nrow = 3, 5
    icam = 0
        
    rW1color = np.unique(fastmeta['RW1'])
    npage = len(rW1color)

    cmap = plt.cm.get_cmap('jet')
    #cmap = sns.color_palette(as_cmap=True)
    cnorm = Normalize(vmin=np.min(fastmeta['ZOBJ']), vmax=np.max(fastmeta['ZOBJ']))
    #cnorm = Normalize(vmin=np.min(rW1color), vmax=np.max(rW1color))
        
    if pdffile:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(pdffile)

    inches_wide = 16
    inches_fullspec = 6
    inches_perline = inches_fullspec / 2.0
    nlinepanels = 4

    nline = len(set(EMFit.linetable['plotgroup']))
    nlinerows = np.ceil(nline / nlinepanels).astype(int)
    nrows = 1 + nlinerows

    height_ratios = np.hstack([1, [0.5]*nlinerows])

    plotsig_default = 300.0 # [km/s]
    meanwaves, deltawaves, sigmas, linenames = [], [], [], []
    for plotgroup in set(EMFit.linetable['plotgroup']):
        I = np.where(plotgroup == EMFit.linetable['plotgroup'])[0]
        linenames.append(EMFit.linetable['nicename'][I[0]])
        meanwaves.append(np.mean(EMFit.linetable['restwave'][I]))
        deltawaves.append((np.max(EMFit.linetable['restwave'][I]) - np.min(EMFit.linetable['restwave'][I])) / 2)
        sigmas.append(plotsig_default)
    srt = np.argsort(meanwaves)
    meanwaves = np.hstack(meanwaves)[srt]
    deltawaves = np.hstack(deltawaves)[srt]
    sigmas = np.hstack(sigmas)[srt]
    linenames = np.hstack(linenames)[srt]
    
    # make the plot!
    for ipage in np.arange(npage):#[:2]:
        zindx = np.where(rW1color[ipage] == fastmeta['RW1'])[0]
        zfastmeta = fastmeta[zindx]
        zfastspec = fastspec[zindx]

        srt = np.argsort(zfastmeta['MR'])
        zfastmeta = zfastmeta[srt]
        zfastspec = zfastspec[srt]

        absmag = np.unique(zfastmeta['MR'])
        nabspage = len(absmag)

        for iabspage in np.arange(nabspage):#[:1]:#[::2]:
            absindx = np.where((absmag[iabspage] == zfastmeta['MR']))[0]
            absfastmeta = zfastmeta[absindx]
            absfastspec = zfastspec[absindx]

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
        
            for iplot, indx in enumerate(zindx[absindx]):
                print(ipage, iabspage, iplot, len(zindx), len(absindx))

                modelwave, continuum, smooth_continuum, emlinemodel, data = rebuild_fastspec_spectrum(
                    fastspec[indx], wave, flux[indx, :], ivar[indx, :], CFit, EMFit)
                #pdb.set_trace()

                redshift = data['zredrock']
                emlineflux = data['flux'][icam] - continuum - smooth_continuum

                modelwave /= (1+redshift) # rest-frame

                label = 'z=[{:.1f}-{:.1f}] (N={})'.format(
                    fastmeta['ZOBJMIN'][indx], fastmeta['ZOBJMAX'][indx],
                    np.sum(fastmeta['ZOBJ'][zindx[absindx]] == fastmeta['ZOBJ'][indx]))
                #label = '[{:.2f},{:.2f}],[{:.1f},{:.1f}]'.format(
                #    fastmeta['ZOBJMIN'][indx], fastmeta['ZOBJMAX'][indx],
                #    fastmeta['GIMIN'][indx], fastmeta['GIMAX'][indx])
                #bigax.plot(modelwave/(1+redshift), emlineflux, color='gray')
                bigax.plot(modelwave, emlinemodel, label=label, color=cmap(cnorm(fastmeta['ZOBJ'][indx])))

                if -np.max(emlinemodel)*0.05 < bigymin:
                    bigymin = -np.max(emlinemodel)*0.05
                if np.max(emlinemodel)*1.1 > bigymax:
                    bigymax = np.max(emlinemodel)*1.1

                # zoom in on individual emission lines
                for iax, (meanwave, deltawave, sig, linename) in enumerate(zip(meanwaves, deltawaves, sigmas, linenames)):
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
            bigax.set_title(
                r'${:.2f}<r-W1<{:.2f}\ {:.1f}<M_{{r}}<{:.1f}$'.format(
                absfastmeta['RW1MIN'][0], absfastmeta['RW1MAX'][0],
                absfastmeta['MRMIN'][0], absfastmeta['MRMAX'][0]
                ))
            #bigax.set_xlabel('Observed-frame Wavelength ($\AA$)')

            plt.subplots_adjust(wspace=0.28, left=0.07, right=0.95, top=0.95, bottom=0.1)
            
            if pdffile:
                pdf.savefig(fig)
                
            plt.close()

    if pdffile:
        log.info('Writing {}'.format(pdffile))
        pdf.close()

    pdb.set_trace()

def qa_fastspec_fullspec(fastspecfile, CFit, EMFit, pdffile=None):
    """Full-spectrum QA."""

    sns, _ = plot_style()

    wave = fitsio.read(fastspecfile, ext='WAVE')
    flux = fitsio.read(fastspecfile, ext='FLUX')
    ivar = fitsio.read(fastspecfile, ext='IVAR')

    fastmeta = Table(fitsio.read(fastspecfile, ext='METADATA'))
    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC'))
    nobj = len(fastmeta)

    ncol, nrow = 3, 5
    icam = 0
        
    zobj = np.unique(fastmeta['ZOBJ'])
    npage = len(zobj)
    
    if pdffile:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(pdffile)

    for ipage in np.arange(npage):#[:1]:
        zindx = np.where(zobj[ipage] == fastmeta['ZOBJ'])[0]
        zfastmeta = fastmeta[zindx]
        zfastspec = fastspec[zindx]

        srt = np.argsort(zfastmeta['MR'])
        zfastmeta = zfastmeta[srt]
        zfastspec = zfastspec[srt]

        absmag = np.unique(zfastmeta['MR'])
        nabspage = len(absmag)

        for iabspage in np.arange(nabspage):#[:2]:

            absindx = np.where((absmag[iabspage] == zfastmeta['MR']))[0]
            absfastmeta = zfastmeta[absindx]
            absfastspec = zfastspec[absindx]

            fig, allax = plt.subplots(nrow, ncol, figsize=(12, 16), sharex=True, sharey=True)
            for iplot, (indx, ax) in enumerate(zip(zindx[absindx], allax.flatten())):
                print(ipage, iabspage, iplot, len(zindx), len(absindx))

                # rebuild the best-fitting spectrum
                modelwave, continuum, smooth_continuum, emlinemodel, data = rebuild_fastspec_spectrum(
                    fastspec[indx], wave, flux[indx, :], ivar[indx, :], CFit, EMFit)

                #sigma, _ = ivar2var(data['ivar'][icam], sigma=True)
                #ax.fill_between(data['wave'][icam], data['flux'][icam]-sigma,
                #                            data['flux'][icam]+sigma, color='skyblue')
                ax.plot(data['wave'][icam], data['flux'][icam], color='skyblue')
                ax.plot(modelwave[::3], (continuum+emlinemodel)[::3], color='firebrick', alpha=0.5)
                ax.plot(modelwave[::3], continuum[::3], color='blue', alpha=0.5)
                #ax.plot(modelwave[::3], (continuum+smooth_continuum)[::3], color='gray', alpha=0.3)
                ax.plot(modelwave[::3], smooth_continuum[::3], color='gray', alpha=0.7)

                ymin, ymax = 1e6, -1e6

                filtflux = median_filter(data['flux'][icam], 51, mode='nearest')
                sigflux = np.std(data['flux'][icam][data['ivar'][icam] > 0])
                if -2 * sigflux < ymin:
                    ymin = -2 * sigflux
                if sigflux * 5 > ymax:
                    ymax = sigflux * 5
                if np.max(filtflux) > ymax:
                    ymax = np.max(filtflux) * 1.4

                ax.text(0.96, 0.06, '\n'.join(( r'${:.1f}<g-i<{:.1f}$'.format(fastmeta['GIMIN'][indx], fastmeta['GIMAX'][indx]),
                                                r'${:.2f}<r-W1<{:.2f}$'.format(fastmeta['RW1MIN'][indx], fastmeta['RW1MAX'][indx]) )),
                                                ha='right', va='bottom', transform=ax.transAxes, fontsize=10,
                                                bbox=dict(boxstyle='round', facecolor='gray', alpha=0.25))
                ax.text(0.04, 0.96,
                        '\n'.join(( 'N={}, S/N={:.1f}'.format(fastmeta['NOBJ'][indx], fastspec['CONTINUUM_SNR_ALL'][indx]), )),
                    ha='left', va='top', transform=ax.transAxes, fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='gray', alpha=0.25))

                ax.set_xlim(modelwave.min(), modelwave.max())
                ax.set_ylim(ymin, ymax)
                ax.set_xticklabels([])
                ax.set_yticklabels([])

                plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.07, right=0.95, top=0.95, bottom=0.1)

                if iplot == ncol*nrow-1:
                    break
            
            fig.text(0.52, 0.968, r'${:.1f}<z<{:.1f}\ {:.1f}<M_{{r}}<{:.1f}$'.format(
                absfastmeta['ZOBJMIN'][0], absfastmeta['ZOBJMAX'][0],
                absfastmeta['MRMIN'][0], absfastmeta['MRMAX'][0]),
                ha='center', va='center', fontsize=22)

            for rem in np.arange(ncol*nrow-iplot-1)+iplot+1:
                allax.flatten()[rem].axis('off')
                
            if pdffile:
                pdf.savefig(fig)
                
            plt.close()

    if pdffile:
        log.info('Writing {}'.format(pdffile))
        pdf.close()



        
def qa_photometry_lrg(phot, spec, meta, bins=None, png_obs=None,
                      png_rest=None, png_rest_bins=None):

    cmap = plt.cm.get_cmap('RdYlBu')    

    def obs(phot, png=None):
        zobslim = (16, 22)
        W1obslim = (16, 21)
        grobslim = (-0.2, 5)
        rzobslim = (0.3, 3)
        zW1obslim = (0, 2.8)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

        ax1.hexbin(phot['RMAG']-phot['ZMAG'], phot['GMAG']-phot['RMAG'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,
                   cmap=cmap)
        ax1.set_xlabel(r'$(r - z)_{\rm obs}$')
        ax1.set_ylabel(r'$(g - r)_{\rm obs}$')
        ax1.set_xlim(rzobslim)
        ax1.set_ylim(grobslim)

        ax2.hexbin(phot['ZMAG']-phot['W1MAG'], phot['RMAG']-phot['ZMAG'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,                    
                   cmap=cmap)

        ax2.set_ylabel(r'$(r - z)_{\rm obs}$')
        ax2.set_xlabel(r'$(z - W1)_{\rm obs}$')
        ax2.set_xlim(zW1obslim)
        ax2.set_ylim(rzobslim)

        ax3.hexbin(phot['ZMAG'], phot['RMAG']-phot['ZMAG'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,
                   cmap=cmap)
        ax3.set_ylabel(r'$(r - z)_{\rm obs}$')
        ax3.set_xlabel(r'$z$')
        ax3.set_xlim(zobslim)
        ax3.set_ylim(rzobslim)

        hb = ax4.hexbin(phot['W1MAG'], phot['ZMAG']-phot['W1MAG'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,
                   cmap=cmap)
        ax4.set_ylabel(r'$(z - W1)_{\rm obs}$')
        ax4.set_xlabel(r'$W1$')
        ax4.set_xlim(W1obslim)
        ax4.set_ylim(zW1obslim)

        cax = fig.add_axes([0.88, 0.12, 0.02, 0.83])
        formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
        fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')

        for aa in (ax1, ax2, ax3, ax4):
            aa.grid(True)

        plt.subplots_adjust(left=0.1, top=0.95, wspace=0.25, hspace=0.32, right=0.85, bottom=0.13)

        if png:
            print('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    def rest(phot, spec, meta, bins=None, png=None):
        zlim, Mrlim, gilim, rW1lim = (0.0, 1.2), (-19, -25), (0.2, 1.6), (-1.4, 1.4)

        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(18, 10))

        ax1.hexbin(meta['Z'], phot['ABSMAG_R'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,
                   cmap=cmap) 
        ax1.set_ylim(Mrlim)
        ax1.set_xlim(zlim)
        ax1.set_xlabel('Redshift')
        ax1.set_ylabel(r'$M_{0.0r}$')

        if bins:
            #ax1.add_patch(Rectangle((bins['z']['min'], bins['absmag']['min']),
            #                         np.ptp(bins['z']['grid'])+bins['z']['del'], 
            #                         np.ptp(bins['absmag']['grid'])+bins['absmag']['del'],
            #                         facecolor='none', edgecolor='k', lw=3))
            dx, dy = bins['z']['del'], bins['Mr']['del']
            [ax1.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['z']['grid'] for yy in bins['Mr']['grid']]            

        ax2.hexbin(meta['Z'], phot['ABSMAG_G']-phot['ABSMAG_I'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,           
                   cmap=plt.cm.get_cmap('RdYlBu'))
        ax2.set_xlim(zlim)
        ax2.set_ylim(gilim)
        ax2.set_xlabel('Redshift')
        ax2.set_ylabel(r'$^{0.0}(g - i)$')

        if bins:
            dx, dy = bins['z']['del'], bins['gi']['del']
            [ax2.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['z']['grid'] for yy in bins['gi']['grid']]

        ax3.hexbin(meta['Z'], phot['ABSMAG_R']-phot['ABSMAG_W1'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,               
                   cmap=plt.cm.get_cmap('RdYlBu'))
        ax3.set_xlabel('Redshift')
        ax3.set_ylabel(r'$^{0.0}(r - W1)$')
        ax3.set_ylim(rW1lim)
        ax3.set_xlim(zlim)

        if bins:
            dx, dy = bins['z']['del'], bins['rW1']['del']
            [ax3.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['z']['grid'] for yy in bins['rW1']['grid']]

        ax4.hexbin(phot['ABSMAG_R'], phot['ABSMAG_G']-phot['ABSMAG_I'], mincnt=1, bins='log', 
                        #C=cat['weight'], reduce_C_function=np.sum,                    
                        cmap=plt.cm.get_cmap('RdYlBu'))
        ax4.set_xlabel(r'$M_{0.0r}$')
        ax4.set_ylabel(r'$^{0.0}(g - i)$')
        ax4.set_xlim(Mrlim)
        ax4.set_ylim(gilim)

        if bins:
            dx, dy = bins['Mr']['del'], bins['gi']['del']
            [ax4.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mr']['grid'] for yy in bins['gi']['grid']]

        ax5.hexbin(phot['ABSMAG_R'], phot['ABSMAG_R']-phot['ABSMAG_W1'], mincnt=1, bins='log', 
                   #C=cat['weight'], reduce_C_function=np.sum,                    
                   cmap=plt.cm.get_cmap('RdYlBu'))
        ax5.set_xlabel(r'$M_{0.0r}$')
        ax5.set_ylabel(r'$^{0.0}(r - W1)$')
        ax5.set_xlim(Mrlim)
        ax5.set_ylim(rW1lim)

        if bins:
            dx, dy = bins['Mr']['del'], bins['rW1']['del']
            [ax5.add_patch(Rectangle((xx, yy), dx, dy, facecolor='none', edgecolor='k'))
             for xx in bins['Mr']['grid'] for yy in bins['rW1']['grid']]    

        hb = ax6.hexbin(phot['ABSMAG_R']-phot['ABSMAG_W1'], phot['ABSMAG_G']-phot['ABSMAG_I'], mincnt=1, bins='log', 
                        #C=cat['weight'], reduce_C_function=np.sum,               
                        cmap=plt.cm.get_cmap('RdYlBu'))
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
            print('Writing {}'.format(png))
            fig.savefig(png)
            plt.close()

    obs(phot, png=png_obs)
    rest(phot, spec, meta, png=png_rest)
    rest(phot, spec, meta, bins, png=png_rest_bins)

