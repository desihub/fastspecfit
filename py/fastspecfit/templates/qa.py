"""
fastspecfit.templates.qa
========================

QA for templates

"""
import pdb

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

cmap = plt.cm.get_cmap('RdYlBu')
    
def qa_photometry_lrg(phot, spec, meta, bins=None, png_obs=None,
                      png_rest=None, png_rest_bins=None):

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

def qa_fastspec(fastspecfile, CFit, EMFit, png=None):
    """Spectral QA.

    """
    import os
    import numpy as np
    import fitsio
    from astropy.table import Table
    from scipy.ndimage import median_filter

    from fastspecfit.util import ivar2var
    from fastspecfit.templates.templates import rebuild_fastspec_spectrum

    sns, _ = plot_style()

    from matplotlib import colors
    col1 = [colors.to_hex(col) for col in ['skyblue', 'darkseagreen', 'tomato']]
    col2 = [colors.to_hex(col) for col in ['navy', 'forestgreen', 'firebrick']]
    col3 = [colors.to_hex(col) for col in ['blue', 'green', 'red']]

    wave = fitsio.read(fastspecfile, ext='WAVE')
    flux = fitsio.read(fastspecfile, ext='FLUX')
    ivar = fitsio.read(fastspecfile, ext='IVAR')

    fastmeta = Table(fitsio.read(fastspecfile, ext='METADATA'))
    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC'))
    nobj = len(fastmeta)

    #fastmeta = fastmeta.group_by(['ZOBJ', 'MR', 'GI', 'RW1'])
    ncol = 3
    nrow = 5
    bbox = dict(boxstyle='round', facecolor='gray', alpha=0.25)
        
    if png:
        from matplotlib.backends.backend_pdf import PdfPages
        pdffile = png.replace('.png', '.pdf')
        pdf = PdfPages(pdffile)

    zobj = np.unique(fastmeta['ZOBJ'])
    nzpage = len(zobj)

    for izpage in np.arange(nzpage)[:1]:
        zindx = np.where(zobj[izpage] == fastmeta['ZOBJ'])[0]
        zfastmeta = fastmeta[zindx]
        zfastspec = fastspec[zindx]

        srt = np.argsort(zfastmeta['MR'])
        zfastmeta = zfastmeta[srt]
        zfastspec = zfastspec[srt]

        absmag = np.unique(zfastmeta['MR'])
        nabspage = len(absmag)

        for iabspage in np.arange(nabspage)[:2]:

            absindx = np.where((absmag[iabspage] == zfastmeta['MR']))[0]
            absfastmeta = zfastmeta[absindx]
            absfastspec = zfastspec[absindx]

            fig, allax = plt.subplots(nrow, ncol, figsize=(12, 16), sharex=True, sharey=True)
            for iplot, (indx, ax) in enumerate(zip(zindx[absindx], allax.flatten())):
                print(izpage, iabspage, iplot, len(zindx), len(absindx))

                # rebuild the best-fitting spectrum
                modelwave, continuum, smooth_continuum, emlinemodel, data = rebuild_fastspec_spectrum(
                    fastspec[indx], wave, flux[indx, :], ivar[indx, :], CFit, EMFit)

                redshift = data['zredrock']

                icam = 0
                wavelims = (modelwave.min(), modelwave.max())
                ymin, ymax = 1e6, -1e6

                sigma, _ = ivar2var(data['ivar'][icam], sigma=True)
                ax.fill_between(data['wave'][icam], data['flux'][icam]-sigma,
                                            data['flux'][icam]+sigma, color='skyblue')
                ax.plot(modelwave, continuum+emlinemodel, color='firebrick', alpha=0.5)
                ax.plot(modelwave, continuum, color='blue', alpha=0.5)
                ax.plot(modelwave, continuum+smooth_continuum, color='gray', alpha=0.3)
                #ax.plot(modelwave, smooth_continuum, color='gray')#col3[icam])#, alpha=0.3, lw=2)#, color='k')

                # get the robust range
                filtflux = median_filter(data['flux'][icam], 51, mode='nearest')
                sigflux = np.std(data['flux'][icam][data['ivar'][icam] > 0])
                if -2 * sigflux < ymin:
                    ymin = -2 * sigflux
                if sigflux * 5 > ymax:
                    ymax = sigflux * 5
                if np.max(filtflux) > ymax:
                    ymax = np.max(filtflux) * 1.4

                txt = '\n'.join((
                    r'${:.1f}<g-i<{:.1f}$'.format(fastmeta['GIMIN'][indx], fastmeta['GIMAX'][indx]),
                    r'${:.2f}<r-W1<{:.2f}$'.format(fastmeta['RW1MIN'][indx], fastmeta['RW1MAX'][indx])
                    ))
                ax.text(0.96, 0.06, txt, ha='right', va='bottom', transform=ax.transAxes, fontsize=10, bbox=bbox)

                txt = '\n'.join((
                    'N={}, S/N={:.1f}'.format(fastmeta['NOBJ'][indx], fastspec['CONTINUUM_SNR_ALL'][indx]),
                    ))
                ax.text(0.04, 0.96, txt, ha='left', va='top', transform=ax.transAxes, fontsize=10, bbox=bbox)

                ax.set_xlim(wavelims)
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
                
            if png:
                pdf.savefig(fig)

    if png:
        log.info('Writing {}'.format(pdffile))
        pdf.close()
        
    pdb.set_trace()
