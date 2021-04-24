"""
fastspecfit.templates.qa
========================

QA for templates

"""
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

import seaborn as sns
sns.set(context='talk', style='ticks', palette='deep', font_scale=1.2)#, rc=rc)
colors = sns.color_palette()

cmap = plt.cm.get_cmap('RdYlBu')
    
def qa_obs(phot, png=None):

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
    
    cax = fig.add_axes([0.88, 0.20, 0.02, 0.68])

    formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
    fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')
                 #label=r'$\log_{10}$ (Number of Galaxies)')

    #cb = plt.colorbar(hb)
    #cb.set_label(r'$\log_{10}$ (Number of Galaxies)')  
    
    for aa in (ax1, ax2, ax3, ax4):
        aa.grid(True)

    plt.subplots_adjust(wspace=0.22, hspace=0.32, right=0.85, bottom=0.15)
    
    if png:
        print('Writing {}'.format(png))
        fig.savefig(png)
    
def qa_rest(phot, spec, meta, bins=None, png=None):

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
    
    cax = fig.add_axes([0.88, 0.13, 0.02, 0.75])
    formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
    fig.colorbar(hb, cax=cax, format=formatter, label='Number of Galaxies')

    for aa in (ax1, ax2, ax3, ax4, ax5, ax6):
        aa.grid(True)
    
    plt.subplots_adjust(wspace=0.35, hspace=0.3, right=0.85)
    
    if png:
        print('Writing {}'.format(png))
        fig.savefig(png)
