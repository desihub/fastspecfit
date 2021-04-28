"""
fastspecfit.templates.templates
===============================

Code for building templates of various target classes. Called by bin/desi-templates.

"""
import pdb # for debugging

import os, time
import numpy as np
import fitsio
from astropy.table import Table

from desiutil.log import get_logger
log = get_logger()

def read_tileinfo(targetclass, remove_vi=False, min_efftime=None,
                  specprod='denali', png=None):
    """Optionally remove tiles which are being visually inspected.
     /global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/

    """
    vilookup = {'lrg': [80605, 80609],
                'elg': [80606, 80608, 86010],
                'bgs': [80613]}
    
    reduxdir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', specprod)
    tileinfo = Table.read(os.path.join(reduxdir, 'tiles-{}.csv'.format(specprod)))
    tileinfo = tileinfo[tileinfo['SURVEY'] != 'unknown']
    tileinfo = tileinfo[np.argsort(tileinfo['TILEID'])]

    itargtiles = [targetclass in program for program in tileinfo['FAPRGRM']]
    targtiles = tileinfo[itargtiles]
    
    efftime = tileinfo['EFFTIME_SPEC'] / 60

    if targetclass and remove_vi:
        ivitiles = np.isin(tileinfo['TILEID'], vilookup[targetclass])
        vitiles = tileinfo[ivitiles]
        log.info('Removing {} VI tiles.'.format(np.sum(ivitiles)))
        tileinfo = tileinfo[np.logical_not(ivitiles)]
    else:
        vitiles = None
        
    if min_efftime:
        ishallowtiles = tileinfo['EFFTIME_SPEC'] / 60 <= min_efftime
        shallowtiles = tileinfo[ishallowtiles]
        log.info('Removing {} tiles with efftime < {:.1f} min.'.format(
            np.sum(ishallowtiles), min_efftime))
        tileinfo = tileinfo[np.logical_not(ishallowtiles)]
    else:
        shallowtiles = None

    if png:
        import matplotlib.pyplot as plt
        from fastspecfit.templates.qa import plot_style
        sns, _ = plot_style()

        xlim = (efftime.min(), efftime.max())
        fig, ax = plt.subplots(figsize=(9, 6))
        _ = ax.hist(tileinfo['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
                    label='All Tiles (N={})'.format(len(tileinfo)))
        _ = ax.hist(targtiles['EFFTIME_SPEC'] / 60, bins=50, range=xlim, alpha=0.9,
                    label='{} Tiles (N={})'.format(targetclass.upper(), len(targtiles)))
        
        if vitiles:
          _ = ax.hist(vitiles['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
                      label='VI Tiles (N={})'.format(len(vitiles)))
        if shallowtiles:
          _ = ax.hist(shallowtiles['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
                      label='Shallow (<{:.0f} min) Tiles (N={})'.format(
                          min_efftime, len(shallowtiles)))
          
        ax.set_xlabel('Effective Time (spec, min)')
        ax.set_ylabel('Number of Tiles')
           
        ax.legend(loc='upper right', fontsize=16)

        plt.subplots_adjust(right=0.95, top=0.95, bottom=0.17)

        print('Writing {}'.format(png))
        fig.savefig(png)
        plt.close()
        
    return tileinfo

def read_fastspecfit(tileinfo, specprod='denali', targetclass='lrg'):
    """Read the fastspecfit output for this production.

    """
    from desitarget.targets import main_cmx_or_sv

    fastspecdir = os.path.join(os.getenv('FASTSPECFIT_DATA'), specprod, 'tiles')
    specfile = os.path.join(fastspecdir, 'merged', 'fastspec-{}-cumulative.fits'.format(specprod))
    photfile = os.path.join(fastspecdir, 'merged', 'fastphot-{}-cumulative.fits'.format(specprod))

    spec = Table(fitsio.read(specfile, 'FASTSPEC'))
    meta = Table(fitsio.read(specfile, 'METADATA'))
    phot = Table(fitsio.read(photfile, 'FASTPHOT'))

    assert(np.all(spec['TARGETID'] == phot['TARGETID']))
    
    log.info('Read {} objects from {}'.format(len(spec), specfile))
    log.info('Read {} objects from {}'.format(len(phot), photfile))
    
    ontiles = np.where(np.isin(meta['TILEID'], tileinfo['TILEID']))[0]
    spec = spec[ontiles]
    meta = meta[ontiles]
    phot = phot[ontiles]
    
    log.info('Keeping {} objects on {}/{} unique tiles.'.format(
        len(ontiles), len(np.unique(meta['TILEID'])), len(tileinfo)))
    
    ngal = len(spec)
    
    # convenience magnitudes and targeting variables
    for band in ('G', 'R', 'Z', 'W1'):
        phot['{}MAG'.format(band)] = np.zeros(ngal, 'f4')
        good = np.where(meta['FLUX_{}'.format(band)] > 0)[0]
        phot['{}MAG'.format(band)][good] = 22.5 - 2.5 * np.log10(meta['FLUX_{}'.format(band)][good])
        
    for band in ('G', 'R', 'Z'):
        phot['{}FIBERMAG'.format(band)] = np.zeros(ngal, 'f4')
        good = np.where(meta['FIBERFLUX_{}'.format(band)] > 0)[0]
        phot['{}FIBERMAG'.format(band)][good] = 22.5 - 2.5 * np.log10(meta['FIBERFLUX_{}'.format(band)][good])

    targs = ['BGS_ANY', 'ELG', 'LRG', 'QSO']
    targcols = ['BGS', 'ELG', 'LRG', 'QSO']
    for targcol in targcols:
        spec[targcol] = np.zeros(ngal, bool)
        phot[targcol] = np.zeros(ngal, bool)
        
    for tile in tileinfo['TILEID']:
        I = np.where(meta['TILEID'] == tile)[0]
        if len(I) == 0:
            continue

        (desicol, bgscol, mwscol), (desimask, bgsmask, mwsmask), survey = main_cmx_or_sv(meta[I])

        for targcol, targ in zip(targcols, targs):
            phot[targcol][I] = meta[desicol][I] & desimask.mask(targ) != 0
            spec[targcol][I] = meta[desicol][I] & desimask.mask(targ) != 0
            
    #for targcol in targcols:
    #    log.info('  {}: {}'.format(targcol, np.sum(phot[targcol])))

    itarg = phot[targetclass.upper().replace('BGS', 'BGS_ANY')]
    log.info('Keeping {} {} targets.'.format(np.sum(itarg), targetclass.upper()))

    phot = phot[itarg]
    spec = spec[itarg]
    meta = meta[itarg]
    
    return phot, spec, meta

def _select_lrgs(iparent, phot, spec, meta, z_minmax=(0.1, 1.1), Mr_minmax=None,
                 gi_minmax=None, rW1_minmax=None, targetclass='LRG', verbose=False,
                 png=None):
    """Select LRGs. This method should be called with:

    select_parent_sample(phot, spec, meta, targetclass='lrg')

    """
    if z_minmax:
        iparent *= (meta['Z'] > z_minmax[0]) * (meta['Z'] < z_minmax[1])
    if Mr_minmax:
        iparent *= (phot['ABSMAG_R'] > Mr_minmax[0]) * (phot['ABSMAG_R'] < Mr_minmax[1])
    if gi_minmax:
        gi = phot['ABSMAG_G'] - phot['ABSMAG_I']
        iparent *= (gi > gi_minmax[0]) * (gi < gi_minmax[1])
    if rW1_minmax:
        rW1 = phot['ABSMAG_R'] - phot['ABSMAG_W1']
        iparent *= (rW1 > rW1_minmax[0]) * (rW1 < rW1_minmax[1])

    if verbose:
        log.info('Selecting {}/{} {}s.'.format(np.sum(iparent), len(meta),
                                               targetclass.upper()))

    if png:
        import matplotlib.pyplot as plt
        from fastspecfit.templates.qa import plot_style
        sns, _ = plot_style()        

        #zlim = (meta['Z'].min(), meta['Z'].max())
        zlim = (0.0, 1.5)
        chi2lim = (-2, 4)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

        ax1.hist(meta['Z'], bins=75, range=zlim,
                 label='All (N={})'.format(len(phot)))
        ax1.hist(meta['Z'][iparent], bins=75, range=zlim, alpha=0.7,
                 label='Parent (N={})'.format(np.sum(iparent)))
        ax1.set_xlim(zlim)
        ax1.set_xlabel('Redshift')
        ax1.set_ylabel('Number of {} Targets'.format(targetclass))
        
        ax2.hist(np.log10(phot['CONTINUUM_CHI2']), bins=75, range=chi2lim,
                 label='All (N={})'.format(len(phot)))
        ax2.hist(np.log10(phot['CONTINUUM_CHI2'][iparent]), bins=75, range=chi2lim, alpha=0.7,
                 label='Parent (N={})'.format(np.sum(iparent)))
        ax2.set_xlim(chi2lim)
        ax2.set_xlabel(r'$\log_{10}\,\chi^{2}_{\nu}$ [fastphot]')
        
        ax1.legend(loc='upper right', fontsize=14)

        plt.subplots_adjust(left=0.14, wspace=0.09, right=0.95, top=0.95, bottom=0.17)

        print('Writing {}'.format(png))
        fig.savefig(png)
        plt.close()

    return iparent

def select_parent_sample(phot, spec, meta, targetclass='lrg', deltachi2_cut=40,
                         continuumchi2_cut=20, smoothcorr_cut=10, return_indices=False,
                         verbose=False, **kwargs):
    """High-level sample selection.

    """
    iparent = (
        #(meta['SPECTYPE'] == 'GALAXY') *
        (meta['DELTACHI2'] > deltachi2_cut) * 
        (meta['FLUX_G'] > 0) * 
        (meta['FLUX_R'] > 0) * 
        (meta['FLUX_Z'] > 0) * 
        (meta['FLUX_W1'] > 0) *
        (phot['CONTINUUM_CHI2'] < continuumchi2_cut) *
        (np.abs(spec['CONTINUUM_SMOOTHCORR_B']) < smoothcorr_cut) *
        (np.abs(spec['CONTINUUM_SMOOTHCORR_R']) < smoothcorr_cut) *
        (np.abs(spec['CONTINUUM_SMOOTHCORR_Z']) < smoothcorr_cut)
    )

    if targetclass == 'lrg':
        _select_lrgs(iparent, phot, spec, meta, verbose=verbose, **kwargs)

    if return_indices:
        return np.where(iparent)[0]
    else:
        return phot[iparent], spec[iparent], meta[iparent]

def lrg_stacking_bins(verbose=False):

    # define the stacking limits and the number of bin *centers*
    zlim, nz = [0.1, 1.1], 10
    Mrlim, nMr = [-24.5, -20], 9
    gilim, ngi = [0.4, 1.4], 5 
    rW1lim, nrW1 = [-1.0, 1.0], 8
    
    dz = (zlim[1] - zlim[0]) / nz
    dMr = (Mrlim[1] - Mrlim[0]) / nMr
    dgi = (gilim[1] - gilim[0]) / ngi
    drW1 = (rW1lim[1] - rW1lim[0]) / nrW1
    
    # build the array of (left) bin *edges*
    zgrid = np.arange(zlim[0], zlim[1], dz)
    Mrgrid = np.arange(Mrlim[0], Mrlim[1], dMr)
    gigrid = np.arange(gilim[0], gilim[1], dgi)
    rW1grid = np.arange(rW1lim[0], rW1lim[1], drW1)
    
    bins = {
        'zobj': {'min': zlim[0], 'max': zlim[1], 'del': dz, 'grid': zgrid},
        'Mr': {'min': Mrlim[0], 'max': Mrlim[1], 'del': dMr, 'grid': Mrgrid}, 
        'gi': {'min': gilim[0], 'max': gilim[1], 'del': dgi, 'grid': gigrid}, 
        'rW1': {'min': rW1lim[0], 'max': rW1lim[1], 'del': drW1, 'grid': rW1grid},
        }
    
    nbins = 1
    for key in bins.keys():
        nbins *= len(bins[key]['grid'])
        if verbose:
            log.debug(len(bins[key]['grid']))
    if verbose:
        log.debug(nbins)

    return bins, nbins

def spectra_in_bins(tileinfo, minperbin=3, targetclass='lrg', specprod='denali',
                    bins=None, CFit=None, verbose=False):
    """Select objects in bins of rest-frame properties.
    
    """
    from astropy.table import Table, Column
    
    fastspecdir = os.path.join(os.getenv('FASTSPECFIT_DATA'), specprod, 'tiles')

    def sample_template(bins):
        sample1 = Table()
        sample1.add_column(Column(name='TILE', dtype=np.int32, length=1))
        sample1.add_column(Column(name='IBIN', dtype=np.int32, length=1))
        sample1.add_column(Column(name='NOBJ', dtype=np.int32, length=1))
        for binkey in bins.keys():
            sample1.add_column(Column(name=binkey.upper(), dtype='f4', length=1)) # mean bin center
            sample1.add_column(Column(name='{}MIN'.format(binkey.upper()), dtype='f4', length=1))
            sample1.add_column(Column(name='{}MAX'.format(binkey.upper()), dtype='f4', length=1))
            
        return sample1
        
    samples, data = [], [] # these lists will be aligned

    dz = bins['zobj']['del']
    dMr = bins['Mr']['del']
    dgi = bins['gi']['del']
    drW1 = bins['rW1']['del']
    
    wave = None

    tiles = tileinfo['TILEID']
    for itile, tile in enumerate(tiles):
        if itile % 10 == 0:
            log.info('Working on tile {}/{}'.format(itile, len(tiles)))
        
        restfile = os.path.join(fastspecdir, 'deredshifted', '{}-{}-restflux.fits'.format(
            targetclass.lower(), tile))
        if not os.path.isfile(restfile): # not all of them exist
            continue

        # don't read the spectra if we're just counting
        allphot = Table(fitsio.read(restfile, ext='FASTPHOT'))
        allspec = Table(fitsio.read(restfile, ext='FASTSPEC'))
        allmeta = Table(fitsio.read(restfile, ext='METADATA'))
        flux = fitsio.read(restfile, ext='FLUX')
        ivar = fitsio.read(restfile, ext='IVAR')
        # the wavelength vector is identical, so just read one
        if wave is None:
            wave = fitsio.read(restfile, ext='WAVE')

        # select the sample of interest on this tile
        ibin = 0
        #binkeys = bins.keys()
        
        for zmin in bins['zobj']['grid']:
            for Mrmin in bins['Mr']['grid']:
                for gimin in bins['gi']['grid']:
                    for rW1min in bins['rW1']['grid']:
                        I = select_parent_sample(allphot, allspec, allmeta, 
                                                 z_minmax=[zmin, zmin+dz],
                                                 Mr_minmax=[Mrmin, Mrmin+dMr],
                                                 gi_minmax=[gimin, gimin+dgi],
                                                 rW1_minmax=[rW1min, rW1min+drW1],
                                                 verbose=False, return_indices=True)
                        nobj = len(I)
                        
                        if nobj >= minperbin:
                            if verbose:
                                log.info('N={:04d}, z: {:.3f} {:.3f}, Mr: {:.2f} {:.2f}, g-i: {:.3f} {:.3f}, r-W1: {:.3f} {:.3f}'.format(
                                    len(I), zmin, zmin+dz, Mrmin, Mrmin+dMr, gimin, 
                                    gimin+dgi, rW1min, rW1min+drW1))
                                
                            _sample = sample_template(bins)
                            _sample['TILE'] = tile
                            _sample['IBIN'] = ibin
                            _sample['NOBJ'] = len(I)
                            _sample['ZOBJ'] = zmin + dz / 2
                            _sample['ZOBJMIN'] = zmin
                            _sample['ZOBJMAX'] = zmin + dz
                            _sample['MR'] = Mrmin + dMr / 2
                            _sample['MRMIN'] = Mrmin
                            _sample['MRMAX'] = Mrmin + dMr
                            _sample['GI'] = gimin + dgi / 2
                            _sample['GIMIN'] = gimin
                            _sample['GIMAX'] = gimin + dgi
                            _sample['RW1'] = rW1min + drW1 / 2
                            _sample['RW1MIN'] = rW1min
                            _sample['RW1MAX'] = rW1min + drW1
                            samples.append(_sample)
                            
                            _data = {}
                            _data['flux'] = flux[I, :]
                            _data['ivar'] = ivar[I, :]
                            _data['fastphot'] = allphot[I]
                            _data['fastspec'] = allspec[I]
                            _data['metadata'] = allmeta[I]

                            # rebuild the best-fitting continuum model fits
                            if CFit is not None:
                                _data['cflux'] = []
                                for iobj in np.arange(nobj):
                                    cflux1, _ = CFit.SSP2data(
                                        CFit.sspflux, CFit.sspwave, 
                                        redshift=allmeta[I][iobj]['Z'],
                                        AV=allphot['CONTINUUM_AV'][I][iobj],
                                        coeff=allphot['CONTINUUM_COEFF'][I][iobj],# * CFit.massnorm,
                                        synthphot=False)
                                    _data['cflux'].append(cflux1)# * (1 + allmeta[I[iobj]]['Z']) # deredshift
                                _data['cflux'] = np.vstack(_data['cflux'])
                            data.append(_data)
                                
                        ibin += 1 # next bin

    # Stack the bin-level statistics table
    samples = Table(np.hstack(samples))
    data = np.array(data)

    # ...and now stack the data in each (unique) bin number.
    samplestack, sampledata = [], {}
    
    for ibin in sorted(set(samples['IBIN'])):
        I = np.where(ibin == samples['IBIN'])[0]
        _samplestack = samples[[I[0]]].copy()
        _samplestack.remove_column('TILE')
        _samplestack['NOBJ'] = np.sum(samples[I]['NOBJ'])
        samplestack.append(_samplestack)

        _sampledata = {}
        _sampledata['flux'] = np.vstack([_data['flux'] for _data in data[I]])
        _sampledata['ivar'] = np.vstack([_data['ivar'] for _data in data[I]])
        _sampledata['fastphot'] = Table(np.hstack([_data['fastphot'] for _data in data[I]]))
        _sampledata['fastspec'] = Table(np.hstack([_data['fastspec'] for _data in data[I]]))
        _sampledata['metadata'] = Table(np.hstack([_data['metadata'] for _data in data[I]]))

        if CFit is not None:
            _sampledata['cflux'] = np.vstack([_data['cflux'] for _data in data[I]])
        
        sampledata.update({str(ibin): _sampledata})
        del _sampledata, _samplestack

    del data, samples
    samplestack = Table(np.hstack(samplestack))
    
    return samplestack, sampledata, wave




