"""
fastspecfit.templates.sample
============================

Code for defining samples and reading the data needed to build templates for the
various target classes. Called by bin/desi-templates.

"""
import pdb # for debugging

import os, time
from copy import copy
import numpy as np
import fitsio
from astropy.table import Table

from desiutil.log import get_logger
log = get_logger()

VITILES_TARGETCLASS = {'lrg': [80605, 80609],
                       'elg': [80606, 80608, 80610],
                       'bgs': [80613]}

def select_tiles(targetclass, remove_vi=True, min_efftime=5.0,
                 specprod='denali', outfile=None, png=None):
    """Select tiles to use. Remove low exposure-time tiles and also 
    remove tiles which are being visually inspected.
    
      /global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/

    """
    reduxdir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', specprod)
    tilestable = Table.read(os.path.join(reduxdir, 'tiles-{}.csv'.format(specprod)))
    tilestable = tilestable[tilestable['SURVEY'] != 'unknown']
    tilestable = tilestable[np.argsort(tilestable['TILEID'])]

    itargtiles = [targetclass in program for program in tilestable['FAPRGRM']]
    targtiles = tilestable[itargtiles]
    
    efftime = tilestable['EFFTIME_SPEC'] / 60

    if targetclass and remove_vi:
        ivitiles = np.isin(tilestable['TILEID'], VITILES_TARGETCLASS[targetclass])
        vitiles = tilestable[ivitiles]
        log.info('Removing {} {} VI tiles: {}'.format(np.sum(ivitiles), targetclass.upper(),
                                                      ', '.join(vitiles['TILEID'].astype(str))))
        tilestable = tilestable[np.logical_not(ivitiles)]
    else:
        vitiles = None
        
    if min_efftime:
        ishallowtiles = tilestable['EFFTIME_SPEC'] / 60 <= min_efftime
        shallowtiles = tilestable[ishallowtiles]
        log.info('Removing {} tiles with efftime < {:.1f} min.'.format(
            np.sum(ishallowtiles), min_efftime))
        tilestable = tilestable[np.logical_not(ishallowtiles)]
    else:
        shallowtiles = None

    if outfile:
        log.info('Writing {} tiles to {}'.format(len(tilestable), outfile))
        tilestable.write(outfile, overwrite=True)

    if png:
        import matplotlib.pyplot as plt
        from fastspecfit.templates.qa import plot_style
        sns, _ = plot_style()

        xlim = (efftime.min(), efftime.max())
        fig, ax = plt.subplots(figsize=(9, 6))
        _ = ax.hist(tilestable['EFFTIME_SPEC'] / 60, bins=50, range=xlim,
                    label='All Tiles (N={})'.format(len(tilestable)))
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

        log.info('Writing {}'.format(png))
        fig.savefig(png)
        plt.close()

    return tilestable

def read_parent_sample(samplefile):
    """Read the output of select_parent_sample

    """
    phot = Table(fitsio.read(samplefile, 'FASTPHOT'))
    spec = Table(fitsio.read(samplefile, 'FASTSPEC'))
    meta = Table(fitsio.read(samplefile, 'METADATA'))
    return phot, spec, meta

def read_fastspecfit(tilestable, specprod='denali', targetclass='lrg'):
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
    
    ontiles = np.where(np.isin(meta['TILEID'], tilestable['TILEID']))[0]
    spec = spec[ontiles]
    meta = meta[ontiles]
    phot = phot[ontiles]
    
    log.info('Keeping {} objects on {}/{} unique tiles.'.format(
        len(ontiles), len(np.unique(meta['TILEID'])), len(tilestable)))
    
    ngal = len(spec)

    # correct for extinction!
    # see https://github.com/desihub/fastspecfit/issues/23
    if specprod == 'denali':
        log.warning('Correcting for MW extinction in denali production!')

        from desiutil.dust import SFDMap, ext_odonnell
        from speclite import filters
        
        RV = 3.1
        SFD = SFDMap(scaling=0.86) # SF11 recalibration of the SFD maps
        ebv = SFD.ebv(meta['RA'], meta['DEC'])

        effwave_north, effwave_south = {}, {}
        for band, nfilt, sfilt in zip(['G', 'R', 'Z', 'W1', 'W2'],
                               ['BASS-g', 'BASS-r', 'MzLS-z', 'wise2010-W1', 'wise2010-W2'],
                               ['decam2014-g', 'decam2014-r', 'decam2014-z', 'wise2010-W1', 'wise2010-W2']):
            effwave_north[band] = filters.load_filters(nfilt).effective_wavelengths.value[0]
            effwave_south[band] = filters.load_filters(sfilt).effective_wavelengths.value[0]

    # correct for extinction!
    # see https://github.com/desihub/fastspecfit/issues/23
    def _get_mw_transmission(good, band):
        mw_transmission = np.ones(len(good))
        isouth = np.where(meta['PHOTSYS'][good] == 'S')[0]
        inorth = np.where(meta['PHOTSYS'][good] == 'N')[0]
        if len(isouth) > 0:
            mw_transmission[isouth] = 10**(-0.4 * ebv[good][isouth] * RV * ext_odonnell(effwave_south[band], Rv=RV))
        if len(inorth) > 0:
            mw_transmission[inorth] = 10**(-0.4 * ebv[good][inorth] * RV * ext_odonnell(effwave_south[band], Rv=RV))
        return mw_transmission

    # convenience magnitudes and targeting variables
    for band in ('G', 'R', 'Z', 'W1'):
        phot['{}MAG'.format(band)] = np.zeros(ngal, 'f4')
        good = np.where(meta['FLUX_{}'.format(band)] > 0)[0]
        phot['{}MAG'.format(band)][good] = 22.5 - 2.5 * np.log10(meta['FLUX_{}'.format(band)][good] / _get_mw_transmission(good, band))
        
    for band in ('G', 'R', 'Z'):
        phot['{}FIBERMAG'.format(band)] = np.zeros(ngal, 'f4')
        good = np.where(meta['FIBERFLUX_{}'.format(band)] > 0)[0]
        phot['{}FIBERMAG'.format(band)][good] = 22.5 - 2.5 * np.log10(meta['FIBERFLUX_{}'.format(band)][good] / _get_mw_transmission(good, band))

    targs = ['BGS_ANY', 'ELG', 'LRG', 'QSO']
    targcols = ['BGS', 'ELG', 'LRG', 'QSO']
    for targcol in targcols:
        spec[targcol] = np.zeros(ngal, bool)
        phot[targcol] = np.zeros(ngal, bool)
        
    for tile in tilestable['TILEID']:
        I = np.where(meta['TILEID'] == tile)[0]
        if len(I) == 0:
            continue

        (desicol, bgscol, mwscol), (desimask, bgsmask, mwsmask), survey = main_cmx_or_sv(meta[I])

        for targcol, targ in zip(targcols, targs):
            phot[targcol][I] = meta[desicol][I] & desimask.mask(targ) != 0
            spec[targcol][I] = meta[desicol][I] & desimask.mask(targ) != 0
            
    #for targcol in targcols:
    #    log.info('  {}: {}'.format(targcol, np.sum(phot[targcol])))

    itarg = phot[targetclass.upper()]
    log.info('Keeping {} {} targets.'.format(np.sum(itarg), targetclass.upper()))

    phot = phot[itarg]
    spec = spec[itarg]
    meta = meta[itarg]
    
    return phot, spec, meta

def _select_lrg(iparent, phot, spec, meta, z_minmax=(0.1, 1.1), Mr_minmax=None,
                gi_minmax=None, rW1_minmax=None):
    """Select LRGs. This method should be called with:

    select_parent_sample(phot, spec, meta, targetclass='lrg')

    """
    iselect = copy(iparent)
    if z_minmax:
        iselect *= (meta['Z'] > z_minmax[0]) * (meta['Z'] < z_minmax[1])
    if Mr_minmax:
        iselect *= (phot['ABSMAG_R'] > Mr_minmax[0]) * (phot['ABSMAG_R'] < Mr_minmax[1])
    if gi_minmax:
        gi = phot['ABSMAG_G'] - phot['ABSMAG_I']
        iselect *= (gi > gi_minmax[0]) * (gi < gi_minmax[1])
    if rW1_minmax:
        rW1 = phot['ABSMAG_R'] - phot['ABSMAG_W1']
        iselect *= (rW1 > rW1_minmax[0]) * (rW1 < rW1_minmax[1])
    return iselect

def _select_elg(iparent, phot, spec, meta, z_minmax=(0.6, 1.5),
                Mg_minmax=None, gr_minmax=None):
    """Select ELGs. This method should be called with:

    select_parent_sample(phot, spec, meta, targetclass='elg')

    """
    iselect = copy(iparent)
    if z_minmax:
        iselect *= (meta['Z'] > z_minmax[0]) * (meta['Z'] < z_minmax[1])
    if Mg_minmax:
        iselect *= (phot['ABSMAG_G'] > Mg_minmax[0]) * (phot['ABSMAG_G'] < Mg_minmax[1])
    if gr_minmax:
        gr = phot['ABSMAG_G'] - phot['ABSMAG_R']
        iselect *= (gr > gr_minmax[0]) * (gr < gr_minmax[1])
    return iselect

def _select_bgs(iparent, phot, spec, meta, z_minmax=(0.05, 0.55),
                Mr_minmax=None, gr_minmax=None):
    """Select BGS. This method should be called with:

    select_parent_sample(phot, spec, meta, targetclass='bgs')

    """
    iselect = copy(iparent)
    if z_minmax:
        iselect *= (meta['Z'] > z_minmax[0]) * (meta['Z'] < z_minmax[1])
    if Mr_minmax:
        iselect *= (phot['ABSMAG_R'] > Mr_minmax[0]) * (phot['ABSMAG_R'] < Mr_minmax[1])
    if gr_minmax:
        gr = phot['ABSMAG_G'] - phot['ABSMAG_R']
        iselect *= (gr > gr_minmax[0]) * (gr < gr_minmax[1])
    return iselect

def select_parent_sample(phot, spec, meta, targetclass='lrg', specprod='denali',
                         deltachi2_cut=40, fastphot_chi2cut=100.0,
                         fastspec_chi2cut=3.0, smoothcorr_cut=10,
                         return_indices=False, verbose=False, png=None, samplefile=None,
                         **kwargs):
    """High-level sample selection.

    """
    iparent = (
        (meta['DELTACHI2'] > deltachi2_cut) * 
        (meta['FLUX_G'] > 0) * 
        (meta['FLUX_R'] > 0) * 
        (meta['FLUX_Z'] > 0) * 
        (meta['FLUX_W1'] > 0) *
        (spec['CONTINUUM_CHI2'] < fastspec_chi2cut) *
        (phot['CONTINUUM_CHI2'] < fastphot_chi2cut) 
        #(np.abs(spec['CONTINUUM_SMOOTHCORR_B']) < smoothcorr_cut) *
        #(np.abs(spec['CONTINUUM_SMOOTHCORR_R']) < smoothcorr_cut) *
        #(np.abs(spec['CONTINUUM_SMOOTHCORR_Z']) < smoothcorr_cut)
    )

    if targetclass == 'lrg':
        iselect = _select_lrg(iparent, phot, spec, meta, **kwargs)
    elif targetclass == 'elg':
        iselect = _select_elg(iparent, phot, spec, meta, **kwargs)
    elif targetclass == 'bgs':
        iselect = _select_bgs(iparent, phot, spec, meta, **kwargs)
    else:
        pass

    if verbose:
        log.info('Selecting a parent sample of {}/{} {}s.'.format(
            np.sum(iselect), len(meta), targetclass.upper()))

    # optionally write out
    if samplefile:
        from astropy.io import fits
        log.info('Writing {} objects to {}'.format(np.sum(iselect), samplefile))

        hduprim = fits.PrimaryHDU()
        hduphot = fits.convenience.table_to_hdu(phot[iselect])
        hduphot.header['EXTNAME'] = 'FASTPHOT'

        hduspec = fits.convenience.table_to_hdu(spec[iselect])
        hduspec.header['EXTNAME'] = 'FASTSPEC'
    
        hdumeta = fits.convenience.table_to_hdu(meta[iselect])
        hdumeta.header['EXTNAME'] = 'METADATA'
        hdumeta.header['SPECPROD'] = (specprod, 'spectroscopic production name')
        hdumeta.header['TARGET'] = (targetclass, 'target class')
        
        hx = fits.HDUList([hduprim, hduphot, hduspec, hdumeta])
        hx.writeto(samplefile, overwrite=True, checksum=True)

    # return
    if return_indices:
        return np.where(iselect)[0]
    else:
        return phot[iselect], spec[iselect], meta[iselect]

def stacking_bins(targetclass='lrg', verbose=False):

    # define the stacking limits and the number of bin *centers*

    if targetclass == 'lrg':
        zlim, nz = [0.1, 1.1], 10
        Mrlim, nMr = [-24.5, -20], 9
        gilim, ngi = [0.4, 1.4], 5 
        rW1lim, nrW1 = [-1.0, 1.25], 9
        
        dz = (zlim[1] - zlim[0]) / nz
        dMr = (Mrlim[1] - Mrlim[0]) / nMr
        dgi = (gilim[1] - gilim[0]) / ngi
        drW1 = (rW1lim[1] - rW1lim[0]) / nrW1
        
        # build the array of (left) bin *edges*
        zgrid = np.arange(zlim[0], zlim[1], dz)
        Mrgrid = np.arange(Mrlim[0], Mrlim[1], dMr)
        gigrid = np.arange(gilim[0], gilim[1], dgi)
        rW1grid = np.arange(rW1lim[0], rW1lim[1], drW1)
        
        bins = {'zobj': {'min': zlim[0], 'max': zlim[1], 'del': dz, 'grid': zgrid},
                'Mr': {'min': Mrlim[0], 'max': Mrlim[1], 'del': dMr, 'grid': Mrgrid}, 
                'gi': {'min': gilim[0], 'max': gilim[1], 'del': dgi, 'grid': gigrid}, 
                'rW1': {'min': rW1lim[0], 'max': rW1lim[1], 'del': drW1, 'grid': rW1grid}
                }
    elif targetclass == 'elg':
        zlim, nz = [0.6, 1.5], 9
        Mglim, nMg = [-24, -19], 10
        grlim, ngr = [-0.2, 0.6], 4
        
        dz = (zlim[1] - zlim[0]) / nz
        dMg = (Mglim[1] - Mglim[0]) / nMg
        dgr = (grlim[1] - grlim[0]) / ngr
        
        # build the array of (left) bin *edges*
        zgrid = np.arange(zlim[0], zlim[1], dz)
        Mggrid = np.arange(Mglim[0], Mglim[1], dMg)
        grgrid = np.arange(grlim[0], grlim[1], dgr)
        
        bins = {'zobj': {'min': zlim[0], 'max': zlim[1], 'del': dz, 'grid': zgrid},
                'Mg': {'min': Mglim[0], 'max': Mglim[1], 'del': dMg, 'grid': Mggrid}, 
                'gr': {'min': grlim[0], 'max': grlim[1], 'del': dgr, 'grid': grgrid}, 
                }
    elif targetclass == 'bgs':
        zlim, nz = [0.05, 0.55], 10
        Mrlim, nMr = [-24.0, -17.0], 7
        grlim, ngr = [0.0, 1.0], 5
        
        dz = (zlim[1] - zlim[0]) / nz
        dMr = (Mrlim[1] - Mrlim[0]) / nMr
        dgr = (grlim[1] - grlim[0]) / ngr
        
        # build the array of (left) bin *edges*
        zgrid = np.arange(zlim[0], zlim[1], dz)
        Mrgrid = np.arange(Mrlim[0], Mrlim[1], dMr)
        grgrid = np.arange(grlim[0], grlim[1], dgr)
        
        bins = {'zobj': {'min': zlim[0], 'max': zlim[1], 'del': dz, 'grid': zgrid},
                'Mr': {'min': Mrlim[0], 'max': Mrlim[1], 'del': dMr, 'grid': Mrgrid}, 
                'gr': {'min': grlim[0], 'max': grlim[1], 'del': dgr, 'grid': grgrid}, 
                }
    else:
        raise NotImplemented
        
    nbins = 1
    for key in bins.keys():
        nbins *= len(bins[key]['grid'])
        if verbose:
            log.info(len(bins[key]['grid']))
    if verbose:
        log.info(nbins)

    return bins, nbins

def spectra_in_bins(tilestable, minperbin=3, targetclass='lrg', specprod='denali',
                    minwave=None, maxwave=None, fastphot_in_bins=True, verbose=False):
    """Select objects in bins of rest-frame properties.

    fastphot_in_bins - also stack the fastphot continuum-fitting results
    
    """
    from astropy.table import Table, Column
    
    fastspecdir = os.path.join(os.getenv('FASTSPECFIT_DATA'), specprod, 'tiles')

    if fastphot_in_bins:
        from fastspecfit.continuum import ContinuumFit            
        CFit = ContinuumFit(minwave=minwave, maxwave=maxwave)
        continuumwave = CFit.sspwave
    else:
        continuumwave = None
        
    def sample_template(bins):
        sample1 = Table()
        sample1.add_column(Column(name='IBIN', dtype=np.int32, length=1))
        sample1.add_column(Column(name='NOBJ', dtype=np.int32, length=1))
        for binkey in bins.keys():
            sample1.add_column(Column(name=binkey.upper(), dtype='f4', length=1)) # mean bin center
            sample1.add_column(Column(name='{}MIN'.format(binkey.upper()), dtype='f4', length=1))
            sample1.add_column(Column(name='{}MAX'.format(binkey.upper()), dtype='f4', length=1))
        return sample1

    def _get_data(I):
        _data = {}
        _data['flux'] = flux[I, :]
        _data['ivar'] = ivar[I, :]
        _data['fastphot'] = allphot[I]
        _data['fastspec'] = allspec[I]
        _data['metadata'] = allmeta[I]

        # rebuild the best-fitting continuum model fits
        if fastphot_in_bins:
            _data['cflux'] = []
            for iobj in np.arange(len(I)):
                cflux1, _ = CFit.SSP2data(
                    CFit.sspflux, continuumwave, 
                    redshift=allmeta[I][iobj]['Z'],
                    AV=allphot['CONTINUUM_AV'][I][iobj],
                    coeff=allphot['CONTINUUM_COEFF'][I][iobj],# * CFit.massnorm,
                    synthphot=False)
                _data['cflux'].append(cflux1)# * (1 + allmeta[I[iobj]]['Z']) # deredshift
            _data['cflux'] = np.vstack(_data['cflux'])
        return _data

    samples, data = [], [] # these lists will be aligned
    
    bins, nbins = stacking_bins(targetclass, verbose=True)

    wave = None

    tiles = tilestable['TILEID']
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

        # this is boneheaded...vectorize!
        if targetclass == 'lrg':
            dz, dMr, dgi, drW1 = bins['zobj']['del'], bins['Mr']['del'], bins['gi']['del'], bins['rW1']['del']
            for zmin in bins['zobj']['grid']:
                for Mrmin in bins['Mr']['grid']:
                    for gimin in bins['gi']['grid']:
                        for rW1min in bins['rW1']['grid']:
                            I = select_parent_sample(allphot, allspec, allmeta, 
                                                     z_minmax=[zmin, zmin+dz],
                                                     Mr_minmax=[Mrmin, Mrmin+dMr],
                                                     gi_minmax=[gimin, gimin+dgi],
                                                     rW1_minmax=[rW1min, rW1min+drW1],
                                                     targetclass=targetclass,
                                                     verbose=False, return_indices=True)
                            if len(I) >= minperbin:
                                _sample = sample_template(bins)
                                _sample['IBIN'] = ibin
                                _sample['NOBJ'] = len(I)
                                for key, mmin, delt in zip(('ZOBJ', 'MR', 'GI', 'RW1'),
                                                           (zmin, Mrmin, gimin, rW1min),
                                                           (dz, dMr, dgi, drW1)):
                                    _sample[key] = mmin + delt / 2
                                    _sample['{}MIN'.format(key)] = mmin
                                    _sample['{}MAX'.format(key)] = mmin + delt
                                samples.append(_sample)
                                _data = _get_data(I)
                                data.append(_data)
                            ibin += 1 # next bin
        elif targetclass == 'elg':
            dz, dMg, dgr = bins['zobj']['del'], bins['Mg']['del'], bins['gr']['del']
            for zmin in bins['zobj']['grid']:
                for Mgmin in bins['Mg']['grid']:
                    for grmin in bins['gr']['grid']:
                        I = select_parent_sample(allphot, allspec, allmeta, 
                                                 z_minmax=[zmin, zmin+dz],
                                                 Mg_minmax=[Mgmin, Mgmin+dMg],
                                                 gr_minmax=[grmin, grmin+dgr],
                                                 targetclass=targetclass,
                                                 verbose=False, return_indices=True)
                        if len(I) >= minperbin:
                            _sample = sample_template(bins)
                            _sample['IBIN'] = ibin
                            _sample['NOBJ'] = len(I)
                            for key, mmin, delt in zip(('ZOBJ', 'MG', 'GR'),
                                                       (zmin, Mgmin, grmin),
                                                       (dz, dMg, dgr)):
                                _sample[key] = mmin + delt / 2
                                _sample['{}MIN'.format(key)] = mmin
                                _sample['{}MAX'.format(key)] = mmin + delt
                            samples.append(_sample)
                            _data = _get_data(I)
                            data.append(_data)
                        ibin += 1 # next bin
        elif targetclass == 'bgs':
            dz, dMr, dgr = bins['zobj']['del'], bins['Mr']['del'], bins['gr']['del']
            for zmin in bins['zobj']['grid']:
                for Mrmin in bins['Mr']['grid']:
                    for grmin in bins['gr']['grid']:
                        I = select_parent_sample(allphot, allspec, allmeta, 
                                                 z_minmax=[zmin, zmin+dz],
                                                 Mr_minmax=[Mrmin, Mrmin+dMr],
                                                 gr_minmax=[grmin, grmin+dgr],
                                                 targetclass=targetclass,
                                                 verbose=False, return_indices=True)
                        if len(I) >= minperbin:
                            _sample = sample_template(bins)
                            _sample['IBIN'] = ibin
                            _sample['NOBJ'] = len(I)
                            for key, mmin, delt in zip(('ZOBJ', 'MR', 'GR'),
                                                       (zmin, Mrmin, grmin),
                                                       (dz, dMr, dgr)):
                                _sample[key] = mmin + delt / 2
                                _sample['{}MIN'.format(key)] = mmin
                                _sample['{}MAX'.format(key)] = mmin + delt
                            samples.append(_sample)
                            _data = _get_data(I)
                            data.append(_data)
                        ibin += 1 # next bin
        else:
            raise NotImplemented

    # Stack the bin-level statistics table
    samples = Table(np.hstack(samples))
    data = np.array(data)

    # ...and now stack the data in each (unique) bin number.
    samplestack, sampledata = [], {}
    nbin = len(set(samples['IBIN']))

    for count, ibin in enumerate(sorted(set(samples['IBIN']))):
        log.info('Stacking spectra from bin {}/{}'.format(count, nbin))
        I = np.where(ibin == samples['IBIN'])[0]
        _samplestack = samples[[I[0]]].copy()
        _samplestack['NOBJ'] = np.sum(samples[I]['NOBJ'])
        samplestack.append(_samplestack)

        _sampledata = {}
        _sampledata['flux'] = np.vstack([_data['flux'] for _data in data[I]])
        _sampledata['ivar'] = np.vstack([_data['ivar'] for _data in data[I]])
        _sampledata['fastphot'] = Table(np.hstack([_data['fastphot'] for _data in data[I]]))
        _sampledata['fastspec'] = Table(np.hstack([_data['fastspec'] for _data in data[I]]))
        _sampledata['metadata'] = Table(np.hstack([_data['metadata'] for _data in data[I]]))

        if fastphot_in_bins:
            _sampledata['cflux'] = np.vstack([_data['cflux'] for _data in data[I]])
        
        sampledata.update({str(ibin): _sampledata})
        del _sampledata, _samplestack

    del data, samples
    samplestack = Table(np.hstack(samplestack))
    
    return samplestack, sampledata, wave, continuumwave




