"""
fastspecfit.templates.sample
============================

Code for defining samples and reading the data needed to build templates for the
various target classes. Called by bin/desi-templates.

"""
import pdb # for debugging

import os
import numpy as np
import fitsio
from astropy.table import Table, Column
    
from desiutil.log import get_logger
log = get_logger()

VITILES_TARGETCLASS = {'lrg': [80605, 80609],
                       'elg': [80606, 80608, 80610],
                       'bgs': [80613]}

SAMPLE_PROPERTIES = {
    'lrg': {'zminmax': (0.1, 1.1), 'normwave': 4500.0, 'absmag': 'Mr', 'color': 'rW1',
            'absmag_band': 'R', 'color_band1': 'R', 'color_band2': 'W1',
            'absmag_label': 'M_{{0.0r}}', 'color_label': '^{{0.0}}(r-W1)'},            
    'elg': {'zminmax': (0.6, 1.5), 'normwave': 3500.0, 'absmag': 'Mg', 'color': 'gr',
            'absmag_band': 'G', 'color_band1': 'G', 'color_band2': 'R',
            'absmag_label': 'M_{{0.0g}}', 'color_label': '^{{0.0}}(g-r)'},
    'bgs': {'zminmax': (0.05, 1.55), 'normwave': 5500.0, 'absmag': 'Mr', 'color': 'gr',
            'absmag_band': 'R', 'color_band1': 'G', 'color_band2': 'R',
            'absmag_label': 'M_{{0.0r}}', 'color_label': '^{{0.0}}(g-r)'},
    }

def select_tiles(targetclass, remove_vi=False, min_efftime=None,
                 specprod='denali', outfile=None, png=None):
    """Select tiles to use. Remove low exposure-time tiles and also 
    remove tiles which are being visually inspected.
    
      /global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/

    """
    reduxdir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'redux', specprod)
    tilestable = Table.read(os.path.join(reduxdir, 'tiles-{}.csv'.format(specprod)))
    tilestable = tilestable[tilestable['SURVEY'] == 'sv1']
    #tilestable = tilestable[tilestable['SURVEY'] != 'unknown']
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
        #log.info('Writing {} tiles to {}'.format(len(targtiles), outfile))
        #targtiles.write(outfile, overwrite=True)
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

    #return targtiles
    return tilestable

def read_parent_sample(samplefile):
    """Read the output of select_parent_sample

    """
    log.info('Reading {}'.format(samplefile))
    phot = Table(fitsio.read(samplefile, 'FASTPHOT'))
    spec = Table(fitsio.read(samplefile, 'FASTSPEC'))
    meta = Table(fitsio.read(samplefile, 'METADATA'))
    return phot, spec, meta

def read_tilestable(tilefile):
    """Read the output of select_tiles

    """
    tilestable = Table.read(tilefile)
    log.info('Read {} tiles from {}'.format(len(tilestable), tilefile))
    return tilestable

def read_fastspecfit(tilestable, fastspecfit_dir=None, specprod='denali',
                     targetclass='lrg'):
    """Read the fastspecfit output for this production.

    """
    from desitarget.targets import main_cmx_or_sv

    if fastspecfit_dir is None:
        fastspecfit_dir = os.path.join(os.getenv('DESI_ROOT'), 'spectro', 'fastspecfit')

    fastspecdir = os.path.join(fastspecfit_dir, specprod, 'tiles')
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

    # for convenience, add observed-frame photometry, rest-frame colors, and
    # targeting variables
    for band in ('G', 'R', 'Z', 'W1'):
        phot['{}MAG'.format(band)] = np.zeros(ngal, 'f4')
        good = np.where(meta['FLUX_{}'.format(band)] > 0)[0]
        phot['{}MAG'.format(band)][good] = 22.5 - 2.5 * np.log10(meta['FLUX_{}'.format(band)][good] / _get_mw_transmission(good, band))
        
    for band in ('G', 'R', 'Z'):
        phot['{}FIBERMAG'.format(band)] = np.zeros(ngal, 'f4')
        good = np.where(meta['FIBERFLUX_{}'.format(band)] > 0)[0]
        phot['{}FIBERMAG'.format(band)][good] = 22.5 - 2.5 * np.log10(meta['FIBERFLUX_{}'.format(band)][good] / _get_mw_transmission(good, band))

    #phot['ABSMAG_GR'] = phot['ABSMAG_G'] - phot['ABSMAG_R']
    #phot['ABSMAG_RZ'] = phot['ABSMAG_R'] - phot['ABSMAG_Z']
    #phot['ABSMAG_RW1'] = phot['ABSMAG_R'] - phot['ABSMAG_W1']

    # targeting...
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

def select_parent_sample(phot, spec, meta, targetclass='lrg', specprod='denali',
                         deltachi2_cut=40, fastphot_chi2cut=100.0,
                         fastspec_chi2cut=3.0, smoothcorr_cut=10,
                         zobj_minmax=None, absmag_minmax=None, color_minmax=None,
                         return_indices=False, verbose=False, png=None, samplefile=None):
    """High-level sample selection.

    """
    iparent = (
        (meta['Z'] > zobj_minmax[0]) *
        (meta['Z'] < zobj_minmax[1]) * 
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

    if zobj_minmax is not None and absmag_minmax is not None and color_minmax is not None:

        props = SAMPLE_PROPERTIES[targetclass]
        absmagcol = 'ABSMAG_{}'.format(props['absmag_band'])
        color = phot['ABSMAG_{}'.format(props['color_band1'])] - phot['ABSMAG_{}'.format(props['color_band2'])]
        
        iselect = iparent * (
            (meta['Z'] > zobj_minmax[0]) * (meta['Z'] < zobj_minmax[1]) *
            (phot[absmagcol] > absmag_minmax[0]) * (phot[absmagcol] < absmag_minmax[1]) *
            (color > color_minmax[0]) * (color < color_minmax[1]) )
    else:
        iselect = iparent

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
        absmaglim, nabsmag = [-24.5, -20], 9 # Mr
        colorlim, ncolor = [-1.0, 1.25], 9   # r-W1
    elif targetclass == 'elg':
        zlim, nz = [0.6, 1.5], 9
        absmaglim, nabsmag = [-24, -19], 10 # Mg
        colorlim, ncolor = [-0.2, 0.6], 4   # g-r
    elif targetclass == 'bgs':
        zlim, nz = [0.05, 0.55], 10
        absmaglim, nabsmag = [-24.0, -17.0], 7 # Mr
        colorlim, ncolor = [0.0, 1.0], 5       # g-r
    else:
        raise NotImplemented
        
    dz = (zlim[1] - zlim[0]) / nz
    dabsmag = (absmaglim[1] - absmaglim[0]) / nabsmag
    dcolor = (colorlim[1] - colorlim[0]) / ncolor

    # build the array of (left) bin *edges*
    zgrid = np.arange(zlim[0], zlim[1], dz)
    absmaggrid = np.arange(absmaglim[0], absmaglim[1], dabsmag)
    colorgrid = np.arange(colorlim[0], colorlim[1], dcolor)

    nbins = len(zgrid) * len(absmaggrid) * len(colorgrid)

    # pack into a table
    bins = Table()
    bins.add_column(Column(name='TARGETCLASS', dtype='U3', length=nbins))
    bins.add_column(Column(name='IBIN', dtype=np.int32, length=nbins))
    bins.add_column(Column(name='ISUBBIN', dtype=np.int16, length=nbins))
    bins.add_column(Column(name='NOBJ', dtype=np.int32, length=nbins))
    bins.add_column(Column(name='SNR', dtype='f4', length=nbins))
    for col in ('ZOBJ', 'ABSMAG', 'COLOR'):
        bins.add_column(Column(name=col, dtype='f4', length=nbins)) # mean bin center
        bins.add_column(Column(name='{}MIN'.format(col), dtype='f4', length=nbins))
        bins.add_column(Column(name='{}MAX'.format(col), dtype='f4', length=nbins))
        
    bins['TARGETCLASS'] = targetclass
    bins['IBIN'] = np.arange(nbins, dtype=np.int32)

    ibin = 0
    for zmin in zgrid:
        for absmagmin in absmaggrid:
            for colormin in colorgrid:
                for col, mmin, delt in zip(('ZOBJ', 'ABSMAG', 'COLOR'),
                                           (zmin, absmagmin, colormin),
                                           (dz, dabsmag, dcolor)):
                    bins[col][ibin] = mmin + delt / 2             # bin center
                    bins['{}MIN'.format(col)][ibin] = mmin        # left edge
                    bins['{}MAX'.format(col)][ibin] = mmin + delt # right edge
                ibin += 1

    if verbose:
        log.info('Number of {} bins = {}'.format(targetclass, bins))

    return bins
