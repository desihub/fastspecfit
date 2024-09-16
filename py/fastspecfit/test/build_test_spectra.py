"""Build mini spectral dataset.

"""
import os, pdb
import numpy as np
from pathlib import Path
import fitsio
from astropy.table import Table
from redrock.external.desi import write_zbest
from desispec.io import write_spectra, read_spectra, photo#, read_tile_spectra

specprod = 'iron'
os.environ['SPECPROD'] = specprod # needed to get write_spectra have the correct dependency
night = 20210324
tile = 80613
petal = 4
targetid = 39633345008634465

dr9dir = os.path.join(os.environ.get('DESI_ROOT'), 'external', 'legacysurvey', 'dr9')
reduxdir = os.path.join(os.environ.get('DESI_ROOT'), 'spectro', 'redux', specprod)
datadir = os.path.join(reduxdir, 'tiles', 'cumulative', str(tile), str(night))

coaddfile = 'coadd-{}-{}-thru{}.fits'.format(petal, tile, night)
redrockfile = coaddfile.replace('coadd-', 'redrock-')

#spec, redrock = read_tile_spectra(tile, night, specprod=specprod, targets=targetid,
#                                  group='cumulative', redrock=True, coadd=True)

coaddfile = os.path.join(datadir, coaddfile)
redrockfile = os.path.join(datadir, redrockfile)

redhdr = fitsio.read_header(redrockfile)
zbest = Table.read(redrockfile, 'REDSHIFTS')
fibermap = Table.read(redrockfile, 'FIBERMAP')
expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
tsnr2 = Table.read(redrockfile, 'TSNR2')

spechdr = fitsio.read_header(coaddfile)

zbest = zbest[np.isin(zbest['TARGETID'], targetid)]
fibermap = fibermap[np.isin(fibermap['TARGETID'], targetid)]
expfibermap = expfibermap[np.isin(expfibermap['TARGETID'], targetid)]
tsnr2 = tsnr2[np.isin(tsnr2['TARGETID'], targetid)]

archetype_version = None
template_version = {redhdr['TEMNAM{:02d}'.format(nn)]: redhdr['TEMVER{:02d}'.format(nn)] for nn in np.arange(10)}

print('Writing {}'.format(os.path.join('data', os.path.basename(redrockfile))))
write_zbest(os.path.join('data', os.path.basename(redrockfile)),
            zbest, fibermap, expfibermap, tsnr2, template_version,
            archetype_version, spec_header=spechdr)

spec = read_spectra(coaddfile).select(targets=targetid)
print('Writing {}'.format(os.path.join('data', os.path.basename(coaddfile))))
write_spectra(os.path.join('data', os.path.basename(coaddfile)), spec)

# write out a little tiles file
tilesfile = os.path.join(reduxdir, 'tiles-{}.csv'.format(specprod))
tiles = Table.read(tilesfile)
tiles = tiles[tiles['TILEID'] == tile]
outfile = os.path.join('data', os.path.basename(tilesfile))
print('Writing {}'.format(outfile))
tiles.write(filename=outfile, format='csv', overwrite=True)

## gather Tractor photometry
#tractor = photo.gather_tractorphot(fibermap, dr9dir=dr9dir)
#tractor.remove_columns(('TARGETID', 'LS_ID'))
#
#region = 'north'
#brick = tractor['BRICKNAME'][0] # need the header, too
#tractorfile = os.path.join(dr9dir, region, 'tractor', brick[:3], 'tractor-{}.fits'.format(brick))
#tractorhdr = fitsio.read_header(tractorfile)
#outtractorfile = os.path.join('data', region, 'tractor', brick[:3], os.path.basename(tractorfile))
#
#print('Writing {}'.format(outtractorfile))
#os.makedirs(os.path.dirname(outtractorfile), exist_ok=True)
#fitsio.write(outtractorfile, tractor.as_array(), header=tractorhdr, clobber=True)
