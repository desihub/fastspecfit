"""Build mini spectral dataset.

"""
import os, pdb
import numpy as np
import fitsio
from astropy.table import Table
from redrock.external.desi import write_zbest
from desispec.io import write_spectra, read_spectra#, read_tile_spectra

specprod = 'fuji'
night = 20210324
tile = 80613
petal = 4
targetid = 39633345008634465

datadir = os.path.join(os.environ.get('DESI_ROOT'), 'spectro', 'redux', specprod, 'tiles', 'cumulative', str(tile), str(night))

coaddfile = 'coadd-{}-{}-thru{}.fits'.format(petal, tile, night)
redrockfile = coaddfile.replace('coadd-', 'redrock-')

#spec, redrock = read_tile_spectra(tile, night, specprod=specprod, targets=targetid, 
#                                  group='cumulative', redrock=True, coadd=True)

coaddfile = os.path.join(datadir, coaddfile)
redrockfile = os.path.join(datadir, redrockfile)

spechdr = fitsio.read_header(coaddfile)
spec = read_spectra(coaddfile).select(targets=targetid)
print('Writing {}'.format(os.path.join('data', os.path.basename(coaddfile))))
write_spectra(os.path.join('data', os.path.basename(coaddfile)), spec)

redhdr = fitsio.read_header(redrockfile)
zbest = Table.read(redrockfile, 'REDSHIFTS')
fibermap = Table.read(redrockfile, 'FIBERMAP')
expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
tsnr2 = Table.read(redrockfile, 'TSNR2')

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









