"""Build miniature template set.

"""
import os, pdb
import numpy as np
import fitsio
from astropy.table import Table
from astropy.io import fits

os.environ['FASTSPECFIT_TEMPLATES'] = os.path.join(os.environ.get('DESI_ROOT'), 'science', 'gqp', 'templates', 'SSP-CKC14z')

templates_dir = os.environ.get('FASTSPECFIT_TEMPLATES', '.')
sspfile = os.path.join(templates_dir, 'SSP_{}_{}_{}_{}.fits'.format(
    isochrone, library, imf, metallicity))
if not os.path.isfile(sspfile):
    raise IOError('SSP templates file not found {}'.format(sspfile))

outfile = os.path.join('data', 'SSP_{}_{}_{}_{}.fits'.format(
    isochrone, library, imf, metallicity))

print('Reading {}'.format(sspfile))
wave, wavehdr = fitsio.read(sspfile, 'WAVE', header=True)
flux, fluxhdr = fitsio.read(sspfile, 'FLUX', header=True)
meta, metahdr = fitsio.read(sspfile, 'METADATA', header=True)
meta = Table(meta)
    
# Trim the wavelengths and select the number/ages of the templates.
# https://www.sdss.org/dr14/spectro/galaxy_mpajhu
minwave, maxwave = np.min(wave), 6e4
keep = np.where((wave >= minwave) * (wave <= maxwave))[0]
wave = wave[keep]

myages = np.array([0.005, 0.025, 0.1, 0.2, 0.6, 0.9, 1.4, 2.5, 5, 10.0, 13.0])*1e9
iage = np.array([np.argmin(np.abs(meta['age']-myage)) for myage in myages])
flux = flux[:, iage][keep, :] # flux[keep, ::5]
meta = meta[iage]

fitsio.write(outfile, flux, 'FLUX', header=fluxhdr, clobber=True)
fitsio.write(outfile, wave, 'WAVE', header=wavehdr)
fitsio.write(outfile, meta.as_array(), 'METADATA', header=metahdr)

