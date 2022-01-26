#!/usr/bin/env python

"""Load the input sample into a database table.

"""
import os
import numpy as np
import fitsio
import django

from astropy.table import Table, hstack
#from astrometry.util.starutil_numpy import radectoxyz

# RA, Dec in degrees: scalars or 1-d arrays.
# returns xyz of shape (N,3)
def radectoxyz(ra_deg, dec_deg):
    ra  = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cosd = np.cos(dec)
    xyz = np.vstack((cosd * np.cos(ra),
                  cosd * np.sin(ra),
                  np.sin(dec))).T
    assert(xyz.shape[1] == 3)
    return xyz

def main():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fastspecfit.webapp.settings")
    django.setup()

    from fastspecfit.webapp.sample.models import Sample

    DATADIR = '/global/cfs/cdirs/desi/spectro/fastspecfit/everest/catalogs'
    #DATADIR = '/global/cfs/cdirs/cosmo/work/legacysurvey/sga/2020'

    fastspecfile = os.path.join(DATADIR, 'fastspec-everest-sv3-bright.fits')

    meta_columns = ['TARGETID', 'RA', 'DEC']
    fastspec_cols = ['CONTINUUM_Z']
       
    meta = Table(fitsio.read(fastspecfile, ext='METADATA', columns=meta_columns))
    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC', columns=fastspec_cols))
    fast = hstack((meta, fastspec))
    print('Read {} rows from {}'.format(len(fast), fastspecfile))

    #sga.rename_column('TYPE', 'TRACTORTYPE')
    #sga['NICE_GROUP_NAME'] = [gname.replace('_GROUP', ' Group') for gname in sga['GROUP_NAME']]

    print(fast.colnames)

    xyz = radectoxyz(fast['RA'], fast['DEC'])

    objs = []
    nextpow = 1024
    for ii, onegal in enumerate(fast):
        if ii == nextpow:
            print('Row', ii)
            nextpow *= 2

        sam = Sample()
        sam.row_index = ii
        sam.ux = xyz[ii, 0]
        sam.uy = xyz[ii, 1]
        sam.uz = xyz[ii, 2]

        for col in fast.colnames:
            val = onegal[col]
            if type(val) == np.str or type(val) == np.str_:
                val.strip()
            setattr(sam, col.lower(), val)

        objs.append(sam)
            
    print('Bulk creating the database.')
    Sample.objects.bulk_create(objs)

if __name__ == '__main__':
    main()
