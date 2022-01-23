#!/usr/bin/env python

"""Load the input sample into a database table.

"""
import os
import numpy as np
import fitsio
import django

from astropy.table import Table, hstack
from astrometry.util.starutil_numpy import radectoxyz

def main():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "SGA.webapp.settings")
    django.setup()

    from SGA.webapp.sample.models import Sample

    DATADIR = '/global/cfs/cdirs/cosmo/data/sga/2020'
    #DATADIR = '/global/cfs/cdirs/cosmo/work/legacysurvey/sga/2020'

    sgafile = os.path.join(DATADIR, 'SGA-2020.fits')

    sga_columns = ['sga_id', 'galaxy', 'morphtype',
                   'ra', 'dec',
                   'ra_leda', 'dec_leda',
                   'd25_leda', 'pa_leda', 'ba_leda', 'pgc',
                   'ra_moment', 'dec_moment', 'd26', 'pa', 'ba', 'sma_moment',
                   'group_id', 'group_name', 'group_ra', 'group_dec', 'group_diameter', 'group_primary',
                   'g_sma50', 'r_sma50', 'z_sma50', 
                   'sma_sb24', 'sma_sb25', 'sma_sb26',
                   'g_mag_sb24', 'g_mag_sb25', 'g_mag_sb26',
                   'r_mag_sb24', 'r_mag_sb25', 'r_mag_sb26',
                   'z_mag_sb24', 'z_mag_sb25', 'z_mag_sb26',
                   'g_cog_params_mtot', 'r_cog_params_mtot', 'z_cog_params_mtot', 
                   ]
                   
    tractor_cols = ['type', 'sersic', 'shape_r', 'shape_e1', 'shape_e2',
                    'flux_g', 'flux_r', 'flux_z', 'flux_ivar_g', 'flux_ivar_r', 'flux_ivar_z']
       
    sga = Table(fitsio.read(sgafile, ext='ELLIPSE', columns=sga_columns))
    sga_tractor = Table(fitsio.read(sgafile, ext='TRACTOR', columns=tractor_cols))
    sga = hstack((sga, sga_tractor))
    print('Read {} rows from {}'.format(len(sga), sgafile))

    sga.rename_column('TYPE', 'TRACTORTYPE')
    sga['NICE_GROUP_NAME'] = [gname.replace('_GROUP', ' Group') for gname in sga['GROUP_NAME']]

    print(sga.colnames)

    xyz = radectoxyz(sga['RA'], sga['DEC'])
    #xyz = radectoxyz(sga['RA_MOMENT'], sga['DEC_MOMENT'])
    #xyz = radectoxyz(sga['RA_LEDA'], sga['DEC_LEDA'])

    objs = []
    nextpow = 1024
    for ii, onegal in enumerate(sga):
        if ii == nextpow:
            print('Row', ii)
            nextpow *= 2

        sam = Sample()
        sam.row_index = ii
        sam.ux = xyz[ii, 0]
        sam.uy = xyz[ii, 1]
        sam.uz = xyz[ii, 2]

        for col in sga.colnames:
            val = onegal[col]
            if type(val) == np.str or type(val) == np.str_:
                val.strip()
            setattr(sam, col.lower(), val)

        objs.append(sam)
            
    print('Bulk creating the database.')
    Sample.objects.bulk_create(objs)

if __name__ == '__main__':
    main()
