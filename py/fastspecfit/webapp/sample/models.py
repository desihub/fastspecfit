#!/usr/bin/env python

"""Each model will be written as a class here, instantiated and populated by
load.py, with each model stored as a table in the database and the fields stored
as columns.

"""
import os
import numpy as np
from django.db.models import (Model, IntegerField, CharField, FloatField, IPAddressField,
                              DateTimeField, ManyToManyField, TextField, BooleanField)

# python manage.py makemigrations SGA
# python manage.py migrate

class Sample(Model):
    """Model to represent a single galaxy.

    """
    # in FITS table
    row_index = IntegerField(default=-1)

    sga_id = IntegerField(null=True)
    galaxy = CharField(max_length=30, default='')
    morphtype = CharField(max_length=30, default='')
    ra_leda = FloatField(null=True)
    dec_leda = FloatField(null=True)
    d25_leda = FloatField(default=0.0)
    pa_leda = FloatField(default=0.0)
    ba_leda = FloatField(default=0.0)
    pgc = IntegerField(null=True)
    
    d26 = FloatField(default=0.0)
    pa = FloatField(default=0.0)
    ba = FloatField(default=0.0)
    sma_moment = FloatField(default=0.0)
    ra = FloatField(null=True)
    dec = FloatField(null=True)
    ra_moment = FloatField(null=True)
    dec_moment = FloatField(null=True)

    group_id = IntegerField(null=True)
    group_name = CharField(max_length=40, default='')
    nice_group_name = CharField(max_length=40, default='')
    group_ra = FloatField(null=True)
    group_dec = FloatField(null=True)
    group_diameter = FloatField(default=0.0)
    group_primary = BooleanField(default=False)

    sma_sb24 = FloatField(null=True)
    sma_sb25 = FloatField(null=True)
    sma_sb26 = FloatField(null=True)

    g_sma50 = FloatField(null=True)
    r_sma50 = FloatField(null=True)
    z_sma50 = FloatField(null=True)

    g_mag_sb24 = FloatField(null=True)
    g_mag_sb25 = FloatField(null=True)
    g_mag_sb26 = FloatField(null=True)

    r_mag_sb24 = FloatField(null=True)
    r_mag_sb25 = FloatField(null=True)
    r_mag_sb26 = FloatField(null=True)

    z_mag_sb24 = FloatField(null=True)
    z_mag_sb25 = FloatField(null=True)
    z_mag_sb26 = FloatField(null=True)

    g_cog_params_mtot = FloatField(null=True)
    r_cog_params_mtot = FloatField(null=True)
    z_cog_params_mtot = FloatField(null=True)

    tractortype = CharField(max_length=3, default='')
    sersic = FloatField(null=True)
    shape_r = FloatField(null=True)
    shape_e1 = FloatField(null=True)
    shape_e2 = FloatField(null=True)

    flux_g = FloatField(null=True)
    flux_r = FloatField(null=True)
    flux_z = FloatField(null=True)

    flux_ivar_g = FloatField(null=True)
    flux_ivar_r = FloatField(null=True)
    flux_ivar_z = FloatField(null=True)

    # radec2xyz, for cone search in the database
    ux = FloatField(default=-2.0)
    uy = FloatField(default=-2.0)
    uz = FloatField(default=-2.0)

    def base_html_dir(self):
        return '/global/cfs/cdirs/cosmo/data/sga/2020/html/'

    def png_base_url(self):
        baseurl = 'https://portal.nersc.gov/project/cosmo/data/sga/2020/html/'
        baseurl += self.ra_slice() + '/' + self.group_name + '/';
        return baseurl

    def data_base_url(self):
        baseurl = 'https://portal.nersc.gov/project/cosmo/data/sga/2020/data/'
        baseurl += self.ra_slice() + '/' + self.group_name + '/';
        return baseurl

    def hyperleda_html(self):
        if self.pgc > -1:
            url = 'http://leda.univ-lyon1.fr/ledacat.cgi?o=PGC{}'.format(self.pgc)
        elif 'DR8' in self.galaxy:
            url = 'http://leda.univ-lyon1.fr/ledacat.cgi?{}&ob=ra'.format(self.galaxy.strip())
        else:
            url = 'http://leda.univ-lyon1.fr/fG.cgi?n=a000&c=o&p={}%20{}&f=0.1&ob=ra'.format(
                self.ra_leda, self.dec_leda)
        return url

    def mosaic_diam(self):
        if self.group_diameter > 30: # NGC0598=M33 is 61 arcmin in diameter!
            mosaic_diam = self.group_diameter * 2 * 0.7 # [arcmin]
        elif self.group_diameter > 14 and self.group_diameter < 30:
            mosaic_diam = self.group_diameter * 2 * 1.0 # [arcmin]
        else:
            mosaic_diam = self.group_diameter * 2 * 1.5 # [arcmin]
        return '{:.3f}'.format(mosaic_diam) # [arcmin]

    def ra_slice(self):
        raslice = '{:06d}'.format(int(self.group_ra*1000))[:3]
        return raslice

    def sga_id_string(self):
        return '{}'.format(self.sga_id)

    def group_ra_string(self):
        return '{:.7f}'.format(self.group_ra)

    def group_dec_string(self):
        return '{:.7f}'.format(self.group_dec)

    def ra_leda_string(self):
        return '{:.7f}'.format(self.ra_leda)
    def dec_leda_string(self):
        return '{:.7f}'.format(self.dec_leda)

    def ra_string(self):
        return '{:.7f}'.format(self.ra)
        #return '{:.7f}'.format(self.ra_moment)
    def dec_string(self):
        return '{:.7f}'.format(self.dec)
        #return '{:.7f}'.format(self.dec_moment)

    def group_id_string(self):
        return '{}'.format(self.group_id)

    def group_diameter_string(self):
        return '{:.3f}'.format(self.group_diameter)
    
    def sersic_string(self):
        return '{:.2f}'.format(self.sersic)

    def shape_r_string(self):
        return '{:.3f}'.format(self.shape_r)

    def pa_leda_string(self):
        return '{:.1f}'.format(self.pa_leda)

    def eps_leda_string(self):
        return '{:.3f}'.format(1-self.ba_leda)

    def r25_leda_string(self):
        return '{:.3f}'.format(self.d25_leda * 60 / 2)

    def pa_string(self):
        return '{:.1f}'.format(self.pa)

    def eps_string(self):
        return '{:.3f}'.format(1-self.ba)

    def sma_moment_string(self):
        return '{:.3f}'.format(self.sma_moment)

    def g_sma50_string(self):
        if self.g_sma50 < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.g_sma50)

    def r_sma50_string(self):
        if self.r_sma50 < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.r_sma50)

    def z_sma50_string(self):
        if self.z_sma50 < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.z_sma50)

    def sma_sb24_string(self):
        if self.sma_sb24 < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.sma_sb24)

    def sma_sb25_string(self):
        if self.sma_sb25 < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.sma_sb25)

    def sma_sb26_string(self):
        if self.sma_sb26 < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.sma_sb26)

    def gmag_sb24(self):
        if self.g_mag_sb24 > 0:
            return '{:.3f}'.format(self.g_mag_sb24)
        else:
            return '...'
    def rmag_sb24(self):
        if self.r_mag_sb24 > 0:
            return '{:.3f}'.format(self.r_mag_sb24)
        else:
            return '...'
    def zmag_sb24(self):
        if self.z_mag_sb24 > 0:
            return '{:.3f}'.format(self.z_mag_sb24)
        else:
            return '...'

    def gmag_sb25(self):
        if self.g_mag_sb25 > 0:
            return '{:.3f}'.format(self.g_mag_sb25)
        else:
            return '...'
    def rmag_sb25(self):
        if self.r_mag_sb25 > 0:
            return '{:.3f}'.format(self.r_mag_sb25)
        else:
            return '...'
    def zmag_sb25(self):
        if self.z_mag_sb25 > 0:
            return '{:.3f}'.format(self.z_mag_sb25)
        else:
            return '...'

    def gmag_sb26(self):
        if self.g_mag_sb26 > 0:
            return '{:.3f}'.format(self.g_mag_sb26)
        else:
            return '...'
    def rmag_sb26(self):
        if self.r_mag_sb26 > 0:
            return '{:.3f}'.format(self.r_mag_sb26)
        else:
            return '...'
    def zmag_sb26(self):
        if self.z_mag_sb26 > 0:
            return '{:.3f}'.format(self.z_mag_sb26)
        else:
            return '...'

    def cog_gtot(self):
        if self.g_cog_params_mtot < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.g_cog_params_mtot)
        
    def cog_rtot(self):
        if self.r_cog_params_mtot < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.r_cog_params_mtot)
        
    def cog_ztot(self):
        if self.z_cog_params_mtot < 0:
            return '...'
        else:
            return '{:.3f}'.format(self.z_cog_params_mtot)
        
    def tractor_pa_string(self):
        pa = 180 - (-np.rad2deg(np.arctan2(self.shape_e2, self.shape_e1) / 2))
        pa = pa % 180
        return '{:.1f}'.format(pa)
    def tractor_eps_string(self):
        ee = np.hypot(self.shape_e1, self.shape_e2)
        ba = (1 - ee) / (1 + ee)
        return '{:.3f}'.format(1-ba)
    def tractor_gmag(self):
        if self.flux_g > 0:
            return '{:.3f}'.format(22.5-2.5*np.log10(self.flux_g))
        else:
            return '...'
    def tractor_rmag(self):
        if self.flux_r > 0:
            return '{:.3f}'.format(22.5-2.5*np.log10(self.flux_r))
        else:
            return '...'
    def tractor_zmag(self):
        if self.flux_z > 0:
            return '{:.3f}'.format(22.5-2.5*np.log10(self.flux_z))
        else:
            return '...'

    def ellipsefile(self):
        ellipsefile = '{}{}-largegalaxy-{}-ellipse-sbprofile.png'.format(self.png_base_url(), self.group_name, self.sga_id_string())
        return ellipsefile
    
    def ellipse_exists(self):
        ellipsefile = os.path.join(self.base_html_dir(), self.ra_slice(), self.group_name, '{}-largegalaxy-{}-ellipse-sbprofile.png'.format(
            self.group_name, self.sga_id_string()))
        return os.path.isfile(ellipsefile)
