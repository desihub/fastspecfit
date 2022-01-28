#!/usr/bin/env python

"""Each model will be written as a class here, instantiated and populated by
load.py, with each model stored as a table in the database and the fields stored
as columns.

"""
import os
import numpy as np
from django.db.models import (Model, IntegerField, CharField, FloatField, IPAddressField,
                              DateTimeField, ManyToManyField, TextField, BooleanField)

# python manage.py makemigrations sample
# python manage.py migrate

class Sample(Model):
    """Model to represent a single object.

    """
    #def __str__(self):
    #    return 'Sample '+self.target_name
    
    # in FITS table
    row_index = IntegerField(default=-1)

    # target(ing) properties
    specprod = CharField(max_length=15, default='')
    targetid = IntegerField(null=True)
    ra = FloatField(null=True)
    dec = FloatField(null=True)
    survey = CharField(max_length=4, default='')
    faprgrm = CharField(max_length=6, default='')
    hpxpixel = IntegerField(null=True)
    #tileid = IntegerField(null=True)
    #fiber = IntegerField(null=True)
    
    # redrock
    z = FloatField(null=True)
    zwarn = IntegerField(null=True)

    # continuum properties
    continuum_z = FloatField(null=True)

    # emission lines
    oii_3726_amp = FloatField(null=True)
    oii_3726_amp_ivar = FloatField(null=True)

    # derived target name, e.g., sv3-bright-80613-TARGETID
    target_name = CharField(max_length=40, default='')

    # radec2xyz, for cone search in the database
    ux = FloatField(default=-2.0)
    uy = FloatField(default=-2.0)
    uz = FloatField(default=-2.0)

    #def str_tileid(self):
    #    return str(self.tileid)

    #def nice_target_name(self):
    #    return self.target_name

    def str_hpxpixel(self):
        return str(self.hpxpixel)

    def str_targetid(self):
        return '{}'.format(self.targetid)

    def str_ra(self):
        return '{:.7f}'.format(self.ra)

    def str_dec(self):
        return '{:.7f}'.format(self.dec)

    def base_html_dir(self):
        return '/global/cfs/cdirs/desi/spectro/fastspecfit/test/{}/html/'.format(self.specprod)
        #return '/global/cfs/cdirs/desi/spectro/fastspecfit/{}/html/'.format(specprod)

    def png_base_url(self):
        # different for cumulative coadds!
        baseurl = 'https://data.desi.lbl.gov/desi/spectro/fastspecfit/test/{}/html/healpix/{}/{}/'.format(self.specprod, self.survey, self.faprgrm)
        baseurl += str(self.hpxpixel//100) +'/'+ self.str_hpxpixel()
        return baseurl

    def data_base_url(self):
        # different for cumulative coadds!
        baseurl = 'https://data.desi.lbl.gov/desi/spectro/fastspecfit/test/{}/healpix/{}/{}/'.format(self.specprod, self.survey, self.faprgrm) # no html subdir
        baseurl += str(self.hpxpixel//100) +'/'+ self.str_hpxpixel()
        return baseurl

    #def ellipsefile(self):
    #    ellipsefile = '{}{}-largegalaxy-{}-ellipse-sbprofile.png'.format(self.png_base_url(), self.group_name, self.sga_id_string())
    #    return ellipsefile
    #
    #def ellipse_exists(self):
    #    ellipsefile = os.path.join(self.base_html_dir(), self.ra_slice(), self.group_name, '{}-largegalaxy-{}-ellipse-sbprofile.png'.format(
    #        self.group_name, self.sga_id_string()))
    #    return os.path.isfile(ellipsefile)
