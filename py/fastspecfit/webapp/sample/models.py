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
    # in FITS table
    row_index = IntegerField(default=-1)

    targetid = IntegerField(null=True)
    ra = FloatField(null=True)
    dec = FloatField(null=True)

    # radec2xyz, for cone search in the database
    ux = FloatField(default=-2.0)
    uy = FloatField(default=-2.0)
    uz = FloatField(default=-2.0)

    specprod = 'everest'
    survey = 'sv3'
    program = 'dark'
    tile = '80605'

    def base_html_dir(self):
        return '/global/cfs/cdirs/desi/spectro/fastspecfit/{}/html/'.format(specprod)

    def png_base_url(self):
        baseurl = 'https://data.desi.lbl.gov/desi/spectro/fastspecfit/{}/html/'.format(specprod)
        baseurl += self.ra_slice() + '/' + self.group_name + '/';
        return baseurl

    def data_base_url(self):
        baseurl = 'https://data.desi.lbl.gov/desi/spectro/fastspecfit/{}/'.format(specprod)
        baseurl += self.ra_slice() + '/' + self.group_name + '/';
        return baseurl

    #def mosaic_diam(self):
    #    if self.group_diameter > 30: # NGC0598=M33 is 61 arcmin in diameter!
    #        mosaic_diam = self.group_diameter * 2 * 0.7 # [arcmin]
    #    elif self.group_diameter > 14 and self.group_diameter < 30:
    #        mosaic_diam = self.group_diameter * 2 * 1.0 # [arcmin]
    #    else:
    #        mosaic_diam = self.group_diameter * 2 * 1.5 # [arcmin]
    #    return '{:.3f}'.format(mosaic_diam) # [arcmin]

    def targetid_string(self):
        return '{}'.format(self.targetid)

    #def ellipsefile(self):
    #    ellipsefile = '{}{}-largegalaxy-{}-ellipse-sbprofile.png'.format(self.png_base_url(), self.group_name, self.sga_id_string())
    #    return ellipsefile
    #
    #def ellipse_exists(self):
    #    ellipsefile = os.path.join(self.base_html_dir(), self.ra_slice(), self.group_name, '{}-largegalaxy-{}-ellipse-sbprofile.png'.format(
    #        self.group_name, self.sga_id_string()))
    #    return os.path.isfile(ellipsefile)
