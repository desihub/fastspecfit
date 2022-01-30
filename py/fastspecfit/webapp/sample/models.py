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

    # derived target name, e.g., sv3-bright-80613-TARGETID
    specprod = CharField(max_length=15, default='')
    target_name = CharField(max_length=40, default='')
    DATADIR = CharField(max_length=150, default='')

    # for col in tt.colnames:
    #     if tt[col].dtype.name == 'float64':
    #         print('{} = FloatField(null=True)'.format(col.lower()))
    #     elif 'int' in tt[col].dtype.name:
    #         print('{} = IntegerField(null=True)'.format(col.lower()))

    # metadata columns
    targetid = IntegerField(null=True)
    ra = FloatField(null=True)
    dec = FloatField(null=True)
    survey = CharField(max_length=4, default='')
    faprgrm = CharField(max_length=6, default='')
    hpxpixel = IntegerField(null=True)
    #tileid = IntegerField(null=True)
    #fiber = IntegerField(null=True)
    desi_target = IntegerField(null=True)
    bgs_target = IntegerField(null=True)
    mws_target = IntegerField(null=True)
    sv1_desi_target = IntegerField(null=True)
    sv1_bgs_target = IntegerField(null=True)
    sv1_mws_target = IntegerField(null=True)
    sv2_desi_target = IntegerField(null=True)
    sv2_bgs_target = IntegerField(null=True)
    sv2_mws_target = IntegerField(null=True)
    sv3_desi_target = IntegerField(null=True)
    sv3_bgs_target = IntegerField(null=True)
    sv3_mws_target = IntegerField(null=True)
    scnd_target = IntegerField(null=True)
    sv1_scnd_target = IntegerField(null=True)
    sv2_scnd_target = IntegerField(null=True)
    sv3_scnd_target = IntegerField(null=True)

    z = FloatField(null=True)
    zwarn = IntegerField(null=True)
    deltachi2 = FloatField(null=True)
    spectype = CharField(max_length=6, default='')
    photsys = CharField(max_length=1, default='')

    mw_transmission_g = FloatField(null=True)
    mw_transmission_r = FloatField(null=True)
    mw_transmission_z = FloatField(null=True)
    mw_transmission_w1 = FloatField(null=True)
    mw_transmission_w2 = FloatField(null=True)
    fiberflux_g = FloatField(null=True)
    fiberflux_r = FloatField(null=True)
    fiberflux_z = FloatField(null=True)
    fibertotflux_g = FloatField(null=True)
    fibertotflux_r = FloatField(null=True)
    fibertotflux_z = FloatField(null=True)
    flux_g = FloatField(null=True)
    flux_r = FloatField(null=True)
    flux_z = FloatField(null=True)
    flux_w1 = FloatField(null=True)
    flux_w2 = FloatField(null=True)
    flux_ivar_g = FloatField(null=True)
    flux_ivar_r = FloatField(null=True)
    flux_ivar_z = FloatField(null=True)
    flux_ivar_w1 = FloatField(null=True)
    flux_ivar_w2 = FloatField(null=True)

    abmag_g = CharField(max_length=17, default='')
    abmag_r = CharField(max_length=17, default='')
    abmag_z = CharField(max_length=17, default='')
    abmag_w1 = CharField(max_length=17, default='')
    abmag_w2 = CharField(max_length=17, default='')

    fiberabmag_g = CharField(max_length=6, default='')
    fiberabmag_r = CharField(max_length=6, default='')
    fiberabmag_z = CharField(max_length=6, default='')

    fibertotabmag_g = CharField(max_length=6, default='')
    fibertotabmag_r = CharField(max_length=6, default='')
    fibertotabmag_z = CharField(max_length=6, default='')

    abmag_synth_g = CharField(max_length=6, default='')
    abmag_synth_r = CharField(max_length=6, default='')
    abmag_synth_z = CharField(max_length=6, default='')
    
    abmag_synth_model_g = CharField(max_length=6, default='')
    abmag_synth_model_r = CharField(max_length=6, default='')
    abmag_synth_model_z = CharField(max_length=6, default='')
    
    # continuum properties
    continuum_z = FloatField(null=True)
    #continuum_coeff = FloatField(null=True)
    continuum_chi2 = FloatField(null=True)
    continuum_age = FloatField(null=True)
    continuum_av = FloatField(null=True)
    continuum_av_ivar = FloatField(null=True)
    continuum_vdisp = FloatField(null=True)
    continuum_vdisp_ivar = FloatField(null=True)
    continuum_snr_b = FloatField(null=True)
    continuum_snr_r = FloatField(null=True)
    continuum_snr_z = FloatField(null=True)
    continuum_smoothcorr_b = FloatField(null=True)
    continuum_smoothcorr_r = FloatField(null=True)
    continuum_smoothcorr_z = FloatField(null=True)
    dn4000 = FloatField(null=True)
    dn4000_ivar = FloatField(null=True)
    dn4000_model = FloatField(null=True)
    flux_synth_g = FloatField(null=True)
    flux_synth_r = FloatField(null=True)
    flux_synth_z = FloatField(null=True)
    flux_synth_model_g = FloatField(null=True)
    flux_synth_model_r = FloatField(null=True)
    flux_synth_model_z = FloatField(null=True)
    balmer_z = FloatField(null=True)
    forbidden_z = FloatField(null=True)
    broad_z = FloatField(null=True)
    balmer_sigma = FloatField(null=True)
    forbidden_sigma = FloatField(null=True)
    broad_sigma = FloatField(null=True)
    mgii_doublet_ratio = FloatField(null=True)
    oii_doublet_ratio = FloatField(null=True)
    sii_doublet_ratio = FloatField(null=True)
    oi_1304_amp = FloatField(null=True)
    oi_1304_amp_ivar = FloatField(null=True)
    oi_1304_flux = FloatField(null=True)
    oi_1304_flux_ivar = FloatField(null=True)
    oi_1304_boxflux = FloatField(null=True)
    oi_1304_vshift = FloatField(null=True)
    oi_1304_sigma = FloatField(null=True)
    oi_1304_cont = FloatField(null=True)
    oi_1304_cont_ivar = FloatField(null=True)
    oi_1304_ew = FloatField(null=True)
    oi_1304_ew_ivar = FloatField(null=True)
    oi_1304_flux_limit = FloatField(null=True)
    oi_1304_ew_limit = FloatField(null=True)
    oi_1304_chi2 = FloatField(null=True)
    oi_1304_npix = IntegerField(null=True)
    siliv_1396_amp = FloatField(null=True)
    siliv_1396_amp_ivar = FloatField(null=True)
    siliv_1396_flux = FloatField(null=True)
    siliv_1396_flux_ivar = FloatField(null=True)
    siliv_1396_boxflux = FloatField(null=True)
    siliv_1396_vshift = FloatField(null=True)
    siliv_1396_sigma = FloatField(null=True)
    siliv_1396_cont = FloatField(null=True)
    siliv_1396_cont_ivar = FloatField(null=True)
    siliv_1396_ew = FloatField(null=True)
    siliv_1396_ew_ivar = FloatField(null=True)
    siliv_1396_flux_limit = FloatField(null=True)
    siliv_1396_ew_limit = FloatField(null=True)
    siliv_1396_chi2 = FloatField(null=True)
    siliv_1396_npix = IntegerField(null=True)
    civ_1549_amp = FloatField(null=True)
    civ_1549_amp_ivar = FloatField(null=True)
    civ_1549_flux = FloatField(null=True)
    civ_1549_flux_ivar = FloatField(null=True)
    civ_1549_boxflux = FloatField(null=True)
    civ_1549_vshift = FloatField(null=True)
    civ_1549_sigma = FloatField(null=True)
    civ_1549_cont = FloatField(null=True)
    civ_1549_cont_ivar = FloatField(null=True)
    civ_1549_ew = FloatField(null=True)
    civ_1549_ew_ivar = FloatField(null=True)
    civ_1549_flux_limit = FloatField(null=True)
    civ_1549_ew_limit = FloatField(null=True)
    civ_1549_chi2 = FloatField(null=True)
    civ_1549_npix = IntegerField(null=True)
    siliii_1892_amp = FloatField(null=True)
    siliii_1892_amp_ivar = FloatField(null=True)
    siliii_1892_flux = FloatField(null=True)
    siliii_1892_flux_ivar = FloatField(null=True)
    siliii_1892_boxflux = FloatField(null=True)
    siliii_1892_vshift = FloatField(null=True)
    siliii_1892_sigma = FloatField(null=True)
    siliii_1892_cont = FloatField(null=True)
    siliii_1892_cont_ivar = FloatField(null=True)
    siliii_1892_ew = FloatField(null=True)
    siliii_1892_ew_ivar = FloatField(null=True)
    siliii_1892_flux_limit = FloatField(null=True)
    siliii_1892_ew_limit = FloatField(null=True)
    siliii_1892_chi2 = FloatField(null=True)
    siliii_1892_npix = IntegerField(null=True)
    ciii_1908_amp = FloatField(null=True)
    ciii_1908_amp_ivar = FloatField(null=True)
    ciii_1908_flux = FloatField(null=True)
    ciii_1908_flux_ivar = FloatField(null=True)
    ciii_1908_boxflux = FloatField(null=True)
    ciii_1908_vshift = FloatField(null=True)
    ciii_1908_sigma = FloatField(null=True)
    ciii_1908_cont = FloatField(null=True)
    ciii_1908_cont_ivar = FloatField(null=True)
    ciii_1908_ew = FloatField(null=True)
    ciii_1908_ew_ivar = FloatField(null=True)
    ciii_1908_flux_limit = FloatField(null=True)
    ciii_1908_ew_limit = FloatField(null=True)
    ciii_1908_chi2 = FloatField(null=True)
    ciii_1908_npix = IntegerField(null=True)
    mgii_2796_amp = FloatField(null=True)
    mgii_2796_amp_ivar = FloatField(null=True)
    mgii_2796_flux = FloatField(null=True)
    mgii_2796_flux_ivar = FloatField(null=True)
    mgii_2796_boxflux = FloatField(null=True)
    mgii_2796_vshift = FloatField(null=True)
    mgii_2796_sigma = FloatField(null=True)
    mgii_2796_cont = FloatField(null=True)
    mgii_2796_cont_ivar = FloatField(null=True)
    mgii_2796_ew = FloatField(null=True)
    mgii_2796_ew_ivar = FloatField(null=True)
    mgii_2796_flux_limit = FloatField(null=True)
    mgii_2796_ew_limit = FloatField(null=True)
    mgii_2796_chi2 = FloatField(null=True)
    mgii_2796_npix = IntegerField(null=True)
    mgii_2803_amp = FloatField(null=True)
    mgii_2803_amp_ivar = FloatField(null=True)
    mgii_2803_flux = FloatField(null=True)
    mgii_2803_flux_ivar = FloatField(null=True)
    mgii_2803_boxflux = FloatField(null=True)
    mgii_2803_vshift = FloatField(null=True)
    mgii_2803_sigma = FloatField(null=True)
    mgii_2803_cont = FloatField(null=True)
    mgii_2803_cont_ivar = FloatField(null=True)
    mgii_2803_ew = FloatField(null=True)
    mgii_2803_ew_ivar = FloatField(null=True)
    mgii_2803_flux_limit = FloatField(null=True)
    mgii_2803_ew_limit = FloatField(null=True)
    mgii_2803_chi2 = FloatField(null=True)
    mgii_2803_npix = IntegerField(null=True)
    nev_3346_amp = FloatField(null=True)
    nev_3346_amp_ivar = FloatField(null=True)
    nev_3346_flux = FloatField(null=True)
    nev_3346_flux_ivar = FloatField(null=True)
    nev_3346_boxflux = FloatField(null=True)
    nev_3346_vshift = FloatField(null=True)
    nev_3346_sigma = FloatField(null=True)
    nev_3346_cont = FloatField(null=True)
    nev_3346_cont_ivar = FloatField(null=True)
    nev_3346_ew = FloatField(null=True)
    nev_3346_ew_ivar = FloatField(null=True)
    nev_3346_flux_limit = FloatField(null=True)
    nev_3346_ew_limit = FloatField(null=True)
    nev_3346_chi2 = FloatField(null=True)
    nev_3346_npix = IntegerField(null=True)
    nev_3426_amp = FloatField(null=True)
    nev_3426_amp_ivar = FloatField(null=True)
    nev_3426_flux = FloatField(null=True)
    nev_3426_flux_ivar = FloatField(null=True)
    nev_3426_boxflux = FloatField(null=True)
    nev_3426_vshift = FloatField(null=True)
    nev_3426_sigma = FloatField(null=True)
    nev_3426_cont = FloatField(null=True)
    nev_3426_cont_ivar = FloatField(null=True)
    nev_3426_ew = FloatField(null=True)
    nev_3426_ew_ivar = FloatField(null=True)
    nev_3426_flux_limit = FloatField(null=True)
    nev_3426_ew_limit = FloatField(null=True)
    nev_3426_chi2 = FloatField(null=True)
    nev_3426_npix = IntegerField(null=True)
    oii_3726_amp = FloatField(null=True)
    oii_3726_amp_ivar = FloatField(null=True)
    oii_3726_flux = FloatField(null=True)
    oii_3726_flux_ivar = FloatField(null=True)
    oii_3726_boxflux = FloatField(null=True)
    oii_3726_vshift = FloatField(null=True)
    oii_3726_sigma = FloatField(null=True)
    oii_3726_cont = FloatField(null=True)
    oii_3726_cont_ivar = FloatField(null=True)
    oii_3726_ew = FloatField(null=True)
    oii_3726_ew_ivar = FloatField(null=True)
    oii_3726_flux_limit = FloatField(null=True)
    oii_3726_ew_limit = FloatField(null=True)
    oii_3726_chi2 = FloatField(null=True)
    oii_3726_npix = IntegerField(null=True)
    oii_3729_amp = FloatField(null=True)
    oii_3729_amp_ivar = FloatField(null=True)
    oii_3729_flux = FloatField(null=True)
    oii_3729_flux_ivar = FloatField(null=True)
    oii_3729_boxflux = FloatField(null=True)
    oii_3729_vshift = FloatField(null=True)
    oii_3729_sigma = FloatField(null=True)
    oii_3729_cont = FloatField(null=True)
    oii_3729_cont_ivar = FloatField(null=True)
    oii_3729_ew = FloatField(null=True)
    oii_3729_ew_ivar = FloatField(null=True)
    oii_3729_flux_limit = FloatField(null=True)
    oii_3729_ew_limit = FloatField(null=True)
    oii_3729_chi2 = FloatField(null=True)
    oii_3729_npix = IntegerField(null=True)
    neiii_3869_amp = FloatField(null=True)
    neiii_3869_amp_ivar = FloatField(null=True)
    neiii_3869_flux = FloatField(null=True)
    neiii_3869_flux_ivar = FloatField(null=True)
    neiii_3869_boxflux = FloatField(null=True)
    neiii_3869_vshift = FloatField(null=True)
    neiii_3869_sigma = FloatField(null=True)
    neiii_3869_cont = FloatField(null=True)
    neiii_3869_cont_ivar = FloatField(null=True)
    neiii_3869_ew = FloatField(null=True)
    neiii_3869_ew_ivar = FloatField(null=True)
    neiii_3869_flux_limit = FloatField(null=True)
    neiii_3869_ew_limit = FloatField(null=True)
    neiii_3869_chi2 = FloatField(null=True)
    neiii_3869_npix = IntegerField(null=True)
    hei_3889_amp = FloatField(null=True)
    hei_3889_amp_ivar = FloatField(null=True)
    hei_3889_flux = FloatField(null=True)
    hei_3889_flux_ivar = FloatField(null=True)
    hei_3889_boxflux = FloatField(null=True)
    hei_3889_vshift = FloatField(null=True)
    hei_3889_sigma = FloatField(null=True)
    hei_3889_cont = FloatField(null=True)
    hei_3889_cont_ivar = FloatField(null=True)
    hei_3889_ew = FloatField(null=True)
    hei_3889_ew_ivar = FloatField(null=True)
    hei_3889_flux_limit = FloatField(null=True)
    hei_3889_ew_limit = FloatField(null=True)
    hei_3889_chi2 = FloatField(null=True)
    hei_3889_npix = IntegerField(null=True)
    hei_broad_3889_amp = FloatField(null=True)
    hei_broad_3889_amp_ivar = FloatField(null=True)
    hei_broad_3889_flux = FloatField(null=True)
    hei_broad_3889_flux_ivar = FloatField(null=True)
    hei_broad_3889_boxflux = FloatField(null=True)
    hei_broad_3889_vshift = FloatField(null=True)
    hei_broad_3889_sigma = FloatField(null=True)
    hei_broad_3889_cont = FloatField(null=True)
    hei_broad_3889_cont_ivar = FloatField(null=True)
    hei_broad_3889_ew = FloatField(null=True)
    hei_broad_3889_ew_ivar = FloatField(null=True)
    hei_broad_3889_flux_limit = FloatField(null=True)
    hei_broad_3889_ew_limit = FloatField(null=True)
    hei_broad_3889_chi2 = FloatField(null=True)
    hei_broad_3889_npix = IntegerField(null=True)
    h6_amp = FloatField(null=True)
    h6_amp_ivar = FloatField(null=True)
    h6_flux = FloatField(null=True)
    h6_flux_ivar = FloatField(null=True)
    h6_boxflux = FloatField(null=True)
    h6_vshift = FloatField(null=True)
    h6_sigma = FloatField(null=True)
    h6_cont = FloatField(null=True)
    h6_cont_ivar = FloatField(null=True)
    h6_ew = FloatField(null=True)
    h6_ew_ivar = FloatField(null=True)
    h6_flux_limit = FloatField(null=True)
    h6_ew_limit = FloatField(null=True)
    h6_chi2 = FloatField(null=True)
    h6_npix = IntegerField(null=True)
    h6_broad_amp = FloatField(null=True)
    h6_broad_amp_ivar = FloatField(null=True)
    h6_broad_flux = FloatField(null=True)
    h6_broad_flux_ivar = FloatField(null=True)
    h6_broad_boxflux = FloatField(null=True)
    h6_broad_vshift = FloatField(null=True)
    h6_broad_sigma = FloatField(null=True)
    h6_broad_cont = FloatField(null=True)
    h6_broad_cont_ivar = FloatField(null=True)
    h6_broad_ew = FloatField(null=True)
    h6_broad_ew_ivar = FloatField(null=True)
    h6_broad_flux_limit = FloatField(null=True)
    h6_broad_ew_limit = FloatField(null=True)
    h6_broad_chi2 = FloatField(null=True)
    h6_broad_npix = IntegerField(null=True)
    hepsilon_amp = FloatField(null=True)
    hepsilon_amp_ivar = FloatField(null=True)
    hepsilon_flux = FloatField(null=True)
    hepsilon_flux_ivar = FloatField(null=True)
    hepsilon_boxflux = FloatField(null=True)
    hepsilon_vshift = FloatField(null=True)
    hepsilon_sigma = FloatField(null=True)
    hepsilon_cont = FloatField(null=True)
    hepsilon_cont_ivar = FloatField(null=True)
    hepsilon_ew = FloatField(null=True)
    hepsilon_ew_ivar = FloatField(null=True)
    hepsilon_flux_limit = FloatField(null=True)
    hepsilon_ew_limit = FloatField(null=True)
    hepsilon_chi2 = FloatField(null=True)
    hepsilon_npix = IntegerField(null=True)
    hepsilon_broad_amp = FloatField(null=True)
    hepsilon_broad_amp_ivar = FloatField(null=True)
    hepsilon_broad_flux = FloatField(null=True)
    hepsilon_broad_flux_ivar = FloatField(null=True)
    hepsilon_broad_boxflux = FloatField(null=True)
    hepsilon_broad_vshift = FloatField(null=True)
    hepsilon_broad_sigma = FloatField(null=True)
    hepsilon_broad_cont = FloatField(null=True)
    hepsilon_broad_cont_ivar = FloatField(null=True)
    hepsilon_broad_ew = FloatField(null=True)
    hepsilon_broad_ew_ivar = FloatField(null=True)
    hepsilon_broad_flux_limit = FloatField(null=True)
    hepsilon_broad_ew_limit = FloatField(null=True)
    hepsilon_broad_chi2 = FloatField(null=True)
    hepsilon_broad_npix = IntegerField(null=True)
    hdelta_amp = FloatField(null=True)
    hdelta_amp_ivar = FloatField(null=True)
    hdelta_flux = FloatField(null=True)
    hdelta_flux_ivar = FloatField(null=True)
    hdelta_boxflux = FloatField(null=True)
    hdelta_vshift = FloatField(null=True)
    hdelta_sigma = FloatField(null=True)
    hdelta_cont = FloatField(null=True)
    hdelta_cont_ivar = FloatField(null=True)
    hdelta_ew = FloatField(null=True)
    hdelta_ew_ivar = FloatField(null=True)
    hdelta_flux_limit = FloatField(null=True)
    hdelta_ew_limit = FloatField(null=True)
    hdelta_chi2 = FloatField(null=True)
    hdelta_npix = IntegerField(null=True)
    hdelta_broad_amp = FloatField(null=True)
    hdelta_broad_amp_ivar = FloatField(null=True)
    hdelta_broad_flux = FloatField(null=True)
    hdelta_broad_flux_ivar = FloatField(null=True)
    hdelta_broad_boxflux = FloatField(null=True)
    hdelta_broad_vshift = FloatField(null=True)
    hdelta_broad_sigma = FloatField(null=True)
    hdelta_broad_cont = FloatField(null=True)
    hdelta_broad_cont_ivar = FloatField(null=True)
    hdelta_broad_ew = FloatField(null=True)
    hdelta_broad_ew_ivar = FloatField(null=True)
    hdelta_broad_flux_limit = FloatField(null=True)
    hdelta_broad_ew_limit = FloatField(null=True)
    hdelta_broad_chi2 = FloatField(null=True)
    hdelta_broad_npix = IntegerField(null=True)
    hgamma_amp = FloatField(null=True)
    hgamma_amp_ivar = FloatField(null=True)
    hgamma_flux = FloatField(null=True)
    hgamma_flux_ivar = FloatField(null=True)
    hgamma_boxflux = FloatField(null=True)
    hgamma_vshift = FloatField(null=True)
    hgamma_sigma = FloatField(null=True)
    hgamma_cont = FloatField(null=True)
    hgamma_cont_ivar = FloatField(null=True)
    hgamma_ew = FloatField(null=True)
    hgamma_ew_ivar = FloatField(null=True)
    hgamma_flux_limit = FloatField(null=True)
    hgamma_ew_limit = FloatField(null=True)
    hgamma_chi2 = FloatField(null=True)
    hgamma_npix = IntegerField(null=True)
    hgamma_broad_amp = FloatField(null=True)
    hgamma_broad_amp_ivar = FloatField(null=True)
    hgamma_broad_flux = FloatField(null=True)
    hgamma_broad_flux_ivar = FloatField(null=True)
    hgamma_broad_boxflux = FloatField(null=True)
    hgamma_broad_vshift = FloatField(null=True)
    hgamma_broad_sigma = FloatField(null=True)
    hgamma_broad_cont = FloatField(null=True)
    hgamma_broad_cont_ivar = FloatField(null=True)
    hgamma_broad_ew = FloatField(null=True)
    hgamma_broad_ew_ivar = FloatField(null=True)
    hgamma_broad_flux_limit = FloatField(null=True)
    hgamma_broad_ew_limit = FloatField(null=True)
    hgamma_broad_chi2 = FloatField(null=True)
    hgamma_broad_npix = IntegerField(null=True)
    oiii_4363_amp = FloatField(null=True)
    oiii_4363_amp_ivar = FloatField(null=True)
    oiii_4363_flux = FloatField(null=True)
    oiii_4363_flux_ivar = FloatField(null=True)
    oiii_4363_boxflux = FloatField(null=True)
    oiii_4363_vshift = FloatField(null=True)
    oiii_4363_sigma = FloatField(null=True)
    oiii_4363_cont = FloatField(null=True)
    oiii_4363_cont_ivar = FloatField(null=True)
    oiii_4363_ew = FloatField(null=True)
    oiii_4363_ew_ivar = FloatField(null=True)
    oiii_4363_flux_limit = FloatField(null=True)
    oiii_4363_ew_limit = FloatField(null=True)
    oiii_4363_chi2 = FloatField(null=True)
    oiii_4363_npix = IntegerField(null=True)
    hei_4471_amp = FloatField(null=True)
    hei_4471_amp_ivar = FloatField(null=True)
    hei_4471_flux = FloatField(null=True)
    hei_4471_flux_ivar = FloatField(null=True)
    hei_4471_boxflux = FloatField(null=True)
    hei_4471_vshift = FloatField(null=True)
    hei_4471_sigma = FloatField(null=True)
    hei_4471_cont = FloatField(null=True)
    hei_4471_cont_ivar = FloatField(null=True)
    hei_4471_ew = FloatField(null=True)
    hei_4471_ew_ivar = FloatField(null=True)
    hei_4471_flux_limit = FloatField(null=True)
    hei_4471_ew_limit = FloatField(null=True)
    hei_4471_chi2 = FloatField(null=True)
    hei_4471_npix = IntegerField(null=True)
    hei_broad_4471_amp = FloatField(null=True)
    hei_broad_4471_amp_ivar = FloatField(null=True)
    hei_broad_4471_flux = FloatField(null=True)
    hei_broad_4471_flux_ivar = FloatField(null=True)
    hei_broad_4471_boxflux = FloatField(null=True)
    hei_broad_4471_vshift = FloatField(null=True)
    hei_broad_4471_sigma = FloatField(null=True)
    hei_broad_4471_cont = FloatField(null=True)
    hei_broad_4471_cont_ivar = FloatField(null=True)
    hei_broad_4471_ew = FloatField(null=True)
    hei_broad_4471_ew_ivar = FloatField(null=True)
    hei_broad_4471_flux_limit = FloatField(null=True)
    hei_broad_4471_ew_limit = FloatField(null=True)
    hei_broad_4471_chi2 = FloatField(null=True)
    hei_broad_4471_npix = IntegerField(null=True)
    heii_4686_amp = FloatField(null=True)
    heii_4686_amp_ivar = FloatField(null=True)
    heii_4686_flux = FloatField(null=True)
    heii_4686_flux_ivar = FloatField(null=True)
    heii_4686_boxflux = FloatField(null=True)
    heii_4686_vshift = FloatField(null=True)
    heii_4686_sigma = FloatField(null=True)
    heii_4686_cont = FloatField(null=True)
    heii_4686_cont_ivar = FloatField(null=True)
    heii_4686_ew = FloatField(null=True)
    heii_4686_ew_ivar = FloatField(null=True)
    heii_4686_flux_limit = FloatField(null=True)
    heii_4686_ew_limit = FloatField(null=True)
    heii_4686_chi2 = FloatField(null=True)
    heii_4686_npix = IntegerField(null=True)
    heii_broad_4686_amp = FloatField(null=True)
    heii_broad_4686_amp_ivar = FloatField(null=True)
    heii_broad_4686_flux = FloatField(null=True)
    heii_broad_4686_flux_ivar = FloatField(null=True)
    heii_broad_4686_boxflux = FloatField(null=True)
    heii_broad_4686_vshift = FloatField(null=True)
    heii_broad_4686_sigma = FloatField(null=True)
    heii_broad_4686_cont = FloatField(null=True)
    heii_broad_4686_cont_ivar = FloatField(null=True)
    heii_broad_4686_ew = FloatField(null=True)
    heii_broad_4686_ew_ivar = FloatField(null=True)
    heii_broad_4686_flux_limit = FloatField(null=True)
    heii_broad_4686_ew_limit = FloatField(null=True)
    heii_broad_4686_chi2 = FloatField(null=True)
    heii_broad_4686_npix = IntegerField(null=True)
    hbeta_amp = FloatField(null=True)
    hbeta_amp_ivar = FloatField(null=True)
    hbeta_flux = FloatField(null=True)
    hbeta_flux_ivar = FloatField(null=True)
    hbeta_boxflux = FloatField(null=True)
    hbeta_vshift = FloatField(null=True)
    hbeta_sigma = FloatField(null=True)
    hbeta_cont = FloatField(null=True)
    hbeta_cont_ivar = FloatField(null=True)
    hbeta_ew = FloatField(null=True)
    hbeta_ew_ivar = FloatField(null=True)
    hbeta_flux_limit = FloatField(null=True)
    hbeta_ew_limit = FloatField(null=True)
    hbeta_chi2 = FloatField(null=True)
    hbeta_npix = IntegerField(null=True)
    hbeta_broad_amp = FloatField(null=True)
    hbeta_broad_amp_ivar = FloatField(null=True)
    hbeta_broad_flux = FloatField(null=True)
    hbeta_broad_flux_ivar = FloatField(null=True)
    hbeta_broad_boxflux = FloatField(null=True)
    hbeta_broad_vshift = FloatField(null=True)
    hbeta_broad_sigma = FloatField(null=True)
    hbeta_broad_cont = FloatField(null=True)
    hbeta_broad_cont_ivar = FloatField(null=True)
    hbeta_broad_ew = FloatField(null=True)
    hbeta_broad_ew_ivar = FloatField(null=True)
    hbeta_broad_flux_limit = FloatField(null=True)
    hbeta_broad_ew_limit = FloatField(null=True)
    hbeta_broad_chi2 = FloatField(null=True)
    hbeta_broad_npix = IntegerField(null=True)
    oiii_4959_amp = FloatField(null=True)
    oiii_4959_amp_ivar = FloatField(null=True)
    oiii_4959_flux = FloatField(null=True)
    oiii_4959_flux_ivar = FloatField(null=True)
    oiii_4959_boxflux = FloatField(null=True)
    oiii_4959_vshift = FloatField(null=True)
    oiii_4959_sigma = FloatField(null=True)
    oiii_4959_cont = FloatField(null=True)
    oiii_4959_cont_ivar = FloatField(null=True)
    oiii_4959_ew = FloatField(null=True)
    oiii_4959_ew_ivar = FloatField(null=True)
    oiii_4959_flux_limit = FloatField(null=True)
    oiii_4959_ew_limit = FloatField(null=True)
    oiii_4959_chi2 = FloatField(null=True)
    oiii_4959_npix = IntegerField(null=True)
    oiii_5007_amp = FloatField(null=True)
    oiii_5007_amp_ivar = FloatField(null=True)
    oiii_5007_flux = FloatField(null=True)
    oiii_5007_flux_ivar = FloatField(null=True)
    oiii_5007_boxflux = FloatField(null=True)
    oiii_5007_vshift = FloatField(null=True)
    oiii_5007_sigma = FloatField(null=True)
    oiii_5007_cont = FloatField(null=True)
    oiii_5007_cont_ivar = FloatField(null=True)
    oiii_5007_ew = FloatField(null=True)
    oiii_5007_ew_ivar = FloatField(null=True)
    oiii_5007_flux_limit = FloatField(null=True)
    oiii_5007_ew_limit = FloatField(null=True)
    oiii_5007_chi2 = FloatField(null=True)
    oiii_5007_npix = IntegerField(null=True)
    nii_5755_amp = FloatField(null=True)
    nii_5755_amp_ivar = FloatField(null=True)
    nii_5755_flux = FloatField(null=True)
    nii_5755_flux_ivar = FloatField(null=True)
    nii_5755_boxflux = FloatField(null=True)
    nii_5755_vshift = FloatField(null=True)
    nii_5755_sigma = FloatField(null=True)
    nii_5755_cont = FloatField(null=True)
    nii_5755_cont_ivar = FloatField(null=True)
    nii_5755_ew = FloatField(null=True)
    nii_5755_ew_ivar = FloatField(null=True)
    nii_5755_flux_limit = FloatField(null=True)
    nii_5755_ew_limit = FloatField(null=True)
    nii_5755_chi2 = FloatField(null=True)
    nii_5755_npix = IntegerField(null=True)
    hei_5876_amp = FloatField(null=True)
    hei_5876_amp_ivar = FloatField(null=True)
    hei_5876_flux = FloatField(null=True)
    hei_5876_flux_ivar = FloatField(null=True)
    hei_5876_boxflux = FloatField(null=True)
    hei_5876_vshift = FloatField(null=True)
    hei_5876_sigma = FloatField(null=True)
    hei_5876_cont = FloatField(null=True)
    hei_5876_cont_ivar = FloatField(null=True)
    hei_5876_ew = FloatField(null=True)
    hei_5876_ew_ivar = FloatField(null=True)
    hei_5876_flux_limit = FloatField(null=True)
    hei_5876_ew_limit = FloatField(null=True)
    hei_5876_chi2 = FloatField(null=True)
    hei_5876_npix = IntegerField(null=True)
    hei_broad_5876_amp = FloatField(null=True)
    hei_broad_5876_amp_ivar = FloatField(null=True)
    hei_broad_5876_flux = FloatField(null=True)
    hei_broad_5876_flux_ivar = FloatField(null=True)
    hei_broad_5876_boxflux = FloatField(null=True)
    hei_broad_5876_vshift = FloatField(null=True)
    hei_broad_5876_sigma = FloatField(null=True)
    hei_broad_5876_cont = FloatField(null=True)
    hei_broad_5876_cont_ivar = FloatField(null=True)
    hei_broad_5876_ew = FloatField(null=True)
    hei_broad_5876_ew_ivar = FloatField(null=True)
    hei_broad_5876_flux_limit = FloatField(null=True)
    hei_broad_5876_ew_limit = FloatField(null=True)
    hei_broad_5876_chi2 = FloatField(null=True)
    hei_broad_5876_npix = IntegerField(null=True)
    oi_6300_amp = FloatField(null=True)
    oi_6300_amp_ivar = FloatField(null=True)
    oi_6300_flux = FloatField(null=True)
    oi_6300_flux_ivar = FloatField(null=True)
    oi_6300_boxflux = FloatField(null=True)
    oi_6300_vshift = FloatField(null=True)
    oi_6300_sigma = FloatField(null=True)
    oi_6300_cont = FloatField(null=True)
    oi_6300_cont_ivar = FloatField(null=True)
    oi_6300_ew = FloatField(null=True)
    oi_6300_ew_ivar = FloatField(null=True)
    oi_6300_flux_limit = FloatField(null=True)
    oi_6300_ew_limit = FloatField(null=True)
    oi_6300_chi2 = FloatField(null=True)
    oi_6300_npix = IntegerField(null=True)
    nii_6548_amp = FloatField(null=True)
    nii_6548_amp_ivar = FloatField(null=True)
    nii_6548_flux = FloatField(null=True)
    nii_6548_flux_ivar = FloatField(null=True)
    nii_6548_boxflux = FloatField(null=True)
    nii_6548_vshift = FloatField(null=True)
    nii_6548_sigma = FloatField(null=True)
    nii_6548_cont = FloatField(null=True)
    nii_6548_cont_ivar = FloatField(null=True)
    nii_6548_ew = FloatField(null=True)
    nii_6548_ew_ivar = FloatField(null=True)
    nii_6548_flux_limit = FloatField(null=True)
    nii_6548_ew_limit = FloatField(null=True)
    nii_6548_chi2 = FloatField(null=True)
    nii_6548_npix = IntegerField(null=True)
    halpha_amp = FloatField(null=True)
    halpha_amp_ivar = FloatField(null=True)
    halpha_flux = FloatField(null=True)
    halpha_flux_ivar = FloatField(null=True)
    halpha_boxflux = FloatField(null=True)
    halpha_vshift = FloatField(null=True)
    halpha_sigma = FloatField(null=True)
    halpha_cont = FloatField(null=True)
    halpha_cont_ivar = FloatField(null=True)
    halpha_ew = FloatField(null=True)
    halpha_ew_ivar = FloatField(null=True)
    halpha_flux_limit = FloatField(null=True)
    halpha_ew_limit = FloatField(null=True)
    halpha_chi2 = FloatField(null=True)
    halpha_npix = IntegerField(null=True)
    halpha_broad_amp = FloatField(null=True)
    halpha_broad_amp_ivar = FloatField(null=True)
    halpha_broad_flux = FloatField(null=True)
    halpha_broad_flux_ivar = FloatField(null=True)
    halpha_broad_boxflux = FloatField(null=True)
    halpha_broad_vshift = FloatField(null=True)
    halpha_broad_sigma = FloatField(null=True)
    halpha_broad_cont = FloatField(null=True)
    halpha_broad_cont_ivar = FloatField(null=True)
    halpha_broad_ew = FloatField(null=True)
    halpha_broad_ew_ivar = FloatField(null=True)
    halpha_broad_flux_limit = FloatField(null=True)
    halpha_broad_ew_limit = FloatField(null=True)
    halpha_broad_chi2 = FloatField(null=True)
    halpha_broad_npix = IntegerField(null=True)
    nii_6584_amp = FloatField(null=True)
    nii_6584_amp_ivar = FloatField(null=True)
    nii_6584_flux = FloatField(null=True)
    nii_6584_flux_ivar = FloatField(null=True)
    nii_6584_boxflux = FloatField(null=True)
    nii_6584_vshift = FloatField(null=True)
    nii_6584_sigma = FloatField(null=True)
    nii_6584_cont = FloatField(null=True)
    nii_6584_cont_ivar = FloatField(null=True)
    nii_6584_ew = FloatField(null=True)
    nii_6584_ew_ivar = FloatField(null=True)
    nii_6584_flux_limit = FloatField(null=True)
    nii_6584_ew_limit = FloatField(null=True)
    nii_6584_chi2 = FloatField(null=True)
    nii_6584_npix = IntegerField(null=True)
    sii_6716_amp = FloatField(null=True)
    sii_6716_amp_ivar = FloatField(null=True)
    sii_6716_flux = FloatField(null=True)
    sii_6716_flux_ivar = FloatField(null=True)
    sii_6716_boxflux = FloatField(null=True)
    sii_6716_vshift = FloatField(null=True)
    sii_6716_sigma = FloatField(null=True)
    sii_6716_cont = FloatField(null=True)
    sii_6716_cont_ivar = FloatField(null=True)
    sii_6716_ew = FloatField(null=True)
    sii_6716_ew_ivar = FloatField(null=True)
    sii_6716_flux_limit = FloatField(null=True)
    sii_6716_ew_limit = FloatField(null=True)
    sii_6716_chi2 = FloatField(null=True)
    sii_6716_npix = IntegerField(null=True)
    sii_6731_amp = FloatField(null=True)
    sii_6731_amp_ivar = FloatField(null=True)
    sii_6731_flux = FloatField(null=True)
    sii_6731_flux_ivar = FloatField(null=True)
    sii_6731_boxflux = FloatField(null=True)
    sii_6731_vshift = FloatField(null=True)
    sii_6731_sigma = FloatField(null=True)
    sii_6731_cont = FloatField(null=True)
    sii_6731_cont_ivar = FloatField(null=True)
    sii_6731_ew = FloatField(null=True)
    sii_6731_ew_ivar = FloatField(null=True)
    sii_6731_flux_limit = FloatField(null=True)
    sii_6731_ew_limit = FloatField(null=True)
    sii_6731_chi2 = FloatField(null=True)
    sii_6731_npix = IntegerField(null=True)
    oii_7320_amp = FloatField(null=True)
    oii_7320_amp_ivar = FloatField(null=True)
    oii_7320_flux = FloatField(null=True)
    oii_7320_flux_ivar = FloatField(null=True)
    oii_7320_boxflux = FloatField(null=True)
    oii_7320_vshift = FloatField(null=True)
    oii_7320_sigma = FloatField(null=True)
    oii_7320_cont = FloatField(null=True)
    oii_7320_cont_ivar = FloatField(null=True)
    oii_7320_ew = FloatField(null=True)
    oii_7320_ew_ivar = FloatField(null=True)
    oii_7320_flux_limit = FloatField(null=True)
    oii_7320_ew_limit = FloatField(null=True)
    oii_7320_chi2 = FloatField(null=True)
    oii_7320_npix = IntegerField(null=True)
    oii_7330_amp = FloatField(null=True)
    oii_7330_amp_ivar = FloatField(null=True)
    oii_7330_flux = FloatField(null=True)
    oii_7330_flux_ivar = FloatField(null=True)
    oii_7330_boxflux = FloatField(null=True)
    oii_7330_vshift = FloatField(null=True)
    oii_7330_sigma = FloatField(null=True)
    oii_7330_cont = FloatField(null=True)
    oii_7330_cont_ivar = FloatField(null=True)
    oii_7330_ew = FloatField(null=True)
    oii_7330_ew_ivar = FloatField(null=True)
    oii_7330_flux_limit = FloatField(null=True)
    oii_7330_ew_limit = FloatField(null=True)
    oii_7330_chi2 = FloatField(null=True)
    oii_7330_npix = IntegerField(null=True)
    siii_9069_amp = FloatField(null=True)
    siii_9069_amp_ivar = FloatField(null=True)
    siii_9069_flux = FloatField(null=True)
    siii_9069_flux_ivar = FloatField(null=True)
    siii_9069_boxflux = FloatField(null=True)
    siii_9069_vshift = FloatField(null=True)
    siii_9069_sigma = FloatField(null=True)
    siii_9069_cont = FloatField(null=True)
    siii_9069_cont_ivar = FloatField(null=True)
    siii_9069_ew = FloatField(null=True)
    siii_9069_ew_ivar = FloatField(null=True)
    siii_9069_flux_limit = FloatField(null=True)
    siii_9069_ew_limit = FloatField(null=True)
    siii_9069_chi2 = FloatField(null=True)
    siii_9069_npix = IntegerField(null=True)
    siii_9532_amp = FloatField(null=True)
    siii_9532_amp_ivar = FloatField(null=True)
    siii_9532_flux = FloatField(null=True)
    siii_9532_flux_ivar = FloatField(null=True)
    siii_9532_boxflux = FloatField(null=True)
    siii_9532_vshift = FloatField(null=True)
    siii_9532_sigma = FloatField(null=True)
    siii_9532_cont = FloatField(null=True)
    siii_9532_cont_ivar = FloatField(null=True)
    siii_9532_ew = FloatField(null=True)
    siii_9532_ew_ivar = FloatField(null=True)
    siii_9532_flux_limit = FloatField(null=True)
    siii_9532_ew_limit = FloatField(null=True)
    siii_9532_chi2 = FloatField(null=True)
    siii_9532_npix = IntegerField(null=True)
            
    # radec2xyz, for cone search in the database
    ux = FloatField(default=-2.0)
    uy = FloatField(default=-2.0)
    uz = FloatField(default=-2.0)

    #def str_abmag_g(self):
    #    return self.abmag_g.replace('+/-', r'&#177;')
    #def str_abmag_r(self):
    #    return self.abmag_r.replace('+/-', r'&#177;')
    #def str_abmag_z(self):
    #    return self.abmag_z.replace('+/-', r'&#177;')
    #def str_abmag_w1(self):
    #    return self.abmag_w1.replace('+/-', r'&#177;')
    #def str_abmag_w2(self):
    #    return self.abmag_w2.replace('+/-', r'&#177;')

    def str_hpxpixel(self):
        return str(self.hpxpixel)

    def str_targetid(self):
        return '{}'.format(self.targetid)

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
