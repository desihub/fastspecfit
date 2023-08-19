#!/usr/bin/env python

"""Load the input sample into a database table.

"""
import os, pdb
import numpy as np
import fitsio
import django
from pkg_resources import resource_filename
from astropy.table import Table, hstack

C_LIGHT = 299792.458 # [km/s]

# change me!
specprod = 'fuji'

DATADIR = '/global/cfs/cdirs/desi/spectro/fastspecfit/{}/catalogs'.format(specprod)

fastspecfile = os.path.join(DATADIR, 'fastspec-{}.fits'.format(specprod))
#fastphotfile = os.path.join(DATADIR, 'fastphot-{}.fits'.format(specprod))

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

    from fastspecfit.webapp.fastmodel.models import FastModel

    meta_columns = [
        'TARGETID',
        'SURVEY',
        'PROGRAM',
        'HEALPIX',
        'TILEID_LIST',
        #'TILEID',
        'RA',
        'DEC',
        'COADD_FIBERSTATUS',
        'CMX_TARGET',
        'DESI_TARGET',
        'BGS_TARGET',
        'MWS_TARGET',
        'SV1_DESI_TARGET',
        'SV1_BGS_TARGET',
        'SV1_MWS_TARGET',
        'SV2_DESI_TARGET',
        'SV2_BGS_TARGET',
        'SV2_MWS_TARGET',
        'SV3_DESI_TARGET',
        'SV3_BGS_TARGET',
        'SV3_MWS_TARGET',
        'SCND_TARGET',
        'SV1_SCND_TARGET',
        'SV2_SCND_TARGET',
        'SV3_SCND_TARGET',
        'Z',
        'ZWARN',
        'DELTACHI2',
        'SPECTYPE',
        'Z_RR',
        'TSNR2_BGS',
        'TSNR2_LRG',
        'TSNR2_ELG',
        'TSNR2_QSO',
        'TSNR2_LYA',
        'PHOTSYS',
        'EBV',
        'MW_TRANSMISSION_G',
        'MW_TRANSMISSION_R',
        'MW_TRANSMISSION_Z',
        'MW_TRANSMISSION_W1',
        'MW_TRANSMISSION_W2',
        'MW_TRANSMISSION_W3',
        'MW_TRANSMISSION_W4',
        'FIBERFLUX_G',
        'FIBERFLUX_R',
        'FIBERFLUX_Z',
        'FIBERTOTFLUX_G',
        'FIBERTOTFLUX_R',
        'FIBERTOTFLUX_Z',
        'FLUX_G',
        'FLUX_R',
        'FLUX_Z',
        'FLUX_W1',
        'FLUX_W2',
        'FLUX_W3',
        'FLUX_W4',
        'FLUX_IVAR_G',
        'FLUX_IVAR_R',
        'FLUX_IVAR_Z',
        'FLUX_IVAR_W1',
        'FLUX_IVAR_W2',
        'FLUX_IVAR_W3',
        'FLUX_IVAR_W4',
        ]
        
    fastspec_cols = [
        #'Z',
        #'COEFF',
        'RCHI2',
        'RCHI2_CONT',
        'RCHI2_PHOT',
        'SNR_B',
        'SNR_R',
        'SNR_Z',
        'SMOOTHCORR_B',
        'SMOOTHCORR_R',
        'SMOOTHCORR_Z',
        'VDISP',
        'VDISP_IVAR',
        'AGE',
        'ZZSUN',
        'LOGMSTAR',
        'SFR',
        #'FAGN',
        'AV',
        'DN4000',
        'DN4000_OBS',
        'DN4000_IVAR',
        'DN4000_MODEL',
        'FLUX_SYNTH_G',
        'FLUX_SYNTH_R',
        'FLUX_SYNTH_Z',
        'FLUX_SYNTH_SPECMODEL_G',
        'FLUX_SYNTH_SPECMODEL_R',
        'FLUX_SYNTH_SPECMODEL_Z',
        'FLUX_SYNTH_PHOTMODEL_G',
        'FLUX_SYNTH_PHOTMODEL_R',
        'FLUX_SYNTH_PHOTMODEL_Z',
        'FLUX_SYNTH_PHOTMODEL_W1',
        'FLUX_SYNTH_PHOTMODEL_W2',
        'FLUX_SYNTH_PHOTMODEL_W3',
        'FLUX_SYNTH_PHOTMODEL_W4',
        'ABSMAG10_DECAM_G',
        'ABSMAG10_IVAR_DECAM_G',
        'KCORR10_DECAM_G',
        'ABSMAG10_DECAM_R',
        'ABSMAG10_IVAR_DECAM_R',
        'KCORR10_DECAM_R',
        'ABSMAG10_DECAM_Z',
        'ABSMAG10_IVAR_DECAM_Z',
        'KCORR10_DECAM_Z',
        'ABSMAG00_U',
        'ABSMAG00_IVAR_U',
        'KCORR00_U',
        'ABSMAG00_B',
        'ABSMAG00_IVAR_B',
        'KCORR00_B',
        'ABSMAG00_V',
        'ABSMAG00_IVAR_V',
        'KCORR00_V',
        'ABSMAG01_SDSS_U',
        'ABSMAG01_IVAR_SDSS_U',
        'KCORR01_SDSS_U',
        'ABSMAG01_SDSS_G',
        'ABSMAG01_IVAR_SDSS_G',
        'KCORR01_SDSS_G',
        'ABSMAG01_SDSS_R',
        'ABSMAG01_IVAR_SDSS_R',
        'KCORR01_SDSS_R',
        'ABSMAG01_SDSS_I',
        'ABSMAG01_IVAR_SDSS_I',
        'KCORR01_SDSS_I',
        'ABSMAG01_SDSS_Z',
        'ABSMAG01_IVAR_SDSS_Z',
        'KCORR01_SDSS_Z',
        'ABSMAG01_W1',
        'ABSMAG01_IVAR_W1',
        'KCORR01_W1',
        #'ABSMAG_W2',
        #'ABSMAG_IVAR_W2',
        #'KCORR_W2',
        'LOGLNU_1500',
        'LOGLNU_2800',
        'LOGL_1450',
        'LOGL_1700',
        'LOGL_3000',
        'LOGL_5100',
        'APERCORR',
        'APERCORR_G',
        'APERCORR_R',
        'APERCORR_Z',
        'RCHI2_LINE',
        'DELTA_LINECHI2',
        'DELTA_LINENDOF',
        'NARROW_Z',
        'BROAD_Z',
        'UV_Z',
        'NARROW_SIGMA',
        'BROAD_SIGMA',
        'UV_SIGMA',
        'NARROW_ZRMS',
        'BROAD_ZRMS',
        'UV_ZRMS',
        'NARROW_SIGMARMS',
        'BROAD_SIGMARMS',
        'UV_SIGMARMS',
        'MGII_DOUBLET_RATIO',
        'OII_DOUBLET_RATIO',
        'SII_DOUBLET_RATIO',
        'LYALPHA_AMP',
        'LYALPHA_AMP_IVAR',
        'LYALPHA_FLUX',
        'LYALPHA_FLUX_IVAR',
        'LYALPHA_BOXFLUX',
        'LYALPHA_VSHIFT',
        'LYALPHA_SIGMA',
        'LYALPHA_CONT',
        'LYALPHA_CONT_IVAR',
        'LYALPHA_EW',
        'LYALPHA_EW_IVAR',
        'LYALPHA_FLUX_LIMIT',
        'LYALPHA_EW_LIMIT',
        'LYALPHA_CHI2',
        'LYALPHA_NPIX',
        'OI_1304_AMP',
        'OI_1304_AMP_IVAR',
        'OI_1304_FLUX',
        'OI_1304_FLUX_IVAR',
        'OI_1304_BOXFLUX',
        'OI_1304_VSHIFT',
        'OI_1304_SIGMA',
        'OI_1304_CONT',
        'OI_1304_CONT_IVAR',
        'OI_1304_EW',
        'OI_1304_EW_IVAR',
        'OI_1304_FLUX_LIMIT',
        'OI_1304_EW_LIMIT',
        'OI_1304_CHI2',
        'OI_1304_NPIX',
        'SILIV_1396_AMP',
        'SILIV_1396_AMP_IVAR',
        'SILIV_1396_FLUX',
        'SILIV_1396_FLUX_IVAR',
        'SILIV_1396_BOXFLUX',
        'SILIV_1396_VSHIFT',
        'SILIV_1396_SIGMA',
        'SILIV_1396_CONT',
        'SILIV_1396_CONT_IVAR',
        'SILIV_1396_EW',
        'SILIV_1396_EW_IVAR',
        'SILIV_1396_FLUX_LIMIT',
        'SILIV_1396_EW_LIMIT',
        'SILIV_1396_CHI2',
        'SILIV_1396_NPIX',
        'CIV_1549_AMP',
        'CIV_1549_AMP_IVAR',
        'CIV_1549_FLUX',
        'CIV_1549_FLUX_IVAR',
        'CIV_1549_BOXFLUX',
        'CIV_1549_VSHIFT',
        'CIV_1549_SIGMA',
        'CIV_1549_CONT',
        'CIV_1549_CONT_IVAR',
        'CIV_1549_EW',
        'CIV_1549_EW_IVAR',
        'CIV_1549_FLUX_LIMIT',
        'CIV_1549_EW_LIMIT',
        'CIV_1549_CHI2',
        'CIV_1549_NPIX',
        'HEII_1640_AMP',
        'HEII_1640_AMP_IVAR',
        'HEII_1640_FLUX',
        'HEII_1640_FLUX_IVAR',
        'HEII_1640_BOXFLUX',
        'HEII_1640_VSHIFT',
        'HEII_1640_SIGMA',
        'HEII_1640_CONT',
        'HEII_1640_CONT_IVAR',
        'HEII_1640_EW',
        'HEII_1640_EW_IVAR',
        'HEII_1640_FLUX_LIMIT',
        'HEII_1640_EW_LIMIT',
        'HEII_1640_CHI2',
        'HEII_1640_NPIX',
        'ALIII_1857_AMP',
        'ALIII_1857_AMP_IVAR',
        'ALIII_1857_FLUX',
        'ALIII_1857_FLUX_IVAR',
        'ALIII_1857_BOXFLUX',
        'ALIII_1857_VSHIFT',
        'ALIII_1857_SIGMA',
        'ALIII_1857_CONT',
        'ALIII_1857_CONT_IVAR',
        'ALIII_1857_EW',
        'ALIII_1857_EW_IVAR',
        'ALIII_1857_FLUX_LIMIT',
        'ALIII_1857_EW_LIMIT',
        'ALIII_1857_CHI2',
        'ALIII_1857_NPIX',
        'SILIII_1892_AMP',
        'SILIII_1892_AMP_IVAR',
        'SILIII_1892_FLUX',
        'SILIII_1892_FLUX_IVAR',
        'SILIII_1892_BOXFLUX',
        'SILIII_1892_VSHIFT',
        'SILIII_1892_SIGMA',
        'SILIII_1892_CONT',
        'SILIII_1892_CONT_IVAR',
        'SILIII_1892_EW',
        'SILIII_1892_EW_IVAR',
        'SILIII_1892_FLUX_LIMIT',
        'SILIII_1892_EW_LIMIT',
        'SILIII_1892_CHI2',
        'SILIII_1892_NPIX',
        'CIII_1908_AMP',
        'CIII_1908_AMP_IVAR',
        'CIII_1908_FLUX',
        'CIII_1908_FLUX_IVAR',
        'CIII_1908_BOXFLUX',
        'CIII_1908_VSHIFT',
        'CIII_1908_SIGMA',
        'CIII_1908_CONT',
        'CIII_1908_CONT_IVAR',
        'CIII_1908_EW',
        'CIII_1908_EW_IVAR',
        'CIII_1908_FLUX_LIMIT',
        'CIII_1908_EW_LIMIT',
        'CIII_1908_CHI2',
        'CIII_1908_NPIX',
        'MGII_2796_AMP',
        'MGII_2796_AMP_IVAR',
        'MGII_2796_FLUX',
        'MGII_2796_FLUX_IVAR',
        'MGII_2796_BOXFLUX',
        'MGII_2796_VSHIFT',
        'MGII_2796_SIGMA',
        'MGII_2796_CONT',
        'MGII_2796_CONT_IVAR',
        'MGII_2796_EW',
        'MGII_2796_EW_IVAR',
        'MGII_2796_FLUX_LIMIT',
        'MGII_2796_EW_LIMIT',
        'MGII_2796_CHI2',
        'MGII_2796_NPIX',
        'MGII_2803_AMP',
        'MGII_2803_AMP_IVAR',
        'MGII_2803_FLUX',
        'MGII_2803_FLUX_IVAR',
        'MGII_2803_BOXFLUX',
        'MGII_2803_VSHIFT',
        'MGII_2803_SIGMA',
        'MGII_2803_CONT',
        'MGII_2803_CONT_IVAR',
        'MGII_2803_EW',
        'MGII_2803_EW_IVAR',
        'MGII_2803_FLUX_LIMIT',
        'MGII_2803_EW_LIMIT',
        'MGII_2803_CHI2',
        'MGII_2803_NPIX',
        'NEV_3346_AMP',
        'NEV_3346_AMP_IVAR',
        'NEV_3346_FLUX',
        'NEV_3346_FLUX_IVAR',
        'NEV_3346_BOXFLUX',
        'NEV_3346_VSHIFT',
        'NEV_3346_SIGMA',
        'NEV_3346_CONT',
        'NEV_3346_CONT_IVAR',
        'NEV_3346_EW',
        'NEV_3346_EW_IVAR',
        'NEV_3346_FLUX_LIMIT',
        'NEV_3346_EW_LIMIT',
        'NEV_3346_CHI2',
        'NEV_3346_NPIX',
        'NEV_3426_AMP',
        'NEV_3426_AMP_IVAR',
        'NEV_3426_FLUX',
        'NEV_3426_FLUX_IVAR',
        'NEV_3426_BOXFLUX',
        'NEV_3426_VSHIFT',
        'NEV_3426_SIGMA',
        'NEV_3426_CONT',
        'NEV_3426_CONT_IVAR',
        'NEV_3426_EW',
        'NEV_3426_EW_IVAR',
        'NEV_3426_FLUX_LIMIT',
        'NEV_3426_EW_LIMIT',
        'NEV_3426_CHI2',
        'NEV_3426_NPIX',
        'OII_3726_AMP',
        'OII_3726_AMP_IVAR',
        'OII_3726_FLUX',
        'OII_3726_FLUX_IVAR',
        'OII_3726_BOXFLUX',
        'OII_3726_VSHIFT',
        'OII_3726_SIGMA',
        'OII_3726_CONT',
        'OII_3726_CONT_IVAR',
        'OII_3726_EW',
        'OII_3726_EW_IVAR',
        'OII_3726_FLUX_LIMIT',
        'OII_3726_EW_LIMIT',
        'OII_3726_CHI2',
        'OII_3726_NPIX',
        'OII_3729_AMP',
        'OII_3729_AMP_IVAR',
        'OII_3729_FLUX',
        'OII_3729_FLUX_IVAR',
        'OII_3729_BOXFLUX',
        'OII_3729_VSHIFT',
        'OII_3729_SIGMA',
        'OII_3729_CONT',
        'OII_3729_CONT_IVAR',
        'OII_3729_EW',
        'OII_3729_EW_IVAR',
        'OII_3729_FLUX_LIMIT',
        'OII_3729_EW_LIMIT',
        'OII_3729_CHI2',
        'OII_3729_NPIX',
        'NEIII_3869_AMP',
        'NEIII_3869_AMP_IVAR',
        'NEIII_3869_FLUX',
        'NEIII_3869_FLUX_IVAR',
        'NEIII_3869_BOXFLUX',
        'NEIII_3869_VSHIFT',
        'NEIII_3869_SIGMA',
        'NEIII_3869_CONT',
        'NEIII_3869_CONT_IVAR',
        'NEIII_3869_EW',
        'NEIII_3869_EW_IVAR',
        'NEIII_3869_FLUX_LIMIT',
        'NEIII_3869_EW_LIMIT',
        'NEIII_3869_CHI2',
        'NEIII_3869_NPIX',
        #'HEI_3889_AMP',
        #'HEI_3889_AMP_IVAR',
        #'HEI_3889_FLUX',
        #'HEI_3889_FLUX_IVAR',
        #'HEI_3889_BOXFLUX',
        #'HEI_3889_VSHIFT',
        #'HEI_3889_SIGMA',
        #'HEI_3889_CONT',
        #'HEI_3889_CONT_IVAR',
        #'HEI_3889_EW',
        #'HEI_3889_EW_IVAR',
        #'HEI_3889_FLUX_LIMIT',
        #'HEI_3889_EW_LIMIT',
        #'HEI_3889_CHI2',
        #'HEI_3889_NPIX',
        #'HEI_BROAD_3889_AMP',
        #'HEI_BROAD_3889_AMP_IVAR',
        #'HEI_BROAD_3889_FLUX',
        #'HEI_BROAD_3889_FLUX_IVAR',
        #'HEI_BROAD_3889_BOXFLUX',
        #'HEI_BROAD_3889_VSHIFT',
        #'HEI_BROAD_3889_SIGMA',
        #'HEI_BROAD_3889_CONT',
        #'HEI_BROAD_3889_CONT_IVAR',
        #'HEI_BROAD_3889_EW',
        #'HEI_BROAD_3889_EW_IVAR',
        #'HEI_BROAD_3889_FLUX_LIMIT',
        #'HEI_BROAD_3889_EW_LIMIT',
        #'HEI_BROAD_3889_CHI2',
        #'HEI_BROAD_3889_NPIX',
        'H6_AMP',
        'H6_AMP_IVAR',
        'H6_FLUX',
        'H6_FLUX_IVAR',
        'H6_BOXFLUX',
        'H6_VSHIFT',
        'H6_SIGMA',
        'H6_CONT',
        'H6_CONT_IVAR',
        'H6_EW',
        'H6_EW_IVAR',
        'H6_FLUX_LIMIT',
        'H6_EW_LIMIT',
        'H6_CHI2',
        'H6_NPIX',
        'H6_BROAD_AMP',
        'H6_BROAD_AMP_IVAR',
        'H6_BROAD_FLUX',
        'H6_BROAD_FLUX_IVAR',
        'H6_BROAD_BOXFLUX',
        'H6_BROAD_VSHIFT',
        'H6_BROAD_SIGMA',
        'H6_BROAD_CONT',
        'H6_BROAD_CONT_IVAR',
        'H6_BROAD_EW',
        'H6_BROAD_EW_IVAR',
        'H6_BROAD_FLUX_LIMIT',
        'H6_BROAD_EW_LIMIT',
        'H6_BROAD_CHI2',
        'H6_BROAD_NPIX',
        'HEPSILON_AMP',
        'HEPSILON_AMP_IVAR',
        'HEPSILON_FLUX',
        'HEPSILON_FLUX_IVAR',
        'HEPSILON_BOXFLUX',
        'HEPSILON_VSHIFT',
        'HEPSILON_SIGMA',
        'HEPSILON_CONT',
        'HEPSILON_CONT_IVAR',
        'HEPSILON_EW',
        'HEPSILON_EW_IVAR',
        'HEPSILON_FLUX_LIMIT',
        'HEPSILON_EW_LIMIT',
        'HEPSILON_CHI2',
        'HEPSILON_NPIX',
        'HEPSILON_BROAD_AMP',
        'HEPSILON_BROAD_AMP_IVAR',
        'HEPSILON_BROAD_FLUX',
        'HEPSILON_BROAD_FLUX_IVAR',
        'HEPSILON_BROAD_BOXFLUX',
        'HEPSILON_BROAD_VSHIFT',
        'HEPSILON_BROAD_SIGMA',
        'HEPSILON_BROAD_CONT',
        'HEPSILON_BROAD_CONT_IVAR',
        'HEPSILON_BROAD_EW',
        'HEPSILON_BROAD_EW_IVAR',
        'HEPSILON_BROAD_FLUX_LIMIT',
        'HEPSILON_BROAD_EW_LIMIT',
        'HEPSILON_BROAD_CHI2',
        'HEPSILON_BROAD_NPIX',
        'HDELTA_AMP',
        'HDELTA_AMP_IVAR',
        'HDELTA_FLUX',
        'HDELTA_FLUX_IVAR',
        'HDELTA_BOXFLUX',
        'HDELTA_VSHIFT',
        'HDELTA_SIGMA',
        'HDELTA_CONT',
        'HDELTA_CONT_IVAR',
        'HDELTA_EW',
        'HDELTA_EW_IVAR',
        'HDELTA_FLUX_LIMIT',
        'HDELTA_EW_LIMIT',
        'HDELTA_CHI2',
        'HDELTA_NPIX',
        'HDELTA_BROAD_AMP',
        'HDELTA_BROAD_AMP_IVAR',
        'HDELTA_BROAD_FLUX',
        'HDELTA_BROAD_FLUX_IVAR',
        'HDELTA_BROAD_BOXFLUX',
        'HDELTA_BROAD_VSHIFT',
        'HDELTA_BROAD_SIGMA',
        'HDELTA_BROAD_CONT',
        'HDELTA_BROAD_CONT_IVAR',
        'HDELTA_BROAD_EW',
        'HDELTA_BROAD_EW_IVAR',
        'HDELTA_BROAD_FLUX_LIMIT',
        'HDELTA_BROAD_EW_LIMIT',
        'HDELTA_BROAD_CHI2',
        'HDELTA_BROAD_NPIX',
        'HGAMMA_AMP',
        'HGAMMA_AMP_IVAR',
        'HGAMMA_FLUX',
        'HGAMMA_FLUX_IVAR',
        'HGAMMA_BOXFLUX',
        'HGAMMA_VSHIFT',
        'HGAMMA_SIGMA',
        'HGAMMA_CONT',
        'HGAMMA_CONT_IVAR',
        'HGAMMA_EW',
        'HGAMMA_EW_IVAR',
        'HGAMMA_FLUX_LIMIT',
        'HGAMMA_EW_LIMIT',
        'HGAMMA_CHI2',
        'HGAMMA_NPIX',
        'HGAMMA_BROAD_AMP',
        'HGAMMA_BROAD_AMP_IVAR',
        'HGAMMA_BROAD_FLUX',
        'HGAMMA_BROAD_FLUX_IVAR',
        'HGAMMA_BROAD_BOXFLUX',
        'HGAMMA_BROAD_VSHIFT',
        'HGAMMA_BROAD_SIGMA',
        'HGAMMA_BROAD_CONT',
        'HGAMMA_BROAD_CONT_IVAR',
        'HGAMMA_BROAD_EW',
        'HGAMMA_BROAD_EW_IVAR',
        'HGAMMA_BROAD_FLUX_LIMIT',
        'HGAMMA_BROAD_EW_LIMIT',
        'HGAMMA_BROAD_CHI2',
        'HGAMMA_BROAD_NPIX',
        'OIII_4363_AMP',
        'OIII_4363_AMP_IVAR',
        'OIII_4363_FLUX',
        'OIII_4363_FLUX_IVAR',
        'OIII_4363_BOXFLUX',
        'OIII_4363_VSHIFT',
        'OIII_4363_SIGMA',
        'OIII_4363_CONT',
        'OIII_4363_CONT_IVAR',
        'OIII_4363_EW',
        'OIII_4363_EW_IVAR',
        'OIII_4363_FLUX_LIMIT',
        'OIII_4363_EW_LIMIT',
        'OIII_4363_CHI2',
        'OIII_4363_NPIX',
        'HEI_4471_AMP',
        'HEI_4471_AMP_IVAR',
        'HEI_4471_FLUX',
        'HEI_4471_FLUX_IVAR',
        'HEI_4471_BOXFLUX',
        'HEI_4471_VSHIFT',
        'HEI_4471_SIGMA',
        'HEI_4471_CONT',
        'HEI_4471_CONT_IVAR',
        'HEI_4471_EW',
        'HEI_4471_EW_IVAR',
        'HEI_4471_FLUX_LIMIT',
        'HEI_4471_EW_LIMIT',
        'HEI_4471_CHI2',
        'HEI_4471_NPIX',
        'HEI_BROAD_4471_AMP',
        'HEI_BROAD_4471_AMP_IVAR',
        'HEI_BROAD_4471_FLUX',
        'HEI_BROAD_4471_FLUX_IVAR',
        'HEI_BROAD_4471_BOXFLUX',
        'HEI_BROAD_4471_VSHIFT',
        'HEI_BROAD_4471_SIGMA',
        'HEI_BROAD_4471_CONT',
        'HEI_BROAD_4471_CONT_IVAR',
        'HEI_BROAD_4471_EW',
        'HEI_BROAD_4471_EW_IVAR',
        'HEI_BROAD_4471_FLUX_LIMIT',
        'HEI_BROAD_4471_EW_LIMIT',
        'HEI_BROAD_4471_CHI2',
        'HEI_BROAD_4471_NPIX',
        'HEII_4686_AMP',
        'HEII_4686_AMP_IVAR',
        'HEII_4686_FLUX',
        'HEII_4686_FLUX_IVAR',
        'HEII_4686_BOXFLUX',
        'HEII_4686_VSHIFT',
        'HEII_4686_SIGMA',
        'HEII_4686_CONT',
        'HEII_4686_CONT_IVAR',
        'HEII_4686_EW',
        'HEII_4686_EW_IVAR',
        'HEII_4686_FLUX_LIMIT',
        'HEII_4686_EW_LIMIT',
        'HEII_4686_CHI2',
        'HEII_4686_NPIX',
        'HEII_BROAD_4686_AMP',
        'HEII_BROAD_4686_AMP_IVAR',
        'HEII_BROAD_4686_FLUX',
        'HEII_BROAD_4686_FLUX_IVAR',
        'HEII_BROAD_4686_BOXFLUX',
        'HEII_BROAD_4686_VSHIFT',
        'HEII_BROAD_4686_SIGMA',
        'HEII_BROAD_4686_CONT',
        'HEII_BROAD_4686_CONT_IVAR',
        'HEII_BROAD_4686_EW',
        'HEII_BROAD_4686_EW_IVAR',
        'HEII_BROAD_4686_FLUX_LIMIT',
        'HEII_BROAD_4686_EW_LIMIT',
        'HEII_BROAD_4686_CHI2',
        'HEII_BROAD_4686_NPIX',
        'HBETA_AMP',
        'HBETA_AMP_IVAR',
        'HBETA_FLUX',
        'HBETA_FLUX_IVAR',
        'HBETA_BOXFLUX',
        'HBETA_VSHIFT',
        'HBETA_SIGMA',
        'HBETA_CONT',
        'HBETA_CONT_IVAR',
        'HBETA_EW',
        'HBETA_EW_IVAR',
        'HBETA_FLUX_LIMIT',
        'HBETA_EW_LIMIT',
        'HBETA_CHI2',
        'HBETA_NPIX',
        'HBETA_BROAD_AMP',
        'HBETA_BROAD_AMP_IVAR',
        'HBETA_BROAD_FLUX',
        'HBETA_BROAD_FLUX_IVAR',
        'HBETA_BROAD_BOXFLUX',
        'HBETA_BROAD_VSHIFT',
        'HBETA_BROAD_SIGMA',
        'HBETA_BROAD_CONT',
        'HBETA_BROAD_CONT_IVAR',
        'HBETA_BROAD_EW',
        'HBETA_BROAD_EW_IVAR',
        'HBETA_BROAD_FLUX_LIMIT',
        'HBETA_BROAD_EW_LIMIT',
        'HBETA_BROAD_CHI2',
        'HBETA_BROAD_NPIX',
        'OIII_4959_AMP',
        'OIII_4959_AMP_IVAR',
        'OIII_4959_FLUX',
        'OIII_4959_FLUX_IVAR',
        'OIII_4959_BOXFLUX',
        'OIII_4959_VSHIFT',
        'OIII_4959_SIGMA',
        'OIII_4959_CONT',
        'OIII_4959_CONT_IVAR',
        'OIII_4959_EW',
        'OIII_4959_EW_IVAR',
        'OIII_4959_FLUX_LIMIT',
        'OIII_4959_EW_LIMIT',
        'OIII_4959_CHI2',
        'OIII_4959_NPIX',
        'OIII_5007_AMP',
        'OIII_5007_AMP_IVAR',
        'OIII_5007_FLUX',
        'OIII_5007_FLUX_IVAR',
        'OIII_5007_BOXFLUX',
        'OIII_5007_VSHIFT',
        'OIII_5007_SIGMA',
        'OIII_5007_CONT',
        'OIII_5007_CONT_IVAR',
        'OIII_5007_EW',
        'OIII_5007_EW_IVAR',
        'OIII_5007_FLUX_LIMIT',
        'OIII_5007_EW_LIMIT',
        'OIII_5007_CHI2',
        'OIII_5007_NPIX',
        'NII_5755_AMP',
        'NII_5755_AMP_IVAR',
        'NII_5755_FLUX',
        'NII_5755_FLUX_IVAR',
        'NII_5755_BOXFLUX',
        'NII_5755_VSHIFT',
        'NII_5755_SIGMA',
        'NII_5755_CONT',
        'NII_5755_CONT_IVAR',
        'NII_5755_EW',
        'NII_5755_EW_IVAR',
        'NII_5755_FLUX_LIMIT',
        'NII_5755_EW_LIMIT',
        'NII_5755_CHI2',
        'NII_5755_NPIX',
        'HEI_5876_AMP',
        'HEI_5876_AMP_IVAR',
        'HEI_5876_FLUX',
        'HEI_5876_FLUX_IVAR',
        'HEI_5876_BOXFLUX',
        'HEI_5876_VSHIFT',
        'HEI_5876_SIGMA',
        'HEI_5876_CONT',
        'HEI_5876_CONT_IVAR',
        'HEI_5876_EW',
        'HEI_5876_EW_IVAR',
        'HEI_5876_FLUX_LIMIT',
        'HEI_5876_EW_LIMIT',
        'HEI_5876_CHI2',
        'HEI_5876_NPIX',
        'HEI_BROAD_5876_AMP',
        'HEI_BROAD_5876_AMP_IVAR',
        'HEI_BROAD_5876_FLUX',
        'HEI_BROAD_5876_FLUX_IVAR',
        'HEI_BROAD_5876_BOXFLUX',
        'HEI_BROAD_5876_VSHIFT',
        'HEI_BROAD_5876_SIGMA',
        'HEI_BROAD_5876_CONT',
        'HEI_BROAD_5876_CONT_IVAR',
        'HEI_BROAD_5876_EW',
        'HEI_BROAD_5876_EW_IVAR',
        'HEI_BROAD_5876_FLUX_LIMIT',
        'HEI_BROAD_5876_EW_LIMIT',
        'HEI_BROAD_5876_CHI2',
        'HEI_BROAD_5876_NPIX',
        'OI_6300_AMP',
        'OI_6300_AMP_IVAR',
        'OI_6300_FLUX',
        'OI_6300_FLUX_IVAR',
        'OI_6300_BOXFLUX',
        'OI_6300_VSHIFT',
        'OI_6300_SIGMA',
        'OI_6300_CONT',
        'OI_6300_CONT_IVAR',
        'OI_6300_EW',
        'OI_6300_EW_IVAR',
        'OI_6300_FLUX_LIMIT',
        'OI_6300_EW_LIMIT',
        'OI_6300_CHI2',
        'OI_6300_NPIX',
        'SIII_6312_AMP',
        'SIII_6312_AMP_IVAR',
        'SIII_6312_FLUX',
        'SIII_6312_FLUX_IVAR',
        'SIII_6312_BOXFLUX',
        'SIII_6312_VSHIFT',
        'SIII_6312_SIGMA',
        'SIII_6312_CONT',
        'SIII_6312_CONT_IVAR',
        'SIII_6312_EW',
        'SIII_6312_EW_IVAR',
        'SIII_6312_FLUX_LIMIT',
        'SIII_6312_EW_LIMIT',
        'SIII_6312_CHI2',
        'SIII_6312_NPIX',
        'NII_6548_AMP',
        'NII_6548_AMP_IVAR',
        'NII_6548_FLUX',
        'NII_6548_FLUX_IVAR',
        'NII_6548_BOXFLUX',
        'NII_6548_VSHIFT',
        'NII_6548_SIGMA',
        'NII_6548_CONT',
        'NII_6548_CONT_IVAR',
        'NII_6548_EW',
        'NII_6548_EW_IVAR',
        'NII_6548_FLUX_LIMIT',
        'NII_6548_EW_LIMIT',
        'NII_6548_CHI2',
        'NII_6548_NPIX',
        'HALPHA_AMP',
        'HALPHA_AMP_IVAR',
        'HALPHA_FLUX',
        'HALPHA_FLUX_IVAR',
        'HALPHA_BOXFLUX',
        'HALPHA_VSHIFT',
        'HALPHA_SIGMA',
        'HALPHA_CONT',
        'HALPHA_CONT_IVAR',
        'HALPHA_EW',
        'HALPHA_EW_IVAR',
        'HALPHA_FLUX_LIMIT',
        'HALPHA_EW_LIMIT',
        'HALPHA_CHI2',
        'HALPHA_NPIX',
        'HALPHA_BROAD_AMP',
        'HALPHA_BROAD_AMP_IVAR',
        'HALPHA_BROAD_FLUX',
        'HALPHA_BROAD_FLUX_IVAR',
        'HALPHA_BROAD_BOXFLUX',
        'HALPHA_BROAD_VSHIFT',
        'HALPHA_BROAD_SIGMA',
        'HALPHA_BROAD_CONT',
        'HALPHA_BROAD_CONT_IVAR',
        'HALPHA_BROAD_EW',
        'HALPHA_BROAD_EW_IVAR',
        'HALPHA_BROAD_FLUX_LIMIT',
        'HALPHA_BROAD_EW_LIMIT',
        'HALPHA_BROAD_CHI2',
        'HALPHA_BROAD_NPIX',
        'NII_6584_AMP',
        'NII_6584_AMP_IVAR',
        'NII_6584_FLUX',
        'NII_6584_FLUX_IVAR',
        'NII_6584_BOXFLUX',
        'NII_6584_VSHIFT',
        'NII_6584_SIGMA',
        'NII_6584_CONT',
        'NII_6584_CONT_IVAR',
        'NII_6584_EW',
        'NII_6584_EW_IVAR',
        'NII_6584_FLUX_LIMIT',
        'NII_6584_EW_LIMIT',
        'NII_6584_CHI2',
        'NII_6584_NPIX',
        'SII_6716_AMP',
        'SII_6716_AMP_IVAR',
        'SII_6716_FLUX',
        'SII_6716_FLUX_IVAR',
        'SII_6716_BOXFLUX',
        'SII_6716_VSHIFT',
        'SII_6716_SIGMA',
        'SII_6716_CONT',
        'SII_6716_CONT_IVAR',
        'SII_6716_EW',
        'SII_6716_EW_IVAR',
        'SII_6716_FLUX_LIMIT',
        'SII_6716_EW_LIMIT',
        'SII_6716_CHI2',
        'SII_6716_NPIX',
        'SII_6731_AMP',
        'SII_6731_AMP_IVAR',
        'SII_6731_FLUX',
        'SII_6731_FLUX_IVAR',
        'SII_6731_BOXFLUX',
        'SII_6731_VSHIFT',
        'SII_6731_SIGMA',
        'SII_6731_CONT',
        'SII_6731_CONT_IVAR',
        'SII_6731_EW',
        'SII_6731_EW_IVAR',
        'SII_6731_FLUX_LIMIT',
        'SII_6731_EW_LIMIT',
        'SII_6731_CHI2',
        'SII_6731_NPIX',
        'OII_7320_AMP',
        'OII_7320_AMP_IVAR',
        'OII_7320_FLUX',
        'OII_7320_FLUX_IVAR',
        'OII_7320_BOXFLUX',
        'OII_7320_VSHIFT',
        'OII_7320_SIGMA',
        'OII_7320_CONT',
        'OII_7320_CONT_IVAR',
        'OII_7320_EW',
        'OII_7320_EW_IVAR',
        'OII_7320_FLUX_LIMIT',
        'OII_7320_EW_LIMIT',
        'OII_7320_CHI2',
        'OII_7320_NPIX',
        'OII_7330_AMP',
        'OII_7330_AMP_IVAR',
        'OII_7330_FLUX',
        'OII_7330_FLUX_IVAR',
        'OII_7330_BOXFLUX',
        'OII_7330_VSHIFT',
        'OII_7330_SIGMA',
        'OII_7330_CONT',
        'OII_7330_CONT_IVAR',
        'OII_7330_EW',
        'OII_7330_EW_IVAR',
        'OII_7330_FLUX_LIMIT',
        'OII_7330_EW_LIMIT',
        'OII_7330_CHI2',
        'OII_7330_NPIX',
        'SIII_9069_AMP',
        'SIII_9069_AMP_IVAR',
        'SIII_9069_FLUX',
        'SIII_9069_FLUX_IVAR',
        'SIII_9069_BOXFLUX',
        'SIII_9069_VSHIFT',
        'SIII_9069_SIGMA',
        'SIII_9069_CONT',
        'SIII_9069_CONT_IVAR',
        'SIII_9069_EW',
        'SIII_9069_EW_IVAR',
        'SIII_9069_FLUX_LIMIT',
        'SIII_9069_EW_LIMIT',
        'SIII_9069_CHI2',
        'SIII_9069_NPIX',
        'SIII_9532_AMP',
        'SIII_9532_AMP_IVAR',
        'SIII_9532_FLUX',
        'SIII_9532_FLUX_IVAR',
        'SIII_9532_BOXFLUX',
        'SIII_9532_VSHIFT',
        'SIII_9532_SIGMA',
        'SIII_9532_CONT',
        'SIII_9532_CONT_IVAR',
        'SIII_9532_EW',
        'SIII_9532_EW_IVAR',
        'SIII_9532_FLUX_LIMIT',
        'SIII_9532_EW_LIMIT',
        'SIII_9532_CHI2',
        'SIII_9532_NPIX'
        ]

    F = fitsio.FITS(fastspecfile)
    nrows = F['METADATA'].get_nrows()

    if False:
        rows = np.arange(10000) + 10000
    elif False:
        rand = np.random.RandomState(seed=1)
        rows = rand.choice(nrows, size=75, replace=False)
        np.sort(rows)
    else:
        rows = np.arange(nrows)

    meta = Table(fitsio.read(fastspecfile, ext='METADATA', columns=meta_columns, rows=rows))

    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC', columns=fastspec_cols, rows=rows))
    #print('Ignoring fastphot!!')
    #if False:
    #    fastphot = Table(fitsio.read(fastphotfile, ext='FASTPHOT', columns=fastphot_cols))
    #    assert(np.all(meta['TARGETID'] == fastphot['TARGETID']))

    # This will be different for healpix vs tile coadds. E.g., sv3-bright-HEALPIX-TARGETID
    meta['TARGET_NAME'] = ['{}-{}-{}-{}'.format(survey, program, healpix, targetid) for
                           survey, program, healpix, targetid in zip(
                               meta['SURVEY'], meta['PROGRAM'], meta['HEALPIX'], meta['TARGETID'])]
    meta['SPECPROD'] = specprod
    #print(meta)
    #print(meta.colnames, fast.colnames)

    # parse the targeting bit names
    desi_bitnames = np.zeros(len(meta), dtype='U300')
    bgs_bitnames = np.zeros(len(meta), dtype='U300')
    mws_bitnames = np.zeros(len(meta), dtype='U300')
    scnd_bitnames = np.zeros(len(meta), dtype='U300')
    cmx_bitnames = np.zeros(len(meta), dtype='U300')
    targetclass = np.zeros(len(meta), dtype='U300')

    def get_targetclass(targetclass, name):
        for cc in ['BGS', 'LRG', 'ELG', 'QSO', 'MWS', 'SCND', 'STD']:
            if cc in name:
                for iobj, tclass in enumerate(targetclass):
                    if tclass == '':
                        targetclass[iobj] = cc
                    else:
                        if not cc in tclass: # only once
                            targetclass[iobj] = ' '.join([tclass, cc])
        return targetclass
    
    for survey, prefix in zip(['CMX', 'SV1', 'SV2', 'SV3', 'MAIN'], ['CMX_', 'SV1_', 'SV2_', 'SV3_', '']):
        I = meta['SURVEY'] == survey.lower()
        if np.sum(I) > 0:
            if survey == 'MAIN':
                from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
            elif survey == 'SV1':
                from desitarget.sv1.sv1_targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
            elif survey == 'SV2':
                from desitarget.sv2.sv2_targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
            elif survey == 'SV3':
                from desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
            elif survey == 'CMX':
                from desitarget.cmx.cmx_targetmask import cmx_mask

            if survey == 'CMX':
                for name in cmx_mask.names():
                    J = np.where(meta['CMX_TARGET'.format(prefix)] & cmx_mask.mask(name) != 0)[0]
                    if len(J) > 0:
                        cmx_bitnames[J] = [' '.join([bit, name]) for bit in cmx_bitnames[J]]
                        #if 'QSO' in name:
                        #    pdb.set_trace()
                        #print(name, targetclass[J])
                        targetclass[J] = get_targetclass(targetclass[J], name)
            else:
                for name in desi_mask.names():
                    J = np.where(meta['{}DESI_TARGET'.format(prefix)] & desi_mask.mask(name) != 0)[0]
                    if len(J) > 0:
                        desi_bitnames[J] = [' '.join([bit, name]) for bit in desi_bitnames[J]]
                        targetclass[J] = get_targetclass(targetclass[J], name)
                        
                for name in bgs_mask.names():
                    J = np.where(meta['{}BGS_TARGET'.format(prefix)] & bgs_mask.mask(name) != 0)[0]
                    if len(J) > 0:
                        bgs_bitnames[J] = [' '.join([bit, name]) for bit in bgs_bitnames[J]]
                        targetclass[J] = get_targetclass(targetclass[J], name)
                        
                for name in mws_mask.names():
                    J = np.where(meta['{}MWS_TARGET'.format(prefix)] & mws_mask.mask(name) != 0)[0]
                    if len(J) > 0:
                        mws_bitnames[J] = [' '.join([bit, name]) for bit in mws_bitnames[J]]
                        targetclass[J] = get_targetclass(targetclass[J], name)
                        
                for name in scnd_mask.names():
                    J = np.where(meta['{}SCND_TARGET'.format(prefix)] & scnd_mask.mask(name) != 0)[0]
                    if len(J) > 0:
                        scnd_bitnames[J] = [' '.join([bit, name]) for bit in scnd_bitnames[J]]
                        targetclass[J] = get_targetclass(targetclass[J], name)

    meta['DESI_BITNAMES'] = desi_bitnames
    meta['BGS_BITNAMES'] = bgs_bitnames
    meta['MWS_BITNAMES'] = mws_bitnames
    meta['SCND_BITNAMES'] = scnd_bitnames
    meta['CMX_BITNAMES'] = cmx_bitnames
    meta['TARGETCLASS'] = targetclass

    # rename a couple columns
    meta.rename_column('Z_RR', 'ZREDROCK')
            
    # join metadata and fastspec fitting results
    data = hstack((meta, fastspec))
    del meta, fastspec

    ## join everything with fastphot fitting results but we need to add a prefix
    #print('Ignoring fastphot!!')
    #if False:
    #    for col in fastphot.colnames:
    #        fastphot.rename_column(col, 'PHOT_{}'.format(col))
    #    fastphot.remove_column('PHOT_TARGETID')
    #    data = hstack((data, fastphot))

    # Parse the photometry into strings with limits. Most of this code is taken
    # from continuum.ContinuumTools.parse_photometry deal with measurements
    def convert_phot(data, prefix, suffix, bands, nsigma=2.0):
        for band in bands:
            maggies = data['{}FLUX{}_{}'.format(prefix, suffix, band)]
            if 'FIBER' in prefix or 'SYNTH' in suffix:
                # e.g., FIBERABMAG_G, FIBERTOTMAG_G, ABMAG_SYNTH_G, and ABMAG_SYNTH_PHOTMODEL_G, and ABMAG_SYNTH_SPECMODEL_G
                data['{}ABMAG{}_{}'.format(prefix, suffix, band)] = np.array('', dtype='U50') 
                good = np.where(maggies > 0)[0]
                if len(good) > 0:
                    data['{}ABMAG{}_{}'.format(prefix, suffix, band)][good] = np.array(list(
                        map(lambda x: '{:.3f}'.format(x), -2.5 * np.log10(1e-9 * maggies[good]))))
                neg = np.where(maggies <= 0)[0]
                if len(neg) > 0:
                    data['{}ABMAG{}_{}'.format(prefix, suffix, band)][neg] = '...'
            else:
                # e.g., ABMAG_G and ABMAG_ERR_G
                data['ABMAG_{}'.format(band)] = np.array('', dtype='U50')
                data['ABMAG_ERR_{}'.format(band)] = np.array('', dtype='U50')
            
                ivarmaggies = data['FLUX_IVAR_{}'.format(band)]
        
                # deal with the uncertainties
                snr = maggies * np.sqrt(ivarmaggies)
        
                # upper limits
                upper = np.where((ivarmaggies > 0) * (snr <= nsigma))[0]
                if len(upper) > 0:
                    data['ABMAG_{}'.format(band)][upper] = np.array(list(map(lambda x: '>{:.3f}'.format(x), -2.5 * np.log10(1e-9 * nsigma / np.sqrt(ivarmaggies[upper])))))
                            
                # significant detections
                good = np.where(snr > nsigma)[0]
                if len(good) > 0:
                    ## approximate the uncertainty as being symmetric in magnitude
                    abmag_ivar = (ivarmaggies[good] * (maggies[good] * 0.4 * np.log(10))**2)
                    abmag_err = 1.0 / np.sqrt(abmag_ivar)
        
                    #errmaggies = 1 / np.sqrt(ivarmaggies[good])
                    #abmag_brighter = errmaggies / (0.4 * np.log(10) * (maggies[good]+errmaggies)) # bright end (flux upper limit)
                    #abmag_fainterr = errmaggies / (0.4 * np.log(10) * (maggies[good]-errmaggies)) # faint end (flux lower limit)
        
                    data['ABMAG_{}'.format(band)][good] = np.array(list(map(lambda x: '{:.3f}'.format(x), -2.5 * np.log10(1e-9 * maggies[good]))))
                    data['ABMAG_ERR_{}'.format(band)][good] = np.array(list(map(lambda x: '&#177;{:.3f}'.format(x), abmag_err)))
                    
        return data

    for prefix, suffix, bands in zip(
            ['', 'FIBER', 'FIBERTOT', '', '', ''],
            ['', '', '', '_SYNTH', '_SYNTH_SPECMODEL', '_SYNTH_PHOTMODEL'],
            [['G', 'R', 'Z', 'W1', 'W2', 'W3', 'W4'],
             ['G', 'R', 'Z'],
             ['G', 'R', 'Z'],
             ['G', 'R', 'Z'],
             ['G', 'R', 'Z'],
             ['G', 'R', 'Z', 'W1', 'W2', 'W3', 'W4']
             ]):
        data = convert_phot(data, prefix, suffix, bands)

    # parse some uncertainties
    data['VDISP_ERR'] = np.zeros(len(data), dtype='U50')
    I = data['VDISP_IVAR'] > 0
    if np.any(I):
        vdisp_err = 1 / np.sqrt(data['VDISP_IVAR'][I])
        data['VDISP_ERR'][I] = np.array(list(map(lambda x: '&#177;{:.0f}'.format(x).strip(), vdisp_err)))

    for suffix in ['NARROW', 'BROAD', 'UV']:
        data['{}_DV'.format(suffix)] = C_LIGHT*(data['{}_Z'.format(suffix)]-data['Z'])
        data['{}_DV_ERR'.format(suffix)] = np.array(list(map(lambda x: '&#177;{:.0f}'.format(x).strip(), C_LIGHT*data['{}_ZRMS'.format(suffix)])))
        
    for suffix in ['NARROW', 'BROAD', 'UV']:
        data['{}_SIGMA_ERR'.format(suffix)] = np.array(list(map(lambda x: '&#177;{:.0f}'.format(x).strip(), data['{}_SIGMARMS'.format(suffix)])))

    # get uncertainties on the emission-line measurements
    def val_plusminus(val, err):
        ret = '{:.3g}&#177;{:.3g}'.format(val, err)
        #if len(ret) > 16:
        #    import pdb ; pdb.set_trace()
        return ret

    linefile = resource_filename('fastspecfit', 'data/emlines.ecsv')    
    linetable = Table.read(linefile, format='ascii.ecsv', guess=False)
    lines = np.array([name.upper() for name in linetable['name'].data])
    #lines1 = np.array([col[:-4] for col in data.colnames if col[-4:] == '_AMP'])
    
    for line, restwave in zip(lines, linetable['restwave']):
        linewave = '{:.4f}'.format(restwave)
        if linewave[-2:] == '00':
            linewave = linewave[:-2]
        if linewave[-1:] == '0':
            linewave = linewave[:-1]
        data['{}_WAVE'.format(line)] = linewave
        
        # S/N, vshift, sigma, and chi2
        data['{}_SNR'.format(line)] = np.array('', dtype='U50')
        data['{}_VSHIFT_STR'.format(line)] = np.array('', dtype='U50')
        data['{}_SIGMA_STR'.format(line)] = np.array('', dtype='U50')
        data['{}_CHI2_STR'.format(line)] = np.array('', dtype='U50')

        snr = data['{}_AMP'.format(line)]*np.sqrt(data['{}_AMP_IVAR'.format(line)])
        good = np.where(snr > 0)[0]
        if len(good) > 0:
            data['{}_SNR'.format(line)][good] = np.array(list(map(lambda x: '{:.3g}'.format(x), snr[good])))
            data['{}_VSHIFT_STR'.format(line)][good] = np.array(list(map(lambda x: '{:.4g}'.format(x), data['{}_VSHIFT'.format(line)][good])))
            data['{}_SIGMA_STR'.format(line)][good] = np.array(list(map(lambda x: '{:.4g}'.format(x), data['{}_SIGMA'.format(line)][good])))
            data['{}_CHI2_STR'.format(line)][good] = np.array(list(map(lambda x: '{:.3g}'.format(x), data['{}_CHI2'.format(line)][good])))
            
        zero = np.where(snr == 0)[0]
        if len(zero) > 0:
            data['{}_SNR'.format(line)][zero] = '0'
            data['{}_VSHIFT_STR'.format(line)][zero] = '0'
            data['{}_SIGMA_STR'.format(line)][zero] = '0'
            data['{}_CHI2_STR'.format(line)][zero] = '0'

        # amp, flux, cont, and EW
        for suffix in ['AMP', 'FLUX', 'CONT', 'EW']:
            data['{}_{}_ERR'.format(line, suffix)] = np.array('', dtype='U50')
            good = np.where(data['{}_{}_IVAR'.format(line, suffix)] > 0)[0]
            if len(good) > 0:
                data['{}_{}_ERR'.format(line, suffix)][good] = np.array(list(map(val_plusminus, data['{}_{}'.format(line, suffix)][good],
                                                                                 1 / np.sqrt(data['{}_{}_IVAR'.format(line, suffix)][good]) )))
            zero = np.where(data['{}_{}_IVAR'.format(line, suffix)] == 0)[0]
            if len(zero) > 0:
                data['{}_{}_ERR'.format(line, suffix)][zero] = '0'
                #print(data['{}_{}_ERR'.format(line, suffix)][good])
                
    print(data.colnames)
    print('Read {} rows from {}'.format(len(data), fastspecfile))
    xyz = radectoxyz(data['RA'], data['DEC'])

    # Load the database in chunks to avoid memory issues.
    chunksize = 50000
    chunkindx = np.arange(len(data))
    chunks = np.array_split(chunkindx, np.ceil(len(data) / chunksize).astype(int))

    for ichunk, chunk in enumerate(chunks):
        objs = []
        for ii in chunk:
            sam = FastModel()
            sam.row_index = ii
            sam.ux = xyz[ii, 0]
            sam.uy = xyz[ii, 1]
            sam.uz = xyz[ii, 2]
            for col in data.colnames:
                val = data[ii][col]
                if type(val) == np.str_:
                    val = val.strip()
                setattr(sam, col.lower(), val)
            objs.append(sam)
            
        print('Loading chunk {} / {} (N={}) into the database.'.format(
            ichunk+1, len(chunks), len(chunk)))
        FastModel.objects.bulk_create(objs)

    del data

#    ##################
#    from fastspecfit.webapp.sample.models import TestSample
#    for ichunk, chunk in enumerate(chunks[7:8]):
#        objs = []
#        for ii in chunk[551:552]:
#            sam = TestSample()
#            for col in data.colnames:
#                if hasattr(sam, col.lower()):
#                    val = data[ii][col]
#                    if type(val) == np.str_:
#                        val = val.strip()
#                    setattr(sam, col.lower(), val)
#            objs.append(sam)
#            
#        print('Loading chunk {} / {} (N={}) into the database.'.format(
#            ichunk+1, len(chunks), len(chunk)))
#        try:
#            TestSample.objects.bulk_create(objs)
#        except:
#            #for field in TestSample._meta.get_fields():
#            #    print([field.name, getattr(field, 'max_length', 'Not Supported'), type(field)])
#            pdb.set_trace()
#
#    print('Success')
#    pdb.set_trace()
#    del data
#
#
#    ##################
#
#    for ichunk, chunk in enumerate(chunks[7:8]):
#        objs = []
#        for ii in chunk[551:552]:
#            sam = Sample()
#            sam.row_index = ii
#            sam.ux = xyz[ii, 0]
#            sam.uy = xyz[ii, 1]
#            sam.uz = xyz[ii, 2]
#            for col in data.colnames:
#                val = data[ii][col]
#                if type(val) == np.str_:
#                    val = val.strip()
#                    #if 'BITNAMES' in col:
#                    if data[col].dtype.num == 4:
#                        print(ii, col, val, len(val), data[col].dtype)
#                    #if len(val) > data[col].dtype.num:
#                    #    print(ii, col, val, len(val), data[col].dtype)
#                setattr(sam, col.lower(), val)
#            objs.append(sam)
#            
#        print('Loading chunk {} / {} (N={}) into the database.'.format(
#            ichunk+1, len(chunks), len(chunk)))
#        try:
#            Sample.objects.bulk_create(objs)
#        except:
#            for field in Sample._meta.get_fields():
#                print([field.name, getattr(field, 'max_length', 'Not Supported'), type(field)])
#            pdb.set_trace()
#            for obj in objs:
#                l = len(obj['your_problematic_field'])
#                if len(l) > 7:
#                    print(f'obj with id {obj["pk"]} has that field at {l} characters long')
#                        
#            pdb.set_trace()
#
#    print('Success')
#    pdb.set_trace()
#    del data

if __name__ == '__main__':
    main()
