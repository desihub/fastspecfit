.. _fastspec datamodel:

===================
Fastspec Data Model
===================

:Summary: Spectroscopic fitting results.
:Naming Convention:
    ``fastspec-SURVEY-PROGRAM-HEALPIX.fits.gz``, where
    ``SURVEY`` is the survey name (e.g., *main* or *sv3*), ``PROGRAM`` is the
    program name (e.g., *bright* or *dark*), and ``HEALPIX`` is the healpixel number.
:Regex: ``fastspec-(cmx|main|special|sv1|sv2|sv3)-(backup|bright|dark|other)-[0-9]+\.fits.gz``
:File Type: FITS

Contents
========

====== ============ ======== ======================
Number EXTNAME      Type     Contents
====== ============ ======== ======================
HDU00_ PRIMARY      IMAGE    Keywords only.
HDU01_ FASTSPEC     BINTABLE Table with spectral fitting results.
HDU02_ METADATA     BINTABLE Table with sample metadata.
HDU03_ MODELS       IMAGE    Model spectra.
====== ============ ======== ======================

FITS Header Units
=================

HDU00
-----

EXTNAME = PRIMARY

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

.. collapse:: Required Header Keywords Table
    :open:

    .. rst-class:: keywords

    ======== ================ ==== ==============================================
    KEY      Example Value    Type Comment
    ======== ================ ==== ==============================================
    SPECPROD iron             str  Spectroscopic production name
    COADDTYP healpix          str  Type of spectral coadd (healpix,cumulative,pernight,perexp)
    CHECKSUM ZHWHZHVFZHVFZHVF str  HDU checksum updated 2022-08-02T19:06:22
    DATASUM  0                str  data unit checksum updated 2022-08-02T19:06:22
    ======== ================ ==== ==============================================

Empty HDU.

HDU01
-----

EXTNAME = FASTSPEC

Spectroscopic fitting results.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

.. collapse:: Required Header Keywords Table
    :open:

    .. rst-class:: keywords

    ======== ================ ==== ==============================================
    KEY      Example Value    Type Comment
    ======== ================ ==== ==============================================
    NAXIS1   3115             int  length of dimension 1
    NAXIS2   338              int  length of dimension 2
    CHECKSUM WLDnaKDlWKDlaKDl str  HDU checksum updated 2022-08-03T14:25:08
    DATASUM  756558178        str  data unit checksum updated 2022-08-03T14:25:08
    ======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: columns

============================ ============ ============================= ============================================
Name                         Type         Units                         Description
============================ ============ ============================= ============================================
                    TARGETID        int64                               Unique target ID.
                 SURVEY [1]_       bytes3                               Survey name.
                PROGRAM [1]_       bytes6                               Program name.
                HEALPIX [1]_        int32                               Healpixel number.
                 TILEID [2]_        int32                               Tile ID number.
                  NIGHT [2]_        int32                               Night of observation.
                  FIBER [2]_        int32                               Fiber number.
                  EXPID [3]_        int32                               Exposure ID number.
                           Z      float64                               Stellar continuum redshift.
                       COEFF float32[168]                               Continuum coefficients.
                       RCHI2      float32                               Reduced chi-squared of the full-spectrum fit.
                  RCHI2_CONT      float32                               Reduced chi-squared of the fit to the stellar continuum.
                  RCHI2_PHOT      float32                               Reduced chi-squared of the fit to the broadband photometry.
                       SNR_B      float32                               Median signal-to-noise ratio per pixel in *b* camera.
                       SNR_R      float32                               Median signal-to-noise ratio per pixel in *r* camera.
                       SNR_Z      float32                               Median signal-to-noise ratio per pixel in *z* camera.
                SMOOTHCORR_B      float32                       percent Mean value of the smooth continuum correction divided by the best-fitting continuum model in the *b* camera.
                SMOOTHCORR_R      float32                       percent Mean value of the smooth continuum correction divided by the best-fitting continuum model in the *r* camera.
                SMOOTHCORR_Z      float32                       percent Mean value of the smooth continuum correction divided by the best-fitting continuum model in the *z* camera.
                       VDISP      float32                        km / s Stellar velocity dispersion.
                  VDISP_IVAR      float32                      s2 / km2 Inverse variance of VDISP.
                          AV      float32                           mag Attenuation of the integrated stellar population.
                         AGE      float32                           Gyr Light-weighted age.
                       ZZSUN      float32                               Logarithmic stellar metallicity relative to solar.
                    LOGMSTAR      float32                          Msun Logarithmic stellar mass (h=1.0, Chabrier+2003 initial mass function).
                         SFR      float32                     Msun / yr Instantaneous star formation rate (h=1.0, Chabrier+2003 initial mass function).
                      DN4000      float32                               Narrow 4000-A break index (from Balogh et al. 1999) measured from the emission-line subtracted data.
                  DN4000_OBS      float32                               Narrow 4000-A break index measured from the observed spectroscopic data.
                 DN4000_IVAR      float32                               Inverse variance of DN4000_OBS or DN4000.
                DN4000_MODEL      float32                               Narrow 4000-A break index measured from the best-fitting continuum model.
                FLUX_SYNTH_G      float32                          nmgy g-band flux synthesized from the data.
                FLUX_SYNTH_R      float32                          nmgy r-band flux synthesized from the data.
                FLUX_SYNTH_Z      float32                          nmgy z-band flux synthesized from the data.
      FLUX_SYNTH_SPECMODEL_G      float32                          nmgy g-band flux synthesized from the best-fitting spectroscopic model.
      FLUX_SYNTH_SPECMODEL_R      float32                          nmgy r-band flux synthesized from the best-fitting spectroscopic model.
      FLUX_SYNTH_SPECMODEL_Z      float32                          nmgy z-band flux synthesized from the best-fitting spectroscopic model.
      FLUX_SYNTH_PHOTMODEL_G      float32                          nmgy g-band flux synthesized from the best-fitting photometric continuum model.
      FLUX_SYNTH_PHOTMODEL_R      float32                          nmgy r-band flux synthesized from the best-fitting photometric continuum model.
      FLUX_SYNTH_PHOTMODEL_Z      float32                          nmgy z-band flux synthesized from the best-fitting photometric continuum model.
     FLUX_SYNTH_PHOTMODEL_W1      float32                          nmgy W1-band flux synthesized from the best-fitting photometric continuum model.
     FLUX_SYNTH_PHOTMODEL_W2      float32                          nmgy W2-band flux synthesized from the best-fitting photometric continuum model.
     FLUX_SYNTH_PHOTMODEL_W3      float32                          nmgy W3-band flux synthesized from the best-fitting photometric continuum model.
     FLUX_SYNTH_PHOTMODEL_W4      float32                          nmgy W4-band flux synthesized from the best-fitting photometric continuum model.
                     KCORR_U      float32                           mag K-correction used to derive ABSMAG_U band-shifted to z=0.0 assuming h=1.0.
                    ABSMAG_U      float32                           mag Absolute magnitude in Johnson/Cousins U-band band-shifted to z=0.0 assuming h=1.0.
               ABSMAG_IVAR_U      float32                      1 / mag2 Inverse variance corresponding to ABSMAG_U.
                     KCORR_B      float32                           mag Like KCORR_U but for Johnson/Cousins B-band.
                    ABSMAG_B      float32                           mag Like ABSMAG_U but for Johnson/Cousins B-band.
               ABSMAG_IVAR_B      float32                      1 / mag2 Like ABSMAG_IVAR_U but for Johnson/Cousins B-band.
                     KCORR_V      float32                           mag Like KCORR_U but for Johnson/Cousins V-band.
                    ABSMAG_V      float32                           mag Like ABSMAG_U but for Johnson/Cousins V-band.
               ABSMAG_IVAR_V      float32                      1 / mag2 Like ABSMAG_IVAR_U but for Johnson/Cousins V-band.
                KCORR_SDSS_U      float32                           mag K-correction used to derive ABSMAG_SDSS_U band-shifted to z=0.1 assuming h=1.0.
               ABSMAG_SDSS_U      float32                           mag Absolute magnitude in SDSS u-band band-shifted to z=0.1 assuming h=1.0.
          ABSMAG_IVAR_SDSS_U      float32                      1 / mag2 Inverse variance corresponding to ABSMAG_SDSS_U.
                KCORR_SDSS_G      float32                           mag Like KCORR_SDSS_U but for SDSS g-band.
               ABSMAG_SDSS_G      float32                           mag Like ABSMAG_SDSS_U but for SDSS g-band.
          ABSMAG_IVAR_SDSS_G      float32                      1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS g-band.
                KCORR_SDSS_R      float32                           mag Like KCORR_SDSS_U but for SDSS r-band.
               ABSMAG_SDSS_R      float32                           mag Like ABSMAG_SDSS_U but for SDSS r-band.
          ABSMAG_IVAR_SDSS_R      float32                      1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS r-band.
                KCORR_SDSS_I      float32                           mag Like KCORR_SDSS_U but for SDSS i-band.
               ABSMAG_SDSS_I      float32                           mag Like ABSMAG_SDSS_U but for SDSS i-band.
          ABSMAG_IVAR_SDSS_I      float32                      1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS i-band.
                KCORR_SDSS_Z      float32                           mag Like KCORR_SDSS_U but for SDSS z-band.
               ABSMAG_SDSS_Z      float32                           mag Like ABSMAG_SDSS_U but for SDSS z-band.
          ABSMAG_IVAR_SDSS_Z      float32                      1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS z-band.
                    KCORR_W1      float32                           mag K-correction used to derive ABSMAG_W1 band-shifted to z=0.0 assuming h=1.0.
                   ABSMAG_W1      float32                           mag Absolute magnitude in WISE W1-band band-shifted to z=0.0 assuming h=1.0.
              ABSMAG_IVAR_W1      float32                      1 / mag2 Inverse variance corresponding to ABSMAG_W1.
                    KCORR_W2      float32                           mag K-correction used to derive ABSMAG_W2 band-shifted to z=0.0 assuming h=1.0.
                   ABSMAG_W2      float32                           mag Absolute magnitude in WISE W2-band band-shifted to z=0.0 assuming h=1.0.
              ABSMAG_IVAR_W2      float32                      1 / mag2 Inverse variance corresponding to ABSMAG_W2.
                 LOGLNU_1500      float32            1e-28 erg / (s Hz) Monochromatic luminosity at 1500 A in the rest-frame.
                 LOGLNU_2800      float32            1e-28 erg / (s Hz) Monochromatic luminosity at 2800 A in the rest-frame.
                   LOGL_5100      float32                    1e+10 Lsun Total luminosity at 5100 A in the rest-frame.
              FOII_3727_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 3728.483 A in the rest-frame.
                 FHBETA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 4862.683 A in the rest-frame.
             FOIII_5007_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 5008.239 A in the rest-frame.
                FHALPHA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 6564.613 A in the rest-frame.
                  RCHI2_LINE      float32                               Reduced chi-squared of an emission-line model fit.
             DELTA_LINERCHI2      float32                               Difference in the reduced chi-squared values between an emission-line model with narrow lines only and a model with both broad and narrow lines.
                     APERCOR      float32                               Median aperture correction factor.
                   APERCOR_G      float32                               Median aperture correction factor measured in the g-band.
                   APERCOR_R      float32                               Median aperture correction factor measured in the r-band.
                   APERCOR_Z      float32                               Median aperture correction factor measured in the z-band.
                    NARROW_Z      float32                               Mean redshift of well-measured narrow rest-frame optical emission lines (defaults to Z).
                 NARROW_ZRMS      float32                               Root-mean-square scatter in NARROW_Z.
                     BROAD_Z      float32                               Mean redshift of well-measured broad rest-frame optical emission lines (defaults to Z).
                  BROAD_ZRMS      float32                               Root-mean-square scatter in BROAD_Z.
                        UV_Z      float32                               Mean redshift of well-measured rest-frame UV emission lines (defaults to Z).
                     UV_ZRMS      float32                               Root-mean-square scatter in UV_Z.
                NARROW_SIGMA      float32                        km / s Mean line-width of well-measured narrow rest-frame optical emission lines.
             NARROW_SIGMARMS      float32                        km / s Root-mean-square scatter in NARROW_SIGMA.
                 BROAD_SIGMA      float32                        km / s Mean line-width of well-measured broad rest-frame optical emission lines.
              BROAD_SIGMARMS      float32                        km / s Root-mean-square scatter in BROAD_SIGMA.
                    UV_SIGMA      float32                        km / s Mean line-width of well-measured rest-frame UV emission lines.
                 UV_SIGMARMS      float32                        km / s Root-mean-square scatter in UV_SIGMA.
          MGII_DOUBLET_RATIO      float32                               MgII 2796 / 2803 doublet line-ratio.
           OII_DOUBLET_RATIO      float32                               [OII] 3726 / 3729 doublet line-ratio.
           SII_DOUBLET_RATIO      float32                               [SII] 6731 / 6716 doublet line-ratio.
                 LYALPHA_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
            LYALPHA_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                LYALPHA_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
           LYALPHA_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
             LYALPHA_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
        LYALPHA_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
              LYALPHA_VSHIFT      float32                        km / s Velocity shift relative to Z.
               LYALPHA_SIGMA      float32                        km / s Gaussian emission-line width.
                LYALPHA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
           LYALPHA_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                  LYALPHA_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
             LYALPHA_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
          LYALPHA_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
            LYALPHA_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                LYALPHA_CHI2      float32                               Reduced chi^2 of the line-fit.
                LYALPHA_NPIX        int32                               Number of pixels attributed to the emission line.
                 OI_1304_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
            OI_1304_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                OI_1304_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
           OI_1304_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
             OI_1304_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
        OI_1304_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
              OI_1304_VSHIFT      float32                        km / s Velocity shift relative to Z.
               OI_1304_SIGMA      float32                        km / s Gaussian emission-line width.
                OI_1304_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
           OI_1304_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                  OI_1304_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
             OI_1304_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
          OI_1304_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
            OI_1304_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                OI_1304_CHI2      float32                               Reduced chi^2 of the line-fit.
                OI_1304_NPIX        int32                               Number of pixels attributed to the emission line.
              SILIV_1396_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
         SILIV_1396_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
             SILIV_1396_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
        SILIV_1396_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
          SILIV_1396_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
     SILIV_1396_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
           SILIV_1396_VSHIFT      float32                        km / s Velocity shift relative to Z.
            SILIV_1396_SIGMA      float32                        km / s Gaussian emission-line width.
             SILIV_1396_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
        SILIV_1396_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
               SILIV_1396_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
          SILIV_1396_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
       SILIV_1396_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
         SILIV_1396_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
             SILIV_1396_CHI2      float32                               Reduced chi^2 of the line-fit.
             SILIV_1396_NPIX        int32                               Number of pixels attributed to the emission line.
                CIV_1549_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           CIV_1549_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               CIV_1549_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          CIV_1549_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            CIV_1549_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       CIV_1549_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             CIV_1549_VSHIFT      float32                        km / s Velocity shift relative to Z.
              CIV_1549_SIGMA      float32                        km / s Gaussian emission-line width.
               CIV_1549_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          CIV_1549_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 CIV_1549_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            CIV_1549_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         CIV_1549_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           CIV_1549_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               CIV_1549_CHI2      float32                               Reduced chi^2 of the line-fit.
               CIV_1549_NPIX        int32                               Number of pixels attributed to the emission line.
             SILIII_1892_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        SILIII_1892_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            SILIII_1892_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       SILIII_1892_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         SILIII_1892_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    SILIII_1892_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
          SILIII_1892_VSHIFT      float32                        km / s Velocity shift relative to Z.
           SILIII_1892_SIGMA      float32                        km / s Gaussian emission-line width.
            SILIII_1892_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       SILIII_1892_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              SILIII_1892_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
         SILIII_1892_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
      SILIII_1892_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        SILIII_1892_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            SILIII_1892_CHI2      float32                               Reduced chi^2 of the line-fit.
            SILIII_1892_NPIX        int32                               Number of pixels attributed to the emission line.
               CIII_1908_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          CIII_1908_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              CIII_1908_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         CIII_1908_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           CIII_1908_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      CIII_1908_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            CIII_1908_VSHIFT      float32                        km / s Velocity shift relative to Z.
             CIII_1908_SIGMA      float32                        km / s Gaussian emission-line width.
              CIII_1908_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         CIII_1908_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                CIII_1908_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           CIII_1908_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        CIII_1908_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          CIII_1908_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              CIII_1908_CHI2      float32                               Reduced chi^2 of the line-fit.
              CIII_1908_NPIX        int32                               Number of pixels attributed to the emission line.
               MGII_2796_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          MGII_2796_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              MGII_2796_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         MGII_2796_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           MGII_2796_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
            MGII_2796_VSHIFT      float32                        km / s Velocity shift relative to Z.
             MGII_2796_SIGMA      float32                        km / s Gaussian emission-line width.
              MGII_2796_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         MGII_2796_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                MGII_2796_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           MGII_2796_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        MGII_2796_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          MGII_2796_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              MGII_2796_CHI2      float32                               Reduced chi^2 of the line-fit.
              MGII_2796_NPIX        int32                               Number of pixels attributed to the emission line.
               MGII_2803_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          MGII_2803_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              MGII_2803_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         MGII_2803_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           MGII_2803_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      MGII_2803_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            MGII_2803_VSHIFT      float32                        km / s Velocity shift relative to Z.
             MGII_2803_SIGMA      float32                        km / s Gaussian emission-line width.
              MGII_2803_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         MGII_2803_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                MGII_2803_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           MGII_2803_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        MGII_2803_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          MGII_2803_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              MGII_2803_CHI2      float32                               Reduced chi^2 of the line-fit.
              MGII_2803_NPIX        int32                               Number of pixels attributed to the emission line.
                NEV_3346_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           NEV_3346_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               NEV_3346_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          NEV_3346_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            NEV_3346_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       NEV_3346_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             NEV_3346_VSHIFT      float32                        km / s Velocity shift relative to Z.
              NEV_3346_SIGMA      float32                        km / s Gaussian emission-line width.
               NEV_3346_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          NEV_3346_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 NEV_3346_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            NEV_3346_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         NEV_3346_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           NEV_3346_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               NEV_3346_CHI2      float32                               Reduced chi^2 of the line-fit.
               NEV_3346_NPIX        int32                               Number of pixels attributed to the emission line.
                NEV_3426_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           NEV_3426_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               NEV_3426_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          NEV_3426_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            NEV_3426_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       NEV_3426_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             NEV_3426_VSHIFT      float32                        km / s Velocity shift relative to Z.
              NEV_3426_SIGMA      float32                        km / s Gaussian emission-line width.
               NEV_3426_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          NEV_3426_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 NEV_3426_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            NEV_3426_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         NEV_3426_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           NEV_3426_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               NEV_3426_CHI2      float32                               Reduced chi^2 of the line-fit.
               NEV_3426_NPIX        int32                               Number of pixels attributed to the emission line.
                OII_3726_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           OII_3726_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               OII_3726_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          OII_3726_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            OII_3726_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       OII_3726_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             OII_3726_VSHIFT      float32                        km / s Velocity shift relative to Z.
              OII_3726_SIGMA      float32                        km / s Gaussian emission-line width.
               OII_3726_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          OII_3726_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 OII_3726_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            OII_3726_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         OII_3726_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           OII_3726_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               OII_3726_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
               OII_3726_NPIX        int32                               Number of pixels attributed to the emission line.
                OII_3729_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           OII_3729_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               OII_3729_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          OII_3729_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            OII_3729_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       OII_3729_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             OII_3729_VSHIFT      float32                        km / s Velocity shift relative to Z.
              OII_3729_SIGMA      float32                        km / s Gaussian emission-line width.
               OII_3729_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          OII_3729_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 OII_3729_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            OII_3729_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         OII_3729_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           OII_3729_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               OII_3729_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
               OII_3729_NPIX        int32                               Number of pixels attributed to the emission line.
              NEIII_3869_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
         NEIII_3869_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
             NEIII_3869_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
        NEIII_3869_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
          NEIII_3869_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
     NEIII_3869_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
           NEIII_3869_VSHIFT      float32                        km / s Velocity shift relative to Z.
            NEIII_3869_SIGMA      float32                        km / s Gaussian emission-line width.
             NEIII_3869_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
        NEIII_3869_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
               NEIII_3869_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
          NEIII_3869_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
       NEIII_3869_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
         NEIII_3869_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
             NEIII_3869_CHI2      float32                               Reduced chi^2 of the line-fit.
             NEIII_3869_NPIX        int32                               Number of pixels attributed to the emission line.
                      H6_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
                 H6_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                     H6_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
                H6_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
                  H6_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
             H6_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
                   H6_VSHIFT      float32                        km / s Velocity shift relative to Z.
                    H6_SIGMA      float32                        km / s Gaussian emission-line width.
                     H6_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
                H6_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                       H6_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
                  H6_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
               H6_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
                 H6_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                     H6_CHI2      float32                               Reduced chi^2 of the line-fit.
                     H6_NPIX        int32                               Number of pixels attributed to the emission line.
                H6_BROAD_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           H6_BROAD_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               H6_BROAD_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          H6_BROAD_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            H6_BROAD_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
             H6_BROAD_VSHIFT      float32                        km / s Velocity shift relative to Z.
              H6_BROAD_SIGMA      float32                        km / s Gaussian emission-line width.
               H6_BROAD_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          H6_BROAD_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 H6_BROAD_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            H6_BROAD_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         H6_BROAD_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           H6_BROAD_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               H6_BROAD_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
               H6_BROAD_NPIX        int32                               Number of pixels attributed to the emission line.
                HEPSILON_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           HEPSILON_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               HEPSILON_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          HEPSILON_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            HEPSILON_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       HEPSILON_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             HEPSILON_VSHIFT      float32                        km / s Velocity shift relative to Z.
              HEPSILON_SIGMA      float32                        km / s Gaussian emission-line width.
               HEPSILON_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          HEPSILON_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 HEPSILON_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            HEPSILON_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         HEPSILON_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           HEPSILON_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               HEPSILON_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
               HEPSILON_NPIX        int32                               Number of pixels attributed to the emission line.
          HEPSILON_BROAD_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
     HEPSILON_BROAD_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
         HEPSILON_BROAD_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
    HEPSILON_BROAD_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
      HEPSILON_BROAD_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
 HEPSILON_BROAD_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
       HEPSILON_BROAD_VSHIFT      float32                        km / s Velocity shift relative to Z.
        HEPSILON_BROAD_SIGMA      float32                        km / s Gaussian emission-line width.
         HEPSILON_BROAD_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
    HEPSILON_BROAD_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
           HEPSILON_BROAD_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
      HEPSILON_BROAD_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
   HEPSILON_BROAD_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
     HEPSILON_BROAD_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
         HEPSILON_BROAD_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
         HEPSILON_BROAD_NPIX        int32                               Number of pixels attributed to the emission line.
                  HDELTA_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
             HDELTA_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                 HDELTA_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
            HDELTA_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
              HDELTA_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         HDELTA_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
               HDELTA_VSHIFT      float32                        km / s Velocity shift relative to Z.
                HDELTA_SIGMA      float32                        km / s Gaussian emission-line width.
                 HDELTA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
            HDELTA_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                   HDELTA_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
              HDELTA_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
           HDELTA_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
             HDELTA_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                 HDELTA_CHI2      float32                               Reduced chi^2 of the line-fit.
                 HDELTA_NPIX        int32                               Number of pixels attributed to the emission line.
            HDELTA_BROAD_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       HDELTA_BROAD_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           HDELTA_BROAD_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      HDELTA_BROAD_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        HDELTA_BROAD_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
   HDELTA_BROAD_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
         HDELTA_BROAD_VSHIFT      float32                        km / s Velocity shift relative to Z.
          HDELTA_BROAD_SIGMA      float32                        km / s Gaussian emission-line width.
           HDELTA_BROAD_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      HDELTA_BROAD_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             HDELTA_BROAD_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
        HDELTA_BROAD_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
     HDELTA_BROAD_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       HDELTA_BROAD_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           HDELTA_BROAD_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
           HDELTA_BROAD_NPIX        int32                               Number of pixels attributed to the emission line.
                  HGAMMA_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
             HGAMMA_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                 HGAMMA_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
            HGAMMA_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
              HGAMMA_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         HGAMMA_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
               HGAMMA_VSHIFT      float32                        km / s Velocity shift relative to Z.
                HGAMMA_SIGMA      float32                        km / s Gaussian emission-line width.
                 HGAMMA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
            HGAMMA_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                   HGAMMA_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
              HGAMMA_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
           HGAMMA_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
             HGAMMA_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                 HGAMMA_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
                 HGAMMA_NPIX        int32                               Number of pixels attributed to the emission line.
            HGAMMA_BROAD_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       HGAMMA_BROAD_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           HGAMMA_BROAD_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      HGAMMA_BROAD_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        HGAMMA_BROAD_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
   HGAMMA_BROAD_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
         HGAMMA_BROAD_VSHIFT      float32                        km / s Velocity shift relative to Z.
          HGAMMA_BROAD_SIGMA      float32                        km / s Gaussian emission-line width.
           HGAMMA_BROAD_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      HGAMMA_BROAD_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             HGAMMA_BROAD_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
        HGAMMA_BROAD_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
     HGAMMA_BROAD_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       HGAMMA_BROAD_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           HGAMMA_BROAD_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
           HGAMMA_BROAD_NPIX        int32                               Number of pixels attributed to the emission line.
               OIII_4363_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          OIII_4363_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              OIII_4363_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         OIII_4363_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           OIII_4363_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      OIII_4363_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            OIII_4363_VSHIFT      float32                        km / s Velocity shift relative to Z.
             OIII_4363_SIGMA      float32                        km / s Gaussian emission-line width.
              OIII_4363_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         OIII_4363_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                OIII_4363_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           OIII_4363_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        OIII_4363_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          OIII_4363_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              OIII_4363_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
              OIII_4363_NPIX        int32                               Number of pixels attributed to the emission line.
                HEI_4471_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           HEI_4471_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               HEI_4471_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          HEI_4471_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            HEI_4471_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       HEI_4471_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             HEI_4471_VSHIFT      float32                        km / s Velocity shift relative to Z.
              HEI_4471_SIGMA      float32                        km / s Gaussian emission-line width.
               HEI_4471_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          HEI_4471_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 HEI_4471_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            HEI_4471_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         HEI_4471_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           HEI_4471_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               HEI_4471_CHI2      float32                               Reduced chi^2 of the line-fit.
               HEI_4471_NPIX        int32                               Number of pixels attributed to the emission line.
               HEII_4686_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HEII_4686_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              HEII_4686_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HEII_4686_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           HEII_4686_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      HEII_4686_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            HEII_4686_VSHIFT      float32                        km / s Velocity shift relative to Z.
             HEII_4686_SIGMA      float32                        km / s Gaussian emission-line width.
              HEII_4686_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HEII_4686_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                HEII_4686_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           HEII_4686_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        HEII_4686_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HEII_4686_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HEII_4686_CHI2      float32                               Reduced chi^2 of the line-fit.
              HEII_4686_NPIX        int32                               Number of pixels attributed to the emission line.
                   HBETA_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
              HBETA_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                  HBETA_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
             HBETA_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
               HBETA_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          HBETA_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
                HBETA_VSHIFT      float32                        km / s Velocity shift relative to Z.
                 HBETA_SIGMA      float32                        km / s Gaussian emission-line width.
                  HBETA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
             HBETA_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                    HBETA_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
               HBETA_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
            HBETA_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
              HBETA_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                  HBETA_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
                  HBETA_NPIX        int32                               Number of pixels attributed to the emission line.
             HBETA_BROAD_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        HBETA_BROAD_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            HBETA_BROAD_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       HBETA_BROAD_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         HBETA_BROAD_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    HBETA_BROAD_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
          HBETA_BROAD_VSHIFT      float32                        km / s Velocity shift relative to Z.
           HBETA_BROAD_SIGMA      float32                        km / s Gaussian emission-line width.
            HBETA_BROAD_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       HBETA_BROAD_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              HBETA_BROAD_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
         HBETA_BROAD_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
      HBETA_BROAD_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        HBETA_BROAD_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            HBETA_BROAD_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
            HBETA_BROAD_NPIX        int32                               Number of pixels attributed to the emission line.
               OIII_4959_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          OIII_4959_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              OIII_4959_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         OIII_4959_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           OIII_4959_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      OIII_4959_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            OIII_4959_VSHIFT      float32                        km / s Velocity shift relative to Z.
             OIII_4959_SIGMA      float32                        km / s Gaussian emission-line width.
              OIII_4959_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         OIII_4959_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                OIII_4959_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           OIII_4959_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        OIII_4959_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          OIII_4959_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              OIII_4959_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
              OIII_4959_NPIX        int32                               Number of pixels attributed to the emission line.
               OIII_5007_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          OIII_5007_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              OIII_5007_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         OIII_5007_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           OIII_5007_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      OIII_5007_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            OIII_5007_VSHIFT      float32                        km / s Velocity shift relative to Z.
             OIII_5007_SIGMA      float32                        km / s Gaussian emission-line width.
              OIII_5007_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         OIII_5007_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                OIII_5007_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           OIII_5007_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        OIII_5007_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          OIII_5007_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              OIII_5007_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
              OIII_5007_NPIX        int32                               Number of pixels attributed to the emission line.
                NII_5755_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           NII_5755_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               NII_5755_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          NII_5755_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            NII_5755_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       NII_5755_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             NII_5755_VSHIFT      float32                        km / s Velocity shift relative to Z.
              NII_5755_SIGMA      float32                        km / s Gaussian emission-line width.
               NII_5755_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          NII_5755_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 NII_5755_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            NII_5755_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         NII_5755_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           NII_5755_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               NII_5755_CHI2      float32                               Reduced chi^2 of the line-fit.
               NII_5755_NPIX        int32                               Number of pixels attributed to the emission line.
                HEI_5876_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           HEI_5876_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               HEI_5876_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          HEI_5876_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            HEI_5876_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       HEI_5876_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             HEI_5876_VSHIFT      float32                        km / s Velocity shift relative to Z.
              HEI_5876_SIGMA      float32                        km / s Gaussian emission-line width.
               HEI_5876_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          HEI_5876_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 HEI_5876_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            HEI_5876_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         HEI_5876_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           HEI_5876_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               HEI_5876_CHI2      float32                               Reduced chi^2 of the line-fit.
               HEI_5876_NPIX        int32                               Number of pixels attributed to the emission line.
                 OI_6300_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
            OI_6300_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                OI_6300_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
           OI_6300_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
             OI_6300_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
        OI_6300_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
              OI_6300_VSHIFT      float32                        km / s Velocity shift relative to Z.
               OI_6300_SIGMA      float32                        km / s Gaussian emission-line width.
                OI_6300_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
           OI_6300_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                  OI_6300_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
             OI_6300_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
          OI_6300_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
            OI_6300_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                OI_6300_CHI2      float32                               Reduced chi^2 of the line-fit.
                OI_6300_NPIX        int32                               Number of pixels attributed to the emission line.
                NII_6548_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           NII_6548_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               NII_6548_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          NII_6548_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            NII_6548_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       NII_6548_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             NII_6548_VSHIFT      float32                        km / s Velocity shift relative to Z.
              NII_6548_SIGMA      float32                        km / s Gaussian emission-line width.
               NII_6548_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          NII_6548_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 NII_6548_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            NII_6548_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         NII_6548_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           NII_6548_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               NII_6548_CHI2      float32                               Reduced chi^2 of the line-fit.
               NII_6548_NPIX        int32                               Number of pixels attributed to the emission line.
                  HALPHA_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
             HALPHA_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                 HALPHA_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
            HALPHA_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
              HALPHA_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         HALPHA_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
               HALPHA_VSHIFT      float32                        km / s Velocity shift relative to Z.
                HALPHA_SIGMA      float32                        km / s Gaussian emission-line width.
                 HALPHA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
            HALPHA_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                   HALPHA_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
              HALPHA_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
           HALPHA_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
             HALPHA_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                 HALPHA_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
                 HALPHA_NPIX        int32                               Number of pixels attributed to the emission line.
            HALPHA_BROAD_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       HALPHA_BROAD_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           HALPHA_BROAD_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      HALPHA_BROAD_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        HALPHA_BROAD_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
   HALPHA_BROAD_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
         HALPHA_BROAD_VSHIFT      float32                        km / s Velocity shift relative to Z.
          HALPHA_BROAD_SIGMA      float32                        km / s Gaussian emission-line width.
           HALPHA_BROAD_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      HALPHA_BROAD_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             HALPHA_BROAD_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
        HALPHA_BROAD_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
     HALPHA_BROAD_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       HALPHA_BROAD_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           HALPHA_BROAD_CHI2      float32                               Reduced chi^2 of the line-fit (default value 1e6).
           HALPHA_BROAD_NPIX        int32                               Number of pixels attributed to the emission line.
                NII_6584_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           NII_6584_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               NII_6584_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          NII_6584_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            NII_6584_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       NII_6584_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             NII_6584_VSHIFT      float32                        km / s Velocity shift relative to Z.
              NII_6584_SIGMA      float32                        km / s Gaussian emission-line width.
               NII_6584_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          NII_6584_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 NII_6584_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            NII_6584_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         NII_6584_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           NII_6584_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               NII_6584_CHI2      float32                               Reduced chi^2 of the line-fit.
               NII_6584_NPIX        int32                               Number of pixels attributed to the emission line.
                SII_6716_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           SII_6716_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               SII_6716_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          SII_6716_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            SII_6716_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       SII_6716_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             SII_6716_VSHIFT      float32                        km / s Velocity shift relative to Z.
              SII_6716_SIGMA      float32                        km / s Gaussian emission-line width.
               SII_6716_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          SII_6716_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 SII_6716_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            SII_6716_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         SII_6716_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           SII_6716_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               SII_6716_CHI2      float32                               Reduced chi^2 of the line-fit.
               SII_6716_NPIX        int32                               Number of pixels attributed to the emission line.
                SII_6731_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           SII_6731_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               SII_6731_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          SII_6731_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            SII_6731_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       SII_6731_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             SII_6731_VSHIFT      float32                        km / s Velocity shift relative to Z.
              SII_6731_SIGMA      float32                        km / s Gaussian emission-line width.
               SII_6731_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          SII_6731_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 SII_6731_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            SII_6731_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         SII_6731_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           SII_6731_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               SII_6731_CHI2      float32                               Reduced chi^2 of the line-fit.
               SII_6731_NPIX        int32                               Number of pixels attributed to the emission line.
                OII_7320_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           OII_7320_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               OII_7320_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          OII_7320_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            OII_7320_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       OII_7320_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             OII_7320_VSHIFT      float32                        km / s Velocity shift relative to Z.
              OII_7320_SIGMA      float32                        km / s Gaussian emission-line width.
               OII_7320_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          OII_7320_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 OII_7320_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            OII_7320_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         OII_7320_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           OII_7320_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               OII_7320_CHI2      float32                               Reduced chi^2 of the line-fit.
               OII_7320_NPIX        int32                               Number of pixels attributed to the emission line.
                OII_7330_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           OII_7330_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               OII_7330_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          OII_7330_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            OII_7330_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       OII_7330_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
             OII_7330_VSHIFT      float32                        km / s Velocity shift relative to Z.
              OII_7330_SIGMA      float32                        km / s Gaussian emission-line width.
               OII_7330_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          OII_7330_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 OII_7330_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
            OII_7330_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
         OII_7330_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           OII_7330_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               OII_7330_CHI2      float32                               Reduced chi^2 of the line-fit.
               OII_7330_NPIX        int32                               Number of pixels attributed to the emission line.
               SIII_9069_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          SIII_9069_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              SIII_9069_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         SIII_9069_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           SIII_9069_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      SIII_9069_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            SIII_9069_VSHIFT      float32                        km / s Velocity shift relative to Z.
             SIII_9069_SIGMA      float32                        km / s Gaussian emission-line width.
              SIII_9069_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         SIII_9069_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                SIII_9069_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           SIII_9069_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        SIII_9069_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          SIII_9069_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              SIII_9069_CHI2      float32                               Reduced chi^2 of the line-fit.
              SIII_9069_NPIX        int32                               Number of pixels attributed to the emission line.
               SIII_9532_AMP      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          SIII_9532_AMP_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              SIII_9532_FLUX      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         SIII_9532_FLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           SIII_9532_BOXFLUX      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      SIII_9532_BOXFLUX_IVAR      float32           1e+34 cm4 s2 / erg2 Inverse variance of boxcar-integrated flux.
            SIII_9532_VSHIFT      float32                        km / s Velocity shift relative to Z.
             SIII_9532_SIGMA      float32                        km / s Gaussian emission-line width.
              SIII_9532_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         SIII_9532_CONT_IVAR      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                SIII_9532_EW      float32                      Angstrom Rest-frame emission-line equivalent width.
           SIII_9532_EW_IVAR      float32                 1 / Angstrom2 Inverse variance of equivalent width.
        SIII_9532_FLUX_LIMIT      float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          SIII_9532_EW_LIMIT      float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              SIII_9532_CHI2      float32                               Reduced chi^2 of the line-fit.
              SIII_9532_NPIX        int32                               Number of pixels attributed to the emission line.
============================ ============ ============================= ============================================

HDU02
-----

EXTNAME = METADATA

Metadata associated with each objected fitted.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

.. collapse:: Required Header Keywords Table
    :open:

    .. rst-class:: keywords

    ======== ================ ==== ==============================================
    KEY      Example Value    Type Comment
    ======== ================ ==== ==============================================
    NAXIS1   339              int  length of dimension 1
    NAXIS2   338              int  length of dimension 2
    CHECKSUM hFY6jCV3hCV3hCV3 str  HDU checksum updated 2022-08-03T14:25:08
    DATASUM  1759692941       str  data unit checksum updated 2022-08-03T14:25:08
    ======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: columns

====================== =========== ========== ==========================================
Name                   Type        Units      Description
====================== =========== ========== ==========================================
              TARGETID   int64                Unique target ID.
           SURVEY [1]_  bytes3                Survey name.
          PROGRAM [1]_  bytes6                Program name.
          HEALPIX [1]_   int32                Healpixel number.
           TILEID [2]_   int32                Tile ID number.
            FIBER [2]_   int32                Fiber number.
            NIGHT [2]_   int32                Night of observation.
            EXPID [3]_   int32                Exposure ID number.
           TILEID_LIST    str5                List of tile IDs that went into healpix coadd.
                    RA float64            deg Right ascension from target catalog.
                   DEC float64            deg Declination from target catalog.
     COADD_FIBERSTATUS   int64                Fiber status bit.
       CMX_TARGET [4]_   int64                Commissioning (CMX) targeting bit.
           DESI_TARGET   int64                DESI targeting bit.
            BGS_TARGET   int64                BGS targeting bit.
            MWS_TARGET   int64                MWS targeting bit.
           SCND_TARGET   int64                Secondary target targeting bit.
  SV1_DESI_TARGET [4]_   int64                SV1 DESI targeting bit.
   SV1_BGS_TARGET [4]_   int64                SV1 BGS targeting bit.
   SV1_MWS_TARGET [4]_   int64                SV1 MWS targeting bit.
  SV2_DESI_TARGET [4]_   int64                SV2 DESI targeting bit.
   SV2_BGS_TARGET [4]_   int64                SV2 BGS targeting bit.
   SV2_MWS_TARGET [4]_   int64                SV2 MWS targeting bit.
  SV3_DESI_TARGET [4]_   int64                SV3 DESI targeting bit.
   SV3_BGS_TARGET [4]_   int64                SV3 BGS targeting bit.
   SV3_MWS_TARGET [4]_   int64                SV3 MWS targeting bit.
  SV1_SCND_TARGET [4]_   int64                SV1 secondary targeting bit.
  SV2_SCND_TARGET [4]_   int64                SV2 secondary targeting bit.
  SV3_SCND_TARGET [4]_   int64                SV3 secondary targeting bit.
                     Z float64                Redshift based on Redrock or QuasarNet (for QSO targets only).
                 ZWARN    int8                Redrock zwarning bit.
             DELTACHI2 float64                Redrock delta-chi-squared.
              SPECTYPE    str6                Redrock spectral classification.
                  Z_RR float64                Redrock redshift.
             TSNR2_BGS float32                Template signal-to-noise ratio squared for BGS targets.
             TSNR2_LRG float32                Template signal-to-noise ratio squared for LRG targets.
             TSNR2_ELG float32                Template signal-to-noise ratio squared for ELG targets.
             TSNR2_QSO float32                Template signal-to-noise ratio squared for QSO targets.
             TSNR2_LYA float32                Template signal-to-noise ratio squared for LYA targets.
               PHOTSYS    str1                Photometric system (*N* or *S*).
               LS_ID     int64                Unique Legacy Surveys identification number.
           FIBERFLUX_G float32           nmgy Fiber g-band flux corrected for Galactic extinction.
           FIBERFLUX_R float32           nmgy Fiber r-band flux corrected for Galactic extinction.
           FIBERFLUX_Z float32           nmgy Fiber z-band flux corrected for Galactic extinction.
        FIBERTOTFLUX_G float32           nmgy Fibertot g-band flux corrected for Galactic extinction.
        FIBERTOTFLUX_R float32           nmgy Fibertot r-band flux corrected for Galactic extinction.
        FIBERTOTFLUX_Z float32           nmgy Fibertot z-band flux corrected for Galactic extinction.
                FLUX_G float32           nmgy Total g-band flux corrected for Galactic extinction.
                FLUX_R float32           nmgy Total r-band flux corrected for Galactic extinction.
                FLUX_Z float32           nmgy Total z-band flux corrected for Galactic extinction.
               FLUX_W1 float32           nmgy Total W1-band flux corrected for Galactic extinction.
               FLUX_W2 float32           nmgy Total W2-band flux corrected for Galactic extinction.
               FLUX_W3 float32           nmgy Total W3-band flux corrected for Galactic extinction.
               FLUX_W4 float32           nmgy Total W4-band flux corrected for Galactic extinction.
           FLUX_IVAR_G float32      1 / nmgy2 Inverse variance of FLUX_G corrected for Galactic extinction.
           FLUX_IVAR_R float32      1 / nmgy2 Inverse variance of FLUX_R corrected for Galactic extinction.
           FLUX_IVAR_Z float32      1 / nmgy2 Inverse variance of FLUX_Z corrected for Galactic extinction.
          FLUX_IVAR_W1 float32      1 / nmgy2 Inverse variance of FLUX_W1 corrected for Galactic extinction.
          FLUX_IVAR_W2 float32      1 / nmgy2 Inverse variance of FLUX_W2 corrected for Galactic extinction.
          FLUX_IVAR_W3 float32      1 / nmgy2 Inverse variance of FLUX_W3 corrected for Galactic extinction.
          FLUX_IVAR_W4 float32      1 / nmgy2 Inverse variance of FLUX_W4 corrected for Galactic extinction.
                   EBV float32            mag Milky Way foreground dust reddening.
     MW_TRANSMISSION_G float32                Milky Way foreground dust transmission factor [0-1] in the g-band.
     MW_TRANSMISSION_R float32                Milky Way foreground dust transmission factor [0-1] in the r-band.
     MW_TRANSMISSION_Z float32                Milky Way foreground dust transmission factor [0-1] in the z-band.
    MW_TRANSMISSION_W1 float32                Milky Way foreground dust transmission factor [0-1] in the W1-band.
    MW_TRANSMISSION_W2 float32                Milky Way foreground dust transmission factor [0-1] in the W2-band.
    MW_TRANSMISSION_W3 float32                Milky Way foreground dust transmission factor [0-1] in the W3-band.
    MW_TRANSMISSION_W4 float32                Milky Way foreground dust transmission factor [0-1] in the W4-band.
====================== =========== ========== ==========================================

HDU03
-----

EXTNAME = MODELS

Best-fitting model spectra (corrected for Galactic extinction).

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

.. collapse:: Required Header Keywords Table
    :open:

    .. rst-class:: keywords

    ======== ============================ ===== ==============================================
    KEY      Example Value Type Comment
    ======== ============================ ===== ==============================================
    NAXIS1   7781                         int   number of pixels
    NAXIS2   3                            int   number of models
    NAXIS3   338                          int   number of objects
    BUNIT    10**-17 erg/(s cm2 Angstrom) str   flux unit
    CUNIT1   Angstrom                     str   wavelength unit
    CTYPE1   WAVE                         str   type of axis
    CRVAL1   3600.0                       float wavelength of pixel CRPIX1 (Angstrom)
    CRPIX1   0                            int   0-indexed pixel number corresponding to CRVAL1
    CDELT1   0.8                          float pixel size (Angstrom)
    DC-FLAG  0                            int   0 = linear wavelength vector
    AIRORVAC vac                          str   wavelengths in vacuum (vac)
    CHECKSUM Q9eTT8bTQ8bTQ8bT             str   HDU checksum updated 2022-08-03T14:25:08
    DATASUM  3074907107                   str   data unit checksum updated 2022-08-03T14:25:08
    ======== ============================ ===== ==============================================

Data: FITS image [int32, 7781x3,338]

.. [1] Column only present when fitting healpix coadds.
       
.. [2] Column only present when fitting cumulative, per-night, or per-expopsure tile-based coadds.
       
.. [3] Column only present when fitting per-exposure tile-based coadds.

.. [4] Column only present in Commissioning and Survey Validation spectroscopic observations.
