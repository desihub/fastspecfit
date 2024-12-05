.. _fastspec datamodel:

===================
Fastspec Data Model
===================

:Summary: Spectrophotometric fitting results.
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
HDU01_ METADATA     BINTABLE Table with sample metadata.
HDU02_ SPECPHOT     BINTABLE Table with spectrophotometric quantities.
HDU03_ FASTSPEC     BINTABLE Table with spectral fitting results.
HDU04_ MODELS       IMAGE    Model continuum, smooth continuum, and emission-line spectra.
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

    ======== ================ ===== ==============================================
    KEY      Example Value    Type  Comment
    ======== ================ ===== ==============================================
    SPECPROD iron             str   Spectroscopic production name.
    COADDTYP healpix          str   Type of spectral coadd (healpix|cumulative|pernight|perexp).
    INPUTZ   F                bool  Input redshifts provided.
    INPUTS   F                bool  Input seeds provided.
    CONSAGE  F                bool  Stellar population synthesis ages constrained by redshift.
    USEQNET  T                bool  QuasarNet redshifts used where applicable.
    NMONTE   50               int   Number of Monte Carlo realizations.
    SEED     1                int   Random seed used for Monte Carlo realizations.
    NOSCORR  F                bool  No smooth continuum correction derived.
    NOPHOTO  F                bool  No fitting to photometry carried out.
    BRDLFIT  T                bool  Broad Balmer emission-line modeling attempted.
    UFLOOR   0.01             float Spectroscopic fractional uncertainty floor.
    SNRBBALM 2.5              float Minimum broad Balmer signal-to-noise ratio.
    CHECKSUM ZHWHZHVFZHVFZHVF str   HDU checksum updated 2022-08-02T19:06:22
    DATASUM  0                str   data unit checksum updated 2022-08-02T19:06:22
    ======== ================ ===== ==============================================

Empty HDU.

HDU01
-----

EXTNAME = METADATA

Metadata associated with each object fitted.

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
           TILEID_LIST    str?                List of tile IDs that went into healpix coadd.
                    RA float64            deg Right ascension from target catalog.
                   DEC float64            deg Declination from target catalog.
     COADD_FIBERSTATUS   int64                Fiber status bit.
       CMX_TARGET [5]_   int64                Commissioning (CMX) targeting bit.
           DESI_TARGET   int64                DESI targeting bit.
            BGS_TARGET   int64                BGS targeting bit.
            MWS_TARGET   int64                MWS targeting bit.
           SCND_TARGET   int64                Secondary target targeting bit.
  SV1_DESI_TARGET [5]_   int64                SV1 DESI targeting bit.
   SV1_BGS_TARGET [5]_   int64                SV1 BGS targeting bit.
   SV1_MWS_TARGET [5]_   int64                SV1 MWS targeting bit.
  SV2_DESI_TARGET [5]_   int64                SV2 DESI targeting bit.
   SV2_BGS_TARGET [5]_   int64                SV2 BGS targeting bit.
   SV2_MWS_TARGET [5]_   int64                SV2 MWS targeting bit.
  SV3_DESI_TARGET [5]_   int64                SV3 DESI targeting bit.
   SV3_BGS_TARGET [5]_   int64                SV3 BGS targeting bit.
   SV3_MWS_TARGET [5]_   int64                SV3 MWS targeting bit.
  SV1_SCND_TARGET [5]_   int64                SV1 secondary targeting bit.
  SV2_SCND_TARGET [5]_   int64                SV2 secondary targeting bit.
  SV3_SCND_TARGET [5]_   int64                SV3 secondary targeting bit.
                     Z float64                Redshift based on Redrock or QuasarNet (for QSO targets only).
                 ZWARN    int8                Redrock zwarning bit.
             DELTACHI2 float64                Redrock delta-chi-squared.
              SPECTYPE    str6                Redrock spectral classification type.
               SUBTYPE   str20                Redrock spectral subtype.
                  Z_RR float64                Redrock redshift.
             TSNR2_BGS float32                Template signal-to-noise ratio squared for BGS targets.
             TSNR2_LRG float32                Like TSNR2_BGS but for LRG targets.
             TSNR2_ELG float32                Like TSNR2_BGS but for ELG targets.
             TSNR2_QSO float32                Like TSNR2_BGS but for QSO targets.
             TSNR2_LYA float32                Like TSNR2_BGS but for LYA targets.
               PHOTSYS    str1                Photometric system (*N* or *S*).
                 LS_ID   int64                Unique Legacy Surveys identification number.
           FIBERFLUX_G float32           nmgy Fiber *g*-band flux corrected for Galactic extinction.
           FIBERFLUX_R float32           nmgy Like FIBERFLUX_G but for the *r*-band.
           FIBERFLUX_Z float32           nmgy Like FIBERFLUX_G but for the *z*-band.
        FIBERTOTFLUX_G float32           nmgy Fibertot *g*-band flux corrected for Galactic extinction.
        FIBERTOTFLUX_R float32           nmgy Like FIBERTOTFLUX_G but for the *r*-band.
        FIBERTOTFLUX_Z float32           nmgy Like FIBERTOTFLUX_G but for the *z*-band.
                FLUX_G float32           nmgy Total *g*-band flux corrected for Galactic extinction.
                FLUX_R float32           nmgy Like FLUX_G but for the *r*-band.
                FLUX_Z float32           nmgy Like FLUX_G but for the *z*-band.
               FLUX_W1 float32           nmgy Like FLUX_G but for the *W1*-band.
               FLUX_W2 float32           nmgy Like FLUX_G but for the *W2*-band.
               FLUX_W3 float32           nmgy Like FLUX_G but for the *W3*-band.
               FLUX_W4 float32           nmgy Like FLUX_G but for the *W4*-band.
           FLUX_IVAR_G float32      1 / nmgy2 Inverse variance of FLUX_G corrected for Galactic extinction.
           FLUX_IVAR_R float32      1 / nmgy2 Like FLUX_IVAR_G but for the *r*-band.
           FLUX_IVAR_Z float32      1 / nmgy2 Like FLUX_IVAR_G but for the *z*-band.
          FLUX_IVAR_W1 float32      1 / nmgy2 Like FLUX_IVAR_G but for the *W1*-band.
          FLUX_IVAR_W2 float32      1 / nmgy2 Like FLUX_IVAR_G but for the *W2*-band.
          FLUX_IVAR_W3 float32      1 / nmgy2 Like FLUX_IVAR_G but for the *W3*-band.
          FLUX_IVAR_W4 float32      1 / nmgy2 Like FLUX_IVAR_G but for the *W4*-band.
                   EBV float32            mag Milky Way foreground color excess.
     MW_TRANSMISSION_G float32                Milky Way foreground dust transmission factor [0-1] in the *g*-band.
     MW_TRANSMISSION_R float32                Like MW_TRANSMISSION_G but for the *r*-band.
     MW_TRANSMISSION_Z float32                Like MW_TRANSMISSION_G but for the *z*-band.
    MW_TRANSMISSION_W1 float32                Like MW_TRANSMISSION_G but for the *W1*-band.
    MW_TRANSMISSION_W2 float32                Like MW_TRANSMISSION_G but for the *W2*-band.
    MW_TRANSMISSION_W3 float32                Like MW_TRANSMISSION_G but for the *W3*-band.
    MW_TRANSMISSION_W4 float32                Like MW_TRANSMISSION_G but for the *W4*-band.
====================== =========== ========== ==========================================

HDU02
-----

EXTNAME = SPECPHOT

Table with spectrophotometric quantities.

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

============================= ============ ============================= ============================================
Name                          Type         Units                         Description
============================= ============ ============================= ============================================
                     TARGETID        int64                               Unique target ID.
                  SURVEY [1]_       bytes3                               Survey name.
                 PROGRAM [1]_       bytes6                               Program name.
                 HEALPIX [1]_        int32                               Healpixel number.
                  TILEID [2]_        int32                               Tile ID number.
                   NIGHT [2]_        int32                               Night of observation.
                   FIBER [2]_        int32                               Fiber number.
                   EXPID [3]_        int32                               Exposure ID number.
                         SEED        int64                               Random seed used for Monte Carlo realizations.
                        COEFF   float32[5]  1e-17 erg / (Angstrom cm2 s) Stellar continuum coefficients.
                        RCHI2      float32                               Reduced chi-squared of the full-spectrophotometric fit.
                   RCHI2_LINE      float32                               Reduced chi-squared of the emission-line model fit.
                   RCHI2_CONT      float32                               Reduced chi-squared of the fit to the stellar continuum.
                   RCHI2_PHOT      float32                               Reduced chi-squared of the fit to the broadband photometry.
                        VDISP      float32                        km / s Stellar velocity dispersion.
                   VDISP_IVAR      float32                      s2 / km2 Inverse variance of VDISP.
                         TAUV      float32                               *V*-band dust optical depth of the integrated stellar population.
                    TAUV_IVAR      float32                               Inverse variance of TAUV.
                          AGE      float32                           Gyr Light-weighted age.
                     AGE_IVAR      float32                      1 / Gyr2 Inverse variance of AGE.
                        ZZSUN      float32                               Logarithmic stellar metallicity relative to solar.
                   ZZSUN_IVAR      float32                               Inverse variance of ZZSUN.
                     LOGMSTAR      float32                          Msun Logarithmic stellar mass (h=1.0, Chabrier+2003 initial mass function).
                LOGMSTAR_IVAR      float32                     1 / Msun2 Inverse variance of LOGMSTAR.
                          SFR      float32                     Msun / yr Instantaneous star formation rate (h=1.0, Chabrier+2003 initial mass function).
                     SFR_IVAR      float32                   yr2 / Msun2 Inverse variance of SFR.
                       DN4000      float32                               Narrow 4000-A break index (from Balogh et al. 1999) measured from the emission-line subtracted spectrum.
                   DN4000_OBS      float32                               Narrow 4000-A break index measured from the observed spectrum.
                  DN4000_IVAR      float32                               Inverse variance of DN4000_OBS and of DN4000.
                 DN4000_MODEL      float32                               Narrow 4000-A break index measured from the best-fitting continuum model.
            DN4000_MODEL_IVAR      float32                               Inverse variance of DN4000_MODEL.
                 FLUX_SYNTH_G      float32                          nmgy *g*-band flux (in the PHOTSYS photometric system) synthesized from the observed spectrum.
                 FLUX_SYNTH_R      float32                          nmgy Like FLUX_SYNTH_G but for the *r*-band.
                 FLUX_SYNTH_Z      float32                          nmgy Like FLUX_SYNTH_G but for the *z*-band.
       FLUX_SYNTH_SPECMODEL_G      float32                          nmgy *g*-band flux (in the PHOTSYS photometric system) synthesized from the best-fitting spectroscopic model.
       FLUX_SYNTH_SPECMODEL_R      float32                          nmgy Like FLUX_SYNTH_SPECMODEL_G but in the *r*-band.
       FLUX_SYNTH_SPECMODEL_Z      float32                          nmgy Like FLUX_SYNTH_SPECMODEL_G but in the *z*-band.
       FLUX_SYNTH_PHOTMODEL_G      float32                          nmgy *g*-band flux (in the PHOTSYS photometric system) synthesized from the best-fitting photometric continuum model.
       FLUX_SYNTH_PHOTMODEL_R      float32                          nmgy Like FLUX_SYNTH_PHOTMODEL_G but in the *r*-band.
       FLUX_SYNTH_PHOTMODEL_Z      float32                          nmgy Like FLUX_SYNTH_PHOTMODEL_G but in the *z*-band.
      FLUX_SYNTH_PHOTMODEL_W1      float32                          nmgy Like FLUX_SYNTH_PHOTMODEL_G but in the *W1*-band.
      FLUX_SYNTH_PHOTMODEL_W2      float32                          nmgy Like FLUX_SYNTH_PHOTMODEL_G but in the *W2*-band.
      FLUX_SYNTH_PHOTMODEL_W3      float32                          nmgy Like FLUX_SYNTH_PHOTMODEL_G but in the *W3*-band.
      FLUX_SYNTH_PHOTMODEL_W4      float32                          nmgy Like FLUX_SYNTH_PHOTMODEL_G but in the *W4*-band.
        ABSMAG10_DECAM_G [4]_      float32                           mag Absolute magnitude in DECam *g*-band band-shifted to z=1.0 assuming h=1.0.
        ABSMAG10_IVAR_DECAM_G      float32                      1 / mag2 Inverse variance corresponding to ABSMAG10_DECAM_G.
       ABSMAG10_SYNTH_DECAM_G      float32                           mag Synthesized absolute magnitude in DECam *g*-band band-shifted to z=1.0 assuming h=1.0.
  ABSMAG10_SYNTH_IVAR_DECAM_G      float32                      1 / mag2 Inverse variance corresponding to ABSMAG10_SYNTH_DECAM_G.
        ABSMAG10_DECAM_R [4]_      float32                           mag Like ABSMAG10_DECAM_G but for DECam *r*-band.
        ABSMAG10_IVAR_DECAM_R      float32                      1 / mag2 Like ABSMAG10_IVAR_DECAM_G but for DECam *r*-band.
       ABSMAG10_SYNTH_DECAM_R      float32                           mag Like ABSMAG10_SYNTH_DECAM_G but for DECam *r*-band.
  ABSMAG10_SYNTH_IVAR_DECAM_R      float32                      1 / mag2 Like ABSMAG10_SYNTH_IVAR_DECAM_G but for DECam *r*-band.
        ABSMAG10_DECAM_Z [4]_      float32                           mag Like ABSMAG10_DECAM_G but for DECam *z*-band.
        ABSMAG10_IVAR_DECAM_Z      float32                      1 / mag2 Like ABSMAG10_IVAR_DECAM_G but for DECam *z*-band.
       ABSMAG10_SYNTH_DECAM_Z      float32                           mag Like ABSMAG10_SYNTH_DECAM_G but for DECam *z*-band.
  ABSMAG10_SYNTH_IVAR_DECAM_Z      float32                      1 / mag2 Like ABSMAG10_SYNTH_IVAR_DECAM_G but for DECam *z*-band.
              ABSMAG00_U [4]_      float32                           mag Absolute magnitude in Johnson/Cousins *U*-band band-shifted to z=0.0 assuming h=1.0.
              ABSMAG00_IVAR_U      float32                      1 / mag2 Inverse variance corresponding to ABSMAG00_U.
             ABSMAG00_SYNTH_U      float32                           mag Synthesized absolute magnitude in Johnson/Cousins *U*-band band-shifted to z=0.0 assuming h=1.0.
        ABSMAG00_SYNTH_IVAR_U      float32                      1 / mag2 Inverse variance corresponding to ABSMAG00_SYNTH_U.
              ABSMAG00_B [4]_      float32                           mag Like ABSMAG00_U but for Johnson/Cousins *B*-band.
              ABSMAG00_IVAR_B      float32                      1 / mag2 Like ABSMAG00_IVAR_U but for Johnson/Cousins *B*-band.
             ABSMAG00_SYNTH_B      float32                           mag Like ABSMAG00_SYNTH_U but for Johnson/Cousins *B*-band.
        ABSMAG00_SYNTH_IVAR_B      float32                      1 / mag2 Like ABSMAG00_SYNTH_IVAR_U but for Johnson/Cousins *B*-band.
              ABSMAG00_V [4]_      float32                           mag Like ABSMAG00_U but for Johnson/Cousins *V*-band.
              ABSMAG00_IVAR_V      float32                      1 / mag2 Like ABSMAG00_IVAR_U but for Johnson/Cousins *V*-band.
             ABSMAG00_SYNTH_V      float32                           mag Like ABSMAG00_SYNTH_U but for Johnson/Cousins *V*-band.
        ABSMAG00_SYNTH_IVAR_V      float32                      1 / mag2 Like ABSMAG00_SYNTH_IVAR_U but for Johnson/Cousins *V*-band.
      ABSMAG00_TWOMASS_J [4]_      float32                           mag Absolute magnitude in 2MASS *J*-band band-shifted to z=0.0 assuming h=1.0.
      ABSMAG00_IVAR_TWOMASS_J      float32                      1 / mag2 Inverse variance corresponding to ABSMAG00_TWOMASS_J.
     ABSMAG00_SYNTH_TWOMASS_J      float32                           mag Synthesized absolute magnitude in 2MASS *J*-band band-shifted to z=0.0 assuming h=1.0.
ABSMAG00_SYNTH_IVAR_TWOMASS_J      float32                      1 / mag2 Inverse variance corresponding to ABSMAG01_SYNTH_TWOMASS_J.
         ABSMAG01_SDSS_U [4]_      float32                           mag Absolute magnitude in SDSS *u*-band band-shifted to z=0.1 assuming h=1.0.
         ABSMAG01_IVAR_SDSS_U      float32                      1 / mag2 Inverse variance corresponding to ABSMAG01_SDSS_U.
        ABSMAG01_SYNTH_SDSS_U      float32                           mag Synthesized absolute magnitude in SDSS *u*-band band-shifted to z=0.1 assuming h=1.0.
   ABSMAG01_SYNTH_IVAR_SDSS_U      float32                      1 / mag2 Inverse variance corresponding to ABSMAG01_SYNTH_SDSS_U.
         ABSMAG01_SDSS_G [4]_      float32                           mag Like ABSMAG01_SDSS_U but for SDSS *g*-band.
         ABSMAG01_IVAR_SDSS_G      float32                      1 / mag2 Like ABSMAG01_IVAR_SDSS_U but for SDSS *g*-band.
        ABSMAG01_SYNTH_SDSS_G      float32                           mag Like ABSMAG01_SYNTH_SDSS_U but for SDSS *g*-band.
   ABSMAG01_SYNTH_IVAR_SDSS_G      float32                      1 / mag2 Like ABSMAG01_SYNTH_IVAR_SDSS_U but for SDSS *g*-band.
         ABSMAG01_SDSS_R [4]_      float32                           mag Like ABSMAG01_SDSS_U but for SDSS *r*-band.
         ABSMAG01_IVAR_SDSS_R      float32                      1 / mag2 Like ABSMAG01_IVAR_SDSS_U but for SDSS *r*-band.
        ABSMAG01_SYNTH_SDSS_R      float32                           mag Like ABSMAG01_SYNTH_SDSS_U but for SDSS *r*-band.
   ABSMAG01_SYNTH_IVAR_SDSS_R      float32                      1 / mag2 Like ABSMAG01_SYNTH_IVAR_SDSS_U but for SDSS *r*-band.
         ABSMAG01_SDSS_I [4]_      float32                           mag Like ABSMAG01_SDSS_U but for SDSS *i*-band.
         ABSMAG01_IVAR_SDSS_I      float32                      1 / mag2 Like ABSMAG01_IVAR_SDSS_U but for SDSS *i*-band.
        ABSMAG01_SYNTH_SDSS_I      float32                           mag Like ABSMAG01_SYNTH_SDSS_U but for SDSS *i*-band.
   ABSMAG01_SYNTH_IVAR_SDSS_I      float32                      1 / mag2 Like ABSMAG01_SYNTH_IVAR_SDSS_U but for SDSS *i*-band.
         ABSMAG01_SDSS_Z [4]_      float32                           mag Like ABSMAG01_SDSS_U but for SDSS *z*-band.
         ABSMAG01_IVAR_SDSS_Z      float32                      1 / mag2 Like ABSMAG01_IVAR_SDSS_U but for SDSS *z*-band.
        ABSMAG01_SYNTH_SDSS_Z      float32                           mag Like ABSMAG01_SYNTH_SDSS_U but for SDSS *z*-band.
   ABSMAG01_SYNTH_IVAR_SDSS_Z      float32                      1 / mag2 Like ABSMAG01_SYNTH_IVAR_SDSS_U but for SDSS *z*-band.
             ABSMAG01_W1 [4]_      float32                           mag Absolute magnitude in WISE *W1*-band band-shifted to z=0.1 assuming h=1.0.
             ABSMAG01_IVAR_W1      float32                      1 / mag2 Inverse variance corresponding to ABSMAG01_W1.
            ABSMAG01_SYNTH_W1      float32                           mag Synthesized absolute magnitude in WISE *W1*-band band-shifted to z=0.1 assuming h=1.0.
       ABSMAG01_SYNTH_IVAR_W1      float32                      1 / mag2 Inverse variance corresponding to ABSMAG01_SYNTH_W1.
              KCORR10_DECAM_G      float32                           mag K-correction used to derive ABSMAG10_DECAM_G band-shifted to z=1.0.
              KCORR10_DECAM_R      float32                           mag Like KCORR10_DECAM_G but for DECam *r*-band.
              KCORR10_DECAM_Z      float32                           mag Like KCORR10_DECAM_G but for DECam *z*-band.
                    KCORR00_U      float32                           mag K-correction used to derive ABSMAG00_U band-shifted to z=0.0.
                    KCORR00_B      float32                           mag Like KCORR00_U but for Johnson/Cousins *B*-band.
                    KCORR00_V      float32                           mag Like KCORR00_U but for Johnson/Cousins *V*-band.
               KCORR01_SDSS_U      float32                           mag K-correction used to derive ABSMAG01_SDSS_U band-shifted to z=0.1.
               KCORR01_SDSS_G      float32                           mag Like KCORR01_SDSS_U but for SDSS *g*-band.
               KCORR01_SDSS_R      float32                           mag Like KCORR01_SDSS_U but for SDSS *r*-band.
               KCORR01_SDSS_I      float32                           mag Like KCORR01_SDSS_U but for SDSS *i*-band.
               KCORR01_SDSS_Z      float32                           mag Like KCORR01_SDSS_U but for SDSS *z*-band.
                   KCORR01_W1      float32                           mag K-correction used to derive ABSMAG01_W1 band-shifted to z=0.0.
                  LOGLNU_1500      float32            1e-28 erg / (Hz s) Monochromatic luminosity at 1500 A in the rest-frame.
             LOGLNU_1500_IVAR      float32           1e+56 Hz2 s2 / erg2 Inverse variance in LOGLNU_1500.
                  LOGLNU_2800      float32            1e-28 erg / (Hz s) Monochromatic luminosity at 2800 A in the rest-frame.
             LOGLNU_2800_IVAR      float32           1e+56 Hz2 s2 / erg2 Inverse variance in LOGLNU_2800.
                    LOGL_1450      float32                    1e-10 Lsun Integrated luminosity at 1450 A in the rest-frame.
               LOGL_1450_IVAR      float32                 1e+20 / Lsun2 Inverse variance in LOGL_1450.
                    LOGL_1700      float32                    1e-10 Lsun Integrated luminosity at 1700 A in the rest-frame.
               LOGL_1700_IVAR      float32                 1e+20 / Lsun2 Inverse variance in LOGL_1700.
                    LOGL_3000      float32                    1e-10 Lsun Integrated luminosity at 3000 A in the rest-frame.
               LOGL_3000_IVAR      float32                 1e+20 / Lsun2 Inverse variance in LOGL_3000.
                    LOGL_5100      float32                    1e-10 Lsun Integrated luminosity at 5100 A in the rest-frame.
               LOGL_5100_IVAR      float32                 1e+20 / Lsun2 Inverse variance in LOGL_5100.
               FLYA_1215_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 1215.67 A in the rest-frame.
          FLYA_1215_CONT_IVAR      float32 1e+34 cm4 Angstrom2 s2 / erg2 Inverse variance in FLYA_1215_CONT.
               FOII_3727_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 3728.483 A in the rest-frame.
          FOII_3727_CONT_IVAR      float32 1e+34 cm4 Angstrom2 s2 / erg2 Inverse variance in FOII_3727_CONT.
                  FHBETA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 4862.683 A in the rest-frame.
             FHBETA_CONT_IVAR      float32 1e+34 cm4 Angstrom2 s2 / erg2 Inverse variance in FHBETA_CONT.
              FOIII_5007_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 5008.239 A in the rest-frame.
         FOIII_5007_CONT_IVAR      float32 1e+34 cm4 Angstrom2 s2 / erg2 Inverse variance in FOIII_5007_CONT.
                 FHALPHA_CONT      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at 6564.613 A in the rest-frame.
            FHALPHA_CONT_IVAR      float32 1e+34 cm4 Angstrom2 s2 / erg2 Inverse variance in FHALPHA_CONT.
============================= ============ ============================= ============================================

HDU03
-----

EXTNAME = FASTSPEC

Table with spectral fitting results.

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

=========================== ============ ============================= ============================================
Name                        Type         Units                         Description
=========================== ============ ============================= ============================================
                   TARGETID        int64                               Unique target ID.
                SURVEY [1]_       bytes3                               Survey name.
               PROGRAM [1]_       bytes6                               Program name.
               HEALPIX [1]_        int32                               Healpixel number.
                TILEID [2]_        int32                               Tile ID number.
                 NIGHT [2]_        int32                               Night of observation.
                 FIBER [2]_        int32                               Fiber number.
                 EXPID [3]_        int32                               Exposure ID number.
                      SNR_B      float32                               Median signal-to-noise ratio per pixel in the *b* camera.
                      SNR_R      float32                               Like SNR_B but in the *r* camera.
                      SNR_Z      float32                               Like SNR_B but in *z* camera.
               SMOOTHCORR_B      float32                       percent Median value of the smooth continuum correction relative to the data in the *b* camera.
               SMOOTHCORR_R      float32                       percent Like SMOOTHCORR_B but in the *r* camera.
               SMOOTHCORR_Z      float32                       percent Like SMOOTHCORR_B but in the *z* camera.
                   APERCORR      float32                               Median aperture correction factor.
                 APERCORR_G      float32                               Aperture correction factor measured in the *g*-band.
                 APERCORR_R      float32                               Like APERCORR_G but measured in the *r*-band.
                 APERCORR_Z      float32                               Like APERCORR_G but measured in the *z*-band.
              INIT_SIGMA_UV      float32                        km / s Initial line width for the rest-frame UV emission lines.
          INIT_SIGMA_NARROW      float32                        km / s Like INIT_SIGMA_UV but for the forbidden plus narrow Balmer emission lines.
          INIT_SIGMA_BALMER      float32                        km / s Like INIT_SIGMA_UV but for broad Balmer emission lines.
             INIT_VSHIFT_UV      float32                        km / s Initial velocity shift for the UV emission lines.
         INIT_VSHIFT_NARROW      float32                        km / s Like INIT_VSHIFT_UV but for narrow (forbidden and narrow Balmer) emission lines.
         INIT_VSHIFT_BALMER      float32                        km / s Like INIT_VSHIFT_UV but for broad Balmer emission lines.
          INIT_BALMER_BROAD         bool                               Boolean flag indicating whether a broad Balmer emission line was initially identified in the spectral range.
             DELTA_LINECHI2      float32                               Chi-squared difference between an emission-line model without and with broad lines.
             DELTA_LINENDOF        int32                               Difference in the degrees of freedom between an emission-line model without and with broad lines.
         MGII_DOUBLET_RATIO      float32                               MgII 2796 / 2803 doublet line-ratio.
    MGII_DOUBLET_RATIO_IVAR      float32                               Inverse variance in MGII_DOUBLET_RATIO.
          OII_DOUBLET_RATIO      float32                               [OII] 3726 / 3729 doublet line-ratio.
     OII_DOUBLET_RATIO_IVAR      float32                               Inverse variance in OII_DOUBLET_RATIO.
         OIII_DOUBLET_RATIO      float32                               [OIII] 5007 / 4959 doublet line-ratio.
    OIII_DOUBLET_RATIO_IVAR      float32                               Inverse variance in OIII_DOUBLET_RATIO.
          SII_DOUBLET_RATIO      float32                               [SII] 6731 / 6716 doublet line-ratio.
     SII_DOUBLET_RATIO_IVAR      float32                               Inverse variance in SII_DOUBLET_RATIO.
          NII_DOUBLET_RATIO      float32                               [NII] 6584 / 6548 doublet line-ratio.
     NII_DOUBLET_RATIO_IVAR      float32                               Inverse variance in NII_DOUBLET_RATIO.
       OIIRED_DOUBLET_RATIO      float32                               [OII] 7320 / 7330 doublet (actually a quadruplet) line-ratio.
  OIIRED_DOUBLET_RATIO_IVAR      float32                               Inverse variance in OIIRED_DOUBLET_RATIO.
     LINENAME_MODELAMP [6]_      float32  1e-17 erg / (Angstrom cm2 s) Model emission-line amplitude (before convolution with the resolution matrix).
          LINENAME_AMP [6]_      float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
     LINENAME_AMP_IVAR [6]_      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance in LINENAME_AMP.
         LINENAME_FLUX [6]_      float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
    LINENAME_FLUX_IVAR [6]_      float32           1e+34 cm4 s2 / erg2 Inverse variance in LINENAME_FLUX.
      LINENAME_BOXFLUX [6]_      float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
 LINENAME_BOXFLUX_IVAR [6]_      float32           1e+34 cm4 s2 / erg2 Inverse variance in LINENAME_BOXFLUX.
       LINENAME_VSHIFT [6]_      float32                        km / s Velocity shift relative to Z.
  LINENAME_VSHIFT_IVAR [6]_      float32                      s2 / km2 Inverse variance in LINENAME_VSHIFT.
        LINENAME_SIGMA [6]_      float32                        km / s Gaussian emission-line width (before convolution with the resolution matrix).
   LINENAME_SIGMA_IVAR [6]_      float32                      s2 / km2 Inverse variance in LINENAME_SIGMA.
         LINENAME_CONT [6]_      float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
    LINENAME_CONT_IVAR [6]_      float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance in LINENAME_CONT.
           LINENAME_EW [6]_      float32                      Angstrom Rest-frame emission-line equivalent width.
      LINENAME_EW_IVAR [6]_      float32                 1 / Angstrom2 Inverse variance in LINENAME_EW.
   LINENAME_FLUX_LIMIT [6]_      float32                 erg / (cm2 s) One-sigma upper-limit on the emission line flux.
         LINENAME_CHI2 [6]_      float32                               Chi-squared of the line-fit.
         LINENAME_NPIX [6]_        int32                               Number of pixels attributed to the emission line.
     LINENAME2_MOMENT1 [7]_      float32                      Angstrom First moment of the continuum-subtracted pixel values within +/- 5-sigma of the Gaussian line-center.
LINENAME2_MOMENT1_IVAR [7]_      float32                 1 / Angstrom2 Inverse variance in LINENAME2_MOMENT1.
     LINENAME2_MOMENT2 [7]_      float32                     Angstrom2 Second moment of the flux relative to LINENAME2_MOMENT1.
LINENAME2_MOMENT2_IVAR [7]_      float32                 1 / Angstrom4 Inverse variance in LINENAME2_MOMENT2.
     LINENAME2_MOMENT3 [7]_      float32                     Angstrom3 Third moment of the flux relative to LINENAME2_MOMENT1.
LINENAME2_MOMENT3_IVAR [7]_      float32                 1 / Angstrom6 Inverse variance in LINENAME2_MOMENT3.
=========================== ============ ============================= ============================================

HDU04
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

.. [2] Column only present when fitting cumulative, per-night, or per-expopsure
       tile-based coadds.

.. [3] Column only present when fitting per-exposure tile-based coadds.

.. [4] Only observed photometry with a minimum signal-to-noise ratio of two is
       used to compute K-corrections. Absolute magnitudes are estimated using
       from the (S/N>2) observed-frame bandpass nearest in wavelength to the
       desired band-shifted rest-frame bandpass. If no observed-frame photometry
       is available, then the absolute magnitude is synthesized from the
       best-fitting model and the corresponding inverse variance is set to zero.

.. [5] Column only present in Commissioning and Survey Validation spectroscopic
       observations.

.. [6] `LINENAME` represents the following modeled emission lines: NV_1240,
       OI_1304, SILIV_1396, CIV_1549, HEII_1640, ALIII_1857, SILIII_1892,
       CIII_1908, MGII_2796, MGII_2803, NEV_3346, NEV_3426, OII_3726, OII_3729,
       NEIII_3869, H6, H6_BROAD, HEPSILON, HEPSILON_BROAD, HDELTA, HDELTA_BROAD,
       HGAMMA, HGAMMA_BROAD, OIII_4363, HEI_4471, HEII_4686, HBETA, HBETA_BROAD,
       OIII_4959, OIII_5007, NII_5755, HEI_5876, OI_6300, SIII_6312, NII_6548,
       HALPHA, HALPHA_BROAD, NII_6584, SII_6716, SII_6731, ARIII_7135, OII_7320,
       OII_7330, SIII_9069, and SIII_9532.

.. [7] `LINENAME2` represents the following emission lines: CIV_1549, MGII_2800
       (using the mean wavelength of the doublet), HBETA, and OIII_5007.
