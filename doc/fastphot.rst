.. _fastphot datamodel:

===================
fastphot Data Model
===================

:Summary: Broadband photometric fitting results.
:Naming Convention:
    ``fastphot-SURVEY-PROGRAM-HEALPIX.fits``, where
    ``SURVEY`` is the survey name (e.g., *main* or *sv3*), ``PROGRAM`` is the
    program name (e.g., *bright* or *dark*), and ``HEALPIX`` is the healpixel number.
:Regex: ``fastphot-(cmx|main|special|sv1|sv2|sv3)-(backup|bright|dark|other)-[0-9]+\.fits``
:File Type: FITS

Contents
========

====== ============ ======== ======================
Number EXTNAME      Type     Contents
====== ============ ======== ======================
HDU00_ PRIMARY      IMAGE    Keywords only.
HDU01_ FASTPHOT     BINTABLE Table with photometric fitting results.
HDU02_ METADATA     BINTABLE Table with sample metadata.
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
    CHECKSUM UKXHaIWHUIWHaIWH str  HDU checksum updated 2022-08-02T00:24:13
    DATASUM  0                str  data unit checksum updated 2022-08-02T00:24:13
    ======== ================ ==== ==============================================

Empty HDU.

HDU01
-----

EXTNAME = FASTPHOT

Broadband photometric fitting results.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

.. collapse:: Required Header Keywords Table
    :open:

    .. rst-class:: keywords

    ======== ================ ==== ==============================================
    KEY      Example Value    Type Comment
    ======== ================ ==== ==============================================
    NAXIS1   335              int  length of dimension 1
    NAXIS2   338              int  length of dimension 2
    CHECKSUM PI3dSH0aPH0aPH0a str  HDU checksum updated 2022-08-02T00:24:13
    DATASUM  72956540         str  data unit checksum updated 2022-08-02T00:24:13
    ======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: columns

======================== ============ ============================ ===================================================
Name                     Type         Units                        Description
======================== ============ ============================ ===================================================
                TARGETID int64                                     Unique target ID.
             SURVEY [1]_ bytes3                                    Survey name.
            PROGRAM [1]_ bytes6                                    Program name.
            HEALPIX [1]_ int32                                     Healpixel number.
             TILEID [2]_ int32                                     Tile ID number.
              NIGHT [2]_ int32                                     Night of observation.
              FIBER [2]_ int32                                     Fiber number.
              EXPID [3]_ int32                                     Exposure ID number.
                       Z float64                                   Stellar continuum redshift.
                   COEFF float32[168]                              Continuum coefficients.
              RCHI2_PHOT float32                                   Reduced chi-squared of the fit to the broadband photometry.
                   VDISP float32                            km / s Stellar velocity dispersion.
                      AV float32                               mag Attenuation of the integrated stellar population.
                     AGE float32                               Gyr Light-weighted age.
                   ZZSUN float32                                   Logarithmic stellar metallicity relative to solar.
                LOGMSTAR float32                              Msun Logarithmic stellar mass (h=1.0, Chabrier+2003 initial mass function).
                     SFR float32                         Msun / yr Instantaneous star formation rate (h=1.0, Chabrier+2003 initial mass function).
            DN4000_MODEL float32                                   Narrow 4000-A break index (from Balogh et al. 1999) measured from the best-fitting continuum model.
  FLUX_SYNTH_PHOTMODEL_G float32                              nmgy g-band flux synthesized from the best-fitting photometric continuum model.
  FLUX_SYNTH_PHOTMODEL_R float32                              nmgy r-band flux synthesized from the best-fitting photometric continuum model.
  FLUX_SYNTH_PHOTMODEL_Z float32                              nmgy z-band flux synthesized from the best-fitting photometric continuum model.
 FLUX_SYNTH_PHOTMODEL_W1 float32                              nmgy W1-band flux synthesized from the best-fitting photometric continuum model.
 FLUX_SYNTH_PHOTMODEL_W2 float32                              nmgy W2-band flux synthesized from the best-fitting photometric continuum model.
 FLUX_SYNTH_PHOTMODEL_W3 float32                              nmgy W3-band flux synthesized from the best-fitting photometric continuum model.
 FLUX_SYNTH_PHOTMODEL_W4 float32                              nmgy W4-band flux synthesized from the best-fitting photometric continuum model.
                 KCORR_U float32                               mag K-correction used to derive ABSMAG_U band-shifted to z=0.0 assuming h=1.0.
                ABSMAG_U float32                               mag Absolute magnitude in Johnson/Cousins U-band band-shifted to z=0.0 assuming h=1.0.
           ABSMAG_IVAR_U float32                          1 / mag2 Inverse variance corresponding to ABSMAG_U.
                 KCORR_B float32                               mag Like KCORR_U but for Johnson/Cousins B-band.
                ABSMAG_B float32                               mag Like ABSMAG_U but for Johnson/Cousins B-band.
           ABSMAG_IVAR_B float32                          1 / mag2 Like ABSMAG_IVAR_U but for Johnson/Cousins B-band.
                 KCORR_V float32                               mag Like KCORR_U but for Johnson/Cousins V-band.
                ABSMAG_V float32                               mag Like ABSMAG_U but for Johnson/Cousins V-band.
           ABSMAG_IVAR_V float32                          1 / mag2 Like ABSMAG_IVAR_U but for Johnson/Cousins V-band.
            KCORR_SDSS_U float32                               mag K-correction used to derive ABSMAG_SDSS_U band-shifted to z=0.1 assuming h=1.0.
           ABSMAG_SDSS_U float32                               mag Absolute magnitude in SDSS u-band band-shifted to z=0.1 assuming h=1.0.
      ABSMAG_IVAR_SDSS_U float32                          1 / mag2 Inverse variance corresponding to ABSMAG_SDSS_U.
            KCORR_SDSS_G float32                               mag Like KCORR_SDSS_U but for SDSS g-band.
           ABSMAG_SDSS_G float32                               mag Like ABSMAG_SDSS_U but for SDSS g-band.
      ABSMAG_IVAR_SDSS_G float32                          1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS g-band.
            KCORR_SDSS_R float32                               mag Like KCORR_SDSS_U but for SDSS r-band.
           ABSMAG_SDSS_R float32                               mag Like ABSMAG_SDSS_U but for SDSS r-band.
      ABSMAG_IVAR_SDSS_R float32                          1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS r-band.
            KCORR_SDSS_I float32                               mag Like KCORR_SDSS_U but for SDSS i-band.
           ABSMAG_SDSS_I float32                               mag Like ABSMAG_SDSS_U but for SDSS i-band.
      ABSMAG_IVAR_SDSS_I float32                          1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS i-band.
            KCORR_SDSS_Z float32                               mag Like KCORR_SDSS_U but for SDSS z-band.
           ABSMAG_SDSS_Z float32                               mag Like ABSMAG_SDSS_U but for SDSS z-band.
      ABSMAG_IVAR_SDSS_Z float32                          1 / mag2 Like ABSMAG_IVAR_SDSS_U but for SDSS z-band.
                KCORR_W1 float32                               mag K-correction used to derive ABSMAG_W1 band-shifted to z=0.0 assuming h=1.0.
               ABSMAG_W1 float32                               mag Absolute magnitude in WISE W1-band band-shifted to z=0.0 assuming h=1.0.
          ABSMAG_IVAR_W1 float32                          1 / mag2 Inverse variance corresponding to ABSMAG_W1.
                KCORR_W2 float32                               mag K-correction used to derive ABSMAG_W2 band-shifted to z=0.0 assuming h=1.0.
               ABSMAG_W2 float32                               mag Absolute magnitude in WISE W2-band band-shifted to z=0.0 assuming h=1.0.
          ABSMAG_IVAR_W2 float32                          1 / mag2 Inverse variance corresponding to ABSMAG_W2.
             LOGLNU_1500 float32                1e-28 erg / (s Hz) Monochromatic luminosity at 1500 A in the rest-frame.
             LOGLNU_2800 float32                1e-28 erg / (s Hz) Monochromatic luminosity at 2800 A in the rest-frame.
               LOGL_5100 float32                        1e+10 Lsun Total luminosity at 5100 A in the rest-frame.
          FOII_3727_CONT float32      1e-17 erg / (Angstrom cm2 s) Continuum flux at 3728.483 A in the rest-frame.
             FHBETA_CONT float32      1e-17 erg / (Angstrom cm2 s) Continuum flux at 4862.683 A in the rest-frame.
         FOIII_5007_CONT float32      1e-17 erg / (Angstrom cm2 s) Continuum flux at 5008.239 A in the rest-frame.
            FHALPHA_CONT float32      1e-17 erg / (Angstrom cm2 s) Continuum flux at 6564.613 A in the rest-frame.
======================== ============ ============================ ===================================================

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
    CHECKSUM iFX6iFW6iFW6iFW6 str  HDU checksum updated 2022-08-02T00:24:13
    DATASUM  1759692941       str  data unit checksum updated 2022-08-02T00:24:13
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
               PHOTSYS  bytes1                Photometric system (*N* or *S*).
                 LS_ID   int64                Unique Legacy Surveys identification number.
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

.. [1] Column only present when fitting healpix coadds.
       
.. [2] Column only present when fitting cumulative, per-night, or per-expopsure tile-based coadds.
       
.. [3] Column only present when fitting per-exposure tile-based coadds.

.. [4] Column only present in Commissioning and Survey Validation spectroscopic observations.
