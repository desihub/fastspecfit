===================
fastspec Data Model
===================

:Summary: Spectroscopic fitting results.
:Naming Convention: ``fastspec-{petal}-{tileid}-{night}.fits``, where
    ``{petal}`` is the petal or spetrograph number (0-9), ``{tileid}`` is the
    tileid and ``{night}`` is the night of the observation.
:Regex: ``fastspec-[0-9]+*+\.fits``
:File Type: FITS, <1 MB

Contents
========

====== ============ ======== ======================
Number EXTNAME      Type     Contents
====== ============ ======== ======================
HDU00_ PRIMARY      IMAGE    Empty
HDU01_ FASTSPEC     BINTABLE Fitting results table.
HDU02_ METADATA     BINTABLE Object metadata table.
====== ============ ======== ======================

FITS Header Units
=================

HDU00
-----

EXTNAME = PRIMARY

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

Checksum not yet implemented.

======== ================ ==== ==============================================
KEY      Example Value    Type Comment
======== ================ ==== ==============================================
CHECKSUM EAnFF7l9EAlEE5l9 str  HDU checksum updated 2018-03-29T22:45:34
DATASUM  0                str  data unit checksum updated 2018-03-29T22:45:34
======== ================ ==== ==============================================

Empty HDU.

HDU01
-----

EXTNAME = FASTSPEC

Fitting results.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================ ==== ==============================================
KEY      Example Value    Type Comment
======== ================ ==== ==============================================
NAXIS1   325              int  length of dimension 1
NAXIS2   3000             int  length of dimension 2
CHECKSUM EAnFF7l9EAlEE5l9 str  HDU checksum updated 2018-03-29T22:45:34
DATASUM  0                str  data unit checksum updated 2018-03-29T22:45:34
======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

========================= =========== ============================= ============================================
Name                      Type        Units                         Description
========================= =========== ============================= ============================================
                 TARGETID       int64                               Unique target ID.
              CONTINUUM_Z     float64                               Stellar continuum redshift.
          CONTINUUM_COEFF float64[11]                               Continuum coefficients.
           CONTINUUM_CHI2     float32                               Reduced chi^2 of the continuum fit.
            CONTINUUM_AGE     float32                           Gyr Light-weighted age.
             CONTINUUM_AV     float32                           mag Intrinsic attenuation.
        CONTINUUM_AV_IVAR     float32                      1 / mag2 Inverse variance of CONTINUUM_AV.
          CONTINUUM_VDISP     float32                        km / s Stellar velocity dispersion.
     CONTINUUM_VDISP_IVAR     float32                      s2 / km2 Inverse variance of CONTINUUM_VDISP.
            CONTINUUM_SNR  float32[3]                               Median signal-to-noise ratio in each camera.
                    D4000     float32                               4000-A break index from the data.
               D4000_IVAR     float32                               Inverse variance of D4000.
              D4000_MODEL     float32                               4000-A break index from the best-fitting continuum model.
             FLUX_SYNTH_G     float32                          nmgy g-band flux synthesized from the data.
             FLUX_SYNTH_R     float32                          nmgy r-band flux synthesized from the data.
             FLUX_SYNTH_Z     float32                          nmgy z-band flux synthesized from the data.
       FLUX_SYNTH_MODEL_G     float32                          nmgy g-band flux synthesized from the best-fitting continuum model.
       FLUX_SYNTH_MODEL_R     float32                          nmgy r-band flux synthesized from the best-fitting continuum model.
       FLUX_SYNTH_MODEL_Z     float32                          nmgy z-band flux synthesized from the best-fitting continuum model.
     LINEVSHIFT_FORBIDDEN     float32                        km / s Velocity shift of the forbidden lines relative to the redrock redshift.
LINEVSHIFT_FORBIDDEN_IVAR     float32                      s2 / km2 Inverse variance of LINEVSHIFT_FORBIDDEN.
        LINEVSHIFT_BALMER     float32                        km / s Velocity shift of the Balmer lines relative to the redrock redshift.
   LINEVSHIFT_BALMER_IVAR     float32                      s2 / km2 Inverse variance of LINEVSHIFT_BALMER.
      LINESIGMA_FORBIDDEN     float32                        km / s Line-width of forbidden lines.
 LINESIGMA_FORBIDDEN_IVAR     float32                      s2 / km2 Inverse variance of LINESIGMA_FORBIDDEN.
         LINESIGMA_BALMER     float32                        km / s Line-width of Balmer lines.
    LINESIGMA_BALMER_IVAR     float32                      s2 / km2 Inverse variance of LINESIGMA_BALMER.
             OII_3726_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        OII_3726_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OII_3726_AMP.
            OII_3726_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       OII_3726_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OII_3726_FLUX
         OII_3726_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    OII_3726_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OII_3726_BOXFLUX
            OII_3726_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       OII_3726_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OII_3726_CONT
              OII_3726_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         OII_3726_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of OII_3726_EW
      OII_3726_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        OII_3726_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            OII_3726_CHI2     float32                               Reduced chi^2 of the line-fit.
            OII_3726_NPIX       int32                               Number of pixels attributed to the emission line.
             OII_3729_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        OII_3729_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OII_3729_AMP.
            OII_3729_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       OII_3729_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OII_3729_FLUX
         OII_3729_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    OII_3729_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OII_3729_BOXFLUX
            OII_3729_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       OII_3729_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OII_3729_CONT
              OII_3729_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         OII_3729_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of OII_3729_EW
      OII_3729_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        OII_3729_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            OII_3729_CHI2     float32                               Reduced chi^2 of the line-fit.
            OII_3729_NPIX       int32                               Number of pixels attributed to the emission line.
========================= =========== ============================= ============================================

HDU02
-----

EXTNAME = METADATA

Metadata associated with each objected fitted.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================ ==== ==============================================
KEY      Example Value    Type Comment
======== ================ ==== ==============================================
NAXIS1   155              int  length of dimension 1
NAXIS2   3000             int  length of dimension 2
SPECPROD daily            str  spectroscopic production name
CHECKSUM EAnFF7l9EAlEE5l9 str  HDU checksum updated 2018-03-29T22:45:34
DATASUM  0                str  data unit checksum updated 2018-03-29T22:45:34
======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

====================== =========== ========== ==========================================
Name                   Type        Units      Description
====================== =========== ========== ==========================================
              TARGETID   int64                Unique target ID.
                    RA float64            deg Right ascension from target catalog.
                   DEC float64            deg Declination from target catalog.
                 FIBER   int32                Fiber ID number.
                TILEID   int32                Tile ID number.
                 NIGHT   int32                Night (not present when fitting deep coadds).
                 EXPID   int32                Exposure ID number (not present when fitting coadds).
               PHOTSYS    str1                Photometric system ('N' or 'S').
       SV1_DESI_TARGET   int64                SV1 DESI targeting bit.
        SV1_BGS_TARGET   int64                SV1 BGS targeting bit.
        SV1_MWS_TARGET   int64                SV1 MWS targeting bit.
           DESI_TARGET   int64                DESI targeting bit.
            BGS_TARGET   int64                BGS targeting bit.
            MWS_TARGET   int64                MWS targeting bit.
                     Z float64                Redrock redshift.
             DELTACHI2 float64                Redrock delta-chi-squared.
              SPECTYPE    str6                Redrock spectral classification.
====================== =========== ========== ==========================================

Notes and Examples
==================


Upcoming changes
================
