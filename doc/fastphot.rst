===================
fastphot Data Model
===================

:Summary: Photometric fitting results.
:Naming Convention: ``fastphot-{petal}-{tileid}-{night}.fits``, where
    ``{petal}`` is the petal or spetrograph number (0-9), ``{tileid}`` is the
    tileid and ``{night}`` is the night of the observation.
:Regex: ``fastphot-[0-9]+*+\.fits``
:File Type: FITS, <1 MB

Contents
========

====== ============ ======== ======================
Number EXTNAME      Type     Contents
====== ============ ======== ======================
HDU00_ PRIMARY      IMAGE    Empty
HDU01_ FASTPHOT     BINTABLE Fitting results table.
HDU02_ METADATA     BINTABLE Object metadata table.
====== ============ ======== ======================

FITS Header Units
=================

HDU00
-----

EXTNAME = PRIMARY

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================ ==== ==============================================
KEY      Example Value    Type Comment
======== ================ ==== ==============================================
CHECKSUM EAnFF7l9EAlEE5l9 str  HDU checksum updated 2018-03-29T22:45:34
DATASUM  0                str  data unit checksum updated 2018-03-29T22:45:34
======== ================ ==== ==============================================

Empty HDU.

HDU01
-----

EXTNAME = FASTPHOT

Fitting results. Checksum not yet implemented.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================ ==== ==============================================
KEY      Example Value    Type Comment
======== ================ ==== ==============================================
NAXIS1   188              int  length of dimension 1
NAXIS2   3000             int  length of dimension 2
CHECKSUM EAnFF7l9EAlEE5l9 str  HDU checksum updated 2018-03-29T22:45:34
DATASUM  0                str  data unit checksum updated 2018-03-29T22:45:34
======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

====================== =========== ========== ==========================================
Name                   Type        Units      Description
====================== =========== ========== ==========================================
              TARGETID   int64                Unique target ID.
       CONTINUUM_COEFF float64[11]            Continuum coefficients.
        CONTINUUM_CHI2 float32                Reduced chi^2 of the continuum fit.
         CONTINUUM_AGE float32            Gyr Light-weighted age.
          CONTINUUM_AV float32            mag Intrinsic attenuation.
     CONTINUUM_AV_IVAR float32     1 / mag^2  Inverse variance of CONTINUUM_AV.
           D4000_MODEL float32                4000-A break index from the best-fitting continuum model.
               KCORR_U float32            mag K-correction used to derive ABSMAG_U.
              ABSMAG_U float32            mag Absolute magnitude in DECam u-band.
         ABSMAG_IVAR_U float32      1 / mag^2 Inverse variance corresponding to ABSMAG_U.
               KCORR_G float32            mag Like KCORR_U but for DECam g-band.
              ABSMAG_G float32            mag Like ABSMAG_U but for DECam g-band.
         ABSMAG_IVAR_G float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam r-band.
               KCORR_R float32            mag Like KCORR_U but for DECam r-band.
              ABSMAG_R float32            mag Like ABSMAG_U but for DECam r-band.
         ABSMAG_IVAR_R float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam i-band.
               KCORR_I float32            mag Like KCORR_U but for DECam i-band.
              ABSMAG_I float32            mag Like ABSMAG_U but for DECam i-band.
         ABSMAG_IVAR_I float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam z-band.
               KCORR_Z float32            mag Like KCORR_U but for DECam z-band.
              ABSMAG_Z float32            mag Like ABSMAG_U but for DECam z-band.
         ABSMAG_IVAR_Z float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam g-band.
              KCORR_W1 float32            mag Like KCORR_U but for DECam W1-band.
             ABSMAG_W1 float32            mag Like ABSMAG_U but for DECam W1-band.
        ABSMAG_IVAR_W1 float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam W1-band.
====================== =========== ========== ==========================================

HDU02
-----

EXTNAME = METADATA

Fitting results. Checksum not yet implemented.

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
        FIBERTOTFLUX_G float32           nmgy Fibertot g-band flux from targeting catalog.
        FIBERTOTFLUX_R float32           nmgy Fibertot r-band flux from targeting catalog.
        FIBERTOTFLUX_Z float32           nmgy Fibertot z-band flux from targeting catalog.
                FLUX_G float32           nmgy Total g-band flux from targeting catalog.
                FLUX_R float32           nmgy Total r-band flux from targeting catalog.
                FLUX_Z float32           nmgy Total z-band flux from targeting catalog.
               FLUX_W1 float32           nmgy Total W1-band flux from targeting catalog.
               FLUX_W2 float32           nmgy Total W2-band flux from targeting catalog.
           FLUX_IVAR_G float32     1 / nmgy^2 Inverse variance of FLUX_G from targeting catalog.
           FLUX_IVAR_R float32     1 / nmgy^2 Inverse variance of FLUX_R from targeting catalog.
           FLUX_IVAR_Z float32     1 / nmgy^2 Inverse variance of FLUX_Z from targeting catalog.
          FLUX_IVAR_W1 float32     1 / nmgy^2 Inverse variance of FLUX_W1 from targeting catalog.
          FLUX_IVAR_W2 float32     1 / nmgy^2 Inverse variance of FLUX_W2 from targeting catalog.
====================== =========== ========== ==========================================

Notes and Examples
==================

If the inverse variance on a given absolutely magnitude is zero it means that
the absolute magnitude was derived from *synthesized* photometry based on the
best-fitting model (i.e., use with care).

Upcoming changes
================
