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

====== ============ ======== =====================
Number EXTNAME      Type     Contents
====== ============ ======== =====================
HDU00_ PRIMARY      IMAGE    Empty
HDU01_ RESULTS      BINTABLE Fitting results table
====== ============ ======== =====================

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

EXTNAME = RESULTS

Fitting results. Checksum not yet implemented.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================ ==== ==============================================
KEY      Example Value    Type Comment
======== ================ ==== ==============================================
NAXIS1   325              int  length of dimension 1
NAXIS2   1225             int  length of dimension 2
CHECKSUM EAnFF7l9EAlEE5l9 str  HDU checksum updated 2018-03-29T22:45:34
DATASUM  0                str  data unit checksum updated 2018-03-29T22:45:34
======== ================ ==== ==============================================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================= =========== ========== ==========================================
Name                              Type        Units      Description
================================= =========== ========== ==========================================
              TARGETID              int64                Unique target ID.
             TARGET_RA            float64            deg Right ascension from fibermap.
            TARGET_DEC            float64            deg Declination from fibermap.
                 NIGHT              int32                Write me.
                TILEID              int32                Write me.
                 FIBER              int32                Write me.
                 EXPID              int32                Write me.
                     Z            float64                Redrock redshift.
             DELTACHI2            float64                Redrock delta-chi-squared.
         PHOTSYS_SOUTH               bool                Photometric system.
       SV1_DESI_TARGET              int64                Write me.
        SV1_BGS_TARGET              int64                Write me.
        SV1_MWS_TARGET              int64                Write me.
           DESI_TARGET              int64                Write me.
            BGS_TARGET              int64                Write me.
            MWS_TARGET              int64                Write me.
  CONTINUUM_PHOT_COEFF            float64[11]            Continuum coefficients.
   CONTINUUM_PHOT_CHI2            float32                Write me.
    CONTINUUM_PHOT_AGE            float32            Gyr Write me.
     CONTINUUM_PHOT_AV            float32            mag Write me.
CONTINUUM_PHOT_AV_IVAR            float32     1 / mag^2  Write me.
      D4000_MODEL_PHOT            float32                Write me.
        FIBERTOTFLUX_G            float32           nmgy Write me.
        FIBERTOTFLUX_R            float32           nmgy Write me.
        FIBERTOTFLUX_Z            float32           nmgy Write me.
                FLUX_G            float32           nmgy Write me.
           FLUX_IVAR_G            float32     1 / nmgy^2 Write me.
                FLUX_R            float32           nmgy Write me.
           FLUX_IVAR_R            float32     1 / nmgy^2 Write me.
                FLUX_Z            float32           nmgy Write me.
           FLUX_IVAR_Z            float32     1 / nmgy^2 Write me.
               FLUX_W1            float32           nmgy Write me.
          FLUX_IVAR_W1            float32     1 / nmgy^2 Write me.
               FLUX_W2            float32           nmgy Write me.
          FLUX_IVAR_W2            float32     1 / nmgy^2 Write me.
               KCORR_U            float32            mag K-correction used to derive ABSMAG_U.
              ABSMAG_U            float32            mag Absolute magnitude in DECam u-band.
         ABSMAG_IVAR_U            float32      1 / mag^2 Inverse variance corresponding to ABSMAG_U.
               KCORR_G            float32            mag Like KCORR_U but for DECam g-band.
              ABSMAG_G            float32            mag Like ABSMAG_U but for DECam g-band.
         ABSMAG_IVAR_G            float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam r-band.
               KCORR_R            float32            mag Like KCORR_U but for DECam r-band.
              ABSMAG_R            float32            mag Like ABSMAG_U but for DECam r-band.
         ABSMAG_IVAR_R            float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam i-band.
               KCORR_I            float32            mag Like KCORR_U but for DECam i-band.
              ABSMAG_I            float32            mag Like ABSMAG_U but for DECam i-band.
         ABSMAG_IVAR_I            float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam z-band.
               KCORR_Z            float32            mag Like KCORR_U but for DECam z-band.
              ABSMAG_Z            float32            mag Like ABSMAG_U but for DECam z-band.
         ABSMAG_IVAR_Z            float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam g-band.
              KCORR_W1            float32            mag Like KCORR_U but for DECam W1-band.
             ABSMAG_W1            float32            mag Like ABSMAG_U but for DECam W1-band.
        ABSMAG_IVAR_W1            float32      1 / mag^2 Like ABSMAG_IVAR_U but for DECam W1-band.
================================= =========== ========== ==========================================

Notes and Examples
==================

If the inverse variance on a given absolutely magnitude is zero it means that
the absolute magnitude was derived from *synthesized* photometry based on the
best-fitting model (i.e., use with care).

Upcoming changes
================
