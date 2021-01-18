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

========================= =========== ============================= ==========================================
Name                      Type        Units                         Description
========================= =========== ============================= ==========================================
                 TARGETID       int64                               Unique target ID.
                TARGET_RA     float64                           deg Right ascension from fibermap.
               TARGET_DEC     float64                           deg Declination from fibermap.
                    NIGHT       int32                               Write me.
                   TILEID       int32                               Write me.
                    FIBER       int32                               Write me.
                    EXPID       int32                               Write me.
                        Z     float64                               Redrock redshift.
                DELTACHI2     float64                               Redrock delta-chi-squared.
            PHOTSYS_SOUTH        bool                               Photometric system.
          SV1_DESI_TARGET       int64                               Write me.
           SV1_BGS_TARGET       int64                               Write me.
           SV1_MWS_TARGET       int64                               Write me.
              DESI_TARGET       int64                               Write me.
               BGS_TARGET       int64                               Write me.
               MWS_TARGET       int64                               Write me.
            CONTINUUM_SNR  float32[3]                               Write me.
              CONTINUUM_Z     float64                               Write me.
          CONTINUUM_COEFF float64[11]                               Write me.
           CONTINUUM_CHI2     float32                               Write me.
            CONTINUUM_AGE     float32                           Gyr Write me.
             CONTINUUM_AV     float32                           mag Write me.
        CONTINUUM_AV_IVAR     float32                      1 / mag2 Write me.
          CONTINUUM_VDISP     float32                        km / s Write me.
     CONTINUUM_VDISP_IVAR     float32                      s2 / km2 Write me.
                    D4000     float32                               Write me.
               D4000_IVAR     float32                               Write me.
              D4000_MODEL     float32                               Write me.
             FLUX_SYNTH_G     float32                          nmgy Write me.
             FLUX_SYNTH_R     float32                          nmgy Write me.
             FLUX_SYNTH_Z     float32                          nmgy Write me.
       FLUX_SYNTH_MODEL_G     float32                          nmgy Write me.
       FLUX_SYNTH_MODEL_R     float32                          nmgy Write me.
       FLUX_SYNTH_MODEL_Z     float32                          nmgy Write me.
     LINEVSHIFT_FORBIDDEN     float32                               Write me.
LINEVSHIFT_FORBIDDEN_IVAR     float32                               Write me.
        LINEVSHIFT_BALMER     float32                               Write me.
   LINEVSHIFT_BALMER_IVAR     float32                               Write me.
      LINESIGMA_FORBIDDEN     float32                        km / s Write me.
 LINESIGMA_FORBIDDEN_IVAR     float32                      s2 / km2 Write me.
         LINESIGMA_BALMER     float32                        km / s Write me.
    LINESIGMA_BALMER_IVAR     float32                      s2 / km2 Write me.
             OII_3726_AMP     float32  1e-17 erg / (Angstrom cm2 s) Write me.
        OII_3726_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Write me.
            OII_3726_FLUX     float32           1e-17 erg / (cm2 s) Write me.
       OII_3726_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Write me.
         OII_3726_BOXFLUX     float32           1e-17 erg / (cm2 s) Write me.
    OII_3726_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Write me.
            OII_3726_CONT     float32  1e-17 erg / (Angstrom cm2 s) Write me.
       OII_3726_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Write me.
              OII_3726_EW     float32                      Angstrom Write me.
         OII_3726_EW_IVAR     float32                 1 / Angstrom2 Write me.
      OII_3726_FLUX_LIMIT     float32                 erg / (cm2 s) Write me.
        OII_3726_EW_LIMIT     float32                      Angstrom Write me.
            OII_3726_CHI2     float32                               Write me.
            OII_3726_NPIX       int32                               Write me.
========================= =========== ============================= ==========================================

Notes and Examples
==================


Upcoming changes
================
