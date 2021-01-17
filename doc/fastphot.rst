=============
fastphot.fits
=============

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

================================= =========== ======== ===========
Name                              Type        Units    Description
================================= =========== ======== ===========
              TARGETID              int64              Write me.
                    RA            float64          deg Write me.
                   DEC            float64          deg Write me.
                 NIGHT              int32              Write me.
                TILEID              int32              Write me.
                 FIBER              int32              Write me.
                 EXPID              int32              Write me.
                     Z            float64              Write me.
             DELTACHI2            float64              Write me.
         PHOTSYS_SOUTH               bool              Write me.
       SV1_DESI_TARGET              int64              Write me.
        SV1_BGS_TARGET              int64              Write me.
        SV1_MWS_TARGET              int64              Write me.
           DESI_TARGET              int64              Write me.
            BGS_TARGET              int64              Write me.
            MWS_TARGET              int64              Write me.
  CONTINUUM_PHOT_COEFF            float64[11]          Write me.
   CONTINUUM_PHOT_CHI2            float32              Write me.
    CONTINUUM_PHOT_AGE            float32          Gyr Write me.
     CONTINUUM_PHOT_AV            float32          mag Write me.
CONTINUUM_PHOT_AV_IVAR            float32     1 / mag2 Write me.
      D4000_MODEL_PHOT            float32              Write me.
        FIBERTOTFLUX_G            float32         nmgy Write me.
        FIBERTOTFLUX_R            float32         nmgy Write me.
        FIBERTOTFLUX_Z            float32         nmgy Write me.
                FLUX_G            float32         nmgy Write me.
           FLUX_IVAR_G            float32      / nmgy2 Write me.
                FLUX_R            float32         nmgy Write me.
           FLUX_IVAR_R            float32      / nmgy2 Write me.
                FLUX_Z            float32         nmgy Write me.
           FLUX_IVAR_Z            float32      / nmgy2 Write me.
               FLUX_W1            float32         nmgy Write me.
          FLUX_IVAR_W1            float32      / nmgy2 Write me.
               FLUX_W2            float32         nmgy Write me.
          FLUX_IVAR_W2            float32      / nmgy2 Write me.
               KCORR_U            float32          mag Write me.
              ABSMAG_U            float32          mag Write me.
         ABSMAG_IVAR_U            float32     1 / mag2 Write me.
               KCORR_G            float32          mag Write me.
              ABSMAG_G            float32          mag Write me.
         ABSMAG_IVAR_G            float32     1 / mag2 Write me.
               KCORR_R            float32          mag Write me.
              ABSMAG_R            float32          mag Write me.
         ABSMAG_IVAR_R            float32     1 / mag2 Write me.
               KCORR_I            float32          mag Write me.
              ABSMAG_I            float32          mag Write me.
         ABSMAG_IVAR_I            float32     1 / mag2 Write me.
               KCORR_Z            float32          mag Write me.
              ABSMAG_Z            float32          mag Write me.
         ABSMAG_IVAR_Z            float32     1 / mag2 Write me.
================================= =========== ======== ===========

Notes and Examples
==================

Upcoming changes
================
