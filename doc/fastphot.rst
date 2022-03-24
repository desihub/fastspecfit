===================
fastphot Data Model
===================

:Summary: Photometric fitting results.
:Naming Convention:
    ``fastphot-{petal}-{tileid}-{all,deep,night,exposures}.fits``, where
    ``{petal}`` is the petal or spetrograph number (0-9), ``{tileid}`` is the
    tileid, and the ``{all}``, ``{deep}``, ``{night}``, or ``{exposures}``
    suffix indicates which type of spectral coadds (and corresponding redshifts)
    were fitted (used).
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

Fitting results.

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
              TARGETID int64                  Unique target ID.
       CONTINUUM_COEFF float64[11]            Continuum coefficients.
        CONTINUUM_CHI2 float32                Reduced chi^2 of the continuum fit.
         CONTINUUM_AGE float32            Gyr Light-weighted age.
          CONTINUUM_AV float32            mag Intrinsic attenuation.
     CONTINUUM_AV_IVAR float32     1 / mag^2  Inverse variance of CONTINUUM_AV.
          DN4000_MODEL float32                Narrow 4000-A break index (from Balogh et al. 1999) measured from the best-fitting continuum model.
    FLUX_SYNTH_MODEL_G float32           nmgy g-band flux synthesized from the best-fitting continuum model.
    FLUX_SYNTH_MODEL_R float32           nmgy r-band flux synthesized from the best-fitting continuum model.
    FLUX_SYNTH_MODEL_Z float32           nmgy z-band flux synthesized from the best-fitting continuum model.
   FLUX_SYNTH_MODEL_W1 float32           nmgy W1-band flux synthesized from the best-fitting continuum model.
   FLUX_SYNTH_MODEL_W2 float32           nmgy W2-band flux synthesized from the best-fitting continuum model.
               KCORR_U float32            mag K-correction used to derive ABSMAG_U.
              ABSMAG_U float32            mag Absolute magnitude in Johnson/Cousins U-band.
         ABSMAG_IVAR_U float32      1 / mag^2 Inverse variance corresponding to ABSMAG_U.
               KCORR_B float32            mag Like KCORR_U but for Johnson/Cousins B-band.
              ABSMAG_B float32            mag Like ABSMAG_U but for Johnson/Cousins B-band.
         ABSMAG_IVAR_B float32      1 / mag^2 Like ABSMAG_IVAR_U but for Johnson/Cousins B-band.
               KCORR_V float32            mag Like KCORR_U but for Johnson/Cousins V-band.
              ABSMAG_V float32            mag Like ABSMAG_U but for Johnson/Cousins V-band.
         ABSMAG_IVAR_V float32      1 / mag^2 Like ABSMAG_IVAR_U but for Johnson/Cousins V-band.
          KCORR_SDSS_U float32            mag Like KCORR_U but for SDSS u-band.
         ABSMAG_SDSS_U float32            mag Like ABSMAG_U but for SDSS u-band.
    ABSMAG_IVAR_SDSS_U float32      1 / mag^2 Like ABSMAG_IVAR_U but for SDSS u-band.
          KCORR_SDSS_G float32            mag Like KCORR_U but for SDSS g-band.
         ABSMAG_SDSS_G float32            mag Like ABSMAG_U but for SDSS g-band.
    ABSMAG_IVAR_SDSS_G float32      1 / mag^2 Like ABSMAG_IVAR_U but for SDSS g-band.
          KCORR_SDSS_R float32            mag Like KCORR_U but for SDSS r-band.
         ABSMAG_SDSS_R float32            mag Like ABSMAG_U but for SDSS r-band.
    ABSMAG_IVAR_SDSS_R float32      1 / mag^2 Like ABSMAG_IVAR_U but for SDSS r-band.
          KCORR_SDSS_I float32            mag Like KCORR_U but for SDSS i-band.
         ABSMAG_SDSS_I float32            mag Like ABSMAG_U but for SDSS i-band.
    ABSMAG_IVAR_SDSS_I float32      1 / mag^2 Like ABSMAG_IVAR_U but for SDSS i-band.
          KCORR_SDSS_Z float32            mag Like KCORR_U but for SDSS z-band.
         ABSMAG_SDSS_Z float32            mag Like ABSMAG_U but for SDSS z-band.
    ABSMAG_IVAR_SDSS_Z float32      1 / mag^2 Like ABSMAG_IVAR_U but for SDSS z-band.
              KCORR_W1 float32            mag Like KCORR_U but for WISE W1-band.
             ABSMAG_W1 float32            mag Like ABSMAG_U but for WISE W1-band.
        ABSMAG_IVAR_W1 float32      1 / mag^2 Like ABSMAG_IVAR_U but for WISE W1-band.
====================== =========== ========== ==========================================

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
SPECPROD fuji             str  spectroscopic production name
COADDTYP healpix          str  type of spectral coadd (*healpix*, *cumulative*, *pernight*, *perexp*)
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
                SURVEY  bytes3                Survey name (e.g., 'SV3'); only present when fitting healpix coadds.
               PROGRAM  bytes6                Program name (e.g., *bright*, *dark*); only present when fitting healpix coadds.
               HEALPIX   int32                Healpixel number (nside=64, nested=True); only present when fitting healpix coadds.
           TILEID_LIST    str5                List of tile IDs that went into healpix coadd.
                TILEID   int32                Tile ID number; only present when fitting tile-level (not healpix) coadds.
                 FIBER   int32                Fiber ID number; only present with TILEID.
                 NIGHT   int32                Night (or *thrunight* for cumulative coadds); only present with TILEID.
                 EXPID   int32                Exposure ID number; only present with TILEID and when fitting per-exposure spectra.
           DESI_TARGET   int64                DESI targeting bit.
            BGS_TARGET   int64                BGS targeting bit.
            MWS_TARGET   int64                MWS targeting bit.
           SCND_TARGET   int64                Secondary target targeting bit.
       SV1_DESI_TARGET   int64                SV1 DESI targeting bit; only present in fuji / Early Data Release.
        SV1_BGS_TARGET   int64                SV1 BGS targeting bit; only present in fuji / Early Data Release.
        SV1_MWS_TARGET   int64                SV1 MWS targeting bit; only present in fuji / Early Data Release.
       SV2_DESI_TARGET   int64                SV2 DESI targeting bit; only present in fuji / Early Data Release.
        SV2_BGS_TARGET   int64                SV2 BGS targeting bit; only present in fuji / Early Data Release.
        SV2_MWS_TARGET   int64                SV2 MWS targeting bit; only present in fuji / Early Data Release.
       SV3_DESI_TARGET   int64                SV3 DESI targeting bit; only present in fuji / Early Data Release.
        SV3_BGS_TARGET   int64                SV3 BGS targeting bit; only present in fuji / Early Data Release.
        SV3_MWS_TARGET   int64                SV3 MWS targeting bit; only present in fuji / Early Data Release.
       SV1_SCND_TARGET   int64                SV1 secondary targeting bit; only present in fuji / Early Data Release.
       SV2_SCND_TARGET   int64                SV2 secondary targeting bit; only present in fuji / Early Data Release.
       SV3_SCND_TARGET   int64                SV3 secondary targeting bit; only present in fuji / Early Data Release.
                     Z float64                Redrock redshift.
                 ZWARN   int64                Redrock zwarning bit.
             DELTACHI2 float64                Redrock delta-chi-squared.
              SPECTYPE  bytes6                Redrock spectral classification.
               PHOTSYS  bytes1                Photometric system (*N* or *S*).
     MW_TRANSMISSION_G float32                Milky Way foreground dust transmission factor [0-1] in the g-band.
     MW_TRANSMISSION_R float32                Milky Way foreground dust transmission factor [0-1] in the r-band.
     MW_TRANSMISSION_Z float32                Milky Way foreground dust transmission factor [0-1] in the z-band.
    MW_TRANSMISSION_W1 float32                Milky Way foreground dust transmission factor [0-1] in the W1-band.
    MW_TRANSMISSION_W2 float32                Milky Way foreground dust transmission factor [0-1] in the W2-band.
    MW_TRANSMISSION_W3 float32                Milky Way foreground dust transmission factor [0-1] in the W3-band.
    MW_TRANSMISSION_W4 float32                Milky Way foreground dust transmission factor [0-1] in the W4-band.
           FIBERFLUX_G float32           nmgy Fiber g-band flux from targeting catalog.
           FIBERFLUX_R float32           nmgy Fiber r-band flux from targeting catalog.
           FIBERFLUX_Z float32           nmgy Fiber z-band flux from targeting catalog.
        FIBERTOTFLUX_G float32           nmgy Fibertot g-band flux from targeting catalog.
        FIBERTOTFLUX_R float32           nmgy Fibertot r-band flux from targeting catalog.
        FIBERTOTFLUX_Z float32           nmgy Fibertot z-band flux from targeting catalog.
                FLUX_G float32           nmgy Total g-band flux from targeting catalog.
                FLUX_R float32           nmgy Total r-band flux from targeting catalog.
                FLUX_Z float32           nmgy Total z-band flux from targeting catalog.
               FLUX_W1 float32           nmgy Total W1-band flux from targeting catalog.
               FLUX_W2 float32           nmgy Total W2-band flux from targeting catalog.
               FLUX_W3 float32           nmgy Total W3-band flux from targeting catalog.
               FLUX_W4 float32           nmgy Total W4-band flux from targeting catalog.
           FLUX_IVAR_G float32     1 / nmgy^2 Inverse variance of FLUX_G from targeting catalog.
           FLUX_IVAR_R float32     1 / nmgy^2 Inverse variance of FLUX_R from targeting catalog.
           FLUX_IVAR_Z float32     1 / nmgy^2 Inverse variance of FLUX_Z from targeting catalog.
          FLUX_IVAR_W1 float32     1 / nmgy^2 Inverse variance of FLUX_W1 from targeting catalog.
          FLUX_IVAR_W2 float32     1 / nmgy^2 Inverse variance of FLUX_W2 from targeting catalog.
          FLUX_IVAR_W3 float32     1 / nmgy^2 Inverse variance of FLUX_W3 from targeting catalog.
          FLUX_IVAR_W4 float32     1 / nmgy^2 Inverse variance of FLUX_W4 from targeting catalog.
====================== =========== ========== ==========================================

Notes and Examples
==================

If the inverse variance on a given absolutely magnitude is zero it means that
the absolute magnitude was derived from *synthesized* photometry based on the
best-fitting model (i.e., use with care).

Similarly, if CONTINUUM_AV_IVAR is zero it means that fitted for the (intrinsic)
dust extinction failed.

In general, one should use the value of CONTINUUM_CHI2 to assess the quality of
the fit to the broadband photometry.

Upcoming changes
================

A basic stellar mass estimate will be added.

