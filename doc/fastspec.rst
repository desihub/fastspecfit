===================
fastspec Data Model
===================

:Summary: Spectroscopic fitting results.
:Naming Convention:
    ``fastspec-{petal}-{tileid}-{all,deep,night,exposures}.fits``, where
    ``{petal}`` is the petal or spetrograph number (0-9), ``{tileid}`` is the
    tileid, and the ``{all}``, ``{deep}``, ``{night}``, or ``{exposures}``
    suffix indicates which type of spectral coadds (and corresponding redshifts)
    were fitted (used).
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
NAXIS2   1225             int  length of dimension 2
SPECPROD blanc            str  spectroscopic production name
COADDTYP deep             str  spectral coadd fitted
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
          CONTINUUM_SNR_B     float32                               Median signal-to-noise ratio per pixel in *b* camera.
          CONTINUUM_SNR_R     float32                               Median signal-to-noise ratio per pixel in *r* camera.
          CONTINUUM_SNR_Z     float32                               Median signal-to-noise ratio per pixel in *z* camera.
   CONTINUUM_SMOOTHCORR_B     float32                               Mean value of the smooth continuum correction divided by the best-fitting continuum model in the *b* camera.
   CONTINUUM_SMOOTHCORR_R     float32                               Mean value of the smooth continuum correction divided by the best-fitting continuum model in the *r* camera.
   CONTINUUM_SMOOTHCORR_Z     float32                               Mean value of the smooth continuum correction divided by the best-fitting continuum model in the *z* camera.
                   DN4000     float32                               Narrow 4000-A break index (from Balogh et al. 1999) measured from the data.
              DN4000_IVAR     float32                               Inverse variance of DN4000.
             DN4000_MODEL     float32                               Narrow 4000-A break index (from Balogh et al. 1999) measured from the best-fitting continuum model.
             FLUX_SYNTH_G     float32                          nmgy g-band flux synthesized from the data.
             FLUX_SYNTH_R     float32                          nmgy r-band flux synthesized from the data.
             FLUX_SYNTH_Z     float32                          nmgy z-band flux synthesized from the data.
       FLUX_SYNTH_MODEL_G     float32                          nmgy g-band flux synthesized from the best-fitting continuum model.
       FLUX_SYNTH_MODEL_R     float32                          nmgy r-band flux synthesized from the best-fitting continuum model.
       FLUX_SYNTH_MODEL_Z     float32                          nmgy z-band flux synthesized from the best-fitting continuum model.
                 BALMER_Z     float32                        km / s Mean redshift of well-measured Balmer emission lines (defaults to the redrock redshift).
              FORBIDDEN_Z     float32                        km / s Mean redshift of well-measured forbidden emission lines (defaults to the redrock redshift).
                  BROAD_Z     float32                        km / s Mean redshift of well-measured broad emission lines (defaults to the redrock redshift).
             BALMER_SIGMA     float32                        km / s Mean line-width of well-measured Balmer emission lines.
          FORBIDDEN_SIGMA     float32                        km / s Mean line-width of well-measured forbidden emission lines.
              BROAD_SIGMA     float32                        km / s Mean line-width of well-measured broad emission lines.
              OI_1304_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
         OI_1304_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
             OI_1304_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
        OI_1304_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
          OI_1304_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
           OI_1304_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
            OI_1304_SIGMA     float32                        km / s Gaussian emission-line width.
             OI_1304_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
        OI_1304_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
               OI_1304_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
          OI_1304_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
       OI_1304_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
         OI_1304_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
             OI_1304_CHI2     float32                               Reduced chi^2 of the line-fit.
             OI_1304_NPIX       int32                               Number of pixels attributed to the emission line.
           SILIV_1396_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
      SILIV_1396_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
          SILIV_1396_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
     SILIV_1396_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
       SILIV_1396_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
        SILIV_1396_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
         SILIV_1396_SIGMA     float32                        km / s Gaussian emission-line width.
          SILIV_1396_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
     SILIV_1396_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
            SILIV_1396_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
       SILIV_1396_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
    SILIV_1396_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
      SILIV_1396_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
          SILIV_1396_CHI2     float32                               Reduced chi^2 of the line-fit.
          SILIV_1396_NPIX       int32                               Number of pixels attributed to the emission line.
             CIV_1549_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        CIV_1549_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            CIV_1549_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       CIV_1549_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         CIV_1549_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          CIV_1549_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           CIV_1549_SIGMA     float32                        km / s Gaussian emission-line width.
            CIV_1549_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       CIV_1549_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              CIV_1549_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         CIV_1549_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      CIV_1549_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        CIV_1549_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            CIV_1549_CHI2     float32                               Reduced chi^2 of the line-fit.
            CIV_1549_NPIX       int32                               Number of pixels attributed to the emission line.
          SILIII_1892_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
     SILIII_1892_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
         SILIII_1892_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
    SILIII_1892_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
      SILIII_1892_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       SILIII_1892_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
        SILIII_1892_SIGMA     float32                        km / s Gaussian emission-line width.
         SILIII_1892_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
    SILIII_1892_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
           SILIII_1892_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
      SILIII_1892_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
   SILIII_1892_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
     SILIII_1892_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
         SILIII_1892_CHI2     float32                               Reduced chi^2 of the line-fit.
         SILIII_1892_NPIX       int32                               Number of pixels attributed to the emission line.
            CIII_1908_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       CIII_1908_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           CIII_1908_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      CIII_1908_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        CIII_1908_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         CIII_1908_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          CIII_1908_SIGMA     float32                        km / s Gaussian emission-line width.
           CIII_1908_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      CIII_1908_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             CIII_1908_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        CIII_1908_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     CIII_1908_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       CIII_1908_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           CIII_1908_CHI2     float32                               Reduced chi^2 of the line-fit.
           CIII_1908_NPIX       int32                               Number of pixels attributed to the emission line.
            MGII_2796_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       MGII_2796_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           MGII_2796_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      MGII_2796_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        MGII_2796_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         MGII_2796_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          MGII_2796_SIGMA     float32                        km / s Gaussian emission-line width.
           MGII_2796_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      MGII_2796_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             MGII_2796_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        MGII_2796_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     MGII_2796_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       MGII_2796_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           MGII_2796_CHI2     float32                               Reduced chi^2 of the line-fit.
           MGII_2796_NPIX       int32                               Number of pixels attributed to the emission line.
            MGII_2803_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       MGII_2803_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           MGII_2803_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      MGII_2803_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        MGII_2803_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         MGII_2803_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          MGII_2803_SIGMA     float32                        km / s Gaussian emission-line width.
           MGII_2803_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      MGII_2803_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             MGII_2803_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        MGII_2803_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     MGII_2803_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       MGII_2803_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           MGII_2803_CHI2     float32                               Reduced chi^2 of the line-fit.
           MGII_2803_NPIX       int32                               Number of pixels attributed to the emission line.
             NEV_3346_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NEV_3346_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            NEV_3346_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NEV_3346_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         NEV_3346_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          NEV_3346_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           NEV_3346_SIGMA     float32                        km / s Gaussian emission-line width.
            NEV_3346_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NEV_3346_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              NEV_3346_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NEV_3346_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      NEV_3346_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NEV_3346_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NEV_3346_CHI2     float32                               Reduced chi^2 of the line-fit.
            NEV_3346_NPIX       int32                               Number of pixels attributed to the emission line.
             NEV_3426_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NEV_3426_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            NEV_3426_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NEV_3426_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         NEV_3426_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          NEV_3426_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           NEV_3426_SIGMA     float32                        km / s Gaussian emission-line width.
            NEV_3426_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NEV_3426_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              NEV_3426_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NEV_3426_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      NEV_3426_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NEV_3426_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NEV_3426_CHI2     float32                               Reduced chi^2 of the line-fit.
            NEV_3426_NPIX       int32                               Number of pixels attributed to the emission line.
             OII_3726_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        OII_3726_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            OII_3726_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       OII_3726_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         OII_3726_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          OII_3726_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           OII_3726_SIGMA     float32                        km / s Gaussian emission-line width.
            OII_3726_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       OII_3726_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              OII_3726_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         OII_3726_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      OII_3726_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        OII_3726_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            OII_3726_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
            OII_3726_NPIX       int32                               Number of pixels attributed to the emission line.
             OII_3729_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        OII_3729_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            OII_3729_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       OII_3729_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         OII_3729_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          OII_3729_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           OII_3729_SIGMA     float32                        km / s Gaussian emission-line width.
            OII_3729_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       OII_3729_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              OII_3729_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         OII_3729_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      OII_3729_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        OII_3729_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            OII_3729_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
            OII_3729_NPIX       int32                               Number of pixels attributed to the emission line.
           NEIII_3869_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
      NEIII_3869_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
          NEIII_3869_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
     NEIII_3869_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
       NEIII_3869_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
        NEIII_3869_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
         NEIII_3869_SIGMA     float32                        km / s Gaussian emission-line width.
          NEIII_3869_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
     NEIII_3869_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
            NEIII_3869_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
       NEIII_3869_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
    NEIII_3869_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
      NEIII_3869_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
          NEIII_3869_CHI2     float32                               Reduced chi^2 of the line-fit.
          NEIII_3869_NPIX       int32                               Number of pixels attributed to the emission line.
             HEI_3889_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        HEI_3889_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            HEI_3889_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       HEI_3889_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         HEI_3889_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          HEI_3889_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           HEI_3889_SIGMA     float32                        km / s Gaussian emission-line width.
            HEI_3889_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       HEI_3889_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              HEI_3889_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         HEI_3889_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      HEI_3889_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        HEI_3889_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            HEI_3889_CHI2     float32                               Reduced chi^2 of the line-fit.
            HEI_3889_NPIX       int32                               Number of pixels attributed to the emission line.
                   H6_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
              H6_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
                  H6_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
             H6_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
               H6_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
                H6_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
                 H6_SIGMA     float32                        km / s Gaussian emission-line width.
                  H6_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
             H6_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                    H6_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
               H6_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
            H6_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
              H6_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
                  H6_CHI2     float32                               Reduced chi^2 of the line-fit.
                  H6_NPIX       int32                               Number of pixels attributed to the emission line.
             HEPSILON_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        HEPSILON_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            HEPSILON_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       HEPSILON_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         HEPSILON_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          HEPSILON_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           HEPSILON_SIGMA     float32                        km / s Gaussian emission-line width.
            HEPSILON_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       HEPSILON_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              HEPSILON_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         HEPSILON_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      HEPSILON_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        HEPSILON_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            HEPSILON_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
            HEPSILON_NPIX       int32                               Number of pixels attributed to the emission line.
               HDELTA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HDELTA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              HDELTA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HDELTA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           HDELTA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
            HDELTA_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
             HDELTA_SIGMA     float32                        km / s Gaussian emission-line width.
              HDELTA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HDELTA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                HDELTA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
           HDELTA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
        HDELTA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HDELTA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HDELTA_CHI2     float32                               Reduced chi^2 of the line-fit.
              HDELTA_NPIX       int32                               Number of pixels attributed to the emission line.
               HGAMMA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HGAMMA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              HGAMMA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HGAMMA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           HGAMMA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
            HGAMMA_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
             HGAMMA_SIGMA     float32                        km / s Gaussian emission-line width.
              HGAMMA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HGAMMA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                HGAMMA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
           HGAMMA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
        HGAMMA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HGAMMA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HGAMMA_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
              HGAMMA_NPIX       int32                               Number of pixels attributed to the emission line.
            OIII_4363_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       OIII_4363_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           OIII_4363_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      OIII_4363_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        OIII_4363_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         OIII_4363_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          OIII_4363_SIGMA     float32                        km / s Gaussian emission-line width.
           OIII_4363_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      OIII_4363_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             OIII_4363_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        OIII_4363_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     OIII_4363_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       OIII_4363_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           OIII_4363_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
           OIII_4363_NPIX       int32                               Number of pixels attributed to the emission line.
             HEI_4471_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        HEI_4471_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            HEI_4471_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       HEI_4471_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         HEI_4471_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          HEI_4471_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           HEI_4471_SIGMA     float32                        km / s Gaussian emission-line width.
            HEI_4471_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       HEI_4471_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              HEI_4471_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         HEI_4471_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      HEI_4471_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        HEI_4471_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            HEI_4471_CHI2     float32                               Reduced chi^2 of the line-fit.
            HEI_4471_NPIX       int32                               Number of pixels attributed to the emission line.
            HEII_4686_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       HEII_4686_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           HEII_4686_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      HEII_4686_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        HEII_4686_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         HEII_4686_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          HEII_4686_SIGMA     float32                        km / s Gaussian emission-line width.
           HEII_4686_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      HEII_4686_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             HEII_4686_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        HEII_4686_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     HEII_4686_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       HEII_4686_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           HEII_4686_CHI2     float32                               Reduced chi^2 of the line-fit.
           HEII_4686_NPIX       int32                               Number of pixels attributed to the emission line.
                HBETA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           HBETA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
               HBETA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          HBETA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
            HBETA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
             HBETA_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
              HBETA_SIGMA     float32                        km / s Gaussian emission-line width.
               HBETA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          HBETA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                 HBETA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
            HBETA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
         HBETA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           HBETA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               HBETA_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
               HBETA_NPIX       int32                               Number of pixels attributed to the emission line.
            OIII_4959_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       OIII_4959_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           OIII_4959_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      OIII_4959_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        OIII_4959_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         OIII_4959_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          OIII_4959_SIGMA     float32                        km / s Gaussian emission-line width.
           OIII_4959_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      OIII_4959_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             OIII_4959_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        OIII_4959_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     OIII_4959_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       OIII_4959_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           OIII_4959_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
           OIII_4959_NPIX       int32                               Number of pixels attributed to the emission line.
            OIII_5007_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       OIII_5007_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           OIII_5007_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      OIII_5007_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        OIII_5007_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         OIII_5007_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          OIII_5007_SIGMA     float32                        km / s Gaussian emission-line width.
           OIII_5007_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      OIII_5007_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             OIII_5007_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        OIII_5007_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     OIII_5007_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       OIII_5007_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           OIII_5007_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
           OIII_5007_NPIX       int32                               Number of pixels attributed to the emission line.
             NII_5755_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NII_5755_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            NII_5755_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NII_5755_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         NII_5755_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          NII_5755_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           NII_5755_SIGMA     float32                        km / s Gaussian emission-line width.
            NII_5755_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NII_5755_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              NII_5755_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NII_5755_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      NII_5755_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NII_5755_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NII_5755_CHI2     float32                               Reduced chi^2 of the line-fit.
            NII_5755_NPIX       int32                               Number of pixels attributed to the emission line.
             HEI_5876_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        HEI_5876_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            HEI_5876_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       HEI_5876_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         HEI_5876_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          HEI_5876_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           HEI_5876_SIGMA     float32                        km / s Gaussian emission-line width.
            HEI_5876_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       HEI_5876_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              HEI_5876_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         HEI_5876_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      HEI_5876_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        HEI_5876_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            HEI_5876_CHI2     float32                               Reduced chi^2 of the line-fit.
            HEI_5876_NPIX       int32                               Number of pixels attributed to the emission line.
              OI_6300_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
         OI_6300_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
             OI_6300_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
        OI_6300_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
          OI_6300_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
           OI_6300_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
            OI_6300_SIGMA     float32                        km / s Gaussian emission-line width.
             OI_6300_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
        OI_6300_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
               OI_6300_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
          OI_6300_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
       OI_6300_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
         OI_6300_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
             OI_6300_CHI2     float32                               Reduced chi^2 of the line-fit.
             OI_6300_NPIX       int32                               Number of pixels attributed to the emission line.
             NII_6548_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NII_6548_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            NII_6548_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NII_6548_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         NII_6548_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          NII_6548_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           NII_6548_SIGMA     float32                        km / s Gaussian emission-line width.
            NII_6548_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NII_6548_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              NII_6548_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NII_6548_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      NII_6548_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NII_6548_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NII_6548_CHI2     float32                               Reduced chi^2 of the line-fit.
            NII_6548_NPIX       int32                               Number of pixels attributed to the emission line.
               HALPHA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HALPHA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
              HALPHA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HALPHA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
           HALPHA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
            HALPHA_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
             HALPHA_SIGMA     float32                        km / s Gaussian emission-line width.
              HALPHA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HALPHA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
                HALPHA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
           HALPHA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
        HALPHA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HALPHA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HALPHA_CHI2     float32                               Reduced chi^2 of the line-fit (default value 1e6).
              HALPHA_NPIX       int32                               Number of pixels attributed to the emission line.
             NII_6584_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NII_6584_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            NII_6584_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NII_6584_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         NII_6584_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          NII_6584_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           NII_6584_SIGMA     float32                        km / s Gaussian emission-line width.
            NII_6584_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NII_6584_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              NII_6584_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NII_6584_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      NII_6584_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NII_6584_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NII_6584_CHI2     float32                               Reduced chi^2 of the line-fit.
            NII_6584_NPIX       int32                               Number of pixels attributed to the emission line.
             SII_6716_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        SII_6716_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            SII_6716_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       SII_6716_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         SII_6716_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          SII_6716_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           SII_6716_SIGMA     float32                        km / s Gaussian emission-line width.
            SII_6716_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       SII_6716_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              SII_6716_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         SII_6716_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      SII_6716_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        SII_6716_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            SII_6716_CHI2     float32                               Reduced chi^2 of the line-fit.
            SII_6716_NPIX       int32                               Number of pixels attributed to the emission line.
             SII_6731_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        SII_6731_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            SII_6731_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       SII_6731_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         SII_6731_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          SII_6731_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           SII_6731_SIGMA     float32                        km / s Gaussian emission-line width.
            SII_6731_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       SII_6731_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              SII_6731_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         SII_6731_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      SII_6731_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        SII_6731_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            SII_6731_CHI2     float32                               Reduced chi^2 of the line-fit.
            SII_6731_NPIX       int32                               Number of pixels attributed to the emission line.
             OII_7320_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        OII_7320_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            OII_7320_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       OII_7320_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         OII_7320_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          OII_7320_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           OII_7320_SIGMA     float32                        km / s Gaussian emission-line width.
            OII_7320_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       OII_7320_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              OII_7320_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         OII_7320_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      OII_7320_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        OII_7320_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            OII_7320_CHI2     float32                               Reduced chi^2 of the line-fit.
            OII_7320_NPIX       int32                               Number of pixels attributed to the emission line.
             OII_7330_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        OII_7330_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
            OII_7330_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       OII_7330_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
         OII_7330_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
          OII_7330_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
           OII_7330_SIGMA     float32                        km / s Gaussian emission-line width.
            OII_7330_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       OII_7330_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
              OII_7330_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         OII_7330_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
      OII_7330_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        OII_7330_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            OII_7330_CHI2     float32                               Reduced chi^2 of the line-fit.
            OII_7330_NPIX       int32                               Number of pixels attributed to the emission line.
            SIII_9069_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       SIII_9069_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           SIII_9069_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      SIII_9069_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        SIII_9069_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         SIII_9069_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          SIII_9069_SIGMA     float32                        km / s Gaussian emission-line width.
           SIII_9069_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      SIII_9069_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             SIII_9069_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        SIII_9069_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     SIII_9069_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       SIII_9069_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           SIII_9069_CHI2     float32                               Reduced chi^2 of the line-fit.
           SIII_9069_NPIX       int32                               Number of pixels attributed to the emission line.
            SIII_9532_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       SIII_9532_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of line-amplitude.
           SIII_9532_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      SIII_9532_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of integrated flux.
        SIII_9532_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
         SIII_9532_VSHIFT     float32                        km / s Velocity shift relative to the redrock redshift.
          SIII_9532_SIGMA     float32                        km / s Gaussian emission-line width.
           SIII_9532_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      SIII_9532_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of continuum flux.
             SIII_9532_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        SIII_9532_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of equivalent width.
     SIII_9532_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       SIII_9532_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           SIII_9532_CHI2     float32                               Reduced chi^2 of the line-fit.
           SIII_9532_NPIX       int32                               Number of pixels attributed to the emission line.
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
                SURVEY  bytes3                Survey name (e.g., 'SV3'); only present when fitting healpix coadds.
               FAPRGRM  bytes6                Program name (*bright* or *dark*); only present when fitting healpix coadds.
              HPXPIXEL   int32                Healpixel number (nside=64); only present when fitting healpix coadds.
                TILEID   int32                Tile ID number; only present when fitting tile-level (not healpix) coadds.
                 FIBER   int32                Fiber ID number; only present with TILEID.
             THRUNIGHT   int32                Last night of coadded data; only present when fitting cumulative coadds.
                 NIGHT   int32                Night; only present when fitting per-night or per-exposure spectra.
                 EXPID   int32                Exposure ID number; only present when fitting per-exposure spectra.
               PHOTSYS  bytes1                Photometric system (*N* or *S*).
     MW_TRANSMISSION_G float32                Milky Way foreground dust transmission factor [0-1] in the g-band.
     MW_TRANSMISSION_R float32                Milky Way foreground dust transmission factor [0-1] in the r-band.
     MW_TRANSMISSION_Z float32                Milky Way foreground dust transmission factor [0-1] in the z-band.
    MW_TRANSMISSION_W1 float32                Milky Way foreground dust transmission factor [0-1] in the W1-band.
    MW_TRANSMISSION_W2 float32                Milky Way foreground dust transmission factor [0-1] in the W2-band.
           DESI_TARGET   int64                DESI targeting bit.
            BGS_TARGET   int64                BGS targeting bit.
            MWS_TARGET   int64                MWS targeting bit.
       SV1_DESI_TARGET   int64                SV1 DESI targeting bit.
        SV1_BGS_TARGET   int64                SV1 BGS targeting bit.
        SV1_MWS_TARGET   int64                SV1 MWS targeting bit.
       SV2_DESI_TARGET   int64                SV2 DESI targeting bit.
        SV2_BGS_TARGET   int64                SV2 BGS targeting bit.
        SV2_MWS_TARGET   int64                SV2 MWS targeting bit.
       SV3_DESI_TARGET   int64                SV3 DESI targeting bit.
        SV3_BGS_TARGET   int64                SV3 BGS targeting bit.
        SV3_MWS_TARGET   int64                SV3 MWS targeting bit.
           SCND_TARGET   int64                Secondary target targeting bit.
       SV1_SCND_TARGET   int64                SV1 secondary targeting bit.
       SV2_SCND_TARGET   int64                SV2 secondary targeting bit.
       SV3_SCND_TARGET   int64                SV3 secondary targeting bit.
                     Z float64                Redrock redshift.
                 ZWARN    int8                Redrock zwarning bit.
             DELTACHI2 float64                Redrock delta-chi-squared.
              SPECTYPE    str6                Redrock spectral classification.
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
           FLUX_IVAR_G float32     1 / nmgy^2 Inverse variance of FLUX_G from targeting catalog.
           FLUX_IVAR_R float32     1 / nmgy^2 Inverse variance of FLUX_R from targeting catalog.
           FLUX_IVAR_Z float32     1 / nmgy^2 Inverse variance of FLUX_Z from targeting catalog.
          FLUX_IVAR_W1 float32     1 / nmgy^2 Inverse variance of FLUX_W1 from targeting catalog.
          FLUX_IVAR_W2 float32     1 / nmgy^2 Inverse variance of FLUX_W2 from targeting catalog.
====================== =========== ========== ==========================================

Notes and Examples
==================


Upcoming changes
================
