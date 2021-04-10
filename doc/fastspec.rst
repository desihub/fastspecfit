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
            OIII_4959_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       OIII_4959_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OIII_4959_AMP.
           OIII_4959_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      OIII_4959_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OIII_4959_FLUX
        OIII_4959_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
   OIII_4959_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OIII_4959_BOXFLUX
           OIII_4959_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      OIII_4959_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OIII_4959_CONT
             OIII_4959_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        OIII_4959_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of OIII_4959_EW
     OIII_4959_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       OIII_4959_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           OIII_4959_CHI2     float32                               Reduced chi^2 of the line-fit.
           OIII_4959_NPIX       int32                               Number of pixels attributed to the emission line.
            OIII_5007_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
       OIII_5007_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OIII_5007_AMP.
           OIII_5007_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
      OIII_5007_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OIII_5007_FLUX
        OIII_5007_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
   OIII_5007_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of OIII_5007_BOXFLUX
           OIII_5007_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
      OIII_5007_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of OIII_5007_CONT
             OIII_5007_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
        OIII_5007_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of OIII_5007_EW
     OIII_5007_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
       OIII_5007_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
           OIII_5007_CHI2     float32                               Reduced chi^2 of the line-fit.
           OIII_5007_NPIX       int32                               Number of pixels attributed to the emission line.
             NII_6548_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NII_6548_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of NII_6548_AMP.
            NII_6548_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NII_6548_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of NII_6548_FLUX
         NII_6548_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    NII_6548_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of NII_6548_BOXFLUX
            NII_6548_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NII_6548_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of NII_6548_CONT
              NII_6548_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NII_6548_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of NII_6548_EW
      NII_6548_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NII_6548_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NII_6548_CHI2     float32                               Reduced chi^2 of the line-fit.
            NII_6548_NPIX       int32                               Number of pixels attributed to the emission line.
             NII_6584_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        NII_6584_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of NII_6584_AMP.
            NII_6584_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       NII_6584_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of NII_6584_FLUX
         NII_6584_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    NII_6584_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of NII_6584_BOXFLUX
            NII_6584_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       NII_6584_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of NII_6584_CONT
              NII_6584_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         NII_6584_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of NII_6584_EW
      NII_6584_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        NII_6584_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            NII_6584_CHI2     float32                               Reduced chi^2 of the line-fit.
            NII_6584_NPIX       int32                               Number of pixels attributed to the emission line.
             SII_6716_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        SII_6716_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of SII_6716_AMP.
            SII_6716_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       SII_6716_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of SII_6716_FLUX
         SII_6716_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    SII_6716_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of SII_6716_BOXFLUX
            SII_6716_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       SII_6716_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of SII_6716_CONT
              SII_6716_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         SII_6716_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of SII_6716_EW
      SII_6716_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        SII_6716_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            SII_6716_CHI2     float32                               Reduced chi^2 of the line-fit.
            SII_6716_NPIX       int32                               Number of pixels attributed to the emission line.
             SII_6731_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        SII_6731_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of SII_6731_AMP.
            SII_6731_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       SII_6731_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of SII_6731_FLUX
         SII_6731_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    SII_6731_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of SII_6731_BOXFLUX
            SII_6731_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       SII_6731_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of SII_6731_CONT
              SII_6731_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         SII_6731_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of SII_6731_EW
      SII_6731_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        SII_6731_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            SII_6731_CHI2     float32                               Reduced chi^2 of the line-fit.
            SII_6731_NPIX       int32                               Number of pixels attributed to the emission line.
             HEPSILON_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
        HEPSILON_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HEPSILON_AMP.
            HEPSILON_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
       HEPSILON_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HEPSILON_FLUX
         HEPSILON_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
    HEPSILON_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HEPSILON_BOXFLUX
            HEPSILON_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
       HEPSILON_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HEPSILON_CONT
              HEPSILON_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
         HEPSILON_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of HEPSILON_EW
      HEPSILON_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
        HEPSILON_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
            HEPSILON_CHI2     float32                               Reduced chi^2 of the line-fit.
            HEPSILON_NPIX       int32                               Number of pixels attributed to the emission line.
               HGAMMA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HGAMMA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HGAMMA_AMP.
              HGAMMA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HGAMMA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HGAMMA_FLUX
           HGAMMA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      HGAMMA_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HGAMMA_BOXFLUX
              HGAMMA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HGAMMA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HGAMMA_CONT
                HGAMMA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
           HGAMMA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of HGAMMA_EW
        HGAMMA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HGAMMA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HGAMMA_CHI2     float32                               Reduced chi^2 of the line-fit.
              HGAMMA_NPIX       int32                               Number of pixels attributed to the emission line.
               HGAMMA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HGAMMA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HGAMMA_AMP.
              HGAMMA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HGAMMA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HGAMMA_FLUX
           HGAMMA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      HGAMMA_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HGAMMA_BOXFLUX
              HGAMMA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HGAMMA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HGAMMA_CONT
                HGAMMA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
           HGAMMA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of HGAMMA_EW
        HGAMMA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HGAMMA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HGAMMA_CHI2     float32                               Reduced chi^2 of the line-fit.
              HGAMMA_NPIX       int32                               Number of pixels attributed to the emission line.
                HBETA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
           HBETA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HBETA_AMP.
               HBETA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
          HBETA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HBETA_FLUX
            HBETA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
       HBETA_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HBETA_BOXFLUX
               HBETA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
          HBETA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HBETA_CONT
                 HBETA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
            HBETA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of HBETA_EW
         HBETA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
           HBETA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
               HBETA_CHI2     float32                               Reduced chi^2 of the line-fit.
               HBETA_NPIX       int32                               Number of pixels attributed to the emission line.
               HALPHA_AMP     float32  1e-17 erg / (Angstrom cm2 s) Emission line amplitude.
          HALPHA_AMP_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HALPHA_AMP.
              HALPHA_FLUX     float32           1e-17 erg / (cm2 s) Gaussian-integrated emission-line flux.
         HALPHA_FLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HALPHA_FLUX
           HALPHA_BOXFLUX     float32           1e-17 erg / (cm2 s) Boxcar-integrated emission-line flux.
      HALPHA_BOXFLUX_IVAR     float32           1e+34 cm4 s2 / erg2 Inverse variance of HALPHA_BOXFLUX
              HALPHA_CONT     float32  1e-17 erg / (Angstrom cm2 s) Continuum flux at line center.
         HALPHA_CONT_IVAR     float32 1e+34 Angstrom2 cm4 s2 / erg2 Inverse variance of HALPHA_CONT
                HALPHA_EW     float32                      Angstrom Rest-frame emission-line equivalent width.
           HALPHA_EW_IVAR     float32                 1 / Angstrom2 Inverse variance of HALPHA_EW
        HALPHA_FLUX_LIMIT     float32                 erg / (cm2 s) One-sigma upper limit on the emission line flux.
          HALPHA_EW_LIMIT     float32                      Angstrom One-sigma upper limit on the emission line equivalent width.
              HALPHA_CHI2     float32                               Reduced chi^2 of the line-fit.
              HALPHA_NPIX       int32                               Number of pixels attributed to the emission line.
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
                 ZWARN    int8                Redrock zwarning bit.
             DELTACHI2 float64                Redrock delta-chi-squared.
              SPECTYPE    str6                Redrock spectral classification.
====================== =========== ========== ==========================================

Notes and Examples
==================


Upcoming changes
================
