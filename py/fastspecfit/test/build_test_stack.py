"""
make_lrg_test_stack.py
======================
Build a single stacked LRG rest-frame spectrum from five sv3/dark targets
for use as a unit-test fixture in the desigal test suite.

The five targets span z ~ 0.41–0.45 and were selected by hand to be
representative LRGs with reliable fastspecfit solutions.  Spectra are
read from the ``loa`` spectroscopic production, shifted to the rest frame,
normalised in a 50 Å window centred on 4500 Å (rest), and combined with
inverse-variance weighting into a single spectrum written via
``write_binned_stacks``.

Usage
-----
Run directly::

    python build_test_stack

Output
------
``stack-LRG.fits`` in the current working directory.

"""
import numpy as np
from astropy.table import Table
from desigal.specutils.stack import stack_spectra, write_binned_stacks

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
targetid = np.array([
    39627763635720044,
    39627829507262695,
    39627921001810574,
    39633290050669693,
    39627775727900698,
], dtype=np.int64)
healpix = np.array([
    26275,
    27247,
    26091,
    11425,
    26276,
], dtype=np.int32)
zz = np.array([
    0.4052904043347254,
    0.411587106145884,
    0.4462853971996436,
    0.44740186572631313,
    0.4263606540076314,
])
targets = Table({
    "TARGETID": targetid,
    "SURVEY":  ["sv3"] * len(targetid),
    "PROGRAM": ["dark"] * len(targetid),
    "HEALPIX": healpix,
    "Z": zz,
})

# ---------------------------------------------------------------------------
# Rest-frame wavelength grid (same conventions as build_stacks)
# ---------------------------------------------------------------------------
SPECPROD       = "loa"
NORM_WINDOW    = [4475.0, 4525.0]   # LRG: centred on ~4500 Å rest
DWAVE          = 0.4                # Å/pixel
OBSWAVE_MIN    = 3600.0
OBSWAVE_MAX    = 9800.0

zmin = targets["Z"].min()
zmax = targets["Z"].max()
restwave = np.arange(OBSWAVE_MIN / (1.0 + zmax), OBSWAVE_MAX / (1.0 + zmin), DWAVE)

# ---------------------------------------------------------------------------
# Read spectra
# ---------------------------------------------------------------------------
from desispec.io.spectra import read_spectra_parallel

spectra = read_spectra_parallel(
    targets,
    prefix="coadd",
    specprod=SPECPROD,
    nproc=1,
    rdspec_kwargs={
        "skip_hdus": ["EXP_FIBERMAP", "SCORES", "EXTRA_CATALOG", "RESOLUTION"]
    },
)
assert np.all(spectra.target_ids() == targets["TARGETID"]), \
    "Spectra order does not match targets table."

# ---------------------------------------------------------------------------
# Stack
# ---------------------------------------------------------------------------
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    (flux, ivar), _ = stack_spectra(
        spectra=spectra,
        redshift=targets["Z"].data,
        stack_redshift=0.0,
        norm_flux_window=NORM_WINDOW,
        norm_method="flux-window",
        resample_method="linear",
        stack_method="ivar-weighted-mean",
        output_wave_grid=restwave,
        n_workers=1,
    )

# ---------------------------------------------------------------------------
# Pack metadata and write
# ---------------------------------------------------------------------------
outmeta = Table(
    {
        "IBIN":    np.array([0], dtype=np.int32),
        "NALLOBJ": np.array([len(targets)], dtype=np.int32),
        "NOBJ":    np.array([len(targets)], dtype=np.int32),
        "STACKID": np.array([0], dtype=np.int32),
        "Z":       np.array([0.0]),          # rest-frame stack; always 0
        "ZMED":    np.array([np.median(targets["Z"])], dtype=np.float32),
    }
)

igood = np.where(np.isfinite(flux) & np.isfinite(ivar))[0]
snr = float(np.median(flux[igood] * np.sqrt(ivar[igood])))
outmeta["SNR"] = np.array([snr], dtype=np.float32)

outfile = "stack-LRG.fits"
write_binned_stacks(
    outfile,
    restwave,
    flux[np.newaxis, :],   # write_binned_stacks expects (nbins, npix)
    ivar[np.newaxis, :],
    resolution=None,
    stackinfo=outmeta,
)
print(f"Wrote {outfile}  (SNR={snr:.2f}, {len(igood)} good pixels)")
