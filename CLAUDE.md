# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

FastSpecFit performs fast spectral energy distribution (SED) modeling and emission-line fitting of [DESI](https://desi.lbl.gov) galaxy and QSO/AGN spectra and broadband photometry. It is a production astronomy pipeline used by DESI collaborators to generate value-added catalogs.

## Commands

### Installation
```bash
pip install -e ".[test]"          # development install with test deps
pip install -e ".[coverage]"      # with coverage tools
pip install -e ".[doc]"           # with Sphinx documentation tools
```

### Running tests
```bash
pytest                            # run all tests
pytest py/fastspecfit/test/test_fastspecfit.py  # single test file
pytest py/fastspecfit/test/test_fastspecfit.py::test_fastphot  # single test
pytest --cov                      # with coverage
```

Note: The main integration tests (`test_fastphot`, `test_fastspec`) download template files (~100 MB) from `data.desi.lbl.gov` on first run and require the `FTEMPLATES_DIR`-equivalent path to be available.

### Building docs
```bash
sphinx-build -W --keep-going -b html doc doc/_build/html
```

Technical notes (LaTeX) are in `doc/technote/`. Build with `make` in that directory.

### Style check
CI runs `flake8` with custom ignores defined in `pyproject.toml` under `[tool.flake8]`.

## Required Environment Variables

These must be set before running `fastspec` or `fastphot` (unless the corresponding CLI flags override them):

| Variable | Purpose |
|---|---|
| `DESI_SPECTRO_REDUX` | Root of DESI spectroscopic reduction outputs |
| `DUST_DIR` | Schlegel, Finkbeiner & Davis dust maps |
| `FPHOTO_DIR` | Legacy Survey DR9 broadband photometry |
| `FTEMPLATES_DIR` | Stellar population synthesis template files |

## CLI Entry Points

Defined in `pyproject.toml` and implemented in `py/fastspecfit/fastspecfit.py`:

- `fastspec` — full spectrophotometric fitting (continuum + emission lines)
- `fastphot` — photometry-only continuum fitting
- `fastqso` — QSO/AGN spectrophotometric fitting (power law + UV Fe + IR torus)
- `stackfit` — fit stacked/coadded spectra
- `fastqa` — generate QA figures from output catalogs
- `mpi-fastspecfit` (bin script) — MPI parallel execution across many files
- `profile-fastspec` (bin script) — profiling tool; supports `--mode fastqso`

## Architecture

### Data Flow

1. **Input**: Redrock redshift catalogs + DESI coadded spectra (FITS files)
2. **`fastspec_one()`** in `fastspecfit.py` — processes a single object:
   - Calls `continuum_specfit()` to fit the stellar continuum and photometry
   - Calls `emline_specfit()` to fit emission lines on the residual spectrum
3. **Output**: Multi-extension FITS with `METADATA`, `SPECPHOT`, `FASTSPEC`/`FASTQSO`, and `MODELS` HDUs

**`fastqso` mode** (`continuum.qso_continuum_fastspec`): fits QSOs/AGN with a three-component model — power law (F_λ ∝ (λ/1450Å)^(α−2) with IGM), UV Fe emission (Vestergaard & Wilkes 2001 template), and IR AGN torus (Nenkova+08, amplitude from energy balance). Dust attenuation uses a steeper power-law extinction curve (slope −1.2 vs. −0.7 for galaxies). Parameters: PL_SLOPE (α), PL_AMPLITUDE, FE_VDISP (σ_Fe), FE_AMPLITUDE, TAUV (τ_V), TORUS_AMPLITUDE. Fitting uses VARPRO: outer grid on α, inner Brent minimization on τ_V, NNLS for amplitudes; all uncertainties from Monte Carlo.

### Single-Copy Global State (`singlecopy.py`)

The `sc_data` singleton (a `Singletons` instance) is initialized once per process and shared across all worker threads/processes. It holds:
- `sc_data.templates` — stellar population synthesis templates (`Templates`)
- `sc_data.emlines` — emission line table (`LineTable`)
- `sc_data.photometry` — filter curves and photometric parameters (`Photometry`)
- `sc_data.cosmology` — tabulated DESI fiducial cosmology (`TabulatedDESI`)
- `sc_data.igm` — Inoue+14 IGM attenuation model (`Inoue14`)

This pattern avoids re-reading large files in multiprocessing workers.

### Key Modules

| Module | Role |
|---|---|
| `fastspecfit.py` | CLI parsing, top-level `fastspec`/`fastphot` drivers, per-object dispatch |
| `continuum.py` | `ContinuumTools` class — stellar SED fitting against SPS templates with dust attenuation and velocity dispersion |
| `emlines.py` | `EMFitTools` class — emission-line fitting orchestration |
| `emline_fit/` | Numba-accelerated emission-line model (Gaussian point evaluation at pixel centers), Jacobian, sparse representation, and parameter mapping |
| `io.py` | `DESISpectra` class for reading DESI spectra; output FITS writing |
| `photometry.py` | `Photometry` class — filter curves (via speclite), photometric bands, dust reddening |
| `templates.py` | `Templates` class — reads SPS template FITS files, manages FFT convolution caching for velocity dispersion broadening |
| `resolution.py` | DESI resolution matrix handling; deconvolution using Koposov/rvspecfit approach; defines `SIGMA0_ANGSTROM` pre-convolution Gaussian sigma |
| `linetable.py` | Reads `data/emlines.ecsv` emission-line parameter table |
| `igm.py` | `Inoue14` IGM attenuation (Ly-α forest) |
| `cosmo.py` | Tabulated DESI fiducial cosmology interpolation |
| `mpi.py` | MPI/multiprocessing utilities for large-scale production runs |
| `qa.py` | Quality assurance figure generation |
| `linemasker.py` | Masking of spectral regions around emission lines |

### Parallelism

- **Multiprocessing**: `--mp N` flag launches N worker processes per MPI rank via `MPPool` in `util.py`
- **MPI**: `mpi-fastspecfit` script (wraps `mpi.py`) distributes work across MPI ranks; requires `mpi4py`
- Workers are initialized with `sc_data` shared state via pool initializers

### Bundled Data Files (`py/fastspecfit/data/`)

- `emlines.ecsv` — emission line table (wavelengths, line types, constraints)
- `legacysurvey-dr9.yaml` / `legacysurvey-dr10.yaml` — photometric filter/band configuration
- `stacked-phot.yaml` — photometric configuration for stacked spectra
- `desi_fiducial_cosmology.dat` — tabulated cosmology table
- `LAFcoeff.txt` / `DLAcoeff.txt` — IGM attenuation coefficients

### Performance-Critical Code

Numba `@jit` decorators are used heavily in `emline_fit/model.py`, `emline_fit/jacobian.py`, `continuum.py`, `resolution.py`, `igm.py`, and `util.py`. Expect JIT compilation overhead on first call. The `Templates` class pre-caches FFTs for velocity dispersion convolution up to `MAX_PRE_VDISP = 500 km/s`.

## Output File Structure

`fastspec` output FITS extensions:
- `METADATA` — targeting metadata, redshifts, photometry from input catalogs
- `SPECPHOT` — spectrophotometric measurements (e.g., synthesized magnitudes)
- `FASTSPEC` — fitted parameters (continuum, emission lines, stellar mass, etc.)
- `MODELS` — best-fit model spectra arrays

`fastphot` output lacks `FASTSPEC` and `MODELS` (photometry-only mode).
