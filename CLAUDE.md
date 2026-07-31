# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

FastSpecFit performs fast spectral energy distribution (SED) modeling and emission-line fitting of [DESI](https://desi.lbl.gov) galaxy spectra and broadband photometry. It is a production astronomy pipeline used by DESI collaborators to generate value-added catalogs.

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
- `stackfit` — fit stacked/coadded spectra
- `fastqa` — generate QA figures from output catalogs
- `mpi-fastspecfit` (bin script) — MPI parallel execution across many files
- `generate-fast-slurm` (bin script) — generate load-balanced Slurm job scripts
- `parse-fast-logs` (bin script) — parse log files and report timing statistics
- `compare-fast-vacs` (bin script) — compare two FastSpecFit VACs
- `profile-fastspec` (bin script) — profiling tool for fastspec performance
- `build-mini-specprod` (bin script) — subset Redrock/coadd/afterburner/Tractor files for a sample of TARGETIDs into a mini `$DESI_ROOT`-like tree, for testing off-NERSC
- `vac-statistics` (bin script) — print summary statistics for a FastSpecFit VAC
- `validate-fast-vac` (bin script) — run fitsverify/NaN checks across a VAC directory, or write per-directory SHA-256 checksums
- `fastbayes` (bin script) — grid-based Bayesian broadband-photometry-only fitting (see `fastspecfit.fastbayes`)
- `fastbayes-qa` (bin script) — regenerate `fastbayes` QA figures from an already-written output file
- `build-bayesian-templates` (bin script) — build the FSPS/python-fsps template grid consumed by `fastbayes`, from a YAML grid-axis config (`data/bayesian-grid.yaml`)
- `build-bayesian-photometry` (bin script) — pre-synthesize broadband photometry for a `build-bayesian-templates` grid on a fixed redshift grid

## Architecture

### Data Flow

1. **Input**: Redrock redshift catalogs + DESI coadded spectra (FITS files)
2. **`fastspec_one()`** in `fastspecfit.py` — processes a single object:
   - Calls `continuum_specfit()` to fit the stellar continuum and photometry
   - Calls `emline_specfit()` to fit emission lines on the residual spectrum
3. **Output**: Multi-extension FITS with `METADATA`, `SPECPHOT`, `FASTSPEC`, and `MODELS` HDUs

### Single-Copy Global State (`singlecopy.py`)

The `sc_data` singleton (a `Singletons` instance) is initialized once per process and shared across all worker threads/processes. It holds:
- `sc_data.templates` — stellar population synthesis templates (`Templates`)
- `sc_data.emlines` — emission line table (`LineTable`)
- `sc_data.photometry` — filter curves and photometric parameters (`Photometry`)
- `sc_data.cosmology` — tabulated DESI fiducial cosmology (`TabulatedDESI`)
- `sc_data.igm` — Inoue+14 IGM attenuation model (`Inoue14`)

This pattern avoids re-reading large files in multiprocessing workers.

`fastbayes.py` has its own parallel singleton, `bg_data` (a `BayesianGrid` instance), holding the pre-synthesized Bayesian grid (metadata, per-axis grid values/dimensions, redshift-interpolated photometry, and a lazy handle to the raw template spectra). It's populated separately from `sc_data` (via `_initialize_fastbayes_workers`) since `fastbayes` never needs the stellar template basis or emission-line tables.

### Key Modules

| Module | Role |
|---|---|
| `fastspecfit.py` | CLI parsing, top-level `fastspec`/`fastphot` drivers, per-object dispatch |
| `continuum.py` | `ContinuumTools` class — stellar SED fitting against SPS templates with dust attenuation and velocity dispersion |
| `emlines.py` | `EMFitTools` class — emission-line fitting orchestration |
| `emline_fit/` | Numba-accelerated emission-line model (Gaussian point evaluation at pixel centers), Jacobian, sparse representation, parameter mapping, and public interface (`interface.py`, `utils.py`) |
| `io.py` | `DESISpectra` class for reading DESI spectra; output FITS writing |
| `photometry.py` | `Photometry` class — filter curves (via speclite), photometric bands, dust reddening |
| `templates.py` | `Templates` class — reads SPS template FITS files, manages FFT convolution caching for velocity dispersion broadening |
| `resolution.py` | DESI resolution matrix handling; deconvolution using Koposov/rvspecfit approach; defines `SIGMA0_ANGSTROM` pre-convolution Gaussian sigma |
| `linetable.py` | Reads `data/emlines.ecsv` emission-line parameter table; `EmlineConstraints` reads `data/emline-constraints.yaml` for kinematic groups, doublet bounds, and fitting-strategy parameters |
| `igm.py` | `Inoue14` IGM attenuation (Ly-α forest) |
| `cosmo.py` | Tabulated DESI fiducial cosmology interpolation |
| `mpi.py` | MPI/multiprocessing utilities for large-scale production runs |
| `qa.py` | Quality assurance figure generation |
| `linemasker.py` | Masking of spectral regions around emission lines |
| `fastbayes.py` | `BayesianGrid` class (`bg_data` singleton) and grid-based Bayesian broadband-photometric SED fitting — closed-form per-template mass amplitude, full discrete posterior via likelihood weighting, sub-grid parabola refinement along each grid axis (axis list read per-grid-file from the templates' `AXES` FITS extension, not hard-coded) |

### Parallelism

- **Multiprocessing**: `--mp N` flag launches N worker processes per MPI rank via `MPPool` in `util.py`
- **MPI**: `mpi-fastspecfit` script (wraps `mpi.py`) distributes work across MPI ranks; requires `mpi4py`
- Workers are initialized with `sc_data` shared state via pool initializers

### Bundled Data Files (`py/fastspecfit/data/`)

- `emlines.ecsv` — emission line table (wavelengths, line types, constraints)
- `emline-constraints.yaml` — YAML-based kinematic groups, doublet bounds, amplitude constraints, and per-group fitting-strategy parameters
- `emlines-hii.ecsv` / `emline-constraints-hii.yaml` — alternate line list/constraints for HII-region fitting; select via the `--emlinesfile`/`--constraintsfile` CLI flags
- `legacysurvey-dr9.yaml` / `legacysurvey-dr10.yaml` — photometric filter/band configuration
- `stacked-phot.yaml` — photometric configuration for stacked spectra
- `desi_fiducial_cosmology.dat` — tabulated cosmology table
- `LAFcoeff.txt` / `DLAcoeff.txt` — IGM attenuation coefficients
- `bayesian-grid.yaml` — free-parameter grid-axis definitions for `bin/build-bayesian-templates` (axis kind/range/spacing, fixed FSPS params); select an alternate design via `--paramsfile`

### Performance-Critical Code

Numba `@jit` decorators are used heavily in `emline_fit/model.py`, `emline_fit/jacobian.py`, `continuum.py`, `resolution.py`, `igm.py`, and `util.py`. Expect JIT compilation overhead on first call. The `Templates` class pre-caches FFTs for velocity dispersion convolution up to `MAX_PRE_VDISP = 500 km/s`.

## Output File Structure

`fastspec` output FITS extensions:
- `METADATA` — targeting metadata, redshifts, photometry from input catalogs
- `SPECPHOT` — spectrophotometric measurements (e.g., synthesized magnitudes)
- `FASTSPEC` — fitted parameters (continuum, emission lines, stellar mass, etc.)
- `MODELS` — best-fit model spectra arrays

`fastphot` output lacks `FASTSPEC` and `MODELS` (photometry-only mode).

`fastbayes` output FITS extensions (distinct schema, written by `fastspecfit.fastbayes.write_fastbayes`):
- `METADATA` — targeting metadata, redshifts, observed/synthesized photometry
- `FASTBAYES` — per-grid-axis refined values + formal uncertainties + posterior mean/mode/percentiles, plus derived `LOGMSTAR`/`SFR`, K-corrections/absolute magnitudes, rest-frame luminosities, and `DN4000_MODEL`

No `WAVE`/`MODELS` extensions: the refined rest-frame model spectrum is cheap and deterministic to rebuild from the stored `FLUX_*`/`FLUX_IVAR_*`/`Z`/`PHOTSYS` columns (`fastspecfit.fastbayes._solve_grid` + `_build_refined_spectrum`), so `fastbayes-qa` always regenerates it on demand from the raw templates file rather than persisting it per object.
