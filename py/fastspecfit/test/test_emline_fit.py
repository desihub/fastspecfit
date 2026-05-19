"""
fastspecfit.test.test_emline_fit
================================
Unit tests for the emline_fit subpackage: utilities, emission-line model,
parameter mapping, and Jacobian.
"""
import numpy as np
import pytest

from fastspecfit.util import C_LIGHT
from fastspecfit.resolution import SIGMA0_ANGSTROM as SIGMA_G


# ── shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def grid():
    """Linear-lambda grid matching the DESI wavelength spacing (0.8 Å/pixel)."""
    dlambda = 0.8  # Angstroms, nominal DESI pixel scale
    edges = np.arange(3600., 9825., dlambda)
    centers = 0.5 * (edges[:-1] + edges[1:])
    log_centers = np.log(centers)
    bin_widths = np.full(len(centers), dlambda)
    return {'log_centers': log_centers, 'bin_widths': bin_widths, 'nbin': len(centers)}


@pytest.fixture(scope='module')
def flux_grid():
    """Linear-lambda grid around Hα at z=0.1 with DESI pixel scale (0.8 Å/pixel)."""
    dlambda = 0.8  # Angstroms
    lam_center = 6563. * 1.1  # Hα observed at z=0.1 (~7219 Å)
    nbin = 400
    edges = lam_center + (np.arange(nbin + 1) - nbin / 2) * dlambda
    centers = 0.5 * (edges[:-1] + edges[1:])
    log_centers = np.log(centers)
    bin_widths = np.full(nbin, dlambda)
    dloglam_center = dlambda / lam_center  # pixel width in log-lambda at line center
    return {'log_centers': log_centers, 'bin_widths': bin_widths,
            'nbin': nbin, 'dloglam_center': dloglam_center}


@pytest.fixture(scope='module')
def single_line():
    """One Hα emission line (6563 Å) at z = 0.1 with σ = 150 km/s."""
    return {
        'rest_wave': 6563.,
        'redshift': 0.1,
        'amplitude': 1.0,
        'vshift': 0.0,
        'sigma': 150.0,
    }


def _line_arrays(sl):
    """Return (line_wavelengths, line_parameters) arrays from a single-line dict."""
    return (np.array([sl['rest_wave']]),
            np.array([sl['amplitude'], sl['vshift'], sl['sigma']]))


def _mapping(nParms, isFree, tiedSources=None, tiedFactors=None,
             doubletTargets=None, doubletSources=None):
    """Construct a ParamsMapping, filling optional args with neutral defaults."""
    from fastspecfit.emline_fit.params_mapping import ParamsMapping
    if tiedSources is None:
        tiedSources = np.full(nParms, -1, dtype=np.int32)
    else:
        tiedSources = np.asarray(tiedSources, dtype=np.int32)
    if tiedFactors is None:
        tiedFactors = np.zeros(nParms, dtype=np.float64)
    else:
        tiedFactors = np.asarray(tiedFactors, dtype=np.float64)
    if doubletTargets is None:
        doubletTargets = np.empty(0, dtype=np.int32)
    else:
        doubletTargets = np.asarray(doubletTargets, dtype=np.int32)
    if doubletSources is None:
        doubletSources = np.empty(0, dtype=np.int32)
    else:
        doubletSources = np.asarray(doubletSources, dtype=np.int32)
    return ParamsMapping(nParms, np.asarray(isFree, dtype=bool),
                         tiedSources, tiedFactors,
                         doubletTargets, doubletSources)


# ── max_buffer_width ──────────────────────────────────────────────────────────

class TestMaxBufferWidth:

    def test_positive(self, grid):
        from fastspecfit.emline_fit.utils import max_buffer_width
        assert max_buffer_width(grid['log_centers'], np.array([100.])) > 0

    def test_wider_sigma_gives_wider_buffer(self, grid):
        from fastspecfit.emline_fit.utils import max_buffer_width
        w_narrow = max_buffer_width(grid['log_centers'], np.array([50.]))
        w_wide   = max_buffer_width(grid['log_centers'], np.array([500.]))
        assert w_wide > w_narrow

    def test_padding(self, grid):
        from fastspecfit.emline_fit.utils import max_buffer_width
        sigmas = np.array([100.])
        w0 = max_buffer_width(grid['log_centers'], sigmas, padding=0)
        w5 = max_buffer_width(grid['log_centers'], sigmas, padding=5)
        assert w5 == w0 + 10  # 2*padding added per call


# ── emline_model ──────────────────────────────────────────────────────────────

class TestEmlineModel:

    def test_zero_amplitude_gives_zeros(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        p = p.copy(); p[0] = 0.
        model = emline_model(lw, p, grid['log_centers'], single_line['redshift'], 0.)
        assert np.all(model == 0.)

    def test_out_of_range_line_gives_zeros(self, grid):
        from fastspecfit.emline_fit.model import emline_model
        lw = np.array([50000.])  # far redward of the grid
        p  = np.array([1., 0., 100.])
        model = emline_model(lw, p, grid['log_centers'], 0., 0.)
        assert np.all(model == 0.)

    def test_nonnegative(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        model = emline_model(lw, p, grid['log_centers'], single_line['redshift'], 0.)
        assert np.all(model >= 0.)

    def test_amplitude_scaling(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        m1 = emline_model(lw, p, grid['log_centers'], single_line['redshift'], 0.)
        p3 = p.copy(); p3[0] *= 3.
        m3 = emline_model(lw, p3, grid['log_centers'], single_line['redshift'], 0.)
        assert np.allclose(3. * m1, m3)

    def test_flux_conservation(self, grid, single_line):
        """Riemann sum F_σ · Δλ matches the analytic line flux sqrt(2π)·A·σ·λ*."""
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        model   = emline_model(lw, p, grid['log_centers'], single_line['redshift'], 0.)
        total   = np.sum(model * grid['bin_widths'])
        sigma   = single_line['sigma'] / C_LIGHT
        lam_obs = single_line['rest_wave'] * (1. + single_line['redshift'])
        expected = np.sqrt(2 * np.pi) * single_line['amplitude'] * sigma * lam_obs
        assert np.isclose(total, expected, rtol=2e-3)

    def test_peak_at_redshifted_wavelength(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        model       = emline_model(lw, p, grid['log_centers'], single_line['redshift'], 0.)
        pk          = np.argmax(model)
        peak_lambda = np.exp(grid['log_centers'][pk])
        expected    = single_line['rest_wave'] * (1. + single_line['redshift'])
        assert abs(peak_lambda - expected) < 5.  # within 5 Å

    def test_vshift_moves_peak_redward(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        m0    = emline_model(lw, p, grid['log_centers'], single_line['redshift'], 0.)
        p_red = p.copy(); p_red[1] = 500.  # 500 km/s positive shift
        m_red = emline_model(lw, p_red, grid['log_centers'], single_line['redshift'], 0.)
        assert np.argmax(m_red) > np.argmax(m0)

    def test_perline_models_match_composite(self, grid, single_line):
        """Sum of per-line profiles equals the composite model."""
        from fastspecfit.emline_fit.model import emline_model, emline_perline_models
        lw, p = _line_arrays(single_line)
        z     = single_line['redshift']
        lc    = grid['log_centers']

        composite          = emline_model(lw, p, lc, z, 0.)
        endpts, profiles   = emline_perline_models(lw, p, lc, z, 0., padding=0)

        reconstructed = np.zeros_like(composite)
        for j in range(len(lw)):
            s, e = endpts[j]
            reconstructed[s:e] += profiles[j, :e - s]

        assert np.allclose(reconstructed, composite)


# ── ParamsMapping ─────────────────────────────────────────────────────────────

class TestParamsMapping:

    def test_all_free_roundtrip(self):
        pm   = _mapping(3, [True, True, True])
        free = np.array([1., 2., 3.])
        assert np.allclose(pm.mapFreeToFull(free), free)

    def test_fixed_parameter_is_zero(self):
        pm   = _mapping(3, [True, False, True])
        full = pm.mapFreeToFull(np.array([5., 7.]))
        assert full[0] == 5. and full[1] == 0. and full[2] == 7.

    def test_fixed_mask(self):
        pm = _mapping(4, [True, False, True, False])
        assert np.array_equal(pm.fixedMask(), [False, True, False, True])

    def test_nFreeParms(self):
        pm = _mapping(5, [True, False, True, True, False])
        assert pm.nFreeParms == 3

    def test_tied_parameter(self):
        tiedSrc = [-1, 0, -1]
        tiedFac = [0., 2.5, 0.]
        pm   = _mapping(3, [True, False, True], tiedSrc, tiedFac)
        free = np.array([4., 3.])
        full = pm.mapFreeToFull(free)
        assert full[0] == 4.
        assert np.isclose(full[1], 4. * 2.5)
        assert full[2] == 3.

    def test_jacobian_matvec_identity(self):
        """For an all-free mapping, J_S @ v == v."""
        from fastspecfit.emline_fit.params_mapping import ParamsMapping
        pm   = _mapping(3, [True, True, True])
        free = np.array([1., 2., 3.])
        J_S  = pm.getJacobian(free)
        w    = np.empty(3)
        ParamsMapping._matvec(J_S, free, w)
        assert np.allclose(w, free)

    def test_jacobian_rmatvec_identity(self):
        """For an all-free mapping, J_S.T @ v == v."""
        from fastspecfit.emline_fit.params_mapping import ParamsMapping
        pm  = _mapping(3, [True, True, True])
        J_S = pm.getJacobian(np.array([1., 2., 3.]))
        v   = np.array([4., 5., 6.])
        w   = np.zeros(3)
        ParamsMapping._add_rmatvec(J_S, v, w)
        assert np.allclose(w, v)


# ── emline_model_jacobian ─────────────────────────────────────────────────────

class TestEmlineModelJacobian:

    def test_amplitude_derivative_exact(self, grid, single_line):
        """d(model)/d(amplitude) == model / amplitude (model is linear in A)."""
        from fastspecfit.emline_fit.model import emline_model
        from fastspecfit.emline_fit.jacobian import emline_model_jacobian
        lw, p = _line_arrays(single_line)
        z     = single_line['redshift']
        lc    = grid['log_centers']
        model      = emline_model(lw, p, lc, z, SIGMA_G)
        endpts, dd = emline_model_jacobian(p, lc, z, lw, SIGMA_G, padding=0)
        s, e = endpts[0]
        analytic = np.zeros_like(model)
        analytic[s:e] = dd[0, :e - s]
        assert np.allclose(analytic, model / p[0], atol=1e-12)

    def test_numerical_jacobian(self, grid, single_line):
        """Analytical Jacobian matches central finite differences for all parameters."""
        from fastspecfit.emline_fit.model import emline_model
        from fastspecfit.emline_fit.jacobian import emline_model_jacobian
        lw, p = _line_arrays(single_line)
        z      = single_line['redshift']
        lc     = grid['log_centers']
        nlines = len(lw)

        endpts, dd = emline_model_jacobian(p, lc, z, lw, SIGMA_G, padding=0)

        for k in range(len(p)):  # amplitude, vshift, sigma
            h   = 1e-4 * max(1., abs(p[k]))
            dp  = p.copy(); dp[k] += h
            dm  = p.copy(); dm[k] -= h
            fd  = (emline_model(lw, dp, lc, z, SIGMA_G) -
                   emline_model(lw, dm, lc, z, SIGMA_G)) / (2 * h)

            row = k * nlines
            s, e = endpts[row]
            analytic = np.zeros_like(fd)
            analytic[s:e] = dd[row, :e - s]

            assert np.allclose(analytic, fd, atol=1e-6, rtol=1e-3), \
                f"Jacobian mismatch for parameter index {k}"

    def test_zero_sigma_no_error(self, grid, single_line):
        """Jacobian does not raise when line_sigma=0 (sigma derivative has a finite limit)."""
        from fastspecfit.emline_fit.jacobian import emline_model_jacobian
        lw, p = _line_arrays(single_line)
        p = p.copy()
        p[2] = 0.  # zero line width
        endpts, dd = emline_model_jacobian(p, grid['log_centers'],
                                           single_line['redshift'], lw, SIGMA_G, padding=0)
        # sigma derivative row should be all finite (no NaN/inf)
        nlines = len(lw)
        s, e = endpts[2 * nlines]
        assert np.all(np.isfinite(dd[2 * nlines, :e - s]))


# ── flux recovery ─────────────────────────────────────────────────────────────

class TestFluxRecovery:
    """Verify that the Gaussian point-evaluation model recovers the correct
    integrated flux, including independence of sub-pixel line center position.

    The key formula: for a pre-convolved Gaussian
      F_σ(λ_j) = A · (σ/σ_eff) · exp(-x_j² / (2·σ_eff²)),
      x_j = log(λ_j) - log(λ*)
      σ_eff = sqrt(σ² + (σ_G/λ*)²)
    the Riemann sum Σ F_σ(λ_j)·Δλ_j ≈ sqrt(2π)·A·σ·λ* (the analytic line flux),
    independently of where the line center falls within the pixel grid.
    """

    def test_model_core_formula(self, flux_grid):
        """emline_model_core pixel values exactly match the analytic F_σ formula."""
        from fastspecfit.emline_fit.model import emline_model_core

        amp, vshift, sigma_kms = 2.5, 0., 80.
        rest_wave, redshift = 6563., 0.1

        vals = np.zeros(flux_grid['nbin'])
        lo, hi = emline_model_core(rest_wave, amp, vshift, sigma_kms,
                                   flux_grid['log_centers'], redshift, SIGMA_G, vals)

        sigma            = sigma_kms / C_LIGHT
        shifted_line     = rest_wave * (1. + redshift + vshift / C_LIGHT)
        log_shifted_line = np.log(shifted_line)
        sigma_G_log      = SIGMA_G / shifted_line
        sigma_eff        = np.sqrt(sigma**2 + sigma_G_log**2)
        amp_eff          = amp * sigma / sigma_eff
        inv_2_sig2       = 0.5 / sigma_eff**2

        expected = np.array([
            amp_eff * np.exp(-(flux_grid['log_centers'][i] - log_shifted_line)**2 * inv_2_sig2)
            for i in range(lo, hi)
        ])
        assert np.allclose(vals[:hi - lo], expected, rtol=1e-12)

    def test_flux_conservation(self, flux_grid):
        """Riemann sum F_σ·Δλ matches sqrt(2π)·A·σ·λ* to within 1% for a broad line."""
        from fastspecfit.emline_fit.model import emline_model

        amp, vshift, sigma_kms = 1.0, 0., 150.
        rest_wave, redshift = 6563., 0.1
        lw = np.array([rest_wave])
        p  = np.array([amp, vshift, sigma_kms])

        model   = emline_model(lw, p, flux_grid['log_centers'], redshift, SIGMA_G)
        total   = np.sum(model * flux_grid['bin_widths'])
        sigma   = sigma_kms / C_LIGHT
        lam_obs = rest_wave * (1. + redshift)
        expected = np.sqrt(2 * np.pi) * amp * sigma * lam_obs
        assert np.isclose(total, expected, rtol=1e-2)

    def test_flux_position_independence(self, flux_grid):
        """Flux sum varies by < 2% RMS as the line center moves across one pixel.

        Scans σ = [5, 15, 50] km/s.  At DESI resolution, σ_G = 0.5 Å dominates
        the effective width for narrow lines, keeping σ_eff ≥ 0.87 pixels and
        ensuring the Riemann sum is accurate regardless of sub-pixel position.
        """
        from fastspecfit.emline_fit.model import emline_model_core

        amp, rest_wave, redshift = 1.0, 6563., 0.1
        n_steps = 25
        dloglam = flux_grid['dloglam_center']

        for sigma_kms in [5., 15., 50.]:
            fluxes = []
            for step in range(n_steps):
                # shift the line center by step/n_steps of one pixel
                dvshift = (step / n_steps) * dloglam * C_LIGHT
                vals = np.zeros(flux_grid['nbin'])
                lo, hi = emline_model_core(rest_wave, amp, dvshift, sigma_kms,
                                           flux_grid['log_centers'], redshift,
                                           SIGMA_G, vals)
                if hi > lo:
                    fluxes.append(
                        np.sum(vals[:hi - lo] * flux_grid['bin_widths'][lo:hi])
                    )

            fluxes = np.array(fluxes)
            rms = np.std(fluxes) / np.mean(fluxes)
            assert rms < 0.02, (
                f"sigma={sigma_kms} km/s: sub-pixel flux RMS variation "
                f"{100*rms:.2f}% exceeds 2% threshold"
            )
