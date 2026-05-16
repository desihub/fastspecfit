"""
fastspecfit.test.test_emline_fit
================================
Unit tests for the emline_fit subpackage: utilities, emission-line model,
parameter mapping, and Jacobian.
"""
import numpy as np
import pytest

from fastspecfit.util import C_LIGHT


# ── shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def grid():
    """Uniform log-lambda bin grid spanning the DESI wavelength range."""
    nbin = 5000
    edges = np.exp(np.linspace(np.log(3600.), np.log(9900.), nbin + 1))
    log_edges = np.log(edges)
    ibin_widths = 1. / np.diff(log_edges)
    return {'log_edges': log_edges, 'ibin_widths': ibin_widths, 'nbin': nbin}


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


# ── norm_cdf ──────────────────────────────────────────────────────────────────

class TestNormCDF:

    def test_midpoint(self):
        from fastspecfit.emline_fit.utils import norm_cdf
        assert np.isclose(norm_cdf(0.), 0.5)

    def test_symmetry(self):
        from fastspecfit.emline_fit.utils import norm_cdf
        for x in (0.5, 1.0, 2.0, 3.0):
            assert np.isclose(norm_cdf(x) + norm_cdf(-x), 1.0, atol=1e-10)

    def test_extremes(self):
        from fastspecfit.emline_fit.utils import norm_cdf
        assert norm_cdf(8.) > 1. - 1e-10
        assert norm_cdf(-8.) < 1e-10

    def test_monotone(self):
        from fastspecfit.emline_fit.utils import norm_cdf
        cdfs = np.array([norm_cdf(x) for x in np.linspace(-5, 5, 50)])
        assert np.all(np.diff(cdfs) > 0)


# ── norm_pdf ──────────────────────────────────────────────────────────────────

class TestNormPDF:

    def test_peak(self):
        from fastspecfit.emline_fit.utils import norm_pdf
        assert np.isclose(norm_pdf(0.), 1. / np.sqrt(2 * np.pi))

    def test_symmetry(self):
        from fastspecfit.emline_fit.utils import norm_pdf
        for x in (0.5, 1.0, 2.0):
            assert np.isclose(norm_pdf(x), norm_pdf(-x))

    def test_positive(self):
        from fastspecfit.emline_fit.utils import norm_pdf
        assert all(norm_pdf(x) > 0 for x in np.linspace(-5, 5, 20))


# ── max_buffer_width ──────────────────────────────────────────────────────────

class TestMaxBufferWidth:

    def test_positive(self, grid):
        from fastspecfit.emline_fit.utils import max_buffer_width
        assert max_buffer_width(grid['log_edges'], np.array([100.])) > 0

    def test_wider_sigma_gives_wider_buffer(self, grid):
        from fastspecfit.emline_fit.utils import max_buffer_width
        w_narrow = max_buffer_width(grid['log_edges'], np.array([50.]))
        w_wide   = max_buffer_width(grid['log_edges'], np.array([500.]))
        assert w_wide > w_narrow

    def test_padding(self, grid):
        from fastspecfit.emline_fit.utils import max_buffer_width
        sigmas = np.array([100.])
        w0 = max_buffer_width(grid['log_edges'], sigmas, padding=0)
        w5 = max_buffer_width(grid['log_edges'], sigmas, padding=5)
        assert w5 == w0 + 10  # 2*padding added per call


# ── emline_model ──────────────────────────────────────────────────────────────

class TestEmlineModel:

    def test_zero_amplitude_gives_zeros(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        p = p.copy(); p[0] = 0.
        model = emline_model(lw, p, grid['log_edges'], single_line['redshift'],
                             grid['ibin_widths'])
        assert np.all(model == 0.)

    def test_out_of_range_line_gives_zeros(self, grid):
        from fastspecfit.emline_fit.model import emline_model
        lw = np.array([50000.])  # far redward of the grid
        p  = np.array([1., 0., 100.])
        model = emline_model(lw, p, grid['log_edges'], 0., grid['ibin_widths'])
        assert np.all(model == 0.)

    def test_nonnegative(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        model = emline_model(lw, p, grid['log_edges'], single_line['redshift'],
                             grid['ibin_widths'])
        assert np.all(model >= 0.)

    def test_amplitude_scaling(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        m1 = emline_model(lw, p, grid['log_edges'], single_line['redshift'],
                          grid['ibin_widths'])
        p3 = p.copy(); p3[0] *= 3.
        m3 = emline_model(lw, p3, grid['log_edges'], single_line['redshift'],
                          grid['ibin_widths'])
        assert np.allclose(3. * m1, m3)

    def test_flux_conservation(self, grid, single_line):
        """Integrated model flux matches the expected Gaussian normalization."""
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        model = emline_model(lw, p, grid['log_edges'], single_line['redshift'],
                             grid['ibin_widths'])
        total    = np.sum(model / grid['ibin_widths'])
        sigma    = single_line['sigma'] / C_LIGHT
        lam_obs  = single_line['rest_wave'] * (1. + single_line['redshift'])
        expected = (np.sqrt(2 * np.pi) * sigma * np.exp(0.5 * sigma**2)
                    * single_line['amplitude'] * lam_obs)
        assert np.isclose(total, expected, rtol=1e-3)

    def test_peak_at_redshifted_wavelength(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        model = emline_model(lw, p, grid['log_edges'], single_line['redshift'],
                             grid['ibin_widths'])
        pk = np.argmax(model)
        log_edges    = grid['log_edges']
        peak_lambda  = np.exp(0.5 * (log_edges[pk] + log_edges[pk + 1]))
        expected_lam = single_line['rest_wave'] * (1. + single_line['redshift'])
        assert abs(peak_lambda - expected_lam) < 5.  # within 5 Å

    def test_vshift_moves_peak_redward(self, grid, single_line):
        from fastspecfit.emline_fit.model import emline_model
        lw, p = _line_arrays(single_line)
        m0    = emline_model(lw, p, grid['log_edges'], single_line['redshift'],
                             grid['ibin_widths'])
        p_red = p.copy(); p_red[1] = 500.  # 500 km/s positive shift
        m_red = emline_model(lw, p_red, grid['log_edges'], single_line['redshift'],
                             grid['ibin_widths'])
        assert np.argmax(m_red) > np.argmax(m0)

    def test_perline_models_match_composite(self, grid, single_line):
        """Sum of per-line profiles equals the composite model."""
        from fastspecfit.emline_fit.model import emline_model, emline_perline_models
        from fastspecfit.emline_fit.utils import max_buffer_width
        lw, p = _line_arrays(single_line)
        z     = single_line['redshift']
        le    = grid['log_edges']
        ibw   = grid['ibin_widths']

        composite = emline_model(lw, p, le, z, ibw)

        sigmas = np.array([p[2]])  # one line
        padding = 0
        endpts, profiles = emline_perline_models(lw, p, le, z, ibw, padding)

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
        model = emline_model(lw, p, grid['log_edges'], z, grid['ibin_widths'])
        endpts, dd = emline_model_jacobian(p, grid['log_edges'], grid['ibin_widths'],
                                          z, lw, padding=0)
        s, e = endpts[0]
        analytic = np.zeros_like(model)
        analytic[s:e] = dd[0, :e - s]
        assert np.allclose(analytic, model / p[0], atol=1e-12)

    def test_numerical_jacobian(self, grid, single_line):
        """Analytical Jacobian matches central finite differences for all parameters."""
        from fastspecfit.emline_fit.model import emline_model
        from fastspecfit.emline_fit.jacobian import emline_model_jacobian
        lw, p = _line_arrays(single_line)
        z   = single_line['redshift']
        le  = grid['log_edges']
        ibw = grid['ibin_widths']
        nlines = len(lw)

        endpts, dd = emline_model_jacobian(p, le, ibw, z, lw, padding=0)

        for k in range(len(p)):  # amplitude, vshift, sigma
            h   = 1e-4 * max(1., abs(p[k]))
            dp  = p.copy(); dp[k] += h
            dm  = p.copy(); dm[k] -= h
            fd  = (emline_model(lw, dp, le, z, ibw) -
                   emline_model(lw, dm, le, z, ibw)) / (2 * h)

            row = k * nlines
            s, e = endpts[row]
            analytic = np.zeros_like(fd)
            analytic[s:e] = dd[row, :e - s]

            assert np.allclose(analytic, fd, atol=1e-6, rtol=1e-3), \
                f"Jacobian mismatch for parameter index {k}"
