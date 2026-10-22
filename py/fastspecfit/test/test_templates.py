"""
fastspecfit.test.test_templates
===============================

"""
import os
import numpy as np
import pytest


# See conftest.py for fixtures used by multiple unit tests.

def test_templates_nofilename(templatedir, template_version, templates):
    from fastspecfit.templates import Templates
    os.environ['FTEMPLATES_DIR'] = str(templatedir.parent)
    temp = Templates()
    assert(os.path.isfile(temp.file))


def test_templates(templates):
    from fastspecfit.templates import Templates
    temp = Templates(template_file=templates)
    assert(os.path.isfile(temp.file))
    assert(hasattr(temp, 'pixkms_bounds'))


# ── vdisp_nominal_kernel ──────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def loaded_T(templates):
    from fastspecfit.templates import Templates
    return Templates(template_file=templates)


class TestVdispNominal:
    """Tests for the vdisp_nominal / vdisp_nominal_kernel distinction."""

    def test_kernel_attribute_value(self, loaded_T):
        """vdisp_nominal_kernel == sqrt(vdisp_nominal² − SIGMA_C3K²)."""
        from fastspecfit.templates import Templates, VDISP_NOMINAL
        expected = float(np.sqrt(max(0., VDISP_NOMINAL**2 - Templates.SIGMA_C3K**2)))
        assert np.isclose(loaded_T.vdisp_nominal_kernel, expected)

    def test_kernel_less_than_nominal(self, loaded_T):
        """Pre-convolution kernel is always less than the reported σ_stars."""
        assert loaded_T.vdisp_nominal_kernel < loaded_T.vdisp_nominal

    def test_flux_nomvdisp_uses_kernel(self, loaded_T):
        """flux_nomvdisp is broadened by vdisp_nominal_kernel."""
        T = loaded_T
        broadened_kernel = T.convolve_vdisp(T.flux, T.vdisp_nominal_kernel)
        assert np.allclose(T.flux_nomvdisp, broadened_kernel)

    def test_reversed_bounds_raises(self, templates):
        """Reversed vdisp_bounds raises ValueError with a clear message."""
        from fastspecfit.templates import Templates
        with pytest.raises(ValueError, match='lo <= hi'):
            Templates(template_file=templates, vdisp_bounds=(500., 50.))


# ── σ–M* relation constants and kernel formula ────────────────────────────────

class TestVdispSigmaRelation:

    def test_default_coefficients(self):
        from fastspecfit.templates import VDISP_SIGMA_RELATION
        assert VDISP_SIGMA_RELATION == (2.30, 0.25)

    def test_lower_bound_is_50(self):
        from fastspecfit.templates import VDISP_BOUNDS
        assert VDISP_BOUNDS[0] == 50.

    def test_formula_at_reference_mass(self):
        """At log M* = 11, σ_stars = 10^a exactly (the b term vanishes)."""
        from fastspecfit.templates import VDISP_SIGMA_RELATION
        a, b = VDISP_SIGMA_RELATION
        vdisp_stars = 10.**(a + b * (11.0 - 11.))
        assert np.isclose(vdisp_stars, 10.**a)

    def test_kernel_roundtrip(self):
        """sqrt(kernel² + SIGMA_C3K²) recovers σ_stars to floating-point precision."""
        from fastspecfit.templates import Templates, VDISP_SIGMA_RELATION, VDISP_BOUNDS
        a, b = VDISP_SIGMA_RELATION
        vdisp_stars = 10.**(a + b * (11.0 - 11.))
        vdisp_kernel = float(np.sqrt(max(0., vdisp_stars**2 - Templates.SIGMA_C3K**2)))
        vdisp_kernel = float(np.clip(vdisp_kernel, *VDISP_BOUNDS))
        assert np.isclose(np.sqrt(vdisp_kernel**2 + Templates.SIGMA_C3K**2), vdisp_stars)

    def test_low_mstar_clamps_to_lower_bound(self):
        """Very low M* (σ_stars < SIGMA_C3K) clamps kernel to vdisp_bounds[0]."""
        from fastspecfit.templates import Templates, VDISP_SIGMA_RELATION, VDISP_BOUNDS
        a, b = VDISP_SIGMA_RELATION
        vdisp_stars = 10.**(a + b * (6.0 - 11.))  # ~11 km/s << SIGMA_C3K ~42 km/s
        vdisp_kernel = float(np.sqrt(max(0., vdisp_stars**2 - Templates.SIGMA_C3K**2)))
        vdisp_kernel = float(np.clip(vdisp_kernel, *VDISP_BOUNDS))
        assert vdisp_kernel == VDISP_BOUNDS[0]


# ── velocity-dispersion convolution at PIXKMS_BOUNDS ──────────────────────────

class TestConvolveVdispEdges:
    """The convolved range is padded, so its edge pixels are not darkened."""

    VDISPS = (75., 250., 500.)

    def _edges(self, T, ncheck=20):
        lo, hi = T.pixkms_bounds
        return np.r_[lo:lo+ncheck, hi-ncheck:hi]

    def test_constant_preserved(self, loaded_T):
        T = loaded_T
        flux = np.ones(T.npix)
        for vdisp in self.VDISPS + (600.,): # >MAX_PRE_VDISP exercises _pad_interp
            conv = T.convolve_vdisp(flux, vdisp)
            assert np.allclose(conv, 1., rtol=0., atol=1e-10)

    def test_powerlaw_preserved(self, loaded_T):
        """A smooth power law is unchanged to well below the old edge dip."""
        T = loaded_T
        flux = (T.wave / 5000.)**(-2.)
        edges = self._edges(T)
        for vdisp in self.VDISPS:
            conv = T.convolve_vdisp(flux, vdisp)
            assert np.allclose(conv[edges] / flux[edges], 1., rtol=0., atol=1e-3)

    def test_pre_matches_direct(self, loaded_T):
        """The precomputed-FFT path agrees with convolve_vdisp everywhere."""
        T = loaded_T
        rows = np.arange(min(3, T.ntemplates))
        flux = T.flux_nolines[rows, :]
        flux_lohi, ft_flux_mid, fft_len = T.conv_pre_select(T.conv_pre_nolines, rows)
        for vdisp in self.VDISPS:
            direct = T.convolve_vdisp(flux, vdisp)
            for ii in range(len(rows)):
                pre = T.convolve_vdisp_from_pre(flux_lohi[ii], ft_flux_mid[ii],
                                                T.npix, fft_len, vdisp)
                assert np.allclose(pre, direct[ii], rtol=1e-8,
                                   atol=1e-10 * np.max(np.abs(direct[ii])))

