"""
fastspecfit.test.test_resolution
=================================
Unit tests for the DESI resolution matrix handling in resolution.py.
"""
import numpy as np
import pytest


# ── shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def diag_matrix():
    """Synthetic normalized Gaussian band matrix in DESI diagonal-storage format."""
    ndiag, npix = 11, 100
    w2 = ndiag // 2
    D = np.zeros((ndiag, npix))
    for i in range(ndiag):
        k = w2 - i
        D[i, :] = np.exp(-0.5 * k**2 / 2.0**2)
    D /= D.sum(axis=0, keepdims=True)
    return D


@pytest.fixture(scope='module')
def res(diag_matrix):
    from fastspecfit.resolution import Resolution
    return Resolution(diag_matrix)


# ── _mat_torows / _mat_tocolumns ──────────────────────────────────────────────

class TestDiagConversions:

    def test_torows_shape_preserved(self, diag_matrix):
        from fastspecfit.resolution import _mat_torows
        result = _mat_torows(diag_matrix)
        assert result.shape == diag_matrix.shape

    def test_tocolumns_shape_preserved(self, diag_matrix):
        from fastspecfit.resolution import _mat_tocolumns
        result = _mat_tocolumns(diag_matrix)
        assert result.shape == diag_matrix.shape

    def test_torows_involution(self, diag_matrix):
        """_mat_torows applied twice is the identity (it is its own inverse)."""
        from fastspecfit.resolution import _mat_torows
        assert np.allclose(_mat_torows(_mat_torows(diag_matrix)), diag_matrix)

    def test_torows_preserves_center_diagonal(self, diag_matrix):
        from fastspecfit.resolution import _mat_torows
        w2 = diag_matrix.shape[0] // 2
        result = _mat_torows(diag_matrix)
        assert np.allclose(result[w2, :], diag_matrix[w2, :])

    def test_tocolumns_preserves_center_diagonal(self, diag_matrix):
        from fastspecfit.resolution import _mat_tocolumns
        w2 = diag_matrix.shape[0] // 2
        result = _mat_tocolumns(diag_matrix)
        assert np.allclose(result[w2, :], diag_matrix[w2, :])


# ── deconvolve_resolution_matrix ──────────────────────────────────────────────

class TestDeconvolve:

    def test_output_shape(self, diag_matrix):
        from fastspecfit.resolution import deconvolve_resolution_matrix
        W = deconvolve_resolution_matrix(diag_matrix)
        assert W.shape == diag_matrix.shape

    def test_nonstandard_width_uses_fallback(self):
        """Width != 11 triggers the non-cached solve path."""
        from fastspecfit.resolution import deconvolve_resolution_matrix
        ndiag, npix = 7, 80
        w2 = ndiag // 2
        D = np.zeros((ndiag, npix))
        for i in range(ndiag):
            k = w2 - i
            D[i, :] = np.exp(-0.5 * k**2 / 2.0**2)
        D /= D.sum(axis=0, keepdims=True)
        W = deconvolve_resolution_matrix(D)
        assert W.shape == (ndiag, npix)


# ── Resolution class ──────────────────────────────────────────────────────────

class TestResolution:

    def test_shape(self, diag_matrix, res):
        npix = diag_matrix.shape[1]
        assert res.shape == (npix, npix)

    def test_ndiag(self, diag_matrix, res):
        assert res.ndiag == diag_matrix.shape[0]

    def test_dot_output_shape(self, res):
        npix = res.shape[0]
        v = np.ones(npix)
        assert res.dot(v).shape == (npix,)

    def test_dot_out_parameter(self, res):
        npix = res.shape[0]
        v   = np.ones(npix)
        out = np.empty(npix)
        result = res.dot(v, out=out)
        assert result is out
        assert np.allclose(result, res.dot(v))

    def test_dot_linearity(self, res):
        """dot(a * v) == a * dot(v)."""
        rng  = np.random.default_rng(42)
        v    = rng.standard_normal(res.shape[0])
        a    = 3.7
        assert np.allclose(res.dot(a * v), a * res.dot(v))

    def test_dot_additivity(self, res):
        """dot(v1 + v2) == dot(v1) + dot(v2)."""
        rng = np.random.default_rng(0)
        v1  = rng.standard_normal(res.shape[0])
        v2  = rng.standard_normal(res.shape[0])
        assert np.allclose(res.dot(v1 + v2), res.dot(v1) + res.dot(v2))

    def test_dot_zero_vector(self, res):
        v = np.zeros(res.shape[0])
        assert np.all(res.dot(v) == 0.)

    def test_rowdata_shape(self, res):
        rows = res.rowdata()
        assert rows.shape == (res.shape[0], res.ndiag)

    def test_rowdata_cached(self, res):
        """rowdata() returns the same object on repeated calls."""
        r1 = res.rowdata()
        r2 = res.rowdata()
        assert r1 is r2
