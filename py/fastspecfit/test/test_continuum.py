"""
fastspecfit.test.test_continuum
================================
Unit tests for continuum.py covering standalone functions and
ContinuumTools methods that do not require the template download or
a fully initialized ContinuumTools instance.

Methods requiring templates (build_stellar_continuum,
continuum_to_photometry, fit_stellar_continuum, etc.) are exercised
indirectly by the integration tests in test_fastspecfit.py.
"""
import numpy as np
import pytest


# ── can_compute_vdisp ─────────────────────────────────────────────────────────

class TestCanComputeVdisp:

    def _call(self, z, wave, **kwargs):
        from fastspecfit.continuum import can_compute_vdisp
        return can_compute_vdisp(z, wave, **kwargs)

    def test_wide_coverage_at_z0(self):
        """Spectrum covering 3500-7000 Å at z=0 reaches the Ca H&K region."""
        wave = np.linspace(3500., 7000., 1000)
        ok, pix = self._call(0., wave)
        assert ok
        assert pix != (0, 0)

    def test_too_red_at_z0(self):
        """Spectrum starting at 5000 Å at z=0 does not reach below 3800 Å."""
        wave = np.linspace(5000., 9000., 1000)
        ok, pix = self._call(0., wave)
        assert not ok
        assert pix == (0, 0)

    def test_redshifted_coverage(self):
        """At z=0.5 a 5400-10500 Å spectrum covers 3600-7000 Å rest-frame."""
        z = 0.5
        wave = np.linspace(5400., 10500., 2000)
        ok, pix = self._call(z, wave)
        assert ok

    def test_pixel_range_within_bounds(self):
        """Returned pixel range must be valid indices into the wave array."""
        wave = np.linspace(3500., 7000., 1000)
        ok, (s, e) = self._call(0., wave)
        assert ok
        assert 0 <= s < e <= len(wave)

    def test_just_missing_blue_edge(self):
        """A spectrum starting just above the minimum rest-frame limit returns False."""
        min_lo = 3800.
        wave = np.linspace(min_lo + 1., 9000., 2000)  # starts above threshold
        ok, pix = self._call(0., wave)
        assert not ok


# ── _younger_than_universe ────────────────────────────────────────────────────

class TestYoungerThanUniverse:

    def _call(self, age_yr, tuniv_gyr, agepad=0.5):
        from fastspecfit.continuum import ContinuumTools
        return ContinuumTools._younger_than_universe(age_yr, tuniv_gyr, agepad)

    def test_all_younger(self):
        """All templates younger than the universe are kept."""
        age = np.array([1e8, 5e8, 1e9])   # yr — all < 5.5 Gyr
        tuniv = 5.0                         # Gyr
        idx = self._call(age, tuniv)
        assert len(idx) == len(age)

    def test_all_older(self):
        """No templates older than the universe are kept."""
        age = np.array([7e9, 10e9, 13e9])  # yr — all > 5.5 Gyr
        tuniv = 5.0
        idx = self._call(age, tuniv)
        assert len(idx) == 0

    def test_mixed_ages(self):
        """Only templates younger than the threshold are returned."""
        tuniv = 5.0   # Gyr
        agepad = 0.5
        threshold = 1e9 * (tuniv + agepad)  # = 5.5e9 yr
        age = np.array([1e9, 3e9, 5e9, 6e9, 13e9])
        idx = self._call(age, tuniv, agepad)
        expected = np.where(age <= threshold)[0]
        assert np.array_equal(idx, expected)


# ── lums_keys / cfluxes_keys ─────────────────────────────────────────────────

class TestKeyLists:

    def test_lums_keys_length_matches(self):
        from fastspecfit.continuum import ContinuumTools
        keys, waves = ContinuumTools.lums_keys()
        assert len(keys) == len(waves)
        assert all(isinstance(k, str) for k in keys)
        assert all(w > 0 for w in waves)

    def test_cfluxes_keys_length_matches(self):
        from fastspecfit.continuum import ContinuumTools
        keys, waves = ContinuumTools.cfluxes_keys()
        assert len(keys) == len(waves)
        assert all(isinstance(k, str) for k in keys)
        assert all(w > 0 for w in waves)


# ── attenuate_nodust ──────────────────────────────────────────────────────────

class TestAttenuateNoDust:

    def _call(self, M, A, zfactors):
        from fastspecfit.continuum import ContinuumTools
        Mc = M.copy().astype(np.float64)
        ContinuumTools.attenuate_nodust(Mc, A, zfactors)
        return Mc

    def test_no_op_when_A_and_zfactors_are_one(self):
        """A=1, zfactors=1: spectrum is unchanged."""
        M = np.array([1., 2., 3., 4., 5.])
        result = self._call(M, np.ones_like(M), np.ones_like(M))
        assert np.allclose(result, M)

    def test_complete_attenuation(self):
        """A=0 zeroes the spectrum regardless of zfactors."""
        M = np.array([1., 2., 3.])
        result = self._call(M, np.zeros_like(M), np.ones_like(M))
        assert np.all(result == 0.)

    def test_scaling_by_zfactors(self):
        """zfactors uniformly scales the output."""
        M = np.array([1., 2., 3.])
        scale = 2.5
        r1 = self._call(M, np.ones_like(M), np.ones_like(M))
        r2 = self._call(M, np.ones_like(M), np.full_like(M, scale))
        assert np.allclose(r2, scale * r1)


# ── attenuate ─────────────────────────────────────────────────────────────────

class TestAttenuate:

    def _call(self, M, A, zfactors, wave, dustflux):
        from fastspecfit.continuum import ContinuumTools
        Mc = M.copy().astype(np.float64)
        ContinuumTools.attenuate(Mc, A, zfactors, wave, dustflux)
        return Mc

    def test_no_op_when_A_and_zfactors_one_no_dust(self):
        """A=1, zfactors=1, dustflux=0: lbol_diff=0 so spectrum is unchanged."""
        M = np.array([1., 2., 3., 4.])
        wave = np.array([3000., 4000., 5000., 6000.], dtype=np.float64)
        result = self._call(M, np.ones_like(M), np.ones_like(M), wave,
                            np.zeros_like(M))
        assert np.allclose(result, M)

    def test_complete_attenuation_no_dust(self):
        """A=0, dustflux=0: all flux is absorbed and no dust emission added."""
        M = np.array([1., 2., 3., 4.])
        wave = np.array([3000., 4000., 5000., 6000.], dtype=np.float64)
        result = self._call(M, np.zeros_like(M), np.ones_like(M), wave,
                            np.zeros_like(M))
        assert np.all(result == 0.)

    def test_dustflux_zero_matches_nodust(self):
        """With dustflux=0, attenuate and attenuate_nodust give identical results."""
        from fastspecfit.continuum import ContinuumTools
        rng = np.random.default_rng(0)
        n = 50
        M     = rng.uniform(0.5, 2., n)
        A     = rng.uniform(0.3, 1., n)
        zf    = rng.uniform(0.5, 1.5, n)
        wave  = np.linspace(3000., 9000., n)

        r_att     = self._call(M, A, zf, wave, np.zeros(n))
        Mc_nodust = M.copy().astype(np.float64)
        ContinuumTools.attenuate_nodust(Mc_nodust, A, zf)
        assert np.allclose(r_att, Mc_nodust)


# ── stellar_continuum_chi2 ────────────────────────────────────────────────────

class TestStellarContinuumChi2:

    def _call(self, resid, ncoeff, vdisp_fitted, **kwargs):
        from fastspecfit.continuum import ContinuumTools
        # self is never accessed inside this method.
        return ContinuumTools.stellar_continuum_chi2(
            None, resid, ncoeff, vdisp_fitted, **kwargs)

    def test_zero_residuals_give_zero_chi2(self):
        resid = np.zeros(20)
        rchi2_spec, rchi2_phot, rchi2_tot = self._call(
            resid, ncoeff=3, vdisp_fitted=False,
            split=10, ndof_spec=10, ndof_phot=10)
        assert rchi2_spec == 0.
        assert rchi2_phot == 0.
        assert rchi2_tot == 0.

    def test_spec_only(self):
        """ndof_phot=0: rchi2_phot=0 and rchi2_tot equals rchi2_spec."""
        resid = np.ones(10)
        rchi2_spec, rchi2_phot, rchi2_tot = self._call(
            resid, ncoeff=2, vdisp_fitted=False,
            split=10, ndof_spec=10, ndof_phot=0)
        assert rchi2_phot == 0.
        assert rchi2_tot == rchi2_spec

    def test_phot_only(self):
        """ndof_spec=0: rchi2_spec=0 and rchi2_tot equals rchi2_phot."""
        resid = np.ones(5)
        rchi2_spec, rchi2_phot, rchi2_tot = self._call(
            resid, ncoeff=2, vdisp_fitted=False,
            split=0, ndof_spec=0, ndof_phot=5)
        assert rchi2_spec == 0.
        assert rchi2_tot == rchi2_phot

    def test_vdisp_fitted_increases_nfree(self):
        """Fitting vdisp adds one degree of freedom, raising rchi2."""
        resid = np.ones(20)
        _, _, rchi2_no  = self._call(resid, ncoeff=3, vdisp_fitted=False,
                                     split=10, ndof_spec=10, ndof_phot=10)
        _, _, rchi2_yes = self._call(resid, ncoeff=3, vdisp_fitted=True,
                                     split=10, ndof_spec=10, ndof_phot=10)
        assert rchi2_yes > rchi2_no

    def test_split_correctly_separates_spec_and_phot(self):
        """Residuals before split go to spec chi2; residuals after go to phot."""
        resid = np.concatenate([np.zeros(10), np.ones(5)])
        rchi2_spec, rchi2_phot, _ = self._call(
            resid, ncoeff=2, vdisp_fitted=False,
            split=10, ndof_spec=10, ndof_phot=5)
        assert rchi2_spec == 0.
        assert rchi2_phot > 0.


# ── smooth_continuum ──────────────────────────────────────────────────────────

class TestSmoothContinuum:

    @pytest.fixture
    def flat_spectrum(self):
        rng = np.random.default_rng(7)
        n = 1000
        wave = np.linspace(4000., 8000., n)
        flux = 5.0 + 0.3 * rng.standard_normal(n)
        ivar = np.ones(n)
        linemask = np.zeros(n, bool)
        camerapix = np.array([[0, n]])
        return wave, flux, ivar, linemask, camerapix

    def test_output_shape(self, flat_spectrum):
        from fastspecfit.continuum import ContinuumTools
        wave, flux, ivar, linemask, camerapix = flat_spectrum
        result = ContinuumTools.smooth_continuum(wave, flux, ivar, linemask,
                                                 camerapix)
        assert result.shape == wave.shape

    def test_linemask_length_mismatch_raises(self, flat_spectrum):
        from fastspecfit.continuum import ContinuumTools
        wave, flux, ivar, linemask, camerapix = flat_spectrum
        bad_mask = np.zeros(len(wave) + 1, bool)
        with pytest.raises(ValueError, match='Linemask'):
            ContinuumTools.smooth_continuum(wave, flux, ivar, bad_mask,
                                            camerapix)

    def test_smooth_of_flat_spectrum_near_mean(self, flat_spectrum):
        """Smooth continuum of a noisy flat spectrum should be near the mean."""
        from fastspecfit.continuum import ContinuumTools
        wave, flux, ivar, linemask, camerapix = flat_spectrum
        result = ContinuumTools.smooth_continuum(wave, flux, ivar, linemask,
                                                 camerapix)
        unmasked = result[result != 0.]
        assert len(unmasked) > 0
        assert abs(np.median(unmasked) - 5.0) < 1.0

    def test_all_masked_ivar_returns_zeros(self, flat_spectrum):
        """With ivar=0 everywhere, the smooth continuum should be all zeros."""
        from fastspecfit.continuum import ContinuumTools
        wave, flux, _, linemask, camerapix = flat_spectrum
        ivar_zero = np.zeros_like(wave)
        flux_zero = np.zeros_like(wave)
        result = ContinuumTools.smooth_continuum(wave, flux_zero, ivar_zero,
                                                 linemask, camerapix)
        assert np.all(result == 0.)
