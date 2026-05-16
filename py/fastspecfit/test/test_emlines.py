"""
fastspecfit.test.test_emlines
==============================
Unit tests for emlines.py covering EMFitTools initialization,
compute_inrange_lines, build_linemodels, the chi2 static method,
and the build_coadded_models module-level function.

Methods requiring real spectra, resolution matrices, or templates
(optimize, bestfit, populate_emtable, emline_specfit, etc.) are
covered indirectly by the integration tests in test_fastspecfit.py.
"""
import numpy as np
import pytest
from astropy.table import Table


# ── shared fixtures ───────────────────────────────────────────────────────────

def _make_emfit(**kwargs):
    from fastspecfit.linetable import LineTable
    from fastspecfit.emlines import EMFitTools
    return EMFitTools(LineTable().table, **kwargs)


@pytest.fixture(scope='module')
def emfit():
    return _make_emfit()


# ── ParamType ─────────────────────────────────────────────────────────────────

class TestParamType:

    def test_values(self):
        from fastspecfit.emlines import ParamType
        assert int(ParamType.AMPLITUDE) == 0
        assert int(ParamType.VSHIFT) == 1
        assert int(ParamType.SIGMA) == 2

    def test_ordering(self):
        from fastspecfit.emlines import ParamType
        assert ParamType.AMPLITUDE < ParamType.VSHIFT < ParamType.SIGMA


# ── EMFitTools initialization ─────────────────────────────────────────────────

class TestEMFitToolsInit:

    def test_param_table_length(self, emfit):
        assert len(emfit.param_table) == 3 * len(emfit.line_table)

    def test_line_map_complete(self, emfit):
        for name in emfit.line_table['name']:
            assert name in emfit.line_map

    def test_line_map_values_are_valid_indices(self, emfit):
        nlines = len(emfit.line_table)
        for idx in emfit.line_map.values():
            assert 0 <= idx < nlines

    def test_doublet_targets_in_doublet_idx(self, emfit):
        """Known doublet targets appear in doublet_idx."""
        known_targets = {'mgii_2796', 'oii_3726', 'sii_6731'}
        target_indices = {emfit.line_map[n] for n in known_targets}
        assert target_indices.issubset(set(emfit.doublet_idx.tolist()))

    def test_isBroad_isNarrow_mutually_exclusive(self, emfit):
        """isBroad and isNarrow are mutually exclusive (broad Balmer lines belong to neither)."""
        assert not np.any(emfit.isBroad & emfit.isNarrow)

    def test_param_types_cover_all_three(self, emfit):
        from fastspecfit.emlines import ParamType
        types = set(int(t) for t in emfit.param_table['type'])
        for pt in ParamType:
            assert int(pt) in types

    def test_param_table_has_modelname(self, emfit):
        assert 'modelname' in emfit.param_table.colnames

    def test_stronglines_flag_reduces_table(self, emfit):
        """stronglines=True produces a smaller line table."""
        emfit_strong = _make_emfit(stronglines=True)
        assert len(emfit_strong.line_table) < len(emfit.line_table)


# ── compute_inrange_lines ─────────────────────────────────────────────────────

class TestComputeInrangeLines:

    def test_result_shape(self):
        emfit = _make_emfit()
        emfit.compute_inrange_lines(0.)
        assert emfit.line_in_range.shape == (len(emfit.line_table),)

    def test_z0_common_optical_lines_in_range(self):
        """At z=0 the key optical lines are within the default 3600-9900 Å window."""
        emfit = _make_emfit()
        emfit.compute_inrange_lines(0.)
        for linename in ('halpha', 'hbeta', 'oiii_5007', 'oii_3729'):
            idx = emfit.line_map[linename]
            assert emfit.line_in_range[idx], f'{linename} expected in range at z=0'

    def test_z5_far_fewer_lines_than_z0(self):
        """At z=5 far fewer lines are in range than at z=0 (most optical lines are shifted out)."""
        emfit_z0 = _make_emfit()
        emfit_z5 = _make_emfit()
        emfit_z0.compute_inrange_lines(0.)
        emfit_z5.compute_inrange_lines(5.)
        assert np.sum(emfit_z5.line_in_range) < np.sum(emfit_z0.line_in_range) // 4

    def test_narrow_wavelims_reduce_in_range_count(self):
        emfit_wide   = _make_emfit()
        emfit_narrow = _make_emfit()
        emfit_wide.compute_inrange_lines(0., wavelims=(3600., 9900.))
        emfit_narrow.compute_inrange_lines(0., wavelims=(5000., 7000.))
        assert np.sum(emfit_narrow.line_in_range) < np.sum(emfit_wide.line_in_range)

    def test_no_lines_outside_wavelims(self):
        """No in-range line falls outside the specified wavelength limits."""
        emfit = _make_emfit()
        wlo, whi = 5000., 7000.
        emfit.compute_inrange_lines(0., wavelims=(wlo, whi))
        in_range_waves = emfit.line_table['restwave'][emfit.line_in_range]
        assert np.all(in_range_waves > wlo)
        assert np.all(in_range_waves < whi)


# ── build_linemodels ──────────────────────────────────────────────────────────

class TestBuildLinemodels:

    @pytest.fixture(scope='class')
    def linemodels(self):
        emfit = _make_emfit()
        emfit.compute_inrange_lines(0.1)
        broad, nobroad = emfit.build_linemodels()
        return emfit, broad, nobroad

    def test_columns_present(self, linemodels):
        _, broad, nobroad = linemodels
        for col in ('free', 'fixed', 'tiedtoparam', 'tiedfactor'):
            assert col in broad.colnames
            assert col in nobroad.colnames

    def test_row_count_equals_param_table(self, linemodels):
        emfit, broad, nobroad = linemodels
        expected = len(emfit.param_table)
        assert len(broad) == expected
        assert len(nobroad) == expected

    def test_free_fixed_mutually_exclusive(self, linemodels):
        _, broad, nobroad = linemodels
        for lm in (broad, nobroad):
            assert not np.any(lm['free'] & lm['fixed'])

    def test_broad_has_at_least_as_many_free_params(self, linemodels):
        """Broad model must free at least as many parameters as narrow-only model."""
        _, broad, nobroad = linemodels
        assert np.sum(broad['free']) >= np.sum(nobroad['free'])

    def test_out_of_range_params_fixed_or_tied(self, linemodels):
        """Every parameter of an out-of-range line is either fixed or tied."""
        emfit, broad, _ = linemodels
        out_of_range = ~emfit.line_in_range
        for line_idx in np.where(out_of_range)[0]:
            for param_idx in emfit.line_table['params'][line_idx]:
                row = broad[param_idx]
                assert row['fixed'] or (row['tiedtoparam'] != -1), \
                    f"Out-of-range param {param_idx} is neither fixed nor tied"

    def test_tiedfactor_zero_for_untied(self, linemodels):
        """tiedfactor must be 0 for parameters that are not tied."""
        _, broad, nobroad = linemodels
        for lm in (broad, nobroad):
            untied = lm['tiedtoparam'] == -1
            assert np.all(lm['tiedfactor'][untied] == 0.)


# ── chi2 (static method) ──────────────────────────────────────────────────────

class TestChi2:

    def _lm(self, n_free):
        """Minimal linemodel Table with n_free free parameters."""
        n = n_free + 2
        free = np.zeros(n, dtype=bool)
        free[:n_free] = True
        return Table({'free': free})

    def _call(self, lm, flux, ivar, model, **kwargs):
        from fastspecfit.emlines import EMFitTools
        wave = np.arange(len(flux), dtype=float)
        return EMFitTools.chi2(lm, wave, flux, ivar, model, **kwargs)

    def test_perfect_fit_zero_chi2(self):
        flux = np.ones(20)
        chi2 = self._call(self._lm(3), flux, np.ones(20), flux)
        assert chi2 == 0.

    def test_nonzero_residuals_positive_chi2(self):
        flux  = np.ones(20)
        model = np.zeros(20)
        chi2 = self._call(self._lm(2), flux, np.ones(20), model)
        assert chi2 > 0.

    def test_dof_zero_returns_zero(self):
        """When free params exceed good pixels, dof <= 0 and chi2 is returned as 0."""
        flux = np.ones(10)
        chi2 = self._call(self._lm(15), flux, np.ones(10), np.zeros(10))
        assert chi2 == 0.

    def test_return_dof_tuple(self):
        flux = np.ones(20)
        ivar = np.ones(20)
        result = self._call(self._lm(3), flux, ivar, flux, return_dof=True)
        assert len(result) == 3
        chi2, dof, nfree = result
        assert nfree == 3
        assert dof == 17
        assert chi2 == 0.

    def test_nfree_patches_reduces_dof(self):
        """Adding nfree_patches decreases dof, increasing chi2 for the same residuals."""
        lm    = self._lm(2)
        flux  = np.ones(20)
        ivar  = np.ones(20)
        model = np.zeros(20)
        chi2_no  = self._call(lm, flux, ivar, model)
        chi2_yes = self._call(lm, flux, ivar, model, nfree_patches=3)
        assert chi2_yes > chi2_no

    def test_continuum_model_applied(self):
        """continuum_model is added to emlineflux_model before computing chi2."""
        lm   = self._lm(2)
        flux = np.ones(20)
        ivar = np.ones(20)
        chi2 = self._call(lm, flux, ivar, np.zeros(20),
                          continuum_model=np.ones(20))
        assert chi2 == 0.

    def test_zero_ivar_pixels_excluded(self):
        """Pixels with ivar=0 do not count toward dof or chi2."""
        lm   = self._lm(2)
        n    = 20
        flux = np.ones(n)
        # Half the pixels masked; same residuals → fewer good pixels
        ivar_full = np.ones(n)
        ivar_half = np.zeros(n)
        ivar_half[:n//2] = 1.
        model = np.zeros(n)
        _, dof_full, _ = self._call(lm, flux, ivar_full, model, return_dof=True)
        _, dof_half, _ = self._call(lm, flux, ivar_half, model, return_dof=True)
        assert dof_half < dof_full


# ── build_coadded_models ──────────────────────────────────────────────────────

class TestBuildCoaddedModels:

    @pytest.fixture
    def inputs(self):
        dwave = 0.8
        n = 200
        wave = 3600. + dwave * np.arange(n)
        data = {'coadd_wave': wave}
        emlinewave        = wave.copy()
        emlineflux_model  = np.ones(n)
        continuum_flux    = np.full(n, 2.)
        smooth_continuum  = np.full(n, 0.5)
        return data, emlinewave, emlineflux_model, continuum_flux, smooth_continuum

    def test_returns_four_values(self, inputs):
        from fastspecfit.emlines import build_coadded_models
        assert len(build_coadded_models(*inputs)) == 4

    def test_output_arrays_same_length(self, inputs):
        from fastspecfit.emlines import build_coadded_models
        wave_out, cont_out, emline_out, _ = build_coadded_models(*inputs)
        assert len(wave_out) == len(cont_out) == len(emline_out)

    def test_wave_spacing_preserved(self, inputs):
        from fastspecfit.emlines import build_coadded_models
        dwave_in = inputs[0]['coadd_wave'][1] - inputs[0]['coadd_wave'][0]
        wave_out, *_ = build_coadded_models(*inputs)
        assert np.allclose(np.diff(wave_out), dwave_in, rtol=1e-3)

    def test_spectra_table_columns(self, inputs):
        from fastspecfit.emlines import build_coadded_models
        _, _, _, spectra = build_coadded_models(*inputs)
        for col in ('CONTINUUM', 'SMOOTHCONTINUUM', 'EMLINEMODEL'):
            assert col in spectra.colnames

    def test_spectra_table_has_wcs_meta(self, inputs):
        from fastspecfit.emlines import build_coadded_models
        _, _, _, spectra = build_coadded_models(*inputs)
        for key in ('CRVAL1', 'CDELT1', 'NAXIS1'):
            assert key in spectra.meta

    def test_flat_continuum_interpolated_correctly(self, inputs):
        """Flat continuum_flux=2 should interpolate to ~2 on the output grid."""
        from fastspecfit.emlines import build_coadded_models
        _, cont_out, _, _ = build_coadded_models(*inputs)
        assert np.allclose(cont_out, 2., atol=1e-4)
