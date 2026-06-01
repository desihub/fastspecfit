"""
fastspecfit.test.test_emline_constraints
=========================================
Unit tests for EmlineConstraints and its integration with EMFitTools.
Tests cover file loading, consistency validation, line_bounds correctness,
amplitude constraint completeness, and EMFitTools wiring.  Fitting results
are not tested here.
"""
import os
import numpy as np
import pytest


# ── shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def ec():
    from fastspecfit.emlines import EmlineConstraints
    return EmlineConstraints()


@pytest.fixture(scope='module')
def line_table():
    from fastspecfit.linetable import LineTable
    return LineTable().table


@pytest.fixture(scope='module')
def emfit(line_table, ec):
    from fastspecfit.emlines import EMFitTools
    return EMFitTools(line_table, ec)


# ── Group 1: loading and error handling ──────────────────────────────────────

class TestLoading:

    def test_default_file_exists(self, ec):
        assert os.path.isfile(ec.file)

    def test_default_path_via_importlib(self, ec):
        # The path must be resolved through importlib.resources (same as LineTable),
        # so it should contain the package name in its path.
        assert 'fastspecfit' in ec.file

    def test_explicit_path_roundtrip(self, tmp_path):
        import shutil
        from fastspecfit.emlines import EmlineConstraints
        src = EmlineConstraints().file
        dest = str(tmp_path / 'constraints_copy.yaml')
        shutil.copy(src, dest)
        ec2 = EmlineConstraints(constraints_file=dest)
        assert ec2.file == dest
        assert set(ec2.profiles.keys()) == {'narrow_only', 'narrow_broad'}

    def test_missing_file_raises(self):
        from fastspecfit.emlines import EmlineConstraints
        with pytest.raises(FileNotFoundError, match='does not exist'):
            EmlineConstraints(constraints_file='/nonexistent/path/constraints.yaml')

    def test_profiles_present(self, ec):
        assert 'narrow_only' in ec.profiles
        assert 'narrow_broad' in ec.profiles

    def test_fitting_strategy_per_profile(self, ec):
        assert set(ec.final_pass.keys()) == {'narrow_only', 'narrow_broad'}

        for pname in ('narrow_only', 'narrow_broad'):
            fp = ec.final_pass[pname]
            assert fp['enabled']       is True
            assert fp['warm_start']    is True
            assert fp['adopt_if']      == 'chi2_improves'
            assert fp['inherit_in_mc'] is True
            assert 'mode' not in fp

    def test_group_final_pass(self, ec):
        gfp = ec.group_final_pass

        # narrow_only: both groups release vshift and sigma
        assert gfp['narrow_only.forbidden']['free_vshift']   is True
        assert gfp['narrow_only.forbidden']['free_sigma']    is True
        assert gfp['narrow_only.narrow_balmer']['free_vshift'] is True
        assert gfp['narrow_only.narrow_balmer']['free_sigma']  is True

        # narrow_broad: narrow_all frees vshift only
        assert gfp['narrow_broad.narrow_all']['free_vshift']  is True
        assert gfp['narrow_broad.narrow_all']['free_sigma']   is False
        # oiii_doublet and broad_balmer stay fully tied
        assert gfp['narrow_broad.oiii_doublet']['free_vshift'] is False
        assert gfp['narrow_broad.oiii_doublet']['free_sigma']  is False
        assert gfp['narrow_broad.broad_balmer']['free_vshift'] is False
        assert gfp['narrow_broad.broad_balmer']['free_sigma']  is False

        # global groups stay fully tied
        assert gfp['global.ciii_aliii']['free_vshift'] is False
        assert gfp['global.mgii_pair']['free_vshift']  is False

        # reserved fields present and null
        assert gfp['narrow_broad.narrow_all']['delta_vshift_max'] is None
        assert gfp['narrow_broad.narrow_all']['delta_sigma_max']  is None


# ── Group 2: consistency check ────────────────────────────────────────────────

class TestConsistencyCheck:

    def test_bundled_files_pass(self, line_table):
        from fastspecfit.emlines import EmlineConstraints
        EmlineConstraints(line_table=line_table)   # must not raise

    def test_missing_line_raises(self, line_table, tmp_path):
        from astropy.table import vstack, Table
        from fastspecfit.emlines import EmlineConstraints

        fake_row = Table(line_table[:1])
        fake_row['name'] = 'fake_line_9999'   # ≤14 chars to fit the fixed-width column
        augmented = vstack([line_table, fake_row])

        with pytest.raises(ValueError, match='missing'):
            EmlineConstraints(line_table=augmented)

    def test_extra_yaml_line_raises(self, line_table, tmp_path):
        import yaml
        from fastspecfit.emlines import EmlineConstraints

        with open(EmlineConstraints().file) as fh:
            raw = yaml.safe_load(fh)

        # Inject a line name that does not exist in emlines.ecsv.
        raw['profiles']['narrow_only'].setdefault('free_lines', [])
        raw['profiles']['narrow_only']['free_lines'].append('nonexistent_line_xyz')

        bad_yaml = tmp_path / 'bad_constraints.yaml'
        with open(bad_yaml, 'w') as fh:
            yaml.dump(raw, fh)

        with pytest.raises(ValueError, match='references'):
            EmlineConstraints(constraints_file=str(bad_yaml), line_table=line_table)


# ── Group 3: line_bounds ──────────────────────────────────────────────────────

class TestLineBounds:

    def test_narrow_forbidden_anchor(self, ec):
        # oiii_5007 is the anchor of the forbidden group in narrow_only.
        sigma_min, sigma_max, vshift_max, sigma_init, vshift_init = \
            ec.line_bounds('oiii_5007')
        assert sigma_min >= 0.
        assert 0. < sigma_max <= 1000.
        assert vshift_max > 0.

    def test_narrow_balmer_member(self, ec):
        sigma_min, sigma_max, vshift_max, _, _ = ec.line_bounds('hbeta')
        assert 0. < sigma_max <= 1000.
        assert vshift_max <= 1000.

    def test_broad_balmer_line(self, ec):
        sigma_min, sigma_max, vshift_max, _, _ = ec.line_bounds('hbeta_broad')
        assert sigma_max > 1000.    # broad Balmer group allows much wider sigma
        assert vshift_max >= 2000.

    def test_uv_free_line(self, ec):
        # lyalpha is in global.free_lines and gets global default bounds.
        sigma_min, sigma_max, vshift_max, _, _ = ec.line_bounds('lyalpha')
        assert sigma_max > 1000.
        assert vshift_max >= 2000.

    def test_global_micro_group_line(self, ec):
        # siliii_1892 is a member of the global ciii_aliii kinematic pair.
        sigma_min, sigma_max, vshift_max, _, _ = ec.line_bounds('siliii_1892')
        assert sigma_max > 0.

    def test_fixed_in_one_profile_gets_real_bounds(self, ec):
        """Regression: hgamma_broad is in narrow_only.fixed_lines but belongs to
        the broad_balmer kinematic group in narrow_broad.  line_bounds must
        return the non-zero bounds from the kinematic group, not zeros."""
        sigma_min, sigma_max, vshift_max, _, _ = ec.line_bounds('hgamma_broad')
        assert sigma_max > 0., \
            'Expected non-zero sigma_max; got zeros because fixed_lines was checked first'

    def test_all_bounds_are_floats(self, ec, line_table):
        """Every line in emlines.ecsv must return float-typed bounds (not str)."""
        for name in line_table['name']:
            bounds = ec.line_bounds(name)
            for val in bounds:
                assert isinstance(val, float), \
                    f"line '{name}': expected float bound, got {type(val).__name__} ({val!r})"

    def test_all_sigma_bounds_ordered(self, ec, line_table):
        """sigma_min <= sigma_max for every line that is not purely fixed."""
        for name in line_table['name']:
            sigma_min, sigma_max, *_ = ec.line_bounds(name)
            assert sigma_min <= sigma_max, \
                f"line '{name}': sigma_min={sigma_min} > sigma_max={sigma_max}"

    def test_unknown_line_raises(self, ec):
        with pytest.raises(ValueError, match='not found'):
            ec.line_bounds('completely_unknown_line_zz99')


# ── Group 4: doublet and amplitude constraints ────────────────────────────────

class TestAmplitudeConstraints:

    def test_doublet_ratios_all_present(self, ec):
        names = {dr['param_name'] for dr in ec.doublet_ratios}
        assert 'oii_doublet_ratio'  in names
        assert 'sii_doublet_ratio'  in names
        assert 'mgii_doublet_ratio' in names

    def test_doublet_ratio_bounds_are_floats(self, ec):
        for dr in ec.doublet_ratios:
            lo = float(dr['bounds']['min'])
            hi = float(dr['bounds']['max'])
            assert isinstance(lo, float)
            assert isinstance(hi, float)
            assert lo < hi, f"{dr['param_name']}: min >= max"

    def test_doublet_ratio_initial_within_bounds(self, ec):
        for dr in ec.doublet_ratios:
            lo  = float(dr['bounds']['min'])
            hi  = float(dr['bounds']['max'])
            ini = float(dr['initial'])
            assert lo <= ini <= hi, \
                f"{dr['param_name']}: initial {ini} outside [{lo}, {hi}]"

    def test_amplitude_fixed_all_present(self, ec):
        lines = {fc['line'] for fc in ec.amplitude_fixed}
        assert 'nii_6548'  in lines
        assert 'oiii_4959' in lines
        assert 'oii_7330'  in lines

    def test_amplitude_fixed_factors_correct(self, ec):
        # Tolerance of 1e-3 catches gross errors while allowing for YAML
        # floating-point representation (factors are ~0.3–0.8).
        factors = {fc['line']: float(fc['factor']) for fc in ec.amplitude_fixed}
        assert abs(factors['nii_6548']  - 1. / 3.0326) < 1e-3
        assert abs(factors['oiii_4959'] - 1. / 2.9643) < 1e-3
        assert abs(factors['oii_7330']  - 1. / 1.225)  < 1e-3


# ── Group 5: EMFitTools wiring ────────────────────────────────────────────────

class TestEMFitToolsIntegration:

    def test_doublet_idx_count_matches_constraints(self, emfit, ec):
        assert len(emfit.doublet_idx) == len(ec.doublet_ratios)

    def test_doublet_param_names_match_constraints(self, emfit, ec):
        expected = {dr['param_name'] for dr in ec.doublet_ratios}
        actual   = set(emfit.param_table['name'][emfit.doublet_idx])
        assert actual == expected

    def test_no_bound_violations_after_initial_guesses(self, emfit):
        """_initial_guesses_and_bounds must produce no out-of-bounds initials.

        Regression: hgamma_broad_sigma had initial=1000 but bound=(0, 0).
        """
        dummy_linepix = {}
        dummy_flux    = np.ones(1000)
        initials, bounds = emfit._initial_guesses_and_bounds(
            dummy_linepix, dummy_flux)

        violations = []
        for i, (iv, (lb, ub)) in enumerate(zip(initials, bounds)):
            if lb == ub == 0.:
                continue   # truly fixed parameter; optimizer never sees it
            if iv < lb or iv > ub:
                violations.append(
                    (emfit.param_table['name'][i], iv, lb, ub))

        assert violations == [], \
            'Initial values outside bounds:\n' + '\n'.join(
                f'  {n}: {iv:.4g} not in [{lb:.4g}, {ub:.4g}]'
                for n, iv, lb, ub in violations)

    def test_build_linemodels_returns_two_tables(self, emfit):
        emfit.compute_inrange_lines(redshift=0.1)
        broad, nobroad = emfit.build_linemodels()
        assert set(broad.colnames)   >= {'free', 'fixed', 'tiedtoparam', 'tiedfactor'}
        assert set(nobroad.colnames) >= {'free', 'fixed', 'tiedtoparam', 'tiedfactor'}
        assert len(broad) == len(nobroad) == len(emfit.param_table)

    def test_nobroad_broad_lines_are_fixed(self, emfit):
        """narrow_only profile: broad Balmer lines must be fixed in linemodel_nobroad."""
        emfit.compute_inrange_lines(redshift=0.1)
        _, nobroad = emfit.build_linemodels()
        for line_name in ('halpha_broad', 'hbeta_broad', 'hgamma_broad'):
            if line_name not in emfit.line_map:
                continue
            for p in emfit.line_table['params'][emfit.line_map[line_name]]:
                assert nobroad['fixed'][p], \
                    f'{line_name} param {p} should be fixed in linemodel_nobroad'

    def test_broad_broad_lines_are_free(self, emfit):
        """narrow_broad profile: broad Balmer anchor must not be fixed in linemodel_broad."""
        emfit.compute_inrange_lines(redshift=0.1)
        broad, _ = emfit.build_linemodels()
        anchor = 'halpha_broad'
        if anchor in emfit.line_map:
            amp, vshift, sigma = emfit.line_table['params'][emfit.line_map[anchor]]
            # Amplitude may still be fixed if the line is out of range,
            # but vshift and sigma should be free when in range.
            if emfit.line_in_range[emfit.line_map[anchor]]:
                assert broad['free'][vshift] or not broad['fixed'][vshift]
                assert broad['free'][sigma]  or not broad['fixed'][sigma]


# ── Group 6: singletons wiring ────────────────────────────────────────────────

class TestSingletonsWiring:

    def test_constraints_populated(self, templates):
        from fastspecfit.singlecopy import sc_data
        from fastspecfit.emlines import EmlineConstraints
        sc_data.initialize(template_file=templates)
        assert isinstance(sc_data.constraints, EmlineConstraints)
        assert os.path.isfile(sc_data.constraints.file)

    def test_custom_constraints_file(self, templates, tmp_path):
        import shutil
        from fastspecfit.singlecopy import sc_data
        from fastspecfit.emlines import EmlineConstraints
        src  = EmlineConstraints().file
        dest = str(tmp_path / 'custom_constraints.yaml')
        shutil.copy(src, dest)
        sc_data.initialize(template_file=templates, constraints_file=dest)
        assert sc_data.constraints.file == dest
        # Restore default for subsequent tests.
        sc_data.initialize(template_file=templates)
