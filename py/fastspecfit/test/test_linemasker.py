"""
fastspecfit.test.test_linemasker
=================================
Unit tests for LineMasker.linepix_and_contpix in linemasker.py.

build_linemask is not tested here: it requires a fully initialized
EMFitTools instance and resolution matrices, and is covered indirectly
by the integration tests in test_fastspecfit.py.
"""
import numpy as np
import pytest
from astropy.table import Table

from fastspecfit.util import C_LIGHT


# ── shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def wave():
    """Observed-frame wavelength grid from 4000 to 8000 Å (3000 pixels)."""
    return np.linspace(4000., 8000., 3000)


@pytest.fixture(scope='module')
def ivar(wave):
    return np.ones_like(wave)


@pytest.fixture(scope='module')
def linetable():
    """Minimal emission-line table with a single OIII 5007 line."""
    return Table({'name': ['oiii_5007'], 'restwave': [5007.0]})


@pytest.fixture(scope='module')
def linesigmas():
    return np.array([150.])  # km/s


# ── helpers ───────────────────────────────────────────────────────────────────

def _call(wave, ivar, linetable, linesigmas, **kwargs):
    """Thin wrapper so each test passes a fresh copy of linesigmas."""
    from fastspecfit.linemasker import LineMasker
    return LineMasker.linepix_and_contpix(
        wave, ivar, linetable, linesigmas.copy(), **kwargs)


# ── linepix_and_contpix tests ─────────────────────────────────────────────────

class TestLinepixAndContpix:

    def test_linepix_contains_line_center(self, wave, ivar, linetable, linesigmas):
        """Pixels bracketing the line rest wavelength should be in linepix."""
        pix = _call(wave, ivar, linetable, linesigmas, redshift=0.)
        assert 'oiii_5007' in pix['linepix']
        lpix = pix['linepix']['oiii_5007']
        assert len(lpix) > 0
        # The observed wavelength of the line center should be bracketed.
        line_wave = 5007.
        assert wave[lpix[0]] <= line_wave <= wave[lpix[-1]]

    def test_contpix_does_not_overlap_linepix(self, wave, ivar, linetable, linesigmas):
        pix = _call(wave, ivar, linetable, linesigmas, redshift=0.)
        lpix = set(pix['linepix']['oiii_5007'].tolist())
        cpix = set(pix['contpix']['oiii_5007'].tolist())
        assert lpix.isdisjoint(cpix)

    def test_keys_match(self, wave, ivar, linetable, linesigmas):
        """linepix and contpix must have the same set of line names."""
        pix = _call(wave, ivar, linetable, linesigmas, redshift=0.)
        assert pix['linepix'].keys() == pix['contpix'].keys()

    def test_linepix_absent_for_out_of_range_line(self, ivar, linesigmas):
        """A line outside the wavelength range should not appear in linepix."""
        wave_short = np.linspace(4000., 6000., 1000)
        ivar_short = np.ones_like(wave_short)
        far_table  = Table({'name': ['oiii_5007'], 'restwave': [9000.]})
        pix = _call(wave_short, ivar_short, far_table, linesigmas, redshift=0.,
                    get_contpix=False)
        assert 'oiii_5007' not in pix['linepix']

    def test_zero_ivar_excluded_from_linepix(self, wave, linetable, linesigmas):
        """Pixels with ivar == 0 must not appear in linepix."""
        ivar_masked = np.ones_like(wave)
        # Mask pixels right at the line center.
        center_idx = np.argmin(np.abs(wave - 5007.))
        ivar_masked[center_idx - 2 : center_idx + 3] = 0.
        pix = _call(wave, ivar_masked, linetable, linesigmas, redshift=0.,
                    get_contpix=False)
        if 'oiii_5007' in pix['linepix']:
            masked_indices = np.where(ivar_masked == 0.)[0]
            assert len(np.intersect1d(pix['linepix']['oiii_5007'],
                                      masked_indices)) == 0

    def test_get_contpix_false_returns_empty_contpix(self, wave, ivar,
                                                      linetable, linesigmas):
        pix = _call(wave, ivar, linetable, linesigmas, redshift=0.,
                    get_contpix=False)
        assert len(pix['contpix']) == 0

    def test_positive_vshift_moves_linepix_redward(self, wave, ivar,
                                                    linetable, linesigmas):
        """A positive velocity shift should move the line mask to higher wavelength."""
        pix0 = _call(wave, ivar, linetable, linesigmas, redshift=0.,
                     linevshifts=np.array([0.]), get_contpix=False)
        pix1 = _call(wave, ivar, linetable, linesigmas, redshift=0.,
                     linevshifts=np.array([2000.]), get_contpix=False)
        center0 = np.mean(wave[pix0['linepix']['oiii_5007']])
        center1 = np.mean(wave[pix1['linepix']['oiii_5007']])
        assert center1 > center0

    def test_linevshifts_length_mismatch_raises(self, wave, ivar,
                                                 linetable, linesigmas):
        """linevshifts with a different length than linesigmas should raise ValueError."""
        from fastspecfit.linemasker import LineMasker
        with pytest.raises(ValueError, match='Mismatch'):
            LineMasker.linepix_and_contpix(
                wave, ivar, linetable, linesigmas.copy(),
                linevshifts=np.array([0., 0.]),  # wrong length
                get_contpix=False)

    def test_get_snr_adds_keys(self, wave, ivar, linetable, linesigmas):
        """get_snr=True should add snr, amp, clocal, cnoise to the result."""
        rng = np.random.default_rng(7)
        flux      = np.ones_like(wave) + 0.1 * rng.standard_normal(len(wave))
        residuals = 0.1 * rng.standard_normal(len(wave))
        pix = _call(wave, ivar, linetable, linesigmas, redshift=0.,
                    flux=flux, residuals=residuals, get_snr=True)
        assert 'oiii_5007' in pix['linepix']
        for key in ('snr', 'amp', 'clocal', 'cnoise'):
            assert 'oiii_5007' in pix[key]

    def test_redshift_moves_linepix(self, wave, ivar, linetable, linesigmas):
        """Nonzero redshift should shift linepix to longer wavelengths."""
        pix0 = _call(wave, ivar, linetable, linesigmas, redshift=0.,
                     get_contpix=False)
        pix1 = _call(wave, ivar, linetable, linesigmas, redshift=0.1,
                     get_contpix=False)
        center0 = np.mean(wave[pix0['linepix']['oiii_5007']])
        center1 = np.mean(wave[pix1['linepix']['oiii_5007']])
        assert center1 > center0
