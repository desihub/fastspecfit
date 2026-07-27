"""
fastspecfit.test.test_cosmo
===========================

"""
import os
import pytest
import numpy as np


@pytest.fixture
def cosmo():
    from fastspecfit.cosmo import TabulatedDESI
    cosmo = TabulatedDESI()
    yield cosmo


@pytest.fixture
def expected():
    zoutofrange = -1.
    redshifts = np.array([1e-4, 0.113423, 1.87988, 8.2323])
    dlums = np.array([2.99823456e-01, 3.68257334e+02, 9.95610914e+03, 5.71941290e+04])
    ages = np.array([9.28173365, 8.26036459, 2.32759459, 0.40124926])
    age_universe = 9.28271137
    yield {'redshifts': redshifts, 'dlums': dlums, 'ages': ages,
           'age_universe': age_universe, 'zoutofrange': zoutofrange}


def test_cosmo(cosmo):
    assert(os.path.isfile(cosmo.file))


def test_dlum(expected, cosmo):
    # check array
    dlums = cosmo.luminosity_distance(expected['redshifts'])
    assert(np.allclose(dlums, expected['dlums']))

    # check scalar
    dlum = cosmo.luminosity_distance(expected['redshifts'][0])
    assert(np.allclose(dlum, expected['dlums'][0]))


def test_age(expected, cosmo):
    ages = cosmo.universe_age(expected['redshifts'])
    assert(np.allclose(ages, expected['ages']))

    age_univ = cosmo.universe_age(0.)
    assert(np.isclose(age_univ, expected['age_universe']))


def test_zoutofrange(expected, cosmo):
    with pytest.raises(ValueError):
        cosmo.comoving_radial_distance(expected['zoutofrange'])

    with pytest.raises(ValueError):
        cosmo.efunc(expected['zoutofrange'])


@pytest.fixture
def flatcosmo():
    from fastspecfit.cosmo import FlatLambdaCDM
    cosmo = FlatLambdaCDM(omega_m=0.3)
    yield cosmo


@pytest.fixture
def flatexpected():
    # independently-derived expectations (analytic flat LambdaCDM, h=1),
    # not just a re-run of the implementation under test
    zoutofrange = -1.
    redshifts = np.array([1e-4, 0.113423, 1.87988, 8.2323])
    dlums = np.array([2.99815691e-01, 3.68733414e+02, 1.00826215e+04, 5.82223181e+04])
    ages = np.array([9.42591104, 8.40324420, 2.39719007, 0.42404841])
    age_universe = 9.42688876 # closed-form: (2/3) t_H / sqrt(Om_L) * asinh(sqrt(Om_L/Om_M))
    yield {'redshifts': redshifts, 'dlums': dlums, 'ages': ages,
           'age_universe': age_universe, 'zoutofrange': zoutofrange}


def test_flatcosmo_efunc(flatcosmo):
    assert(np.isclose(flatcosmo.efunc(0.), 1.))
    # efunc is analytic, so unlike TabulatedDESI it isn't restricted to
    # the interpolation range used for comoving_radial_distance
    assert(np.isclose(flatcosmo.efunc(1000.),
                       np.sqrt(flatcosmo.omega_m * 1001.**3 + flatcosmo.omega_l)))


def test_flatcosmo_dlum(flatexpected, flatcosmo):
    # check array
    dlums = flatcosmo.luminosity_distance(flatexpected['redshifts'])
    assert(np.allclose(dlums, flatexpected['dlums'], rtol=1e-5))

    # check scalar
    dlum = flatcosmo.luminosity_distance(flatexpected['redshifts'][0])
    assert(np.allclose(dlum, flatexpected['dlums'][0], rtol=1e-5))


def test_flatcosmo_age(flatexpected, flatcosmo):
    ages = flatcosmo.universe_age(flatexpected['redshifts'])
    assert(np.allclose(ages, flatexpected['ages'], rtol=1e-5))

    age_univ = flatcosmo.universe_age(0.)
    assert(np.isclose(age_univ, flatexpected['age_universe'], rtol=1e-5))


def test_flatcosmo_zoutofrange(flatexpected, flatcosmo):
    with pytest.raises(ValueError):
        flatcosmo.comoving_radial_distance(flatexpected['zoutofrange'])


class _FakeArgs(object):
    """Minimal stand-in for an argparse.Namespace, for build_cosmology()."""
    def __init__(self, omega_m=0.3):
        self.omega_m = omega_m


def test_build_cosmology_default():
    from fastspecfit.cosmo import build_cosmology
    assert(build_cosmology(None, _FakeArgs()) is None)


def test_build_cosmology_flatlcdm():
    from fastspecfit.cosmo import build_cosmology, FlatLambdaCDM
    cosmo = build_cosmology('flatLCDM', _FakeArgs(omega_m=0.25))
    assert(isinstance(cosmo, FlatLambdaCDM))
    assert(np.isclose(cosmo.omega_m, 0.25))


def test_build_cosmology_unrecognized():
    from fastspecfit.cosmo import build_cosmology
    with pytest.raises(ValueError):
        build_cosmology('not-a-real-model', _FakeArgs())
