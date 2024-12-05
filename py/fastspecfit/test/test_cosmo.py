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
