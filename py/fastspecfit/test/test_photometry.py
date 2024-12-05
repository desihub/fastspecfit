"""
fastspecfit.test.test_photometry
================================

"""
import os
import pytest
import numpy as np


@pytest.fixture
def phot():
    from fastspecfit.photometry import Photometry
    phot = Photometry()
    assert(os.path.isfile(phot.fphotofile))
    yield phot


@pytest.fixture
def phot_stacked():
    from fastspecfit.photometry import Photometry
    phot_stacked = Photometry(fitstack=True)
    yield phot_stacked


@pytest.fixture
def phot_ignore():
    from fastspecfit.photometry import Photometry
    phot_ignore = Photometry(ignore_photometry=True)
    yield phot_ignore


@pytest.fixture
def data(phot):
    from fastspecfit.cosmo import TabulatedDESI
    cosmo = TabulatedDESI()

    zwave = np.logspace(np.log10(50.), np.log10(40e4), 10000)
    zflux = np.ones_like(zwave)
    zivar = np.zeros_like(zwave) + 100.

    redshift = 0.2
    dmod = cosmo.distance_modulus(redshift)
    photsys = 'S'

    absmag_filters = np.array(['decam2014-g', 'decam2014-r', 'decam2014-z',
                               'bessell-U', 'bessell-B', 'bessell-V', 'twomass-J',
                               'sdss2010-u', 'sdss2010-g', 'sdss2010-r',
                               'sdss2010-i', 'sdss2010-z', 'wise2010-W1'])

    nanomaggies = np.ones(len(phot.filters[photsys]))
    nanomaggies_ivar = np.ones_like(nanomaggies) + 100.

    data = {'zwave': zwave, 'zflux': zflux, 'zivar': zivar,
            'redshift': redshift, 'dmod': dmod, 'photsys': photsys,
            'absmag_filters': absmag_filters, 'nanomaggies': nanomaggies,
            'nanomaggies_ivar': nanomaggies_ivar}

    yield data


@pytest.fixture
def expected(phot, data):

    synth_absmag = np.array([-16.19285649, -16.81623895, -17.58668895, -17.0631772 ,
                             -17.49934411, -17.98466151, -19.74184508, -16.85210932,
                             -17.43688202, -18.02971177, -18.45005832, -18.82373475,
                             -21.71275126])
    rest_nanomaggies = np.array([0.6403109 ,  1.13695657,  2.31165958, 1.42731286,
                                 2.13297587,  3.33513103, 16.82621093, 1.17514654,
                                 2.01372938,  3.47642615,  5.12001964, 7.22340889,
                                 103.35911934])
    kcorr = np.array([-1.30718301e+00, -6.83789219e-01,  8.66694134e-02, -4.36832738e-01,
                      -6.69176268e-04, -1.38803532e-01,  8.47950031e-01, -6.47904968e-01,
                      -6.31308290e-02, -9.37527202e-02, -4.43827450e-01, -7.01527417e-02,
                      -9.00627082e-03])
    dn4000 = 1.078368181
    dn4000_ivar = 4298.47424

    expected = {'synth_absmag': synth_absmag, 'rest_nanomaggies': rest_nanomaggies,
                'kcorr': kcorr, 'dn4000': dn4000, 'dn4000_ivar': dn4000_ivar}

    yield expected


#def test_photometry(self):
#    from fastspecfit.photometry import Photometry
#    phot_stacked = Photometry(fitstack=True)
#    phot_ignore = Photometry(ignore_photometry=True)
#
#    self.assertTrue(os.path.isfile(self.phot.fphotofile))
#    self.assertTrue(np.all(self.phot.absmag_filters.names == self.absmag_filters))


def test_kcorr_and_absmag(phot, data, expected):
    synth_absmag, synth_maggies_rest = phot.synth_absmag(
        data['redshift'], data['dmod'], data['zwave'],
        data['zflux'])
    assert(np.all(np.isclose(synth_absmag, expected['synth_absmag'], 1e-4)))
    assert(np.all(np.isclose(1e9 * synth_maggies_rest, expected['rest_nanomaggies'], 1e-4)))

    kcorr, absmag, ivarabsmag, synth_mmaggies_obs = phot.kcorr_and_absmag(
        data['nanomaggies'], data['nanomaggies_ivar'], data['redshift'],
        data['dmod'], data['photsys'], data['zwave'], data['zflux'],
        synth_absmag, synth_maggies_rest)
    assert(np.all(np.isclose(kcorr, expected['kcorr'], 1e-4)))

    # test redshift=0.
    synth_absmag, synth_maggies_rest = phot.synth_absmag(
        0., data['dmod'], data['zwave'], data['zflux'])
    assert(np.all(synth_absmag == 0.))
    assert(np.all(synth_maggies_rest == 0.))

    kcorr, absmag, ivarabsmag, synth_maggies_obs = phot.kcorr_and_absmag(
        data['nanomaggies'], data['nanomaggies_ivar'], 0.,
        data['dmod'], data['photsys'], data['zwave'],
        data['zflux'], synth_absmag, synth_maggies_rest)
    assert(np.all(kcorr == 0.))
    assert(np.all(absmag == 0.))
    assert(np.all(ivarabsmag == 0.))
    assert(np.all(synth_maggies_obs == 0.))


def test_dn4000(phot, data, expected):
    dn4000_1, dn4000_ivar_1 = phot.get_dn4000(
        data['zwave'], data['zflux'], redshift=data['redshift'],
        rest=False)
    assert(np.isclose(dn4000_1, expected['dn4000'], 1e-4))
    assert(np.isclose(dn4000_ivar_1, 0.))

    dn4000_2, dn4000_ivar_2 = phot.get_dn4000(
        data['zwave'] / (1. + data['redshift']), data['zflux'] * (1. + data['redshift']),
        flam_ivar=data['zivar'] / (1. + data['redshift'])**2, redshift=None, rest=True)
    assert(np.isclose(dn4000_2, expected['dn4000'], 1e-4))
    assert(np.isclose(dn4000_ivar_2, expected['dn4000_ivar'], 1e-4))

    I = data['zwave'] / (1. + data['redshift']) > 3850.
    dn4000_3, dn4000_ivar_3 = phot.get_dn4000(
        data['zwave'][I] / (1. + data['redshift']), data['zflux'][I] * (1. + data['redshift']),
        flam_ivar=data['zivar'][I] / (1. + data['redshift'])**2, redshift=None, rest=True)
    assert(dn4000_3 == 0.)
    assert(dn4000_ivar_3 == 0.)

