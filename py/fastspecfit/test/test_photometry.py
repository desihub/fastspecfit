"""
fastspecfit.test.test_photometry
================================

"""
import os
import numpy as np
import unittest

class TestUtil(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from fastspecfit.photometry import Photometry
        from fastspecfit.cosmo import TabulatedDESI
        cosmo = TabulatedDESI()
        cls.phot = Photometry()

        cls.zwave = np.logspace(np.log10(50.), np.log10(40e4), 10000)
        cls.zflux = np.ones_like(cls.zwave)
        cls.zivar = np.zeros_like(cls.zwave) + 100.

        cls.absmag_filters = np.array(['decam2014-g', 'decam2014-r', 'decam2014-z',
                                       'bessell-U', 'bessell-B', 'bessell-V', 'twomass-J',
                                       'sdss2010-u', 'sdss2010-g', 'sdss2010-r',
                                       'sdss2010-i', 'sdss2010-z', 'wise2010-W1'])

        cls.redshift = 0.2
        cls.dmod = cosmo.distance_modulus(cls.redshift)
        cls.photsys = 'S'

        nfilt = len(cls.phot.filters[cls.photsys])
        cls.nanomaggies = np.ones(nfilt)
        cls.nanomaggies_ivar = np.ones_like(cls.nanomaggies) + 100.

        cls.synth_absmag = np.array([-16.19285649, -16.81623895, -17.58668895, -17.0631772 ,
                                     -17.49934411, -17.98466151, -19.74184508, -16.85210932,
                                     -17.43688202, -18.02971177, -18.45005832, -18.82373475,
                                     -21.71275126])
        cls.rest_nanomaggies = np.array([  0.6403109 ,   1.13695657,   2.31165958,   1.42731286,
                                           2.13297587,   3.33513103,  16.82621093,   1.17514654,
                                           2.01372938,   3.47642615,   5.12001964,   7.22340889,
                                           103.35911934])
        cls.kcorr = np.array([-1.30718301e+00, -6.83789219e-01,  8.66694134e-02, -4.36832738e-01,
                              -6.69176268e-04, -1.38803532e-01,  8.47950031e-01, -6.47904968e-01,
                              -6.31308290e-02, -9.37527202e-02, -4.43827450e-01, -7.01527417e-02,
                              -9.00627082e-03])

        cls.dn4000 = 1.078368181
        cls.dn4000_ivar = 4298.47424


    @classmethod
    def tearDownClass(cls):
        pass


    def test_photometry(self):
        from fastspecfit.photometry import Photometry
        phot_stacked = Photometry(fitstack=True)
        phot_ignore = Photometry(ignore_photometry=True)

        self.assertTrue(os.path.isfile(self.phot.fphotofile))
        self.assertTrue(np.all(self.phot.absmag_filters.names == self.absmag_filters))


    def test_kcorr_and_absmag(self):
        synth_absmag, synth_maggies_rest = self.phot.synth_absmag(
            self.redshift, self.dmod, self.zwave, self.zflux)
        self.assertTrue(np.all(np.isclose(synth_absmag, self.synth_absmag, 1e-4)))
        self.assertTrue(np.all(np.isclose(1e9 * synth_maggies_rest, self.rest_nanomaggies, 1e-4)))

        kcorr, absmag, ivarabsmag, synth_mmaggies_obs = self.phot.kcorr_and_absmag(
            self.nanomaggies, self.nanomaggies_ivar, self.redshift,
            self.dmod, self.photsys, self.zwave, self.zflux, synth_absmag,
            synth_maggies_rest)
        self.assertTrue(np.all(np.isclose(kcorr, self.kcorr, 1e-4)))

        # test redshift=0.
        synth_absmag, synth_maggies_rest = self.phot.synth_absmag(
            0., self.dmod, self.zwave, self.zflux)
        self.assertTrue(np.all(synth_absmag == 0.))
        self.assertTrue(np.all(synth_maggies_rest == 0.))

        kcorr, absmag, ivarabsmag, synth_maggies_obs = self.phot.kcorr_and_absmag(
            self.nanomaggies, self.nanomaggies_ivar, 0.,
            self.dmod, self.photsys, self.zwave, self.zflux, synth_absmag,
            synth_maggies_rest)
        self.assertTrue(np.all(kcorr == 0.))
        self.assertTrue(np.all(absmag == 0.))
        self.assertTrue(np.all(ivarabsmag == 0.))
        self.assertTrue(np.all(synth_maggies_obs == 0.))


    def test_dn4000(self):
        dn4000_1, dn4000_ivar_1 = self.phot.get_dn4000(
            self.zwave, self.zflux, redshift=self.redshift, rest=False)
        self.assertTrue(np.isclose(dn4000_1, self.dn4000, 1e-4))
        self.assertTrue(np.isclose(dn4000_ivar_1, 0.))

        dn4000_2, dn4000_ivar_2 = self.phot.get_dn4000(
            self.zwave / (1. + self.redshift), self.zflux * (1. + self.redshift),
            flam_ivar=self.zivar / (1. + self.redshift)**2, redshift=None, rest=True)
        self.assertTrue(np.isclose(dn4000_2, self.dn4000, 1e-4))
        self.assertTrue(np.isclose(dn4000_ivar_2, self.dn4000_ivar, 1e-4))

        I = self.zwave / (1. + self.redshift) > 3850.
        dn4000_3, dn4000_ivar_3 = self.phot.get_dn4000(
            self.zwave[I] / (1. + self.redshift), self.zflux[I] * (1. + self.redshift),
            flam_ivar=self.zivar[I] / (1. + self.redshift)**2, redshift=None, rest=True)
        self.assertTrue(dn4000_3 == 0.)
        self.assertTrue(dn4000_ivar_3 == 0.)

        #import pdb ; pdb.set_trace()


