"""
fastspecfit.test.test_cosmo
===========================

"""
import os
import numpy as np
import unittest

class TestCosmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from fastspecfit.cosmo import TabulatedDESI
        cls.cosmo = TabulatedDESI()
        cls.redshifts = np.array([1e-4, 0.113423, 1.87988, 8.2323])
        cls.dlums = np.array([2.99823456e-01, 3.68257334e+02, 9.95610914e+03, 5.71941290e+04])
        cls.ages = np.array([9.28173365, 8.26036459, 2.32759459, 0.40124926])


    @classmethod
    def tearDownClass(cls):
        pass


    def test_cosmo(self):
        self.assertTrue(os.path.isfile(self.cosmo.file))


    def test_dlum(self):
        dlums = self.cosmo.luminosity_distance(self.redshifts)
        self.assertTrue(np.allclose(dlums, self.dlums))


    def test_age(self):
        ages = self.cosmo.universe_age(self.redshifts)
        self.assertTrue(np.allclose(ages, self.ages))
