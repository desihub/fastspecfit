"""
fastspecfit.test.test_util
===========================

"""
import os
import numpy as np
import unittest

class TestUtil(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filtername = 'decam2014-g'
        cls.ebv = np.array([0.07, 0.458])
        cls.mwdust = np.array([0.81295031, 0.25796548])


    @classmethod
    def tearDownClass(cls):
        pass


    def test_ZWarningMask(self):
        from fastspecfit.util import ZWarningMask
        self.assertTrue(ZWarningMask.NODATA == 2**9)
        self.assertTrue(len(ZWarningMask.flags()) == 12)


    def test_mwdust(self):
        from fastspecfit.util import mwdust_transmission
        mwdust = mwdust_transmission(self.ebv, self.filtername)
        self.assertTrue(np.allclose(mwdust, self.mwdust))


    def test_var_tofrom_ivar(self):
        from fastspecfit.util import var2ivar, ivar2var

        self.assertTrue(var2ivar(10.) == 0.1)
        self.assertTrue(var2ivar(10., sigma=True) == 0.01)
        self.assertTrue(var2ivar(0.) == 0.)
        self.assertTrue(var2ivar(0., sigma=True) == 0.)

        vals = np.array([0., 10., 0., 3., 100.])
        vals2 = np.array([0., 1e-4, 0.1])

        var, mask = ivar2var(vals)
        self.assertTrue(mask[0] == False & mask[2] == False)

        var, mask = ivar2var(vals, sigma=True)
        self.assertTrue(var[4] == 0.1)

        with self.assertRaises(ValueError):
            ivar2var(np.array([0.]))

        var, mask = ivar2var(np.array([0.]), allmasked_ok=True)
        self.assertTrue(mask[0] == False)

        var, mask = ivar2var(vals2)
        self.assertTrue(np.all(var[[0, 1]] == 0.))


    def test_minima(self):
        from fastspecfit.util import find_minima, minfit
        xx = np.arange(-7, 8, 0.5)
        x0 = 1.6
        yy = 2 * (xx-x0)**2

        imin = find_minima(yy)
        a, b, c, warn = minfit(xx, yy)
        self.assertTrue(xx[imin] == 1.5)
        self.assertTrue(np.isclose(a, x0))


    def test_stats(self):
        from fastspecfit.util import sigmaclip, quantile, median
        xx = np.arange(35)
        xgood, mask = sigmaclip(xx)
        self.assertTrue(np.all(mask))

        q1, q2 = quantile(xx, [0.1, 0.9])
        self.assertTrue(np.isclose(q1, 3.4) & np.isclose(q2, 30.6))

        with self.assertRaises(ValueError):
            quantile(xx, [10., 90.])

        self.assertTrue(np.isclose(median(xx), 17.))
