"""
fastspecfit.test.test_util
===========================

"""
import os
import numpy as np
import pytest

@pytest.fixture
def data():
    filtername = 'decam2014-g'
    ebv = np.array([0.07, 0.458])
    yield {'filtername': filtername, 'ebv': ebv}


@pytest.fixture
def expected():
    mwdust = np.array([0.81295031, 0.25796548])
    yield {'mwdust': mwdust}


def test_ZWarningMask():
    from fastspecfit.util import ZWarningMask
    assert(ZWarningMask.NODATA == 2**9)
    assert(len(ZWarningMask.flags()) == 12)


def test_mwdust(data, expected):
    from fastspecfit.util import mwdust_transmission
    mwdust = mwdust_transmission(data['ebv'], data['filtername'])
    assert(np.allclose(mwdust, expected['mwdust']))


def test_var_tofrom_ivar():
    from fastspecfit.util import var2ivar, ivar2var

    assert(var2ivar(10.) == 0.1)
    assert(var2ivar(10., sigma=True) == 0.01)
    assert(var2ivar(0.) == 0.)
    assert(var2ivar(0., sigma=True) == 0.)

    vals = np.array([0., 10., 0., 3., 100.])
    vals2 = np.array([0., 1e-4, 0.1])

    var, mask = ivar2var(vals)
    assert(mask[0] == False & mask[2] == False)

    var, mask = ivar2var(vals, sigma=True)
    assert(var[4] == 0.1)

    with pytest.raises(ValueError):
        ivar2var(np.array([0.]))

    var, mask = ivar2var(np.array([0.]), allmasked_ok=True)
    assert(mask[0] == False)

    var, mask = ivar2var(vals2)
    assert(np.all(var[[0, 1]] == 0.))


def test_minima():
    from fastspecfit.util import find_minima, minfit
    xx = np.arange(-7, 8, 0.5)
    x0 = 1.6
    yy = 2 * (xx-x0)**2

    imin = find_minima(yy)
    a, b, c, warn = minfit(xx, yy)
    assert(xx[imin] == 1.5)
    assert(np.isclose(a, x0))


def test_stats():
    from fastspecfit.util import sigmaclip, quantile, median
    xx = np.arange(35)
    xgood, mask = sigmaclip(xx)
    assert(np.all(mask))

    q1, q2 = quantile(xx, np.asarray([0.1, 0.9]))
    assert(np.isclose(q1, 3.4) & np.isclose(q2, 30.6))

    with pytest.raises(ValueError):
        quantile(xx, np.asarray([10., 90.]))

    assert(np.isclose(median(xx), 17.))
