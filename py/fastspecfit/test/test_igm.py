"""
fastspecfit.test.test_igm
=========================

"""
import pytest
import numpy as np

@pytest.fixture
def data():
    z = 5.8
    lobs = np.arange(3600, 9800, 0.8)
    tau_indx = [4000, 5000]
    yield {'z': z, 'lobs': lobs, 'tau_indx': tau_indx}


@pytest.fixture
def expected():
    tau_vals = np.array([0.09645156, 0.08515676])
    yield {'tau_vals': tau_vals}


def test_igm(data, expected):
    from fastspecfit.igm import Inoue14
    igm = Inoue14()
    tau = igm.full_IGM(data['z'], data['lobs'])
    assert(np.all(np.isclose(tau[data['tau_indx']], expected['tau_vals'])))
