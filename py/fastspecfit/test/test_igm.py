"""
fastspecfit.test.test_igm
=========================

"""
import pytest
import numpy as np

@pytest.fixture
def setup():
    z = 5.8
    lobs = np.arange(3600, 9800, 0.8)
    tau_indx = [4000, 5000]
    tau_vals = np.array([0.09645156, 0.08515676])
    yield {'z': z, 'lobs': lobs, 'tau_indx': tau_indx, 'tau_vals': tau_vals}


def test_igm(setup):
    from fastspecfit.igm import Inoue14
    igm = Inoue14()
    tau = igm.full_IGM(setup['z'], setup['lobs'])
    assert(np.all(np.isclose(tau[setup['tau_indx']], setup['tau_vals'])))
