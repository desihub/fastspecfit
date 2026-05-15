"""
fastspecfit.igm
===============

Tools for handling intergalactic medium (IGM) attenuation.

The :class:`Inoue14` implementation is adapted from Gabriel Brammer's
``eazy-py`` package and is used here under the MIT License.
Copyright (c) 2016-2022 Gabriel Brammer.

"""
import numpy as np
from numba import jit


class Inoue14(object):
    r"""IGM absorption model from Inoue et al. (2014).

    Parameters
    ----------
    scale_tau : float, optional
        Scaling factor applied to the IGM optical depth :math:`\tau`.
        The transmission is :math:`f_\mathrm{igm} = e^{-\mathrm{scale\_tau}\,\tau}`.
        Defaults to 1.

    Attributes
    ----------
    scale_tau : float
        Optical depth scaling factor.
    igm_params : tuple
        Pre-computed LAF and DLA absorption coefficients.
    reference : str
        Citation string for the IGM model.

    """
    igm_params = None

    def __init__(self, scale_tau=1.):

        self.scale_tau = scale_tau
        self.igm_params = self._load_data()
        self.reference = 'Inoue+2014'

    @staticmethod
    def _load_data():
        from importlib import resources
        LAF_file = resources.files('fastspecfit').joinpath('data/LAFcoeff.txt')
        DLA_file = resources.files('fastspecfit').joinpath('data/DLAcoeff.txt')

        data = np.loadtxt(LAF_file, unpack=True)
        _, lam, ALAF1, ALAF2, ALAF3 = data

        cALAF1 = ALAF1 / lam**1.2
        cALAF2 = ALAF2 / lam**3.7
        cALAF3 = ALAF3 / lam**5.5

        data = np.loadtxt(DLA_file, unpack=True)
        _, __, ADLA1, ADLA2 = data
        cADLA1 = ADLA1 / lam**2
        cADLA2 = ADLA2 / lam**3

        return (
            lam,
            cALAF1, cALAF2, cALAF3,
            cADLA1, cADLA2
        )


    def full_IGM(self, z, lobs):
        r"""Compute the full IGM transmission at observed wavelengths.

        Parameters
        ----------
        z : float
            Source redshift.
        lobs : :class:`numpy.ndarray`
            Observed-frame wavelengths in Angstroms.

        Returns
        -------
        :class:`numpy.ndarray`
            IGM transmission :math:`e^{-\tau}` at each wavelength.

        """
        return self._full_IGM(z, lobs,
                              self.scale_tau,
                              self.igm_params)


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _full_IGM(z, lobs, scale_tau, igm_params):

        lam, cALAF1, cALAF2, cALAF3, cADLA1, cADLA2 = igm_params

        tau_LS = \
            _tLSLAF(z, lobs, lam, cALAF1, cALAF2, cALAF3) + \
            _tLSDLA(z, lobs, lam, cADLA1, cADLA2)

        tau_LC = _tLCLAF(z, lobs) + _tLCDLA(z, lobs)

        ### Upturn at short wavelengths, low-z
        ### (add to tau_LC + tau_LS below)
        #k = 1./100
        #l0 = 600-6/k
        #clip = lobs/(1+z) < 600.
        #tau_clip = 100*(1-1./(1+np.exp(-k*(lobs/(1+z)-l0))))

        return np.exp(-scale_tau * (tau_LC + tau_LS))


@jit(nopython=True, fastmath=True, nogil=True, cache=True)
def _tLSLAF(zS, lobs, lam,
            cALAF1, cALAF2, cALAF3):
    """
    Lyman series, Lyman-alpha forest
    """
    z1LAF = 1. + 1.2
    z2LAF = 1. + 4.7
    z = 1. + zS

    r = np.empty_like(lobs)

    for i in range(len(lobs)):
        acc = 0.
        for j in range(len(lam)):
            if lobs[i] < lam[j] * z:
                if   lobs[i] < lam[j] * z1LAF:
                    acc += cALAF1[j] * lobs[i]**1.2
                elif lobs[i] < lam[j] * z2LAF:
                    acc += cALAF2[j] * lobs[i]**3.7
                else:
                    acc += cALAF3[j] * lobs[i]**5.5
        r[i] = acc

    return r


@jit(nopython=True, fastmath=True, nogil=True, cache=True)
def _tLSDLA(zS, lobs, lam,
            cADLA1, cADLA2):
    """
    Lyman Series, DLA
    """
    z1DLA = 1. + 2.
    z = 1. + zS

    r = np.empty_like(lobs)

    for i in range(len(lobs)):
        acc = 0.
        for j in range(len(lam)):
            if lobs[i] < lam[j] * z:
                if lobs[i] < lam[j] * z1DLA:
                    acc += cADLA1[j] * lobs[i]**2
                else:
                    acc += cADLA2[j] * lobs[i]**3
        r[i] = acc

    return r


@jit(nopython=True, fastmath=True, nogil=True, cache=True)
def _tLCDLA(zS, lobs):
    """
    Lyman continuum, DLA
    """
    z1DLA = 1. + 2.
    lamL = 911.8
    z = 1. + zS

    r = np.zeros_like(lobs)

    if z < z1DLA:
        for i in range(len(lobs)):
            if lobs[i]/lamL < z:
                r[i] = \
                    0.2113  * z**2 - \
                    0.07661 * z**2.3 * (lobs[i]/lamL) ** -0.3 - \
                    0.1347  * (lobs[i]/lamL) ** 2
    else:
        for i in range(len(lobs)):
            if lobs[i]/lamL < z:
                if lobs[i]/lamL < z1DLA:
                    r[i] =\
                        0.6340 + \
                        0.04696 * z**3 - \
                        0.01779 * z**3.3 * (lobs[i]/lamL) ** -0.3 - \
                        0.1347  * (lobs[i]/lamL) ** 2 - \
                        0.2905  * (lobs[i]/lamL) ** -0.3
                else:
                    r[i] = \
                        0.04696 * z**3 - \
                        0.01779 * z**3.3 * (lobs[i]/lamL) ** -0.3 - \
                        0.02916 * (lobs[i]/lamL) ** 3

    return r


@jit(nopython=True, fastmath=True, nogil=True, cache=True)
def _tLCLAF(zS, lobs):
    """
    Lyman continuum, LAF
    """
    z1LAF = 1. + 1.2
    z2LAF = 1. + 4.7
    lamL = 911.8
    z = 1. + zS

    r = np.zeros_like(lobs)

    if z < z1LAF:
        for i in range(len(lobs)):
            if lobs[i]/lamL < z:
                r[i] = \
                    0.3248 * ((lobs[i]/lamL) ** 1.2 - \
                              z**-0.9 * (lobs[i]/lamL) ** 2.1)
    elif z < z2LAF:
        for i in range(len(lobs)):
            if lobs[i]/lamL < z:
                if lobs[i]/lamL < z1LAF:
                    r[i] = \
                        0.02545 * z**1.6 * (lobs[i]/lamL) ** 2.1 + \
                        0.3248  * (lobs[i]/lamL) ** 1.2 - \
                        0.2496  * (lobs[i]/lamL) ** 2.1
                else:
                    r[i] = \
                        0.02545 * (z**1.6 * (lobs[i]/lamL) ** 2.1 - \
                                   (lobs[i]/lamL) ** 3.7)
    else:
        for i in range(len(lobs)):
            if lobs[i]/lamL < z:
                if lobs[i]/lamL < z1LAF:
                    r[i] = \
                        5.221e-4 * z**3.4 * (lobs[i]/lamL) ** 2.1 + \
                        0.3248   * (lobs[i]/lamL) ** 1.2 - \
                        0.03140  * (lobs[i]/lamL) ** 2.1
                elif lobs[i]/lamL < z2LAF:
                    r[i] = \
                        5.221e-4 * z**3.4 * (lobs[i]/lamL) ** 2.1 + \
                        0.2182   * (lobs[i]/lamL) ** 2.1 - \
                        0.02545  * (lobs[i]/lamL) ** 3.7
                else:
                    r[i] = \
                        5.221e-4 * (z**3.4 * (lobs[i]/lamL) ** 2.1 - \
                                    (lobs[i]/lamL) ** 5.5)

    return r
