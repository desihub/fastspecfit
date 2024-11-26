"""
fastspecfit.cosmo
=================

Cosmology utilities.

"""
import numpy as np

class TabulatedDESI(object):
    """
    Class to load tabulated z->E(z) and z->comoving_radial_distance(z) relations within DESI fiducial cosmology
    (in LSS/data/desi_fiducial_cosmology.dat) and perform the (linear) interpolations at any z.

    >>> cosmo = TabulatedDESI()
    >>> distance = cosmo.comoving_radial_distance([0.1, 0.2])
    >>> efunc = cosmo.efunc(0.3)

    The cosmology is defined in https://github.com/abacusorg/AbacusSummit/blob/master/Cosmologies/abacus_cosm000/CLASS.ini
    and the tabulated file was obtained using https://github.com/adematti/cosmoprimo/blob/main/cosmoprimo/fiducial.py.

    Notes
    -----
    Redshift interpolation range is [0, 100].

    """
    def __init__(self):
        from importlib import resources
        cosmofile = resources.files('fastspecfit').joinpath('data/desi_fiducial_cosmology.dat')
        self.file = cosmofile

        self._z, self._efunc, self._comoving_radial_distance = np.loadtxt(
            cosmofile, comments='#', usecols=None, unpack=True)

        self.H0 = 100.
        self.h = self.H0 / 100.
        self.hubble_time = 3.08567758e19 / 3.15576e16 / self.H0 # Hubble time [Gyr]


    def efunc(self, z):
        r"""Return :math:`E(z)`, where the Hubble parameter is defined as :math:`H(z) = H_{0} E(z)`, unitless."""
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any():
            raise ValueError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._efunc, left=None, right=None)


    def comoving_radial_distance(self, z):
        r"""Return comoving radial distance, in :math:`\mathrm{Mpc}/h`."""
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any():
            raise ValueError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._comoving_radial_distance, left=None, right=None)


    # public interface
    def luminosity_distance(self, z):
        r"""Return luminosity distance, in :math:`\mathrm{Mpc}/h`."""
        return self.comoving_radial_distance(z) * (1. + z)


    def distance_modulus(self, z):
        """Return the distance modulus at the given redshift (Hogg Eq. 24)."""
        return 5. * np.log10(self.luminosity_distance(z)) + 25.


    def universe_age(self, z):
        """Return the age of the universe at the given redshift.

        """
        from scipy.integrate import quad

        def _agefunc(z):
            return 1. / self.efunc(z) / (1. + z)


        if np.isscalar(z):
            integ, _ =  quad(_agefunc, z, self._z[-1])
            return integ * self.hubble_time
        else:
            age = []
            for _z in z:
                integ, _ =  quad(_agefunc, _z, self._z[-1])
                age.append(integ * self.hubble_time)
            return np.array(age)
