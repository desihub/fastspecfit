"""
fastspecfit.cosmo
=================

Cosmology utilities.

"""
import numpy as np

class TabulatedDESI(object):
    """Tabulated DESI fiducial cosmology for fast redshift interpolation.

    Loads tabulated :math:`E(z)` and comoving radial distance as a function
    of redshift and performs linear interpolation. The cosmology parameters
    are defined by the AbacusSummit baseline (Planck 2018 ΛCDM).

    Attributes
    ----------
    file : str
        Path to the tabulated cosmology file.
    H0 : float
        Hubble constant in km/s/Mpc.
    h : float
        Dimensionless Hubble parameter.
    hubble_time : float
        Hubble time in Gyr.

    Notes
    -----
    Redshift interpolation range is [0, 100]. Cosmology defined at
    https://github.com/abacusorg/AbacusSummit; tabulated file generated
    with https://github.com/adematti/cosmoprimo.

    Examples
    --------
    >>> cosmo = TabulatedDESI()
    >>> distance = cosmo.comoving_radial_distance([0.1, 0.2])
    >>> efunc = cosmo.efunc(0.3)

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
        r"""Return :math:`E(z) = H(z) / H_0`.

        Parameters
        ----------
        z : float or array-like
            Redshift.

        Returns
        -------
        float or :class:`numpy.ndarray`
            Dimensionless Hubble parameter at the given redshift(s).

        """
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any():
            raise ValueError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._efunc, left=None, right=None)


    def comoving_radial_distance(self, z):
        r"""Return comoving radial distance.

        Parameters
        ----------
        z : float or array-like
            Redshift.

        Returns
        -------
        float or :class:`numpy.ndarray`
            Comoving radial distance in :math:`\mathrm{Mpc}/h`.

        """
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any():
            raise ValueError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._comoving_radial_distance, left=None, right=None)


    # public interface
    def luminosity_distance(self, z):
        r"""Return luminosity distance.

        Parameters
        ----------
        z : float or array-like
            Redshift.

        Returns
        -------
        float or :class:`numpy.ndarray`
            Luminosity distance in :math:`\mathrm{Mpc}/h`.

        """
        return self.comoving_radial_distance(z) * (1. + z)


    def distance_modulus(self, z):
        """Return the distance modulus.

        Parameters
        ----------
        z : float or array-like
            Redshift.

        Returns
        -------
        float or :class:`numpy.ndarray`
            Distance modulus in magnitudes (Hogg 1999, Eq. 24).

        """
        return 5. * np.log10(self.luminosity_distance(z)) + 25.


    def universe_age(self, z):
        """Return the age of the universe at the given redshift.

        Parameters
        ----------
        z : float or array-like
            Redshift.

        Returns
        -------
        float or :class:`numpy.ndarray`
            Age of the universe in Gyr.

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
