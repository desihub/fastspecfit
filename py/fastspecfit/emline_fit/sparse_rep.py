import numpy as np
from scipy.sparse.linalg import LinearOperator

from numba import jit

from .params_mapping import ParamsMapping


class EMLineJacobian(LinearOperator):
    """Sparse Jacobian of objective function for emission-line
    fitting. For each camera's pixel range, the Jacobian is a matrix
    product

      W * M * J_I * J_S

    where
      J_S is the Jacobian of the parameter expansion
      J_I is the ideal Jacobian for the Gaussian mixture
      M is the camera's resolution matrix
      W is a diagonal matrix of per-observation weights

    Note that we precompute W*M*J_I for each camera
    with the external mulWMJ() function.  This product
    has one contiguous run of nonzero entries per column,
    while J_S has either one or two nonzero entries
    per row.

    When we use patch pedestals, the above Jacobian is
    extended with additional rows corresponding to the
    two free parameters (slope, intercept) of each patch.

    The components of this Jacobian are computed in params_mapping.py
    for J_S and jacobian.py for everything else.

    """

    def __init__(self, shape, nLineFreeParms, camerapix, jacs,
                 J_S, J_P=None):
        """
        Parameters
        ----------
        shape : :class:`tuple` [2]
          Shape of Jacobian matrix.
        nLineFreeParms : :class:`int`
          Number of *line*-related free parameters in set.
        camerapix : :class:`np.ndarray` [# cameras]
          Array of start and end wavelength bin indices per camera.
        jacs : :class:`tuple` [# cameras]
          Tuple of partial Jacobians (W * M * J_I) for each camera.
        J_S : :class:`np.ndarray` [# line params x # nLineFreeParms]
          Parameter expansion Jacobian.
        J_P : :class:`np.ndarray` [# patch params x # wavelength bins] or :None:
          Jacobian of patch pedestals, if used.

        """
        self.camerapix = camerapix
        self.jacs      = tuple(jacs)
        self.J_S       = J_S
        self.J_P       = J_P
        self.nLineFreeParms = nLineFreeParms

        if J_P is not None:
            nPatchParms = J_P[0].shape[0]

        dtype = jacs[0][1].dtype
        nParms = jacs[0][1].shape[0]  # num line params in full set
        self.vFull = np.empty(nParms, dtype=dtype)

        super().__init__(dtype, shape)


    def _matvec(self, v):
        """
        Compute left matrix-vector product J * v, where we are J.

        Parameters
        ----------
        v : :class:`np.ndarray` [# of free parameters]

        Returns
        -------
        w : :class:`np.ndarray [# of wavelength bins]

        """

        nBins = self.shape[0]
        w = np.empty(nBins, dtype=v.dtype)

        vLines = v[:self.nLineFreeParms]

        # return result in self.vFull
        ParamsMapping._matvec(self.J_S, vLines, self.vFull)

        for campix, jac in zip(self.camerapix, self.jacs):
            s, e = campix

            # write result to w[s:e]
            self._matvec_J(jac, self.vFull, w[s:e])

        # add contribution of patch pedestals, if any
        if self.J_P is not None:
            vPatches = v[self.nLineFreeParms:]
            _matvec_J_add(self.J_P, vPatches, w)

        return w


    def _matmat(self, M):
        """
        Compute left matrix-matrix product J * M, where we are J.
        This exists mainly to avoid having to ravel v in the
        more frequently called _matvec().  M has an arbitrary
        second dimension N.

        Parameters
        ----------
        M : :class:`np.ndarray` [# of free parameters x N]

        Returns
        -------
        w : :class:`np.ndarray [N x # of wavelength bins]

        """

        nBins = self.shape[0]
        nVecs = M.shape[1]
        R = np.empty((nVecs, nBins), dtype=M.dtype)  # transpose of result

        for i in range(nVecs):
            w = R[i,:]

            v = M[:,i].ravel()

            # return result in self.vFull
            vLines = v[:self.nLineFreeParms]
            ParamsMapping._matvec(self.J_S, vLines, self.vFull)

            for campix, jac in zip(self.camerapix, self.jacs):
                s, e = campix

                # write result to w[s:e]
                self._matvec_J(jac, self.vFull, w[s:e])

            # add contribution of patch pedestals, if any
            if self.J_P is not None:
                vPatches = v[self.nLineFreeParms:]
                _matvec_J_add(self.J_P, vPatches, w)

        return R.T


    def _rmatvec(self, v):
        """
        Compute right matrix-vector product v * J.T, where we are J.

        Parameters
        ----------
        v : :class:`np.ndarray` [# of wavelength bins]

        Returns
        -------
        w : :class:`np.ndarray [# of free parameters]

        """

        nFreeParms = self.shape[1]
        w = np.zeros(nFreeParms, dtype=v.dtype)

        wLines = w[:self.nLineFreeParms]
        for campix, jac in zip(self.camerapix, self.jacs):
            s, e = campix

            # return result in self.vFull
            self._rmatvec_J(jac, v[s:e], self.vFull)

            # add result to w
            ParamsMapping._add_rmatvec(self.J_S, self.vFull, wLines)

        if self.J_P is not None:
            wPatches = w[self.nLineFreeParms:]
            self._rmatvec_J(self.J_P, v, wPatches)

        return w


    #
    # Multiply ideal Jacobian J * v, writing result to w.
    #
    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _matvec_J(J, v, w):
        """
        Multiply partial Jacobian J with v,
        returning result in supplied w.

        Parameters
        ----------
        J : :class:`np.ndarray` [# wavelength bins x # line parameters]
        v : :class:`np.ndarray` [# of line parameters]
        w : :class:`np.ndarray` [# wavength bins] (output param)

        """

        nbins = len(w)
        for j in range(nbins):
            w[j] = 0.

        _matvec_J_add(J, v, w)


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _rmatvec_J(J, v, w):
        """
        Multiply v with partial Jacobian J.T,
        returning result in supplied w.

        Parameters
        ----------
        J : :class:`np.ndarray` [# wavelength bins x # line parameters]
        v : :class:`np.ndarray` [# wavength bins]
        w : :class:`np.ndarray` [# of line parameters] (output param)

        """

        endpts, values = J
        nvars = endpts.shape[0]

        for i in range(nvars):
            s, e = endpts[i]
            vals = values[i]   # row i of transpose

            acc = 0.
            for j in range(e - s):
                acc += vals[j] * v[j + s]
            w[i] = acc


@jit(nopython=True, nogil=True, cache=True)
def _matvec_J_add(J, v, w):
    """
    Multiply partial Jacobian J with v,
    *adding* result to supplied w.

    Parameters
    ----------
    J : :class:`np.ndarray` [# wavelength bins x # line parameters]
    v : :class:`np.ndarray` [# of line parameters]
    w : :class:`np.ndarray` [# wavength bins] (output param)

    Notes
    -----
    This function is standalone, rather than a static method of
    the EMLineJacobian class, only because Numba cannot call a static
    class method from JIT'd code.

    """

    endpts, values = J
    nvars = endpts.shape[0]

    for i in range(nvars):
        s, e = endpts[i]
        vals = values[i]    # column i

        for j in range(e - s):
            w[j + s] += vals[j] * v[i]
