import numpy as np
from scipy.sparse.linalg import LinearOperator

from numba import jit

from .params_mapping import ParamsMapping


class EMLineJacobian(LinearOperator):
    """Sparse Jacobian of the emission-line fitting objective function.

    For each camera's pixel range, the Jacobian is a matrix product
    ``W * M * J_I * J_S``, where ``J_S`` is the Jacobian of the
    parameter expansion, ``J_I`` is the ideal Jacobian for the Gaussian
    mixture, ``M`` is the camera's resolution matrix, and ``W`` is a
    diagonal matrix of per-observation weights.

    The product ``W * M * J_I`` is precomputed for each camera.  When
    patch pedestals are used, the Jacobian is extended with additional
    columns for the slope and intercept of each patch.

    Parameters
    ----------
    shape : tuple
        Shape ``(nrows, ncols)`` of the Jacobian matrix.
    nLineFreeParms : int
        Number of free line parameters.
    camerapix : :class:`numpy.ndarray`
        Start and end wavelength bin indices for each camera.
    jacs : tuple
        Precomputed partial Jacobians ``(W * M * J_I)`` for each camera.
    J_S : :class:`numpy.ndarray`
        Parameter expansion Jacobian mapping free to full line parameters.
    J_P : :class:`numpy.ndarray` or None, optional
        Sparse Jacobian of patch pedestals, or ``None`` if not used.

    """

    def __init__(self, shape, nLineFreeParms, camerapix, jacs,
                 J_S, J_P=None):
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
        """Compute the left matrix-vector product ``J @ v``.

        Parameters
        ----------
        v : :class:`numpy.ndarray`
            Vector of free parameters.

        Returns
        -------
        w : :class:`numpy.ndarray`
            Result vector over all wavelength bins.

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
        """Compute the left matrix-matrix product ``J @ M``.

        Parameters
        ----------
        M : :class:`numpy.ndarray`, shape (nfreeparms, N)
            Matrix whose columns are free-parameter vectors.

        Returns
        -------
        result : :class:`numpy.ndarray`, shape (nbins, N)
            Product of the Jacobian with each column of ``M``.

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
        """Compute the right matrix-vector product ``J.T @ v``.

        Parameters
        ----------
        v : :class:`numpy.ndarray`
            Vector over all wavelength bins.

        Returns
        -------
        w : :class:`numpy.ndarray`
            Result vector of free parameters.

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
        """Multiply sparse Jacobian ``J`` by ``v``, writing the result to ``w``.

        Parameters
        ----------
        J : tuple
            Sparse Jacobian in ``(endpts, values)`` form.
        v : :class:`numpy.ndarray`
            Input vector of line parameters.
        w : :class:`numpy.ndarray`
            Output vector of wavelength bins (overwritten).

        """

        nbins = len(w)
        for j in range(nbins):
            w[j] = 0.

        _matvec_J_add(J, v, w)


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _rmatvec_J(J, v, w):
        """Multiply ``v`` by the transpose of sparse Jacobian ``J``, writing the result to ``w``.

        Parameters
        ----------
        J : tuple
            Sparse Jacobian in ``(endpts, values)`` form.
        v : :class:`numpy.ndarray`
            Input vector of wavelength bins.
        w : :class:`numpy.ndarray`
            Output vector of line parameters (overwritten).

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
    """Multiply sparse Jacobian ``J`` by ``v``, adding the result to ``w``.

    Parameters
    ----------
    J : tuple
        Sparse Jacobian in ``(endpts, values)`` form.
    v : :class:`numpy.ndarray`
        Input vector of line parameters.
    w : :class:`numpy.ndarray`
        Output vector of wavelength bins (accumulated in-place).

    Notes
    -----
    This function is module-level rather than a static method of
    :class:`EMLineJacobian` because Numba cannot call a static class
    method from JIT-compiled code.

    """

    endpts, values = J
    nvars = endpts.shape[0]

    for i in range(nvars):
        s, e = endpts[i]
        vals = values[i]    # column i

        for j in range(e - s):
            w[j + s] += vals[j] * v[i]
