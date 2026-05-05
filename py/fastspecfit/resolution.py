import numpy as np
from numba import jit
import scipy


def resolution_mat_torows(mat):
    """Code from rvspecfit by Sergey Koposov"""
    w = mat.shape[0]
    w2 = w // 2
    return np.array([np.roll(mat[_], _ - w2) for _ in range(w)])[::-1]


def resolution_mat_tocolumns(mat):
    """Code from rvspecfit by Sergey Koposov"""
    w = mat.shape[0]
    w2 = w // 2
    return np.array([np.roll(mat[::-1][_], w2 - _) for _ in range(w)])


def deconvolve_resolution_matrix(mat0,
                                 sigma0_angstrom=0.5,
                                 pix_size_angstrom=0.8):
    """Code from rvspecfit by Sergey Koposov"""
    width, npix = mat0.shape
    sig_pix = sigma0_angstrom / pix_size_angstrom
    xs = np.arange(width)
    gau_mat = np.array([
        1. / np.sqrt(2 * np.pi) / sig_pix *
        np.exp(-0.5 * ((xs - i) / sig_pix)**2) for i in range(len(xs))
    ])
    w2 = width // 2
    mat_rows = resolution_mat_torows(mat0)
    # Now this stores rows of the banded matrix rather than columns
    for i in range(w2):
        mat_rows[:w2 - i - 1, i] = 0
        # mat_rows[:, i] = mat_rows[:, i] / mat_rows[:, i].sum()
        j = npix - 1 - i
        mat_rows[w2 + 1 + i:, j] = 0
        # mat_rows[:, j] = mat_rows[:, j] / mat_rows[:, j].sum()

    mat_rows1 = scipy.linalg.solve(gau_mat, mat_rows)
    # shift backward
    mat = resolution_mat_tocolumns(mat_rows1)
    # this a hack no deconvolution of edges
    # mat[:, :w2] = mat0[:, :w2]
    # mat[:, -w2:] = mat0[:, -w2:]
    return mat


class Resolution(object):
    """Resolution matrix, in the style of desispec.resolution.Resolution

    The base implementation is analogous to scipy.sparse.dia_matrix,
    storing a [ndiag x dim] 2D array giving the values on diagonals
    [ndiag//2 ... -ndiag//2] of a dim x dim matrix.  We reimplement
    the dot() method in Numba, special-casing for this specific
    pattern of diagonal offsets; the result is faster than using
    dia_matrix, even though the latter's dot method is in C.

    We also provide conversion to a 2D array giving the nonzero
    entries in each *row* of the matrix, which is used by emline
    fitting code.

    """

    def __init__(self, data0):
        """
        Parameters
        ----------
        data : :class:`np.ndarray` [ndiag x dim]
           Array of values for diagonals [ndiag//2 .. -ndiag//2];
           note that some entries at the ends of rows are ignored.
           All other diagonals are assumed to be zero.

        """
        data = deconvolve_resolution_matrix(data0)
        ndiag, dim = data.shape
        self.shape = (dim, dim)
        self.ndiag = ndiag
        self.data  = data

        self.rows  = None # compute on first use


    def rowdata(self):
        """
        Compute sparse row matrix from diagonal representation.

        Returns
        -------
        2D array of size [dim x ndiag] with nonzero entries
        for each row of the matrix.

        """
        if self.rows is None:
            self.rows = self._dia_to_rows(self.data)
        return self.rows


    def dot(self, v, out=None):
        """
        Compute our dot product with vector v, storing
        result in out (which is created if None).

        Parameters
        ----------
        v : :class:`np.ndarray` [dim]
           vector of values to dot with us.
        out : :class:`np.ndarray` [dim] or None
           array to hold output; created if None.

        Returns
        -------
        array of length [dim] holding output.

        """
        if out is None:
            out = np.empty(self.shape[0], dtype=v.dtype)

        return self._matvec(self.data, v, out)


    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True, cache=True)
    def _matvec(D, v, out):
        """
        Compute matrix-vector product of a Resolution matrix A and a vector v.

        A has shape dim x dim, and its nonzero entries are given by
        a matrix D of size ndiag x dim, providing values on diagonals
        [ndiag//2 .. -ndiag//2], inclusive. v is an array of length dim.

        Return result in array out of size dim.

        """
        out[:] = 0.

        ndiag, dim = D.shape
        for i in range(ndiag):
            k = ndiag//2 - i   # diagonal offset

            i_start = np.maximum(0, -k)
            j_start = np.maximum(0,  k)
            j_end   = np.minimum(dim + k, dim)

            for n in range(j_end - j_start):
                out[i_start + n] += D[i, j_start + n] * v[j_start + n]

        return out


    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True, cache=True)
    def _dia_to_rows(D):
        """
        Convert a diagonally sparse matrix M in the form
        stored by DESI into a sparse row representation.

        Input M is represented as a 2D array D of size ndiag x nrow,
        whose rows are M's diagonals:
             M[i,j] = D[ndiag//2 - (j - i), j]
        ndiag is assumed to be odd, and entries in D that would be
        outside the bounds of M are ignored.
        """

        ndiag, dim = D.shape
        hdiag = ndiag//2

        A = np.empty((dim, ndiag), dtype=D.dtype)

        for i in range(dim):
            # min and max column for row
            jmin = np.maximum(i - hdiag,       0)
            jmax = np.minimum(i + hdiag, dim - 1)

            for j in range(jmin, jmax + 1):
                A[i, j - (i - hdiag)] = D[hdiag + i - j, j]

        return A
