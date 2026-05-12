"""
fastspecfit.resolution
======================

Class and functions for handling the DESI resolution matrix.

"""
import numpy as np
from numba import jit
import scipy.linalg


# Constants: must match emline_model_core's sigma_eff parameterization
SIGMA0_ANGSTROM   = 0.5   # fiducial pre-convolution Gaussian sigma [Angstrom]
PIX_SIZE_ANGSTROM = 0.8   # DESI pixel size [Angstrom]
_NDIAG            = 11    # DESI standard resolution matrix half-bandwidth


def _make_gaussian_matrix(width, sigma0_angstrom, pix_size_angstrom):
    """
    Build the width x width Gaussian RBF matrix M used in the deconvolution
    W = M^{-1} R (Koposov, rvspecfit).

    """
    sig_pix = sigma0_angstrom / pix_size_angstrom
    xs = np.arange(width)
    return np.array([
        1.0 / np.sqrt(2 * np.pi) / sig_pix *
        np.exp(-0.5 * ((xs - i) / sig_pix)**2)
        for i in range(width)
    ])


# Pre-factorize the constant Gaussian matrix once at module load.
# Amortizes O(width^3) LU cost across all Resolution instantiations.
_GAU_MAT_LU = scipy.linalg.lu_factor(
    _make_gaussian_matrix(_NDIAG, SIGMA0_ANGSTROM, PIX_SIZE_ANGSTROM)
)


@jit(nopython=True, nogil=True, cache=True)
def _mat_torows(mat):
    """
    Convert DESI diagonal-storage matrix to row representation.
    Numba equivalent of Koposov's resolution_mat_torows, replacing
    np.roll + list comprehension with direct index arithmetic.

    For row i of the output:
        result[i, j] = mat[w-1-i, (j - ((w-1-i) - w2)) % npix]

    """
    w, npix = mat.shape
    w2 = w // 2
    result = np.empty((w, npix), dtype=mat.dtype)
    for i in range(w):
        r     = w - 1 - i
        shift = r - w2
        for j in range(npix):
            result[i, j] = mat[r, (j - shift) % npix]
    return result


@jit(nopython=True, nogil=True, cache=True)
def _mat_tocolumns(mat):
    """
    Convert row representation back to DESI diagonal-storage matrix.
    Numba equivalent of Koposov's resolution_mat_tocolumns.

    For row r of the output:
        result[r, j] = mat[r, (j - (r - w2)) % npix]

    """
    w, npix = mat.shape
    w2 = w // 2
    result = np.empty((w, npix), dtype=mat.dtype)
    for r in range(w):
        shift = r - w2
        for j in range(npix):
            result[r, j] = mat[r, (j - shift) % npix]
    return result


def deconvolve_resolution_matrix(mat0,
                                 sigma0_angstrom=SIGMA0_ANGSTROM,
                                 pix_size_angstrom=PIX_SIZE_ANGSTROM):
    """
    Compute W = M^{-1} R: the resolution matrix for spectra pre-convolved
    with a Gaussian of width sigma0_angstrom (Koposov, rvspecfit).

    Parameters
    ----------
    mat0 : np.ndarray [ndiag x npix]
        Raw DESI resolution matrix in diagonal storage format.
    sigma0_angstrom : float
        Fiducial Gaussian sigma [Angstrom]. Must match the sigma_eff
        parameterization used in emline_model_core.
    pix_size_angstrom : float
        Pixel size [Angstrom].

    Returns
    -------
    np.ndarray [ndiag x npix] — deconvolved resolution matrix W.

    """
    width, npix = mat0.shape
    w2 = width // 2

    mat_rows = _mat_torows(mat0)

    # Zero boundary columns where the banded structure breaks down.
    # Edge pixels are masked during fitting so no further correction needed.
    for i in range(w2):
        mat_rows[:w2 - i - 1, i]            = 0.
        mat_rows[w2 + 1 + i:, npix - 1 - i] = 0.

    # Solve M @ W_rows = R_rows using pre-factorized LU when possible.
    if width == _NDIAG:
        mat_rows1 = scipy.linalg.lu_solve(_GAU_MAT_LU, mat_rows)
    else:
        # Fallback for non-standard matrix widths
        gau_mat = _make_gaussian_matrix(width, sigma0_angstrom, pix_size_angstrom)
        mat_rows1 = scipy.linalg.solve(gau_mat, mat_rows)

    return _mat_tocolumns(mat_rows1)


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
        data0 : :class:`np.ndarray` [ndiag x dim]
           Array of values for diagonals [ndiag//2 .. -ndiag//2];
           note that some entries at the ends of rows are ignored.
           All other diagonals are assumed to be zero.

        """
        data = deconvolve_resolution_matrix(data0)
        ndiag, dim  = data.shape
        self.shape  = (dim, dim)
        self.ndiag  = ndiag
        self.data   = data
        self.rows   = None # compute on first use


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
