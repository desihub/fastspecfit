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
    """Build the ``width x width`` Gaussian RBF matrix M for deconvolution W = M^{-1} R."""
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
    """Convert a DESI diagonal-storage matrix to row representation.

    Parameters
    ----------
    mat : :class:`numpy.ndarray`, shape (ndiag, npix)
        DESI resolution matrix in diagonal storage format.

    Returns
    -------
    result : :class:`numpy.ndarray`, shape (ndiag, npix)
        Row representation where ``result[i, j]`` gives the matrix value
        at row ``i``, column ``j``.

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
    """Convert a row-representation matrix back to DESI diagonal-storage format.

    Parameters
    ----------
    mat : :class:`numpy.ndarray`, shape (ndiag, npix)
        Row representation of the matrix.

    Returns
    -------
    result : :class:`numpy.ndarray`, shape (ndiag, npix)
        DESI diagonal-storage format.

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
    """Deconvolve a DESI resolution matrix: compute W = M^{-1} R.

    Produces the resolution matrix appropriate for spectra pre-convolved
    with a Gaussian of width ``sigma0_angstrom`` (Koposov / rvspecfit
    approach).

    Parameters
    ----------
    mat0 : :class:`numpy.ndarray`, shape (ndiag, npix)
        Raw DESI resolution matrix in diagonal storage format.
    sigma0_angstrom : :class:`float`, optional
        Fiducial Gaussian sigma in Angstroms. Must match the ``sigma_eff``
        parameterization used in the emission-line model core. Default is
        :data:`SIGMA0_ANGSTROM`.
    pix_size_angstrom : :class:`float`, optional
        Pixel size in Angstroms. Default is :data:`PIX_SIZE_ANGSTROM`.

    Returns
    -------
    W : :class:`numpy.ndarray`, shape (ndiag, npix)
        Deconvolved resolution matrix in diagonal storage format.

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
    """Resolution matrix modeled after ``desispec.resolution.Resolution``.

    Stores data as a ``(ndiag, dim)`` array of diagonal values, analogous
    to :class:`scipy.sparse.dia_matrix` for diagonals
    ``[ndiag//2, ..., -ndiag//2]`` of a ``dim x dim`` matrix. The
    :meth:`dot` method is implemented in Numba for speed. Also provides
    :meth:`rowdata` for a sparse row representation used by the
    emission-line fitting code.

    Parameters
    ----------
    data0 : :class:`numpy.ndarray`, shape (ndiag, dim)
        Raw DESI resolution matrix in diagonal storage format. Entries at
        the ends of rows that fall outside the ``dim x dim`` matrix are
        ignored; all other diagonals are assumed zero.

    """
    def __init__(self, data0):
        data = deconvolve_resolution_matrix(data0)
        ndiag, dim  = data.shape
        self.shape  = (dim, dim)
        self.ndiag  = ndiag
        self.data   = data
        self.rows   = None # compute on first use


    def rowdata(self):
        """Return the sparse row representation of the resolution matrix.

        Returns
        -------
        rows : :class:`numpy.ndarray`, shape (dim, ndiag)
            Nonzero entries for each row of the resolution matrix.

        """
        if self.rows is None:
            self.rows = self._dia_to_rows(self.data)
        return self.rows


    def dot(self, v, out=None):
        """Multiply the resolution matrix by vector ``v``.

        Parameters
        ----------
        v : :class:`numpy.ndarray`, shape (dim,)
            Input vector.
        out : :class:`numpy.ndarray` or None, optional
            Pre-allocated output array of length ``dim``; allocated if
            ``None``.

        Returns
        -------
        out : :class:`numpy.ndarray`, shape (dim,)
            Result of the matrix-vector product.

        """
        if out is None:
            out = np.empty(self.shape[0], dtype=v.dtype)

        return self._matvec(self.data, v, out)


    def matmat(self, V, out=None):
        """Apply the resolution matrix to all rows of ``V`` (ntemplates × dim).

        Parameters
        ----------
        V : :class:`numpy.ndarray`, shape (ntemplates, dim)
        out : :class:`numpy.ndarray` or None, optional
            Pre-allocated output of the same shape; allocated if ``None``.

        Returns
        -------
        out : :class:`numpy.ndarray`, shape (ntemplates, dim)

        """
        if out is None:
            out = np.empty_like(V)

        return self._matmat(self.data, V, out)


    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True, cache=True)
    def _matvec(D, v, out):
        """Compute the resolution matrix-vector product, writing result into ``out``."""
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
    def _matmat(D, V, out):
        """Apply the resolution matrix to all rows of ``V``, writing results into ``out``."""
        out[:] = 0.

        ntemplates = V.shape[0]
        ndiag, dim = D.shape
        for t in range(ntemplates):
            for i in range(ndiag):
                k = ndiag//2 - i

                i_start = np.maximum(0, -k)
                j_start = np.maximum(0,  k)
                j_end   = np.minimum(dim + k, dim)

                for n in range(j_end - j_start):
                    out[t, i_start + n] += D[i, j_start + n] * V[t, j_start + n]

        return out


    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True, cache=True)
    def _dia_to_rows(D):
        """Convert a DESI diagonal-storage matrix to a sparse row representation."""
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
