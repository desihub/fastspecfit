#
# Sparse representations of resolution matrices and
# Jacobian of objective for emission line fitting
# (Ideal Jacobian rep is generated in EMLines_jacobian())
#

import numpy as np
from scipy.sparse.linalg import LinearOperator

from numba import jit

from .params_mapping import ParamsMapping

    
#
# resolution matrix
# A resolution matrix M of size nrow x nrow is stored as a 2D array A of
# size nrow x ndiag, where ndiag is the number of nonzero diagonals
# (which must be odd).  The rows of A are M's rows, but with only the
# nonzero entries stored.  The nonzero entries on row i run from
# j = i - diag//2 to i + diag//2, so
#            M[i,j] = A[i, j - (i - diag//2)]
#
class ResMatrix(object):

    def __init__(self, D):
        self.data = self._from_dia_matrix(D)
    
    def ndiag(self):
        return self.data.shape[1]
    
    def matvec(self, v, w):
        self._matvec(self.data, v, w)

    # for compatibility
    def dot(self, v):
        w = np.empty(self.data.shape[0])
        self.matvec(v, w)
        return w
    
    #
    # _from_dia_matrix()
    # Convert a diagonally sparse matrix M in the form
    # stored by DESI into a sparse row rerpesentation.
    #
    # Input M is represented as a 2D array D of size ndiag x nrow,
    # whose rows are M's diagonals:
    #            M[i,j] = D[ndiag//2 - (j - i), j]
    # ndiag is assumed to be odd, and entries in D that would be
    # outside the bounds of M are ignored.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _from_dia_matrix(D):
        ndiag, nrow = D.shape
        hdiag = ndiag//2

        A = np.empty((nrow, ndiag), dtype=D.dtype)
        
        for i in range(nrow):
            # min and max column for row
            jmin = np.maximum(i - hdiag,        0)
            jmax = np.minimum(i + hdiag, nrow - 1)

            for j in range(jmin, jmax + 1):
                A[i, j - i + hdiag] = D[hdiag + i - j, j]
                
        return A
    
    #
    # _matvec()
    # Compute the matrix-vector product M * v, where
    # M is a row-sparse matrix with a limited
    # number of diagonals created by
    # dia_to_row_matrix().
    #
    # w is an output parameter
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _matvec(M, v, w):
        nrow, ndiag = M.shape
        hdiag = ndiag//2
        
        for i in range(nrow):
            jmin = np.maximum(i - hdiag,    0)
            jmax = np.minimum(i + hdiag, nrow - 1)
            
            acc = 0.
            for j in range(jmin, jmax + 1):
                acc += M[i, j - i + hdiag] * v[j]
        
            w[i] = acc


#################################################################

#
# Sparse Jacobian of objective function.  For
# each camera's pixel range, the Jacobian
# is a matrix product
#
#    W * M * J_I * J_S
#
# where
#  J_S is the Jacobian of the parameter expansion
#  J_I is the ideal Jacobian
#  M is the camera's resolution matrix
#  w is a diagonal matrix of per-observation weights
#
# Note that we precompute W*M*J_I for each camera 
# with the external mulWMJ() function.  This product
# has one contiguous run of nonzero entries per column,
# while J_S has either one or two nonzero entries
# per row.
#
class EMLineJacobian(LinearOperator):
    
    #
    # CONSTRUCTOR ARGS:
    #   shape of Jacobian
    #   number of *line*-related free parameters in set
    #   array of start and end obs bin indices
    #     for each camera
    #   partial Jacobian jac = (W * M  J_I)
    #     for each camera
    #   parameter expansion Jacobian J_S
    #
    def __init__(self, shape, nLineFreeParms, camerapix, jacs, J_S,
                 J_P = None):
        
        self.camerapix = camerapix
        self.jacs      = tuple(jacs)
        self.J_S       = J_S
        self.J_P       = J_P
        self.nLineFreeParms = nLineFreeParms
        
        if J_P is not None:
            nPatchParms = J_P[0].shape[0]
        
        dtype = jacs[0][1].dtype
        nParms = jacs[0][1].shape[0] # num line params in full set
        self.vFull = np.empty(nParms, dtype=dtype)
        
        super().__init__(dtype, shape)


    #
    # Compute left matrix-vector product J * v
    # |v| = number of free parameters
    #
    def _matvec(self, v):
        
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


    #
    # Compute left matrix-matrix product J * M
    # We do this just to avoid needing to ravel
    # v in the more commonly used _matvec()
    #
    def _matmat(self, M):

        nBins = self.shape[0]
        nVecs = M.shape[1]
        R = np.empty((nVecs, nBins), dtype=M.dtype) # transpose of result

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
    
        
    #
    # Compute right matrix product product v * J^T
    # |v| = number of observable bins
    #
    def _rmatvec(self, v):

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
    # col_norms()
    # Return the column norms of the Jacobian
    #
    def col_norms(self):
        
        norms = np.empty(self.shape[1], dtype=self.jacs[0][1].dtype)
        
        v = np.zeros(self.shape[0])
        for i in range(self.shape[1]):

            v[i] = 1.
            
            w = self._matvec(v)
            norms[i] = np.sqrt(np.sum(w * w))
            
            v[i] = 0.
            
        return norms
    

    #
    # Multiply ideal Jacobian J * v, writing result to w.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _matvec_J(J, v, w):

        nbins = len(w)
        for j in range(nbins):
            w[j] = 0.
            
        _matvec_J_add(J, v, w)
        
    #
    # Multiply ideal Jacobian v * J^T, writing result to w.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _rmatvec_J(J, v, w):
    
        endpts, values = J
        nvars = endpts.shape[0]
    
        for i in range(nvars):
            s, e = endpts[i]
            vals = values[i]   # row i of transpose
            
            acc = 0.
            for j in range(e - s):
                acc += vals[j] * v[j + s]
            w[i] = acc

        
    #
    # for debugging -- compute singular values of Jacobian
    # and estimate its condition number.
    #
    def estimateConditionNumber(self):

        from scipi.sparse.linalg import svds

        try:
            svs = svds(self, return_singular_vectors=False,
                       k=nFreeParms - 1, which="LM")
            sv0 = svds(self, return_singular_vectors=False,
                   k=1, which="SM")[0]
            cond = svs[-1] / sv0
            print(np.hstack((sv0, svs)))       
            print(f"cond(J) = {cond:.3e}")
        
        except:
            print("Failed to compute Jacobian condition number")


#
# Multiply ideal Jacobian J * v, *adding* result to w.  We have to
# make this a standalone function, rather than a static method of the
# EMLine_Jacobian class, because Numba cannot handle calling a static
# class method from JIT'd code.
#
@jit(nopython=True, fastmath=False, nogil=True)
def _matvec_J_add(J, v, w):
    
    endpts, values = J
    nvars = endpts.shape[0]
    
    for i in range(nvars):
        s, e = endpts[i]
        vals = values[i]    # column i
        
        for j in range(e - s):
            w[j + s] += vals[j] * v[i]  
