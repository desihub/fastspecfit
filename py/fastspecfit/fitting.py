"""
fastspecfit.fitting
===================

Methods and classes for optimization routines.

Includes sparse representations of resolution matrices and Jacobian of objective
for emission line fitting. (NB: The ideal Jacobian rep is generated in
EMLines_jacobian().)

"""
import numpy as np
from numba import jit
from scipy.sparse.linalg import LinearOperator

class ParamsMapping(object):
    """Compute a mapping from the free parameters of a spectrum fitting problem to
    the full set of parameters for the EMLine model, as well the Jacobian of
    this mapping.  We capture most of the complexity of the mapping once and do
    the minimum necessary updating to it for each new set of freeParameters.

    Written by Jeremy Buhler (https://github.com/jdbuhler) - June 2024.

    """
    def __init__(self, fixedParameters, freeParms,
                 tiedParms, tiedSources, tiedFactors,
                 doubletRatios, doubletSources):
        
        self.nParms     = len(fixedParameters)
        self.nFreeParms = len(freeParms)
        
        # permutation mapping each free parameter in full list
        # to its location in free list
        pFree = np.empty(self.nParms, dtype=np.int32)
        pFree[freeParms] = np.arange(self.nFreeParms, dtype=np.int32)
        
        self._precomputeMapping(fixedParameters, freeParms,
                                tiedParms, tiedSources, tiedFactors,
                                doubletRatios, doubletSources,
                                pFree)

        self._precomputeJacobian()

        

    # return Boolean mask with "True" for each fixed parameter
    def fixedMask(self):
        return self.isFixed

    
    #
    # mapFreeToFull()
    # Given a vector of free parameters, return the corresponding
    # list of full parameters, accounting for fixed, tied, and
    # doublet features.
    #
    def mapFreeToFull(self, freeParms):

        return self._mapFreeToFull(freeParms,
                                   self.nParms,
                                   self.sources,
                                   self.factors,
                                   self.doubletPatches)

    #
    # _mapFreeToFull()
    # Given a vector of free parameters, return the corresponding
    # list of full parameters, accounting for fixed, tied, and
    # doublet features.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _mapFreeToFull(freeParms, nParms, sources, factors, doubletPatches):
        
        for j, src_j_free in doubletPatches:
            factors[j] = freeParms[src_j_free]

        fullParms = np.empty(nParms, dtype=freeParms.dtype)

        for j, src_j_free in enumerate(sources):
            fullParms[j] = factors[j] # copy fixed value
            if src_j_free != -1:
                fullParms[j] *= freeParms[src_j_free]
        
        return fullParms
        
    
    #
    # getJacobian()
    # Given a vector v of free parameters, return the Jacobian
    # of the transformation from free to full at v.  The Jacobian
    # is a sparse matrix represented as an array of nonzero entries
    # (jacElts) and their values (jacFactors)
    #
    def getJacobian(self, freeParms):
        for j, src_j_free in self.jacDoubletPatches:
            self.jacFactors[j] = freeParms[src_j_free]

        return ((self.nParms, self.nFreeParms), self.jacElts, self.jacFactors)
    

    #
    # Multiply parameter Jacobian J_S * v, writing result to w.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _matvec(J_S, v, w):

        shape, elts, factors = J_S
        
        for j in range(shape[0]): # total params
            w[j] = 0.
        
        for i, (dst, src) in enumerate(elts):
            w[dst] += factors[i] * v[src]
            


    #
    # Multiply parameter Jacobian v * J_S^T, *adding* result to w.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _add_rmatvec(J_S, v, w):

        _, elts, factors = J_S
        
        for i, (dst, src) in enumerate(elts):
            w[src] += factors[i] * v[dst]

    
    ###########################################################
    
    #
    # Precompute all the transformations from free parameters
    # to full parameters that do not require knowledge of the
    # free parameter values.
    #
    def _precomputeMapping(self, fixedParameters, freeParms,
                           tiedParms, tiedSources, tiedFactors,
                           doubletRatios, doubletSources,
                           p):
        
        # by default, assume parameters are fixed and that
        # they take on the values in fixedParameters
        sources = np.full(self.nParms, -1, dtype=np.int32)
        factors = fixedParameters.copy()

        for j in freeParms:
            sources[j] = p[j]
            factors[j] = 1.

        for j, src_j, factor in zip(tiedParms, tiedSources, tiedFactors):
            if src_j not in freeParms:
                #print(f"SOURCE {src_j} tied to {j} is not free!")
                # if source is fixed, so is target, and it's in fixedParameters
                pass
            else:
                sources[j] = p[src_j]
                factors[j] = factor

        doubletPatches = []
        for (j, src_j) in zip(doubletRatios, doubletSources):
            #if j not in freeParms:
            #    print(f"ratio {j} in doublet with {src_j} is not free!")
            #if src_j not in freeParms:
            #    print(f"amplitude {src_j} in doublet with {j} is not free!")

            if j not in freeParms or src_j not in freeParms:
                continue

            # j's factor should be v[ p[j] ], so that its value
            # becomes v[ p[j] ] * v[ p[src_j] ]. We will patch it
            # dynamically at mapping time.
            
            sources[j] = p[src_j]
            doubletPatches.append((j, p[j])) # record where to grab factor
        
        self.sources = sources
        self.factors = factors
        self.doubletPatches = np.array(doubletPatches)

        # record fixed parameter mask
        self.isFixed = (self.sources == -1)

        
    #
    # Precompute as much of the Jacobian of the transformation
    # from free parameters to full parameters as does not require
    # knowledge of the free parameter values.  We rely on the
    # precomputation for the mapping to avoid recalculating all
    # the source/factor information.
    #
    def _precomputeJacobian(self):
        
        isLive  = np.logical_not(self.isFixed)
        liveIdx = np.where(isLive)[0]
        nLive = len(liveIdx)
        
        # jacElts compresses the source/factor arrays to a sparse array
        # of coefficients for just the live (non-fixed)
        # parameters, plus second coeffs for each ratio param of a
        # doublet pair.  We need not create entries for fixed parametesr
        # because their derivatives w/r to the free parameters are zero.
        nDoublets = len(self.doubletPatches)
        nElts = nLive + nDoublets
        jacElts = np.empty((nElts,2), dtype=np.int32)
        jacFactors = np.empty(nElts, dtype=self.factors.dtype)
        
        jacElts[:nLive,0] = liveIdx
        jacElts[:nLive,1] = self.sources[liveIdx]
        jacFactors[:nLive] = self.factors[liveIdx]  # makes a copy
        
        # Create doublet patches for Jacobian w/r to compressed array,
        # and record new patch locations
        
        # offsets of orig parms in live array
        liveOffsets = np.cumsum(isLive) - isLive
        
        jacDoubletPatches = np.empty((2*len(self.doubletPatches), 2), dtype=np.int32)
        
        for i, (j, p_j) in enumerate(self.doubletPatches):
            
            p_src_j = self.sources[j]

            # jacElts already has coefficient (j, p[src_j]).
            # its factor should be patched from v[ p[j] ]
            jacDoubletPatches[2*i,  :] = (liveOffsets[j], p_j)
            
            # add a second coefficient (j, p[j]).
            # its factor should be patched from v[ p[src_j] ]
            jacElts[nLive + i,:] = (j, p_j)
            jacDoubletPatches[2*i+1,:] = (nLive + i, p_src_j)
        
        self.jacElts = jacElts
        self.jacFactors = jacFactors
        self.jacDoubletPatches = jacDoubletPatches


class ResMatrix(object):
    """Representation of a resolution matrix.

    A resolution matrix M of size nrow x nrow is stored as a 2D array A of size
    nrow x ndiag, where ndiag is the number of nonzero diagonals (which must be
    odd).  The rows of A are M's rows, but with only the nonzero entries stored.
    The nonzero entries on row i run from j = i - diag//2 to i + diag//2, so

      M[i,j] = A[i, j - (i - diag//2)]

    """
    def __init__(self, D):
        self.data = self._from_dia_matrix(D)

    def matvec(self, v, w):
        self._matvec(self.data, v, w)
        
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
    #   number of parameters in full set
    #   array of start and end obs bin indices
    #     for each camera
    #   partial Jacobian jac = (W * M  J_I)
    #     for each camera
    #   parameter expansion Jacobian J_S
    #
    def __init__(self, shape, nParms, camerapix, jacs, J_S):
        
        self.camerapix = camerapix
        self.jacs      = tuple(jacs)
        self.J_S       = J_S
        
        dtype = jacs[0][1].dtype

        self.vFull = np.empty(nParms, dtype=dtype)
        
        super().__init__(dtype, shape)


    #
    # Compute left matrix-vector product J * v
    # |v| = number of free parameters
    #
    def _matvec(self, v):
        
        nBins = self.shape[0]
        w = np.empty(nBins, dtype=v.dtype)
        
        # return result in self.vFull
        ParamsMapping._matvec(self.J_S, v, self.vFull)
        
        for campix, jac in zip(self.camerapix, self.jacs):
            s, e = campix
            
            # write result to w[s:e]
            self._matvec_J(jac, self.vFull, w[s:e])
        
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
            
            # return result in self.vFull
            ParamsMapping._matvec(self.J_S, M[:,i].ravel(), self.vFull)

            for campix, jac in zip(self.camerapix, self.jacs):
                s, e = campix

                # write result to w[s:e]
                self._matvec_J(jac, self.vFull, w[s:e])
                
        return R.T
    
        
    #
    # Compute right matrix product product v * J^T
    # |v| = number of observable bins
    #
    def _rmatvec(self, v):

        nFreeParms = self.shape[1]
        w = np.zeros(nFreeParms, dtype=v.dtype)
        
        for campix, jac in zip(self.camerapix, self.jacs):
            s, e = campix
            
            # return result in self.vFull
            self._rmatvec_J(jac, v[s:e], self.vFull)
            
            # add result to w
            ParamsMapping._add_rmatvec(self.J_S, self.vFull, w)
            
        return w

    #
    # Multiply ideal Jacobian J * v, writing result to w.
    #
    @staticmethod
    @jit(nopython=True, fastmath=False, nogil=True)
    def _matvec_J(J, v, w):
    
        endpts, values = J
        nvars = endpts.shape[0]
        nbins = len(w)
        
        for j in range(nbins):
            w[j] = 0.
        
        for i in range(nvars):
            s, e = endpts[i]
            vals = values[i]    # column i
            
            for j in range(e - s):
                w[j + s] += vals[j] * v[i]  
    
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
