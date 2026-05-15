import numpy as np

from numba import jit


class ParamsMapping(object):
    """Map free line parameters to the full EMLine model parameter set.

    Precomputes and caches the mapping from free (unconstrained) parameters
    to the full parameter vector, including fixed, tied, and doublet
    relationships, as well as the Jacobian of that mapping.

    Parameters
    ----------
    nParms : int
        Total number of parameters (free + fixed/tied).
    isFree : :class:`numpy.ndarray` of bool
        Boolean mask indicating which parameters are free.
    tiedSources : :class:`numpy.ndarray` of int
        Index of the source parameter for each tied parameter (-1 if not tied).
    tiedFactors : :class:`numpy.ndarray` of float
        Multiplicative factor from source to target for each tied parameter.
    doubletTargets : :class:`numpy.ndarray` of int
        Indices of parameters that are doublet ratio targets.
    doubletSources : :class:`numpy.ndarray` of int
        Source parameter indices for each doublet ratio.

    """

    def __init__(self, nParms,
                 isFree,
                 tiedSources, tiedFactors,
                 doubletTargets, doubletSources):

        self.nParms     = nParms
        self.nFreeParms = np.sum(isFree)

        # permutation mapping each free parameter in full list
        # to its location in free list
        pFree = np.empty(self.nParms, dtype=np.int32)
        pFree[isFree] = np.arange(self.nFreeParms, dtype=np.int32)

        self._precomputeMapping(isFree, pFree,
                                tiedSources, tiedFactors,
                                doubletTargets, doubletSources)

        self._precomputeJacobian()


    def fixedMask(self):
        """Return a Boolean mask with ``True`` for each fixed parameter.

        Returns
        -------
        mask : :class:`numpy.ndarray` of bool
            ``True`` where parameters are fixed (not free or tied).

        """
        return self.isFixed


    def mapFreeToFull(self, freeParms, out=None, patchDoublets=True):
        """Map a vector of free parameters to the full parameter set.

        Parameters
        ----------
        freeParms : :class:`numpy.ndarray`
            Values of the free parameters.
        out : :class:`numpy.ndarray` or None, optional
            Pre-allocated output array; if ``None``, one is allocated.
        patchDoublets : bool, optional
            If ``True``, replace doublet ratio parameters with their actual
            values. Defaults to ``True``.

        Returns
        -------
        fullParms : :class:`numpy.ndarray`
            Values of all parameters (free, fixed, tied, and doublet).

        """

        return self._mapFreeToFull(freeParms,
                                   self.nParms,
                                   self.sources,
                                   self.factors,
                                   self.doubletPatches,
                                   out,
                                   patchDoublets)


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _mapFreeToFull(freeParms, nParms, sources, factors,
                       doubletPatches, fullParms, patchDoublets):

        for j, src_j_free in doubletPatches:
            factors[j] = freeParms[src_j_free] if patchDoublets else 1.

        if fullParms is None:
            fullParms = np.empty(nParms, dtype=freeParms.dtype)

        for j, src_j_free in enumerate(sources):
            fullParms[j] = factors[j]  # copy fixed zeros
            if src_j_free != -1:
                fullParms[j] *= freeParms[src_j_free]

        return fullParms


    def getJacobian(self, freeParms):
        """Return the Jacobian of the free-to-full parameter mapping.

        Parameters
        ----------
        freeParms : :class:`numpy.ndarray`
            Current values of the free parameters.

        Returns
        -------
        J_S : tuple
            Sparse Jacobian as ``(shape, elts, factors)``, where ``shape``
            is ``(nParms, nFreeParms)``, ``elts`` contains the row and column
            indices of each nonzero element, and ``factors`` contains their
            values.

        """

        for j, src_j_free in self.jacDoubletPatches:
            self.jacFactors[j] = freeParms[src_j_free]

        return ((self.nParms, self.nFreeParms), self.jacElts, self.jacFactors)


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _matvec(J_S, v, w):
        """Multiply sparse parameter Jacobian ``J_S`` by ``v``, writing the result to ``w``.

        Parameters
        ----------
        J_S : tuple
            Sparse Jacobian from :meth:`getJacobian`.
        v : :class:`numpy.ndarray`
            Input vector of free parameters.
        w : :class:`numpy.ndarray`
            Output vector of full parameters (overwritten).

        """

        shape, elts, factors = J_S

        for j in range(shape[0]):  # total params
            w[j] = 0.

        for i, (dst, src) in enumerate(elts):
            w[dst] += factors[i] * v[src]


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _add_rmatvec(J_S, v, w):
        """Multiply ``v`` by the transpose of sparse Jacobian ``J_S``, adding the result to ``w``.

        Parameters
        ----------
        J_S : tuple
            Sparse Jacobian from :meth:`getJacobian`.
        v : :class:`numpy.ndarray`
            Input vector of full parameters.
        w : :class:`numpy.ndarray`
            Output vector of free parameters (accumulated in-place).

        """

        _, elts, factors = J_S

        for i, (dst, src) in enumerate(elts):
            w[src] += factors[i] * v[dst]


    ###########################################################


    def _precomputeMapping(self, isFree, pFree,
                           tiedSources, tiedFactors,
                           doubletTargets, doubletSources):
        """Precompute source/factor arrays and doublet patches for :meth:`mapFreeToFull`."""

        # by default, assume parameters are fixed
        sources = np.full(self.nParms, -1, dtype=np.int32)
        factors = np.zeros(self.nParms, dtype=np.float64)

        # record all free parameters
        sources[isFree] = pFree[isFree]
        factors[isFree] = 1.

        for j, src_j in enumerate(tiedSources):
            if src_j != -1 and isFree[src_j]:
                # j is tied to src_j with factor tiedFactors[j]
                sources[j] = pFree[src_j]
                factors[j] = tiedFactors[j]

        doubletPatches = []
        for j, src_j in zip(doubletTargets, doubletSources):
            if isFree[j] and isFree[src_j]:
                # j's factor should be v[ p[src_j] ], so that its value
                # becomes v[ p[j] ] * v[ p[src_j] ]. We will patch it
                # dynamically at mapping time.

                doubletPatches.append((j, pFree[src_j]))

        self.sources = sources
        self.factors = factors

        # if there are no patches, we need to create a null array of
        # the right type and shape to appease Numba
        if len(doubletPatches) == 0:
            self.doubletPatches = np.empty((0,2), dtype=np.int32)
        else:
            self.doubletPatches = np.array(doubletPatches, dtype=np.int32)

        # record fixed parameter mask
        self.isFixed = (self.sources == -1)


    def _precomputeJacobian(self):
        """Precompute sparse Jacobian structure and doublet patches for :meth:`getJacobian`."""

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

        for i, (j, p_src_j) in enumerate(self.doubletPatches):

            p_j = self.sources[j]

            # jacElts already has coefficient (j, p[j]).
            # its factor should be patched from v[ p[src_j] ]
            jacDoubletPatches[2*i, :] = (liveOffsets[j], p_src_j)

            # add a second coefficient (j, p[src_j]).
            # its factor should be patched from v[ p[j] ]
            jacElts[nLive + i,:] = (j, p_src_j)
            jacDoubletPatches[2*i+1,:] = (nLive + i, p_j)

        self.jacElts = jacElts
        self.jacFactors = jacFactors
        self.jacDoubletPatches = jacDoubletPatches
