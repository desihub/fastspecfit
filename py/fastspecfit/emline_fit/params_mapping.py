import numpy as np

from numba import jit


class ParamsMapping(object):
    """Compute a mapping from the free (line) parameters of a spectrum
    fitting problem to the full set of parameters for the EMLine
    model, as well the Jacobian of this mapping.  We capture most of
    the complexity of the mapping once and do the minimum necessary
    work to update it for each new set of free parameters.

    """

    def __init__(self, nParms,
                 isFree,
                 tiedSources, tiedFactors,
                 doubletTargets, doubletSources):
        """
        Parameters
        ----------
        nParms : :class:`int`
          Total number of parameters (free + fixed/tied).
        isFree : :class:`np.ndarray` of `bool` [nParms]
          Is each parameter free or fixed/tied?
        tiedSources : :class:`np.ndarray` of `int` [nParms]
          For parameters tied to a source param, idx of source
          (-1 if not tied).
        tiedFactors : :class:`np.ndarray` of `np.float64` [nParms]
          For parameters tied to a source param, multiplier from
          source to target.
        doubletTargets: :class:`np.ndarray` of `int`
          Parameters that are doublet ratios vs. a source.
        doubletSources: :class:`np.ndarray` of `int`
          Source parameters for all doublet ratios.

        """

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
        """
        Return Boolean mask with "True" for each fixed parameter.
        """
        return self.isFixed


    def mapFreeToFull(self, freeParms, out=None, patchDoublets=True):
        """
        Given a vector of free parameters, return the corresponding
        list of full parameters, accounting for fixed, tied, and
        (if patchDoublets is true) doublet features.

        Parameters
        ----------
        freeParms : :class:`np.ndarray` of `np.float64`
          Values of free parameters.
        out : :class: `np.ndarray` of `np.float64` or :None:
          Vector to hold values of full parameter set; if none,
          allocate one.
        patchDoublets : :class: `bool`
          Replace doublet target ratios by their actual values?

        Returns
        -------
        A `np.ndarray` with the values of all parameters.

        Notes
        ----
        Accelerated via Numba implementation below.
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
        """
        Given a vector v of free parameters, return the Jacobian
        of the transformation from free to full parameters at v.

        Parameters
        ----------
        freeParms : :class:`np.ndarray` of `np.float64`
          Vector of free parameter values

        Returns
        -------
        Sparse Jacobian as a tuple (shape, elts, factors), where
          - shape is shape of Jacobian [# parms x # free parms]
          - elts gives the 2D indices of nonzero elements of Jacobian
          - factors gives the values of the nonzero elements

        """

        for j, src_j_free in self.jacDoubletPatches:
            self.jacFactors[j] = freeParms[src_j_free]

        return ((self.nParms, self.nFreeParms), self.jacElts, self.jacFactors)


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _matvec(J_S, v, w):
        """
        Multiply sparse parameter Jacobian J_S * v, writing result to w.

        Parameters
        ----------
        J_S : :class: `tuple`:
          Sparse Jacobian computed by getJacobian()
        v : :class: `np.ndarray` [# free parameters]
        w : :class: `np.ndarray` [# total parameters] (output)

        """

        shape, elts, factors = J_S

        for j in range(shape[0]):  # total params
            w[j] = 0.

        for i, (dst, src) in enumerate(elts):
            w[dst] += factors[i] * v[src]


    @staticmethod
    @jit(nopython=True, nogil=True, cache=True)
    def _add_rmatvec(J_S, v, w):
        """
        Multiply v with sparse parameter Jacobian J_S.T,
        *adding* result to w.

        Parameters
        ----------
        J_S : :class: `tuple`:
          Sparse Jacobian computed by getJacobian()
        v : :class: `np.ndarray` [# total parameters]
        w : :class: `np.ndarray` [# free parameters] (output)

        """

        _, elts, factors = J_S

        for i, (dst, src) in enumerate(elts):
            w[src] += factors[i] * v[dst]


    ###########################################################


    def _precomputeMapping(self, isFree, pFree,
                           tiedSources, tiedFactors,
                           doubletTargets, doubletSources):
        """
        Precompute and store all the transformations from free
        parameters to full parameters that do not require knowledge of
        the free parameter values.  Create "patches" for doublets so
        that we can quickly compute the values of all doublet targets
        once their source values and ratios are known.

        isFree : :class:`np.ndarray` of `bool` [nParms]
          Is each parameter free or fixed/tied?
        pFree: :class:`np.ndarray` of `int` [nParms]
          If a parameter in the full list isFree, which
          index in the free parameter list does it come from?
          (If not free, entry is undefined)
        tiedSources : :class:`np.ndarray` of `int` [nParms]
          For parameters tied to a source param, idx of source
          (-1 if not tied).
        tiedFactors : :class:`np.ndarray` of `np.float64` [nParms]
          For parameters tied to a source param, multiplier from
          source to target.
        doubletTargets: :class:`np.ndarray` of `int`
          Parameters that are doublet ratios vs. a source.
        doubletSources: :class:`np.ndarray` of `int`
          Source parameters for all doublet ratios.

        """

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
        """
        Precompute and store as much of the Jacobian of the
        transformation from free parameters to full parameters as does
        not require knowledge of the free parameter values.  We rely
        on the precomputation for the mapping to avoid recalculating
        all the source/factor information.  We transform doublet
        patches for the mapping into a larger set of patches for the
        Jacobian.

        """

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
