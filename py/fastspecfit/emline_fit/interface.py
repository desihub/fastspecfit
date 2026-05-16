import numpy as np
from math import erf, erfc

from numba import jit

from .params_mapping import ParamsMapping
from .sparse_rep import EMLineJacobian

from .utils import (
    MAX_SDEV,
    norm_pdf,
    norm_cdf,
    max_buffer_width,
)

from .model import (
    emline_model,
    emline_perline_models,
)

from .jacobian import (
    emline_model_jacobian,
    patch_jacobian,
)


class EMLine_Objective(object):
    """Objective function for emission-line least-squares fitting.

    Parameters
    ----------
    obs_bin_centers : :class:`numpy.ndarray`
        Center wavelength of each observed wavelength bin.
    obs_fluxes : :class:`numpy.ndarray`
        Observed flux in each wavelength bin.
    obs_weights : :class:`numpy.ndarray`
        Weighting factors for each wavelength bin in the residual.
    redshift : float
        Redshift of the observed spectrum.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
        Resolution matrices for each camera.
    camerapix : :class:`numpy.ndarray` of int
        Start and end wavelength bin indices for each camera.
    params_mapping : :class:`~fastspecfit.emline_fit.params_mapping.ParamsMapping`
        Mapping from free to full line parameters.
    continuum_patches : dict or None, optional
        Patch pedestal parameters, or ``None`` if not used.

    """

    def __init__(self,
                 obs_bin_centers,
                 obs_fluxes,
                 obs_weights,
                 redshift,
                 line_wavelengths,
                 resolution_matrices,
                 camerapix,
                 params_mapping,
                 continuum_patches=None):

        self.dtype = obs_fluxes.dtype

        self.obs_fluxes = obs_fluxes
        self.obs_weights = obs_weights
        self.redshift = redshift
        self.line_wavelengths = line_wavelengths
        self.resolution_matrices = resolution_matrices
        self.camerapix = camerapix
        self.params_mapping = params_mapping
        self.obs_bin_centers = obs_bin_centers

        if continuum_patches is not None:
            self.nPatches = len(continuum_patches)

            self.patch_endpts = continuum_patches['endpts'].value
            self.patch_pivotwave = continuum_patches['pivotwave'].value

            self.J_P = patch_jacobian(obs_bin_centers,
                                      obs_weights,
                                      self.patch_endpts,
                                      self.patch_pivotwave)
        else:
            self.nPatches = 0
            self.J_P = None

        self.log_obs_bin_edges, self.ibin_widths = \
            _prepare_bins(obs_bin_centers, camerapix)

        # temporary storage to prevent allocation in params_mapping
        # on every call to objective/jacobian
        self.line_parameters = \
            np.empty(self.params_mapping.nParms, dtype=self.dtype)


    def objective(self, free_parameters):
        """Compute the weighted residuals between the emission-line model and observations.

        Parameters
        ----------
        free_parameters : :class:`numpy.ndarray`
            Values of all free parameters.

        Returns
        -------
        residuals : :class:`numpy.ndarray`
            Weighted residuals for each wavelength bin.

        """

        nLineParameters = len(free_parameters) - 2 * self.nPatches
        line_free_parameters = free_parameters[:nLineParameters]

        #
        # expand line free parameters into complete
        # line parameter array, handling tied params
        # and doublets
        #
        line_parameters = self.params_mapping.mapFreeToFull(line_free_parameters,
                                                            out=self.line_parameters)

        model_fluxes = np.empty_like(self.obs_fluxes, dtype=self.dtype)

        _build_model_core(line_parameters,
                          self.line_wavelengths,
                          self.redshift,
                          self.log_obs_bin_edges,
                          self.ibin_widths,
                          self.resolution_matrices,
                          self.camerapix,
                          model_fluxes)

        if self.nPatches > 0:
            slopes     = free_parameters[nLineParameters:nLineParameters + self.nPatches]
            intercepts = free_parameters[nLineParameters+self.nPatches:]

            # add patch pedestals to line model
            _add_patches(self.obs_bin_centers, model_fluxes,
                         self.patch_endpts,
                         self.patch_pivotwave,
                         slopes,
                         intercepts)

        # turn model fluxes into residuals in-place to avoid
        # unwanted memory allocation
        residuals  = model_fluxes
        residuals -= self.obs_fluxes
        residuals *= self.obs_weights

        return residuals


    def jacobian(self, free_parameters):
        """Compute the sparse Jacobian of the objective function.

        Parameters
        ----------
        free_parameters : :class:`numpy.ndarray`
            Values of all free parameters.

        Returns
        -------
        J : :class:`~fastspecfit.emline_fit.sparse_rep.EMLineJacobian`
            Sparse Jacobian at the given parameter values.

        """

        nLineFreeParms = len(free_parameters) - 2 * self.nPatches
        line_free_parameters = free_parameters[:nLineFreeParms]

        #
        # expand free paramters into complete
        # parameter array, handling tied params
        # and doublets
        #

        line_parameters = self.params_mapping.mapFreeToFull(line_free_parameters,
                                                            out=self.line_parameters)

        J_S = self.params_mapping.getJacobian(line_free_parameters)

        jacs = []
        for icam, campix in enumerate(self.camerapix):
            s, e = campix

            # Actual inverse bin widths are in ibin_widths[s+1:e+2].
            # Setup guarantees that at least one more array entry
            # exists to either side of this range, so we can pass
            # those in as dummies.
            ibw = self.ibin_widths[s:e+3]

            idealJac = \
                emline_model_jacobian(line_parameters,
                                      self.log_obs_bin_edges[s+icam:e+icam+1],
                                      ibw,
                                      self.redshift,
                                      self.line_wavelengths,
                                      self.resolution_matrices[icam].ndiag)

            # ignore any columns corresponding to fixed parameters
            endpts = idealJac[0]
            endpts[self.params_mapping.fixedMask(), :] = 0

            jacs.append( mulWMJ(self.obs_weights[s:e],
                                self.resolution_matrices[icam].rowdata(),
                                idealJac) )

        nBins = np.sum(np.diff(self.camerapix))
        nFreeParms = len(free_parameters)
        J =  EMLineJacobian((nBins, nFreeParms), nLineFreeParms,
                            self.camerapix, jacs, J_S,
                            self.J_P)

        return J


def build_model(redshift,
                line_parameters,
                line_wavelengths,
                obs_bin_centers,
                resolution_matrices,
                camerapix,
                continuum_patches=None):
    """Compute the resolution-convolved emission-line model fluxes.

    Parameters
    ----------
    redshift : float
        Redshift of the observed spectrum.
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    obs_bin_centers : :class:`numpy.ndarray`
        Center wavelength of each observed wavelength bin.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
        Resolution matrices for each camera.
    camerapix : :class:`numpy.ndarray` of int
        Start and end wavelength bin indices for each camera.
    continuum_patches : dict or None, optional
        Patch pedestal parameters, or ``None`` if not used.

    Returns
    -------
    model_fluxes : :class:`numpy.ndarray`
        Modeled flux in each observed wavelength bin.

    """

    log_obs_bin_edges, ibin_widths = _prepare_bins(obs_bin_centers, camerapix)

    model_fluxes = np.empty_like(obs_bin_centers, dtype=obs_bin_centers.dtype)

    _build_model_core(line_parameters,
                      line_wavelengths,
                      redshift,
                      log_obs_bin_edges,
                      ibin_widths,
                      resolution_matrices,
                      camerapix,
                      model_fluxes)

    # suppress negative pixels arising from resolution matrix
    model_fluxes[model_fluxes < 0.] = 0.

    if continuum_patches is not None:
        # add patch pedestals to model fluxes
        _add_patches(obs_bin_centers, model_fluxes,
                     continuum_patches['endpts'].value,
                     continuum_patches['pivotwave'].value,
                     continuum_patches['slope'].value,
                     continuum_patches['intercept'].value)

    return model_fluxes



class MultiLines(object):
    """Per-line emission-line flux models across one or more cameras.

    Builds and stores individual resolution-convolved line profiles so
    that the caller can retrieve a sparse model for any single line.

    Parameters
    ----------
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    obs_bin_centers : :class:`numpy.ndarray`
        Center wavelength of each observed wavelength bin.
    redshift : float
        Redshift of the observed spectrum.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
        Resolution matrices for each camera.
    camerapix : :class:`numpy.ndarray` of int
        Start and end wavelength bin indices for each camera.

    """

    def __init__(self,
                 line_parameters,
                 obs_bin_centers,
                 redshift,
                 line_wavelengths,
                 resolution_matrices,
                 camerapix):

        @jit(nopython=True, nogil=True, cache=True)
        def _suppress_negative_fluxes(endpts, M):
            """Clamp negative per-line fluxes to zero."""
            for i in range(M.shape[0]):
                s, e = endpts[i]
                for j in range(e-s):
                    M[i,j] = np.maximum(M[i,j], 0.)

        if resolution_matrices is None:
            # create trivial diagonal resolution matrices
            from fastspecfit.resolution import Resolution
            rm = [ Resolution(np.ones((1, e - s))) for (s, e) in camerapix ]
            resolution_matrices = tuple(rm)

        self.line_models = []
        _build_multimodel_core(line_parameters,
                               obs_bin_centers,
                               redshift,
                               line_wavelengths,
                               resolution_matrices,
                               camerapix,
                               lambda m: self.line_models.append(m))

        self.line_models = tuple(self.line_models)

        for endpts, M in self.line_models:
            _suppress_negative_fluxes(endpts, M)


    def getLine(self, line):
        """Return the model for one emission line as a sparse array.

        Parameters
        ----------
        line : int
            Index of the line to return.

        Returns
        -------
        range : tuple of int
            ``(s, e)`` giving the half-open bin range ``[s, e)`` in the
            combined observed wavelength array with nonzero flux.
        data : :class:`numpy.ndarray`
            Flux values for bins ``obs_bin_centers[s:e]``.

        """

        s = 1000000000
        e =-1

        # Determine which cameras' multiline models contain
        # bins with nonzero flux for this line, and compute
        # the lowest and highest such bins in the combined
        # observed flux array.
        live_models = []
        for i, line_model in enumerate(self.line_models):

            ls, le = line_model[0][line]

            if ls < le: # line has nonzero flux bins
                s = np.minimum(s, ls)
                e = np.maximum(e, le)
                live_models.append(i)

        if len(live_models) == 0:
            # line has no nonzero flux bins
            return (0, 0), np.empty((0), dtype=np.float64)
        elif len(live_models) == 1:
            # line appears on only one camera
            i = live_models[0]
            return (s, e), self.line_models[i][1][line][:e-s]
        else:
            # build combind array across multiple cameras,
            # allowing for arbitrary overlap.  This array
            # represents the range obs_bins[s:e]
            data = np.zeros(e-s, self.line_models[0][1].dtype)

            # copy the live bins for each camera to the
            # combind array.
            for i in live_models:
                ls, le = self.line_models[i][0][line]
                ldata  = self.line_models[i][1][line]

                data[ls-s:le-s] = ldata[:le-ls]

            return (s, e), data


def find_peak_amplitudes(line_parameters,
                         obs_bin_centers,
                         redshift,
                         line_wavelengths,
                         resolution_matrices,
                         camerapix):
    """Return the maximum flux contributed by each line to any observed bin.

    Parameters
    ----------
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    obs_bin_centers : :class:`numpy.ndarray`
        Center wavelength of each observed wavelength bin.
    redshift : float
        Redshift of the observed spectrum.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
        Resolution matrices for each camera.
    camerapix : :class:`numpy.ndarray` of int
        Start and end wavelength bin indices for each camera.

    Returns
    -------
    max_amps : :class:`numpy.ndarray`
        Maximum flux in any observed bin for each line.

    """

    @jit(nopython=True, nogil=True, cache=True)
    def _update_line_maxima(max_amps, line_models):
        """Update per-line maximum amplitudes from a camera's line models."""
        endpts, vals = line_models

        # find the highest flux for each peak; if it's
        # bigger than any seen so far, update global max
        for i in range(vals.shape[0]):
            ps, pe = endpts[i]
            if pe > ps:
                max_amps[i] = np.maximum(max_amps[i],
                                         np.max(vals[i,:pe-ps]))

    max_amps = np.zeros_like(line_wavelengths, dtype=line_parameters.dtype)

    _build_multimodel_core(line_parameters,
                           obs_bin_centers,
                           redshift,
                           line_wavelengths,
                           resolution_matrices,
                           camerapix,
                           lambda m: _update_line_maxima(max_amps, m))

    return max_amps


def find_peak_amplitudes_and_fluxes(line_parameters,
                                    obs_bin_centers,
                                    redshift,
                                    line_wavelengths,
                                    resolution_matrices,
                                    camerapix):
    """Return the peak amplitude and integrated flux for each emission line.

    Parameters
    ----------
    line_parameters : :class:`numpy.ndarray`
        Concatenated array of amplitudes, velocity shifts, and sigmas for
        all lines.
    obs_bin_centers : :class:`numpy.ndarray`
        Center wavelength of each observed wavelength bin.
    redshift : float
        Redshift of the observed spectrum.
    line_wavelengths : :class:`numpy.ndarray`
        Nominal rest-frame wavelengths of all fitted lines in Angstroms.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
        Resolution matrices for each camera.
    camerapix : :class:`numpy.ndarray` of int
        Start and end wavelength bin indices for each camera.

    Returns
    -------
    max_amps : :class:`numpy.ndarray`
        Maximum flux in any observed bin for each line.
    line_fluxes : :class:`numpy.ndarray`
        Wavelength-integrated flux of the convolved line profile for each line.

    """
    @jit(nopython=True, nogil=True, cache=True)
    def _update_line_maxima_and_fluxes(max_amps, line_fluxes, line_models, wave_weight):
        endpts, vals = line_models

        for i in range(vals.shape[0]):
            ps, pe = endpts[i]
            if pe > ps:
                max_amps[i] = np.maximum(max_amps[i], np.max(vals[i, :pe-ps]))
                for k in range(pe - ps):
                    line_fluxes[i] += vals[i, k] * wave_weight[ps + k]

    nlines = len(line_wavelengths)
    nbins  = len(obs_bin_centers)

    # Build per-pixel integration weights accounting for wavelength overlap
    wave_weight = np.empty(nbins)
    for s, e in camerapix:
        wave_weight[s:e] = np.gradient(obs_bin_centers[s:e])

    # For each camera pair, find wavelength overlap and halve weights there.
    for i, (s1, e1) in enumerate(camerapix):
        for j, (s2, e2) in enumerate(camerapix):
            if j <= i:
                continue
            wmin = max(obs_bin_centers[s1:e1].min(), obs_bin_centers[s2:e2].min())
            wmax = min(obs_bin_centers[s1:e1].max(), obs_bin_centers[s2:e2].max())
            if wmax > wmin:
                mask1 = np.zeros(nbins, dtype=bool)
                mask2 = np.zeros(nbins, dtype=bool)
                mask1[s1:e1] = (obs_bin_centers[s1:e1] >= wmin) & \
                               (obs_bin_centers[s1:e1] <= wmax)
                mask2[s2:e2] = (obs_bin_centers[s2:e2] >= wmin) & \
                               (obs_bin_centers[s2:e2] <= wmax)
                wave_weight[mask1] /= 2.
                wave_weight[mask2] /= 2.

    max_amps    = np.zeros(nlines, dtype=line_parameters.dtype)
    line_fluxes = np.zeros(nlines, dtype=line_parameters.dtype)

    _build_multimodel_core(line_parameters,
                           obs_bin_centers,
                           redshift,
                           line_wavelengths,
                           resolution_matrices,
                           camerapix,
                           lambda m: _update_line_maxima_and_fluxes(max_amps, line_fluxes, m, wave_weight))

    return max_amps, line_fluxes


def _build_model_core(line_parameters,
                      line_wavelengths,
                      redshift,
                      log_obs_bin_edges,
                      ibin_widths,
                      resolution_matrices,
                      camerapix,
                      model_fluxes):
    """Compute resolution-convolved combined model fluxes, writing into ``model_fluxes``."""

    for icam, campix in enumerate(camerapix):

        # start and end for obs fluxes of camera icam
        s, e = campix

        # Actual inverse bin widths are in ibin_widths[s+1:e+2].
        # Setup guarantees that at least one more array entry
        # exists to either side of this range, so we can pass
        # those in as dummies.
        ibw = ibin_widths[s:e+3]

        mf = emline_model(line_wavelengths,
                          line_parameters,
                          log_obs_bin_edges[s+icam:e+icam+1],
                          redshift,
                          ibw)

        # convolve model with resolution matrix and store in
        # this camera's subrange of model_fluxes
        resolution_matrices[icam].dot(mf, model_fluxes[s:e])


def _build_multimodel_core(line_parameters,
                           obs_bin_centers,
                           redshift,
                           line_wavelengths,
                           resolution_matrices,
                           camerapix,
                           consumer_fun):
    """Compute sparse per-line models for each camera and pass each to ``consumer_fun``."""

    log_obs_bin_edges, ibin_widths = _prepare_bins(obs_bin_centers,
                                                   camerapix)

    for icam, campix in enumerate(camerapix):

        # start and end for obs fluxes of camera icam
        s, e = campix

        # Actual inverse bin widths are in ibin_widths[s+1:e+2].
        # Setup guarantees that at least one more array entry
        # exists to either side of this range, so we can pass
        # those in as dummies.
        ibw = ibin_widths[s:e+3]

        # compute model waveform for each spectral line
        line_models = emline_perline_models(line_wavelengths,
                                            line_parameters,
                                            log_obs_bin_edges[s+icam:e+icam+1],
                                            redshift,
                                            ibw,
                                            resolution_matrices[icam].ndiag)

        # convolve each line's waveform with resolution matrix
        endpts, M = mulWMJ(np.ones(e - s),
                           resolution_matrices[icam].rowdata(),
                           line_models)

        # adjust endpoints to reflect camera range
        endpts += s

        consumer_fun((endpts, M))


@jit(nopython=True, nogil=True, cache=True)
def _add_patches(obs_bin_centers,
                 model_fluxes,
                 patch_endpts,
                 patch_pivotwaves,
                 slopes,
                 intercepts):
    """Add affine patch pedestals to ``model_fluxes`` in-place."""

    # add patch pedestals to line model
    nPatches = len(slopes)

    for ipatch in range(nPatches):
        s, e      = patch_endpts[ipatch]
        pivotwave = patch_pivotwaves[ipatch]

        slope     = slopes[ipatch]
        intercept = intercepts[ipatch]

        for j in range(s,e):
            model_fluxes[j] += \
                slope * (obs_bin_centers[j] - pivotwave) + intercept


@jit(nopython=True, nogil=True, cache=True)
def mulWMJ(w, M, Jsp):
    """Compute the sparse matrix product ``P = W @ M @ J``.

    ``W`` is a diagonal weight matrix, ``M`` is a resolution matrix, and
    ``J`` is a column-sparse matrix with one contiguous nonzero range per
    column.  The result is written back into ``Jsp`` in-place.

    Parameters
    ----------
    w : :class:`numpy.ndarray`
        Diagonal of the weight matrix ``W``.
    M : :class:`numpy.ndarray`, shape (nbins, ndiag)
        Resolution matrix in sparse row form.
    Jsp : tuple
        Column-sparse matrix ``(endpts, J)``, where ``endpts[j] = (s, e)``
        gives the half-open nonzero range for column ``j`` and
        ``J[j, :e-s]`` holds the values.

    Returns
    -------
    result : tuple
        Product ``W @ M @ J`` in the same column-sparse form as ``Jsp``.

    Notes
    -----
    The input ``Jsp`` must have sufficient padding so that each column of
    ``J`` can be overwritten with the corresponding column of ``P``.

    """
    endpts, J = Jsp

    nbins, ndiag = M.shape
    ncol, maxColSize = J.shape

    hdiag = ndiag//2

    # temporary buffer for each column of WMJ
    buf = np.empty(maxColSize, dtype=J.dtype)

    for j in range(ncol):
        # boundaries of nonzero entries
        # in jth column of J
        s, e = endpts[j]

        if s == e: # no nonzero values in column j
            continue

        # boundaries of entries in jth column of P
        # impacted by matrix multiply
        imin = np.maximum(s - hdiag, 0)
        imax = np.minimum(e + hdiag, nbins) # one past last impacted entry

        for i in range(imin, imax):

            # boundaries of interval of k where both
            # M[i, k] and J[k, j] are nonzero.
            kmin = np.maximum(i - hdiag,     s)
            kmax = np.minimum(i + hdiag, e - 1)

            acc = 0.
            for k in range(kmin, kmax + 1):
                acc += M[i, k - i + hdiag] * J[j, k - s]
            buf[i - imin] = acc * w[i]

        # write result back to J and update endpts for col
        newS = np.maximum(imin, 0)
        newE = np.minimum(imax, nbins)
        J[j, :newE - newS] = buf[:newE - newS]
        endpts[j] = np.array([newS, newE])

    return (endpts, J)


@jit(nopython=True, nogil=True, cache=True)
def _prepare_bins(centers, camerapix):
    """Convert observed bin centers to log bin edges and inverse bin widths.

    Parameters
    ----------
    centers : :class:`numpy.ndarray`
        Center wavelength of each observed wavelength bin.
    camerapix : :class:`numpy.ndarray` of int
        Start and end wavelength bin indices for each camera.

    Returns
    -------
    log_obs_bin_edges : :class:`numpy.ndarray`
        Natural log of each bin edge wavelength.  Edges are placed
        halfway between adjacent centers, with extrapolation at the ends.
        One extra edge per camera gap is included.
    ibin_widths : :class:`numpy.ndarray`
        Inverse bin widths, zero-padded by one entry on each end.

    """

    ncameras = camerapix.shape[0]
    edges = np.empty(len(centers) + ncameras, dtype=centers.dtype)
    ibin_widths = np.empty(len(centers) + 2,  dtype=centers.dtype)

    for icam, campix in enumerate(camerapix):

        s, e = campix
        icenters = centers[s:e]

        #- interior edges are just points half way between bin centers
        int_edges = 0.5 * (icenters[:-1] + icenters[1:])

        #- exterior edges are extrapolation of interior bin sizes
        edge_l = icenters[ 0] - (icenters[ 1] - int_edges[ 0])
        edge_r = icenters[-1] + (icenters[-1] - int_edges[-1])

        edges[s + icam]              = edge_l
        edges[s + icam + 1:e + icam] = int_edges
        edges[e + icam]              = edge_r

        # add 1 to indices i ibin_widths to skip dummy at 0
        ibin_widths[s+1:e+1] = 1. / np.diff(edges[s+icam : e+icam+1])

    # dummies before and after widths are needed
    # for corner cases in edge -> bin computation
    ibin_widths[0]  = 0.
    ibin_widths[-1] = 0.

    return (np.log(edges), ibin_widths)
