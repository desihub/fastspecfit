import numpy as np
from math import erf, erfc

from numba import jit

from fastspecfit.resolution import Resolution

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
    """
    Objective function for emission-line fitting.
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
        """
        Parameters
        ----------
        obs_bin_centers : :class:`np.ndarray` [# obs wavelength bins]
          Center wavelength of each observed wavelength bin.
        obs_fluxes : :class:`np.ndarray` [# obs wavelength bins]
          Observed flux in each wavelength bin.
        obs_weights : :class:`np.ndarray` [# obs wavelength bins]
          Weighting factors for each wavelength bin in residual.
        redshift : :class:`np.float64`
          Red shift of observed spectrum.
        line_wavelengths : :class:`np.ndarray` [# lines]
          Array of nominal wavelengths for all fitted lines.
        resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
          Resolution matrices for each camera.
        params_mapping : :class:`.params_mapping.ParamsMapping`
          Mapping from free to full parameters.
        continuum_patches : patch dictionary or :None:
          If using patch pedestals, a dictionary of patch parameters;
          else None.

        """

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
        """
        Objective function for least-squares optimization.

        Build the emline model as described above and compute the weighted
        vector of residuals between the modeled fluxes and the observations.

        Parameters
        ----------
        free_parameters : :class:`np.ndarray`
          Values of all free parameters.

        Returns
        -------
        `np.ndarray` of residuals for each wavelength bin.

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
        """
        Jacobian of objective function for emission line fitting
        optimization. The result of the detailed calculation is converted
        to a sparse matrix, since it is extremely sparse, to speed up
        subsequent matrix-vector multiplies in the optimizer.

        Parameters
        ----------
        free_parameters : :class:`np.ndarray`
          Values of all free parameters.

        Returns
        -------
        Sparse Jacobian at given parameter value as
        defined in sparse_rep.py.

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


##################################################################################


def build_model(redshift,
                line_parameters,
                line_wavelengths,
                obs_bin_centers,
                resolution_matrices,
                camerapix,
                continuum_patches=None):
    """
    Compatibility entry point to compute modeled fluxes.

    Parameters
    ----------
    redshift : :class:`np.float64`
      Red shift of observed spectrum.
    line_parameters : :class:`np.ndarray`
      Array of all fitted line parameters.
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    obs_bin_centers : :class:`np.ndarray` [# obs wavelength bins]
      Center wavelength of each observed wavelength bin.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
      Resolution matrices for each camera.
    camerapix : :class:`np.ndarray` of `int` [# cameras x 2]
      Pixels corresponding to each camera in obs wavelength array.
    continuum_patches : patch dictionary or :None:
      If using patch pedestals, a dictionary of patch parameters;
      else None.

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
    """
    Construct models for each of a list of lines across one or
    more cameras.  Return an object that lets caller obtain
    a rendered model for each individual line as a sparse array.
    """

    def __init__(self,
                 line_parameters,
                 obs_bin_centers,
                 redshift,
                 line_wavelengths,
                 resolution_matrices,
                 camerapix):
        """
        Given fitted parameters for all emission lines, compute
        models for each line for each camera.

        Parameters
        ----------
        line_parameters : :class:`np.ndarray`
          Array of all fitted line parameters.
        obs_bin_centers : :class:`np.ndarray` [# obs wavelength bins]
          Center wavelength of each observed wavelength bin.
        redshift : :class:`np.float64`
          Red shift of observed spectrum.
        line_wavelengths : :class:`np.ndarray` [# lines]
          Array of nominal wavelengths for all fitted lines.
        resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
          Resolution matrices for each camera.
        camerapix : :class:`np.ndarray` of `int` [# cameras x 2]
          Pixels corresponding to each camera in obs wavelength array.

        """

        @jit(nopython=True, nogil=True, cache=True)
        def _suppress_negative_fluxes(endpts, M):
            """
            suppress negative fluxes arising from resolution matrix
            """
            for i in range(M.shape[0]):
                s, e = endpts[i]
                for j in range(e-s):
                    M[i,j] = np.maximum(M[i,j], 0.)

        if resolution_matrices is None:
            # create trivial diagonal resolution matrices
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
        """
        Return a model for one emission line.

        Parameters
        ----------
        line : :class:`int`
          Index of line to return.

        Returns
        -------
        Sparse line model as tuple ((s, e), data), where the range
        [s,e) in the combined bin array across all cameras (as in
        obs_bin_centers) contains all bins with nonzero fluxes for
        that line, and data is an array of length e - s that contains
        the bin values.

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
    """
    Given fitted parameters for all emission lines, report for
    each line the largest flux that it contributes to any observed
    bin.

    Parameters
    ----------
    parameters : :class:`np.ndarray`
      Array of all fitted line parameters.
    bin_data : :class:`tuple`
      Preprocessed bin data in form returned by _prepare_bins
    redshift : :class:`np.float64`
      Red shift of observed spectrum.
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
      Resolution matrices for each camera.
    camerapix : :class:`np.ndarray` of `int` [# cameras x 2]
      Pixels corresponding to each camera in obs wavelength array.

    Returns
    -------
    Array with maximum amplitude observed for each line.

    """

    @jit(nopython=True, nogil=True, cache=True)
    def _update_line_maxima(max_amps, line_models):
        """
        Given an array of line waveforms and an array of maximum
        amplitudes observed for each line so far, if the waveform contains
        a higher amplitude than the previous maximum, update the maximum
        in place.
        """
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


##########################################################################


def _build_model_core(line_parameters,
                      line_wavelengths,
                      redshift,
                      log_obs_bin_edges,
                      ibin_widths,
                      resolution_matrices,
                      camerapix,
                      model_fluxes):
    """
    Core loop for computing a combined model flux from a set of
    spectral emission lines.

    Parameters
    ----------
    line_parameters : :class:`np.ndarray`
      Combined array of amplitudes, vshifts, and sigmas
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    redshift : :class:`np.float64`
      Red shift of observed spectrum.
    log_obs_bin_edges : :class:`np.ndarray` [# obs wavelength bins + 1]
      Natural logs of observed wavelength bin edges
    ibin_widths : :class:`np.ndarray` [# obs wavelength bins]
      Inverse widths of each observed wavelength bin.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
      Resolution matrices for each camera.
    camerapix : :class:`np.ndarray` of `int` [# cameras x 2]
      Pixels corresponding to each camera in obs wavelength array.
    model_fluxes : :class:`np.ndarray` [# obs wavelength bins] (output)
      Returns computed model flux for each wavelength bin.

    """

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
    """
    Core loop for computing array of individual line flux models sparsely
    from a set of spectral emission lines.  Result is not returned but
    passed to a consumer function supplied by the caller.

    Parameters
    ----------
    line_parameters : :class:`np.ndarray`
      Combined array of amplitudes, vshifts, and sigmas
    obs_bin_centers : :class:`np.ndarray` [# obs wavelength bins]
      Center wavelength of each observed wavelength bin.
    redshift : :class:`np.float64`
      Redshift of observed spectrum.
    line_wavelengths : :class:`np.ndarray` [# lines]
      Array of nominal wavelengths for all fitted lines.
    resolution_matrices : tuple of :class:`fastspecfit.resolution.Resolution`
      Resolution matrices for each camera.
    camerapix : :class:`np.ndarray` of `int` [# cameras x 2]
      Pixels corresponding to each camera in obs wavelength array.
    consumer_fun :
      Function to which we pass computed sparse array of line models
      for each camera.

    """

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
    """
    Add patch pedestals for patches with given endpts, slopes,
    intercepts, and pivots to an array of model fluxes.

    Parameters
    ----------
    obs_bin_centers : :class:`np.ndarray` [# obs wavelength bins]
      Center wavelength of each obesrved bin.
    model_fluxes : :class:`np.ndarray` [# obs wavelength bins] (output)
      Fluxes computed from base emission line model.  Modified
      by addition of patch pedestal contribution to each bin.
    patch_endpts: :class:`np.ndarray` of `int` [2 * # patches]
      Endpoint indices of each patch in wavelength bin array
    patch_pivotwave: :class:`np.ndarray` [# patches]
      Offset for each patch pedestal.
    slopes : :class:`np.ndarray` [# patches]
      Slope of each patch pedestal.
    itercepts: :class:`np.ndarray` [# patches]
      Intercept of each patch pedestal.

    """

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


###################################################################


@jit(nopython=True, nogil=True, cache=True)
def mulWMJ(w, M, Jsp):
    """Compute the sparse matrix product P = WMJ, where
    W is a diagonal weight matrix
    M is a resolution matrix
    J is a column-sparse matrix giving the nonzero
      values in one contiguous range per column

    Parameters
    ----------
    w : :class:`np.ndarray` [# obs wavelength bins]
      Diagonal of matrix W.
    M : :class:`np.ndarray` [# obs wavelength bins x # diags]
      Resolution matrix, represented in sparse ROW form
      giving nonzero entries in each row.
    Jsp : :class:`tuple`
      column-sparse matrix (endpts, J), where endpts = (s,e)
      gives a half-open range of indices [s, e) with values
      for this column, and P[:e-s] contains these values.

    Returns
    -------
    Product WMJ in column-sparse form as tuple (endpts, P)
    where endpts = (s,e) gives a half-open range of indices
    [s, e) with values for this column, and P[:e-s] contains these
    values.

    Note
    ----
    We assume that the column-sparse input Jsp has enough space
    allocated that we can overwrite each column of Jsp with the
    corresponding column of P.

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
    """
    Convert bin centers to the info needed by the optimizer,
    returned as an opaque tuple.

    Parameters
    ----------
    centers : :class:`np.ndarray` [# obs wavelength bins]
      Center wavelength of each obs bin
    camerapix : :class:`np.ndarray` of `int` [# cameras x 2]
      pixels corresponding to each camera in obs wavelength array

    Returns
    -------
    Tuple containing
       - array of log bin edges.  Edges are placed halfway between centers,
         with extrapolation at the ends.  We return the natural log of
         each edge's wavelength, since that is what model and Jacobian
          building need.
      - array of inverse widths for each wavelength bin. The array is
        zero-padded by one cell on the left and right to accomodate
        edge-to-bin computations.

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
