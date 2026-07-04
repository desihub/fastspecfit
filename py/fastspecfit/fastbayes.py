"""
fastspecfit.fastbayes
=====================

Grid-based Bayesian broadband-photometric SED fitting.

Given a pre-synthesized photometry grid (see ``bin/build-bayesian-templates``
and ``bin/build-bayesian-photometry``), fit each object's observed broadband
photometry by interpolating the grid to the object's (Redrock) redshift,
solving for the chi2-minimizing stellar-mass amplitude of every grid template
in closed form, and building up the posterior probability distribution of
every grid parameter by weighting each template by its likelihood. Point
estimates (mean, mode, and 16/50/84th percentiles) are reported for every
parameter from this discrete-grid posterior.

Because the underlying grid has finite resolution, the maximum-likelihood
point is additionally refined: since the grid is a full factorial design
(every combination of its 7 axes exists), the templates neighboring the
discrete chi2 minimum are known exactly via index arithmetic, with no need
for scattered-data interpolation. A local parabola fit to chi2 along each
axis independently (:func:`fastspecfit.util.minfit`, the same utility used
elsewhere in this package to refine a velocity-dispersion grid minimum)
gives a sub-grid-resolution parameter estimate and a formal (delta-chi2=1)
uncertainty for each of the 7 grid axes. The same local fractional offsets
feed an N-linear interpolation over the (up to 2**7) neighboring templates,
which is used to build a refined model spectrum, re-solve the mass
amplitude, and recompute chi2 for that refined model -- so the reported
CHI2 corresponds to the same model as the reported parameters, not the
unrefined discrete grid point. K-corrections, absolute magnitudes,
rest-frame luminosities, and the model Dn(4000) index are computed from
that refined rest-frame spectrum, which is also written to a MODELS
extension for later QA/analysis without needing to re-read the templates
file.

Axes built log-uniform (age, zzsun, umin, qpah) are refined and reported in
log10 space (as LOGAGE, LOGZZSUN, LOGUMIN, LOGQPAH) so that their formal
uncertainties are symmetric; axes built linear-uniform (tau, dustn, gamma)
are refined and reported in linear space.

IGM attenuation and cosmological dimming for the *grid* photometry are
already baked into the photometry grid (see ``bin/build-bayesian-
photometry``); this module only needs the IGM model and cosmology directly
to redshift the one refined rest-frame spectrum per object for the K-
correction/absolute-magnitude calculation.

"""
import os, sys, time, logging
import numpy as np
import fitsio
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.util import MPPool, fsftime, minfit, C_LIGHT
from fastspecfit.photometry import Photometry
from fastspecfit.singlecopy import sc_data

# Grid axes in the exact nested-loop order used by bin/build-bayesian-templates
# (age outermost ... qpah innermost). This order is what lets the flattened
# grid's N-D structure be recovered via np.unravel_index/np.ravel_multi_index.
GRID_AXIS_COLUMNS = ('age', 'zzsun', 'tau', 'dustn', 'umin', 'gamma', 'qpah')

# Axes built log-uniform: refined/reported in log10 space so that the formal
# (delta-chi2=1) uncertainty is symmetric. The remaining axes (tau, dustn,
# gamma) were built linear-uniform and are refined/reported in linear space.
LOG_AXES = frozenset(('age', 'zzsun', 'umin', 'qpah'))

# Output parameter name for each grid axis, and its inverse.
_AXIS_OUTNAME = {
    'age': 'LOGAGE', 'zzsun': 'LOGZZSUN', 'tau': 'TAU', 'dustn': 'DUSTN',
    'umin': 'LOGUMIN', 'gamma': 'GAMMA', 'qpah': 'LOGQPAH',
}
_OUTNAME_TO_AXIS = {v: k for k, v in _AXIS_OUTNAME.items()}

# All reported parameters: the 7 grid axes, plus LOGMASS and SFR, which are
# derived per-object from the closed-form amplitude solve.
PARAM_NAMES = tuple(_AXIS_OUTNAME[col] for col in GRID_AXIS_COLUMNS) + ('LOGMASS', 'SFR')

# Rest-frame luminosity output keys and reference wavelengths (Angstrom),
# matching the LOGL_*/LOGLNU_* columns in fastspecfit's own SPECPHOT schema
# (see fastspecfit.continuum.ContinuumTools.lums_keys).
LUM_KEYS = ('LOGL_1450', 'LOGLNU_1500', 'LOGL_1700', 'LOGLNU_2800', 'LOGL_3000', 'LOGL_5100')
LUM_WAVES = (1450., 1500., 1700., 2800., 3000., 5100.)

_PC_CM = 3.0856775814913673e18 # [cm]
_LSUN = 3.846e33 # [erg/s]


class BayesianGrid(object):
    """Pre-synthesized Bayesian grid photometry, shared read-only across MPPool workers.

    Populated by :func:`_initialize_fastbayes_workers` once per worker
    process, mirroring the ``sc_data`` singleton pattern in
    ``fastspecfit.singlecopy``.

    """
    def __init__(self):
        self.file = None
        self._templates_fits = None
        self._wave = None


    def load(self, gridfile):
        """Load a ``bayesian-photometry-*.fits`` grid file (idempotent)."""
        if self.file == gridfile:
            return

        T = fitsio.FITS(gridfile)
        prihdr = T[0].read_header()

        self.gridnumber = prihdr.get('GRIDNUM')
        self.imf = prihdr.get('IMF')
        self.fphotofile = prihdr.get('FPHOTO')
        self.templates_file = prihdr.get('TEMPFILE')
        self.logspace = bool(prihdr.get('LOGSPACE'))

        photsys_hdr = prihdr.get('PHOTSYS')
        self.photsys_keys = [''] if photsys_hdr in (None, 'NONE') else photsys_hdr.split(',')

        self.redshift = T['REDSHIFT'].read().astype('f8')
        self.zgrid_interp = np.log10(1. + self.redshift) if self.logspace else self.redshift

        self.meta = Table(T['METADATA'].read())
        self.ntemplate = len(self.meta)

        # Recover the grid's N-D axis structure. The grid is a full
        # factorial/Cartesian-product design (every combination of the 7
        # axes exists), so the per-axis grid values and the flattened
        # index<->multi-index mapping are recoverable directly from the
        # metadata table -- no separate shape bookkeeping needs to be
        # stored in the FITS file.
        self.axis_values = {col: np.unique(np.asarray(self.meta[col], dtype='f8')) for col in GRID_AXIS_COLUMNS}
        self.dims = tuple(len(self.axis_values[col]) for col in GRID_AXIS_COLUMNS)
        if int(np.prod(self.dims)) != self.ntemplate:
            errmsg = (f'Grid file {gridfile} metadata is not a full factorial grid over '
                      f'{GRID_AXIS_COLUMNS} (dims product {int(np.prod(self.dims))} != ntemplate {self.ntemplate}).')
            log.critical(errmsg)
            raise ValueError(errmsg)

        filt = T['FILTERS'].read()
        self.bands = np.array([b.strip() for b in filt['band']])
        self.nband = len(self.bands)

        self.maggies = {}
        for key in self.photsys_keys:
            extname = 'MAGGIES' + ('' if key == '' else f'_{key}')
            self.maggies[key] = T[extname].read() # [ntemplate, nz, nband]

        # File is set last so a failed/partial load is retried rather than
        # treated as already-cached.
        self.file = gridfile
        log.debug(f'Cached Bayesian grid {self.file}')


    def templates_fits(self):
        """Lazily open (and cache) a handle to the raw templates FITS file.

        Individual FLUX rows are read on demand (tens of KB each) so that
        the full multi-GB FLUX array is never loaded into memory.

        """
        if self._templates_fits is None:
            if not self.templates_file:
                errmsg = f'Grid file {self.file} has no TEMPFILE header keyword.'
                log.critical(errmsg)
                raise ValueError(errmsg)
            self._templates_fits = fitsio.FITS(self.templates_file)
        return self._templates_fits


    def template_wave(self):
        """Native rest-frame wavelength array shared by every template."""
        if self._wave is None:
            self._wave = self.templates_fits()['WAVE'].read()
        return self._wave


    def template_flux_row(self, itemplate):
        """Read one raw rest-frame template spectrum by flat grid index."""
        return self.templates_fits()['FLUX'][int(itemplate), :].ravel()


    def interpolate_at_z(self, photsys, redshift):
        """Linearly interpolate ``MAGGIES(template, z, band)`` to a single redshift.

        Parameters
        ----------
        photsys : :class:`str`
            Photometric-system key (e.g. ``'N'``, ``'S'``, or ``''``).
        redshift : :class:`float`
            Object redshift.

        Returns
        -------
        :class:`numpy.ndarray`
            Model maggies per solar mass formed, shape (ntemplate, nband).

        """
        if photsys not in self.maggies:
            errmsg = (f"photsys '{photsys}' not present in grid file {self.file} "
                      f"(available: {list(self.maggies.keys())}).")
            log.critical(errmsg)
            raise KeyError(errmsg)

        maggies = self.maggies[photsys] # [ntemplate, nz, nband]

        zvar = np.log10(1. + redshift) if self.logspace else redshift
        zvar = np.clip(zvar, self.zgrid_interp[0], self.zgrid_interp[-1])

        idx = int(np.clip(np.searchsorted(self.zgrid_interp, zvar), 1, len(self.zgrid_interp) - 1))
        z0, z1 = self.zgrid_interp[idx - 1], self.zgrid_interp[idx]
        frac = 0. if z1 == z0 else (zvar - z0) / (z1 - z0)

        return maggies[:, idx - 1, :] + frac * (maggies[:, idx, :] - maggies[:, idx - 1, :])


# global structures with single-copy data, initially empty
bg_data = BayesianGrid()
_igm = None
_cosmo = None


def _initialize_fastbayes_workers(fphotofile=None, gridfile=None):
    """MPPool initializer: populate ``sc_data.photometry``, ``bg_data``, the
    IGM model, and the cosmology in each worker.

    ``sc_data.photometry`` is populated directly (rather than via
    ``sc_data.initialize()``) so this mode never loads the stellar template
    basis or emission-line tables, neither of which it needs.

    """
    global _igm, _cosmo

    sc_data.photometry = Photometry(fphotofile=fphotofile)
    bg_data.load(gridfile)

    if list(bg_data.bands) != list(sc_data.photometry.bands):
        errmsg = ('Band mismatch between the photometry configuration and the Bayesian grid file; '
                 'they must be built from the same (or a compatible) fphoto configuration file.')
        log.critical(errmsg)
        raise ValueError(errmsg)

    if _igm is None:
        from fastspecfit.igm import Inoue14
        _igm = Inoue14()
    if _cosmo is None:
        from fastspecfit.cosmo import TabulatedDESI
        _cosmo = TabulatedDESI()


def _weighted_percentile(values, weights, percentiles):
    """Weighted percentiles of ``values`` via linear interpolation of the weighted CDF."""
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w)
    cw /= cw[-1]
    return np.interp(np.asarray(percentiles) / 100., cw, v)


def _weighted_mode(values, weights):
    """Weighted mode: the (exact) grid value with the largest total weight."""
    uniq, inv = np.unique(values, return_inverse=True)
    binweights = np.bincount(inv, weights=weights)
    return uniq[np.argmax(binweights)]


def _axis_posterior_values(col):
    """Per-template values of grid-axis column ``col``, in its fit coordinate
    (log10 for :data:`LOG_AXES`, linear otherwise)."""
    vals = np.asarray(bg_data.meta[col], dtype='f8')
    return np.log10(vals) if col in LOG_AXES else vals


def _refine_grid_axes(ibest, chi2):
    """Locally refine each grid axis's ML value via a 3-point parabola fit.

    For each of the 7 grid axes, holds the other 6 fixed at their best-fit
    (discrete argmin chi2) index and fits a parabola to the already-computed
    chi2 at the immediate neighboring grid points along that axis alone
    (:func:`fastspecfit.util.minfit`). Because the grid is a full factorial
    design, "neighboring template" is exact index arithmetic, not a
    scattered-data interpolation/triangulation problem.

    Parameters
    ----------
    ibest : :class:`int`
        Flat index of the discrete chi2 minimum.
    chi2 : :class:`numpy.ndarray`
        chi2 for every grid template, shape (ntemplate,).

    Returns
    -------
    fit_value : :class:`dict`
        Per grid-axis-column refined value, in its fit coordinate.
    fit_ivar : :class:`dict`
        Per grid-axis-column formal (delta-chi2=1) inverse variance, in the
        same coordinate; 0 where refinement was not possible (grid edge, or
        a failed/degenerate parabola fit).
    frac : :class:`dict`
        Per grid-axis-column fractional offset in [0, 1) toward the
        relevant neighbor (0 when unrefined); feeds the N-linear
        interpolation weights over neighboring templates.
    frac_dir : :class:`dict`
        Per grid-axis-column direction of the offset (+1 or -1 grid steps
        relative to the discrete best-fit index); meaningless when the
        corresponding ``frac`` is 0.

    """
    multi_index = list(np.unravel_index(ibest, bg_data.dims))

    fit_value, fit_ivar, frac, frac_dir = {}, {}, {}, {}

    for axis_pos, col in enumerate(GRID_AXIS_COLUMNS):
        n = bg_data.dims[axis_pos]
        i0 = multi_index[axis_pos]
        grid_vals = bg_data.axis_values[col]
        log_axis = col in LOG_AXES

        center_coord = np.log10(grid_vals[i0]) if log_axis else grid_vals[i0]

        if i0 <= 0 or i0 >= n - 1:
            # at a grid edge -- cannot bracket a minimum, no refinement
            fit_value[col] = center_coord
            fit_ivar[col] = 0.
            frac[col] = 0.
            frac_dir[col] = 1
            continue

        idx = list(multi_index)
        xs, ys = [], []
        for di in (-1, 0, 1):
            idx[axis_pos] = i0 + di
            flat = np.ravel_multi_index(tuple(idx), bg_data.dims)
            v = grid_vals[i0 + di]
            xs.append(np.log10(v) if log_axis else v)
            ys.append(chi2[flat])

        x0, xerr, y0, zwarn = minfit(np.array(xs), np.array(ys))

        if zwarn != 0 or not (xs[0] < x0 < xs[2]):
            fit_value[col] = center_coord
            fit_ivar[col] = 0.
            frac[col] = 0.
            frac_dir[col] = 1
            continue

        fit_value[col] = x0
        fit_ivar[col] = 1. / xerr**2 if xerr > 0. else 0.

        if x0 >= xs[1]:
            f = (x0 - xs[1]) / (xs[2] - xs[1])
            frac_dir[col] = +1
        else:
            f = (xs[1] - x0) / (xs[1] - xs[0])
            frac_dir[col] = -1
        # Snap negligible offsets to exactly zero so that a refined point
        # sitting (to floating-point precision) exactly on a grid point
        # collapses to that single corner, rather than picking up spurious
        # near-zero-weight neighbors (and their disk reads) from rounding
        # noise in the parabola fit.
        frac[col] = 0. if f < 1e-9 else f

    return fit_value, fit_ivar, frac, frac_dir


def _corner_weights(ibest, frac, frac_dir):
    """Build N-linear interpolation weights over the grid corners bracketing
    the locally refined point.

    Parameters
    ----------
    ibest : :class:`int`
        Flat index of the discrete chi2 minimum.
    frac : :class:`dict`
        Per-axis fractional offset, from :func:`_refine_grid_axes`.
    frac_dir : :class:`dict`
        Per-axis offset direction, from :func:`_refine_grid_axes`.

    Returns
    -------
    indices : :class:`numpy.ndarray`
        Flat grid indices of the contributing corner templates (up to
        2**7 = 128; fewer whenever some axes were not refined).
    weights : :class:`numpy.ndarray`
        Corresponding non-negative weights, summing to 1.

    """
    multi_index = list(np.unravel_index(ibest, bg_data.dims))
    naxes = len(GRID_AXIS_COLUMNS)

    indices, weights = [], []
    for corner in range(2**naxes):
        idx = list(multi_index)
        w = 1.
        for axis_pos, col in enumerate(GRID_AXIS_COLUMNS):
            t = frac[col]
            if (corner >> axis_pos) & 1:
                w *= t
                idx[axis_pos] = multi_index[axis_pos] + frac_dir[col]
            else:
                w *= (1. - t)
        if w <= 0.:
            continue
        indices.append(np.ravel_multi_index(tuple(idx), bg_data.dims))
        weights.append(w)

    indices = np.array(indices, dtype='i8')
    weights = np.array(weights, dtype='f8')
    weights /= weights.sum()

    return indices, weights


def get_fastbayes_dtype(phot, topk=0):
    """Build the output dtype for the per-object Bayesian-fitting results.

    Parameters
    ----------
    phot : :class:`fastspecfit.photometry.Photometry`
        Photometry configuration, used to size the ABSMAG/KCORR columns.
    topk : :class:`int`, optional
        Number of top-weight grid templates to store per object (sparse
        joint posterior). Disabled when ``0`` (default).

    Returns
    -------
    :class:`numpy.dtype`

    """
    cols = []
    for pname in PARAM_NAMES:
        cols += [(pname, 'f4'), (f'{pname}_MEAN', 'f4'), (f'{pname}_MODE', 'f4'),
                (f'{pname}_P16', 'f4'), (f'{pname}_P50', 'f4'), (f'{pname}_P84', 'f4')]

    # Formal delta-chi2=1 uncertainties for the 7 native grid axes only;
    # LOGMASS/SFR are derived quantities without a directly fit chi2(x)
    # curve, so no formal uncertainty is reported for them yet.
    for col in GRID_AXIS_COLUMNS:
        cols += [(f'{_AXIS_OUTNAME[col]}_IVAR', 'f4')]

    cols += [('CHI2', 'f4'), ('NDOF', 'i2')]

    # K-corrections, absolute magnitudes, rest-frame luminosities, and the
    # model Dn(4000) index, computed from the refined rest-frame spectrum
    # (point estimates only for now -- no Monte Carlo/uncertainty framework
    # analogous to fastspecfit's nmonte machinery exists here yet).
    for band, shift in zip(phot.absmag_bands, phot.band_shift):
        band = band.upper()
        shift = int(10 * shift)
        cols += [(f'ABSMAG{shift:02d}_{band}', 'f4'), (f'KCORR{shift:02d}_{band}', 'f4')]
    for key in LUM_KEYS:
        cols += [(key, 'f4')]
    cols += [('DN4000_MODEL', 'f4')]

    if topk > 0:
        cols += [('TOPK_INDEX', 'i4', (topk,)), ('TOPK_WEIGHT', 'f4', (topk,))]

    return np.dtype(cols)


def fastbayes_one(iobj, data, meta, fastbayes_dtype, topk=0, uncertainty_floor=0.01,
                  qa=False, qadir='.', coadd_type='healpix'):
    """Fit one object's broadband photometry against the Bayesian grid.

    Parameters
    ----------
    iobj : :class:`int`
        Index of the object in the input list, used for log messages.
    data : :class:`dict`
        Per-object data dictionary from
        :meth:`fastspecfit.io.DESISpectra.read` (``fastphot=True``).
    meta : :class:`astropy.table.Row`
        Metadata row for this object.
    fastbayes_dtype : :class:`numpy.dtype`
        Output dtype, from :func:`get_fastbayes_dtype`.
    topk : :class:`int`, optional
        Number of top-weight grid templates to store (0 disables).
    uncertainty_floor : :class:`float`, optional
        Minimum fractional photometric uncertainty added in quadrature.
    qa : :class:`bool`, optional
        If ``True``, generate a QA figure for this object. Generated
        inline (rather than as a separate post-processing pass) because
        the full per-template posterior weights only exist in memory
        during fitting. Default is ``False``.
    qadir : :class:`str`, optional
        Output directory for the QA figure. Default is ``'.'``.
    coadd_type : :class:`str`, optional
        Coadd type, used to build the QA target label/filename. Default
        is ``'healpix'``.

    Returns
    -------
    meta : :class:`astropy.table.Row`
        Updated metadata row with observed photometry filled in.
    result : :class:`numpy.ndarray`
        Bayesian-fitting output row.
    modelspectrum : :class:`numpy.ndarray`
        Refined rest-frame maximum-likelihood spectrum (erg/s/cm2/A at
        10 pc for the fitted stellar mass), on the shared template
        wavelength grid; all zeros when the fit was skipped.

    """
    from fastspecfit.io import one_spectrum

    phot = sc_data.photometry
    npix = len(bg_data.template_wave())

    log.info(f'Bayesian fitting object {iobj} [{phot.uniqueid_col.lower()} '
            f'{data["uniqueid"]}, z={data["redshift"]:.6f}].')

    one_spectrum(data, meta, uncertainty_floor=uncertainty_floor, fastphot=True, synthphot=True)

    # Copy parsed photometry from the 'data' dictionary to the 'meta' table
    # (mirroring fastspec_one's convention).
    flux = data['photometry']['nanomaggies']
    fluxivar = data['photometry']['nanomaggies_ivar']
    for iband, band in enumerate(phot.bands):
        meta[f'FLUX_{band.upper()}'] = flux[iband]
        meta[f'FLUX_IVAR_{band.upper()}'] = fluxivar[iband]

    result = np.zeros(1, dtype=fastbayes_dtype)[0]
    modelspectrum = np.zeros(npix, dtype='f4')

    redshift = data['redshift']
    photsys = data['photsys']

    if redshift > bg_data.redshift[-1]:
        log.warning(f'Object {iobj} [{phot.uniqueid_col.lower()} {data["uniqueid"]}] redshift '
                   f'{redshift:.6f} exceeds the grid maximum {bg_data.redshift[-1]:.6f}; skipping the fit.')
        return meta, result, modelspectrum

    flam = np.asarray(data['photometry']['flam'], dtype='f8')
    flam_ivar = np.asarray(data['photometry']['flam_ivar'], dtype='f8') * phot.bands_to_fit
    lambda_eff = np.asarray(data['photometry']['lambda_eff'], dtype='f8')

    model_maggies = bg_data.interpolate_at_z(photsys, redshift) # [ntemplate, nband]
    model_flam = Photometry.get_photflam(model_maggies, lambda_eff) # [ntemplate, nband]

    # Closed-form, non-negative, chi2-minimizing mass amplitude per template.
    numer = (model_flam * (flam_ivar * flam)[np.newaxis, :]).sum(axis=1)
    denom = (model_flam**2 * flam_ivar[np.newaxis, :]).sum(axis=1)
    amplitude = np.divide(numer, denom, out=np.zeros_like(numer), where=denom > 0.)
    amplitude = np.clip(amplitude, 0., None)

    resid = flam[np.newaxis, :] - amplitude[:, np.newaxis] * model_flam
    chi2 = (flam_ivar[np.newaxis, :] * resid**2).sum(axis=1)

    chi2min = np.min(chi2)
    weight = np.exp(-0.5 * (chi2 - chi2min))
    weight /= weight.sum()

    ibest = np.argmin(chi2)
    logmass = np.log10(np.clip(amplitude, 1e-30, None))

    # The grid stores sfr per solar mass formed (like the flux templates
    # themselves), so scale by the fitted mass amplitude to get each
    # template's actual SFR estimate for this object.
    sfr_per_template = amplitude * np.asarray(bg_data.meta['sfr'], dtype='f8') # [ntemplate], Msun/yr

    # --- local refinement of the maximum-likelihood point -----------------
    fit_value, fit_ivar, frac, frac_dir = _refine_grid_axes(ibest, chi2)
    corner_idx, corner_weight = _corner_weights(ibest, frac, frac_dir)

    refined_model_maggies = (corner_weight[:, np.newaxis] * model_maggies[corner_idx, :]).sum(axis=0)
    refined_model_flam = Photometry.get_photflam(refined_model_maggies, lambda_eff)

    refined_numer = np.sum(refined_model_flam * flam_ivar * flam)
    refined_denom = np.sum(refined_model_flam**2 * flam_ivar)
    refined_amplitude = refined_numer / refined_denom if refined_denom > 0. else 0.
    refined_amplitude = max(refined_amplitude, 0.)

    refined_resid = flam - refined_amplitude * refined_model_flam
    chi2_refined = np.sum(flam_ivar * refined_resid**2)

    refined_logmass = np.log10(max(refined_amplitude, 1e-30))
    refined_sfr_per_msun = np.sum(corner_weight * np.asarray(bg_data.meta['sfr'], dtype='f8')[corner_idx])
    refined_sfr = refined_amplitude * refined_sfr_per_msun

    # --- refined rest-frame spectrum (N-linear combination of the
    # neighboring raw template spectra), scaled by the refined mass -------
    restwave = bg_data.template_wave()
    restflux = np.zeros(npix, dtype='f8')
    for idx, w in zip(corner_idx, corner_weight):
        restflux += w * bg_data.template_flux_row(idx)
    restflux *= refined_amplitude # [erg/s/cm2/A at 10pc, actual stellar mass]

    # --- K-corrections, absolute magnitudes, rest-frame luminosities, and
    # the model Dn(4000) index, all derived from the refined rest-frame
    # spectrum -------------------------------------------------------------
    zwave = restwave * (1. + redshift)
    igm_trans = _igm.full_IGM(redshift, zwave)
    dlum = _cosmo.luminosity_distance(redshift) # [Mpc]
    zfactor = igm_trans * (10. / (1e6 * dlum))**2 / (1. + redshift)
    zflux = restflux * zfactor # [erg/s/cm2/A, observed frame]

    dmod = _cosmo.distance_modulus(redshift)
    synth_absmag, synth_maggies_rest = phot.synth_absmag(redshift, dmod, zwave, zflux)
    kcorr, absmag, ivarabsmag, _ = phot.kcorr_and_absmag(
        flux, fluxivar, redshift, dmod, photsys, zwave, zflux,
        synth_absmag, synth_maggies_rest)

    dn4000_model, _ = Photometry.get_dn4000(restwave, restflux, rest=True)

    fourpi_10pc2 = 4. * np.pi * (10. * _PC_CM)**2 # [cm2]
    lums = {}
    for cwave, key in zip(LUM_WAVES, LUM_KEYS):
        cflux = np.interp(cwave, restwave, restflux) * fourpi_10pc2 # [erg/s/A]
        if 'LOGL_' in key:
            val = cflux * cwave / _LSUN / 1e10 # [1e10 Lsun]
        else:
            val = cflux * cwave**2 / (C_LIGHT * 1e13) / 1e28 # [1e-28 erg/s/Hz]
        lums[key] = np.log10(val) if val > 0. else 0.

    # --- fill the output row ------------------------------------------------
    derived = {'LOGMASS': refined_logmass, 'SFR': refined_sfr}
    posterior = {'LOGMASS': logmass, 'SFR': sfr_per_template}

    for col in GRID_AXIS_COLUMNS:
        pname = _AXIS_OUTNAME[col]
        result[pname] = fit_value[col]
        result[f'{pname}_IVAR'] = fit_ivar[col]

    posterior_arrays = {}
    for pname in PARAM_NAMES:
        if pname in derived:
            result[pname] = derived[pname]
            vals = posterior[pname]
        else:
            vals = _axis_posterior_values(_OUTNAME_TO_AXIS[pname])
        posterior_arrays[pname] = (vals, weight)
        result[f'{pname}_MEAN'] = np.sum(weight * vals)
        result[f'{pname}_MODE'] = _weighted_mode(vals, weight)
        p16, p50, p84 = _weighted_percentile(vals, weight, (16., 50., 84.))
        result[f'{pname}_P16'] = p16
        result[f'{pname}_P50'] = p50
        result[f'{pname}_P84'] = p84

    # dof = number of fitted bands minus the one continuous free parameter
    # (the per-template mass amplitude) solved for above.
    nbands_used = int(np.sum(flam_ivar > 0.))
    result['CHI2'] = chi2_refined
    result['NDOF'] = max(nbands_used - 1, 0)

    for iband, (band, shift) in enumerate(zip(phot.absmag_bands, phot.band_shift)):
        band = band.upper()
        shift = int(10 * shift)
        result[f'ABSMAG{shift:02d}_{band}'] = absmag[iband]
        result[f'KCORR{shift:02d}_{band}'] = kcorr[iband]

    for key in LUM_KEYS:
        result[key] = lums[key]
    result['DN4000_MODEL'] = dn4000_model

    if topk > 0:
        idx = np.argsort(weight)[::-1][:topk]
        result['TOPK_INDEX'][:len(idx)] = idx.astype('i4')
        result['TOPK_WEIGHT'][:len(idx)] = weight[idx].astype('f4')

    if qa:
        _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                          zwave, zflux, redshift, coadd_type=coadd_type, outdir=qadir)

    return meta, result, restflux.astype('f4')


def _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                      zwave, zflux, redshift, coadd_type='healpix', outdir='.'):
    """Generate a QA figure for one Bayesian-fit object.

    Reuses :func:`fastspecfit.qa._target_label` and
    :func:`fastspecfit.qa._fetch_cutout`, which are generic utilities with
    no fastspec/fastphot-specific coupling. Called from
    :func:`fastbayes_one` rather than as a separate post-processing pass,
    since ``posterior_arrays`` (the full per-template posterior weights)
    only exists in memory during fitting.

    Parameters
    ----------
    data : :class:`dict`
        Per-object data dictionary.
    meta : :class:`astropy.table.Row`
        Metadata row for this object.
    result : :class:`numpy.ndarray`
        Bayesian-fitting output row, from :func:`fastbayes_one`.
    posterior_arrays : :class:`dict`
        Mapping of parameter name to ``(values, weights)`` per-template
        posterior arrays.
    restwave, restflux : :class:`numpy.ndarray`
        Refined rest-frame maximum-likelihood spectrum.
    zwave, zflux : :class:`numpy.ndarray`
        The same spectrum redshifted (and IGM/distance-attenuated) to the
        observed frame.
    redshift : :class:`float`
    coadd_type : :class:`str`, optional
    outdir : :class:`str`, optional

    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from fastspecfit.io import get_qa_filename
    from fastspecfit.qa import _target_label, _fetch_cutout

    phot = sc_data.photometry

    pngfile = get_qa_filename(meta, coadd_type, outprefix='fastbayes', outdir=outdir)
    figdir = os.path.dirname(pngfile)
    if figdir and not os.path.isdir(figdir):
        os.makedirs(figdir, exist_ok=True)

    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(5, 9)

    # image cutout, only if this photometry configuration has viewer info
    if hasattr(phot, 'viewer_layer') and hasattr(phot, 'viewer_pixscale'):
        img, wcs, _, _ = _fetch_cutout(meta, figdir, pngfile, phot.viewer_layer, phot.viewer_pixscale)
        cutax = fig.add_subplot(gs[0:2, 6:9], projection=wcs)
        cutax.imshow(img, origin='lower')
        cutax.set_xlabel('RA')
        cutax.set_ylabel('Dec')

    # SED panel: observed photometry vs. the refined maximum-likelihood model
    sedax = fig.add_subplot(gs[0:2, 0:6])
    lambda_obs = np.asarray(data['photometry']['lambda_eff']) / 1e4 # [micron]
    flam = np.asarray(data['photometry']['flam'])
    flam_ivar = np.asarray(data['photometry']['flam_ivar'])
    good = flam_ivar > 0.
    ferr = np.zeros_like(flam)
    ferr[good] = 1. / np.sqrt(flam_ivar[good])

    sedax.plot(zwave / 1e4, zflux, color='firebrick', alpha=0.8, lw=1, label='Best-fit model', zorder=1)
    sedax.errorbar(lambda_obs[good], flam[good], yerr=ferr[good], fmt='o', color='k',
                   markersize=6, label='Observed', zorder=2)
    sedax.set_xscale('log')
    sedax.set_xlabel(r'Observed-frame wavelength ($\mu$m)')
    sedax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    sedax.legend(fontsize=9, loc='best')
    sedax.set_title(' / '.join(_target_label(meta, coadd_type)) + f'  (z={redshift:.4f})', fontsize=10)

    # posterior panel: weighted 1D marginal histogram per parameter
    nparam = len(PARAM_NAMES)
    ncols = 3
    nrows = int(np.ceil(nparam / ncols))
    post_gs = gs[2:5, 0:9].subgridspec(nrows, ncols, hspace=0.5, wspace=0.3)
    for i, pname in enumerate(PARAM_NAMES):
        ax = fig.add_subplot(post_gs[i // ncols, i % ncols])
        vals, w = posterior_arrays[pname]
        lo, hi = np.min(vals), np.max(vals)
        if hi > lo:
            ax.hist(vals, bins=20, range=(lo, hi), weights=w, color='gray', edgecolor='k', alpha=0.8)
        ax.axvline(result[pname], color='C0', lw=1.5)
        ax.set_xlabel(pname, fontsize=9)
        ax.tick_params(labelleft=False, labelsize=7)

    fig.tight_layout()
    fig.savefig(pngfile)
    plt.close(fig)
    log.info(f'Wrote {pngfile}')


def write_fastbayes(meta, results, modelwave, modelspectra, outfile, gridfile, fphotofile, topk=0):
    """Write the Bayesian-fitting output to a multi-extension FITS file.

    Parameters
    ----------
    meta : :class:`astropy.table.Table`
        Output metadata table (see :func:`fastspecfit.io.create_output_meta`).
    results : :class:`astropy.table.Table`
        Output fitting-results table.
    modelwave : :class:`numpy.ndarray`
        Shared rest-frame wavelength array for ``modelspectra``.
    modelspectra : :class:`numpy.ndarray`
        Refined rest-frame maximum-likelihood spectra, shape (nobj, nwave).
    outfile : :class:`str`
        Full path of the output FITS file (``.gz`` triggers gzip
        compression, as in :func:`fastspecfit.io.write_fastspecfit`).
    gridfile : :class:`str`
        Full path of the Bayesian grid file used for fitting.
    fphotofile : :class:`str`
        Full path of the photometry configuration file used.
    topk : :class:`int`, optional
        Number of top-weight grid templates stored per object, if any.

    """
    import gzip, shutil
    from astropy.io import fits

    outdir = os.path.dirname(os.path.abspath(os.path.expanduser(os.path.expandvars(outfile))))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    if outfile.endswith('.gz'):
        tmpfile = outfile[:-3] + '.tmp'
    else:
        tmpfile = outfile + '.tmp'

    hduprim = fits.PrimaryHDU()
    hduprim.header['GRIDFILE'] = os.path.abspath(str(gridfile))
    hduprim.header['FPHOTO'] = os.path.abspath(str(fphotofile)) if fphotofile else ''
    hduprim.header['TOPK'] = topk

    hdumeta = fits.convenience.table_to_hdu(meta)
    hdumeta.header['EXTNAME'] = 'METADATA'

    hduresults = fits.convenience.table_to_hdu(results)
    hduresults.header['EXTNAME'] = 'FASTBAYES'

    hduwave = fits.ImageHDU(modelwave.astype('f4'))
    hduwave.header['EXTNAME'] = 'WAVE'
    hduwave.header['BUNIT'] = 'Angstrom'
    hduwave.header['AIRORVAC'] = ('vac', 'vacuum wavelengths')

    hdumodels = fits.ImageHDU(modelspectra)
    hdumodels.header['EXTNAME'] = 'MODELS'
    hdumodels.header['BUNIT'] = 'erg/(s cm2 Angstrom)'
    hdumodels.header['COMMENT'] = 'rest-frame refined maximum-likelihood model spectra'

    hx = fits.HDUList([hduprim, hdumeta, hduresults, hduwave, hdumodels])

    nobj = len(meta)
    if nobj == 1:
        log.info(f'Writing 1 object to {outfile}')
    else:
        log.info(f'Writing {nobj:,d} objects to {outfile}')

    hx.writeto(tmpfile, overwrite=True, checksum=True)

    # compress via another tempfile if needed, otherwise just rename; either
    # way the final `outfile` only ever appears atomically
    if outfile.endswith('.gz'):
        tmpfilegz = outfile[:-3] + '.tmp.gz'
        with open(tmpfile, 'rb') as f_in:
            with gzip.open(tmpfilegz, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.rename(tmpfilegz, outfile)
        os.remove(tmpfile)
    else:
        os.rename(tmpfile, outfile)


def parse(options=None):
    """Parse input arguments to the ``fastbayes`` script.

    Parameters
    ----------
    options : list of str or None, optional
        Command-line argument strings. If ``None``, reads from ``sys.argv``.

    Returns
    -------
    args : :class:`argparse.Namespace`

    """
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('redrockfiles', nargs='+', help='Full path to input redrock file(s).')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename (required).')
    parser.add_argument('--gridfile', type=str, required=True,
                        help='Full path to the Bayesian photometry grid file '
                        '(output of bin/build-bayesian-photometry).')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to process in each file, zero-indexed.')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--input-redshifts', type=str, default=None,
                        help='Comma-separated list of input redshifts corresponding to the (required) --targetids input.')
    parser.add_argument('--zmin', type=float, default=None, help='Override the default minimum redshift required for modeling.')
    parser.add_argument('--topk', type=int, default=0,
                        help='Number of top-weight grid templates to save per object (sparse joint posterior); 0 disables.')
    parser.add_argument('--qa', action='store_true', help='Generate a QA figure for each object.')
    parser.add_argument('--qadir', type=str, default='.', help='Output directory for QA figures.')
    parser.add_argument('--use-quasarnet', dest='use_quasarnet', default=False, action='store_true',
                        help='Use QuasarNet to improve QSO redshifts.')
    parser.add_argument('--fphotodir', type=str, default=None, help='Top-level location of the source photometry.')
    parser.add_argument('--fphotofile', type=str, default=None,
                        help='Photometric configuration file (default: taken from the grid file header).')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--redux_dir', type=str, default=None, help='Optional full path $DESI_SPECTRO_REDUX.')
    parser.add_argument('--specproddir', type=str, default=None, help='Optional directory name for the spectroscopic production.')
    parser.add_argument('--uncertainty-floor', type=float, default=0.01,
                        help='Minimum fractional uncertainty to add in quadrature to the formal inverse variance.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    log.info('fastbayes {}'.format(' '.join(options)))

    return parser.parse_args(options)


def fastbayes(args=None, mp_pool=None):
    """Main fastbayes engine: read, fit, and write Bayesian grid photometric results.

    Parameters
    ----------
    args : :class:`argparse.Namespace` or list of str or None, optional
        Pre-parsed arguments or raw argument list. If ``None``, reads from
        ``sys.argv``.
    mp_pool : :class:`fastspecfit.util.MPPool` or None, optional
        Pre-built worker pool; a new one is created (and closed) when
        ``None``.

    Returns
    -------
    :class:`int`
        Exit code (0 on success).

    """
    from astropy.table import vstack
    from fastspecfit.io import DESISpectra, create_output_meta, create_output_table

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    envlist = []
    if args.redux_dir is None:
        envlist.append('DESI_SPECTRO_REDUX')
    if args.mapdir is None:
        envlist.append('DUST_DIR')
    if args.fphotodir is None:
        envlist.append('FPHOTO_DIR')
    for env in envlist:
        if env not in os.environ:
            errmsg = f'Mandatory environment variable {env} missing.'
            log.critical(errmsg)
            raise KeyError(errmsg)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    targetids = None
    input_redshifts = None
    if args.targetids is not None:
        targetids = [int(x) for x in args.targetids.split(',')]
        if args.input_redshifts is not None:
            input_redshifts = [float(x) for x in args.input_redshifts.split(',')]
            if len(input_redshifts) != len(targetids):
                errmsg = 'targetids and input_redshifts must have the same number of elements.'
                log.critical(errmsg)
                raise ValueError(errmsg)

    gridfile = os.path.expandvars(args.gridfile)
    if not os.path.isfile(gridfile):
        errmsg = f'Bayesian grid file {gridfile} not found.'
        log.critical(errmsg)
        raise IOError(errmsg)

    fphotofile = args.fphotofile
    if fphotofile is None:
        with fitsio.FITS(gridfile) as F:
            fphotofile = F[0].read_header().get('FPHOTO')

    init_argdict = {'fphotofile': fphotofile, 'gridfile': gridfile}

    t0 = time.time()
    _initialize_fastbayes_workers(**init_argdict)
    log.info(fsftime('bg_data_init', time.time() - t0))

    _own_pool = False
    if mp_pool is None:
        mp_pool = MPPool(args.mp, initializer=_initialize_fastbayes_workers, init_argdict=init_argdict)
        _own_pool = True

    phot = sc_data.photometry

    Spec = DESISpectra(phot=phot, cosmo=_cosmo, fphotodir=args.fphotodir,
                       mapdir=args.mapdir, redux_dir=args.redux_dir)

    Spec.gather_metadata(args.redrockfiles, firsttarget=args.firsttarget,
                         targetids=targetids, input_redshifts=input_redshifts,
                         ntargets=args.ntargets, zmin=args.zmin,
                         redrockfile_prefix=args.redrockfile_prefix,
                         use_quasarnet=args.use_quasarnet, specprod_dir=args.specproddir)
    if len(Spec.specfiles) == 0:
        return 0

    data, meta = Spec.read(phot, fastphot=True)

    nobj = len(meta)
    fastbayes_dtype = get_fastbayes_dtype(phot, topk=args.topk)

    fitargs = [{
        'iobj': iobj,
        'data': data[iobj],
        'meta': meta[iobj],
        'fastbayes_dtype': fastbayes_dtype,
        'topk': args.topk,
        'uncertainty_floor': args.uncertainty_floor,
        'qa': args.qa,
        'qadir': args.qadir,
        'coadd_type': Spec.coadd_type,
    } for iobj in range(nobj)]

    t0 = time.time()
    out = mp_pool.starmap(fastbayes_one, fitargs)
    out = list(zip(*out))

    outmeta = create_output_meta(vstack(out[0]), phot=phot, fastphot=True)

    units = {'LOGAGE': 'dex(Gyr)'}
    units.update({f'LOGAGE_{stat}': 'dex(Gyr)' for stat in ('MEAN', 'MODE', 'P16', 'P50', 'P84')})
    results = create_output_table(out[1], outmeta, units)

    modelwave = bg_data.template_wave()
    modelspectra = np.array(out[2])

    if _own_pool:
        mp_pool.close()

    log.info(fsftime('fastbayes_all', time.time() - t0, context=f'nobj={nobj}'))

    write_fastbayes(outmeta, results, modelwave, modelspectra, outfile=args.outfile,
                    gridfile=gridfile, fphotofile=phot.fphotofile, topk=args.topk)

    return 0
