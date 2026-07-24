"""
fastspecfit.fastbayes
=====================

Grid-based Bayesian broadband-photometric SED fitting.

Given a pre-synthesized photometry grid (see ``bin/build-bayesian-templates``
and ``bin/build-bayesian-photometry``), fit each object's observed broadband
photometry by interpolating the grid to the object's (Redrock) redshift,
solving for the chi2-minimizing stellar-mass amplitude of every grid template
in closed form, and building up the posterior probability distribution of
every grid parameter by weighting each template by its likelihood. The
discrete maximum-likelihood point is then refined to sub-grid precision via
a local parabola fit along each of the 7 grid axes and an N-linear
interpolation over the neighboring templates, so that the reported
parameters, CHI2, and model spectrum are all mutually consistent.

See ``doc/technote/fastbayes.tex`` for the full grid design, statistical
methodology, sub-grid refinement algorithm, and derived-quantity
definitions.

"""
import os, sys, time, logging
import numpy as np
import fitsio
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.util import MPPool, fsftime, ZWarningMask, C_LIGHT
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

# All reported parameters: the 7 grid axes, plus LOGMSTAR and LOGSFR, which
# are derived per-object from the closed-form amplitude solve.
PARAM_NAMES = tuple(_AXIS_OUTNAME[col] for col in GRID_AXIS_COLUMNS) + ('LOGMSTAR', 'LOGSFR')

# Human-friendly axis labels for the QA posterior-histogram grid.
_PARAM_LABELS = {
    'LOGAGE': r'$\log_{10}({\rm Age/Gyr})$',
    'LOGZZSUN': r'$\log_{10}(Z/Z_{\odot})$',
    'TAU': r'$\tau$',
    'DUSTN': r'$n_{\rm dust}$',
    'LOGUMIN': r'$\log_{10}(U_{\rm min})$',
    'GAMMA': r'$\gamma$',
    'LOGQPAH': r'$\log_{10}(Q_{\rm PAH})$',
    'LOGMSTAR': r'$\log_{10}(M/M_{\odot})$',
    'LOGSFR': r'$\log_{10}({\rm SFR}/[M_{\odot}/{\rm yr}])$',
}

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
        self._axis_cache = {}


    @staticmethod
    def get_templates_filename(gridnumber, imf, templatedir=None):
        """Build the raw templates FITS filename from ``FTEMPLATES_DIR`` (or an
        explicit override), grid number, and IMF, mirroring
        :func:`fastspecfit.templates.Templates.get_templates_filename`.

        Returns ``None`` (rather than raising) when ``templatedir`` is not
        given and ``$FTEMPLATES_DIR`` is unset, so that QA-only workflows
        (:func:`fastbayes_qa`), which never need the raw templates file, are
        not forced to have it configured; :meth:`templates_fits` raises a
        clear error if raw template access is actually attempted without it.

        """
        if templatedir is None:
            templatedir = os.environ.get('FTEMPLATES_DIR')
            if templatedir is None:
                return None
            templatedir = os.path.join(os.path.expandvars(templatedir), 'bayesian')
        else:
            templatedir = os.path.expandvars(templatedir)
        return os.path.join(templatedir, str(gridnumber), f'bayesian-templates-{imf}-{gridnumber}.fits')


    def load(self, gridfile, templatedir=None):
        """Load a ``bayesian-photometry-*.fits`` grid file (idempotent)."""
        if self.file == gridfile:
            return

        T = fitsio.FITS(gridfile)
        prihdr = T[0].read_header()

        self.gridnumber = prihdr.get('GRIDNUM')
        self.imf = prihdr.get('IMF')
        self.fphotofile = prihdr.get('FPHOTO')
        self.templates_file = self.get_templates_filename(self.gridnumber, self.imf, templatedir=templatedir)
        self.logspace = bool(prihdr.get('LOGSPACE'))

        photsys_hdr = prihdr.get('PHOTSYS')
        self.photsys_keys = [''] if photsys_hdr in (None, 'NONE') else photsys_hdr.split(',')

        self.redshift = T['REDSHIFT'].read().astype('f8')
        self.zgrid_interp = np.log10(1. + self.redshift) if self.logspace else self.redshift

        self.meta = Table(T['METADATA'].read())
        self.ntemplate = len(self.meta)
        self.sfr = np.asarray(self.meta['sfr'], dtype='f8') # [ntemplate], Msun/yr per Msun formed
        self._axis_cache = {}

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


    def axis_posterior_cache(self, col):
        """Lazily cache per-grid-axis posterior statistics for grid-axis column
        ``col``. These depend only on ``self.meta`` (identical for every
        object fit against this grid), so they are computed once per worker
        rather than once per object.

        Returns
        -------
        vals : :class:`numpy.ndarray`
            Per-template values of ``col`` in its fit coordinate (log10 for
            :data:`LOG_AXES`, linear otherwise).
        order : :class:`numpy.ndarray`
            ``np.argsort(vals)``.
        sorted_vals : :class:`numpy.ndarray`
            ``vals[order]``.
        uniq : :class:`numpy.ndarray`
            Unique values of ``vals``, from ``np.unique``.
        inv : :class:`numpy.ndarray`
            Inverse-index mapping from ``np.unique(vals, return_inverse=True)``.

        """
        cached = self._axis_cache.get(col)
        if cached is None:
            vals = np.asarray(self.meta[col], dtype='f8')
            if col in LOG_AXES:
                vals = np.log10(vals)
            order = np.argsort(vals)
            sorted_vals = vals[order]
            uniq, inv = np.unique(vals, return_inverse=True)
            cached = (vals, order, sorted_vals, uniq, inv)
            self._axis_cache[col] = cached
        return cached


# global structures with single-copy data, initially empty
bg_data = BayesianGrid()
_igm = None
_cosmo = None


def _initialize_fastbayes_workers(fphotofile=None, gridfile=None, templatedir=None,
                                  require_templates=True, mapdir=None, cutout_unreachable=None):
    """MPPool initializer: populate ``sc_data.photometry``, ``bg_data``, the
    IGM model, and the cosmology in each worker.

    ``sc_data.photometry`` is populated directly (rather than via
    ``sc_data.initialize()``) so this mode never loads the stellar template
    basis or emission-line tables, neither of which it needs; the Milky Way
    dust-map directory is set the same way, via
    :meth:`~fastspecfit.singlecopy.Singletons.set_mapdir`, since
    :func:`fastspecfit.io.one_spectrum` (called from :func:`fastbayes_one`)
    still needs :attr:`sc_data.sfdmap` to be usable.

    Parameters
    ----------
    mapdir : :class:`str` or None, optional
        Directory containing the Milky Way dust maps; defaults to
        ``$DUST_DIR/maps`` when ``None``. Not needed by
        :func:`fastbayes_qa`/:func:`fastbayes_qa_one`, which never call
        :func:`fastspecfit.io.one_spectrum`.
    require_templates : :class:`bool`, optional
        If ``True`` (default), require the raw templates FITS file (used by
        :meth:`BayesianGrid.template_flux_row` to build the refined rest-frame
        spectrum during fitting) to be present. QA regeneration from an
        already-written FASTBAYES output file
        (:func:`fastbayes_qa`/:func:`fastbayes_qa_one`) never reads raw
        template rows -- the refined spectrum is read back from that file's
        ``MODELS`` extension instead -- so it passes ``False`` here to avoid
        requiring the (multi-GB) raw templates file to be present just to
        regenerate plots.
    cutout_unreachable : :class:`bool` or None, optional
        Seeds this worker's copy of :data:`fastspecfit.qa._cutout_unreachable`
        (same sticky per-process mechanism ``fastqa`` uses to detect an
        unreachable Legacy Survey viewer once, up front, rather than letting
        every object separately discover the outage via
        :func:`fastspecfit.qa._fetch_cutout`'s retry/backoff loop). Left
        untouched (``None``, the default) for the plain fitting path, which
        only imports ``fastspecfit.qa`` at all when ``--qa`` is requested.

    """
    global _igm, _cosmo

    sc_data.photometry = Photometry(fphotofile=fphotofile)
    sc_data.set_mapdir(mapdir)
    bg_data.load(gridfile, templatedir=templatedir)

    if cutout_unreachable is not None:
        import fastspecfit.qa as qa_module
        qa_module._cutout_unreachable = cutout_unreachable

    if require_templates and not os.path.isfile(bg_data.templates_file):
        errmsg = (f'Bayesian templates file {bg_data.templates_file} not found; '
                  'check $FTEMPLATES_DIR or --templatedir.')
        log.critical(errmsg)
        raise IOError(errmsg)

    if list(bg_data.bands) != list(sc_data.photometry.bands):
        errmsg = (f"Band mismatch between the photometry configuration {sc_data.photometry.fphotofile} "
                  f"and the Bayesian grid file {gridfile}, which was built with fphoto file "
                  f"'{bg_data.fphotofile}'; they must be built from the same (or a compatible) "
                  "fphoto configuration file.")
        log.critical(errmsg)
        raise ValueError(errmsg)

    if _igm is None:
        from fastspecfit.igm import Inoue14
        _igm = Inoue14()
    if _cosmo is None:
        from fastspecfit.cosmo import TabulatedDESI
        _cosmo = TabulatedDESI()


def _weighted_percentile(values, weights, percentiles, order=None, sorted_values=None):
    """Weighted percentiles of ``values`` via linear interpolation of the weighted CDF.

    ``order``/``sorted_values`` (``np.argsort(values)`` and ``values[order]``)
    may be passed in when already known -- e.g. because ``values`` is
    grid-only data, identical for every object -- to skip the O(n log n)
    argsort.

    """
    if order is None:
        order = np.argsort(values)
        sorted_values = values[order]
    w = weights[order]
    cw = np.cumsum(w)
    cw /= cw[-1]
    return np.interp(np.asarray(percentiles) / 100., cw, sorted_values)


def _weighted_mode(values, weights, uniq=None, inv=None):
    """Weighted mode: the (exact) grid value with the largest total weight.

    ``uniq``/``inv`` (from ``np.unique(values, return_inverse=True)``) may be
    passed in when already known -- e.g. because ``values`` is grid-only
    data, identical for every object -- to skip the redundant call.

    """
    if uniq is None:
        uniq, inv = np.unique(values, return_inverse=True)
    binweights = np.bincount(inv, weights=weights)
    return uniq[np.argmax(binweights)]


def _axis_posterior_values(col):
    """Per-template values of grid-axis column ``col``, in its fit coordinate
    (log10 for :data:`LOG_AXES`, linear otherwise)."""
    return bg_data.axis_posterior_cache(col)[0]


def _parabola_vertex3(x, y):
    """Closed-form vertex of the unique parabola through 3 ``(x, y)`` points.

    Equivalent to :func:`fastspecfit.util.minfit` specialized to exactly 3
    points: with 3 points and a 2nd-degree polynomial the fit is an exact
    interpolation (not a least-squares approximation), so the coefficients
    have a direct algebraic solution (Newton's divided-difference form)
    instead of paying for a generic least-squares/SVD solve
    (:func:`numpy.polyfit`) on every call. ``x`` is assumed strictly
    increasing, guaranteed by its grid-neighbor construction in
    :func:`_refine_grid_axes`.

    Returns
    -------
    x0, xerr, y0, zwarn : Same meaning and edge-case semantics as
        :func:`fastspecfit.util.minfit`.

    """
    x0v, x1v, x2v = x
    y0v, y1v, y2v = y

    f01 = (y1v - y0v) / (x1v - x0v)
    f12 = (y2v - y1v) / (x2v - x1v)
    a = (f12 - f01) / (x2v - x0v)
    b = f01 - a * (x0v + x1v)
    c = y0v - a * x0v**2 - b * x0v

    if a == 0.:
        return -1, -1, -1, ZWarningMask.BAD_MINFIT

    x0 = -b / (2. * a)
    y0 = -(b**2) / (4. * a) + c

    zwarn = 0
    if (x0 <= x0v) or (x2v <= x0):
        zwarn |= ZWarningMask.BAD_MINFIT
    if y0 <= 0.:
        zwarn |= ZWarningMask.BAD_MINFIT

    if a > 0.:
        xerr = 1. / np.sqrt(a)
    else:
        xerr = 1. / np.sqrt(-a)
        zwarn |= ZWarningMask.BAD_MINFIT

    return x0, xerr, y0, zwarn


def _solve_grid(flam, flam_ivar, lambda_eff, photsys, redshift):
    """Closed-form chi2-minimizing amplitude solve over the full grid, plus
    sub-grid refinement of the maximum-likelihood point.

    Factored out of :func:`fastbayes_one` so that it can also be called by
    :func:`fastbayes_qa_one` to regenerate QA figures from an
    already-written FASTBAYES output file, without repeating the fit's
    per-object I/O (only ``flam``/``flam_ivar``/``lambda_eff``/``photsys``/
    ``redshift`` are needed, all recoverable from the output file's
    ``METADATA`` table) -- this keeps the two callers' math from being able
    to drift apart.

    Parameters
    ----------
    flam, flam_ivar, lambda_eff : :class:`numpy.ndarray`
        Observed flux density, inverse variance, and effective wavelength,
        one entry per band.
    photsys : :class:`str`
        Photometric-system key.
    redshift : :class:`float`
        Object redshift.

    Returns
    -------
    :class:`dict`
        ``chi2``, ``weight``, ``amplitude`` (all shape (ntemplate,)),
        ``ibest`` (flat index of the discrete chi2 minimum), ``fit_value``/
        ``fit_ivar`` (per grid-axis-column refined value/formal ivar, from
        :func:`_refine_grid_axes`), ``corner_idx``/``corner_weight`` (N-linear
        interpolation corners, from :func:`_corner_weights`),
        ``refined_model_maggies`` (shape (nband,)), ``refined_amplitude``,
        and ``chi2_refined``.

    """
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

    return {
        'chi2': chi2, 'weight': weight, 'amplitude': amplitude, 'ibest': ibest,
        'fit_value': fit_value, 'fit_ivar': fit_ivar,
        'corner_idx': corner_idx, 'corner_weight': corner_weight,
        'refined_model_maggies': refined_model_maggies,
        'refined_amplitude': refined_amplitude, 'chi2_refined': chi2_refined,
    }


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

        x0, xerr, y0, zwarn = _parabola_vertex3(xs, ys)

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
    multi_index = np.array(np.unravel_index(ibest, bg_data.dims))
    naxes = len(GRID_AXIS_COLUMNS)
    t = np.array([frac[col] for col in GRID_AXIS_COLUMNS])
    direction = np.array([frac_dir[col] for col in GRID_AXIS_COLUMNS])

    # bits[corner, axis_pos] is the axis_pos-th bit of `corner`, for every
    # corner of the 2**naxes hypercube at once.
    bits = (np.arange(2**naxes)[:, np.newaxis] >> np.arange(naxes)[np.newaxis, :]) & 1

    per_axis_weight = np.where(bits == 1, t[np.newaxis, :], 1. - t[np.newaxis, :])
    corner_weight = np.prod(per_axis_weight, axis=1)

    keep = corner_weight > 0.
    corner_multi_index = multi_index[np.newaxis, :] + bits[keep] * direction[np.newaxis, :]

    indices = np.ravel_multi_index(tuple(corner_multi_index.T), bg_data.dims).astype('i8')
    weights = corner_weight[keep]
    weights = weights / weights.sum()

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

    # Formal delta-chi2=1 uncertainties for the 7 native grid axes.
    for col in GRID_AXIS_COLUMNS:
        cols += [(f'{_AXIS_OUTNAME[col]}_IVAR', 'f4')]

    # LOGMSTAR/LOGSFR are derived quantities without a directly fit chi2(x)
    # curve, so their uncertainty is instead estimated from the weighted
    # variance of their discrete-grid posterior.
    cols += [('LOGMSTAR_IVAR', 'f4'), ('LOGSFR_IVAR', 'f4')]

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

    soln = _solve_grid(flam, flam_ivar, lambda_eff, photsys, redshift)
    chi2 = soln['chi2']
    weight = soln['weight']
    amplitude = soln['amplitude']
    ibest = soln['ibest']
    fit_value = soln['fit_value']
    fit_ivar = soln['fit_ivar']
    corner_idx = soln['corner_idx']
    corner_weight = soln['corner_weight']
    refined_model_maggies = soln['refined_model_maggies']
    refined_amplitude = soln['refined_amplitude']
    chi2_refined = soln['chi2_refined']

    logmstar = np.log10(np.clip(amplitude, 1e-30, None))

    # The grid stores sfr per solar mass formed (like the flux templates
    # themselves), so scale by the fitted mass amplitude to get each
    # template's actual SFR estimate for this object. Guard against
    # zero/negative SFR (e.g., old, passively evolving SSPs) before taking
    # the log, given LOGSFR's large dynamic range.
    sfr_per_template = amplitude * bg_data.sfr # [ntemplate], Msun/yr
    logsfr_per_template = np.log10(np.clip(sfr_per_template, 1e-30, None))

    refined_logmstar = np.log10(max(refined_amplitude, 1e-30))
    refined_sfr_per_msun = np.sum(corner_weight * bg_data.sfr[corner_idx])
    refined_sfr = refined_amplitude * refined_sfr_per_msun
    refined_logsfr = np.log10(max(refined_sfr, 1e-30))

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
    derived = {'LOGMSTAR': refined_logmstar, 'LOGSFR': refined_logsfr}
    posterior = {'LOGMSTAR': logmstar, 'LOGSFR': logsfr_per_template}

    for col in GRID_AXIS_COLUMNS:
        pname = _AXIS_OUTNAME[col]
        result[pname] = fit_value[col]
        result[f'{pname}_IVAR'] = fit_ivar[col]

    posterior_arrays = {}
    for pname in PARAM_NAMES:
        order = sorted_vals = uniq = inv = None
        if pname in derived:
            result[pname] = derived[pname]
            vals = posterior[pname]
        else:
            # Grid-only data, identical for every object -- fetch the
            # cached argsort/unique instead of recomputing them here.
            vals, order, sorted_vals, uniq, inv = bg_data.axis_posterior_cache(_OUTNAME_TO_AXIS[pname])
        posterior_arrays[pname] = (vals, weight)
        mean = np.sum(weight * vals)
        result[f'{pname}_MEAN'] = mean
        result[f'{pname}_MODE'] = _weighted_mode(vals, weight, uniq=uniq, inv=inv)
        p16, p50, p84 = _weighted_percentile(vals, weight, (16., 50., 84.), order=order, sorted_values=sorted_vals)
        result[f'{pname}_P16'] = p16
        result[f'{pname}_P50'] = p50
        result[f'{pname}_P84'] = p84

        # LOGMSTAR/LOGSFR have no direct grid axis (and hence no delta-chi2=1
        # parabola fit), so estimate their formal uncertainty from the
        # weighted variance of their discrete-grid posterior instead.
        if pname in derived:
            var = np.sum(weight * (vals - mean)**2)
            result[f'{pname}_IVAR'] = 1. / var if var > 0. else 0.

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
        synth_maggies = refined_model_maggies * refined_amplitude # [nband], observed-frame maggies
        _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                          zwave, zflux, synth_maggies, redshift, coadd_type=coadd_type, outdir=qadir)

    return meta, result, restflux.astype('f4')


def fastbayes_qa_one(iobj, meta, result, restwave, restflux, qadir='.', coadd_type='healpix'):
    """Regenerate the QA figure for one object from an already-written
    FASTBAYES output file, without repeating the fit.

    Rebuilds everything :func:`_fastbayes_qa_one` needs purely from data
    already present in the output file's ``METADATA``/``FASTBAYES``/``WAVE``/
    ``MODELS`` extensions: the observed photometry (``FLUX_*``/
    ``FLUX_IVAR_*`` columns), redshift (``Z``) and photometric system
    (``PHOTSYS``), and the refined rest-frame spectrum (``restwave``/
    ``restflux``, read back rather than recomputed). The only work redone is
    the cheap, vectorized full-grid solve (:func:`_solve_grid`) needed to
    reconstruct the per-template posterior weights, which are not stored in
    the output file -- no raw DESI spectra or individual template rows are
    read.

    Parameters
    ----------
    iobj : :class:`int`
        Index of the object in the input list, used for log messages.
    meta : :class:`astropy.table.Row`
        ``METADATA`` row for this object, from the FASTBAYES output file.
    result : :class:`numpy.void`
        ``FASTBAYES`` results row for this object, from the FASTBAYES
        output file.
    restwave : :class:`numpy.ndarray`
        Shared rest-frame wavelength array (the output file's ``WAVE``
        extension).
    restflux : :class:`numpy.ndarray`
        This object's refined rest-frame spectrum (one row of the output
        file's ``MODELS`` extension).
    qadir : :class:`str`, optional
        Output directory for the QA figure. Default is ``'.'``.
    coadd_type : :class:`str`, optional
        Coadd type, used to build the QA target label/filename. Not stored
        in the output file, so it must be supplied by the caller. Default
        is ``'healpix'``.

    """
    phot = sc_data.photometry

    redshift = float(meta['Z'])
    if redshift <= 0.:
        # mirror the (non-persisted) floor applied to the local `redshift`
        # value in fastspecfit.io.DESISpectra.read
        redshift = 1e-8
    photsys = str(meta['PHOTSYS']).strip()

    log.info(f'Regenerating QA for object {iobj} [{phot.uniqueid_col.lower()} '
            f'{meta[phot.uniqueid_col]}, z={redshift:.6f}].')

    if redshift > bg_data.redshift[-1]:
        log.warning(f'Object {iobj} [{phot.uniqueid_col.lower()} {meta[phot.uniqueid_col]}] redshift '
                   f'{redshift:.6f} exceeds the grid maximum {bg_data.redshift[-1]:.6f}; skipping QA.')
        return

    nanomaggies = np.array([meta[f'FLUX_{band.upper()}'] for band in phot.bands], dtype='f8')
    nanomaggies_ivar = np.array([meta[f'FLUX_IVAR_{band.upper()}'] for band in phot.bands], dtype='f8')
    lambda_eff = phot.filters[photsys].effective_wavelengths.value

    phot_tbl = Photometry.parse_photometry(
        phot.bands, maggies=nanomaggies, ivarmaggies=nanomaggies_ivar,
        lambda_eff=lambda_eff, min_uncertainty=phot.min_uncertainty)
    flam = np.asarray(phot_tbl['flam'], dtype='f8')
    flam_ivar = np.asarray(phot_tbl['flam_ivar'], dtype='f8') * phot.bands_to_fit

    soln = _solve_grid(flam, flam_ivar, lambda_eff, photsys, redshift)
    weight = soln['weight']
    amplitude = soln['amplitude']
    refined_model_maggies = soln['refined_model_maggies']
    refined_amplitude = soln['refined_amplitude']

    logmstar = np.log10(np.clip(amplitude, 1e-30, None))
    sfr_per_template = amplitude * bg_data.sfr # [ntemplate], Msun/yr
    logsfr_per_template = np.log10(np.clip(sfr_per_template, 1e-30, None))
    posterior = {'LOGMSTAR': logmstar, 'LOGSFR': logsfr_per_template}

    posterior_arrays = {}
    for pname in PARAM_NAMES:
        if pname in posterior:
            vals = posterior[pname]
        else:
            vals = bg_data.axis_posterior_cache(_OUTNAME_TO_AXIS[pname])[0]
        posterior_arrays[pname] = (vals, weight)

    zwave = restwave * (1. + redshift)
    igm_trans = _igm.full_IGM(redshift, zwave)
    dlum = _cosmo.luminosity_distance(redshift) # [Mpc]
    zfactor = igm_trans * (10. / (1e6 * dlum))**2 / (1. + redshift)
    zflux = restflux * zfactor # [erg/s/cm2/A, observed frame]

    synth_maggies = refined_model_maggies * refined_amplitude # [nband], observed-frame maggies

    data = {'photometry': {
        'nanomaggies': nanomaggies,
        'nanomaggies_ivar': nanomaggies_ivar,
        'lambda_eff': lambda_eff,
    }}

    _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                      zwave, zflux, synth_maggies, redshift, coadd_type=coadd_type, outdir=qadir)


def _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                      zwave, zflux, synth_maggies, redshift, coadd_type='healpix', outdir='.'):
    """Generate a QA figure for one Bayesian-fit object.

    Reuses :func:`fastspecfit.qa._target_label` and
    :func:`fastspecfit.qa._fetch_cutout`, which are generic utilities with
    no fastspec/fastphot-specific coupling. Called either inline from
    :func:`fastbayes_one` during fitting (``--qa``), or from
    :func:`fastbayes_qa_one` to regenerate QA from an already-written
    FASTBAYES output file without repeating the fit -- both callers build
    the same ``posterior_arrays``/``zwave``/``zflux``/``synth_maggies``
    inputs, just from different sources.

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
    synth_maggies : :class:`numpy.ndarray`
        Observed-frame photometry (AB maggies) synthesized from the refined
        maximum-likelihood model in each band.
    redshift : :class:`float`
    coadd_type : :class:`str`, optional
    outdir : :class:`str`, optional

    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib import colors
    from matplotlib.patches import Circle
    from matplotlib.lines import Line2D
    import seaborn as sns

    from fastspecfit.io import get_qa_filename
    from fastspecfit.qa import _target_label, _fetch_cutout

    sns.set(context='talk', style='ticks', font_scale=1.1)

    phot = sc_data.photometry
    phot_wavelims = (0.1, 35.) # [micron]

    photcol1 = colors.to_hex('darkorange')
    fontsize1, fontsize2 = 14, 20
    legxpos, legypos2, legfntsz1, legfntsz = 0.98, 0.05, 16, 18
    bbox = dict(boxstyle='round', facecolor='lightgray', alpha=0.15)
    bbox2 = dict(boxstyle='round', facecolor='lightgray', alpha=0.7)

    @ticker.FuncFormatter
    def major_formatter(x, pos):
        if (x >= 0.01) and (x < 0.1):
            return f'{x:.2f}'
        elif (x >= 0.1) and (x < 1):
            return f'{x:.1f}'
        else:
            return f'{x:.0f}'

    def _fmt(val, ivar, fmt):
        if ivar > 0.:
            return fmt.format(val) + r'\pm' + fmt.format(1. / np.sqrt(ivar))
        return fmt.format(val)

    pngfile = get_qa_filename(meta, coadd_type, outprefix='fastbayes', outdir=outdir)
    figdir = os.path.dirname(pngfile)
    if figdir and not os.path.isdir(figdir):
        os.makedirs(figdir, exist_ok=True)

    # 7 rows x 8 cols: rows 0:3 are the SED/cutout block (matching the
    # column proportions of fastqa's fastphot layout: sedax 5/8, cutax
    # 3/8, so labels/legends sized for that layout still fit), row 3 is
    # a (shrunk) blank gap, and rows 4:7 hold the posterior-histogram grid.
    # Explicit margins mirror fastqa's own fastphot subplots_adjust (tight
    # left/bottom, generous right/top for the fig.text labels/legends).
    fig = plt.figure(figsize=(18, 15))
    gs = fig.add_gridspec(7, 8, height_ratios=[1, 1, 1, 0.2, 1.1, 1.1, 1.1],
                          left=0.09, right=0.88, top=0.9, bottom=0.06,
                          hspace=0.6, wspace=0.3)

    sedax = fig.add_subplot(gs[0:3, 0:5])

    # image cutout, only if this photometry configuration has viewer info
    have_cutout = hasattr(phot, 'viewer_layer') and hasattr(phot, 'viewer_pixscale')
    if have_cutout:
        img, wcs, _, _ = _fetch_cutout(meta, figdir, pngfile, phot.viewer_layer, phot.viewer_pixscale)
        cutax = fig.add_subplot(gs[0:2, 5:8], projection=wcs)
        cutax.imshow(img, origin='lower')
        cutax.set_xlabel('RA [J2000]')
        cutax.set_ylabel('Dec [J2000]')
        cutax.invert_yaxis()
        cutax.coords[1].set_ticks_position('r')
        cutax.coords[1].set_ticklabel_position('r')
        cutax.coords[1].set_axislabel_position('r')
        sgn = '+' if meta['DEC'] > 0 else ''
        cutax.text(0.04, 0.95, '$(\\alpha,\\delta)$=({:.7f}, {}{:.6f})'.format(
                   meta['RA'], sgn, meta['DEC']), ha='left', va='top', color='k',
                   fontsize=fontsize1, bbox=bbox2, transform=cutax.transAxes)
        sz = img.shape
        pixscale = phot.viewer_pixscale
        cutax.add_artist(Circle((sz[1]/2, sz[0]/2), radius=1.5/2/pixscale,
                                facecolor='none', edgecolor='yellow', ls='-', alpha=0.8))
        cutax.add_artist(Circle((sz[1]/2, sz[0]/2), radius=10/2/pixscale,
                                facecolor='none', edgecolor='yellow', ls='--', alpha=0.8))
        handles = [Line2D([0], [0], color='yellow', lw=2, ls='-', label='1.5 arcsec'),
                  Line2D([0], [0], color='yellow', lw=2, ls='--', label='10 arcsec')]
        cutax.legend(handles=handles, loc='lower left', fontsize=fontsize1, facecolor='lightgray')
    else:
        cutax = fig.add_subplot(gs[0:2, 5:8])
        cutax.axis('off')

    # --- SED panel: observed photometry (filled markers, upper limits where
    # undetected) vs. the refined maximum-likelihood model (curve) and its
    # synthesized photometry (open markers), all in AB mag ------------------
    nanomaggies = np.asarray(data['photometry']['nanomaggies'], dtype='f8')
    nanomaggies_ivar = np.asarray(data['photometry']['nanomaggies_ivar'], dtype='f8')
    lambda_eff = np.asarray(data['photometry']['lambda_eff'], dtype='f8')

    phot_tbl = Photometry.parse_photometry(
        phot.bands, maggies=nanomaggies, ivarmaggies=nanomaggies_ivar,
        lambda_eff=lambda_eff, min_uncertainty=phot.min_uncertainty, get_abmag=True)

    abmag = np.asarray(phot_tbl['abmag'])
    abmag_limit = np.asarray(phot_tbl['abmag_limit'])
    abmag_ivar = np.asarray(phot_tbl['abmag_ivar'])
    yerr = np.vstack((np.asarray(phot_tbl['abmag_brighterr']), np.asarray(phot_tbl['abmag_fainterr'])))

    # best-fit model curve, converted from erg/s/cm2/A to AB mag; restrict to
    # the plotted wavelength range so that, e.g., far-IR dust emission well
    # beyond phot_wavelims doesn't skew the y-axis limits below.
    factor = 10**(0.4 * 48.6) * zwave**2 / (C_LIGHT * 1e13) # [erg/s/cm2/A --> maggies]
    zwave_um = zwave / 1e4
    mgood = (zflux > 0.) & (zwave_um >= phot_wavelims[0]) & (zwave_um <= phot_wavelims[1])
    sedmodel_abmag = np.full_like(zflux, np.nan)
    sedmodel_abmag[mgood] = -2.5 * np.log10(zflux[mgood] * factor[mgood])
    sedax.plot(zwave_um[mgood], sedmodel_abmag[mgood], color='grey', alpha=0.9, zorder=1)

    # synthesized photometry (open diamonds)
    synth_good = synth_maggies > 0.
    synth_abmag = np.full_like(synth_maggies, np.nan, dtype='f8')
    synth_abmag[synth_good] = -2.5 * np.log10(synth_maggies[synth_good])
    sedax.scatter(lambda_eff[synth_good] / 1e4, synth_abmag[synth_good], marker='D',
                 s=450, color='k', facecolor='none', linewidth=2, alpha=1.0, zorder=2)

    # observed photometry: filled markers (fitted bands) or hollow markers
    # (bands not used in the fit) for detections, upper limits (lolims)
    # for non-detections
    markersize = 14

    def _plot_obs(idx, facecolor, alpha):
        idx = np.asarray(idx)
        if len(idx) == 0:
            return
        good = idx[(abmag[idx] > 0.) & (abmag_limit[idx] == 0.)]
        upper = idx[abmag_limit[idx] > 0.]
        if len(good) > 0:
            sedax.errorbar(lambda_eff[good] / 1e4, abmag[good], yerr=yerr[:, good], fmt='o',
                           markersize=markersize, markeredgewidth=1, markeredgecolor='k',
                           markerfacecolor=facecolor, elinewidth=3, ecolor=photcol1,
                           capsize=4, alpha=alpha)
        if len(upper) > 0:
            sedax.errorbar(lambda_eff[upper] / 1e4, abmag_limit[upper], lolims=True, yerr=0.75,
                           fmt='o', markersize=markersize, markeredgewidth=3, markeredgecolor='k',
                           markerfacecolor=facecolor, elinewidth=3, ecolor=photcol1,
                           capsize=4, alpha=alpha)

    _plot_obs(np.where(phot.bands_to_fit)[0], photcol1, 1.0)
    _plot_obs(np.where(~phot.bands_to_fit)[0], 'none', 0.7)

    # y-axis limits (AB mag; brighter/smaller values at the top)
    dm = 1.5
    candidates = []
    if np.any(abmag_ivar > 0.):
        candidates.append(abmag[abmag_ivar > 0.])
    if np.any(abmag_limit > 0.):
        candidates.append(abmag_limit[abmag_limit > 0.])
    if np.any(mgood):
        candidates.append(sedmodel_abmag[mgood])
    if candidates:
        allmags = np.concatenate(candidates)
        sed_ymin = np.nanmax(allmags) + dm
        sed_ymax = np.nanmin(allmags) - dm
    else:
        sed_ymin, sed_ymax = 30., 20.

    sedax.set_xlim(phot_wavelims[0], phot_wavelims[1])
    sedax.set_ylim(sed_ymin, sed_ymax)
    sedax.set_xscale('log')
    sedax.set_ylabel('AB mag')
    sedax.set_xlabel(r'Observed-frame Wavelength ($\mu$m)')
    sedax.xaxis.set_major_formatter(major_formatter)
    obsticks = np.array([0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 20.0])
    obsticks = obsticks[(obsticks >= phot_wavelims[0]) * (obsticks <= phot_wavelims[1])]
    sedax.set_xticks(obsticks)

    sedax_twin = sedax.twiny()
    sedax_twin.set_xlim(phot_wavelims[0] / (1. + redshift), phot_wavelims[1] / (1. + redshift))
    sedax_twin.set_xscale('log')
    sedax_twin.xaxis.set_major_formatter(major_formatter)
    restticks = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10.0, 15.0, 20.0])
    restticks = restticks[(restticks >= phot_wavelims[0] / (1. + redshift)) *
                         (restticks <= phot_wavelims[1] / (1. + redshift))]
    sedax_twin.set_xticks(restticks)

    # reduced chi2 (top-left, no box)
    rchi2 = result['CHI2'] / result['NDOF'] if result['NDOF'] > 0 else result['CHI2']
    sedax.text(0.02, 0.94, r'$\chi^{2}_{\nu,\mathrm{phot}}=$' + r'${:.2f}$'.format(rchi2),
              ha='left', va='top', transform=sedax.transAxes, fontsize=legfntsz)

    # basic fitting info (bottom-right, shaded box)
    txt = '\n'.join((
        r'$\log_{{10}}(Z/Z_{{\odot}})={}$'.format(_fmt(result['LOGZZSUN'], result['LOGZZSUN_IVAR'], '{:.2f}')),
        r'$\tau={}$'.format(_fmt(result['TAU'], result['TAU_IVAR'], '{:.2f}')),
        r'$\log_{{10}}(\mathrm{{SFR}}/[M_{{\odot}}/\mathrm{{yr}}])={}$'.format(
            _fmt(result['LOGSFR'], result['LOGSFR_IVAR'], '{:.2f}')),
        r'$\log_{{10}}(\mathrm{{Age}}/\mathrm{{Gyr}})={}$'.format(
            _fmt(result['LOGAGE'], result['LOGAGE_IVAR'], '{:.2f}')),
        r'$\log_{{10}}(M/M_{{\odot}})={}$'.format(_fmt(result['LOGMSTAR'], result['LOGMSTAR_IVAR'], '{:.2f}')),
    ))
    sedax.text(legxpos, legypos2, txt, ha='right', va='bottom', transform=sedax.transAxes,
              fontsize=legfntsz1, bbox=bbox, linespacing=1.4)

    # target label above the cutout, and rest-frame wavelength label above the SED
    cpos = cutax.get_position()
    spos = sedax.get_position()
    fig.text(cpos.x0, cpos.y1 + 0.02, '\n'.join(_target_label(meta, coadd_type)),
             ha='left', va='bottom', fontsize=fontsize2, linespacing=1.4)
    fig.text((spos.x0 + spos.x1) / 2., spos.y1 + 0.05, r'Rest-frame Wavelength ($\mu$m)',
             ha='center', va='bottom', fontsize=fontsize2)

    # z / Dn(4000) and absolute-magnitude boxes below the cutout (below its
    # RA/Dec axis label and tick labels, not just its image edge). Anchoring
    # the left box at the column's left edge (ha='left') and the right box
    # at the column's right edge (ha='right') maximizes the gap between them
    # for a given cutout column width, rather than splitting at a fixed
    # fraction that can overlap if either box's text runs long.
    ytext = cpos.y0 - 0.06
    txt = [r'$z={:.7f}$'.format(redshift), '',
          r'$D_{{n}}(4000)_{{\mathrm{{model}}}}={:.3f}$'.format(result['DN4000_MODEL'])]
    fig.text(cpos.x0, ytext, '\n'.join(txt), ha='left', va='top', fontsize=legfntsz,
             bbox=bbox, linespacing=1.4)

    gindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1. + phot.band_shift) - 4300))
    rindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1. + phot.band_shift) - 5600))
    zindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1. + phot.band_shift) - 8100))
    absmag_gband, absmag_rband, absmag_zband = phot.absmag_bands[gindx], phot.absmag_bands[rindx], phot.absmag_bands[zindx]
    shift_gband, shift_rband, shift_zband = phot.band_shift[gindx], phot.band_shift[rindx], phot.band_shift[zindx]

    def _absmag_col(band, shift):
        return f'ABSMAG{int(10 * shift):02d}_{band.upper()}'

    txt = [r'$M_{{{}{}}}={:.2f}$'.format(
        str(shift_rband), absmag_rband.lower().replace('decam_', '').replace('sdss_', ''),
        result[_absmag_col(absmag_rband, shift_rband)])]
    if gindx != rindx:
        gr = result[_absmag_col(absmag_gband, shift_gband)] - result[_absmag_col(absmag_rband, shift_rband)]
        txt += [r'$M_{{{}{}}}-M_{{{}{}}}={:.3f}$'.format(
            str(shift_gband), absmag_gband.lower(), str(shift_rband), absmag_rband.lower(),
            gr).replace('decam_', '').replace('sdss_', '')]
    if zindx != rindx:
        rz = result[_absmag_col(absmag_rband, shift_rband)] - result[_absmag_col(absmag_zband, shift_zband)]
        txt += [r'$M_{{{}{}}}-M_{{{}{}}}={:.3f}$'.format(
            str(shift_rband), absmag_rband.lower(), str(shift_zband), absmag_zband.lower(),
            rz).replace('decam_', '').replace('sdss_', '')]
    fig.text(cpos.x1, ytext, '\n'.join(txt), ha='right', va='top', fontsize=legfntsz,
             bbox=bbox, linespacing=1.4)

    # --- posterior panel: weighted 1D marginal histogram per parameter -----
    nparam = len(PARAM_NAMES)
    ncols = 3
    nrows = int(np.ceil(nparam / ncols))
    post_gs = gs[4:7, 0:8].subgridspec(nrows, ncols, hspace=0.5, wspace=0.1)
    for i, pname in enumerate(PARAM_NAMES):
        ax = fig.add_subplot(post_gs[i // ncols, i % ncols])
        vals, w = posterior_arrays[pname]
        uniq, inv = np.unique(vals, return_inverse=True)

        if len(uniq) <= 5:
            # native grid axis: only a handful of discrete values exist, so
            # a bar per actual grid value avoids mostly-empty uniform bins
            binweight = np.bincount(inv, weights=w, minlength=len(uniq))
            width = np.min(np.diff(uniq)) * 0.8 if len(uniq) > 1 else 1.
            ax.bar(uniq, binweight, width=width, color='gray', edgecolor='k', alpha=0.8)
        else:
            # derived quantity (LOGMSTAR/LOGSFR): continuous across the
            # grid, so zoom the range to where the posterior weight
            # actually is rather than the full (often much wider) range
            lo, hi = _weighted_percentile(vals, w, (0.5, 99.5))
            if hi > lo:
                pad = 0.05 * (hi - lo)
                ax.hist(vals, bins=30, range=(lo - pad, hi + pad), weights=w,
                        color='gray', edgecolor='k', alpha=0.8)

        ax.axvline(result[pname], color='C0', lw=1.5)
        ax.set_xlabel(_PARAM_LABELS.get(pname, pname), fontsize=10)
        ax.tick_params(labelleft=False, labelsize=9)

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
    import gzip, shutil, warnings
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyUserWarning

    outdir = os.path.dirname(os.path.abspath(os.path.expanduser(os.path.expandvars(outfile))))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    if outfile.endswith('.gz'):
        tmpfile = outfile[:-3] + '.tmp'
    else:
        tmpfile = outfile + '.tmp'

    hduprim = fits.PrimaryHDU()
    hduprim.header['GRIDFILE'] = os.path.abspath(str(gridfile))
    # Basename only (not the full path): grids are often built on one
    # machine and fit/QA'd on another (e.g. laptop -> NERSC), so an absolute
    # path recorded here would not generally be reloadable. This is purely a
    # diagnostic/provenance record -- fastbayes_qa requires --fphotofile to
    # be passed explicitly rather than reading it back from here.
    hduprim.header['FPHOTO'] = os.path.basename(str(fphotofile)) if fphotofile else ''
    hduprim.header['TOPK'] = topk

    with warnings.catch_warnings():
        # e.g., 'dex(Gyr)' cannot be represented in native FITS format. This
        # is raised by fits.convenience.table_to_hdu as AstropyUserWarning
        # (not units.UnitsWarning, which is a different category and does
        # not actually match this warning); scope the filter to this
        # specific, expected message so other AstropyUserWarnings still
        # surface normally.
        warnings.filterwarnings('ignore', message="The unit .* could not be saved in native FITS format",
                               category=AstropyUserWarning)
        hdumeta = fits.convenience.table_to_hdu(meta)
        hduresults = fits.convenience.table_to_hdu(results)
    hdumeta.header['EXTNAME'] = 'METADATA'
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
    parser.add_argument('--templatedir', type=str, default=None,
                        help='Top-level location of the raw Bayesian templates file '
                        '(default: $FTEMPLATES_DIR/bayesian).')
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
                        help='Photometric configuration file (default: bundled fastspecfit configuration).')
    parser.add_argument('--redrockfile-prefix', type=str, default='redrock-', help='Prefix of the input Redrock file name(s).')
    parser.add_argument('--mapdir', type=str, default=None, help='Optional directory name for the dust maps.')
    parser.add_argument('--redux_dir', type=str, default=None, help='Optional full path $DESI_SPECTRO_REDUX.')
    parser.add_argument('--specprod', type=str, default=None, help="""Optional override of the on-disk spectroscopic
        production directory name under --redux_dir, when it differs from the SPECPROD dependency recorded in the
        Redrock/coadd file headers (e.g., a relocated or "mini" production tree).""")
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
    if args.templatedir is None:
        envlist.append('FTEMPLATES_DIR')
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

    init_argdict = {'fphotofile': args.fphotofile, 'gridfile': gridfile,
                    'templatedir': args.templatedir, 'mapdir': args.mapdir}

    t0 = time.time()
    _initialize_fastbayes_workers(**init_argdict)
    log.info(fsftime('bg_data_init', time.time() - t0))

    _own_pool = False
    if mp_pool is None:
        mp_pool = MPPool(args.mp, initializer=_initialize_fastbayes_workers, init_argdict=init_argdict)
        _own_pool = True

    phot = sc_data.photometry

    Spec = DESISpectra(phot=phot, cosmo=_cosmo, fphotodir=args.fphotodir,
                       redux_dir=args.redux_dir)

    Spec.gather_metadata(args.redrockfiles, firsttarget=args.firsttarget,
                         targetids=targetids, input_redshifts=input_redshifts,
                         ntargets=args.ntargets, zmin=args.zmin,
                         redrockfile_prefix=args.redrockfile_prefix,
                         use_quasarnet=args.use_quasarnet, specprod=args.specprod)
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


def parse_qa(options=None):
    """Parse input arguments to the ``fastbayes-qa`` script.

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

    parser.add_argument('fastbayesfile', help='Full path to an existing FASTBAYES output FITS file.')
    parser.add_argument('--gridfile', type=str, required=True,
                        help='Full path to the Bayesian photometry grid file used to produce '
                        '--fastbayesfile (not read back from its header, since that path may not '
                        'resolve on this machine, e.g. if the fit was run elsewhere -- pass the '
                        'same --gridfile used at fit time).')
    parser.add_argument('--fphotofile', type=str, default=None,
                        help='Photometric configuration file used to produce --fastbayesfile '
                        '(default: bundled fastspecfit configuration; pass the same --fphotofile '
                        'used at fit time if a custom one was used).')
    parser.add_argument('--qadir', type=str, default='.', help='Output directory for QA figures.')
    parser.add_argument('--coadd-type', dest='coadd_type', type=str, default='healpix',
                        choices=['healpix', 'uniqpix', 'cumulative', 'pernight', 'perexp', 'custom', 'stacked'],
                        help='Coadd type, used to build the QA target label/filename (not stored in the '
                        'FASTBAYES output file, so it cannot be inferred).')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of objects to process.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to process, zero-indexed.')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    log.info('fastbayes-qa {}'.format(' '.join(options)))

    return parser.parse_args(options)


def fastbayes_qa(args=None, mp_pool=None):
    """Regenerate QA figures for a FASTBAYES output file, without repeating the fit.

    Reads back an already-written FASTBAYES output file's ``METADATA``,
    ``FASTBAYES``, ``WAVE``, and ``MODELS`` extensions, then calls
    :func:`fastbayes_qa_one` for each requested object. The only work redone
    is the cheap, vectorized :func:`_solve_grid` solve needed to rebuild the
    per-template posterior weights (the one thing not stored in the output
    file); no raw DESI spectra or individual raw template rows are read, so
    this does not need ``$FTEMPLATES_DIR`` (or the original redrock/spectra
    files) to be available.

    ``--gridfile``/``--fphotofile`` must be passed explicitly (matching
    whatever was used to build the grid and run the fit) rather than being
    read back from the output file's header: the header's ``GRIDFILE`` is an
    absolute path that may not exist on this machine (e.g. the fit ran on a
    different machine than this QA regeneration), and ``FPHOTO`` deliberately
    stores only a basename for the same reason (see commit ``ad3df20``) --
    neither is reliably reloadable without user input.

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
    from astropy.table import Table

    if isinstance(args, (list, tuple, type(None))):
        args = parse_qa(args)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    fastbayesfile = os.path.expandvars(args.fastbayesfile)
    if not os.path.isfile(fastbayesfile):
        errmsg = f'FASTBAYES output file {fastbayesfile} not found.'
        log.critical(errmsg)
        raise IOError(errmsg)

    gridfile = os.path.expandvars(args.gridfile)
    if not os.path.isfile(gridfile):
        errmsg = (f'Bayesian grid file {gridfile} not found. Pass the same --gridfile that was '
                  'used to produce this FASTBAYES output file.')
        log.critical(errmsg)
        raise IOError(errmsg)

    fphotofile = args.fphotofile

    F = fitsio.FITS(fastbayesfile)
    meta = Table(F['METADATA'].read())
    results = F['FASTBAYES'].read()
    restwave = F['WAVE'].read()
    modelspectra = F['MODELS'].read()

    nobj = len(meta)
    keep = np.arange(nobj)

    if args.targetids is not None:
        targetids = [int(x) for x in args.targetids.split(',')]
        keep = np.where(np.isin(meta['TARGETID'], targetids))[0]
        if len(keep) == 0:
            log.warning('No objects match the requested --targetids.')
            return 0
    else:
        firsttarget = args.firsttarget
        lasttarget = firsttarget + args.ntargets if args.ntargets is not None else nobj
        keep = keep[firsttarget:lasttarget]

    if len(keep) == 0:
        log.warning('No objects to process.')
        return 0

    # Probe the Legacy Survey viewer once for this call and seed the
    # (per-process) cutout-reachability flag, both here in the main process
    # (covers the args.mp<=1 case) and via the pool initializer below (covers
    # each multiprocessing worker) -- same mechanism as fastqa().
    from fastspecfit.qa import _probe_cutout_host
    cutout_unreachable = not _probe_cutout_host()

    init_argdict = {'fphotofile': fphotofile, 'gridfile': gridfile, 'require_templates': False,
                    'cutout_unreachable': cutout_unreachable}

    t0 = time.time()
    _initialize_fastbayes_workers(**init_argdict)
    log.info(fsftime('bg_data_init', time.time() - t0))

    _own_pool = False
    if mp_pool is None:
        mp_pool = MPPool(args.mp, initializer=_initialize_fastbayes_workers, init_argdict=init_argdict)
        _own_pool = True

    qaargs = [{
        'iobj': iobj,
        'meta': meta[iobj],
        'result': results[iobj],
        'restwave': restwave,
        'restflux': modelspectra[iobj],
        'qadir': args.qadir,
        'coadd_type': args.coadd_type,
    } for iobj in keep]

    t0 = time.time()
    # Forces evaluation: MPPool.starmap returns a lazy itertools.starmap
    # generator when running serially (nworkers=1), which would otherwise
    # never actually execute fastbayes_qa_one (same idiom as qa.desiqa_one).
    for _ in mp_pool.starmap(fastbayes_qa_one, qaargs):
        pass

    if _own_pool:
        mp_pool.close()

    log.info(fsftime('fastbayes_qa_all', time.time() - t0, context=f'nobj={len(qaargs)}'))

    return 0
