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
a local parabola fit along each grid axis (whose names, order, and log-ness
are read per-grid-file from the templates' AXES FITS extension, see
:meth:`BayesianGrid.load`) and an N-linear interpolation over the
neighboring templates, so that the reported parameters, CHI2, and model
spectrum are all mutually consistent.

See ``doc/technote/fastbayes.tex`` for the full grid design, statistical
methodology, sub-grid refinement algorithm, and derived-quantity
definitions.

"""
import os, sys, time, logging
from collections import OrderedDict
import numpy as np
import fitsio
from astropy.table import Table

from fastspecfit.logger import log
from fastspecfit.util import MPPool, fsftime, ZWarningMask, C_LIGHT, FLUXNORM
from fastspecfit.photometry import Photometry
from fastspecfit.singlecopy import sc_data

# The grid's free axes -- names, nested-loop order, log10-refinement flags,
# output column names, and QA labels -- are no longer hard-coded here: they
# are read per-grid-file from the AXES FITS extension (written by
# bin/build-bayesian-templates from its grid-parameter YAML) into instance
# attributes on BayesianGrid (axis_columns/log_axes/axis_outname/
# outname_to_axis/axis_labels/param_names; see BayesianGrid.load below) --
# this is what lets different grid files use entirely different axis
# combinations without any code change here. Axis order is what lets the
# flattened grid's N-D structure be recovered via
# np.unravel_index/np.ravel_multi_index. A fixed (N=1) grid axis is already
# handled correctly (no refinement, weight fully on the single point) by the
# generic per-axis edge case in _refine_grid_axes/_corner_weights below, so
# no special-casing is needed to support one.

# Derived (non-axis) output parameters, computed per-object from the
# closed-form amplitude solve rather than read off a grid axis.
#
# LOGMSTAR is the *surviving* stellar mass (matching fastspec/fastphot's
# own LOGMSTAR convention -- current mass, not total mass ever formed):
# the grid's per-template `mstar` column is FSPS's surviving-mass fraction
# per solar mass formed (mass loss/return already applied), so the fitted
# mass amplitude (Msun formed) times bg_data.mstar gives each template's
# surviving mass. (The total mass *formed*, log10(amplitude) alone, has no
# output column -- it's used internally to scale LOGL_*/LOGLNU_*, which
# are extensive quantities tied to the templates' `_permass` normalization,
# not to the surviving-mass fraction.)
#
# SFR/SSFR (unlike LOGMSTAR) are reported in linear space, matching
# fastspecfit.continuum's SFR/SFR_IVAR convention: bg_data.sfr is exactly
# zero for every template older than the ~100 Myr window baked into the
# templates' precomputed sfr column (passively-evolving population), and a
# log transform cannot represent that exactly, forcing an arbitrary floor
# that collapses every quiescent template onto one identical value and
# produces a spuriously tiny (or exactly zero) log-space posterior
# variance -- reporting SFR/SSFR linearly instead lets 0 be an exact,
# well-behaved output with a meaningful (possibly zero) formal uncertainty.
# SSFR = SFR/Mstar is exactly bg_data.ssfr = bg_data.sfr/bg_data.mstar
# (both per Msun formed): the fitted mass amplitude cancels since SFR and
# Mstar both scale linearly with it, so SSFR is amplitude-independent like
# a grid axis, and its posterior/uncertainty come from bg_data.ssfr's own
# discrete posterior rather than from (incorrectly) differencing the
# separately-refined SFR and LOGMSTAR.
#
# T50 (half-mass assembly time) and MSTARAGE (mass-weighted age, lookback
# convention) are likewise amplitude-independent, purely geometric
# functions of each template's Dense Basis age/t33/t67 -- see
# BayesianGrid.load's t50/mstarage precomputation.
DERIVED_PARAM_NAMES = ('LOGMSTAR', 'SFR', 'SSFR', 'T50', 'MSTARAGE')

# Human-friendly labels for the derived (non-axis) parameters; per-grid
# axis labels come from BayesianGrid.axis_labels instead (see above).
_DERIVED_PARAM_LABELS = {
    'LOGMSTAR': r'$\log_{10}(M_\star/M_{\odot})$',
    'SFR': r'${\rm SFR}\ [M_{\odot}/{\rm yr}]$',
    'SSFR': r'${\rm sSFR}\ [{\rm Gyr}^{-1}]$',
    'T50': r'$t_{50}\ [{\rm Gyr}]$',
    'MSTARAGE': r'$\langle{\rm Age}\rangle_M\ [{\rm Gyr}]$',
}

# Preferred-order candidates for the QA figure's compact text summary
# (fastbayes_qa_one): (output param name, LaTeX template with one {}
# placeholder for the formatted value+error, value format spec). Entries
# whose param name isn't present in this grid's bg_data.param_names (e.g. a
# custom grid without a 'fagn'/'tau' axis) are simply skipped, rather than
# hard-requiring this exact axis set to exist. LOGZZSUN/LOGAGE are displayed
# here (and in the posterior panel below) as linear Z/Zsun and Age -- see
# _QA_LOG_TO_LINEAR -- even though the stored/output columns stay in log10
# space. All entries share the same fixed-decimal format spec so that a
# value and its error always print at the same number of decimal places
# (a per-entry '{:.3g}' spec instead made SFR's error round to a wildly
# different decimal place than its value, e.g. "1.65+/-0.000866").
_QA_SUMMARY_CANDIDATES = (
    ('LOGZZSUN', r'$Z/Z_{{\odot}}={}$', '{:.2f}'),
    ('LOGAGE', r'${{\rm Age}}={}\ {{\rm Gyr}}$', '{:.2f}'),
    ('LOGMSTAR', r'$\log_{{10}}(M_\star/M_{{\odot}})={}$', '{:.2f}'),
    ('SFR', r'$\mathrm{{SFR}}={}\ M_{{\odot}}/\mathrm{{yr}}$', '{:.2f}'),
    ('TAUV', r'$\tau_V={}$', '{:.2f}'),
    ('FAGN', r'$f_{{\rm AGN}}={}$', '{:.2f}'),
)

# Grid-axis output columns whose stored value/error are in log10 space but
# are shown in the QA text summary and posterior panel as their linear
# equivalent (Age in Gyr, Z/Zsun) -- see _log_to_linear for the (division-
# free) error propagation.
_QA_LOG_TO_LINEAR = frozenset(('LOGAGE', 'LOGZZSUN'))
_QA_LINEAR_LABELS = {
    'LOGAGE': r'${\rm Age}\ [{\rm Gyr}]$',
    'LOGZZSUN': r'$Z/Z_{\odot}$',
}

# Grid axes that are still legitimate FASTBAYES output columns (refined and
# reported like any other axis) but excluded from the QA posterior-panel
# grid: t33frac/t67frac are internal Dense Basis SFH bookkeeping, not
# physically meaningful on their own to someone reading the QA figure --
# T50/MSTARAGE (see DERIVED_PARAM_NAMES) are the intended human-readable
# substitutes.
_QA_HIDDEN_PARAMS = frozenset(('T33FRAC', 'T67FRAC'))

# Above this many distinct posterior values, a QA posterior panel bins the
# values into a histogram rather than drawing one bar per value; native
# grid axes never exceed this (the largest in the bundled default grid is
# TAUV's 10), so a larger unique count means an amplitude-scaled quantity
# (LOGMSTAR, SFR, SSFR) or a quantity derived from a combination of axes
# (T50, MSTARAGE, which only depend on the 3 SFH axes but still range over
# up to 8*5*5=200 distinct values).
_DISCRETE_UNIQ_MAX = 30

# Natural minimum posterior-histogram bin width [Gyr] for QA panels of
# these output parameters (see _posterior_binwidth); LOGMSTAR's analogous
# 0.05 dex width is handled generically for any log-refined quantity.
_GYR_BINWIDTH_PARAMS = frozenset(('T50', 'MSTARAGE'))

# Rest-frame luminosity output keys and reference wavelengths (Angstrom).
# The first six match the LOGL_*/LOGLNU_* columns in fastspecfit's own
# SPECPHOT schema (see fastspecfit.continuum.ContinuumTools.lums_keys); the
# last three (3.4, 12, and 22 micron, roughly WISE W1/W3/W4) are fastbayes-
# specific -- monochromatic IR luminosities usable as SFR indicators given
# adequate IR photometric coverage -- and have no main-pipeline counterpart
# since the fastspec/fastphot SPS templates have no dust re-emission.
LUM_KEYS = ('LOGL_1450', 'LOGLNU_1500', 'LOGL_1700', 'LOGLNU_2800', 'LOGL_3000', 'LOGL_5100',
           'LOGL_3MU', 'LOGL_12MU', 'LOGL_22MU')
LUM_WAVES = (1450., 1500., 1700., 2800., 3000., 5100., 34000., 120000., 220000.)

_PC_CM = 3.0856775814913673e18 # [cm]
_LSUN = 3.846e33 # [erg/s]

# Default threshold (see BayesianGrid._ensure_flux_cache) below which the raw
# templates' full FLUX array is preloaded once per worker rather than read
# one row at a time on demand: the current default grid (3600 templates) is
# ~375 MB, comparable to the MAGGIES arrays that are always fully preloaded,
# so this is generous headroom for the default case while still falling back
# to lazy per-row reads for much larger custom grids.
DEFAULT_FLUX_CACHE_MB = 2048.

# Bound on the number of individual raw template rows cached (in the lazy,
# above-threshold path) per worker process, so a custom grid too large to
# preload outright still can't grow this cache without bound.
_FLUX_LRU_MAXSIZE = 5000


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
        self.flux_cache_mb = DEFAULT_FLUX_CACHE_MB
        self._flux_cache = None # full (ntemplate, npix) FLUX array, if small enough to preload
        self._flux_lru = None   # bounded per-row cache over lazy fitsio reads, otherwise


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
        return os.path.join(templatedir, f'bayesian-templates-{imf}-{gridnumber}.fits')


    def load(self, gridfile, templatedir=None, templatesfile=None,
            flux_cache_mb=DEFAULT_FLUX_CACHE_MB):
        """Load a ``bayesian-photometry-*.fits`` grid file (idempotent).

        Parameters
        ----------
        templatesfile : :class:`str` or None, optional
            Full path to the raw Bayesian templates FITS file, overriding
            the ``GRIDNUM``/``IMF``/``templatedir`` reconstruction. Not
            read back from ``gridfile``'s header -- that path may not
            resolve on this machine (e.g. switching between NERSC and a
            laptop) -- so pass it explicitly here (or via ``--templatesfile``)
            when the templates live somewhere non-standard.
        flux_cache_mb : :class:`float`, optional
            Threshold (see :meth:`_ensure_flux_cache`) below which the raw
            templates' full FLUX array is preloaded once per worker rather
            than read one row at a time on demand.

        Raises
        ------
        KeyError
            If ``gridfile`` predates the ``FPHOTO_FILE`` dependency
            keyword (i.e. was built with an older
            ``bin/build-bayesian-photometry``) -- regenerate the grid
            file rather than relying on a fallback.

        """
        from desiutil.depend import getdep

        self.flux_cache_mb = flux_cache_mb

        if self.file == gridfile:
            return

        T = fitsio.FITS(gridfile)
        prihdr = T[0].read_header()

        self.gridnumber = prihdr.get('GRIDNUM')
        self.imf = prihdr.get('IMF')
        self.fphotofile = getdep(prihdr, 'FPHOTO_FILE')
        if templatesfile is not None:
            self.templates_file = os.path.expandvars(templatesfile)
        else:
            self.templates_file = self.get_templates_filename(self.gridnumber, self.imf, templatedir=templatedir)
        self.logspace = bool(prihdr.get('LOGSPACE'))

        photsys_hdr = prihdr.get('PHOTSYS')
        self.photsys_keys = [''] if photsys_hdr in (None, 'NONE') else photsys_hdr.split(',')

        self.redshift = T['REDSHIFT'].read().astype('f8')
        self.zgrid_interp = np.log10(1. + self.redshift) if self.logspace else self.redshift

        self.meta = Table(T['METADATA'].read())
        self.ntemplate = len(self.meta)

        # Each grid template is a single continuity-SFH age bin (sfh=3), so
        # unlike the old single-burst-SSP (sfh=0) design, FSPS's "current"
        # SFR is physically meaningful here -- but bin/build-bayesian-templates
        # doesn't use FSPS's sp.sfr directly (that would normalize by the
        # bin's own width, not by a fixed averaging window); instead it
        # precomputes this column via the same trailing-100-Myr windowed-
        # overlap formula used for the main (tabular-SFH) templates in
        # fastspecfit.continuum._get_sps_properties, specialized to a single
        # age-bin template. Just read it back directly.
        self.sfr = np.asarray(self.meta['sfr'], dtype='f8') # [ntemplate], Msun/yr per Msun formed

        # Surviving stellar mass fraction per solar mass formed (FSPS's
        # sp.stellar_mass, mass loss/return already applied) -- scaled by
        # the fitted mass amplitude in fastbayes_one/fastbayes_qa_one to get
        # each template's actual (surviving) stellar mass estimate, matching
        # fastspec/fastphot's own LOGMSTAR convention. SSFR = SFR/Mstar is
        # precomputed here too since it's likewise amplitude-independent
        # (both scale linearly with the fitted amplitude, which cancels).
        self.mstar = np.asarray(self.meta['mstar'], dtype='f8') # [ntemplate], dimensionless
        self.ssfr = self.sfr / self.mstar # [ntemplate], 1/yr

        # T50 (half-mass assembly time) and MSTARAGE (mass-weighted age,
        # lookback convention) are purely geometric functions of each
        # template's Dense Basis tertile boundaries (t33, t67) and total age
        # -- see bin/build-bayesian-templates's module docstring for the SFH
        # construction. T50 = (t33+t67)/2 because SFR is constant within the
        # middle tertile segment, so cumulative mass grows linearly there and
        # the 50%-mass point is exactly that segment's midpoint. MSTARAGE
        # averages the three equal-mass (1/3 each) segments' midpoint
        # lookback times: age - mean([t33/2, (t33+t67)/2, (t67+age)/2])
        # = (5*age - 2*t33 - 2*t67) / 6.
        t33_grid = np.asarray(self.meta['t33'], dtype='f8') # [Gyr]
        t67_grid = np.asarray(self.meta['t67'], dtype='f8') # [Gyr]
        age_grid = np.asarray(self.meta['age'], dtype='f8')  # [Gyr]
        self.t50 = (t33_grid + t67_grid) / 2. # [Gyr]
        self.mstarage = (5. * age_grid - 2. * t33_grid - 2. * t67_grid) / 6. # [Gyr]

        # Rest-frame Dn(4000) and monochromatic luminosities-per-solar-mass-
        # formed, precomputed once per template by bin/build-bayesian-templates
        # (they depend only on the template's intrinsic rest-frame spectral
        # shape, not on redshift or observed photometry) -- cached here so
        # per-object LOGL_*/DN4000_MODEL uncertainties (fastbayes_one) are a
        # cheap vectorized weighted-variance calculation, with no extra raw
        # spectrum reads.
        self.dn4000_model_permass = np.asarray(self.meta['dn4000_model'], dtype='f8')
        self.lum_permass = {key: np.asarray(self.meta[f'{key.lower()}_permass'], dtype='f8') for key in LUM_KEYS}

        self._axis_cache = {}
        self._flux_cache = None
        self._flux_lru = None

        # Per-grid-file axis list/order/log-ness/output names/labels, read
        # from the AXES extension (written by bin/build-bayesian-templates
        # from its grid-parameter YAML) rather than hard-coded here -- this
        # is what lets different grid files use entirely different axis
        # combinations without any fastbayes.py code change.
        axes_tbl = T['AXES'].read()
        self.axis_columns = tuple(n.strip() for n in axes_tbl['name'])
        self.log_axes = frozenset(n for n, islog in zip(self.axis_columns, axes_tbl['log']) if bool(islog))
        self.axis_outname = {n: out.strip() for n, out in zip(self.axis_columns, axes_tbl['outname'])}
        self.outname_to_axis = {v: k for k, v in self.axis_outname.items()}
        self.axis_labels = {n: lbl.strip() for n, lbl in zip(self.axis_columns, axes_tbl['label'])}
        self.param_names = tuple(self.axis_outname[col] for col in self.axis_columns) + DERIVED_PARAM_NAMES
        self.param_labels = {self.axis_outname[col]: self.axis_labels[col] for col in self.axis_columns}
        self.param_labels.update(_DERIVED_PARAM_LABELS)

        # Recover the grid's N-D axis structure. The grid is a full
        # factorial/Cartesian-product design (every combination of the
        # axes exists), so the per-axis grid values and the flattened
        # index<->multi-index mapping are recoverable directly from the
        # metadata table -- no separate shape bookkeeping needs to be
        # stored in the FITS file.
        self.axis_values = {col: np.unique(np.asarray(self.meta[col], dtype='f8')) for col in self.axis_columns}
        self.dims = tuple(len(self.axis_values[col]) for col in self.axis_columns)
        if int(np.prod(self.dims)) != self.ntemplate:
            errmsg = (f'Grid file {gridfile} metadata is not a full factorial grid over '
                      f'{self.axis_columns} (dims product {int(np.prod(self.dims))} != ntemplate {self.ntemplate}).')
            log.critical(errmsg)
            raise ValueError(errmsg)

        # Output names of single-grid-point (N=1) axes, e.g. UMIN/GAMMA in
        # the bundled default grid: there is exactly one possible value, so
        # `vals` (the per-template axis value fed into the weighted
        # MEAN/MODE/percentile/ERR computation in fastbayes_one/
        # fastbayes_qa_one) is bit-identical across every template and the
        # true weighted variance is exactly zero -- computing it anyway
        # returns floating-point roundoff noise (e.g. ~1e-18) instead, which
        # both callers special-case away using this set.
        self.fixed_outnames = frozenset(
            self.axis_outname[col] for col, n in zip(self.axis_columns, self.dims) if n == 1)

        filt = T['FILTERS'].read()
        self.bands = np.array([b.strip() for b in filt['band']])
        self.nband = len(self.bands)

        self.maggies = {}
        for key in self.photsys_keys:
            extname = 'MAGGIES' + ('' if key == '' else f'_{key}')
            self.maggies[key] = T[extname].read() # [ntemplate, nz, nband]

        # Rest-frame (band-shifted) photometry, tabulated the same way as
        # MAGGIES above by bin/build-bayesian-photometry (see
        # synthesize_photometry_grid), and used by fastbayes_one to get
        # ABSMAG*_SYNTH_ERR_*/KCORR*_ERR_* for free (a vectorized weighted-
        # variance over the whole grid, no per-object filter convolution).
        # Optional: grid files predating this feature simply lack the
        # extension, in which case those output columns are omitted rather
        # than erroring -- no fallback/rebuild is forced on older grids.
        extnames = [hdu.get_extname() for hdu in T]
        self.has_restmaggies = 'RESTMAGGIES' in extnames
        if self.has_restmaggies:
            self.restmaggies = T['RESTMAGGIES'].read() # [ntemplate, nz, nabs]
            restfilt = T['RESTFILTERS'].read()
            self.absmag_bands_grid = np.array([b.strip() for b in restfilt['band']])
        else:
            self.restmaggies = None
            self.absmag_bands_grid = None
            log.debug(f'Grid file {gridfile} has no RESTMAGGIES extension; '
                     'ABSMAG*_SYNTH_ERR_*/KCORR*_ERR_* will be omitted.')

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
                errmsg = (f'Could not resolve a raw templates file for grid {self.file} '
                          '(missing GRIDNUM/IMF header keyword(s), or $FTEMPLATES_DIR is unset; '
                          'check --templatedir or pass --templatesfile explicitly).')
                log.critical(errmsg)
                raise ValueError(errmsg)
            self._templates_fits = fitsio.FITS(self.templates_file)
        return self._templates_fits


    def template_wave(self):
        """Native rest-frame wavelength array shared by every template."""
        if self._wave is None:
            self._wave = self.templates_fits()['WAVE'].read()
        return self._wave


    def _ensure_flux_cache(self):
        """Decide, on first raw-template access, whether to preload the
        entire FLUX array or fall back to a bounded per-row LRU cache over
        lazy ``fitsio`` reads.

        The decision is made lazily (not in :meth:`load`) because the raw
        templates file is not always required up front -- e.g.
        :func:`fastbayes_qa` with ``--ndraw 0`` never touches it -- so
        opening it here, on first actual access, keeps that case working
        without a raw templates file at all.

        """
        if self._flux_cache is not None or self._flux_lru is not None:
            return

        npix = len(self.template_wave())
        nbytes = self.ntemplate * npix * 4 # FLUX is stored as f4

        if nbytes <= self.flux_cache_mb * 1e6:
            self._flux_cache = self.templates_fits()['FLUX'].read()
            log.debug(f'Preloaded the full {nbytes/1e6:.0f} MB FLUX array for {self.file}.')
        else:
            self._flux_lru = OrderedDict()
            log.debug(f'FLUX array ({nbytes/1e6:.0f} MB) exceeds --flux-cache-mb='
                     f'{self.flux_cache_mb:.0f}; using a bounded per-row LRU cache instead.')


    def template_flux_row(self, itemplate):
        """Read one raw rest-frame template spectrum by flat grid index.

        Served from the fully preloaded FLUX array, or from a bounded
        per-worker LRU cache over lazy per-row ``fitsio`` reads, depending
        on which mode :meth:`_ensure_flux_cache` picked for this grid; either
        way, repeated corner hits across many objects within one worker
        never re-read the same row from disk more than necessary.

        """
        self._ensure_flux_cache()
        itemplate = int(itemplate)

        if self._flux_cache is not None:
            return self._flux_cache[itemplate]

        row = self._flux_lru.get(itemplate)
        if row is not None:
            self._flux_lru.move_to_end(itemplate)
            return row

        row = self.templates_fits()['FLUX'][itemplate, :].ravel()
        self._flux_lru[itemplate] = row
        if len(self._flux_lru) > _FLUX_LRU_MAXSIZE:
            self._flux_lru.popitem(last=False)
        return row


    def _z_interp_index(self, redshift):
        """Shared linear-interpolation index/fraction into the grid's
        redshift axis, used by both :meth:`interpolate_at_z` and
        :meth:`interpolate_restmaggies_at_z`.

        Returns
        -------
        idx : :class:`int`
        frac : :class:`float`

        """
        zvar = np.log10(1. + redshift) if self.logspace else redshift
        zvar = np.clip(zvar, self.zgrid_interp[0], self.zgrid_interp[-1])

        idx = int(np.clip(np.searchsorted(self.zgrid_interp, zvar), 1, len(self.zgrid_interp) - 1))
        z0, z1 = self.zgrid_interp[idx - 1], self.zgrid_interp[idx]
        frac = 0. if z1 == z0 else (zvar - z0) / (z1 - z0)

        return idx, frac


    def interpolate_at_z(self, photsys, redshift):
        """Interpolate ``MAGGIES(template, z, band)`` to a single redshift.

        The raw tabulated ``MAGGIES`` include the luminosity-distance
        dimming factor ``(10/(1e6*dlum(z)))**2/(1+z)``, which diverges as
        ``~1/z**2`` toward ``z=0`` -- far too steep to linearly interpolate
        across the grid's z nodes without large error close to ``zmin``
        (e.g. for a nearby-galaxy grid with ``zmin`` pushed well below the
        grid's node spacing). So the divergent factor is divided out of
        each bracketing node before interpolating (leaving only the smooth
        IGM/band-shift shape, safe to interpolate) and the *exact* factor
        for ``redshift`` -- not a grid-interpolated one -- is multiplied
        back in via :func:`_distfactor`.

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
        idx, frac = self._z_interp_index(redshift)
        z0, z1 = self.redshift[idx - 1], self.redshift[idx]

        shape0 = maggies[:, idx - 1, :] / _distfactor(z0)
        shape1 = maggies[:, idx, :] / _distfactor(z1)

        return (shape0 + frac * (shape1 - shape0)) * _distfactor(redshift)


    def interpolate_restmaggies_at_z(self, redshift):
        """Interpolate ``RESTMAGGIES(template, z, absmag_band)`` to a single
        redshift; mirrors :meth:`interpolate_at_z` (including the exact,
        non-interpolated re-application of the divergent distance-dimming
        factor via :func:`_distfactor`) but has no photsys split, since
        ``Photometry.filters_out`` (the rest-frame band-shifted output
        filters) does not depend on photsys.

        Parameters
        ----------
        redshift : :class:`float`
            Object redshift.

        Returns
        -------
        :class:`numpy.ndarray`
            Model rest-frame maggies per solar mass formed, shape
            (ntemplate, nabs).

        """
        idx, frac = self._z_interp_index(redshift)
        z0, z1 = self.redshift[idx - 1], self.redshift[idx]

        shape0 = self.restmaggies[:, idx - 1, :] / _distfactor(z0, rest=True)
        shape1 = self.restmaggies[:, idx, :] / _distfactor(z1, rest=True)

        return (shape0 + frac * (shape1 - shape0)) * _distfactor(redshift, rest=True)


    def axis_posterior_cache(self, col):
        """Lazily cache per-grid-axis posterior statistics for grid-axis column
        ``col``. These depend only on ``self.meta`` (identical for every
        object fit against this grid), so they are computed once per worker
        rather than once per object.

        Returns
        -------
        vals : :class:`numpy.ndarray`
            Per-template values of ``col`` in its fit coordinate (log10 for
            axes in :attr:`log_axes`, linear otherwise).
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
            if col in self.log_axes:
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
                                  templatesfile=None, require_templates=True, mapdir=None,
                                  cutout_unreachable=None, flux_cache_mb=DEFAULT_FLUX_CACHE_MB):
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
    templatesfile : :class:`str` or None, optional
        Full path to the raw Bayesian templates FITS file, overriding
        ``templatedir`` reconstruction. See :meth:`BayesianGrid.load`.
    mapdir : :class:`str` or None, optional
        Directory containing the Milky Way dust maps; defaults to
        ``$DUST_DIR/maps`` when ``None``. Not needed by
        :func:`fastbayes_qa`/:func:`fastbayes_qa_one`, which never call
        :func:`fastspecfit.io.one_spectrum`.
    require_templates : :class:`bool`, optional
        If ``True`` (default), require the raw templates FITS file (used by
        :meth:`BayesianGrid.template_flux_row` to build the refined
        rest-frame spectrum) to be present. Both :func:`fastbayes` (during
        fitting) and :func:`fastbayes_qa` (regenerating the spectrum via
        :func:`_build_refined_spectrum`, since it is not stored in the
        output file) always need it, so neither passes ``False`` here in
        practice; this remains available for callers that only need, e.g.,
        the grid's tabulated per-template quantities (``sfr``,
        ``dn4000_model_permass``, ``lum_permass``) without ever touching a
        raw spectrum.
    cutout_unreachable : :class:`bool` or None, optional
        Seeds this worker's copy of :data:`fastspecfit.qa._cutout_unreachable`
        (same sticky per-process mechanism ``fastqa`` uses to detect an
        unreachable Legacy Survey viewer once, up front, rather than letting
        every object separately discover the outage via
        :func:`fastspecfit.qa._fetch_cutout`'s retry/backoff loop). Left
        untouched (``None``, the default) by the plain fitting path
        (:func:`fastbayes`), which never generates QA figures and so never
        imports ``fastspecfit.qa`` at all; only :func:`fastbayes_qa` passes
        this.
    flux_cache_mb : :class:`float`, optional
        Threshold passed straight through to :meth:`BayesianGrid.load`; see
        there for details.

    """
    global _igm, _cosmo

    sc_data.photometry = Photometry(fphotofile=fphotofile)
    sc_data.set_mapdir(mapdir)
    bg_data.load(gridfile, templatedir=templatedir, templatesfile=templatesfile,
                flux_cache_mb=flux_cache_mb)

    if cutout_unreachable is not None:
        import fastspecfit.qa as qa_module
        qa_module._cutout_unreachable = cutout_unreachable

    if require_templates and not os.path.isfile(bg_data.templates_file):
        errmsg = (f'Bayesian templates file {bg_data.templates_file} not found; '
                  'check $FTEMPLATES_DIR, --templatedir, or --templatesfile.')
        log.critical(errmsg)
        raise IOError(errmsg)

    if list(bg_data.bands) != list(sc_data.photometry.bands):
        errmsg = (f"Band mismatch between the photometry configuration {sc_data.photometry.fphotofile} "
                  f"and the Bayesian grid file {gridfile}, which was built with fphoto file "
                  f"'{bg_data.fphotofile}'; they must be built from the same (or a compatible) "
                  "fphoto configuration file.")
        log.critical(errmsg)
        raise ValueError(errmsg)

    if bg_data.has_restmaggies and list(bg_data.absmag_bands_grid) != list(sc_data.photometry.absmag_bands):
        errmsg = (f"absmag_bands mismatch between the photometry configuration {sc_data.photometry.fphotofile} "
                  f"and the Bayesian grid file {gridfile}'s RESTMAGGIES/RESTFILTERS extensions; "
                  "they must be built from the same (or a compatible) fphoto configuration file.")
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
    (log10 for axes in ``bg_data.log_axes``, linear otherwise)."""
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


def _distfactor(redshift, rest=False):
    """Exact (non-interpolated) luminosity-distance dimming factor, matching
    the grid-build convention in ``bin/build-bayesian-photometry``'s
    ``synthesize_photometry_grid`` (``distfactor``/``rest_zfactor``'s
    ``dlum``-dependent term).

    """
    dlum = _cosmo.luminosity_distance(redshift) # [Mpc]
    df = (10. / (1e6 * dlum))**2
    return df if rest else df / (1. + redshift)


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
        ``ibest`` (flat index of the discrete chi2 minimum), ``fit_value``
        (per grid-axis-column refined value, from :func:`_refine_grid_axes`),
        ``corner_idx``/``corner_weight`` (N-linear interpolation corners,
        from :func:`_corner_weights`), ``refined_model_maggies`` (shape
        (nband,)), ``refined_amplitude``, ``chi2_refined``, and
        ``model_maggies`` (shape (ntemplate, nband) -- the full per-template,
        z-interpolated observed-band photometry already computed here, reused
        by :func:`fastbayes_one` to get ``KCORR*_ERR_*`` for free from
        ``BayesianGrid.restmaggies`` without any extra filter convolution).

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
    fit_value, frac, frac_dir = _refine_grid_axes(ibest, chi2)
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
        'fit_value': fit_value,
        'corner_idx': corner_idx, 'corner_weight': corner_weight,
        'refined_model_maggies': refined_model_maggies,
        'refined_amplitude': refined_amplitude, 'chi2_refined': chi2_refined,
        'model_maggies': model_maggies,
    }


def _build_refined_spectrum(corner_idx, corner_weight, refined_amplitude):
    """N-linearly interpolate the refined rest-frame spectrum from the
    neighboring raw grid templates (:func:`_solve_grid`'s ``corner_idx``/
    ``corner_weight``), scaled by the refined mass amplitude.

    Factored out so that both :func:`fastbayes_one` (during fitting) and
    :func:`fastbayes_qa_one` (regenerating QA figures from an
    already-written output file, which no longer stores the spectrum
    itself) rebuild it identically, from the same handful of
    :meth:`BayesianGrid.template_flux_row` reads.

    """
    npix = len(bg_data.template_wave())
    restflux = np.zeros(npix, dtype='f8')
    for idx, w in zip(corner_idx, corner_weight):
        restflux += w * bg_data.template_flux_row(idx)
    restflux *= refined_amplitude # [erg/s/cm2/A at 10pc, actual stellar mass]
    return restflux


def _nearest_observed_band(lambda_obs, lambda_out, maggies, ivarmaggies, redshift, snrmin=2.):
    """Nearest-qualifying-observed-band index per rest-frame output band.

    Duplicates (deliberately, rather than importing) the band-selection
    logic inside :meth:`fastspecfit.photometry.Photometry.kcorr_and_absmag`,
    so that :func:`fastbayes_one` can reuse the *same* ``oband`` selection
    to pull each grid template's already-tabulated observed-band photometry
    (``BayesianGrid``'s ``model_maggies``) for the per-template ``KCORR``
    needed by the weighted-posterior-variance estimate -- without changing
    that shared method's signature/behavior (used as-is by
    ``fastspecfit.continuum`` too).

    Parameters
    ----------
    lambda_obs, lambda_out : :class:`numpy.ndarray`
        Effective wavelengths (Angstroms) of the observed and rest-frame
        (band-shifted) output filter sets, respectively.
    maggies, ivarmaggies : :class:`numpy.ndarray`
        Observed photometry (maggies) and its inverse variance, one entry
        per observed band.
    redshift : :class:`float`
    snrmin : :class:`float`, optional
        Minimum S/N in an observed bandpass for it to qualify. Default 2.

    Returns
    -------
    oband : :class:`numpy.ndarray`
        Index into the observed bands, one per output (rest-frame) band.
    good : :class:`numpy.ndarray`
        Boolean mask, ``True`` where ``oband``'s pick actually meets
        ``snrmin`` (mirrors ``kcorr_and_absmag``'s ``I``).

    """
    nabs = len(lambda_out)
    oband = np.empty(nabs, dtype=np.int16)
    for jj in range(nabs):
        lambdadist = np.abs(lambda_obs / (1. + redshift) - lambda_out[jj])
        oband[jj] = np.argmin(lambdadist + (maggies * np.sqrt(ivarmaggies) < snrmin) * 1e10)
    good = (maggies[oband] * np.sqrt(ivarmaggies[oband]) > snrmin)
    return oband, good


def _refine_grid_axes(ibest, chi2):
    """Locally refine each grid axis's ML value via a 3-point parabola fit.

    For each grid axis, holds all other axes fixed at their best-fit
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

    fit_value, frac, frac_dir = {}, {}, {}

    for axis_pos, col in enumerate(bg_data.axis_columns):
        n = bg_data.dims[axis_pos]
        i0 = multi_index[axis_pos]
        grid_vals = bg_data.axis_values[col]
        log_axis = col in bg_data.log_axes

        center_coord = np.log10(grid_vals[i0]) if log_axis else grid_vals[i0]

        if i0 <= 0 or i0 >= n - 1:
            # at a grid edge -- cannot bracket a minimum, no refinement
            fit_value[col] = center_coord
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
            frac[col] = 0.
            frac_dir[col] = 1
            continue

        fit_value[col] = x0

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

    return fit_value, frac, frac_dir


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
        2**naxes; fewer whenever some axes were not refined).
    weights : :class:`numpy.ndarray`
        Corresponding non-negative weights, summing to 1.

    """
    multi_index = np.array(np.unravel_index(ibest, bg_data.dims))
    naxes = len(bg_data.axis_columns)
    t = np.array([frac[col] for col in bg_data.axis_columns])
    direction = np.array([frac_dir[col] for col in bg_data.axis_columns])

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
    cols = [('CHI2', 'f4'), ('NDOF', 'i2')]

    # One block per parameter: refined maximum-likelihood value and its
    # formal uncertainty (the sqrt of the weighted 2nd moment of the
    # marginalized discrete posterior), then descriptive statistics (mean,
    # mode, and the 25th/50th/75th percentiles) of that same posterior.
    for pname in bg_data.param_names:
        cols += [(pname, 'f4'), (f'{pname}_ERR', 'f4'), (f'{pname}_MEAN', 'f4'),
                (f'{pname}_MODE', 'f4'), (f'{pname}_P25', 'f4'), (f'{pname}_P50', 'f4'),
                (f'{pname}_P75', 'f4')]

    # K-corrections and absolute magnitudes, computed from the refined
    # rest-frame spectrum blended with the real observed photometry.
    # ABSMAG_SYNTH is the purely synthetic (model-only, no K-correction
    # blending) absolute magnitude. KCORR_ERR/ABSMAG_SYNTH_ERR (weighted
    # std over the full discrete posterior, from BayesianGrid.restmaggies --
    # see fastbayes_one) are only declared when the loaded grid actually has
    # a RESTMAGGIES extension; older grid files simply omit these columns
    # rather than erroring. ABSMAG_IVAR is fixed up using the same two
    # quantities (see fastbayes_one for exactly how, to avoid double-
    # counting): observational-noise-only when a qualifying observed band
    # was used for K-correction (as before), now combined in quadrature with
    # KCORR_ERR's own model uncertainty; ABSMAG_SYNTH_ERR when there wasn't
    # (absmag falls back to synth_absmag exactly, previously always zero).
    band_shifts = [(band.upper(), int(10 * shift))
                  for band, shift in zip(phot.absmag_bands, phot.band_shift)]
    cols += [(f'KCORR{shift:02d}_{band}', 'f4') for band, shift in band_shifts]
    cols += [(f'ABSMAG{shift:02d}_{band}', 'f4') for band, shift in band_shifts]
    cols += [(f'ABSMAG{shift:02d}_SYNTH_{band}', 'f4') for band, shift in band_shifts]
    cols += [(f'ABSMAG{shift:02d}_IVAR_{band}', 'f4') for band, shift in band_shifts]
    if bg_data.has_restmaggies:
        cols += [(f'KCORR{shift:02d}_ERR_{band}', 'f4') for band, shift in band_shifts]
        cols += [(f'ABSMAG{shift:02d}_SYNTH_ERR_{band}', 'f4') for band, shift in band_shifts]

    # Rest-frame luminosities and the model Dn(4000) index, computed from the
    # refined rest-frame spectrum, plus their formal uncertainty (sqrt of the
    # weighted 2nd moment of their marginalized discrete posterior, computed
    # from per-template quantities tabulated once at grid-build time -- see
    # BayesianGrid.lum_permass/dn4000_model_permass).
    for key in LUM_KEYS:
        cols += [(key, 'f4'), (f'{key}_ERR', 'f4')]
    cols += [('DN4000_MODEL', 'f4'), ('DN4000_MODEL_ERR', 'f4')]

    if topk > 0:
        cols += [('TOPK_INDEX', 'i4', (topk,)), ('TOPK_WEIGHT', 'f4', (topk,))]

    return np.dtype(cols)


def fastbayes_one(iobj, data, meta, fastbayes_dtype, topk=0):
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

    Returns
    -------
    meta : :class:`astropy.table.Row`
        Updated metadata row with observed photometry filled in.
    result : :class:`numpy.ndarray`
        Bayesian-fitting output row.

    """
    from fastspecfit.io import one_spectrum

    phot = sc_data.photometry

    log.info(f'Bayesian fitting object {iobj} [{phot.uniqueid_col.lower()} '
            f'{data["uniqueid"]}, z={data["redshift"]:.6f}].')

    one_spectrum(data, meta, fastphot=True, synthphot=True)

    # Copy parsed photometry from the 'data' dictionary to the 'meta' table
    # (mirroring fastspec_one's convention).
    flux = data['photometry']['nanomaggies']
    fluxivar = data['photometry']['nanomaggies_ivar']
    for iband, band in enumerate(phot.bands):
        meta[f'FLUX_{band.upper()}'] = flux[iband]
        meta[f'FLUX_IVAR_{band.upper()}'] = fluxivar[iband]

    result = np.zeros(1, dtype=fastbayes_dtype)[0]

    redshift = data['redshift']
    photsys = data['photsys']

    if redshift > bg_data.redshift[-1]:
        log.warning(f'Object {iobj} [{phot.uniqueid_col.lower()} {data["uniqueid"]}] redshift '
                   f'{redshift:.6f} exceeds the grid maximum {bg_data.redshift[-1]:.6f}; skipping the fit.')
        return meta, result

    flam = np.asarray(data['photometry']['flam'], dtype='f8')
    flam_ivar = np.asarray(data['photometry']['flam_ivar'], dtype='f8') * phot.bands_to_fit
    lambda_eff = np.asarray(data['photometry']['lambda_eff'], dtype='f8')

    # With <=1 usable band the closed-form amplitude solve in _solve_grid
    # can drive chi2 to *exactly* zero for every template in the grid (with
    # 0 bands because denom==0 everywhere forces amplitude=0; with exactly
    # 1 band because a single free amplitude always exactly zeros a
    # single-point residual, for any template shape) -- np.argmin then
    # silently returns the grid's arbitrary first/lowest-edge template
    # instead of a real fit, poisoning every output column (LOGMSTAR floors
    # to log10(1e-30)=-30, every grid-axis column snaps to its lowest grid
    # value, SSFR/T50/MSTARAGE report that arbitrary template's intrinsic
    # values, and CHI2 comes out exactly 0 as if it were a perfect fit).
    # Skip the fit entirely and leave `result` at its all-zeros default
    # instead, mirroring both the redshift-out-of-grid-range return above
    # and continuum.py's analogous `ndof_phot == 0` branch.
    nbands_used = int(np.sum(flam_ivar > 0.))
    if nbands_used <= 1:
        log.warning(f'Object {iobj} [{phot.uniqueid_col.lower()} {data["uniqueid"]}] has only '
                   f'{nbands_used} usable photometric band(s); skipping the fit.')
        return meta, result

    soln = _solve_grid(flam, flam_ivar, lambda_eff, photsys, redshift)
    chi2 = soln['chi2']
    weight = soln['weight']
    amplitude = soln['amplitude']
    ibest = soln['ibest']
    fit_value = soln['fit_value']
    corner_idx = soln['corner_idx']
    corner_weight = soln['corner_weight']
    refined_model_maggies = soln['refined_model_maggies']
    refined_amplitude = soln['refined_amplitude']
    chi2_refined = soln['chi2_refined']
    model_maggies = soln['model_maggies']

    # log10(Msun formed) -- internal only (no output column): used below to
    # scale LOGL_*/LOGLNU_*, which are extensive quantities tied to the
    # templates' `_permass` (per Msun *formed*) normalization, not to the
    # surviving-mass fraction that LOGMSTAR uses.
    amplitude_log = np.log10(np.clip(amplitude, 1e-30, None))

    # Surviving stellar mass: the grid's per-template `mstar` column is
    # FSPS's surviving-mass fraction per solar mass formed (mass loss/return
    # already applied), so scale by the fitted mass amplitude to get each
    # template's actual surviving-mass estimate -- matching fastspec/
    # fastphot's own LOGMSTAR convention (current mass, not total mass ever
    # formed).
    mstar_per_template = amplitude * bg_data.mstar # [ntemplate], Msun (surviving)
    logmstar = np.log10(np.clip(mstar_per_template, 1e-30, None))

    # The grid stores sfr per solar mass formed (like the flux templates
    # themselves), so scale by the fitted mass amplitude to get each
    # template's actual SFR estimate for this object. Reported linearly
    # (not logged): bg_data.sfr is exactly zero for every template older
    # than the ~100 Myr window baked into the templates' precomputed sfr
    # column, and SFR=0 is a valid, physically meaningful value (passively
    # evolving population) that a log transform cannot represent.
    sfr_per_template = amplitude * bg_data.sfr # [ntemplate], Msun/yr

    refined_mstar_permass = np.sum(corner_weight * bg_data.mstar[corner_idx])
    refined_logmstar = np.log10(np.clip(refined_amplitude * refined_mstar_permass, 1e-30, None))

    refined_sfr_per_msun = np.sum(corner_weight * bg_data.sfr[corner_idx])
    refined_sfr = refined_amplitude * refined_sfr_per_msun

    # SSFR = SFR/Mstar = bg_data.ssfr = bg_data.sfr/bg_data.mstar: the
    # fitted mass amplitude cancels exactly (both SFR and Mstar scale
    # linearly with it), so SSFR is amplitude-independent like a grid axis,
    # and its refined value is the same N-linear interpolation over the
    # corner simplex used for every other intrinsic per-template quantity
    # (not a division of the separately-refined SFR and LOGMSTAR, which
    # would propagate their individual refinement/rounding error
    # incorrectly).
    refined_ssfr = np.sum(corner_weight * bg_data.ssfr[corner_idx]) # [1/yr]

    # T50/MSTARAGE are likewise amplitude-independent (see
    # BayesianGrid.load), so refine them the same way.
    refined_t50 = np.sum(corner_weight * bg_data.t50[corner_idx]) # [Gyr]
    refined_mstarage = np.sum(corner_weight * bg_data.mstarage[corner_idx]) # [Gyr]

    # --- refined rest-frame spectrum (N-linear combination of the
    # neighboring raw template spectra), scaled by the refined mass -------
    restwave = bg_data.template_wave()
    restflux = _build_refined_spectrum(corner_idx, corner_weight, refined_amplitude)

    # Synthetic observed-frame photometry from the interpolated grid itself
    # (refined_model_maggies/refined_amplitude, already computed above by
    # _solve_grid) -- consistent with what CHI2 was actually computed from,
    # unlike a fresh speclite resynthesis of the resampled spectrum. Written
    # into the METADATA row (FLUX_SYNTH_*, pre-declared by the fastbayes()
    # driver), alongside the real FLUX_*/FLUX_IVAR_* columns.
    synth_maggies_grid = refined_model_maggies * refined_amplitude # [nband], observed-frame maggies
    for iband, band in enumerate(phot.bands):
        meta[f'FLUX_SYNTH_{band.upper()}'] = 1e9 * synth_maggies_grid[iband] # [nanomaggies]

    # --- K-corrections, absolute magnitudes, rest-frame luminosities, and
    # the model Dn(4000) index, all derived from the refined rest-frame
    # spectrum -------------------------------------------------------------
    zwave = restwave * (1. + redshift)
    igm_trans = _igm.full_IGM(redshift, zwave)
    dlum = _cosmo.luminosity_distance(redshift) # [Mpc]
    zfactor = igm_trans * (10. / (1e6 * dlum))**2 / (1. + redshift)
    zflux = restflux * zfactor # [erg/s/cm2/A, observed frame]

    dmod = _cosmo.distance_modulus(redshift)
    # Photometry.synth_absmag/kcorr_and_absmag internally divide by FLUXNORM
    # (see fastspecfit.util), a convention that only makes sense because
    # fastspecfit.continuum multiplies its own model flux by FLUXNORM before
    # ever calling them; restflux/zflux here are genuine physical flux
    # (erg/s/cm2/A), so FLUXNORM must be applied at the call site instead, or
    # every ABSMAG*_SYNTH_* (and the ABSMAG*_* fallback when no observed band
    # qualifies) comes out exactly 2.5*log10(FLUXNORM) = 42.5 mag too faint.
    zflux_fluxnorm = zflux * FLUXNORM
    synth_absmag, synth_maggies_rest = phot.synth_absmag(redshift, dmod, zwave, zflux_fluxnorm)
    kcorr, absmag, ivarabsmag, _ = phot.kcorr_and_absmag(
        flux, fluxivar, redshift, dmod, photsys, zwave, zflux_fluxnorm,
        synth_absmag, synth_maggies_rest)

    # --- ABSMAG_SYNTH_ERR/KCORR_ERR from the full discrete posterior, using
    # BayesianGrid.restmaggies (pre-synthesized by bin/build-bayesian-photometry,
    # same as model_maggies above) -- a vectorized weighted-variance over
    # every template, no per-object filter convolution or re-solving of the
    # grid, unlike fastspecfit.continuum's --nmonte. Omitted (no columns in
    # fastbayes_dtype at all) when the loaded grid predates RESTMAGGIES.
    if bg_data.has_restmaggies:
        restmaggies_at_z = bg_data.interpolate_restmaggies_at_z(redshift) # [ntemplate, nabs]

        # KCORR is amplitude-independent (it cancels in the ratio), so its
        # per-template value needs no mass scaling at all.
        oband, good_oband = _nearest_observed_band(
            phot.filters[photsys].effective_wavelengths.value,
            phot.filters_out.effective_wavelengths.value, flux * 1e-9,
            (fluxivar / 1e-9**2) * phot.bands_to_fit, redshift)
        with np.errstate(divide='ignore', invalid='ignore'):
            kcorr_per_template = 2.5 * np.log10(restmaggies_at_z / model_maggies[:, oband]) # [ntemplate, nabs]
        kcorr_mean = np.sum(weight[:, np.newaxis] * kcorr_per_template, axis=0)
        kcorr_err = np.sqrt(np.sum(weight[:, np.newaxis] * (kcorr_per_template - kcorr_mean)**2, axis=0))

        # Floor amplitude (matching the amplitude_log/logmstar convention
        # above) rather than leaving exact-zero-amplitude templates (routine
        # -- the non-negative closed-form solve in _solve_grid clips many
        # poorly-fitting templates to amplitude=0) to produce log10(0)=-inf:
        # left unfloored, weight[i]*(+inf) is 0*inf=nan for any such
        # template whose weight underflows to exactly 0, which corrupts the
        # *entire* weighted sum (np.sum does not skip NaN) rather than
        # contributing the intended ~0.
        amplitude_floor = np.clip(amplitude, 1e-30, None)
        with np.errstate(divide='ignore', invalid='ignore'):
            synth_absmag_per_template = -2.5 * np.log10(restmaggies_at_z * amplitude_floor[:, np.newaxis]) - dmod
        synth_absmag_mean = np.sum(weight[:, np.newaxis] * synth_absmag_per_template, axis=0)
        absmag_synth_err = np.sqrt(np.sum(weight[:, np.newaxis] * (synth_absmag_per_template - synth_absmag_mean)**2, axis=0))

        # Fix ABSMAG_IVAR to actually reflect the total uncertainty in
        # ABSMAG, combining (in quadrature, i.e. as variances) the source of
        # uncertainty relevant to whichever formula produced it -- never
        # both, to avoid double-counting: kcorr_err only where a qualifying
        # observed band was used (absmag = -2.5*log10(maggies_obs) - dmod -
        # kcorr there, so kcorr's own model uncertainty is missing from the
        # observational-noise-only ivarabsmag computed by kcorr_and_absmag
        # above); absmag_synth_err only in the fallback case (absmag =
        # synth_absmag exactly there, previously reported with zero
        # uncertainty).
        absmag_var = np.zeros_like(ivarabsmag)
        with np.errstate(divide='ignore', invalid='ignore'):
            absmag_var[good_oband] = 1. / ivarabsmag[good_oband] + kcorr_err[good_oband]**2
        absmag_var[~good_oband] = absmag_synth_err[~good_oband]**2
        with np.errstate(divide='ignore', invalid='ignore'):
            ivarabsmag = np.where(absmag_var > 0., 1. / absmag_var, 0.)

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
    derived = {'LOGMSTAR': refined_logmstar, 'SFR': refined_sfr, 'SSFR': refined_ssfr,
              'T50': refined_t50, 'MSTARAGE': refined_mstarage}
    posterior = {'LOGMSTAR': logmstar, 'SFR': sfr_per_template, 'SSFR': bg_data.ssfr,
                'T50': bg_data.t50, 'MSTARAGE': bg_data.mstarage}

    for col in bg_data.axis_columns:
        pname = bg_data.axis_outname[col]
        result[pname] = fit_value[col]

    for pname in bg_data.param_names:
        order = sorted_vals = uniq = inv = None
        if pname in derived:
            result[pname] = derived[pname]
            vals = posterior[pname]
        else:
            # Grid-only data, identical for every object -- fetch the
            # cached argsort/unique instead of recomputing them here.
            vals, order, sorted_vals, uniq, inv = bg_data.axis_posterior_cache(bg_data.outname_to_axis[pname])

        if pname in bg_data.fixed_outnames:
            # Single-grid-point axis (e.g. UMIN/GAMMA): `vals` is
            # bit-identical for every template, so the weighted computation
            # below would return floating-point roundoff noise (e.g.
            # ERR ~ 1e-18) instead of the true, exact zero -- report the
            # fixed value directly and skip it.
            fixed_value = vals[0]
            result[f'{pname}_MEAN'] = fixed_value
            result[f'{pname}_MODE'] = fixed_value
            result[f'{pname}_P25'] = fixed_value
            result[f'{pname}_P50'] = fixed_value
            result[f'{pname}_P75'] = fixed_value
            result[f'{pname}_ERR'] = 0.
            continue

        mean = np.sum(weight * vals)
        result[f'{pname}_MEAN'] = mean
        result[f'{pname}_MODE'] = _weighted_mode(vals, weight, uniq=uniq, inv=inv)
        p25, p50, p75 = _weighted_percentile(vals, weight, (25., 50., 75.), order=order, sorted_values=sorted_vals)
        result[f'{pname}_P25'] = p25
        result[f'{pname}_P50'] = p50
        result[f'{pname}_P75'] = p75

        # Formal uncertainty for every parameter: the sqrt of the weighted
        # variance of its marginalized discrete posterior (no division, so
        # no overflow guard is needed even when the posterior collapses to
        # a single dominant template).
        var = np.sum(weight * (vals - mean)**2)
        result[f'{pname}_ERR'] = np.sqrt(var)

    # dof = number of fitted bands minus the one continuous free parameter
    # (the per-template mass amplitude) solved for above (nbands_used
    # computed above, before the guard that skips the fit entirely).
    result['CHI2'] = chi2_refined
    result['NDOF'] = max(nbands_used - 1, 0)

    for iband, (band, shift) in enumerate(zip(phot.absmag_bands, phot.band_shift)):
        band = band.upper()
        shift = int(10 * shift)
        result[f'KCORR{shift:02d}_{band}'] = kcorr[iband]
        result[f'ABSMAG{shift:02d}_{band}'] = absmag[iband]
        result[f'ABSMAG{shift:02d}_SYNTH_{band}'] = synth_absmag[iband]
        result[f'ABSMAG{shift:02d}_IVAR_{band}'] = ivarabsmag[iband]
        if bg_data.has_restmaggies:
            result[f'KCORR{shift:02d}_ERR_{band}'] = kcorr_err[iband]
            result[f'ABSMAG{shift:02d}_SYNTH_ERR_{band}'] = absmag_synth_err[iband]

    for key in LUM_KEYS:
        result[key] = lums[key]
    result['DN4000_MODEL'] = dn4000_model

    # Formal uncertainty on DN4000_MODEL/LOGL_* from the weighted 2nd moment
    # of their marginalized discrete posterior, over the *full* grid (no
    # top-K truncation needed: both depend only on tabulated per-template
    # scalars, not on a rebuilt spectrum, so this is a cheap vectorized
    # calculation, not an extra disk read). DN4000_MODEL is intensive (the
    # fitted mass amplitude cancels in the ratio, like a grid axis); LOGL_*/
    # LOGLNU_* are extensive (linear in amplitude), so log10(amplitude *
    # L_permass) = amplitude_log + log10(L_permass) -- deliberately
    # amplitude_log (Msun formed), not logmstar (surviving mass): L_permass
    # is normalized per Msun *formed*, matching the flux templates
    # themselves.
    dn4000_mean = np.sum(weight * bg_data.dn4000_model_permass)
    result['DN4000_MODEL_ERR'] = np.sqrt(np.sum(weight * (bg_data.dn4000_model_permass - dn4000_mean)**2))
    for key in LUM_KEYS:
        logl_per_template = amplitude_log + np.log10(np.clip(bg_data.lum_permass[key], 1e-30, None))
        logl_mean = np.sum(weight * logl_per_template)
        result[f'{key}_ERR'] = np.sqrt(np.sum(weight * (logl_per_template - logl_mean)**2))

    if topk > 0:
        idx = np.argsort(weight)[::-1][:topk]
        result['TOPK_INDEX'][:len(idx)] = idx.astype('i4')
        result['TOPK_WEIGHT'][:len(idx)] = weight[idx].astype('f4')

    return meta, result


def fastbayes_qa_one(iobj, meta, result, qadir='.', coadd_type='healpix', ndraw=50):
    """Regenerate the QA figure for one object from an already-written
    FASTBAYES output file, without repeating the fit.

    Rebuilds everything :func:`_fastbayes_qa_one` needs purely from data
    already present in the output file's ``METADATA``/``FASTBAYES``
    extensions: the observed photometry (``FLUX_*``/``FLUX_IVAR_*``
    columns), redshift (``Z``) and photometric system (``PHOTSYS``). The
    refined rest-frame spectrum is not stored in the output file (see
    :func:`write_fastbayes`) and is instead rebuilt here, from the same
    cheap, vectorized full-grid solve (:func:`_solve_grid`) that also
    reconstructs the per-template posterior weights, via
    :func:`_build_refined_spectrum` -- so this always touches the raw
    templates file (unlike the original design, where ``--ndraw 0`` could
    skip it entirely).

    Always produces a figure, even for an object with no real fit (redshift
    beyond the grid, or <=1 usable photometric band -- see ``fit_ok`` in
    :func:`_fastbayes_qa_one`): the cutout, observed photometry, and
    redshift are always genuine data independent of whether a model could
    be solved for, so only the model curve/synthesized photometry, the
    Bayesian-fit text summary, and the posterior-panel histograms are
    omitted in that case.

    Parameters
    ----------
    iobj : :class:`int`
        Index of the object in the input list, used for log messages.
    meta : :class:`astropy.table.Row`
        ``METADATA`` row for this object, from the FASTBAYES output file.
    result : :class:`numpy.void`
        ``FASTBAYES`` results row for this object, from the FASTBAYES
        output file.
    qadir : :class:`str`, optional
        Output directory for the QA figure. Default is ``'.'``.
    coadd_type : :class:`str`, optional
        Coadd type, used to build the QA target label/filename. Not stored
        in the output file, so it must be supplied by the caller. Default
        is ``'healpix'``.
    ndraw : :class:`int`, optional
        Number of additional grid templates to draw (without replacement,
        probability-weighted by the discrete posterior) and overplot on the
        SED panel as a visual sense of the model uncertainty, alongside the
        one refined maximum-likelihood curve. ``0`` disables this (default
        ``50``). Draws are capped at the number of templates with strictly
        positive posterior weight (typically a small fraction of the grid
        -- no attempt is made to force exactly ``ndraw`` draws via
        replacement).

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

    nanomaggies = np.array([meta[f'FLUX_{band.upper()}'] for band in phot.bands], dtype='f8')
    nanomaggies_ivar = np.array([meta[f'FLUX_IVAR_{band.upper()}'] for band in phot.bands], dtype='f8')
    lambda_eff = phot.filters[photsys].effective_wavelengths.value

    phot_tbl = Photometry.parse_photometry(
        phot.bands, maggies=nanomaggies, ivarmaggies=nanomaggies_ivar,
        lambda_eff=lambda_eff, min_uncertainty=phot.min_uncertainty)
    flam = np.asarray(phot_tbl['flam'], dtype='f8')
    flam_ivar = np.asarray(phot_tbl['flam_ivar'], dtype='f8') * phot.bands_to_fit

    restwave = restflux = zwave = zflux = synth_maggies = posterior_arrays = None
    family_zflux = np.empty((0, 0))

    # Same two degenerate conditions guarded against in fastbayes_one (see
    # there for why <=1 usable band makes the grid solve meaningless): QA
    # still gets a figure below (cutout, observed photometry, redshift,
    # placeholder posterior panels) -- it just has no model/fit to show.
    nbands_used = int(np.sum(flam_ivar > 0.))
    fit_ok = (redshift <= bg_data.redshift[-1]) and (nbands_used > 1)

    if not fit_ok:
        if redshift > bg_data.redshift[-1]:
            log.warning(f'Object {iobj} [{phot.uniqueid_col.lower()} {meta[phot.uniqueid_col]}] redshift '
                       f'{redshift:.6f} exceeds the grid maximum {bg_data.redshift[-1]:.6f}; '
                       f'QA will show the cutout and observed photometry only.')
        else:
            log.warning(f'Object {iobj} [{phot.uniqueid_col.lower()} {meta[phot.uniqueid_col]}] has only '
                       f'{nbands_used} usable photometric band(s); QA will show the cutout and '
                       f'observed photometry only.')
    else:
        soln = _solve_grid(flam, flam_ivar, lambda_eff, photsys, redshift)
        weight = soln['weight']
        amplitude = soln['amplitude']
        corner_idx = soln['corner_idx']
        corner_weight = soln['corner_weight']
        refined_model_maggies = soln['refined_model_maggies']
        refined_amplitude = soln['refined_amplitude']

        restwave = bg_data.template_wave()
        restflux = _build_refined_spectrum(corner_idx, corner_weight, refined_amplitude)

        logmstar = np.log10(np.clip(amplitude * bg_data.mstar, 1e-30, None)) # [ntemplate]; surviving mass
        sfr_per_template = amplitude * bg_data.sfr # [ntemplate], Msun/yr
        posterior = {'LOGMSTAR': logmstar, 'SFR': sfr_per_template, 'SSFR': bg_data.ssfr,
                    'T50': bg_data.t50, 'MSTARAGE': bg_data.mstarage}

        posterior_arrays = {}
        for pname in bg_data.param_names:
            if pname in posterior:
                vals = posterior[pname]
            else:
                vals = bg_data.axis_posterior_cache(bg_data.outname_to_axis[pname])[0]
            posterior_arrays[pname] = (vals, weight)

        zwave = restwave * (1. + redshift)
        igm_trans = _igm.full_IGM(redshift, zwave)
        dlum = _cosmo.luminosity_distance(redshift) # [Mpc]
        zfactor = igm_trans * (10. / (1e6 * dlum))**2 / (1. + redshift)
        zflux = restflux * zfactor # [erg/s/cm2/A, observed frame]

        synth_maggies = refined_model_maggies * refined_amplitude # [nband], observed-frame maggies

        # Additional templates drawn (without replacement) with probability
        # proportional to the discrete posterior weight, purely for a visual
        # sense of the model uncertainty in the SED panel -- distinct from the
        # one refined (N-linearly interpolated) maximum-likelihood spectrum
        # above. Seeded per-object so repeated QA regeneration is reproducible.
        if ndraw > 0:
            nonzero = np.flatnonzero(weight > 0.)
            ndraw_actual = min(ndraw, len(nonzero)) # don't force replacement/padding
            if ndraw_actual > 0:
                rng = np.random.default_rng(int(meta[phot.uniqueid_col]))
                p = weight[nonzero] / weight[nonzero].sum()
                draw_idx = rng.choice(nonzero, size=ndraw_actual, replace=False, p=p)
                family_zflux = np.array([bg_data.template_flux_row(idx) * amplitude[idx]
                                         for idx in draw_idx]) * zfactor

    data = {'photometry': {
        'nanomaggies': nanomaggies,
        'nanomaggies_ivar': nanomaggies_ivar,
        'lambda_eff': lambda_eff,
    }}

    _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                      zwave, zflux, synth_maggies, redshift, family_zflux=family_zflux,
                      fit_ok=fit_ok, coadd_type=coadd_type, outdir=qadir)


def _posterior_binwidth(pname):
    """Natural minimum posterior-histogram bin width for output parameter
    ``pname`` (see the posterior panel in :func:`_fastbayes_qa_one`), or
    ``None`` if no natural width applies (fall back to a fixed bin count
    instead). A narrow posterior sliced into a fixed bin *count* reads as
    noisy scatter rather than a clean, physically meaningful histogram."""
    if pname in _QA_LOG_TO_LINEAR:
        # displayed/binned linearly (see _QA_LOG_TO_LINEAR), so a fixed
        # dex-space width doesn't apply here -- fall back to a fixed bin
        # count instead.
        return None
    if pname == 'LOGMSTAR' or bg_data.outname_to_axis.get(pname) in bg_data.log_axes:
        return 0.05 # [dex]
    if pname in _GYR_BINWIDTH_PARAMS:
        return 0.3 # [Gyr]
    return None


def _highlight_ml_bar(patches, ml_value):
    """Recolor whichever bar/bin of a posterior panel contains the
    reported (ML) value, so it stays identifiable even when the panel is
    zoomed to the broader marginalized posterior rather than centered on
    it. No-op if ``patches`` is ``None`` (nothing was plotted)."""
    if patches is None:
        return
    for patch in patches:
        x0 = patch.get_x()
        if x0 <= ml_value <= x0 + patch.get_width():
            patch.set_facecolor('C0')
            patch.set_alpha(1.0)
            return


def _sfr_like_hist(ax, vals, w, ml_value):
    """Histogram an SFR/SSFR-like posterior for a QA panel.

    SFR and SSFR are reported and weighted in linear space so that an
    exact zero (a passively evolving population) is represented exactly,
    but the nonzero population can span several decades -- for which a
    linear-space histogram crushes all but the largest values into one
    bin. The strictly-positive subset is binned in log space instead, at
    a fixed 0.3~dex minimum width, and the axis uses a symmetric-log
    (``symlog``) scale so the exact-zero population still has a
    well-defined position near the origin alongside the log-spaced
    nonzero bins (a true log scale cannot represent zero at all).

    Returns
    -------
    :class:`matplotlib.container.BarContainer` or None
        The histogram's bar patches (for :func:`_highlight_ml_bar`), or
        ``None`` if the entire posterior is exactly zero, in which case
        nothing is plotted.

    """
    positive = vals > 0.
    if not np.any(positive):
        return None

    vals_pos, w_pos = vals[positive], w[positive]
    logvals = np.log10(vals_pos)
    lo, hi = _weighted_percentile(logvals, w_pos, (0.5, 99.5))
    if ml_value > 0.:
        lo, hi = min(lo, np.log10(ml_value)), max(hi, np.log10(ml_value))
    pad = 0.05 * max(hi - lo, 1e-6)

    binwidth = 0.3 # [dex]
    lo_edge = binwidth * np.floor((lo - pad) / binwidth)
    hi_edge = binwidth * np.ceil((hi + pad) / binwidth)
    hi_edge = max(hi_edge, lo_edge + binwidth) # at least one bin
    edges = 10.**np.arange(lo_edge, hi_edge + binwidth, binwidth)

    ax.set_xscale('symlog', linthresh=edges[0])
    _, _, patches = ax.hist(vals, bins=np.concatenate(([0.], edges)), weights=w,
                            color='gray', edgecolor='k', alpha=0.8)
    return patches


def _fastbayes_qa_one(data, meta, result, posterior_arrays, restwave, restflux,
                      zwave, zflux, synth_maggies, redshift, family_zflux=None,
                      fit_ok=True, coadd_type='healpix', outdir='.'):
    """Generate a QA figure for one Bayesian-fit object.

    Reuses :func:`fastspecfit.qa._target_label` and
    :func:`fastspecfit.qa._fetch_cutout`, which are generic utilities with
    no fastspec/fastphot-specific coupling. Called from
    :func:`fastbayes_qa_one` to regenerate QA from an already-written
    FASTBAYES output file (``bin/fastbayes-qa``) -- QA generation is no
    longer available inline from :func:`fastbayes_one` during fitting.

    Parameters
    ----------
    data : :class:`dict`
        Per-object data dictionary.
    meta : :class:`astropy.table.Row`
        Metadata row for this object.
    result : :class:`numpy.ndarray`
        Bayesian-fitting output row, from :func:`fastbayes_one`.
    posterior_arrays : :class:`dict` or None
        Mapping of parameter name to ``(values, weights)`` per-template
        posterior arrays. ``None`` when ``fit_ok`` is ``False``.
    restwave, restflux : :class:`numpy.ndarray` or None
        Refined rest-frame maximum-likelihood spectrum (``restwave`` is not
        otherwise used by this function). Both ``None`` when ``fit_ok`` is
        ``False``.
    zwave, zflux : :class:`numpy.ndarray` or None
        The same spectrum redshifted (and IGM/distance-attenuated) to the
        observed frame. ``None`` when ``fit_ok`` is ``False``.
    synth_maggies : :class:`numpy.ndarray` or None
        Observed-frame photometry (AB maggies) synthesized from the refined
        maximum-likelihood model in each band. ``None`` when ``fit_ok`` is
        ``False``.
    redshift : :class:`float`
    family_zflux : :class:`numpy.ndarray` or None, optional
        Shape ``(ndraw, npix)`` observed-frame spectra of additional
        templates drawn probability-weighted from the discrete posterior
        (:func:`fastbayes_qa_one`), overplotted faintly behind the one
        refined maximum-likelihood curve as a visual sense of the model
        uncertainty. Empty or ``None`` (default) plots none.
    fit_ok : :class:`bool`, optional
        Whether a real grid fit was actually performed for this object
        (default ``True``). ``False`` means the redshift was beyond the
        grid, or too few photometric bands were usable to constrain a fit
        (see :func:`fastbayes_one`/:func:`fastbayes_qa_one`) -- the figure
        is still generated (cutout, observed photometry, redshift), but the
        model curve/synthesized photometry, the Bayesian-fit text summary,
        and the posterior-panel histograms are omitted (empty placeholder
        panels only) since there is nothing genuine to show for them.
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
    legfntsz = 18
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

    def _fmt(val, err, fmt):
        if err > 0.:
            return fmt.format(val) + r'\pm' + fmt.format(err)
        return fmt.format(val)

    def _log_to_linear(val, err):
        # d(10**x)/dx = ln(10) * 10**x, so this never divides by val/err
        # (which, unlike a naive value_linear/value_log ratio, is safe even
        # when the log-space value is exactly 0, e.g. Age = 1 Gyr or
        # Z/Zsun = 1).
        linear_val = 10.**val
        linear_err = np.log(10.) * linear_val * err
        return linear_val, linear_err

    pngfile = get_qa_filename(meta, coadd_type, outprefix='fastbayes', outdir=outdir,
                              uniqueid_col=phot.uniqueid_col)
    figdir = os.path.dirname(pngfile)
    if figdir and not os.path.isdir(figdir):
        os.makedirs(figdir, exist_ok=True)

    # 5 rows x 8 cols: rows 0:3 are the SED/cutout block (matching the
    # column proportions of fastqa's fastphot layout: sedax 5/8, cutax
    # 3/8, so labels/legends sized for that layout still fit), row 3 is
    # a blank gap (also hosting the z/Dn4000/absmag boxes below the
    # cutout), and row 4 holds the (2-row x 3-col, 6-parameter) posterior grid.
    # Explicit margins mirror fastqa's own fastphot subplots_adjust (tight
    # left/bottom, generous right/top for the fig.text labels/legends).
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(5, 8, height_ratios=[1, 1, 1, 0.7, 3.5],
                          left=0.09, right=0.92, top=0.9, bottom=0.07,
                          hspace=0.1, wspace=0.25)

    sedax = fig.add_subplot(gs[0:3, 0:5])

    # image cutout, only if this photometry configuration has viewer info
    cutgs = gs[0:2, 5:8]
    have_cutout = hasattr(phot, 'viewer_layer') and hasattr(phot, 'viewer_pixscale')
    if have_cutout:
        img, wcs, _, _ = _fetch_cutout(meta, figdir, pngfile, phot.viewer_layer, phot.viewer_pixscale)
        cutax = fig.add_subplot(cutgs, projection=wcs)
        # WCS projection forces equal aspect, and this cell is wider than it
        # is tall, so matplotlib shrinks the axes box to a square; anchoring
        # 'NW' keeps that box flush against sedax (and the cell's top edge)
        # instead of centering it, which would leave a gap on both sides.
        cutax.set_anchor('NW')
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
        cutax = fig.add_subplot(cutgs)
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
    # beyond phot_wavelims doesn't skew the y-axis limits below. Skipped
    # entirely when fit_ok is False (no model to show -- see fit_ok in the
    # docstring above): the SED panel then shows only real, model-
    # independent data (observed photometry, both axes' wavelength scales).
    if fit_ok:
        factor = 10**(0.4 * 48.6) * zwave**2 / (C_LIGHT * 1e13) # [erg/s/cm2/A --> maggies]
        zwave_um = zwave / 1e4
        mgood = (zflux > 0.) & (zwave_um >= phot_wavelims[0]) & (zwave_um <= phot_wavelims[1])
        sedmodel_abmag = np.full_like(zflux, np.nan)
        sedmodel_abmag[mgood] = -2.5 * np.log10(zflux[mgood] * factor[mgood])
    else:
        mgood = np.array([], dtype=bool)

    # y-axis limits (AB mag; brighter/smaller values at the top) -- set
    # *before* any errorbar(lolims=True) calls below: matplotlib bakes in
    # the upper-limit caret's orientation using the axis's inversion state
    # at call time, so setting (and thus inverting) the y-axis afterward
    # leaves the caret pointing the wrong way (toward the point instead of
    # away from it, into the fainter/undetected region).
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
        sed_ymin = min(np.nanmax(allmags) + dm, 32.) # never fainter than AB=32
        sed_ymax = np.nanmin(allmags) - dm
    else:
        sed_ymin, sed_ymax = 30., 20.

    sedax.set_xlim(phot_wavelims[0], phot_wavelims[1])
    sedax.set_ylim(sed_ymin, sed_ymax)

    if fit_ok:
        # additional posterior-weighted draws, plotted faintly behind the one
        # refined maximum-likelihood curve as a visual sense of the model
        # uncertainty (same grey hue, distinguished only by opacity/linewidth so
        # as not to clash with the photometry's orange color scheme below)
        if family_zflux is not None and len(family_zflux) > 0:
            for fam_flux in family_zflux:
                fam_mgood = (fam_flux > 0.) & (zwave_um >= phot_wavelims[0]) & (zwave_um <= phot_wavelims[1])
                if np.any(fam_mgood):
                    fam_abmag = -2.5 * np.log10(fam_flux[fam_mgood] * factor[fam_mgood])
                    sedax.plot(zwave_um[fam_mgood], fam_abmag, color='grey', alpha=0.12, lw=0.5, zorder=0)

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

    # reduced chi2 (top-left, no box) -- omitted when fit_ok is False: CHI2/
    # NDOF are both exactly 0 in that case (no fit was attempted), and
    # "chi2=0.00" would misleadingly read as a perfect fit rather than "no
    # fit".
    if fit_ok:
        rchi2 = result['CHI2'] / result['NDOF'] if result['NDOF'] > 0 else result['CHI2']
        sedax.text(0.02, 0.94, r'$\chi^{2}_{\nu,\mathrm{phot}}=$' + r'${:.2f}$'.format(rchi2),
                  ha='left', va='top', transform=sedax.transAxes, fontsize=legfntsz)

    # target label above the cutout, and rest-frame wavelength label above the SED.
    # Use the cutout's *gridspec cell* position, not cutax.get_position(): the
    # WCS-projected axes has equal aspect and gets shrunk to a (now
    # left-anchored) square within that wider cell, so the actual axes box is
    # narrower than the allocated column -- anchoring labels/boxes to the full
    # cell instead keeps them aligned with the RA/Dec axis labels, which sit
    # outside the shrunk box near the cell's edges.
    cpos = cutgs.get_position(fig)
    spos = sedax.get_position()
    fig.text(cpos.x0, cpos.y1 + 0.02, '\n'.join(_target_label(meta, coadd_type, uniqueid_col=phot.uniqueid_col)),
             ha='left', va='bottom', fontsize=fontsize2, linespacing=1.4)
    fig.text((spos.x0 + spos.x1) / 2., spos.y1 + 0.05, r'Rest-frame Wavelength ($\mu$m)',
             ha='center', va='bottom', fontsize=fontsize2)

    # Two shaded info boxes below the cutout (below its RA/Dec axis label and
    # tick labels, not just its image edge): observed/spectrum-derived
    # quantities on the left (z, Dn(4000), rest-frame absolute magnitude and
    # colors), Bayesian fitting results on the right (age, metallicity,
    # stellar mass, SFR -- in the same order as the posterior panel below).
    # Anchoring the left box at the column's left edge (ha='left') and the
    # right box at the column's right edge (ha='right') maximizes the gap
    # between them for a given cutout column width, rather than splitting at
    # a fixed fraction that can overlap if either box's text runs long.
    ytext = cpos.y0 - 0.08

    gindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1. + phot.band_shift) - 4300))
    rindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1. + phot.band_shift) - 5600))
    zindx = np.argmin(np.abs(phot.absmag_filters.effective_wavelengths.value / (1. + phot.band_shift) - 8100))
    absmag_gband, absmag_rband, absmag_zband = phot.absmag_bands[gindx], phot.absmag_bands[rindx], phot.absmag_bands[zindx]
    shift_gband, shift_rband, shift_zband = phot.band_shift[gindx], phot.band_shift[rindx], phot.band_shift[zindx]

    def _absmag_col(band, shift):
        return f'ABSMAG{int(10 * shift):02d}_{band.upper()}'

    # DN4000_MODEL/absolute-magnitude/color lines all come from the refined
    # rest-frame model, so they're meaningless (still at their all-zeros
    # init) when fit_ok is False -- z is real data either way, so it's the
    # only line kept in that case.
    txt = [r'$z={:.7f}$'.format(redshift)]
    if fit_ok:
        txt += [r'$D_{{n}}(4000)_{{\mathrm{{model}}}}={:.3f}$'.format(result['DN4000_MODEL']), '']
        txt += [r'$M_{{{}{}}}={:.2f}$'.format(
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
    fig.text(cpos.x0, ytext, '\n'.join(txt), ha='left', va='top', fontsize=fontsize1,
             bbox=bbox, linespacing=1.6)

    # Iterate in bg_data.param_names order (i.e. the AXES order the
    # posterior panel below also follows), not _QA_SUMMARY_CANDIDATES'
    # own tuple order, so the two stay visually aligned. A short placeholder
    # note stands in for the whole box when fit_ok is False, rather than
    # printing e.g. LOGMSTAR's all-zeros init as if it were a measurement.
    if fit_ok:
        _qa_summary = {pname: (template, fmt) for pname, template, fmt in _QA_SUMMARY_CANDIDATES}
        txt = []
        for pname in bg_data.param_names:
            if pname not in _qa_summary:
                continue
            template, fmt = _qa_summary[pname]
            val, err = result[pname], result[f'{pname}_ERR']
            if pname in _QA_LOG_TO_LINEAR:
                val, err = _log_to_linear(val, err)
            txt.append(template.format(_fmt(val, err, fmt)))
    else:
        txt = ['No fit', '(insufficient', 'photometry or', 'redshift beyond', 'the grid)']
    fig.text(cpos.x1+0.04, ytext, '\n'.join(txt), ha='right', va='top', fontsize=fontsize1,
             bbox=bbox, linespacing=1.6)

    # --- posterior panel: weighted 1D marginal histogram of every *free*
    # fitted parameter (bg_data.param_names minus bg_data.fixed_outnames --
    # a fixed axis like GAMMA/UMIN has exactly one possible value, so its
    # panel is never informative -- and minus _QA_HIDDEN_PARAMS, raw Dense
    # Basis SFH bookkeeping axes not physically meaningful on their own;
    # T50/MSTARAGE are the human-readable substitutes), plus LOGMSTAR/SFR/
    # SSFR/T50/MSTARAGE. A standalone gridspec (not a subgridspec of `gs`)
    # so its right edge can extend past gs's own right margin to cpos.x1 +
    # 0.04, matching the right-hand gray box/Dec label; top/bottom match
    # row 4 of `gs` exactly. ncols is fixed at 4 (a layout choice); nrows
    # adapts to however many free parameters this grid actually has.
    free_param_names = [pname for pname in bg_data.param_names
                        if pname not in bg_data.fixed_outnames and pname not in _QA_HIDDEN_PARAMS]
    ncols = 5
    nrows = -(-len(free_param_names) // ncols) # ceil division
    row4pos = gs[4, 0:8].get_position(fig)
    post_gs = fig.add_gridspec(nrows, ncols, left=spos.x0-0.02, right=cpos.x1 + 0.04,
                               bottom=row4pos.y0, top=row4pos.y1-0.02, wspace=0.1, hspace=0.6)
    for i, pname in enumerate(free_param_names):
        row, col = divmod(i, ncols)
        ax = fig.add_subplot(post_gs[row, col])

        if not fit_ok:
            # Placeholder-only panel: no posterior to show (see fit_ok in
            # the docstring above), so just the empty frame and axis label,
            # keeping the same panel layout/count as a real fit's figure.
            ax.tick_params(labelbottom=False)
            ax.set_xlabel(_QA_LINEAR_LABELS.get(pname, bg_data.param_labels.get(pname, pname)), fontsize=16)
            ax.tick_params(labelleft=False, labelsize=13)
            continue

        vals, w = posterior_arrays[pname]
        ml_value = result[pname]

        # LOGAGE/LOGZZSUN are plotted here (and summarized in the text box
        # above) as linear Age/Z-Zsun rather than their stored log10 form;
        # SSFR is rescaled yr^-1 -> Gyr^-1. See _log_to_linear/
        # _QA_LOG_TO_LINEAR and _DERIVED_PARAM_LABELS['SSFR'].
        if pname in _QA_LOG_TO_LINEAR:
            vals, ml_value = 10.**vals, 10.**ml_value
        elif pname == 'SSFR':
            vals, ml_value = vals * 1e9, ml_value * 1e9

        uniq, inv = np.unique(vals, return_inverse=True)
        patches = None

        if len(uniq) <= _DISCRETE_UNIQ_MAX:
            # discrete axis, or any quantity whose posterior happens to
            # collapse onto a handful of values: one bar per actual value
            # avoids mostly-empty uniform bins from arbitrary bin edges
            # landing between (or splitting) the real, finite set of
            # values. A single distinct value (a fixed axis like GAMMA/
            # UMIN, or a posterior that collapses onto one value) has
            # nothing to show beyond the vertical ML line below.
            if len(uniq) > 1:
                binweight = np.bincount(inv, weights=w, minlength=len(uniq))
                # Per-point width (0.8x the narrower of its two neighbor
                # gaps) rather than one global width: axes that are
                # log-uniformly spaced in their stored (log) form -- e.g.
                # LOGAGE/LOGZZSUN converted to linear above -- are not
                # evenly spaced once linear, so a single global min-diff
                # width would make the widely-separated high-value bars
                # a barely visible sliver. Reduces to the original global
                # width for genuinely uniform axes.
                diffs = np.diff(uniq)
                width = 0.8 * np.minimum(np.r_[diffs[0], diffs], np.r_[diffs, diffs[-1]])
                patches = ax.bar(uniq, binweight, width=width, color='gray', edgecolor='k', alpha=0.8)
        elif pname in ('SFR', 'SSFR'):
            patches = _sfr_like_hist(ax, vals, w, ml_value)
        else:
            # many distinct values (LOGMSTAR, T50, MSTARAGE, or any axis
            # with a large point count): bin at a natural minimum width
            # (_posterior_binwidth), zoomed to where the posterior weight
            # actually is. The zoom window always includes the reported
            # (ML) value, even if it falls in a low-weight tail outside
            # the 0.5-99.5 percentile band.
            lo, hi = _weighted_percentile(vals, w, (0.5, 99.5))
            lo, hi = min(lo, ml_value), max(hi, ml_value)
            if hi > lo:
                pad = 0.05 * (hi - lo)
                binwidth = _posterior_binwidth(pname)
                if binwidth is not None:
                    lo_edge = binwidth * np.floor((lo - pad) / binwidth)
                    hi_edge = binwidth * np.ceil((hi + pad) / binwidth)
                    hi_edge = max(hi_edge, lo_edge + binwidth) # at least one bin
                    bins = np.arange(lo_edge, hi_edge + binwidth, binwidth)
                else:
                    bins = np.linspace(lo - pad, hi + pad, 31)
                _, _, patches = ax.hist(vals, bins=bins, weights=w, color='gray', edgecolor='k', alpha=0.8)

        # Recolor whichever bar/bin contains the ML value, so it stays
        # identifiable even in a panel zoomed/shaped by the broader
        # marginalized posterior rather than centered on it.
        _highlight_ml_bar(patches, ml_value)

        # An SFR/SSFR posterior with no nonzero contribution at all (every
        # template passively evolving) has no meaningful spread to show --
        # suppress the numeric x-tick labels (matplotlib's default
        # auto-ranged ticks on an otherwise-empty axis are pure noise),
        # while keeping the axis label so it's still clear what the panel is.
        if pname in ('SFR', 'SSFR') and not np.any(vals > 0.):
            ax.tick_params(labelbottom=False)

        ax.axvline(ml_value, color='C0', lw=1.5)
        ax.set_xlabel(_QA_LINEAR_LABELS.get(pname, bg_data.param_labels.get(pname, pname)), fontsize=16)
        ax.tick_params(labelleft=False, labelsize=13)

    fig.savefig(pngfile)
    plt.close(fig)
    log.info(f'Wrote {pngfile}')


def write_fastbayes(meta, results, outfile, gridfile, fphotofile,
                    template_file=None, topk=0):
    """Write the Bayesian-fitting output to a multi-extension FITS file.

    The refined rest-frame model spectrum is not stored: it is cheaply and
    deterministically rebuildable from the stored ``FLUX_*``/``FLUX_IVAR_*``/
    ``Z``/``PHOTSYS`` columns via :func:`_solve_grid` and
    :func:`_build_refined_spectrum` (see :func:`fastbayes_qa_one`), so
    persisting it per-object would only be a large, redundant disk cost.

    Parameters
    ----------
    meta : :class:`astropy.table.Table`
        Output metadata table (see :func:`fastspecfit.io.create_output_meta`).
    results : :class:`astropy.table.Table`
        Output fitting-results table.
    outfile : :class:`str`
        Full path of the output FITS file (``.gz`` triggers gzip
        compression, as in :func:`fastspecfit.io.write_fastspecfit`).
    gridfile : :class:`str`
        Full path of the Bayesian grid file used for fitting.
    fphotofile : :class:`str`
        Full path of the photometry configuration file used.
    template_file : :class:`str` or None, optional
        Full path of the raw Bayesian templates FITS file actually used
        for fitting (e.g. ``bg_data.templates_file``).
    topk : :class:`int`, optional
        Number of top-weight grid templates stored per object, if any.

    """
    import gzip, shutil, warnings
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyUserWarning
    from desiutil.depend import add_dependencies, possible_dependencies, setdep

    outdir = os.path.dirname(os.path.abspath(os.path.expanduser(os.path.expandvars(outfile))))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    if outfile.endswith('.gz'):
        tmpfile = outfile[:-3] + '.tmp'
    else:
        tmpfile = outfile + '.tmp'

    hduprim = fits.PrimaryHDU()
    add_dependencies(hduprim.header, module_names=possible_dependencies+['fastspecfit'],
                     envvar_names=('DESI_SPECTRO_REDUX', 'DUST_DIR', 'FTEMPLATES_DIR', 'FPHOTO_DIR'))
    # Recorded purely as diagnostic/provenance (mirroring write_fastspecfit):
    # these are absolute (GRIDFILE, FPHOTO_FILE) or basename-only
    # (FTEMPLATES_FILE) paths that may not exist on this machine, e.g. if the
    # fit ran elsewhere (laptop -> NERSC). fastbayes_qa requires
    # --gridfile/--fphotofile to be passed explicitly rather than reading
    # them back from here.
    setdep(hduprim.header, 'GRIDFILE', os.path.abspath(str(gridfile)))
    if fphotofile:
        setdep(hduprim.header, 'FPHOTO_FILE', str(fphotofile))
    if template_file:
        setdep(hduprim.header, 'FTEMPLATES_FILE', os.path.basename(str(template_file)))
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

    hx = fits.HDUList([hduprim, hdumeta, hduresults])

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
                        '(default: $FTEMPLATES_DIR/bayesian); ignored if --templatesfile is given.')
    parser.add_argument('--templatesfile', type=str, default=None,
                        help='Full path to the raw Bayesian templates FITS file, overriding '
                        '--templatedir/$FTEMPLATES_DIR reconstruction from the --gridfile header '
                        '(useful when the templates live somewhere non-standard, e.g. switching '
                        'between NERSC and a laptop).')
    parser.add_argument('--flux-cache-mb', dest='flux_cache_mb', type=float, default=DEFAULT_FLUX_CACHE_MB,
                        help='Preload the full raw-templates FLUX array once per worker when it is '
                        'smaller than this many MB; otherwise fall back to a bounded per-row LRU '
                        'cache over lazy per-object reads.')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads.')
    parser.add_argument('--checkpoint-size', dest='checkpoint_size', type=int, default=None,
                        help='Fit and write results in batches of this many objects at a time, '
                        'each to its own file under --checkpoint-dir, resuming from any batches '
                        'already present from a previous (e.g. interrupted) run, then merge every '
                        'batch into --outfile at the end. Recommended for very large single-shot '
                        'samples (e.g. hundreds of thousands of objects processed outside the usual '
                        'per-healpix-file workflow), where holding every object\'s result in memory '
                        'at once would otherwise be the memory bottleneck. Disabled (default) for '
                        'the usual per-file scale.')
    parser.add_argument('--checkpoint-dir', dest='checkpoint_dir', type=str, default=None,
                        help='Directory for --checkpoint-size batch files (default: a directory '
                        'named after --outfile, alongside it). Ignored if --checkpoint-size is not given.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to process in each file, zero-indexed.')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--input-redshifts', type=str, default=None,
                        help='Comma-separated list of input redshifts corresponding to the (required) --targetids input.')
    parser.add_argument('--zmin', type=float, default=None, help='Override the default minimum redshift required for modeling.')
    parser.add_argument('--topk', type=int, default=0,
                        help='Number of top-weight grid templates to save per object (sparse joint posterior); 0 disables.')
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
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    log.info('fastbayes {}'.format(' '.join(options)))

    return parser.parse_args(options)


def _assemble_fastbayes_results(out, phot, units):
    """Assemble the ``(outmeta, results)`` output tables from a list of
    per-object ``(meta, result)`` tuples returned by :func:`fastbayes_one`.

    Factored out of :func:`fastbayes` so that the checkpointed code path
    (:func:`_run_fastbayes_chunks`) can call it once per chunk, identically
    to how the non-checkpointed path calls it once for the whole catalog.

    Parameters
    ----------
    out : list of (:class:`astropy.table.Row`, :class:`numpy.ndarray`)
        Per-object ``(meta, result)`` tuples, in the order returned by
        :func:`fastbayes_one`.
    phot : :class:`fastspecfit.photometry.Photometry`
        Photometry configuration, needed to size/label the output columns.
    units : :class:`dict`
        Mapping from output column name to astropy unit, passed straight
        through to :func:`fastspecfit.io.create_output_table`.

    Returns
    -------
    outmeta : :class:`astropy.table.Table`
    results : :class:`astropy.table.Table`

    """
    from astropy.table import vstack
    from fastspecfit.io import create_output_meta, create_output_table

    out = list(zip(*out))

    allmeta = vstack(out[0])
    outmeta = create_output_meta(allmeta, phot=phot, fastphot=True)
    for band in phot.bands:
        col = f'FLUX_SYNTH_{band.upper()}'
        outmeta[col] = allmeta[col]
        outmeta[col].unit = 'nanomaggies'

    results = create_output_table(out[1], outmeta, units)

    return outmeta, results


def _checkpoint_is_complete(chunk_file, expected_nobj):
    """Whether ``chunk_file`` already holds a complete, readable batch of
    ``expected_nobj`` objects, so :func:`_run_fastbayes_chunks` can skip
    re-fitting it (resume-after-crash support).

    :func:`write_fastbayes` writes to a temporary file and only renames it
    into place once complete, so an interrupted run cannot leave a partial
    ``chunk_file`` behind -- the row-count check here is mostly a safety net
    against, e.g., resuming with a different ``--checkpoint-size`` than the
    original run used.

    """
    if not os.path.isfile(chunk_file):
        return False
    try:
        with fitsio.FITS(chunk_file) as F:
            return F['METADATA'].get_nrows() == expected_nobj
    except OSError:
        return False


def _run_fastbayes_chunks(data, meta, phot, fastbayes_dtype, units, args, mp_pool, gridfile):
    """Fit ``data``/``meta`` in resumable batches of ``args.checkpoint_size``
    objects, so that a very large single-shot sample (e.g. hundreds of
    thousands of objects run in one process, outside the usual
    per-healpix-file workflow) never needs every object's result held in
    memory at once, and an interrupted run can pick back up without redoing
    completed batches.

    Each batch is fit and written (via :func:`write_fastbayes`, unchanged)
    to its own self-contained FASTBAYES file under ``args.checkpoint_dir``;
    a batch whose file already exists and is complete
    (:func:`_checkpoint_is_complete`) is skipped entirely. Once every batch
    exists, their ``METADATA``/``FASTBAYES`` tables are concatenated to
    produce the final result -- cheap now that neither extension holds a
    per-object model spectrum.

    Parameters
    ----------
    data : list of :class:`dict`
    meta : :class:`astropy.table.Table`
        Full-catalog per-object data/metadata, from
        :meth:`fastspecfit.io.DESISpectra.read`.
    phot : :class:`fastspecfit.photometry.Photometry`
    fastbayes_dtype : :class:`numpy.dtype`
    units : :class:`dict`
    args : :class:`argparse.Namespace`
        Parsed ``fastbayes`` arguments (uses ``checkpoint_size``,
        ``checkpoint_dir``, ``outfile``, ``topk``).
    mp_pool : :class:`fastspecfit.util.MPPool`
        Worker pool, reused across every batch (the grid is loaded once,
        not re-initialized per batch).
    gridfile : :class:`str`
        Full path of the Bayesian grid file used for fitting.

    Returns
    -------
    outmeta : :class:`astropy.table.Table`
    results : :class:`astropy.table.Table`

    """
    from astropy.table import vstack

    nobj = len(meta)
    checkpoint_size = args.checkpoint_size
    checkpoint_dir = args.checkpoint_dir
    if checkpoint_dir is None:
        outbase = os.path.splitext(os.path.basename(args.outfile.removesuffix('.gz')))[0]
        checkpoint_dir = os.path.join(os.path.dirname(os.path.abspath(args.outfile)),
                                      f'{outbase}-checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)

    nchunk = int(np.ceil(nobj / checkpoint_size))
    log.info(f'Checkpointing {nobj:,d} objects in {nchunk} batches of up to '
            f'{checkpoint_size:,d} under {checkpoint_dir}.')

    chunk_files = []
    for ichunk in range(nchunk):
        start = ichunk * checkpoint_size
        stop = min(start + checkpoint_size, nobj)
        chunk_file = os.path.join(checkpoint_dir, f'chunk-{ichunk:05d}.fits')
        chunk_files.append(chunk_file)

        if _checkpoint_is_complete(chunk_file, stop - start):
            log.info(f'Batch {ichunk+1}/{nchunk} ({chunk_file}) already complete; skipping.')
            continue

        log.info(f'Fitting batch {ichunk+1}/{nchunk}: objects {start}:{stop} ({stop-start:,d} total).')

        chunk_fitargs = [{
            'iobj': iobj,
            'data': data[iobj],
            'meta': meta[iobj],
            'fastbayes_dtype': fastbayes_dtype,
            'topk': args.topk,
        } for iobj in range(start, stop)]

        t0 = time.time()
        out = mp_pool.starmap(fastbayes_one, chunk_fitargs)
        chunk_outmeta, chunk_results = _assemble_fastbayes_results(out, phot, units)
        log.info(fsftime('fastbayes_chunk', time.time() - t0, context=f'nobj={stop-start}'))

        write_fastbayes(chunk_outmeta, chunk_results, outfile=chunk_file, gridfile=gridfile,
                        fphotofile=phot.fphotofile, template_file=bg_data.templates_file,
                        topk=args.topk)

    log.info(f'Merging {nchunk} batches into {args.outfile}.')
    allmeta, allresults = [], []
    for chunk_file in chunk_files:
        with fitsio.FITS(chunk_file) as F:
            allmeta.append(Table(F['METADATA'].read()))
            allresults.append(Table(F['FASTBAYES'].read()))

    return vstack(allmeta), vstack(allresults)


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
    from fastspecfit.io import DESISpectra

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    envlist = []
    if args.redux_dir is None:
        envlist.append('DESI_SPECTRO_REDUX')
    if args.mapdir is None:
        envlist.append('DUST_DIR')
    if args.fphotodir is None:
        envlist.append('FPHOTO_DIR')
    if args.templatesfile is None and args.templatedir is None:
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
                    'templatedir': args.templatedir, 'templatesfile': args.templatesfile,
                    'mapdir': args.mapdir, 'flux_cache_mb': args.flux_cache_mb}

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

    # Pre-declare the grid-based synthetic-photometry columns (populated by
    # fastbayes_one) on the whole metadata table before it's sliced into
    # per-object rows and dispatched to the worker pool.
    for band in phot.bands:
        meta[f'FLUX_SYNTH_{band.upper()}'] = np.zeros(nobj, dtype='f4')

    fastbayes_dtype = get_fastbayes_dtype(phot, topk=args.topk)

    units = {'LOGAGE': 'dex(Gyr)'}
    units.update({f'LOGAGE_{stat}': 'dex(Gyr)' for stat in ('ERR', 'MEAN', 'MODE', 'P25', 'P50', 'P75')})

    if args.checkpoint_size is not None and args.checkpoint_size > 0 and nobj > args.checkpoint_size:
        outmeta, results = _run_fastbayes_chunks(
            data, meta, phot, fastbayes_dtype, units, args, mp_pool, gridfile)
    else:
        fitargs = [{
            'iobj': iobj,
            'data': data[iobj],
            'meta': meta[iobj],
            'fastbayes_dtype': fastbayes_dtype,
            'topk': args.topk,
        } for iobj in range(nobj)]

        t0 = time.time()
        out = mp_pool.starmap(fastbayes_one, fitargs)
        outmeta, results = _assemble_fastbayes_results(out, phot, units)
        log.info(fsftime('fastbayes_all', time.time() - t0, context=f'nobj={nobj}'))

    if _own_pool:
        mp_pool.close()

    write_fastbayes(outmeta, results, outfile=args.outfile,
                    gridfile=gridfile, fphotofile=phot.fphotofile,
                    template_file=bg_data.templates_file, topk=args.topk)

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
                        choices=['healpix', 'uniqpix', 'cumulative', 'pernight', 'perexp', 'custom',
                                'stacked', 'external'],
                        help='Coadd type, used to build the QA target label/filename (not stored in the '
                        'FASTBAYES output file, so it cannot be inferred). Use "external" for non-DESI '
                        'samples (e.g., built from a custom photometry file), identified by the '
                        'photometry configuration\'s uniqueid column rather than SURVEY/PROGRAM/HEALPIX/'
                        'TARGETID.')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of TARGETIDs to process.')
    parser.add_argument('--uniqueids', type=str, default=None,
                        help='Comma-separated list of integer IDs to process, matched against the '
                        'photometry configuration\'s uniqueid column instead of TARGETID. Use for '
                        'external/non-DESI samples (e.g., --coadd-type=external); ignored if '
                        '--targetids is also given.')
    parser.add_argument('-n', '--ntargets', type=int, help='Number of objects to process.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to process, zero-indexed.')
    parser.add_argument('--ndraw', type=int, default=50,
                        help='Number of additional grid templates to draw (probability-weighted by '
                        'the discrete posterior, without replacement) and overplot on the SED panel '
                        'as a visual sense of the model uncertainty, alongside the one refined '
                        'maximum-likelihood spectrum (always rebuilt from --templatedir/'
                        '--templatesfile below; 0 just disables the extra overplotted draws).')
    parser.add_argument('--templatedir', type=str, default=None,
                        help='Top-level location of the raw Bayesian templates file '
                        '(default: $FTEMPLATES_DIR/bayesian); ignored if --templatesfile is given. '
                        'Always required, to rebuild the refined rest-frame spectrum.')
    parser.add_argument('--templatesfile', type=str, default=None,
                        help='Full path to the raw Bayesian templates FITS file, overriding '
                        '--templatedir/$FTEMPLATES_DIR reconstruction. Always required, to rebuild '
                        'the refined rest-frame spectrum.')
    parser.add_argument('--flux-cache-mb', dest='flux_cache_mb', type=float, default=DEFAULT_FLUX_CACHE_MB,
                        help='Preload the full raw-templates FLUX array once per worker when it is '
                        'smaller than this many MB; otherwise fall back to a bounded per-row LRU '
                        'cache over lazy per-object reads.')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing threads.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose (for debugging purposes).')

    if options is None:
        options = sys.argv[1:]

    log.info('fastbayes-qa {}'.format(' '.join(options)))

    return parser.parse_args(options)


def fastbayes_qa(args=None, mp_pool=None):
    """Regenerate QA figures for a FASTBAYES output file, without repeating the fit.

    Reads back an already-written FASTBAYES output file's ``METADATA`` and
    ``FASTBAYES`` extensions, then calls :func:`fastbayes_qa_one` for each
    requested object. No raw DESI spectra are read, so this does not need
    the original redrock/spectra files to be available. It does, however,
    always need the raw Bayesian templates file
    (``$FTEMPLATES_DIR``/``--templatedir``/``--templatesfile``): the refined
    rest-frame spectrum is not stored in the output file (see
    :func:`write_fastbayes`) and is always rebuilt from the raw templates,
    via :func:`_build_refined_spectrum` (cheap in practice, since
    :meth:`BayesianGrid.template_flux_row` preloads the whole FLUX array for
    grids under ``--flux-cache-mb``).

    ``--gridfile``/``--fphotofile`` must be passed explicitly (matching
    whatever was used to build the grid and run the fit) rather than being
    read back from the output file's header: the header's ``GRIDFILE`` and
    ``FPHOTO_FILE`` dependency entries (see :func:`write_fastbayes`) are
    absolute paths that may not exist on this machine (e.g. the fit ran on a
    different machine than this QA regeneration) -- they are purely
    diagnostic/provenance and not reliably reloadable without user input.

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

    nobj = len(meta)
    keep = np.arange(nobj)

    if args.targetids is not None:
        targetids = [int(x) for x in args.targetids.split(',')]
        keep = np.where(np.isin(meta['TARGETID'], targetids))[0]
        if len(keep) == 0:
            log.warning('No objects match the requested --targetids.')
            return 0
    elif args.uniqueids is not None:
        # For external/non-DESI samples, matched against the photometry
        # configuration's uniqueid column instead of the literal TARGETID
        # column; --targetids (above) is left untouched so default DESI
        # behavior is unaffected.
        uniqueid_col = Photometry(fphotofile=fphotofile).uniqueid_col
        uniqueids = [int(x) for x in args.uniqueids.split(',')]
        keep = np.where(np.isin(meta[uniqueid_col], uniqueids))[0]
        if len(keep) == 0:
            log.warning('No objects match the requested --uniqueids.')
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

    init_argdict = {'fphotofile': fphotofile, 'gridfile': gridfile,
                    'templatedir': args.templatedir, 'templatesfile': args.templatesfile,
                    'cutout_unreachable': cutout_unreachable,
                    'flux_cache_mb': args.flux_cache_mb}

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
        'qadir': args.qadir,
        'coadd_type': args.coadd_type,
        'ndraw': args.ndraw,
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
