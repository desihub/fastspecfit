"""
fastspecfit.test.test_io
========================
Unit tests for pure, file-I/O-free functions in io.py:
  - select()         — filter metadata/specphot tables by healpix / tile / night
  - get_qa_filename() — build QA PNG filenames from metadata

The heavy I/O methods (DESISpectra.read, write_fastspecfit, etc.) are
covered by the integration tests in test_fastspecfit.py.
"""
import numpy as np
import pytest
from astropy.table import Table


# ── shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def healpix_tables():
    """Small metadata + specphot Tables with healpix-coadd columns."""
    meta = Table({
        'TARGETID': np.array([1, 2, 3, 4, 5]),
        'HEALPIX':  np.array([100, 200, 100, 300, 200]),
        'SURVEY':   ['main', 'main', 'sv3', 'main', 'sv3'],
        'PROGRAM':  ['bright', 'dark', 'bright', 'dark', 'bright'],
    })
    specphot = Table({'FLUX': np.ones(5)})
    return meta, specphot


@pytest.fixture(scope='module')
def tile_tables():
    """Small metadata + specphot Tables with tile-coadd columns."""
    meta = Table({
        'TARGETID': np.array([10, 20, 30, 40]),
        'TILEID':   np.array([1000, 1000, 2000, 2000]),
        'NIGHT':    np.array([20230101, 20230101, 20230202, 20230303]),
    })
    specphot = Table({'FLUX': np.ones(4)})
    return meta, specphot


@pytest.fixture(scope='module')
def healpix_row():
    """Single astropy Row with healpix columns for get_qa_filename tests."""
    t = Table({
        'SURVEY':   ['main'],
        'PROGRAM':  ['bright'],
        'HEALPIX':  [42],
        'TARGETID': [999],
    })
    return t[0]


@pytest.fixture(scope='module')
def tile_row():
    """Single astropy Row with tile/night columns."""
    t = Table({
        'TILEID':   [1234],
        'NIGHT':    [20230101],
        'EXPID':    [56789],
        'TARGETID': [888],
    })
    return t[0]


@pytest.fixture(scope='module')
def stacked_row():
    """Single astropy Row with STACKID column."""
    t = Table({'STACKID': [7]})
    return t[0]


# ── select: healpix mode ──────────────────────────────────────────────────────

class TestSelectHealpix:

    def _call(self, meta, specphot, **kwargs):
        from fastspecfit.io import select
        return select(meta, specphot, coadd_type='healpix', **kwargs)

    def test_no_filter_keeps_all(self, healpix_tables):
        meta, sp = healpix_tables
        out_meta, out_sp, out_ff = self._call(meta, sp)
        assert len(out_meta) == len(meta)

    def test_single_healpixel_filter(self, healpix_tables):
        meta, sp = healpix_tables
        out_meta, out_sp, _ = self._call(meta, sp, healpixels=['100'])
        assert len(out_meta) == 2
        assert np.all(out_meta['HEALPIX'] == 100)

    def test_multiple_healpixels(self, healpix_tables):
        meta, sp = healpix_tables
        out_meta, out_sp, _ = self._call(meta, sp, healpixels=['100', '300'])
        assert len(out_meta) == 3
        assert set(out_meta['HEALPIX'].tolist()) == {100, 300}

    def test_nonexistent_healpixel_returns_empty(self, healpix_tables):
        meta, sp = healpix_tables
        out_meta, out_sp, _ = self._call(meta, sp, healpixels=['9999'])
        assert len(out_meta) == 0

    def test_return_index_gives_array(self, healpix_tables):
        meta, sp = healpix_tables
        idx = self._call(meta, sp, healpixels=['200'], return_index=True)
        assert isinstance(idx, np.ndarray)
        assert np.all(meta['HEALPIX'][idx] == 200)

    def test_fastfit_filtered_in_sync(self, healpix_tables):
        meta, sp = healpix_tables
        fastfit = Table({'VAL': np.arange(len(meta))})
        out_meta, out_sp, out_ff = self._call(meta, sp, fastfit=fastfit,
                                              healpixels=['100'])
        assert len(out_ff) == len(out_meta)

    def test_specphot_filtered_in_sync(self, healpix_tables):
        meta, sp = healpix_tables
        out_meta, out_sp, _ = self._call(meta, sp, healpixels=['300'])
        assert len(out_sp) == len(out_meta)


# ── select: tile mode ─────────────────────────────────────────────────────────

class TestSelectTile:

    def _call(self, meta, specphot, **kwargs):
        from fastspecfit.io import select
        return select(meta, specphot, coadd_type='cumulative', **kwargs)

    def test_no_filter_keeps_all(self, tile_tables):
        meta, sp = tile_tables
        out_meta, out_sp, _ = self._call(meta, sp)
        assert len(out_meta) == len(meta)

    def test_tile_filter_only(self, tile_tables):
        meta, sp = tile_tables
        out_meta, _, _ = self._call(meta, sp, tiles=['1000'])
        assert len(out_meta) == 2
        assert np.all(out_meta['TILEID'] == 1000)

    def test_night_filter_only(self, tile_tables):
        meta, sp = tile_tables
        out_meta, _, _ = self._call(meta, sp, nights=['20230101'])
        assert len(out_meta) == 2
        assert np.all(out_meta['NIGHT'] == 20230101)

    def test_tile_and_night_filter(self, tile_tables):
        meta, sp = tile_tables
        out_meta, _, _ = self._call(meta, sp, tiles=['1000'], nights=['20230101'])
        assert len(out_meta) == 2

    def test_tile_and_night_no_intersection(self, tile_tables):
        """Tile 2000 + night 20230101 don't overlap → empty result."""
        meta, sp = tile_tables
        out_meta, _, _ = self._call(meta, sp, tiles=['2000'], nights=['20230101'])
        assert len(out_meta) == 0

    def test_return_index_tile_filter(self, tile_tables):
        meta, sp = tile_tables
        idx = self._call(meta, sp, tiles=['2000'], return_index=True)
        assert np.all(meta['TILEID'][idx] == 2000)


# ── get_qa_filename ───────────────────────────────────────────────────────────

class TestGetQaFilename:

    def _call(self, metadata, coadd_type, **kwargs):
        from fastspecfit.io import get_qa_filename
        return get_qa_filename(metadata, coadd_type, **kwargs)

    def test_healpix_default_prefix(self, healpix_row):
        fn = self._call(healpix_row, 'healpix')
        assert fn.endswith('.png')
        assert 'fastspec' in fn
        assert 'main' in fn
        assert 'bright' in fn
        assert '42' in fn
        assert '999' in fn

    def test_custom_coadd_type_same_as_healpix(self, healpix_row):
        fn_hp  = self._call(healpix_row, 'healpix')
        fn_cus = self._call(healpix_row, 'custom')
        assert fn_hp == fn_cus

    def test_fastphot_prefix(self, healpix_row):
        fn = self._call(healpix_row, 'healpix', fastphot=True)
        assert fn.startswith('./fastphot')

    def test_custom_outprefix(self, healpix_row):
        fn = self._call(healpix_row, 'healpix', outprefix='myprefix')
        assert 'myprefix' in fn

    def test_custom_outdir(self, healpix_row, tmp_path):
        fn = self._call(healpix_row, 'healpix', outdir=str(tmp_path))
        assert fn.startswith(str(tmp_path))

    def test_cumulative(self, tile_row):
        fn = self._call(tile_row, 'cumulative')
        assert 'cumulative' in fn
        assert '1234' in fn
        assert '888' in fn

    def test_pernight(self, tile_row):
        fn = self._call(tile_row, 'pernight')
        assert '20230101' in fn
        assert '888' in fn

    def test_perexp(self, tile_row):
        fn = self._call(tile_row, 'perexp')
        assert '56789' in fn
        assert '888' in fn

    def test_stacked(self, stacked_row):
        fn = self._call(stacked_row, 'stacked')
        assert 'stacked' in fn
        assert '7' in fn

    def test_unknown_coadd_type_raises(self, healpix_row):
        with pytest.raises(ValueError):
            self._call(healpix_row, 'bogus_type')

    def test_table_input_returns_list(self, healpix_row):
        """Passing a multi-row Table returns a list of filenames."""
        t = Table({
            'SURVEY':   ['main', 'sv3'],
            'PROGRAM':  ['dark', 'bright'],
            'HEALPIX':  [10, 20],
            'TARGETID': [1, 2],
        })
        result = self._call(t, 'healpix')
        assert isinstance(result, list)
        assert len(result) == 2
