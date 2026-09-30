"""
fastspecfit.test.test_photometry
================================

"""
import os
import pytest
import numpy as np
import fitsio
from astropy.table import Table


@pytest.fixture
def phot():
    from fastspecfit.photometry import Photometry
    phot = Photometry()
    assert(os.path.isfile(phot.fphotofile))
    yield phot


@pytest.fixture
def phot_stacked():
    from fastspecfit.photometry import Photometry
    phot_stacked = Photometry(fitstack=True)
    yield phot_stacked


@pytest.fixture
def phot_ignore():
    from fastspecfit.photometry import Photometry
    phot_ignore = Photometry(ignore_photometry=True)
    yield phot_ignore


@pytest.fixture
def data(phot):
    from fastspecfit.cosmo import TabulatedDESI
    cosmo = TabulatedDESI()

    zwave = np.logspace(np.log10(50.), np.log10(40e4), 10000)
    zflux = np.ones_like(zwave)
    zivar = np.zeros_like(zwave) + 100.

    redshift = 0.2
    dmod = cosmo.distance_modulus(redshift)
    photsys = 'S'

    absmag_filters = np.array(['decam2014-g', 'decam2014-r', 'decam2014-z',
                               'bessell-U', 'bessell-B', 'bessell-V', 'twomass-J',
                               'sdss2010-u', 'sdss2010-g', 'sdss2010-r',
                               'sdss2010-i', 'sdss2010-z', 'wise2010-W1'])

    nanomaggies = np.ones(len(phot.filters[photsys]))
    nanomaggies_ivar = np.ones_like(nanomaggies) + 100.

    data = {'zwave': zwave, 'zflux': zflux, 'zivar': zivar,
            'redshift': redshift, 'dmod': dmod, 'photsys': photsys,
            'absmag_filters': absmag_filters, 'nanomaggies': nanomaggies,
            'nanomaggies_ivar': nanomaggies_ivar}

    yield data


@pytest.fixture
def expected(phot, data):

    synth_absmag = np.array([-16.19285649, -16.81623895, -17.58668895, -17.0631772 ,
                             -17.49934411, -17.98466151, -19.74184508, -16.85210932,
                             -17.43688202, -18.02971177, -18.45005832, -18.82373475,
                             -21.71275126])
    rest_nanomaggies = np.array([0.6403109 ,  1.13695657,  2.31165958, 1.42731286,
                                 2.13297587,  3.33513103, 16.82621093, 1.17514654,
                                 2.01372938,  3.47642615,  5.12001964, 7.22340889,
                                 103.35911934])
    kcorr = np.array([-1.30718301e+00, -6.83789219e-01,  8.66694134e-02, -4.36832738e-01,
                      -6.69176268e-04, -1.38803532e-01,  8.47950031e-01, -6.47904968e-01,
                      -6.31308290e-02, -9.37527202e-02, -4.43827450e-01, -7.01527417e-02,
                      -9.00627082e-03])
    dn4000 = 1.078368181
    dn4000_ivar = 4298.47424

    expected = {'synth_absmag': synth_absmag, 'rest_nanomaggies': rest_nanomaggies,
                'kcorr': kcorr, 'dn4000': dn4000, 'dn4000_ivar': dn4000_ivar}

    yield expected


#def test_photometry(self):
#    from fastspecfit.photometry import Photometry
#    phot_stacked = Photometry(fitstack=True)
#    phot_ignore = Photometry(ignore_photometry=True)
#
#    self.assertTrue(os.path.isfile(self.phot.fphotofile))
#    self.assertTrue(np.all(self.phot.absmag_filters.names == self.absmag_filters))


def test_kcorr_and_absmag(phot, data, expected):
    synth_absmag, synth_maggies_rest = phot.synth_absmag(
        data['redshift'], data['dmod'], data['zwave'],
        data['zflux'])
    assert(np.all(np.isclose(synth_absmag, expected['synth_absmag'], 1e-4)))
    assert(np.all(np.isclose(1e9 * synth_maggies_rest, expected['rest_nanomaggies'], 1e-4)))

    kcorr, absmag, ivarabsmag, synth_mmaggies_obs = phot.kcorr_and_absmag(
        data['nanomaggies'], data['nanomaggies_ivar'], data['redshift'],
        data['dmod'], data['photsys'], data['zwave'], data['zflux'],
        synth_absmag, synth_maggies_rest)
    assert(np.all(np.isclose(kcorr, expected['kcorr'], 1e-4)))

    # test redshift=0.
    synth_absmag, synth_maggies_rest = phot.synth_absmag(
        0., data['dmod'], data['zwave'], data['zflux'])
    assert(np.all(synth_absmag == 0.))
    assert(np.all(synth_maggies_rest == 0.))

    kcorr, absmag, ivarabsmag, synth_maggies_obs = phot.kcorr_and_absmag(
        data['nanomaggies'], data['nanomaggies_ivar'], 0.,
        data['dmod'], data['photsys'], data['zwave'],
        data['zflux'], synth_absmag, synth_maggies_rest)
    assert(np.all(kcorr == 0.))
    assert(np.all(absmag == 0.))
    assert(np.all(ivarabsmag == 0.))
    assert(np.all(synth_maggies_obs == 0.))


def test_dn4000(phot, data, expected):
    dn4000_1, dn4000_ivar_1 = phot.get_dn4000(
        data['zwave'], data['zflux'], redshift=data['redshift'],
        rest=False)
    assert(np.isclose(dn4000_1, expected['dn4000'], 1e-4))
    assert(np.isclose(dn4000_ivar_1, 0.))

    dn4000_2, dn4000_ivar_2 = phot.get_dn4000(
        data['zwave'] / (1. + data['redshift']), data['zflux'] * (1. + data['redshift']),
        flam_ivar=data['zivar'] / (1. + data['redshift'])**2, redshift=None, rest=True)
    assert(np.isclose(dn4000_2, expected['dn4000'], 1e-4))
    assert(np.isclose(dn4000_ivar_2, expected['dn4000_ivar'], 1e-4))

    I = data['zwave'] / (1. + data['redshift']) > 3850.
    dn4000_3, dn4000_ivar_3 = phot.get_dn4000(
        data['zwave'][I] / (1. + data['redshift']), data['zflux'][I] * (1. + data['redshift']),
        flam_ivar=data['zivar'][I] / (1. + data['redshift'])**2, redshift=None, rest=True)
    assert(dn4000_3 == 0.)
    assert(dn4000_ivar_3 == 0.)


# ---------------------------------------------------------------------------
# gather_tractorphot
# ---------------------------------------------------------------------------

# Two synthetic Tractor sources sharing one brick, used across the tests
# below. OBJ_A and OBJ_B have different RELEASE values so that a mismatch
# between the two is unambiguous.
BRICK = '1501p020'
OBJ_A = dict(OBJID=100, BRICKID=555, RELEASE=9010, RA=150.00, DEC=2.00,
            BRICK_PRIMARY=True, FLUX_R=10.0)
OBJ_B = dict(OBJID=200, BRICKID=555, RELEASE=9011, RA=150.01, DEC=2.01,
            BRICK_PRIMARY=True, FLUX_R=20.0)


def _write_fake_tractor_brick(tractorfile, rows):
    """Write a minimal but structurally-valid Tractor brick FITS file."""
    from fastspecfit.photometry import tractorphot_datamodel

    datamodel = tractorphot_datamodel(datarelease='dr9')
    tractor = Table(np.hstack(np.repeat(datamodel, len(rows))))
    for irow, row in enumerate(rows):
        tractor['BRICKNAME'][irow] = BRICK
        for key, val in row.items():
            tractor[key][irow] = val

    os.makedirs(os.path.dirname(tractorfile), exist_ok=True)
    fitsio.write(tractorfile, tractor.as_array(), clobber=True)


@pytest.fixture
def legacysurveydir(tmp_path):
    # Naming the directory 'dr9' lets gather_tractorphot infer the data
    # release without a warning.
    legacysurveydir = tmp_path / 'dr9'
    tractorfile = legacysurveydir / 'south' / 'tractor' / BRICK[:3] / f'tractor-{BRICK}.fits'
    _write_fake_tractor_brick(str(tractorfile), [OBJ_A, OBJ_B])
    return str(legacysurveydir)


def test_gather_tractorphot_direct_match(legacysurveydir):
    from desitarget.targets import encode_targetid
    from fastspecfit.photometry import gather_tractorphot

    # TARGETID is consistent with the (correct) BRICKID/BRICK_OBJID/RELEASE
    # metadata, so no RELEASE mismatch and no repair pass.
    targetid = encode_targetid(objid=OBJ_A['OBJID'], brickid=OBJ_A['BRICKID'], release=OBJ_A['RELEASE'])
    cat = Table({
        'TARGETID': [targetid], 'TARGET_RA': [OBJ_A['RA']], 'TARGET_DEC': [OBJ_A['DEC']],
        'BRICKNAME': [BRICK], 'BRICKID': [OBJ_A['BRICKID']], 'BRICK_OBJID': [OBJ_A['OBJID']],
        'RELEASE': [OBJ_A['RELEASE']], 'PHOTSYS': ['S'],
    })
    out = gather_tractorphot(cat, legacysurveydir=legacysurveydir)
    assert out['OBJID'][0] == OBJ_A['OBJID']
    assert out['FLUX_R'][0] == OBJ_A['FLUX_R']
    assert out['LS_ID'][0] != 0


def test_gather_tractorphot_positional_match(legacysurveydir):
    from desitarget.targets import encode_targetid
    from fastspecfit.photometry import gather_tractorphot

    # No BRICKID/BRICK_OBJID supplied, so matching falls back to RA/Dec;
    # TARGETID encodes the true (matching) RELEASE so no repair is triggered.
    targetid = encode_targetid(objid=OBJ_B['OBJID'], brickid=OBJ_B['BRICKID'], release=OBJ_B['RELEASE'])
    cat = Table({
        'TARGETID': [targetid], 'TARGET_RA': [OBJ_B['RA']], 'TARGET_DEC': [OBJ_B['DEC']],
    })
    out = gather_tractorphot(cat, legacysurveydir=legacysurveydir)
    assert out['OBJID'][0] == OBJ_B['OBJID']
    assert out['FLUX_R'][0] == OBJ_B['FLUX_R']


def test_gather_tractorphot_repair_path(legacysurveydir):
    """Regression test: a RELEASE mismatch must trigger a positional
    re-match whose result is merged back into the output (out[bug] = bugout).
    """
    from desitarget.targets import encode_targetid
    from fastspecfit.photometry import gather_tractorphot

    # TARGETID encodes the true match (OBJ_B), but the catalog's BRICKID/
    # BRICK_OBJID/RELEASE columns are stale and point at OBJ_A instead. The
    # direct-match pass will therefore land on the wrong source (OBJ_A), the
    # RELEASE check will catch the discrepancy, and the repair pass should
    # recover OBJ_B via positional matching on the (correct) RA/Dec.
    targetid = encode_targetid(objid=OBJ_B['OBJID'], brickid=OBJ_B['BRICKID'], release=OBJ_B['RELEASE'])
    cat = Table({
        'TARGETID': [targetid], 'TARGET_RA': [OBJ_B['RA']], 'TARGET_DEC': [OBJ_B['DEC']],
        'BRICKNAME': [BRICK], 'BRICKID': [OBJ_A['BRICKID']], 'BRICK_OBJID': [OBJ_A['OBJID']],
        'RELEASE': [OBJ_A['RELEASE']], 'PHOTSYS': ['S'],
    })
    out = gather_tractorphot(cat, legacysurveydir=legacysurveydir)
    assert out['OBJID'][0] == OBJ_B['OBJID']
    assert out['RELEASE'][0] == OBJ_B['RELEASE']
    assert out['FLUX_R'][0] == OBJ_B['FLUX_R']


def test_gather_tractorphot_columns_subset(legacysurveydir):
    from desitarget.targets import encode_targetid
    from fastspecfit.photometry import gather_tractorphot

    targetid = encode_targetid(objid=OBJ_A['OBJID'], brickid=OBJ_A['BRICKID'], release=OBJ_A['RELEASE'])
    cat = Table({
        'TARGETID': [targetid], 'TARGET_RA': [OBJ_A['RA']], 'TARGET_DEC': [OBJ_A['DEC']],
        'BRICKNAME': [BRICK], 'BRICKID': [OBJ_A['BRICKID']], 'BRICK_OBJID': [OBJ_A['OBJID']],
        'RELEASE': [OBJ_A['RELEASE']], 'PHOTSYS': ['S'],
    })
    out = gather_tractorphot(cat, legacysurveydir=legacysurveydir, columns=['TARGETID', 'FLUX_R'])
    assert out.colnames == ['TARGETID', 'FLUX_R']


def test_gather_tractorphot_missing_column():
    from fastspecfit.photometry import gather_tractorphot

    cat = Table({'TARGETID': [1]})  # missing TARGET_RA, TARGET_DEC
    with pytest.raises(ValueError):
        gather_tractorphot(cat)


def test_gather_tractorphot_bad_legacysurveydir(tmp_path):
    from fastspecfit.photometry import gather_tractorphot

    cat = Table({'TARGETID': [1], 'TARGET_RA': [150.0], 'TARGET_DEC': [2.0]})
    with pytest.raises(IOError):
        gather_tractorphot(cat, legacysurveydir=str(tmp_path / 'nonexistent'))


def test_gather_tractorphot_bad_photsys(legacysurveydir):
    from fastspecfit.photometry import gather_tractorphot

    cat = Table({
        'TARGETID': [1], 'TARGET_RA': [OBJ_A['RA']], 'TARGET_DEC': [OBJ_A['DEC']],
        'BRICKNAME': [BRICK], 'BRICKID': [OBJ_A['BRICKID']], 'BRICK_OBJID': [OBJ_A['OBJID']],
        'RELEASE': [OBJ_A['RELEASE']], 'PHOTSYS': ['X'],
    })
    with pytest.raises(ValueError):
        gather_tractorphot(cat, legacysurveydir=legacysurveydir)

