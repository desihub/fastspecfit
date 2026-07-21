"""
fastspecfit.test.test_fastspecfit
==================================

"""
import os
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Template tests
# ---------------------------------------------------------------------------

@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_template_dt_column(templates):
    """Templates >=2.1.0 must carry a dt column for 100-Myr SFR averaging."""
    from fastspecfit.templates import Templates
    T = Templates(template_file=templates)
    assert 'dt' in T.info.colnames, 'dt column missing'
    assert np.all(T.info['dt'] > 0), 'dt values must be positive'


# ---------------------------------------------------------------------------
# Fitting integration tests
# ---------------------------------------------------------------------------

@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_fastphot(fastphot_output):
    """Test that fastphot produces a valid output file."""
    import fitsio
    assert os.path.exists(fastphot_output)
    fits = fitsio.FITS(fastphot_output)
    for hdu in fits:
        if hdu.has_data():
            assert(hdu.get_extname() in ['METADATA', 'SPECPHOT'])


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_stackfit(stackfit_output):
    """Test that stackfit produces a valid output file."""
    import fitsio
    assert os.path.exists(stackfit_output)
    fits = fitsio.FITS(stackfit_output)
    for hdu in fits:
        if hdu.has_data():
            assert(hdu.get_extname() in ['METADATA', 'SPECPHOT', 'FASTSPEC', 'MODELS'])


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_fastspec(fastspec_output):
    """Test that fastspec produces a valid output file."""
    import fitsio
    assert os.path.exists(fastspec_output)
    fits = fitsio.FITS(fastspec_output)
    for hdu in fits:
        if hdu.has_data():
            assert(hdu.get_extname() in ['METADATA', 'SPECPHOT', 'FASTSPEC', 'MODELS'])


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_fastspec_hii_morelines(fastspec_hii_output):
    """The HII-region maximal line list (63 lines) exceeds the FITS 999-column
    limit, triggering the FASTSPEC/MORELINES HDU split (see issue #236)."""
    import fitsio
    from fastspecfit.io import read_fastspecfit, MAX_FASTSPEC_COLUMNS, MORELINES_IDCOLS

    fits = fitsio.FITS(fastspec_hii_output)
    extnames = [fits[i].get_extname() for i in range(1, len(fits))]
    assert 'MORELINES' in extnames

    fastspec_cols = fits['FASTSPEC'].get_colnames()
    morelines_cols = fits['MORELINES'].get_colnames()
    assert len(fastspec_cols) <= MAX_FASTSPEC_COLUMNS
    assert 'TARGETID' in morelines_cols

    # a known strong line stays in FASTSPEC; a known weak/auroral line moves
    # to MORELINES
    assert 'OIII_5007_AMP' in fastspec_cols
    assert 'HEI_4471_AMP' in morelines_cols
    assert 'HEI_4471_AMP' not in fastspec_cols

    # reading transparently merges FASTSPEC+MORELINES back into one table
    _, _, fastfit, _, _ = read_fastspecfit(fastspec_hii_output)
    idcols_in_morelines = [c for c in morelines_cols if c in MORELINES_IDCOLS]
    expected_ncols = len(fastspec_cols) + len(morelines_cols) - len(idcols_in_morelines)
    assert len(fastfit.colnames) == expected_ncols
    assert 'HEI_4471_AMP' in fastfit.colnames
    assert 'OIII_5007_AMP' in fastfit.colnames
    assert np.all(np.isfinite(fastfit['HEI_4471_AMP']))


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_sfr_values(fastspec_output, fastphot_output):
    """SFR in SPECPHOT must be finite and non-negative for all objects."""
    import fitsio
    for outfile in [fastspec_output, fastphot_output]:
        data = fitsio.read(outfile, ext='SPECPHOT')
        assert 'SFR' in data.dtype.names
        assert np.all(data['SFR'] >= 0), f'negative SFR in {outfile}'
        assert np.all(np.isfinite(data['SFR'])), f'non-finite SFR in {outfile}'


# ---------------------------------------------------------------------------
# QA integration tests
# ---------------------------------------------------------------------------

@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_qa_fastphot(fastphot_output, filenames, templates, outdir):
    """Test that fastqa runs and produces output for a fastphot file."""
    from pathlib import Path
    from fastspecfit.qa import fastqa, parse as qa_parse
    qa_outdir = str(outdir / 'qa_fastphot')
    cmd = (f'fastqa {fastphot_output} '
           f'--redrockfiles {filenames["redrockfile"]} '
           f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} '
           f'--templates {templates} --outdir {qa_outdir} --overwrite')
    fastqa(args=qa_parse(options=cmd.split()[1:]))
    assert len(list(Path(qa_outdir).glob('*.png'))) > 0


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_qa_constraintsfile(fastphot_output, filenames, templates, outdir):
    """Test that fastqa accepts a custom --constraintsfile argument."""
    from pathlib import Path
    from importlib import resources
    from fastspecfit.qa import fastqa, parse as qa_parse
    constraintsfile = str(resources.files('fastspecfit').joinpath('data/emline-constraints.yaml'))
    qa_outdir = str(outdir / 'qa_constraintsfile')
    cmd = (f'fastqa {fastphot_output} '
           f'--redrockfiles {filenames["redrockfile"]} '
           f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} '
           f'--templates {templates} --constraintsfile {constraintsfile} '
           f'--outdir {qa_outdir} --overwrite')
    fastqa(args=qa_parse(options=cmd.split()[1:]))
    assert len(list(Path(qa_outdir).glob('*.png'))) > 0


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_qa_fastspec(fastspec_output, filenames, templates, outdir):
    """Test that fastqa runs and produces output for a fastspec file."""
    from pathlib import Path
    from fastspecfit.qa import fastqa, parse as qa_parse
    qa_outdir = str(outdir / 'qa_fastspec')
    cmd = (f'fastqa {fastspec_output} '
           f'--redrockfiles {filenames["redrockfile"]} '
           f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} '
           f'--templates {templates} --outdir {qa_outdir} --overwrite')
    fastqa(args=qa_parse(options=cmd.split()[1:]))
    assert len(list(Path(qa_outdir).glob('*.png'))) > 0


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_fastspec_suprime(fastspec_suprime_output, suprime_filenames):
    """Test fastspec with an external photometric catalog (Suprime/HSC mode)."""
    import fitsio
    assert os.path.exists(fastspec_suprime_output)
    fits = fitsio.FITS(fastspec_suprime_output)
    for hdu in fits:
        if hdu.has_data():
            assert hdu.get_extname() in ['METADATA', 'SPECPHOT', 'FASTSPEC', 'MODELS']
    meta = fitsio.read(fastspec_suprime_output, ext='METADATA')
    assert meta['TARGETID'][0] == suprime_filenames['targetid']


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_qa_stackfit(stackfit_output, filenames, templates, outdir):
    """Test that fastqa runs and produces output for a stackfit file."""
    from pathlib import Path
    from fastspecfit.qa import fastqa, parse as qa_parse
    qa_outdir = str(outdir / 'qa_stackfit')
    cmd = (f'fastqa {stackfit_output} '
           f'--redrockfiles {filenames["stackfile"]} '
           f'--templates {templates} --outdir {qa_outdir} --overwrite')
    fastqa(args=qa_parse(options=cmd.split()[1:]))
    assert len(list(Path(qa_outdir).glob('*.png'))) > 0
