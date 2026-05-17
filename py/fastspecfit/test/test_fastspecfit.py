"""
fastspecfit.test.test_fastspecfit
==================================

"""
import os
import pytest


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
