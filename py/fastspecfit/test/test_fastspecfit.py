"""
fastspecfit.test.test_fastspecfit
=================================

"""
import os
import pytest


@pytest.fixture
def filenames(outdir):
    from importlib import resources

    redux_dir = resources.files('fastspecfit').joinpath('test/data')
    specproddir = resources.files('fastspecfit').joinpath('test/data')
    mapdir = resources.files('fastspecfit').joinpath('test/data')
    fphotodir = resources.files('fastspecfit').joinpath('test/data')
    redrockfile = resources.files('fastspecfit').joinpath('test/data/redrock-4-80613-thru20210324.fits')
    stackfile = resources.files('fastspecfit').joinpath('test/data/stack-LRG.fits')
    fastspec_outfile = os.path.join(outdir, 'fastspec.fits')
    fastphot_outfile = os.path.join(outdir, 'fastphot.fits')
    stackfit_outfile = os.path.join(outdir, 'stackfit.fits')

    filenames = {'redux_dir': redux_dir, 'specproddir': specproddir,
                 'mapdir': mapdir, 'fphotodir': fphotodir,
                 'redrockfile': redrockfile, 'stackfile': stackfile,
                 'fastspec_outfile': fastspec_outfile,
                 'fastphot_outfile': fastphot_outfile,
                 'stackfit_outfile': stackfit_outfile}

    yield filenames


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_fastphot(filenames, templates):
    """Test fastphot."""
    import fitsio
    from fastspecfit.fastspecfit import fastphot, parse

    outfile = filenames["fastphot_outfile"]

    cmd = f'fastphot {filenames["redrockfile"]} -o {outfile} ' + \
        f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} ' + \
        f'--redux_dir {filenames["redux_dir"]} ' + \
        f'--specproddir {filenames["specproddir"]} --templates {templates}'

    args = parse(options=cmd.split()[1:])
    fastphot(args=args)

    assert(os.path.exists(outfile))

    fits = fitsio.FITS(outfile)
    for hdu in fits:
        if hdu.has_data(): # skip zeroth extension
            assert(hdu.get_extname() in ['METADATA', 'SPECPHOT'])


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_stackfit(filenames, templates):
    """Test stackfit."""
    import fitsio
    from fastspecfit.fastspecfit import stackfit, parse

    outfile = filenames["stackfit_outfile"]

    cmd = f'stackfit {filenames["stackfile"]} -o {outfile} --templates {templates}'

    args = parse(options=cmd.split()[1:])
    stackfit(args=args)

    assert(os.path.exists(outfile))

    fits = fitsio.FITS(outfile)
    for hdu in fits:
        if hdu.has_data():
            assert(hdu.get_extname() in ['METADATA', 'SPECPHOT', 'FASTSPEC', 'MODELS'])


@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_fastspec(filenames, templates):
    """Test fastspec."""
    import fitsio
    from fastspecfit.fastspecfit import fastspec, parse

    outfile = filenames["fastspec_outfile"]

    cmd = f'fastspec {filenames["redrockfile"]} -o {outfile} ' + \
        f'--redux_dir {filenames["redux_dir"]} ' + \
        f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} ' + \
        f'--specproddir {filenames["specproddir"]} --templates {templates}'

    args = parse(options=cmd.split()[1:])
    fastspec(args=args)

    assert(os.path.exists(outfile))

    fits = fitsio.FITS(outfile)
    for hdu in fits:
        if hdu.has_data(): # skip zeroth extension
            assert(hdu.get_extname() in ['METADATA', 'SPECPHOT', 'FASTSPEC', 'MODELS'])
