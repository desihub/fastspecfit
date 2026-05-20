"""
fastspecfit.test.conftest
=========================

"""
import os
import pathlib
import pytest
from urllib.request import urlretrieve

# Use the non-interactive Agg backend so QA tests can render figures in
# headless CI environments without a display.
os.environ.setdefault('MPLBACKEND', 'Agg')


@pytest.fixture(scope='session')
def template_version():
    yield '2.0.0'


@pytest.fixture(scope='session')
def outdir(tmp_path_factory):
    outdir = tmp_path_factory.mktemp('data')
    yield outdir


@pytest.fixture(scope='session')
def templatedir(outdir, template_version):
    cache_dir = os.environ.get('FTEMPLATES_CACHE_DIR')
    if cache_dir:
        templatedir = pathlib.Path(cache_dir) / template_version
        templatedir.mkdir(parents=True, exist_ok=True)
    else:
        templatedir = outdir / template_version
        templatedir.mkdir()
    yield templatedir


@pytest.fixture(scope='session')
def templates(templatedir, template_version):
    templates_file = f'ftemplates-chabrier-{template_version}.fits'
    templates = os.path.join(templatedir, templates_file)

    url = f"https://data.desi.lbl.gov/public/external/templates/fastspecfit/2.0.0/{templates_file}"
    if not os.path.isfile(templates):
        urlretrieve(url, templates)
    yield templates

    # Skip cleanup when using a persistent cache directory.
    if not os.environ.get('FTEMPLATES_CACHE_DIR') and os.path.isfile(templates):
        os.remove(templates)


@pytest.fixture(scope='session')
def filenames(outdir):
    from importlib import resources
    redux_dir    = resources.files('fastspecfit').joinpath('test/data')
    specproddir  = resources.files('fastspecfit').joinpath('test/data')
    mapdir       = resources.files('fastspecfit').joinpath('test/data')
    fphotodir    = resources.files('fastspecfit').joinpath('test/data')
    redrockfile  = resources.files('fastspecfit').joinpath('test/data/redrock-4-80613-thru20210324.fits')
    stackfile    = resources.files('fastspecfit').joinpath('test/data/stack-LRG.fits')
    yield {
        'redux_dir':        redux_dir,
        'specproddir':      specproddir,
        'mapdir':           mapdir,
        'fphotodir':        fphotodir,
        'redrockfile':      redrockfile,
        'stackfile':        stackfile,
        'fastphot_outfile': os.path.join(outdir, 'fastphot.fits'),
        'fastspec_outfile': os.path.join(outdir, 'fastspec.fits'),
        'stackfit_outfile': os.path.join(outdir, 'stackfit.fits'),
    }


@pytest.fixture(scope='session')
def fastphot_output(filenames, templates):
    from fastspecfit.fastspecfit import fastphot, parse
    outfile = filenames['fastphot_outfile']
    cmd = (f'fastphot {filenames["redrockfile"]} -o {outfile} '
           f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} '
           f'--redux_dir {filenames["redux_dir"]} '
           f'--specproddir {filenames["specproddir"]} --templates {templates}')
    fastphot(args=parse(options=cmd.split()[1:]))
    yield outfile


@pytest.fixture(scope='session')
def fastspec_output(filenames, templates):
    from fastspecfit.fastspecfit import fastspec, parse
    outfile = filenames['fastspec_outfile']
    cmd = (f'fastspec {filenames["redrockfile"]} -o {outfile} '
           f'--redux_dir {filenames["redux_dir"]} '
           f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} '
           f'--specproddir {filenames["specproddir"]} --templates {templates}')
    fastspec(args=parse(options=cmd.split()[1:]))
    yield outfile


@pytest.fixture(scope='session')
def stackfit_output(filenames, templates):
    from fastspecfit.fastspecfit import stackfit, parse
    outfile = filenames['stackfit_outfile']
    cmd = f'stackfit {filenames["stackfile"]} -o {outfile} --templates {templates}'
    stackfit(args=parse(options=cmd.split()[1:]))
    yield outfile
