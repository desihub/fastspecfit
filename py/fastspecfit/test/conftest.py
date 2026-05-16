"""
fastspecfit.test.conftest
=========================

"""
import os
import pathlib
import pytest
from urllib.request import urlretrieve

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
