"""
fastspecfit.test.test_templates
===============================

"""
import os
import pytest
import numpy as np
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

    # Optional cleanup
    if os.path.isfile(templates):
        os.remove(templates)


def test_templates_nofilename(outdir, template_version, templates):
    from fastspecfit.templates import Templates
    os.environ['FTEMPLATES_DIR'] = str(outdir)
    temp = Templates()
    assert(os.path.isfile(temp.file))


def test_templates(templates):
    from fastspecfit.templates import Templates
    temp = Templates(template_file=templates)
    assert(os.path.isfile(temp.file))
    assert(hasattr(temp, 'pixkms_bounds'))
