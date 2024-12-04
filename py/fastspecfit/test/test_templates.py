"""
fastspecfit.test.test_templates
===============================

"""
import os


# See conftest.py for fixtures used by multiple unit tests.

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
