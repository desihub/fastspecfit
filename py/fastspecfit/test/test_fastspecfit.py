"""
fastspecfit.test.test_fastspecfit
=================================

Test fastspecfit.fastspecfit.fastspec

# ToOs -- no broadband photometry
fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/19/20210501/redrock-4-19-thru20210501.fits -o fastspec.fits --targetids 43977515642913220,43977515642913243

# secondary targets with good targetphot photometry
fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/80895/20210403/redrock-7-80895-thru20210403.fits -o fastspec.fits --targetids 39632936277902799,39632931181824699,39632931173436143,39632936273709365,39632936273708359

# secondary targets with no good targetphot photometry
fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/80856/20210318/redrock-9-80856-thru20210318.fits -o fastspec.fits --targetids 6432023904256,6448025174016

"""
import os
import tempfile
import pytest
import numpy as np
from urllib.request import urlretrieve


def setUpClass(cls):
    os.environ['DESI_SPECTRO_REDUX'] = str(resources.files('fastspecfit').joinpath('test/data'))

    #cls.templates = '/Users/ioannis/work/desi/users/ioannis/fastspecfit/templates/2.0.0/ftemplates-chabrier-2.0.0.fits'
    #cls.outdir = tempfile.mkdtemp()
    #cls.templates = os.path.join(cls.outdir, 'ftemplates-chabrier-2.0.0.fits')
    #if os.path.isfile(cls.templates):
    #    os.remove(cls.templates)
    ##url = "https://portal.nersc.gov/project/cosmo/temp/ioannis/tmp/ftemplates-chabrier-2.0.0.fits"
    #url = "https://data.desi.lbl.gov/public/external/templates/fastspecfit/2.0.0/ftemplates-chabrier-2.0.0.fits"
    #urlretrieve(url, cls.templates)

    cls.fastspec_outfile = os.path.join(cls.outdir, 'fastspec.fits')
    cls.fastphot_outfile = os.path.join(cls.outdir, 'fastphot.fits')


@pytest.fixture
def filenames(outdir):
    from importlib import resources

    specproddir = resources.files('fastspecfit').joinpath('test/data')
    mapdir = resources.files('fastspecfit').joinpath('test/data')
    fphotodir = resources.files('fastspecfit').joinpath('test/data')
    redrockfile = resources.files('fastspecfit').joinpath('test/data/redrock-4-80613-thru20210324.fits')
    fastspec_outfile = os.path.join(outdir, 'fastspec.fits')
    fastphot_outfile = os.path.join(outdir, 'fastphot.fits')

    filenames = {'specproddir': specproddir, 'mapdir': mapdir, 'fphotodir': fphotodir,
                 'redrockfile': redrockfile, 'fastspec_outfile': fastspec_outfile,
                 'fastphot_outfile': fastphot_outfile, }

    yield filenames


def test_fastphot(filenames, templates):
    """Test fastphot."""
    import fitsio
    from fastspecfit.fastspecfit import fastphot, parse

    outfile = filenames["fastphot_outfile"]

    cmd = f'fastphot {filenames["redrockfile"]} -o {outfile} ' + \
        f'--mapdir {filenames["mapdir"]} --fphotodir {filenames["fphotodir"]} ' + \
        f'--specproddir {filenames["specproddir"]} --templates {templates}'

    args = parse(options=cmd.split()[1:])
    fastphot(args=args)

    assert(os.path.exists(outfile))

    #fits = fitsio.FITS(outfile)
    #for hdu in fits:
    #    if hdu.has_data(): # skip zeroth extension
    #        assert(hdu.get_extname() in ['METADATA', 'SPECPHOT'])


#def test_fastspec(self):
#    """Test fastspec."""
#    import fitsio
#    from fastspecfit.fastspecfit import fastspec, parse
#
#    cmd = f'fastspec {self.redrockfile} -o {self.fastspec_outfile} --mapdir {self.mapdir} ' + \
#        f'--fphotodir {self.fphotodir} --specproddir {self.specproddir} --templates {self.templates}'
#    args = parse(options=cmd.split()[1:])
#    fastspec(args=args)
#
#    self.assertTrue(os.path.exists(self.fastspec_outfile))
#
#    fits = fitsio.FITS(self.fastspec_outfile)
#    for hdu in fits:
#        if hdu.has_data(): # skip zeroth extension
#            self.assertTrue(hdu.get_extname() in ['METADATA', 'SPECPHOT', 'FASTSPEC', 'MODELS'])
