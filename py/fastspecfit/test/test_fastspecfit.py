"""
fastspecfit.test.test_continuum
===============================

Test fastspecfit.fastspecfit.fastspec

# ToOs -- no broadband photometry
fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/19/20210501/redrock-4-19-thru20210501.fits -o fastspec.fits --targetids 43977515642913220,43977515642913243

# secondary targets with good targetphot photometry
fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/80895/20210403/redrock-7-80895-thru20210403.fits -o fastspec.fits --targetids 39632936277902799,39632931181824699,39632931173436143,39632936273709365,39632936273708359

# secondary targets with no good targetphot photometry
fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/80856/20210318/redrock-9-80856-thru20210318.fits -o fastspec.fits --targetids 6432023904256,6448025174016

"""
import pdb
import unittest, os, shutil, tempfile, subprocess
import numpy as np
from unittest.mock import patch, call
from pkg_resources import resource_filename

class TestFastspec(unittest.TestCase):
    """Test fastspecfit.fastspecfit.fastspec"""
    @classmethod
    def setUpClass(cls):
        os.environ['DESI_ROOT'] = resource_filename('fastspecfit.test', 'data')
        cls.mapdir = resource_filename('fastspecfit.test', 'data')
        cls.dr9dir = resource_filename('fastspecfit.test', 'data')
        cls.redrockfile = resource_filename('fastspecfit.test', 'data/redrock-4-80613-thru20210324.fits')

        cls.outdir = tempfile.mkdtemp()
        cls.templates = os.path.join(cls.outdir, 'ftemplates-chabrier-1.0.0.fits')
        cmd = 'wget -O {} https://portal.nersc.gov/project/cosmo/temp/ioannis/templates/ftemplates-chabrier-1.0.0.fits'.format(
            cls.templates)
        err = subprocess.call(cmd.split())
        cls.cwd = os.getcwd()
        cls.fastspec_outfile = os.path.join(cls.outdir, 'fastspec.fits')
        cls.fastphot_outfile = os.path.join(cls.outdir, 'fastphot.fits')

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_ContinuumTools(self):
        """Test the ContinuumTools class."""
        from fastspecfit.continuum import ContinuumTools
        CTools = ContinuumTools(mapdir=self.mapdir, templates=self.templates)

        # expected attributes
        self.assertTrue(CTools.imf in ['salpeter', 'chabrier', 'kroupa'])

    def test_fastphot(self):
        """Test fastphot."""
        import fitsio
        from fastspecfit.fastspecfit import fastphot, parse

        cmd = 'fastphot {} -o {} --mapdir {} --dr9dir {} --templates {}'.format(
            self.redrockfile, self.fastphot_outfile, self.mapdir, self.dr9dir, self.templates)
        args = parse(options=cmd.split()[1:])
        fastphot(args=args)

        self.assertTrue(os.path.exists(self.fastphot_outfile))

        fits = fitsio.FITS(self.fastphot_outfile)
        for hdu in fits:
            if hdu.has_data(): # skip zeroth extension
                self.assertTrue(hdu.get_extname() in ['METADATA', 'FASTPHOT'])

    def test_fastspec(self):
        """Test fastspec."""
        import fitsio
        from fastspecfit.fastspecfit import fastspec, parse
    
        cmd = 'fastspec {} -o {} --mapdir {} --dr9dir {} --templates {}'.format(
            self.redrockfile, self.fastspec_outfile, self.mapdir, self.dr9dir, self.templates)
        args = parse(options=cmd.split()[1:])
        fastspec(args=args)
    
        self.assertTrue(os.path.exists(self.fastspec_outfile))
    
        fits = fitsio.FITS(self.fastspec_outfile)
        for hdu in fits:
            if hdu.has_data(): # skip zeroth extension
                self.assertTrue(hdu.get_extname() in ['METADATA', 'FASTSPEC', 'MODELS'])

if __name__ == '__main__':
    unittest.main()
