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
        cls.templates = os.path.join(resource_filename('fastspecfit.test', 'data'),
                                        'SSP_Padova_CKC14z_Kroupa_Z0.0190.fits')
        cls.mapdir = resource_filename('fastspecfit.test', 'data')
        cls.dr9dir = resource_filename('fastspecfit.test', 'data')
        cls.redrockfile = resource_filename('fastspecfit.test', 'data/redrock-4-80613-thru20210324.fits')
        cls.cwd = os.getcwd()
        cls.outdir = tempfile.mkdtemp()
        cls.fastspec_outfile = os.path.join(cls.outdir, 'fastspec.fits')
        cls.fastphot_outfile = os.path.join(cls.outdir, 'fastphot.fits')

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_ContinuumFit(self):
        """Test the ContinuumTools class."""
        from fastspecfit.continuum import ContinuumFit
        CFit = ContinuumFit(mapdir=self.mapdir, templates=self.templates)

        # expected attributes
        self.assertTrue(CFit.metallicity in ['Z0.0190'])
        self.assertTrue(CFit.library in ['CKC14z'])
        self.assertTrue(CFit.isochrone in ['Padova'])
        self.assertTrue(CFit.imf in ['Kroupa'])

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
