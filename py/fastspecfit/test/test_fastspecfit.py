"""
fastspecfit.test.test_continuum
===============================

Test fastspecfit.fastspecfit.fastspec

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
        os.environ['FASTSPECFIT_TEMPLATES'] = resource_filename('fastspecfit.test', 'data')
        cls.mapdir = resource_filename('fastspecfit.test', 'data')
        cls.redrockfile = resource_filename('fastspecfit.test', 'data/redrock-4-80613-thru20210324.fits')
        cls.cwd = os.getcwd()
        cls.outdir = tempfile.mkdtemp()
        cls.fastspec_outfile = os.path.join(cls.outdir, 'fastspec.fits')
        cls.fastphot_outfile = os.path.join(cls.outdir, 'fastphot.fits')

        #self.ra = np.array([84.56347552,  88.25858593,
        #                    85.18114653,  84.04246538, 83.22215524])

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_ContinuumFit(self):
        """Test the ContinuumTools class."""
        from fastspecfit.continuum import ContinuumFit
        CFit = ContinuumFit(mapdir=self.mapdir)

        # expected attributes
        self.assertTrue(CFit.metallicity in ['Z0.0190'])
        self.assertTrue(CFit.library in ['CKC14z'])
        self.assertTrue(CFit.isochrone in ['Padova'])
        self.assertTrue(CFit.imf in ['Kroupa'])

    def test_fastphot(self):
        """Test fastphot."""
        import fitsio
        from fastspecfit.fastspecfit import fastphot, parse

        cmd = 'fastphot {} -o {} --mapdir {}'.format(self.redrockfile, self.fastphot_outfile, self.mapdir)
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
    
        cmd = 'fastspec {} -o {} --mapdir {}'.format(self.redrockfile, self.fastspec_outfile, self.mapdir)
        args = parse(options=cmd.split()[1:])
        fastspec(args=args)
    
        self.assertTrue(os.path.exists(self.fastspec_outfile))
    
        fits = fitsio.FITS(self.fastspec_outfile)
        for hdu in fits:
            if hdu.has_data(): # skip zeroth extension
                self.assertTrue(hdu.get_extname() in ['METADATA', 'FASTSPEC'])

if __name__ == '__main__':
    unittest.main()
