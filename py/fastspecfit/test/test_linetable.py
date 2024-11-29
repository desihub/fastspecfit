"""
fastspecfit.test.test_linetable
===============================

"""
import os
import numpy as np
import unittest

class TestUtil(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass


    @classmethod
    def tearDownClass(cls):
        pass


    def test_linetable(self):
        import astropy
        from fastspecfit.linetable import LineTable
        emline_table = LineTable(emlines_file=None)
        self.assertTrue(os.path.isfile(emline_table.file))
        self.assertTrue(type(emline_table.table) == astropy.table.Table)

        with self.assertRaises(IOError):
            LineTable(emlines_file='doesnotexist.ecsv')
