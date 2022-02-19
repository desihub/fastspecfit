# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
test top-level fastspecfit functions

"""
import unittest
import re
from fastspecfit import __version__ as theVersion

class TestTopLevel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.versionre = re.compile(r'([0-9]+!)?([0-9]+)(\.[0-9]+)*((a|b|rc|\.post|\.dev)[0-9]+)?')

    @classmethod
    def tearDownClass(cls):
        pass

    def test_version(self):
        """Ensure the version conforms to PEP386/PEP440.

        """
        self.assertRegex(theVersion,self.versionre)

if __name__ == '__main__':
    unittest.main()
