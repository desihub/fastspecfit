import unittest

def test_suite():
    """Returns unittest.TestSuite of fastspecfit tests for use by setup.py"""

    from os.path import dirname
    fastspecfit_dir = dirname(dirname(__file__))
    print(fastspecfit_dir)
    return unittest.defaultTestLoader.discover(fastspecfit_dir,
        top_level_dir=dirname(fastspecfit_dir))
