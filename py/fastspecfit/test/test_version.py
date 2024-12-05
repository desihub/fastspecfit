"""
fastspecfit.test.test_version
=============================

"""
import re
from fastspecfit import __version__ as theVersion


def test_version():
    """Ensure the version conforms to PEP386/PEP440.

    """
    assert(re.search(r'([0-9]+!)?([0-9]+)(\.[0-9]+)*((a|b|rc|\.post|\.dev)[0-9]+)?', theVersion))
