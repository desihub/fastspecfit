#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file.

import sys
from setuptools import setup

# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.

API_HELP = """
Note: Generating api.rst files is no longer done using 'python setup.py api'. Instead
you will need to run:

    desi_api_file

which is part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

MODULE_HELP = """
Note: Generating Module files is no longer done using 'python setup.py api'. Instead
you will need to run:

    desiInstall

or

    desi_module_file

depending on your exact situation.  desiInstall is preferred.  Both commands are
part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

VERSION_HELP = """
Note: Generating version strings is no longer done using 'python setup.py version'. Instead
you will need to run:

    desi_update_version [-t TAG] desiutil

which is part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    pytest

If you don't already have pytest installed, you can install it with:

    pip install pytest
"""

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py {0}'. Instead you will need to run:

    sphinx-build -W --keep-going -b html doc doc/_build/html

If you don't already have Sphinx installed, you can install it with:

    pip install Sphinx
"""

message = {'api': API_HELP,
           'module_file': MODULE_HELP,
           'test': TEST_HELP,
           'version': VERSION_HELP,
           'build_docs': DOCS_HELP.format('build_docs'),
           'build_sphinx': DOCS_HELP.format('build_sphinx'), }

for m in message:
    if m in sys.argv:
        print(message[m])
        sys.exit(1)

setup()






#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard imports

import glob
import os, re
import sys
from os.path import abspath, basename, exists, isdir, isfile, join

from setuptools import setup, find_packages

# DESI support code.

#from desiutil.setup import DesiTest, DesiVersion, get_version

def find_version_directory(productname):
    """Return the name of a directory containing version information.

    Looks for files in the following places:

    * py/`productname`/_version.py
    * `productname`/_version.py

    Parameters
    ----------
    productname : :class:`str`
        The name of the package.

    Returns
    -------
    :class:`str`
        Name of a directory that can or does contain version information.

    Raises
    ------
    IOError
        If no valid directory can be found.
    """
    setup_dir = abspath('.')
    if isdir(join(setup_dir, 'py', productname)):
        version_dir = join(setup_dir, 'py', productname)
    elif isdir(join(setup_dir, productname)):
        version_dir = join(setup_dir, productname)
    else:
        raise IOError("Could not find a directory containing version information!")
    return version_dir

def get_version(productname):
    """Get the value of ``__version__`` without having to import the module.

    Parameters
    ----------
    productname : :class:`str`
        The name of the package.

    Returns
    -------
    :class:`str`
        The value of ``__version__``.
    """
    ver = 'unknown'
    try:
        version_dir = find_version_directory(productname)
    except IOError:
        return ver
    version_file = join(version_dir, '_version.py')
    if not isfile(version_file):
        update_version(productname)
    with open(version_file, "r") as f:
        for line in f.readlines():
            mo = re.match("__version__ = '(.*)'", line)
            if mo:
                ver = mo.group(1)
    return ver

# Begin setup

setup_keywords = dict()

# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.

setup_keywords['name'] = 'fastspecfit'
setup_keywords['description'] = 'Fast fitting of DESI spectroscopy and photometry.'
setup_keywords['author'] = 'John Moustakas & the DESI Collaboration'
setup_keywords['author_email'] = 'jmoustakas@siena.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/desihub/fastspecfit'

# END OF SETTINGS THAT NEED TO BE CHANGED.

setup_keywords['version'] = get_version(setup_keywords['name'])

# Use README.rst as long_description.

setup_keywords['long_description'] = ''
if os.path.exists('README.rst'):
    with open('README.rst') as readme:
        setup_keywords['long_description'] = readme.read()

# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.

# Treat everything in bin/ except *.rst as a script to be installed.

if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>3.0.0)']
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'':'py'}
setup_keywords['test_suite']='{name}.test.test_suite'.format(**setup_keywords)

# Add internal data directories.
setup_keywords['package_data'] = {'fastspecfit': ['data/*', 'test/data/*']}

# Run setup command.
setup(**setup_keywords)
