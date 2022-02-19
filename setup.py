#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard imports

import glob
import os
import sys

from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages

# DESI support code.

from desiutil.setup import DesiTest, DesiVersion, get_version

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
