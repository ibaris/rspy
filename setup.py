# -*- coding: utf-8 -*-
"""
Created on 31.03.2019 by Ismail Baris

For COPYING and LICENSE details, please refer to the LICENSE file
(c) 2017- Ismail Baris

"""
import os
from glob import glob

import numpy
from setuptools import find_packages

try:
    from Cython.Distutils import build_ext

except ImportError:
    USING_CYTHON = False

else:
    USING_CYTHON = True

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension


# ------------------------------------------------------------------------------------------------------------
# Environmental Functions
# ------------------------------------------------------------------------------------------------------------
def get_version():
    """
    Function to get the version.

    Returns
    -------
    out : str
    """
    version = dict()

    with open("rspy/__version__.py") as fp:
        exec (fp.read(), version)

    return version['__version__']


def get_packages():
    """
    Function get necessary modules.

    Returns
    -------
    out : setuptools.find_packages
    """
    return find_packages(exclude=['docs'])


# Open and Read Requirements --------------------------------------------------------------------------------
with open('requirements.txt') as f:
    required = f.read().splitlines()

# ------------------------------------------------------------------------------------------------------------
# Build Cython Extensions
# ------------------------------------------------------------------------------------------------------------
ext_modules = []
cmdclass = {}

ext = 'pyx' if USING_CYTHON else 'c'  # Define extension. If cython is present then use the pyx files.
sources = glob('rspy/**/*/*.{0}'.format(ext))  # Define source directories.

for item in sources:
    ext_modules.append(Extension(item.split('.')[0].replace(os.path.sep, '.'),
                                 [item],
                                 include_dirs=['rspy/bin', os.path.dirname(item), '.']))

if USING_CYTHON:
    cmdclass.update({'build_ext': build_ext})

setup(name='rspy',
      author="Ismail Baris",
      maintainer='Ismail Baris',
      version=get_version(),
      description='Fundamental Formulas for Radar and Angle Management',

      # Build Extensions -----------------------------------------------------------------------------------
      packages=get_packages(),
      cmdclass=cmdclass,
      include_dirs=[numpy.get_include()],
      ext_modules=ext_modules,

      # General Package Information ------------------------------------------------------------------------
      # ~ license='APACHE 2',
      url='https://github.com/ibaris/rspy',
      long_description='Fundamental Formulas for Radar and Angle Management',

      keywords=["radar", "remote-sensing", "optics", "integration",
                "microwave", "estimation", "physics", "radiative transfer"],

      # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Developers",
          "Natural Language :: English",
          "License :: Other/Proprietary License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Programming Language :: Python :: 2",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.3",
          "Programming Language :: Python :: 3.4"],

      # Include Package Data -------------------------------------------------------------------------------
      package_data={'rspy/bin/bin_unit': ['*.pxd', '*.c']},
      include_package_data=True,
      zip_safe=False,

      # Package Requirement --------------------------------------------------------------------------------
      install_requires=required,

      # Package Testing ------------------------------------------------------------------------------------
      setup_requires=['pytest-runner'],
      tests_require=['pytest'])
