#!/usr/bin/env python

# System imports
from distutils.core import setup, Extension
import os

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

PACKAGE = 'gmes'
VERSION = open('VERSION').read().strip()
base = os.getcwd() + '/'

# _pointwise_material module
_pointwise_material = Extension(name = 'gmes._pointwise_material',
                                sources = ['src/pointwise_material.i',
                                           'src/pointwise_material.cc',
                                           'src/pointwise_cpml.cc',
                                           'src/pointwise_dielectric.cc',
                                           'src/pointwise_dummy.cc',
                                           'src/pointwise_zero.cc',
                                           'src/pointwise_one.cc',
                                           'src/pointwise_upml.cc',
                                           'src/constants.cc'],
                                depends = ['src/pointwise_material.hh',
                                           'src/pointwise_cpml.hh',
                                           'src/pointwise_dielectric.hh',
                                           'src/pointwise_dummy.hh',
                                           'src/pointwise_zero.hh',
                                           'src/pointwise_one.hh',
                                           'src/pointwise_upml.hh',
                                           'src/constants.hh'],
                                include_dirs = [numpy_include],
                                swig_opts = ['-python', '-c++'],
                                language = 'c++',
                                extra_compile_args=[])

# _constants module
_constants = Extension(name = 'gmes._constants',
                       sources = ['src/constants.i',
                                  'src/constants.cc'],
                       depends = ['src/constants.hh'],
                       include_dirs = [numpy_include],
                       swig_opts = ['-python', '-c++'],
                       language = 'c++',
                       extra_compile_args=[])

setup(name = PACKAGE,
      version = VERSION,
      description = "GIST Maxwell's Equations Solver",
      long_description = """
      GMES is a Python package to solve the Maxwell's equations
      using the explicit Finite-Difference Time-Domain method.""",
      author = 'Kyungwon Chun',
      author_email = 'kwchun@gist.ac.kr',
      maintainer = 'Kyungwon Chun',
      maintainer_email = 'kwchun@gist.ac.kr',
      url = 'http://sourceforge.net/projects/gmes',
      license = 'http://www.gnu.org/licenses/gpl.html',
      packages = [PACKAGE],
#      py_modules = ['gmes.dielectric',
#                    'gmes.fdtd',
#                    'gmes.geometry',
#                    'gmes.material',
#                    'gmes.show',
#                    'gmes.sources',
#                    'gmes.pointwise_material'],
      ext_modules = [_pointwise_material,
                     _constants])
