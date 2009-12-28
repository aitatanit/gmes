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

# _pw_material module
_pw_material = Extension(name = 'gmes._pw_material',
                              sources = ['src/pw_material.i',
                                         'src/pw_material.cc',
                                         'src/pw_cpml.cc',
                                         'src/pw_dielectric.cc',
                                         'src/pw_dummy.cc',
                                         'src/pw_const.cc',
                                         'src/pw_upml.cc',
                                         'src/pw_drude.cc'],
                              depends = ['src/pw_material.hh',
                                         'src/pw_cpml.hh',
                                         'src/pw_dielectric.hh',
                                         'src/pw_dummy.hh',
                                         'src/pw_const.hh',
                                         'src/pw_upml.hh',
                                         'src/pw_drude.hh',
                                         'src/constants.hh'],
                              include_dirs = [numpy_include],
                              swig_opts = ['-c++', '-outdir', 'gmes'],
                              #language = 'c++',
                              extra_compile_args=[])

# _constants module
_constants = Extension(name = 'gmes._constants',
                       sources = ['src/constants.i',
                                  'src/constants.cc'],
                       depends = ['src/constants.hh'],
                       include_dirs = [numpy_include],
                       swig_opts = ['-c++', '-outdir', 'gmes'],
                       #language = 'c++',
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
      ext_modules = [_pw_material, _constants])
