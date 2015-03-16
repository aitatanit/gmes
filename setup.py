#!/usr/bin/env python
# -*- coding: utf-8 -*-

# System imports
from os import getcwd
from glob import glob
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory. This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

PACKAGE = 'gmes'
VERSION = open('VERSION').read().strip()
base = getcwd() + '/'

pw_src_lst = glob('src/pw_*.cc')
pw_src_lst.extend(glob('src/pw_*.i'))
pw_dep_lst = glob('src/pw_*.hh')

# pw_material module
pw_material = Extension(name = 'gmes._pw_material',
                        sources = pw_src_lst,
                        depends = pw_dep_lst,
                        include_dirs = [numpy_include],
                        swig_opts = ['-c++', '-outdir', 'gmes'],
                        language = 'c++',
                        extra_compile_args=['-std=c++0x'])

# constant module
constant = Extension(name = 'gmes._constant',
                     sources = ['src/constant.i', 'src/constant.cc'],
                     depends = ['src/constant.hh'],
                     include_dirs = [numpy_include],
                     swig_opts = ['-c++', '-outdir', 'gmes'],
                     language = 'c++',
                     extra_compile_args=['-std=c++0x'])

# pygeom module
pygeom = Extension(name = 'gmes.pygeom',
                   sources = ['src/pygeom.pyx'],
                   include_dirs = [numpy_include])

# material module
material = Extension(name = 'gmes.material',
                     sources = ['src/material.pyx'])

setup(name = PACKAGE,
      version = VERSION,
      description = "GIST Maxwell's Equations Solver",
      long_description = """
      GMES is a free Python package to solve the Maxwell's 
      equations using the explicit Finite-Difference Time-Domain 
      method.""",
      author = 'Kyungwon Chun',
      author_email = 'kwchun@gist.ac.kr',
      maintainer = 'Kyungwon Chun',
      maintainer_email = 'kwchun@gist.ac.kr',
      url = 'http://sourceforge.net/projects/gmes',
      classifiers = [
        'Development Status :: 3 - Alpha',
        'User Interface :: Console/Terminal',
        'Intended Audience :: Science/Research',
        'License :: GNU General Public License version 3.0 (GPLv3)',
        'Operating System :: OS Portable (Source code to work with many OS platforms)',
        'Programming Language :: C++, Python',
        'Topic :: Physics, Simulations',
        'Translations: English, Korean'],
      license = 'http://www.gnu.org/licenses/gpl.html',
      packages = [PACKAGE],
      ext_modules = [pw_material, constant, pygeom, material],
      cmdclass = {'build_ext': build_ext})
