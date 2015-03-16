#!/usr/bin/env pythonssssssssssssssssssssssssssssssssssssssssssssssssssss
# -*- coding: utf-8 -*-

##    gmes - GIST Maxwell's Equations Solver
##    Copyright (C) 2007-2012  Kyungwon Chun
##
##    This library is free software; you can redistribute it and/or
##    modify it under the terms of the GNU Library General Public
##    License as published by the Free Software Foundation; either
##    version 3 of the License, or (at your option) any later version.
##
##    This library is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##    Library General Public License for more details.
##
##    You should have received a copy of the GNU Library General Public
##    License along with this library; if not, write to the Free
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##
##    Kyungwon Chun
##    kwchun@gist.ac.kr

"""A Python implementation of the explicit FDTD method

GMES is a Python implementation of the explicit finite-difference
time-domain (FDTD) method. It is designed to simulate the photonic
device in 1, 2, and 3-d Cartesian coordinates.

Modules:
    fdtd --- Provide various simulation classes suitable for 1, 2, and 3-d.
    geometry --- Provide coordinate and geometric primitives.
    show --- Real-time display classes
    constant --- Physical and simulation constants
    source --- Define the input sources
    pw_source --- Source update mechanism
    material --- Define the propagating medium
    pw_material --- Provide the update mechanism 

"""

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')

from fdtd import *
from geometry import *
from constant import *
from source import *
from material import *

import fdtd, geometry, show, constant, source, material
import pw_material, pw_source

# List here only the objects we want to be publicly available
_module = ['fdtd', 'geometry', 'show', 'constant', 'source', 'pw_source', 'material', 'pw_material']
_class = ['TimeStep', 'FDTD', 'TExFDTD', 'TEyFDTD', 'TEzFDTD', 'TMxFDTD', 'TMyFDTD', 'TMzFDTD', 'TEMxFDTD', 'TEMyFDTD', 'TEMzFDTD', 
          'Cartesian', 'DefaultMedium', 'Cone', 'Cylinder', 'Block', 'Ellipsoid', 'Sphere', 'Shell', 
          'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', 'Jx', 'Jy', 'Jz', 'Mx', 'My', 'Mz', 'X', 'Y', 'Z', 'PlusX', 'MinusX', 'PlusY', 'MinusY', 'PlusZ', 'MinusZ', 
          'Continuous', 'Bandpass', 'DifferentiatedGaussian', 'PointSource', 'TotalFieldScatteredField', 'GaussianBeam', 
          'Dummy', 'Const', 'Dielectric', 'Upml', 'Cpml', 'DrudePole', 'LorentzPole', 'CriticalPoint', 'DcpAde', 'DcpPlrc', 'DcpRc', 'Drude', 'Lorentz', 'Dm2']
_constant = ['pi', 'c0', 'mu0', 'eps0', 'Z0', 'PETA', 'TERA', 'GIGA', 'MEGA', 'KILO', 'MILLI', 'MICRO', 'NANO', 'PICO', 'FEMTO', 'ATTO',
             'inf']
__all__ = []
__all__.extend(_module)
__all__.extend(_class)
__all__.extend(_constant)
