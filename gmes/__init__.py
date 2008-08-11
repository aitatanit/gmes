#!/usr/bin/env python

##    gmes - GIST Maxwell's Equations Solver
##    Copyright (C) 2007  Kyungwon Chun
##
##    This library is free software; you can redistribute it and/or
##    modify it under the terms of the GNU Library General Public
##    License as published by the Free Software Foundation; either
##    version 2 of the License, or (at your option) any later version.
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
device in 1-d, 2-d, and 3-d Cartesian coordinates.

Available modules:
fdtd --- Provide various simulation classes suitable for 1, 2, and 3-d cases.
geometric --- Provide coordinate and geometric primitives
show --- Real-time display classes
constants --- Physical and simulation constants
source --- Define the input sources
pointwise_source --- Source update mechanism
material --- Define the propagating medium
pointwise_material --- Provide the update mechanism 

"""

# List here only the objects we want to be publicly available
__all__ = ['fdtd', 'geometric', 'show', 'constants', 
           'source', 'pointwise_source', 
           'material', 'pointwise_material'
           ]
