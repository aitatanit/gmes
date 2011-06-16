#!/usr/bin/env python

"""Models two-dimensional TMz cylindrical-wave propagation in air.

This script models two-dimensional TMz cylindrical-wave propagation in
air. A single Ez component located at the center of the space 
oscillates sinusoidally. A simple on-time visualization display will 
show the Ez, Hx, and Hy fields of the outgoing wave distributed within
the grid. You can compare the spatial-symmetry properties of these 
fields with respect to the center of the space where the excitation is
applied.

"""

from time import time

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

start = time()

SIZE = (10,10,0)

space = geometry.Cartesian(size=SIZE, resolution=20, parallel=True)
geom_list = (geometry.DefaultMedium(material=material.Dielectric()), 
             geometry.Boundary(material=material.UPML(), thickness=0.5, size=SIZE))
src_list = (source.Dipole(src_time=source.Continuous(freq=0.8), 
                          component=constants.Ez, pos=(0,0,0)),)

my_fdtd = fdtd.TMzFDTD(space, geom_list, src_list)

print time() - start

my_fdtd.show_ez(constants.Z, 0)
my_fdtd.show_hx(constants.Z, 0)
my_fdtd.show_hy(constants.Z, 0)

start = time()
my_fdtd.step_utill_n(100)
print time() - start
