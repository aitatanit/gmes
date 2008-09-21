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

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *


SIZE = (5,5,0)

space = geometry.Cartesian(size=SIZE, resolution=20, parallel=True)
geom_list = (geometry.DefaultMaterial(material=material.Dielectric()), 
             geometry.Boundary(material=material.UPML(), thickness=0.5, size=SIZE))
src_list = (source.Dipole(src_time=source.Continuous(freq=0.8), 
                          component=constants.Ez, pos=(0,0,0)),)

my_fdtd = fdtd.TMzFDTD(space, geom_list, src_list)

try:
    import psyco
    psyco.full()
except ImportError:
    pass

my_fdtd.show_ez(constants.Z, 0)
my_fdtd.show_hx(constants.Z, 0)
my_fdtd.show_hy(constants.Z, 0)

while True:
    my_fdtd.step()
    if space.my_id == 0:
        print int(my_fdtd.time_step.n)
