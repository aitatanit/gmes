#!/usr/bin/env python

"""Launch a dipole source in two dimensional free space."""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *


SIZE = (5,0,5)

space = geometric.Cartesian(size=SIZE, resolution=20, parallel=True)
geom_list = (geometric.DefaultMaterial(material=material.Dielectric()), 
             geometric.Boundary(material=material.CPML(), thickness=0.5, size=SIZE))
src_list = (source.Dipole(src_time=source.Continuous(freq=0.8), 
                          component=constants.Ey, pos=(0,0,2)),)

my_fdtd = fdtd.TMyFDTD(space, geom_list, src_list)

try:
    import psyco
    psyco.full()
    psyco.log()
    psyco.profile()
except ImportError:
    pass

my_fdtd.show_ey(constants.Y, 0)

while True:
    my_fdtd.step()
    if space.my_id == 0:
        print int(my_fdtd.time_step.n)
