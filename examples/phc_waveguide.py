#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Shows a Ez field in a photonic crystal waveguide.

A simple example showing the Ez field in a two-dimensional 
photonic crystal waveguide.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

space = Cartesian(size=(16,8,0), resolution=20)
geom_list = [DefaultMedium(material=Dielectric())]
geom_list.extend([Cylinder(material=Dielectric(8.9),
                           radius=0.38,
                           center=(x,y,0))
                  for x in xrange(-8, 9) 
                  for y in xrange(-4, 5)
                  if y != 0])
geom_list.append(Shell(material=Cpml()))
src_list = [PointSource(src_time=Continuous(freq=0.43),
                        component=Ez,
                        center=(-7,0,0))]
my_fdtd = TMzFDTD(space, geom_list, src_list)
my_fdtd.init()
my_fdtd.show_permittivity(Ez, Z, 0)
my_fdtd.show_field(Ez, Z, 0)
my_fdtd.step_until_t(200)
