#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Launch a planwave in vacuum with a scatterer.

A simple example which showing how to launch a planewave in vacuum 
with a cylindrical dielectric scatterer.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

SIZE = (5,5,0)

space = Cartesian(size=SIZE, resolution=20)
geom_list = [DefaultMedium(Dielectric()),
             Cylinder(Dielectric(3),
                      center=(0,0,0),
                      radius=1,
                      axis=(0,0,1)),
             Shell(material=Cpml(), thickness=0.5)]
src_list = [TotalFieldScatteredField(src_time=Continuous(freq=0.8),
                                     center=(0,0,0),
                                     size=(3,3,1),
                                     direction=(1,-1,0),
                                     polarization=(0,0,1))]

my_fdtd = TMzFDTD(space, geom_list, src_list)

my_fdtd.init()

my_fdtd.show_field(Ez, Z, 0)
# my_fdtd.write_field(Ez, (-5,-5,0), (5,5,0))
my_fdtd.step_until_t(200)
