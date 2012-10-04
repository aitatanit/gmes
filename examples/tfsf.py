#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A simple example which showing how to use launch a planewave.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

SIZE = (5,5,0)

space = Cartesian(size=SIZE, resolution=20)
geom_list = [DefaultMedium(material=Dielectric()),
             Shell(material=Cpml())]
src_list = [TotalFieldScatteredField(src_time=Continuous(freq=0.8),
                                     center=(0,0,0),
                                     size=(1,1,1),
                                     direction=(0,-1,0),
                                     polarization=(0,0,1))]

my_fdtd = TMzFDTD(space, geom_list, src_list)

my_fdtd.init()

# my_fdtd.show_field(Ez, Z, 0)
# my_fdtd.show_field(Hx, Z, 0)
# my_fdtd.show_field(Hy, Z, 0)
my_fdtd.step_until_t(10)
# my_fdtd.write_field(Ez, (-5,-5,0), (5,5,0))
