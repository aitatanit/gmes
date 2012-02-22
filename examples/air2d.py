#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Models two-dimensional TMz cylindrical-wave propagation in air.

This script models two-dimensional TMz cylindrical-wave propagation in air.
A single Ez component located at the center of the space oscillates 
sinusoidally. A simple on-time visualization display will show the Ez, Hx, 
and Hy fields of the outgoing wave distributed within the grid. You can 
compare the spatial-symmetry properties of these fields with respect to the
center of the space where the excitation is applied.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

SIZE = (10,10,0)

space = geometry.Cartesian(size=SIZE, resolution=20)
geom_list = (geometry.DefaultMedium(material=material.Dielectric()),
             geometry.Boundary(material=material.UPML(), thickness=0.5))
src_list = (source.PointSource(src_time=source.Continuous(freq=0.8), 
                               component=constant.Ez, pos=(0,0,0)),)

my_fdtd = fdtd.TMzFDTD(space, geom_list, src_list)

my_fdtd.init()

my_fdtd.show_ez(constant.Z, 0)
my_fdtd.show_hx(constant.Z, 0)
my_fdtd.show_hy(constant.Z, 0)

my_fdtd.step_until_t(100)
