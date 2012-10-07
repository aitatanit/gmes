#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Shows a man-shaped structure.

This script shows how to set up a geometrical structure in the 
calculation domain. You can see a man-shaped structure in three-
dimensional space constructed using the script.

"""

import os, sys, datetime
new_path = os.path.abspath('../')
sys.path.append(new_path)

from datetime import datetime
print os.uname()
print 'python version:', sys.version
start_time = datetime.now()

from numpy import cross
from gmes import *

space = Cartesian(size=(6,6,6), resolution=20)
body = Block(Dielectric(1), size=(1,1,2))
head = Sphere(Dielectric(2), center=(0,0,-1.5), radius=0.5)
hat = Cone(Dielectric(3), center=(0,0.2,-2.15), axis=(0,0.2,-1), 
           radius=0.7, height=0.5)
leg1 = Cylinder(Dielectric(4), center=(0,0.5,2), axis=(0,-0.2,-1),
                radius=0.2, height=2)
leg2 = Cylinder(Dielectric(5), center=(0,-0.5,2), axis=(0,0.2,-1),
                radius=0.2, height=2)
arm1 = Ellipsoid(Dielectric(6), center=(0, 1.3, -0.2),
                 e1=(1, 0, 0), e2=(0, 1, -1), 
                 e3=cross((1, 0, 0), (0, 1, -1)),
                 size=(0.5, 0.5, 1.5))
arm2 = Ellipsoid(Dielectric(7), center=(0, -1.3, -0.2),
                 e1=(1, 0, 0), e2=(0, 1, 1), 
                 e3=cross((1, 0, 0), (0, 1, 1)),
                 size=(0.5, 0.5, 1.5))
geom_list = (DefaultMedium(Dielectric(10)),
             body, head, hat, leg1, leg2, arm1, arm2)

my_fdtd = fdtd.FDTD(space, geom_list, ())
my_fdtd.init()

my_fdtd.show_permittivity(Ex, X, 0)
my_fdtd.show_permittivity(Ex, Y, 0)
my_fdtd.show_permittivity(Ex, Z, 0)
