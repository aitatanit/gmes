#!/usr/local/bin/python2.5-mpi
# -*- coding: utf-8 -*-

import os, sys, datetime
new_path = os.path.abspath('../')
sys.path.append(new_path)

from datetime import datetime
print os.uname()
print 'python version:', sys.version
start_time = datetime.now()
print 'starting initialization:', start_time

from sys import argv
from math import pi, sin, cos
from numpy import inf
from gmes import material, geometry, fdtd, source, constant

ORDINAL = float(argv[1])
NJOBS = float(argv[2])
x_size, y_size = 4, 0.5
SIZE = (x_size, y_size, 0)
angle = 0
wl= 10 + (50 - 10) * ORDINAL / NJOBS
k0 = 2 * pi / wl
air = material.Dielectric()
gold = material.GoldRc(a=20e-9)
host = geometry.DefaultMedium(material=air)
cylinder = geometry.Cylinder(center=(0, 0, 0),
                             axis=(1, 0, 0),
                             radius=1000,
                             height=1,
                             material=gold)
boundary = geometry.Boundary(material=material.CPML(),
                             thickness=1,
                             minus_y=False, plus_y=False)

space = geometry.Cartesian(size=SIZE, resolution=100, parallel=True)
geom_list = (host, cylinder, boundary)
source_list = (source.GaussianBeam(
    src_time=source.Continuous(freq=1/wl,
                               width=50),
    directivity=constant.MinusX,
    center=(0.6, 0, 0),
    size=(0, y_size + 1, 1),
    direction=(-1 * cos(angle), sin(angle), 0),
    polarization=(0, 0, 1)),)

my_fdtd = fdtd.TMzFDTD(space,
                       geom_list,
                       source_list,
                       bloch=(0, k0 * sin(angle), 0))

# directory = os.path.dirname(__file__) + '/../data'
# my_fdtd.set_probe((0.7, 0, 0), directory + '/r_wl=%f' % wl)
# my_fdtd.set_probe((-0.7, 0, 0), directory + '/t_wl=%f' % wl)

# if os.uname()[1] == 'magi':
#     my_fdtd.show_ez(constant.Z, 0)

end_time = datetime.now()
print 'ending initialization:', end_time
print 'elasped time:', end_time - start_time

start_time = datetime.now()
print 'starting update:', start_time
my_fdtd.step_until_t(200)

print 'n:', my_fdtd.time_step.n
print 't:', my_fdtd.time_step.t
end_time = datetime.now()
print 'ending update:', end_time
print 'elasped time:', end_time - start_time
print my_fdtd.time_step.n
