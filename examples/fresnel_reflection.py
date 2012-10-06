#!/usr/local/bin/python2.5-mpi
# -*- coding: utf-8 -*-

""" Transmittance and reflectance through a thin gold layer.

This script is to obtain the transmittance and reflectance of 
TE polarized light through a thin gold layer.

"""

from __future__ import division

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
from gmes import *

x_size, y_size = 4, 0.5
SIZE = (x_size, y_size, 0)
angle = 0
wl= 10
k0 = 2 * pi / wl

a=20e-9
dp = DrudePole(omega=13.1839e15 * a / c0,
               gamma=0.109173e15 * a / c0)
cp1 = CriticalPoint(amp=0.273222,
                    phi=-1.18299,
                    omega=3.88123e15 * a / c0,
                    gamma=0.452006e15 * a / c0)
cp2 = CriticalPoint(amp=3.04155,
                    phi=-1.09115,
                    omega=4.20737 * a / c0,
                    gamma=2.35409 * a / c0)
gold = DcpPlrc(eps_inf=1.11683, mu_inf=1, dps=(dp,), cps=(cp1, cp2))

space = Cartesian(size=SIZE, resolution=100, parallel=True)
geom_list = [DefaultMedium(Dielectric()),
             Cylinder(center=(0, 0, 0),
                      axis=(1, 0, 0),
                      radius=1000,
                      height=1,
                      material=gold),
             Shell(material=Cpml(),
                   minus_y=False, plus_y=False)]
source_list = [GaussianBeam(
        src_time=Continuous(freq=1/wl, width=50),
        directivity=MinusX,
        center=(0.6, 0, 0),
        size=(0, y_size + 1, 1),
        direction=(-1 * cos(angle), sin(angle), 0),
        polarization=(0, 0, 1))]

my_fdtd = TMzFDTD(space,
                  geom_list,
                  source_list,
                  bloch=(0, k0 * sin(angle), 0))

my_fdtd.init()
my_fdtd.set_probe((0.7, 0, 0), 'r_wl=%f' % wl)
my_fdtd.set_probe((-0.7, 0, 0), 't_wl=%f' % wl)
my_fdtd.step_until_t(200)
