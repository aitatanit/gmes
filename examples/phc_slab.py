#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from numpy import *
from time import *

from gmes import *

FREQ = 0.3
#FREQ = 0.42
RADIUS = 0.35
#RADIUS = 0.48
SLAB_CORE = .5 / sqrt(11.8336)
#SLAB_CORE = .5 / sqrt(13)
SLAB_THICK = 3 * SLAB_CORE
PML_THICK = .5
ZLENGTH = 3 * SLAB_THICK + 2 * PML_THICK
SIZE = (15, 15, 0)

AIR = material.Dielectric(1)
SiO2 = material.Dielectric(2.1316)
Si = material.Dielectric(11.8336)
#Si = material.Dielectric(13)

def make_hole(center):
    return geometric.Cylinder(material=AIR, axis=(0,0,1), radius=RADIUS,
                              height=SLAB_THICK, center=center)

def make_crystals(x_size, y_size):
    a1 = array((cos(pi/3), sin(pi/3)))
    a2 = array((cos(pi/3), -sin(pi/3)))
    crystals = []
    
    for i in arange(-ceil(x_size), ceil(x_size)):
        for j in arange(-ceil(y_size), ceil(y_size)):
            center = tuple(i * a1 + j * a2) + (0,)
            if fabs(center[0]) <= .5 * x_size + RADIUS and fabs(center[1]) <= .5 * y_size + RADIUS:
                crystals.append(make_hole(center))
    return tuple(crystals)

def fill_hole(center):
    filler = (geometric.Cylinder(material=SiO2, axis=(0,0,1), radius=RADIUS,
                                 height=SLAB_THICK, center=center),
              geometric.Cylinder(material=Si, axis=(0,0,1), radius=RADIUS,
                                 height=SLAB_CORE, center=center))

    return filler

def make_line_defect(length):
    line_defect = ()
    for i in arange(-.5 * length - 1, .5 * length + 2):
        line_defect += fill_hole((i, 0, 0))

    return tuple(line_defect)
    
geom_list = (geometric.DefaultMaterial(material=AIR),
             geometric.Block(material=SiO2, size=(SIZE[0], SIZE[1], SLAB_THICK)),
             geometric.Block(material=Si, size=(SIZE[0], SIZE[1], SLAB_CORE))) + \
             make_crystals(*SIZE[:2]) + make_line_defect(SIZE[0]) + \
             (geometric.Boundary(material=material.UPML(), thickness=PML_THICK, size=SIZE),)

space = geometric.Cartesian(size=SIZE, resolution=(25, 25, 10), parallel=True)

# src_list = (source.Dipole(src_time=source.Continuous(freq=FREQ),
#                          component=constants.Hz, pos=(-.5 * SIZE[0] + 1, 0, 0)),)

src_list = (source.Transparent(direction=constants.PlusX,
                               center=(-.5 * SIZE[0] + 1, 0, 0),
                               size=(0, RADIUS, 0),
                               freq=FREQ,
                               polarization=constants.Y),)

print 'communicator shape:', space.cart_comm.dims

start_time = time()
my_fdtd = fdtd.FDTD(space, geom_list, src_list)
if space.my_id == 0:
    print 'initialization:', time() - start_time

try:
    import psyco
    psyco.full()
except ImportError:
    pass

# shape = my_fdtd.material_ex.shape[:2]
# eps = empty(shape, float)
# for idx in ndindex(*shape):
#     eps[idx] = my_fdtd.material_ex[idx[0], idx[1], 0].epsilon
    
# import pylab as p

# p.imshow(eps)
# p.colorbar()
# p.xlabel('y')
# p.ylabel('x')
# p.show()

my_fdtd.show_hz(constants.Z, 0)

start_time = time()
while True:
    my_fdtd.step()
    if space.my_id == 0:
        print int(my_fdtd.time_step.n)
if space.my_id == 0:
    print '100 loop:', time() - start_time
