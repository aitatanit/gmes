#!/usr/bin/env python

"""Simulate a photonic crystal slab waveguide.

The photonic crystal slab consists of a silicon-on-insulator substrate
with a triangular array of holes. The whole waveguide is in the air. 
This structure is presented at 

'N. Moll and G.-L. Bona, "Comparison of three-dimensional photonic 
crystal slab waveguides with two-dimensional photonic crystal 
waveguides: Efficient butt coupling into these photonic crystal 
waveguides," J. Appl. Phys., vol. 93, no. 9, pp. 4986-4991, 2003.'

A transparent source is located at the one end of the waveguide. You 
should be careful to execute this script, because it consumes lots of 
memory.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from numpy import *

from gmes import *


# Define simulation parameters.

FREQ = 0.3
RADIUS = 0.35
SLAB_CORE = .5 / sqrt(11.8336)
SLAB_THICK = 3 * SLAB_CORE
PML_THICK = .5
ZLENGTH = 3 * SLAB_THICK + 2 * PML_THICK
SIZE = (15, 15, 0)

AIR = material.Dielectric(1)
SiO2 = material.Dielectric(2.1316)
Si = material.Dielectric(11.8336)


def make_hole(center):
    """Punch a air hole on the slab."""
    return geometry.Cylinder(material=AIR, axis=(0,0,1), radius=RADIUS,
                              height=SLAB_THICK, center=center)


def make_crystals(x_size, y_size):
    """Punch holes on the lattice positions."""
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
    """Remove a air hole."""
    filler = (geometry.Cylinder(material=SiO2, axis=(0,0,1), radius=RADIUS,
                                 height=SLAB_THICK, center=center),
              geometry.Cylinder(material=Si, axis=(0,0,1), radius=RADIUS,
                                 height=SLAB_CORE, center=center))

    return filler


def make_line_defect(length):
    """Remove air holes to form a line defect."""
    line_defect = ()
    for i in arange(-.5 * length - 1, .5 * length + 2):
        line_defect += fill_hole((i, 0, 0))

    return tuple(line_defect)
    
    
geom_list = (geometry.DefaultMaterial(material=AIR),
             geometry.Block(material=SiO2, size=(SIZE[0], SIZE[1], SLAB_THICK)),
             geometry.Block(material=Si, size=(SIZE[0], SIZE[1], SLAB_CORE))) + \
             make_crystals(*SIZE[:2]) + make_line_defect(SIZE[0]) + \
             (geometry.Boundary(material=material.UPML(), thickness=PML_THICK, size=SIZE),)

space = geometry.Cartesian(size=SIZE, resolution=(25, 25, 10), parallel=True)

# src_list = (source.Dipole(src_time=source.Continuous(freq=FREQ),
#                          component=constants.Hz, pos=(-.5 * SIZE[0] + 1, 0, 0)),)

src_list = (source.Transparent(direction=constants.PlusX,
                               center=(-.5 * SIZE[0] + 1, 0, 0),
                               size=(0, RADIUS, 0),
                               freq=FREQ,
                               polarization=constants.Y),)

my_fdtd = fdtd.FDTD(space, geom_list, src_list)

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

while True:
    my_fdtd.step()
    if space.my_id == 0:
        print int(my_fdtd.time_step.n)
