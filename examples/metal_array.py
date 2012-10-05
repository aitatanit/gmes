#!/usr/bin/env python-mpi
# -*- coding: utf-8 -*-

"""A plasmon waveguide consisting of six silver nanospheres in the
air.

This script models a plasmon waveguide consisting of six silver 
nanospheres in the air. A dipole source oscilating along the array
is used for field excitation. A simple on-time visualization 
display will show the Ey fields of the propagating longitudinal
mode. This script was set to use the parallel environment when it 
is executed using a MPI-enabled Python interpreter, like 
'python-mpi' in the mpi4py package. This script requires about 1.1
GB of memory.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

class Silver(DcpPlrc):
    """This silver permittivity value in the range of 200-1,000 nm.
    
    This parameters has a fitness value of 0.0266134 to the real 
    permittivity data in the range of 200-1,000 nm of 
    
    P. B. Johnson and R. W. Christy, "Optical constants of the 
    noble metals,"  Phys. Rev. B 6, 4370 (1972).
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
                
        """
        dp1 = DrudePole(omega=1.38737e16 * a / c0, 
                        gamma=2.07331e13 * a / c0)
        cp1 = CriticalPoint(amp=1.3735,
                            phi=-0.504658,
                            omega=7.59914e15 * a / c0, 
                            gamma=4.28431e15 * a / c0)
        cp2 = CriticalPoint(amp=0.304478,
                            phi=-1.48944,
                            omega=6.15009e15 * a / c0, 
                            gamma=6.59262e14 * a / c0)
        DcpPlrc.__init__(self, eps_inf=0.89583, mu_inf=1, sigma=0,
                        dps=(dp1,), cps=(cp1,cp2))

space = Cartesian(size=(2, 8, 2), resolution=40, parallel=True)
geom_list = [DefaultMedium(Dielectric())]
for y in range(-2, 4):
    geom_list.append(Sphere(Silver(75 * NANO),
                                   radius=1.0 / 3,
                                   center=(0, y, 0)))
geom_list.append(Shell(Cpml(), thickness=0.5))
src_list = [PointSource(Continuous(freq=0.207),
                        center=(0, -3, 0),
                        component=Jy)]
my_fdtd = FDTD(space, geom_list, src_list, courant_ratio=0.5)
my_fdtd.init()
my_fdtd.show_field(Ey, Z, 0, (-1e-5, 1e-5))
my_fdtd.step_until_t(50)
