#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Shows an Ez field in a dielectric slab waveguide.

A simple example showing the Ez field in a dielectric slab 
waveguide. This is a GMES version of the script in Fig.12 of

A. F. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J. D. 
Joannopoulos, and S. G. Johnson, "Meep: A flexible free-software 
package for electromagnetic simulations by the FDTD method," 
Comput. Phys. Commun. 181, 687-702 (2010).

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes import *

space = Cartesian(size=(16,8,0), resolution=10)
geom_list = [DefaultMedium(material=Dielectric()),
             Block(material=Dielectric(12),
                   size=(inf, 1, inf)),
             Shell(material=Cpml())]
src_list = [PointSource(src_time=Continuous(freq=0.15),
                        component=Ez,
                        center=(-7,0,0))]
my_fdtd = TMzFDTD(space, geom_list, src_list)
my_fdtd.init()
my_fdtd.show_field(Ez, Z, 0)
my_fdtd.step_until_t(200)
