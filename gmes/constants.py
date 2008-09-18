#!/usr/bin/env python

from math import sqrt

# physical constants

# Vaccum permittivity
# const double epsilon0 = 8.854187817e-12; # in SI unit (Ampere Second / Meter / Vol)
epsilon0 = 1.0 # normalized

# Vacuum permeability
# const double mu0 = 1.2566370614359173e-6 # in SI unit (Second Volt / Ampere / Meter)
mu0 = 1.0 #  Normalized

# speed of light in vacuum
c0 = 1 / sqrt(epsilon0 * mu0)

# vacuum impedance
Z0 = sqrt(mu0 / epsilon0)


# decimal factors
PETA = 1e15
TERA = 1e12
GIGA = 1e9
MEGA = 1e6
KILO = 1e3
MILLI = 1e-3
MICRO = 1e-6
NANO = 1e-9
PICO = 1e-12
FEMTO = 1e-15
ATTO = 1e-18

# The following classes should be used as singletons.

class Component:
    def __init__(self):
        raise SyntaxError, 'component classes should be used as a class'
    
    
class Electric(Component):
    pass


class Ex(Electric):
    tag = 0

    
class Ey(Electric):
    tag = 1
    
    
class Ez(Electric):
    tag = 2
    
class Magnetic(Component):
    pass


class Hx(Magnetic):
    tag = 3
    
    
class Hy(Magnetic):
    tag = 4
    
    
class Hz(Magnetic):
    tag = 5


class Directional:
    def __init__(self):
        raise SyntaxError, 'directional classes should be used as a class'
    
    
class X(Directional): pass
class Y(Directional): pass
class Z(Directional): pass
class PlusX(X): pass
class MinusX(X): pass
class PlusY(Y): pass
class MinusY(Y): pass
class PlusZ(Z): pass
class MinusZ(Z): pass
