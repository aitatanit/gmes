#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from gmes.constants import *


class DummyElectricCmplx(object):
    def __init__(self, idx, epsilon_r=1):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_r) * epsilon0
        
    def update(self, efield, hfield1, hfield2, d1, d2, dt, t):
        pass
    
    
class DummyExCmplx(DummyElectricCmplx): pass
    

class DummyEyCmplx(DummyElectricCmplx): pass
    

class DummyEzCmplx(DummyElectricCmplx): pass
    
    
class DummyMagneticCmplx(object):
    def __init__(self, idx, mu_r=1):
        self.i, self.j, self.k = array(idx, int)
        self.mu = float(mu_r) * mu0
        
    def update(self, hfield, efield1, efield2, d1, d2, dt, t):
        pass
    
    
class DummyHxCmplx(DummyMagneticCmplx): pass
    
    
class DummyHyCmplx(DummyMagneticCmplx): pass
    
    
class DummyHzCmplx(DummyMagneticCmplx): pass
    
    
if __name__ == '__main__':
    pass
 