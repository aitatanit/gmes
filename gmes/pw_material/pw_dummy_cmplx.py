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
    def __init__(self, idx, epsilon_r):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_r) * epsilon0
        
        
class DummyExCmplx(DummyElectricCmplx):
    def __init__(self, idx, epsilon_r=1):
        DummyElectricCmplx.__init__(self, idx, epsilon_r)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        pass
    

class DummyEyCmplx(DummyElectricCmplx):
    def __init__(self, idx, epsilon_r=1):
        DummyElectricCmplx.__init__(self, idx, epsilon_r)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        pass
    

class DummyEzCmplx(DummyElectricCmplx):
    def __init__(self, idx, epsilon_r=1):
        DummyElectricCmplx.__init__(self, idx, epsilon_r)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        pass
    
    
class DummyMagneticCmplx(object):
    def __init__(self, idx, mu_r):
        self.i, self.j, self.k = array(idx, int)
        self.mu = float(mu_r) * mu0
        
        
class DummyHxCmplx(DummyMagneticCmplx):
    def __init__(self, idx, mu_r=1):
        DummyMagneticCmplx.__init__(self, idx, mu_r)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        pass
    
    
class DummyHyCmplx(DummyMagneticCmplx):
    def __init__(self, idx, mu_r=1):
        DummyMagneticCmplx.__init__(self, idx, mu_r)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        pass
    
    
class DummyHzCmplx(DummyMagneticCmplx):
    def __init__(self, idx, mu_r=1):
        DummyMagneticCmplx.__init__(self, idx, mu_r)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        pass
    
    
if __name__ == '__main__':
    pass
 