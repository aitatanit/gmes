#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from gmes.constants import *


class ZeroElectricCmplx(object):
    def __init__(self, idx, epsilon_r):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_r) * epsilon0
        
        
class ZeroExCmplx(ZeroElectricCmplx):
    def __init__(self, idx, epsilon_r=1):
        ZeroElectricCmplx.__init__(self, idx, epsilon_r)

    def update(self, ex, hz, hy, dt, dy, dz):
        ex[self.i, self.j, self.k] = 0
        
        
class ZeroEyCmplx(ZeroElectricCmplx):
    def __init__(self, idx, epsilon_r=1):
        ZeroElectricCmplx.__init__(self, idx, epsilon_r)

    def update(self, ey, hx, hz, dt, dz, dx):
        ey[self.i, self.j, self.k] = 0
        
        
class ZeroEzCmplx(ZeroElectricCmplx):
    def __init__(self, idx, epsilon_r=1):
        ZeroElectricCmplx.__init__(self, idx, epsilon_r)

    def update(self, ez, hy, hx, dt, dx, dy):
        ez[self.i, self.j, self.k] = 0
        
        
class ZeroMagneticCmplx(object):
    def __init__(self, idx, mu_r):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(mu_r) * mu0
        
        
class ZeroHxCmplx(ZeroMagneticCmplx):
    def __init__(self, idx, mu_r=1):
        ZeroMagneticCmplx.__init__(self, idx, mu_r)

    def update(self, hx, ez, ey, dt, dy, dz):
        hx[self.i, self.j, self.k] = 0
        
        
class ZeroHyCmplx(ZeroMagneticCmplx):
    def __init__(self, idx, mu_r=1):
        ZeroMagneticCmplx.__init__(self, idx, mu_r)

    def update(self, hy, ex, ez, dt, dz, dx):
        hy[self.i, self.j, self.k] = 0
        

class ZeroHzCmplx(ZeroMagneticCmplx):
    def __init__(self, idx, mu_r=1):
        ZeroMagneticCmplx.__init__(self, idx, mu_r)

    def update(self, hz, ey, ex, dt, dx, dy):
        hz[self.i, self.j, self.k] = 0
            
        
if __name__ == '__main__':
    pass
 