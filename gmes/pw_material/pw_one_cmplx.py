#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from gmes.constants import *


class OneElectricCmplx(object):
    def __init__(self, idx, epsilon_r=1):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_r) * epsilon0
        
    def update(self, efield, hfield1, hfield2, d1, d2, dt, t):
        efield[self.i, self.j, self.k] = 1
        
            
class OneExCmplx(OneElectricCmplx): pass
        
        
class OneEyCmplx(OneElectricCmplx): pass
        
        
class OneEzCmplx(OneElectricCmplx): pass
        
        
class OneMagneticCmplx(object):
    def __init__(self, idx, mu_r=1):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(mu_r) * mu0
        
    def update(self, hfield, efield1, efield2, d1, d2, dt, t):
        hfield[self.i, self.j, self.k] = 0
        
        
class OneHxCmplx(OneMagneticCmplx): pass
        
        
class OneHyCmplx(OneMagneticCmplx): pass
        

class OneHzCmplx(OneMagneticCmplx): pass
            
        
if __name__ == '__main__':
    pass
 