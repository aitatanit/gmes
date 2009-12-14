#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from gmes.constants import *


class DielectricElectricCmplx(object):
    def __init__(self, idx, epsilon_r):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_r) * epsilon0
        
        
class DielectricExCmplx(DielectricElectricCmplx):
    def update(self, ex, hz, hy, dy, dz, dt, t):
        i, j, k = self.i, self.j, self.k
        
        ex[i,j,k] += dt / self.epsilon * ((hz[i+1,j+1,k] - hz[i+1,j,k]) / dy -
                                          (hy[i+1,j,k+1] - hy[i+1,j,k]) / dz)


class DielectricEyCmplx(DielectricElectricCmplx):
    def update(self, ey, hx, hz, dz, dx, dt, t):
        i, j, k = self.i, self.j, self.k
        
        ey[i,j,k] += dt / self.epsilon * ((hx[i,j+1,k+1] - hx[i,j+1,k]) / dz -
                                          (hz[i+1,j+1,k] - hz[i,j+1,k]) / dx)
        

class DielectricEzCmplx(DielectricElectricCmplx):
    def update(self, ez, hy, hx, dx, dy, dt, t):
        i, j, k = self.i, self.j, self.k

        ez[i,j,k] += dt / self.epsilon * ((hy[i+1,j,k+1] - hy[i,j,k+1]) / dx - 
                                          (hx[i,j+1,k+1] - hx[i,j,k+1]) / dy)
            
            
class DielectricMagneticCmplx(object):
    def __init__(self, idx, mu_r):
        self.i, self.j, self.k = array(idx, int)
        self.mu = float(mu_r) * mu0
        
        
class DielectricHxCmplx(DielectricMagneticCmplx):    
    def update(self, hx, ez, ey, dy, dz, dt, t):
        i, j, k = self.i, self.j, self.k

        hx[i,j,k] -= dt / self.mu * ((ez[i,j,k-1] - ez[i,j-1,k-1]) / dy - 
                                     (ey[i,j-1,k] - ey[i,j-1,k-1]) / dz)


class DielectricHyCmplx(DielectricMagneticCmplx):     
    def update(self, hy, ex, ez, dz, dx, dt, t):
        i, j, k = self.i, self.j, self.k

        hy[i,j,k] -= dt / self.mu * ((ex[i-1,j,k] - ex[i-1,j,k-1]) / dz -
                                     (ez[i,j,k-1] - ez[i-1,j,k-1]) / dx)


class DielectricHzCmplx(DielectricMagneticCmplx):    
    def update(self, hz, ey, ex, dx, dy, dt, t):
        i, j, k = self.i, self.j, self.k

        hz[i,j,k] -= dt / self.mu * ((ey[i,j-1,k] - ey[i-1,j-1,k]) / dx -
                                     (ex[i-1,j,k] - ex[i-1,j-1,k]) / dy)
        

if __name__ == '__main__':
    pass
    