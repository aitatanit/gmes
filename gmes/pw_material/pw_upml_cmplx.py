#!/usr/bin/env python

# This implementation is based on the following article.
# S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,"
# Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
# Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
# 2005, pp. 273-328.

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from gmes.constants import *


class UpmlElectricCmplx(object):
    def __init__(self, idx, epsilon_r, c1, c2, c3, c4, c5, c6):
        self.i, self.j, self.k = array(idx, int)
        
        self.epsilon = float(epsilon_r) * epsilon0
        
        self.c1, self.c2, self.c3 = float(c1), float(c2), float(c3)
        self.c4, self.c5, self.c6 = float(c4), float(c5), float(c6)
        
        self.d = 0


class UpmlExCmplx(UpmlElectricCmplx):
    def update(self, ex, hz, hy, dy, dz, dt, t):
        i, j, k = self.i, self.j, self.k
        
        dstore = self.d
        self.d = self.c1 * self.d + self.c2 * ((hz[i+1,j+1,k] - hz[i+1,j,k]) / dy - 
                                               (hy[i+1,j,k+1] - hy[i+1,j,k]) / dz)
        ex[i,j,k] = self.c3 * ex[i,j,k] + \
        self.c4 * (self.c5 * self.d - self.c6 * dstore) / self.epsilon
    
        
class UpmlEyCmplx(UpmlElectricCmplx):
    def update(self, ey, hx, hz, dz, dx, dt, t):
        i, j, k = self.i, self.j, self.k
        
        dstore = self.d
        self.d = self.c1 * self.d + self.c2 * ((hx[i,j+1,k+1] - hx[i,j+1,k]) / dz - 
                                               (hz[i+1,j+1,k] - hz[i,j+1,k]) / dx)        
        ey[i,j,k] = self.c3 * ey[i,j,k] + \
        self.c4 * (self.c5 * self.d - self.c6 * dstore) / self.epsilon
    
    
class UpmlEzCmplx(UpmlElectricCmplx):
    def update(self, ez, hy, hx, dx, dy, dt, t):
        i, j, k = self.i, self.j, self.k
        
        dstore = self.d
        self.d = self.c1 * self.d + self.c2 * ((hy[i+1,j,k+1] - hy[i,j,k+1]) / dx - 
                                               (hx[i,j+1,k+1] - hx[i,j,k+1]) / dy)
        ez[i,j,k] = self.c3 * ez[i,j,k] + \
        self.c4 * (self.c5 * self.d - self.c6 * dstore) / self.epsilon
    
    
class UpmlMagneticCmplx(object):
    def __init__(self, idx, mu_r, c1, c2, c3, c4, c5, c6):
        self.i, self.j, self.k = array(idx, int)
        
        self.mu = float(mu_r) * mu0
        
        self.c1, self.c2, self.c3 = float(c1), float(c2), float(c3)
        self.c4, self.c5, self.c6 = float(c4), float(c5), float(c6)
        
        self.b = 0
        
        
class UpmlHxCmplx(UpmlMagneticCmplx):
    def update(self, hx, ez, ey, dy, dz, dt, t):
        i, j, k = self.i, self.j, self.k
        
        bstore = self.b
        self.b = self.c1 * self.b - self.c2 * ((ez[i,j,k-1] - ez[i,j-1,k-1]) / dy - 
                                               (ey[i,j-1,k] - ey[i,j-1,k-1]) / dz)
        hx[i,j,k] = self.c3 * hx[i,j,k] + \
        self.c4 * (self.c5 * self.b - self.c6 * bstore) / self.mu
    
    
class UpmlHyCmplx(UpmlMagneticCmplx):
    def update(self, hy, ex, ez, dz, dx, dt, t):
        i, j, k = self.i, self.j, self.k
        
        bstore = self.b
        self.b = self.c1 * self.b - self.c2 * ((ex[i-1,j,k] - ex[i-1,j,k-1]) / dz - 
                                               (ez[i,j,k-1] - ez[i-1,j,k-1]) / dx)
        hy[i,j,k] = self.c3 * hy[i,j,k] + \
        self.c4 * (self.c5 * self.b - self.c6 * bstore) / self.mu
    
    
class UpmlHzCmplx(UpmlMagneticCmplx):
    def update(self, hz, ey, ex, dx, dy, dt, t):
        i, j, k = self.i, self.j, self.k
        
        bstore = self.b
        self.b = self.c1 * self.b - self.c2 * ((ey[i,j-1,k] - ey[i-1,j-1,k]) / dx - 
                                               (ex[i-1,j,k] - ex[i-1,j-1,k]) / dy)
        hz[i,j,k] = self.c3 * hz[i,j,k] + \
        self.c4 * (self.c5 * self.b - self.c6 * bstore) / self.mu
    
    
if __name__ == '__main__':
    pass
