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

class CpmlElectricCmplx(object):
    def __init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2, psi1=0, psi2=0):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_r) * epsilon0
        self.b1, self.b2 = float(b1), float(b2)
        self.c1, self.c2 = float(c1), float(c2)
        self.kappa1, self.kappa2 = float(kappa1), float(kappa2)
        self.psi1, self.psi2 = float(psi1), float(psi2)


class CpmlExCmplx(CpmlElectricCmplx):
    def __init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2):
        CpmlElectricCmplx.__init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2)
        
    def update(self, ex, hz, hy, dy, dz, dt, t):
        i, j, k = self.i, self.j, self.k
        
        self.psi1 = self.b1 * self.psi1 + self.c1 * (hz[i+1,j+1,k] - hz[i+1,j,k]) / dy;
        self.psi2 = self.b2 * self.psi2 + self.c2 * (hy[i+1,j,k+1] - hy[i+1,j,k]) / dz;

        ex[i,j,k] += dt / self.epsilon * ((hz[i+1,j+1,k] - hz[i+1,j,k]) / dy / self.kappa1 - 
                                          (hy[i+1,j,k+1] - hy[i+1,j,k]) / dz / self.kappa2 + 
                                          self.psi1 - self.psi2)
        
        
class CpmlEyCmplx(CpmlElectricCmplx):
    def __init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2):
        CpmlElectricCmplx.__init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2)
        
    def update(self, ey, hx, hz, dz, dx, dt, t):
        i, j, k = self.i, self.j, self.k
        
        self.psi1 = self.b1 * self.psi1 + self.c1 * (hx[i,j+1,k+1] - hx[i,j+1,k]) / dz
        self.psi2 = self.b2 * self.psi2 + self.c2 * (hz[i+1,j+1,k] - hz[i,j+1,k]) / dx

        ey[i,j,k] += dt / self.epsilon * ((hx[i,j+1,k+1] - hx[i,j+1,k]) / dz / self.kappa1 - 
                                          (hz[i+1,j+1,k] - hz[i,j+1,k]) / dx / self.kappa2 + 
                                          self.psi1 - self.psi2)
    
    
class CpmlEzCmplx(CpmlElectricCmplx):
    def __init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2):
        CpmlElectricCmplx.__init__(self, idx, epsilon_r, b1, b2, c1, c2, kappa1, kappa2)
        
    def update(self, ez, hy, hx, dx, dy, dt, t):
        i, j, k = self.i, self.j, self.k
        
        self.psi1 = self.b1 * self.psi1 + self.c1 * (hy[i+1,j,k+1] - hy[i,j,k+1]) / dx
        self.psi2 = self.b2 * self.psi2 + self.c2 * (hx[i,j+1,k+1] - hx[i,j,k+1]) / dy

        ez[i,j,k] += dt / self.epsilon * ((hy[i+1,j,k+1] - hy[i,j,k+1]) / dx / self.kappa1 - 
                                          (hx[i,j+1,k+1] - hx[i,j,k+1]) / dy / self.kappa2 + 
                                          self.psi1 - self.psi2)


class CpmlMagneticCmplx(object):
    def __init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2, psi1=0, psi2=0):
        self.i, self.j, self.k = array(idx, int)
        self.mu = float(mu_r) * mu0
        self.b1, self.b2 = float(b1), float(b2)
        self.c1, self.c2 = float(c1), float(c2)
        self.kappa1, self.kappa2 = float(kappa1), float(kappa2)
        self.psi1, self.psi2 = float(psi1), float(psi2)
        
        
class CpmlHxCmplx(CpmlMagneticCmplx):
    def __init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2):
        CpmlMagneticCmplx.__init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2)
        
    def update(self, hx, ez, ey, dy, dz, dt, t):
        i, j, k = self.i, self.j, self.k
        
        self.psi1 = self.b1 * self.psi1 + self.c1 * (ez[i,j,k-1] - ez[i,j-1,k-1]) / dy
        self.psi2 = self.b2 * self.psi2 + self.c2 * (ey[i,j-1,k] - ey[i,j-1,k-1]) / dz

        hx[i,j,k] -= dt / self.mu * ((ez[i,j,k-1] - ez[i,j-1,k-1]) / dy / self.kappa1 - 
                                     (ey[i,j-1,k] - ey[i,j-1,k-1]) / dz / self.kappa2 + 
                                     self.psi1 - self.psi2)
        

class CpmlHyCmplx(CpmlMagneticCmplx):
    def __init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2):
        CpmlMagneticCmplx.__init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2)
        
    def update(self, hy, ex, ez, dz, dx, dt, t):
        i, j, k = self.i, self.j, self.k
        
        self.psi1 = self.b1 * self.psi1 + self.c1 * (ex[i-1,j,k] - ex[i-1,j,k-1]) / dz
        self.psi2 = self.b2 * self.psi2 + self.c2 * (ez[i,j,k-1] - ez[i-1,j,k-1]) / dx

        hy[i,j,k] -= dt / self.mu * ((ex[i-1,j,k] - ex[i-1,j,k-1]) / dz / self.kappa1 - 
                                     (ez[i,j,k-1] - ez[i-1,j,k-1]) / dx / self.kappa2 + 
                                     self.psi1 - self.psi2)
        

class CpmlHzCmplx(CpmlMagneticCmplx):
    def __init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2):
        CpmlMagneticCmplx.__init__(self, idx, mu_r, b1, b2, c1, c2, kappa1, kappa2)
        
    def update(self, hz, ey, ex, dx, dy, dt, t):
        i, j, k = self.i, self.j, self.k
        
        self.psi1 = self.b1 * self.psi1 + self.c1 * (ey[i,j-1,k] - ey[i-1,j-1,k]) / dx
        self.psi2 = self.b2 * self.psi2 + self.c2 * (ex[i-1,j,k] - ex[i-1,j-1,k]) / dy

        hz[i,j,k] -= dt / self.mu * ((ey[i,j-1,k] - ey[i-1,j-1,k]) / dx / self.kappa1 - 
                                     (ex[i-1,j,k] - ex[i-1,j-1,k]) / dy / self.kappa2 + 
                                     self.psi1 - self.psi2)
        
        
if __name__ == '__main__':
    pass
