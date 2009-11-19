#!/usr/bin/env python

# This implementation is based on the following article.
# M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
# Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006.

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from pw_dielectric_cmplx import *
from gmes.constants import *


class DrudeElectricCmplx(object):
    def __init__(self, idx, epsilon_inf, omega_p, gamma_p):
        self.i, self.j, self.k = array(idx, int)
        self.epsilon = float(epsilon_inf) * epsilon0
        self.omega_p = array(omega_p, float)
        self.gamma_p = array(gamma_p, float)
        self.q_new = zeros(self.omega_p.size, complex)
        self.q_old = zeros(self.omega_p.size, complex)
        
        
class DrudeExCmplx(DrudeElectricCmplx):
    def __init__(self, idx, epsilon_inf, omega_p, gamma_p):
        DrudeElectricCmplx.__init__(self, idx, epsilon_inf, omega_p, gamma_p)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        i, j, k = self.i, self.j, self.k
        
        q_tmp = (4 * self.q_new + (self.gamma_p * dt - 2) * self.q_old - 
                 (2 * dt**2 * epsilon0 * self.omega_p**2) * ex[i,j,k]) / (self.gamma_p * dt + 2)
        
        self.q_old[:] = self.q_new[:]
        self.q_new[:] = q_tmp[:]

        q_diff_sum = self.q_new.sum() - self.q_old.sum()
    
        ex[i,j,k] += (dt * ((hz[i+1,j+1,k] - hz[i+1,j,k]) / dy - 
                            (hy[i+1,j,k+1] - hy[i+1,j,k]) / dz) + q_diff_sum) / self.epsilon
    
    
class DrudeEyCmplx(DrudeElectricCmplx):
    def __init__(self, idx, epsilon_inf, omega_p, gamma_p):
        DrudeElectricCmplx.__init__(self, idx, epsilon_inf, omega_p, gamma_p)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        i, j, k = self.i, self.j, self.k
        
        q_tmp = (4 * self.q_new + (self.gamma_p * dt - 2) * self.q_old - 
                 (2 * dt**2 * epsilon0 * self.omega_p**2) * ey[i,j,k]) / (self.gamma_p * dt + 2)
                 
        self.q_old[:] = self.q_new[:]
        self.q_new[:] = q_tmp[:]
        
        q_diff_sum = self.q_new.sum() - self.q_old.sum()
        
        ey[i,j,k] += (dt * ((hx[i,j+1,k+1] - hx[i,j+1,k]) / dz - 
                            (hz[i+1,j+1,k] - hz[i,j+1,k]) / dx) + q_diff_sum) / self.epsilon
    
    
class DrudeEzCmplx(DrudeElectricCmplx):
    def __init__(self, idx, epsilon_inf, omega_p, gamma_p):
        DrudeElectricCmplx.__init__(self, idx, epsilon_inf, omega_p, gamma_p)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        i, j, k = self.i, self.j, self.k
        
        q_tmp = (4 * self.q_new + (self.gamma_p * dt - 2) * self.q_old - 
                 (2 * dt**2 * epsilon0 * self.omega_p**2) * ez[i,j,k]) / (self.gamma_p * dt + 2)
    
        self.q_old[:] = self.q_new[:]
        self.q_new[:] = q_tmp[:]

        q_diff_sum = self.q_new.sum() - self.q_old.sum()

        ez[i,j,k] += (dt * ((hy[i+1,j,k+1] - hy[i,j,k+1]) / dx - 
                            (hx[i,j+1,k+1] - hx[i,j,k+1]) / dy) + q_diff_sum) / self.epsilon
    
    
class DrudeHxCmplx(DielectricHxCmplx):
    def __init__(self, idx, mu_r=1):
        DielectricHxCmplx.__init__(self, idx, mu_r)


class DrudeHyCmplx(DielectricHyCmplx):
    def __init__(self, idx, mu_r=1):
        DielectricHyCmplx.__init__(self, idx, mu_r)
        
        
class DrudeHzCmplx(DielectricHzCmplx):
    def __init__(self, idx, mu_r=1):
        DielectricHzCmplx.__init__(self, idx, mu_r)
        
        
if __name__ == '__main__':
    pass
