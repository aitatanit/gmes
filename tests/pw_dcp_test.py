#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np

from gmes.material import DcpPlrc, DcpAde, DrudePole, CriticalPoint
from gmes.constant import c0
from gmes.geometry import Cartesian


class Gold(DcpAde):
    """This gold permittivity value in the range of 200-1,000 nm.

    This parameters has a fitness value of 0.0907699 to the real 
    permittivity data in the range of 200-1,000 nm of 

    P. B. Johnson and R. W. Christy, "Optical constants of the
    noble metals," Phys. Rev. B 6, 4370 (1972).

    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.31839e16 * a / c0, 
                        gamma=1.09173e14 * a / c0)
        cp1 = CriticalPoint(amp=0.273222,
                            phi=-1.18299,
                            omega=3.88123e15 * a / c0, 
                            gamma=4.52006e14 * a / c0)
        cp2 = CriticalPoint(amp=3.04155,
                            phi=-1.09115,
                            omega=4.20737e15 * a / c0, 
                            gamma=2.35409e15 * a / c0)
        DcpAde.__init__(self, eps_inf=1.11683, mu_inf=1, sigma=0, 
                        dps=(dp1,), cps=(cp1,cp2))


class Silver(DcpPlrc):
    """This silver permittivity value in the range of 200-1,000 nm.
    
    This parameters has a fitness value of 0.0266134 to the real 
    permittivity data in the range of 200-1,000 nm of 
    
    P. B. Johnson and R. W. Christy, "Optical constants of the 
    noble metals,"  Phys. Rev. B 6, 4370 (1972).
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.38737e16 * a / c0, 
                        gamma=2.07331e13 * a / c0)
        cp1 = CriticalPoint(amp=1.3735,
                            phi=-0.504658,
                            omega=7.59914e15 * a / c0, 
                            gamma=4.28431e15 * a / c0)
        cp2 = CriticalPoint(amp=0.304478,
                            phi=-1.48944,
                            omega=6.15009e15 * a / c0, 
                            gamma=6.59262e14 * a / c0)
        DcpPlrc.__init__(self, eps_inf=0.89583, mu_inf=1, sigma=0,
                         dps=(dp1,), cps=(cp1,cp2))


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.spc = Cartesian((0, 0, 0))
        self.spc.dt = 1
        
        self.dcp_ade = Gold(a=1e-6)
        self.dcp_ade.init(self.spc)
        
        self.dcp_plrc = Silver(a=1e-6)
        self.dcp_plrc.init(self.spc)
        
    def testAdeExReal(self):
        sample = self.dcp_ade.get_pw_material_ex(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_ade.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = self.spc.dt
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ex[idx], 0)

    def testAdeEyReal(self):
        sample = self.dcp_ade.get_pw_material_ey(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_ade.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ey[idx], 0)
       
    def testAdeEzReal(self):
        sample = self.dcp_ade.get_pw_material_ez(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_ade.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ez[idx], 0)
    
    def testAdeHxReal(self):
        sample = self.dcp_ade.get_pw_material_hx(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_ade.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hx[idx], 0)
        
    def testAdeHyReal(self):
        sample = self.dcp_ade.get_pw_material_hy(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_ade.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hy[idx], 0)
       
    def testAdeHzReal(self):
        sample = self.dcp_ade.get_pw_material_hz(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_ade.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hz[idx], 0)

    def testAdeExCmplx(self):
        sample = self.dcp_ade.get_pw_material_ex(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_ade.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ex[idx], 0)
        
    def testAdeEyCmplx(self):
        sample = self.dcp_ade.get_pw_material_ey(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_ade.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ey[idx], 0)
        
    def testAdeEzCmplx(self):
        sample = self.dcp_ade.get_pw_material_ez(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_ade.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ez[idx], 0)
    
    def testAdeHxCmplx(self):
        sample = self.dcp_ade.get_pw_material_hx(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_ade.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hx[idx], 0)
        
    def testAdeHyCmplx(self):
        sample = self.dcp_ade.get_pw_material_hy(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_ade.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hy[idx], 0)
       
    def testAdeHzCmplx(self):
        sample = self.dcp_ade.get_pw_material_hz(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_ade.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hz[idx], 0)
        
    def testPlrcExReal(self):
        sample = self.dcp_plrc.get_pw_material_ex(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_plrc.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = self.spc.dt
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ex[idx], 0)

    def testPlrcEyReal(self):
        sample = self.dcp_plrc.get_pw_material_ey(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_plrc.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ey[idx], 0)
       
    def testPlrcEzReal(self):
        sample = self.dcp_plrc.get_pw_material_ez(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_plrc.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ez[idx], 0)
    
    def testPlrcHxReal(self):
        sample = self.dcp_plrc.get_pw_material_hx(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_plrc.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hx[idx], 0)
        
    def testPlrcHyReal(self):
        sample = self.dcp_plrc.get_pw_material_hy(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_plrc.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hy[idx], 0)
       
    def testPlrcHzReal(self):
        sample = self.dcp_plrc.get_pw_material_hz(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_plrc.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hz[idx], 0)

    def testPlrcExCmplx(self):
        sample = self.dcp_plrc.get_pw_material_ex(self.idx, (0,0,0), cmplx=True)

        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_plrc.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ex[idx], 0j)
        
    def testPlrcEyCmplx(self):
        sample = self.dcp_plrc.get_pw_material_ey(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_plrc.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ey[idx], 0j)
        
    def testPlrcEzCmplx(self):
        sample = self.dcp_plrc.get_pw_material_ez(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.dcp_plrc.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ez[idx], 0j)
    
    def testPlrcHxCmplx(self):
        sample = self.dcp_plrc.get_pw_material_hx(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_plrc.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hx[idx], 0j)
        
    def testPlrcHyCmplx(self):
        sample = self.dcp_plrc.get_pw_material_hy(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_plrc.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hy[idx], 0j)
       
    def testPlrcHzCmplx(self):
        sample = self.dcp_plrc.get_pw_material_hz(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.dcp_plrc.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hz[idx], 0j)

        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
