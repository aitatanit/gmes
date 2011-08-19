#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np

from gmes.material import GoldPlrc, GoldAde
from gmes.geometry import Cartesian


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.spc = Cartesian((0, 0, 0))
        self.spc.dt = 1
        
        self.dcp_ade = GoldAde(a=1)
        self.dcp_ade.init(self.spc)
        
        self.dcp_plrc = GoldPlrc(a=1e-6)
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
    
