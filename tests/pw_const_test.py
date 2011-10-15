#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np
from random import random

from gmes.material import Const
from gmes.geometry import Cartesian

    
class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)

        self.spc = Cartesian((0, 0, 0))
        self.spc.dt = 1
        
        self.const_real = Const(value=random(), 
                                eps_inf=random(), mu_inf=random())
        self.const_real.init(self.spc)

        self.const_cmplx = Const(value=random() + 1j * random(), 
                                 eps_inf=random(), mu_inf=random())
        self.const_cmplx.init(self.spc)

    def testExReal(self):
        sample = \
            self.const_real.get_pw_material_ex(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.const_real.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = self.spc.dt
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(ex[idx], self.const_real.value)
            else:
                self.assertEqual(ex[idx], 0)

    def testEyReal(self):
        sample = \
            self.const_real.get_pw_material_ey(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.const_real.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(ey[idx], self.const_real.value)
            else:
                self.assertEqual(ey[idx], 0)

    def testEzReal(self):
        sample = \
            self.const_real.get_pw_material_ez(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.const_real.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(ez[idx], self.const_real.value)
            else:
                self.assertEqual(ez[idx], 0)

    def testHxReal(self):
        sample = \
            self.const_real.get_pw_material_hx(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.const_real.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(hx[idx], self.const_real.value)
            else:
                self.assertEqual(hx[idx], 0)

    def testHyReal(self):
        sample = \
            self.const_real.get_pw_material_hy(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.const_real.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(hy[idx], self.const_real.value)
            else:
                self.assertEqual(hy[idx], 0)

    def testHzReal(self):
        sample = \
            self.const_real.get_pw_material_hz(self.idx, (0,0,0), cmplx=False)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.const_real.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(hz[idx], self.const_real.value)
            else:
                self.assertEqual(hz[idx], 0)

    def testExCmplx(self):
        sample = \
            self.const_cmplx.get_pw_material_ex(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.const_cmplx.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(ex[idx], self.const_cmplx.value)
            else:
                self.assertEqual(ex[idx], 0j)

    def testEyCmplx(self):
        sample = \
            self.const_cmplx.get_pw_material_ey(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.const_cmplx.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(ey[idx], self.const_cmplx.value)
            else:
                self.assertEqual(ey[idx], 0j)

    def testEzCmplx(self):
        sample = \
            self.const_cmplx.get_pw_material_ez(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps_inf(idx), self.const_cmplx.eps_inf)
            else:
                self.assertEqual(sample.get_eps_inf(idx), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(ez[idx], self.const_cmplx.value)
            else:
                self.assertEqual(ez[idx], 0j)

    def testHxCmplx(self):
        sample = \
            self.const_cmplx.get_pw_material_hx(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.const_cmplx.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(hx[idx], self.const_cmplx.value)
            else:
                self.assertEqual(hx[idx], 0j)

    def testHyCmplx(self):
        sample = \
            self.const_cmplx.get_pw_material_hy(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.const_cmplx.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(hy[idx], self.const_cmplx.value)
            else:
                self.assertEqual(hy[idx], 0j)

    def testHzCmplx(self):
        sample = \
            self.const_cmplx.get_pw_material_hz(self.idx, (0,0,0), cmplx=True)
        
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu_inf(idx), self.const_cmplx.mu_inf)
            else:
                self.assertEqual(sample.get_mu_inf(idx), 0)

        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(hz[idx], self.const_cmplx.value)
            else:
                self.assertEqual(hz[idx], 0j)


if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
