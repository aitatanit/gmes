#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np
from random import random

from gmes.pw_material import CpmlElectricParamReal, CpmlElectricParamCmplx
from gmes.pw_material import CpmlMagneticParamReal, CpmlMagneticParamCmplx
from gmes.pw_material import CpmlExReal, CpmlExCmplx
from gmes.pw_material import CpmlEyReal, CpmlEyCmplx
from gmes.pw_material import CpmlEzReal, CpmlEzCmplx
from gmes.pw_material import CpmlHxReal, CpmlHxCmplx
from gmes.pw_material import CpmlHyReal, CpmlHyCmplx
from gmes.pw_material import CpmlHzReal, CpmlHzCmplx


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.e_param_real = CpmlElectricParamReal()
        self.e_param_real.eps = random()
        self.e_param_real.b1 = self.e_param_real.b2 = 1
        self.e_param_real.c1 = self.e_param_real.c2 = 1
        self.e_param_real.kappa1 = self.e_param_real.kappa2 = 1
        
        self.e_param_cmplx = CpmlElectricParamCmplx()
        self.e_param_cmplx.eps = random()
        self.e_param_cmplx.b1 = self.e_param_cmplx.b2 = 1
        self.e_param_cmplx.c1 = self.e_param_cmplx.c2 = 1
        self.e_param_cmplx.kappa1 = self.e_param_cmplx.kappa2 = 1

        self.h_param_real = CpmlMagneticParamReal()
        self.h_param_real.mu = random()
        self.h_param_real.b1 = self.h_param_real.b2 = 1
        self.h_param_real.c1 = self.h_param_real.c2 = 1
        self.h_param_real.kappa1 = self.h_param_real.kappa2 = 1

        self.h_param_cmplx = CpmlMagneticParamCmplx()
        self.h_param_cmplx.mu = random()
        self.h_param_cmplx.b1 = self.h_param_cmplx.b2 = 1
        self.h_param_cmplx.c1 = self.h_param_cmplx.c2 = 1
        self.h_param_cmplx.kappa1 = self.h_param_cmplx.kappa2 = 1

    def testExReal(self):
        sample = CpmlExReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)

    def testEyReal(self):
        sample = CpmlEyReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
       
    def testEzReal(self):
        sample = CpmlEzReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testHxReal(self):
        sample = CpmlHxReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testHyReal(self):
        sample = CpmlHyReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testHzReal(self):
        sample = CpmlHzReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)

    def testExCmplx(self):
        sample = CpmlExCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)
        
    def testEyCmplx(self):
        sample = CpmlEyCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
        
    def testEzCmplx(self):
        sample = CpmlEzCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testHxCmplx(self):
        sample = CpmlHxCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testHyCmplx(self):
        sample = CpmlHyCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testHzCmplx(self):
        sample = CpmlHzCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
