#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np
from random import random

from gmes.pw_material import DielectricElectricParamReal
from gmes.pw_material import DielectricElectricParamCmplx
from gmes.pw_material import DielectricMagneticParamReal
from gmes.pw_material import DielectricMagneticParamCmplx
from gmes.pw_material import DielectricExReal, DielectricExCmplx
from gmes.pw_material import DielectricEyReal, DielectricEyCmplx
from gmes.pw_material import DielectricEzReal, DielectricEzCmplx
from gmes.pw_material import DielectricHxReal, DielectricHxCmplx
from gmes.pw_material import DielectricHyReal, DielectricHyCmplx
from gmes.pw_material import DielectricHzReal, DielectricHzCmplx
    

class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.e_param_real = DielectricElectricParamReal()
        self.e_param_real.eps = random()

        self.e_param_cmplx = DielectricElectricParamCmplx()
        self.e_param_cmplx.eps = random()

        self.h_param_real = DielectricMagneticParamReal()
        self.h_param_real.mu = random()

        self.h_param_cmplx = DielectricMagneticParamCmplx()
        self.h_param_cmplx.mu = random()
        
    def testExReal(self):
        sample = DielectricExReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)

    def testEyReal(self):
        sample = DielectricEyReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
       
    def testEzReal(self):
        sample = DielectricEzReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testHxReal(self):
        sample = DielectricHxReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testHyReal(self):
        sample = DielectricHyReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testHzReal(self):
        sample = DielectricHzReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)

    def testExCmplx(self):
        sample = DielectricExCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)
        
    def testEyCmplx(self):
        sample = DielectricEyCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
        
    def testEzCmplx(self):
        sample = DielectricEzCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testHxCmplx(self):
        sample = DielectricHxCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testHyCmplx(self):
        sample = DielectricHyCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testHzCmplx(self):
        sample = DielectricHzCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
