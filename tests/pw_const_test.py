#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np
from random import random

from gmes.pw_material import ConstElectricParamReal, ConstElectricParamCmplx
from gmes.pw_material import ConstMagneticParamReal, ConstMagneticParamCmplx
from gmes.pw_material import ConstExReal, ConstExCmplx
from gmes.pw_material import ConstEyReal, ConstEyCmplx
from gmes.pw_material import ConstEzReal, ConstEzCmplx
from gmes.pw_material import ConstHxReal, ConstHxCmplx
from gmes.pw_material import ConstHyReal, ConstHyCmplx
from gmes.pw_material import ConstHzReal, ConstHzCmplx

    
class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.e_param_real = ConstElectricParamReal()
        self.e_param_real.eps = random()
        self.e_param_real.value = random()

        self.e_param_cmplx = ConstElectricParamCmplx()
        self.e_param_cmplx.eps = random()
        self.e_param_cmplx.value = random() + 1j * random()

        self.h_param_real = ConstMagneticParamReal()
        self.h_param_real.mu = random()
        self.h_param_real.value = random()

        self.h_param_cmplx = ConstMagneticParamCmplx()
        self.h_param_cmplx.mu = random()
        self.h_param_cmplx.value = random() + 1j * random()

    def testExReal(self):
        sample = ConstExReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ex = np.random.random((3,3,3))
        hz = hy = np.empty((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], self.e_param_real.value)
        
    def testEyReal(self):
        sample = ConstEyReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ey = np.random.random((3,3,3))
        hx = hz = np.empty((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], self.e_param_real.value)
       
    def testEzReal(self):
        sample = ConstEzReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ez = np.random.random((3,3,3))
        hy = hx = np.empty((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], self.e_param_real.value)
    
    def testHxReal(self):
        sample = ConstHxReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hx = np.random.random((3,3,3))
        ez = ey = np.empty((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], self.h_param_real.value)
        
    def testHyReal(self):
        sample = ConstHyReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hy = np.random.random((3,3,3))
        ex = ez = np.empty((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], self.h_param_real.value)
       
    def testHzReal(self):
        sample = ConstHzReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hz = np.random.random((3,3,3))
        ey = ex = np.empty((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], self.h_param_real.value)

    def testExCmplx(self):
        sample = ConstExCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ex = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        hz = hy = np.empty((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], self.e_param_cmplx.value)
        
    def testEyCmplx(self):
        sample = ConstEyCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ey = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        hx = hz = np.empty((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], self.e_param_cmplx.value)
        
    def testEzCmplx(self):
        sample = ConstEzCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_cmplx.eps)
        
        ez = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        hy = hx = np.empty((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], self.e_param_cmplx.value)
    
    def testHxCmplx(self):
        sample = ConstHxCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hx = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        ez = ey = np.empty((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], self.h_param_cmplx.value)
        
    def testHyCmplx(self):
        sample = ConstHyCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hy = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        ex = ez = np.empty((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], self.h_param_cmplx.value)
       
    def testHzCmplx(self):
        sample = ConstHzCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_cmplx.mu)
        
        hz = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        ey = ex = np.empty((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], self.h_param_cmplx.value)


if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
