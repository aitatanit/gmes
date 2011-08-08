#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np
from random import random

from gmes.pw_material import DummyElectricParamReal, DummyElectricParamCmplx
from gmes.pw_material import DummyMagneticParamReal, DummyMagneticParamCmplx
from gmes.pw_material import DummyExReal, DummyExCmplx
from gmes.pw_material import DummyEyReal, DummyEyCmplx
from gmes.pw_material import DummyEzReal, DummyEzCmplx
from gmes.pw_material import DummyHxReal, DummyHxCmplx
from gmes.pw_material import DummyHyReal, DummyHyCmplx
from gmes.pw_material import DummyHzReal, DummyHzCmplx


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.e_param_real = DummyElectricParamReal()
        self.e_param_real.eps = random()

        self.e_param_cmplx = DummyElectricParamCmplx()
        self.e_param_cmplx.eps = random()

        self.h_param_real = DummyMagneticParamReal()
        self.h_param_real.mu = random()

        self.h_param_cmplx = DummyMagneticParamCmplx()
        self.h_param_cmplx.mu = random()
        
    def testExReal(self):
        sample = DummyExReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)

        ex = np.random.random((3,3,3))
        hz = hy = np.empty((3,3,3))
        dy = dz = dt = 1
        n = 0
        value = ex[self.idx]
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], value)
        
    def testEyReal(self):
        sample = DummyEyReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ey = np.random.random((3,3,3))
        hx = hz = np.empty((3,3,3))
        dz = dx = dt = 1
        n = 0
        value = ey[self.idx]
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], value)

    def testEzReal(self):
        sample = DummyEzReal()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)
        
        ez = np.random.random((3,3,3))
        hy = hx = np.empty((3,3,3))
        dx = dy = dt = 1
        n = 0
        value = ez[self.idx]
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], value)
        
    def testHxReal(self):
        sample = DummyHxReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hx = np.random.random((3,3,3))
        ez = ey = np.empty((3,3,3))
        dy = dz = dt = 1
        n = 0
        value = hx[self.idx]
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], value)
        
    def testHyReal(self):
        sample = DummyHyReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)

        hy = np.random.random((3,3,3))
        ex = ez = np.empty((3,3,3))
        dz = dx = dt = 1
        n = 0
        value = hy[self.idx]
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], value)

    def testHzReal(self):
        sample = DummyHzReal()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)
        
        hz = np.random.random((3,3,3))
        ey = ex = np.empty((3,3,3))
        dx = dy = dt = 1
        n = 0
        value = hz[self.idx]
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], value)

    def testExCmplx(self):
        sample = DummyExCmplx()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)

        ex = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        hz = hy = np.empty((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        value = ex[self.idx]
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], value)
        
    def testEyCmplx(self):
        sample = DummyEyCmplx()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)

        ey = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        hx = hz = np.empty((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        value = ey[self.idx]
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], value)

    def testEzCmplx(self):
        sample = DummyEzCmplx()
        sample.attach(self.idx, self.e_param_real)
        self.assertEqual(sample.get_eps(self.idx), self.e_param_real.eps)

        ez = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        hy = hx = np.empty((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        value = ez[self.idx]
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], value)
        
    def testHxCmplx(self):
        sample = DummyHxCmplx()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)

        hx = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        ez = ey = np.empty((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        value = hx[self.idx]
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], value)
        
    def testHyCmplx(self):
        sample = DummyHyCmplx()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)

        hy = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        ex = ez = np.empty((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        value = hy[self.idx]
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], value)

    def testHzCmplx(self):
        sample = DummyHzCmplx()
        sample.attach(self.idx, self.h_param_real)
        self.assertEqual(sample.get_mu(self.idx), self.h_param_real.mu)

        hz = np.random.random((3,3,3)) + 1j * np.random.random((3,3,3))
        ey = ex = np.empty((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        value = hz[self.idx]
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], value)


if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
