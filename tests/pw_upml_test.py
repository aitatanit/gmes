#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np
from random import random

from gmes.pw_material import UpmlElectricParamReal, UpmlElectricParamCmplx
from gmes.pw_material import UpmlMagneticParamReal, UpmlMagneticParamCmplx
from gmes.pw_material import UpmlExReal, UpmlExCmplx
from gmes.pw_material import UpmlEyReal, UpmlEyCmplx
from gmes.pw_material import UpmlEzReal, UpmlEzCmplx
from gmes.pw_material import UpmlHxReal, UpmlHxCmplx
from gmes.pw_material import UpmlHyReal, UpmlHyCmplx
from gmes.pw_material import UpmlHzReal, UpmlHzCmplx


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.e_param_real = UpmlElectricParamReal()
        self.e_param_real.eps = random()
        self.e_param_real.c1 = self.e_param_real.c2 = 1
        self.e_param_real.c3 = self.e_param_real.c4 = 1
        self.e_param_real.c5 = self.e_param_real.c6 = 1
        
        self.e_param_cmplx = UpmlElectricParamCmplx()
        self.e_param_cmplx.eps = random()
        self.e_param_cmplx.c1 = self.e_param_cmplx.c2 = 1
        self.e_param_cmplx.c3 = self.e_param_cmplx.c4 = 1
        self.e_param_cmplx.c5 = self.e_param_cmplx.c6 = 1

        self.h_param_real = UpmlMagneticParamReal()
        self.h_param_real.mu = random()
        self.h_param_real.c1 = self.h_param_real.c2 = 1
        self.h_param_real.c3 = self.h_param_real.c4 = 1
        self.h_param_real.c5 = self.h_param_real.c6 = 1

        self.h_param_cmplx = UpmlMagneticParamCmplx()
        self.h_param_cmplx.mu = random()
        self.h_param_cmplx.c1 = self.h_param_cmplx.c2 = 1
        self.h_param_cmplx.c3 = self.h_param_cmplx.c4 = 1
        self.h_param_cmplx.c5 = self.h_param_cmplx.c6 = 1

    def testExReal(self):
        sample = UpmlExReal()
        sample.attach(self.idx, self.e_param_real)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps(idx), self.e_param_real.eps)
            else:
                self.assertEqual(sample.get_eps(idx), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ex[idx], 0)

    def testEyReal(self):
        sample = UpmlEyReal()
        sample.attach(self.idx, self.e_param_real)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps(idx), self.e_param_real.eps)
            else:
                self.assertEqual(sample.get_eps(idx), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ey[idx], 0)

    def testEzReal(self):
        sample = UpmlEzReal()
        sample.attach(self.idx, self.e_param_real)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps(idx), self.e_param_real.eps)
            else:
                self.assertEqual(sample.get_eps(idx), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ez[idx], 0)
    
    def testHxReal(self):
        sample = UpmlHxReal()
        sample.attach(self.idx, self.h_param_real)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu(idx), self.h_param_real.mu)
            else:
                self.assertEqual(sample.get_mu(idx), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hx[idx], 0)
        
    def testHyReal(self):
        sample = UpmlHyReal()
        sample.attach(self.idx, self.h_param_real)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu(idx), self.h_param_real.mu)
            else:
                self.assertEqual(sample.get_mu(idx), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hy[idx], 0)
       
    def testHzReal(self):
        sample = UpmlHzReal()
        sample.attach(self.idx, self.h_param_real)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu(idx), self.h_param_real.mu)
            else:
                self.assertEqual(sample.get_mu(idx), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hz[idx], 0)

    def testExCmplx(self):
        sample = UpmlExCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps(idx), self.e_param_cmplx.eps)
            else:
                self.assertEqual(sample.get_eps(idx), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ex[idx], 0)
        
    def testEyCmplx(self):
        sample = UpmlEyCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps(idx), self.e_param_cmplx.eps)
            else:
                self.assertEqual(sample.get_eps(idx), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ey[idx], 0)
        
    def testEzCmplx(self):
        sample = UpmlEzCmplx()
        sample.attach(self.idx, self.e_param_cmplx)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_eps(idx), self.e_param_cmplx.eps)
            else:
                self.assertEqual(sample.get_eps(idx), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(ez[idx], 0)
    
    def testHxCmplx(self):
        sample = UpmlHxCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu(idx), self.h_param_cmplx.mu)
            else:
                self.assertEqual(sample.get_mu(idx), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hx[idx], 0)
        
    def testHyCmplx(self):
        sample = UpmlHyCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu(idx), self.h_param_cmplx.mu)
            else:
                self.assertEqual(sample.get_mu(idx), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hy[idx], 0)
       
    def testHzCmplx(self):
        sample = UpmlHzCmplx()
        sample.attach(self.idx, self.h_param_cmplx)
        for idx in np.ndindex(3, 3, 3):
            if idx == self.idx:
                self.assertEqual(sample.get_mu(idx), self.h_param_cmplx.mu)
            else:
                self.assertEqual(sample.get_mu(idx), 0)

        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        for idx in np.ndindex(3, 3, 3):
            self.assertEqual(hz[idx], 0)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
