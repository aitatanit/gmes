#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np

from gmes.material import Gold2
from gmes.geometry import Cartesian
from gmes.pw_material import DrudeExReal, DrudeExCmplx
from gmes.pw_material import DrudeEyReal, DrudeEyCmplx
from gmes.pw_material import DrudeEzReal, DrudeEzCmplx
from gmes.pw_material import DrudeHxReal, DrudeHxCmplx
from gmes.pw_material import DrudeHyReal, DrudeHyCmplx
from gmes.pw_material import DrudeHzReal, DrudeHzCmplx


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1, 1, 1)
        self.idx2 = (1, 1, 1)
        while self.idx == self.idx2:
            self.idx2 = tuple(np.random.random_integers(0, 2) for x in range(3))

        self.spc = Cartesian((0, 0, 0))
        self.spc.dt = 1
        
        self.drude = Gold2(a=1)
        self.drude.init(self.spc)
        
    def testExReal(self):
        sample = DrudeExReal()
        sample_param = self.drude.get_ex_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = self.spc.dt
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)

    def testEyReal(self):
        sample = DrudeEyReal()
        sample_param = self.drude.get_ey_param(self.idx, (0,0,0), cmplx=False)

        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
       
    def testEzReal(self):
        sample = DrudeEzReal()
        sample_param = self.drude.get_ez_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testHxReal(self):
        sample = DrudeHxReal()
        sample_param = self.drude.get_hx_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testHyReal(self):
        sample = DrudeHyReal()
        sample_param = self.drude.get_hy_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testHzReal(self):
        sample = DrudeHzReal()
        sample_param = self.drude.get_hz_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)

    def testExCmplx(self):
        sample = DrudeExCmplx()
        sample_param = self.drude.get_ex_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)
        
    def testEyCmplx(self):
        sample = DrudeEyCmplx()
        sample_param = self.drude.get_ey_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
        
    def testEzCmplx(self):
        sample = DrudeEzCmplx()
        sample_param = self.drude.get_ez_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testHxCmplx(self):
        sample = DrudeHxCmplx()
        sample_param = self.drude.get_hx_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testHyCmplx(self):
        sample = DrudeHyCmplx()
        sample_param = self.drude.get_hy_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testHzCmplx(self):
        sample = DrudeHzCmplx()
        sample_param = self.drude.get_hz_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)


if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    
