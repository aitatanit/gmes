#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
import numpy as np

from gmes.material import GoldPlrc, GoldAde
from gmes.geometry import Cartesian

# ADE implementation
from gmes.pw_material import DcpAdeExReal, DcpAdeExCmplx
from gmes.pw_material import DcpAdeEyReal, DcpAdeEyCmplx
from gmes.pw_material import DcpAdeEzReal, DcpAdeEzCmplx
from gmes.pw_material import DcpAdeHxReal, DcpAdeHxCmplx
from gmes.pw_material import DcpAdeHyReal, DcpAdeHyCmplx
from gmes.pw_material import DcpAdeHzReal, DcpAdeHzCmplx

# PLRC implementation
from gmes.pw_material import DcpPlrcExReal, DcpPlrcExCmplx
from gmes.pw_material import DcpPlrcEyReal, DcpPlrcEyCmplx
from gmes.pw_material import DcpPlrcEzReal, DcpPlrcEzCmplx
from gmes.pw_material import DcpPlrcHxReal, DcpPlrcHxCmplx
from gmes.pw_material import DcpPlrcHyReal, DcpPlrcHyCmplx
from gmes.pw_material import DcpPlrcHzReal, DcpPlrcHzCmplx


class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        self.idx2 = (1, 1, 1)
        while self.idx == self.idx2:
            self.idx2 = tuple(np.random.random_integers(0, 2) for x in range(3))
        
        self.spc = Cartesian((0, 0, 0))
        self.spc.dt = 1
        
        self.dcp_ade = GoldAde(a=1)
        self.dcp_ade.init(self.spc)
        
        self.dcp_plrc = GoldPlrc(a=1e-6)
        self.dcp_plrc.init(self.spc)
        
    def testAdeExReal(self):
        sample = DcpAdeExReal()
        sample_param = self.dcp_ade.get_ex_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = self.spc.dt
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)

    def testAdeEyReal(self):
        sample = DcpAdeEyReal()
        sample_param = self.dcp_ade.get_ey_param(self.idx, (0,0,0), cmplx=False)

        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
       
    def testAdeEzReal(self):
        sample = DcpAdeEzReal()
        sample_param = self.dcp_ade.get_ez_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testAdeHxReal(self):
        sample = DcpAdeHxReal()
        sample_param = self.dcp_ade.get_hx_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testAdeHyReal(self):
        sample = DcpAdeHyReal()
        sample_param = self.dcp_ade.get_hy_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testAdeHzReal(self):
        sample = DcpAdeHzReal()
        sample_param = self.dcp_ade.get_hz_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)

    def testAdeExCmplx(self):
        sample = DcpAdeExCmplx()
        sample_param = self.dcp_ade.get_ex_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)
        
    def testAdeEyCmplx(self):
        sample = DcpAdeEyCmplx()
        sample_param = self.dcp_ade.get_ey_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
        
    def testAdeEzCmplx(self):
        sample = DcpAdeEzCmplx()
        sample_param = self.dcp_ade.get_ez_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testAdeHxCmplx(self):
        sample = DcpAdeHxCmplx()
        sample_param = self.dcp_ade.get_hx_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testAdeHyCmplx(self):
        sample = DcpAdeHyCmplx()
        sample_param = self.dcp_ade.get_hy_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testAdeHzCmplx(self):
        sample = DcpAdeHzCmplx()
        sample_param = self.dcp_ade.get_hz_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hz = ey = ex = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)
        
    def testPlrcExReal(self):
        sample = DcpPlrcExReal()
        sample_param = self.dcp_plrc.get_ex_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ex = hz = hy = np.zeros((3,3,3))
        dy = dz = dt = self.spc.dt
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)

    def testPlrcEyReal(self):
        sample = DcpPlrcEyReal()
        sample_param = self.dcp_plrc.get_ey_param(self.idx, (0,0,0), cmplx=False)

        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ey = hx = hz = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
       
    def testPlrcEzReal(self):
        sample = DcpPlrcEzReal()
        sample_param = self.dcp_plrc.get_ez_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ez = hy = hx = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testPlrcHxReal(self):
        sample = DcpPlrcHxReal()
        sample_param = self.dcp_plrc.get_hx_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hx = ez = ey = np.zeros((3,3,3))
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testPlrcHyReal(self):
        sample = DcpPlrcHyReal()
        sample_param = self.dcp_plrc.get_hy_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hy = ex = ez = np.zeros((3,3,3))
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testPlrcHzReal(self):
        sample = DcpPlrcHzReal()
        sample_param = self.dcp_plrc.get_hz_param(self.idx, (0,0,0), cmplx=False)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hz = ey = ex = np.zeros((3,3,3))
        dx = dy = dt = 1
        n = 0
        sample.update_all(hz, ey, ex, dx, dy, dt, n)
        self.assertEqual(hz[self.idx], 0)

    def testPlrcExCmplx(self):
        sample = DcpPlrcExCmplx()
        sample_param = self.dcp_plrc.get_ex_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ex = hz = hy = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(ex, hz, hy, dy, dz, dt, n)
        self.assertEqual(ex[self.idx], 0)
        
    def testPlrcEyCmplx(self):
        sample = DcpPlrcEyCmplx()
        sample_param = self.dcp_plrc.get_ey_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ey = hx = hz = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(ey, hx, hz, dz, dx, dt, n)
        self.assertEqual(ey[self.idx], 0)
        
    def testPlrcEzCmplx(self):
        sample = DcpPlrcEzCmplx()
        sample_param = self.dcp_plrc.get_ez_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_eps(self.idx), sample_param.eps)
        self.assertEqual(sample.get_eps(self.idx2), 0)

        ez = hy = hx = np.zeros((3,3,3), complex)
        dx = dy = dt = 1
        n = 0
        sample.update_all(ez, hy, hx, dx, dy, dt, n)
        self.assertEqual(ez[self.idx], 0)
    
    def testPlrcHxCmplx(self):
        sample = DcpPlrcHxCmplx()
        sample_param = self.dcp_plrc.get_hx_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hx = ez = ey = np.zeros((3,3,3), complex)
        dy = dz = dt = 1
        n = 0
        sample.update_all(hx, ez, ey, dy, dz, dt, n)
        self.assertEqual(hx[self.idx], 0)
        
    def testPlrcHyCmplx(self):
        sample = DcpPlrcHyCmplx()
        sample_param = self.dcp_plrc.get_hy_param(self.idx, (0,0,0), cmplx=True)
        sample.attach(self.idx, sample_param)
        self.assertEqual(sample.get_mu(self.idx), sample_param.mu)
        self.assertEqual(sample.get_mu(self.idx2), 0)

        hy = ex = ez = np.zeros((3,3,3), complex)
        dz = dx = dt = 1
        n = 0
        sample.update_all(hy, ex, ez, dz, dx, dt, n)
        self.assertEqual(hy[self.idx], 0)
       
    def testPlrcHzCmplx(self):
        sample = DcpPlrcHzCmplx()
        sample_param = self.dcp_plrc.get_hz_param(self.idx, (0,0,0), cmplx=True)
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
    
