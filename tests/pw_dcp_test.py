#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
from numpy import *
from gmes.material import *
from gmes.geometry import *

class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.realA = ones((3,3,3), float)
        self.realB = ones((3,3,3), float)
        self.realC = ones((3,3,3), float)
        
        self.cmplxA = ones((3,3,3), complex)
        self.cmplxB = ones((3,3,3), complex)
        self.cmplxC = ones((3,3,3), complex)
        
        self.diff = 1
        self.n = 0
        
        self.spc = Cartesian((0,0,0))
        self.spc.dt = 1
        
        self.dcp = DCP(dps=(DrudePole(omega=1,gamma=1),), cps=(CriticalPoint(amp=1,phi=1,omega=1,gamma=1),))
        self.dcp.init(self.spc)
        
    def testEx(self):
        sampleReal = self.dcp.get_pw_material_ex(self.idx, (0,0,0), cmplx=False)
        sampleCmplx = self.dcp.get_pw_material_ex(self.idx, (0,0,0), cmplx=True)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.n)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sampleReal = self.dcp.get_pw_material_ey(self.idx, (0,0,0), cmplx=False)
        sampleCmplx = self.dcp.get_pw_material_ey(self.idx, (0,0,0), cmplx=True)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.n)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sampleReal = self.dcp.get_pw_material_ez(self.idx, (0,0,0), cmplx=False)
        sampleCmplx = self.dcp.get_pw_material_ez(self.idx, (0,0,0), cmplx=True)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.n)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sampleReal = self.dcp.get_pw_material_hx(self.idx, (0,0,0), cmplx=False)
        sampleCmplx = self.dcp.get_pw_material_hx(self.idx, (0,0,0), cmplx=True)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.n)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sampleReal = self.dcp.get_pw_material_hy(self.idx, (0,0,0), cmplx=False)
        sampleCmplx = self.dcp.get_pw_material_hy(self.idx, (0,0,0), cmplx=True)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.n)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sampleReal = self.dcp.get_pw_material_hz(self.idx, (0,0,0), cmplx=False)
        sampleCmplx = self.dcp.get_pw_material_hz(self.idx, (0,0,0), cmplx=True)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.n)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    