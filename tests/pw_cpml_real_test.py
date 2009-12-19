#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
from gmes.pw_material import * 

    
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
        self.t = 0
        self.epsilon_r = 1
        self.mu_r = 1
        self.b1 = 1
        self.b2 = 1
        self.c1 = 1
        self.c2 = 1
        self.kappa1 = 1
        self.kappa2 = 1
        
    def testEx(self):
        sampleReal = CpmlExReal(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleCmplx = CpmlExCmplx(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sampleReal = CpmlEyReal(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleCmplx = CpmlEyCmplx(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sampleReal = CpmlEzReal(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleCmplx = CpmlEzCmplx(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sampleReal = CpmlHxReal(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleCmplx = CpmlHxCmplx(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sampleReal = CpmlHyReal(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleCmplx = CpmlHyCmplx(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sampleReal = CpmlHzReal(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleCmplx = CpmlHzCmplx(self.idx, self.epsilon_r, self.b1, self.b2, 
                             self.c1, self.c2, self.kappa1, self.kappa2)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    