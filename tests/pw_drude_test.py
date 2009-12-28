#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
from numpy import *
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
        self.epsilon_inf = 1
        self.omega_p = self.gamma_p = (1,)
        
    def testEx(self):
        sampleReal = DrudeExReal(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sampleCmplx = DrudeExCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sampleReal = DrudeEyReal(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sampleCmplx = DrudeEyCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sampleReal = DrudeEzReal(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sampleCmplx = DrudeEzCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sampleReal = DrudeHxReal(self.idx, self.epsilon_inf)
        sampleCmplx = DrudeHxCmplx(self.idx, self.epsilon_inf)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sampleReal = DrudeHyReal(self.idx, self.epsilon_inf)
        sampleCmplx = DrudeHyCmplx(self.idx, self.epsilon_inf)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sampleReal = DrudeHzReal(self.idx, self.epsilon_inf)
        sampleCmplx = DrudeHzCmplx(self.idx, self.epsilon_inf)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff, self.t)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    