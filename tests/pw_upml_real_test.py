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
        self.epsilon = self.mu = 1
        self.c = ones(6)
        
    def testEx(self):
        sampleReal = UpmlExReal(self.idx, self.epsilon, *self.c)
        sampleCmplx = UpmlExCmplx(self.idx, self.epsilon, *self.c)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sampleReal = UpmlEyReal(self.idx, self.epsilon, *self.c)
        sampleCmplx = UpmlEyCmplx(self.idx, self.epsilon, *self.c)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sampleReal = UpmlEzReal(self.idx, self.epsilon, *self.c)
        sampleCmplx = UpmlEzCmplx(self.idx, self.epsilon, *self.c)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sampleReal = UpmlHxReal(self.idx, self.mu, *self.c)
        sampleCmplx = UpmlHxCmplx(self.idx, self.mu, *self.c)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sampleReal = UpmlHyReal(self.idx, self.mu, *self.c)
        sampleCmplx = UpmlHyCmplx(self.idx, self.mu, *self.c)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sampleReal = UpmlHzReal(self.idx, self.mu, *self.c)
        sampleCmplx = UpmlHzCmplx(self.idx, self.mu, *self.c)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    