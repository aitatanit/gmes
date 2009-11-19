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
        
    def testEx(self):
        sampleReal = DielectricExReal(self.idx)
        sampleCmplx = DielectricExCmplx(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sampleReal = DielectricEyReal(self.idx)
        sampleCmplx = DielectricEyCmplx(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sampleReal = DielectricEzReal(self.idx)
        sampleCmplx = DielectricEzCmplx(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sampleReal = DielectricHxReal(self.idx)
        sampleCmplx = DielectricHxCmplx(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sampleReal = DielectricHyReal(self.idx)
        sampleCmplx = DielectricHyCmplx(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sampleReal = DielectricHzReal(self.idx)
        sampleCmplx = DielectricHzCmplx(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                           self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    