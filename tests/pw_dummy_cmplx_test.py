#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
from gmes.pw_material import * 

    
class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.cmplxA = array(random.random((3,3,3)), complex)
        self.cmplxB = empty((3,3,3))
        self.cmplxC = empty((3,3,3))
        
        self.copyA = array(self.cmplxA)
        
        self.diff = 1
        self.t = 0
        
    def testEx(self):
        sampleCmplx = DummyExCmplx(self.idx)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.cmplxA == self.copyA).all(), True)
                         
    def testEy(self):
        sampleCmplx = DummyEyCmplx(self.idx)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.cmplxA == self.copyA).all(), True)
        
    def testEz(self):
        sampleCmplx = DummyEzCmplx(self.idx)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.cmplxA == self.copyA).all(), True)
        
    def testHx(self):
        sampleCmplx = DummyHxCmplx(self.idx)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.cmplxA == self.copyA).all(), True)
        
    def testHy(self):
        sampleCmplx = DummyHyCmplx(self.idx)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.cmplxA == self.copyA).all(), True)
        
    def testHz(self):
        sampleCmplx = DummyHzCmplx(self.idx)
        sampleCmplx.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.cmplxA == self.copyA).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    