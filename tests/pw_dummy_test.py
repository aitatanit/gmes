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
        
        self.realA = random.random((3,3,3))
        self.realB = empty((3,3,3))
        self.realC = empty((3,3,3))
        
        self.copyA = array(self.realA)
        
        self.diff = 1
        self.n = 0
        
    def testEx(self):
        sampleReal = DummyExReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
                         
    def testEy(self):
        sampleReal = DummyEyReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testEz(self):
        sampleReal = DummyEzReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHx(self):
        sampleReal = DummyHxReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHy(self):
        sampleReal = DummyHyReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHz(self):
        sampleReal = DummyHzReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.n)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    