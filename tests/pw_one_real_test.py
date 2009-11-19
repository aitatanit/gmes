#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import unittest
from gmes.pw_material import * 

    
class TestSequence(unittest.TestCase):
    def setUp(self):
        self.idx = (1,1,1)
        
        self.realA = random.random((3,3,3))
        self.realB = empty((3,3,3))
        self.realC = empty((3,3,3))
        
        self.copyA = array(self.realA)
        self.copyA[self.idx] = 1
        
        self.diff = 1
        
    def testEx(self):
        sampleReal = OneExReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
                         
    def testEy(self):
        sampleReal = OneEyReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testEz(self):
        sampleReal = OneEzReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHx(self):
        sampleReal = OneHxReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHy(self):
        sampleReal = OneHyReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHz(self):
        sampleReal = OneHzReal(self.idx)
        sampleReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    