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
        
        self.epsilon = random.random()
        self.value = random.random()
        self.diff = 1
        self.t = 0
        
        self.copyA = array(self.realA)
        self.copyA[self.idx] = self.value
        
    def testAttributes(self):
        ex = ConstExReal(self.idx, self.epsilon, self.value)
        self.assertEqual(self.idx == (ex.i, ex.j, ex.k), True)
        self.assertEqual(self.epsilon == ex.epsilon, True)
        self.assertEqual(self.value == ex.value, True)
        
        ey = ConstEyReal(self.idx, self.epsilon, self.value)
        self.assertEqual(self.idx == (ey.i, ey.j, ey.k), True)
        self.assertEqual(self.epsilon == ey.epsilon, True)
        self.assertEqual(self.value == ey.value, True)
        
        ez = ConstEyReal(self.idx, self.epsilon, self.value)
        self.assertEqual(self.idx == (ez.i, ez.j, ez.k), True)
        self.assertEqual(self.epsilon == ez.epsilon, True)
        self.assertEqual(self.value == ez.value, True)
        
        hx = ConstExReal(self.idx, self.epsilon, self.value)
        self.assertEqual(self.idx == (hx.i, hx.j, hx.k), True)
        self.assertEqual(self.epsilon == hx.epsilon, True)
        self.assertEqual(self.value == hx.value, True)
        
        hy = ConstEyReal(self.idx, self.epsilon, self.value)
        self.assertEqual(self.idx == (hy.i, hy.j, hy.k), True)
        self.assertEqual(self.epsilon == hy.epsilon, True)
        self.assertEqual(self.value == hy.value, True)
        
        hz = ConstEyReal(self.idx, self.epsilon, self.value)
        self.assertEqual(self.idx == (hz.i, hz.j, hz.k), True)
        self.assertEqual(self.epsilon == hz.epsilon, True)
        self.assertEqual(self.value == hz.value, True)
        
    def testEx(self):
        ConstReal = ConstExReal(self.idx, self.epsilon, self.value)
        ConstReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testEy(self):
        ConstReal = ConstEyReal(self.idx, self.epsilon, self.value)
        ConstReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testEz(self):
        ConstReal = ConstEzReal(self.idx, self.epsilon, self.value)
        ConstReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHx(self):
        ConstReal = ConstHxReal(self.idx, self.epsilon, self.value)
        ConstReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHy(self):
        ConstReal = ConstHyReal(self.idx, self.epsilon, self.value)
        ConstReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
    def testHz(self):
        ConstReal = ConstHzReal(self.idx, self.epsilon, self.value)
        ConstReal.update(self.realA, self.realB, self.realC, 
                          self.diff, self.diff, self.diff, self.t)
        
        self.assertEqual((self.realA == self.copyA).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    