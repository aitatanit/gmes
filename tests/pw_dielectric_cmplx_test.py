#!/usr/bin/env python

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from gmes.pw_material import *
import unittest


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
        sample = DielectricExCmplx(self.idx)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sample = DielectricEyCmplx(self.idx)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sample = DielectricEzCmplx(self.idx)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sample = DielectricHxCmplx(self.idx)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sample = DielectricHyCmplx(self.idx)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sample = DielectricHzCmplx(self.idx)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    