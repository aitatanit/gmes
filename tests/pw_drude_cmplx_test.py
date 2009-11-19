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
        self.epsilon_inf = 1
        self.omega_p = self.gamma_p = (1,)
        
    def testEx(self):
        sample1 = DrudeExCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sample1.update(self.realA, self.realB, self.realC, 
                       self.diff, self.diff, self.diff)
        sample2 = DrudeExCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sample2.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                       self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEy(self):
        sample1 = DrudeEyCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sample1.update(self.realA, self.realB, self.realC, 
                       self.diff, self.diff, self.diff)
        sample2 = DrudeEyCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sample2.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                       self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testEz(self):
        sample1 = DrudeEzCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sample1.update(self.realA, self.realB, self.realC, 
                       self.diff, self.diff, self.diff)
        sample2 = DrudeEzCmplx(self.idx, self.epsilon_inf, self.omega_p, self.gamma_p)
        sample2.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                       self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHx(self):
        sample = DrudeHxCmplx(self.idx, self.epsilon_inf)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHy(self):
        sample = DrudeHyCmplx(self.idx, self.epsilon_inf)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
    def testHz(self):
        sample = DrudeHzCmplx(self.idx, self.epsilon_inf)
        sample.update(self.realA, self.realB, self.realC, 
                      self.diff, self.diff, self.diff)
        sample.update(self.cmplxA, self.cmplxB, self.cmplxC, 
                      self.diff, self.diff, self.diff)
        self.assertEqual((self.realA == self.cmplxA.real).all(), True)
        
        
if __name__ == '__main__':
    unittest.main(argv=('', '-v'))
    