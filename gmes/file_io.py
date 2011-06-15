#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')
    
from sys import modules
if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')
import pylab

#from tables import openFile

# GMES modules
from pw_material import MaterialElectricReal, MaterialElectricCmplx
from pw_material import MaterialMagneticReal, MaterialMagneticCmplx

class Probe(object):
    def __init__(self, filename, pw_material):
        self.pw_mat = pw_material

        if isinstance(self.pw_mat, 
                      (MaterialElectricReal, MaterialElectricCmplx)):
            self.epsilon = self.pw_mat.epsilon
            
        if isinstance(self.pw_mat, 
                      (MaterialMagneticReal, MaterialMagneticCmplx)):
            self.mu = self.pw_mat.mu

        self.f = open(filename, 'w')
        
    def __del__(self):
        self.f.close()
        
    def update(self, field1, field2, field3, d1, d2, dt, n, i, j, k):
        idx = i, j, k
        self.pw_mat.update(field1, field2, field3, d1, d2, dt, n, *idx)
        self.f.write(str(n) + ' ' + str(field1[idx]) + '\n')
        
        
def write_hdf5(data, name, low_index, high_index):
    h5file = openFile(name + '.h5', mode='w')
    group = h5file.createGroup('/')
    h5file.createArray(group, name, 
                       data[low_index[0]:high_index[0], 
                            low_index[1]:high_index[1], 
                            low_index[2]:high_index[2]])
    h5file.close()

def snapshot(data, filename, title):
    pylab.title(title)
    pylab.imshow(data, origin="lower")
    pylab.savefig(filename)
    
