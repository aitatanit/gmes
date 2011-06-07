#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from sys import modules
if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')
import pylab

#from tables import openFile

class Probe(object):
    def __init__(self, filename, pw_material):
        self.pw_mat = pw_material
        self.f = open(filename, 'w')
        
    def __del__(self):
        self.f.close()
        
    def update(self, field1, field2, field3, d1, d2, dt, n):
        self.pw_mat.update(field1, field2, field3, d1, d2, dt, n)
        idx = self.pw_mat.i, self.pw_mat.j, self.pw_mat.k
        self.f.write(str(n) + ' ' + str(field1[idx]) + '\n')
        
        
def write_hdf5(data, name, low_index, high_index):
    h5file = openFile(name + '.h5', mode='w')
    group = h5file.createGroup('/')
    h5file.createArray(group, name, data[low_index[0]:high_index[0], low_index[1]:high_index[1], low_index[2]:high_index[2]])
    h5file.close()

def snapshot(data, filename, title):
    pylab.title(title)
    pylab.imshow(data, origin="lower")
    pylab.savefig(filename)
    
