#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr
from os.path import exists

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass
    
from sys import modules
if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')
import pylab

# from tables import openFile

# GMES modules
from pw_material import MaterialElectricReal, MaterialElectricCmplx
from pw_material import MaterialMagneticReal, MaterialMagneticCmplx


class Probe(object):
    def __init__(self, idx, field, filename):
        """
        idx: index of probing point. type: tuple-3
        field: field to probe. type: numpy.array
        filename: recording file name. type: str

        """
        self.idx = tuple(idx)

        self.field = field

        f_name = str(filename)
        if exists(f_name):
            stderr.write('Warning: ' + f_name + ' already exists.\n')
        try:
            self.f = open(f_name, 'w')
        except IOError:
            self.f = None
            print('Warning: Can\'t open file ' + f_name + '.\n')

    def __del__(self):
        self.f.close()
    
    def write_header(self, p, dt):
        """Write some meta-data on the header of the recording file.
        
        p: space coordinates. type: tuple-3
        dt: time-step. type: float
        
        """
        self.f.write('# location=' + str(p) + '\n')
        self.f.write('# dt=' + str(dt) + '\n')

    def write(self, n):
        self.f.write(str(n) + ' ' + str(self.field[self.idx]) + '\n')
       
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
    
