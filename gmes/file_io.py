#!/usr/bin/env python

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

from tables import openFile
import pylab


def write_hdf5(data, name, low_index, high_index):
    h5file = openFile(name + '.h5', mode='w')
    group = h5file.createGroup('/')
    h5file.createArray(group, name, data[low_index[0]:high_index[0], low_index[1]:high_index[1], low_index[2]:high_index[2]])
    h5file.close()

def snapshot(data, filename, title):
    pylab.title(title)
    pylab.imshow(data, origin="lower")
    pylab.savefig(filename)