#!/usr/bin/env python
# file: pmltest2d02.py
# author: Huioon Kim

"""Perform a {C,U}PML test with altering PML parameters."""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from numpy import array, arange
from sys import stdout
from gmes import *
from pmltest2d01 import *

# general settings
save_fname = 'ref_20080825.dat'

# common settings #1
res = 20
def_mat = geometric.DefaultMaterial(material = material.Dielectric())
src = [source.Dipole(src_time = source.Continuous(freq = 1), component = constants.Ez, pos = (0, 0, 0))]

tst_size = (2, 2, 0)
tst_space = geometric.Cartesian(size = tst_size, resolution = res)

pml_thickness = 0.5

probe_ez_idx1_x = (tst_size[0] / 2) - (pml_thickness + 0.1) # X component value of Ez index to probe at upper corner of test space
probe_ez_idx1_y = (tst_size[1] / 2) - (pml_thickness + 0.1) # Y component value of Ez index to probe at upper corner of test sapce

probe_ez_idx2_x = 0 # X component value of Ez index to probe at right edge of test space
probe_ez_idx2_y = (tst_size[1] / 2) - (pml_thickness + 0.1) # Y component value of Ez index to probe at right edge of test space

probe_ez_idx3_x = -probe_ez_idx1_x # X component value of Ez index to probe at lower corner of test space
probe_ez_idx3_y = -probe_ez_idx1_y # Y component value of Ez index to probe at lower corner of test space

tst_probe_ez_idx1 = tst_space.space_to_ez_index( \
        (probe_ez_idx1_x, probe_ez_idx1_y, 0) \
        ) # Ez index to probe in test space (upper corner of test space)
tst_probe_ez_idx2 = tst_space.space_to_ez_index( \
        (probe_ez_idx2_x, probe_ez_idx2_y, 0) \
        ) # Ez index to probe in test space (right edge of test space)
tst_probe_ez_idx3 = tst_space.space_to_ez_index(( \
        probe_ez_idx3_x, probe_ez_idx3_y, 0) \
        ) # Ez index to probe in test space (lower corner of test space)

tst_prob_ez_idxs = [tst_probe_ez_idx1, tst_probe_ez_idx2, tst_probe_ez_idx3]

# generate test fdtd objects
print "---------------------"
print "FDTD generation start"
print "---------------------"
print

tst_fdtd_list = []

for a_max in arange(0, 0.26, 0.01):
    print "-- Generate FDTD with CPML whose alpha_max value is %s." % a_max
    print

    cpml_boundary = geometric.Boundary(material = material.CPML(alpha_max = a_max), thickness = pml_thickness, size = tst_size)
    cpml_tst_geoms = [def_mat, cpml_boundary]
    tst_fdtd_list.append(create_fdtd(tst_space, cpml_tst_geoms, src))
    print

#for k_max in range(0, 45):
#    print "-- Generate FDTD with CPML whose kappa max value is %s." % k_max
#    print
#
#    cpml_boundary = geometric.Boundary(material = material.CPML(kappa_max = k_max), thickness = pml_thickness, size = tst_size)
#    cpml_tst_geoms = [def_mat, cpml_boundary]
#    tst_fdtd_list.append(create_fdtd(tst_space, cpml_tst_geoms, src))
#    print

print "-------------------"
print "FDTD generation end"
print "-------------------"
print

print "---------------------------------"
print "CPML test values acqusition start"
print "---------------------------------"
print

ref_prob_ez_vals = load_vals(save_fname)

tst_prob_ez_vals_list = []

count = 1
for tst_fdtd in tst_fdtd_list:
    print "-- test values acquisition (%s/%s)." % (count, len(tst_fdtd_list))
    print

    tst_prob_ez_vals_list.append(acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals[0])))
    count += 1

print "------------------------------------"
print "CPML test values acqusition complete"
print "------------------------------------"
print

print "Now, the result graph is drawn..."
print

max_error_list = []

for tst_prob_ez_vals in tst_prob_ez_vals_list:
    errors = abs((array(ref_prob_ez_vals[0]) - array(tst_prob_ez_vals[0])) / max(ref_prob_ez_vals[0]))
    max_error_list.append(max(errors))

#print "Max errors:"
#print max_error_list

import pylab

pylab.plot(max_error_list)
#pylab.semilogy(plot_vals)
pylab.show()
