#!/usr/bin/env python
# file: pmltest2d04.py
# author: Huioon Kim

"""Perform a CPML and UPML test relating the alpha maximum value
with the sigma ratio, i.e, sigma maximum divided by sigma optimum value
using a low frequency source."""

import sys
sys.path.append('../')

from gmes import *
from numpy import array, arange, log10
from sys import stdout
from pmltest2d01 import *
from pmltest2d03 import *

# general settings
acquisition = True
ref_save_fname = 'ref_20080830.dat'
tst_save_fname = 'tst_20080831_lfs_sm_am.dat' # 'lfs', 'sm', 'am'  mean 'low frequency source', 'sigma max', and 'alpha max', respectively.

# common settings #1
res = 20
def_mat = geometric.DefaultMaterial(material = material.Dielectric())

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

ref_prob_ez_vals = load_vals(ref_save_fname)

a_max_range = arange(0, 0.21, 0.01)
s_ratio_range = arange(2.0, 4.1, 0.1)

if acquisition == True:
    tst_prob_ez_vals_list = []

    count = 1
    for a_max in a_max_range:
        temp_list = []

        for s_ratio in s_ratio_range:
            counter_str = "[%s / %s]" % (count, len(a_max_range) * len(s_ratio_range))

            print "-" * len(counter_str)
            print counter_str
            print "-" * len(counter_str)
            print

            print "-- Generate FDTD with CPML whose alpha max value is %s and sigma max ratio is %s." % (a_max, s_ratio)
            print

            cpml_boundary = geometric.Boundary(material = material.CPML( \
                    alpha_max = a_max, sigma_max_ratio = s_ratio \
                    ), thickness = pml_thickness, size = tst_size)
            cpml_tst_geoms = [def_mat, cpml_boundary]
            tst_fdtd = create_fdtd(tst_space, cpml_tst_geoms)
            print

            print "-- Generation complete."
            print
            
            print "-- Acquire test values."
            print
            
            temp_list.append(acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals[0])))

            print "-- Acquisition complete."
            print

            count += 1

        tst_prob_ez_vals_list.append(temp_list)

    save_vals(tst_prob_ez_vals_list, tst_save_fname)

print "Now, the result graph is drawn..."
print

tst_prob_ez_vals_list = load_vals(tst_save_fname)

max_error_list = []

for item in tst_prob_ez_vals_list:
    temp_list = []

    for tst_prob_ez_vals in item:
        errors = abs((array(ref_prob_ez_vals[0]) - array(tst_prob_ez_vals[0])) / max(ref_prob_ez_vals[0]))
        temp_list.append(max(errors))

    max_error_list.append(temp_list)

import pylab

#pylab.contour(max_error_list)
pylab.imshow(10*log10(array(max_error_list)), origin='lower')
pylab.title('Maximum relative error')
pylab.xlabel(r'$\sigma_\mathrm{max}$' + '/' + r'$\sigma_\mathrm{opt}$')
pylab.ylabel(r'$\alpha_\mathrm{max}$')
pylab.colorbar()
pylab.show()

