#!/usr/bin/env python
# file: pmltest2d04_mpi.py
# author: Huioon Kim

"""Perform a CPML and UPML test relating the kappa max value
with the alpha max value using a low frequency source."""

import sys
sys.path.append('../')

from numpy import *
from sys import stdout

import mpi

from gmes import *

from pmltest2d01 import *
from pmltest2d03 import *

# mpi parameters
myid = mpi.rank
numprocs = mpi.size

# general settings
acquisition = False
ref_save_fname = 'ref_20080826.dat'
tst_save_fname = 'tst_20080827.dat'

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

def relocate(size_list, collection):
    accum_size = add.accumulate(size_list)
    relocated = []
    
    for i in range(size_list[0]):
        for j in range(len(size_list)):
            if i == size_list[j]:
                break

            if j == 0:
                relocated.append(collection[i])
            else:
                relocated.append(collection[i + accum_size[j-1]])
                    
    return relocated

if acquisition == True:
    tst_prob_ez_vals_list = []

    count = 1
    a_max_list = arange(myid, 3, 1 * numprocs)
    k_max_list = arange(1, 3, 1)

    for a_max in a_max_list:
        temp_list = []

        for k_max in k_max_list:
            print 'node %d: [%d / %d]' % (myid, count, len(a_max_list) * len(k_max_list))

            cpml_boundary = geometric.Boundary(material = material.CPML(kappa_max = k_max, alpha_max = a_max), thickness = pml_thickness, size = tst_size)
            cpml_tst_geoms = [def_mat, cpml_boundary]
            tst_fdtd = create_fdtd(tst_space, cpml_tst_geoms)
            temp_list.append(acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals[0])))

            count += 1

        tst_prob_ez_vals_list.append(temp_list)

    size_list = mpi.gather([len(tst_prob_ez_vals_list)])
    collection = mpi.gather(tst_prob_ez_vals_list)
    
    if myid == 0:        
        relocated_collection = relocate(size_list, collection)
        save_vals(relocated_collection, tst_save_fname)

if myid == 0:
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

    max_error_list = 10 * log10(array(max_error_list).T)

    import pylab

    #pylab.plot(max_error_list)
    #pylab.contour(max_error_list)
    pylab.xlabel(r'$\alpha_\mathrm{max}$')
    pylab.ylabel(r'$\kappa_\mathrm{max}$')
    pylab.imshow(max_error_list)
    pylab.colorbar()
    pylab.show()

