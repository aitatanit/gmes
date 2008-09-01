#!/usr/bin/env python
# file: pmltest2d05.py
# author: Huioon Kim

"""Perform a CPML and UPML test relating the alpha maximum value
with the sigma ratio, i.e, sigma maximum divided by sigma optimum value
using a low frequency source, and parallelize using MPI."""

# import statements
import sys
sys.path.append('../')

from time import localtime
from numpy import add, arange, array, log10

try:
    import mpi
    myid = mpi.rank
    numprocs = mpi.size
except ImportError:
    myid = 0
    numprocs = 1

from gmes import *

from pmltest2d01 import *
from pmltest2d03 import *

# general settings
acquisition = True
ref_save_fname = 'ref_20080830.dat'
tst_save_fname = 'tst_20080901_lfs_sm_am.dat' # 'lfs', 'sm', 'am'  mean 'low frequency source', 'sigma max ratio', and 'alpha max', respectively.

# simulation parameters
A_MAX_START = 0.01
A_MAX_STOP = 0.11
A_MAX_STEP = 0.01

S_RATIO_START = 3.0
S_RATIO_STOP = 4.1
S_RATIO_STEP = 0.1

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

    a_max_range = arange(A_MAX_START + A_MAX_STEP * myid, A_MAX_STOP, A_MAX_STEP * numprocs)
    s_ratio_range = arange(S_RATIO_START, S_RATIO_STOP, S_RATIO_STEP)

    count = 1
    for a_max in a_max_range:
        temp_list = []

        for s_ratio in s_ratio_range:
            print 'node %d: [%d / %d]' % (myid, count, len(a_max_range) * len(s_ratio_range))
            print

            cpml_boundary = geometric.Boundary(material = material.CPML( \
                    alpha_max = a_max, sigma_max_ratio = s_ratio \
                    ), thickness = pml_thickness, size = tst_size)
            cpml_tst_geoms = [def_mat, cpml_boundary]
            tst_fdtd = create_fdtd(tst_space, cpml_tst_geoms)
            temp_list.append(acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals[0])))

            count += 1

        tst_prob_ez_vals_list.append(temp_list)

    save_vals(tst_prob_ez_vals_list, tst_save_fname)

    if numprocs > 0:
        if myid == 0:
            print "Collecting data from nodes."
        size_list = mpi.gather([len(tst_prob_ez_vals_list)])
        collection = mpi.gather(tst_prob_ez_vals_list)
    
        if myid == 0:
            relocated_collection = relocate(size_list, collection)
            save_vals(relocated_collection, tst_save_fname)
    else:
        save_vals(tst_prob_ez_vals_list, tst_save_fname)

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

    import pylab

    pylab.imshow( \
            10 * log10(array(max_error_list)), \
            origin='lower', aspect='auto', \
            extent=(S_RATIO_START, S_RATIO_STOP - S_RATIO_STEP, A_MAX_START, A_MAX_STOP - A_MAX_STEP))
    pylab.title('Maximum relative error')
#    pylab.xlabel(r'$\sigma_\mathrm{max}$' + '/' + r'$\sigma_\mathrm{opt}$')
#    pylab.ylabel(r'$\alpha_\mathrm{max}$')
    pylab.xlabel(r'$\sigma_{max}$' + '/' + r'$\sigma_{opt}$')
    pylab.ylabel(r'$\alpha_{max}$')
    pylab.colorbar()
    pylab.show()

