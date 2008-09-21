# file: pmltest2d04.py
# author: Huioon Kim

"""CPML and UPML test script.

Perform a CPML and UPML test relating the kappa max value
with the alpha max value using a low frequency source,
and parallelize using MPI."""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from numpy import *

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
tst_save_fname = 'tst_20080903_lfs_km_am.dat' # 'lfs', 'km', 'am'  mean 'low frequency source', 'kappa max', and 'alpha max', respectively.

# simulation parameters
START_ALPHA = 0.
END_ALPHA = 11.
DELTA_ALPHA = 1.

START_KAPPA = 1.
END_KAPPA = 11.
DELTA_KAPPA = 1.

# common settings #1
res = 20
def_mat = geometry.DefaultMaterial(material = material.Dielectric())

tst_size = (2, 2, 0)
tst_space = geometry.Cartesian(size = tst_size, resolution = res)

pml_thickness = 0.5

# X and Y component value of Ez index to probe at upper corner of test space
probe_ez_idx1_x = (tst_size[0] / 2) - (pml_thickness + 0.1)
probe_ez_idx1_y = (tst_size[1] / 2) - (pml_thickness + 0.1)

# X and Y component value of Ez index to probe at right edge of test space
probe_ez_idx2_x = 0
probe_ez_idx2_y = (tst_size[1] / 2) - (pml_thickness + 0.1)

 # X and Y component value of Ez index to probe at lower corner of test space
probe_ez_idx3_x = -probe_ez_idx1_x
probe_ez_idx3_y = -probe_ez_idx1_y

# Ez index to probe in test space (upper corner of test space)
tst_probe_ez_idx1 = tst_space.space_to_ez_index( \
        (probe_ez_idx1_x, probe_ez_idx1_y, 0) \
        )

# Ez index to probe in test space (right edge of test space)
tst_probe_ez_idx2 = tst_space.space_to_ez_index( \
        (probe_ez_idx2_x, probe_ez_idx2_y, 0) \
        ) 

# Ez index to probe in test space (lower corner of test space)
tst_probe_ez_idx3 = tst_space.space_to_ez_index(( \
        probe_ez_idx3_x, probe_ez_idx3_y, 0) \
        ) 

tst_prob_ez_idxs = (tst_probe_ez_idx1, tst_probe_ez_idx2, tst_probe_ez_idx3)

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
    a_max_list = arange(START_ALPHA + DELTA_ALPHA * myid, END_ALPHA, DELTA_ALPHA * numprocs)
    k_max_list = arange(START_KAPPA, END_KAPPA, DELTA_KAPPA)

    for a_max in a_max_list:
        temp_list = []

        for k_max in k_max_list:
            print 'node %d: [%d / %d]' % (myid, count, len(a_max_list) * len(k_max_list))

            cpml_boundary = geometry.Boundary(material = material.CPML(kappa_max = k_max, alpha_max = a_max), thickness = pml_thickness, size = tst_size)
            cpml_tst_geoms = [def_mat, cpml_boundary]
            tst_fdtd = create_fdtd(tst_space, cpml_tst_geoms, verbose=False)
            temp_list.append(acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals[0]), verbose=False))

            count += 1

        tst_prob_ez_vals_list.append(temp_list)

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
    print "Drawing graph."
    print

    tst_prob_ez_vals_list = load_vals(tst_save_fname)
    #print array(tst_prob_ez_vals_list).shape
    max_error_list = []
    
    for item in tst_prob_ez_vals_list:
        temp_list = []
        
        for tst_prob_ez_vals in item:
            errors = abs((array(ref_prob_ez_vals[0]) - array(tst_prob_ez_vals[0])) / max(ref_prob_ez_vals[0]))
#            print 'DEBUG:', len(errors)
            temp_list.append(max(errors))

        max_error_list.append(temp_list)

    max_error_list = 10 * log10(array(max_error_list).T)

    import pylab

    pylab.title('Maximum relative error')
#    pylab.xlabel(r'$\alpha_\mathrm{max}$')
#    pylab.ylabel(r'$\kappa_\mathrm{max}$')
    pylab.xlabel(r'$\alpha_{max}$')
    pylab.ylabel(r'$\kappa_{max}$')
    pylab.imshow(max_error_list, origin='lower', aspect='auto', \
            extent=(START_ALPHA, END_ALPHA - DELTA_ALPHA, START_KAPPA, END_KAPPA - DELTA_KAPPA))
    pylab.colorbar()
    pylab.show()

