# file: pmltest2d01.py
# author: Huioon Kim

"""Perform a CPML and UPML test.

Class
    AcqMode -- Contain data acquisition modes

Functions
    create_fdtd -- Create and return a new FDTD object.
    acquire_ez_vals -- Acquire Ez values from updating a FDTD object.
    save_vals -- Write data into a file
    load_vals -- Read data from a file

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

from numpy import array
from sys import stdout
from gmes import *

class AcqMode:
    """Contain data acquisition modes, reference or test."""

    REFERENCE = 0
    TEST = 1

def create_fdtd(space, geoms, src, verbose = True):
    """Create and return a new FDTD object using passing parameters."""

    return fdtd.TMzFDTD(space, geoms, src, verbose)

def acquire_ez_vals(fdtd, prob_ez_idxs, mode, stepnum_for_tst = None, verbose = True):
    """Acquire Ez values from updating passing FDTD object. Do differently according to the mode."""

    #fdtd.show_ez(constants.Z(), 0)

    if verbose == True:
        try:
            print "-- The Ez index of the position of source: %s" % str(fdtd.space.space_to_ez_index(fdtd.src_list[0].pos))
            print
        except IndexError:
            pass

        print "-- The space size:", tuple([int(item) for item in fdtd.space.half_size * 2])
        print

    if mode == AcqMode.REFERENCE:
        cent_rmost_ez_idx = fdtd.space.space_to_ez_index( \
                (0, fdtd.space.half_size[1] - 0.1, 0) \
                ) # The Ez index 0.1 space unit away from center-rightmost of the space

        if verbose == True:
            print "-- The Ez index 0.1 space unit away from center-rightmost of the space:", cent_rmost_ez_idx
            print

    result_ez_vals = []

    for idx in prob_ez_idxs:
        if verbose == True:
            print "-- The Ez index to probe:", idx
        result_ez_vals.append([])
    if verbose == True:
        print

    if verbose == True:
        print "-- FDTD update start..."
        print

        print "Time step:"

    if mode == AcqMode.REFERENCE:
        while 1:
            fdtd.step()

            for i, idx in zip(range(len(prob_ez_idxs)), prob_ez_idxs):
                result_ez_vals[i].append(fdtd.ez[idx])

            if verbose == True:
                print "\r[%s]" % int(fdtd.time_step.n),
                stdout.flush()

            if fdtd.ez[cent_rmost_ez_idx] != 0.0:
                break

        if verbose == True:
            print "\n"
    elif mode == AcqMode.TEST:
        for x in range(stepnum_for_tst):
            fdtd.step()

            for i, idx in zip(range(len(prob_ez_idxs)), prob_ez_idxs):
                result_ez_vals[i].append(fdtd.ez[idx])

            if verbose == True:
                print "\r[%s]" % int(fdtd.time_step.n),
                stdout.flush()

        if verbose == True:
            print "\n"

    if verbose == True:
        print "-- FDTD update end..."
        print

    return result_ez_vals

def save_vals(vals, fname, verbose = True):
    """Write passing values into a file whose name is passing through the second parameter using pickling."""

    import cPickle

    vals_file = open(fname, 'w')
    cPickle.dump(vals, vals_file)
    vals_file.close()

    if verbose == True:
        print "The values are successfully saved in %s file." % fname
        print

def load_vals(fname, verbose = True):
    """Read values from a file whose name is passing through the parameter using unpickling."""

    import cPickle

    vals_file = open(fname)
    vals = cPickle.load(vals_file)
    vals_file.close()
    
    if verbose == True:
        print "Values are successfully loaded from %s file." % fname
        print

    return vals

if __name__ == "__main__":
    # general settings
    debug_mode = False
    ref_acq = False
    save_fname = 'ref_20080905.dat'

    # common values #1
    res = 20
    def_mat = geometric.DefaultMaterial(material = material.Dielectric())
    src = [source.Dipole(src_time = source.Continuous(freq = 1), component = constants.Ez, pos = (0, 0, 0))]

    # reference settings
    if ref_acq == True:
        ref_size = (50, 50, 0)
        ref_space = geometric.Cartesian(size = ref_size, resolution = res)
        ref_geoms = [def_mat]

    # test settings
    tst_size = (2, 2, 0)
    tst_space = geometric.Cartesian(size = tst_size, resolution = res)

    pml_thickness = 0.5
    cpml_boundary = geometric.Boundary(material = material.CPML(), thickness = pml_thickness, size = tst_size)
    #upml_boundary = geometric.Boundary(material = material.UPML(), thickness = pml_thickness, size = tst_size)

    cpml_tst_geoms = [def_mat, cpml_boundary]
    #upml_tst_geoms = [def_mat, upml_boundary]

    # common values #2
    probe_ez_idx1_x = (tst_size[0] / 2) - (pml_thickness + 0.1) # X component value of Ez index to probe at upper corner of test space
    probe_ez_idx1_y = (tst_size[1] / 2) - (pml_thickness + 0.1) # Y component value of Ez index to probe at upper corner of test sapce

    probe_ez_idx2_x = 0 # X component value of Ez index to probe at right edge of test space
    probe_ez_idx2_y = (tst_size[1] / 2) - (pml_thickness + 0.1) # Y component value of Ez index to probe at right edge of test space

    probe_ez_idx3_x = -probe_ez_idx1_x # X component value of Ez index to probe at lower corner of test space
    probe_ez_idx3_y = -probe_ez_idx1_y # Y component value of Ez index to probe at lower corner of test space

    # reference values acquisition
    if ref_acq == True:
        print "---------------------------------"
        print "Reference values acqusition start"
        print "---------------------------------"
        print

        ref_probe_ez_idx1 = ref_space.space_to_ez_index( \
                (probe_ez_idx1_x, probe_ez_idx1_y, 0) \
                ) # Ez index to probe in reference space (upper corner of test space)
        ref_probe_ez_idx2 = ref_space.space_to_ez_index( \
                (probe_ez_idx2_x, probe_ez_idx2_y, 0) \
                ) # Ez index to probe in reference space (right edge of test space)
        ref_probe_ez_idx3 = ref_space.space_to_ez_index(( \
                probe_ez_idx3_x, probe_ez_idx3_y, 0) \
                ) # Ez index to probe in reference space (lower corner of test space)

        ref_prob_ez_idxs = [ref_probe_ez_idx1, ref_probe_ez_idx2, ref_probe_ez_idx3]
        
        ref_fdtd = create_fdtd(ref_space, ref_geoms, src)
        print

        ref_ez_vals = acquire_ez_vals(ref_fdtd, ref_prob_ez_idxs, AcqMode.REFERENCE)

        save_vals(ref_ez_vals, save_fname)

        print "------------------------------------"
        print "Reference values acqusition complete"
        print "------------------------------------"
        print

    # test values acquistion
    print "---------------------------------"
    print "CPML test values acqusition start"
    print "---------------------------------"
    print

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

    tst_fdtd = create_fdtd(tst_space, cpml_tst_geoms, src)
    print
    
    tst_ez_vals = acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_ez_vals[0]))
    
    print "------------------------------------"
    print "CPML test values acqusition complete"
    print "------------------------------------"
    print

    if ref_acq != True:
        ref_ez_vals = load_vals(save_fname)

    #######DEBUG MODE MESSAGE BLOCK#######
    if debug_mode == True:
        print "!!!!!!!!!!!START OF DEBUG MESSAGES!!!!!!!!!!!"
        print "Length of ref_prob_ez_vals:", len(ref_ez_vals)
        print
        print "Length of tst_prob_ez_vals:", len(tst_ez_vals)
        print
        print "Length of ref_prob_ez_vals[0]:", len(ref_ez_vals[0])
        print
        print "Length of tst_prob_ez_vals[0]:", len(tst_ez_vals[0])
        print
        print "ref_prob_ez_vals:", array(ref_ez_vals)
        print
        print "tst_prob_ez_vals:", array(tst_ez_vals)
        print "!!!!!!!!!!!END OF DEBUG MESSAGES!!!!!!!!!!!"
        print
    #######DEBUG MODE MESSAGE BLOCK#######

    print "Now, the result graph is drawn..."

    plot_vals = abs((array(ref_ez_vals[1]) - array(tst_ez_vals[1])) / max(ref_ez_vals[1]))

    import pylab

    pylab.plot(plot_vals)
    pylab.title("right edge point")
    pylab.semilogy(plot_vals)
    pylab.show()

