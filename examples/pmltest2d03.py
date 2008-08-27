#!/usr/bin/env python
# file: pmltest2d03.py
# author: Huioon Kim

"""Perfrom a CPML and UPML test with a low frequency source.

- Class
LowFreqSrc: Contain the update mathod which implement a low frequency source properties.

- Functions
create_fdtd: Create and return a new FDTD object which have low frequency source.
"""

import sys
sys.path.append('../')

from gmes import *
from numpy import array, exp
from sys import stdout
from pmltest2d01 import *

class LowFreqSrc:
    """Generate low frequency source.

    __init__: constructor
    update: update Ez value according to time sequence.
    """

    def __init__(self, idx):
        """Assign initial values."""

        self.idx = idx
        self.tw = 0.00795349
        self.t0 = 4 * self.tw
        self.t = 0.0
        self.epsilon = 1

    def update(self, ez, hy, hx, dt, dx, dy):
        """Update Ez values according to time sequence."""

        i, j, k = self.idx
        j_src = -2 * (self.t - self.t0) / self.tw * exp(-((self.t - self.t0) / self.tw) ** 2)
        ez[self.idx] += dt / self.epsilon * ((hy[i+1,j,k+1] - hy[i,j,k+1]) / dx - (hx[i,j+1,k+1] - hx[i,j,k+1]) / dy - j_src)
        self.t += dt

def create_fdtd(space, geoms):
    """Create and return a new FDTD object which have low frequency source using passing parameters."""

    src_idx = space.space_to_ez_index((0, 0, 0))
    my_fdtd = fdtd.TMzFDTD(space, geoms, [])
    my_fdtd.material_ez[src_idx] = LowFreqSrc(src_idx)

    return my_fdtd

if __name__ == "__main__":
    # general settings
    debug_mode = False
    ref_acq = True
    save_fname = 'ref_20080826.dat'

    # common values #1
    res = 20
    def_mat = geometric.DefaultMaterial(material = material.Dielectric())

    # reference settings
    if ref_acq == True:
        ref_size = (50, 50, 0)
        ref_space = geometric.Cartesian(size = ref_size, resolution = res, courant_ratio = 0.99)
        ref_geoms = [def_mat]

    # test settings
    tst_size = (2, 2, 0)
    tst_space = geometric.Cartesian(size = tst_size, resolution = res, courant_ratio = 0.99)

    pml_thickness = 0.5
    cpml_boundary = geometric.Boundary(material = material.CPML(), thickness = pml_thickness, size = tst_size)
    upml_boundary = geometric.Boundary(material = material.UPML(), thickness = pml_thickness, size = tst_size)

    cpml_tst_geoms = [def_mat, cpml_boundary]
    upml_tst_geoms = [def_mat, upml_boundary]

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

        ref_fdtd = create_fdtd(ref_space, ref_geoms)
        print

        ref_prob_ez_vals = acquire_ez_vals(ref_fdtd, ref_prob_ez_idxs, AcqMode.REFERENCE)

        save_vals(ref_prob_ez_vals, save_fname)

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

    tst_fdtd = create_fdtd(tst_space, cpml_tst_geoms)
    print

    #########################TEMPORARY BLOCK#########################
#    import cPickle
#
#    ref_ez_vals1_file = open('ref_ez_vals1.dat')
#    ref_ez_vals2_file = open('ref_ez_vals2.dat')
#    ref_ez_vals3_file = open('ref_ez_vals3.dat')
#
#    ref_prob_ez_vals1 = cPickle.load(ref_ez_vals1_file)
#    ref_prob_ez_vals2 = cPickle.load(ref_ez_vals2_file)
#    ref_prob_ez_vals3 = cPickle.load(ref_ez_vals3_file)
#
#    ref_ez_vals1_file.close()
#    ref_ez_vals2_file.close()
#    ref_ez_vals3_file.close()
#
#    tst_prob_ez_vals = acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals1))
    #########################TEMPORARY BLOCK#########################

    if ref_acq != True:
        ref_prob_ez_vals = load_vals(save_fname)

    tst_prob_ez_vals = acquire_ez_vals(tst_fdtd, tst_prob_ez_idxs, AcqMode.TEST, len(ref_prob_ez_vals[0]))

    print "------------------------------------"
    print "CPML test values acqusition complete"
    print "------------------------------------"
    print

    #######DEBUG MODE MESSAGE BLOCK#######
    if debug_mode == True:
        print "!!!!!!!!!!!START OF DEBUG MESSAGES!!!!!!!!!!!"
        print "Length of ref_prob_ez_vals:", len(ref_prob_ez_vals)
        print
        print "Length of tst_prob_ez_vals:", len(tst_prob_ez_vals)
        print
        print "Length of ref_prob_ez_vals[0]:", len(ref_prob_ez_vals[0])
        print
        print "Length of tst_prob_ez_vals[0]:", len(tst_prob_ez_vals[0])
        print
        print "ref_prob_ez_vals:", array(ref_prob_ez_vals)
        print
        print "tst_prob_ez_vals:", array(tst_prob_ez_vals)
        print "!!!!!!!!!!!END OF DEBUG MESSAGES!!!!!!!!!!!"
        print
    #######DEBUG MODE MESSAGE BLOCK#######

    print "Now, the result graph is drawn..."

    plot_vals = abs((array(ref_prob_ez_vals[1]) - array(tst_prob_ez_vals[1])) / max(ref_prob_ez_vals[1]))

    import pylab

    pylab.plot(plot_vals)
    pylab.title("right edge point")
    pylab.semilogy(plot_vals)
    pylab.show()

