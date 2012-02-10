#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script can handle parallel execution and data-collection. The fields 
in the waveguide are recorded following probe lines. Each field data is 
pickled and saved in separate files.

"""

import os, sys
new_path = os.path.abspath('../')
sys.path.append(new_path)

import cPickle, time
from numpy import *
from gmes import *

# Define simulation parameters.

FREQ = 0.4015
FWIDTH = 0.0455
RADIUS = 0.2
PML_THICK = .5
SIZE = (15, 13, 0)
RESOLUTION = 25
TIME_STEP_T = 800

START_PNT_SPC_1 = (-5, -1, 0)
END_PNT_SPC_1 = (-5, 1, 0)

START_PNT_SPC_2 = (6, -1, 0)
END_PNT_SPC_2 = (6, 1, 0)

FILENAME_EZ_SUMS_1 = 'ez_sums1_t' + `TIME_STEP_T` + '.dat'
FILENAME_HX_SUMS_1 = 'hx_sums1_t' + `TIME_STEP_T` + '.dat'
FILENAME_HY_SUMS_1 = 'hy_sums1_t' + `TIME_STEP_T` + '.dat'

FILENAME_EZ_SUMS_2 = 'ez_sums2_t' + `TIME_STEP_T` + '.dat'
FILENAME_HX_SUMS_2 = 'hx_sums2_t' + `TIME_STEP_T` + '.dat'
FILENAME_HY_SUMS_2 = 'hy_sums2_t' + `TIME_STEP_T` + '.dat'

AIR = material.Dielectric(1)
Al2O3 = material.Dielectric(8.9)


def make_rod(center):
    """Stand a rod."""

    return geometry.Cylinder(material=Al2O3, axis=(0,0,1), radius=RADIUS,
            height=10, center=center)


def make_crystals(x_size, y_size):
    """Stand rods on the lattice positions."""

    a1 = array((1,0), float)
    a2 = array((0,1), float)
    crystals = []

    for i in arange(-ceil(x_size), ceil(x_size)):
        for j in arange(-ceil(y_size), ceil(y_size)):
            center = tuple(i * a1 + j * a2) + (0,)
            if fabs(center[0]) <= .5 * x_size + RADIUS \
                    and fabs(center[1]) <= .5 * y_size + RADIUS:
                        crystals.append(make_rod(center))
    return tuple(crystals)


def pullout_rod(center):
    """Remove a dielectric rod."""

    return geometry.Cylinder(material=AIR, axis=(0,0,1), radius=RADIUS,
            height=10, center=center)


def make_line_defect(length):
    """Remove dielectric rods to form a line defect."""

    line_defect = []
    for i in arange(-ceil(.5 * length + 1), ceil(.5 * length + 2)):
        line_defect.append(pullout_rod((i, 0, 0)))

    return tuple(line_defect)


def make_highpass_bragg_grating(period):
    l = .5 / .28
    bragg = []

    for i in arange(-(.5 * period - .5) * l, (.5 * period + .5) * l, l):
        rod = geometry.Cylinder(material=Al2O3, axis=(0,0,1), radius=0.0623302*l,
                height=10, center=(i, 0, 0))
        bragg.append(rod)

    return tuple(bragg)


def make_lowpass_bragg_grating(period):
    l = .5 / .48
    bragg = []

    for i in arange(-(.5 * period - .5) * l, (.5 * period + .5) * l, l):
        rod = geometry.Cylinder(material=Al2O3, axis=(0,0,1), radius=0.0362054*l,
                                height=10, center=(i, 0, 0))
        bragg.append(rod)

    return tuple(bragg)


if __name__ == "__main__":
    start = time.time()
    geom_list = (geometry.DefaultMedium(material=AIR),) + \
            make_crystals(*SIZE[:2]) + make_line_defect(SIZE[0]) + \
            (geometry.Boundary(material=material.UPML(), thickness=PML_THICK, size=SIZE),)

    space = geometry.Cartesian(size=SIZE, resolution=RESOLUTION)

    ez1_id = 0 % space.numprocs
    hx1_id = 1 % space.numprocs
    hy1_id = 2 % space.numprocs
    ez2_id = 3 % space.numprocs
    hx2_id = 4 % space.numprocs
    hy2_id = 5 % space.numprocs

    src_list = (source.PointSource(src_time=source.Bandpass(freq=FREQ, fwidth=FWIDTH),
                                   component=constant.Ez, 
                                   pos=(-.5 * SIZE[0] + 1, 0, 0)),)

    my_fdtd = fdtd.TMzFDTD(space, geom_list, src_list)
    my_fdtd.init()
    my_fdtd.show_permittivity_ez(constant.Z, 0)
    my_fdtd.show_ez(constant.Z, 0)

    ez_sums1 = []
    hx_sums1 = []
    hy_sums1 = []

    start_ez_idx1 = space.space_to_ez_index(*START_PNT_SPC_1)
    start_hx_idx1 = space.space_to_hx_index(*START_PNT_SPC_1)
    start_hy_idx1 = space.space_to_hy_index(*START_PNT_SPC_1)

    end_ez_idx1 = space.space_to_ez_index(*END_PNT_SPC_1)
    end_hx_idx1 = space.space_to_hx_index(*END_PNT_SPC_1)
    end_hy_idx1 = space.space_to_hy_index(*END_PNT_SPC_1)

    ez_sums2 = []
    hx_sums2 = []
    hy_sums2 = []

    start_ez_idx2 = space.space_to_ez_index(*START_PNT_SPC_2)
    start_hx_idx2 = space.space_to_hx_index(*START_PNT_SPC_2)
    start_hy_idx2 = space.space_to_hy_index(*START_PNT_SPC_2)

    end_ez_idx2 = space.space_to_ez_index(*END_PNT_SPC_2)
    end_hx_idx2 = space.space_to_hx_index(*END_PNT_SPC_2)
    end_hy_idx2 = space.space_to_hy_index(*END_PNT_SPC_2)

    while my_fdtd.time_step.t < TIME_STEP_T:
        if space.my_id == 0 and my_fdtd.time_step.n % 100 == 0:
            print my_fdtd.time_step.t

        ez_sum1, hx_sum1, hy_sum1 = 0, 0, 0
        ez_sum2, hx_sum2, hy_sum2 = 0, 0, 0

        for i in range(start_ez_idx1[1], end_ez_idx1[1] + 1):
            idx1 = start_ez_idx1[0], i, start_ez_idx1[2]
            if geometry.in_range(idx1, my_fdtd.ez.shape, constant.Ez):
                ez_sum1 += my_fdtd.ez[idx1]

        for i in range(start_hx_idx1[1], end_hx_idx1[1] + 1):
            idx1 = start_hx_idx1[0], i, start_hx_idx1[2]
            if geometry.in_range(idx1, my_fdtd.hx.shape, constant.Hx):
                hx_sum1 += my_fdtd.hx[idx1]

        for i in range(start_hy_idx1[1], end_hx_idx1[1] + 1):
            idx1 = start_hy_idx1[0], i, start_hy_idx1[2]
            if geometry.in_range(idx1, my_fdtd.hy.shape, constant.Hy):
                hy_sum1 += my_fdtd.hy[idx1]

        for i in range(start_ez_idx2[1], end_ez_idx2[1] + 1):
            idx2 = start_ez_idx2[0], i, start_ez_idx2[2]
            if geometry.in_range(idx2, my_fdtd.ez.shape, constant.Ez):
                ez_sum2 += my_fdtd.ez[idx2]

        for i in range(start_hx_idx2[1], end_hx_idx2[1] + 1):
            idx2 = start_hx_idx2[0], i, start_hx_idx2[2]
            if geometry.in_range(idx2, my_fdtd.hx.shape, constant.Hx):
                hx_sum2 += my_fdtd.hx[idx2]

        for i in range(start_hy_idx2[1], end_hx_idx2[1] + 1):
            idx2 = start_hy_idx2[0], i, start_hy_idx2[2]
            if geometry.in_range(idx2, my_fdtd.hy.shape, constant.Hy):
                hy_sum2 += my_fdtd.hy[idx2]

        ez_sum1 = space.cart_comm.reduce(ez_sum1, root=ez1_id)
        hx_sum1 = space.cart_comm.reduce(hx_sum1, root=hx1_id)
        hy_sum1 = space.cart_comm.reduce(hy_sum1, root=hy1_id)

        ez_sum2 = space.cart_comm.reduce(ez_sum2, root=ez2_id)
        hx_sum2 = space.cart_comm.reduce(hx_sum2, root=hx2_id)
        hy_sum2 = space.cart_comm.reduce(hy_sum2, root=hy2_id)

        ez_sums1.append(ez_sum1)
        hx_sums1.append(hx_sum1)
        hy_sums1.append(hy_sum1)

        ez_sums2.append(ez_sum2)
        hx_sums2.append(hx_sum2)
        hy_sums2.append(hy_sum2)

        my_fdtd.step()

    if space.my_id == ez1_id:
        ez_sums1_file = open(FILENAME_EZ_SUMS_1, 'w')
        cPickle.dump(ez_sums1, ez_sums1_file)
        ez_sums1_file.close()

    if space.my_id == hx1_id:
        hx_sums1_file = open(FILENAME_HX_SUMS_1, 'w')
        cPickle.dump(hx_sums1, hx_sums1_file)
        hx_sums1_file.close()

    if space.my_id == hy1_id:
        hy_sums1_file = open(FILENAME_HY_SUMS_1, 'w')
        cPickle.dump(hy_sums1, hy_sums1_file)
        hy_sums1_file.close()

    if space.my_id == ez2_id:
        ez_sums2_file = open(FILENAME_EZ_SUMS_2, 'w')
        cPickle.dump(ez_sums2, ez_sums2_file)
        ez_sums2_file.close()

    if space.my_id == hx2_id:
        hx_sums2_file = open(FILENAME_HX_SUMS_2, 'w')
        cPickle.dump(hx_sums2, hx_sums2_file)
        hx_sums2_file.close()

    if space.my_id == hy2_id:
        hy_sums2_file = open(FILENAME_HY_SUMS_2, 'w')
        cPickle.dump(hy_sums2, hy_sums2_file)
        hy_sums2_file.close()

    if space.my_id == 0:
        print 'dt =', space.dt
        print time.time() - start, 'sec.'
