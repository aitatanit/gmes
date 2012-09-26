#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass

from copy import deepcopy
from math import sqrt
from cmath import exp as cexp
from numpy import ndindex, arange, inf, array
from datetime import datetime, timedelta

import numpy as np

# GMES modules
from geometry import GeomBoxTree, in_range, DefaultMedium
from file_io import Probe
#from file_io import write_hdf5, snapshot
from show import ShowLine, ShowPlane, Snapshot
from material import Dummy
from pygeom import GeomBox
from constant import *


class TimeStep(object):
    """Store the current time-step and time.
    
    Attributes:
    dt -- time-step size
    n -- current time-step
    t -- current time
    
    """
    def __init__(self, dt, n=0.0, t=0.0):
        """Constructor.
        
        Keyword arguments:
        dt -- time-step size (no default)
        n -- initial value of the time-step (default 0.0)
        t -- initial value of the time (default 0.0)
        
        """
        self.n = float(n)
        self.t = float(t)
        self.dt = float(dt)
        
    def half_step_up(self):
        """Increase n and t for the electric or magnetic field update.
        
        """
        self.n += 0.5
        self.t = self.n * self.dt


class FDTD(object):
    """three dimensional finite-difference time-domain class
    
    Attributes:
    space -- geometry.Cartesian instance
    cmplx -- Boolean of whether field is complex. Determined by the 
        space.period.
    dr -- space differentials: dx, dy, dz
    dt -- time-step size
    courant_ratio -- the ratio of dt to Courant stability bound
    bloch -- Bloch wave vector
    e_field_component
    h_field_component -- component list of the
    time_step -- an instance of the TimeStep class
    ex, ey, ez -- arrays for the electric fields
    hx, hy, hz -- arrays for the magnetic fields
    worker -- list of the update methods
    chatter -- list of the syncronizaion methods

    """
    def __init__(self, space=None, geom_list=None, src_list=None,
                 courant_ratio=.99, dt=None, bloch=None, verbose=True):
        """Constructor.
        
        Keyword arguments:
        space -- an instance which represents the coordinate system 
            (default None)
        geom_list -- a list which represents the geometric structure
            (default None)
        src_list -- a list of source instances
            (default None)
        courant_ratio -- the ratio of dt to Courant stability bound 
            (default 0.99)
        dt -- time-step size. If None is given, dt is calculated using space 
            differentials and courant_ratio. (default None)
        bloch -- Bloch wave vector (default None)
        verbose -- whether it prints the details (default True)

        """
        self._init_field_compnt()
        
        self._updater = {Ex: self.update_ex, Ey: self.update_ey,
                         Ez: self.update_ez, Hx: self.update_hx,
                         Hy: self.update_hy, Hz: self.update_hz}
        
        self._chatter = {Ex: self.talk_with_ex_neighbors,
                         Ey: self.talk_with_ey_neighbors,
                         Ez: self.talk_with_ez_neighbors,
                         Hx: self.talk_with_hx_neighbors,
                         Hy: self.talk_with_hy_neighbors,
                         Hz: self.talk_with_hz_neighbors}
        
        self.e_recorder = []
        self.h_recorder = []

        self.verbose = bool(verbose)

        self.space = space
                
        self._fig_id = int(self.space.my_id)
            
        self.dx, self.dy, self.dz = self.space.dr

        default_medium = (i for i in geom_list 
                          if isinstance(i, DefaultMedium)).next()
        eps_inf = default_medium.material.eps_inf
        mu_inf = default_medium.material.mu_inf
        dt_limit = self._dt_limit(space, eps_inf, mu_inf)

        if dt is None:
            self.courant_ratio = float(courant_ratio)
            time_step_size = self.courant_ratio * dt_limit
        else:
            time_step_size = float(dt)
            self.courant_ratio = time_step_size / dt_limit

        # Some codes in geometry.py and source.py use dt.
        self.space.dt = time_step_size

        self.time_step = TimeStep(time_step_size)
        
        if self.verbose:
            self.space.display_info()

        if self.verbose:
            print 'dt:', time_step_size
            print 'courant ratio:', self.courant_ratio
            
        if self.verbose:
            print 'Initializing the geometry list...',
            
        self.geom_list = geom_list

        for go in self.geom_list:
            go.init(self.space)
            
        if self.verbose:
            print 'done.'
            
        if self.verbose:
            print 'Generating geometric binary search tree...',
            
        self.geom_tree = GeomBoxTree(self.geom_list)

        if self.verbose:
            print 'done.'
            
        if self.verbose:
            print 'The geometric tree follows...'
            self.geom_tree.display_info()
                
        if bloch is None:
            self.cmplx = False
            self.bloch = None
        else:
            self.cmplx = True

            self.bloch = np.array(bloch, np.double)

            if self.verbose:
                print 'Bloch wave vector is', self.bloch

            # Calculate accurate a Bloch wave vector using equation 4.14a 
            # at p.112 of 'A. Taflove and S. C. Hagness, Computational 
            # Electrodynamics: The Finite-Difference Time-Domain Method, 
            # Third Edition, 3rd ed. Artech House Publishers, 2005'.
            ds = np.array((self.dx, self.dy, self.dz))
            default_medium = (i for i in geom_list 
                                if isinstance(i, DefaultMedium)).next()
            eps_inf = default_medium.material.eps_inf
            mu_inf = default_medium.material.mu_inf
            c = 1 / sqrt(eps_inf * mu_inf)
            S = c * time_step_size / ds
            ref_n = sqrt(eps_inf)
            numeric_bloch = 2 * ref_n / ds * np.arcsin(np.sin(self.bloch * S * ds / 2) / S)

            if self.verbose:
                print 'numerical Bloch wave vector is', numeric_bloch
            
        if self.verbose:
            print 'Initializing source...',
            
        self.src_list = src_list
        for so in self.src_list:
            so.init(self.geom_tree, self.space, self.cmplx)
            
        if self.verbose:
            print 'done.'
            
        if self.verbose:
            print 'The source list information follows...'
            for so in self.src_list:
                so.display_info()

        self.pw_material = {}

    def init(self):
        """Initialize sources.

        """
        st = datetime.now()
        
        if self.verbose:
            print 'Allocating memory for the electromagnetic fields...',
            
        # storage for the electromagnetic field 
        self.ex = self.space.get_ex_storage(self.e_field_compnt, self.cmplx)
        self.ey = self.space.get_ey_storage(self.e_field_compnt, self.cmplx)
        self.ez = self.space.get_ez_storage(self.e_field_compnt, self.cmplx)
        self.hx = self.space.get_hx_storage(self.h_field_compnt, self.cmplx)
        self.hy = self.space.get_hy_storage(self.h_field_compnt, self.cmplx)
        self.hz = self.space.get_hz_storage(self.h_field_compnt, self.cmplx)
        
        self.field = {Ex: self.ex, Ey: self.ey, Ez: self.ez,
                      Hx: self.hx, Hy: self.hy, Hz: self.hz}

        if self.verbose:
            print 'done.'
            
        if self.verbose:
            for comp in self.field:
                print comp.__name__, 'field:', self.field[comp].dtype, self.field[comp].shape 

        if self.verbose:
            print 'Mapping the piecewise material.',
            print 'This will take some times...'

        self.init_material()
        
        if self.verbose:
            print 'Mapping the pointwise source...',
            
        self.pw_source = {}

        self.init_source()

        et = datetime.now()
        print 'Elapsed time:', (et - st)
        
    def _print_pw_obj(self, pw_obj):
        """Print information of the piecewise material and source.

        """
        for o in pw_obj.itervalues():
            print o.name(), 'at', o.idx_size(), 'point(s).',
        if len(pw_obj):
            print
        else:
            print None

    def _init_field_compnt(self):
        """Set the significant electromagnetic field components.
        
        """
        self.e_field_compnt = (Ex, Ey, Ez)
        self.h_field_compnt = (Hx, Hy, Hz)

    def _dt_limit(self, space, eps_inf, mu_inf):
        """Courant stability bound of a time-step.

        Equation 4.54b at p.131 of A. Taflove and S. C. Hagness, Computational
        Electrodynamics: The Finite-Difference Time-Domain Method, Third 
        Edition, 3rd ed. Artech House Publishers, 2005.

        """
        # Pick out meaningful space-differential(s).
        # 0 for dx, 1 for dy, and 2 for dz.
        ds = {Ex: (1, 2), Ey: (2, 0), Ez: (0, 1),
              Hx: (1, 2), Hy: (2, 0), Hz: (0, 1)}
        
        e_ds = set()
        for i in self.e_field_compnt:
            e_ds.update(ds[i])
        h_ds = set()
        for i in self.h_field_compnt:
            h_ds.update(ds[i])

        non_inf = e_ds.intersection(h_ds)

        dr = array(space.dr, np.double)
        for i in range(3):
            if i not in non_inf: dr[i] = inf
        
        c = 1 / sqrt(eps_inf * mu_inf)
        return 1 / c / sqrt(sum(dr**-2))
        
    def _step_aux_fdtd(self):
        for src in self.src_list:
            src.step()
        
    def __deepcopy__(self, memo={}):
        """The classes generated by swig do not have __deepcopy__ method.
        Thus, the tricky constructor call follows.
        
        """
        newcopy = self.__class__(space=self.space, 
                                 geom_list=self.geom_list, 
                                 src_list=self.src_list, 
                                 courant_ratio=self.courant_ratio, 
                                 dt=self.time_step.dt, 
                                 wavevector=self.wavevector, 
                                 verbose=self.verbose)

        newcopy.ex = np.array(self.ex)
        newcopy.ey = np.array(self.ey)
        newcopy.ez = np.array(self.ez)
        newcopy.hx = np.array(self.hx)
        newcopy.hy = np.array(self.hy)
        newcopy.hz = np.array(self.hz)
        
        newcopy.time_step = deepcopy(self.time_step)
        return newcopy
    	
    def init_material_ex(self):
        """Set up the update mechanism for Ex field.
        
        Set up the update mechanism for Ex field and stores the result
        at self.pw_material[Ex].
        
        """
        self.pw_material[Ex] = {}
        shape = self.ex.shape
        for idx in ndindex(shape):
            spc = self.space.ex_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[1] == shape[1] - 1 or idx[2] == shape[2] - 1:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = mat_obj.get_pw_material_ex(idx, spc, underneath, self.cmplx)
            
            if self.pw_material[Ex].has_key(type(pw_obj)):
                self.pw_material[Ex][type(pw_obj)].merge(pw_obj)
            else:
                self.pw_material[Ex][type(pw_obj)] = pw_obj

    def init_material_ey(self):
        """Set up the update mechanism for Ey field.
        
        Set up the update mechanism for Ey field and stores the result
        at self.pw_material[Ey].
        
        """
        self.pw_material[Ey] = {}
        shape = self.ey.shape
        for idx in ndindex(shape):
            spc = self.space.ey_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[2] == shape[2] - 1 or idx[0] == shape[0] - 1:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = mat_obj.get_pw_material_ey(idx, spc, underneath, self.cmplx)

            if self.pw_material[Ey].has_key(type(pw_obj)):
                self.pw_material[Ey][type(pw_obj)].merge(pw_obj)
            else:
                self.pw_material[Ey][type(pw_obj)] = pw_obj

    def init_material_ez(self):
        """Set up the update mechanism for Ez field.
        
        Set up the update mechanism for Ez field and stores the result
        at self.pw_material[Ez].
        
        """
        self.pw_material[Ez] = {}
        shape = self.ez.shape
        for idx in ndindex(shape):
            spc = self.space.ez_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[0] == shape[0] - 1 or idx[1] == shape[1] - 1:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = mat_obj.get_pw_material_ez(idx, spc, underneath, self.cmplx)

            if self.pw_material[Ez].has_key(type(pw_obj)):
                self.pw_material[Ez][type(pw_obj)].merge(pw_obj)
            else:
                self.pw_material[Ez][type(pw_obj)] = pw_obj

    def init_material_hx(self):
        """Set up the update mechanism for Hx field.
        
        Set up the update mechanism for Hx field and stores the result
        at self.pw_material[Hx].
        
        """
        self.pw_material[Hx] = {}
        shape = self.hx.shape
        for idx in ndindex(shape):
            spc = self.space.hx_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[1] == 0 or idx[2] == 0:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = mat_obj.get_pw_material_hx(idx, spc, underneath, self.cmplx)

            if self.pw_material[Hx].has_key(type(pw_obj)):
                self.pw_material[Hx][type(pw_obj)].merge(pw_obj)
            else:
                self.pw_material[Hx][type(pw_obj)] = pw_obj

    def init_material_hy(self):
        """Set up the update mechanism for Hy field.
        
        Set up the update mechanism for Hy field and stores the result
        at self.pw_material[Hy].
        
        """
        self.pw_material[Hy] = {}
        shape = self.hy.shape
        for idx in ndindex(shape):
            spc = self.space.hy_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[2] == 0 or idx[0] == 0:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = mat_obj.get_pw_material_hy(idx, spc, underneath, self.cmplx)

            if self.pw_material[Hy].has_key(type(pw_obj)):
                self.pw_material[Hy][type(pw_obj)].merge(pw_obj)
            else:
                self.pw_material[Hy][type(pw_obj)] = pw_obj

    def init_material_hz(self):
        """Set up the update mechanism for Hz field.
        
        Set up the update mechanism for Hz field and stores the result
        at self.pw_material[Hz].
        
        """
        self.pw_material[Hz] = {}
        shape = self.hz.shape
        for idx in ndindex(shape):
            spc = self.space.hz_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[0] == 0 or idx[1] == 0:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = mat_obj.get_pw_material_hz(idx, spc, underneath, self.cmplx)

            if self.pw_material[Hz].has_key(type(pw_obj)):
                self.pw_material[Hz][type(pw_obj)].merge(pw_obj)
            else:
                self.pw_material[Hz][type(pw_obj)] = pw_obj

    def init_material(self):
        init_mat_func = {Ex: self.init_material_ex,
                         Ey: self.init_material_ey,
                         Ez: self.init_material_ez,
                         Hx: self.init_material_hx,
                         Hy: self.init_material_hy,
                         Hz: self.init_material_hz}
        
        for comp in self.e_field_compnt:
            if self.verbose:
                print 'Mapping materials for', comp.__name__, 'field'
            init_mat_func[comp]()
            
            if self.verbose:
                self._print_pw_obj(self.pw_material[comp])

        for comp in self.h_field_compnt:
            if self.verbose:
                print 'Mapping materials for', comp.__name__, 'field'
            init_mat_func[comp]()

            if self.verbose:
                self._print_pw_obj(self.pw_material[comp])

    def init_source_ex(self):
        self.pw_source[Ex] = {}
        for so in self.src_list:
            pw_src = so.get_pw_source_ex(self.ex, self.space, self.geom_tree)

            if pw_src is None:
                continue

            if self.pw_source[Ex].has_key(type(pw_src)):
                self.pw_source[Ex][type(pw_src)].merge(pw_src)
            else:
                self.pw_source[Ex][type(pw_src)] = pw_src
            
    def init_source_ey(self):
        self.pw_source[Ey] = {}
        for so in self.src_list:
            pw_src = so.get_pw_source_ey(self.ey, self.space, self.geom_tree)
            
            if pw_src is None:
                continue

            if self.pw_source[Ey].has_key(type(pw_src)):
                self.pw_source[Ey][type(pw_src)].merge(pw_src)
            else:
                self.pw_source[Ey][type(pw_src)] = pw_src

    def init_source_ez(self):
        self.pw_source[Ez] = {}
        for so in self.src_list:
            pw_src = so.get_pw_source_ez(self.ez, self.space, self.geom_tree)

            if pw_src is None:
                continue

            if self.pw_source[Ez].has_key(type(pw_src)):
                self.pw_source[Ez][type(pw_src)].merge(pw_src)
            else:
                self.pw_source[Ez][type(pw_src)] = pw_src

    def init_source_hx(self):
        self.pw_source[Hx] = {}
        for so in self.src_list:
            pw_src = so.get_pw_source_hx(self.hx, self.space, self.geom_tree)

            if pw_src is None:
                continue

            if self.pw_source[Hx].has_key(type(pw_src)):
                self.pw_source[Hx][type(pw_src)].merge(pw_src)
            else:
                self.pw_source[Hx][type(pw_src)] = pw_src
            
    def init_source_hy(self):
        self.pw_source[Hy] = {}
        for so in self.src_list:
            pw_src = so.get_pw_source_hy(self.hy, self.space, self.geom_tree)

            if pw_src is None:
                continue

            if self.pw_source[Hy].has_key(type(pw_src)):
                self.pw_source[Hy][type(pw_src)].merge(pw_src)
            else:
                self.pw_source[Hy][type(pw_src)] = pw_src
            
    def init_source_hz(self):
        self.pw_source[Hz] = {}
        for so in self.src_list:
            pw_src = so.get_pw_source_hz(self.hz, self.space, self.geom_tree)

            if pw_src is None:
                continue

            if self.pw_source[Hz].has_key(type(pw_src)):
                self.pw_source[Hz][type(pw_src)].merge(pw_src)
            else:
                self.pw_source[Hz][type(pw_src)] = pw_src
            
    def init_source(self):
        init_src_func = {Ex: self.init_source_ex,
                         Ey: self.init_source_ey,
                         Ez: self.init_source_ez,
                         Hx: self.init_source_hx,
                         Hy: self.init_source_hy,
                         Hz: self.init_source_hz}

        for comp in self.e_field_compnt:
            if self.verbose:
                print 'Mapping sources for', comp.__name__, 'field.'
            init_src_func[comp]()
            if self.verbose:
                self._print_pw_obj(self.pw_source[comp])
 
        for comp in self.h_field_compnt:
            if self.verbose:
                print 'Mapping sources for', comp.__name__, 'field.'
            init_src_func[comp]()
            if self.verbose:
                self._print_pw_obj(self.pw_source[comp])
        
    def set_probe(self, p, prefix=None):
        """
        p: space coordinates. type: tuple-3
        prefix: prefix of the recording file name. type: str
        
        """
        spc2idx = {Ex: self.space.space_to_ex_index,
                   Ey: self.space.space_to_ey_index,
                   Ez: self.space.space_to_ez_index,
                   Hx: self.space.space_to_hx_index,
                   Hy: self.space.space_to_hy_index,
                   Hz: self.space.space_to_hz_index}
        
        idx2spc = {Ex: self.space.ex_index_to_space,
                   Ey: self.space.ey_index_to_space,
                   Ez: self.space.ez_index_to_space,
                   Hx: self.space.hx_index_to_space,
                   Hy: self.space.hy_index_to_space,
                   Hz: self.space.hz_index_to_space}

        validity = {Ex: 
                    (lambda idx: in_range(idx, self.ex.shape, Ex)),
                    Ey: 
                    (lambda idx: in_range(idx, self.ey.shape, Ey)),
                    Ez: 
                    (lambda idx: in_range(idx, self.ez.shape, Ez)),
                    Hx:
                    (lambda idx: in_range(idx, self.hx.shape, Hx)),
                    Hy: 
                    (lambda idx: in_range(idx, self.hy.shape, Hy)),
                    Hz: 
                    (lambda idx: in_range(idx, self.hz.shape, Hz))}

        postfix = {Ex: '_ex.dat', Ey: '_ey.dat',
                   Ez: '_ez.dat', Hx: '_hx.dat',
                   Hy: '_hy.dat', Hz: '_hz.dat'}
        
        for comp in self.e_field_compnt:
            idx = spc2idx[comp](*p)
            if validity[comp](idx):
                if prefix:
                    filename = prefix + postfix[comp]
                else:
                    filename = postfix[comp]
                recorder = Probe(idx, self.field[comp], filename)
                loc = idx2spc[comp](*idx)
                recorder.write_header(loc, self.time_step.dt)
                self.e_recorder.append(recorder)

        for comp in self.h_field_compnt:
            idx = spc2idx[comp](*p)
            if validity[comp](idx):
                if prefix:
                    filename = prefix + postfix[comp]
                else:
                    filename = postfix[comp]
                recorder = Probe(idx, self.field[comp], filename)
                loc = idx2spc[comp](*idx)
                recorder.write_header(loc, self.time_step.dt)
                self.h_recorder.append(recorder)

    def update_ex(self):
        for pw_obj in self.pw_material[Ex].itervalues():
            pw_obj.update_all(self.ex, self.hz, self.hy, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.pw_source[Ex].itervalues():
            pw_obj.update_all(self.ex, self.hz, self.hy, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)
        
    def update_ey(self):
        for pw_obj in self.pw_material[Ey].itervalues():
            pw_obj.update_all(self.ey, self.hx, self.hz, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)
		
        for pw_obj in self.pw_source[Ey].itervalues():
            pw_obj.update_all(self.ey, self.hx, self.hz, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)

    def update_ez(self):
        for pw_obj in self.pw_material[Ez].itervalues():
            pw_obj.update_all(self.ez, self.hy, self.hx, self.dx, self.dy,
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.pw_source[Ez].itervalues():
            pw_obj.update_all(self.ez, self.hy, self.hx, self.dx, self.dy,
                              self.time_step.dt, self.time_step.n)
        
    def update_hx(self):
        for pw_obj in self.pw_material[Hx].itervalues():
            pw_obj.update_all(self.hx, self.ez, self.ey, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.pw_source[Hx].itervalues():
            pw_obj.update_all(self.hx, self.ez, self.ey, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)
		
    def update_hy(self):
        for pw_obj in self.pw_material[Hy].itervalues():
            pw_obj.update_all(self.hy, self.ex, self.ez, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.pw_source[Hy].itervalues():
            pw_obj.update_all(self.hy, self.ex, self.ez, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)
		
    def update_hz(self):
        for pw_obj in self.pw_material[Hz].itervalues():
            pw_obj.update_all(self.hz, self.ey, self.ex, self.dx, self.dy, 
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.pw_source[Hz].itervalues():
            pw_obj.update_all(self.hz, self.ey, self.ex, self.dx, self.dy, 
                              self.time_step.dt, self.time_step.n)

    def talk_with_ex_neighbors(self):
        """Synchronize ex data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # Send ex field data to -y direction and receive from +y direction.
        src, dest = self.space.cart_comm.Shift(1, -1)
        if dest == -1 or src == -1:
            return

        if self.cmplx:
            dest_spc = self.space.ex_index_to_space(0, self.ex.shape[1] - 1, 0)[1]
        
            src_spc = self.space.ex_index_to_space(0, 0, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Ex.tag,
                                                    None, src, Ex.tag)
            
            phase_shift = cexp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.ex[:, -1, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ex[:, 0, :], dest, Ex.tag,
                                      None, src, Ex.tag)
        
        # Send ex field data to -z direction and receive from +z direction.
        src, dest = self.space.cart_comm.Shift(2, -1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.ex_index_to_space(0, 0, self.ex.shape[2] - 1)[2]
        
            src_spc = self.space.ex_index_to_space(0, 0, 0)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Ex.tag,
                                                    None, src, Ex.tag)
        
            phase_shift = cexp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.ex[:, :, -1] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ex[:, :, 0], dest, Ex.tag,
                                      None, src, Ex.tag)
        
    def talk_with_ey_neighbors(self):
        """Synchronize ey data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # Send ey field data to -z direction and receive from +z direction.
        src, dest = self.space.cart_comm.Shift(2, -1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.ey_index_to_space(0, 0, self.ey.shape[2] - 1)[2]
        
            src_spc = self.space.ey_index_to_space(0, 0, 0)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Ey.tag,
                                                    None, src, Ey.tag) 
        
            phase_shift = cexp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 1
            
        self.ey[:, :, -1] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ey[:, :, 0], dest, Ey.tag,
                                      None, src, Ey.tag)
        
        # Send ey field data to -x direction and receive from +x direction.
        src, dest = self.space.cart_comm.Shift(0, -1)
        if dest == -1 or src == -1:
            return

        if self.cmplx:
            dest_spc = self.space.ey_index_to_space(self.ey.shape[0] - 1, 0, 0)[0]
        
            src_spc = self.space.ey_index_to_space(0, 0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Ey.tag,
                                                    None, src, Ey.tag)
            
            phase_shift = cexp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 1
            
        self.ey[-1, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ey[0, :, :], dest, Ey.tag,
                                      None, src, Ey.tag)
        
    def talk_with_ez_neighbors(self):
        """Synchronize ez data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # Send ez field data to -x direction and receive from +x direction.
        src, dest = self.space.cart_comm.Shift(0, -1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.ez_index_to_space(self.ez.shape[0] - 1, 0, 0)[0]
        
            src_spc = self.space.ez_index_to_space(0, 0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Ez.tag,
                                                    None, src, Ez.tag) 
            
            phase_shift = cexp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.ez[-1, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ez[0, :, :], dest, Ez.tag,
                                      None, src, Ez.tag)
        
        # Send ez field data to -y direction and receive from +y direction.
        src, dest = self.space.cart_comm.Shift(1, -1)
        if dest == -1 or src == -1:
            return

        if self.cmplx:
            dest_spc = self.space.ez_index_to_space(0, self.ez.shape[1] - 1, 0)[1]
        
            src_spc = self.space.ez_index_to_space(0, 0, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Ez.tag,
                                                    None, src, Ez.tag) 
            
            phase_shift = cexp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 1

        self.ez[:, -1, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ez[:, 0, :], dest, Ez.tag,
                                      None, src, Ez.tag)
        
    def talk_with_hx_neighbors(self):
        """Synchronize hx data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # Send hx field data to +y direction and receive from -y direction.
        src, dest = self.space.cart_comm.Shift(1, 1)
        if dest == -1 or src == -1:
            return

        if self.cmplx:
            dest_spc = self.space.hx_index_to_space(0, 0, 0)[1]
        
            src_spc = self.space.hx_index_to_space(0, self.hx.shape[1] - 1, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Hx.tag,
                                                    None, src, Hx.tag)
        
            phase_shift = cexp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.hx[:, 0, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hx[:, -1, :], dest, Hx.tag,
                                      None, src, Hx.tag)
            
        # Send hx field data to +z direction and receive from -z direction.    
        src, dest = self.space.cart_comm.Shift(2, 1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.hx_index_to_space(0, 0, 0)[2]
        
            src_spc = self.space.hx_index_to_space(0, 0, self.hx.shape[2] - 1)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Hx.tag,
                                                    None, src, Hx.tag)
        
            phase_shift = cexp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.hx[:, :, 0] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hx[:, :, -1], dest, Hx.tag,
                                      None, src, Hx.tag)
            
    def talk_with_hy_neighbors(self):
        """Synchronize hy data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # Send hy field data to +z direction and receive from -z direction.
        src, dest = self.space.cart_comm.Shift(2, 1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.hy_index_to_space(0, 0, 0)[2]
        
            src_spc = self.space.hy_index_to_space(0, 0, self.hy.shape[2] - 1)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Hy.tag,
                                                    None, src, Hy.tag)
        
            phase_shift = cexp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.hy[:, :, 0] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hy[:, :, -1], dest, Hy.tag,
                                      None, src, Hy.tag)
            
        # Send hy field data to +x direction and receive from -x direction.
        src, dest = self.space.cart_comm.Shift(0, 1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.hy_index_to_space(0, 0, 0)[0]
        
            src_spc = self.space.hy_index_to_space(self.hy.shape[0] - 1, 0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Hy.tag,
                                                    None, src, Hy.tag)
        
            phase_shift = cexp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.hy[0, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hy[-1, :, :], dest, Hy.tag,
                                      None, src, Hy.tag)
            
    def talk_with_hz_neighbors(self):
        """Synchronize hz data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # Send hz field data to +x direction and receive from -x direction.
        src, dest = self.space.cart_comm.Shift(0, 1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.hz_index_to_space(0, 0, 0)[0]
        
            src_spc = self.space.hz_index_to_space(self.hz.shape[0] - 1, 0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Hz.tag,
                                                    None, src, Hz.tag)
        
            phase_shift = cexp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.hz[0, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hz[-1, :, :], dest, Hz.tag,
                                      None, src, Hz.tag)
            
        # Send hz field data to +y direction and receive from -y direction.
        src, dest = self.space.cart_comm.Shift(1, 1)
        if dest == -1 or src == -1:
            return
        
        if self.cmplx:
            dest_spc = self.space.hz_index_to_space(0, 0, 0)[1]
        
            src_spc = self.space.hz_index_to_space(0, self.hz.shape[1] - 1, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, Hz.tag,
                                                    None, src, Hz.tag)
        
            phase_shift = cexp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 1
        
        self.hz[:, 0, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hz[:, -1, :], dest, Hz.tag,
                                      None, src, Hz.tag)
        
    def step(self):
        self.time_step.half_step_up()

        for comp in self.h_field_compnt:
            self._chatter[comp]()
            
        for comp in self.e_field_compnt:
            self._updater[comp]()

        for probe in self.e_recorder:
            probe.write(self.time_step.n)
            
        self.time_step.half_step_up()

        self._step_aux_fdtd()

        for comp in self.e_field_compnt:
            self._chatter[comp]()
        
        for comp in self.h_field_compnt:
            self._updater[comp]()

        for probe in self.h_recorder:
            probe.write(self.time_step.n)

    def step_while_zero(self, component, point, modulus=inf):
        """Run self.step() while the field value at the given point is 0.
        
        Keyword arguments:
        component: filed component to check.
            one of gmes.constant.{Ex, Ey, Ez, Hx, Hy, Hz}.
        point: coordinates of the location to check the field. 
            tuple of three scalars
        modulus: print n and t at every modulus steps.

        """
        spc_to_idx = {Ex: self.space.space_to_ex_index,
                      Ey: self.space.space_to_ey_index,
                      Ez: self.space.space_to_ez_index,
                      Hx: self.space.space_to_hx_index,
                      Hy: self.space.space_to_hy_index,
                      Hz: self.space.space_to_hz_index}

        idx = spc_to_idx[component](*point)

        if in_range(idx, self.field[component].shape, component):
            hot_node = self.space.my_id
        else:
            hot_node = None
            
        hot_node = self.space.bcast(hot_node, hot_node)
            
        flag = True
        while flag:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t
            if self.space.my_id == hot_node and self.field[component][idx] != 0:
                flag = False
            flag = self.space.cart_comm.bcast(flag, hot_node)

    def step_until_n(self, n=0, modulus=inf):
        """Run self.step() until time step reaches n.

        """
        st = datetime.now()

        if self.time_step.n < n:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t
                
            et = datetime.now()
            estimated_t = (n - 1) * (et - st).seconds
            print 'Estimated time of completion:', timedelta(seconds=estimated_t)
        
        while self.time_step.n < n:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t

        et = datetime.now()
        print 'Elapsed time:', (et - st)

    def step_until_t(self, t=0, modulus=inf):
        """Run self.step() until time reaches t.

        """
        st = datetime.now()
        sn = self.time_step.n

        if self.time_step.t < t:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t

            et = datetime.now()
            num_of_steps = (t - self.time_step.t) / self.time_step.dt
            estimated_t = num_of_steps * (et - st).seconds
            print 'Estimated time of completion:', timedelta(seconds=estimated_t)
        
        while self.time_step.t < t:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t
                
        et = datetime.now()
        en = self.time_step.n
        print 'Elapsed time:', (et - st),
        print '(%d timesteps)' % (en - sn)

    def show_field_line(self, comp, start, end, vrange=(-1,1), interval=2500):
        """Show the real value of the feild along the line.

        comp: field component
        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        title = {Ex: 'Ex field', Ey: 'Ey field', Ez: 'Ez field',
                 Hx: 'Hx field', Hy: 'Hy field', Hz: 'Hz field'}

        showcase = ShowLine(self, comp, start, end, vrange, interval, 
                            title[comp], self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
    
    def show_field(self, comp, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ex on the plone.

        comp: field component
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        title = {Ex: 'Ex field', Ey: 'Ey field', Ez: 'Ez field',
                 Hx: 'Hx field', Hy: 'Hy field', Hz: 'Hz field'}

        showcase = ShowLine(self, comp, start, end, vrange, interval, 
                            title[comp], self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
    
    def show_field(self, comp, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ex on the plone.

        comp: field component
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        title = {Ex: 'Ex field', Ey: 'Ey field', Ez: 'Ez field',
                 Hx: 'Hx field', Hy: 'Hy field', Hz: 'Hz field'}
        
        showcase = ShowPlane(self, comp, axis, cut, vrange, 
                             interval, title[comp], self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase

    def show_permittivity(self, comp, axis, cut, vrange=None):
        """Show permittivity for the ex on the plane.

        Keyword arguments:
        comp: field component
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. A tuple of length two.
        
        """
        title = {Ex: 'Permittivity for Ex field', 
                 Ey: 'Permittivity for Ey field', 
                 Ez: 'Permittivity for Ez field',
                 Hx: 'Permittivity for Hx field', 
                 Hy: 'Permittivity for Hy field', 
                 Hz: 'Permittivity for Hz field'}
        
        showcase = Snapshot(self, comp, axis, cut, vrange, 
                            title[comp], self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def write_field(self, comp, low=(-inf,-inf,-inf), high=(inf,inf,inf), prefix=None, postfix=None):
        """Dump the field values on a file.

        comp: Component of the field
        low: Coordinate of lower boundary
        high: Coordinate of higher boundary
        prefix: A prefix of the output filename. It should not contain
                parentheses.
        postfix: A postfix of the output filename. It should not contain
                 parentheses.

        The current version saves the field in npy format. The 
        format will be changed to hdf5 for the metadata.

        """
        spc2idx = {Ex: self.space.space_to_ex_index,
                   Ey: self.space.space_to_ey_index,
                   Ez: self.space.space_to_ez_index,
                   Hx: self.space.space_to_hx_index,
                   Hy: self.space.space_to_hy_index,
                   Hz: self.space.space_to_hz_index}
        
        idx2spc = {Ex: self.space.ex_index_to_space,
                   Ey: self.space.ey_index_to_space,
                   Ez: self.space.ez_index_to_space,
                   Hx: self.space.hx_index_to_space,
                   Hy: self.space.hy_index_to_space,
                   Hz: self.space.hz_index_to_space}

        low = array(low, np.double)
        high = array(high, np.double)
        
        for i, diff in enumerate(high - low):
            if 0 <= diff < self.space.dr[i]:
                high[i] += self.space.dr[i] / 2
                low[i] -= self.space.dr[i] / 2
            elif diff < 0:
                print 'ERROR'

        dump_volume = GeomBox(low, high)

        if issubclass(comp, Electric):
            low_bndry_idx = (0, 0, 0)
            if comp is Ex:
                high_bndry_idx = (self.field[comp].shape[0] - 1, 
                                  self.field[comp].shape[1] - 2,
                                  self.field[comp].shape[2] - 2)
            elif comp is Ey:
                high_bndry_idx = (self.field[comp].shape[0] - 2, 
                                  self.field[comp].shape[1] - 1,
                                  self.field[comp].shape[2] - 2)
            elif comp is Ez:
                high_bndry_idx = (self.field[comp].shape[0] - 2, 
                                  self.field[comp].shape[1] - 2,
                                  self.field[comp].shape[2] - 1)
        elif issubclass(comp, Magnetic):
            high_bndry_idx = [i - 1 for i in self.field[comp].shape]
            if comp is Hx:
                low_bndry_idx = (0, 1, 1)
            elif comp is Hy:
                low_bndry_idx = (1, 0, 1)
            elif comp is Hz:
                low_bndry_idx = (1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)

        low_bndry_pnt = idx2spc[comp](*low_bndry_idx)
        high_bndry_pnt = idx2spc[comp](*high_bndry_idx)

        field_volume = GeomBox(low_bndry_pnt, high_bndry_pnt)
        
        if field_volume.overlap(dump_volume):
            dump_volume.intersection(field_volume)
        else:
            return None
                
        low_idx = spc2idx[comp](*(dump_volume.low))
        high_idx = map(lambda i: i + 1, spc2idx[comp](*(dump_volume.high)))

        topo = self.space.cart_comm.topo[2]
        name = '%s_t=%f_(%d,%d,%d)' % (comp.str(), self.time_step.t, topo[0], topo[1], topo[2])
        if prefix is not None:
            name = prefix + name
        if postfix is not None:
            name = name + postfix
            
        np.save(name, self.field[comp][low_idx[0]: high_idx[0],
                                       low_idx[1]: high_idx[1],
                                       low_idx[2]: high_idx[2]])

    def write_field_all(self, low=(-inf,-inf,-inf), high=(inf,inf,inf), prefix=None, postfix=None):
        """Write the all current fields.
        
        """
        for comp in self.e_field_compnt:
            self.write_field(comp, low, high, prefix, postfix)
        for comp in self.h_field_compnt:
            self.write_field(comp, low, high, prefix, postfix)
        
    def snapshot_field(self, comp, axis, cut):
        """Take a graphical snapshot of a field.

        comp: field component
        axis: normal axis to the snapshot plane
        cut: coordinate of axis to take a snoptshot

        """
        spc2idx = {Ex: self.space.space_to_ex_index,
                   Ey: self.space.space_to_ey_index,
                   Ez: self.space.space_to_ez_index,
                   Hx: self.space.space_to_hx_index,
                   Hy: self.space.space_to_hy_index,
                   Hz: self.space.space_to_hz_index}
        
        if axis is X:
            cut_idx = spc2idx[comp](cut, 0, 0)[0]
            data = self.field[comp][cut_idx, :, :]
        elif axis is Y:
            cut_idx = spc2idx[comp](0, cut, 0)[1]
            data = self.field[comp][:, cut_idx, :]
        elif axis is Z:
            cut_idx = spc2idx[comp](0, 0, cut)[2]
            data = self.field[comp][:, :, cut_idx]
        else:
            raise TypeError

        filename = 't=' + str(self.time_step[1] * space.dt)
        snapshot(data, filename, comp)
        

class TExFDTD(FDTD):
    """2-D fdtd which has transverse-electric mode with respect to x.
    
    Assume that the structure and incident wave are uniform in 
    the x direction. TExFDTD updates only Ey, Ez, and Hx field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ey, Ez)
        self.h_field_compnt = (Hx,)

        
class TEyFDTD(FDTD):
    """2-D FDTD which has transverse-electric mode with respect to y.
    
    Assume that the structure and incident wave are uniform in the
    y direction. TEyFDTD updates only Ez, Ex, and Hy field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ex, Ez)
        self.h_field_compnt = (Hy,)


class TEzFDTD(FDTD):
    """2-D FDTD which has transverse-electric mode with respect to z.

    Assume that the structure and incident wave are uniform in the
    z direction. TEzFDTD updates only Ex, Ey, and Hz field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ex, Ey)
        self.h_field_compnt = (Hz,)


class TMxFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to x.

    Assume that the structure and incident wave are uniform in the
    x direction. TMxFDTD updates only Hy, Hz, and Ex field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ex,)
        self.h_field_compnt = (Hy, Hz)

        
class TMyFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to y

    Assume that the structure and incident wave are uniform in the
    y direction. TMyFDTD updates only Hz, Hx, and Ey field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ey,)
        self.h_field_compnt = (Hz, Hx)


class TMzFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to z
    
    Assume that the structure and incident wave are uniform in the
    z direction. TMzFDTD updates only Hx, Hy, and Ez field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ez,)
        self.h_field_compnt = (Hx, Hy)


class TEMxFDTD(FDTD):
    """y-polarized and x-directed one dimensional fdtd class

    Assume that the structure and incident wave are uniform in 
    transverse direction. TEMxFDTD updates only Ey and Hz field 
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ey,)
        self.h_field_compnt = (Hz,)

        
class TEMyFDTD(FDTD):
    """z-polarized and y-directed one dimensional fdtd class

    Assume that the structure and incident wave are uniform in 
    transverse direction. TEMyFDTD updates only Ez and Hx field
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ez,)
        self.h_field_compnt = (Hx,)

        
class TEMzFDTD(FDTD):
    """x-polarized and z-directed one dimensional fdtd class
    
    Assume that the structure and incident wave are uniform in 
    transverse direction. TEMzFDTD updates only Ex and Hy field
    components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (Ex,)
        self.h_field_compnt = (Hy,)


if __name__ == '__main__':
    from math import sin
    
    from numpy import inf
    
    from geometry import Cylinder, Cartesian
    from material import Dielectric
    
    low = Dielectric(index=1)
    hi = Dielectric(index=3)
    width_hi = low.eps_inf / (low.eps_inf + hi.eps_inf)
    space = Cartesian(size=[1, 1, 1])
    geom_list = [DefaultMedium(material=low), 
                 Cylinder(material=hi, 
                          axis=[1, 0, 0], 
                          radius=inf, 
                          height=width_hi)]
    
    fdtd = FDTD(space=space, geometry=geom_list)
    
    while True:
        fdtd.step()
        fdtd.ex[7, 7, 7] = sin(a.n)
        print a.n
