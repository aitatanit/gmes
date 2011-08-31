#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')

from copy import deepcopy
from math import sqrt
from cmath import exp

import numpy as np
from numpy import ndindex, arange, inf, array

# GMES modules
from geometry import GeomBoxTree, in_range, DefaultMedium
from file_io import Probe
#from file_io import write_hdf5, snapshot
from show import ShowLine, ShowPlane, Snapshot
from material import Dummy
import constant as const


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
    dx, dy, dz -- space differentials
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
        
        self._updater = {const.Ex: self.update_ex, const.Ey: self.update_ey,
                         const.Ez: self.update_ez, const.Hx: self.update_hx,
                         const.Hy: self.update_hy, const.Hz: self.update_hz}
        
        self._chatter = {const.Ex: self.talk_with_ex_neighbors,
                         const.Ey: self.talk_with_ey_neighbors,
                         const.Ez: self.talk_with_ez_neighbors,
                         const.Hx: self.talk_with_hx_neighbors,
                         const.Hy: self.talk_with_hy_neighbors,
                         const.Hz: self.talk_with_hz_neighbors}
        
        self.e_recorder = []
        self.h_recorder = []

        self.verbose = bool(verbose)

        self.space = space
                
        self._fig_id = int(self.space.my_id)
            
        self.dx = float(space.dx)
        self.dy = float(space.dy)
        self.dz = float(space.dz)

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
        
        if verbose:
            self.space.display_info()

        if verbose:
            print 'dt:', time_step_size
            print 'courant ratio:', self.courant_ratio
            
        if verbose:
            print 'Initializing the geometry list...',
            
        self.geom_list = deepcopy(geom_list)
        for geom_obj in self.geom_list:
            geom_obj.init(self.space)
            
        if verbose:
            print 'done.'
            
        if verbose:
            print 'Generating geometric binary search tree...',
            
        self.geom_tree = GeomBoxTree(self.geom_list)

        if verbose:
            print 'done.'
            
        if verbose:
            print 'The geometric tree follows...'
            self.geom_tree.display_info()
                
        if bloch is None:
            self.cmplx = False
            self.bloch = None
        else:
            self.cmplx = True

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
            self.bloch = np.array(bloch, float)
            self.bloch = 2 * ref_n / ds \
                * np.arcsin(np.sin(self.bloch * S * ds / 2) / S)
            
        if verbose:
            print 'Bloch wave vector is', self.bloch
            
        if verbose:
            print 'Initializing source...',
            
        self.src_list = deepcopy(src_list)
        for so in self.src_list:
            so.init(self.geom_tree, self.space, self.cmplx)
            
        if verbose:
            print 'done.'
            
        if verbose:
            print 'The source list information follows...'
            for so in self.src_list:
                so.display_info()
                
        if verbose:
            print 'Allocating memory for the electric & magnetic fields...',
            
        # storage for the electric & magnetic field 
        self.ex = space.get_ex_storage(self.e_field_compnt, self.cmplx)
        self.ey = space.get_ey_storage(self.e_field_compnt, self.cmplx)
        self.ez = space.get_ez_storage(self.e_field_compnt, self.cmplx)
        self.hx = space.get_hx_storage(self.h_field_compnt, self.cmplx)
        self.hy = space.get_hy_storage(self.h_field_compnt, self.cmplx)
        self.hz = space.get_hz_storage(self.h_field_compnt, self.cmplx)
        
        if verbose:
            print 'done.'
            
        if verbose:
            print 'ex field:', self.ex.dtype, self.ex.shape 
            print 'ey field:', self.ey.dtype, self.ey.shape 
            print 'ez field:', self.ez.dtype, self.ez.shape
            print 'hx field:', self.hx.dtype, self.hx.shape
            print 'hy field:', self.hy.dtype, self.hy.shape
            print 'hz field:', self.hz.dtype, self.hz.shape
        
        if verbose:
            print 'Mapping the pointwise material.',
            print 'This will take some tiems...'

        self.material_ex = {}
        self.material_ey = {}
        self.material_ez = {}
        self.material_hx = {}
        self.material_hy = {}
        self.material_hz = {}

        self.init_material()
        
        if verbose:
            print 'done.'

        # pw_material information for electric & magnetic fields
        if verbose:
            print 'ex material:',
            self._print_pw_obj(self.material_ex)

            print 'ey material:',
            self._print_pw_obj(self.material_ey)

            print 'ez material:',
            self._print_pw_obj(self.material_ez)

            print 'hx material:',
            self._print_pw_obj(self.material_hx)

            print 'hy material:',
            self._print_pw_obj(self.material_hy)

            print 'hz material:',
            self._print_pw_obj(self.material_hz)

        if verbose:
            print 'Mapping the pointwise source...',
            
        self.source_ex = {}
        self.source_ey = {}
        self.source_ez = {}
        self.source_hx = {}
        self.source_hy = {}
        self.source_hz = {}

        self.init_source()

        if verbose:
            print 'done.'
    
        # pw_source information for electric & magnetic fields
        if verbose:
            print 'ex source:',
            self._print_pw_obj(self.source_ex)

            print 'ey source:',
            self._print_pw_obj(self.source_ey)

            print 'ez source:',
            self._print_pw_obj(self.source_ez)

            print 'hx source:',
            self._print_pw_obj(self.source_hx)

            print 'hy source:',
            self._print_pw_obj(self.source_hy)

            print 'hz source:',
            self._print_pw_obj(self.source_hz)

    def _print_pw_obj(self, pw_obj):
        print_count = 0
        for o in pw_obj.itervalues():
            print type(o), 'with', o.idx_size(), 'point(s).'
            print_count += 1
        if print_count == 0: print

    def _init_field_compnt(self):
        """Set the significant electromagnetic field components.
        
        """
        self.e_field_compnt = (const.Ex, const.Ey, const.Ez)
        self.h_field_compnt = (const.Hx, const.Hy, const.Hz)

    def _dt_limit(self, space, eps_inf, mu_inf):
        """Courant stability bound of a time-step.

        Equation 4.54b at p.131 of A. Taflove and S. C. Hagness, Computational
        Electrodynamics: The Finite-Difference Time-Domain Method, Third 
        Edition, 3rd ed. Artech House Publishers, 2005.

        """
        # Pick out meaningful space-differential(s).
        # 0 for dx, 1 for dy, and 2 for dz.
        ds = {const.Ex: (1, 2), const.Ey: (2, 0), const.Ez: (0, 1),
              const.Hx: (1, 2), const.Hy: (2, 0), const.Hz: (0, 1)}
        
        e_ds = set()
        for i in self.e_field_compnt:
            e_ds.update(ds[i])
        h_ds = set()
        for i in self.h_field_compnt:
            h_ds.update(ds[i])

        non_inf = e_ds.intersection(h_ds)

        dr = array((space.dx, space.dy, space.dz), float)
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
        at self.material_ex.
        
        """
        shape = self.ex.shape
        for idx in ndindex(shape):
            spc = self.space.ex_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[1] == shape[1] - 1 or idx[2] == shape[2] - 1:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = \
                mat_obj.get_pw_material_ex(idx, spc, underneath, self.cmplx)
            
            if self.material_ex.has_key(type(pw_obj)):
                self.material_ex[type(pw_obj)].merge(pw_obj)
            else:
                self.material_ex[type(pw_obj)] = pw_obj

    def init_material_ey(self):
        """Set up the update mechanism for Ey field.
        
        Set up the update mechanism for Ey field and stores the result
        at self.material_ey.
        
        """
        shape = self.ey.shape
        for idx in ndindex(shape):
            spc = self.space.ey_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[2] == shape[2] - 1 or idx[0] == shape[0] - 1:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = \
                mat_obj.get_pw_material_ey(idx, spc, underneath, self.cmplx)

            if self.material_ey.has_key(type(pw_obj)):
                self.material_ey[type(pw_obj)].merge(pw_obj)
            else:
                self.material_ey[type(pw_obj)] = pw_obj

    def init_material_ez(self):
        """Set up the update mechanism for Ez field.
        
        Set up the update mechanism for Ez field and stores the result
        at self.material_ez.
        
        """
        shape = self.ez.shape
        for idx in ndindex(shape):
            spc = self.space.ez_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[0] == shape[0] - 1 or idx[1] == shape[1] - 1:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = \
                mat_obj.get_pw_material_ez(idx, spc, underneath, self.cmplx)

            if self.material_ez.has_key(type(pw_obj)):
                self.material_ez[type(pw_obj)].merge(pw_obj)
            else:
                self.material_ez[type(pw_obj)] = pw_obj

    def init_material_hx(self):
        """Set up the update mechanism for Hx field.
        
        Set up the update mechanism for Hx field and stores the result
        at self.material_hx.
        
        """
        shape = self.hx.shape
        for idx in ndindex(shape):
            spc = self.space.hx_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[1] == 0 or idx[2] == 0:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = \
                mat_obj.get_pw_material_hx(idx, spc, underneath, self.cmplx)

            if self.material_hx.has_key(type(pw_obj)):
                self.material_hx[type(pw_obj)].merge(pw_obj)
            else:
                self.material_hx[type(pw_obj)] = pw_obj

    def init_material_hy(self):
        """Set up the update mechanism for Hy field.
        
        Set up the update mechanism for Hy field and stores the result
        at self.material_hy.
        
        """
        shape = self.hy.shape
        for idx in ndindex(shape):
            spc = self.space.hy_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[2] == 0 or idx[0] == 0:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = \
                mat_obj.get_pw_material_hy(idx, spc, underneath, self.cmplx)

            if self.material_hy.has_key(type(pw_obj)):
                self.material_hy[type(pw_obj)].merge(pw_obj)
            else:
                self.material_hy[type(pw_obj)] = pw_obj

    def init_material_hz(self):
        """Set up the update mechanism for Hz field.
        
        Set up the update mechanism for Hz field and stores the result
        at self.material_hz.
        
        """
        shape = self.hz.shape
        for idx in ndindex(shape):
            spc = self.space.hz_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(spc)
            if idx[0] == 0 or idx[1] == 0:
                mat_obj = Dummy(mat_obj.eps_inf, mat_obj.mu_inf)
            pw_obj = \
                mat_obj.get_pw_material_hz(idx, spc, underneath, self.cmplx)

            if self.material_hz.has_key(type(pw_obj)):
                self.material_hz[type(pw_obj)].merge(pw_obj)
            else:
                self.material_hz[type(pw_obj)] = pw_obj

    def init_material(self):
        init_mat_func = {const.Ex: self.init_material_ex,
                         const.Ey: self.init_material_ey,
                         const.Ez: self.init_material_ez,
                         const.Hx: self.init_material_hx,
                         const.Hy: self.init_material_hy,
                         const.Hz: self.init_material_hz}
        
        for comp in self.e_field_compnt:
            init_mat_func[comp]()

        for comp in self.h_field_compnt:
            init_mat_func[comp]()

    def init_source_ex(self):
        for so in self.src_list:
            pw_src = so.get_pw_source_ex(self.ex, self.space)

            if pw_src is None:
                continue

            if self.source_ex.has_key(type(pw_src)):
                self.source_ex[type(pw_src)].merge(pw_src)
            else:
                self.source_ex[type(pw_src)] = pw_src
            
    def init_source_ey(self):
        for so in self.src_list:
            pw_src = so.get_pw_source_ey(self.ey, self.space)
            
            if pw_src is None:
                continue

            if self.source_ey.has_key(type(pw_src)):
                self.source_ey[type(pw_src)].merge(pw_src)
            else:
                self.source_ey[type(pw_src)] = pw_src

    def init_source_ez(self):
        for so in self.src_list:
            pw_src = so.get_pw_source_ez(self.ez, self.space)
            
            if pw_src is None:
                continue

            if self.source_ez.has_key(type(pw_src)):
                self.source_ez[type(pw_src)].merge(pw_src)
            else:
                self.source_ez[type(pw_src)] = pw_src

    def init_source_hx(self):
        for so in self.src_list:
            pw_src = so.get_pw_source_hx(self.hx, self.space)

            if pw_src is None:
                continue

            if self.source_hx.has_key(type(pw_src)):
                self.source_hx[type(pw_src)].merge(pw_src)
            else:
                self.source_hx[type(pw_src)] = pw_src
            
    def init_source_hy(self):
        for so in self.src_list:
            pw_src = so.get_pw_source_hy(self.hy, self.space)

            if pw_src is None:
                continue

            if self.source_hy.has_key(type(pw_src)):
                self.source_hy[type(pw_src)].merge(pw_src)
            else:
                self.source_hy[type(pw_src)] = pw_src
            
    def init_source_hz(self):
        for so in self.src_list:
            pw_src = so.get_pw_source_hz(self.hz, self.space)

            if pw_src is None:
                continue

            if self.source_hz.has_key(type(pw_src)):
                self.source_hz[type(pw_src)].merge(pw_src)
            else:
                self.source_hz[type(pw_src)] = pw_src
            
    def init_source(self):
        init_src_func = {const.Ex: self.init_source_ex,
                         const.Ey: self.init_source_ey,
                         const.Ez: self.init_source_ez,
                         const.Hx: self.init_source_hx,
                         const.Hy: self.init_source_hy,
                         const.Hz: self.init_source_hz}

        for comp in self.e_field_compnt:
            init_src_func[comp]()
            
        for comp in self.h_field_compnt:
            init_src_func[comp]()
        
    def set_probe(self, p, prefix):
        """
        p: space coordinates. type: tuple-3
        prefix: prefix of the recording file name. type: str
        
        """
        spc2idx = {const.Ex: self.space.space_to_ex_index,
                   const.Ey: self.space.space_to_ey_index,
                   const.Ez: self.space.space_to_ez_index,
                   const.Hx: self.space.space_to_hx_index,
                   const.Hy: self.space.space_to_hy_index,
                   const.Hz: self.space.space_to_hz_index}
        
        idx2spc = {const.Ex: self.space.ex_index_to_space,
                   const.Ey: self.space.ey_index_to_space,
                   const.Ez: self.space.ez_index_to_space,
                   const.Hx: self.space.hx_index_to_space,
                   const.Hy: self.space.hy_index_to_space,
                   const.Hz: self.space.hz_index_to_space}

        validity = {const.Ex: 
                    (lambda idx: in_range(idx, self.ex.shape, const.Ex)),
                    const.Ey: 
                    (lambda idx: in_range(idx, self.ey.shape, const.Ey)),
                    const.Ez: 
                    (lambda idx: in_range(idx, self.ez.shape, const.Ez)),
                    const.Hx:
                    (lambda idx: in_range(idx, self.hx.shape, const.Hx)),
                    const.Hy: 
                    (lambda idx: in_range(idx, self.hy.shape, const.Hy)),
                    const.Hz: 
                    (lambda idx: in_range(idx, self.hz.shape, const.Hz))}

        field = {const.Ex: self.ex, const.Ey: self.ey, const.Ez: self.ez,
                 const.Hx: self.hy, const.Hy: self.hy, const.Hz: self.hz}

        postfix = {const.Ex: '_ex.dat', const.Ey: '_ey.dat',
                   const.Ez: '_ez.dat', const.Hx: '_hx.dat',
                   const.Hy: '_hy.dat', const.Hz: '_hz.dat'}
        
        for comp in self.e_field_compnt:
            idx = spc2idx[comp](*p)
            if validity[comp](idx):
                recorder = Probe(idx, field[comp], postfix[comp])
                loc = idx2spc[comp](*idx)
                recorder.write_header(loc, self.time_step.dt)
                self.e_recorder.append(recorder)

        for comp in self.h_field_compnt:
            idx = spc2idx[comp](*p)
            if validity[comp](idx):
                recorder = Probe(idx, field[comp], postfix[comp])
                loc = idx2spc[comp](*idx)
                recorder.write_header(loc, self.time_step.dt)
                self.h_recorder.append(recorder)

    def update_ex(self):
        for pw_obj in self.material_ex.itervalues():
            pw_obj.update_all(self.ex, self.hz, self.hy, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.source_ex.itervalues():
            pw_obj.update_all(self.ex, self.hz, self.hy, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)
        
    def update_ey(self):
        for pw_obj in self.material_ey.itervalues():
            pw_obj.update_all(self.ey, self.hx, self.hz, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)
		
        for pw_obj in self.source_ey.itervalues():
            pw_obj.update_all(self.ey, self.hx, self.hz, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)

    def update_ez(self):
        for pw_obj in self.material_ez.itervalues():
            pw_obj.update_all(self.ez, self.hy, self.hx, self.dx, self.dy,
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.source_ez.itervalues():
            pw_obj.update_all(self.ez, self.hy, self.hx, self.dx, self.dy,
                              self.time_step.dt, self.time_step.n)
        
    def update_hx(self):
        for pw_obj in self.material_hx.itervalues():
            pw_obj.update_all(self.hx, self.ez, self.ey, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.source_hx.itervalues():
            pw_obj.update_all(self.hx, self.ez, self.ey, self.dy, self.dz, 
                              self.time_step.dt, self.time_step.n)
		
    def update_hy(self):
        for pw_obj in self.material_hy.itervalues():
            pw_obj.update_all(self.hy, self.ex, self.ez, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.source_hy.itervalues():
            pw_obj.update_all(self.hy, self.ex, self.ez, self.dz, self.dx,
                              self.time_step.dt, self.time_step.n)
		
    def update_hz(self):
        for pw_obj in self.material_hz.itervalues():
            pw_obj.update_all(self.hz, self.ey, self.ex, self.dx, self.dy, 
                              self.time_step.dt, self.time_step.n)

        for pw_obj in self.source_hz.itervalues():
            pw_obj.update_all(self.hz, self.ey, self.ex, self.dx, self.dy, 
                              self.time_step.dt, self.time_step.n)

    def talk_with_ex_neighbors(self):
        """Synchronize ex data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # send ex field data to -y direction and receive from +y direction.
        src, dest = self.space.cart_comm.Shift(1, -1)
        
        if self.cmplx:
            dest_spc = self.space.ex_index_to_space(0, 
                                                    self.ex.shape[1] - 1, 0)[1]
        
            src_spc = self.space.ex_index_to_space(0, 0, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Ex.tag,
                                                    None, src, const.Ex.tag)
            
            phase_shift = exp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.ex[:, -1, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ex[:, 0, :], dest, const.Ex.tag,
                                      None, src, const.Ex.tag)
        
        # send ex field data to -z direction and receive from +z direction.
        src, dest = self.space.cart_comm.Shift(2, -1)
        
        if self.cmplx:
            dest_spc = self.space.ex_index_to_space(0, 0, 
                                                    self.ex.shape[2] - 1)[2]
        
            src_spc = self.space.ex_index_to_space(0, 0, 0)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Ex.tag,
                                                    None, src, const.Ex.tag)
        
            phase_shift = exp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.ex[:, :, -1] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ex[:, :, 0], dest, const.Ex.tag,
                                      None, src, const.Ex.tag)
        
    def talk_with_ey_neighbors(self):
        """Synchronize ey data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # send ey field data to -z direction and receive from +z direction.
        src, dest = self.space.cart_comm.Shift(2, -1)
        
        if self.cmplx:
            dest_spc = self.space.ey_index_to_space(0, 0, 
                                                    self.ey.shape[2] - 1)[2]
        
            src_spc = self.space.ey_index_to_space(0, 0, 0)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Ey.tag,
                                                    None, src, const.Ey.tag) 
        
            phase_shift = exp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 0
            
        self.ey[:, :, -1] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ey[:, :, 0], dest, const.Ey.tag,
                                      None, src, const.Ey.tag)
        
        # send ey field data to -x direction and receive from +x direction.
        src, dest = self.space.cart_comm.Shift(0, -1)

        if self.cmplx:
            dest_spc = self.space.ey_index_to_space(self.ey.shape[0] - 1, 
                                                    0, 0)[0]
        
            src_spc = self.space.ey_index_to_space(0, 0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Ey.tag,
                                                    None, src, const.Ey.tag)
            
            phase_shift = exp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 0
            
        self.ey[-1, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ey[0, :, :], dest, const.Ey.tag,
                                      None, src, const.Ey.tag)
        
    def talk_with_ez_neighbors(self):
        """Synchronize ez data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # send ez field data to -x direction and receive from +x direction.
        src, dest = self.space.cart_comm.Shift(0, -1)
        
        if self.cmplx:
            dest_spc = self.space.ez_index_to_space(self.ez.shape[0] - 1, 
                                                    0, 0)[0]
        
            src_spc = self.space.ez_index_to_space(0, 0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Ez.tag,
                                                    None, src, const.Ez.tag) 
            
            phase_shift = exp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.ez[-1, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ez[0, :, :], dest, const.Ez.tag,
                                      None, src, const.Ez.tag)
        
        # send ez field data to -y direction and receive from +y direction.
        src, dest = self.space.cart_comm.Shift(1, -1)

        if self.cmplx:
            dest_spc = self.space.ez_index_to_space(0, 
                                                    self.ez.shape[1] - 1, 0)[1]
        
            src_spc = self.space.ez_index_to_space(0, 0, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Ez.tag,
                                                    None, src, const.Ez.tag) 
            
            phase_shift = exp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 0

        self.ez[:, -1, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.ez[:, 0, :], dest, const.Ez.tag,
                                      None, src, const.Ez.tag)
        
    def talk_with_hx_neighbors(self):
        """Synchronize hx data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # send hx field data to +y direction and receive from -y direction.
        src, dest = self.space.cart_comm.Shift(1, 1)

        if self.cmplx:
            dest_spc = self.space.hx_index_to_space(0, 0, 0)[1]
        
            src_spc = self.space.hx_index_to_space(0, 
                                                   self.hx.shape[1] - 1, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Hx.tag,
                                                    None, src, const.Hx.tag)
        
            phase_shift = exp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.hx[:, 0, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hx[:, -1, :], dest, const.Hx.tag,
                                      None, src, const.Hx.tag)
            
        # send hx field data to +z direction and receive from -z direction.    
        src, dest = self.space.cart_comm.Shift(2, 1)
        
        if self.cmplx:
            dest_spc = self.space.hx_index_to_space(0, 0, 0)[2]
        
            src_spc = self.space.hx_index_to_space(0, 0, 
                                                   self.hx.shape[2] - 1)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Hx.tag,
                                                    None, src, const.Hx.tag)
        
            phase_shift = exp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.hx[:, :, 0] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hx[:, :, -1], dest, const.Hx.tag,
                                      None, src, const.Hx.tag)
            
    def talk_with_hy_neighbors(self):
        """Synchronize hy data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # send hy field data to +z direction and receive from -z direction.
        src, dest = self.space.cart_comm.Shift(2, 1)
        
        if self.cmplx:
            dest_spc = self.space.hy_index_to_space(0, 0, 0)[2]
        
            src_spc = self.space.hy_index_to_space(0, 0, 
                                                   self.hy.shape[2] - 1)[2]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Hy.tag,
                                                    None, src, const.Hy.tag)
        
            phase_shift = exp(1j * self.bloch[2] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.hy[:, :, 0] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hy[:, :, -1], dest, const.Hy.tag,
                                      None, src, const.Hy.tag)
            
        # send hy field data to +x direction and receive from -x direction.
        src, dest = self.space.cart_comm.Shift(0, 1)
        
        if self.cmplx:
            dest_spc = self.space.hy_index_to_space(0, 0, 0)[0]
        
            src_spc = self.space.hy_index_to_space(self.hy.shape[0] - 1, 
                                                   0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Hy.tag,
                                                    None, src, const.Hy.tag)
        
            phase_shift = exp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.hy[0, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hy[-1, :, :], dest, const.Hy.tag,
                                      None, src, const.Hy.tag)
            
    def talk_with_hz_neighbors(self):
        """Synchronize hz data.
        
        This method uses the object serialization interface of MPI4Python.
        
        """
        # send hz field data to +x direction and receive from -x direction.
        src, dest = self.space.cart_comm.Shift(0, 1)
        
        if self.cmplx:
            dest_spc = self.space.hz_index_to_space(0, 0, 0)[0]
        
            src_spc = self.space.hz_index_to_space(self.hz.shape[0] - 1, 
                                                   0, 0)[0]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Hz.tag,
                                                    None, src, const.Hz.tag)
        
            phase_shift = exp(1j * self.bloch[0] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.hz[0, :, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hz[-1, :, :], dest, const.Hz.tag,
                                      None, src, const.Hz.tag)
            
        # send hz field data to +y direction and receive from -y direction.
        src, dest = self.space.cart_comm.Shift(1, 1)
        
        if self.cmplx:
            dest_spc = self.space.hz_index_to_space(0, 0, 0)[1]
        
            src_spc = self.space.hz_index_to_space(0, 
                                                   self.hz.shape[1] - 1, 0)[1]
            src_spc = self.space.cart_comm.sendrecv(src_spc, dest, const.Hz.tag,
                                                    None, src, const.Hz.tag)
        
            phase_shift = exp(1j * self.bloch[1] * (dest_spc - src_spc))
        else:
            phase_shift = 0
        
        self.hz[:, 0, :] = phase_shift * \
        self.space.cart_comm.sendrecv(self.hz[:, -1, :], dest, const.Hz.tag,
                                      None, src, const.Hz.tag)
        
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
        spc_to_idx = {const.Ex: self.space.space_to_ex_index,
                      const.Ey: self.space.space_to_ey_index,
                      const.Ez: self.space.space_to_ez_index,
                      const.Hx: self.space.space_to_hx_index,
                      const.Hy: self.space.space_to_hy_index,
                      const.Hz: self.space.space_to_hz_index}

        field = {const.Ex: self.ex, const.Ey: self.ey,
                 const.Ez: self.ez, const.Hx: self.hx,
                 const.Hy: self.hy, const.Hz: self.hz}

        idx = spc_to_idx[component](*point)

        if in_range(idx, field[component].shape, component):
            hot_node = self.space.my_id
        else:
            hot_node = None
            
        hot_node = self.space.bcast(hot_node, hot_node)
            
        flag = True
        while flag:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t
            if self.space.my_id == hot_node and field[component][idx] != 0:
                flag = False
            flag = self.space.cart_comm.bcast(flag, hot_node)

    def step_until_n(self, n=0, modulus=inf):
        """Run self.step() until time step reaches n.

        """
        while self.time_step.n < n:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t

    def step_until_t(self, t=0, modulus=inf):
        """Run self.step() until time reaches t.

        """
        while self.time_step.t < t:
            self.step()
            if self.time_step.n % modulus == 0:
                print 'n:', self.time_step.n, 't:', self.time_step.t
    
    def show_line_ex(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the ex along the line.


        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Ex, start, end, vrange, interval, 
                            'Ex field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
		
    def show_line_ey(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the ey along the line.

        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Ey, start, end, vrange, interval, 
                            'Ey field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
		
    def show_line_ez(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the ez along the line.

        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Ez, start, end, vrange, interval, 
                            'Ez field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
		
    def show_line_hx(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the hx along the line.

        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Hx, start, end, vrange, interval, 
                            'Hx field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
		
    def show_line_hy(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the hy along the line.

        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Hy, start, end, vrange, interval, 
                            'Hy field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
		
    def show_line_hz(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the hz along the line.

        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Hz, start, end, vrange, interval, 
                            'Hz field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase

    def show_ex(self, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ex on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Ex, axis, cut, vrange, 
                             interval, 'Ex field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_ey(self, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ey on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Ey, axis, cut, vrange, 
                             interval, 'Ey field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_ez(self, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ez on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Ez, axis, cut, vrange, 
                             interval, 'Ez field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_hx(self, axis, cut, vrange=(-1, 1), msecs=2500):
        """Show the real value of the hx on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Hx, axis, cut, vrange, 
                             interval, 'Hx field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_hy(self, axis, cut, vrange=(-1, 1), msecs=2500):
        """Show the real value of the hy on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Hy, axis, cut, vrange, 
                             interval, 'Hy field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_hz(self, axis, cut, vrange=(-1, 1), msecs=2500):
        """Show the real value of the hz on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Hz, axis, cut, vrange, 
                             interval, 'Hz field', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_permittivity_ex(self, axis, cut, vrange=None):
        """Show permittivity for the ex on the plane.

        Keyword arguments:
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. A tuple of length two.
        
        """
        showcase = Snapshot(self, const.Ex, axis, cut, vrange, 
                            'Permittivity for Ex', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_permittivity_ey(self, axis, cut, vrange=None):
        """Show permittivity for the ey on the plane.

        Keyword arguments:
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. A tuple of length two.
        
        """
        showcase = Snapshot(self, const.Ey, axis, cut, vrange, 
                            'Permittivity for Ey', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
            
    def show_permittivity_ez(self, axis, cut, vrange=None):
        """Show permittivity for the ez on the plane.

        Keyword arguments:
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. A tuple of length two.
        
        """
        showcase = Snapshot(self, const.Ez, axis, cut, vrange, 
                            'Permittivity for Ez', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase

    def show_permeability_hx(self, axis, cut, vrange=None):
        """Show permeability for the hx on the plane.

        Keyword arguments:
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. A tuple of length two.
        
        """
        showcase = Snapshot(self, const.Hx, axis, cut, vrange, 
                            'Permittivity for Hx', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_permeability_hy(self, axis, cut, vrange=None):
        """Show permeability for the hy on the plane.

        Keyword arguments:
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. A tuple of length two.
        
        """
        showcase = Snapshot(self, const.Hy, axis, cut, vrange, 
                            'Permittivity for Hy', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
            
    def show_permeability_hz(self, axis, cut, vrange=None):
        """Show permeability for the hz on the plane.

        Keyword arguments:
        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constant.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis.
        vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Hz, axis, cut, vrange, 
                            'Permittivity for Hz', self._fig_id)
        self._fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def write_ex(self, low=None, high=None, prefix=None, postfix=None):
        if low is None:
            low_idx = (0, 0, 0)
        else:
            low_idx = self.space.space_to_ex_index(low)
            
        if low is None:
            high_idx = self.ex.shape
        else:
            high_idx = self.space.space_to_ex_index(high)
        
        high_idx = [i + 1 for i in high_idx]
        
        name = ''
        if prefix is not None:
            name = prefix + name
        if postfix is not None:
            name = name + postfix
            
        write_hdf5(self.ex, name, low_idx, high_idx)
    	
    def write_ey(self):
        pass
    	
    def write_ez(self):
        pass

    def write_hx(self):
        pass
    	
    def write_hy(self):
        pass
    	
    def write_hz(self):
        pass
        
    def snapshot_ex(self, axis, cut):
        if axis is const.X:
            cut_idx = self.space.space_to_index(cut, 0, 0)[0]
            data = self.ex[cut_idx, :, :]
        elif axis is const.Y:
            cut_idx = self.space.space_to_index(0, cut, 0)[1]
            data = self.ex[:, cut_idx, :]
        elif axis is const.Z:
            cut_idx = self.space.space_to_index(0, 0, cut)[2]
            data = self.ex[:, :, cut_idx]
        else:
            pass
        
        filename = 't=' + str(self.time_step[1] * space.dt)
        snapshot(data, filename, const.Ex)
        
    def snapshot_ey(self, axis=const.Z, cut=0, vrange=(-.1, .1), size=(400, 400)):
        pass
    
    def snapshot_ez(self, axis=const.Z, cut=0, vrange=(-.1, .1), size=(400, 400)):
        pass
    
    def snapshot_hx(self, axis=const.Z, cut=0, vrange=(-.1, .1), size=(400, 400)):
        pass
        
    def snapshot_hy(self, axis=const.Z, cut=0, vrange=(-.1, .1), size=(400, 400)):
        pass
        
    def snapshot_hz(self, axis=const.Z, cut=0, vrange=(-.1, .1), size=(400, 400)):
        pass
        

class TExFDTD(FDTD):
    """2-D fdtd which has transverse-electric mode with respect to x.
    
    Assume that the structure and incident wave are uniform in the x 
    direction. TExFDTD updates only Ey, Ez, and Hx field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ey, const.Ez)
        self.h_field_compnt = (const.Hx,)

        
class TEyFDTD(FDTD):
    """2-D FDTD which has transverse-electric mode with respect to y.
    
    Assume that the structure and incident wave are uniform in the y direction.
    TEyFDTD updates only Ez, Ex, and Hy field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ex, const.Ez)
        self.h_field_compnt = (const.Hy,)


class TEzFDTD(FDTD):
    """2-D FDTD which has transverse-electric mode with respect to z.

    Assume that the structure and incident wave are uniform in the z direction.
    TEzFDTD updates only Ex, Ey, and Hz field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ex, const.Ey)
        self.h_field_compnt = (const.Hz,)


class TMxFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to x.

    Assume that the structure and incident wave are uniform in the x direction.
    TMxFDTD updates only Hy, Hz, and Ex field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ex,)
        self.h_field_compnt = (const.Hy, const.Hz)

        
class TMyFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to y

    Assume that the structure and incident wave are uniform in the y direction.
    TMyFDTD updates only Hz, Hx, and Ey field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ey,)
        self.h_field_compnt = (const.Hz, const.Hx)


class TMzFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to z
    
    Assume that the structure and incident wave are uniform in the z direction.
    TMzFDTD updates only Hx, Hy, and Ez field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ez,)
        self.h_field_compnt = (const.Hx, const.Hy)


class TEMxFDTD(FDTD):
    """y-polarized and x-directed one dimensional fdtd class

    Assume that the structure and incident wave are uniform in transverse 
    direction. TEMxFDTD updates only Ey and Hz field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ey,)
        self.h_field_compnt = (const.Hz,)

        
class TEMyFDTD(FDTD):
    """z-polarized and y-directed one dimensional fdtd class

    Assume that the structure and incident wave are uniform in transverse 
    direction. TEMyFDTD updates only Ez and Hx field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ez,)
        self.h_field_compnt = (const.Hx,)

        
class TEMzFDTD(FDTD):
    """x-polarized and z-directed one dimensional fdtd class
    
    Assume that the structure and incident wave are uniform in transverse 
    direction. TEMzFDTD updates only Ex and Hy field components.
    
    """
    def _init_field_compnt(self):
        self.e_field_compnt = (const.Ex,)
        self.h_field_compnt = (const.Hy,)


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
