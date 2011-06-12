#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')

try:
    from threading import Thread, Lock
except ImportError:
    stderr.write('No module named threading. Using dummy_threading instead.\n')
    from dummy_threading import Thread, Lock

from copy import deepcopy
from math import sqrt
from cmath import exp

import numpy as np
from numpy import ndindex, arange

# GMES modules
from geometry import GeomBoxTree, in_range, DefaultMedium
from file_io import Probe
#from file_io import write_hdf5, snapshot
from show import ShowLine, ShowPlane, Snapshot
from material import Dummy
import constants as const


class TimeStep(object):
    """Store the current time-step and time.
    
    Attributes:
        n -- current time-step
        t -- current time
        
    """
    def __init__(self, n=0.0, t=0.0):
        self.n = n
        self.t = t
        
        
class FDTD(object):
    """three dimensional finite-difference time-domain class
    
    Attributes:
        space -- geometry.Cartesian instance
        cmplx -- Boolean of whether field is complex. Determined by the 
            space.period.
        
    """
    def __init__(self, space=None, geom_list=None, src_list=None,
                 courant_ratio=.99, dt=None, bloch=None, verbose=True):
        """
        Argumetns:
        space -- an instance which represents the coordinate system.
        geom_list -- a list which represents the geometric structure.
        src_list -- a list of source instances.
        courant_ratio -- the ratio of dt to Courant stability bound
                         default: 0.99
        dt -- the time differential
              If None is given, dt is calculated using space differentials 
              and courant_ratio.
              default: None
        bloch -- Bloch wave vector.
        verbose --
        
        Attributes:
        lock_ex, lock_ey, lock_ez
        lock_hx, lock_hy, lock_hz
        space
        fig_id
        dx, dy, dz
        dt
        courant_ratio: the ratio of dt to Courant stability bound
        bloch: Bloch wave vector

        """
        self.verbose = verbose

        self.lock_ex, self.lock_ey = Lock(), Lock()
        self.lock_ez, self.lock_hx = Lock(), Lock()
        self.lock_hy, self.lock_hz = Lock(), Lock()
       	
        self.space = space
                
        self.fig_id = self.space.my_id
            
        self.dx, self.dy, self.dz = space.dx, space.dy, space.dz

        default_medium = (i for i in geom_list 
                            if isinstance(i, DefaultMedium)).next()
        eps = default_medium.material.epsilon
        mu = default_medium.material.mu
        dt_limit = self.dt_limit(space, eps, mu)

        if dt is None:
            self.courant_ratio = float(courant_ratio)
            self.dt = self.courant_ratio * dt_limit
        else:
            self.dt = float(dt)
            self.courant_ratio = self.dt / dt_limit

        # Some codes in geometry.py and source.py use space.dt.
        self.space.dt = self.dt

        self.time_step = TimeStep()
        
        if verbose:
            self.space.display_info()

        if verbose:
            print 'dt:', self.dt
            print 'courant ratio:', self.courant_ratio
            
        if verbose:
            print "Initializing the geometry list...",
            
        self.geom_list = deepcopy(geom_list)
        for geom_obj in self.geom_list:
            geom_obj.init(self.space)
            
        if verbose:
            print "done."
            
        if verbose:
            print "Generating geometric binary search tree...",
            
        self.geom_tree = GeomBoxTree(self.geom_list)

        if verbose:
            print "done."
            
        if verbose:
            print "The geometric tree follows..."
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
            eps = default_medium.material.epsilon
            mu = default_medium.material.mu
            c = 1 / sqrt(eps * mu)
            S = c * self.dt / ds
            ref_n = sqrt(eps)
            self.bloch = np.array(bloch, float)
            self.bloch = 2 * ref_n / ds \
                * np.arcsin(np.sin(self.bloch * S * ds / 2) / S)
            
        if verbose:
            print "Bloch wave vector is", self.bloch
            
        if verbose:
            print "Initializing source...",
            
        self.src_list = deepcopy(src_list)
        for so in self.src_list:
            so.init(self.geom_tree, self.space, self.cmplx)
            
        if verbose:
            print "done."
            
        if verbose:
            print "The source list information follows..."
            for so in self.src_list:
                so.display_info()
                
        if verbose:
            print "Allocating memory for the electric & magnetic fields...",
            
        # storage for the electric & magnetic field 
        self.ex = space.get_ex_storage(self.cmplx)
        self.ey = space.get_ey_storage(self.cmplx)
        self.ez = space.get_ez_storage(self.cmplx)
        self.hx = space.get_hx_storage(self.cmplx)
        self.hy = space.get_hy_storage(self.cmplx)
        self.hz = space.get_hz_storage(self.cmplx)
        
        if verbose:
            print "done."
            
        if verbose:
            print "ex field:", self.ex.dtype, self.ex.shape 
            print "ey field:", self.ey.dtype, self.ey.shape 
            print "ez field:", self.ez.dtype, self.ez.shape
            print "hx field:", self.hx.dtype, self.hx.shape
            print "hy field:", self.hy.dtype, self.hy.shape
            print "hz field:", self.hz.dtype, self.hz.shape
        
        if verbose:
            print "Mapping the pointwise material...",
            
        self.material_ex = self.material_ey = self.material_ez = None
        self.material_hx = self.material_hy = self.material_hz = None
        
        self.init_field_components()

        self.init_material()

        if verbose:
            print "done."

        # medium information for electric & magnetic fields
        if verbose:
            print "ex material:",
            if self.material_ex is None: print None
            else: print self.material_ex.dtype, self.material_ex.shape
                
            print "ey material:", 
            if self.material_ey is None: print None
            else: print self.material_ey.dtype, self.material_ey.shape
                
            print "ez material:", 
            if self.material_ez is None: print None
            else: print self.material_ez.dtype, self.material_ez.shape
            
            print "hx material:", 
            if self.material_hx is None: print None
            else: print self.material_hx.dtype, self.material_hx.shape
                
            print "hy material:", 
            if self.material_hy is None: print None
            else: print self.material_hy.dtype, self.material_hy.shape
                
            print "hz material:", 
            if self.material_hz is None: print None
            else: print self.material_hz.dtype, self.material_hz.shape
		
        if verbose:
            print "done."
            
        if verbose:
            print "Mapping the pointwise source...",
            
        self.init_source()

        if verbose:
            print "done."
    
    def init_field_components(self):
        self.e_field_components = (const.Ex, const.Ey, const.Ez)
        self.h_field_components = (const.Hx, const.Hy, const.Hz)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.

        Equation 4.54b at p.131 of A. Taflove and S. C. Hagness, Computational
        Electrodynamics: The Finite-Difference Time-Domain Method, Third 
        Edition, 3rd ed. Artech House Publishers, 2005.

        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / c / sqrt(space.dx**-2 + space.dy**-2 + space.dz**-2)
        
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
                                 dt=self.dt, 
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
        self.lock_ex.acquire()
        
        self.material_ex = self.space.get_material_ex_storage()
        shape = self.ex.shape
        for idx in ndindex(shape):
            coords = self.space.ex_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(coords)
            if idx[1] == shape[1] - 1 or idx[2] == shape[2] - 1:
                mat_obj = Dummy(mat_obj.epsilon, mat_obj.mu)
            self.material_ex[idx] = \
            mat_obj.get_pw_material_ex(idx, coords, underneath, self.cmplx)
        
        self.lock_ex.release()
        
    def init_material_ey(self):
        """Set up the update mechanism for Ey field.
        
        Set up the update mechanism for Ey field and stores the result
        at self.material_ey.
        
        """
        self.lock_ey.acquire()
        
        self.material_ey = self.space.get_material_ey_storage()
        shape = self.ey.shape
        for idx in ndindex(shape):
            coords = self.space.ey_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(coords)
            if idx[2] == shape[2] - 1 or idx[0] == shape[0] - 1:
                mat_obj = Dummy(mat_obj.epsilon, mat_obj.mu)
            self.material_ey[idx] = \
            mat_obj.get_pw_material_ey(idx, coords, underneath, self.cmplx)
            
        self.lock_ey.release()
        
    def init_material_ez(self):
        """Set up the update mechanism for Ez field.
        
        Set up the update mechanism for Ez field and stores the result
        at self.material_ez.
        
        """
        self.lock_ez.acquire()
        
        self.material_ez = self.space.get_material_ez_storage()
        shape = self.ez.shape
        for idx in ndindex(shape):
            coords = self.space.ez_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(coords)
            if idx[0] == shape[0] - 1 or idx[1] == shape[1] - 1:
                mat_obj = Dummy(mat_obj.epsilon, mat_obj.mu)
            self.material_ez[idx] = \
            mat_obj.get_pw_material_ez(idx, coords, underneath, self.cmplx)
            
        self.lock_ez.release()
        
    def init_material_hx(self):
        """Set up the update mechanism for Hx field.
        
        Set up the update mechanism for Hx field and stores the result
        at self.material_hx.
        
        """
        self.lock_hx.acquire()
        
        self.material_hx = self.space.get_material_hx_storage()
        shape = self.hx.shape
        for idx in ndindex(shape):
            coords = self.space.hx_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(coords)
            if idx[1] == 0 or idx[2] == 0:
                mat_obj = Dummy(mat_obj.epsilon, mat_obj.mu)
            self.material_hx[idx] = \
            mat_obj.get_pw_material_hx(idx, coords, underneath, self.cmplx)
                
        self.lock_hx.release()
        
    def init_material_hy(self):
        """Set up the update mechanism for Hy field.
        
        Set up the update mechanism for Hy field and stores the result
        at self.material_hy.
        
        """
        self.lock_hy.acquire()
        
        self.material_hy = self.space.get_material_hy_storage()
        shape = self.hy.shape
        for idx in ndindex(shape):
            coords = self.space.hy_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(coords)
            if idx[2] == 0 or idx[0] == 0:
                mat_obj = Dummy(mat_obj.epsilon, mat_obj.mu)
            self.material_hy[idx] = \
            mat_obj.get_pw_material_hy(idx, coords, underneath, self.cmplx)
            
        self.lock_hy.release()
        
    def init_material_hz(self):
        """Set up the update mechanism for Hz field.
        
        Set up the update mechanism for Hz field and stores the result
        at self.material_hz.
        
        """
        self.lock_hz.acquire()
        
        self.material_hz = self.space.get_material_hz_storage()
        shape = self.hz.shape
        for idx in ndindex(shape):
            coords = self.space.hz_index_to_space(*idx)
            mat_obj, underneath = self.geom_tree.material_of_point(coords)
            if idx[0] == 0 or idx[1] == 0:
                mat_obj = Dummy(mat_obj.epsilon, mat_obj.mu)
            self.material_hz[idx] = \
            mat_obj.get_pw_material_hz(idx, coords, underneath, self.cmplx)
            
        self.lock_hz.release()
        
    def init_material(self):
        init_mat_func = {const.Ex: self.init_material_ex,
                         const.Ey: self.init_material_ey,
                         const.Ez: self.init_material_ez,
                         const.Hx: self.init_material_hx,
                         const.Hy: self.init_material_hy,
                         const.Hz: self.init_material_hz}
        
        for comp in self.e_field_components:
            init_mat_func[comp]()

        for comp in self.h_field_components:
            init_mat_func[comp]()

    def init_source_ex(self):
        for so in self.src_list:
            so.set_pointwise_source_ex(self.material_ex, self.space)
            
    def init_source_ey(self):
        for so in self.src_list:
            so.set_pointwise_source_ey(self.material_ey, self.space)
			
    def init_source_ez(self):
        for so in self.src_list:
            so.set_pointwise_source_ez(self.material_ez, self.space)
			
    def init_source_hx(self):
        for so in self.src_list:
            so.set_pointwise_source_hx(self.material_hx, self.space)
			
    def init_source_hy(self):
        for so in self.src_list:
            so.set_pointwise_source_hy(self.material_hy, self.space)
			
    def init_source_hz(self):
        for so in self.src_list:
            so.set_pointwise_source_hz(self.material_hz, self.space)
			
    def init_source(self):
        init_src_func = {const.Ex: self.init_source_ex,
                         const.Ey: self.init_source_ey,
                         const.Ez: self.init_source_ez,
                         const.Hx: self.init_source_hx,
                         const.Hy: self.init_source_hy,
                         const.Hz: self.init_source_hz}

        for comp in self.e_field_components:
            init_src_func[comp]()
            
        for comp in self.h_field_components:
            init_src_func[comp]()
            
    def set_probe(self, x, y, z, prefix):
        if self.material_ex is not None:
            idx = self.space.space_to_ex_index(x, y, z)
            if in_range(idx, self.material_ex, const.Ex):
                self.material_ex[idx] = Probe(prefix + '_ex.dat', 
                                              self.material_ex[idx])
                loc = self.space.ex_index_to_space(*idx)
                self.material_ex[idx].f.write('# location=' + str(loc) + '\n')
                self.material_ex[idx].f.write('# dt=' + str(self.dt) + '\n')
        
        if self.material_ey is not None:
            idx = self.space.space_to_ey_index(x, y, z)
            if in_range(idx, self.material_ey, const.Ey):
                self.material_ey[idx] = Probe(prefix + '_ey.dat', 
                                              self.material_ey[idx])
                loc = self.space.ey_index_to_space(*idx)
                self.material_ey[idx].f.write('# location=' + str(loc) + '\n')
                self.material_ey[idx].f.write('# dt=' + str(self.dt) + '\n')
        
        if self.material_ez is not None:
            idx = self.space.space_to_ez_index(x, y, z)
            if in_range(idx, self.material_ez, const.Ez):
                self.material_ez[idx] = Probe(prefix + '_ez.dat', 
                                              self.material_ez[idx])
                loc = self.space.ez_index_to_space(*idx)
                self.material_ez[idx].f.write('# location=' + str(loc) + '\n')
                self.material_ez[idx].f.write('# dt=' + str(self.dt) + '\n')
        
        if self.material_hx is not None:
            idx = self.space.space_to_hx_index(x, y, z)
            if in_range(idx, self.material_hx, const.Hx):
                self.material_hx[idx] = Probe(prefix + '_hx.dat', 
                                              self.material_hx[idx])
                loc = self.space.hx_index_to_space(*idx)
                self.material_hx[idx].f.write('# location=' + str(loc) + '\n')
                self.material_hx[idx].f.write('# dt=' + str(self.dt) + '\n')
        
        if self.material_hy is not None:
            idx = self.space.space_to_hy_index(x, y, z)
            if in_range(idx, self.material_hy, const.Hy):
                self.material_hy[idx] = Probe(prefix + '_hy.dat', 
                                              self.material_hy[idx])
                loc = self.space.hy_index_to_space(*idx)
                self.material_hy[idx].f.write('# location=' + str(loc) + '\n')
                self.material_hy[idx].f.write('# dt=' + str(self.dt) + '\n')
        
        if self.material_hz is not None:
            idx = self.space.space_to_hz_index(x, y, z)
            if in_range(idx, self.material_hz, const.Hz):
                self.material_hz[idx] = Probe(prefix + '_hz.dat', 
                                              self.material_hz[idx])
                loc = self.space.hz_index_to_space(*idx)
                self.material_hz[idx].f.write('# location=' + str(loc) + '\n')
                self.material_hz[idx].f.write('# dt=' + str(self.dt) + '\n')
            
    def update_ex(self):
        self.lock_ex.acquire()
        for mo in self.material_ex.flat:
            mo.update(self.ex, self.hz, self.hy, self.dy, self.dz, self.dt, 
                      self.time_step.n)
        self.lock_ex.release()
        
    def update_ey(self):
        self.lock_ey.acquire()
        for mo in self.material_ey.flat:
            mo.update(self.ey, self.hx, self.hz, self.dz, self.dx, self.dt, 
                      self.time_step.n)
        self.lock_ey.release()
		
    def update_ez(self):
        self.lock_ez.acquire()
        for mo in self.material_ez.flat:
            mo.update(self.ez, self.hy, self.hx, self.dx, self.dy, self.dt, 
                      self.time_step.n)
        self.lock_ez.release()
		
    def update_hx(self):
        self.lock_hx.acquire()
        for mo in self.material_hx.flat:
            mo.update(self.hx, self.ez, self.ey, self.dy, self.dz, self.dt, 
                      self.time_step.n)
        self.lock_hx.release()
		
    def update_hy(self):
        self.lock_hy.acquire()
        for mo in self.material_hy.flat:
            mo.update(self.hy, self.ex, self.ez, self.dz, self.dx, self.dt, 
                      self.time_step.n)
        self.lock_hy.release()
		
    def update_hz(self):
        self.lock_hz.acquire()
        for mo in self.material_hz.flat:
            mo.update(self.hz, self.ey, self.ex, self.dx, self.dy, self.dt, 
                      self.time_step.n)
        self.lock_hz.release()

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
        
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hx_neighbors),
                                  Thread(target=self.talk_with_hy_neighbors), 
                                  Thread(target=self.talk_with_hz_neighbors))

        self._e_chatter_threads = (Thread(target=self.talk_with_ex_neighbors),
                                  Thread(target=self.talk_with_ey_neighbors), 
                                  Thread(target=self.talk_with_ez_neighbors))

        self._h_worker_threads = (Thread(target=self.update_hx),
                                 Thread(target=self.update_hy),
                                 Thread(target=self.update_hz))

        self._e_worker_threads = (Thread(target=self.update_ex),
                                 Thread(target=self.update_ey),
                                 Thread(target=self.update_ez))

    def step(self):
        self._init_threads()

        self.time_step.n += .5
        self.time_step.t = self.time_step.n * self.dt
        
        for chatter in self._h_chatter_threads:
            chatter.start()
                   
        for chatter in self._h_chatter_threads:
            chatter.join()
        
        for worker in self._e_worker_threads:
            worker.start()
           
        for worker in self._e_worker_threads:
            worker.join()

        self.time_step.n += .5
        self.time_step.t = self.time_step.n * self.dt

        self._step_aux_fdtd()
        
        for chatter in self._e_chatter_threads:
            chatter.start()
                   
        for chatter in self._e_chatter_threads:
            chatter.join()
        
        for worker in self._h_worker_threads:
            worker.start()
           
        for worker in self._h_worker_threads:
            worker.join()
        
    def step_while_zero(self, component, point):
        """Run self.step() while the field value at the given point is 0.
        
        Arguments:
            component: filed component to check.
                one of gmes.constants.{Ex, Ey, Ez, Hx, Hy, Hz}.
            point: coordinates of the location to check the field. 
                tuple of three scalars

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

        if in_range(idx, field[component], component):
            hot_node = self.space.my_id
        else:
            hot_node = None
            
        hot_node = self.space.bcast(hot_node, hot_node)
            
        flag = True
        while flag:
            self.step()
            if self.space.my_id == hot_node and field[component][idx] != 0:
                flag = False
            flag = self.space.cart_comm.bcast(flag, hot_node)

    def step_until_n(self, n=0):
        """Run self.step() untill time step reaches n.

        """
        while self.time_step.n < n:
            self.step()

    def step_untill_t(self, t=0):
        """Run self.step() untill time reaches t.

        """
        while self.time_step.t < t:
            self.step()
    
    def show_line_ex(self, start, end, vrange=(-1, 1), interval=2500):
        """Show the real value of the ex along the line.


        start: The start point of the probing line.
        end: The end point of the probing line.
        vrange: Plot range of the y axis.
        interval: Refresh rate of the plot in milliseconds.

        """
        showcase = ShowLine(self, const.Ex, start, end, vrange, interval, 
                            'Ex field', self.fig_id)
        self.fig_id += self.space.numprocs
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
                            'Ey field', self.fig_id)
        self.fig_id += self.space.numprocs
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
                            'Ez field', self.fig_id)
        self.fig_id += self.space.numprocs
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
                            'Hx field', self.fig_id)
        self.fig_id += self.space.numprocs
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
                            'Hy field', self.fig_id)
        self.fig_id += self.space.numprocs
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
                            'Hz field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase

    def show_ex(self, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ex on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constants.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Ex, axis, cut, vrange, 
                             interval, 'Ex field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_ey(self, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ey on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constants.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Ey, axis, cut, vrange, 
                             interval, 'Ey field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_ez(self, axis, cut, vrange=(-1, 1), interval=2500):
        """Show the real value of the ez on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constants.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Ez, axis, cut, vrange, 
                             interval, 'Ez field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_hx(self, axis, cut, vrange=(-1, 1), msecs=2500):
        """Show the real value of the hx on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constants.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Hx, axis, cut, vrange, 
                             interval, 'Hx field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_hy(self, axis, cut, vrange=(-1, 1), msecs=2500):
        """Show the real value of the hy on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constants.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Hy, axis, cut, vrange, 
                             interval, 'Hy field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_hz(self, axis, cut, vrange=(-1, 1), msecs=2500):
        """Show the real value of the hz on the plone.

        axis: Specify the normal axis to the show plane.
            This should be one of the gmes.constants.Directional.
        cut: A scalar value which specifies the cut position on the 
            axis. 
        vrange: Specify the colorbar range.
        inerval: Refresh rates in millisecond.

        """
        showcase = ShowPlane(self, const.Hz, axis, cut, vrange, 
                             interval, 'Hz field', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_permittivity_ex(self, axis, cut, vrange=None):
        """Show permittivity for the ex on the plane.

        Arguments:
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis.
            vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Ex, axis, cut, vrange, 
                             'Permittivity for Ex', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_permittivity_ey(self, axis, cut, vrange=None):
        """Show permittivity for the ey on the plane.

        Arguments:
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis.
            vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Ey, axis, cut, vrange, 
                             'Permittivity for Ey', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
            
    def show_permittivity_ez(self, axis, cut, vrange=None):
        """Show permittivity for the ez on the plane.

        Arguments:
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis.
            vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Ez, axis, cut, vrange, 
                             'Permittivity for Ez', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase

    def show_permeability_hx(self, axis, cut, vrange=None):
        """Show permeability for the hx on the plane.

        Arguments:
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis.
            vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Hx, axis, cut, vrange, 
                             'Permittivity for Hx', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
        
    def show_permeability_hy(self, axis, cut, vrange=None):
        """Show permeability for the hy on the plane.

        Arguments:
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis.
            vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Hy, axis, cut, vrange, 
                             'Permittivity for Hy', self.fig_id)
        self.fig_id += self.space.numprocs
        showcase.start()
        return showcase
            
    def show_permeability_hz(self, axis, cut, vrange=None):
        """Show permeability for the hz on the plane.

        Arguments:
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis.
            vrange: Specify the colorbar range. a tuple of length two.
        
        """
        showcase = Snapshot(self, const.Hz, axis, cut, vrange, 
                             'Permittivity for Hz', self.fig_id)
        self.fig_id += self.space.numprocs
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
    def init_field_components(self):
        self.e_field_components = (const.Ey, const.Ez)
        self.h_field_components = (const.Hx,)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.

        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / c / sqrt(space.dy**-2 + space.dz**-2)
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hx_neighbors),)

        self._e_chatter_threads = (Thread(target=self.talk_with_ey_neighbors), 
                                  Thread(target=self.talk_with_ez_neighbors))

        self._h_worker_threads = (Thread(target=self.update_hx),)

        self._e_worker_threads = (Thread(target=self.update_ey),
                                 Thread(target=self.update_ez))

        
class TEyFDTD(FDTD):
    """2-D FDTD which has transverse-electric mode with respect to y.
    
    Assume that the structure and incident wave are uniform in the y direction.
    TEyFDTD updates only Ez, Ex, and Hy field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ex, const.Ez)
        self.h_field_components = (const.Hy,)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.
        
        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / c /sqrt(space.dx**-2 + space.dz**-2)
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hy_neighbors),)

        self._e_chatter_threads = (Thread(target=self.talk_with_ex_neighbors),
                                  Thread(target=self.talk_with_ez_neighbors))

        self._h_worker_threads = (Thread(target=self.update_hy),)

        self._e_worker_threads = (Thread(target=self.update_ex),
                                 Thread(target=self.update_ez))


class TEzFDTD(FDTD):
    """2-D FDTD which has transverse-electric mode with respect to z.

    Assume that the structure and incident wave are uniform in the z direction.
    TEzFDTD updates only Ex, Ey, and Hz field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ex, const.Ey)
        self.h_field_components = (const.Hz,)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.
        
        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / c / sqrt(space.dx**-2 + space.dy**-2)
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hz_neighbors),)

        self._e_chatter_threads = (Thread(target=self.talk_with_ex_neighbors),
                                  Thread(target=self.talk_with_ey_neighbors))

        self._h_worker_threads = (Thread(target=self.update_hz),)

        self._e_worker_threads = (Thread(target=self.update_ex),
                                 Thread(target=self.update_ey))

        
class TMxFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to x.

    Assume that the structure and incident wave are uniform in the x direction.
    TMxFDTD updates only Hy, Hz, and Ex field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ex,)
        self.h_field_components = (const.Hy, const.Hz)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.

        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / c / sqrt(space.dy**-2 + space.dz**-2)
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hy_neighbors), 
                                  Thread(target=self.talk_with_hz_neighbors))

        self._e_chatter_threads = (Thread(target=self.talk_with_ex_neighbors),)

        self._h_worker_threads = (Thread(target=self.update_hy),
                                 Thread(target=self.update_hz))

        self._e_worker_threads = (Thread(target=self.update_ex),)

        
class TMyFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to y

    Assume that the structure and incident wave are uniform in the y direction.
    TMyFDTD updates only Hz, Hx, and Ey field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ey,)
        self.h_field_components = (const.Hz, const.Hx)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.

        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / sqrt(space.dx**-2 + space.dz**-2)
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hx_neighbors),
                                  Thread(target=self.talk_with_hz_neighbors))

        self._e_chatter_threads = (Thread(target=self.talk_with_ey_neighbors),)

        self._h_worker_threads = (Thread(target=self.update_hx),
                                 Thread(target=self.update_hz))

        self._e_worker_threads = (Thread(target=self.update_ey),)

        
class TMzFDTD(FDTD):
    """2-D FDTD which has transverse-magnetic mode with respect to z
    
    Assume that the structure and incident wave are uniform in the z direction.
    TMzFDTD updates only Hx, Hy, and Ez field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ez,)
        self.h_field_components = (const.Hx, const.Hy)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.

        """
        c = 1 / sqrt(epsilon * mu)
        return 1 / c / sqrt(space.dx**-2 + space.dy**-2)
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hx_neighbors),
                                  Thread(target=self.talk_with_hy_neighbors))

        self._e_chatter_threads = (Thread(target=self.talk_with_ez_neighbors),)

        self._h_worker_threads = (Thread(target=self.update_hx),
                                 Thread(target=self.update_hy))

        self._e_worker_threads = (Thread(target=self.update_ez),)


class TEMxFDTD(FDTD):
    """y-polarized and x-directed one dimensional fdtd class

    Assume that the structure and incident wave are uniform in transverse direction.
    TEMxFDTD updates only Ey and Hz field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ey,)
        self.h_field_components = (const.Hz,)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.
        
        """
        c = 1 / sqrt(epsilon * mu)
        return space.dx / c
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hz_neighbors),)

        self._e_chatter_threads = (Thread(target=self.talk_with_ey_neighbors),)

        self._h_worker_threads = (Thread(target=self.update_hz),)

        self._e_worker_threads = (Thread(target=self.update_ey),)

        
class TEMyFDTD(FDTD):
    """z-polarized and y-directed one dimensional fdtd class

    Assume that the structure and incident wave are uniform in transverse direction.
    TEMyFDTD updates only Ez and Hx field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ez,)
        self.h_field_components = (const.Hx,)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.
        
        """
        c = 1 / sqrt(epsilon * mu)
        return space.dy / c
    
    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hx_neighbors),)

        self._e_chatter_threads = (Thread(target=self.talk_with_ez_neighbors),)

        self._h_worker_threads = (Thread(target=self.update_hx),)

        self._e_worker_threads = (Thread(target=self.update_ez),)

        
class TEMzFDTD(FDTD):
    """x-polarized and z-directed one dimensional fdtd class
    
    Assume that the structure and incident wave are uniform in transverse 
    direction. TEMzFDTD updates only Ex and Hy field components.
    
    """
    def init_field_components(self):
        self.e_field_components = (const.Ex,)
        self.h_field_components = (const.Hy,)

    def dt_limit(self, space, epsilon, mu):
        """Courant stability bound of a time step.
        
        """
        c = 1 / sqrt(epsilon * mu)
        return space.dz / c

    def _init_threads(self):
        self._h_chatter_threads = (Thread(target=self.talk_with_hy_neighbors),)

        self._e_chatter_threads = (Thread(target=self.talk_with_ex_neighbors),)

        self._h_worker_threads = (Thread(target=self.update_hy),)

        self._e_worker_threads = (Thread(target=self.update_ex),)


if __name__ == '__main__':
    from math import sin
    
    from numpy import inf
    
    from geometry import Cylinder, Cartesian
    from material import Dielectric
    
    low = Dielectric(index=1)
    hi = Dielectric(index=3)
    width_hi = low.epsilon / (low.epsilon + hi.epsilon)
    space = Cartesian(size=[1, 1, 1])
    geom_list = [DefaultMedium(material=low), 
                 Cylinder(material=hi, 
                          axis=[1, 0, 0], 
                          radius=inf, 
                          height=width_hi)]
    
    a = FDTD(space=space, geometry=geom_list)
    
    while True:
        a.step()
        a.ex[7, 7, 7] = sin(a.n)
        print a.n
