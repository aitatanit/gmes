#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This code is based on libctl 3.0.2.

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass

from copy import deepcopy

import numpy as np
from numpy import empty, zeros, inf, dot, array

# GMES modules
import constant as const
from pygeom import *


class AuxiCartComm(object):
    """Auxiliary MPI Cartesian communicator for the absence of MPI implementation.
    
    Make an instance with default parameters only.

    Attributes:
    dims -- size of mpi Cartesian communicator
    ndims -- dimensionality of this Cartesian topology
        
    """
    def __init__(self, dims=(1,1,1), periods=None, reorder=(0,0,0)):
        """Constructor. 
        
        Keyword arguments:
        dims -- dimensions of the communicator (default (1,1,1))
        periods -- type: 3-int tuple, default (1,1,1)
        reorder -- default (0,0,0)

        """
        self.rank = 0
        self.dim = 3
        if periods:
            cyclic = tuple(map(int, periods))
        else:
            cyclic = (0, 0, 0)
        self.topo = ((1,1,1), cyclic, (0,0,0))
        
    def Get_cart_rank(self, coords):
        return 0
    
    def Get_coords(self, rank):
        """Get local or remote grid coordinate.
        
        Keyword arguments:
        rank -- local rank
            
        """
        return 0, 0, 0
    
    def Get_dim(self):
        return self.dim
    
    def Get_topo(self):
        return self.topo
    
    def Get_size(self):
        return 1

    def Shift(self, direction, disp):
        """Get source/destination with specified shift.
        
        Keyword arguments:
        direction -- 0 <= Dimension to move < ndims
        displacement -- steps to take in that dimension
            
        """
        if self.topo[1][direction]:
            return 0, 0
        else:
            return -1, -1
    
    def sendrecv(self, sendbuf, dest=0, sendtag=0,
                 recvbuf=None, source=0, recvtag=0, status=None):
        """Mimic Sendrecv method.
        
        All arguments except message are ignored.
        
        """
        if dest == -1 or source == -1:
            return None
        else:
            return sendbuf
        
    def reduce(self, value, root=0, op=None):
        """Mimic reduce method.
        
        """
        return value
    
    def bcast(self, obj=None, root=0):
        """Mimic bcast method.

        """
        return obj


class Cartesian(object):
    """Define the calculation space with Cartesian coordinates.
    
    Attributes:
    half_size -- the half size of whole calculation volume
    res -- number of sections of one unit length
    dx, dy, dz -- the space differentials
    dt -- the time differential
    whole_field_size -- the total array size for the each component of 
        the electromagnetic field except the communication buffers
    my_id -- mpi rank of this node
    numprocs -- the number of mpi nodes
    cart_comm -- mpi Cartesian communicator
    my_cart_idx -- the coordinates of this node in mpi Cartesian communicator
    general_field_size -- the general array size for the each component of 
        the electromagnetic field except the communication buffers
    my_field_size -- the specific array size for the each component of
        the electromagnetic field of this node except the communication buffers
            
    """
    def __init__(self, size, resolution=15, parallel=False):
        """Constructor

        Keyword arguments:
        size -- a length three sequence consists of non-negative numbers
        resolution -- number of sections of one unit. scalar or 3-tuple
            (default 15)
        parallel -- whether space be divided into segments (default False)    

        """
        try:
            if len(resolution) == 3:
                self.res = array(resolution, np.double)
        except TypeError:
            self.res = array((resolution,)*3, np.double)

        self.dr = array(1 / self.res, np.double)

        self.half_size = 0.5 * array(size, np.double)
        
        for i, v in enumerate(self.dr):
            if self.half_size[i] == 0:
                self.half_size[i] = .5 * v
            
        # the size of the whole field arrays 
        self.whole_field_size = \
            array((2 * self.half_size * self.res).round(), np.int)
        
        self.my_id = 0
        self.numprocs = 1
        self.cart_comm = AuxiCartComm((1,1,1), (1,1,1))
        try:
            if parallel:
                from mpi4py import MPI
                self.my_id = MPI.COMM_WORLD.rank
                self.numprocs = MPI.COMM_WORLD.size
                self.cart_comm = MPI.COMM_WORLD.Create_cart(self.find_best_deploy(), (1, 1, 1))
        except ImportError:
            pass

        self.my_cart_idx = self.cart_comm.topo[2]
        
        # Usually the my_field_size is general_field_size,
        # except the last node in each dimension.
        self.general_field_size = \
            self.whole_field_size / self.cart_comm.topo[0]
        
        # my_field_size may be different than general_field_size at the last 
        # node in each dimension.
        self.my_field_size = self.get_my_field_size()
    
    def bcast(self, obj=None, root=None):
        """Same with the Broadcast but, it handles for unknown root among 
        the nodes.
        
        """
        size = self.cart_comm.Get_size()
        if size == 1:
            return obj

        if root is None:
            obj = self.cart_comm.recv(source = MPI.ANY_SOURCE)
        else:
            for dest in xrange(size):
                if dest != root:
                    self.cart_comm.send(obj, dest)
            
        return obj

    def get_my_field_size(self):
        """Return the field size of this node.
        
        This method depends on 
        self.general_field_size
        self.whole_field_size
        self.my_cart_idx
            
        """
        field_size = empty(3, np.int)
        dims = self.cart_comm.topo[0]

        for i in xrange(3):
            # At the last node of that dimension.
            if self.my_cart_idx[i] == dims[i] - 1:
                field_size[i] = \
                self.whole_field_size[i] - (self.my_cart_idx[i] * 
                                            self.general_field_size[i])
            else:
                field_size[i] = self.general_field_size[i]
        
        return field_size
    
    def find_best_deploy(self):
        """Return the minimum load deploy of the nodes.
        
        This method depends on
        self.numprocs
            
        """
        best_partition = ()
        min_load = inf

        factors = [i for i in xrange(1, self.numprocs + 1) if self.numprocs % i == 0]
        for l in factors: 
            for m in factors:
                if l * m > self.numprocs:
                    break
                n = self.numprocs / (l * m)
                if self.numprocs % n == 0:
                        tmp_load = self.load_metric(l, m, n)
                        if tmp_load < min_load:
                            best_partition = l, m, n
                            min_load = tmp_load

        return best_partition

    def load_metric(self, l, m, n):
        """Estimate the load on a node.
    
        Keyword arguments:
        l, m, n -- the number of node in each direction
        
        This method depends on
        self.whole_field_size
        
        """
        l, m, n = float(l), float(m), float(n)
        
        # network load ratio compared to CPU
        R = 1000
        
        cpu_load = (self.whole_field_size[0] * self.whole_field_size[1] * 
                    self.whole_field_size[2]) / (l * m * n)
        net_load = 4 * R * (self.whole_field_size[0] / l * self.whole_field_size[1] / m + 
                            self.whole_field_size[1] / m * self.whole_field_size[2] / n + 
                            self.whole_field_size[2] / n * self.whole_field_size[0] / l)
        
        return cpu_load + net_load

    def _get_em_field_storage(self, shape, cmplx):
        if cmplx:
            return zeros(shape, complex)
        else:
            return zeros(shape, np.double)

    def get_ex_storage(self, field_compnt, cmplx=False):
        """Return an initialized array for Ex field component.
        
        """
        if const.Ex in field_compnt:
            shape = (self.my_field_size[0], self.my_field_size[1] + 1,
                     self.my_field_size[2] + 1)
        else:
            shape = (1, 1, 1)
        
        return self._get_em_field_storage(shape, cmplx)
        
    def get_ey_storage(self, field_compnt, cmplx=False):
        """Return an initialized array for Ey field component.
        
        """
        if const.Ey in field_compnt:
            shape = (self.my_field_size[0] + 1, self.my_field_size[1],
                     self.my_field_size[2] + 1)
        else:
            shape = (1, 1, 1)
        
        return self._get_em_field_storage(shape, cmplx)

    def get_ez_storage(self, field_compnt, cmplx=False):
        """Return an initialized array for Ez field component.
        
        """
        if const.Ez in field_compnt:
            shape = (self.my_field_size[0] + 1, self.my_field_size[1] + 1,
                     self.my_field_size[2])
        else:
            shape = (1, 1, 1)

        return self._get_em_field_storage(shape, cmplx)
        
    def get_hx_storage(self, field_compnt, cmplx=False):
        """Return an initialized array for Hx field component.
        
        """
        if const.Hx in field_compnt:
            shape = (self.my_field_size[0], self.my_field_size[1] + 1,
                     self.my_field_size[2] + 1)
        else:
            shape = (1, 1, 1)

        return self._get_em_field_storage(shape, cmplx)
        
    def get_hy_storage(self, field_compnt, cmplx=False):  
        """Return an initialized array for Hy field component.
        
        """
        if const.Hy in field_compnt:
            shape = (self.my_field_size[0] + 1, self.my_field_size[1],
                     self.my_field_size[2] + 1)
        else:
            shape = (1, 1, 1)
        
        return self._get_em_field_storage(shape, cmplx)
        
    def get_hz_storage(self, field_compnt, cmplx=False):
        """Return an initialized array for Hz field component.
        
        """
        if const.Hz in field_compnt:
            shape = (self.my_field_size[0] + 1, self.my_field_size[1] + 1,
                     self.my_field_size[2])
        else:
            shape = (1, 1, 1)
        
        return self._get_em_field_storage(shape, cmplx)

    def ex_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Ex mesh point.
        
        Keyword arguments:
        i, j, k -- array index
        
        """
        idx = array((i, j, k), np.int)  
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        spc_0 = (global_idx[0] + .5) * self.dr[0] - self.half_size[0]
        spc_1 = global_idx[1] * self.dr[1] - self.half_size[1]
        spc_2 = global_idx[2] * self.dr[2] - self.half_size[2]
        
        return spc_0, spc_1, spc_2

    def spc_to_exact_ex_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Ex mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Keyword arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)

        global_idx = empty(3, np.double)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] - .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1]
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2]
        
        idx = empty(3, np.double)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                  
        return tuple(idx)
    
    def space_to_ex_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Ex mesh point of 
        the given (global) space coordinate. The return index could be 
        out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        spc = array((x,y,z), np.double)

        global_idx = empty(3, np.int)
        global_idx[0] = (spc[0] + self.half_size[0]) / self.dr[0]
        global_idx[1] = (spc[1] + self.half_size[1]) / self.dr[1] + .5
        global_idx[2] = (spc[2] + self.half_size[2]) / self.dr[2] + .5

        idx = empty(3, np.int)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                  
        return tuple(idx)
    
    def ey_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Ey mesh point.
        
        Keyword arguments:
        i, j, k -- array index
        
        """
        idx = array((i,j,k), np.int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = global_idx[0] * self.dr[0] - self.half_size[0]
        coords_1 = (global_idx[1] + .5) * self.dr[1] - self.half_size[1]
        coords_2 = global_idx[2] * self.dr[2] - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_ey_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Ey mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.double)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0]
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] - .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2]
    
        idx = empty(3, np.double)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                  
        return tuple(idx)
        
    def space_to_ey_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Ey mesh point of 
        the given (global) space coordinate. The return index could be 
        out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1]
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] + .5
    
        idx = empty(3, np.int)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                  
        return tuple(idx)
    
    def ez_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Ez mesh point.
        
        Keyword arguments:
        i, j, k -- array index
        
        """
        idx = array((i, j, k), np.int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = global_idx[0] * self.dr[0] - self.half_size[0]
        coords_1 = global_idx[1] * self.dr[1] - self.half_size[1]
        coords_2 = (global_idx[2] + .5) * self.dr[2] - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_ez_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Ez mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x, y, z), np.double)
            
        global_idx = empty(3, np.double)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0]
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1]
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] - .5
        
        idx = empty(3, np.double)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                    
        return tuple(idx)
        
    def space_to_ez_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Ez mesh point of 
        the given (global) space coordinate. The return index could be 
        out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2]
        
        idx = empty(3, np.int)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                    
        return tuple(idx)

    def hx_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Hx mesh point.
        
        Keyword arguments:
        i, j, k -- array index
        
        """
        idx = array((i, j, k), np.int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = global_idx[0] * self.dr[0] - self.half_size[0]
        coords_1 = (global_idx[1] - .5) * self.dr[1] - self.half_size[1]
        coords_2 = (global_idx[2] - .5) * self.dr[2] - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_hx_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Hx mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.double)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0]
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] + .5

        idx = global_idx - self.my_cart_idx * self.general_field_size
        if self.whole_field_size[0] == 1:
            idx[0] = 0
        if self.whole_field_size[1] == 1:
            idx[1] = 1
        if self.whole_field_size[2] == 1:
            idx[2] = 1
            
        return tuple(idx)
    
    def space_to_hx_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Hx mesh point of 
        the given (global) space coordinate. The return index could be 
        out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] + 1
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] + 1

        idx = global_idx - self.my_cart_idx * self.general_field_size
        if self.whole_field_size[0] == 1:
            idx[0] = 0
        if self.whole_field_size[1] == 1:
            idx[1] = 1
        if self.whole_field_size[2] == 1:
            idx[2] = 1
            
        return tuple(idx)

    def hy_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Hy mesh point.
        
        Keyword arguments:
        i, j, k -- array index
        
        """
        idx = array((i,j,k), np.int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = (global_idx[0] - .5) * self.dr[0] - self.half_size[0]
        coords_1 = global_idx[1] * self.dr[1] - self.half_size[1]
        coords_2 = (global_idx[2] - .5) * self.dr[2] - self.half_size[2]
        
        return coords_0, coords_1, coords_2
        
    def spc_to_exact_hy_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Hy mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.double)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1]
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] + .5

        idx = global_idx - self.my_cart_idx * self.general_field_size
        if self.whole_field_size[0] == 1:
            idx[0] = 1
        if self.whole_field_size[1] == 1:
            idx[1] = 0
        if self.whole_field_size[2] == 1:
            idx[2] = 1
        
        return tuple(idx)
            
    def space_to_hy_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Hy mesh point of 
        the given (global) space coordinate. The return index could be
        out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + 1
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] + 1

        idx = global_idx - self.my_cart_idx * self.general_field_size
        if self.whole_field_size[0] == 1:
            idx[0] = 1
        if self.whole_field_size[1] == 1:
            idx[1] = 0
        if self.whole_field_size[2] == 1:
            idx[2] = 1
        
        return tuple(idx)
    
    def hz_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Hz mesh point.
        
        Keyword arguments:
        i, j, k -- array index
        
        """
        idx = array((i,j,k), np.int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = (global_idx[0] - .5) * self.dr[0] - self.half_size[0]
        coords_1 = (global_idx[1] - .5) * self.dr[1] - self.half_size[1]
        coords_2 = global_idx[2] * self.dr[2] - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_hz_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Hz mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Keyword arguments:
        x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.double)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2]

        idx = global_idx - self.my_cart_idx * self.general_field_size
        if self.whole_field_size[0] == 1:
            idx[0] = 1
        if self.whole_field_size[1] == 1:
            idx[1] = 1
        if self.whole_field_size[2] == 1:
            idx[2] = 0
        
        return tuple(idx)
            
    def space_to_hz_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Hy mesh point of 
        the given (global) space coordinate. The return index could be 
        out-of-range.
        
        Keyword arguments:
         x, y, z -- (global) space coordinate
        
        """
        coords = array((x,y,z), np.double)
            
        global_idx = empty(3, np.int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dr[0] + 1
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dr[1] + 1
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dr[2] + .5

        idx = global_idx - self.my_cart_idx * self.general_field_size
        if self.whole_field_size[0] == 1:
            idx[0] = 1
        if self.whole_field_size[1] == 1:
            idx[1] = 1
        if self.whole_field_size[2] == 1:
            idx[2] = 0
        
        return tuple(idx)
    
    def display_info(self, indent=0):
        print " " * indent, "Cartesian space"
        
        print " " * indent, "MPI topology:",
        print self.my_id, 'of', self.numprocs

        print " " * indent,
        print "size:", 2 * self.half_size,
        print "resolution:", self.res
        
        print " " * indent,
        print "dx:", self.dr[0], "dy:", self.dr[1], "dz:", self.dr[2]
        
        print " " * indent,
        print "number of participating nodes:", self.numprocs


def in_range(idx, shape, component):
    """Perform bounds checking.
    
    
    Keyword arguments:
        idx -- index of an array
        shape -- shape of the array to be checked
        component -- specify field component
        
    """
    if component is const.Ex:
        if idx[0] < 0 or idx[0] >= shape[0]:
            return False
        if idx[1] < 0 or idx[1] >= shape[1] - 1:
            return False
        if idx[2] < 0 or idx[2] >= shape[2] - 1:
            return False
        
    elif component is const.Ey:
        if idx[0] < 0 or idx[0] >= shape[0] - 1:
            return False
        if idx[1] < 0 or idx[1] >= shape[1]:
            return False
        if idx[2] < 0 or idx[2] >= shape[2] - 1:
            return False
          
    elif component is const.Ez:
        if idx[0] < 0 or idx[0] >= shape[0] - 1:
            return False
        if idx[1] < 0 or idx[1] >= shape[1] - 1:
            return False
        if idx[2] < 0 or idx[2] >= shape[2]:
            return False
        
    elif component is const.Hx:
        if idx[0] < 0 or idx[0] >= shape[0]:
            return False
        if idx[1] <= 0 or idx[1] >= shape[1]:
            return False
        if idx[2] <= 0 or idx[2] >= shape[2]:
            return False
        
    elif component is const.Hy:
        if idx[0] <= 0 or idx[0] >= shape[0]:
            return False
        if idx[1] < 0 or idx[1] >= shape[1]:
            return False
        if idx[2] <= 0 or idx[2] >= shape[2]:
            return False
        
    elif component is const.Hz:
        if idx[0] <= 0 or idx[0] >= shape[0]:
            return False
        if idx[1] <= 0 or idx[1] >= shape[1]:
            return False
        if idx[2] < 0 or idx[2] >= shape[2]:
            return False
        
    else:
        raise ValueError
       
    return True


if __name__ == '__main__':
    from material import Dielectric
    
    geom_list = [DefaultMedium(material=Dielectric()),
                 Cone(0, (1, 0, 0), 1, 1, Dielectric(), (0, 0, 2)),
                 Cone(0, (1, 0, 0), 1, 1, Dielectric(), (0, 0, -2))]
    t = GeomBoxTree(geom_list)
    t.display_info()
    space = Cartesian(size=(5, 5, 5))
    ex = space.get_ex_storage()
    print "ex shape:", ex.shape
    print "ex:", space.ex_index_to_space(0, 0, 0)
    print "ex:", space.ex_index_to_space(74, 75, 75)
    print "ex:", space.space_to_ex_index(0, 0, 0)
    print "ey:", space.space_to_ey_index(0, 0, 0)
    print "ez:", space.space_to_ez_index(0, 0, 0)
    print "hx:", space.space_to_hx_index(0, 0, 0)
    print "hy:", space.space_to_hy_index(0, 0, 0)
    print "hz:", space.space_to_hz_index(0, 0, 0)
