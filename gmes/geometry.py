#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

# this code is based on libctl 3.0.2.

from copy import deepcopy

import numpy as np
from numpy import empty, zeros, inf, dot
from scipy.linalg import norm

try:
    from mpi4py import MPI
except ImportError:
    pass

import constants as const
from material import Compound, PML


class AuxiCartComm(object):
    """Auxiliary MPI Cartesian communicator for the absence of MPI implementation.
    
    Attributes:
        dims -- size of mpi Cartesian communicator
        ndims -- dimensionality of this Cartesian topology
        
    """
    def __init__(self, dims=(1,1,1), periods=(1,1,1), reorder=(0,0,0)):
        self.rank = 0
        self.dim = 3
        cyclic = tuple(map(int, periods))
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
        return sendbuf
        
    def Reduce(self, value, root=0, op=None):
        """Mimic Reduce method.
        
        """
        return value
    
    
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
        """
        Arguments:
            size -- a length three sequence consists of non-negative numbers
            resolution -- number of sections of one unit (scalar or length 3 
                sequence)
                default: 15
            parallel -- whether space be divided into segments.
            
        """
        try:
            if len(resolution) == 3:
                self.res = np.array(resolution, float)
        except TypeError:
            self.res = np.array((resolution,)*3, float)
        
        self.dx, self.dy, self.dz = 1 / self.res

        # local spatial differentials to calculate dt
        dx, dy, dz = self.dx, self.dy, self.dz
        
        self.half_size = 0.5 * np.array(size, float)
        
        if self.half_size[0] == 0:
            self.half_size[0] = .5 * self.dx
            
        if self.half_size[1] == 0:
            self.half_size[1] = .5 * self.dy
            
        if self.half_size[2] == 0:
            self.half_size[2] = .5 * self.dz
            
        # the size of the whole field arrays 
        self.whole_field_size = np.array((2 * self.half_size * self.res).round(), int)
        
        try:
            if not parallel:
                raise StandardError
            
            self.my_id = MPI.COMM_WORLD.rank
            self.numprocs = MPI.COMM_WORLD.size
            self.cart_comm = \
            MPI.COMM_WORLD.Create_cart(self.find_best_deploy(), (1, 1, 1))
        except StandardError:
            self.my_id = 0
            self.numprocs = 1
            self.cart_comm = AuxiCartComm()
            
        self.my_cart_idx = self.cart_comm.topo[2]
        
        # usually the my_field_size is general_field_size,
        # except the last node in each dimension.
        self.general_field_size = \
        self.whole_field_size / self.cart_comm.topo[0]
        
        # my_field_size may be different than general_field_size at the last 
        # node in each dimension.
        self.my_field_size = self.get_my_field_size()
        
    def get_my_field_size(self):
        """Return the field size of this node.
        
        This method depends on 
            self.general_field_size
            self.whole_field_size
            self.my_cart_idx
            
        """
        field_size = empty(3, int)
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
    
        Arguments:
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
        
    def get_ex_storage(self, cmplx=False):
        """Return an initialized array for Ex field component.
        
        """
        ex_shape = (self.my_field_size[0], self.my_field_size[1] + 1,
                    self.my_field_size[2] + 1)
        
        if cmplx:
            return zeros(ex_shape, complex)
        else:
            return zeros(ex_shape, float)
        
    def get_ey_storage(self, cmplx=False):
        """Return an initialized array for Ey field component.
        
        """
        ey_shape = (self.my_field_size[0] + 1, self.my_field_size[1],
                    self.my_field_size[2] + 1)
        
        if cmplx:
            return zeros(ey_shape, complex)
        else:
            return zeros(ey_shape, float)
        
    def get_ez_storage(self, cmplx=False):
        """Return an initialized array for Ez field component.
        
        """
        ez_shape = (self.my_field_size[0] + 1, self.my_field_size[1] + 1,
                    self.my_field_size[2])
        
        if cmplx:
            return zeros(ez_shape, complex)
        else:
            return zeros(ez_shape, float)
        
    def get_hx_storage(self, cmplx=False):
        """Return an initialized array for Hx field component.
        
        """
        hx_shape = (self.my_field_size[0], self.my_field_size[1] + 1,
                    self.my_field_size[2] + 1)
        
        if cmplx:
            return zeros(hx_shape, complex)
        else:
            return zeros(hx_shape, float)
        
    def get_hy_storage(self, cmplx=False):  
        """Return an initialized array for Hy field component.
        
        """
        hy_shape = (self.my_field_size[0] + 1, self.my_field_size[1],
                    self.my_field_size[2] + 1)
        
        if cmplx:
            return zeros(hy_shape, complex)
        else:
            return zeros(hy_shape, float)
        
    def get_hz_storage(self, cmplx=False):
        """Return an initialized array for Hz field component.
        
        """
        hz_shape = (self.my_field_size[0] + 1, self.my_field_size[1] + 1,
                    self.my_field_size[2])
        
        if cmplx:
            return zeros(hz_shape, complex)
        else:
            return zeros(hz_shape, float)
        
    def get_material_ex_storage(self):
        """Return an array for Ex pointwise materials.
        
        """
        ex_shape = (self.my_field_size[0], self.my_field_size[1] + 1,
                    self.my_field_size[2] + 1)
        return empty(ex_shape, object)
        
    def get_material_ey_storage(self):
        """Return an array for Ey pointwise materials.
        
        """
        ey_shape = (self.my_field_size[0] + 1, self.my_field_size[1],
                    self.my_field_size[2] + 1)
        return empty(ey_shape, object)
        
    def get_material_ez_storage(self):
        """Return an array for Ez pointwise materials.
        
        """
        ez_shape = (self.my_field_size[0] + 1, self.my_field_size[1] + 1,
                    self.my_field_size[2])
        return empty(ez_shape, object)
    
    def get_material_hx_storage(self):
        """Return an array for Hx pointwise materials.
        
        """
        hx_shape = (self.my_field_size[0], self.my_field_size[1] + 1,
                    self.my_field_size[2] + 1)
        return empty(hx_shape, object)
        
    def get_material_hy_storage(self):
        """Return an array for Hy pointwise materials.
        
        """
        hy_shape = (self.my_field_size[0] + 1, self.my_field_size[1],
                    self.my_field_size[2] + 1)
        return empty(hy_shape, object)
        
    def get_material_hz_storage(self):
        """Return an array for Hz pointwise materials.
        
        """
        hz_shape = (self.my_field_size[0] + 1, self.my_field_size[1] + 1,
                    self.my_field_size[2])
        return empty(hz_shape, object)
    
    def ex_index_to_space(self, i, j, k):
        """Return space coordinate of the given index.
        
        This method returns the (global) space coordinates corresponding to 
        the given (local) index of Ex mesh point.
        
        Arguments:
            i, j, k -- array index
        
        """
        idx = np.array((i, j, k), int)  
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        spc_0 = (global_idx[0] + .5) * self.dx - self.half_size[0]
        spc_1 = global_idx[1] * self.dy - self.half_size[1]
        spc_2 = global_idx[2] * self.dz - self.half_size[2]
        
        return spc_0, spc_1, spc_2

    def spc_to_exact_ex_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Ex mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)

        global_idx = empty(3, float)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx - .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz
        
        idx = empty(3, float)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                  
        return tuple(idx)
    
    def space_to_ex_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Ex mesh point of the
        given (global) space coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        spc = np.array((x,y,z), float)

        global_idx = empty(3, int)
        global_idx[0] = (spc[0] + self.half_size[0]) / self.dx
        global_idx[1] = (spc[1] + self.half_size[1]) / self.dy + .5
        global_idx[2] = (spc[2] + self.half_size[2]) / self.dz + .5

        idx = empty(3, int)
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
        
        Arguments:
            i, j, k -- array index
        
        """
        idx = np.array((i,j,k), int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = global_idx[0] * self.dx - self.half_size[0]
        coords_1 = (global_idx[1] + .5) * self.dy - self.half_size[1]
        coords_2 = global_idx[2] * self.dz - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_ey_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Ey mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, float)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy - .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz
    
        idx = empty(3, float)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                  
        return tuple(idx)
        
    def space_to_ey_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Ey mesh point of the
        given (global) space coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz + .5
    
        idx = empty(3, int)
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
        
        Arguments
            i, j, k -- array index
        
        """
        idx = np.array((i, j, k), int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = global_idx[0] * self.dx - self.half_size[0]
        coords_1 = global_idx[1] * self.dy - self.half_size[1]
        coords_2 = (global_idx[2] + .5) * self.dz - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_ez_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Ez mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x, y, z), float)
            
        global_idx = empty(3, float)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz - .5
        
        idx = empty(3, float)
        for i in xrange(3):
            if self.whole_field_size[i] == 1:
                idx[i] = 0
            else:
                idx[i] = global_idx[i] - (self.my_cart_idx[i] * 
                                          self.general_field_size[i])
                    
        return tuple(idx)
        
    def space_to_ez_index(self, x, y, z):
        """Return the nearest mesh point of the given space coordinate.
        
        This method returns the (local) index of the nearest Ez mesh point of the
        given (global) space coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz
        
        idx = empty(3, int)
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
        
        Arguments:
            i, j, k -- array index
        
        """
        idx = np.array((i, j, k), int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = global_idx[0] * self.dx - self.half_size[0]
        coords_1 = (global_idx[1] - .5) * self.dy - self.half_size[1]
        coords_2 = (global_idx[2] - .5) * self.dz - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_hx_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Hx mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, float)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz + .5

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
        
        This method returns the (local) index of the nearest Hx mesh point of the
        given (global) space coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy + 1
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz + 1

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
        
        Arguments
            i, j, k -- array index
        
        """
        idx = np.array((i,j,k), int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = (global_idx[0] - .5) * self.dx - self.half_size[0]
        coords_1 = global_idx[1] * self.dy - self.half_size[1]
        coords_2 = (global_idx[2] - .5) * self.dz - self.half_size[2]
        
        return coords_0, coords_1, coords_2
        
    def spc_to_exact_hy_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Hy mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, float)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz + .5

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
        
        This method returns the (local) index of the nearest Hy mesh point of the
        given (global) space coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + 1
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz + 1

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
        
        Arguments:
            i, j, k -- array index
        
        """
        idx = np.array((i,j,k), int)
            
        global_idx = idx + self.general_field_size * self.my_cart_idx
        
        coords_0 = (global_idx[0] - .5) * self.dx - self.half_size[0]
        coords_1 = (global_idx[1] - .5) * self.dy - self.half_size[1]
        coords_2 = global_idx[2] * self.dz - self.half_size[2]
        
        return coords_0, coords_1, coords_2

    def spc_to_exact_hz_idx(self, x, y, z):
        """Return the exact mesh point of the given space coordinate.
        
        This method returns the (local) position, in index dimension, 
        of the nearest Hz mesh point of the given (global) space 
        coordinate. The return index could be out-of-range.
        
        Arguments:
            x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, float)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + .5
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy + .5
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz

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
        
        This method returns the (local) index of the nearest Hy mesh point of the
        given (global) space coordinate. The return index could be out-of-range.
        
        Arguments:
             x, y, z -- (global) space coordinate
        
        """
        coords = np.array((x,y,z), float)
            
        global_idx = empty(3, int)
        global_idx[0] = (coords[0] + self.half_size[0]) / self.dx + 1
        global_idx[1] = (coords[1] + self.half_size[1]) / self.dy + 1
        global_idx[2] = (coords[2] + self.half_size[2]) / self.dz + .5

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
        
        print " " * indent,
        print "size:", 2 * self.half_size,
        print "resolution:", self.res
        
        print " " * indent,
        print "dx:", self.dx, "dy:", self.dy, "dz:", self.dz
        
        print " " * indent,
        print "number of nodes participating:", self.numprocs


####################################################################
#                                                                  #
#                      Fast geometry routines                      #
#                                                                  #
# Using the geometry list is way too slow, especially when there   #
# are lots of objects to test.                                     #
#                                                                  #
# The basic idea here is twofold.  (1) Compute bounding boxes for  #
# each geometric object, for which inclusion tests can be computed #
# quickly.  (2) Build a tree that recursively breaks down the unit #
# cell in half, allowing us to perform searches in logarithmic     #
# time.                                                            #
#                                                                  #
####################################################################

class GeomBox(object):
    """A bounding box of a geometric object.
    
    Attributes:
        low -- the coordinates of the lowest vertex
        high -- the coordinates of the highest vertex
    
    """
    def __init__(self, low=None, high=None):
        self.low = np.array(low, float)
        self.high = np.array(high, float)
        
    def union(self, box):
        """Enlarge the box to include the given box.
        
        """
        self.low = map(min, self.low, box.low)
        self.high = map(max, self.high, box.high)
        
    def intersection(self, box):
        """Reduce the box to intersect volume with the given box.
        
        """
        self.low = map(max, self.low, box.low)
        self.high = map(min, self.high, box.high)        
        
    def add_point(self, point):
        """Enlarge the box to include the given point.
        
        """
        self.low = map(min, self.low, point)
        self.high = map(max, self.high, point)
           
    @staticmethod
    def __between__(x, low, high):
        """Return truth of low <= x <= high.
        
        """
        return low <= x <= high
    
    def in_box(self, point):
        """Check whether the given point is in the box.
        
        """
        componential = map(self.__between__, point, self.low, self.high)
        return reduce(lambda x, y: x and y, componential)
    
    def intersect(self, box):
        """Check whether the given box intersect with the box.
        
        """
        componential = (self.__between__(box.low[0], self.low[0], self.high[0]) or
                        self.__between__(box.high[0], self.low[0], self.high[0]) or
                        self.__between__(self.low[0], box.low[0], box.high[0])), \
                       (self.__between__(box.low[1], self.low[1], self.high[1]) or
                        self.__between__(box.high[1], self.low[1], self.high[1]) or
                        self.__between__(self.low[1], box.low[1], box.high[1])), \
                       (self.__between__(box.low[2], self.low[2], self.high[2]) or
                        self.__between__(box.high[2], self.low[2], self.high[2]) or
                        self.__between__(self.low[2], box.low[2], box.high[2]))
        return reduce(lambda x, y: x and y, componential)
    
    def divide(self, axis, x):
        high1 = deepcopy(self.high)
        high1[axis] = x
        
        low2 = deepcopy(self.low)
        low2[axis] = x
        
        return GeomBox(self.low, high1), GeomBox(low2, self.high)
    
    def display_info(self, indent=0):
        print " " * indent, "geom box:",
        print "low:", self.low, "high:", self.high
        
    def __str__(self):
        return "low: " + self.low.__str__() + " high: " + self.high.__str__()
    
    
class GeomBoxNode(object):
    """Node class which makes up a binary search tree.
    
    Attributes:
        box -- a bounding box enclosing the volume of this node
        t1 -- left branch from this node
        t2 -- right branch from this node
        geom_list -- a geometric object list overlapping the volume of this node.
        depth -- depth from the root of this binary search tree
        
    """
    def __init__(self, box, geom_list, depth):
        self.box = box
        self.t1, self.t2 = None, None
        self.geom_list = geom_list
        self.depth = depth
        
        
class GeomBoxTree(object):
    """A tree for the fast inclusion test of geometric objects within them.
    
    The tree recursively partitions the unit cell, allowing us to perform 
    binary searches for the object containing a give point.
    
    Attributes:
        root -- root node of the binary search tree
        
    """
    def __init__(self, geom_list=None):
        box = GeomBox((-inf, -inf, -inf), (inf, inf, inf))
        self.root = GeomBoxNode(box, geom_list, 0)
        self.branch_out(self.root)
    
    @staticmethod    
    def find_best_partition(node, divideAxis):
        """
        Find the best place to "cut" along the axis divideAxis in 
        order to maximally divide the objects between the partitions.
        Upon return, n1 and n2 are the number of objects below and 
        above the partition, respectively.
        
        """
        small = 1e-6 # only 1e-6 works
        bestPartition = None
        
        n1 = n2 = len(node.geom_list)
        
        # Search for the best partition, by checking all possible partitions 
        # either just above the high end of an object or just below the low 
        # end of an object. 
        
        for i in node.geom_list:
            curPartition = i.box.high[divideAxis] + small
            curN1 = curN2 = 0
            for j in node.geom_list:
                if j.box.low[divideAxis] <= curPartition:
                    curN1 += 1
                if j.box.high[divideAxis] >= curPartition:
                    curN2 += 1
            if max(curN1, curN2) < max(n1, n2):
                bestPartition = curPartition
                n1 = curN1
                n2 = curN2
                
        for i in node.geom_list:
            curPartition = i.box.low[divideAxis] - small
            curN1 = curN2 = 0
            for j in node.geom_list:
                if j.box.low[divideAxis] <= curPartition:
                    curN1 += 1
                if j.box.high[divideAxis] >= curPartition:
                    curN2 += 1
            if max(curN1, curN2) < max(n1, n2):
                bestPartition = curPartition
                n1 = curN1
                n2 = curN2
        
        return bestPartition, n1, n2
    
    @staticmethod
    def divide_geom_box_tree(node):
        """Divide box in two, along the axis that maximally partitions the boxes.
        
        """
        # Try partitioning along each dimension, counting the
        # number of objects in the partitioned boxes and finding
        # the best partition.
        best = 0
        division = []
        for i in xrange(3):
            partition, n1, n2 = GeomBoxTree.find_best_partition(node, i)
            division.append((partition, n1, n2))
            if max(division[i][1], division[i][2]) < max(division[best][1], division[best][2]):
                best = i
        
        # Don't do anything if division makes the worst case worse or if
        # it fails to improve the best case:       
        if division[best][0] is None:
            return None, None
        
        box1, box2 = node.box.divide(best, division[best][0])
        b1GeomList = []
        b2GeomList = []
        
        for i in node.geom_list:
            if box1.intersect(i.box):
                b1GeomList.append(i)
            if box2.intersect(i.box):
                b2GeomList.append(i)
        
        return GeomBoxNode(box1, b1GeomList, node.depth + 1), GeomBoxNode(box2, b2GeomList, node.depth + 1)
    
    @staticmethod
    def branch_out(node):
        node.t1, node.t2 = GeomBoxTree.divide_geom_box_tree(node)
        
        if not (node.t1 or node.t2):
            return
    
        GeomBoxTree.branch_out(node.t1)
        GeomBoxTree.branch_out(node.t2)
    
    @staticmethod
    def tree_search(node, point):
        if node.box.in_box(point) == False: 
            return None
        else:
            if not (node.t1 and node.t2):
                return node
            else:
                if node.t1.box.in_box(point):
                    return GeomBoxTree.tree_search(node.t1, point)
        
                if node.t2.box.in_box(point):
                    return GeomBoxTree.tree_search(node.t2, point)
    
    def object_of_point(self, point):
        leaf = self.tree_search(self.root, point)
        geom_obj, idx = find_object(point, leaf.geom_list)
        
        # epsilon and mu of compound material
        # refer to the underneath material.
        if isinstance(geom_obj.material, Compound):
            aux_geom_list = leaf.geom_list[:idx]
            underneath_obj, trash = find_object(point, aux_geom_list)
        else:
            underneath_obj = None
            
#            geom_obj.material.epsilon = geom_obj2.material.epsilon
#            geom_obj.material.mu = geom_obj2.material.mu
        
        return geom_obj, underneath_obj
        
    def material_of_point(self, point):
        geom_obj, underneath_obj = self.object_of_point(point)

        if underneath_obj: 
            underneath_material = underneath_obj.material
        else:
            underneath_material = None
            
        return geom_obj.material, underneath_material
        
    def display_info(self, node=None, indent=0):
        if not node: node = self.root
        
        print " " * indent, "depth:", node.depth, node.box
        for i in node.geom_list:
            print " " * (indent + 5), "bounding box:", i.box
            i.display_info(indent + 5)
            
        print " "

        if node.t1: self.display_info(node.t1, indent + 5)
        if node.t2: self.display_info(node.t2, indent + 5)
 
        
####################################################################
#                                                                  #
#                       Geometric primitives                       #
#                                                                  #
####################################################################

# TODO: Redefine the hierarchy of GeometricObject as like the one of Meep. 

class GeometricObject(object):
    """Base class for geometric object types.
    
    This class and its descendants are used to specify the solid 
    geometric objects that form the structure being simulated. One 
    normally does not create objects of type geometric-object directly,
    however; instead, you use one of the subclasses. Recall that 
    subclasses inherit the properties of their superclass, so these 
    subclasses automatically have the material property (which must be 
    specified, since they have no default values). In a two- or one-
    dimensional calculation, only the intersections of the objects with
    the simulation plane or line are considered. 
    
    Attributes:
        material -- Filling up material.
        box -- bounding box enclosing this geometric object
    
    """
    def __init__(self, material):
        """
        
        Arguments:
            material -- The material that the object is made of. No default.
            
            """ 
        self.material = material
        self.box = self.geom_box()
        
    def init(self, space):
        self.material.init(space)
        
    def geom_box(self):
        """Return a bounding box enclosing this geometric object.
        
        The derived classes should override this method.
        """
        
        raise NotImplementedError
        
    
    def in_object(self, point):
        """Return whether or not the point is inside.
        
        Return whether or not the point (in the lattice basis) is 
        inside this geometric object. This method additionally requires
        that fixObject has been called on this object (if the lattice 
        basis is non-orthogonal).
        
        """ 
        raise NotImplementedError
    
    def display_info(self, indent=0):
        """Display some information about this geometric object.
        
        """
        print " " * indent, "geometric object"
        print " " * indent, "center", self.center
        if self.material:
            self.material.display_info(indent + 5)
       
            
class DefaultMedium(GeometricObject):
    """A geometric object expanding the whole space.
    
    """
    def __init__(self, material):
        GeometricObject.__init__(self, material)
        
    def in_object(self, point):
        """
        Override GeometricObject.in_object.
        
        """
        return self.geom_box().in_box(point)

    def geom_box(self):
        """
        Override GeometriObject.geom_box.
        
        """
        return GeomBox((-inf, -inf, -inf), (inf, inf, inf))

    def display_info(self, indent=0):
        """
        Override GeometricObject.display_info.
        
        """
        print " " * indent, "default material"
        if self.material:
            self.material.display_info(indent + 5)
    
    
class Cone(GeometricObject):
    """Form a cone or possibly a truncated cone. 
    
    Attributes:
        center -- coordinates of the center of this geometric object
        axis -- unit vector of axis
        height -- length of axis
        box -- bounding box
        
    """
    def __init__(self, material, radius2=0, axis=(1, 0, 0),
                 radius=1, height=1, center=(0, 0, 0)):
        """
        
        Arguments:
        radius2 -- Radius of the tip of the cone (i.e. the end of the 
            cone pointed to by the axis vector). Defaults to zero 
            (a "sharp" cone).
            
        """
        if radius < 0:
            msg = "radius must be non-negative."
            raise ValueError(msg)
        else:
            self.radius = float(radius) # low side radius
            
        if radius2 < 0:
            msg = "radius2 must be non-negative."
            raise ValueError(msg)
        else:
            self.radius2 = float(radius2) # high side radius

        self.center = np.array(center, float)
        self.axis = np.array(axis, float) / norm(axis)
        self.height = float(height)
        
        GeometricObject.__init__(self, material)
        
    def in_object(self, x):
        """
        Override GeometricObject.in_object.
        
        """
        x = np.array(x, float)
        r = x - self.center
        proj = dot(self.axis, r)
        if np.abs(proj) <= .5 * self.height:
            if self.radius2 == self.radius == inf:
                return True
            radius = self.radius
            radius += (proj / self.height + .5) * (self.radius2 - radius)
            truth = radius != 0 and norm(r - proj * self.axis) <= np.abs(radius)
            return truth
        else:
            return False
           
    def display_info(self, indent=0):
        """
        Override GeometricObject.display_info.
        
        """
        print " " * indent, "cone"
        print " " * indent,
        print "center", self.center,
        print "radius", self.radius,
        print "height" , self.height,
        print "axis", self.axis,
        print "radius2", self.radius2
        if self.material:
            self.material.display_info(indent + 5)
        
    def geom_box(self):
        """
        Override GeometricObject.geom_box.
        
        """
        tmpBox1 = GeomBox(low=self.center, high=self.center)

        h = .5 * self.height
        r = np.sqrt(1 - self.axis * self.axis)

        # set tmpBox2 to center of object
        tmpBox2 = deepcopy(tmpBox1)
        
        # bounding box for -h*axis cylinder end
        tmpBox1.low -= h * self.axis + r * self.radius
        tmpBox1.high -= h * self.axis - r * self.radius
        
        # bounding box for +h*axis cylinder end
        tmpBox2.low += h * self.axis - r * self.radius
        tmpBox2.high += h * self.axis + r * self.radius
        
        tmpBox1.union(tmpBox2)
        
        return tmpBox1
       
        
class Cylinder(Cone):
    """Form a cylinder.
    
    """
    def __init__(self, material, axis=(0, 0, 1),
                 radius=1, height=1, center=(0, 0, 0)):
        """
        Arguments:
            axis -- Direction of the cylinder's axis; the length of 
                this vector is ignored. Defaults to point parallel to 
                the z axis i.e., (0,0,1).
            radius -- Radius of the cylinder's cross-section. Default is 1.
            height -- Length of the cylinder along its axis. Default is 1. 
            material -- The material that the object is made of. 
                No default.
            center -- Center point of the object. Default is (0,0,0). 
        
        """
        Cone.__init__(self, material, radius, axis, radius, height, center)
        
    def display_info(self, indent=0):
        """
        Override Cone.display_info.
        
        """
        print " " * indent, "cylinder"
        print " " * indent,
        print "center", self.center,
        print "radius", self.radius,
        print "height", self.height,
        print "axis", self.axis
        if self.material:
            self.material.display_info(indent + 5)


class _Block(GeometricObject):
    """Base class for Block and Ellipsoid class.
    
    """
    def __init__(self, material,
                 e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1),
                 size=(0, 0, 0), center=(0, 0, 0)):
        self.center = np.array(center, float)
        
        self.e1 = np.array(e1, float) / norm(e1)
        self.e2 = np.array(e2, float) / norm(e2)
        self.e3 = np.array(e3, float) / norm(e3)
        self.size = np.array(size, float)
        
        self.projection_matrix = np.array([self.e1, self.e2, self.e3])
        
        GeometricObject.__init__(self, material)
        
    def in_object(self, x):
        """
        Override GeometricObject.in_object.
        
        """
        x = np.array(x, float)
        r = x - self.center
        proj = dot(self.projection_matrix, r)

        return (np.abs(proj) <= .5 * self.size).all()
        
    def geom_box(self):
        """
        Override GeometricObject.geom_box.
        
        """
        tmpBox = GeomBox(low=self.center, high=self.center)
        # enlarge the box to be big enough to contain all 8 corners
        # of the block.
        s1 = self.size[0] * self.e1
        s2 = self.size[1] * self.e2
        s3 = self.size[2] * self.e3
        
        corner = self.center + (-0.5 * (s1 + s2 + s3))
        
        tmpBox.add_point(corner)
        tmpBox.add_point(corner + s1)
        tmpBox.add_point(corner + s2)
        tmpBox.add_point(corner + s3)
        tmpBox.add_point(corner + s1 + s2)
        tmpBox.add_point(corner + s1 + s3)
        tmpBox.add_point(corner + s3 + s2)
        tmpBox.add_point(corner + s1 + s2 + s3)
        
        return tmpBox
        
        
class Block(_Block):
    """Form a parallelpiped (i.e., a brick, possibly with non-orthogonal axes).
    
    """
    def __init__(self, material,
                 e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1),
                 size=(1, 1, 1), center=(0, 0, 0)):
        """
        
        Arguments:
            size -- The lengths of the block edges along each of its 
                three axes. Default is (1, 1, 1).
            e1, e2, e3 -- The directions of the axes of the block; the 
                lengths of these vectors are ignored. Must be linearly 
                independent. They default to the three Cartesian axis.
                
        """
        _Block.__init__(self, material, e1, e2, e3, size, center)
        
    def display_info(self, indent=0):
        print " " * indent, "block"
        print " " * indent,
        print "center", self.center,
        print "size", self.size,
        print "axes", self.e1, self.e2, self.e3
        if self.material:
            self.material.display_info(indent + 5)


class Ellipsoid(_Block):
    """Form an ellipsoid.
    
    """
    def __init__(self, material,
                 e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1),
                 size=(1, 1, 1), center=(0, 0, 0)):
        _Block.__init__(self, material, e1, e2, e3, size, center)

        self.inverse_semi_axes = 2 / np.array(size, float)

    def in_object(self, x):
        x = np.array(x, float)
        r = x - self.center
        proj = dot(self.projection_matrix, r)
        q = proj * self.inverse_semi_axes
        return sum(q * q) <= 1

    def display_info(self, indent=0):
        print " " * indent, "ellipsoid"
        print " " * indent,
        print "center", self.center,
        print "size", self.size,
        print "axis", self.e1, self.e2, self.e3
        if self.material:
            self.material.display_info(indent + 5)


class Sphere(Ellipsoid):
    """Form a sphere.
    
    Attributes:
        radius -- Radius of the sphere.
    
    """
    def __init__(self, material, radius=1, center=(0, 0, 0)):
        """
        
        Arguments:
            radius -- Radius of the sphere. Default is 1. 
            material -- The material that the object is made of. 
                No default.
            center -- Center point of the object. Default is (0,0,0).
        """
        if radius < 0:
            msg = "radius must be non-negative."
            raise ValueError(msg)
        else:
            self.radius = radius
            
        Ellipsoid.__init__(self, material, (1, 0, 0), (0, 1, 0), (0, 0, 1),
                           (2 * radius, 2 * radius, 2 * radius), center)

    def display_info(self, indent=0):
        print " " * indent, "sphere"
        print " " * indent,
        print "center", self.center,
        print "radius", self.radius
        if self.material:
            self.material.display_info(indent + 5)


class Boundary(GeometricObject):
    """Form a boundary.
     
    """
    def __init__(self, material, thickness=None, size=None,
                 plus_x=True, minus_x=True,
                 plus_y=True, minus_y=True,
                 plus_z=True, minus_z=True):
        """
        
         Arguments:
             material --
             thickness -- The spatial thickness of the Boundary layer 
                 (which extends from the boundary towards the inside of
                 the computational cell). Default value.
             size --
             plus_x -- Specify whether the high of the boundary in 
                 direction x is set. Default is True.
             minus_x -- Specify whether the low of the boundary in 
                 direction x is set. Default is True.
             plus_y -- Specify whether the high of the boundary in 
                 direction y is set. Default is True.
             minus_y -- Specify whether the low of the boundary in 
                 direction y is set. Default is True.
             plus_z -- Specify whether the high of the boundary in 
                 direction z is set. Default is True.
             minus_z -- Specify whether the low of the boundary in 
                 direction z is set. Default is True.
        
        """
        self.d = float(thickness)
        
        self.half_size = .5 * np.array(size, float)
        
        self.box_list = []
        
        self.minus_x, self.plus_x = minus_x, plus_x
        self.minus_y, self.plus_y = minus_y, plus_y
        self.minus_z, self.plus_z = minus_z, plus_z

        # do someting for the PML derived class?
        if isinstance(material, PML):
            pass
        
        GeometricObject.__init__(self, material)
                
    def init(self, space):
        if 2 * self.half_size[0] < space.dx:
            self.half_size[0] = 0.5 * space.dx
        
        if 2 * self.half_size[1] < space.dy:
            self.half_size[1] = 0.5 * space.dy
            
        if 2 * self.half_size[2] < space.dz:
            self.half_size[2] = 0.5 * space.dz
            
        if 2 * self.half_size[0] > space.dx:
            if self.plus_x:
                low = (self.half_size[0] - self.d, -self.half_size[1],
                       -self.half_size[2])
                high = self.half_size[:]
                self.box_list.append(GeomBox(low, high))
            if self.minus_x:
                low = -self.half_size[:]
                high = (-self.half_size[0] + self.d, self.half_size[1],
                        self.half_size[2])
                self.box_list.append(GeomBox(low, high))
                
        if 2 * self.half_size[1] > space.dy:
            if self.plus_y:
                low = (-self.half_size[0], self.half_size[1] - self.d,
                       -self.half_size[2])
                high = self.half_size[:]
                self.box_list.append(GeomBox(low, high))
            if self.minus_y:
                low = -self.half_size[:]
                high = (self.half_size[0], -self.half_size[1] + self.d,
                        self.half_size[2])
                self.box_list.append(GeomBox(low, high))

        if 2 * self.half_size[2] > space.dz:
            if self.plus_z:
                low = (-self.half_size[0], -self.half_size[1],
                       self.half_size[2] - self.d)
                high = self.half_size[:]
                self.box_list.append(GeomBox(low, high))
            if self.minus_z:
                low = -self.half_size[:]
                high = (self.half_size[0], self.half_size[1],
                        -self.half_size[2] + self.d)
                self.box_list.append(GeomBox(low, high))
        
        self.material.init(space, self.d)
        
    def in_object(self, point):
        for box in self.box_list:
            if box.in_box(point):
                return True
        return False
        
    def geom_box(self):
        return GeomBox(-self.half_size, self.half_size)
        
    def display_info(self, indent=0):
        print " " * indent, "boundary"
        print " " * indent,
        print "+x:", self.plus_x, "-x:", self.minus_x,
        print "+y:", self.plus_y, "-y:", self.minus_y,
        print "+z:", self.plus_z, "-z:", self.minus_z
        if self.material:
            self.material.display_info(indent + 5)
    
    
def find_object(point, geom_list):
    """Find the last object including point in geom_list.
    
    find_object returns (object, array index). If no object includes 
    the given point it returns (geom_list[0], 0).
    """
    for i in range(len(geom_list) - 1, -1, -1):
        if geom_list[i].in_object(point):
            break
        
    return geom_list[i], i


def in_range(idx, numpy_array, component):
    """Perform bounds checking.
    
    
    Keyword arguments:
        idx -- index of an array
        numpy_array -- the array to be checked
        component -- specify field component
        
    """
    if component is const.Ex:
        if idx[0] < 0 or idx[0] >= numpy_array.shape[0]:
            return False
        if idx[1] < 0 or idx[1] >= numpy_array.shape[1] - 1:
            return False
        if idx[2] < 0 or idx[2] >= numpy_array.shape[2] - 1:
            return False
        
    elif component is const.Ey:
        if idx[0] < 0 or idx[0] >= numpy_array.shape[0] - 1:
            return False
        if idx[1] < 0 or idx[1] >= numpy_array.shape[1]:
            return False
        if idx[2] < 0 or idx[2] >= numpy_array.shape[2] - 1:
            return False
          
    elif component is const.Ez:
        if idx[0] < 0 or idx[0] >= numpy_array.shape[0] - 1:
            return False
        if idx[1] < 0 or idx[1] >= numpy_array.shape[1] - 1:
            return False
        if idx[2] < 0 or idx[2] >= numpy_array.shape[2]:
            return False
        
    elif component is const.Hx:
        if idx[0] < 0 or idx[0] >= numpy_array.shape[0]:
            return False
        if idx[1] <= 0 or idx[1] >= numpy_array.shape[1]:
            return False
        if idx[2] <= 0 or idx[2] >= numpy_array.shape[2]:
            return False
        
    elif component is const.Hy:
        if idx[0] <= 0 or idx[0] >= numpy_array.shape[0]:
            return False
        if idx[1] < 0 or idx[1] >= numpy_array.shape[1]:
            return False
        if idx[2] <= 0 or idx[2] >= numpy_array.shape[2]:
            return False
        
    elif component is const.Hz:
        if idx[0] <= 0 or idx[0] >= numpy_array.shape[0]:
            return False
        if idx[1] <= 0 or idx[1] >= numpy_array.shape[1]:
            return False
        if idx[2] < 0 or idx[2] >= numpy_array.shape[2]:
            return False
        
    else:
        raise ValueError
       
    return True


if __name__ == '__main__':
    from material import *
    
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
    
