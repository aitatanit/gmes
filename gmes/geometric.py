#!/usr/bin/env python

# this code is based on the libctl 3.0.2.

from copy import deepcopy
from math import ceil
from numpy import array, zeros, empty, dot, inf, sqrt, fabs
from scipy.linalg import norm

from constants import c0


class Cartesian:
    """Define the whole calculation space.
    size: a sequence with length three consists of
    non-negative numbers.
    a: physical length of a unit in the coordinate. (scalar)
    resolution: number of sections of one unit (scalar or length 3 sequence)
    """
    
    def __init__(self, size, a=1, resolution=15, courant=.5, dt=None):
        self.a = float(a)
        self.half_size = 0.5 * array(size, float)
        try:
            if len(resolution) == 3:
                self.res = array(resolution, float)
        except TypeError:
            self.res = array((resolution,)*3, float)
            
        self.courant = float(courant)
        
        self.dx, self.dy, self.dz = self.a / self.res

        # local spatial differentials to calculate dt
        dx, dy, dz = self.a / self.res
        
        if self.half_size[0] == 0 or self.half_size[0] == inf:
            self.half_size[0] = .5 * self.dx
            dx = inf
            
        if self.half_size[1] == 0 or self.half_size[1] == inf:
            self.half_size[1] = .5 * self.dy
            dy = inf
            
        if self.half_size[2] == 0 or self.half_size[2] == inf:
            self.half_size[2] = .5 * self.dz
            dz = inf
            
        if dt is None:
            self.dt = self.courant / c0 / sqrt(dx**-2 + dy**-2 + dz**-2)
        else:
            self.dt = float(dt)
            self.courant = self.dt * c0 * sqrt(dx**-2 + dy**-2 + dz**-2)

        self.field_size = []
        for i in (int(x + .5) for x in 2 * self.half_size * self.res):
            if i == 0:
                self.field_size.append(1)
            else:
                self.field_size.append(i)

        self.ex_shape = self.field_size[0], self.field_size[1] + 1, self.field_size[2] + 1
        self.ey_shape = self.field_size[0] + 1, self.field_size[1], self.field_size[2] + 1
        self.ez_shape = self.field_size[0] + 1, self.field_size[1] + 1, self.field_size[2]
        self.hx_shape = self.ex_shape
        self.hy_shape = self.ey_shape
        self.hz_shape = self.ez_shape
        
    def get_ex_storage(self):
        """Return an initialized array for Ex field component.
        """
        return zeros(self.ex_shape, float)
    
    def get_ey_storage(self):
        return zeros(self.ey_shape, float)
    
    def get_ez_storage(self):
        return zeros(self.ez_shape, float)
    
    def get_hx_storage(self):
        return zeros(self.hx_shape, float)
    
    def get_hy_storage(self):    
        return zeros(self.hy_shape, float)
    
    def get_hz_storage(self):
        return zeros(self.hz_shape, float)

    def get_material_ex_storage(self):
        return empty(self.ex_shape, object)
    
    def get_material_ey_storage(self):
        return empty(self.ey_shape, object)
    
    def get_material_ez_storage(self):
        return empty(self.ez_shape, object)
    
    def get_material_hx_storage(self):
        return empty(self.hx_shape, object)
    
    def get_material_hy_storage(self):    
        return empty(self.hy_shape, object)
    
    def get_material_hz_storage(self):
        return empty(self.hz_shape, object)

    def ex_index_to_space(self, idx):
        """Return space coordinates corresponding to the given index of Ex mesh point. 
        """
        
        x0 = (idx[0] + .5) * self.dx - self.half_size[0]
        x1 = idx[1] * self.dy - self.half_size[1]
        x2 = idx[2] * self.dz - self.half_size[2]
        return x0, x1, x2
        
    def space_to_ex_index(self, x):
        """Return the index of the nearest Ex mesh point.
        
        x: space coordinates (length 3)
        """
        
        idx0 = int((x[0] + self.half_size[0]) / self.dx)
        if idx0 < 0: idx0 = 0
        elif idx0 > self.ex_shape[0] - 1: idx0 = self.ex_shape[0] - 1
        
        idx1 = int((x[1] + self.half_size[1]) / self.dy + .5)
        if idx1 < 0: idx1 = 0
        elif idx1 > self.ex_shape[1] - 2: idx1 = self.ex_shape[1] - 2

        idx2 = int((x[2] + self.half_size[2]) / self.dz + .5)
        if idx2 < 0: idx2 = 0
        elif idx2 > self.ex_shape[2] - 2: idx2 = self.ex_shape[2] - 2

        return idx0, idx1, idx2

    def ey_index_to_space(self, idx):
        x0 = idx[0] * self.dx - self.half_size[0]
        x1 = (idx[1] + .5) * self.dy - self.half_size[1]
        x2 = idx[2] * self.dz - self.half_size[2]
        return x0, x1, x2
    
    def space_to_ey_index(self, x):
        idx0 = int((x[0] + self.half_size[0]) / self.dx + .5)
        if idx0 < 0: idx0 = 0
        elif idx0 > self.ey_shape[0] - 2: idx0 = self.ey_shape[0] - 2

        idx1 = int((x[1] + self.half_size[1]) / self.dy)
        if idx1 < 0: idx1 = 0
        elif idx1 > self.ey_shape[1] - 1: idx1 = self.ey_shape[1] - 1

        idx2 = int((x[2] + self.half_size[2]) / self.dz + .5)
        if idx2 < 0: idx2 = 0
        elif idx2 > self.ey_shape[2] - 2: idx2 = self.ey_shape[2] - 2

        return idx0, idx1, idx2
    
    def ez_index_to_space(self, idx):
        x0 = idx[0] * self.dx - self.half_size[0]
        x1 = idx[1] * self.dy - self.half_size[1]
        x2 = (idx[2] + .5) * self.dz - self.half_size[2]
        return x0, x1, x2
    
    def space_to_ez_index(self, x):
        idx0 = int((x[0] + self.half_size[0]) / self.dx + .5)
        if idx0 < 0: idx0 = 0
        elif idx0 > self.ez_shape[0] - 2: idx0 = self.ez_shape[0] - 2

        idx1 = int((x[1] + self.half_size[1]) / self.dy + .5)
        if idx1 < 0: idx1 = 0
        elif idx1 > self.ez_shape[1] - 2: idx1 = self.ez_shape[1] - 2

        idx2 = int((x[2] + self.half_size[2]) / self.dz)
        if idx2 < 0: idx2 = 0
        elif idx2 > self.ez_shape[2] - 1: idx2 = self.ez_shape[2] - 1

        return idx0, idx1, idx2

    def hx_index_to_space(self, idx):
        x0 = idx[0] * self.dx - self.half_size[0]
        x1 = (idx[1] - .5) * self.dy - self.half_size[1]
        x2 = (idx[2] - .5) * self.dz - self.half_size[2]
        return x0, x1, x2

    def space_to_hx_index(self, x):
        idx0 = int((x[0] + self.half_size[0]) / self.dx + .5)
        if idx0 < 0: idx0 = 0
        elif idx0 > self.hx_shape[0] - 1: idx0 = self.hx_shape[0] - 1
        
        idx1 = int((x[1] + self.half_size[1]) / self.dy + 1)
        if idx1 < 1: idx1 = 1
        elif idx1 > self.hx_shape[1] - 1: idx1 = self.hx_shape[1] - 1

        idx2 = int((x[2] + self.half_size[2]) / self.dz + 1)
        if idx2 < 1: idx2 = 1
        elif idx2 > self.hx_shape[2] - 1: idx2 = self.hx_shape[2] - 1

        return idx0, idx1, idx2

    def hy_index_to_space(self, idx):
        x0 = (idx[0] - .5) * self.dx - self.half_size[0]
        x1 = idx[1] * self.dy - self.half_size[1]
        x2 = (idx[2] - .5) * self.dz - self.half_size[2]
        return x0, x1, x2
        
    def space_to_hy_index(self, x):
        idx0 = int((x[0] + self.half_size[0]) / self.dx + 1)
        if idx0 < 1: idx0 = 1
        elif idx0 > self.hy_shape[0] - 1: idx0 = self.hy_shape[0] - 1
        
        idx1 = int((x[1] + self.half_size[1]) / self.dy + .5)
        if idx1 < 0: idx1 = 0
        elif idx1 > self.hy_shape[1] - 1: idx1 = self.hy_shape[1] - 1

        idx2 = int((x[2] + self.half_size[2]) / self.dz + 1)
        if idx2 < 1: idx2 = 1
        elif idx2 > self.hy_shape[2] - 1: idx2 = self.hy_shape[2] - 1

        return idx0, idx1, idx2

    def hz_index_to_space(self, idx):
        x0 = (idx[0] - .5) * self.dx - self.half_size[0]
        x1 = (idx[1] - .5) * self.dy - self.half_size[1]
        x2 = idx[2] * self.dz - self.half_size[2]
        return x0, x1, x2
        
    def space_to_hz_index(self, x):
        idx0 = int((x[0] + self.half_size[0]) / self.dx + 1)
        if idx0 < 1: idx0 = 1
        elif idx0 > self.hz_shape[0] - 1: idx0 = self.hz_shape[0] - 1

        idx1 = int((x[1] + self.half_size[1]) / self.dy + 1)
        if idx1 < 1: idx1 = 1
        elif idx1 > self.hz_shape[1] - 1: idx1 = self.hz_shape[1] - 1

        idx2 = int((x[2] + self.half_size[2]) / self.dz + .5)
        if idx2 < 0: idx2 = 0
        elif idx2 > self.hz_shape[2] - 1: idx2 = self.hz_shape[2] - 1
        
        return idx0, idx1, idx2
    
    def display_info(self, indent=0):
        print " " * indent, "Cartesian space"
        print " " * indent,
        print "size:", 2 * self.half_size,
        print "a:", self.a,
        print "resolution:", self.res
        print "dx:", self.dx, "dy:", self.dy, "dz:", self.dz
        print "dt:", self.dt
        
    
class GeomBox:
    def __init__(self, low=None, high=None):
        self.low = array(low, float)
        self.high = array(high, float)
        
    def union(self, box):
        self.low = map(min, self.low, box.low)
        self.high = map(max, self.high, box.high)
        
    def intersection(self, box):
        self.low = map(max, self.low, box.low)
        self.high = map(min, self.high, box.high)        
        
    def add_point(self, x):
        self.low = map(min, self.low, x)
        self.high = map(max, self.high, x)
           
    @staticmethod
    def __between__(x, low, high):
        return low <= x <= high
    
    def in_box(self, x):
        componential = map(self.__between__, x, self.low, self.high)
        return reduce(lambda x, y: x and y, componential)
    
    def intersect(self, box):
        """
        Return whether or not the given box intersect with this.
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
    
    
class GeomBoxNode:
    def __init__(self, box, geom_list, depth):
        self.box = box
        self.t1, self.t2 = None, None
        self.geom_list = geom_list
        self.depth = depth
        
        
class GeomBoxTree:
    def __init__(self, geom_list=None):
        box = GeomBox((-inf, -inf, -inf), (inf, inf, inf))
        self.root = GeomBoxNode(box, geom_list, 0)
        self.branch_out(self.root)
    
    @staticmethod    
    def find_best_partition(node, divideAxis):
        """ Find the best place to "cut" along the axis divideAxis in 
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
        """divide_geom_box_tree divide box in two, along the axis
        that maximally partitions the boxes. 
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
        
        return GeomBoxNode(box1, b1GeomList, node.depth+1), GeomBoxNode(box2, b2GeomList, node.depth+1)
    
    @staticmethod
    def branch_out(node):
        node.t1, node.t2 = GeomBoxTree.divide_geom_box_tree(node)
        
        if not (node.t1 or node.t2):
            return
    
        GeomBoxTree.branch_out(node.t1)
        GeomBoxTree.branch_out(node.t2)
    
    @staticmethod
    def tree_search(node, x):
        if node.box.in_box(x) == False: 
            return None
        else:
            if not (node.t1 and node.t2):
                return node
            else:
                if node.t1.box.in_box(x):
                    return GeomBoxTree.tree_search(node.t1,x)
        
                if node.t2.box.in_box(x):
                    return GeomBoxTree.tree_search(node.t2,x)
    
    def object_of_point(self, x):
        leaf = self.tree_search(self.root,x)
        for go in reversed(leaf.geom_list):
            if go.in_object(x):
                return go
        
#    def materialOfPoint(self, x):
#        go = self.object_of_point(x)
#        if (go is not None): return go.material
        
    def display_info(self, node=None, indent=0):
        if not node: node = self.root
        
        print " " * indent, "depth:", node.depth, node.box
        for i in node.geom_list:
            print " " * (indent + 5), "bounding box:", i.box
            i.display_info(indent + 5)
            
        print " "

        if node.t1: self.display_info(node.t1,indent + 5)
        if node.t2: self.display_info(node.t2,indent + 5)
        
        
class GeometricObject:
    def __init__(self, material=None, center=(0,0,0)):
        self.material = material
        self.center = array(center, float)
    
    def init(self, space):
        pass
        
    def in_object(self, x):
        """Return whether or not the x (in the lattice basis) is inside
        this geometric object. This method additionally requires that 
        fixObject has been called on this object (if the lattice basis is
        non-orthogonal).
        """
        return False
    
    def display_info(self, indent=0):
        """Display some information about this geometric object.
        """
        print " " * indent, "geometric object"
        print " " * indent, "center", self.center
        if self.material:
            self.material.display_info(indent + 5)
       
            
class DefaultMaterial(GeometricObject):
    def __init__(self, material=None):
        self.material = material
        self.box = self.geom_box()
        
    def in_object(self, x):
        return self.geom_box().in_box(x)

    def geom_box(self):
        return GeomBox((-inf,-inf,-inf), (inf,inf,inf))

    def display_info(self, indent=0):
        print " " * indent, "default material"
        if self.material:
            self.material.display_info(indent + 5)
    
    
class Cone(GeometricObject):
    def __init__(self, radius2=0, axis=(1,0,0), radius=1, height=1, material=None, center=(0,0,0)):
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

        GeometricObject.__init__(self, material, center)

        self.axis = array(axis, float) / norm(axis)
        self.height = float(height)
        
        self.box = self.geom_box()
        
    def in_object(self, x):
        x = array(x, float)
        r = x - self.center
        proj = dot(self.axis, r)
        if fabs(proj) <= .5 * self.height:
            if self.radius2 == self.radius == inf:
                return True
            radius = self.radius
            radius += (proj / self.height + .5) * (self.radius2 - radius)
            truth = radius != 0 and \
                    norm(r - proj * self.axis) <= fabs(radius)
            return truth
        else:
            return False
           
    def display_info(self, indent=0):
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
        tmpBox1 = GeomBox(low=self.center, high=self.center)

        h = .5 * self.height
        r = sqrt(1 - self.axis * self.axis)
        
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
    def __init__(self, axis=(1,0,0), radius=1, height=1, material=None, center=(0,0,0)):
        Cone.__init__(self, radius, axis, radius, height, material, center)
        
    def display_info(self, indent=0):
        print " " * indent, "cylinder"
        print " " * indent, 
        print "center", self.center,
        print "radius", self.radius,
        print "height" , self.height,
        print "axis", self.axis
        if self.material:
            self.material.display_info(indent + 5)


class _Block(GeometricObject):
    def __init__(self, e1=(1,0,0), e2=(0,1,0), e3=(0,0,1), size=(0,0,0), material=None, center=(0,0,0)):
        GeometricObject.__init__(self, material, center)
        self.e1 = array(e1, float) / norm(e1)
        self.e2 = array(e2, float) / norm(e2)
        self.e3 = array(e3, float) / norm(e3)
        self.size = array(size, float)
        
        self.projection_matrix = array([self.e1, self.e2, self.e3])
        
        self.box = self.geom_box()
        
    def in_object(self, x):
        x = array(x, float)
        r = x - self.center
        proj = dot(self.projection_matrix, r)

        return (fabs(proj) <= .5 * self.size).all()
        
    def geom_box(self):
        tmpBox = GeomBox(low=self.center, high=self.center)
        # enlarge the box to be big enough to contain all 8 corners of the block.
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
    def __init__(self, e1=(1,0,0), e2=(0,1,0), e3=(0,0,1), size=(0,0,0), material=None, center=(0,0,0)):
        _Block.__init__(self, e1, e2, e3, size, material, center)
        
    def display_info(self, indent=0):
        print " " * indent, "block"
        print " " * indent,
        print "center", self.center,
        print "size", self.size,
        print "axes", self.e1, self.e2, self.e3
        if self.material:
            self.material.display_info(indent + 5)


class Ellipsoid(_Block):
    def __init__(self, e1=(1,0,0), e2=(0,1,0), e3=(0,0,1), size=(1,1,1), material=None, center=(0,0,0)):
        _Block.__init__(self, e1, e2, e3, size, material, center)

        self.inverse_semi_axes = 2 / array(size, float)

    def in_object(self, x):
        x = array(x, float)
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
    def __init__(self, radius=1, material=None, center=(0,0,0)):
        if radius < 0:
            msg = "radius must be non-negative."
            raise ValueError(msg)
        else:
            self.radius = radius
            
        Ellipsoid.__init__(self, (1,0,0), (0,1,0), (0,0,1), (radius,radius,radius), material, center)

    def display_info(self, indent=0):
        print " " * indent, "sphere"
        print " " * indent,
        print "center", self.center,
        print "radius", self.radius
        if self.material:
            self.material.display_info(indent + 5)


class Boundary:
    def __init__(self, material=None, thickness=None, size=None, plusX=True, minusX=True, plusY=True, minusY=True, plusZ=True, minusZ=True):
        self.d = float(thickness)
        self.size = array(size, float)
        half_size = []
        for i in self.size:
            if i <= 2 * self.d:
                i = inf
            else:
                i *= 0.5
            half_size.append(i)
        self.half_size = array(half_size)

        self.boxList = []
        if size[0]:
            if plusX:
                low = -self.half_size
                low[0] = self.half_size[0] - self.d
                high = self.half_size
                self.boxList.append(GeomBox(low, high))
            if minusX:
                low = -self.half_size
                high = deepcopy(self.half_size)
                high[0] = -self.half_size[0] + self.d
                self.boxList.append(GeomBox(low, high))
        if size[1]:
            if plusY:
                low = -self.half_size
                low[1] = self.half_size[1] - self.d
                high = self.half_size
                self.boxList.append(GeomBox(low, high))
            if minusX:
                low = -self.half_size
                high = deepcopy(self.half_size)
                high[1] = -self.half_size[1] + self.d
                self.boxList.append(GeomBox(low, high))   
        if size[2]:   
            if plusZ:
                low = -self.half_size
                low[2] = self.half_size[2] - self.d
                high = self.half_size
                self.boxList.append(GeomBox(low, high))
            if minusZ:
                low = -self.half_size
                high = deepcopy(self.half_size)
                high[2] = -self.half_size[2] + self.d
                self.boxList.append(GeomBox(low, high)) 
                    
        self.plusX = plusX
        self.minusX = minusX
        self.plusY = plusY
        self.minusY = minusY
        self.plusZ = plusZ
        self.minusZ = minusZ
         
        if material.__class__ .__name__ == 'CPML':
            pass
        
        elif material.__class__.__name__ == 'UPML':
            pass
        
        self.material = material
        
        self.box = self.geom_box()

    def init(self, space):
        self.material.init(self.d, space)
        
    def in_object(self, x):
        for box in self.boxList:
            if box.in_box(x):
                return True
        return False
    
    def geom_box(self):
        return GeomBox(-self.half_size, self.half_size)
        
    def display_info(self, indent=0):
        print " " * indent, "boundary"
        print " " * indent,
        print "+x:", self.plusX, "-x:", self.minusX,
        print "+y:", self.plusY, "-y:", self.minusY,
        print "+z:", self.plusZ, "-z:", self.minusZ
        if self.material:
            self.material.display_info(indent + 5)
    
    
def find_object(x, geom_list):
    for go in reversed(geom_list):
        if go.in_object(x):
            return go


if __name__ == '__main__':
    from material import Dielectric   
    
    geom_list = [DefaultMaterial(material=Dielectric()), Cone(0, (1,0,0), 1, 1, Dielectric(), (0,0,2)), Cone(0, (1,0,0), 1, 1, Dielectric(), (0,0,-2))]
    t = GeomBoxTree(geom_list)
    t.display_info()
    space = Cartesian(size=(5,5,5))
    ex = space.get_ex_storage()
    print "ex shape:", ex.shape
    print "ex:", space.ex_index_to_space((0,0,0))
    print "ex:", space.ex_index_to_space((74,75,75))
    print "ex:", space.space_to_ex_index((0,0,0))
    print "ey:", space.space_to_ey_index((0,0,0))
    print "ez:", space.space_to_ez_index((0,0,0))
    print "hx:", space.space_to_hx_index((0,0,0))
    print "hy:", space.space_to_hy_index((0,0,0))
    print "hz:", space.space_to_hz_index((0,0,0))

    
