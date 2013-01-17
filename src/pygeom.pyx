# -*- coding: utf-8 -*-
# cython: boundscheck=False
# cython: wraparound=False

from __future__ import division
from copy import deepcopy
from scipy import sqrt
import numpy as np
from material import Compound, Pml

cimport numpy as np
np.import_array()
cimport cython


cdef double norm(object p):
    return sqrt(p[0]**2 + p[1]**2 + p[2]**2)


cdef class Material(object):
    """A base class for material types.
    
    """
    cdef public double eps_inf, mu_inf

    def __init__(self, eps_inf=1, mu_inf=1):
        self.eps_inf = float(eps_inf)
        self.mu_inf = float(mu_inf)

    def __getstate__(self):
        d = {}
        d['eps_inf'] = self.eps_inf
        d['mu_inf'] = self.mu_inf
        return d

    def __reduce__(self):
        return self.__class__, (), self.__getstate__()

    def __setstate__(self, d):
        self.eps_inf = d['eps_inf']
        self.mu_inf = d['mu_inf']

    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        raise NotImplementedError

    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        """Return an ElectricParam structure of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        """Return an ElectricParam structure of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        """Return an ElectricParam structure of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        """Return a MagneticParam structure of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        """Return a MagneticParam structure of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        """Return a MagneticParam structure of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError

    def init(self, space, param=None):
        raise NotImplementedError


####################################################################
#                                                                  #
#                      Fast geometry routines                      #
#                                                                  #
# Using the geometry list is way too slow, especially when there   #
# are lots of objects to test.                                     #
#                                                                  #
# The basic idea here is twofold. (1) Compute bounding boxes for   #
# each geometric object, for which inclusion tests can be computed #
# quickly. (2) Build a tree that recursively breaks down the unit  #
# cell in half, allowing us to perform searches in logarithmic     #
# time.                                                            #
#                                                                  #
####################################################################

cdef tuple find_object(tuple point, tuple geom_list):
    """Find the last object including point in geom_list.
    
    find_object returns (object, array index). If no object includes 
    the given point it returns (geom_list[0], 0).

    """
    cdef int i
    
    i = len(geom_list) - 1
    while i > 0 and geom_list[i].in_object(point) is False:
        i -= 1
        
    return geom_list[i], i


cdef class GeomBox(object):
    """A bounding box of a geometric object.
    
    Attributes:
    low -- the coordinates of the lowest vertex
    high -- the coordinates of the highest vertex
    
    """
    cdef public np.ndarray low, high

    def __init__(self, low=(0,0,0), high=(0,0,0)):
        self.low  = np.array(low, np.double)
        self.high = np.array(high, np.double)
        
    def __getstate__(self):
        d = {}
        d['low'] = self.low
        d['high'] = self.high
        return d

    def __reduce__(self):        
        return self.__class__, (), self.__getstate__()

    def __setstate__(self, d):
        self.low.setfield(d['low'], np.double)
        self.high.setfield(d['high'], np.double)

    def union(self, box):
        """Enlarge the box to include the given box.
        
        """
        self.low.setfield(map(min, self.low, box.low), np.double)
        self.high.setfield(map(max, self.high, box.high), np.double)
        
    def intersection(self, box):
        """Reduce the box to intersect volume with the given box.
        
        """
        self.low.setfield(map(max, self.low, box.low), np.double)
        self.high.setfield(map(min, self.high, box.high), np.double)
        
    def add_point(self, point):
        """Enlarge the box to include the given point.
        
        """
        self.low.setfield(map(min, self.low, point), np.double)
        self.high.setfield(map(max, self.high, point), np.double)
    
    cdef inline bint between(self, double x, double low, double high):
        """Return truth of low <= x <= high.
        
        """
        return low <= x <= high
    
    cpdef bint in_box(self, tuple point):
        """Check whether the given point is in this box.
        
        """
        cdef bint truth

        truth = (self.between(point[0], self.low[0], self.high[0]) and
                 self.between(point[1], self.low[1], self.high[1]) and
                 self.between(point[2], self.low[2], self.high[2]))

        return truth
    
    def overlap(self, box):
        """Check whether the given box intersect with the box.
        
        """
        truth = ((self.between(box.low[0], self.low[0], self.high[0]) or
                  self.between(box.high[0], self.low[0], self.high[0]) or
                  self.between(self.low[0], box.low[0], box.high[0])) and
                 (self.between(box.low[1], self.low[1], self.high[1]) or
                  self.between(box.high[1], self.low[1], self.high[1]) or
                  self.between(self.low[1], box.low[1], box.high[1])) and
                 (self.between(box.low[2], self.low[2], self.high[2]) or
                  self.between(box.high[2], self.low[2], self.high[2]) or
                  self.between(self.low[2], box.low[2], box.high[2])))

        return truth
    
    def divide(self, axis, x):
        high1 = list(self.high)
        high1[axis] = x
        
        low2 = list(self.low)
        low2[axis] = x
        
        return GeomBox(self.low, high1), GeomBox(low2, self.high)
    
    def display_info(self, indent=0):
        print " " * indent, "geom box:",
        print "low:", self.low, "high:", self.high
        
    def __str__(self):
        return "low: " + self.low.__str__() + " high: " + self.high.__str__()
    
    
cdef class GeomBoxNode(object):
    """Node class which makes up a binary search tree.
    
    Attributes:
    box -- a bounding box enclosing the volume of this node
    t1 -- left branch from this node
    t2 -- right branch from this node
    geom_list -- a geometric object list overlapping the volume of this node.
    depth -- depth from the root of this binary search tree
        
    """
    cdef public GeomBox box
    cdef public GeomBoxNode t1, t2
    cdef public tuple geom_list
    cdef public int depth

    def __init__(self, GeomBox box, object geom_list, int depth):
        self.box = box
        self.t1, self.t2 = None, None
        self.geom_list = tuple(geom_list)
        self.depth = depth
        
    def __getstate__(self):
        d = {}
        d['box'] = self.box
        d['t1'] = self.t1
        d['t2'] = self.t2
        d['geom_list'] = self.geom_list
        d['depth'] = self.depth
        return d

    def __reduce__(self):        
        return self.__class__, (None, (), 0), self.__getstate__()

    def __setstate__(self, d):
        self.box = deepcopy(d['box'])
        self.t1 = deepcopy(d['t1'])
        self.t2 = deepcopy(d['t2'])
        self.geom_list = deepcopy(d['geom_list'])
        self.depth = d['depth']


cdef class GeomBoxTree(object):
    """A tree for the fast inclusion test of geometric objects within them.
    
    The tree recursively partitions the unit cell, allowing us to perform 
    binary searches for the object containing a give point.
    
    Attributes:
    root -- root node of the binary search tree
        
    """
    cdef public GeomBoxNode root

    def __init__(self, geom_list):
        box = GeomBox((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf))
        self.root = GeomBoxNode(box, geom_list, 0)
        self.branch_out(self.root)
    
    def __getstate__(self):
        d = {}
        d['root'] = self.root
        return d

    def __reduce__(self):        
        return self.__class__, ((),), self.__getstate__()

    def __setstate__(self, d):
        self.root = deepcopy(d['root'])
        
    def find_best_partition(self, node, divide_axis):
        """
        Find the best place to "cut" along the axis divide_axis in 
        order to maximally divide the objects between the partitions.
        Upon return, n1 and n2 are the number of objects below and 
        above the partition, respectively.
        
        """
        small = 1e-7
        best_partition = None
        
        n1 = n2 = len(node.geom_list)
        
        # Search for the best partition, by checking all possible partitions 
        # either just above the high end of an object or just below the low 
        # end of an object. 
        
        for i in node.geom_list:
            curPartition = i.box.high[divide_axis] + small
            curN1 = curN2 = 0
            for j in node.geom_list:
                if j.box.low[divide_axis] <= curPartition:
                    curN1 += 1
                if j.box.high[divide_axis] >= curPartition:
                    curN2 += 1
            if max(curN1, curN2) < max(n1, n2):
                best_partition = curPartition
                n1 = curN1
                n2 = curN2
                
        for i in node.geom_list:
            curPartition = i.box.low[divide_axis] - small
            curN1 = curN2 = 0
            for j in node.geom_list:
                if j.box.low[divide_axis] <= curPartition:
                    curN1 += 1
                if j.box.high[divide_axis] >= curPartition:
                    curN2 += 1
            if max(curN1, curN2) < max(n1, n2):
                best_partition = curPartition
                n1 = curN1
                n2 = curN2
        
        return best_partition, n1, n2
    
    def divide_geom_box_tree(self, node):
        """Divide box in two, along the axis that maximally partitions the boxes.
        
        """
        # Try partitioning along each dimension, counting the
        # number of objects in the partitioned boxes and finding
        # the best partition.
        best = 0
        division = []
        for i in range(3):
            partition, n1, n2 = self.find_best_partition(node, i)
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
            if box1.overlap(i.box):
                b1GeomList.append(i)
            if box2.overlap(i.box):
                b2GeomList.append(i)
        
        return GeomBoxNode(box1, b1GeomList, node.depth + 1), GeomBoxNode(box2, b2GeomList, node.depth + 1)
    
    def branch_out(self, node):
        node.t1, node.t2 = self.divide_geom_box_tree(node)
        
        if node.t1 or node.t2:
            self.branch_out(node.t1)
            self.branch_out(node.t2)
    
    cdef GeomBoxNode tree_search(self, GeomBoxNode node, tuple point):
        if node.box.in_box(point) == False: 
            return None
        else:
            if not (node.t1 and node.t2):
                return node
            else:
                if node.t1.box.in_box(point):
                    return self.tree_search(node.t1, point)
        
                if node.t2.box.in_box(point):
                    return self.tree_search(node.t2, point)
    
    cpdef tuple object_of_point(self, tuple point):
        cdef GeomBoxNode leaf
        cdef GeometricObject geom_obj
        cdef int idx
        
        leaf = self.tree_search(self.root, point)
        geom_obj, idx = find_object(point, leaf.geom_list)
        
        # eps_inf and mu_inf of compound material
        # refer to the underneath material.
        if isinstance(geom_obj.material, Compound):
            aux_geom_list = leaf.geom_list[:idx]
            underneath_obj, trash = find_object(point, aux_geom_list)
        else:
            underneath_obj = None
            
        return geom_obj, underneath_obj
        
    cpdef tuple material_of_point(self, tuple point):
        cdef GeometricObject geom_obj, underneath_obj
        cdef Material underneath_material

        geom_obj, underneath_obj = self.object_of_point(point)

        if underneath_obj: 
            underneath_material = underneath_obj.material
        else:
            underneath_material = None
            
        return geom_obj.material, underneath_material
        
    def display_info(self, node=None, indent=0):
        if not node: node = self.root
        
        print " " * indent, "depth:", node.depth, node.box
        
        if node.t1 is None or node.t2 is None:
            for i in node.geom_list:
                print " " * (indent + 5), "bounding box:", i.box
                i.display_info(indent + 5)
            
        if node.t1: self.display_info(node.t1, indent + 5)
        if node.t2: self.display_info(node.t2, indent + 5)


####################################################################
#                                                                  #
#                       Geometric primitives                       #
#                                                                  #
####################################################################

cdef class GeometricObject(object):
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
    cdef public Material material
    cdef public GeomBox box

    def __init__(self, material):
        """
        
        Keyword arguments:
        material -- The material that the object is made of. No default.
            
        """ 
        self.material = material

    def __getstate__(self):
        d = {}
        d['material'] = self.material
        d['box'] = self.box
        return d

    def __reduce__(self):
        return self.__class__, (None,), self.__getstate__()

    def __setstate__(self, d):
        self.material = deepcopy(d['material'])
        self.box = deepcopy(d['box'])

    def init(self, space):
        self.material.init(space)
        self.box = self.geom_box()

    def geom_box(self):
        """Return a bounding box enclosing this geometric object.
        
        The derived classes should override this method.

        """        
        raise NotImplementedError
    
    cpdef bint in_object(self, tuple point):
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
        print ' ' * indent, 'geometric object'
        print ' ' * indent, 'center:', self.center
        if self.material:
            self.material.display_info(indent + 5)
       

cdef class DefaultMedium(GeometricObject):
    """A geometric object expanding the whole space.
    
    """
    def __init__(self, material):
        GeometricObject.__init__(self, material)

    cpdef bint in_object(self, tuple point):
        """
        Override GeometricObject.in_object.
        
        """
        
        return self.box.in_box(point)

    def geom_box(self):
        """
        Override GeometriObject.geom_box.
        
        """
        return GeomBox((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf))

    def display_info(self, indent=0):
        """
        Override GeometricObject.display_info.
        
        """
        print ' ' * indent, 'default medium'
        if self.material:
            self.material.display_info(indent + 5)
    
    
cdef class Cone(GeometricObject):
    """Form a cone or possibly a truncated cone. 
    
    Attributes:
    center -- coordinates of the center of this geometric object
    axis -- unit vector of axis
    height -- length of axis
    box -- bounding box
        
    """
    cdef public double radius, radius2, height
    cdef public np.ndarray center, axis

    def __init__(self, material, object center=(0,0,0), double radius2=0, object axis=(1,0,0), double radius=1, double height=1):
        """
        
        Keyword arguments:
        radius2 -- Radius of the tip of the cone (i.e. the end of the 
            cone pointed to by the axis vector). Defaults to zero 
            (a "sharp" cone).
            
        """
        GeometricObject.__init__(self, material)

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

        self.center = np.array(center, np.double)
        self.axis = np.array(axis, np.double) / norm(axis)
        self.height = float(height)

    def __getstate__(self):
        d = GeometricObject.__getstate__(self)
        d['radius'] = self.radius
        d['radius2'] = self.radius2
        d['height'] = self.height
        d['center'] = self.center
        d['axis'] = self.axis
        return d

    def __setstate__(self, d):
        GeometricObject.__setstate__(self, d)
        self.radius = d['radius']
        self.radius2 = d['radius2']
        self.height = d['height']
        self.center.setfield(d['center'], np.double)
        self.axis.setfield(d['axis'], np.double)

    cpdef bint in_object(self, tuple point):
        """Check whether the given point is in this Cone.
        
        """
        cdef np.ndarray p, r
        cdef double proj, radius
        cdef bint truth

        p = np.array(point, np.double)
        r = p - self.center
        proj = np.dot(self.axis, r)
        if abs(proj) <= .5 * self.height:
            if self.radius2 == self.radius == np.inf:
                truth = True
            else:
                radius = self.radius
                radius += (proj / self.height + .5) * (self.radius2 - radius)
                truth = norm(r - proj * self.axis) <= abs(radius)
        else:
            truth = False

        return truth

    def display_info(self, indent=0):
        """
        Override GeometricObject.display_info.
        
        """
        print ' ' * indent, 'cone'
        print ' ' * indent,
        print 'center:', self.center,
        print 'radius:', self.radius,
        print 'height:' , self.height,
        print 'axis:', self.axis,
        print 'radius2:', self.radius2
        if self.material:
            self.material.display_info(indent + 5)
        
    def geom_box(self):
        """
        Override GeometricObject.geom_box.
        
        """
        h = .5 * self.height
        
        # radial is a vector consisting of Cartesian unit vectors
        # which are orthogonal to the self.axis.
        r = np.ones((3,), np.double)
        radial = r - self.axis * r

        # bounding box for -h*axis cylinder end
        tmpBox1 = GeomBox(low=self.center, high=self.center)
        tmpBox1.low -= h * self.axis + self.radius * radial
        tmpBox1.high -= h * self.axis - self.radius * radial
        
        # bounding box for +h*axis cylinder end
        tmpBox2 = GeomBox(low=self.center, high=self.center)
        tmpBox2.low += h * self.axis - self.radius2 * radial
        tmpBox2.high += h * self.axis + self.radius2 * radial
        
        tmpBox1.union(tmpBox2)
        
        return tmpBox1
       
    
class Cylinder(Cone):
    """Form a cylinder.
    
    """	
    def __init__(self, material, center=(0, 0, 0), axis=(0, 0, 1),
                 radius=1, height=1):
        """
        Keyword arguments:
            material -- The material that the object is made of. 
                No default.
            center -- Center point of the object. Default is (0,0,0). 
            axis -- Direction of the cylinder's axis; the length of 
                this vector is ignored. Defaults to point parallel to 
                the z axis i.e., (0,0,1).
            radius -- Radius of the cylinder's cross-section. Default is 1.
            height -- Length of the cylinder along its axis. Default is 1. 
        
        """
        Cone.__init__(self, material, center, radius, axis, radius, height)

    def display_info(self, indent=0):
        """Display information of this cylinder.
        
        """
        print ' ' * indent, 'cylinder'
        print ' ' * indent,
        print 'center:', self.center,
        print 'radius:', self.radius,
        print 'height:', self.height,
        print 'axis:', self.axis
        if self.material:
            self.material.display_info(indent + 5)


cdef class Block(GeometricObject):
    """Form a parallelpiped(i.e., a brick, possibly with non-orthogonal axes.)
    
    """
    cdef public np.ndarray center, e1, e2, e3, size, projection_matrix

    def __init__(self, material, center=(0, 0, 0), 
                 e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1), 
                 size=(1, 1, 1)):
        """
        Keyword arguments:
            center -- center location. Default is (0, 0, 0).
            e1, e2, e3 -- The directions of the axes of the block; the 
                lengths of these vectors are ignored. Must be linearly 
                independent. They default to the three Cartesian axis.
            size -- The lengths of the block edges along each of its 
                three axes. Default is (1, 1, 1).
            
        """
        GeometricObject.__init__(self, material)

        self.center = np.array(center, np.double)
        
        self.e1 = np.array(e1, np.double) / norm(e1)
        self.e2 = np.array(e2, np.double) / norm(e2)
        self.e3 = np.array(e3, np.double) / norm(e3)
        self.size = np.array(size, np.double)
        
        self.projection_matrix = np.array([self.e1, self.e2, self.e3])

    def __getstate__(self):
        d = GeometricObject.__getstate__(self)
        d['center'] = self.center
        d['e1'] = self.e1
        d['e2'] = self.e2
        d['e3'] = self.e3
        d['size'] = self.size
        d['pm'] = self.projection_matrix
        return d

    def __setstate__(self, d):
        GeometricObject.__setstate__(self, d)
        self.center.setfield(d['center'], np.double)
        self.e1.setfield(d['e1'], np.double)
        self.e2.setfield(d['e2'], np.double)
        self.e3.setfield(d['e3'], np.double)
        self.size.setfield(d['size'], np.double)
        self.projection_matrix.setfield(d['pm'], np.double)
        
    cpdef bint in_object(self, tuple point):
        """Check whether the given point is in this block.

        """
        cdef np.ndarray p, r, proj
        cdef bint truth

        p = np.array(point, np.double)
        r = p - self.center
        proj = np.dot(self.projection_matrix, r)
        truth = (np.abs(proj) <= .5 * self.size).all()

        return truth
        
    def geom_box(self):
        """Return a GeomBox for this block.

        """
        tmpBox = GeomBox(low=self.center, high=self.center)
        # enlarge the box to be big enough to contain all 8 corners
        # of the block.
        s1 = self.size[0] * self.e1
        s2 = self.size[1] * self.e2
        s3 = self.size[2] * self.e3
        
        corner = self.center - 0.5 * (s1 + s2 + s3)
        
        tmpBox.add_point(corner)
        tmpBox.add_point(corner + s1)
        tmpBox.add_point(corner + s2)
        tmpBox.add_point(corner + s3)
        tmpBox.add_point(corner + s1 + s2)
        tmpBox.add_point(corner + s1 + s3)
        tmpBox.add_point(corner + s3 + s2)
        tmpBox.add_point(corner + s1 + s2 + s3)

        return tmpBox
        
    def display_info(self, indent=0):
        """Display information of this block.

        """
        print ' ' * indent, 'block'
        print ' ' * indent,
        print 'center:', self.center,
        print 'size:', self.size,
        print 'axes:', self.e1, self.e2, self.e3
        if self.material:
            self.material.display_info(indent + 5)


cdef class Ellipsoid(Block):
    """Form an ellipsoid.
    
    """
    cdef public np.ndarray inverse_semi_axes

    def __init__(self, material, center=(0, 0, 0),
                 e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1),
                 size=(1, 1, 1)):
        """

        """
        Block.__init__(self, material, center, e1, e2, e3, size)
        self.inverse_semi_axes = 2 / np.array(size, np.double)

    def __getstate__(self):
        d = Block.__getstate__(self)
        d['isa'] = self.inverse_semi_axes
        return d

    def __setstate__(self, d):
        Block.__setstate__(self, d)
        self.inverse_semi_axes = d['isa']

    cpdef bint in_object(self, tuple point):
        """Check whether the given point is in this ellipsoid.

        """
        cdef np.ndarray p, r, q, proj
        cdef bint truth

        p = np.array(point, np.double)
        r = p - self.center
        proj = np.dot(self.projection_matrix, r)
        q = proj * self.inverse_semi_axes
        truth = sum(q * q) <= 1

        return truth

    def display_info(self, indent=0):
        """Display information of this ellipsoid.

        """
        print ' ' * indent, 'ellipsoid'
        print ' ' * indent,
        print 'center:', self.center,
        print 'size:', self.size,
        print 'axis:', self.e1, self.e2, self.e3
        if self.material:
            self.material.display_info(indent + 5)


cdef class Sphere(GeometricObject):
    """Form a sphere.
    
    Attributes:
    radius -- Radius of the sphere.
    
    """
    cdef public double radius
    cdef public np.ndarray center

    def __init__(self, material, center=(0, 0, 0), radius=1):
        """
        
        Keyword arguments:
        radius -- Radius of the sphere. Default is 1. 
        material -- The material that the object is made of. 
            No default.
        center -- Center point of the object. Default is (0,0,0).

        """
        GeometricObject.__init__(self, material)

        if radius < 0:
            msg = "radius must be non-negative."
            raise ValueError(msg)
        else:
            self.radius = float(radius)

        self.center = np.array(center, np.double)


    def __getstate__(self):
        d = GeometricObject.__getstate__(self)
        d['radius'] = self.radius
        d['center'] = self.center
        return d

    def __setstate__(self, d):
        GeometricObject.__setstate__(self, d)
        self.radius = d['radius']
        self.center.setfield(d['center'], np.double)

    def geom_box(self):
        """Return GeomBox for the sphere.

        """
        box = GeomBox(low=self.center, high=self.center)

        box.low -= self.radius
        box.high += self.radius

        return box

    cpdef bint in_object(self, tuple point):
        """Check whether the given point is in the sphere.

        """
        cdef np.ndarray p, r
        cdef bint truth

        p = np.array(point, np.double)
        r = p - self.center
        truth = norm(r) <= self.radius;

        return truth

    def display_info(self, indent=0):
        """Display information of the sphere.

        """
        print ' ' * indent, 'sphere'
        print ' ' * indent,
        print 'center:', self.center,
        print 'radius:', self.radius
        if self.material:
            self.material.display_info(indent + 5)


cdef class Shell(GeometricObject):
    """Form a boundary.
     
    """
    cdef public double d
    cdef public np.ndarray center, half_size
    cdef public list box_list
    cdef public bint minus_x, plus_x, minus_y, plus_y, minus_z, plus_z
    cdef public bint boundary

    def __init__(self, material, center=(0,0,0), size=None,
                 thickness=1,
                 plus_x=True, minus_x=True,
                 plus_y=True, minus_y=True,
                 plus_z=True, minus_z=True):
        """
        Keyword arguments:
            material -- The filling material
            center -- coordinates of the center of this geometric object. Default is (0,0,0).
            size -- size of ot the shell. Default is None.
            thickness -- The spatial thickness of the Shell (the 
                distance between inner and outer surface)
                Default value is 1.
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
        GeometricObject.__init__(self, material)
        
        self.center = np.array(center, np.double)
        self.d = float(thickness)
        
        self.box_list = []
        
        self.minus_x, self.plus_x = minus_x, plus_x
        self.minus_y, self.plus_y = minus_y, plus_y
        self.minus_z, self.plus_z = minus_z, plus_z
        
        if size is None:
            self.half_size = np.zeros((3,), np.double)
            self.boundary = True
        else:
            self.half_size = np.array(map(lambda i: 0.5 * i, size),
                                      np.double)
            self.boundary = False
            
        # do someting for the PML derived class?
        if isinstance(material, Pml):
            pass
        
    def __getstate__(self):
        d = GeometricObject.__getstate__(self)
        d['center'] = self.center
        d['half_size'] = self.half_size
        d['d'] = self.d
        d['box_list'] = self.box_list
        d['minus_x'] = self.minus_x
        d['plus_x'] = self.plus_x
        d['minus_y'] = self.minus_y
        d['plus_y'] = self.plus_y
        d['minus_z'] = self.minus_z
        d['plus_z'] = self.plus_z
        d['boundary'] = self.boundary
        return d
       
    def __setstate__(self, d):
        GeometricObject.__setstate__(self, d)
        self.center.setfield(d['center'], np.double)
        self.half_size.setfield(d['half_size'], np.double)
        self.d = d['d']
        self.box_list = deepcopy(d['box_list'])
        self.minus_x = d['minus_x']
        self.plus_x = d['plus_x']
        self.minus_y = d['minus_y']
        self.plus_y = d['plus_y']
        self.minus_z = d['minus_z']
        self.plus_z = d['plus_z']
        self.boundary = d['boundary']

    def init(self, space):
        if self.boundary:
            self.half_size.setfield(space.half_size, np.double)

        for i in range(3):
            if 2 * self.half_size[i] < space.dr[i]:
                self.half_size[i] = 0.5 * space.dr[i]

        if 2 * self.half_size[0] <= space.dr[0]:
            self.plus_x = False
            self.minus_x = False

        if 2 * self.half_size[1] <= space.dr[1]:
            self.plus_y = False
            self.minus_y = False

        if 2 * self.half_size[2] <= space.dr[2]:
            self.plus_z = False
            self.minus_z = False

        if self.plus_x:
            low = (self.center[0] + self.half_size[0] - self.d,
                   self.center[1] - self.half_size[1],
                   self.center[2] - self.half_size[2])
            high = self.center + self.half_size
            self.box_list.append(GeomBox(low, high))

        if self.minus_x:
            low = self.center - self.half_size
            high = (self.center[0] - self.half_size[0] + self.d, 
                    self.center[1] + self.half_size[1],
                    self.center[2] + self.half_size[2])
            self.box_list.append(GeomBox(low, high))

        if self.plus_y:
            low = (self.center[0] - self.half_size[0], 
                   self.center[1] + self.half_size[1] - self.d,
                   self.center[2] - self.half_size[2])
            high = self.center + self.half_size
            self.box_list.append(GeomBox(low, high))

        if self.minus_y:
            low = self.center - self.half_size
            high = (self.center[0] + self.half_size[0], 
                    self.center[1] - self.half_size[1] + self.d,
                    self.center[2] + self.half_size[2])
            self.box_list.append(GeomBox(low, high))

        if self.plus_z:
            low = (self.center[0] - self.half_size[0], 
                   self.center[1] - self.half_size[1],
                   self.center[2] + self.half_size[2] - self.d)
            high = self.center + self.half_size
            self.box_list.append(GeomBox(low, high))

        if self.minus_z:
            low = self.center - self.half_size
            high = (self.center[0] + self.half_size[0],
                    self.center[1] + self.half_size[1],
                    self.center[2] - self.half_size[2] + self.d)
            self.box_list.append(GeomBox(low, high))
        
        self.material.init(space, 
                           (self.center, self.half_size, self.d))
        self.box = self.geom_box()

    cpdef bint in_object(self, tuple point):
        cdef GeomBox box

        for box in self.box_list:
            if box.in_box(point):
                return True
        return False
        
    def geom_box(self):
        return GeomBox(-self.half_size, self.half_size)
        
    def display_info(self, indent=0):
        print ' ' * indent, 'shell'
        print ' ' * indent, 'center:', self.center
        print ' ' * indent, 'size:', 2 * self.half_size
        print ' ' * indent,
        print '+x:', self.plus_x, '-x:', self.minus_x,
        print '+y:', self.plus_y, '-y:', self.minus_y,
        print '+z:', self.plus_z, '-z:', self.minus_z
        if self.material:
            self.material.display_info(indent + 5)
