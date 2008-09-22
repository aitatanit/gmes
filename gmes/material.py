#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *

from pointwise_material import *
import constants as const


class Material(object):
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "material type object"

    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            underneath -- underneath material object of the target point.
        """
        raise NotImplementedError
    
    def get_pointwise_mateiral_ey(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            underneath -- underneath material object of the target point.
        """
        raise NotImplementedError
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            underneath -- underneath material object of the target point.
        """
        raise NotImplementedError
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            underneath -- underneath material object of the target point.
        """
        raise NotImplementedError
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            underneath -- underneath material object of the target point.
        """
        raise NotImplementedError
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            underneath -- underneath material object of the target point.
        """
        raise NotImplementedError
    
    
class Dummy(Material):
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "dummy object"
        
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        pw_obj = DummyEx(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_mateiral_ey(self, idx, coords, underneath=None):
        pw_obj = DummyEy(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        pw_obj = DummyEz(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        pw_obj = DummyHx(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        pw_obj = DummyHy(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        pw_obj = DummyHz(idx, self.mu_r)
        return pw_obj
    
    
class Dielectric(Material):
    def __init__(self, epsilon_r=1, mu_r=1):
        """Representation of dielectric medium.
        
        Arguments:
            epsilon_r -- relative permittivity
            mu_r -- relative permeability
        
        """
        self.epsilon_r = epsilon_r
        self.mu_r = mu_r

    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "dielectric"
        print " " * indent, 
        print "relative permittivity:", self.epsilon_r,
        print "relative permeability:", self.mu_r

    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        pw_obj = DielectricEx(idx, self.epsilon_r)
        return pw_obj

    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        pw_obj = DielectricEy(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        pw_obj = DielectricEz(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        pw_obj = DielectricHx(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        pw_obj = DielectricHy(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        pw_obj = DielectricHz(idx, self.mu_r)
        return pw_obj

class Compound(object):
    """Represent the compound material.
    
    Compound material means that its electric permittivity and 
    magnetic permeability refers to the underneath material. 
    
    """
    pass

class PML(Material, Compound):
    """Base class of PML materials.
    
    Attributes:
        d -- thickness of PML medium
        half_size --
        dt -- time differential
        dw -- tuple of space differentials
        sigma_opt --
        initialized -- initialization semaphore 
        
    """
    def init(self, thickness, space):
        """
        The thickness of PML layer is provided from the boundary instance 
        which contain the pml. Also, the differential of space and time 
        should get from the space instance.
        
        """
        self.d = float(thickness)
        
        half_size = []
        for i in space.half_size:
            if i <= self.d: i = inf
            half_size.append(i)
        self.half_size = array(half_size, float)
        
        self.dt = space.dt
        self.dw = array((space.dx, space.dy, space.dz), float)
        self.sigma_max = self.sigma_max_ratio * self.get_sigma_opt()
        
        self.initialized = True
        
    def get_sigma_opt(self):
        """Calculate the optimal value of sigma.
        
        """
        eta = sqrt(self.effective_mu_r / self.effective_epsilon_r) * const.Z0
        return .8 * (self.m + 1) / (eta * self.dw)
    
    def sigma(self, w, component):
        """Polynomial grading of sigma.
        
        """
        if w <= self.d - self.half_size[component]:
            return self.sigma_max[component] * (1 - (self.half_size[component] + w) / self.d)**self.m
        elif self.half_size[component] - self.d <= w:
            return self.sigma_max[component] * (1 - (self.half_size[component] - w) / self.d)**self.m
        else:
            return 0.0
        
    def kappa(self, w, component):
        """Polynomial grading of kappa.
        
        """
        if w <= self.d - self.half_size[component]:
            return 1 + (self.kappa_max - 1) * (1 - (self.half_size[component] + w) / self.d)**self.m
        elif self.half_size[component] - self.d <= w:
            return 1 + (self.kappa_max - 1) * (1 - (self.half_size[component] - w) / self.d)**self.m
        else:
            return 1.0
        
        
class UPML(PML):
    """Form Uniaxial Perfectly Matched Layer (UPML).
    
    This class implements Uniaxial PML represented in 'A. Taflove and S. C. Hagness, 
    Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
    3rd ed., Artech House, Inc., 2005.
    
    Attributes:
        m --
        kappa_max --
        effective_epsilon_r -- the effective relative permittivity of incident mode impinging on the PML boundary
        effective_mu_r -- the effective relative permeability of incident mode impinging on the PML boundary
        sigma_max_ratio --

    """    
    def __init__(self, effective_epsilon_r=1, effective_mu_r=1, m=3.5, kappa_max=1, sigma_max_ratio=.75):
        self.initialized = False
        
        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.effective_epsilon_r = float(effective_epsilon_r)
        self.effective_mu_r = float(effective_epsilon_r)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
    def display_info(self, indent=0):
        """Display the parameter values.

        Override PML.display_info.
        
        """
        print " " * indent, "UPML"
        print " " * indent, 
        print "relative effective permittivity:", self.effective_epsilon_r,
        print "relative effective permeability:", self.effective_mu_r,
        print "sigma_max:", self.sigma_max,
        print "m:", self.m,
        print "kappa_max:", self.kappa_max,
        
    def c1(self, w, component):
        numerator = 2 * const.epsilon0 * self.kappa(w, component) - self.sigma(w, component) * self.dt
        denominator = 2 * const.epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c2(self, w, component):
        numerator = 2 * const.epsilon0 * self.dt
        denominator = 2 * const.epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c3(self, w, component):
        numerator = 2 * const.epsilon0 * self.kappa(w, component) - self.sigma(w, component) * self.dt
        denominator = 2 * const.epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c4(self, w, component):
        denominator = 2 * const.epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return 1 / denominator
    
    def c5(self, w, component):
        numerator = 2 * const.epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator
    
    def c6(self, w, component):
        numerator = 2 * const.epsilon0 * self.kappa(w, component) - self.sigma(w, component) * self.dt
        return numerator
    
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        c1 = self.c1(coords[1], 1)
        c2 = self.c2(coords[1], 1)
        c3 = self.c3(coords[2], 2)
        c4 = self.c4(coords[2], 2)
        c5 = self.c5(coords[0], 0)
        c6 = self.c6(coords[0], 0)
        pw_obj = UPMLEx(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        c1 = self.c1(coords[2], 2)
        c2 = self.c2(coords[2], 2)
        c3 = self.c3(coords[0], 0)
        c4 = self.c4(coords[0], 0)
        c5 = self.c5(coords[1], 1)
        c6 = self.c6(coords[1], 1)
        pw_obj = UPMLEy(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        c1 = self.c1(coords[0], 0)
        c2 = self.c2(coords[0], 0)
        c3 = self.c3(coords[1], 1)
        c4 = self.c4(coords[1], 1)
        c5 = self.c5(coords[2], 2)
        c6 = self.c6(coords[2], 2)
        pw_obj = UPMLEz(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        c1 = self.c1(coords[1], 1)
        c2 = self.c2(coords[1], 1)
        c3 = self.c3(coords[2], 2)
        c4 = self.c4(coords[2], 2)
        c5 = self.c5(coords[0], 0)
        c6 = self.c6(coords[0], 0)
        pw_obj = UPMLHx(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        c1 = self.c1(coords[2], 2)
        c2 = self.c2(coords[2], 2)
        c3 = self.c3(coords[0], 0)
        c4 = self.c4(coords[0], 0)
        c5 = self.c5(coords[1], 1)
        c6 = self.c6(coords[1], 1)
        pw_obj = UPMLHy(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        c1 = self.c1(coords[0], 0)
        c2 = self.c2(coords[0], 0)
        c3 = self.c3(coords[1], 1)
        c4 = self.c4(coords[1], 1)
        c5 = self.c5(coords[2], 2)
        c6 = self.c6(coords[2], 2)
        pw_obj = UPMLHz(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    
class CPML(PML):
    """Form Complex Frequency Shifted (CFS) Perfectly Matched Layer (PML).
    
    This class implements CFS PML represented in 'A. Taflove and S. C. Hagness, 
    Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
    3rd ed., Artech House, Inc., 2005.
    
    Attributes:
        m --
        kappa_max --
        a_max --
        effective_epsilon_r -- the effective relative permittivity of incident mode impinging on the PML boundary
        effective_mu_r -- the effective relative permeability of incident mode impinging on the PML boundary
        sigma_max_ratio --
    
    """
    def __init__(self, effective_epsilon_r=1, effective_mu_r=1, m=3, kappa_max=15, m_a=1, a_max=0.2, sigma_max_ratio=1):
        self.initialized = False
        
        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.m_a = float(m_a)
        self.a_max = float(a_max)
        self.effective_epsilon_r = float(effective_epsilon_r)
        self.effective_mu_r = float(effective_mu_r)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
    def display_info(self, indent=0):
        """Display the parameter values.

        Override PML.display_info.

        """
        print " " * indent, "CPML"
        print " " * indent, 
        print "effective relative permittivity:", self.effective_epsilon_r,
        print "effective relative permeability:", self.effective_mu_r,
        
        print " " * indent,
        print "sigma_max:", self.sigma_max,
        print "m:", self.m,
        print "kappa_max:", self.kappa_max,
        print "m_a:", self.m_a,
        print "a_max:", self.a_max

    def a(self, w, component):
        if w <= self.d - self.half_size[component]:
            return self.a_max * ((self.half_size[component] + w) / self.d)**self.m_a
        elif self.half_size[component] - self.d <= w:
            return self.a_max * ((self.half_size[component] - w) / self.d)**self.m_a
        else:
            return 0.0
       
    def b(self, w, component):
        exponent = -(self.sigma(w, component) / self.kappa(w, component) + 
                     self.a(w, component)) * self.dt / const.epsilon0
        return exp(exponent)
        
    def c(self, w, component):
        numerator = self.sigma(w, component) * (self.b(w, component) - 1)
        if numerator == 0:
            return 0.0
        denominator = (self.sigma(w, component) + 
                       self.kappa(w, component) * self.a(w, component)) * self.kappa(w, component)
        return numerator / denominator
    
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        by = self.b(coords[1], 1)
        bz = self.b(coords[2], 2)
        cy = self.c(coords[1], 1)
        cz = self.c(coords[2], 2)
        kappay = self.kappa(coords[1], 1)
        kappaz = self.kappa(coords[2], 2)
        pw_obj = CPMLEx(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        bz = self.b(coords[2], 2)
        bx = self.b(coords[0], 0)
        cz = self.c(coords[2], 2)
        cx = self.c(coords[0], 0)
        kappaz = self.kappa(coords[2], 2)
        kappax = self.kappa(coords[0], 0)
        pw_obj = CPMLEy(idx, underneath.epsilon_r, bz, bx, cz, cx, kappaz, kappax)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        bx = self.b(coords[0], 0)
        by = self.b(coords[1], 1)
        cx = self.c(coords[0], 0)
        cy = self.c(coords[1], 1)
        kappax = self.kappa(coords[0], 0)
        kappay = self.kappa(coords[1], 1)
        pw_obj = CPMLEz(idx, underneath.epsilon_r, bx, by, cx, cy, kappax, kappay)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        by = self.b(coords[1], 1)
        bz = self.b(coords[2], 2)
        cy = self.c(coords[1], 1)
        cz = self.c(coords[2], 2)
        kappay = self.kappa(coords[1], 1)
        kappaz = self.kappa(coords[2], 2)
        pw_obj = CPMLHx(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        bz = self.b(coords[2], 2)
        bx = self.b(coords[0], 0)
        cz = self.c(coords[2], 2)
        cx = self.c(coords[0], 0)
        kappaz = self.kappa(coords[2], 2)
        kappax = self.kappa(coords[0], 0)
        pw_obj = CPMLHy(idx, underneath.mu_r, bz, bx, cz, cx, kappaz, kappax)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        bx = self.b(coords[0], 0)
        by = self.b(coords[1], 1)
        cx = self.c(coords[0], 0)
        cy = self.c(coords[1], 1)
        kappax = self.kappa(coords[0], 0)
        kappay = self.kappa(coords[1], 1)
        pw_obj = CPMLHz(idx, underneath.mu_r, bx, by, cx, cy, kappax, kappay)
        return pw_obj
        
