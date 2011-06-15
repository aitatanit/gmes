#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')

from numpy import *

from pw_material import *
import constants as const


class Material(object):
    """A super class for material types.
    
    """
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "material type object"

    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def init(self, space, param=None):
        pass
        
        
class Dummy(Material):
    """A dummy material type which dosen't update the field component.
    
    """
    def __init__(self, epsilon=1, mu=1):
        self.epsilon = float(epsilon)
        self.mu = float(mu)
        
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "dummy object"
        print " " * indent, 
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyExCmplx(self.epsilon)
        else:
            pw_obj = DummyExReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyEyCmplx(self.epsilon)
        else:
            pw_obj = DummyEyReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyEzCmplx(self.epsilon)
        else:
            pw_obj = DummyEzReal(self.epsilon)
        
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyHxCmplx(self.mu)
        else:
            pw_obj = DummyHxReal(self.mu)
        
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyHyCmplx(self.mu)
        else:
            pw_obj = DummyHyReal(self.mu)
        
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyHzCmplx(self.mu)
        else:
            pw_obj = DummyHzReal(self.mu)
        
        return pw_obj

    
class Zero(Material):
    """A material type which sets the field value zero.
    
    """
    def __init__(self, epsilon=1, mu=1):
        self.epsilon = float(epsilon)
        self.mu = float(mu)
        
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "zero object"
        print " " * indent, 
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ZeroExCmplx(self.epsilon)
        else:
            pw_obj = ZeroExReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ZeroEyCmplx(self.epsilon)
        else:
            pw_obj = ZeroEyReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ZeroEzCmplx(self.epsilon)
        else:
            pw_obj = ZeroEzReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ZeroHxCmplx(self.mu)
        else:
            pw_obj = ZeroHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ZeroHyCmplx(self.mu)
        else:
            pw_obj = ZeroHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ZeroHzCmplx(self.mu)
        else:
            pw_obj = ZeroHzReal(self.mu)
            
        return pw_obj
    

class One(Material):
    """A material type which sets the field value one.
    
    """
    def __init__(self, epsilon=1, mu=1):
        self.epsilon = float(epsilon)
        self.mu = float(mu)
        
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "one object"
        print " " * indent, 
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = OneExCmplx(self.epsilon)
        else:
            pw_obj = OneExReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = OneEyCmplx(self.epsilon)
        else:
            pw_obj = OneEyReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = OneEzCmplx(self.epsilon)
        else:
            pw_obj = OneEzReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = OneHxCmplx(self.mu)
        else:
            pw_obj = OneHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = OneHyCmplx(self.mu)
        else:
            pw_obj = OneHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = OneHzCmplx(self.mu)
        else:
            pw_obj = OneHzReal(self.mu)
            
        return pw_obj
    
            
class Dielectric(Material):
    def __init__(self, epsilon=1, mu=1):
        """Representation of dielectric medium.
        
        Arguments:
            epsilon -- permittivity
            mu -- permeability
        
        """
        self.epsilon = float(epsilon)
        self.mu = float(mu)

    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "dielectric"
        print " " * indent, 
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu

    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricExCmplx(self.epsilon)
        else:
            pw_obj = DielectricExReal(self.epsilon)
            
        return pw_obj

    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricEyCmplx(self.epsilon)
        else:
            pw_obj = DielectricEyReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricEzCmplx(self.epsilon)
        else:
            pw_obj = DielectricEzReal(self.epsilon)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricHxCmplx(self.mu)
        else:
            pw_obj = DielectricHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricHyCmplx(self.mu)
        else:
            pw_obj = DielectricHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricHzCmplx(self.mu)
        else:
            pw_obj = DielectricHzReal(self.mu)
            
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
        half_size -- a tuple that represents the half of the outer volume
        dt -- time differential
        dw -- tuple of space differentials
        sigma_opt -- optimal conductivity
        initialized -- initialization semaphore 
        
    """
    def init(self, space, thickness):
        """
        The thickness of PML layer is provided from the boundary instance 
        which contain the PML. Also, the differential of space and time 
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
        """Calculate the optimal value of conductivity.
        
        """
        eta = sqrt(self.effective_mu / self.effective_epsilon)
        return .8 * (self.m + 1) / (eta * self.dw)
    
    def sigma(self, w, component):
        """Polynomial grading of conductivity.
        
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
    
    This class implements CFS PML represented in
    S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,"
    Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
    Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
    2005, pp. 273-328.
    
    Attributes:
        m -- 
        kappa_max -- maximum of kappa
        effective_epsilon -- the effective permittivity of incident mode impinging on the PML boundary
        effective_mu -- the effective permeability of incident mode impinging on the PML boundary
        sigma_max_ratio -- the ratio between sigma_max and sigma_opt

    """    
    def __init__(self, effective_epsilon=1, effective_mu=1, m=3.5, kappa_max=1, sigma_max_ratio=.75):
        self.initialized = False
        
        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.effective_epsilon = float(effective_epsilon)
        self.effective_mu = float(effective_epsilon)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
        self.epsilon = self.effective_epsilon
        self.mu = self.effective_mu
        
    def display_info(self, indent=0):
        """Display the parameter values.

        Override PML.display_info.
        
        """
        print " " * indent, "UPML"
        print " " * indent, 
        print "effective permittivity:", self.effective_epsilon,
        print "effective permeability:", self.effective_mu,
        print "sigma_max:", self.sigma_max,
        print "m:", self.m,
        print "kappa_max:", self.kappa_max,
        
    def c1(self, w, component):
        numerator = 2 * self.kappa(w, component) \
            - self.sigma(w, component) * self.dt
        denominator = 2 * self.kappa(w, component) \
            + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c2(self, w, component):
        numerator = 2 * self.dt
        denominator = 2 * self.kappa(w, component) \
            + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c3(self, w, component):
        numerator = 2 * self.kappa(w, component) \
            - self.sigma(w, component) * self.dt
        denominator = 2 * self.kappa(w, component) \
            + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c4(self, w, component):
        denominator = 2 * self.kappa(w, component) \
            + self.sigma(w, component) * self.dt
        return 1 / denominator
    
    def c5(self, w, component):
        numerator = 2 * self.kappa(w, component) \
            + self.sigma(w, component) * self.dt
        return numerator
    
    def c6(self, w, component):
        numerator = 2 * self.kappa(w, component) \
            - self.sigma(w, component) * self.dt
        return numerator
    
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        c1 = self.c1(coords[1], 1)
        c2 = self.c2(coords[1], 1)
        c3 = self.c3(coords[2], 2)
        c4 = self.c4(coords[2], 2)
        c5 = self.c5(coords[0], 0)
        c6 = self.c6(coords[0], 0)
        
        if cmplx:
            pw_obj = UpmlExCmplx(underneath.epsilon, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlExReal(underneath.epsilon, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        c1 = self.c1(coords[2], 2)
        c2 = self.c2(coords[2], 2)
        c3 = self.c3(coords[0], 0)
        c4 = self.c4(coords[0], 0)
        c5 = self.c5(coords[1], 1)
        c6 = self.c6(coords[1], 1)
        
        if cmplx:
            pw_obj = UpmlEyCmplx(underneath.epsilon, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlEyReal(underneath.epsilon, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        c1 = self.c1(coords[0], 0)
        c2 = self.c2(coords[0], 0)
        c3 = self.c3(coords[1], 1)
        c4 = self.c4(coords[1], 1)
        c5 = self.c5(coords[2], 2)
        c6 = self.c6(coords[2], 2)
        
        if cmplx:
            pw_obj = UpmlEzCmplx(underneath.epsilon, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlEzReal(underneath.epsilon, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        c1 = self.c1(coords[1], 1)
        c2 = self.c2(coords[1], 1)
        c3 = self.c3(coords[2], 2)
        c4 = self.c4(coords[2], 2)
        c5 = self.c5(coords[0], 0)
        c6 = self.c6(coords[0], 0)
        
        if cmplx:
            pw_obj = UpmlHxCmplx(underneath.mu, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlHxReal(underneath.mu, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        c1 = self.c1(coords[2], 2)
        c2 = self.c2(coords[2], 2)
        c3 = self.c3(coords[0], 0)
        c4 = self.c4(coords[0], 0)
        c5 = self.c5(coords[1], 1)
        c6 = self.c6(coords[1], 1)
        
        if cmplx:
            pw_obj = UpmlHyCmplx(underneath.mu, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlHyReal(underneath.mu, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        c1 = self.c1(coords[0], 0)
        c2 = self.c2(coords[0], 0)
        c3 = self.c3(coords[1], 1)
        c4 = self.c4(coords[1], 1)
        c5 = self.c5(coords[2], 2)
        c6 = self.c6(coords[2], 2)
        
        if cmplx:
            pw_obj = UpmlHzCmplx(underneath.mu, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlHzReal(underneath.mu, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    
class CPML(PML):
    """Form Complex Frequency Shifted (CFS) Perfectly Matched Layer (PML).
    
    This class implements CFS PML represented in
    S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,
    Computational Electrodynamics: The Finite-Difference Time-Domain Method," 
    Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
    2005, pp. 273-328.
    
    Attributes:
        m -- default 3
        kappa_max -- default 15
        a_max -- default 0. CPML works like UPML when a_max = 0.
        effective_epsilon -- the effective permittivity of incident mode impinging on the PML boundary
        effective_mu -- the effective permeability of incident mode impinging on the PML boundary
        sigma_max_ratio -- default 1
    
    """
    def __init__(self, effective_epsilon=1, effective_mu=1, m=3, kappa_max=2, m_a=1, a_max=0, sigma_max_ratio=2):
        self.initialized = False
        
        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.m_a = float(m_a)
        self.a_max = float(a_max)
        self.effective_epsilon = float(effective_epsilon)
        self.effective_mu = float(effective_mu)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
        self.epsilon = self.effective_epsilon
        self.mu = self.effective_mu
        
    def display_info(self, indent=0):
        """Display the parameter values.

        Override PML.display_info.

        """
        print " " * indent, "CPML"
        print " " * indent, 
        print "effective permittivity:", self.effective_epsilon,
        print "effective permeability:", self.effective_mu
        
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
        exponent = -(self.sigma(w, component) / self.kappa(w, component) + self.a(w, component)) * self.dt
        return exp(exponent)
        
    def c(self, w, component):
        numerator = self.sigma(w, component) * (self.b(w, component) - 1)
        if numerator == 0:
            return 0.0
        denominator = (self.sigma(w, component) + 
                       self.kappa(w, component) * self.a(w, component)) * self.kappa(w, component)
        return numerator / denominator
    
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        by = self.b(coords[1], 1)
        bz = self.b(coords[2], 2)
        cy = self.c(coords[1], 1)
        cz = self.c(coords[2], 2)
        kappay = self.kappa(coords[1], 1)
        kappaz = self.kappa(coords[2], 2)
        
        if cmplx:
            pw_obj = CpmlExCmplx(underneath.epsilon, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlExReal(underneath.epsilon, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        bz = self.b(coords[2], 2)
        bx = self.b(coords[0], 0)
        cz = self.c(coords[2], 2)
        cx = self.c(coords[0], 0)
        kappaz = self.kappa(coords[2], 2)
        kappax = self.kappa(coords[0], 0)
        
        if cmplx:
            pw_obj = CpmlEyCmplx(underneath.epsilon, bz, bx, cz, cx, kappaz, kappax)
        else:
            pw_obj = CpmlEyReal(underneath.epsilon, bz, bx, cz, cx, kappaz, kappax)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        bx = self.b(coords[0], 0)
        by = self.b(coords[1], 1)
        cx = self.c(coords[0], 0)
        cy = self.c(coords[1], 1)
        kappax = self.kappa(coords[0], 0)
        kappay = self.kappa(coords[1], 1)
        
        if cmplx:
            pw_obj = CpmlEzCmplx(underneath.epsilon, bx, by, cx, cy, kappax, kappay)
        else:
            pw_obj = CpmlEzReal(underneath.epsilon, bx, by, cx, cy, kappax, kappay)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        by = self.b(coords[1], 1)
        bz = self.b(coords[2], 2)
        cy = self.c(coords[1], 1)
        cz = self.c(coords[2], 2)
        kappay = self.kappa(coords[1], 1)
        kappaz = self.kappa(coords[2], 2)
        
        if cmplx:
            pw_obj = CpmlHxCmplx(underneath.mu, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlHxReal(underneath.mu, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        bz = self.b(coords[2], 2)
        bx = self.b(coords[0], 0)
        cz = self.c(coords[2], 2)
        cx = self.c(coords[0], 0)
        kappaz = self.kappa(coords[2], 2)
        kappax = self.kappa(coords[0], 0)
        
        if cmplx:
            pw_obj = CpmlHyCmplx(underneath.mu, bz, bx, cz, cx, kappaz, kappax)
        else:
            pw_obj = CpmlHyReal(underneath.mu, bz, bx, cz, cx, kappaz, kappax)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        bx = self.b(coords[0], 0)
        by = self.b(coords[1], 1)
        cx = self.c(coords[0], 0)
        cy = self.c(coords[1], 1)
        kappax = self.kappa(coords[0], 0)
        kappay = self.kappa(coords[1], 1)
        
        if cmplx:
            pw_obj = CpmlHzCmplx(underneath.mu, bx, by, cx, cy, kappax, kappay)
        else:
            pw_obj = CpmlHzReal(underneath.mu, bx, by, cx, cy, kappax, kappay)
            
        return pw_obj
        

class DrudePole(object):
    def __init__(self, omega, gamma):
        """
        DrudePole() -> a new Drude pole
        omega: a plasma frequency
        gamma: a relaxation frequency
        
        """
        self.omega = float(omega)
        self.gamma = float(gamma)
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude pole"
        print " " * indent,
        print "plasma frequency:", self.omega,
        print "relaxation frequency:", self.gamma
        

class LorentzPole(object):
    def __init__(self, amp, omega, gamma):
        """
        LorentzPole() -> a new Lorentz pole
        amp: amplitude
        omega: energy of the gap
        gamma: broadening
        
        """
        self.amp = float(amp)
        self.omega = float(omega)
        self.gamma = float(gamma)
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Lorentz pole"
        print " " * indent,
        print "amplitude:", self.amp,
        print "energy of the gap:", self.omega,
        print "broadening:", self.gamma
        
                
class CriticalPoint(object):
    def __init__(self, amp, phi, omega, gamma):
        """
        CriticalPoint() -> a new critical point
        amp: amplitude 
        phi: phase
        omega: energy of the gap
        gamma: broadening
         
        """
        self.amp = float(amp)
        self.phi = float(phi)
        self.omega = float(omega)
        self.gamma = float(gamma)
    
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "critical point"
        print " " * indent,
        print "amplitude:", self.amp,
        print "phase:", self.phi,
        print "energy of the gap:", self.omega,
        print "broadening:", self.gamma
        
        
class DCP(Dielectric):
    """
    The auxiliary differential equation implementation of Drude-critical points model 
    based on the following references.
    * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the optical 
      properties of gold," The Journal of Chemical Physics, vol. 125, no. 16, 
      pp. 164705-3, Oct. 2006.
    * A. Taflove and S. C. Hagness, Computational Electrodynamics: The Finite-
      Difference Time-Domain Method, Third Edition, 3rd ed. Artech House Publishers, 
      2005.
    * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the
      optical properties of gold," J. Chem. Phys. 125, 164705, 2001.
      
    """
    def __init__(self, epsilon=1, mu=1, sigma=0, dps=(), cps=()):
        """
        epsilon: The (frequency-independent) isotropic relative permittivity. Default is 1.
        mu: The (frequency-independent) isotropic relative permeability. Default is 1.
        sigma: The (frequency-independent) isotropic conductivity. Default is 0.
        dps: list of Drude poles. Default is ().
        cps: list of critical points. Default is ().
        
        """
        Dielectric.__init__(self, epsilon=epsilon, mu=mu)
        self.sigma = float(sigma) # instant conductivity
        self.dps = tuple(dps) # tuple of Drude poles
        self.cps = tuple(cps) # tuple of critical points
        
    def init(self, space, param=None):
        self.dt = space.dt
        
        # parameters for the ADE of the Drude model
        self.a = empty((len(self.dps),3), float)
        for i in xrange(len(self.dps)):
            pole = self.dps[i]
            denom = self.dt * pole.gamma + 2.
            self.a[i,0] = (self.dt * pole.gamma - 2) / denom
            self.a[i,1] = 4 / denom
            self.a[i,2] = 2 * (self.dt * pole.omega)**2 / denom
        
        # parameters for the ADE of critical points model
        self.b = empty((len(self.cps),4), float)
        for i in xrange(len(self.cps)):
            pnt = self.cps[i]
            denom = self.dt * pnt.gamma + 1.
            self.b[i,0] = (self.dt * pnt.gamma - 1) / denom
            self.b[i,1] = (2 - self.dt**2 * (pnt.gamma**2 + pnt.omega**2)) / denom
            self.b[i,2] = self.dt * sin(pnt.phi) * pnt.amp * pnt.omega / denom
            self.b[i,3] = 2 * self.dt**2 * pnt.amp * pnt.omega * \
            (cos(pnt.phi) * pnt.omega - sin(pnt.phi) * pnt.gamma) / denom
            
        # parameters for the electric field update equations.
        self.c = empty(4, float)
        denom = self.dt * self.sigma + 2. * (self.epsilon - sum(self.b[:,2]))
        self.c[0] = 2 * self.dt / denom
        self.c[1] = -2 / denom
        self.c[2] = -2 * sum(self.b[:,2]) / denom
        self.c[3] = (2 * (self.epsilon - sum(self.a[:,2]) - sum(self.b[:,3])) - self.dt * self.sigma) / denom
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude-critical points dispersive media"
        print " " * indent,
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu,
        print "conductivity:", self.sigma
        
        print " " * indent, "Drude pole(s):"
        for i in self.dps:
            i.display_info(indent+4)
        print " " * indent, "critical point(s):"
        for i in self.cps:
            i.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPExCmplx(self.epsilon, self.a, self.b, self.c)
        else:
            pw_obj = DCPExReal(self.epsilon, self.a, self.b, self.c)
                
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPEyCmplx(self.epsilon, self.a, self.b, self.c)
        else:
            pw_obj = DCPEyReal(self.epsilon, self.a, self.b, self.c)
                
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPEzCmplx(self.epsilon, self.a, self.b, self.c)
        else:
            pw_obj = DCPEzReal(self.epsilon, self.a, self.b, self.c)
        
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPHxCmplx(self.mu)
        else:
            pw_obj = DCPHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPHyCmplx(self.mu)
        else:
            pw_obj = DCPHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPHzCmplx(self.mu)
        else:
            pw_obj = DCPHzReal(self.mu)
            
        return pw_obj
    

class DCPPLRC(Dielectric):
    """
    The piecewise-linear recursive-convolution implementation of 
    Drude-critical points model based on the following references.
    * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the optical 
      properties of gold," The Journal of Chemical Physics, vol. 125, no. 16, 
      pp. 164705-3, Oct. 2006.
    * A. Taflove and S. C. Hagness, Computational Electrodynamics: The Finite-
      Difference Time-Domain Method, Third Edition, 3rd ed. Artech House Publishers, 
      2005.

    """
    def __init__(self, epsilon=1, mu=1, sigma=0, dps=(), cps=()):
        """
        epsilon: The (frequency-independent) isotropic relative permittivity. Default is 1.
        mu: The (frequency-independent) isotropic relative permeability. Default is 1.
        sigma: The (frequency-independent) isotropic conductivity. Default is 0.
        dps: list of Drude poles. Default is ().
        cps: list of critical points. Default is ().
        
        """
        Dielectric.__init__(self, epsilon=epsilon, mu=mu)
        self.sigma = float(sigma) # instant conductivity
        self.dps = tuple(dps) # tuple of Drude poles
        self.cps = tuple(cps) # tuple of critical points
        
    def init(self, space, param=None):
        self.dt = space.dt

        # parameters of the recursion relation for the Drude pole recursive accumulator.
        self.a = empty((len(self.dps),3), float)
        for i in xrange(len(self.dps)):
            pole = self.dps[i]
            self.a[i,0] = self.delta_chi_dp_0(pole) - self.delta_xi_dp_0(pole)
            self.a[i,1] = self.delta_xi_dp_0(pole)
            self.a[i,2] = exp(-pole.gamma * self.dt)
        
        # parameters of the recursion relation for the critical point recursive accumulator.
        self.b = empty((len(self.cps),3), complex)
        for i in xrange(len(self.cps)):
            pnt = self.cps[i]
            self.b[i,0] = self.delta_chi_cp_0(pnt) - self.delta_xi_cp_0(pnt)
            self.b[i,1] = self.delta_xi_cp_0(pnt)
            self.b[i,2] = exp(-self.dt * (pnt.gamma + 1j * pnt.omega))
            
        # parameters for the electric field update equations.
        chi_0 = sum(map(self.chi_dp_0, self.dps) + map(self.chi_cp_0, self.cps)).real
        xi_0 = sum(map(self.xi_dp_0, self.dps) + map(self.xi_cp_0, self.cps)).real
        
        self.c = empty(3, float)
        denom = self.epsilon - xi_0 + chi_0
        self.c[0] = (self.epsilon - xi_0) / denom
        self.c[1] = self.dt / denom
        self.c[2] = 1 / denom
        
    def chi_dp_0(self, dp):
        omega = dp.omega
        gamma = dp.gamma
        gdt = dp.gamma * self.dt
        
        return (omega / gamma)**2 * (gdt + exp(-gdt) - 1)

    def xi_dp_0(self, dp):
        omega = dp.omega
        gamma = dp.gamma
        gdt = dp.gamma * self.dt
        
        return (omega / gamma)**2 * (gdt / 2 - (1 - (1 + gdt) * exp(-gdt)) / gdt)
    
    def delta_chi_dp_0(self, dp):
        omega = dp.omega
        gamma = dp.gamma
        gdt = dp.gamma * self.dt
        
        return -(omega / gamma * (1 - exp(-gdt)))**2
    
    def delta_xi_dp_0(self, dp):
        gdt = dp.gamma * self.dt
        delta_chi_dp_0 = self.delta_chi_dp_0(dp)
        
        return delta_chi_dp_0 * (1 / (1 - exp(gdt)) + 1 / gdt) 
    
    def chi_cp_0(self, cp):
        go = cp.gamma + 1j * cp.omega
        return 2j * cp.amp * cp.omega * exp(1j * cp.phi) * (1 - exp(-self.dt * go)) / go
    
    def xi_cp_0(self, cp):
        dtgo = self.dt * (cp.gamma + 1j * cp.omega)
        chi_cp_0 = self.chi_cp_0(cp)
        
        return chi_cp_0 * (1 / (1 - exp(dtgo)) + 1 / dtgo)
    
    def delta_chi_cp_0(self, cp):
        go = cp.gamma + 1j * cp.omega
        return 2j * cp.amp * cp.omega * exp(1j * cp.phi) * (1 - exp(-self.dt * go))**2 / go
    
    def delta_xi_cp_0(self, cp):
        dtgo = self.dt * (cp.gamma + 1j * cp.omega)
        delta_chi_cp_0 = self.delta_chi_cp_0(cp)
        
        return delta_chi_cp_0 * (1 / (1 - exp(dtgo)) + 1 / dtgo)
    
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude-critical points dispersive media"
        print " " * indent,
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu,
        print "conductivity:", self.sigma
        
        print " " * indent, "Drude pole(s):"
        for i in self.dps:
            i.display_info(indent+4)
        print " " * indent, "critical point(s):"
        for i in self.cps:
            i.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPPLRCExCmplx(self.epsilon, self.a, self.b, self.c)
        else:
            pw_obj = DCPPLRCExReal(self.epsilon, self.a, self.b, self.c)
                
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPPLRCEyCmplx(self.epsilon, self.a, self.b, self.c)
        else:
            pw_obj = DCPPLRCEyReal(self.epsilon, self.a, self.b, self.c)
                
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPPLRCEzCmplx(self.epsilon, self.a, self.b, self.c)
        else:
            pw_obj = DCPPLRCEzReal(self.epsilon, self.a, self.b, self.c)
        
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPPLRCHxCmplx(self.mu)
        else:
            pw_obj = DCPPLRCHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPPLRCHyCmplx(self.mu)
        else:
            pw_obj = DCPPLRCHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DCPPLRCHzCmplx(self.mu)
        else:
            pw_obj = DCPPLRCHzReal(self.mu)
            
        return pw_obj
    
    
class DCPRC(DCPPLRC):
    """
    The recursive convolution implementation of Drude-critical points model
    based on the following articles.
    * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the optical 
      properties of gold," The Journal of Chemical Physics, vol. 125, no. 16, 
      pp. 164705-3, Oct. 2006.
    * A. Vial, "Implementation of the critical points model in the recursive 
      convolution method for modelling dispersive media with the finite-difference 
      time domain method," Journal of Optics A: Pure and Applied Optics, vol. 9, 
      Jul. 2007, pp. 745-748.
    
    """
    def xi_dp_0(self, dp):
        return .0
    
    def xi_cp_0(self, cp):
        return .0

    def delta_xi_dp_0(self, dp):
        return .0
    
    def delta_xi_cp_0(self, cp):
        return .0
        
        
class Drude(Dielectric):
    """
    The auxiliary differential equation implementation of the Drude model based on the 
    following article.
    * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
    Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006.
    
    """
    def __init__(self, epsilon=1, mu=1, sigma=0, dps=()):
        """
        Arguments:
            epsilon: The (frequency-independent) isotropic relative permittivity. Default is 1.
            mu: The (frequency-independent) isotropic relative permeability. Default is 1.
            sigma: The (frequency-independent) isotropic conductivity. Default is 0.
            dps: list of Drude poles. Default is().
            
        """
        Dielectric.__init__(self, epsilon=epsilon, mu=mu)
        self.sigma = float(sigma)
        self.dps = tuple(dps)
        
    def init(self, space, param=None):
        self.dt = space.dt
        
        # parameters for the ADE of the Drude model.
        self.a = empty((len(self.dps),3), float)
        for i in xrange(len(self.dps)):
            pole = self.dps[i]
            denom = self.dt * pole.gamma + 2.
            self.a[i,0] = (self.dt * pole.gamma - 2) / denom
            self.a[i,1] = 4 / denom
            self.a[i,2] = 2 * (self.dt * pole.omega)**2 / denom
        
        # parameters for the electric field update equations.
        self.c = empty(3, float)
        denom = 2. * self.epsilon + self.dt * self.sigma 
        self.c[0] = 2 * self.dt / denom
        self.c[1] = -2 / denom
        self.c[2] = (2 * self.epsilon - self.dt * self.sigma) / denom
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude dispersion media"
        print " " * indent, 
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu,
        print "conductivity:", self.sigma
        
        print " "* indent, "Drude pole(s):"
        for p in self.dps:
            p.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeExCmplx(self.epsilon, self.a, self.c)
        else:
            pw_obj = DrudeExReal(self.epsilon, self.a, self.c)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeEyCmplx(self.epsilon, self.a, self.c)
        else:
            pw_obj = DrudeEyReal(self.epsilon, self.a, self.c)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeEzCmplx(self.epsilon, self.a, self.c)
        else:
            pw_obj = DrudeEzReal(self.epsilon, self.a, self.c)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeHxCmplx(self.mu)
        else:
            pw_obj = DrudeHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeHyCmplx(self.mu)
        else:
            pw_obj = DrudeHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeHzCmplx(self.mu)
        else:
            pw_obj = DrudeHzReal(self.mu)
            
        return pw_obj
        

class Lorentz(Dielectric):
    """
    The auxiliary differential equation implementation of the Lorentz model.
    
    """
    def __init__(self, epsilon=1, mu=1, sigma=0, lps=()):
        """
        Arguments:
            epsilon: The (frequency-independent) isotropic relative permittivity. Default is 1.
            mu: The (frequency-independent) isotropic relative permeability. Default is 1.
            sigma: The (frequency-independent) isotropic conductivity. Default is 0.
            lps: list of Lorentz poles. Default is().
            
        """
        Dielectric.__init__(self, epsilon=epsilon, mu=mu)
        self.sigma = float(sigma)
        self.lps = tuple(lps)
        
    def init(self, space, param=None):
        self.dt = space.dt
        
        # parameters for the ADE of the Drude model.
        self.a = empty((len(self.lps),3), float)
        for i in xrange(len(self.lps)):
            pole = self.lps[i]
            denom = self.dt * pole.gamma + 2.
            self.a[i,0] = (self.dt * pole.gamma - 2) / denom
            self.a[i,1] = (4 - 2 * (self.dt * pole.omega)**2) / denom
            self.a[i,2] = 2 * pole.amp * (self.dt * pole.omega)**2 / denom
        
        # parameters for the electric field update equations.
        self.c = empty(3, float)
        denom = 2. * self.epsilon + self.dt * self.sigma 
        self.c[0] = 2 * self.dt / denom
        self.c[1] = -2 / denom
        self.c[2] = (2 * self.epsilon - self.dt * self.sigma) / denom
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Lorentz dispersion media"
        print " " * indent, 
        print "permittivity:", self.epsilon,
        print "permeability:", self.mu,
        print "conductivity:", self.sigma
        
        print " "* indent, "Lorentz pole(s):"
        for p in self.lps:
            p.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzExCmplx(self.epsilon, self.a, self.c)
        else:
            pw_obj = LorentzExReal(self.epsilon, self.a, self.c)
            
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzEyCmplx(self.epsilon, self.a, self.c)
        else:
            pw_obj = LorentzEyReal(self.epsilon, self.a, self.c)
            
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzEzCmplx(self.epsilon, self.a, self.c)
        else:
            pw_obj = LorentzEzReal(self.epsilon, self.a, self.c)
            
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzHxCmplx(self.mu)
        else:
            pw_obj = LorentzHxReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzHyCmplx(self.mu)
        else:
            pw_obj = LorentzHyReal(self.mu)
            
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzHzCmplx(self.mu)
        else:
            pw_obj = LorentzHzReal(self.mu)
            
        return pw_obj
    
            
class Gold(DCP):
    """
    The parameters are from the following article.
    * A. Vial and T. Laroche, "Comparison of gold and silver dispersion laws
      suitable for FDTD simulations," Appl. Phys. B, 93, 139-143, 2008.
    
    These parameters represents the permittivity of gold in 200-1,000 nm range.
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.3202e16 * a / const.c0, 
                        gamma=1.0805e14 * a / const.c0)
        cp1 = CriticalPoint(amp=0.26698, 
                            phi=-1.2371, 
                            omega=3.8711e15 * a / const.c0, 
                            gamma=4.4642e14 * a / const.c0)
        cp2 = CriticalPoint(amp=3.0834, 
                            phi=-1.0968, 
                            omega=4.1684e15 * a / const.c0, 
                            gamma=2.3555e15 * a / const.c0)
        DCP.__init__(self, epsilon=1.1431, mu=1, sigma=0, dps=(dp1,), cps=(cp1,cp2))
        

class GoldPLRC(DCPPLRC):
    """
    The parameters are from the following article.
    * A. Vial and T. Laroche, "Comparison of gold and silver dispersion laws
      suitable for FDTD simulations," Appl. Phys. B, 93, 139-143, 2008.
    
    These parameters represents the permittivity of gold in 200-1,000 nm range.
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.3202e16 * a / const.c0, 
                        gamma=1.0805e14 * a / const.c0)
        cp1 = CriticalPoint(amp=0.26698, 
                            phi=-1.2371, 
                            omega=3.8711e15 * a / const.c0, 
                            gamma=4.4642e14 * a / const.c0)
        cp2 = CriticalPoint(amp=3.0834, 
                            phi=-1.0968, 
                            omega=4.1684e15 * a / const.c0, 
                            gamma=2.3555e15 * a / const.c0)
        DCPPLRC.__init__(self, epsilon=1.1431, mu=1, sigma=0, dps=(dp1,), cps=(cp1,cp2))
        

class GoldRC(DCPRC):
    """
    The parameters are from the following article.
    * A. Vial and T. Laroche, "Comparison of gold and silver dispersion laws
      suitable for FDTD simulations," Appl. Phys. B, 93, 139-143, 2008.
    
    These parameters represents the permittivity of gold in 200-1,000 nm range.
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.3202e16 * a / const.c0, 
                        gamma=1.0805e14 * a / const.c0)
        cp1 = CriticalPoint(amp=0.26698, 
                            phi=-1.2371, 
                            omega=3.8711e15 * a / const.c0, 
                            gamma=4.4642e14 * a / const.c0)
        cp2 = CriticalPoint(amp=3.0834, 
                            phi=-1.0968, 
                            omega=4.1684e15 * a / const.c0, 
                            gamma=2.3555e15 * a / const.c0)
        DCPRC.__init__(self, epsilon=1.1431, mu=1, sigma=0, dps=(dp1,), cps=(cp1,cp2))
                        
class Gold2(Drude):
    """
    The parameters are from the following article.
    * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
    Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006.

    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.196e16 * a / const.c0, 
                        gamma=8.052e13 * a / const.c0)
        Drude.__init__(self, epsilon=1, mu=1, sigma=0, dps=(dp1,))
        

class Meep(Lorentz):
    """
    The parameters are from the 'Meep Tutorial/Material dispersion'.

    """
    def __init__(self):
        """
        
        """
        lp1 = LorentzPole(omega=1.1,
                          gamma=1e-5,
                          amp=0.5)
        lp2 = LorentzPole(omega=0.5,
                          gamma=0.1,
                          amp=2e-5)
        Lorentz.__init__(self, epsilon=2.25, mu=1, sigma=0, lps=(lp1,lp2))
        
        
class Silver(DCP):
    """
    The parameters are from the following article.
    * A. Vial and T. Laroche, "Comparison of gold and silver dispersion laws
      suitable for FDTD simulations," Appl. Phys. B, 93, 139-143, 2008.
    
    These parameters represents the permittivity of silver in 200-1,000 nm range.
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=1.3861e16 * a / const.c0, 
                        gamma=4.5841e13 * a / const.c0)
        cp1 = CriticalPoint(amp=1.0171, 
                            phi=-0.93935, 
                            omega=6.6327e15 * a / const.c0, 
                            gamma=1.6666e15 * a / const.c0)
        cp2 = CriticalPoint(amp=15.797, 
                            phi=1.8087, 
                            omega=9.2726e17 * a / const.c0, 
                            gamma=2.3716e17 * a / const.c0)
        DCP.__init__(self, epsilon=15.833, mu=1, sigma=0, dps=(dp1,), cps=(cp1,cp2))
        

class Aluminum(DCP):
    """
    The parameters are from the following article.
    * A. Vial and T. Laroche, "Comparison of gold and silver dispersion laws
      suitable for FDTD simulations," Appl. Phys. B, 93, 139-143, 2008.
    
    These parameters represents the permittivity of silver in 200-1,000 nm range.
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=2.0598e16 * a / const.c0, 
                        gamma=2.2876e14 * a / const.c0)
        cp1 = CriticalPoint(amp=5.2306, 
                            phi=-0.51202, 
                            omega=2.2694e15 * a / const.c0, 
                            gamma=3.2867e14 * a / const.c0)
        cp2 = CriticalPoint(amp=5.2704, 
                            phi=0.42503, 
                            omega=2.4668e15 * a / const.c0, 
                            gamma=1.7731e15 * a / const.c0)
        DCP.__init__(self, epsilon=1.0000, mu=1, sigma=0, dps=(dp1,), cps=(cp1,cp2))
        

class Chromium(DCP):
    """
    The parameters are from the following article.
    * A. Vial and T. Laroche, "Comparison of gold and silver dispersion laws
      suitable for FDTD simulations," Appl. Phys. B, 93, 139-143, 2008.
    
    These parameters represents the permittivity of silver in 200-1,000 nm range.
    
    """
    def __init__(self, a):
        """
        a: lattice constant in meters.
        
        """
        dp1 = DrudePole(omega=8.8128e15 * a / const.c0, 
                        gamma=3.8828e14 * a / const.c0)
        cp1 = CriticalPoint(amp=33.086, 
                            phi=-0.25722, 
                            omega=1.7398e15 * a / const.c0, 
                            gamma=1.6329e15 * a / const.c0)
        cp2 = CriticalPoint(amp=1.6592, 
                            phi=0.83533, 
                            omega=3.7925e15 * a / const.c0, 
                            gamma=7.3567e14 * a / const.c0)
        DCP.__init__(self, epsilon=1.1297, mu=1, sigma=0, dps=(dp1,), cps=(cp1,cp2))
        
