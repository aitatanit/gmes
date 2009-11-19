#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

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

    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pointwise_mateiral_ey(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        """Return PointwiseMaterial object of the given point.
        
        Arguments:
            idx -- (local) array index of the target point
            coords -- (global) space coordinate of the target point
            complex -- whether the EM field has complex value. Default is False.
            underneath -- underneath material object of the target point.
            
        """
        raise NotImplementedError
    
    def init(self, space, param=None):
        self.cmplx = space.cmplx
        
        
class Dummy(Material):
    """A dummy material type which dosen't update the field component.
    
    """
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "dummy object"
        
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DummyExCmplx(idx, self.epsilon_r)
        else:
            pw_obj = DummyExReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_mateiral_ey(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DummyEyCmplx(idx, self.epsilon_r)
        else:
            pw_obj = DummyEyReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DummyEzCmplx(idx, self.epsilon_r)
        else:
            pw_obj = DummyEzReal(idx, self.epsilon_r)
        
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DummyHxCmplx(idx, self.mu_r)
        else:
            pw_obj = DummyHxReal(idx, self.mu_r)
        
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DummyHyCmplx(idx, self.mu_r)
        else:
            pw_obj = DummyHyReal(idx, self.mu_r)
        
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DummyHzCmplx(idx, self.mu_r)
        else:
            pw_obj = DummyHzReal(idx, self.mu_r)
        
        return pw_obj

    
class Zero(Material):
    """A material type which sets the field value zero.
    
    """
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "zero object"
        
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = ZeroExCmplx(idx, self.epsilon_r)
        else:
            pw_obj = ZeroExReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_mateiral_ey(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = ZeroEyCmplx(idx, self.epsilon_r)
        else:
            pw_obj = ZeroEyReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = ZeroEzCmplx(idx, self.epsilon_r)
        else:
            pw_obj = ZeroEzReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = ZeroHxCmplx(idx, self.mu_r)
        else:
            pw_obj = ZeroHxReal(idx, self.mu_r)
            
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = ZeroHyCmplx(idx, self.mu_r)
        else:
            pw_obj = ZeroHyReal(idx, self.mu_r)
            
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = ZeroHzCmplx(idx, self.mu_r)
        else:
            pw_obj = ZeroHzReal(idx, self.mu_r)
            
        return pw_obj
    

class One(Material):
    """A material type which sets the field value one.
    
    """
    def display_info(self, indentby=0):
        """Display the parameter values.
        
        """
        print " " * indentby, "one object"
        
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = OneExCmplx(idx, self.epsilon_r)
        else:
            pw_obj = OneExReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_mateiral_ey(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = OneEyCmplx(idx, self.epsilon_r)
        else:
            pw_obj = OneEyReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = OneEzCmplx(idx, self.epsilon_r)
        else:
            pw_obj = OneEzReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = OneHxCmplx(idx, self.mu_r)
        else:
            pw_obj = OneHxReal(idx, self.mu_r)
            
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = OneHyCmplx(idx, self.mu_r)
        else:
            pw_obj = OneHyReal(idx, self.mu_r)
            
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = OneHzCmplx(idx, self.mu_r)
        else:
            pw_obj = OneHzReal(idx, self.mu_r)
            
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
        if self.cmplx:
            pw_obj = DielectricExCmplx(idx, self.epsilon_r)
        else:
            pw_obj = DielectricExReal(idx, self.epsilon_r)
            
        return pw_obj

    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DielectricEyCmplx(idx, self.epsilon_r)
        else:
            pw_obj = DielectricEyReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DielectricEzCmplx(idx, self.epsilon_r)
        else:
            pw_obj = DielectricEzReal(idx, self.epsilon_r)
            
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DielectricHxCmplx(idx, self.mu_r)
        else:
            pw_obj = DielectricHxReal(idx, self.mu_r)
            
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DielectricHyCmplx(idx, self.mu_r)
        else:
            pw_obj = DielectricHyReal(idx, self.mu_r)
            
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DielectricHzCmplx(idx, self.mu_r)
        else:
            pw_obj = DielectricHzReal(idx, self.mu_r)
            
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
        self.cmplx = space.cmplx
        
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
        eta = sqrt(self.effective_mu_r / self.effective_epsilon_r) * const.Z0
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
        effective_epsilon_r -- the effective relative permittivity of incident mode impinging on the PML boundary
        effective_mu_r -- the effective relative permeability of incident mode impinging on the PML boundary
        sigma_max_ratio -- the ratio between sigma_max and sigma_opt

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
        
        if self.cmplx:
            pw_obj = UpmlExCmplx(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlExReal(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        c1 = self.c1(coords[2], 2)
        c2 = self.c2(coords[2], 2)
        c3 = self.c3(coords[0], 0)
        c4 = self.c4(coords[0], 0)
        c5 = self.c5(coords[1], 1)
        c6 = self.c6(coords[1], 1)
        
        if self.cmplx:
            pw_obj = UpmlEyCmplx(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlEyReal(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        c1 = self.c1(coords[0], 0)
        c2 = self.c2(coords[0], 0)
        c3 = self.c3(coords[1], 1)
        c4 = self.c4(coords[1], 1)
        c5 = self.c5(coords[2], 2)
        c6 = self.c6(coords[2], 2)
        
        if self.cmplx:
            pw_obj = UpmlEzCmplx(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlEzReal(idx, underneath.epsilon_r, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        c1 = self.c1(coords[1], 1)
        c2 = self.c2(coords[1], 1)
        c3 = self.c3(coords[2], 2)
        c4 = self.c4(coords[2], 2)
        c5 = self.c5(coords[0], 0)
        c6 = self.c6(coords[0], 0)
        
        if self.cmplx:
            pw_obj = UpmlHxCmplx(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlHxReal(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        c1 = self.c1(coords[2], 2)
        c2 = self.c2(coords[2], 2)
        c3 = self.c3(coords[0], 0)
        c4 = self.c4(coords[0], 0)
        c5 = self.c5(coords[1], 1)
        c6 = self.c6(coords[1], 1)
        
        if self.cmplx:
            pw_obj = UpmlHyCmplx(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlHyReal(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        c1 = self.c1(coords[0], 0)
        c2 = self.c2(coords[0], 0)
        c3 = self.c3(coords[1], 1)
        c4 = self.c4(coords[1], 1)
        c5 = self.c5(coords[2], 2)
        c6 = self.c6(coords[2], 2)
        
        if self.cmplx:
            pw_obj = UpmlHzCmplx(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
        else:
            pw_obj = UpmlHzReal(idx, underneath.mu_r, c1, c2, c3, c4, c5, c6)
            
        return pw_obj
    
    
class CPML(PML):
    """Form Complex Frequency Shifted (CFS) Perfectly Matched Layer (PML).
    
    This class implements CFS PML represented in
    S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,"
    Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
    Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
    2005, pp. 273-328.
    
    Attributes:
        m -- default 3
        kappa_max -- default 15
        a_max -- default 0. CPML works like UPML when a_max = 0.
        effective_epsilon_r -- the effective relative permittivity of incident mode impinging on the PML boundary
        effective_mu_r -- the effective relative permeability of incident mode impinging on the PML boundary
        sigma_max_ratio -- default 1
    
    """
    def __init__(self, effective_epsilon_r=1, effective_mu_r=1, m=3, kappa_max=2, m_a=1, a_max=0, sigma_max_ratio=2):
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
        
        if self.cmplx:
            pw_obj = CpmlExCmplx(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlExReal(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        bz = self.b(coords[2], 2)
        bx = self.b(coords[0], 0)
        cz = self.c(coords[2], 2)
        cx = self.c(coords[0], 0)
        kappaz = self.kappa(coords[2], 2)
        kappax = self.kappa(coords[0], 0)
        
        if self.cmplx:
            pw_obj = CpmlEyCmplx(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlEyReal(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        bx = self.b(coords[0], 0)
        by = self.b(coords[1], 1)
        cx = self.c(coords[0], 0)
        cy = self.c(coords[1], 1)
        kappax = self.kappa(coords[0], 0)
        kappay = self.kappa(coords[1], 1)
        
        if self.cmplx:
            pw_obj = CpmlEzCmplx(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlEzReal(idx, underneath.epsilon_r, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        by = self.b(coords[1], 1)
        bz = self.b(coords[2], 2)
        cy = self.c(coords[1], 1)
        cz = self.c(coords[2], 2)
        kappay = self.kappa(coords[1], 1)
        kappaz = self.kappa(coords[2], 2)
        
        if self.cmplx:
            pw_obj = CpmlHxCmplx(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlHxReal(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        bz = self.b(coords[2], 2)
        bx = self.b(coords[0], 0)
        cz = self.c(coords[2], 2)
        cx = self.c(coords[0], 0)
        kappaz = self.kappa(coords[2], 2)
        kappax = self.kappa(coords[0], 0)
        
        if self.cmplx:
            pw_obj = CpmlHyCmplx(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlHyReal(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        bx = self.b(coords[0], 0)
        by = self.b(coords[1], 1)
        cx = self.c(coords[0], 0)
        cy = self.c(coords[1], 1)
        kappax = self.kappa(coords[0], 0)
        kappay = self.kappa(coords[1], 1)
        
        if self.cmplx:
            pw_obj = CpmlHzCmplx(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
        else:
            pw_obj = CpmlHzReal(idx, underneath.mu_r, by, bz, cy, cz, kappay, kappaz)
            
        return pw_obj
        
        
class Drude(Dielectric):
    """
    
    """
    def __init__(self, epsilon_inf, omega_p, gamma_p, mu_r=1):
        """
        Arguments:
            epsilon_inf: relative permittivity at infinite frequency
            omega_p: the pole resonant frequency 
            gamma_p: the inverse of the pole relaxation time
            mu_r: relative magnetic permeability
        """
        Dielectric.__init__(self, epsilon_r=epsilon_inf, mu_r=mu_r)
        self.omega_p = array(omega_p, float)
        self.gamma_p = array(gamma_p, float)
    
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude dispersion media"
        print " " * indent, 
        print "infinite permittivity:", self.epsilon_r,
        print "plasma frequency:", self.omega_p
        print "relaxation frequency:", self.gamma_p
        print "relative permeability:", self.mu_r
        
    def get_pointwise_material_ex(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DrudeExCmplx(idx, self.epsilon_r, self.omega_p, self.gamma_p)
        else:
            pw_obj = DrudeExReal(idx, self.epsilon_r, self.omega_p, self.gamma_p)
            
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DrudeEyCmplx(idx, self.epsilon_r, self.omega_p, self.gamma_p)
        else:
            pw_obj = DrudeEyReal(idx, self.epsilon_r, self.omega_p, self.gamma_p)
            
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DrudeEzCmplx(idx, self.epsilon_r, self.omega_p, self.gamma_p)
        else:
            pw_obj = DrudeEzReal(idx, self.epsilon_r, self.omega_p, self.gamma_p)
            
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DrudeHxCmplx(idx, self.mu_r, self.omega_p, self.gamma_p)
        else:
            pw_obj = DrudeHxReal(idx, self.mu_r, self.omega_p, self.gamma_p)
            
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DrudeHyCmplx(idx, self.mu_r, self.omega_p, self.gamma_p)
        else:
            pw_obj = DrudeHyReal(idx, self.mu_r, self.omega_p, self.gamma_p)
            
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, coords, underneath=None):
        if self.cmplx:
            pw_obj = DrudeHzCmplx(idx, self.mu_r, self.omega_p, self.gamma_p)
        else:
            pw_obj = DrudeHzReal(idx, self.mu_r, self.omega_p, self.gamma_p)
            
        return pw_obj
        