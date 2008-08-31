#!/usr/bin/env python

from math import sqrt, exp
from numpy import array, inf

from pointwise_material import *

from constants import Z0, epsilon0


class Material:
    def display_info(self, indentby=0):
        """Display the parameter values
        """
        print " " * indentby, "material type object"

    def get_pointwise_material_ex(self, idx, co):
        raise NotImplementedError
    
    def get_pointwise_mateiral_ey(self, idx, co):
        raise NotImplementedError
    
    def get_pointwise_material_ez(self, idx, co):
        raise NotImplementedError
    
    def get_pointwise_material_hx(self, idx, co):
        raise NotImplementedError
    
    def get_pointwise_material_hy(self, idx, co):
        raise NotImplementedError
    
    def get_pointwise_material_hz(self, idx, co):
        raise NotImplementedError
    
    
class Dummy(Material):
    def display_info(self, indentby=0):
        """
        """
        print " " * indentby, "dummy object"
        
    def get_pointwise_material_ex(self, idx, co):
        pw_obj = DummyEx(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_mateiral_ey(self, idx, co):
        pw_obj = DummyEy(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, co):
        pw_obj = DummyEz(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, co):
        pw_obj = DummyHx(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, co):
        pw_obj = DummyHy(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, co):
        pw_obj = DummyHz(idx, self.mu_r)
        return pw_obj
    
    
class Dielectric(Material):
    def __init__(self, epsilon_r=1, mu_r=1):
        """Representation of dielectric medium.
        
        epsilon_r: relative permittivity
        mu_r: relative permeability
        """
        self.epsilon_r = epsilon_r
        self.mu_r = mu_r

    def display_info(self, indent=0):
        """Display the parameter values
        """
        print " " * indent, "dielectric"
        print " " * indent, 
        print "relative permittivity:", self.epsilon_r,
        print "relative permeability:", self.mu_r

    def get_pointwise_material_ex(self, idx, co):
        pw_obj = DielectricEx(idx, self.epsilon_r)
        return pw_obj

    def get_pointwise_material_ey(self, idx, co):
        pw_obj = DielectricEy(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, co):
        pw_obj = DielectricEz(idx, self.epsilon_r)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, co):
        pw_obj = DielectricHx(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, co):
        pw_obj = DielectricHy(idx, self.mu_r)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, co):
        pw_obj = DielectricHz(idx, self.mu_r)
        return pw_obj


class PML(Material):
    def init(self, thickness=None, space=None):
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
        """ Calculate the optimal value of sigma.
        """
        eta = sqrt(self.mu_r / self.epsilon_r) * Z0
        return .8 * (self.m + 1) / eta / self.dw
    

class UPML(PML):
    """The implementation class of Uniaxial PML
    
    This class implements Uniaxial PML represented in 'A. Taflove and S. C. Hagness, 
    Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
    3rd ed., Artech House, Inc., 2005.
    """
    def __init__(self, epsilon_r=1, mu_r=1, m=3, kappa_max=15, sigma_max_ratio=1):
        self.initialized = False
        
        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.epsilon_r = float(epsilon_r)
        self.mu_r = float(mu_r)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
    def display_info(self, indent=0):
        print " " * indent, "UPML"
        print " " * indent, 
        print "relative permittivity:", self.epsilon_r,
        print "relative permeability:", self.mu_r,
        print "sigma_max:", self.sigma_max,
        print "m:", self.m,
        print "kappa_max:", self.kappa_max,
        
    def sigma(self, w, component):
        if w <= self.d - self.half_size[component]:
            return self.sigma_max[component] * (1 - (w + self.half_size[component]) / self.d)**self.m
        elif self.half_size[component] - self.d <= w:
            return self.sigma_max[component] * (1 + (w - self.half_size[component]) / self.d)**self.m
        else:
            return 0.0
        
    def kappa(self, w, component):
        if w <= self.d - self.half_size[component]:
            return 1 + (self.kappa_max - 1) * (1 - (w + self.half_size[component]) / self.d)**self.m
        elif self.half_size[component] - self.d <= w:
            return 1 + (self.kappa_max - 1) * (1 + (w - self.half_size[component]) / self.d)**self.m
        else:
            return 1.0
       
    def c1(self, w, component):
        numerator = 2 * epsilon0 * self.kappa(w, component) - self.sigma(w, component) * self.dt
        denominator = 2 * epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c2(self, w, component):
        numerator = 2 * epsilon0 * self.dt
        denominator = 2 * epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c3(self, w, component):
        numerator = 2 * epsilon0 * self.kappa(w, component) - self.sigma(w, component) * self.dt
        denominator = 2 * epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator / denominator
    
    def c4(self, w, component):
        denominator = 2 * epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return 1 / denominator
    
    def c5(self, w, component):
        numerator = 2 * epsilon0 * self.kappa(w, component) + self.sigma(w, component) * self.dt
        return numerator
    
    def c6(self, w, component):
        numerator = 2 * epsilon0 * self.kappa(w, component) - self.sigma(w, component) * self.dt
        return numerator
    
    def get_pointwise_material_ex(self, idx, co):
        c1 = self.c1(co[1], 1)
        c2 = self.c2(co[1], 1)
        c3 = self.c3(co[2], 2)
        c4 = self.c4(co[2], 2)
        c5 = self.c5(co[0], 0)
        c6 = self.c6(co[0], 0)
        pw_obj = UPMLEx(idx, self.epsilon_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, co):
        c1 = self.c1(co[2], 2)
        c2 = self.c2(co[2], 2)
        c3 = self.c3(co[0], 0)
        c4 = self.c4(co[0], 0)
        c5 = self.c5(co[1], 1)
        c6 = self.c6(co[1], 1)
        pw_obj = UPMLEy(idx, self.epsilon_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, co):
        c1 = self.c1(co[0], 0)
        c2 = self.c2(co[0], 0)
        c3 = self.c3(co[1], 1)
        c4 = self.c4(co[1], 1)
        c5 = self.c5(co[2], 2)
        c6 = self.c6(co[2], 2)
        pw_obj = UPMLEz(idx, self.epsilon_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, co):
        c1 = self.c1(co[1], 1)
        c2 = self.c2(co[1], 1)
        c3 = self.c3(co[2], 2)
        c4 = self.c4(co[2], 2)
        c5 = self.c5(co[0], 0)
        c6 = self.c6(co[0], 0)
        pw_obj = UPMLHx(idx, self.mu_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, co):
        c1 = self.c1(co[2], 2)
        c2 = self.c2(co[2], 2)
        c3 = self.c3(co[0], 0)
        c4 = self.c4(co[0], 0)
        c5 = self.c5(co[1], 1)
        c6 = self.c6(co[1], 1)
        pw_obj = UPMLHy(idx, self.mu_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, co):
        c1 = self.c1(co[0], 0)
        c2 = self.c2(co[0], 0)
        c3 = self.c3(co[1], 1)
        c4 = self.c4(co[1], 1)
        c5 = self.c5(co[2], 2)
        c6 = self.c6(co[2], 2)
        pw_obj = UPMLHz(idx, self.mu_r, c1, c2, c3, c4, c5, c6)
        return pw_obj
    
    
class CPML(PML):
    """The implementation class of CFS PML
    
    This class implements CFS PML represented in 'A. Taflove and S. C. Hagness, 
    Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
    3rd ed., Artech House, Inc., 2005.
    """
    def __init__(self, epsilon_r=1, mu_r=1, m=3, kappa_max=15, m_alpha=1, alpha_max=0.2, sigma_max_ratio=1):
        """
        epsilon_r: relative permittivity
        mu_r: relative permeability
        m:
        kappa_max:
        m_alpha:
        alpha_max:
        grid:
        """
        
        self.initialized = False
        
        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.m_alpha = float(m_alpha)
        self.alpha_max = float(alpha_max)
        self.epsilon_r = float(epsilon_r)
        self.mu_r = float(mu_r)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
    def display_info(self, indent=0):
        """Display the parameter values
        """
        print " " * indent, "CPML"
        print " " * indent, 
        print "relative permittivity:", self.epsilon_r,
        print "relative permeability:", self.mu_r,
        print "sigma_max:", self.sigma_max,
        print "m:", self.m,
        print "kappa_max:", self.kappa_max,
        print "m_alpha:", self.m_alpha,
        print "alpha_max:", self.alpha_max

    def sigma(self, w, component):
        if w <= self.d - self.half_size[component]:
            return self.sigma_max[component] * (1 - (w + self.half_size[component]) / self.d)**self.m
        elif self.half_size[component] - self.d <= w:
            return self.sigma_max[component] * (1 + (w - self.half_size[component]) / self.d)**self.m
        else:
            return 0.0
        
    def kappa(self, w, component):
        if w <= self.d - self.half_size[component]:
            return 1 + (self.kappa_max - 1) * (1 - (w + self.half_size[component]) / self.d)**self.m
        elif self.half_size[component] - self.d <= w:
            return 1 + (self.kappa_max - 1) * (1 + (w - self.half_size[component]) / self.d)**self.m
        else:
            return 1.0
        
    def alpha(self, w, component):
        if w <= self.d - self.half_size[component]:
            return self.alpha_max * ((w + self.half_size[component]) / self.d)**self.m_alpha
        elif self.half_size[component] - self.d <= w:
            return self.alpha_max * ((self.half_size[component] - w) / self.d)**self.m_alpha
        else:
            return 0.0
       
    def b(self, w, component):
        exponent = -(self.sigma(w, component) / self.kappa(w, component) + self.alpha(w, component)) \
                * self.dt / epsilon0
        return exp(exponent)
        
    def c(self, w, component):
        numerator = self.sigma(w, component) * (self.b(w, component) - 1.)
        if numerator == 0:
            return 0.0
        denominator = (self.sigma(w, component) + \
                       self.kappa(w, component) * self.alpha(w, component)) * self.kappa(w, component)
        return numerator / denominator
           
    def get_pointwise_material_ex(self, idx, co):
        by = self.b(co[1], 1)
        bz = self.b(co[2], 2)
        cy = self.c(co[1], 1)
        cz = self.c(co[2], 2)
        kappay = self.kappa(co[1], 1)
        kappaz = self.kappa(co[2], 2)
        pw_obj = CPMLEx(idx, self.epsilon_r, by, bz, cy, cz, kappay, kappaz)
        return pw_obj
    
    def get_pointwise_material_ey(self, idx, co):
        bz = self.b(co[2], 2)
        bx = self.b(co[0], 0)
        cz = self.c(co[2], 2)
        cx = self.c(co[0], 0)
        kappaz = self.kappa(co[2], 2)
        kappax = self.kappa(co[0], 0)
        pw_obj = CPMLEy(idx, self.epsilon_r, bz, bx, cz, cx, kappaz, kappax)
        return pw_obj
    
    def get_pointwise_material_ez(self, idx, co):
        bx = self.b(co[0], 0)
        by = self.b(co[1], 1)
        cx = self.c(co[0], 0)
        cy = self.c(co[1], 1)
        kappax = self.kappa(co[0], 0)
        kappay = self.kappa(co[1], 1)
        pw_obj = CPMLEz(idx, self.epsilon_r, bx, by, cx, cy, kappax, kappay)
        return pw_obj
    
    def get_pointwise_material_hx(self, idx, co):
        by = self.b(co[1], 1)
        bz = self.b(co[2], 2)
        cy = self.c(co[1], 1)
        cz = self.c(co[2], 2)
        kappay = self.kappa(co[1], 1)
        kappaz = self.kappa(co[2], 2)
        pw_obj = CPMLHx(idx, self.mu_r, by, bz, cy, cz, kappay, kappaz)
        return pw_obj
    
    def get_pointwise_material_hy(self, idx, co):
        bz = self.b(co[2], 2)
        bx = self.b(co[0], 0)
        cz = self.c(co[2], 2)
        cx = self.c(co[0], 0)
        kappaz = self.kappa(co[2], 2)
        kappax = self.kappa(co[0], 0)
        pw_obj = CPMLHy(idx, self.mu_r, bz, bx, cz, cx, kappaz, kappax)
        return pw_obj
    
    def get_pointwise_material_hz(self, idx, co):
        bx = self.b(co[0], 0)
        by = self.b(co[1], 1)
        cx = self.c(co[0], 0)
        cy = self.c(co[1], 1)
        kappax = self.kappa(co[0], 0)
        kappay = self.kappa(co[1], 1)
        pw_obj = CPMLHz(idx, self.mu_r, bx, by, cx, cy, kappax, kappay)
        return pw_obj
    
    