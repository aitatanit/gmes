# -*- coding: utf-8 -*-
# cython: boundscheck=False
# cython: wraparound=False

from __future__ import division
from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass

from cmath import exp as cexp
from collections import Sequence
from copy import deepcopy
from math import sqrt, sin, cos, tanh, exp, pi
from numpy import array, inf, empty, zeros
import numpy as np

from pygeom import Material
from pw_material import *
from constant import c0


class Dummy(Material):
    """A dummy material type which dosen't update the field component.
    
    """
    def __init__(self, eps_inf=1, mu_inf=1):
        Material.__init__(self, eps_inf, mu_inf)
        
    def init(self, space, param=None):
        pass

    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "dummy object"
        print " " * indent, 
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyExCmplx()
            pw_param = DummyElectricParamCmplx()
        else:
            pw_obj = DummyExReal()
            pw_param = DummyElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyEyCmplx()
            pw_param = DummyElectricParamCmplx()
        else:
            pw_obj = DummyEyReal()
            pw_param = DummyElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyEzCmplx()
            pw_param = DummyElectricParamCmplx()
        else:
            pw_obj = DummyEzReal()
            pw_param = DummyElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyHxCmplx()
            pw_param = DummyMagneticParamCmplx()
        else:
            pw_obj = DummyHxReal()
            pw_param = DummyMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyHyCmplx()
            pw_param = DummyMagneticParamCmplx()
        else:
            pw_obj = DummyHyReal()
            pw_param = DummyMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DummyHzCmplx()
            pw_param = DummyMagneticParamCmplx()
        else:
            pw_obj = DummyHzReal()
            pw_param = DummyMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj


class Const(Material):
    """A material type which sets the field to the given value.
    
    """
    def __init__(self, value=0, eps_inf=1, mu_inf=1):
        """Arguments:
            value -- field value
            eps_inf -- frequency independent permittivity
            mu_inf -- frequency independent permeability
        
        """
        Material.__init__(self, eps_inf, mu_inf)

        if type(value) is complex:
            self.value = value
        else:
            self.value = float(value)

    def __getstate__(self):
        d = Material.__getstate__(self)
        d['value'] = self.value
        return d

    def __setstate__(self, d):
        Material.__setstate(self, d)
        self.value = d['value']
        
    def init(self, space, param=None):
        pass

    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "const object"
        print " " * indent, 
        print "value:", self.value,
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ConstExCmplx()
            pw_param = ConstElectricParamCmplx()
        else:
            pw_obj = ConstExReal()
            pw_param = ConstElectricParamReal()
        
        pw_param.value = self.value
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ConstEyCmplx()
            pw_param = ConstElectricParamCmplx()
        else:
            pw_obj = ConstEyReal()
            pw_param = ConstElectricParamReal()
        
        pw_param.value = self.value
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ConstEzCmplx()
            pw_param = ConstElectricParamCmplx()
        else:
            pw_obj = ConstEzReal()
            pw_param = ConstElectricParamReal()
        
        pw_param.value = self.value
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ConstHxCmplx()
            pw_param = ConstMagneticParamCmplx()
        else:
            pw_obj = ConstHxReal()
            pw_param = ConstMagneticParamReal()

        pw_param.value = self.value
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ConstHyCmplx()
            pw_param = ConstMagneticParamCmplx()
        else:
            pw_obj = ConstHyReal()
            pw_param = ConstMagneticParamReal()

        pw_param.value = self.value
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = ConstHzCmplx()
            pw_param = ConstMagneticParamCmplx()
        else:
            pw_obj = ConstHzReal()
            pw_param = ConstMagneticParamReal()

        pw_param.value = self.value
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj


class Dielectric(Material):
    """Representation of non-dispersive isotropic dielectric medium.
        
    """
    def __init__(self, eps_inf=1, mu_inf=1):
        """Arguments:
            eps_inf -- frequency independent permittivity
            mu_inf -- frequency independent permeability
        
        """
        Material.__init__(self, eps_inf, mu_inf)

    def init(self, space, param=None):
        pass

    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "dielectric"
        print " " * indent, 
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf

    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricExCmplx()
            pw_param = DielectricElectricParamCmplx()
        else:
            pw_obj = DielectricExReal()
            pw_param = DielectricElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricEyCmplx()
            pw_param = DielectricElectricParamCmplx()
        else:
            pw_obj = DielectricEyReal()
            pw_param = DielectricElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricEzCmplx()
            pw_param = DielectricElectricParamCmplx()
        else:
            pw_obj = DielectricEzReal()
            pw_param = DielectricElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricHxCmplx()
            pw_param = DielectricMagneticParamCmplx()
        else:
            pw_obj = DielectricHxReal()
            pw_param = DielectricMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricHyCmplx()
            pw_param = DielectricMagneticParamCmplx()
        else:
            pw_obj = DielectricHyReal()
            pw_param = DielectricMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DielectricHzCmplx()
            pw_param = DielectricMagneticParamCmplx()
        else:
            pw_obj = DielectricHzReal()
            pw_param = DielectricMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj


class Compound(object):
    """Represent the compound material.
    
    Compound material means that its electric permittivity and 
    magnetic permeability refers to the underneath material. 
    
    """
    pass


class Pml(Material, Compound):
    """Base class of PML materials.
    
    Attributes:
        d -- thickness of PML medium
        half_size -- a tuple that represents the half of the outer volume
        dt -- time differential
        dw -- tuple of space differentials
        sigma_opt -- optimal conductivity
        initialized -- initialization semaphore 
        
    """
    def __init__(self, eps_inf, mu_inf):
        Material.__init__(self, eps_inf, mu_inf)
        self.initialized = False

    def __getstate__(self):
        d = Material.__getstate__(self)
        d['initialized'] = self.initialized

        if self.initialized:
            d['center'] = self.center
            d['half_size'] = self.half_size
            d['d'] = self.d
            d['dt'] = self.dt
            d['dw'] = self.dw
            d['sigma_max'] = self.sigma_max

        return d
    
    def __setstate__(self, d):
        Material.__setstate__(self, d)
        
        if d['initialized']:
            self.center.setfield(d['center'])
            self.half_size.setfield(d['half_size'])
            self.d = d['d']
            self.dt = d['dt']
            self.dw = d['dw'].copy()
            self.sigma_max.setfield(d['sigma_max'])
        
    def init(self, space, param):
        """
        The thickness of PML layer and size of the Shell instance
        which contain the PML, are required. Also, the 
        differential of space and time should get from the space 
        instance.
        
        param: (center, half_size, thickness) of the the Shell.

        """
        self.center = np.array(param[0], np.double)
        self.half_size = array(param[1], np.double)
        self.d = param[2]

        for i in range(3):
            if self.half_size[i] < self.d:
                self.half_size[i] = np.inf
        
        self.dt = space.dt
        self.dw = array(space.dr, np.double)
        self.sigma_max = self.sigma_max_ratio * self.get_sigma_opt()
        
        self.initialized = True
        
    def get_sigma_opt(self):
        """Calculate the optimal value of conductivity.
        
        """
        eta = sqrt(self.mu_inf / self.eps_inf)
        return 0.8 * (self.m + 1) / (eta * self.dw)
    
    def sigma(self, w, component):
        """Polynomial grading of conductivity.
        
        """
        w -= self.center[component]
        half_size = self.half_size[component]
        
        if w <= self.d - half_size:
            return self.sigma_max[component] * (1 - (half_size + w) / self.d)**self.m
        elif half_size - self.d <= w:
            return self.sigma_max[component] * (1 - (half_size - w) / self.d)**self.m
        else:
            return 0
        
    def kappa(self, w, component):
        """Polynomial grading of kappa.
        
        """
        w -= self.center[component]
        half_size = self.half_size[component]
        if w <= self.d - half_size:
            return 1 + (self.kappa_max - 1) * (1 - (half_size + w) / self.d)**self.m
        elif half_size - self.d <= w:
            return 1 + (self.kappa_max - 1) * (1 - (half_size - w) / self.d)**self.m
        else:
            return 1
        
        
class Upml(Pml):
    """Form Uniaxial Perfectly Matched Layer (UPML).
    
    This class implements UPML represented in
    
    S. D. Gedney, "An anisotropic perfectly matched layer-
    absorbing medium for the truncation of FDTD lattices," IEEE 
    Trans. Antennas Propag. 44, 1630-1639 (1996).
    
    Attributes:
        eps_inf -- the permittivity for incident mode impinging on the PML boundary with infinite frequency. default 1
        mu_inf -- the permeability for incident mode impinging on the PML boundary with infinite frequency. default 1
        m -- degree of kappa. default 3.6
        kappa_max -- maximum of kappa. default 4.6
        sigma_max_ratio -- the ratio between sigma_max and sigma_opt. default 0.745

    """    
    def __init__(self, eps_inf=1, mu_inf=1, m=3.6, kappa_max=4.6, sigma_max_ratio=.745):
        Pml.__init__(self, eps_inf, mu_inf)

        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.sigma_max_ratio = float(sigma_max_ratio)

    def __getstate__(self):
        d = Pml.__getstate__(self)
        d['m'] = self.m
        d['kappa_max'] = self.kappa_max
        d['sigma_max_ratio'] = self.sigma_max_ratio
        return d

    def __setstate__(self, d):
        Pml.__setstate__(self, d)
        self.m = d['m']
        self.kappa_max = d['kappa_max']
        self.sigma_max_ratio = d['sigma_max_ratio']

    def display_info(self, indent=0):
        """Display the parameter values.

        Override PML.display_info.
        
        """
        print ' ' * indent, 'UPML'
        print ' ' * indent, 
        print 'frequency independent permittivity:', self.eps_inf,
        print 'frequency independent permeability:', self.mu_inf

        print ' ' * indent,
        print 'sigma_max:', self.sigma_max,
        print 'm:', self.m,
        print 'kappa_max:', self.kappa_max
        
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
        if cmplx:
            pw_obj = UpmlExCmplx()
            pw_param = UpmlElectricParamCmplx()
        else:
            pw_obj = UpmlExReal()
            pw_param = UpmlElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.c1 = self.c1(coords[1], 1)
        pw_param.c2 = self.c2(coords[1], 1)
        pw_param.c3 = self.c3(coords[2], 2)
        pw_param.c4 = self.c4(coords[2], 2)
        pw_param.c5 = self.c5(coords[0], 0)
        pw_param.c6 = self.c6(coords[0], 0)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = UpmlEyCmplx()
            pw_param = UpmlElectricParamCmplx()
        else:
            pw_obj = UpmlEyReal()
            pw_param = UpmlElectricParamReal()
         
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.c1 = self.c1(coords[2], 2)
        pw_param.c2 = self.c2(coords[2], 2)
        pw_param.c3 = self.c3(coords[0], 0)
        pw_param.c4 = self.c4(coords[0], 0)
        pw_param.c5 = self.c5(coords[1], 1)
        pw_param.c6 = self.c6(coords[1], 1)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = UpmlEzCmplx()
            pw_param = UpmlElectricParamCmplx()
        else:
            pw_obj = UpmlEzReal()
            pw_param = UpmlElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.c1 = self.c1(coords[0], 0)
        pw_param.c2 = self.c2(coords[0], 0)
        pw_param.c3 = self.c3(coords[1], 1)
        pw_param.c4 = self.c4(coords[1], 1)
        pw_param.c5 = self.c5(coords[2], 2)
        pw_param.c6 = self.c6(coords[2], 2)

        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = UpmlHxCmplx()
            pw_param = UpmlMagneticParamCmplx()
        else:
            pw_obj = UpmlHxReal()
            pw_param = UpmlMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_param.c1 = self.c1(coords[1], 1)
        pw_param.c2 = self.c2(coords[1], 1)
        pw_param.c3 = self.c3(coords[2], 2)
        pw_param.c4 = self.c4(coords[2], 2)
        pw_param.c5 = self.c5(coords[0], 0)
        pw_param.c6 = self.c6(coords[0], 0)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = UpmlHyCmplx()
            pw_param = UpmlMagneticParamCmplx()
        else:
            pw_obj = UpmlHyReal()
            pw_param = UpmlMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
            
        pw_param.c1 = self.c1(coords[2], 2)
        pw_param.c2 = self.c2(coords[2], 2)
        pw_param.c3 = self.c3(coords[0], 0)
        pw_param.c4 = self.c4(coords[0], 0)
        pw_param.c5 = self.c5(coords[1], 1)
        pw_param.c6 = self.c6(coords[1], 1)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = UpmlHzCmplx()
            pw_param = UpmlMagneticParamCmplx()
        else:
            pw_obj = UpmlHzReal()
            pw_param = UpmlMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_param.c1 = self.c1(coords[0], 0)
        pw_param.c2 = self.c2(coords[0], 0)
        pw_param.c3 = self.c3(coords[1], 1)
        pw_param.c4 = self.c4(coords[1], 1)
        pw_param.c5 = self.c5(coords[2], 2)
        pw_param.c6 = self.c6(coords[2], 2)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    
class Cpml(Pml):
    """Form Complex Frequency Shifted (CFS) Perfectly Matched Layer (PML).
    
    This class implements CFS PML represented in
    S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,
    Computational Electrodynamics: The Finite-Difference Time-Domain Method," 
    Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
    2005, pp. 273-328.
    
    Attributes:
        eps_inf -- the permittivity for incident mode impinging on the PML boundary with infinite frequency. default 1
        mu_inf -- the permeability of incident mode impinging on the PML boundary with infinite frequency. default 1
        m -- degree of kappa. default 3.4
        kappa_max -- maximum of kappa. default 1
        m_a -- degree of a. default 4.8
        a_max -- maximum of a. default 0.8
        sigma_max_ratio -- default 0.65
    
    """
    def __init__(self, eps_inf=1, mu_inf=1, m=3.4, kappa_max=1, m_a=4.8, a_max=0.8, sigma_max_ratio=0.65):
        Pml.__init__(self, eps_inf, mu_inf)

        self.m = float(m)
        self.kappa_max = float(kappa_max)
        self.m_a = float(m_a)
        self.a_max = float(a_max)
        self.sigma_max_ratio = float(sigma_max_ratio)
        
    def __getstate__(self):
        d = Pml.__getstate__(self)
        d['m'] = self.m
        d['kappa_max'] = self.kappa_max
        d['m_a'] = self.m_a
        d['a_max'] = self.a_max
        d['sigma_max_ratio'] = self.sigma_max_ratio
        return d

    def __setstate__(self, d):
        Pml.__setstate__(self, d)
        self.m = d['m']
        self.kappa_max = d['kappa_max']
        self.m_a = d['m_a']
        self.a_max = d['a_max']
        self.sigma_max_ratio = d['sigma_max_ratio']
        
    def display_info(self, indent=0):
        """Display the parameter values.

        Override PML.display_info.

        """
        print ' ' * indent, 'CPML'
        print ' ' * indent, 
        print 'frequency independent permittivity:', self.eps_inf,
        print 'frequency independent permeability:', self.mu_inf
        
        print ' ' * indent,
        print 'sigma_max:', self.sigma_max,
        print 'm:', self.m,
        print 'kappa_max:', self.kappa_max,
        print 'm_a:', self.m_a,
        print 'a_max:', self.a_max

    def a(self, w, component):
        w -= self.center[component]
        half_size = self.half_size[component]

        if w <= self.d - half_size:
            return self.a_max * ((half_size + w) / self.d)**self.m_a
        elif half_size - self.d <= w:
            return self.a_max * ((half_size - w) / self.d)**self.m_a
        else:
            return 0
        
    def b(self, w, component):
        exponent = -(self.sigma(w, component) / self.kappa(w, component) + self.a(w, component)) * self.dt
        return exp(exponent)
        
    def c(self, w, component):
        sigma = self.sigma(w, component)
        kappa = self.kappa(w, component)
        numerator = sigma * (self.b(w, component) - 1)
        denominator = (sigma + kappa * self.a(w, component)) * kappa
        if denominator:
            return numerator / denominator
        else:
            return 0
    
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = CpmlExCmplx()
            pw_param = CpmlElectricParamCmplx()
        else:
            pw_obj = CpmlExReal()
            pw_param = CpmlElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.b1 = self.b(coords[1], 1)
        pw_param.b2 = self.b(coords[2], 2)
        pw_param.c1 = self.c(coords[1], 1)
        pw_param.c2 = self.c(coords[2], 2)
        pw_param.kappa1 = self.kappa(coords[1], 1)
        pw_param.kappa2 = self.kappa(coords[2], 2)

        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = CpmlEyCmplx()
            pw_param = CpmlElectricParamCmplx()
        else:
            pw_obj = CpmlEyReal()
            pw_param = CpmlElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.b1 = self.b(coords[2], 2)
        pw_param.b2 = self.b(coords[0], 0)
        pw_param.c1 = self.c(coords[2], 2)
        pw_param.c2 = self.c(coords[0], 0)
        pw_param.kappa1 = self.kappa(coords[2], 2)
        pw_param.kappa2 = self.kappa(coords[0], 0)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = CpmlEzCmplx()
            pw_param = CpmlElectricParamCmplx()
        else:
            pw_obj = CpmlEzReal()
            pw_param = CpmlElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.b1 = self.b(coords[0], 0)
        pw_param.b2 = self.b(coords[1], 1)
        pw_param.c1 = self.c(coords[0], 0)
        pw_param.c2 = self.c(coords[1], 1)
        pw_param.kappa1 = self.kappa(coords[0], 0)
        pw_param.kappa2 = self.kappa(coords[1], 1)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = CpmlHxCmplx()
            pw_param = CpmlMagneticParamCmplx()
        else:
            pw_obj = CpmlHxReal()
            pw_param = CpmlMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_param.b1 = self.b(coords[1], 1)
        pw_param.b2 = self.b(coords[2], 2)
        pw_param.c1 = self.c(coords[1], 1)
        pw_param.c2 = self.c(coords[2], 2)
        pw_param.kappa1 = self.kappa(coords[1], 1)
        pw_param.kappa2 = self.kappa(coords[2], 2)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = CpmlHyCmplx()
            pw_param = CpmlMagneticParamCmplx()
        else:
            pw_obj = CpmlHyReal()
            pw_param = CpmlMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_param.b1 = self.b(coords[2], 2)
        pw_param.b2 = self.b(coords[0], 0)
        pw_param.c1 = self.c(coords[2], 2)
        pw_param.c2 = self.c(coords[0], 0)
        pw_param.kappa1 = self.kappa(coords[2], 2)
        pw_param.kappa2 = self.kappa(coords[0], 0)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = CpmlHzCmplx()
            pw_param = CpmlMagneticParamCmplx()
        else:
            pw_obj = CpmlHzReal()
            pw_param = CpmlMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_param.b1 = self.b(coords[0], 0)
        pw_param.b2 = self.b(coords[1], 1)
        pw_param.c1 = self.c(coords[0], 0)
        pw_param.c2 = self.c(coords[1], 1)
        pw_param.kappa1 = self.kappa(coords[0], 0)
        pw_param.kappa2 = self.kappa(coords[1], 1)
        pw_obj.attach(idx, pw_param)
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
        

class DcpAde(Dielectric):
    def __init__(self, eps_inf=1, mu_inf=1, sigma=0, dps=(), cps=()):
        """
        eps_inf: The (frequency-independent) relative permittivity. Default is 1.
        mu_inf: The (frequency-independent) relative permeability. Default is 1.
        sigma: The (frequency-independent) isotropic conductivity. Default is 0.
        dps: list of Drude poles. Default is ().
        cps: list of critical points. Default is ().
        
        """
        Dielectric.__init__(self, eps_inf, mu_inf)
        self.sigma = float(sigma) # instant conductivity
        self.dps = tuple(dps) # tuple of Drude poles
        self.cps = tuple(cps) # tuple of critical points
        self.initialized = False
        
    def __getstate__(self):
        d = Dielectric.__getstate__(self)
        d['sigma'] = self.sigma
        d['dps'] = self.dps
        d['cps'] = self.cps
        d['initialized'] = self.initialized

        if self.initialized:
            d['dt'] = self.dt
            d['a'] = self.a
            d['b'] = self.b
            d['c'] = self.c
            
        return d

    def __setstate__(self, d):
        Dielectric.__setstate__(self, d)
        
        self.sigma = d['sigma']
        self.dps = deepcopy(d['dps'])
        self.cps = deepcopy(d['cps'])
        self.initialized = d['initialized']

        if self.initialized:
            self.dt = d['dt']
            self.a = d['a'].copy()
            self.b = d['b'].copy()
            self.c = d['c'].copy()
            
    def init(self, space, param=None):
        self.dt = space.dt
        
        # parameters for the ADE of the Drude model
        self.a = empty((len(self.dps),3), np.double)
        for i in xrange(len(self.dps)):
            pole = self.dps[i]
            denom = float(self.dt * pole.gamma + 2)
            self.a[i,0] = (self.dt * pole.gamma - 2) / denom
            self.a[i,1] = 4 / denom
            self.a[i,2] = 0.5 * (self.dt * pole.omega)**2 / denom

        # parameters for the ADE of critical points model
        self.b = empty((len(self.cps),5), np.double)
        for i in xrange(len(self.cps)):
            pnt = self.cps[i]
            denom = (self.dt * pnt.gamma + 2)**2 + (self.dt * pnt.omega)**2
            self.b[i,0] = -((self.dt * pnt.gamma - 2)**2 + (self.dt * pnt.omega)**2) / denom
            self.b[i,1] = 2 * (4 - self.dt**2 * (pnt.gamma**2 + pnt.omega**2)) / denom
            self.b[i,2] = 2 * self.dt * pnt.amp * pnt.omega * (self.dt * cos(pnt.phi) * pnt.omega + sin(pnt.phi) * (2 - self.dt * pnt.gamma)) / denom
            self.b[i,3] = 4 * self.dt**2 * pnt.amp * pnt.omega * (cos(pnt.phi) * pnt.omega - sin(pnt.phi) * pnt.gamma) / denom
            self.b[i,4] = 2 * self.dt * pnt.amp * pnt.omega * (self.dt * cos(pnt.phi) * pnt.omega - sin(pnt.phi) * (2 + self.dt * pnt.gamma)) / denom
            
        # parameters for the electric field update equations.
        self.c = empty(4, np.double)
        denom = 0.5 * self.dt * self.sigma + sum(self.a[:,2]) + sum(self.b[:,4]) + self.eps_inf
        self.c[0] = self.dt / denom
        self.c[1] = 1 / denom
        self.c[2] = -(sum(self.a[:,2]) + sum(self.b[:,2])) / denom
        self.c[3] = -(0.5 * self.dt * self.sigma + 2 * sum(self.a[:,2]) + sum(self.b[:,3]) - self.eps_inf) / denom

        self.initialized = True
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude-critical points dispersive media"
        print " " * indent,
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf,
        print "conductivity:", self.sigma
        
        print " " * indent, "Drude pole(s):"
        for i in self.dps:
            i.display_info(indent+4)
        print " " * indent, "critical point(s):"
        for i in self.cps:
            i.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpAdeExCmplx()
            pw_param = DcpAdeElectricParamCmplx()
        else:
            pw_obj = DcpAdeExReal()
            pw_param = DcpAdeElectricParamReal()
                
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.b, self.c)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpAdeEyCmplx()
            pw_param = DcpAdeElectricParamCmplx()
        else:
            pw_obj = DcpAdeEyReal()
            pw_param = DcpAdeElectricParamReal()
                
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.b, self.c)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpAdeEzCmplx()
            pw_param = DcpAdeElectricParamCmplx()
        else:
            pw_obj = DcpAdeEzReal()
            pw_param = DcpAdeElectricParamReal()
                
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.b, self.c)
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpAdeHxCmplx()
            pw_param = DcpAdeMagneticParamCmplx()
        else:
            pw_obj = DcpAdeHxReal()
            pw_param = DcpAdeMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpAdeHyCmplx()
            pw_param = DcpAdeMagneticParamCmplx()
        else:
            pw_obj = DcpAdeHyReal()
            pw_param = DcpAdeMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpAdeHzCmplx()
            pw_param = DcpAdeMagneticParamCmplx()
        else:
            pw_obj = DcpAdeHzReal()
            pw_param = DcpAdeMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj


class DcpPlrc(Dielectric):
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
    def __init__(self, eps_inf=1, mu_inf=1, sigma=0, dps=(), cps=()):
        """
        eps_inf: The (frequency-independent) relative permittivity. Default is 1.
        mu_inf: The (frequency-independent) relative permeability. Default is 1.
        sigma: The (frequency-independent) isotropic conductivity. Default is 0.
        dps: list of Drude poles. Default is ().
        cps: list of critical points. Default is ().
        
        """
        Dielectric.__init__(self, eps_inf=eps_inf, mu_inf=mu_inf)
        self.sigma = float(sigma) # instant conductivity
        self.dps = tuple(dps) # tuple of Drude poles
        self.cps = tuple(cps) # tuple of critical points
        self.initialized = False
        
    def __getstate__(self):
        d = Dielectric.__getstate__(self)
        d['sigma'] = self.sigma
        d['dps'] = self.dps
        d['cps'] = self.cps
        d['initialized'] = self.initialized

        if self.initialized:
            d['dt'] = self.dt
            d['a'] = self.a
            d['b'] = self.b
            d['c'] = self.c

        return d
    
    def __setstate__(self, d):
        Dielectric.__setstate__(self, d)
        self.sigma = d['sigma']
        self.dps = deepcopy(d['dps'])
        self.cps = deepcopy(d['cps'])
        self.initialized = d['initialized']

        if self.initialized:
            self.dt = d['dt']
            self.a = d['a'].copy()
            self.b = d['b'].copy()
            self.c = d['c'].copy()
        
    def init(self, space, param=None):
        self.dt = space.dt

        # parameters of the recursion relation for the Drude pole recursive accumulator.
        self.a = empty((len(self.dps), 3), np.double)
        for i, pole in enumerate(self.dps):
            self.a[i,0] = self.delta_chi_dp_0(pole) - self.delta_xi_dp_0(pole)
            self.a[i,1] = self.delta_xi_dp_0(pole)
            self.a[i,2] = exp(-pole.gamma * self.dt)
        
        # parameters of the recursion relation for the critical point recursive accumulator.
        self.b = empty((len(self.cps), 3), complex)
        for i, pnt in enumerate(self.cps):
            self.b[i,0] = self.delta_chi_cp_0(pnt) - self.delta_xi_cp_0(pnt)
            self.b[i,1] = self.delta_xi_cp_0(pnt)
            self.b[i,2] = cexp(-self.dt * (pnt.gamma + 1j * pnt.omega))
            
        # parameters for the electric field update equations.
        chi_0 = (sum(map(self.chi_dp_0, self.dps) + 
                     map(self.chi_cp_0, self.cps)) + 0j).real

        xi_0 = (sum(map(self.xi_dp_0, self.dps) + 
                    map(self.xi_cp_0, self.cps)) + 0j).real
        
        self.c = empty(3, np.double)
        denom = self.eps_inf - xi_0 + chi_0
        self.c[0] = self.dt / denom
        self.c[1] = (self.eps_inf - xi_0) / denom
        self.c[2] = 1 / denom

        self.initialized = True
        
    def chi_dp_0(self, dp):
        omega = dp.omega
        gamma = dp.gamma
        gdt = dp.gamma * self.dt
        return (omega / gamma)**2 * (exp(-gdt) + gdt - 1)

    def xi_dp_0(self, dp):
        chi_dp_0 = self.chi_dp_0(dp)
        omega = dp.omega
        gamma = dp.gamma
        gdt = dp.gamma * self.dt
        return chi_dp_0 * (1 / (1 - exp(gdt)) + 1 / gdt) + (omega / gamma)**2 * (gdt / 2 / tanh(gdt / 2) - 1)

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
        return 2j * cp.amp * cp.omega * cexp(1j * cp.phi) * (1 - cexp(-self.dt * go)) / go
    
    def xi_cp_0(self, cp):
        dtgo = self.dt * (cp.gamma + 1j * cp.omega)
        chi_cp_0 = self.chi_cp_0(cp)
        return chi_cp_0 * (1 / (1 - cexp(dtgo)) + 1 / dtgo)
    
    def delta_chi_cp_0(self, cp):
        go = cp.gamma + 1j * cp.omega
        return 2j * cp.amp * cp.omega * cexp(1j * cp.phi) * (1 - cexp(-self.dt * go))**2 / go
    
    def delta_xi_cp_0(self, cp):
        dtgo = self.dt * (cp.gamma + 1j * cp.omega)
        delta_chi_cp_0 = self.delta_chi_cp_0(cp)
        return delta_chi_cp_0 * (1 / (1 - cexp(dtgo)) + 1 / dtgo)
    
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude-critical points dispersive media"
        print " " * indent,
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf,
        print "conductivity:", self.sigma
        
        print " " * indent, "Drude pole(s):"
        for i in self.dps:
            i.display_info(indent+4)
        print " " * indent, "critical point(s):"
        for i in self.cps:
            i.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpPlrcExCmplx()
            pw_param = DcpPlrcElectricParamCmplx()
        else:
            pw_obj = DcpPlrcExReal()
            pw_param = DcpPlrcElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.b, self.c)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpPlrcEyCmplx()
            pw_param = DcpPlrcElectricParamCmplx()
        else:
            pw_obj = DcpPlrcEyReal()
            pw_param = DcpPlrcElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.b, self.c)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpPlrcEzCmplx()
            pw_param = DcpPlrcElectricParamCmplx()
        else:
            pw_obj = DcpPlrcEzReal()
            pw_param = DcpPlrcElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.b, self.c)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpPlrcHxCmplx()
            pw_param = DcpPlrcMagneticParamCmplx()
        else:
            pw_obj = DcpPlrcHxReal()
            pw_param = DcpPlrcMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpPlrcHyCmplx()
            pw_param = DcpPlrcMagneticParamCmplx()
        else:
            pw_obj = DcpPlrcHyReal()
            pw_param = DcpPlrcMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DcpPlrcHzCmplx()
            pw_param = DcpPlrcMagneticParamCmplx()
        else:
            pw_obj = DcpPlrcHzReal()
            pw_param = DcpPlrcMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    
class DcpRc(DcpPlrc):
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
        return 0
    
    def xi_cp_0(self, cp):
        return 0j

    def delta_xi_dp_0(self, dp):
        return 0
    
    def delta_xi_cp_0(self, cp):
        return 0j
        
        
class Drude(Dielectric):
    """
    The auxiliary differential equation implementation of the Drude model based on the 
    following article.
    * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
    Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006.
    
    """
    def __init__(self, eps_inf=1, mu_inf=1, sigma=0, dps=()):
        """
        Arguments:
            eps_inf: The (frequency-independent) relative permittivity. Default is 1.
            mu_inf: The (frequency-independent) relative permeability. Default is 1.
            sigma: The (frequency-independent) isotropic conductivity. Default is 0.
            dps: list of Drude poles. Default is ().
            
        """
        Dielectric.__init__(self, eps_inf=eps_inf, mu_inf=mu_inf)
        self.sigma = float(sigma)
        self.dps = tuple(dps)
        self.initialized = False

    def __getstate__(self):
        d = Dielectric.__getstate__(self)
        d['sigma'] = self.sigma
        d['dps'] = self.dps
        d['initialized'] = self.initialized

        if self.initialized:
            d['dt'] = self.dt
            d['a'] = self.a
            d['c'] = self.c

        return d
    
    def __setstate__(self, d):
        Dielectric.__setstate__(self, d)
        self.sigma = d['sigma']
        self.dps = deepcopy(d['dps'])
        self.initialized = d['initialized']

        if self.initialized:
            self.dt = d['dt']
            self.a = d['a'].copy()
            self.c = d['c'].copy()

    def init(self, space, param=None):
        self.dt = space.dt
        
        # parameters for the ADE of the Drude model.
        self.a = empty((len(self.dps), 3), np.double)
        for i in xrange(len(self.dps)):
            pole = self.dps[i]
            denom = self.dt * pole.gamma + 2.
            self.a[i,0] = (self.dt * pole.gamma - 2) / denom
            self.a[i,1] = 4 / denom
            self.a[i,2] = 2 * (self.dt * pole.omega)**2 / denom
        
        # parameters for the electric field update equations.
        self.c = empty(3, np.double)
        denom = 2. * self.eps_inf + self.dt * self.sigma 
        self.c[0] = 2 * self.dt / denom
        self.c[1] = -2 / denom
        self.c[2] = (2 * self.eps_inf - self.dt * self.sigma) / denom

        self.initialized = True
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Drude dispersion media"
        print " " * indent, 
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf,
        print "conductivity:", self.sigma
        
        print " "* indent, "Drude pole(s):"
        for p in self.dps:
            p.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeExCmplx()
            pw_param = DrudeElectricParamCmplx()
        else:
            pw_obj = DrudeExReal()
            pw_param = DrudeElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.c)

        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeEyCmplx()
            pw_param = DrudeElectricParamCmplx()
        else:
            pw_obj = DrudeEyReal()
            pw_param = DrudeElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.c)

        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeEzCmplx()
            pw_param = DrudeElectricParamCmplx()
        else:
            pw_obj = DrudeEzReal()
            pw_param = DrudeElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.c)

        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeHxCmplx()
            pw_param = DrudeMagneticParamCmplx()
        else:
            pw_obj = DrudeHxReal()
            pw_param = DrudeMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeHyCmplx()
            pw_param = DrudeMagneticParamCmplx()
        else:
            pw_obj = DrudeHyReal()
            pw_param = DrudeMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = DrudeHzCmplx()
            pw_param = DrudeMagneticParamCmplx()
        else:
            pw_obj = DrudeHzReal()
            pw_param = DrudeMagneticParamReal()

        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj


class Lorentz(Dielectric):
    """
    The auxiliary differential equation implementation of the Lorentz model.
    
    """
    def __init__(self, eps_inf=1, mu_inf=1, sigma=0, lps=()):
        """
        Arguments:
            eps_inf: The (frequency-independent) relative permittivity. Default is 1.
            mu_inf: The (frequency-independent) relative permeability. Default is 1.
            sigma: The (frequency-independent) conductivity. Default is 0.
            lps: list of Lorentz poles. Default is ().
            
        """
        Dielectric.__init__(self, eps_inf, mu_inf)
        self.sigma = float(sigma)
        self.lps = tuple(lps)
        self.initialized = False

    def __getstate__(self):
        d = Dielectric.__setstate__(self)
        d['sigma'] = self.sigma
        d['lps'] = self.lps
        d['initialized'] = self.initialized

        if self.initialized:
            d['dt'] = self.dt
            d['a'] = self.a
            d['c'] = self.c

    def __setstate__(self, d):
        self.sigma = d['sigma']
        self.lps = deepcopy(d['lps'])
        self.initialized = d['initialized']

        if self.initialized:
            self.dt = d['dt']
            self.a = d['a'].copy()
            self.c = d['c'].coyp()
            
    def init(self, space, param=None):
        self.dt = space.dt
        
        # parameters for the ADE of the Drude model.
        self.a = empty((len(self.lps),3), np.double)
        for i in xrange(len(self.lps)):
            pole = self.lps[i]
            denom = self.dt * pole.gamma + 2.
            self.a[i,0] = (self.dt * pole.gamma - 2) / denom
            self.a[i,1] = (4 - 2 * (self.dt * pole.omega)**2) / denom
            self.a[i,2] = 2 * pole.amp * (self.dt * pole.omega)**2 / denom
        
        # parameters for the electric field update equations.
        self.c = empty(3, np.double)
        denom = 2. * self.eps_inf + self.dt * self.sigma 
        self.c[0] = 2 * self.dt / denom
        self.c[1] = -2 / denom
        self.c[2] = (2 * self.eps_inf - self.dt * self.sigma) / denom

        self.initialized = True
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "Lorentz dispersion media"
        print " " * indent, 
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf,
        print "conductivity:", self.sigma
        
        print " "* indent, "Lorentz pole(s):"
        for p in self.lps:
            p.display_info(indent+4)
        
    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzExCmplx()
            pw_param = LorentzElectricParamCmplx()
        else:
            pw_obj = LorentzExReal()
            pw_param = LorentzElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.c)
        pw_obj.attach(idx, pw_param)
        return pw_obj
    
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzEyCmplx()
            pw_param = LorentzElectricParamCmplx()
        else:
            pw_obj = LorentzEyReal()
            pw_param = LorentzElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.c)
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzEzCmplx()
            pw_param = LorentzElectricParamCmplx()
        else:
            pw_obj = LorentzEzReal()
            pw_param = LorentzElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.set(self.a, self.c)
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzHxCmplx()
            pw_param = LorentzMagneticParamCmplx()
        else:
            pw_obj = LorentzHxReal()
            pw_param = LorentzMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzHyCmplx()
            pw_param = LorentzMagneticParamCmplx()
        else:
            pw_obj = LorentzHyReal()
            pw_param = LorentzMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            pw_obj = LorentzHzCmplx()
            pw_param = LorentzMagneticParamCmplx()
        else:
            pw_obj = LorentzHzReal()
            pw_param = LorentzMagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj


class Dm2(Dielectric):
    """
    The predictor-corrector implementation of the density matrix of
    a two-level medium.
    
    """
    def __init__(self, eps_inf=1, mu_inf=1, omega=(1,), n_atom=(1,), rho30=-1, gamma=1, t1=1, t2=1, hbar=1, rtol=10e-5):
        """
        Arguments:
            eps_inf: float, optional
                The (frequency-independent) relative permittivity. Defaults to 1.
            mu_inf: float, optional
                The (frequency-independent) relative permeability. Defaults to 1.
            omega: tuple, optional
                The angular resonance frequencies of atomic transition from the ground level 
                to the excited level. Defaults to (1,).
            n_atom: tuple, optional
                The densities of polarizable atoms. Defaults to (1,).
            rho30: float, optional
                The initial population difference in the system. Defaults to -1.
            gamma: float, optional
                The dipole coupling coefficient. Defaults to 1.
            t1: float, optional
                The excited-state lifetime. Defaults to 1.
            t2: float, optional
                The dephasing time. Defaults to 1.
            hbar: float, optional
                The normalized reduced Planck constant. Defaults to 1.
            rtol: float, optional
                relative tolerance for solution. Default is 10e-5.

        """
        Dielectric.__init__(self, eps_inf, mu_inf)
        if isinstance(omega, (Sequence, np.ndarray)) :
            self.omega = tuple(map(float, omega))
        else:
            self.omega = (float(omega),)

        if isinstance(n_atom, (Sequence, np.ndarray)):
            self.n_atom = tuple(map(float, n_atom))
        else:
            self.n_atom = (float(n_atom),)

        self.rho30 = float(rho30)
        self.gamma = float(gamma)
        self.t1 = float(t1)
        self.t2 = float(t2)
        self.hbar = float(hbar)
        self.rtol = float(rtol)

    def __getstate__(self):
        d = Dielectric.__setstate__(self)
        d['omega'] = deepcopoy(self.omega)
        d['n_atom'] = deepcopy(self.n_atom)
        d['rho30'] = self.rho30
        d['gamma'] = self.gamma
        d['t1'] = self.t1
        d['t2'] = self.t2
        d['hbar'] = self.hbar
        d['rtol'] = self.rtol
        d['initialized'] = self.initialized

        if self.initialized:
            d['dt'] = self.dt

    def __setstate__(self, d):
        self.omega = deepcopy(d['omega'])
        self.n_atom = deepcopy(d['n_atom'])
        self.rho30 = d['rho30']
        self.gamma = d['gamma']
        self.t1 = d['t1']
        self.t2 = d['t2']
        self.hbar = d['hbar']
        self.rtol = d['rtol']
        self.initialized = d['initialized']

        if self.initialized:
            self.dt = d['dt']
            
    def init(self, space, param=None):
        self.dt = space.dt
        self.initialized = True
        
    def display_info(self, indent=0):
        """Display the parameter values.
        
        """
        print " " * indent, "2-level media"
        print " " * indent, 
        print "frequency independent permittivity:", self.eps_inf,
        print "frequency independent permeability:", self.mu_inf,
        print "angular frequencies of atomic transition resonance energy:", self.omega,
        print "densities of polarizable atoms:", self.n_atom,
        print "initial populaiton difference:", self.rho30,
        print "dipole coupling coefficient:", self.gamma,
        print "excited-state lifetime:", self.t1,
        print "dephasing time:", self.t2
        print "normzlied reduced Planck constant:", self.hbar
        print "relative tolerance:", self.rtol

    def get_pw_material_ex(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            raise ValueError('Dm2 class supports real fields only')
        else:
            pw_obj = Dm2ExReal()
            pw_param = Dm2ElectricParamReal()

        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
            
        pw_param.set(self.omega, self.n_atom)
        pw_param.rho30 = self.rho30
        pw_param.gamma = self.gamma
        pw_param.t1 = self.t1
        pw_param.t2 = self.t2
        pw_param.hbar = self.hbar
        pw_param.rtol = self.rtol

        pw_obj.attach(idx, pw_param)
        return pw_obj
        
    def get_pw_material_ey(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            raise ValueError('Dm2 class supports real fields only')
        else:
            pw_obj = Dm2EyReal()
            pw_param = Dm2ElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.omega = self.omega
        pw_param.rho30 = self.rho30
        pw_param.n_atom = self.n_atom
        pw_param.gamma = self.gamma
        pw_param.t1 = self.t1
        pw_param.t2 = self.t2
        pw_param.hbar = self.hbar
        pw_param.rtol = self.rtol

        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_ez(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            raise ValueError('Dm2 class supports real fields only')
        else:
            pw_obj = Dm2EzReal()
            pw_param = Dm2ElectricParamReal()
            
        if underneath is None:
            pw_param.eps_inf = self.eps_inf
        else:
            pw_param.eps_inf = underneath.eps_inf
        
        pw_param.omega = self.omega
        pw_param.rho30 = self.rho30
        pw_param.n_atom = self.n_atom
        pw_param.gamma = self.gamma
        pw_param.t1 = self.t1
        pw_param.t2 = self.t2
        pw_param.hbar = self.hbar
        pw_param.rtol = self.rtol

        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hx(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            raise ValueError('Dm2 class supports real fields only')
        else:
            pw_obj = Dm2HxReal()
            pw_param = Dm2MagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)    
        return pw_obj

    def get_pw_material_hy(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            raise ValueError('Dm2 class supports real fields only')
        else:
            pw_obj = Dm2HyReal()
            pw_param = Dm2MagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj

    def get_pw_material_hz(self, idx, coords, underneath=None, cmplx=False):
        if cmplx:
            raise ValueError('Dm2 class supports real fields only')
        else:
            pw_obj = Dm2HzReal()
            pw_param = Dm2MagneticParamReal()
            
        if underneath is None:
            pw_param.mu_inf = self.mu_inf
        else:
            pw_param.mu_inf = underneath.mu_inf
        
        pw_obj.attach(idx, pw_param)
        return pw_obj
