#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass
    
from copy import deepcopy
from math import sqrt, pi, sin, cos, exp
from cmath import exp as cexp

import numpy as np
from numpy import inf, cross, dot, ndindex
from numpy.linalg import norm
from scipy.optimize import bisect

import constant as const
from geometry import Cartesian, DefaultMedium, Shell, in_range
from fdtd import TEMzFDTD
from material import Dielectric, Cpml

# for a point source
from pw_source import PointSourceParam
from pw_source import PointSourceEx, PointSourceEy, PointSourceEz
from pw_source import PointSourceHx, PointSourceHy, PointSourceHz

# for a total-field/scattered-field and Gaussian beam source
from pw_source import TransparentElectricParam, TransparentMagneticParam
from pw_source import TransparentEx, TransparentEy, TransparentEz
from pw_source import TransparentHx, TransparentHy, TransparentHz

# 
# SrcTime: Continuous, Bandpass
# Src: PointSource, GaussianBeam, TotalFieldScatteredField
#

class SrcTime(object):
    """Time-dependent part of a source.
    
    """
    def init(self, cmplx):
        raise NotImplementedError

    def oscillator(self, time):
        raise NotImplementedError

    def display_info(self, indent=0):
        raise NotImplementedError
    

class Src(object):
    """Space-dependent part of a source.
    
    """
    def display_info(self, indent=0):
        raise NotImplementedError
    
    def init(self, geom_tree, space, cmplx):
        raise NotImplementedError

    def step(self):
        raise NotImplementedError

    def get_pw_source_ex(self, ex_field, space, geom_tree):
        raise NotImplementedError
    
    def get_pw_source_ey(self, ey_field, space, geom_tree):
        raise NotImplementedError

    def get_pw_source_ez(self, ez_field, space, geom_tree):
        raise NotImplementedError

    def get_pw_source_hx(self, hx_field, space, geom_tree):
        raise NotImplementedError

    def get_pw_source_hy(self, hy_field, space, geom_tree):
        raise NotImplementedError
    
    def get_pw_source_hz(self, hz_field, space, geom_tree):
        raise NotImplementedError


class Continuous(SrcTime):
    """Continuous (CW) source with (optional) slow turn-on and/or turn-off.
    
    """
    def __init__(self, freq, phase=0, start=0, end=inf, width=None):
        self.freq = float(freq)
        self.phase = float(phase)
        self.start = float(start)
        self.end = float(end)
        
        if width is None:
            self.width = 5 / self.freq
        else:
            self.width = float(width)

    def init(self, cmplx):
        self.cmplx = cmplx
        
    def display_info(self, indent=0):
        print ' ' * indent, 'continuous source:'
        print ' ' * indent,
        print 'frequency:', self.freq,
        print 'initial phase advance:', self.phase,
        print 'start time:', self.start,
        print 'end time:', self.end,
        print 'raising duration:', self.width
                
    def oscillator(self, time):
        ts = time - self.start
        te = self.end - time
        
        if ts < 0 or te < 0:
            return 0
        
        # Use Hanning window function to reduce the transition.
        # D. T. Prescott and N. V. Shuley, "Reducing solution time in
        # monochromatic FDTD waveguide simulations", IEEE Trans. Microwave 
        # Theory Tech., vol. 42, no. 8, pp. 1582-1584, 8. 1994.
        if ts < self.width:
            env = sin(0.5 * pi * ts / self.width)**2
        elif te < self.width:
            env = sin(0.5 * pi * te / self.width)**2
        else:
            env = 1
        
        osc = env * cexp(2j * pi * self.freq * time + 1j * self.phase)
        if self.cmplx:
            return osc
        else:
            return osc.real


class Bandpass(SrcTime):
    """a pulse source with Gaussian-envelope
    
    """
    def __init__(self, freq, fwidth, s=10, phase=0):
        self.freq = float(freq)
        self.phase = float(phase)
        self.fwidth = float(fwidth)
        self.width = 1 / self.fwidth
        self.peak_time = self.width * s
        self.cutoff = 2 * self.width * s
        
        # Makes the last_source_time as small as possible.
        while exp(-0.5 * (self.cutoff / self.width)**2) == 0:
            self.cutoff *= 0.9
        
    def init(self, cmplx):
        self.cmplx = cmplx
        
    def display_info(self, indent=0):
        print " " * indent, "bandpass source"
        print " " * indent,
        print "center frequency:", self.freq,
        print "bandwidth:", self.fwidth,
        print "peak time:", self.peak_time,
        print "cutoff:", self.cutoff
        
    def oscillator(self, time):
        tt = time - self.peak_time
        if (abs(tt) > self.cutoff): 
            return 0

        # correction factor so that current amplitude (= d(oscillator)/dt) 
        # is ~1 near the peak of the Gaussian.
        cfactor = 1.0 / (-2j * pi * self.freq)
        
        osc = cfactor * exp(-0.5 * (tt / self.width)**2) \
            * cexp(2j * pi * self.freq * time + 1j * self.phase)
        if self.cmplx:
            return osc
        else:
            return osc.real
        
        
class DifferentiatedGaussian(SrcTime):
    """a differentiated Gaussian pulse
    
    -2((t-t0)/tw)exp(-((t-t0)/tw)**2)

    """
    def __init__(self, tw, t0):
        """
        tw: half-width of pulse
        t0: dealy time
        
        """
        self.tw = float(tw)
        self.t0 = float(t0)
        
    def init(self, cmplx):
        self.cmplx = bool(cmplx)
        
    def oscillator(self, time):
        exponent = -((time - self.t0) / self.tw)**2
        osc = -2 * (time - self.t0) / self.tw * cexp(exponent)
        if self.cmplx:
            return osc
        else:
            return osc.real
        
    def display_info(self, indent=0):
        print ' ' * indent,
        print 'differentiated Gaussian pulse'
        print ' ' * indent,
        print 'half-width:', self.tw,
        print 'delay time:', self.t0


class PointSource(Src):
    def __init__(self, src_time, center, component, amp=1, filename=None):
        self.center = np.array(center, np.double)
        self.comp = component
        self.src_time = src_time
        self.amp = float(amp)
        if filename:
            self.filename = str(filename)
        else:
            self.filename = None
        
    def init(self, geom_tree, space, cmplx):
        self.geom_tree = geom_tree
        self.src_time.init(cmplx)
        
    def step(self):
        pass
    
    def display_info(self, indent=0):
        print ' ' * indent, 'point source:'
        print ' ' * indent, 'center:', self.center
        print ' ' * indent, 'component:', self.comp.str()
        print ' ' * indent, 'maximum amp.:', self.amp
        print ' ' * indent, 'source recording:', self.filename
        
        self.src_time.display_info(4)
        
    def get_pw_source_ex(self, ex_field, space, geom_tree):
        pw_src = None

        if self.comp is const.Ex or self.comp is const.Jx:
            idx = space.space_to_ex_index(*self.center)
            if in_range(idx, ex_field.shape, const.Ex):
                mat_obj, underneath = geom_tree.material_of_point(tuple(self.center))
                pw_src_param = PointSourceParam(self.src_time, 
                                                self.amp,
                                                self.comp,
                                                mat_obj.eps_inf,
                                                mat_obj.mu_inf,
                                                self.filename)
                pw_src = PointSourceEx()
                pw_src.attach(idx, pw_src_param)

        return pw_src

    def get_pw_source_ey(self, ey_field, space, geom_tree):
        pw_src = None

        if self.comp is const.Ey or self.comp is const.Jy:
            idx = space.space_to_ey_index(*self.center)
            if in_range(idx, ey_field.shape, const.Ey):
                mat_obj, underneath = geom_tree.material_of_point(tuple(self.center))
                pw_src_param = PointSourceParam(self.src_time, 
                                                self.amp,
                                                self.comp,
                                                mat_obj.eps_inf,
                                                mat_obj.mu_inf,
                                                self.filename)
                pw_src = PointSourceEy()
                pw_src.attach(idx, pw_src_param)

        return pw_src

    def get_pw_source_ez(self, ez_field, space, geom_tree):
        pw_src = None

        if self.comp is const.Ez or self.comp is const.Jz:
            idx = space.space_to_ez_index(*self.center)
            if in_range(idx, ez_field.shape, const.Ez):
                mat_obj, underneath = geom_tree.material_of_point(tuple(self.center))
                pw_src_param = PointSourceParam(self.src_time, 
                                                self.amp,
                                                self.comp,
                                                mat_obj.eps_inf,
                                                mat_obj.mu_inf,
                                                self.filename)
                pw_src = PointSourceEz()
                pw_src.attach(idx, pw_src_param)

        return pw_src

    def get_pw_source_hx(self, hx_field, space, geom_tree):
        pw_src = None

        if self.comp is const.Hx or self.comp is const.Mx:
            idx = space.space_to_hx_index(*self.center)
            if in_range(idx, hx_field.shape, const.Hx):
                mat_obj, underneath = geom_tree.material_of_point(tuple(self.center))
                pw_src_param = PointSourceParam(self.src_time, 
                                                self.amp,
                                                self.comp,
                                                mat_obj.eps_inf,
                                                mat_obj.mu_inf,
                                                self.filename)
                pw_src = PointSourceHx()
                pw_src.attach(idx, pw_src_param)

        return pw_src

    def get_pw_source_hy(self, hy_field, space, geom_tree):
        pw_src = None

        if self.comp is const.Hy or self.comp is const.My:
            idx = space.space_to_hy_index(*self.center)
            if in_range(idx, hy_field.shape, const.Hy):
                mat_obj, underneath = geom_tree.material_of_point(tuple(self.center))
                pw_src_param = PointSourceParam(self.src_time, 
                                                self.amp,
                                                self.comp,
                                                mat_obj.eps_inf,
                                                mat_obj.mu_inf,
                                                self.filename)
                pw_src = PointSourceHy()
                pw_src.attach(idx, pw_src_param)

        return pw_src

    def get_pw_source_hz(self, hz_field, space, geom_tree):
        pw_src = None

        if self.comp is const.Hz or self.comp is const.Mz:
            idx = space.space_to_hz_index(*self.center)
            if in_range(idx, hz_field.shape, const.Hz):
                mat_obj, underneath = geom_tree.material_of_point(tuple(self.center))
                pw_src_param = PointSourceParam(self.src_time, 
                                                self.amp,
                                                self.comp,
                                                mat_obj.eps_inf,
                                                mat_obj.mu_inf,
                                                self.filename)
                pw_src = PointSourceHz()
                pw_src.attach(idx, pw_src_param)

        return pw_src


class TotalFieldScatteredField(Src):
    """Set a total and scattered field zone to launch a plane wave.
    
    """
    def __init__(self, src_time, center, size, direction, polarization, amp=1):
        """Constructor
        
        Arguments:
        center -- center of the incidence interface. The beam axis crosses
                  this point.
           type: a tuple with three real numbers.
        size --  size of the incidence interface plane.
           type: a tuple with three real numbers.
        direction -- propagation direction of the beam.
           type: a tuple with three real numbers.
        freq -- oscillating frequency of the beam.
           type: a real number
        polarization -- electric field direction of the beam. 
           type: a tuple with three real numbers.
        amp -- amplitude of the plane wave. The default is 1.
           type: a real number

        """
        if isinstance(src_time, SrcTime):
            self.src_time = src_time
        else:
            raise TypeError, 'src_time must be an instance of SrcTime.'
        
        self.k = np.array(direction, np.double) / norm(direction)
        self.center = np.array(center, np.double)
        self.size = np.array(size, np.double)
        
        self.half_size = .5 * self.size
        self.e_direction = np.array(polarization, np.double) / norm(polarization)
        
        # direction of h field
        self.h_direction = cross(self.k, self.e_direction)
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
        
        self.on_axis_k = self._axis_in_k()
        
    def init(self, geom_tree, space, cmplx):
        self.geom_tree = geom_tree
        self.src_time.init(cmplx)
        
        self.aux_fdtd = self._get_aux_fdtd(space, geom_tree, cmplx)
        self.aux_fdtd.init()

    def step(self):
        self.aux_fdtd.step()
        
    def display_info(self, indent=0):
        print " " * indent, "plane-wave source:"
        print " " * indent, "propagation direction:", self.k
        print " " * indent, "center:", self.center
        print " " * indent, "source plane size:", self.size 
        print " " * indent, "polarization direction:", self.e_direction
        print " " * indent, "amplitude:", self.amp
        
        self.src_time.display_info(4)
        
    def mode_function(self, x, y, z):
        return 1.0
        
    def _dist_from_center(self, point):
        """Calculate distance from the interface plane center.
        
        Arguments:
        point -- location in the space coordinate
            
        """
        return norm(self.center - point)
    
    def _metric_from_center_along_beam_axis(self, point):
        """Calculate projected distance from center along the beam axis.

        Returns positive value when the point is located in
        the k direction to the center.
        
        Keyword arguments:
        point -- location in the space coordinate
            
        """
        return dot(self.k, point - self.center)

    def _dist_from_beam_axis(self, x, y, z):
        """Calculate distance from the beam axis.
        
        Keyword arguments:
        point -- location in the space coordinate 
            
        """   
        return norm(cross(self.k, (x, y, z) - self.center))
    
    def _axis_in_k(self):
        """Return the biggest component direction of k.
        
        """
        dot_with_axis = {}
        
        dot_with_axis[dot(const.PlusX.vector, self.k)] = const.PlusX 
        dot_with_axis[dot(const.PlusY.vector, self.k)] = const.PlusY 
        dot_with_axis[dot(const.PlusZ.vector, self.k)] = const.PlusZ 
        dot_with_axis[dot(const.MinusX.vector, self.k)] = const.MinusX 
        dot_with_axis[dot(const.MinusY.vector, self.k)] = const.MinusY 
        dot_with_axis[dot(const.MinusZ.vector, self.k)] = const.MinusZ 
        
        return dot_with_axis[max(dot_with_axis)]

    def _get_wave_number(self, k, eps_inf, mu_inf, space):
        """Calculate the wave number for auxiliary fdtd using Newton's method.
        
        Keyword arguments:
        k -- normalized wave vector
        eps_inf -- permittivity which fills the auxiliary fdtd
        mu_inf -- permeability which fills the auxiliary fdtd
        space -- Cartesian instance
            
        """
        ds = np.array(space.dr)
        dt = space.dt
        v = 1 / sqrt(eps_inf * mu_inf)
        omega = 2 * pi * self.src_time.freq
        wave_number = omega / v
        wave_vector = wave_number * np.array(k)
        
        zeta = bisect(self._3d_dispersion_relation, 0, 2,
                      (v, omega, ds, dt, wave_vector))

        return zeta * wave_number

    def _3d_dispersion_relation(self, zeta, v, omega, ds, dt, k):
        """
        Keyword arguments:
        zeta: a scalar factor which is yet to be determined.
        v: the phase speed of the wave in the default medium.
        omega: the angular frequency of the input wave.
        ds: the space-cell size, (dx, dy, dz)
        dt: the time step
        k: the true wave vector, (kx, ky, kz)

        Equation 5.65 at p.214 of 'A. Taflove and S. C. Hagness, Computational
        Electrodynamics: The Finite-Difference Time-Domain Method, Third 
        Edition, 3rd ed. Artech House Publishers, 2005'.

        """
        lhs = (sin(0.5 * dt * omega) / v / dt)**2
        rhs = sum((np.sin(0.5 * zeta * np.array(ds) * k) / ds)**2)

        return lhs - rhs

    def _1d_dispersion_relation(self, ds, zeta, v, omega, dt, k):
        """
        Keyword arguments:
        ds: an 1D cell-size which is yet to be determined
        zeta: the scalar factor which relates the true and numerical 
              wavenumber
        v: the phase speed of the input wave in the default medium
        omega: the angular frequency of the input wave
        dt: the time step
        k: the true wavenumber

        Equation 5.67 at p.215 of A. Taflove and S. C. Hagness, Computational
        Electrodynamics: The Finite-Difference Time-Domain Method, Third 
        Edition, 3rd ed. Artech House Publishers, 2005.

        """
        lhs = sin(0.5 * omega * dt) / v / dt
        if ds == 0:
            rhs = 0.5 * k * zeta
        else:
            rhs = sin(0.5 * k * zeta * ds) / ds
        return lhs - rhs

    def _get_aux_fdtd(self, space, geom_tree, cmplx):
        """Returns a TEMz FDTD for a reference of a plane wave.
        
        The space-cell size of the aux_fdtd is calculated using the matched
        numerical dispersion technique. This method assumes that dx=dy=dz.

        """
        default_medium = geom_tree.object_of_point((inf, inf, inf))[0]
        eps_inf = default_medium.material.eps_inf
        mu_inf = default_medium.material.mu_inf
        v = 1 / sqrt(eps_inf * mu_inf)
        
        ds = tuple(space.dr)
        dt = space.dt
        omega = 2 * pi * self.src_time.freq
        wave_vector = omega * self.k / v
        zeta = bisect(self._3d_dispersion_relation, 0, 2,
                      (v, omega, ds, dt, wave_vector))
        wave_number = omega / v
        delta_1d = bisect(self._1d_dispersion_relation, 0, 2 * max(ds),
                          (zeta, v, omega, dt, wave_number))
        
        pml_thickness = 50 * delta_1d
        
        # Find the furthest distance, max_dist from the longitudinal
        # axis of the incomming wave
        #
        # FIXME: When self.size contains numpy.inf vertices contains nan
        #        and the following algorithm does not work.
        vertices = []
        for x in (0.5 * self.size[0], -0.5 * self.size[0]):
            for y in (0.5 * self.size[1], -0.5 * self.size[1]):
                for z in (0.5 * self.size[2], -0.5 * self.size[2]):
                    vertices.append(self.center + (x, y, z))

        dist = map(abs, map(self._metric_from_center_along_beam_axis, vertices))
        max_dist = max(dist)

        longitudinal_size = 2 * (max_dist + pml_thickness + 2 * delta_1d)
        aux_size = (0, 0, longitudinal_size)
        
        mat_objs =  self.geom_tree.material_of_point((inf, inf, inf))[0]
        
        aux_space = Cartesian(size=aux_size,
                              resolution=1/delta_1d,
                              parallel=False)
        aux_geom_list = (DefaultMedium(material=mat_objs),
                         Shell(material=Cpml(),
                               thickness=pml_thickness))
        src_pnt = (0, 0, -max_dist - delta_1d)
        aux_src_list = (PointSource(src_time=deepcopy(self.src_time),
                                    component=const.Ex,
                                    center=src_pnt),)
        
        if cmplx:
            aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, aux_src_list,
                                dt=space.dt, bloch=(0,0,0), verbose=False)
        else:
            aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, aux_src_list,
                                dt=space.dt, bloch=None, verbose=False)
        
        return aux_fdtd

    def _get_pw_source(self, space, component, cosine, field,
                       low_idx, high_idx, source, samp_i2s, face):
        """
        Keyword arguments:
        space - the Coordinate object given as a FDTD argument
        component - Specify the field component
        cosine - the cosine of the field vector and the given component
        field - numpy array of the electromagnetic field
        low_idx - the low end index of the source boundary
        high_idx - the high end index of the source boundary
        source - the pointwise source class
        samp_i2s - the corresponding index_to_space function
        face - which side of the interface
            
        """
        aux_ds = {const.PlusX: space.dr[0], const.MinusX: space.dr[0],  
                  const.PlusY: space.dr[1], const.MinusY: space.dr[1],
                  const.PlusZ: space.dr[2], const.MinusZ: space.dr[2]}
        
        idx_to_spc = {const.Ex: space.ex_index_to_space,
                      const.Ey: space.ey_index_to_space,
                      const.Ez: space.ez_index_to_space,
                      const.Hx: space.hx_index_to_space,
                      const.Hy: space.hy_index_to_space,
                      const.Hz: space.hz_index_to_space}
        
        low_idx_array = np.array(low_idx)
        high_idx_array = np.array(high_idx)

        pw_src = source()
        for i, j, k in ndindex(tuple(high_idx_array - low_idx_array)):
            idx = tuple((i, j, k) + low_idx_array)
            if in_range(idx, field.shape, component):
                pnt = idx_to_spc[component](*idx)
                
                mat_obj, underneath = self.geom_tree.material_of_point(pnt)
                
                if underneath is None:
                    eps_inf = mat_obj.eps_inf
                    mu_inf = mat_obj.mu_inf
                else:
                    eps_inf = underneath.eps_inf
                    mu_inf = underneath.mu_inf

                amp = cosine * self.amp * self.mode_function(*pnt)

                samp_pnt = \
                (0, 0, self._metric_from_center_along_beam_axis(samp_i2s(*idx)))
                
                if isinstance(pw_src, (TransparentEx, TransparentEy, 
                                       TransparentEz)):
                    pw_src_param = TransparentElectricParam(eps_inf, amp, 
                                                            self.aux_fdtd, 
                                                            samp_pnt, face)
                    pw_src.attach(idx, pw_src_param)
                if isinstance(pw_src, (TransparentHx, TransparentHy, 
                                       TransparentHz)):
                    pw_src_param = TransparentMagneticParam(mu_inf, amp, 
                                                            self.aux_fdtd, 
                                                            samp_pnt, face)
                    pw_src.attach(idx, pw_src_param)
                
        return pw_src

    def get_pw_source_ex(self, ex_field, space, geom_tree):
        pw_src = TransparentEx()

        cosine = dot(self.h_direction, (0, 0, 1))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_ex_minus_y(ex_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)

            pw_src_tmp = self._get_pw_source_ex_plus_y(ex_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            
        cosine = dot(self.h_direction, (0, 1, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_ex_minus_z(ex_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_ex_plus_z(ex_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)

        if len(pw_src._param) == 0:
            return None
        else:
            return pw_src

    def _get_pw_source_ex_minus_y(self, ex_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ex_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j, k)
            
            pw_src = self._get_pw_source(space, const.Ex, cosine, ex_field, 
                                         low_idx, high_idx, TransparentEx, 
                                         hz_i2s, const.MinusY)
        else:
            pw_src = None
            
        return pw_src
        
    def _get_pw_source_ex_plus_y(self, ex_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
            pw_src = self._get_pw_source(space, const.Ex, cosine, ex_field,
                                         low_idx, high_idx, TransparentEx, 
                                         hz_i2s, const.PlusY)
        
        else:
            pw_src = None

        return pw_src
        
    def _get_pw_source_ex_minus_z(self, ex_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
            
            hy_i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k)
            
            pw_src = self._get_pw_source(space, const.Ex, cosine, ex_field,
                                         low_idx, high_idx, TransparentEx, 
                                         hy_i2s, const.MinusZ)
        else:
            pw_src = None

        return pw_src
    
    def _get_pw_source_ex_plus_z(self, ex_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
            
            i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
            pw_src = self._get_pw_source(space, const.Ex, cosine, ex_field,
                                         low_idx, high_idx, TransparentEx, 
                                         i2s, const.PlusZ)
        else:
            pw_src = None

        return pw_src
        
    def get_pw_source_ey(self, ey_field, space, geom_tree):
        pw_src = TransparentEy()

        cosine = dot(self.h_direction, (1, 0, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_ey_minus_z(ey_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_ey_plus_z(ey_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            
        cosine = dot(self.h_direction, (0, 0, 1))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_ey_minus_x(ey_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_ey_plus_x(ey_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
        
        if len(pw_src._param) == 0:
            return None
        else:
            return pw_src

    def _get_pw_source_ey_minus_z(self, ey_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k)
            
            pw_src = self._get_pw_source(space, const.Ey, cosine, ey_field,
                                         low_idx, high_idx, TransparentEy, 
                                         hx_i2s, const.MinusZ)
        else:
            pw_src = None
            
        return pw_src

    def _get_pw_source_ey_plus_z(self, ey_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)
            
            pw_src = self._get_pw_source(space, const.Ey, cosine, ey_field,
                                         low_idx, high_idx, TransparentEy,
                                         hx_i2s, const.PlusZ)
        else:
            pw_src = None

        return pw_src

    def _get_pw_source_ey_minus_x(self, ey_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i, j + 1, k)
            
            pw_src = self._get_pw_source(space, const.Ey, cosine, ey_field, 
                                     low_idx, high_idx, TransparentEy,
                                     hz_i2s, const.MinusX)
        else:
            pw_src = None
            
        return pw_src

    def _get_pw_source_ey_plus_x(self, ey_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
            pw_src = self._get_pw_source(space, const.Ey, cosine, ey_field, 
                                         low_idx, high_idx, TransparentEy,
                                         hz_i2s, const.PlusX)
        else:
            pw_src = None
            
        return pw_src

    def get_pw_source_ez(self, ez_field, space, geom_tree):
        pw_src = TransparentEz()

        cosine = dot(self.h_direction, (0, 1, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_ez_minus_x(ez_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_ez_plus_x(ez_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            
        cosine = dot(self.h_direction, (1, 0, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_ez_minus_y(ez_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_ez_plus_y(ez_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
        
        if len(pw_src._param) == 0:
            return None
        else:
            return pw_src

    def _get_pw_source_ez_minus_x(self, ez_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ez_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hy_i2s = lambda i, j, k: space.hy_index_to_space(i, j, k + 1)
            
            pw_src = self._get_pw_source(space, const.Ez, cosine, ez_field, 
                                         low_idx, high_idx, TransparentEz, 
                                         hy_i2s, const.MinusX)
        else: 
            pw_src = None

        return pw_src

    def _get_pw_source_ez_plus_x(self, ez_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hy_i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
            pw_src = self._get_pw_source(space, const.Ez, cosine, ez_field,
                                         low_idx, high_idx, TransparentEz,
                                         hy_i2s, const.PlusX)
        else: 
            pw_src = None

        return pw_src

    def _get_pw_source_ez_minus_y(self, ez_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ez_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j, k + 1)

            pw_src = self._get_pw_source(space, const.Ez, cosine, ez_field,
                                         low_idx, high_idx, TransparentEz, 
                                         hx_i2s, const.MinusY)
        else: 
            pw_src = None

        return pw_src

    def _get_pw_source_ez_plus_y(self, ez_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)

            pw_src = self._get_pw_source(space, const.Ez, cosine, ez_field,  
                                         low_idx, high_idx, TransparentEz, 
                                         hx_i2s, const.PlusY)
        else: 
            pw_src = None

        return pw_src

    def get_pw_source_hx(self, hx_field, space, geom_tree):
        pw_src = TransparentHx()

        cosine = dot(self.e_direction, (0, 0, 1))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_hx_minus_y(hx_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_hx_plus_y(hx_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            
        cosine = dot(self.e_direction, (0, 1, 0))
        if cosine != 0:   
            pw_src_tmp = self._get_pw_source_hx_minus_z(hx_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_hx_plus_z(hx_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
        
        if len(pw_src._param) == 0:
            return None
        else:
            return pw_src

    def _get_pw_source_hx_minus_y(self, hx_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ez_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
            
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)    
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i, j, k - 1)
            
            pw_src = self._get_pw_source(space, const.Hx, cosine, hx_field,
                                         low_idx, high_idx, TransparentHx, 
                                         ez_i2s, const.MinusY)
        else: 
            pw_src = None

        return pw_src

    def _get_pw_source_hx_plus_y(self, hx_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
            
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)    
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i, j - 1, k - 1)
            
            pw_src = self._get_pw_source(space, const.Hx, cosine, hx_field, 
                                         low_idx, high_idx, TransparentHx,
                                         ez_i2s, const.PlusY)
    
        else: 
            pw_src = None

        return pw_src

    def _get_pw_source_hx_minus_z(self, hx_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ey_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
            
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i, j - 1, k - 1)
            
            pw_src = self._get_pw_source(space, const.Hx, cosine, hx_field,
                                         low_idx, high_idx, TransparentHx,
                                         ey_i2s, const.MinusZ)
        else: 
            pw_src = None

        return pw_src

    def _get_pw_source_hx_plus_z(self, hx_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            ey_low_idx = space.space_to_ey_index(*low)
            ey_high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
            
            low_idx = (ey_low_idx[0], ey_low_idx[1] + 1, ey_low_idx[2] + 1)
            high_idx = (ey_high_idx[0], ey_high_idx[1] + 1, ey_high_idx[2] + 1)
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i, j - 1, k)
            
            pw_src = self._get_pw_source(space, const.Hx, cosine, hx_field,
                                         low_idx, high_idx, TransparentHx,
                                         ey_i2s, const.PlusZ)
        else: 
            pw_src = None

        return pw_src

    def get_pw_source_hy(self, hy_field, space, geom_tree):
        pw_src = TransparentHy()

        cosine = dot(self.e_direction, (1, 0, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_hy_minus_z(hy_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_hy_plus_z(hy_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)

        cosine = dot(self.e_direction, (0, 0, 1))
        if cosine != 0:   
            pw_src_tmp = self._get_pw_source_hy_minus_x(hy_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_hy_plus_x(hy_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)

        if len(pw_src._param) == 0:
            return None
        else:
            return pw_src

    def _get_pw_source_hy_minus_z(self, hy_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            
            ex_i2s = lambda i, j, k: space.ex_index_to_space(i - 1, j, k)
            
            pw_src = self._get_pw_source(space, const.Hy, cosine, hy_field,
                                         low_idx, high_idx, TransparentHy,
                                         ex_i2s, const.MinusZ)
        else:
            pw_src = None
        
        return pw_src

    def _get_pw_source_hy_plus_z(self, hy_field, space, cosine):
        if 2 * space.half_size[2] > space.dr[2]:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            
            ex_i2s = lambda i, j, k: space.ex_index_to_space(i - 1, j, k - 1)
            
            pw_src = self._get_pw_source(space, const.Hy, cosine, hy_field,
                                         low_idx, high_idx, TransparentHy,
                                         ex_i2s, const.PlusZ)
        else:
            pw_src = None
        
        return pw_src

    def _get_pw_source_hy_minus_x(self, hy_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            ez_low_idx = space.space_to_ez_index(*low)
            ez_high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
                
            low_idx = (ez_low_idx[0], ez_low_idx[1], ez_low_idx[2] + 1)
            high_idx = (ez_high_idx[0], ez_high_idx[1], ez_high_idx[2] + 1)
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i, j, k - 1)
            
            pw_src = self._get_pw_source(space, const.Hy, cosine, hy_field,
                                         low_idx, high_idx, TransparentHy,
                                         ez_i2s, const.MinusX)
        else:
            pw_src = None
        
        return pw_src

    def _get_pw_source_hy_plus_x(self, hy_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            ez_low_idx = space.space_to_ez_index(*low)
            ez_high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
                
            low_idx = (ez_low_idx[0] + 1, ez_low_idx[1], ez_low_idx[2] + 1)
            high_idx = (ez_high_idx[0] + 1, ez_high_idx[1], ez_high_idx[2] + 1)
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i - 1, j, k - 1)
            
            pw_src = self._get_pw_source(space, const.Hy, cosine, hy_field,
                                         low_idx, high_idx, TransparentHy,
                                         ez_i2s, const.PlusX)
        else:
            pw_src = None
        
        return pw_src

    def get_pw_source_hz(self, hz_field, space, geom_tree):
        pw_src = TransparentHz()

        cosine = dot(self.e_direction, (0, 1, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_hz_minus_x(hz_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_hz_plus_x(hz_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            
        cosine = dot(self.e_direction, (1, 0, 0))
        if cosine != 0:
            pw_src_tmp = self._get_pw_source_hz_minus_y(hz_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
            pw_src_tmp = self._get_pw_source_hz_plus_y(hz_field, space, cosine)
            if pw_src_tmp != None: 
                pw_src.merge(pw_src_tmp)
        
        if len(pw_src._param) == 0:
            return None
        else:
            return pw_src

    def _get_pw_source_hz_minus_x(self, hz_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ey_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
                
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i, j - 1, k)
            
            pw_src = self._get_pw_source(space, const.Hz, cosine, hz_field,
                                         low_idx, high_idx, TransparentHz,
                                         ey_i2s, const.MinusX)
        else:
            pw_src = None
        
        return pw_src

    def _get_pw_source_hz_plus_x(self, hz_field, space, cosine):
        if 2 * space.half_size[0] > space.dr[0]:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i - 1, j - 1, k)
            
            pw_src = self._get_pw_source(space, const.Hz, cosine, hz_field,
                                         low_idx, high_idx, TransparentHz,
                                         ey_i2s, const.PlusX)
        else:
            pw_src = None
        
        return pw_src

    def _get_pw_source_hz_minus_y(self, hz_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            
            ex_i2s = lambda i, j, k: space.ex_index_to_space(i, j - 1, k)
            
            pw_src = self._get_pw_source(space, const.Hz, cosine, hz_field,
                                         low_idx, high_idx, TransparentHz,
                                         ex_i2s, const.MinusY)
        else:
            pw_src = None
        
        return pw_src

    def _get_pw_source_hz_plus_y(self, hz_field, space, cosine):
        if 2 * space.half_size[1] > space.dr[1]:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
            ex_2s = lambda i, j, k: space.ex_index_to_space(i - 1, j - 1, k)
            
            pw_src = self._get_pw_source(space, const.Hz, cosine, hz_field,
                                         low_idx, high_idx, TransparentHz,
                                         ex_i2s, const.PlusY)
        else:
            pw_src = None
        
        return pw_src
        

class GaussianBeam(TotalFieldScatteredField):
    """Launch a transparent Gaussian beam.
    
    It works as a guided mode with Gaussian profile is launched through the 
    incidence interface. The incidence interface is transparent, thus the 
    scattered wave can penetrate through the interface plane.
    
    """
    def __init__(self, src_time, directivity, center, size, direction, 
                 polarization, waist=inf, amp=1):
        """
        
        Keyword arguments:
        directivity -- directivity of the incidence interface.
           type: a child class of constant.Directional.
        center -- center of the incidence interface. The beam axis crosses
                  this point.
           type: a tuple with three real numbers.
        size --  size of the incidence interface plane.
           type: a tuple with three real numbers.
        direction -- propagation direction of the beam.
           type: a tuple with three real numbers.
        polarization -- electric field direction of the beam.
           type: a tuple with three real numbers.
        waist -- the Gaussian beam radius. The default is inf.
           type: a tuple with three real numbers.
        amp -- amplitude of the plane wave. The default is 1.
           type: a tuple with three real numbers.

        """
        TotalFieldScatteredField.__init__(self, src_time, center, size, 
                                          direction, polarization, amp)
        
        if issubclass(directivity, const.Directional):
            self.directivity = directivity
        else:
            raise TypeError, 'directivity must be a Directional type.'
        
        # spot size of Gaussian beam
        self.waist = float(waist)

    def init(self, geom_tree, space, cmplx):
        self.geom_tree = geom_tree
        self.src_time.init(cmplx)
        
        aux_fdtd = self._get_aux_fdtd(space, geom_tree, cmplx)
        aux_fdtd.init()
        raising = aux_fdtd.src_list[0].src_time.width
        dist = 2 * aux_fdtd.space.half_size[2]
        default_medium = (i for i in aux_fdtd.geom_list 
                            if isinstance(i, DefaultMedium)).next()
        eps_inf = default_medium.material.eps_inf
        mu_inf = default_medium.material.mu_inf
        v_p = 1 / sqrt(eps_inf * mu_inf)
        passby = raising + dist / v_p

        aux_fdtd.step_until_t(2 * passby)
        
        self.aux_fdtd = _GaussianBeamSrcTime(aux_fdtd)

    def display_info(self, indent=0):
        print ' ' * indent, 'Gaussian beam source:'
        print ' ' * indent, 
        print 'propagation direction:', self.k,
        print 'center:', self.center
        print 'source plane size:', self.size 
        print 'polarization direction:', self.e_direction
        print 'beam waist:', self.waist
        print 'maximum amp.:', self.amp
        
        self.src_time.display_info(indent + 4)
    
    def mode_function(self, x, y, z):
        r = self._dist_from_beam_axis(x, y, z)
        return exp(-(r / self.waist)**2)
        
    def get_pw_source_ex(self, ex_field, space, geom_tree):
        if self.directivity is const.PlusY:
            cosine = dot(self.h_direction, (0, 0, 1))
            return self._get_pw_source_ex_minus_y(ex_field, space, cosine)
            
        elif self.directivity is const.MinusY:
            cosine = dot(self.h_direction, (0, 0, 1))
            return self._get_pw_source_ex_plus_y(ex_field, space, cosine)

        elif self.directivity is const.PlusZ:
            cosine = dot(self.h_direction, (0, 1, 0))
            return self._get_pw_source_ex_minus_z(ex_field, space, cosine)

        elif self.directivity is const.MinusZ:
            cosine = dot(self.h_direction, (0, 1, 0))
            return self._get_pw_source_ex_plus_z(ex_field, space, cosine)

        else:
            return None

    def get_pw_source_ey(self, ey_field, space, geom_tree):
        if self.directivity is const.PlusZ:
            cosine = dot(self.h_direction, (1, 0, 0))
            return self._get_pw_source_ey_minus_z(ey_field, space, cosine)
            
        elif self.directivity is const.MinusZ:
            cosine = dot(self.h_direction, (1, 0, 0))
            return self._get_pw_source_ey_plus_z(ey_field, space, cosine)

        elif self.directivity is const.PlusX:
            cosine = dot(self.h_direction, (0, 0, 1))
            return self._get_pw_source_ey_minus_x(ey_field, space, cosine)

        elif self.directivity is const.MinusX:
            cosine = dot(self.h_direction, (0, 0, 1))
            return self._get_pw_source_ey_plus_x(ey_field, space, cosine)

        else:
            return None

    def get_pw_source_ez(self, ez_field, space, geom_tree):
        if self.directivity is const.PlusX:
            cosine = dot(self.h_direction, (0, 1, 0))
            return self._get_pw_source_ez_minus_x(ez_field, space, cosine)
            
        elif self.directivity is const.MinusX:
            cosine = dot(self.h_direction, (0, 1, 0))
            return self._get_pw_source_ez_plus_x(ez_field, space, cosine)

        elif self.directivity is const.PlusY:
            cosine = dot(self.h_direction, (1, 0, 0))
            return self._get_pw_source_ez_minus_y(ez_field, space, cosine)

        elif self.directivity is const.MinusY:
            cosine = dot(self.h_direction, (1, 0, 0))
            return self._get_pw_source_ez_plus_y(ez_field, space, cosine)

        else:
            return None
        
    def get_pw_source_hx(self, hx_field, space, geom_tree):
        if self.directivity is const.PlusY:
            cosine = dot(self.e_direction, (0, 0, 1))
            return self._get_pw_source_hx_minus_y(hx_field, space, cosine)

        elif self.directivity is const.MinusY:
            cosine = dot(self.e_direction, (0, 0, 1))
            return self._get_pw_source_hx_plus_y(hx_field, space, cosine)

        elif self.directivity is const.PlusZ:
            cosine = dot(self.e_direction, (0, 1, 0))
            return self._get_pw_source_hx_minus_z(hx_field, space, cosine)

        elif self.directivity is const.MinusZ:
            cosine = dot(self.e_direction, (0, 1, 0))
            return self._get_pw_source_hx_plus_z(hx_field, space, cosine)

        else:
            return None
        
    def get_pw_source_hy(self, hy_field, space, geom_tree):
        if self.directivity is const.PlusZ:
            cosine = dot(self.e_direction, (1, 0, 0))
            return self._get_pw_source_hy_minus_z(hy_field, space, cosine)

        elif self.directivity is const.MinusZ:
            cosine = dot(self.e_direction, (1, 0, 0))
            return self._get_pw_source_hy_plus_z(hy_field, space, cosine)

        elif self.directivity is const.PlusX:
            cosine = dot(self.e_direction, (0, 0, 1))
            return self._get_pw_source_hy_minus_x(hy_field, space, cosine)

        elif self.directivity is const.MinusX:
            cosine = dot(self.e_direction, (0, 0, 1))
            return self._get_pw_source_hy_plus_x(hy_field, space, cosine)

        else:
            return None
        
    def get_pw_source_hz(self, hz_field, space, geom_tree):
        if self.directivity is const.PlusX:
            cosine = dot(self.e_direction, (0, 1, 0))
            return self._get_pw_source_hz_minus_x(hz_field, space, cosine)

        elif self.directivity is const.MinusX:
            cosine = dot(self.e_direction, (0, 1, 0))
            return self._get_pw_source_hz_plus_x(hz_field, space, cosine)

        elif self.directivity is const.PlusY:
            cosine = dot(self.e_direction, (1, 0, 0))
            return self._get_pw_source_hz_minus_y(hz_field, space, cosine)

        elif self.directivity is const.MinusY:
            cosine = dot(self.e_direction, (1, 0, 0))
            return self._get_pw_source_hz_plus_y(hz_field, space, cosine)

        else:
            return None


class _GaussianBeamSrcTime(object):
    class EX(object):
        def __init__(self, outer):
            self.outer = outer

        def __getitem__(self, idx):
            return self.outer.envelope() * self.outer.aux_fdtd.ex[idx]

    class HY(object):
        def __init__(self, outer):
            self.outer = outer

        def __getitem__(self, idx):
            return self.outer.envelope() * self.outer.aux_fdtd.hy[idx]

    def __init__(self, aux_fdtd):
        self.aux_fdtd = aux_fdtd
        self.space = self.aux_fdtd.space
        self.ex = self.EX(self)
        self.hy = self.HY(self)
        
        self.n = 0
        self.t = 0

    def step(self):
        self.aux_fdtd.step()
        self.n += 1
        self.t = self.n * self.aux_fdtd.time_step.dt
        
    def envelope(self):
        width = self.aux_fdtd.src_list[0].src_time.width
        if self.t < width:
            env = sin(0.5 * pi * self.t / width)**2
        else:
            env = 1
        return env
