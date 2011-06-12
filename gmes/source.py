#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')
    
from copy import deepcopy
from math import sqrt, pi, sin, cos, exp
import cmath as cm

import numpy as np
from numpy import inf, cross, dot, ndindex
from numpy.linalg import norm
from scipy.optimize import bisect

import constants as const
from geometry import Cartesian, DefaultMedium, Boundary, in_range
from fdtd import TEMzFDTD
from material import Dielectric, CPML


#
# SrcTime: Continuous, Bandpass
# Src: Dipole, GaussianBeam, TotalFieldScatteredField
#

class SrcTime(object):
    """Time-dependent part of a source.
    
    """
    def display_info(self, indent=0):
        pass
    

class Src(object):
    """Space-dependent part of a source.
    
    """
    def display_info(self, indent=0):
        pass
    
        
class Continuous(SrcTime):
    """Continuous (CW) source with (optional) slow turn-on and/or turn-off.
    
    """
    def __init__(self, freq, phase=0, start=0, end=inf, width=None):
        self.freq = float(freq)
        self.phase = float(phase)
        self.start = float(start)
        self.end = float(end)
        
        if width is None:
            self.width = 3 / self.freq
        else:
            self.width = float(width)

    def init(self, cmplx):
        self.cmplx = cmplx
        
    def display_info(self, indent=0):
        print " " * indent, "continuous source"
        print " " * indent,
        print "frequency:", self.freq,
        print "initial phase advance:", self.phase,
        print "start time:", self.start,
        print "end time:", self.end,
        print "raising duration:", self.width
                
    def dipole(self, time):
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
        
        osc = env * cm.exp(-2j * pi * self.freq * time - self.phase)
        if self.cmplx:
            return osc
        else:
            return osc.real


class Bandpass(SrcTime):
    """a pulse source with Gaussian-envelope
    
    """
    def __init__(self, freq, fwidth, s=5, phase=0):
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
        
    def dipole(self, time):
        tt = time - self.peak_time
        if (abs(tt) > self.cutoff): 
            return 0

        # correction factor so that current amplitude (= d(dipole)/dt) is
        # ~ 1 near the peak of the Gaussian.
        cfactor = 1.0 / (-2j * pi * self.freq)
        
        osc = cfactor * exp(-0.5 * (tt / self.width)**2) \
            * cm.exp(-2j * pi * self.freq * time - self.phase)
        if self.cmplx:
            return osc
        else:
            return osc.real
        

from pw_source import DipoleEx, DipoleEy, DipoleEz
from pw_source import DipoleHx, DipoleHy, DipoleHz

        
class Dipole(Src):
    def __init__(self, src_time, pos, component, amp=1, filename=None):
        self.pos = np.array(pos, float)
        self.comp = component
        self.src_time = src_time
        self.amp = float(amp)
        self.filename = filename
        
    def init(self, geom_tree, space, cmplx):
        self.src_time.init(cmplx)
    
    def step(self):
        pass
    
    def display_info(self, indent=0):
        print " " * indent, "Hertzian dipole source:"
        print " " * indent, "center:", self.pos
        print " " * indent, "polarization direction:", self.comp
        print " " * indent, "maximum amp.:", self.amp
        print " " * indent, "source recording: ", self.filename
        
        self.src_time.display_info(4)
        
    def set_pointwise_source_ex(self, material_ex, space):
        if self.comp is const.Ex:
            idx = space.space_to_ex_index(*self.pos)
            if in_range(idx, material_ex, const.Ex):
                material_ex[idx] = DipoleEx(material_ex[idx], self.src_time, 
                                            self.amp, self.filename)
                if self.filename is not None:
                    loc = space.ex_index_to_space(*idx)
                    material_ex[idx].file.write('# location=' + str(loc) + '\n')
                    
    def set_pointwise_source_ey(self, material_ey, space):
        if self.comp is const.Ey:
            idx = space.space_to_ey_index(*self.pos)
            if in_range(idx, material_ey, const.Ey):
                material_ey[idx] = DipoleEy(material_ey[idx], self.src_time, 
                                            self.amp, self.filename)
                if self.filename is not None:
                    loc = space.ey_index_to_space(*idx)
                    material_ey[idx].file.write('# location=' + str(loc) + '\n')
                    
    def set_pointwise_source_ez(self, material_ez, space):
        if self.comp is const.Ez:
            idx = space.space_to_ez_index(*self.pos)
            if in_range(idx, material_ez, const.Ez):
                material_ez[idx] = DipoleEz(material_ez[idx], self.src_time, 
                                            self.amp, self.filename)
                if self.filename is not None:
                    loc = space.ez_index_to_space(*idx)
                    material_ez[idx].file.write('# location=' + str(loc) + '\n')
                    
    def set_pointwise_source_hx(self, material_hx, space):
        if self.comp is const.Hx:
            idx = space.space_to_hx_index(*self.pos)
            if in_range(idx, material_hx, const.Hx):
                material_hx[idx] = DipoleHx(material_hx[idx], self.src_time, 
                                            self.amp, self.filename)
                if self.filename is not None:
                    loc = space.hx_index_to_space(*idx)
                    material_hx[idx].file.write('# location=' + str(loc) + '\n')
                    
    def set_pointwise_source_hy(self, material_hy, space):
        if self.comp is const.Hy:
            idx = space.space_to_hy_index(*self.pos)
            if in_range(idx, material_hy, const.Hy):
                material_hy[idx] = DipoleHy(material_hy[idx], self.src_time, 
                                            self.amp, self.filename)
                if self.filename is not None:
                    loc = space.hy_index_to_space(*idx)
                    material_hy[idx].file.write('# location=' + str(loc) + '\n')
                    
    def set_pointwise_source_hz(self, material_hz, space):
        if self.comp is const.Hz:
            idx = space.space_to_hz_index(*self.pos)
            if in_range(idx, material_hz, const.Hz):
                material_hz[idx] = DipoleHz(material_hz[idx], self.src_time, 
                                            self.amp, self.filename)
                if self.filename is not None:
                    loc = space.hz_index_to_space(*idx)
                    material_hz[idx].file.write('# location=' + str(loc) + '\n')


from pw_source import TransparentElectric, TransparentMagnetic
from pw_source import TransparentEx, TransparentEy, TransparentEz
from pw_source import TransparentHx, TransparentHy, TransparentHz


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
        
        self.k = np.array(direction, float) / norm(direction)
        self.center = np.array(center, float)
        self.size = np.array(size, float)
        
        self.half_size = .5 * self.size
        self.e_direction = np.array(polarization, float) / norm(polarization)
        
        # direction of h field
        self.h_direction = cross(self.k, self.e_direction)
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
        
        self.on_axis_k = self._axis_in_k()
        
    def init(self, geom_tree, space, cmplx):
        self.geom_tree = geom_tree
        self.src_time.init(cmplx)
        
        self.aux_fdtd = self._get_aux_fdtd(space, geom_tree, cmplx)
        
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
        
        Arguments:
            point -- location in the space coordinate
            
        """
        return dot(self.k, point - self.center)

    def _dist_from_beam_axis(self, x, y, z):
        """Calculate distance from the beam axis.
        
        Arguments:
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

    def _get_wave_number(self, k, epsilon, mu, space):
        """Calculate the wave number for auxiliary fdtd using Newton's method.
        
        Arguments:
            k -- normalized wave vector
            epsilon -- permittivity which fills the auxiliary fdtd
            mu -- permeability which fills the auxiliary fdtd
            space -- Cartesian instance
            
        """
        ds = np.array((space.dx, space.dy, space.dz))
        dt = space.dt
        v = 1 / sqrt(epsilon * mu)
        omega = 2 * pi * self.src_time.freq
        wave_number = omega / v
        wave_vector = wave_number * np.array(k)
        
        zeta = bisect(self._3d_dispersion_relation, 0, 2,
                      (v, omega, ds, dt, wave_vector))

        return zeta * wave_number

    def _3d_dispersion_relation(self, zeta, v, omega, ds, dt, k):
        """
        Arguments:
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
        Arguments:
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
        eps = default_medium.material.epsilon
        mu = default_medium.material.mu
        v = 1 / sqrt(eps * mu)
        
        ds = (space.dx, space.dy, space.dz)
        dt = space.dt
        omega = 2 * pi * self.src_time.freq
        wave_vector = omega * self.k / v
        zeta = bisect(self._3d_dispersion_relation, 0, 2,
                      (v, omega, ds, dt, wave_vector))
        wave_number = omega / v
        delta_1d = bisect(self._1d_dispersion_relation, 0, 2 * max(ds),
                          (zeta, v, omega, dt, wave_number))
        
        pml_thickness = 10 * delta_1d
        
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

        longitudinal_size = 2 * (max_dist + pml_thickness + delta_1d)
        aux_size = (0, 0, longitudinal_size)

        mat_objs =  self.geom_tree.material_of_point((inf, inf, inf))[0]
        
        aux_space = Cartesian(size=aux_size,
                              resolution=1/delta_1d,
                              parallel=False)
        aux_geom_list = (DefaultMedium(material=mat_objs),
                         Boundary(material=CPML(kappa_max=2.0,
                                                sigma_max_ratio=2.0),
                                  thickness=pml_thickness,
                                  size=aux_size,
                                  minus_z=False))
        src_pnt = aux_space.ex_index_to_space(0, 0, 0)
        aux_src_list = (Dipole(src_time=deepcopy(self.src_time),
                               component=const.Ex,
                               pos=src_pnt),)
        
        if cmplx:
            aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, aux_src_list,
                                dt=space.dt, bloch=(0,0,0), verbose=False)
        else:
            aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, aux_src_list,
                                dt=space.dt, bloch=None, verbose=False)

        return aux_fdtd

    def _set_pw_source(self, space, component, cosine, material, 
                       low_idx, high_idx, source, samp_i2s, face):
        """
        Arguments:
            space - the Coordinate object given as a FDTD argument.
            component - Specify the field component.
            cosine - the cosine of the field vector and the given component.
            material - pointwise material object array.
            low_idx - the low end index of the source boundary
            high_idx - the high end index of the source boundary
            source - the pointwise source class
            samp_i2s - the corresponding index_to_space function
            face - which side of the interface
            
        """
        aux_ds = {const.PlusX: space.dx, const.MinusX: space.dx,  
                  const.PlusY: space.dy, const.MinusY: space.dy,
                  const.PlusZ: space.dz, const.MinusZ: space.dz}
        
        idx_to_spc = {const.Ex: space.ex_index_to_space,
                      const.Ey: space.ey_index_to_space,
                      const.Ez: space.ez_index_to_space,
                      const.Hx: space.hx_index_to_space,
                      const.Hy: space.hy_index_to_space,
                      const.Hz: space.hz_index_to_space}
        
        low_idx_array = np.array(low_idx)
        high_idx_array = np.array(high_idx)
        
        for i, j, k in ndindex(*(high_idx_array - low_idx_array)):
            idx = tuple((i, j, k) + low_idx_array)
            if in_range(idx, material, component):
                pnt = idx_to_spc[component](*idx)
                
                mat_objs = self.geom_tree.material_of_point(pnt)
                eps = mat_objs[0].epsilon
                mu = mat_objs[0].mu

                amp = cosine * self.amp * self.mode_function(*pnt)

                sample_pnt = \
                (0, 0, self._metric_from_center_along_beam_axis(samp_i2s(*idx)))
                
                if issubclass(source, TransparentElectric):
                    material[idx] = source(material[idx], eps, amp,
                                           self.aux_fdtd, sample_pnt, face)
                if issubclass(source, TransparentMagnetic):
                    material[idx] = source(material[idx], mu, amp,
                                           self.aux_fdtd, sample_pnt, face)

    def set_pointwise_source_ex(self, material_ex, space):
        cosine = dot(self.h_direction, (0, 0, 1))
        if cosine != 0:
            self._set_pw_source_ex_minus_y(material_ex, space, cosine)
            self._set_pw_source_ex_plus_y(material_ex, space, cosine)
            
        cosine = dot(self.h_direction, (0, 1, 0))
        if cosine != 0:
            self._set_pw_source_ex_minus_z(material_ex, space, cosine)
            self._set_pw_source_ex_plus_z(material_ex, space, cosine)
        
    def _set_pw_source_ex_minus_y(self, material_ex, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ex_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j, k)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex, 
                                low_idx, high_idx, TransparentEx, 
                                hz_i2s, const.MinusY)
        
    def _set_pw_source_ex_plus_y(self, material_ex, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex,
                                low_idx, high_idx, TransparentEx, 
                                hz_i2s, const.PlusY)
        
    def _set_pw_source_ex_minus_z(self, material_ex, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
            
            hy_i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex,
                                low_idx, high_idx, TransparentEx, 
                                hy_i2s, const.MinusZ)
    
    def _set_pw_source_ex_plus_z(self, material_ex, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
            
            i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex,
                                low_idx, high_idx, TransparentEx, 
                                i2s, const.PlusZ)
        
    def set_pointwise_source_ey(self, material_ey, space):
        cosine = dot(self.h_direction, (1, 0, 0))
        if cosine != 0:
            self._set_pw_source_ey_minus_z(material_ey, space, cosine)
            self._set_pw_source_ey_plus_z(material_ey, space, cosine)
            
        cosine = dot(self.h_direction, (0, 0, 1))
        if cosine != 0:
            self._set_pw_source_ey_minus_x(material_ey, space, cosine)
            self._set_pw_source_ey_plus_x(material_ey, space, cosine)
        
    def _set_pw_source_ey_minus_z(self, material_ey, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,
                                low_idx, high_idx, TransparentEy, 
                                hx_i2s, const.MinusZ)
    
    def _set_pw_source_ey_plus_z(self, material_ey, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,
                                low_idx, high_idx, TransparentEy,
                                hx_i2s, const.PlusZ)
    
    def _set_pw_source_ey_minus_x(self, material_ey, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i, j + 1, k)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,  
                                low_idx, high_idx, TransparentEy,
                                hz_i2s, const.MinusX)
    
    def _set_pw_source_ey_plus_x(self, material_ey, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,  
                                low_idx, high_idx, TransparentEy,
                                hz_i2s, const.PlusX)
        
    def set_pointwise_source_ez(self, material_ez, space):
        cosine = dot(self.h_direction, (0, 1, 0))
        if cosine != 0:
            self._set_pw_source_ez_minus_x(material_ez, space, cosine)
            self._set_pw_source_ez_plus_x(material_ez, space, cosine)
            
        cosine = dot(self.h_direction, (1, 0, 0))
        if cosine != 0:
            self._set_pw_source_ez_minus_y(material_ez, space, cosine)
            self._set_pw_source_ez_plus_y(material_ez, space, cosine)
        
    def _set_pw_source_ez_minus_x(self, material_ez, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ez_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hy_i2s = lambda i, j, k: space.hy_index_to_space(i, j, k + 1)
            
            self._set_pw_source(space, const.Ez, cosine, material_ez,  
                                low_idx, high_idx, TransparentEz, 
                                hy_i2s, const.MinusX)

    def _set_pw_source_ez_plus_x(self, material_ez, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hy_i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
            self._set_pw_source(space, const.Ez, cosine, material_ez,
                                low_idx, high_idx, TransparentEz,
                                hy_i2s, const.PlusX)
    
    def _set_pw_source_ez_minus_y(self, material_ez, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ez_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j, k + 1)

            self._set_pw_source(space, const.Ez, cosine, material_ez,
                                low_idx, high_idx, TransparentEz, 
                                hx_i2s, const.MinusY)
    
    def _set_pw_source_ez_plus_y(self, material_ez, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(*low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))

            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)

            self._set_pw_source(space, const.Ez, cosine, material_ez,  
                                low_idx, high_idx, TransparentEz, 
                                hx_i2s, const.PlusY)
        
    def set_pointwise_source_hx(self, material_hx, space):
        cosine = dot(self.e_direction, (0, 0, 1))
        if cosine != 0:
            self._set_pw_source_hx_minus_y(material_hx, space, cosine)
            self._set_pw_source_hx_plus_y(material_hx, space, cosine)
            
        cosine = dot(self.e_direction, (0, 1, 0))
        if cosine != 0:   
            self._set_pw_source_hx_minus_z(material_hx, space, cosine)
            self._set_pw_source_hx_plus_z(material_hx, space, cosine)
        
    def _set_pw_source_hx_minus_y(self, material_hx, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ez_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
            
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)    
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i, j, k - 1)
            
            self._set_pw_source(space, const.Hx, cosine, material_hx,
                                low_idx, high_idx, TransparentHx, 
                                ez_i2s, const.MinusY)
    
    def _set_pw_source_hx_plus_y(self, material_hx, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
            
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)    
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i, j - 1, k - 1)
            
            self._set_pw_source(space, const.Hx, cosine, material_hx,  
                                low_idx, high_idx, TransparentHx,
                                ez_i2s, const.PlusY)
    
    def _set_pw_source_hx_minus_z(self, material_hx, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ey_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
            
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i, j - 1, k - 1)
            
            self._set_pw_source(space, const.Hx, cosine, material_hx,
                                low_idx, high_idx, TransparentHx,
                                ey_i2s, const.MinusZ)
    
    def _set_pw_source_hx_plus_z(self, material_hx, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            ey_low_idx = space.space_to_ey_index(*low)
            ey_high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
            
            low_idx = (ey_low_idx[0], ey_low_idx[1] + 1, ey_low_idx[2] + 1)
            high_idx = (ey_high_idx[0], ey_high_idx[1] + 1, ey_high_idx[2] + 1)
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i, j - 1, k)
            
            self._set_pw_source(space, const.Hx, cosine, material_hx,
                                low_idx, high_idx, TransparentHx,
                                ey_i2s, const.PlusZ)
        
    def set_pointwise_source_hy(self, material_hy, space):
        cosine = dot(self.e_direction, (1, 0, 0))
        if cosine != 0:
            self._set_pw_source_hy_minus_z(material_hy, space, cosine)
            self._set_pw_source_hy_plus_z(material_hy, space, cosine)
         
        cosine = dot(self.e_direction, (0, 0, 1))
        if cosine != 0:   
            self._set_pw_source_hy_minus_x(material_hy, space, cosine)
            self._set_pw_source_hy_plus_x(material_hy, space, cosine)
                
    def _set_pw_source_hy_minus_z(self, material_hy, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            
            ex_i2s = lambda i, j, k: space.ex_index_to_space(i - 1, j, k)
            
            self._set_pw_source(space, const.Hy, cosine, material_hy,
                                low_idx, high_idx, TransparentHy,
                                ex_i2s, const.MinusZ)
    
    def _set_pw_source_hy_plus_z(self, material_hy, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            
            ex_i2s = lambda i, j, k: space.ex_index_to_space(i - 1, j, k - 1)
            
            self._set_pw_source(space, const.Hy, cosine, material_hy,
                                low_idx, high_idx, TransparentHy,
                                ex_i2s, const.PlusZ)
    
    def _set_pw_source_hy_minus_x(self, material_hy, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            ez_low_idx = space.space_to_ez_index(*low)
            ez_high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
                
            low_idx = (ez_low_idx[0], ez_low_idx[1], ez_low_idx[2] + 1)
            high_idx = (ez_high_idx[0], ez_high_idx[1], ez_high_idx[2] + 1)
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i, j, k - 1)
            
            self._set_pw_source(space, const.Hy, cosine, material_hy,
                                low_idx, high_idx, TransparentHy,
                                ez_i2s, const.MinusX)
    
    def _set_pw_source_hy_plus_x(self, material_hy, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            ez_low_idx = space.space_to_ez_index(*low)
            ez_high_idx = map(lambda x: x + 1, space.space_to_ez_index(*high))
                
            low_idx = (ez_low_idx[0] + 1, ez_low_idx[1], ez_low_idx[2] + 1)
            high_idx = (ez_high_idx[0] + 1, ez_high_idx[1], ez_high_idx[2] + 1)
            
            ez_i2s = lambda i, j, k: space.ez_index_to_space(i - 1, j, k - 1)
            
            self._set_pw_source(space, const.Hy, cosine, material_hy,
                                low_idx, high_idx, TransparentHy,
                                ez_i2s, const.PlusX)
        
    def set_pointwise_source_hz(self, material_hz, space):
        cosine = dot(self.e_direction, (0, 1, 0))
        if cosine != 0:
            self._set_pw_source_hz_minus_x(material_hz, space, cosine)
            self._set_pw_source_hz_plus_x(material_hz, space, cosine)
            
        cosine = dot(self.e_direction, (1, 0, 0))
        if cosine != 0:
            self._set_pw_source_hz_minus_y(material_hz, space, cosine)
            self._set_pw_source_hz_plus_y(material_hz, space, cosine)
        
    def _set_pw_source_hz_minus_x(self, material_hz, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ey_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
                
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i, j - 1, k)
            
            self._set_pw_source(space, const.Hz, cosine, material_hz,
                                low_idx, high_idx, TransparentHz,
                                ey_i2s, const.MinusX)
    
    def _set_pw_source_hz_plus_x(self, material_hz, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
            ey_i2s = lambda i, j, k: space.ey_index_to_space(i - 1, j - 1, k)
            
            self._set_pw_source(space, const.Hz, cosine, material_hz,
                                low_idx, high_idx, TransparentHz,
                                ey_i2s, const.PlusX)
    
    def _set_pw_source_hz_minus_y(self, material_hz, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            
            ex_i2s = lambda i, j, k: space.ex_index_to_space(i, j - 1, k)
            
            self._set_pw_source(space, const.Hz, cosine, material_hz,
                                low_idx, high_idx, TransparentHz,
                                ex_i2s, const.MinusY)
    
    def _set_pw_source_hz_plus_y(self, material_hz, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(*low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(*high))
                
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
            ex_2s = lambda i, j, k: space.ex_index_to_space(i - 1, j - 1, k)
            
            self._set_pw_source(space, const.Hz, cosine, material_hz,
                                low_idx, high_idx, TransparentHz,
                                ex_i2s, const.PlusY)
        
        
class GaussianBeam(TotalFieldScatteredField):
    """Launch a transparent Gaussian beam.
    
    It works as a guided mode with Gaussian profile is launched through the 
    incidence interface. The incidence interface is transparent, thus the 
    scattered wave can penetrate through the interface plane.
    
    """
    def __init__(self, src_time, directivity, center, size, direction, 
                 polarization, waist=inf, amp=1):
        """
        
        Arguments:
            directivity -- directivity of the incidence interface.
               type: a child class of constants.Directional.
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
        raising = aux_fdtd.src_list[0].src_time.width
        dist = 2 * aux_fdtd.space.half_size[2]
        default_medium = (i for i in aux_fdtd.geom_list 
                            if isinstance(i, DefaultMedium)).next()
        eps = default_medium.material.epsilon
        mu = default_medium.material.mu
        v_p = 1 / sqrt(eps * mu)
        passby = raising + dist / v_p

        while aux_fdtd.time_step.t < 2 * passby:
            aux_fdtd.step()
        
        self.aux_fdtd = _GaussianBeamSrcTime(aux_fdtd)

    def display_info(self, indent=0):
        print " " * indent, "Gaussian beam source:"
        print " " * indent, "propagation direction:", self.k
        print " " * indent, "center:", self.center
        print " " * indent, "source plane size:", self.size 
        print " " * indent, "polarization direction:", self.e_direction
        print " " * indent, "beam waist:", self.waist
        print " " * indent, "maximum amp.:", self.amp
        
        self.src_time.display_info(indent + 4)
    
    def mode_function(self, x, y, z):
        r = self._dist_from_beam_axis(x, y, z)
        return exp(-(r / self.waist)**2)
        
    def set_pointwise_source_ex(self, material_ex, space):
        if self.directivity is const.PlusY:
            cosine = dot(self.h_direction, (0, 0, 1))
            self._set_pw_source_ex_minus_y(material_ex, space, cosine)
            
        elif self.directivity is const.MinusY:
            cosine = dot(self.h_direction, (0, 0, 1))
            self._set_pw_source_ex_plus_y(material_ex, space, cosine)
            
        elif self.directivity is const.PlusZ:
            cosine = dot(self.h_direction, (0, 1, 0))
            self._set_pw_source_ex_minus_z(material_ex, space, cosine)
            
        elif self.directivity is const.MinusZ:
            cosine = dot(self.h_direction, (0, 1, 0))
            self._set_pw_source_ex_plus_z(material_ex, space, cosine)
            
        else:
            return None
        
    def set_pointwise_source_ey(self, material_ey, space):
        if self.directivity is const.PlusZ:
            cosine = dot(self.h_direction, (1, 0, 0))
            self._set_pw_source_ey_minus_z(material_ey, space, cosine)
            
        elif self.directivity is const.MinusZ:
            cosine = dot(self.h_direction, (1, 0, 0))
            self._set_pw_source_ey_plus_z(material_ey, space, cosine)
            
        elif self.directivity is const.PlusX:
            cosine = dot(self.h_direction, (0, 0, 1))
            self._set_pw_source_ey_minus_x(material_ey, space, cosine)
            
        elif self.directivity is const.MinusX:
            cosine = dot(self.h_direction, (0, 0, 1))
            self._set_pw_source_ey_plus_x(material_ey, space, cosine)
            
        else:
            return None
        
    def set_pointwise_source_ez(self, material_ez, space):
        if self.directivity is const.PlusX:
            cosine = dot(self.h_direction, (0, 1, 0))
            self._set_pw_source_ez_minus_x(material_ez, space, cosine)
            
        elif self.directivity is const.MinusX:
            cosine = dot(self.h_direction, (0, 1, 0))
            self._set_pw_source_ez_plus_x(material_ez, space, cosine)
            
        elif self.directivity is const.PlusY:
            cosine = dot(self.h_direction, (1, 0, 0))
            self._set_pw_source_ez_minus_y(material_ez, space, cosine)
            
        elif self.directivity is const.MinusY:
            cosine = dot(self.h_direction, (1, 0, 0))
            self._set_pw_source_ez_plus_y(material_ez, space, cosine)
            
        else:
            return None
        
    def set_pointwise_source_hx(self, material_hx, space):
        if self.directivity is const.PlusY:
            cosine = dot(self.e_direction, (0, 0, 1))
            self._set_pw_source_hx_minus_y(material_hx, space, cosine)
            
        elif self.directivity is const.MinusY:
            cosine = dot(self.e_direction, (0, 0, 1))
            self._set_pw_source_hx_plus_y(material_hx, space, cosine)
            
        elif self.directivity is const.PlusZ:
            cosine = dot(self.e_direction, (0, 1, 0))
            self._set_pw_source_hx_minus_z(material_hx, space, cosine)
            
        elif self.directivity is const.MinusZ:
            cosine = dot(self.e_direction, (0, 1, 0))
            self._set_pw_source_hx_plus_z(material_hx, space, cosine)
            
        else:
            return None
        
    def set_pointwise_source_hy(self, material_hy, space):
        if self.directivity is const.PlusZ:
            cosine = dot(self.e_direction, (1, 0, 0))
            self._set_pw_source_hy_minus_z(material_hy, space, cosine)
            
        elif self.directivity is const.MinusZ:
            cosine = dot(self.e_direction, (1, 0, 0))
            self._set_pw_source_hy_plus_z(material_hy, space, cosine)
            
        elif self.directivity is const.PlusX:
            cosine = dot(self.e_direction, (0, 0, 1))
            self._set_pw_source_hy_minus_x(material_hy, space, cosine)
            
        elif self.directivity is const.MinusX:
            cosine = dot(self.e_direction, (0, 0, 1))
            self._set_pw_source_hy_plus_x(material_hy, space, cosine)
            
        else:
            return None
        
    def set_pointwise_source_hz(self, material_hz, space):
        if self.directivity is const.PlusX:
            cosine = dot(self.e_direction, (0, 1, 0))
            self._set_pw_source_hz_minus_x(material_hz, space, cosine)
            
        elif self.directivity is const.MinusX:
            cosine = dot(self.e_direction, (0, 1, 0))
            self._set_pw_source_hz_plus_x(material_hz, space, cosine)
            
        elif self.directivity is const.PlusY:
            cosine = dot(self.e_direction, (1, 0, 0))
            self._set_pw_source_hz_minus_y(material_hz, space, cosine)
            
        elif self.directivity is const.MinusY:
            cosine = dot(self.e_direction, (1, 0, 0))
            self._set_pw_source_hz_plus_y(material_hz, space, cosine)
            
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
        self.t = self.n * self.aux_fdtd.dt
        
    def envelope(self):
        width = self.aux_fdtd.src_list[0].src_time.width
        if self.t < width:
            env = sin(0.5 * pi * self.t / width)**2
        else:
            env = 1
        return env
