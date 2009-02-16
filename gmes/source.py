#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from copy import deepcopy
from numpy import *
from numpy.linalg import norm

from pointwise_source import *
import constants as const

from geometry import Cartesian, DefaultMaterial, Boundary, in_range
from fdtd import TEMzFDTD
from material import Dielectric, CPML


class SrcTime(object):
    """Time-dependent part of a source.
    
    """
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
        
    def dipole(self, time):
        ts = time - self.start
        te = self.end - time
        
        if ts < 0 or te < 0:
            return 0.0
        
        if ts < self.width:
            env = sin(.5 * pi * ts / self.width)**2
        elif te < self.width:
            env = sin(.5 * pi * te / self.width)**2
        else:
            env = 1
            
        return env * cos(2 * pi * self.freq * time - self.phase)
    
    def display_info(self, indent=0):
        print " " * indent, "continuous source"
        print " " * indent,
        print "frequency:", self.freq,
        print "initial phase advance:", self.phase,
        print "start time:", self.start,
        print "end time:", self.end,
        print "raising duration:", self.width


class Bandpass(SrcTime):
    """Gaussian-envelope source.
    
    """
    def __init__(self, freq, fwidth):
        self.freq = float(freq)
        self.fwidth = float(fwidth)
        width = 1 / self.fwidth
        s = 5
        self.peak = width * s
        self.cutoff = 2 * width * s
        
        while exp(-.5 * self.cutoff * self.cutoff / width / width) == 0:
            self.cutoff *= .9
        
        if self.peak - self.cutoff < 0:
            self.peak += self.cutoff - self.peak
    
    def dipole(self, time):
        tt = time - self.peak
        if (abs(tt) > self.cutoff): return 0.0

        return exp(-.5 * (tt * self.fwidth)**2) * cos(2 * pi * self.freq * time)
    
    def display_info(self, indent=0):
        print " " * indent, "bandpass source"
        print " " * indent,
        print "center frequency:", self.freq,
        print "bandwidth:", self.fwidth,
        print "peak time:", self.peak
        print "cutoff:", self.cutoff
        
        
class Dipole(object):
    def __init__(self, src_time, pos, component, amp=1):
        self.pos = array(pos, float)
        self.comp = component
        self.src_time = src_time
        self.amp = float(amp)
    
    def init(self, geom_tree, space):
        pass
    
    def set_pointwise_source_ex(self, material_ex, space):
        if self.comp is const.Ex:
            idx = space.space_to_ex_index(self.pos)
            if in_range(idx, material_ex, const.Ex):
                material_ex[idx] = DipoleEx(material_ex[idx], self.src_time, 
                                            space.dt, self.amp)
            
    def set_pointwise_source_ey(self, material_ey, space):
        if self.comp is const.Ey:
            idx = space.space_to_ey_index(self.pos)
            if in_range(idx, material_ey, const.Ey):
                material_ey[idx] = DipoleEy(material_ey[idx], self.src_time, 
                                            space.dt, self.amp)
   
    def set_pointwise_source_ez(self, material_ez, space):
        if self.comp is const.Ez:
            idx = space.space_to_ez_index(self.pos)
            if in_range(idx, material_ez, const.Ez):
                material_ez[idx] = DipoleEz(material_ez[idx], self.src_time, 
                                            space.dt, self.amp)
   
    def set_pointwise_source_hx(self, material_hx, space):
        if self.comp is const.Hx:
            idx = space.space_to_hx_index(self.pos)
            if in_range(idx, material_hx, const.Hx):
                material_hx[idx] = DipoleHx(material_hx[idx], self.src_time, 
                                            space.dt, self.amp)
   
    def set_pointwise_source_hy(self, material_hy, space):
        if self.comp is const.Hy:
            idx = space.space_to_hy_index(self.pos)
            if in_range(idx, material_hy, const.Hy):
                material_hy[idx] = DipoleHy(material_hy[idx], self.src_time, 
                                            space.dt, self.amp)
            
    def set_pointwise_source_hz(self, material_hz, space):
        if self.comp is const.Hz:
            idx = space.space_to_hz_index(self.pos)
            if in_range(idx, material_hz, const.Hz):
                material_hz[idx] = DipoleHz(material_hz[idx], self.src_time, 
                                            space.dt, self.amp)


class TotalFieldScatteredField(object):
    def __init__(self, theta, phi, psi, low, high):
        self.theta = float(theta)
        self.phi = float(theta)
        self.psi = float(psi)
        self.low = array(low, float)
        self.high = array(high, float)

    def generate_aux_fdtd(self, space):
        pass
    
    
class GaussianBeam(object):
    """Launch a transparent Gaussian beam.
    
    It works as a guided mode with Gaussian profile is launched through the incidence interface.
    The incidence interface is transparent, thus the scattered wave can penetrate through the 
    interface plane.
    
    """
    def __init__(self, src_time, directivity, center, size, direction, polarization, width=inf, amp=1):
        """
        
        Arguments:
            directivity -- directivity of the incidence interface.
            center -- center of the incidence interface. The beam axis crosses this point.
            size --  size of the incidence interface plane.
            direction -- propagation direction of the beam.
            freq -- oscillating frequency of the beam.
            polarization -- electric field direction of the beam.
            width -- the Gaussian beam radius. The default is inf.
            amp -- amplitude of the plane wave. The default is 1.
            
        """
        if isinstance(src_time, SrcTime):
            self.src_time = src_time
        else:
            raise TypeError, 'src_time must be an instance of SrcTime.'
           
        if issubclass(directivity, const.Directional):
            self.directivity = directivity
        else:
            raise TypeError, 'directivity must be a Directional type.' 

        self.k = array(direction, float) / norm(direction)
        self.center = array(center, float)
        self.size = array(size, float)
        
        self.half_size = .5 * self.size
        self.e_direction = array(polarization, float) / norm(polarization)
        
        # direction of h field
        self.h_direction = cross(self.k, self.e_direction)
        
        # spot size of Gaussian beam
        self.width = float(width)
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
        
    def init(self, geom_tree, space):
        self.geom_tree = geom_tree
        
    def _get_wave_number(self, k, epsilon_r, mu_r, space):
        """Return the numerical wave number for auxiliary fdtd.
        
        Arguments:
            k -- normalized wave vector
            epsilon_r -- relative permittivity which fills the auxiliary fdtd
            mu_r -- relative permeability which fills the auxiliary fdtd
            space -- Cartesian instance
            
        """
        k = array(k, float)
        dx = array((space.dx, space.dy, space.dz))
        dt = space.dt
        k1_scalar = 0
        k2_scalar = 1
        
        while k2_scalar - k1_scalar > 0.1:
            k1_scalar = k2_scalar
            
            numerator = 2 * (sum(((sin(.5 * k1_scalar * k * dx) / dx)**2)) -
                             epsilon_r * mu_r * 
                             (sin(pi * self.src_time.freq * dt) / dt)**2)
            denominator = sum(k1_scalar * sin(k1_scalar * k * dx) / dx)
            
            k2_scalar = k1_scalar - numerator / denominator
            
        return k2_scalar
    
    def _get_aux_fdtd(self, epsilon_r, mu_r, dz, dt):
        border = self.center + self.size
        aux_long_size = \
        2 * self._dist_from_center_along_aux_fdtd(border) + 30 * dz
        src_pnt = (0, 0, -self._dist_from_center_along_aux_fdtd(border) - 5 * dz)        
        aux_size = (0 , 0, aux_long_size)

        aux_space = Cartesian(size=aux_size, 
                              resolution=1./dz, 
                              dt=dt, 
                              parallel=False)
        aux_geom_list = (DefaultMaterial(material=Dielectric(epsilon_r, mu_r)),
                         Boundary(material=CPML(epsilon_r, mu_r), 
                                  thickness=10*dz, 
                                  size=aux_size))
        aux_src_list = (Dipole(src_time=deepcopy(self.src_time), 
                               component=const.Ex, 
                               pos=src_pnt),)
        aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, 
                            aux_src_list, 
                            verbose=False)
        
        return aux_fdtd
        
    def set_pointwise_source_ex(self, material_ex, space):
        cosine = dot(self.e_direction, (1, 0, 0))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        low_idx = space.space_to_ex_index(low)
        
        if self.directivity is const.PlusY:
            TransparentEx = TransparentPlusYEx
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        elif self.directivity is const.MinusY:
            TransparentEx = TransparentMinusYEx
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        elif self.directivity is const.PlusZ:
            TransparentEx = TransparentPlusZEx
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        elif self.directivity is const.MinusZ:
            TransparentEx = TransparentMinusZEx
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i, j, k), material_ex, const.Ex):
                        point = space.ex_index_to_space(i, j, k)
                        
                        r = self._dist_from_beam_axis(point)
                        amp = cosine * self.amp * exp(-(r / self.width)**2)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                                      aux_dz, space.dt)
                        
                        samp_pnt = \
                        (0, 0, self._dist_from_center_along_aux_fdtd(point))
                        
                        # v_in_axis / v_in_k
                        v_ratio = \
                        self._get_wave_number(self.k, epsilon_r, mu_r, space) /\
                        self._get_wave_number(in_axis_k, epsilon_r, mu_r, space)
                                                
                        material_ex[i,j,k] = TransparentEx(material_ex[i,j,k], 
                                                           epsilon_r,
                                                           amp, 
                                                           aux_fdtd, 
                                                           samp_pnt,
                                                           v_ratio)
        
    def set_pointwise_source_ey(self, material_ey, space):
        cosine = dot(self.e_direction, (0, 1, 0))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        low_idx = space.space_to_ey_index(low)
        
        if self.directivity is const.PlusZ:
            TransparentEy = TransparentPlusZEy
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        elif self.directivity is const.MinusZ:
            TransparentEy = TransparentMinusZEy
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        elif self.directivity is const.PlusX:
            TransparentEy = TransparentPlusXEy
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        elif self.directivity is const.MinusX:
            TransparentEy = TransparentMinusXEy
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i, j, k), material_ey, const.Ey):
                        point = space.ey_index_to_space(i, j, k)
                        
                        r = self._dist_from_beam_axis(point)
                        amp = cosine * self.amp * exp(-(r / self.width)**2)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                                      aux_dz, space.dt)
                        
                        samp_pnt = \
                        (0, 0, self._dist_from_center_along_aux_fdtd(point))
                        
                        # v_in_axis / v_in_k
                        v_ratio = \
                        self._get_wave_number(self.k, epsilon_r, mu_r, space) /\
                        self._get_wave_number(in_axis_k, epsilon_r, mu_r, space)
                        
                        material_ey[i,j,k] = TransparentEy(material_ey[i,j,k], 
                                                           epsilon_r, 
                                                           amp, 
                                                           aux_fdtd, 
                                                           samp_pnt,
                                                           v_ratio)

    def set_pointwise_source_ez(self, material_ez, space):
        cosine = dot(self.e_direction, (0, 0, 1))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        low_idx = space.space_to_ez_index(low)
        
        if self.directivity is const.PlusX:
            TransparentEz = TransparentPlusXEz
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        elif self.directivity is const.MinusX:
            TransparentEz = TransparentMinusXEz
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        elif self.directivity is const.PlusY:
            TransparentEz = TransparentPlusYEz
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        elif self.directivity is const.MinusY:
            TransparentEz = TransparentMinusYEz
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i, j, k), material_ez, const.Ez):
                        point = space.ez_index_to_space(i, j, k)
                        
                        r = self._dist_from_beam_axis(point)
                        amp = cosine * self.amp * exp(-(r / self.width)**2)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                                      aux_dz, space.dt)
                        
                        samp_pnt = \
                        (0, 0, self._dist_from_center_along_aux_fdtd(point))
                        
                        # v_in_axis / v_in_k
                        v_ratio = \
                        self._get_wave_number(self.k, epsilon_r, mu_r, space) /\
                        self._get_wave_number(in_axis_k, epsilon_r, mu_r, space)
                                                
                        material_ez[i,j,k] = TransparentEz(material_ez[i,j,k], 
                                                           epsilon_r, 
                                                           amp, 
                                                           aux_fdtd, 
                                                           samp_pnt,
                                                           v_ratio)

    def set_pointwise_source_hx(self, material_hx, space):
        cosine = dot(self.h_direction, (1, 0, 0))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size

        if self.directivity is const.PlusY:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            TransparentHx = TransparentPlusYHx
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        elif self.directivity is const.MinusY:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            TransparentHx = TransparentMinusYHx
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        elif self.directivity is const.PlusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            TransparentHx = TransparentPlusZHx
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        elif self.directivity is const.MinusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            TransparentHx = TransparentMinusZHx
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i, j, k), material_hx, const.Hx):
                        point = space.hx_index_to_space(i, j, k)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                                      aux_dz, space.dt)
                        
                        r = self._dist_from_beam_axis(point)
                        amp = cosine * self.amp * exp(-(r / self.width)**2) 
                        
                        samp_pnt = \
                        (0, 0, self._dist_from_center_along_aux_fdtd(point))
                        
                        # v_in_axis / v_in_k
                        v_ratio = \
                        self._get_wave_number(self.k, epsilon_r, mu_r, space) /\
                        self._get_wave_number(in_axis_k, epsilon_r, mu_r, space)
                                                
                        material_hx[i,j,k] = TransparentHx(material_hx[i,j,k], 
                                                           mu_r, 
                                                           amp, 
                                                           aux_fdtd, 
                                                           samp_pnt,
                                                           v_ratio)

    def set_pointwise_source_hy(self, material_hy, space):
        cosine = dot(self.h_direction, (0, 1, 0))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.directivity is const.PlusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            
            TransparentHy = TransparentPlusZHy
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        elif self.directivity is const.MinusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            TransparentHy = TransparentMinusZHy
            aux_dz = space.dz
            in_axis_k = const.PlusZ.vector
            
        elif self.directivity is const.PlusX:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            TransparentHy = TransparentPlusXHy
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        elif self.directivity is const.MinusX:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            TransparentHy = TransparentMinusXHy
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i, j, k), material_hy, const.Hy):
                        point = space.hy_index_to_space(i, j, k)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                                      aux_dz, space.dt)
                        
                        r = self._dist_from_beam_axis(point)
                        amp = cosine * self.amp * exp(-(r / self.width)**2)
                        
                        samp_pnt = \
                        (0, 0, self._dist_from_center_along_aux_fdtd(point))
                        
                        # v_in_axis / v_in_k
                        v_ratio = \
                        self._get_wave_number(self.k, epsilon_r, mu_r, space) /\
                        self._get_wave_number(in_axis_k, epsilon_r, mu_r, space)
                                                
                        material_hy[i,j,k] = TransparentHy(material_hy[i,j,k], 
                                                           mu_r, 
                                                           amp, 
                                                           aux_fdtd, 
                                                           samp_pnt,
                                                           v_ratio)
                        
    def set_pointwise_source_hz(self, material_hz, space):
        cosine = dot(self.h_direction, (0, 0, 1))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.directivity is const.PlusX:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            TransparentHz = TransparentPlusXHz
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        elif self.directivity is const.MinusX:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            TransparentHz = TransparentMinusXHz
            aux_dz = space.dx
            in_axis_k = const.PlusX.vector
            
        elif self.directivity is const.PlusY:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
    
            TransparentHz = TransparentPlusYHz
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        elif self.directivity is const.MinusY:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            TransparentHz = TransparentMinusYHz
            aux_dz = space.dy
            in_axis_k = const.PlusY.vector
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i, j, k), material_hz, const.Hz):
                        point = space.hz_index_to_space(i, j, k)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                                      aux_dz, space.dt)
                        
                        r = self._dist_from_beam_axis(point)
                        amp = cosine * self.amp * exp(-(r / self.width)**2)
                        
                        samp_pnt = \
                        (0, 0, self._dist_from_center_along_aux_fdtd(point))
                        
                        # v_in_axis / v_in_k
                        v_ratio = \
                        self._get_wave_number(self.k, epsilon_r, mu_r, space) /\
                        self._get_wave_number(in_axis_k, epsilon_r, mu_r, space)
                                                
                        material_hz[i,j,k] = TransparentHz(material_hz[i,j,k],
                                                           mu_r, 
                                                           amp,
                                                           aux_fdtd, 
                                                           samp_pnt,
                                                           v_ratio)
                        
    def _dist_from_beam_axis(self, point):
        """Calculate distance from the beam axis.
        
        Arguments:
            point -- location in the space coordinate 
            
        """   
        return norm(cross(self.k, self.center - point))

    def _dist_from_center(self, point):
        """Calculate distance from the interface plane center.
        
        Arguments:
            point -- location in the space coordinate
            
        """
        return norm(self.center - point)
    
    def _dist_from_center_along_aux_fdtd(self, point):
        """Calculate projected distance from center along the aux_fdtd.
        
        Arguments:
            point -- location iin the space coordinate
            
        """
        return dot(self.k, point - self.center)
    