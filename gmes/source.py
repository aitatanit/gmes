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

#
# SrcTime: Continuous, Bandpass
# Src: Dipole, GaussianBeam
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
            return 0.0
        
        if ts < self.width:
            env = sin(.5 * pi * ts / self.width)**2
        elif te < self.width:
            env = sin(.5 * pi * te / self.width)**2
        else:
            env = 1
            
        return env * cos(2 * pi * self.freq * time - self.phase)


class Bandpass(SrcTime):
    """a pulse source with Gaussian-envelope
    
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
    
    def display_info(self, indent=0):
        print " " * indent, "bandpass source"
        print " " * indent,
        print "center frequency:", self.freq,
        print "bandwidth:", self.fwidth,
        print "peak time:", self.peak
        print "cutoff:", self.cutoff
            
    def dipole(self, time):
        tt = time - self.peak
        if (abs(tt) > self.cutoff): return 0.0

        return exp(-.5 * (tt * self.fwidth)**2) * cos(2 * pi * self.freq * time)
        
        
class Dipole(Src):
    def __init__(self, src_time, pos, component, amp=1):
        self.pos = array(pos, float)
        self.comp = component
        self.src_time = src_time
        self.amp = float(amp)
    
    def init(self, geom_tree, space):
        pass
    
    def display_info(self, indent=0):
        print " " * indent, "Hertzian dipole source:"
        print " " * indent, "center:", self.pos
        print " " * indent, "polarization direction:", self.comp
        print " " * indent, "maximum amp.:", self.amp
        
        self.src_time.display_info(4)
        
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


class TotalFieldScatteredField(Src):
    """Set a total and scattered field zone to launch a plane wave.
    
    """
    def __init__(self, src_time, center, size, direction, polarization, amp=1):
        """Constructor
        
        Arguments:
            center -- center of the incidence interface. The beam axis crosses this point.
            size --  size of the incidence interface plane.
            direction -- propagation direction of the beam.
            freq -- oscillating frequency of the beam.
            polarization -- electric field direction of the beam.
            amp -- amplitude of the plane wave. The default is 1.
            
        """
        if isinstance(src_time, SrcTime):
            self.src_time = src_time
        else:
            raise TypeError, 'src_time must be an instance of SrcTime.'

        self.k = array(direction, float) / norm(direction)
        self.center = array(center, float)
        self.size = array(size, float)
        
        self.half_size = .5 * self.size
        self.e_direction = array(polarization, float) / norm(polarization)
        
        # direction of h field
        self.h_direction = cross(self.k, self.e_direction)
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
    
        self.on_axis_k = self._axis_in_k()
         
    def init(self, geom_tree, space):
        self.geom_tree = geom_tree
        
    def display_info(self, indent=0):
        print " " * indent, "plane-wave source:"
        print " " * indent, "propagation direction:", self.k
        print " " * indent, "center:", self.center
        print " " * indent, "source plane size:", self.size 
        print " " * indent, "polarization direction:", self.e_direction
        print " " * indent, "amp.:", self.amp
        
        self.src_time.display_info(4)

    def _dist_from_center(self, point):
        """Calculate distance from the interface plane center.
        
        Arguments:
            point -- location in the space coordinate
            
        """
        return norm(self.center - point)
    
    def  _dist_from_center_along_beam_axis(self, point):
        """Calculate projected distance from center along the beam axis.
        
        This method returns positive value when the point is located in
        the k direction relative to the center.
        
        Arguments:
            point -- location iin the space coordinate
            
        """
        return dot(self.k, point - self.center)
    
    def _axis_in_k(self):
        dot_with_axis = {}
        
        dot_with_axis[const.PlusX] = dot(const.PlusX.vector, self.k)
        dot_with_axis[const.PlusY] = dot(const.PlusY.vector, self.k)
        dot_with_axis[const.PlusY] = dot(const.PlusY.vector, self.k)
        dot_with_axis[const.MinusX] = dot(const.MinusX.vector, self.k)
        dot_with_axis[const.MinusY] = dot(const.MinusY.vector, self.k)
        dot_with_axis[const.MinusZ] = dot(const.MinusZ.vector, self.k)
        
        max_item =  max(dot_with_axis.items(), key=lambda item:item[1])
        return max_item[0]
    
    def _get_wave_number(self, k, epsilon_r, mu_r, space, error=1.e-3):
        """Return the numerical wave number for auxiliary fdtd.
        
        Arguments:
            k -- normalized wave vector
            epsilon_r -- relative permittivity which fills the auxiliary fdtd
            mu_r -- relative permeability which fills the auxiliary fdtd
            space -- Cartesian instance
            
        """
        ds = array((space.dx, space.dy, space.dz))
        dt = space.dt
        k1_scalar = inf
        k2_scalar = 2 * pi * self.src_time.freq
        
        while abs(k2_scalar - k1_scalar) > error:
            k1_scalar = k2_scalar
            
            f = (sum(((sin(.5 * k1_scalar * k * ds) / ds)**2)) - 
                         sqrt(epsilon_r * mu_r) * 
                         (sin(pi * self.src_time.freq * dt) / dt / const.c0)**2)
            f_prime = .5 * sum(k * sin(k1_scalar * k * ds) / ds)
            
            k2_scalar = k1_scalar - f / f_prime
            
        return k2_scalar

    def _get_aux_fdtd(self, epsilon_r, mu_r, dz, dt):
        border = self.center + .5 * self.size
        aux_long_size = \
        2 * abs(self._dist_from_center_along_beam_axis(border)) + 30 * dz
        src_pnt = (0, 0, -abs(self._dist_from_center_along_beam_axis(border)) - 5 * dz)
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
    
    def _set_pointwise_source(self, space, component, cosine, low_idx, high_idx, source, material):
        aux_ds = {const.PlusX: space.dx, const.MinusX: space.dx,  
                  const.PlusY: space.dy, const.MinusY: space.dy,
                  const.PlusZ: space.dz, const.MinusZ: space.dz}
        
        idx_to_spc = {const.Ex: space.ex_index_to_space,
                      const.Ey: space.ey_index_to_space,
                      const.Ez: space.ez_index_to_space,
                      const.Hx: space.hx_index_to_space,
                      const.Hy: space.hy_index_to_space,
                      const.Hz: space.hz_index_to_space}
        
        low_idx_array = array(low_idx)
        high_idx_array = array(high_idx)
        
        for i, j, k in ndindex(*(high_idx_array - low_idx_array)):
            i, j, k = (i, j, k) + low_idx_array
            
            if in_range((i, j, k), material, component):
                point = idx_to_spc[component](i, j, k)
                
                amp = cosine * self.amp
                
                mat_objs = self.geom_tree.material_of_point(point)
                epsilon_r = mat_objs[0].epsilon_r
                mu_r = mat_objs[0].mu_r
                aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                              aux_ds[self.on_axis_k], space.dt)
                
                samp_pnt = \
                (0, 0, self._dist_from_center_along_beam_axis(point))
                
                # v_in_axis / v_in_k
                v_ratio = \
                self._get_wave_number(self.k, epsilon_r, mu_r, space) / \
                self._get_wave_number(self.on_axis_k.vector, epsilon_r, mu_r, space)
                                        
                material[i,j,k] = source(material[i,j,k], epsilon_r,
                                         amp, aux_fdtd, samp_pnt,
                                         v_ratio)
    
    def set_pointwise_source_ex(self, material_ex, space):
        cosine = dot(self.e_direction, (1, 0, 0))
        
        if cosine == 0:
            return None
        
        # -y interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, -1, 1)
        
        low_idx = space.space_to_ex_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        
        pw_source = TransparentMinusYEx

        self._set_pointwise_source(space, const.Ex, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_ex)
        
        # +y interface
        low = self.center - self.half_size * (1, -1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ex_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        
        pw_source = TransparentPlusYEx

        self._set_pointwise_source(space, const.Ex, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ex)
        
        # -z interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, 1, -1)
        
        low_idx = space.space_to_ex_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        
        pw_source = TransparentMinusZEx

        self._set_pointwise_source(space, const.Ex, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ex)
        
        # +z interface
        low = self.center - self.half_size * (1, 1, -1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ex_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        
        pw_source = TransparentPlusZEx

        self._set_pointwise_source(space, const.Ex, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ex)
    
    def set_pointwise_source_ey(self, material_ey, space):
        cosine = dot(self.e_direction, (0, 1, 0))
        
        if cosine == 0:
            return None
        
        # -z interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, 1, -1)
        
        low_idx = space.space_to_ey_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        
        pw_source = TransparentMinusZEy

        self._set_pointwise_source(space, const.Ey, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ey)
        
        # +z interface
        low = self.center - self.half_size * (1, 1, -1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ey_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        
        pw_source = TransparentPlusZEy

        self._set_pointwise_source(space, const.Ey, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_ey)
        
        # -x interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (-1, 1, 1)
        
        low_idx = space.space_to_ey_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        
        pw_source = TransparentMinusXEy

        self._set_pointwise_source(space, const.Ey, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_ey)
        
        # +x interface
        low = self.center - self.half_size * (-1, 1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ey_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        
        pw_source = TransparentPlusXEy

        self._set_pointwise_source(space, const.Ey, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ey)
    
    def set_pointwise_source_ez(self, material_ez, space):
        cosine = dot(self.e_direction, (0, 0, 1))
        
        if cosine == 0:
            return None
        
        # -x interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (-1, 1, 1)
        
        low_idx = space.space_to_ez_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        
        pw_source = TransparentMinusXEz

        self._set_pointwise_source(space, const.Ez, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ez)
        
        # +x interface
        low = self.center - self.half_size * (-1, 1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ez_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        
        pw_source = TransparentPlusXEz

        self._set_pointwise_source(space, const.Ez, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ez)
        
        # -y interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, -1, 1)
        
        low_idx = space.space_to_ez_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        
        pw_source = TransparentMinusYEz

        self._set_pointwise_source(space, const.Ez, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ez)
        
        # +y interface
        low = self.center - self.half_size * (1, -1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ez_index(low)  
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        
        pw_source = TransparentPlusYEz

        self._set_pointwise_source(space, const.Ez, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_ez)
    
    def set_pointwise_source_hx(self, material_hx, space):
        cosine = dot(self.h_direction, (1, 0, 0))
        
        if cosine == 0:
            return None
        
        # -y interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, -1, 1)
        
        low_idx = space.space_to_ez_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        
        low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)    
        high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
        
        pw_source = TransparentMinusYHx
        
        self._set_pointwise_source(space, const.Hx, cosine,  
                                   low_idx, high_idx, 
                                   pw_source, material_hx)
        
        # +y interface
        low = self.center - self.half_size * (1, -1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ez_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        
        low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
        high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)    
            
        pw_source = TransparentPlusYHx
        
        self._set_pointwise_source(space, const.Hx, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hx)
        
        # -z interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, 1, -1)
        
        low_idx = space.space_to_ey_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        
        low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
        high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
        
        pw_source = TransparentMinusZHx
        
        self._set_pointwise_source(space, const.Hx, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hx)
        
        # +z interface
        low = self.center - self.half_size * (1, 1, -1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ey_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            
        low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
        high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
        
        pw_source = TransparentPlusZHx
            
        self._set_pointwise_source(space, const.Hx, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hx)
    
    def set_pointwise_source_hy(self, material_hy, space):
        cosine = dot(self.h_direction, (0, 1, 0))
        
        if cosine == 0:
            return None
        
        # -z interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, 1, -1)
        
        low_idx = space.space_to_ex_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            
        low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
        high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            
        pw_source = TransparentMinusZHy
        
        self._set_pointwise_source(space, const.Hy, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hy)
        
        # +z interface
        low = self.center - self.half_size * (1, 1, -1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ex_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            
        low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
        high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            
        pw_source = TransparentPlusZHy
        
        self._set_pointwise_source(space, const.Hy, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hy)
            
        # -x interface
        low = self.center - self.half_size * (-1, 1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ez_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            
        low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
        high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            
        pw_source = TransparentMinusXHy
        
        self._set_pointwise_source(space, const.Hy, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hy)
        
        # +x interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (-1, 1, 1)
        
        low_idx = space.space_to_ex_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            
        low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
        high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            
        pw_source = TransparentPlusXHy
        
        self._set_pointwise_source(space, const.Hy, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hy)
    
    def set_pointwise_source_hz(self, material_hz, space):
        cosine = dot(self.h_direction, (0, 0, 1))
        
        if cosine == 0:
            return None
        
        # -x interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (-1, 1, 1)
        
        low_idx = space.space_to_ey_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            
        low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
        high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            
        pw_source = TransparentMinusXHz
        
        self._set_pointwise_source(space, const.Hz, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hz)
        
        # +x interface
        low = self.center - self.half_size * (-1, 1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ey_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            
        low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
        high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
        pw_source = TransparentPlusXHz
        
        self._set_pointwise_source(space, const.Hz, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hz)
        
        # -y interface
        low = self.center - self.half_size
        high = self.center + self.half_size * (1, -1, 1)
        
        low_idx = space.space_to_ex_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            
        low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
        high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
        
        pw_source = TransparentMinusYHz
        
        self._set_pointwise_source(space, const.Hz, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hz)
            
        # +y interface
        low = self.center - self.half_size * (1, -1, 1)
        high = self.center + self.half_size
        
        low_idx = space.space_to_ex_index(low)
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            
        low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
        high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
        pw_source = TransparentPlusYHz
        
        self._set_pointwise_source(space, const.Hz, cosine,  
                                   low_idx, high_idx,
                                   pw_source, material_hz)
        
class GaussianBeam( TotalFieldScatteredField):
    """Launch a transparent Gaussian beam.
    
    It works as a guided mode with Gaussian profile is launched through the incidence interface.
    The incidence interface is transparent, thus the scattered wave can penetrate through the 
    interface plane.
    
    """
    def __init__(self, src_time, directivity, center, size, direction, polarization, waist=inf, amp=1):
        """
        
        Arguments:
            directivity -- directivity of the incidence interface.
            center -- center of the incidence interface. The beam axis crosses this point.
            size --  size of the incidence interface plane.
            direction -- propagation direction of the beam.
            freq -- oscillating frequency of the beam.
            polarization -- electric field direction of the beam.
            waist -- the Gaussian beam radius. The default is inf.
            amp -- amplitude of the plane wave. The default is 1.
            
        """
        TotalFieldScatteredField.__init__(self, src_time, center, size, direction, polarization, amp)
        
        if issubclass(directivity, const.Directional):
            self.directivity = directivity
        else:
            raise TypeError, 'directivity must be a Directional type.'
        
        # spot size of Gaussian beam
        self.waist = float(waist)
        
    def display_info(self, indent=0):
        print " " * indent, "Gaussian beam source:"
        print " " * indent, "propagation direction:", self.k
        print " " * indent, "center:", self.center
        print " " * indent, "source plane size:", self.size 
        print " " * indent, "polarization direction:", self.e_direction
        print " " * indent, "beam waist:", self.waist
        print " " * indent, "maximum amp.:", self.amp
        
        self.src_time.display_info(4)
            
    def _dist_from_beam_axis(self, point):
        """Calculate distance from the beam axis.
        
        Arguments:
            point -- location in the space coordinate 
            
        """   
        return norm(cross(self.k, point - self.center))
    
    def _set_pointwise_source(self, space, component, cosine, material, low_idx, high_idx, source, idx_to_spc):
        aux_ds = {const.PlusX: space.dx, const.MinusX: space.dx,  
                  const.PlusY: space.dy, const.MinusY: space.dy,
                  const.PlusZ: space.dz, const.MinusZ: space.dz}
        
        low_idx_array = array(low_idx)
        high_idx_array = array(high_idx)
        
        for i, j, k in ndindex(*(high_idx_array - low_idx_array)):
            i, j, k = (i, j, k) + low_idx_array
            
            if in_range((i, j, k), material, component):
                point = idx_to_spc(i, j, k)
                
                mat_objs = self.geom_tree.material_of_point(point)
                epsilon_r = mat_objs[0].epsilon_r
                mu_r = mat_objs[0].mu_r
                aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, 
                                              aux_ds[self.directivity], space.dt)
                
                r = self._dist_from_beam_axis(point)
                amp = cosine * self.amp * exp(-(r / self.waist)**2)
                
                samp_pnt = \
                (0, 0, self._dist_from_center_along_beam_axis(point))
                
                # v_in_axis / v_in_k
                v_ratio = \
                self._get_wave_number(self.k, epsilon_r, mu_r, space) / \
                self._get_wave_number(self.directivity.vector, epsilon_r, mu_r, space)

                material[i,j,k] = source(material[i,j,k], mu_r, amp, 
                                         aux_fdtd, samp_pnt, v_ratio)
                
    def set_pointwise_source_ex(self, material_ex, space):
        cosine = dot(self.e_direction, (1, 0, 0))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        low_idx = space.space_to_ex_index(low)
        
        if self.directivity is const.PlusY:
            pw_source = TransparentMinusYEx
            idx_to_spc = lambda i, j, k: space.hz_index_to_space(i + 1, j, k)
            
        elif self.directivity is const.MinusY:
            pw_source = TransparentPlusYEx
            idx_to_spc = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
        elif self.directivity is const.PlusZ:
            pw_source = TransparentMinusZEx
            idx_to_spc = lambda i, j, k: space.hy_index_to_space(i + 1, j, k)
            
        elif self.directivity is const.MinusZ:
            pw_source = TransparentPlusZEx
            idx_to_spc = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
        else:
            return None
            
        self._set_pointwise_source(space, const.Ex, cosine, material_ex, 
                                   low_idx, high_idx,
                                   pw_source, idx_to_spc)
        
    def set_pointwise_source_ey(self, material_ey, space):
        cosine = dot(self.e_direction, (0, 1, 0))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        low_idx = space.space_to_ey_index(low)
        
        if self.directivity is const.PlusZ:
            TransparentEy = TransparentMinusZEy
            idx_to_spc = lambda i, j, k: space.hx_index_to_space(i, j + 1, k)
            
        elif self.directivity is const.MinusZ:
            TransparentEy = TransparentPlusZEy
            idx_to_spc = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)
            
        elif self.directivity is const.PlusX:
            TransparentEy = TransparentMinusXEy
            idx_to_spc = lambda i, j, k: space.hz_index_to_space(i, j + 1, k)
            
        elif self.directivity is const.MinusX:
            TransparentEy = TransparentPlusXEy
            idx_to_spc = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
        else:
            return None
            
        self._set_pointwise_source(space, const.Ey, cosine, material_ey,
                                   low_idx, high_idx,
                                   TransparentEy, idx_to_spc)
        
    def set_pointwise_source_ez(self, material_ez, space):
        cosine = dot(self.e_direction, (0, 0, 1))
        
        if cosine == 0:
            return None
        
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        low_idx = space.space_to_ez_index(low)
        
        if self.directivity is const.PlusX:
            TransparentEz = TransparentMinusXEz
            idx_to_spc = lambda i, j, k: space.hy_index_to_space(i, j, k + 1)
            
        elif self.directivity is const.MinusX:
            TransparentEz = TransparentPlusXEz
            idx_to_spc = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
        elif self.directivity is const.PlusY:
            TransparentEz = TransparentMinusYEz
            idx_to_spc = lambda i, j, k: space.hx_index_to_space(i, j, k + 1)
            
        elif self.directivity is const.MinusY:
            TransparentEz = TransparentPlusYEz
            idx_to_spc = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)
            
        else:
            return None
        
        self._set_pointwise_source(space, const.Ez, cosine, material_ez,
                                   low_idx, high_idx, 
                                   TransparentEz, idx_to_spc)
        
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
            
            TransparentHx = TransparentMinusYHx
            idx_to_spc = lambda i, j, k: space.ez_index_to_space(i, j - 1, k - 1)
            
        elif self.directivity is const.MinusY:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            TransparentHx = TransparentPlusYHx
            idx_to_spc = lambda i, j, k: space.ez_index_to_space(i, j, k - 1)
            
        elif self.directivity is const.PlusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            TransparentHx = TransparentMinusZHx
            idx_to_spc = lambda i, j, k: space.ey_index_to_space(i, j - 1, k - 1)
            
        elif self.directivity is const.MinusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            TransparentHx = TransparentPlusZHx
            idx_to_spc = lambda i, j, k: space.ey_index_to_space(i, j - 1, k)
            
        else:
            return None

        self._set_pointwise_source(space, const.Hx, cosine, material_hx,
                                   low_idx, high_idx,
                                   TransparentHx, idx_to_spc)
        
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
            
            TransparentHy = TransparentMinusZHy
            idx_to_spc = lambda i, j, k: space.ex_index_to_space(i - 1, j, k - 1)
            
        elif self.directivity is const.MinusZ:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            TransparentHy = TransparentPlusZHy
            idx_to_spc = lambda i, j, k: space.ex_index_to_space(i - 1, j, k)
            
        elif self.directivity is const.PlusX:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            TransparentHy = TransparentMinusXHy
            idx_to_spc = lambda i, j, k: space.ez_index_to_space(i - 1, j, k - 1)
            
        elif self.directivity is const.MinusX:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            TransparentHy = TransparentPlusXHy
            idx_to_spc = lambda i, j, k: space.ez_index_to_space(i, j, k - 1)
            
        else:
            return None
        
        self._set_pointwise_source(space, const.Hy, cosine, material_hy,
                                   low_idx, high_idx,
                                   TransparentHy, idx_to_spc)
        
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
            
            TransparentHz = TransparentMinusXHz
            idx_to_spc = lambda i, j, k: space.ey_index_to_space(i - 1, j - 1, k)
            
        elif self.directivity is const.MinusX:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            TransparentHz = TransparentPlusXHz
            idx_to_spc = lambda i, j, k: space.ey_index_to_space(i, j - 1, k)
            
        elif self.directivity is const.PlusY:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
    
            TransparentHz = TransparentMinusYHz
            idx_to_spc = lambda i, j, k: space.ex_index_to_space(i - 1, j - 1, k)
            
        elif self.directivity is const.MinusY:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            TransparentHz = TransparentPlusYHz
            idx_to_spc = lambda i, j, k: space.ex_index_to_space(i - 1, j, k)
            
        else:
            return None
        
        self._set_pointwise_source(space, const.Hz, cosine, material_hz,
                                   low_idx, high_idx, 
                                   TransparentHz, idx_to_spc)
        