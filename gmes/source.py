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
    
    def step(self):
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
        
        self.aux_fdtd = self._get_aux_fdtd(space)
        
    def step(self):
        self.aux_fdtd.step()
        
    def display_info(self, indent=0):
        print " " * indent, "plane-wave source:"
        print " " * indent, "propagation direction:", self.k
        print " " * indent, "center:", self.center
        print " " * indent, "source plane size:", self.size 
        print " " * indent, "polarization direction:", self.e_direction
        print " " * indent, "amp.:", self.amp
        
        self.src_time.display_info(4)
        
    def mode_function(self, x, y, z):
        return 1.0
        
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
        dot_with_axis[dot(const.PlusY.vector, self.k)] = const.PlusY 
        dot_with_axis[dot(const.MinusX.vector, self.k)] = const.MinusX 
        dot_with_axis[dot(const.MinusY.vector, self.k)] = const.MinusY 
        dot_with_axis[dot(const.MinusZ.vector, self.k)] = const.MinusZ 
        
        return dot_with_axis[max(dot_with_axis)]
        
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

    def _get_aux_fdtd(self, space):
        """Returns a TEMz FDTD for a reference of a plane wave. 
        
        """
        aux_ds = {const.PlusX: space.dx, const.MinusX: space.dx,  
                  const.PlusY: space.dy, const.MinusY: space.dy,
                  const.PlusZ: space.dz, const.MinusZ: space.dz}
        
        dz = aux_ds[self.on_axis_k]
        
        pml_thickness = 10 * dz
        border = self.center + .5 * self.size
        
        longitudinal_size = \
        2 * (abs(self._dist_from_center_along_beam_axis(border)) + pml_thickness + dz)
        aux_size = (0, 0, longitudinal_size)
        
        mat_objs =  self.geom_tree.material_of_point((inf, inf, inf))
        
        aux_space = Cartesian(size=aux_size, 
                              resolution=1./dz, 
                              dt=space.dt, 
                              parallel=False)
        aux_geom_list = (DefaultMaterial(material=mat_objs[0]),
                         Boundary(material=CPML(), 
                                  thickness=pml_thickness, 
                                  size=aux_size,
                                  minus_z=False))
        src_pnt = aux_space.ex_index_to_space(0, 0, 0)
        aux_src_list = (Dipole(src_time=deepcopy(self.src_time), 
                               component=const.Ex, 
                               pos=src_pnt),)
        aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, 
                            aux_src_list, 
                            verbose=False)
        
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
        
        low_idx_array = array(low_idx)
        high_idx_array = array(high_idx)
        
        for i, j, k in ndindex(*(high_idx_array - low_idx_array)):
            i, j, k = (i, j, k) + low_idx_array
            
            if in_range((i, j, k), material, component):
                point = idx_to_spc[component](i, j, k)
                
                mat_objs = self.geom_tree.material_of_point(point)
                epsilon_r = mat_objs[0].epsilon_r
                mu_r = mat_objs[0].mu_r

                amp = cosine * self.amp * self.mode_function(*point)

                samp_pnt = (0, 0, 
                            self._dist_from_center_along_beam_axis(samp_i2s(i, j, k)))
                
                # v_in_axis / v_in_k
                v_ratio = \
                self._get_wave_number(self.k, epsilon_r, mu_r, space) / \
                self._get_wave_number(self.on_axis_k.vector, epsilon_r, mu_r, space)
                
                if issubclass(source, TransparentElectric):
                    material[i, j, k] = source(material[i, j, k], epsilon_r,
                                               amp, self.aux_fdtd, samp_pnt,
                                               v_ratio, face)
                if issubclass(source, TransparentMagnetic):
                    material[i, j, k] = source(material[i, j, k], mu_r, 
                                               amp, self.aux_fdtd, samp_pnt, 
                                               v_ratio, face)
        
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
            
            low_idx = space.space_to_ex_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j, k)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex, 
                                low_idx, high_idx, TransparentEx, 
                                hz_i2s, const.MinusY)
        
    def _set_pw_source_ex_plus_y(self, material_ex, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i + 1, j + 1, k)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex,
                                low_idx, high_idx, TransparentEx, 
                                hz_i2s, const.PlusY)
        
    def _set_pw_source_ex_minus_z(self, material_ex, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, 1, -1)
            
            low_idx = space.space_to_ex_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
    
            hy_i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k)
            
            self._set_pw_source(space, const.Ex, cosine, material_ex,
                                low_idx, high_idx, TransparentEx, 
                                hy_i2s, const.MinusZ)
    
    def _set_pw_source_ex_plus_z(self, material_ex, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ex_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            
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
            
            low_idx = space.space_to_ey_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,
                                low_idx, high_idx, TransparentEy, 
                                hx_i2s, const.MinusZ)
    
    def _set_pw_source_ey_plus_z(self, material_ey, space, cosine):
        if 2 * space.half_size[2] > space.dz:
            low = self.center - self.half_size * (1, 1, -1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j + 1, k + 1)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,
                                low_idx, high_idx, TransparentEy,
                                hx_i2s, const.PlusZ)
    
    def _set_pw_source_ey_minus_x(self, material_ey, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size
            high = self.center + self.half_size * (-1, 1, 1)
            
            low_idx = space.space_to_ey_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
    
            hz_i2s = lambda i, j, k: space.hz_index_to_space(i, j + 1, k)
            
            self._set_pw_source(space, const.Ey, cosine, material_ey,  
                                low_idx, high_idx, TransparentEy,
                                hz_i2s, const.MinusX)
    
    def _set_pw_source_ey_plus_x(self, material_ey, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ey_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
    
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
            
            low_idx = space.space_to_ez_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
    
            hy_i2s = lambda i, j, k: space.hy_index_to_space(i, j, k + 1)
            
            self._set_pw_source(space, const.Ez, cosine, material_ez,  
                                low_idx, high_idx, TransparentEz, 
                                hy_i2s, const.MinusX)
    
    def _set_pw_source_ez_plus_x(self, material_ez, space, cosine):
        if 2 * space.half_size[0] > space.dx:
            low = self.center - self.half_size * (-1, 1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
    
            hy_i2s = lambda i, j, k: space.hy_index_to_space(i + 1, j, k + 1)
            
            self._set_pw_source(space, const.Ez, cosine, material_ez,
                                low_idx, high_idx, TransparentEz,
                                hy_i2s, const.PlusX)
    
    def _set_pw_source_ez_minus_y(self, material_ez, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size
            high = self.center + self.half_size * (1, -1, 1)
            
            low_idx = space.space_to_ez_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
    
            hx_i2s = lambda i, j, k: space.hx_index_to_space(i, j, k + 1)
            
            self._set_pw_source(space, const.Ez, cosine, material_ez,
                                low_idx, high_idx, TransparentEz, 
                                hx_i2s, const.MinusY)
    
    def _set_pw_source_ez_plus_y(self, material_ez, space, cosine):
        if 2 * space.half_size[1] > space.dy:
            low = self.center - self.half_size * (1, -1, 1)
            high = self.center + self.half_size
            
            low_idx = space.space_to_ez_index(low)  
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
    
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
            
            low_idx = space.space_to_ez_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            
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
            
            low_idx = space.space_to_ez_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            
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
            
            low_idx = space.space_to_ey_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            
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
            
            ey_low_idx = space.space_to_ey_index(low)
            ey_high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            
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
            
            low_idx = space.space_to_ex_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
                
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
            
            low_idx = space.space_to_ex_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
                
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
            
            ez_low_idx = space.space_to_ez_index(low)
            ez_high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
                
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
            
            ez_low_idx = space.space_to_ez_index(low)
            ez_high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
                
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
            
            low_idx = space.space_to_ey_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
                
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
            
            low_idx = space.space_to_ey_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
                
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
            
            low_idx = space.space_to_ex_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
                
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
            
            low_idx = space.space_to_ex_index(low)
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
                
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            
            ex_2s = lambda i, j, k: space.ex_index_to_space(i - 1, j - 1, k)
            
            self._set_pw_source(space, const.Hz, cosine, material_hz,
                                low_idx, high_idx, TransparentHz,
                                ex_i2s, const.PlusY)    
        
        
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
        