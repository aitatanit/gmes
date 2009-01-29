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
    def __init__(self, pos, component, src_time, amp=1):
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
                material_ex[idx] = DipoleEx(material_ex[idx], self.src_time, space.dt, self.amp)
            
    def set_pointwise_source_ey(self, material_ey, space):
        if self.comp is const.Ey:
            idx = space.space_to_ey_index(self.pos)
            if in_range(idx, material_ey, const.Ey):
                material_ey[idx] = DipoleEy(material_ey[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_ez(self, material_ez, space):
        if self.comp is const.Ez:
            idx = space.space_to_ez_index(self.pos)
            if in_range(idx, material_ez, const.Ez):
                material_ez[idx] = DipoleEz(material_ez[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_hx(self, material_hx, space):
        if self.comp is const.Hx:
            idx = space.space_to_hx_index(self.pos)
            if in_range(idx, material_hx, const.Hx):
                material_hx[idx] = DipoleHx(material_hx[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_hy(self, material_hy, space):
        if self.comp is const.Hy:
            idx = space.space_to_hy_index(self.pos)
            if in_range(idx, material_hy, const.Hy):
                material_hy[idx] = DipoleHy(material_hy[idx], self.src_time, space.dt, self.amp)
            
    def set_pointwise_source_hz(self, material_hz, space):
        if self.comp is const.Hz:
            idx = space.space_to_hz_index(self.pos)
            if in_range(idx, material_hz, const.Hz):
                material_hz[idx] = DipoleHz(material_hz[idx], self.src_time, space.dt, self.amp)


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
    
    """
    
    def __init__(self, direction, center, size, freq, polarization, width=inf, amp=1):
        """
        
        Arguments:
            direction -- propagation direction of the beam.
            center --
            size --
            freq --
            polarization --
            width -- the beam radius. The default is inf.
            amp -- amplitude of the plane wave. The default is 1.
            
        """
        
        # launching plane parameters
        self.direction = direction
        self.center = array(center, float)
        self.size = array(size, float)
        self.half_size = .5 * self.size
        self.freq = float(freq)
        self.polarization = polarization
        
        # spot size of Gaussian beam
        self.width = float(width)
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
        
    def init(self, geom_tree, space):
        self.geom_tree = geom_tree
        
    def _get_aux_fdtd(self, epsilon_r, mu_r, space):
        # two 10 meshes for the ABC,
        # 3 ex point (1 for the source) and 4 hy points for the free space
        aux_size = array((0 , 0, 24), float) / space.res[self.direction.tag]
        aux_space = Cartesian(size=aux_size, resolution=space.res[self.direction.tag], dt=space.dt, parallel=False)
        aux_geom_list = (DefaultMaterial(material=Dielectric(epsilon_r, mu_r)),
                         Boundary(material=CPML(epsilon_r, mu_r), thickness=10. / aux_space.res[2], size=aux_size))
        aux_src_list = (Dipole(src_time=Continuous(freq=self.freq), component=const.Ex, pos=(0,0,0)),)
        aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, aux_src_list, verbose=False)
        
        return aux_fdtd
        
    def set_pointwise_source_ex(self, material_ex, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        low_idx = space.space_to_ex_index(low)
        
        if self.direction is const.PlusY and self.polarization is const.X:
            TransparentEx = TransparentPlusYEx
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.MinusY and self.polarization is const.X:
            TransparentEx = TransparentMinusYEx
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.PlusZ and self.polarization is const.X:
            TransparentEx = TransparentPlusZEx
            distance_metric = self._distance_from_axis_in_z
            
        elif self.direction is const.MinusZ and self.polarization is const.X:
            TransparentEx = TransparentMinusZEx
            distance_metric = self._distance_from_axis_in_z
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_ex, const.Ex):
                        point = space.ex_index_to_space(i, j, k)
                        
                        r = distance_metric(point)
                        amp = self.amp * exp(-(r / self.width)**2)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, space)
                        
                        material_ex[i,j,k] = TransparentEx(material_ex[i,j,k], epsilon_r, amp, aux_fdtd)
        
    def set_pointwise_source_ey(self, material_ey, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        low_idx = space.space_to_ey_index(low)
        
        if self.direction is const.PlusZ and self.polarization is const.Y:
            TransparentEy = TransparentPlusZEy
            distance_metric = self._distance_from_axis_in_z
            
        elif self.direction is const.MinusZ and self.polarization is const.Y:
            TransparentEy = TransparentMinusZEy
            distance_metric = self._distance_from_axis_in_z
            
        elif self.direction is const.PlusX and self.polarization is const.Y:
            TransparentEy = TransparentPlusXEy
            distance_metric = self._distance_from_axis_in_x
            
        elif self.direction is const.MinusX and self.polarization is const.Y:
            TransparentEy = TransparentMinusXEy
            distance_metric = self._distance_from_axis_in_x
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_ey, const.Ey):
                        point = space.ey_index_to_space(i, j, k)
                        
                        r = distance_metric(point)
                        amp = self.amp * exp(-(r / self.width)**2)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, space)
                        
                        material_ey[i,j,k] = TransparentEy(material_ey[i,j,k], epsilon_r, amp, aux_fdtd)

    def set_pointwise_source_ez(self, material_ez, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        low_idx = space.space_to_ez_index(low)
        
        if self.direction is const.PlusY and self.polarization is const.Z:
            TransparentEz = TransparentPlusYEz
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.MinusY and self.polarization is const.Z:
            TransparentEz = TransparentMinusYEz
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.PlusX and self.polarization is const.Z:
            TransparentEz = TransparentPlusXEz
            distance_metric = self._distance_from_axis_in_x
            
        elif self.direction is const.MinusX and self.polarization is const.Z:
            TransparentEz = TransparentMinusXEz
            distance_metric = self._distance_from_axis_in_x
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_ez, const.Ez):
                        point = space.ez_index_to_space(i, j, k)
                        
                        r = distance_metric(point)
                        amp = self.amp * exp(-(r / self.width)**2)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, space)
                        
                        material_ez[i,j,k] = TransparentEz(material_ez[i,j,k], epsilon_r, amp, aux_fdtd)

    def set_pointwise_source_hx(self, material_hx, space):
        high = self.center + self.half_size
        low = self.center - self.half_size

        if self.direction is const.PlusY and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            amp_tmp = self.amp
            TansparentHx = TransparentPlusYHx
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.MinusY and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            amp_tmp = -self.amp
            TransparentHx = TransparentMinusYHx
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.PlusZ and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            amp_tmp = -self.amp
            TransparentHx = TransparentPlusZHx
            distance_metric = self._distance_from_axis_in_z
            
        elif self.direction is const.MinusZ and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            amp_tmp = self.amp
            TransparentHx = TransparentMinusZHx
            distance_metric = self._distance_from_axis_in_z
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_hx, const.Hx):
                        point = space.hx_index_to_space(i, j, k)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, space)
                        
                        r = distance_metric(point)
                        amp = amp_tmp * exp(-(r / self.width)**2) 

                        material_hx[i,j,k] = TransparentHx(material_hx[i,j,k], mu_r, amp, aux_fdtd)

    def set_pointwise_source_hy(self, material_hy, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.direction is const.PlusZ and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            
            amp_tmp = self.amp
            TransparentHy = TransparentPlusZHy
            distance_metric = self._distance_from_axis_in_z
            
        elif self.direction is const.MinusZ and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            amp_tmp = -self.amp
            TransparentHy = TransparentMinusZHy
            distance_metric = self._distance_from_axis_in_z
            
        elif self.direction is const.PlusX and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            amp_tmp = -self.amp
            TransparentHy = TransparentPlusXHy
            distance_metric = self._distance_from_axis_in_x
            
        elif self.direction is const.MinusX and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            amp_tmp = -self.amp
            TransparentHy = TransparentMinusXHy
            distance_metric = self._distance_from_axis_in_x
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_hy, const.Hy):
                        point = space.hy_index_to_space(i,j,k)
                        
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, space)
                        
                        r = distance_metric(point)
                        amp = amp_tmp * exp(-(r / self.width)**2)

                        material_hy[i,j,k] = TransparentHy(material_hy[i,j,k], mu_r, amp, aux_fdtd)
                        
    def set_pointwise_source_hz(self, material_hz, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.direction is const.PlusY and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
    
            amp_tmp = -self.amp
            TransparentHz = TransparentPlusYHz
            distance_metric = self._distance_from_axis_in_y
            
        elif self.direction is const.MinusY and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            amp_tmp = self.amp
            TransparentHz = TransparentMinusYHz
            distance_metric = self._distance_from_axis_in_y
                                   
        elif self.direction is const.PlusX and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            amp_tmp = self.amp
            TransparentHz = TransparentPlusXHz
            distance_metric = self._distance_from_axis_in_x
            
        elif self.direction is const.MinusX and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            amp_tmp = -self.amp
            TransparentHz = TransparentMinusXHz
            distance_metric = self._distance_from_axis_in_x
            
        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_hz, const.Hz):
                        point = space.hz_index_to_space(i, j, k)
                         
                        mat_objs = self.geom_tree.material_of_point(point)
                        epsilon_r = mat_objs[0].epsilon_r
                        mu_r = mat_objs[0].mu_r
                        aux_fdtd = self._get_aux_fdtd(epsilon_r, mu_r, space)
                        
                        r = distance_metric(point)
                        amp = amp_tmp * exp(-(r / self.width)**2)

                        material_hz[i,j,k] = TransparentHz(material_hz[i,j,k], mu_r, amp, self.aux_fdtd)
             
    def _distance_from_axis(self, point, axis_component):
        """Calculate distance from beam axis.
        
        Arguments:
            point -- point location in space coordinate
            axis_component -- constants.{X, Y, Z} 
        """
        
        if axis_component is const.X:
            return norm(point[1:] - self.center[1:])
        elif axis_component is const.Y:
            return norm(point[0:3:2] - self.center[0:3:2])
        elif axis_component is const.Z:
            return norm(point[:2] - self.center[:2])
        else:
            return None
        
    def _distance_from_axis_in_x(self, point):
        return self._distance_from_axis(point, const.X)
    
    def _distance_from_axis_in_y(self, point):
        return self._distance_from_axis(point, const.Y)
    
    def _distance_from_axis_in_z(self, point):
        return self._distance_from_axis(point, const.Z)
