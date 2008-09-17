#!/usr/bin/env python

from copy import deepcopy
from math import cos, sin, exp, pi
from numpy import inf, array

from pointwise_source import *
import constants as const

from geometric import Cartesian, DefaultMaterial, Boundary
from fdtd import TEMzFDTD
from material import Dielectric, UPML


class SrcTime:
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
        
        
class Dipole:
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
            material_ex[idx] = DipoleEx(material_ex[idx], self.src_time, space.dt, self.amp)
            
    def set_pointwise_source_ey(self, material_ey, space):
        if self.comp is const.Ey:
            idx = space.space_to_ey_index(self.pos)
            material_ey[idx] = DipoleEy(material_ey[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_ez(self, material_ez, space):
        if self.comp is const.Ez:
            idx = space.space_to_ez_index(self.pos)
            material_ez[idx] = DipoleEz(material_ez[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_hx(self, material_hx, space):
        if self.comp is const.Hx:
            idx = space.space_to_hx_index(self.pos)
            material_hx[idx] = DipoleHx(material_hx[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_hy(self, material_hy, space):
        if self.comp is const.Hy:
            idx = space.space_to_hy_index(self.pos)
            material_hy[idx] = DipoleHy(material_hy[idx], self.src_time, space.dt, self.amp)
            
    def set_pointwise_source_hz(self, material_hz, space):
        if self.comp is const.Hz:
            idx = space.space_to_hz_index(self.pos)
            material_hz[idx] = DipoleHz(material_hz[idx], self.src_time, space.dt, self.amp)


class TotalFieldScatteredField:
    def __init__(self, theta, phi, psi, low, high):
        self.theta = float(theta)
        self.phi = float(theta)
        self.psi = float(psi)
        self.low = array(low, float)
        self.high = array(high, float)

    def generate_aux_fdtd(self, space):
        pass
    
    
class Transparent:
    def __init__(self, direction, center, size, freq, polarization, amp=1,
                 ex_mode_file=None, ey_mode_file=None, ez_mode_file=None,
                 hx_mode_file=None, hy_mode_file=None, hz_mode_file=None):
        # launching plane parameters
        self.direction = direction
        self.center = array(center, float)
        self.size = array(size, float)
        self.half_size = .5 * self.size
        self.freq = float(freq)
        self.polarization = polarization
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
        
    def init(self, geom_tree, space):
        go = geom_tree.object_of_point(self.center)
        self.epsilon_r = go.material.epsilon_r
        self.mu_r = go.material.mu_r
        self.aux_fdtd = self.get_aux_fdtd(space)
        
    def get_aux_fdtd(self, space):
        # two 10 meshes for the ABC,
        # 1 ex point and 2 hy points for the free space
        aux_size = array((0 , 0, 21), float) / space.res
        aux_space = Cartesian(size=aux_size, resolution=space.res)
        aux_geom_list = (DefaultMaterial(material=Dielectric(self.epsilon_r, self.mu_r)),
                         Boundary(material=UPML(self.epsilon_r, self.mu_r), thickness=10 / space.res[2], size=aux_size))
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
        
        elif self.direction is const.MinusY and self.polarization is const.X:
            TransparentEx = TransparentMinusYEx
            
        elif self.direction is const.PlusZ and self.polarization is const.X:
            TransparentEx = TransparentPlusZEx
            
        elif self.direction is const.MinusZ and self.polarization is const.X:
            TransparentEx = TransparentMinusZEx

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    material_ex[i,j,k] = TransparentEx(material_ex[i,j,k], self.epsilon_r, self.amp, deepcopy(self.aux_fdtd))
        
    def set_pointwise_source_ey(self, material_ey, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        low_idx = space.space_to_ey_index(low)
        
        if self.direction is const.PlusZ and self.polarization is const.Y:
            TransparentEy = TransparentPlusZEy

        elif self.direction is const.MinusZ and self.polarization is const.Y:
            TransparentEy = TransparentMinusZEy
            
        elif self.direction is const.PlusX and self.polarization is const.Y:
            TransparentEy = TransparentPlusXEy
            
        elif self.direction is const.MinusX and self.polarization is const.Y:
            TransparentEy = TransparentMinusXEy

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    material_ey[i,j,k] = TransparentEy(material_ey[i,j,k], self.epsilon_r, self.amp, deepcopy(self.aux_fdtd))

    def set_pointwise_source_ez(self, material_ez, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        low_idx = space.space_to_ez_index(low)
        
        if self.direction is const.PlusY and self.polarization is const.Z:
            TransparentEz = TransparentPlusYEz
            
        elif self.direction is const.MinusY and self.polarization is const.Z:
            TransparentEz = TransparentMinusYEz
            
        elif self.direction is const.PlusX and self.polarization is const.Z:
            TransparentEz = TransparentPlusXEz
            
        elif self.direction is const.MinusX and self.polarization is const.Z:
            TransparentEz = TransparentMinusXEz

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    material_ez[i,j,k] = TransparentEz(material_ez[i,j,k], self.epsilon_r, self.amp, deepcopy(self.aux_fdtd))

    def set_pointwise_source_hx(self, material_hx, space):
        high = self.center + self.half_size
        low = self.center - self.half_size

        if self.direction is const.PlusY and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            amp = self.amp
            TansparentHx = TransparentPlusYHx
            
        elif self.direction is const.MinusY and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            amp = -self.amp
            TransparentHx = TransparentMinusYHx
            
        elif self.direction is const.PlusZ and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            amp = -self.amp
            TransparentHx = TransparentPlusZHx

        elif self.direction is const.MinusZ and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            amp = self.amp
            TransparentHx = TransparentMinusZHx

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    material_hx[i,j,k] = TransparentHx(material_hx[i,j,k], self.mu_r, amp, deepcopy(self.aux_fdtd))
        
    def set_pointwise_source_hy(self, material_hy, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.direction is const.PlusZ and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            
            amp = self.amp
            TransparentHy = TransparentPlusZHy

        elif self.direction is const.MinusZ and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            amp = -self.amp
            TransparentHy = TransparentMinusZHy

        elif self.direction is const.PlusX and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            amp = -self.amp
            TransparentHy = TransparentPlusXHy

        elif self.direction is const.MinusX and self.polarization is const.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            amp = -self.amp
            TransparentHy = TransparentMinusXHy

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    material_hy[i,j,k] = TransparentHy(material_hy[i,j,k], self.mu_r, amp, deepcopy(self.aux_fdtd))
                    
    def set_pointwise_source_hz(self, material_hz, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.direction is const.PlusY and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
    
            amp = -self.amp
            TransparentHz = TransparentPlusYHz

        elif self.direction is const.MinusY and self.polarization is const.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            amp = self.amp
            TransparentHz = TransparentMinusYHz
            
        elif self.direction is const.PlusX and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
    
            amp = self.amp
            TransparentHz = TransparentPlusXHz

        elif self.direction is const.MinusX and self.polarization is const.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            amp = -self.amp
            TransparentHz = TransparentMinusXHz

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    material_hz[i,j,k] = TransparentHz(material_hz[i,j,k], self.mu_r, amp, deepcopy(self.aux_fdtd))
                    