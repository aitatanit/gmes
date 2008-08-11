#!/usr/bin/env python

from copy import deepcopy
from math import cos, sin, exp, pi
from numpy.core import inf, array

from pointwise_source import *
import constants

from geometric import Cartesian, DefaultMaterial, Boundary
from fdtd import TEMzFDTD
from material import Dielectric, UPML


class SrcTime:
    """Time-dependent part of a source.
    """

    
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
        self.pos = array(pos, 'd')
        self.comp = component
        self.src_time = src_time
        self.amp = float(amp)
    
    def init(self, geom_tree, space):
        pass
    
    def set_pointwise_source_ex(self, material_ex, space):
        if self.comp is constants.Ex:
            idx = space.space_to_ex_index(self.pos)
            material_ex[idx] = DipoleEx(material_ex[idx], self.src_time, space.dt, self.amp)
            
    def set_pointwise_source_ey(self, material_ey, space):
        if self.comp is constants.Ey:
            idx = space.space_to_ey_index(self.pos)
            material_ey[idx] = DipoleEy(material_ey[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_ez(self, material_ez, space):
        if self.comp is constants.Ez:
            idx = space.space_to_ez_index(self.pos)
            material_ez[idx] = DipoleEz(material_ez[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_hx(self, material_hx, space):
        if self.comp is constants.Hx:
            idx = space.space_to_hx_index(self.pos)
            material_hx[idx] = DipoleHx(material_hx[idx], self.src_time, space.dt, self.amp)
   
    def set_pointwise_source_hy(self, material_hy, space):
        if self.comp is constants.Hy:
            idx = space.space_to_hy_index(self.pos)
            material_hy[idx] = DipoleHy(material_hy[idx], self.src_time, space.dt, self.amp)
            
    def set_pointwise_source_hz(self, material_hz, space):
        if self.comp is constants.Hz:
            idx = space.space_to_hz_index(self.pos)
            material_hz[idx] = DipoleHz(material_hz[idx], self.src_time, space.dt, self.amp)


class TotalFieldScatteredField:
    def __init__(self, theta, phi, psi, low, high):
        self.theta = float(theta)
        self.phi = float(theta)
        self.psi = float(psi)
        self.low = array(low, 'd')
        self.high = array(high, 'd')

    def generate_aux_fdtd(self, space):
        pass
    
    
class Transparent:
    def __init__(self, direction, center, size, freq, polarization, amp=1,
                 ex_mode_file=None, ey_mode_file=None, ez_mode_file=None,
                 hx_mode_file=None, hy_mode_file=None, hz_mode_file=None):
        # launching plane parameters
        self.direction = direction
        self.center = array(center, 'd')
        self.size = array(size, 'd')
        self.half_size = .5 * self.size
        self.freq = float(freq)
        self.polarization = polarization
        
        # maximum amplitude of stimulus
        self.amp = float(amp)
        
    def init(self, geom_tree, space):
        go = geom_tree.object_of_point(self.center)
        self.epsilon_r = go.material.epsilon_r
        self.mu_r = go.material.mu_r
        
    def get_aux_fdtd(self, space):
        # two 10 meshes for the ABC,
        # 1 ex point and 2 hy points for the free space
        aux_size = array((0 , 0, 21), 'd') / space.res
        aux_space = Cartesian(size=aux_size, resolution=space.res)
        aux_geom_list = (DefaultMaterial(material=Dielectric(self.epsilon_r, self.mu_r)),
                         Boundary(material=UPML(), thickness=10 / space.res[2], size=aux_size))
        aux_src_list = (Dipole(src_time=Continuous(freq=self.freq), component=constants.Ex, pos=(0,0,0)),)
        aux_fdtd = TEMzFDTD(aux_space, aux_geom_list, aux_src_list, verbose=False)
        
        return aux_fdtd
        
    def set_pointwise_source_ex(self, material_ex, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
        low_idx = space.space_to_ex_index(low)
        
        if self.direction is constants.PlusY and self.polarization is constants.X:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ex[(i,j,k)] = TransparentPlusYEx(material_ex[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusY and self.polarization is constants.X:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ex[(i,j,k)] = TransparentMinusYEx(material_ex[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
                        
        elif self.direction is constants.PlusZ and self.polarization is constants.X:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ex[(i,j,k)] = TransparentPlusZEx(material_ex[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusZ and self.polarization is constants.X:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ex[(i,j,k)] = TransparentMinusZEx(material_ex[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))

    def set_pointwise_source_ey(self, material_ey, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
        low_idx = space.space_to_ey_index(low)
        
        if self.direction is constants.PlusZ and self.polarization is constants.Y:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ey[(i,j,k)] = TransparentPlusZEy(material_ey[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusZ and self.polarization is constants.Y:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ey[(i,j,k)] = TransparentMinusZEy(material_ey[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
                        
        elif self.direction is constants.PlusX and self.polarization is constants.Y:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ey[(i,j,k)] = TransparentPlusXEy(material_ey[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusX and self.polarization is constants.Y:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ey[(i,j,k)] = TransparentMinusXEy(material_ey[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))

    def set_pointwise_source_ez(self, material_ez, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
        low_idx = space.space_to_ez_index(low)
        
        if self.direction is constants.PlusY and self.polarization is constants.Z:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ez[(i,j,k)] = TransparentPlusYEz(material_ez[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusY and self.polarization is constants.Z:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ez[(i,j,k)] = TransparentMinusYEz(material_ez[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
                        
        elif self.direction is constants.PlusX and self.polarization is constants.Z:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ez[(i,j,k)] = TransparentPlusXEz(material_ez[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusX and self.polarization is constants.Z:
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_ez[(i,j,k)] = TransparentMinusXEz(material_ez[(i,j,k)], self.epsilon_r, self.amp, self.get_aux_fdtd(space))

    def set_pointwise_source_hx(self, material_hx, space):
        high = self.center + self.half_size
        low = self.center - self.half_size

        if self.direction is constants.PlusY and self.polarization is constants.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hx[(i,j,k)] = TransparentPlusYHx(material_hy[(i,j,k)], self.mu_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusY and self.polarization is constants.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hx[(i,j,k)] = TransparentMinusYHx(material_hy[(i,j,k)], self.mu_r, -self.amp, self.get_aux_fdtd(space))
                                
        elif self.direction is constants.PlusZ and self.polarization is constants.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hx[(i,j,k)] = TransparentPlusZHx(material_hx[(i,j,k)], self.mu_r, -self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusZ and self.polarization is constants.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2] + 1)
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hx[(i,j,k)] = TransparentMinusZHx(material_hx[(i,j,k)], self.mu_r, self.amp, self.get_aux_fdtd(space))

    def set_pointwise_source_hy(self, material_hy, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.direction is constants.PlusZ and self.polarization is constants.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hy[(i,j,k)] = TransparentPlusZHy(material_hy[(i,j,k)], self.mu_r, self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusZ and self.polarization is constants.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hy[(i,j,k)] = TransparentMinusZHy(material_hy[(i,j,k)], self.mu_r, -self.amp, self.get_aux_fdtd(space))

        elif self.direction is constants.PlusX and self.polarization is constants.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ez_index(high))
            low_idx = space.space_to_ez_index(low)
            
            high_idx = (high_idx[0], high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0], low_idx[1], low_idx[2] + 1)
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hy[(i,j,k)] = TransparentPlusXHy(material_hy[(i,j,k)], self.mu_r, -self.amp, self.get_aux_fdtd(space))
        
        elif self.direction is constants.MinusX and self.polarization is constants.Z:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2] + 1)
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2] + 1)
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hy[(i,j,k)] = TransparentMinusXHy(material_hy[(i,j,k)], self.mu_r, self.amp, self.get_aux_fdtd(space))
                                                
    def set_pointwise_source_hz(self, material_hz, space):
        high = self.center + self.half_size
        low = self.center - self.half_size
        
        if self.direction is constants.PlusY and self.polarization is constants.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1], high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1], low_idx[2])
    
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hz[(i,j,k)] = TransparentPlusYHz(material_hz[(i,j,k)], self.mu_r, -self.amp, self.get_aux_fdtd(space))
            
        elif self.direction is constants.MinusY and self.polarization is constants.X:
            high_idx = map(lambda x: x + 1, space.space_to_ex_index(high))
            low_idx = space.space_to_ex_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hz[(i,j,k)] = TransparentMinusYHz(material_hz[(i,j,k)], self.mu_r, self.amp, self.get_aux_fdtd(space))
                                                 
        elif self.direction is constants.PlusX and self.polarization is constants.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0], high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0], low_idx[1] + 1, low_idx[2])
    
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hz[(i,j,k)] = TransparentPlusXHz(material_hz[(i,j,k)], self.mu_r, self.amp, self.get_aux_fdtd(space))
            
        elif self.direction is constants.MinusX and self.polarization is constants.Y:
            high_idx = map(lambda x: x + 1, space.space_to_ey_index(high))
            low_idx = space.space_to_ey_index(low)
            
            high_idx = (high_idx[0] + 1, high_idx[1] + 1, high_idx[2])
            low_idx = (low_idx[0] + 1, low_idx[1] + 1, low_idx[2])
            
            for i in xrange(low_idx[0], high_idx[0]):
                for j in xrange(low_idx[1], high_idx[1]):
                    for k in xrange(low_idx[2], high_idx[2]):
                        material_hz[(i,j,k)] = TransparentMinusXHz(material_hz[(i,j,k)], self.mu_r, -self.amp, self.get_aux_fdtd(space))
                        
                        