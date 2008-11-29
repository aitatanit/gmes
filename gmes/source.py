#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from copy import deepcopy
from numpy import *

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
    
    Attributes:
        freq
        phase -- The initial phase advance of the source.
        start -- The starting time for the source.
        end -- The end time for the source.
        width -- The temporal width of the smoothing.
    
    """
    def __init__(self, frequency, phase=0, start=0, end=inf, width=None):
        """
        
        Arguments:
            frequency
            phase -- Specify the initial phase advance of the source.
                Default is 0.
            start -- The start time for the source. Default is 0 
                (turn on at t = 0). 
            end -- The end time for the source; default is inf 
                (never turn off). 
            width -- The temporal width of the smoothing. 
                Default is None (automatically determines).
        
        """
        self.freq = float(frequency)
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
            return None
        
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
    
    Attributes:
        freq -- -- The center frequency in units of c.
        fwidth -- The frequency width.
        peak
        cutoff -- How many widths the source decays for before we turn it
            off and returns None - this applies for both turn-on and 
            turn-off of the pulse. 
    
    """
    def __init__(self, frequency, fwidth):
        """
        
        Arguments:
            frequency -- The center frequency in units of c.
                No default. 
            fwidth -- The frequency width. No default.
            
            """
        self.freq = float(frequency)
        self.fwidth = float(fwidth)
        width = 1 / self.fwidth
        s = 5
        self.peak = width * s
        self.cutoff = 2 * width * s
        
        while exp(-.5 * self.cutoff * self.cutoff / width / width) == 0:
            self.cutoff *= .9
        
        if self.peak - self.cutoff < 0:
            self.peak = self.cutoff
    
    def dipole(self, time):
        tt = time - self.peak
        if (abs(tt) > self.cutoff): return None

        return exp(-.5 * (tt * self.fwidth)**2) * cos(2 * pi * self.freq * time)
    
    def display_info(self, indent=0):
        print " " * indent, "bandpass source"
        print " " * indent,
        print "center frequency:", self.freq,
        print "bandwidth:", self.fwidth,
        print "peak time:", self.peak
        print "cutoff:", self.cutoff
        
        
class Dipole(object):
    """A hard source acting on a single point.
    
    This source acts as a hard source just when the src_time does not 
    return None.
    
    Attributes:
        pos -- The location of the center of the source in the space.
        comp -- The direction and type of the field component.
        src_time -- Specify the time-dependence of the source. 
        amp -- An overall amplitude multiplying the the source.
    
    """
    def __init__(self, position, component, src_time, amplitude=1):
        """ 
        
        Arguments:
            position -- The location of the source in the space. 
                No default.
            component -- Specify the direction and type of the 
                component: e.g. constants.{Ex, Ey, Ez} for an electric 
                field, and constants.{Hx, Hy, Hz} for a magnetic field. 
                Note that field pointing in an arbitrary direction are 
                specified simply as multiple sources with the 
                appropriate amplitudes for each component. No default.
            src_time -- Specify the time-dependence of the source. 
                No default.
            amplitude -- An overall amplitude multiplying the the 
                source. Default is 1.
        
        """
        
        self.pos = array(position, float)
        self.comp = component
        self.src_time = src_time
        self.amp = float(amplitude)
    
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
    
    
class Transparent(object):
    """
    
    Attributes:
        direction -- The propagation direction.
        center -- The location of the center of the source in the 
            space.
        size -- The size (in the form of length 3 tuple) of the source.
        half_size
        freq
        polarization -- The electric field direction.
    
    """
    def __init__(self, direction, center, size, frequency, polarization, amplitude=1):#,
#                 ex_mode_file=None, ey_mode_file=None, ez_mode_file=None,
#                 hx_mode_file=None, hy_mode_file=None, hz_mode_file=None):
        """
        
        Arguments:
            direction -- Specify the propagation direction: e.g. 
                constants.{PlusX, MinusX, PlusY, MinusY, PlusZ, 
                MinusZ}. No default.
            center -- The location of the center of the source in the 
                space. No default.
            size -- The size (in the form of length 3 tuple) of the 
                source. Note that the component in the propagation 
                direction is ignored. No default value.
            frequency -- No default value.
            polarization -- Specify the electric field direction: e.g.
                constants.{X, Y, Z}. No default value.
            amplitude -- An overall amplitude multiplying the the 
                source. Default is 1.
#            ex_mode_file -- Default is None.
#            ey_mode_file -- Default is None.
#            ez_mode_file -- Default is None.
#            hx_mode_file -- Default is None.
#            hy_mode_file -- Default is None.
#            hz_mode_file -- Default is None.
            
        """
        # launching plane parameters
        self.direction = direction
        self.center = array(center, float)
        self.size = array(size, float)
        self.half_size = .5 * self.size
        self.freq = float(frequency)
        self.polarization = polarization
        
        # maximum amplitude of stimulus
        self.amp = float(amplitude)
        
    def init(self, geom_tree, space):
        mat_objs = geom_tree.material_of_point(self.center)
        self.epsilon_r = mat_objs[0].epsilon_r
        self.mu_r = mat_objs[0].mu_r
        self.aux_fdtd = self.get_aux_fdtd(space)
        
    def get_aux_fdtd(self, space):
        # two 10 meshes for the ABC,
        # 1 Ex point and 2 Hy points for the free space
        aux_size = array((0 , 0, 21), float) / space.res[self.direction.tag]
        aux_space = Cartesian(size=aux_size, resolution=space.res[self.direction.tag], dt=space.dt, parallel=False)
        aux_geom_list = (DefaultMaterial(material=Dielectric(self.epsilon_r, self.mu_r)),
                         Boundary(material=CPML(self.epsilon_r, self.mu_r), thickness=10. / aux_space.res[2], size=aux_size))
        aux_src_list = (Dipole(src_time=Continuous(frequency=self.freq), component=const.Ex, position=(0,0,0)),)
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
                    if in_range((i,j,k), material_ex, const.Ex):
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
                    if in_range((i,j,k), material_ey, const.Ey):
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
                    if in_range((i,j,k), material_ez, const.Ez):
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
            TransparentHx = TransparentPlusYHx

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
                    if in_range((i,j,k), material_hx, const.Hx):
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
            
            amp = self.amp
            TransparentHy = TransparentMinusXHy

        else:
            return None
        
        for i in xrange(low_idx[0], high_idx[0]):
            for j in xrange(low_idx[1], high_idx[1]):
                for k in xrange(low_idx[2], high_idx[2]):
                    if in_range((i,j,k), material_hy, const.Hy):
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
                    if in_range((i,j,k), material_hz, const.Hz):
                        material_hz[i,j,k] = TransparentHz(material_hz[i,j,k], self.mu_r, amp, deepcopy(self.aux_fdtd))
                        
