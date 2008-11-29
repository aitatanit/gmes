#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
import constants as const


class DipoleElectric(object):
    def __init__(self, pw_material, src_time=None, dt=None, amp=1):
        self.pw_material = pw_material
        self.idx = pw_material.i, pw_material.j, pw_material.k
        self.epsilon = pw_material.epsilon
        self.src_time = src_time
        self.amp = float(amp)
        self.t = 0.0
        self.n = 0
        
    def update(self, efield, hfield1, hfield2, dt, space_diff1, space_diff2):
        src_t = self.src_time.dipole(self.t)
        
        if src_t is None:
            self.pw_material.update(efield, hfield1, hfield2, dt, space_diff1, space_diff2)
        else:
            efield[self.idx] = self.amp * src_t

        self.n += 1
        self.t = self.n * dt


class DipoleEx(DipoleElectric): pass
    
    
class DipoleEy(DipoleElectric): pass
    
    
class DipoleEz(DipoleElectric): pass
    
    
class DipoleMagnetic(object):
    def __init__(self, pw_material, src_time=None, dt=None, amp=1):
        self.pw_material = pw_material
        self.idx = pw_material.i, pw_material.j, pw_material.k
        self.mu = pw_material.mu
        self.src_time = src_time
        self.amp = float(amp)
        self.t = .5 * dt
        self.n = 0.5

    def update(self, hfield, efield1, efield2, dt, space_diff1, space_diff2):
        src_t = self.src_time.dipole(self.t)
        
        if src_t is None:
            self.pw_material.update(hfield, efield1, efield2, dt, space_diff1, space_diff2)
        else:
            hfield[self.idx] = self.amp * self.src_t

        self.n += 1
        self.t = self.n * dt
        
        
class DipoleHx(DipoleMagnetic): pass


class DipoleHy(DipoleMagnetic): pass


class DipoleHz(DipoleMagnetic): pass


class _SrcTime(object):
    """Time-dependent part of a source.
    """


class _Continuous(_SrcTime):
    """_Continuous (CW) source with (optional) slow turn-on and/or turn-off.
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
        
        
class TransparentElectric(object):
    def __init__(self, pw_material, epsilon_r, amp, aux_fdtd):
        self.idx = pw_material.i, pw_material.j, pw_material.k
        self.pw_material = pw_material
        self.epsilon = epsilon_r * const.epsilon0
        self.amp = amp
        self.aux_fdtd = aux_fdtd
        

class TransparentPlusYEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        ex[self.idx] -= dt / (self.epsilon * dy) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        
        
class TransparentPlusYEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        ez[self.idx] += dt / (self.epsilon * dy) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        

class TransparentMinusYEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        ex[self.idx] += dt / (self.epsilon * dy) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        
        
class TransparentMinusYEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        ez[self.idx] -= dt / (self.epsilon * dy) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        

class TransparentPlusZEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        ex[self.idx] += dt / (self.epsilon * dz) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        
        
class TransparentPlusZEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        ey[self.idx] -= dt / (self.epsilon * dz) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        

class TransparentMinusZEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        ex[self.idx] -= dt / (self.epsilon * dz) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        
        
class TransparentMinusZEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        ey[self.idx] += dt / (self.epsilon * dz) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()

                                                
class TransparentPlusXEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        ey[self.idx] += dt / (self.epsilon * dx) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        

class TransparentPlusXEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        ez[self.idx] -= dt / (self.epsilon * dx) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        
        
class TransparentMinusXEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        ey[self.idx] -= dt / (self.epsilon * dx) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        

class TransparentMinusXEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        ez[self.idx] += dt / (self.epsilon * dx) * self.amp * self.aux_fdtd.hy[1,0,12]
        self.aux_fdtd.step()
        
        
class TransparentMagnetic(object):
    def __init__(self, pw_material, mu_r, amp, aux_fdtd):
        self.idx = pw_material.i, pw_material.j, pw_material.k
        self.pw_material = pw_material
        self.mu = mu_r * const.mu0
        self.amp = amp
        self.aux_fdtd = aux_fdtd


class TransparentPlusYHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        hz[self.idx] -= dt / (self.mu * dy) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        
class TransparentPlusYHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        hx[self.idx] += dt / (self.mu * dy) * self.amp * self.aux_fdtd.ex[0,0,12]
        

class TransparentMinusYHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        hz[self.idx] += dt / (self.mu * dy) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        
class TransparentMinusYHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        hx[self.idx] -= dt / (self.mu * dy) * self.amp * self.aux_fdtd.ex[0,0,12]
        

class TransparentPlusZHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        hy[self.idx] += dt / (self.mu * dz) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        
class TransparentPlusZHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        hx[self.idx] -= dt / (self.mu * dz) * self.amp * self.aux_fdtd.ex[0,0,12]
        

class TransparentMinusZHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        hy[self.idx] -= dt / (self.mu * dz) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        
class TransparentMinusZHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        hx[self.idx] += dt / (self.mu * dz) * self.amp * self.aux_fdtd.ex[0,0,12]
        
                                
class TransparentPlusXHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        hz[self.idx] += dt / (self.mu * dx) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        
class TransparentPlusXHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        hy[self.idx] -= dt / (self.mu * dx) * self.amp * self.aux_fdtd.ex[0,0,12]
        
                
class TransparentMinusXHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        hz[self.idx] -= dt / (self.mu * dx) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        
class TransparentMinusXHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        hy[self.idx] += dt / (self.mu * dx) * self.amp * self.aux_fdtd.ex[0,0,12]
        
        