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
    def __init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        self.idx = pw_material.i, pw_material.j, pw_material.k
        self.pw_material = pw_material
        self.epsilon = epsilon_r * const.epsilon0
        self.amp = amp
        self.aux_fdtd = aux_fdtd
        
        samp_idx = aux_fdtd.space.spc_to_exact_hy_idx(samp_pnt)
        self.samp_idx0 = tuple(floor(samp_idx))
        self.samp_idx1 = tuple(floor(samp_idx) + (0, 0, 1))
        
        self.r1 = samp_idx[2] - floor(samp_idx[2])
        self.r0 = 1 - self.r1

        self.corrective = corrective
        
class TransparentMinusYEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ex[self.idx] -= dt / (self.corrective * self.epsilon * dy) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        
        
class TransparentMinusYEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ez[self.idx] += dt / (self.corrective * self.epsilon * dy) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        

class TransparentPlusYEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ex[self.idx] += dt / (self.corrective * self.epsilon * dy) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        
        
class TransparentPlusYEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ez[self.idx] -= dt / (self.corrective * self.epsilon * dy) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        

class TransparentMinusZEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ex[self.idx] += dt / (self.corrective * self.epsilon * dz) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        
        
class TransparentMinusZEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ey[self.idx] -= dt / (self.corrective * self.epsilon * dz) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        

class TransparentPlusZEx(TransparentElectric):
    def __init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ex, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ex, hz, hy, dt, dy, dz):
        self.pw_material.update(ex, hz, hy, dt, dy, dz)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ex[self.idx] -= dt / (self.corrective * self.epsilon * dz) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        
        
class TransparentPlusZEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ey[self.idx] += dt / (self.corrective * self.epsilon * dz) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()

                                                
class TransparentMinusXEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ey[self.idx] += dt / (self.corrective * self.epsilon * dx) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        

class TransparentMinusXEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)

        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
    
        ez[self.idx] -= dt / (self.corrective * self.epsilon * dx) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()

        
class TransparentPlusXEy(TransparentElectric):
    def __init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ey, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ey, hx, hz, dt, dz, dx):
        self.pw_material.update(ey, hx, hz, dt, dz, dx)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ey[self.idx] -= dt / (self.corrective * self.epsilon * dx) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        

class TransparentPlusXEz(TransparentElectric):
    def __init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentElectric.__init__(self, pw_material_ez, epsilon_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, ez, hy, hx, dt, dx, dy):
        self.pw_material.update(ez, hy, hx, dt, dx, dy)
        
        aux_fdtd_hy = self.r0 * self.aux_fdtd.hy[self.samp_idx0] + self.r1 * self.aux_fdtd.hy[self.samp_idx1]
        
        ez[self.idx] += dt / (self.corrective * self.epsilon * dx) * self.amp * aux_fdtd_hy
        self.aux_fdtd.step()
        
        
class TransparentMagnetic(object):
    def __init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        self.idx = pw_material.i, pw_material.j, pw_material.k
        self.pw_material = pw_material
        self.mu = mu_r * const.mu0
        self.amp = amp
        self.aux_fdtd = aux_fdtd
        self.samp_idx = aux_fdtd.space.space_to_ex_index(samp_pnt)

        samp_idx = aux_fdtd.space.spc_to_exact_ex_idx(samp_pnt)
        self.samp_idx0 = tuple(floor(samp_idx))
        self.samp_idx1 = tuple(floor(samp_idx) + (0, 0, 1))
        
        self.r1 = samp_idx[2] - floor(samp_idx[2])
        self.r0 = 1 - self.r1
        
        self.corrective = corrective

class TransparentMinusYHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hz[self.idx] -= dt / (self.corrective * self.mu * dy) * self.amp * aux_fdtd_ex
        
        
class TransparentMinusYHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hx[self.idx] += dt / (self.corrective * self.mu * dy) * self.amp * aux_fdtd_ex
        

class TransparentPlusYHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hz[self.idx] += dt / (self.corrective * self.mu * dy) * self.amp * aux_fdtd_ex
        
        
class TransparentPlusYHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hx[self.idx] -= dt / (self.corrective * self.mu * dy) * self.amp * aux_fdtd_ex
        
        
class TransparentMinusZHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hy[self.idx] += dt / (self.corrective * self.mu * dz) * self.amp * aux_fdtd_ex
        
        
class TransparentMinusZHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hx[self.idx] -= dt / (self.corrective * self.mu * dz) * self.amp * aux_fdtd_ex
        

class TransparentPlusZHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hy[self.idx] -= dt / (self.corrective * self.mu * dz) * self.amp * aux_fdtd_ex
        
        
class TransparentPlusZHx(TransparentMagnetic):
    def __init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hx, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hx, ez, ey, dt, dy, dz):
        self.pw_material.update(hx, ez, ey, dt, dy, dz)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hx[self.idx] += dt / (self.corrective * self.mu * dz) * self.amp * aux_fdtd_ex
        
                                
class TransparentMinusXHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hz[self.idx] += dt / (self.corrective * self.mu * dx) * self.amp * aux_fdtd_ex
        
        
class TransparentMinusXHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hy[self.idx] -= dt / (self.corrective * self.mu * dx) * self.amp * aux_fdtd_ex
        
                
class TransparentPlusXHz(TransparentMagnetic):
    def __init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hz, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hz, ey, ex, dt, dx, dy):
        self.pw_material.update(hz, ey, ex, dt, dx, dy)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hz[self.idx] -= dt / (self.corrective * self.mu * dx) * self.amp * aux_fdtd_ex
        
        
class TransparentPlusXHy(TransparentMagnetic):
    def __init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective):
        TransparentMagnetic.__init__(self, pw_material_hy, mu_r, amp, aux_fdtd, samp_pnt, corrective)
        
    def update(self, hy, ex, ez, dt, dz, dx):
        self.pw_material.update(hy, ex, ez, dt, dz, dx)
        self.aux_fdtd.step()
        
        aux_fdtd_ex = self.r0 * self.aux_fdtd.ex[self.samp_idx0] + self.r1 * self.aux_fdtd.ex[self.samp_idx1]
        
        hy[self.idx] += dt / (self.corrective * self.mu * dx) * self.amp * aux_fdtd_ex
        
        