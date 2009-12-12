#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from numpy import *
from constants import *


class DipoleElectric(object):
    def __init__(self, pw_material, src_time=None, dt=None, amp=1):
        self.pw_material = pw_material
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.epsilon = pw_material.epsilon
        self.src_time = src_time
        self.amp = float(amp)
        
    def update(self, efield, hfield1, hfield2, space_diff1, space_diff2, dt, t):
        src_t = self.src_time.dipole(t)
        
        if src_t is None:
            self.pw_material.update(efield, hfield1, hfield2, space_diff1, space_diff2, dt, t)
        else:
            efield[self.i, self.j, self.k] = self.amp * src_t


class DipoleEx(DipoleElectric): pass
    
    
class DipoleEy(DipoleElectric): pass
    
    
class DipoleEz(DipoleElectric): pass
    
    
class DipoleMagnetic(object):
    def __init__(self, pw_material, src_time=None, dt=None, amp=1):
        self.pw_material = pw_material
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.mu = pw_material.mu
        self.src_time = src_time
        self.amp = float(amp)
        
    def update(self, hfield, efield1, efield2, space_diff1, space_diff2, dt, t):
        src_t = self.src_time.dipole(t)
        
        if src_t is None:
            self.pw_material.update(hfield, efield1, efield2, space_diff1, space_diff2, dt, t)
        else:
            hfield[self.i, self.j, self.k] = self.amp * src_t

        
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
        
        
class TransparentElectric(object):
    def __init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.epsilon = epsilon_r * epsilon0
        
        self.aux_fdtd = aux_fdtd
            
        self.face_list = [directional]
        self.amp = {directional: amp}
        
        samp_idx = aux_fdtd.space.spc_to_exact_hy_idx(samp_pnt)
        self.samp_idx0 = {directional: tuple(floor(samp_idx))}
        self.samp_idx1 = {directional: tuple(floor(samp_idx) + (0, 0, 1))}
        
        r1_value = samp_idx[2] - floor(samp_idx[2])
        self.r1 = {directional: r1_value}
        self.r0 = {directional: 1 - r1_value}

        self.corrective = {directional: corrective}
        
        if isinstance(pw_material, TransparentElectric):
            self.pw_material = pw_material.pw_material
            self += pw_material
        else:
            self.pw_material = pw_material
        
    def display_info(self, indent=0):
        print " " * indent, "array index:", self.i, self.j, self.k
        print " " * indent, "electric permittivity:", self.epsilon
        
        for face in self.face_list:
            print " " * indent, 
            print "face:", face,
            print "amp:", self.amp[face],
            print "sampling index 0:", self.samp_idx0[face],
            print "sampling index 1:", self.samp_idx1[face],
            print "sampling ratio 0:", self.r0[face],
            print "sampling ratio 1:", self.r1[face],
            print "corrective:", self.corrective[face]
            
    def __iadd__(self, rhs):
        self.face_list.extend(rhs.face_list)
        
        for face in rhs.face_list:
            self.amp[face] = rhs.amp[face]
            self.samp_idx0[face] = rhs.samp_idx0[face]
            self.samp_idx1[face] = rhs.samp_idx1[face]
            self.r1[face] = rhs.r1[face]
            self.r0[face] = rhs.r0[face]
            self.corrective[face] = rhs.corrective[face]
        
        return self
        
        
class TransparentEx(TransparentElectric):
    def __init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        TransparentElectric.__init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional)
        
        self._consist_cond = {MinusY: self._consistency_minus_y,
                              MinusZ: self._consistency_minus_z,
                              PlusY: self._consistency_plus_y,
                              PlusZ: self._consistency_plus_z}
        
    def display_info(self, indent=0):
        print " " * indent, "TransparentEx"
        TransparentElectric.display_info(self, indent)
         
    def update(self, ex, hz, hy, dy, dz, dt, t):
        self.pw_material.update(ex, hz, hy, dy, dz, dt, t)
        
        for face in self.face_list:
            self._consist_cond[face](ex, hz, hy, dy, dz, dt, face)
        
    def _consistency_minus_y(self, ex, hz, hy, dy, dz, dt, face):
        incidnet_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
                
        idx = self.i, self.j, self.k
        ex[idx] -= (dt / (self.corrective[face] * self.epsilon * dy) *
                    self.amp[face] * incidnet_hz)

    def _consistency_plus_y(self, ex, hz, hy, dy, dz, dt, face):
        incidnet_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ex[idx] += (dt / (self.corrective[face] * self.epsilon * dy) *
                    self.amp[face] * incidnet_hz)
    
    def _consistency_minus_z(self, ex, hz, hy, dy, dz, dt, face):
        incidnet_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ex[idx] += (dt / (self.corrective[face] * self.epsilon * dz) *
                    self.amp[face] * incidnet_hy)

    def _consistency_plus_z(self, ex, hz, hy, dy, dz, dt, face):
        incidnet_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ex[idx] -= (dt / (self.corrective[face] * self.epsilon * dz) *
                    self.amp[face] * incidnet_hy)
        
        
class TransparentEy(TransparentElectric):
    def __init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        TransparentElectric.__init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional)
        
        self._consist_cond = {MinusZ: self._consistency_minus_z,
                              MinusX: self._consistency_minus_x,
                              PlusZ: self._consistency_plus_z,
                              PlusX: self._consistency_plus_x}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentEy"
        TransparentElectric.display_info(self, indent)
        
    def update(self, ey, hx, hz, dz, dx, dt, t):
        self.pw_material.update(ey, hx, hz, dz, dx, dt, t)
        
        for face in self.face_list:
            self._consist_cond[face](ey, hx, hz, dz, dx, dt, face)
        
    def _consistency_minus_z(self, ey, hx, hz, dz, dx, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] -= (dt / (self.corrective[face] * self.epsilon * dz) *
                    self.amp[face] * incident_hx)

    def _consistency_minus_x(self, ey, hx, hz, dz, dx, dt, face):
        incident_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] += (dt / (self.corrective[face] * self.epsilon * dx) *
                    self.amp[face] * incident_hz)
    
    def _consistency_plus_z(self, ey, hx, hz, dz, dx, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] += (dt / (self.corrective[face] * self.epsilon * dz) *
                    self.amp[face] * incident_hx)

    def _consistency_plus_x(self, ey, hx, hz, dz, dx, dt, face):
        incident_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] -= (dt / (self.corrective[face] * self.epsilon * dx) *
                    self.amp[face] * incident_hz)
        

class TransparentEz(TransparentElectric):
    def __init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        TransparentElectric.__init__(self, pw_material, epsilon_r, amp, aux_fdtd, samp_pnt, corrective, directional)
        
        self._consist_cond = {MinusX: self._consistency_minus_x,
                              MinusY: self._consistency_minus_y,
                              PlusX: self._consistency_plus_x,
                              PlusY: self._consistency_plus_y}
        
    def display_info(self, indent=0):
        print " " * indent, "TransparentEz"
        TransparentElectric.display_info(self, indent)
                
    def update(self, ez, hy, hx, dx, dy, dt, t):
        self.pw_material.update(ez, hy, hx, dx, dy, dt, t)
        
        for face in self.face_list:
            self._consist_cond[face](ez, hy, hx, dx, dy, dt, face)
        
    def _consistency_minus_x(self, ez, hy, hx, dx, dy, dt, face):
        incident_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] -= (dt / (self.corrective[face] * self.epsilon * dx) *
                    self.amp[face] * incident_hy)

    def _consistency_minus_y(self, ez, hy, hx, dx, dy, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] += (dt / (self.corrective[face] * self.epsilon * dy) *
                    self.amp[face] * incident_hx)
    
    def _consistency_plus_x(self, ez, hy, hx, dx, dy, dt, face):
        incident_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] += (dt / (self.corrective[face] * self.epsilon * dx) *
                    self.amp[face] * incident_hy)

    def _consistency_plus_y(self, ez, hy, hx, dx, dy, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] -= (dt / (self.corrective[face] * self.epsilon * dy) *
                    self.amp[face] * incident_hx)


class TransparentMagnetic(object):
    def __init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.mu = mu_r * epsilon0
        
        self.aux_fdtd = aux_fdtd
        
        self.face_list = [directional]
        self.amp = {directional: amp}
        
        samp_idx = aux_fdtd.space.spc_to_exact_ex_idx(samp_pnt)
        self.samp_idx0 = {directional: tuple(floor(samp_idx))}
        self.samp_idx1 = {directional: tuple(floor(samp_idx) + (0, 0, 1))}
        
        r1_value = samp_idx[2] - floor(samp_idx[2])
        self.r1 = {directional: r1_value}
        self.r0 = {directional: 1 - r1_value}

        self.corrective = {directional: corrective}
        
        if isinstance(pw_material, TransparentElectric):
            self.pw_material = pw_material.pw_material
            self += pw_material
        else:
            self.pw_material = pw_material
        
    def display_info(self, indent=0):
        print " " * indent, "array index:", self.i, self.j, self.k
        print " " * indent, "magnetic permeability:", self.mu 
        
        for face in self.face_list:
            print " " * indent, 
            print "face:", face,
            print "amp:", self.amp[face],
            print "sampling index 0:", self.samp_idx0[face],
            print "sampling index 1:", self.samp_idx1[face],
            print "sampling ratio 0:", self.r0[face],
            print "sampling ratio 1:", self.r1[face],
            print "corrective:", self.corrective[face]
            
    def __iadd__(self, rhs):
        self.face_list.extend(rhs.face_list)
        
        for face in rhs.face_list:
            self.amp[face] = rhs.amp[face]
            self.samp_idx0[face] = rhs.samp_idx0[face]
            self.samp_idx1[face] = rhs.samp_idx1[face]
            self.r1[face] = rhs.r1[face]
            self.r0[face] = rhs.r0[face]
            self.corrective[face] = rhs.corrective[face]
            
        return self
   

class TransparentHx(TransparentMagnetic):
    def __init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        TransparentMagnetic.__init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional)
        
        self._consist_cond = {MinusY: self._consistency_minus_y,
                              MinusZ: self._consistency_minus_z,
                              PlusY: self._consistency_plus_y,
                              PlusZ: self._consistency_plus_z}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentHx"
        TransparentMagnetic.display_info(self, indent)
        
    def update(self, hx, ez, ey, dy, dz, dt, t):
        self.pw_material.update(hx, ez, ey, dy, dz, dt, t)
        
        for face in self.face_list:
            self._consist_cond[face](hx, ez, ey, dy, dz, dt, face)
            
    def _consistency_minus_y(self, hx, ez, ey, dt, dy, dz, face):
        incidnet_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hx[idx] += (dt / (self.corrective[face] * self.mu * dy) * 
                    self.amp[face] * incidnet_ez)

    def _consistency_plus_y(self, hx, ez, ey, dy, dz, dt, face):
        incidnet_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hx[idx] -= (dt / (self.corrective[face] * self.mu * dy) * 
                    self.amp[face] * incidnet_ez)
    
    def _consistency_minus_z(self, hx, ez, ey, dy, dz, dt, face):
        incidnet_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hx[idx] -= (dt / (self.corrective[face] * self.mu * dz) * 
                    self.amp[face] * incidnet_ey)

    def _consistency_plus_z(self, hx, ez, ey, dy, dz, dt, face):
        incidnet_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hx[idx] += (dt / (self.corrective[face] * self.mu * dz) * 
                    self.amp[face] * incidnet_ey)
        

class TransparentHy(TransparentMagnetic):
    def __init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        TransparentMagnetic.__init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional)
        
        self._consist_cond = {MinusZ: self._consistency_minus_z,
                              MinusX: self._consistency_minus_x,
                              PlusZ: self._consistency_plus_z,
                              PlusX: self._consistency_plus_x}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentHy"
        TransparentMagnetic.display_info(self, indent)
        
    def update(self, hy, ex, ez, dz, dx, dt, t):
        self.pw_material.update(hy, ex, ez, dz, dx, dt, t)
        
        for face in self.face_list:
            self._consist_cond[face](hy, ex, ez, dz, dx, dt, face)
        
        
    def _consistency_minus_z(self, hy, ex, ez, dz, dx, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] += (dt / (self.corrective[face] * self.mu * dz) * 
                    self.amp[face] * incident_ex)

    def _consistency_minus_x(self, hy, ex, ez, dz, dx, dt, face):
        incident_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] -= (dt / (self.corrective[face] * self.mu * dx) * 
                    self.amp[face] * incident_ez)
    
    def _consistency_plus_z(self, hy, ex, ez, dz, dx, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] -= (dt / (self.corrective[face] * self.mu * dz) * 
                    self.amp[face] * incident_ex)

    def _consistency_plus_x(self, hy, ex, ez, dz, dx, dt, face):
        incident_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] += (dt / (self.corrective[face] * self.mu * dx) * 
                    self.amp[face] * incident_ez)
        
        
class TransparentHz(TransparentMagnetic):
    def __init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional):
        TransparentMagnetic.__init__(self, pw_material, mu_r, amp, aux_fdtd, samp_pnt, corrective, directional)
        
        self._consist_cond = {MinusX: self._consistency_minus_x,
                              MinusY: self._consistency_minus_y,
                              PlusX: self._consistency_plus_x,
                              PlusY: self._consistency_plus_y}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentHz"
        TransparentMagnetic.display_info(self, indent)
                
    def update(self, hz, ey, ex, dx, dy, dt, t):
        self.pw_material.update(hz, ey, ex, dx, dy, dt, t)
        
        for face in self.face_list:
            self._consist_cond[face](hz, ey, ex, dx, dy, dt, face)
        
    def _consistency_minus_x(self, hz, ey, ex, dx, dy, dt, face):
        incident_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
    
        idx = self.i, self.j, self.k
        hz[idx] += (dt / (self.corrective[face] * self.mu * dx) * 
                    self.amp[face] * incident_ey)

    def _consistency_minus_y(self, hz, ey, ex, dx, dy, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hz[idx] -= (dt / (self.corrective[face] * self.mu * dy) * 
                    self.amp[face] * incident_ex)
    
    def _consistency_plus_(self, hz, ey, ex, dx, dy, dt, face):
        incident_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hz[idx] -= (dt / (self.corrective[face] * self.mu * dx) * 
                    self.amp[face] * incident_ey)

    def _consistency_plus_y(self, hz, ey, ex, dx, dy, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hz[idx] += (dt / (self.corrective[face] * self.mu * dy) * 
                    self.amp[face] * incident_ex)
        