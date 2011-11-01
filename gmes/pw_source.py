#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass

from math import cos, sin, pi, floor

import numpy as np
from numpy import inf

# GMES modules
import constant as const


class PwSourceParam(object):
    pass


class PwSource(object):
    def __init__(self):
        self._param = {}

    def attach(self, idx, parameter):
        key = tuple(idx)
        if self._param.has_key(key):
            stderr.write('Overwriting the existing index.\n')
        self._param[tuple(idx)] = parameter

    def merge(self, ps):
        self._param.update(ps._param)

    def idx_size(self):
        return len(self._param)

    def update_all(self, inplace_field, in_field1, in_field2, d1, d2, dt, n):
        for idx, param in self._param.iteritems():
            self._update(inplace_field, in_field1, in_field2, d1, d2, dt, n,
                        idx, param)

    def _update(self, inplace_field, in_field1, in_field2, d1, d2, dt, n, 
               idx, param):
        raise NotImplementedError


class DipoleParam(PwSourceParam):
    def __init__(self, src_time=None, amp=1):
        self.src_time = src_time
        self.amp = float(amp)


class DipoleElectric(PwSource):
    def _update(self, e, h1, h2, dr1, dr2, dt, n, idx, param):
        src_t = param.amp * param.src_time.dipole(dt * n)
        e[idx] = src_t
        
        
class DipoleEx(DipoleElectric): pass


class DipoleEy(DipoleElectric): pass


class DipoleEz(DipoleElectric): pass


class DipoleMagnetic(PwSource):
    def _update(self, h, e1, e2, dr1, dr2, dt, n, idx, param):
        src_t = param.amp * param.src_time.dipole(dt * n)
        h[idx] = src_t

        
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


class TransparentParam(PwSourceParam):
    def __init__(self, amp, aux_fdtd, directional):
        self.aux_fdtd = aux_fdtd
            
        self.face_list = [directional]
        self.amp = {directional: amp}
        

class TransparentElectricParam(TransparentParam):
    def __init__(self, eps_inf, amp, aux_fdtd, samp_pnt, directional):
        TransparentParam.__init__(self, amp, aux_fdtd, directional)

        self.eps_inf = float(eps_inf)
        
        samp_idx = aux_fdtd.space.spc_to_exact_hy_idx(*samp_pnt)
        self.samp_idx0 = {directional: tuple(np.floor(samp_idx))}
        self.samp_idx1 = {directional: tuple(np.floor(samp_idx) + (0, 0, 1))}
        
        r1_value = samp_idx[2] - floor(samp_idx[2])
        self.r1 = {directional: r1_value}
        self.r0 = {directional: 1 - r1_value}


class TransparentMagneticParam(TransparentParam):
    def __init__(self, mu_inf, amp, aux_fdtd, samp_pnt, directional):
        TransparentParam.__init__(self, amp, aux_fdtd, directional)

        self.mu_inf = float(mu_inf)

        samp_idx = aux_fdtd.space.spc_to_exact_ex_idx(*samp_pnt)
        self.samp_idx0 = {directional: tuple(np.floor(samp_idx))}
        self.samp_idx1 = {directional: tuple(np.floor(samp_idx) + (0, 0, 1))}
        
        r1_value = samp_idx[2] - floor(samp_idx[2])
        self.r1 = {directional: r1_value}
        self.r0 = {directional: 1 - r1_value}


class TransparentEx(PwSource):
    def __init__(self):
        PwSource.__init__(self)
        self._consist_cond = {const.MinusY: self._consistency_minus_y,
                              const.MinusZ: self._consistency_minus_z,
                              const.PlusY: self._consistency_plus_y,
                              const.PlusZ: self._consistency_plus_z}
        
    def _update(self, ex, hz, hy, dy, dz, dt, n, idx, param):
        for face in param.face_list:
            self._consist_cond[face](ex, hz, hy, dy, dz, dt, face, 
                                     idx, param)
        
    def _consistency_minus_y(self, ex, hz, hy, dy, dz, dt, face, idx, param):
        incident_hz = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] +
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
                
        ex[idx] -= dt / (param.eps_inf * dy) * param.amp[face] * incident_hz

    def _consistency_plus_y(self, ex, hz, hy, dy, dz, dt, face, idx, param):
        incident_hz = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ex[idx] += dt / (param.eps_inf * dy) * param.amp[face] * incident_hz
    
    def _consistency_minus_z(self, ex, hz, hy, dy, dz, dt, face, idx, param):
        incident_hy = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ex[idx] += dt / (param.eps_inf * dz) * param.amp[face] * incident_hy

    def _consistency_plus_z(self, ex, hz, hy, dy, dz, dt, face, idx, param):
        incident_hy = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] +
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ex[idx] -= dt / (param.eps_inf * dz) * param.amp[face] * incident_hy
        
        
class TransparentEy(PwSource):
    def __init__(self):
        PwSource.__init__(self)
        self._consist_cond = {const.MinusZ: self._consistency_minus_z,
                              const.MinusX: self._consistency_minus_x,
                              const.PlusZ: self._consistency_plus_z,
                              const.PlusX: self._consistency_plus_x}
    
    def _update(self, ey, hx, hz, dz, dx, dt, n, idx, param):
        for face in param.face_list:
            self._consist_cond[face](ey, hx, hz, dz, dx, dt, face, idx, param)
        
    def _consistency_minus_z(self, ey, hx, hz, dz, dx, dt, face, idx, param):
        incident_hx = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ey[idx] -= dt / (param.eps_inf * dz) * param.amp[face] * incident_hx

    def _consistency_minus_x(self, ey, hx, hz, dz, dx, dt, face, idx, param):
        incident_hz = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ey[idx] += dt / (param.eps_inf * dx) * param.amp[face] * incident_hz
    
    def _consistency_plus_z(self, ey, hx, hz, dz, dx, dt, face, idx, param):
        incident_hx = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ey[idx] += dt / (param.eps_inf * dz) * param.amp[face] * incident_hx

    def _consistency_plus_x(self, ey, hx, hz, dz, dx, dt, face, idx, param):
        incident_hz = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ey[idx] -= dt / (param.eps_inf * dx) * param.amp[face] * incident_hz
        

class TransparentEz(PwSource):
    def __init__(self):
        PwSource.__init__(self)
        self._consist_cond = {const.MinusX: self._consistency_minus_x,
                              const.MinusY: self._consistency_minus_y,
                              const.PlusX: self._consistency_plus_x,
                              const.PlusY: self._consistency_plus_y}
        
    def _update(self, ez, hy, hx, dx, dy, dt, n, idx, param):
        for face in param.face_list:
            self._consist_cond[face](ez, hy, hx, dx, dy, dt, face, idx, param)

    def _consistency_minus_x(self, ez, hy, hx, dx, dy, dt, face, idx, param):
        incident_hy = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ez[idx] -= dt / (param.eps_inf * dx) * param.amp[face] * incident_hy

    def _consistency_minus_y(self, ez, hy, hx, dx, dy, dt, face, idx, param):
        incident_hx = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ez[idx] += dt / (param.eps_inf * dy) * param.amp[face] * incident_hx

    def _consistency_plus_x(self, ez, hy, hx, dx, dy, dt, face, idx, param):
        incident_hy = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ez[idx] += dt / (param.eps_inf * dx) * param.amp[face] * incident_hy

    def _consistency_plus_y(self, ez, hy, hx, dx, dy, dt, face, idx, param):
        incident_hx = (param.r0[face] * param.aux_fdtd.hy[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.hy[param.samp_idx1[face]])
        
        ez[idx] -= dt / (param.eps_inf * dy) * param.amp[face] * incident_hx

        
class TransparentHx(PwSource):
    def __init__(self):
        PwSource.__init__(self)
        self._consist_cond = {const.MinusY: self._consistency_minus_y,
                              const.MinusZ: self._consistency_minus_z,
                              const.PlusY: self._consistency_plus_y,
                              const.PlusZ: self._consistency_plus_z}
    
    def _update(self, hx, ez, ey, dy, dz, dt, n, idx, param):
        for face in param.face_list:
            self._consist_cond[face](hx, ez, ey, dy, dz, dt, face, idx, param)
        
    def _consistency_minus_y(self, hx, ez, ey, dy, dz, dt, face, idx, param):
        incident_ez = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] +
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])

        hx[idx] += dt / (param.mu_inf * dy) * param.amp[face] * incident_ez

    def _consistency_plus_y(self, hx, ez, ey, dy, dz, dt, face, idx, param):
        incident_ez = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])

        hx[idx] -= dt / (param.mu_inf * dy) * param.amp[face] * incident_ez

    def _consistency_minus_z(self, hx, ez, ey, dy, dz, dt, face, idx, param):
        incident_ey = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hx[idx] -= dt / (param.mu_inf * dz) * param.amp[face] * incident_ey

    def _consistency_plus_z(self, hx, ez, ey, dy, dz, dt, face, idx, param):
        incident_ey = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] +
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hx[idx] += dt / (param.mu_inf * dz) * param.amp[face] * incident_ey
        

class TransparentHy(PwSource):
    def __init__(self):
        PwSource.__init__(self)
        self._consist_cond = {const.MinusZ: self._consistency_minus_z,
                              const.MinusX: self._consistency_minus_x,
                              const.PlusZ: self._consistency_plus_z,
                              const.PlusX: self._consistency_plus_x}
    
    def _update(self, hy, ex, ez, dz, dx, dt, n, idx, param):
        for face in param.face_list:
            self._consist_cond[face](hy, ex, ez, dz, dx, dt, face, idx, param)
        
    def _consistency_minus_z(self, hy, ex, ez, dz, dx, dt, face, idx, param):
        incident_ex = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hy[idx] += dt / (param.mu_inf * dz) * param.amp[face] * incident_ex

    def _consistency_minus_x(self, hy, ex, ez, dz, dx, dt, face, idx, param):
        incident_ez = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hy[idx] -= dt / (param.mu_inf * dx) * param.amp[face] * incident_ez
    
    def _consistency_plus_z(self, hy, ex, ez, dz, dx, dt, face, idx, param):
        incident_ex = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hy[idx] -= dt / (param.mu_inf * dz) * param.amp[face] * incident_ex

    def _consistency_plus_x(self, hy, ex, ez, dz, dx, dt, face, idx, param):
        incident_ez = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hy[idx] += dt / (param.mu_inf * dx) * param.amp[face] * incident_ez
        

class TransparentHz(PwSource):
    def __init__(self):
        PwSource.__init__(self)
        self._consist_cond = {const.MinusX: self._consistency_minus_x,
                              const.MinusY: self._consistency_minus_y,
                              const.PlusX: self._consistency_plus_x,
                              const.PlusY: self._consistency_plus_y}
    
    def _update(self, hz, ey, ex, dx, dy, dt, n, idx, param):
        for face in param.face_list:
            self._consist_cond[face](hz, ey, ex, dx, dy, dt, face, idx, param)
        
    def _consistency_minus_x(self, hz, ey, ex, dx, dy, dt, face, idx, param):
        incident_ey = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
    
        hz[idx] += dt / (param.mu_inf * dx) * param.amp[face] * incident_ey

    def _consistency_minus_y(self, hz, ey, ex, dx, dy, dt, face, idx, param):
        incident_ex = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hz[idx] -= dt / (param.mu_inf * dy) * param.amp[face] * incident_ex
    
    def _consistency_plus_x(self, hz, ey, ex, dx, dy, dt, face, idx, param):
        incident_ey = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hz[idx] -= dt / (param.mu_inf * dx) * param.amp[face] * incident_ey

    def _consistency_plus_y(self, hz, ey, ex, dx, dy, dt, face, idx, param):
        incident_ex = (param.r0[face] * param.aux_fdtd.ex[param.samp_idx0[face]] + 
                       param.r1[face] * param.aux_fdtd.ex[param.samp_idx1[face]])
        
        hz[idx] += dt / (param.mu_inf * dy) * param.amp[face] * incident_ex
