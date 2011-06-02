#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from math import cos, sin, pi, floor

import numpy as np
from numpy import inf

import constants as const


class DipoleElectric(object):
    def __init__(self, pw_material, src_time=None, amp=1, filename=None):
        self.pw_material = pw_material
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.epsilon = pw_material.epsilon
        self.src_time = src_time
        self.amp = float(amp)
        if filename:
            self.file = open(filename, 'w')
        else:
            self.file = None
        
    def __del__(self):
        if self.file:
            self.file.close()
        
    def update(self, efield, hfield1, hfield2, space_diff1, space_diff2, dt, n):
        src_t = self.amp * self.src_time.dipole(dt*n)
        efield[self.i, self.j, self.k] = src_t
        if self.file is not None:
            self.file.write(str(n) + ' ' + str(src_t) + '\n')
        
        
class DipoleEx(DipoleElectric): pass


class DipoleEy(DipoleElectric): pass


class DipoleEz(DipoleElectric): pass


class DipoleMagnetic(object):
    def __init__(self, pw_material, src_time=None,amp=1, filename=None):
        self.pw_material = pw_material
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.mu = pw_material.mu
        self.src_time = src_time
        self.amp = float(amp)
        if filename:
            self.file = open(filename, 'w')
        else:
            self.file = None
    
    def __del__(self):
        if self.file:
            self.file.close()
        
    def update(self, hfield, efield1, efield2, space_diff1, space_diff2, dt, n):
        src_t = self.amp * self.src_time.dipole(dt*n)
        hfield[self.i, self.j, self.k] = src_t
        if self.file is not None:
            self.file.write(str(n) + ' ' + str(src_t)+'\n')
        
        
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
    def __init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional):
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.epsilon = epsilon
        
        self.aux_fdtd = aux_fdtd
            
        self.face_list = [directional]
        self.amp = {directional: amp}
        
        samp_idx = aux_fdtd.space.spc_to_exact_hy_idx(*samp_pnt)
        self.samp_idx0 = {directional: tuple(np.floor(samp_idx))}
        self.samp_idx1 = {directional: tuple(np.floor(samp_idx) + (0, 0, 1))}
        
        r1_value = samp_idx[2] - floor(samp_idx[2])
        self.r1 = {directional: r1_value}
        self.r0 = {directional: 1 - r1_value}

        if isinstance(pw_material, TransparentElectric):
            self.pw_material = pw_material.pw_material
            
            self.face_list.extend(pw_material.face_list)
            self.amp.update(pw_material.amp)
            self.samp_idx0.update(pw_material.samp_idx0)
            self.samp_idx1.update(pw_material.samp_idx1)
            self.r1.update(pw_material.r1)
            self.r0.update(pw_material.r0)
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
            

class TransparentEx(TransparentElectric):
    def __init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional):
        TransparentElectric.__init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional)
        
        self._consist_cond = {const.MinusY: self._consistency_minus_y,
                              const.MinusZ: self._consistency_minus_z,
                              const.PlusY: self._consistency_plus_y,
                              const.PlusZ: self._consistency_plus_z}
        
    def display_info(self, indent=0):
        print " " * indent, "TransparentEx"
        TransparentElectric.display_info(self, indent)
         
    def update(self, ex, hz, hy, dy, dz, dt, n):
        self.pw_material.update(ex, hz, hy, dy, dz, dt, n)
        
        for face in self.face_list:
            self._consist_cond[face](ex, hz, hy, dy, dz, dt, face)
        
    def _consistency_minus_y(self, ex, hz, hy, dy, dz, dt, face):
        incident_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
                
        idx = self.i, self.j, self.k
        ex[idx] -= dt / (self.epsilon * dy) * self.amp[face] * incident_hz

    def _consistency_plus_y(self, ex, hz, hy, dy, dz, dt, face):
        incident_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ex[idx] += dt / (self.epsilon * dy) * self.amp[face] * incident_hz
    
    def _consistency_minus_z(self, ex, hz, hy, dy, dz, dt, face):
        incident_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ex[idx] += dt / (self.epsilon * dz) * self.amp[face] * incident_hy

    def _consistency_plus_z(self, ex, hz, hy, dy, dz, dt, face):
        incident_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ex[idx] -= dt / (self.epsilon * dz) * self.amp[face] * incident_hy
        
        
class TransparentEy(TransparentElectric):
    def __init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional):
        TransparentElectric.__init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional)
        
        self._consist_cond = {const.MinusZ: self._consistency_minus_z,
                              const.MinusX: self._consistency_minus_x,
                              const.PlusZ: self._consistency_plus_z,
                              const.PlusX: self._consistency_plus_x}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentEy"
        TransparentElectric.display_info(self, indent)
        
    def update(self, ey, hx, hz, dz, dx, dt, n):
        self.pw_material.update(ey, hx, hz, dz, dx, dt, n)
        
        for face in self.face_list:
            self._consist_cond[face](ey, hx, hz, dz, dx, dt, face)
        
    def _consistency_minus_z(self, ey, hx, hz, dz, dx, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] -= dt / (self.epsilon * dz) * self.amp[face] * incident_hx

    def _consistency_minus_x(self, ey, hx, hz, dz, dx, dt, face):
        incident_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] += dt / (self.epsilon * dx) * self.amp[face] * incident_hz
    
    def _consistency_plus_z(self, ey, hx, hz, dz, dx, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] += dt / (self.epsilon * dz) * self.amp[face] * incident_hx

    def _consistency_plus_x(self, ey, hx, hz, dz, dx, dt, face):
        incident_hz = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ey[idx] -= dt / (self.epsilon * dx) * self.amp[face] * incident_hz
        

class TransparentEz(TransparentElectric):
    def __init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional):
        TransparentElectric.__init__(self, pw_material, epsilon, amp, aux_fdtd, samp_pnt, directional)
        
        self._consist_cond = {const.MinusX: self._consistency_minus_x,
                              const.MinusY: self._consistency_minus_y,
                              const.PlusX: self._consistency_plus_x,
                              const.PlusY: self._consistency_plus_y}
        
    def display_info(self, indent=0):
        print " " * indent, "TransparentEz"
        TransparentElectric.display_info(self, indent)
                
    def update(self, ez, hy, hx, dx, dy, dt, n):
        self.pw_material.update(ez, hy, hx, dx, dy, dt, n)

        for face in self.face_list:
            self._consist_cond[face](ez, hy, hx, dx, dy, dt, face)

    def _consistency_minus_x(self, ez, hy, hx, dx, dy, dt, face):
        incident_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] -= dt / (self.epsilon * dx) * self.amp[face] * incident_hy

    def _consistency_minus_y(self, ez, hy, hx, dx, dy, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] += dt / (self.epsilon * dy) * self.amp[face] * incident_hx

    def _consistency_plus_x(self, ez, hy, hx, dx, dy, dt, face):
        incident_hy = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] += dt / (self.epsilon * dx) * self.amp[face] * incident_hy

        # # DEBUG CODE STARTS
        # if idx == (55,6,0):
        #     f = open("DEBUG_ez_amp_hy.dat", 'a')
        #     f.write(str(self.amp[face]) + '\n')
        #     f.close()
        #     f = open("DEBUG_ez_incident_hy.dat", 'a')
        #     f.write(str(incident_hy) + '\n')
        #     f.close()
        # # DEBUG CODE ENDS

    def _consistency_plus_y(self, ez, hy, hx, dx, dy, dt, face):
        incident_hx = (self.r0[face] * self.aux_fdtd.hy[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.hy[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        ez[idx] -= dt / (self.epsilon * dy) * self.amp[face] * incident_hx

        
class TransparentMagnetic(object):
    def __init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional):
        self.i = pw_material.i
        self.j = pw_material.j
        self.k = pw_material.k
        self.mu = mu
        
        self.aux_fdtd = aux_fdtd
        
        self.face_list = [directional]
        self.amp = {directional: amp}
        
        samp_idx = aux_fdtd.space.spc_to_exact_ex_idx(*samp_pnt)
        self.samp_idx0 = {directional: tuple(np.floor(samp_idx))}
        self.samp_idx1 = {directional: tuple(np.floor(samp_idx) + (0, 0, 1))}
        
        r1_value = samp_idx[2] - floor(samp_idx[2])
        self.r1 = {directional: r1_value}
        self.r0 = {directional: 1 - r1_value}

        if isinstance(pw_material, TransparentMagnetic):
            self.pw_material = pw_material.pw_material
            
            self.face_list.extend(pw_material.face_list)
            self.amp.update(pw_material.amp)
            self.samp_idx0.update(pw_material.samp_idx0)
            self.samp_idx1.update(pw_material.samp_idx1)
            self.r1.update(pw_material.r1)
            self.r0.update(pw_material.r0)
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
   

class TransparentHx(TransparentMagnetic):
    def __init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional):
        TransparentMagnetic.__init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional)
        
        self._consist_cond = {const.MinusY: self._consistency_minus_y,
                              const.MinusZ: self._consistency_minus_z,
                              const.PlusY: self._consistency_plus_y,
                              const.PlusZ: self._consistency_plus_z}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentHx"
        TransparentMagnetic.display_info(self, indent)
        
    def update(self, hx, ez, ey, dy, dz, dt, n):
        self.pw_material.update(hx, ez, ey, dy, dz, dt, n)

        for face in self.face_list:
            self._consist_cond[face](hx, ez, ey, dy, dz, dt, face)
        
    def _consistency_minus_y(self, hx, ez, ey, dy, dz, dt, face):
        incident_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])

        idx = self.i, self.j, self.k
        hx[idx] += dt / (self.mu * dy) * self.amp[face] * incident_ez

    def _consistency_plus_y(self, hx, ez, ey, dy, dz, dt, face):
        incident_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])

        idx = self.i, self.j, self.k
        hx[idx] -= dt / (self.mu * dy) * self.amp[face] * incident_ez

    def _consistency_minus_z(self, hx, ez, ey, dy, dz, dt, face):
        incident_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hx[idx] -= dt / (self.mu * dz) * self.amp[face] * incident_ey

    def _consistency_plus_z(self, hx, ez, ey, dy, dz, dt, face):
        incident_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] +
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hx[idx] += dt / (self.mu * dz) * self.amp[face] * incident_ey
        

class TransparentHy(TransparentMagnetic):
    def __init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional):
        TransparentMagnetic.__init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional)
        
        self._consist_cond = {const.MinusZ: self._consistency_minus_z,
                              const.MinusX: self._consistency_minus_x,
                              const.PlusZ: self._consistency_plus_z,
                              const.PlusX: self._consistency_plus_x}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentHy"
        TransparentMagnetic.display_info(self, indent)
        
    def update(self, hy, ex, ez, dz, dx, dt, n):
        self.pw_material.update(hy, ex, ez, dz, dx, dt, n)

        for face in self.face_list:
            self._consist_cond[face](hy, ex, ez, dz, dx, dt, face)
        
    def _consistency_minus_z(self, hy, ex, ez, dz, dx, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] += dt / (self.mu * dz) * self.amp[face] * incident_ex

    def _consistency_minus_x(self, hy, ex, ez, dz, dx, dt, face):
        incident_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] -= dt / (self.mu * dx) * self.amp[face] * incident_ez
    
    def _consistency_plus_z(self, hy, ex, ez, dz, dx, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] -= dt / (self.mu * dz) * self.amp[face] * incident_ex

    def _consistency_plus_x(self, hy, ex, ez, dz, dx, dt, face):
        incident_ez = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hy[idx] += dt / (self.mu * dx) * self.amp[face] * incident_ez
        
        # # DEBUG CODE STARTS
        # if idx == (56,6,1):
        #     f = open("DEBUG_hy_amp_ez.dat", 'a')
        #     f.write(str(self.amp[face]) + '\n')
        #     f.close()
        #     f = open("DEBUG_hy_incident_ez.dat", 'a')
        #     f.write(str(incident_ez) + '\n')
        #     f.close()
        # # DEBUG CODE ENDS

class TransparentHz(TransparentMagnetic):
    def __init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional):
        TransparentMagnetic.__init__(self, pw_material, mu, amp, aux_fdtd, samp_pnt, directional)
        
        self._consist_cond = {const.MinusX: self._consistency_minus_x,
                              const.MinusY: self._consistency_minus_y,
                              const.PlusX: self._consistency_plus_x,
                              const.PlusY: self._consistency_plus_y}
    
    def display_info(self, indent=0):
        print " " * indent, "TransparentHz"
        TransparentMagnetic.display_info(self, indent)
                
    def update(self, hz, ey, ex, dx, dy, dt, n):
        self.pw_material.update(hz, ey, ex, dx, dy, dt, n)
        
        for face in self.face_list:
            self._consist_cond[face](hz, ey, ex, dx, dy, dt, face)
        
    def _consistency_minus_x(self, hz, ey, ex, dx, dy, dt, face):
        incident_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
    
        idx = self.i, self.j, self.k
        hz[idx] += dt / (self.mu * dx) * self.amp[face] * incident_ey

    def _consistency_minus_y(self, hz, ey, ex, dx, dy, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hz[idx] -= dt / (self.mu * dy) * self.amp[face] * incident_ex
    
    def _consistency_plus_(self, hz, ey, ex, dx, dy, dt, face):
        incident_ey = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hz[idx] -= dt / (self.mu * dx) * self.amp[face] * incident_ey

    def _consistency_plus_y(self, hz, ey, ex, dx, dy, dt, face):
        incident_ex = (self.r0[face] * self.aux_fdtd.ex[self.samp_idx0[face]] + 
                       self.r1[face] * self.aux_fdtd.ex[self.samp_idx1[face]])
        
        idx = self.i, self.j, self.k
        hz[idx] += dt / (self.mu * dy) * self.amp[face] * incident_ex
        
