#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    pass

try:
    from threading import Thread
except ImportError:
    stderr.write('No module named threading. Using dummy_threading instead.\n')
    from dummy_threading import Thread

from numpy import array, empty, arange, ndindex, linspace
import numpy as np

from sys import modules
if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')

from matplotlib.pyplot import new_figure_manager, cm, show

# GMES modules
from constant import *
from geometry import in_range


class ShowLine(Thread):
    """Animated 1-D on-time display. 
    
    """
    def __init__(self, fdtd, component, start, end, vrange, interval, title, fig_id):
        """Constructor.

        Argumetns:
            fdtd: a FDTD instance.
            component: Specify electric or magnetic field component. 
                This should be one of the gmes.constant.Component. 
            start: The start point of the probing line.
            end: The end point of the probing line.
            vrange: Plot range of the y axis.
            interval: Refresh rate of the plot in milliseconds.
            title: Title string of the figure.
            fig_id:

        """
        Thread.__init__(self)

        self.vrange = array(vrange, np.double)
        self.interval = int(interval)
        self.time_step = fdtd.time_step
        self.note_form = 'time: %f'
        self.title = str(title)
        self.id = fig_id

        comp = component

        spc2idx = {Ex: fdtd.space.space_to_ex_index,
                   Ey: fdtd.space.space_to_ey_index,
                   Ez: fdtd.space.space_to_ez_index,
                   Hx: fdtd.space.space_to_hx_index,
                   Hy: fdtd.space.space_to_hy_index,
                   Hz: fdtd.space.space_to_hz_index}
        
        idx2spc = {Ex: fdtd.space.ex_index_to_space,
                   Ey: fdtd.space.ey_index_to_space,
                   Ez: fdtd.space.ez_index_to_space,
                   Hx: fdtd.space.hx_index_to_space,
                   Hy: fdtd.space.hy_index_to_space,
                   Hz: fdtd.space.hz_index_to_space}
            
        field = fdtd.field[comp].real

        if issubclass(comp, Electric):
            start_bndry_idx = (0, 0, 0)            
            if comp is Ex:
                end_bndry_idx = (field.shape[0] - 1, field.shape[1] - 2,
                                    field.shape[2] - 2)
            elif comp is Ey:
                end_bndry_idx = (field.shape[0] - 2, field.shape[1] - 1,
                                    field.shape[2] - 2)
            elif comp is Ez:
                end_bndry_idx = (field.shape[0] - 2, field.shape[1] - 2,
                                    field.shape[2] - 1)
        elif issubclass(comp, Magnetic):
            end_bndry_idx = [i - 1 for i in field.shape]
            if comp is Hx:
                start_bndry_idx = (0, 1, 1)
            elif comp is Hy:
                start_bndry_idx = (1, 0, 1)
            elif comp is Hz:
                start_bndry_idx = (1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)

        start_idx = array(spc2idx[comp](*start), int)
        end_idx = array(spc2idx[comp](*end), int)
        label = 'x', 'y', 'z'
        for i, v in enumerate(end_idx - start_idx):
            if v > 0:
                for j in (j for j in xrange(3) if j != i):
                    tmp1_idx = array(start_bndry_idx, int)
                    tmp1_idx[j] = start_idx[j]
                    if in_range(tmp1_idx, field.shape, comp) is False:
                        return None
                    tmp2_idx = array(end_bndry_idx, int)
                    tmp2_idx[j] = end_idx[j]
                    if in_range(tmp2_idx, field.shape, comp) is False:
                        return None
                    
                if in_range(start_idx, field.shape, comp) is False:
                    if start_idx[i] > end_bndry_idx[i]:
                        return None
                    start_idx[i] = start_bndry_idx[i]
                    if in_range(start_idx, field.shape, comp) is False:
                        return None
                if in_range(end_idx, field.shape, comp) is False:
                    if end_idx[i] < start_bndry_idx[i]:
                        return None
                    end_idx[i] = end_bndry_idx[i]
                    if in_range(end_idx, field.shape, comp) is False:
                        return None
                
                if i == 0:
                    self.ydata = field[start_idx[0]:(end_idx[0] + 1), 
                                       start_idx[1], start_idx[2]]
                elif i == 1:
                    self.ydata = field[start_idx[0], 
                                       start_idx[1]:(end_idx[1] + 1), 
                                       start_idx[2]]
                elif i == 2:
                    self.ydata = field[start_idx[0], start_idx[1], 
                                       start_idx[2]:(end_idx[2] + 1)]

                local_start = idx2spc[comp](*start_idx)
                local_end = idx2spc[comp](*end_idx)
                self.xdata = linspace(local_start[i], local_end[i], 
                                      end_idx[i] - start_idx[i] + 1)
                self.xlabel = label[i]

                break

        self.ylabel = 'displacement'        
        self.window_title = 'GMES' + ' ' + str(fdtd.space.cart_comm.topo[2])
        
    def animate(self):
        self.line.set_ydata(self.ydata)
        self.line.recache()
        self.time_note.set_text(self.note_form % self.time_step.t)
        self.manager.canvas.draw()
        if self.manager.window is not None:
            self.manager.window.after(self.interval, self.animate)

    def run(self):
        self.manager = new_figure_manager(self.id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.line, = ax.plot(self.xdata, self.ydata)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_xlim(self.xdata[0], self.xdata[-1])
        ax.set_ylim(self.vrange[0], self.vrange[1])
        ax.grid(True)
        ax.set_title(self.title)
        self.time_note = self.manager.canvas.figure.text(.6, .92, self.note_form % self.time_step.t)
        self.manager.window.title(self.window_title)
        self.animate()
        self.manager.show()
        self.manager.window.mainloop()


class ShowPlane(Thread):
    """Animated 2-D on-time display.
    
    """
    def __init__(self, fdtd, component, axis, cut, vrange, interval, title, fig_id):
        """Constructor.

        Arguments:
            fdtd: a FDTD instance
            component: Specify electric or magnetic field component. 
                This should be one of the gmes.constant.Component. 
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constant.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis. 
            vrange: Specify the colorbar range.
            inerval: Refresh rates in millisecond.
            title: title string of the figure.
            fig_id:

        """
        Thread.__init__(self)
        
        self.vrange = array(vrange, np.double)
        self.time_step = fdtd.time_step
        self.title = str(title)
        self.interval = int(interval)
        self.id = fig_id
        
        comp = component
        
        spc2idx = {Ex: fdtd.space.space_to_ex_index,
                   Ey: fdtd.space.space_to_ey_index,
                   Ez: fdtd.space.space_to_ez_index,
                   Hx: fdtd.space.space_to_hx_index,
                   Hy: fdtd.space.space_to_hy_index,
                   Hz: fdtd.space.space_to_hz_index}
        
        idx2spc = {Ex: fdtd.space.ex_index_to_space,
                   Ey: fdtd.space.ey_index_to_space,
                   Ez: fdtd.space.ez_index_to_space,
                   Hx: fdtd.space.hx_index_to_space,
                   Hy: fdtd.space.hy_index_to_space,
                   Hz: fdtd.space.hz_index_to_space}
        
        field = fdtd.field[comp].real

        if issubclass(comp, Electric):
            start_bndry_idx = (0, 0, 0)            
            if comp is Ex:
                end_bndry_idx = (field.shape[0] - 1, 
                                 field.shape[1] - 2,
                                 field.shape[2] - 2)
            elif comp is Ey:
                end_bndry_idx = (field.shape[0] - 2, 
                                 field.shape[1] - 1,
                                 field.shape[2] - 2)
            elif comp is Ez:
                end_bndry_idx = (field.shape[0] - 2, 
                                 field.shape[1] - 2,
                                 field.shape[2] - 1)
        elif issubclass(comp, Magnetic):
            end_bndry_idx = [i - 1 for i in field.shape]
            if comp is Hx:
                start_bndry_idx = (0, 1, 1)
            elif comp is Hy:
                start_bndry_idx = (1, 0, 1)
            elif comp is Hz:
                start_bndry_idx = (1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)

        start_bndry_spc = idx2spc[comp](*start_bndry_idx)
        end_bndry_spc = idx2spc[comp](*end_bndry_idx)

        axis2int = {X:0, Y:1, Z:2}
        axis_int = axis2int[axis]
        cut_spc = array(end_bndry_spc, np.double)
        cut_spc[axis_int] = cut
        cut_idx = spc2idx[comp](*cut_spc)
        if in_range(cut_idx, field.shape, comp) is False:
            return None

        if axis is X:
            self.xlabel = 'z'
            self.ylabel = 'y'
            self.extent = (start_bndry_spc[2], end_bndry_spc[2], 
                           end_bndry_spc[1], start_bndry_spc[1])
        elif axis is Y:
            self.xlabel = 'z'
            self.ylabel = 'x'
            self.extent = (start_bndry_spc[2], end_bndry_spc[2], 
                           end_bndry_spc[0], start_bndry_spc[0])
        elif axis is Z:
            self.xlabel = 'y'
            self.ylabel = 'x'
            self.extent = (start_bndry_spc[1], end_bndry_spc[1], 
                           end_bndry_spc[0], start_bndry_spc[0])

        if axis is X:
            self.data = field[cut_idx[0], 
                              start_bndry_idx[1]:end_bndry_idx[1], 
                              start_bndry_idx[2]:end_bndry_idx[2]]
        elif axis is Y:
            self.data = field[start_bndry_idx[0]:end_bndry_idx[0], 
                              cut_idx[1], 
                              start_bndry_idx[2]:end_bndry_idx[2]]
        elif axis is Z:
            self.data = field[start_bndry_idx[0]:end_bndry_idx[0], 
                              start_bndry_idx[1]:end_bndry_idx[1], 
                              cut_idx[2]]
        else:
            msg = "axis must be gmes.constant.Directional."
            raise ValueError(msg)

        self.window_title = 'GMES' + ' ' + str(fdtd.space.cart_comm.topo[2])
        self.note_form = 'time: %f'
        
    def animate(self):
        self.im.set_data(self.data)
        self.time_note.set_text(self.note_form % self.time_step.t)
        self.manager.canvas.draw()
        if self.manager.window is not None:
            self.manager.window.after(self.interval, self.animate)
        
    def run(self):
        self.manager = new_figure_manager(self.id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, aspect='auto', 
                            vmin=self.vrange[0], vmax=self.vrange[1],
                            cmap=cm.RdBu)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        self.manager.canvas.figure.colorbar(self.im)
        self.time_note = self.manager.canvas.figure.text(.6, .92, self.note_form % self.time_step.t)
        self.manager.window.title(self.window_title)
        self.animate()
        self.manager.show()
        self.manager.window.mainloop()


class Snapshot(Thread):
    """A snapshot of 2-D display. 
    
    In this moment, Snapshot is only used to show the structures. 

    """
    def __init__(self, fdtd, component, axis, cut, vrange, title, fig_id):
        """Constructor.

        Arguments:
            fdtd: a FDTD instance
            component: Specify electric or magnetic field component. 
                This should be one of the gmes.constant.Component.
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constant.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis. 
            vrange: Specify the colorbar range. a tuple of length two.
            title: title string of the figure.
            fig_id:

        """
        Thread.__init__(self)
        
        comp = component

        spc2idx = {Ex: fdtd.space.space_to_ex_index,
                   Ey: fdtd.space.space_to_ey_index,
                   Ez: fdtd.space.space_to_ez_index,
                   Hx: fdtd.space.space_to_hx_index,
                   Hy: fdtd.space.space_to_hy_index,
                   Hz: fdtd.space.space_to_hz_index}
        
        idx2spc = {Ex: fdtd.space.ex_index_to_space,
                   Ey: fdtd.space.ey_index_to_space,
                   Ez: fdtd.space.ez_index_to_space,
                   Hx: fdtd.space.hx_index_to_space,
                   Hy: fdtd.space.hy_index_to_space,
                   Hz: fdtd.space.hz_index_to_space}

        material = fdtd.pw_material[comp]
        field = fdtd.field[comp]

        if issubclass(comp, Electric):
            start_bndry_idx = (0, 0, 0)            
            if comp is Ex:
                end_bndry_idx = (field.shape[0] - 1, 
                                 field.shape[1] - 2,
                                 field.shape[2] - 2)
            elif comp is Ey:
                end_bndry_idx = (field.shape[0] - 2, 
                                 field.shape[1] - 1,
                                 field.shape[2] - 2)
            elif comp is Ez:
                end_bndry_idx = (field.shape[0] - 2, 
                                 field.shape[1] - 2,
                                 field.shape[2] - 1)
        elif issubclass(comp, Magnetic):
            end_bndry_idx = [i - 1 for i in field.shape]
            if comp is Hx:
                start_bndry_idx = idx2spc[comp](0, 1, 1)
            elif comp is Hy:
                start_bndry_idx = idx2spc[comp](1, 0, 1)
            elif comp is Hz:
                start_bndry_idx = idx2spc[comp](1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)
    
        start_bndry_spc = idx2spc[comp](*start_bndry_idx)
        end_bndry_spc = idx2spc[comp](*end_bndry_idx)

        axis2int = {X:0, Y:1, Z:2}
        axis_int = axis2int[axis]
        cut_spc = array(end_bndry_spc, np.double)
        cut_spc[axis_int] = cut
        cut_idx = spc2idx[comp](*cut_spc)
        if in_range(cut_idx, field.shape, comp) is False:
            return None

        if axis is X:
            self.xlabel = 'z'
            self.ylabel = 'y'
            self.extent = (start_bndry_spc[2], end_bndry_spc[2], 
                           end_bndry_spc[1], start_bndry_spc[1])
        elif axis is Y:
            self.xlabel = 'z'
            self.ylabel = 'x'
            self.extent = (start_bndry_spc[2], end_bndry_spc[2], 
                           end_bndry_spc[0], start_bndry_spc[0])
        elif axis is Z:
            self.xlabel = 'y'
            self.ylabel = 'x'
            self.extent = (start_bndry_spc[1], end_bndry_spc[1], 
                           end_bndry_spc[0], start_bndry_spc[0])

        data_shape_3d = array(end_bndry_idx) - array(start_bndry_idx) + 1
        data_shape_2d = [v for i, v in enumerate(data_shape_3d) 
                         if i != axis_int]
        self.data = empty(data_shape_2d, np.double)

        for idx in ndindex(*data_shape_2d):
            mat_idx = empty(3, np.double)
            
            if axis_int == 0:
                mat_idx[1] = idx[0] + start_bndry_idx[1]
                mat_idx[2] = idx[1] + start_bndry_idx[2]
            elif axis_int == 1:
                mat_idx[0] = idx[0] + start_bndry_idx[0]
                mat_idx[2] = idx[1] + start_bndry_idx[2]
            elif axis_int == 2:
                mat_idx[0] = idx[0] + start_bndry_idx[0]
                mat_idx[1] = idx[1] + start_bndry_idx[1]
            
            mat_idx[axis_int] = cut_idx[axis_int]
            for pw_mat in material.itervalues():
                if issubclass(comp, Electric):
                    value = pw_mat.get_eps_inf(tuple(mat_idx))
                elif issubclass(comp, Magnetic): 
                    value = pw_mat.get_mu_inf(tuple(mat_idx))
                if value != 0:
                    self.data[idx] = value

        self.window_title = 'GMES' + ' ' + str(fdtd.space.cart_comm.topo[2])

        if vrange is None:
            data_min = self.data.min()
            data_max = self.data.max()
            if data_min == data_max:
                data_min = 0.5 * data_min
                data_max = 1.5 * data_max
            self.vrange = array((data_min, data_max), np.double)
        else:
            self.vrange = array(vrange, np.double)

        self.title = str(title)

        self.id = fig_id
        
    def run(self):
        self.manager = new_figure_manager(self.id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, 
                            aspect='auto', 
                            vmin=self.vrange[0], 
                            vmax=self.vrange[1], 
                            cmap=cm.bone)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        self.manager.canvas.figure.colorbar(self.im)
        self.manager.window.title(self.window_title)
        
        self.manager.canvas.draw()
        self.manager.show()
        self.manager.window.mainloop()
