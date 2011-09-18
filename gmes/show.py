#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import stderr

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except ImportError:
    stderr.write('No module named psyco. Execution speed might be slow.\n')

try:
    from threading import Thread
except ImportError:
    stderr.write('No module named threading. Using dummy_threading instead.\n')
    from dummy_threading import Thread

from numpy import array, empty, arange, ndindex, linspace

from sys import modules
if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')

from matplotlib.pyplot import new_figure_manager, cm, show

# GMES modules
import constant as const
from geometry import in_range


class ShowLine(Thread):
    """Animated 1-D on-time display. 
    
    """
    def __init__(self, fdtd, component, start, end, vrange, 
                 interval, title, fig_id):
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

        self.vrange = array(vrange, float)
        self.interval = int(interval)
        self.time_step = fdtd.time_step
        self.note_form = 'time: %f'
        self.title = str(title)
        self.id = fig_id

        comp = component

        spc2idx = {const.Ex: fdtd.space.space_to_ex_index,
                   const.Ey: fdtd.space.space_to_ey_index,
                   const.Ez: fdtd.space.space_to_ez_index,
                   const.Hx: fdtd.space.space_to_hx_index,
                   const.Hy: fdtd.space.space_to_hy_index,
                   const.Hz: fdtd.space.space_to_hz_index}
        
        idx2spc = {const.Ex: fdtd.space.ex_index_to_space,
                   const.Ey: fdtd.space.ey_index_to_space,
                   const.Ez: fdtd.space.ez_index_to_space,
                   const.Hx: fdtd.space.hx_index_to_space,
                   const.Hy: fdtd.space.hy_index_to_space,
                   const.Hz: fdtd.space.hz_index_to_space}
            
        fields = {const.Ex: fdtd.ex, const.Ey: fdtd.ey, const.Ez: fdtd.ez,
                  const.Hx: fdtd.hx, const.Hy: fdtd.hy, const.Hz: fdtd.hz}

        field = fields[comp].real

        if issubclass(comp, const.Electric):
            start_boundary_idx = (0, 0, 0)            
            if comp is const.Ex:
                end_boundary_idx = (field.shape[0] - 1, field.shape[1] - 2,
                                    field.shape[2] - 2)
            elif comp is const.Ey:
                end_boundary_idx = (field.shape[0] - 2, field.shape[1] - 1,
                                    field.shape[2] - 2)
            elif comp is const.Ez:
                end_boundary_idx = (field.shape[0] - 2, field.shape[1] - 2,
                                    field.shape[2] - 1)
        elif issubclass(comp, const.Magnetic):
            end_boundary_idx = [i - 1 for i in field.shape]
            if comp is const.Hx:
                start_boundary_idx = idx2spc[comp](0, 1, 1)
            elif comp is const.Hy:
                start_boundary_idx = idx2spc[comp](1, 0, 1)
            elif comp is const.Hz:
                start_boundary_idx = idx2spc[comp](1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)

        start_idx = array(spc2idx[comp](*start), int)
        end_idx = array(spc2idx[comp](*end), int)
        label = 'x', 'y', 'z'
        for i, v in enumerate(end_idx - start_idx):
            if v > 0:
                for j in (j for j in xrange(3) if j != i):
                    tmp1_idx = array(start_boundary_idx, int)
                    tmp1_idx[j] = start_idx[j]
                    if in_range(tmp1_idx, field.shape, comp) is False:
                        return None
                    tmp2_idx = array(end_boundary_idx, int)
                    tmp2_idx[j] = end_idx[j]
                    if in_range(tmp2_idx, field.shape, comp) is False:
                        return None
                    
                if in_range(start_idx, field.shape, comp) is False:
                    if start_idx[i] > end_boundary_idx[i]:
                        return None
                    start_idx[i] = start_boundary_idx[i]
                    if in_range(start_idx, field.shape, comp) is False:
                        return None
                if in_range(end_idx, field.shape, comp) is False:
                    if end_idx[i] < start_boundary_idx[i]:
                        return None
                    end_idx[i] = end_boundary_idx[i]
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
        self.ax = self.manager.canvas.figure.add_subplot(111)
        self.line, = self.ax.plot(self.xdata, self.ydata)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_xlim(self.xdata[0], self.xdata[-1])
        self.ax.set_ylim(self.vrange[0], self.vrange[1])
        self.ax.grid(True)
        self.ax.set_title(self.title)
        self.time_note = \
            self.manager.canvas.figure.text(.6, .92, self.note_form % 
                                             self.time_step.t)
        self.manager.window.title(self.window_title)
        
        self.animate()
        self.manager.show()
        self.manager.window.mainloop()


class ShowPlane(Thread):
    """Animated 2-D on-time display.
    
    """
    def __init__(self, fdtd, component, axis, cut, vrange, interval, 
                 title, fig_id):
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

        self.vrange = array(vrange, float)
        self.time_step = fdtd.time_step
        self.title = str(title)
        self.interval = int(interval)
        self.id = fig_id

        comp = component
        
        spc2idx = {const.Ex: fdtd.space.space_to_ex_index,
                   const.Ey: fdtd.space.space_to_ey_index,
                   const.Ez: fdtd.space.space_to_ez_index,
                   const.Hx: fdtd.space.space_to_hx_index,
                   const.Hy: fdtd.space.space_to_hy_index,
                   const.Hz: fdtd.space.space_to_hz_index}
        
        idx2spc = {const.Ex: fdtd.space.ex_index_to_space,
                   const.Ey: fdtd.space.ey_index_to_space,
                   const.Ez: fdtd.space.ez_index_to_space,
                   const.Hx: fdtd.space.hx_index_to_space,
                   const.Hy: fdtd.space.hy_index_to_space,
                   const.Hz: fdtd.space.hz_index_to_space}
        
        fields = {const.Ex: fdtd.ex, const.Ey: fdtd.ey, const.Ez: fdtd.ez,
                  const.Hx: fdtd.hx, const.Hy: fdtd.hy, const.Hz: fdtd.hz}

        field = fields[comp].real

        if issubclass(comp, const.Electric):
            start_boundary_idx = (0, 0, 0)            
            if comp is const.Ex:
                end_boundary_idx = (field.shape[0] - 1, field.shape[1] - 2,
                                    field.shape[2] - 2)
            elif comp is const.Ey:
                end_boundary_idx = (field.shape[0] - 2, field.shape[1] - 1,
                                    field.shape[2] - 2)
            elif comp is const.Ez:
                end_boundary_idx = (field.shape[0] - 2, field.shape[1] - 2,
                                    field.shape[2] - 1)
        elif issubclass(comp, const.Magnetic):
            end_boundary_idx = [i - 1 for i in field.shape]
            if comp is const.Hx:
                start_boundary_idx = idx2spc[comp](0, 1, 1)
            elif comp is const.Hy:
                start_boundary_idx = idx2spc[comp](1, 0, 1)
            elif comp is const.Hz:
                start_boundary_idx = idx2spc[comp](1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)

        start_boundary_spc = idx2spc[comp](*start_boundary_idx)
        end_boundary_spc = idx2spc[comp](*end_boundary_idx)

        direct2int = {const.X:0, const.Y:1, const.Z:2}
        axis_int = direct2int[axis]
        self.extent = (start_boundary_spc[(axis_int + 2) % 3], 
                       end_boundary_spc[(axis_int + 2) % 3], 
                       end_boundary_spc[(axis_int + 1) % 3], 
                       start_boundary_spc[(axis_int + 1) % 3])
        cut_spc = array(end_boundary_spc, float)
        cut_spc[axis_int] = cut
        cut_idx = spc2idx[comp](*cut_spc)
        if in_range(cut_idx, field.shape, comp) is False:
            return None
        label = 'x', 'y', 'z'
        self.xlabel = label[(axis_int + 1) % 3]
        self.ylabel = label[(axis_int + 2) % 3]
        if axis is const.X:
            self.data = field[cut_idx[0], 
                              start_boundary_idx[1]:end_boundary_idx[1], 
                              start_boundary_idx[2]:end_boundary_idx[2]]
        elif axis is const.Y:
            self.data = field[start_boundary_idx[0]:end_boundary_idx[0], 
                              cut_idx[1], 
                              start_boundary_idx[2]:end_boundary_idx[2]]
        elif axis is const.Z:
            self.data = field[start_boundary_idx[0]:end_boundary_idx[0], 
                              start_boundary_idx[1]:end_boundary_idx[1], 
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
        self.time_note = \
            self.manager.canvas.figure.text(.6, .92, self.note_form % 
                                             self.time_step.t, animated=True)
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

        spc2idx = {const.Ex: fdtd.space.space_to_ex_index,
                   const.Ey: fdtd.space.space_to_ey_index,
                   const.Ez: fdtd.space.space_to_ez_index,
                   const.Hx: fdtd.space.space_to_hx_index,
                   const.Hy: fdtd.space.space_to_hy_index,
                   const.Hz: fdtd.space.space_to_hz_index}
        
        idx2spc = {const.Ex: fdtd.space.ex_index_to_space,
                   const.Ey: fdtd.space.ey_index_to_space,
                   const.Ez: fdtd.space.ez_index_to_space,
                   const.Hx: fdtd.space.hx_index_to_space,
                   const.Hy: fdtd.space.hy_index_to_space,
                   const.Hz: fdtd.space.hz_index_to_space}

        materials = {const.Ex: fdtd.material_ex, 
                     const.Ey: fdtd.material_ey, 
                     const.Ez: fdtd.material_ez,
                     const.Hx: fdtd.material_hx, 
                     const.Hy: fdtd.material_hy,
                     const.Hz: fdtd.material_hz}

        fields = {const.Ex: fdtd.ex, const.Ey: fdtd.ey, const.Ez: fdtd.ez,
                  const.Hx: fdtd.hx, const.Hy: fdtd.hy, const.Hz: fdtd.hz}

        material = materials[comp]
        field = fields[comp]

        if issubclass(comp, const.Electric):
            start_boundary_idx = (0, 0, 0)            
            if comp is const.Ex:
                end_boundary_idx = (field.shape[0] - 1, field.shape[1] - 2,
                                    field.shape[2] - 2)
            elif comp is const.Ey:
                end_boundary_idx = (field.shape[0] - 2, field.shape[1] - 1,
                                    field.shape[2] - 2)
            elif comp is const.Ez:
                end_boundary_idx = (field.shape[0] - 2, field.shape[1] - 2,
                                    field.shape[2] - 1)
        elif issubclass(comp, const.Magnetic):
            end_boundary_idx = [i - 1 for i in field.shape]
            if comp is const.Hx:
                start_boundary_idx = idx2spc[comp](0, 1, 1)
            elif comp is const.Hy:
                start_boundary_idx = idx2spc[comp](1, 0, 1)
            elif comp is const.Hz:
                start_boundary_idx = idx2spc[comp](1, 1, 0)
 	else:
            msg = "component should be of class constant.Component."
            raise ValueError(msg)
    
        start_boundary_spc = idx2spc[comp](*start_boundary_idx)
        end_boundary_spc = idx2spc[comp](*end_boundary_idx)

        direct2int = {const.X:0, const.Y:1, const.Z:2}
        axis_int = direct2int[axis]
        self.extent = (start_boundary_spc[(axis_int + 2) % 3], 
                       end_boundary_spc[(axis_int + 2) % 3], 
                       end_boundary_spc[(axis_int + 1) % 3], 
                       start_boundary_spc[(axis_int + 1) % 3])
        cut_spc = array(end_boundary_spc, float)
        cut_spc[axis_int] = cut
        cut_idx = spc2idx[comp](*cut_spc)
        if in_range(cut_idx, field.shape, comp) is False:
            return None
        label = 'x', 'y', 'z'
        self.xlabel = label[(axis_int + 1) % 3]
        self.ylabel = label[(axis_int + 2) % 3]

        data_shape_3d = array(end_boundary_idx) - array(start_boundary_idx) + 1
        data_shape_2d = [v for i, v in enumerate(data_shape_3d) 
                         if i != axis_int]
        self.data = empty(data_shape_2d, float)

        for idx in ndindex(*data_shape_2d):
            mat_idx = empty(3, float)
            
            if axis_int == 0:
                mat_idx[1] = idx[0] + start_boundary_idx[1]
                mat_idx[2] = idx[1] + start_boundary_idx[2]
            elif axis_int == 1:
                mat_idx[0] = idx[0] + start_boundary_idx[0]
                mat_idx[2] = idx[1] + start_boundary_idx[2]
            elif axis_int == 2:
                mat_idx[0] = idx[0] + start_boundary_idx[0]
                mat_idx[1] = idx[1] + start_boundary_idx[1]
            
            mat_idx[axis_int] = cut_idx[axis_int]
            for pw_mat in material.itervalues():
                if issubclass(comp, const.Electric):
                    value = pw_mat.get_eps_inf(tuple(mat_idx))
                elif issubclass(comp, const.Magnetic): 
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
            self.vrange = array((data_min, data_max), float)
        else:
            self.vrange = array(vrange, float)

        self.title = str(title)

        self.id = fig_id
        
    def run(self):
        self.manager = new_figure_manager(self.id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, aspect='auto', 
                            vmin=self.vrange[0], vmax=self.vrange[1], 
                            cmap=cm.bone)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        self.manager.canvas.figure.colorbar(self.im)
        self.manager.window.title(self.window_title)
        
        self.manager.canvas.draw()
        self.manager.show()
        self.manager.window.mainloop()
