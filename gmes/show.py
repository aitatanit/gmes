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

from numpy import array, empty, arange, ndindex

from sys import modules
if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')

from matplotlib.pyplot import new_figure_manager, cm, show

# GMES modules
import constants as const
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
                This should be one of the gmes.constants.Component. 
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
        self.fig_id = fig_id

        if component is const.Ex:
            field = fdtd.ex.real
            spc_to_idx = fdtd.space.space_to_ex_index
            idx_to_spc = fdtd.space.ex_index_to_space
            tmp_start_idx = (0, 0, 0)
            tmp_end_idx = (field.shape[0] - 1, field.shape[1] - 2, 
                           field.shape[2] - 2)
        elif component is const.Ey:
            field = fdtd.ey.real
            spc_to_idx = fdtd.space.space_to_ey_index
            idx_to_spc = fdtd.space.ey_index_to_space
            tmp_start_idx = (0, 0, 0)
            tmp_end_idx = (field.shape[0] - 2, field.shape[1] - 1, 
                           field.shape[2] - 2)
        elif component is const.Ez:
            field = fdtd.ez.real
            spc_to_idx = fdtd.space.space_to_ez_index
            idx_to_spc = fdtd.space.ez_index_to_space
            tmp_start_idx = (0, 0, 0)
            tmp_end_idx = (field.shape[0] - 2, field.shape[1] - 2, 
                           field.shape[2] - 1)
        elif component is const.Hx:
            field = fdtd.hx.real
            spc_to_idx = fdtd.space.space_to_hx_index
            idx_to_spc = fdtd.space.hx_index_to_space
            tmp_start_idx = idx_to_spc(0, 1, 1)
            tmp_end_idx = [i - 1 for i in field.shape]
        elif component is const.Hy:
            field = fdtd.hy.real
            spc_to_idx = fdtd.space.space_to_hy_index
            idx_to_spc = fdtd.space.hy_index_to_space
            tmp_start_idx = idx_to_spc(1, 0, 1)
            tmp_end_idx = [i - 1 for i in field.shape]
        elif component is const.Hz:
            field = fdtd.hz.real
            spc_to_idx = fdtd.space.space_to_hz_index
            idx_to_spc = fdtd.space.hz_index_to_space
            tmp_start_idx = idx_to_spc(1, 1, 0)
            tmp_end_idx = [i - 1 for i in field.shape]
        else:
            msg = "component should be of class constants.Component."
            raise ValueError(msg)

        global_start_idx = spc_to_idx(*start)
        global_end_idx = [i + 1 for i in spc_to_idx(*end)]
        
        if global_end_idx[0] - global_start_idx[0] > 1:
            start_idx = (tmp_start_idx[0], global_start_idx[1], 
                         global_start_idx[2])
            end_idx = tmp_end_idx[0], global_end_idx[1], global_end_idx[2]
            if in_range(start_idx, field, component) is False:
                return None
            self.ydata = field[start_idx[0]:end_idx[0], start_idx[1], 
                                start_idx[2]]
            
        elif global_end_idx[1] - global_start_idx[1] > 1:
            start_idx = (global_start_idx[0], tmp_start_idx[1], 
                         global_start_idx[2])
            end_idx = global_end_idx[0], tmp_end_idx[1], global_end_idx[2]
            if in_range(start_idx, field, component) is False:
                return None
            self.ydata = field[start_idx[0], start_idx[1]:end_idx[1], 
                                start_idx[2]]
            
        elif global_end_idx[2] - global_start_idx[2] > 1:
            start_idx = (global_start_idx[0], global_start_idx[1], 
                         tmp_start_idx[2])
            end_idx = global_end_idx[0], global_end_idx[1], tmp_end_idx[2]
            if in_range(start_idx, field, component) is False:
                return None
            self.ydata = field[start_idx[0], start_idx[1], 
                                start_idx[2]:end_idx[2]]
        
        start2 = idx_to_spc(*start_idx)
        end2 = idx_to_spc(*end_idx)
        domain_idx = map(lambda x, y: x - y, end_idx, start_idx)
        for i in range(3):
            if domain_idx[i] != 1 and i == 0:
                step = fdtd.space.dx
                self.xlabel = 'x'
                break
            if domain_idx[i] != 1 and i == 1:
                step = fdtd.space.dy
                self.xlabel = 'y'
                break
            if domain_idx[i] != 1 and i == 2:
                step = fdtd.space.dz
                self.xlabel = 'z'
                break
				
        self.xdata = arange(start2[i], end2[i], step)
        
        if len(self.xdata) > len(self.ydata):
            self.xdata = self.xdata[:-1]
			
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
        self.manager = new_figure_manager(self.fig_id)
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
                This should be one of the gmes.constants.Component. 
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
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
        self.fig_id = fig_id

        if component is const.Ex:
            field = fdtd.ex.real
            spc_to_idx = fdtd.space.space_to_ex_index
            idx_to_spc = fdtd.space.ex_index_to_space
            tmp_cut_coords = idx_to_spc(0, 0, 0)
            
        elif component is const.Ey:
            field = fdtd.ey.real
            spc_to_idx = fdtd.space.space_to_ey_index
            idx_to_spc = fdtd.space.ey_index_to_space
            tmp_cut_coords = idx_to_spc(0, 0, 0)
            
        elif component is const.Ez:
            field = fdtd.ez.real
            spc_to_idx = fdtd.space.space_to_ez_index
            idx_to_spc = fdtd.space.ez_index_to_space
            tmp_cut_coords = idx_to_spc(0, 0, 0)
            
        elif component is const.Hx:
            field = fdtd.hx.real
            spc_to_idx = fdtd.space.space_to_hx_index
            idx_to_spc = fdtd.space.hx_index_to_space
            tmp_cut_coords = idx_to_spc(0, 1, 1)
            
        elif component is const.Hy:
            field = fdtd.hy.real
            spc_to_idx = fdtd.space.space_to_hy_index
            idx_to_spc = fdtd.space.hy_index_to_space
            tmp_cut_coords = idx_to_spc(1, 0, 1)
            
        elif component is const.Hz:
            field = fdtd.hz.real
            spc_to_idx = fdtd.space.space_to_hz_index
            idx_to_spc = fdtd.space.hz_index_to_space
            tmp_cut_coords = idx_to_spc(1, 1, 0)
            
        if axis is const.X:
            high_idx = [i - 1 for i in field.shape]
            high = idx_to_spc(*high_idx)
            self.extent = (low[2], high[2], high[1], low[1])
            
            cut_idx = spc_to_idx(cut, tmp_cut_coords[1], tmp_cut_coords[2])
            if in_range(cut_idx, field, component) is False:
                return None
            self.data = field[cut_idx[0], :, :]
            
            self.xlabel, self.ylabel = 'z', 'y'
            
        elif axis is const.Y:
            low = idx_to_spc(0, 0, 0)
            high_idx = [i - 1 for i in field.shape]
            high = idx_to_spc(*high_idx)
            self.extent = (low[2], high[2], high[0], low[0])
            
            cut_idx = spc_to_idx(tmp_cut_coords[0], cut, tmp_cut_coords[2])
            if in_range(cut_idx, field, component) is False:
                return None
            self.data = field[:, cut_idx[1], :]
            
            self.xlabel, self.ylabel = 'z', 'x'
            
        elif axis is const.Z:
            low = idx_to_spc(0, 0, 0)
            high_idx = [i - 1 for i in field.shape]
            high = idx_to_spc(*high_idx)
            self.extent = (low[1], high[1], high[0], low[0])
            
            cut_idx = spc_to_idx(tmp_cut_coords[0], tmp_cut_coords[1], cut)
            if in_range(cut_idx, field, component) is False:
                return None
            self.data = field[:, :, cut_idx[2]]
            
            self.xlabel, self.ylabel = 'y', 'x'
            
        else:
            msg = "axis must be gmes.constants.Directional."
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
        self.manager = new_figure_manager(self.fig_id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, aspect='auto', 
                            vmin=self.vrange[0], vmax=self.vrange[1])
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
                This should be one of the gmes.constants.Component.
            axis: Specify the normal axis to the show plane.
                This should be one of the gmes.constants.Directional.
            cut: A scalar value which specifies the cut position on the 
                axis. 
            vrange: Specify the colorbar range. a tuple of length two.
            title: title string of the figure.
            fig_id:

        """
        Thread.__init__(self)

        if component is const.Ex:
            material = fdtd.material_ex
            spc_to_idx = fdtd.space.space_to_ex_index
            idx_to_spc = fdtd.space.ex_index_to_space
            tmp_cut_coords = idx_to_spc(0, 0, 0)
            
        elif component is const.Ey:
            material = fdtd.material_ey
            spc_to_idx = fdtd.space.space_to_ey_index
            idx_to_spc = fdtd.space.ey_index_to_space
            tmp_cut_coords = idx_to_spc(0, 0, 0)
            
        elif component is const.Ez:
            material = fdtd.material_ez
            spc_to_idx = fdtd.space.space_to_ez_index
            idx_to_spc = fdtd.space.ez_index_to_space
            tmp_cut_coords = idx_to_spc(0, 0, 0)
        
        elif component is const.Hx:
            material = fdtd.material_hx
            spc_to_idx = fdtd.space.space_to_hx_index
            idx_to_spc = fdtd.space.hx_index_to_space
            tmp_cut_coords = idx_to_spc([i - 1 for i in material.shape])
            
        elif component is const.Hy:
            material = fdtd.material_hy
            spc_to_idx = fdtd.space.space_to_hy_index
            idx_to_spc = fdtd.space.hy_index_to_space
            tmp_cut_coords = idx_to_spc([i - 1 for i in material.shape])
            
        elif component is const.Hz:
            material = fdtd.material_hz
            spc_to_idx = fdtd.space.space_to_hz_index
            idx_to_spc = fdtd.space.hz_index_to_space
            tmp_cut_coords = idx_to_spc([i - 1 for i in material.shape])
            
        if axis is const.X:
            high_idx = [i - 1 for i in material.shape]
            high = idx_to_spc(*high_idx)
            self.extent = (low[2], high[2], high[1], low[1])
            
            cut_idx = spc_to_idx(cut, tmp_cut_coords[1], tmp_cut_coords[2])
            if in_range(cut_idx, material, component) is False:
                return None
            
            self.data = empty((material.shape[1], material.shape[2]), float) 
            if issubclass(component, const.Electric):          
                for idx in ndindex(*self.data.shape):
                    material_idx = cut_idx[0], idx[0], idx[1]
                    self.data[idx] = material[material_idx].epsilon
            elif issubclass(component, const.Magnetic): 
                for idx in ndindex(*self.data.shape):
                    material_idx = cut_idx[0], idx[0], idx[1]
                    self.data[idx] = material[material_idx].mu
                    
            self.xlabel, self.ylabel = 'z', 'y'
            
        elif axis is const.Y:
            low = idx_to_spc(0, 0, 0)
            high_idx = [i - 1 for i in material.shape]
            high = idx_to_spc(*high_idx)
            self.extent = (low[2], high[2], high[0], low[0])
            
            cut_idx = spc_to_idx(tmp_cut_coords[0], cut, tmp_cut_coords[2])
            if in_range(cut_idx, material, component) is False:
                return None
            
            self.data = empty((material.shape[0], material.shape[2]), float)
            if issubclass(component, const.Electric):                
                for idx in ndindex(self.data.shape):
                    material_idx = idx[0], cut_idx[1], idx[1]
                    self.data[idx] = material[material_idx].epsilon
            elif issubclass(component, const.Magnetic):
                for idx in ndindex(self.data.shape):
                    material_idx = idx[0], cut_idx[1], idx[1]
                    self.data[idx] = material[material_idx].mu
                    
            self.xlabel, self.ylabel = 'z', 'x'
            
        elif axis is const.Z:
            low = idx_to_spc(0, 0, 0)
            high_idx = [i - 1 for i in material.shape]
            high = idx_to_spc(*high_idx)
            self.extent = (low[1], high[1], high[0], low[0])
            
            cut_idx = spc_to_idx(tmp_cut_coords[0], tmp_cut_coords[1], cut)
            if in_range(cut_idx, material, component) is False:
                return None
            
            self.data = empty((material.shape[0], material.shape[1]), float)
            if issubclass(component, const.Electric):
                for idx in ndindex(self.data.shape):
                    material_idx = idx[0], idx[1], cut_idx[2]
                    self.data[idx] = material[material_idx].epsilon
            elif issubclass(component, const.Magnetic):
                for idx in ndindex(self.data.shape):
                    material_idx = idx[0], idx[1], cut_idx[2]
                    self.data[idx] = material[material_idx].mu
                    
            self.xlabel, self.ylabel = 'y', 'x'
            
        else:
            msg = "axis must be gmes.constants.Directional."
            raise ValueError(msg)

        self.window_title = 'GMES' + ' ' + str(fdtd.space.cart_comm.topo[2])

        if vrange is None:
            self.vrange = array((self.data.min(), self.data.max()))
        else:
            self.vrange = array(vrange, float)

        self.title = str(title)

        self.fig_id = fig_id
        
    def run(self):
        self.manager = new_figure_manager(self.fig_id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, aspect='auto', 
                            vmin=self.vrange[0], vmax=self.vrange[1], 
                            cmap=cm.gray)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        self.manager.canvas.figure.colorbar(self.im)
        self.manager.window.title(self.window_title)
        
        self.manager.canvas.draw()
        self.manager.show()
        self.manager.window.mainloop()
