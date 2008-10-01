#!/usr/bin/env python

try:
    import psyco
    psyco.profile()
    from psyco.classes import *
except:
    pass

from sys import modules
from threading import Thread

from numpy import array
from pylab import get_current_fig_manager, new_figure_manager
from pylab import figure, show

if not 'matplotlib.backends' in modules:
    import matplotlib 
    matplotlib.use('TkAgg')


class ShowLine(Thread):
    """Animated 1-D on-time display. 
    
    """
    def __init__(self, xdata, ydata, yrange, time_step, xlabel='', ylabel='', title='', window_title='GMES', msecs=2500, fig_id=None):
        Thread.__init__(self)
        self.xdata, self.ydata = xdata, ydata
        self.yrange = array(yrange, float)
        self.msecs= msecs
        self.xlabel, self.ylabel = xlabel, ylabel
        self.title = title
        self.window_title = window_title
        self.time_step = time_step
        self.note_form = 'time: %f'
        if fig_id == None:
            self.id = 0
        else:
            self.id = fig_id
        
    def animate(self):
        #self.manager.canvas.restore_region(self.background)
        self.plt.set_ydata(self.ydata)
        self.time_note.set_text(self.note_form % self.time_step.t)
        self.manager.canvas.draw()
        #self.ax.draw()
        #self.ax.draw_artist(self.plt)
        #self.manager.canvas.figure.draw_animated()
        #self.manager.canvas.blit(self.ax.bbox)
        if self.manager.window is not None:
            self.manager.window.after(self.msecs, self.animate)

    def run(self):
        self.manager = new_figure_manager(self.id)
        self.ax = self.manager.canvas.figure.add_subplot(111)
        self.plt, = self.ax.plot(self.xdata, self.ydata)#, animated=True)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        #self.ax.set_xlim(self.xdata[0], self.xdata[1])
        #self.ax.set_ylim(self.yrange[0], self.yrange[1])
        #self.ax.autoscale_view(tight=False, scalex=False, scaley=True) 
        self.ax.axis((self.xdata[0], self.xdata[-1], self.yrange[0], self.yrange[1]))
        self.ax.grid(True)
        self.ax.set_title(self.title)
        #self.manager.canvas.figure.colorbar(self.plt)
        self.time_note = self.manager.canvas.figure.text(.6, .92, self.note_form % self.time_step.t, animated=True)
        self.manager.window.title(self.window_title)
        
        #self.manager.canvas.draw()
        #self.background = self.manager.canvas.copy_from_bbox(self.ax.bbox)
        
        self.animate()
        self.manager.show()
        self.manager.window.mainloop()
        
class Snapshot(Thread):
    """A snapshot of 2-D display. 
    
    """
    def __init__(self, data, extent, range, time_step, xlabel='', ylabel='', title='', window_title='GMES', fig_id=None):
        Thread.__init__(self)
        self.data = data
        self.range = array(range, float)
        self.extent = extent
        self.xlabel, self.ylabel = xlabel, ylabel
        self.title = title
        self.window_title = window_title
        self.time_step = time_step
        self.note_form = 'time: %f'
        if fig_id == None:
            self.id = 0
        else:
            self.id = fig_id
    
    def run(self):
        self.manager = new_figure_manager(self.id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, aspect='auto', vmin=self.range[0], vmax=self.range[1])
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        self.manager.canvas.figure.colorbar(self.im)
        self.time_note = self.manager.canvas.figure.text(.6, .92, self.note_form % self.time_step.t)
        self.manager.window.title(self.window_title)
        self.time_note.set_text(self.note_form % self.time_step.t)
        
        self.manager.canvas.draw()       
        self.manager.show()
        self.manager.window.mainloop()
        
        
class ShowPlane(Thread):
    """Animated 2-D on-time display. 
    
    """
    def __init__(self, data, extent, range, time_step, xlabel='', ylabel='', title='', window_title='GMES', msecs=2500, fig_id=None):
        Thread.__init__(self)
        self.data = data
        self.range = array(range, float)
        self.msecs = msecs
        self.extent = extent
        self.xlabel, self.ylabel = xlabel, ylabel
        self.title = title
        self.window_title = window_title
        self.time_step = time_step
        self.note_form = 'time: %f'
        if fig_id == None:
            self.id = 0
        else:
            self.id = fig_id
        
    def animate(self):
        self.im.set_data(self.data)
        self.time_note.set_text(self.note_form % self.time_step.t)
        self.manager.canvas.draw()
        if self.manager.window is not None:
            self.manager.window.after(self.msecs, self.animate)
        
    def run(self):
        self.manager = new_figure_manager(self.id)
        ax = self.manager.canvas.figure.add_subplot(111)
        self.im = ax.imshow(self.data, extent=self.extent, aspect='auto', vmin=self.range[0], vmax=self.range[1])
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        self.manager.canvas.figure.colorbar(self.im)
        self.time_note = self.manager.canvas.figure.text(.6, .92, self.note_form % self.time_step.t, animated=True)
        self.manager.window.title(self.window_title)
        
        self.animate()
        self.manager.show()
        self.manager.window.mainloop()
