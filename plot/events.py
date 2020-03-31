"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw


# ------------------------------------------------------------------------------


@traced
@logged
class IndexTracker(object):
  """
  Original scroller.
  
  """
  def __init__(self, ax, A, istep=1, **kwargs):
    """
    """
    from fullwavepy.generic.array import slice_array
    from fullwavepy.plot.generic import plot_array
    
    scoord = kw('scoord', 'y', kwargs)
    self.istep = istep
    self.imshow_kwargs = {'cmap': 'twilight'}
    self.ax = ax
    self.A = A
    a = slice_array(A, **kwargs)
    plot_array(a, **kwargs)
    #print('djkaf', plt.gcf().axes[0].images)
    self.im = ax.imshow(a.T, **self.imshow_kwargs)
    self.ind = kw('svalue', 0, kwargs)
    self.update()

    if scoord == 'x':
      svalue_max = A.shape[0]                           
    elif scoord == 'y':
      svalue_max = A.shape[1]
    elif scoord == 'z':
      svalue_max = A.shape[2]
    self.svalue_max = svalue_max


  def onscroll(self, event):
    if event.button == 'up':
      self.ind = (self.ind + self.istep) % self.svalue_max
    else:
      self.ind = (self.ind - self.istep) % self.svalue_max
    self.update()

  def update(self, **kwargs):
    from fullwavepy.generic.array import slice_array    
    kwargs['svalue'] = self.ind
    a = slice_array(self.A, **kwargs)
    self.im.set_data(a.T)
    self.ax.set_ylabel('slice %s' % self.ind)
    self.im.axes.figure.canvas.draw()
  
  
# ------------------------------------------------------------------------------
  
  
@traced
@logged
class IndexTrackerAll(object):
  """
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, fig, A, istep=1, **kwargs):
    """
    """
    from matplotlib.gridspec import GridSpec
    
    self.fig = fig
    self.A = A
    self.istep = istep
    
    self.i = kw('scoords', np.array(A.shape) // 2, kwargs)
    self.imax = A.shape

    self.imshow_kwargs = {'cmap': kw('cmap', 'twilight', kwargs),
                          'vmin': 1500, 'vmax': 6500}    
    self.line_kwargs = {'lw': 1, 'ls': '--', 'c': 'k'}
    gs = GridSpec(3, 1) #, left=0.05, right=0.48, wspace=0.05)
    
    self.ax = list(np.zeros(3))
    self.im = list(np.zeros(3))
    
    
    for j in range(len(self.ax)):
      self.ax[j] = self.fig.add_subplot(gs[j, :])
      a = self.A.take(indices=self.i[j], axis=j)
      self.im[j] = self.ax[j].imshow(a.T, **self.imshow_kwargs)
      
      from mpl_toolkits.axes_grid1 import make_axes_locatable
      divider = make_axes_locatable(self.ax[j])
      cax = divider.append_axes("right", size="5%", pad=0.05)      
      self.ax[j].figure.colorbar(self.im[j], cax=cax)
      #xy = [0,1,2]
      #xy.remove(j)
      #x1 = range(self.imax[xy[0]])
      #y1 = np.ones(len(x1)) * self.i[xy[1]]
      #print('creating horiz', xy[0], xy[1])
      #self.ax[jplot(x1, y1, **self.line_kwargs)
      #x2 = range(self.imax[xy[1]])
      #y2 = np.ones(len(x2)) * self.i[xy[0]]
      #self.ax[jplot(y2, x2, **self.line_kwargs)    
      #print('creating vert', xy[1], xy[0])
      #
      #print('lines of j ', j, self.ax[j].lines)    
      
      self.ax[j].set_title('slice %s' % self.i[j])
      if j == 0:
        self.ax[j].set_xlabel('cross-line node')
        self.ax[j].set_ylabel('depth node')
      elif j == 1:
        self.ax[j].set_xlabel('in-line node')
        self.ax[j].set_ylabel('depth node')
      elif j == 2:
        self.ax[j].set_xlabel('in-line node')
        self.ax[j].set_ylabel('cross-line node')        
        self.ax[j].invert_yaxis()
    
      self.ax[j].set_aspect('equal')

  # -----------------------------------------------------------------------------

  def onscroll(self, event):
    for j in range(len(self.i)):
      if event.button == 'up':
        self.i[j] = (self.i[j] + self.istep) % self.imax[j]
      else:
        self.i[j] = (self.i[j] - self.istep) % self.imax[j]
    self.update(event=event)

  # -----------------------------------------------------------------------------

  def update(self, **kwargs):
    event = kw('event', None, kwargs)
    
    for j in range(len(self.im)):
      if event.inaxes == self.im[j].axes:
        a = self.A.take(indices=self.i[j], axis=j)
        self.im[j].set_data(a.T)

        self.ax[j].set_title('slice %s' % self.i[j])
        self.im[j].axes.figure.canvas.draw()     
     
        #if j == 2:
          #lines = 0
        #else:
          #lines = 1
        
        # CHANGE LINES AT REMAINING AXES
        #remaning_axes = [0,1,2]
        #remaning_axes.remove(j)
        #for axes_index in remaning_axes:
        #  xy = [0,1,2]
        #  xy.remove(axes_index)
        #  #x1 = range(self.imax[xy[0]])
        #  #y1 = np.ones(len(x1)) * self.i[xy[1]]
        #  #self.ax[axes_index].lines[2].set_data(x1, y1)
        #  x2 = range(self.imax[xy[1]])
        #  y2 = np.ones(len(x2)) * self.i[xy[0]]
        #  self.ax[axes_index].lines[1].set_data(y2, x2)
        #break


# ------------------------------------------------------------------------------
        
    
@traced
@logged
class IterationScroller(object):
  """
  Scroll through subsequent iterations.
  
  """
  def __init__(self): # fig, iterfile_list, istep=1, **kwargs):
    """
    """
    from matplotlib.gridspec import GridSpec
    print('eee')
   #self.fig = fig
   #self.A = A
   #self.istep = istep
   #
   #self.i = kw('scoords', np.array(A.shape) // 2, kwargs)
   #self.imax = A.shape
   #
   #self.imshow_kwargs = {'cmap': kw('cmap', 'twilight', kwargs),
   #                      'vmin': 1500, 'vmax': 6500}    
   #self.line_kwargs = {'lw': 1, 'ls': '--', 'c': 'k'}
   #gs = GridSpec(3, 1) #, left=0.05, right=0.48, wspace=0.05)
   #
   #self.ax = list(np.zeros(3))
   #self.im = list(np.zeros(3))
   #
   #
   #for j in range(len(self.ax)):
   #  self.ax[j] = self.fig.add_subplot(gs[j, :])
   #  a = self.A.take(indices=self.i[j], axis=j)
   #  self.im[j] = self.ax[j].imshow(a.T, **self.imshow_kwargs)
   #  
   #  from mpl_toolkits.axes_grid1 import make_axes_locatable
   #  divider = make_axes_locatable(self.ax[j])
   #  cax = divider.append_axes("right", size="5%", pad=0.05)      
   #  self.ax[j].figure.colorbar(self.im[j], cax=cax)
   #  
   #  self.ax[j].set_title('slice %s' % self.i[j])
   #  if j == 0:
   #    self.ax[j].set_xlabel('cross-line node')
   #    self.ax[j].set_ylabel('depth node')
   #  elif j == 1:
   #    self.ax[j].set_xlabel('in-line node')
   #    self.ax[j].set_ylabel('depth node')
   #  elif j == 2:
   #    self.ax[j].set_xlabel('in-line node')
   #    self.ax[j].set_ylabel('cross-line node')        
   #    self.ax[j].invert_yaxis()
   #
   #  self.ax[j].set_aspect('equal')

  # -----------------------------------------------------------------------------

  def onscroll(self, event):
    for j in range(len(self.i)):
      if event.button == 'up':
        self.i[j] = (self.i[j] + self.istep) % self.imax[j]
      else:
        self.i[j] = (self.i[j] - self.istep) % self.imax[j]
    self.update(event=event)

  # -----------------------------------------------------------------------------

  def update(self, **kwargs):
    event = kw('event', None, kwargs)
    
    for j in range(len(self.im)):
      if event.inaxes == self.im[j].axes:
        a = self.A.take(indices=self.i[j], axis=j)
        self.im[j].set_data(a.T)

        self.ax[j].set_title('slice %s' % self.i[j])
        self.im[j].axes.figure.canvas.draw()     
     
        #if j == 2:
          #lines = 0
        #else:
          #lines = 1
        
        # CHANGE LINES AT REMAINING AXES
        #remaning_axes = [0,1,2]
        #remaning_axes.remove(j)
        #for axes_index in remaning_axes:
        #  xy = [0,1,2]
        #  xy.remove(axes_index)
        #  #x1 = range(self.imax[xy[0]])
        #  #y1 = np.ones(len(x1)) * self.i[xy[1]]
        #  #self.ax[axes_index].lines[2].set_data(x1, y1)
        #  x2 = range(self.imax[xy[1]])
        #  y2 = np.ones(len(x2)) * self.i[xy[0]]
        #  self.ax[axes_index].lines[1].set_data(y2, x2)
        #break


# ------------------------------------------------------------------------------
        
    