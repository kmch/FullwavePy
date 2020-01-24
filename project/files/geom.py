"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import AsciiProjFile


# -------------------------------------------------------------------------------


@traced
@logged
class SRFile(AsciiProjFile):
  """
  Sources/receivers files.
  
  Notes
  -----
  Creating a metaclass for pgy/geo seems 
  like an overkill.
  
  """  

  # -----------------------------------------------------------------------------
  
  def create(self, dictio, **kwargs):
    """
    dictio : dictionary
      Dictionary with items:
        ID: [x, y, z]
    """
    from fullwavepy.ioapi.fw3d import save_pgy
    #from fullwavepy.ioapi.segy import save_geo
    
    if self.proj.io == 'fw3d':
      save_pgy(self.fname, dictio, **kwargs)
  
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    
    """
    #if not hasattr(self, 'd'): #IT'S A TINY FILE SO IT'S SAFER TO READ IT EVERYTIME
    io = self.proj.io
    if io == 'sgy':
      from fullwavepy.ioapi.segy import read_geo
      if 'dx' in dir(self.proj): kwargs['dx'] = self.proj.dx
      kwargs['x0'] = self.proj.box[0]
      kwargs['y0'] = self.proj.box[2]
      kwargs['z0'] = self.proj.box[4]
      sr = read_geo(self.fname, **kwargs)
      
    elif io == 'fw3d':
      from fullwavepy.ioapi.fw3d import read_pgy
      sr = read_pgy(self.fname, **kwargs)
    
    else:
      raise ValueError('Unknown io: ' + io)
    
    self.d = sr
    
    return self.d
  
  # -----------------------------------------------------------------------------  
  
  def plot_slice(self, scoord, ax=None, **kwargs):
    """
    Add slice later.
    """
    from fullwavepy.plot.oned import slice_points
    annotate = kw('annotate', False, kwargs)
    annoffset = kw('annoffset', 0, kwargs)
    alpha = kw('alpha', 0.7, kwargs)
    marker = kw('marker', '.', kwargs)
    markersize = kw('markersize', 5, kwargs)
    markeredgecolor = kw('markeredgecolor', 'k', kwargs)
    markerfacecolor = kw('markerfacecolor', 'none', kwargs) # EMPTY MARKERS
    
    if ax is None:
      ax = plt.gca()
    
    d = self.read(**kwargs)
        
    points2d = slice_points(d.values(), scoord)
    for i, key in enumerate(d.keys()):
      x, y = points2d[i]
      if annotate: # NOTE: clip_on IS REQUIRED
        ax.annotate(key, (x+annoffset, y+annoffset), clip_on=True) 
    
    ax.plot([i[0] for i in points2d], [i[1] for i in points2d], '.',
            alpha=alpha, marker=marker, markersize=markersize, markeredgecolor=markeredgecolor,
            markerfacecolor=markerfacecolor)
  
  # -----------------------------------------------------------------------------   
  
  def plot_3slices(self, fig, **kwargs):
    d = self.read(**dict(kwargs, unit='node'))
    
    s3 = kw('slice', 'y', kwargs) #FIXME: THIS MUST BE MERGED WITH arr3d
    s1, s2 = [i for i in ['x', 'y', 'z'] if i != s3]
    s = [s1, s2, s3]
    
    for i in range(3):
      self.plot_slice(s[i], fig.axes[i])

  # -----------------------------------------------------------------------------
  
  def plot(self, fig=None, **kwargs): # FIXME
    """
    Add slice later.
    """
    if fig is None:
      fig = plt.figure()
      gs = fig.add_gridspec(2,2)
      fig.add_subplot(gs[0,0])
      fig.add_subplot(gs[0,1])
      fig.add_subplot(gs[1,:])
    
    self.plot_3slices(fig, **kwargs)
    
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


@traced
@logged
class SourcesFile(SRFile):
  """

  """  
  
  # -----------------------------------------------------------------------------

  def __init__(self, proj, path, **kwargs):
    """
    
    """
    super().__init__(proj, path, **kwargs)
    if proj.io == 'sgy':
      suffix = '-Sources.geo'
    elif proj.io == 'fw3d':
      suffix = '-PointSources.pgy'
    self.name = proj.name + suffix
    self.fname = self.path + self.name # self.path => trailing /  

  # -----------------------------------------------------------------------------  
  
  def plot_slice(self, *args, **kwargs):
    kwargs['marker'] = kw('marker', '*', kwargs)
    kwargs['markersize'] = kw('markersize', 10, kwargs)
    kwargs['markeredgecolor'] = kw('markeredgecolor', 'k', kwargs)
    kwargs['markerfacecolor'] = kw('markerfacecolor', 'w', kwargs)
    super().plot_slice(*args, **kwargs)

# -------------------------------------------------------------------------------


@traced
@logged
class ReceiversFile(SRFile):
  """

  """  
  
  # -----------------------------------------------------------------------------

  def __init__(self, proj, path, **kwargs):
    """
    
    """
    super().__init__(proj, path, **kwargs)
    if proj.io == 'sgy':
      suffix = '-Receivers.geo'
    elif proj.io == 'fw3d':
      suffix = '-PointReceivers.pgy'
    
    self.name = proj.name + suffix
    self.fname = self.path + self.name # self.path => trailing /  
    
  # -----------------------------------------------------------------------------  
  
  def plot_slice(self, *args, **kwargs):
    kwargs['alpha'] = kw('alpha', 1, kwargs)
    kwargs['markers'] = kw('marker', '.', kwargs)
    kwargs['markersize'] = kw('markersize', 2, kwargs)
    kwargs['markeredgecolor'] = kw('markeredgecolor', 'none', kwargs)
    kwargs['markerfacecolor'] = kw('markerfacecolor', 'k', kwargs)
    super().plot_slice(*args, **kwargs)    
    
  # -----------------------------------------------------------------------------  

  def plot(self, **kwargs):
    kwargs['annotate'] = False
    kwargs['s'] = 1e-2
    kwargs['c'] = 'gray'
    kwargs['alpha'] = .0001
    super().plot(**kwargs)
  
  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------

