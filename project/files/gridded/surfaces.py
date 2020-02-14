"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import widgets
from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.ioapi.fw3d import VtrFile
from fullwavepy.project.lists.basic import ShotFileList, TimestepFileList
from fullwavepy.project.files.gridded.models import ModelFileVtr


# -------------------------------------------------------------------------------


@traced
@logged
class SurfaceFile(ModelFileVtr):
  """
  Free surface or a model
  interface.
  
  """
  pass


# -------------------------------------------------------------------------------  


@traced
@logged
class TopographyFile(SurfaceFile):
  """
  Topography of the topmost rock layer:
  seafloor and/or land surface.
  
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, proj, path, **kwargs):
    suffix = 'Topography'
    super().__init__(suffix, proj, path, **kwargs)
  
  # -----------------------------------------------------------------------------  
  
  def plot(self, array=None, **kwargs):
    from fullwavepy.plot.misc import plot_square 
    
    kwargs['cbar'] = kw('cbar', True, kwargs)
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    kwargs['shade'] = kw('shade', True, kwargs)
    kwargs['cmap'] = kw('cmap', [], kwargs)
    kwargs['extent'] = kw('extent', [-8e4, 8e4, 4e4, -4e4], kwargs)
    pad = kw('pad', 10*self.proj.dx, kwargs)
    full = kw('full', False, kwargs)
    xlim = kw('xlim', (self.proj.box[0]-pad, self.proj.box[1]+pad), kwargs)
    ylim = kw('ylim', (self.proj.box[2]-pad, self.proj.box[3]+pad), kwargs)
    if full:
      xlim = None
      ylim = None
    
    plot_square(self.proj.box[0], self.proj.box[1], 
                self.proj.box[2], self.proj.box[3]) 
      #kwargs['slice_at'] = 'z'
      
      #array.plot(**kwargs)  
      #plt.gca().set_aspect('equal')
      #plt.xlim(xlim)
      #plt.ylim(ylim)
      #plt.gca().invert_yaxis()
      
    
    #plot_square(self.proj.box[0], self.proj.box[1], 
                #self.proj.box[2], self.proj.box[3])
    

    ## SOURCES AND RECEIVERS NOTE: LOOP OVER S, R
    #path = '/home/kmc3817/heavy_PhD/meta_data/'
    #sources, receivers = SR_Read_All_PROTEUS(**kwargs)    
    #
    #sx, sy = [], []
    #for key in sources:
    #
    #  x = sources[key][0]
    #  y = sources[key][1]
    #  sx.append(x)
    #  sy.append(y)   
    #
    #plt.scatter(sx, sy, s=1*s_factor, c='lightgray')
    #
    #sx, sy = [], []
    #for key in receivers:
    #  x = receivers[key][0]
    #  y = receivers[key][1]
    #  
    #  sx.append(x)
    #  sy.append(y)
    #  plt.annotate(key, (x + shift_labels, y + shift_labels), 
    #               clip_on=True) # NOTE: clip_on IS REQUIRED
    #plt.scatter(sx, sy, s=10*s_factor, c='orange')    
    #
    #plt.xlim(xlim)
    #plt.ylim(ylim)    
    return None

  # -----------------------------------------------------------------------------
    
  def plotly(self, array=None, fig=None, **kwargs):
    """
    x0 etc. are now hardwired for Santorini bathymetry.
    """
    import plotly.graph_objects as go
    
    stride = kw('stride', 20, kwargs)
    zmin = kw('zmin', -30, kwargs)
    zmax = kw('zmax', 30, kwargs)
    ncontours = kw('ncontours', 40, kwargs)
    cmap = kw('cmap', 'Earth', kwargs)
    box = kw('box', True, kwargs)
    srcs = kw('scrs', True, kwargs)
    recs = kw('recs', True, kwargs)
    
    
    if array is None:
      array = self.read(**kwargs)
    
    if fig is None:
      fig = go.Figure()

    fig.add_trace(go.Contour(z=array[::stride, ::stride, 0].T,
                             zmin=zmin, zmax=zmax,
                             colorscale=cmap, ncontours=ncontours,
                             x0=-8e4, dx=50*stride, 
                             y0=-4e4, dy=50*stride, 
                             ))
    if box:
      self.proj.pbox.plotly(fig, color='red')
    if recs:
      self.proj.i.r.plotly(fig, size=3, **kwargs)
    if srcs:
      self.proj.i.s.plotly(fig, size=5, **kwargs)
    
    return fig
  
  # -----------------------------------------------------------------------------



# -------------------------------------------------------------------------------  


@traced
@logged
class FsFile(SurfaceFile):
  """
  Free surface
  
  """
  def __init__(self, proj, path, **kwargs):
    suffix = 'FreeSurf'
    super().__init__(suffix, proj, path, **kwargs)
  
  def run(self, **kwargs):
    cmd = self.proj.exe['fsprep'] + " " + self.proj.name
    o, e = bash(cmd, path=self.proj.inp.path, **kwargs)

   
  
# -------------------------------------------------------------------------------  

