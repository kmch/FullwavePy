"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import widgets
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.system import bash
from fullwavepy.ioapi.fw3d import VtrFile
from fullwavepy.ndat.manifs import SurfZ
from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.project.files.gridded.generic import GridProjFile, ExtenGridProjFile


@traced
@logged
class SurfZFile(ArrayProjFile, VtrFile):
  """
  File storing a model boundary (e.g. free surface) 
  or interface.
  
  Notes
  -----
  At the moment it assumes the surface to be 
  a function z = f(x,y), i.e. it doesn't allow
  overhangs etc. 

  This implies it has a shape (nx, ny, 1) because
  there is only one z value per (x,y).
  
  Such an array can be stored in a vtr file and 
  that is what we do here.
  
  """
  def __init__(self, suffix, proj, path, **kwargs):
    """
    """
    self.name = proj.name + '-' + suffix + '.vtr'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)

  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    from fullwavepy.ioapi.generic import read_any
    A = np.array(read_any(self.fname))
    assert len(A.shape) == 3
    assert A.shape[-1] == 1
    Az = A[:,:,0]
    self.array = SurfZ(Az)
    return self.array
    
  # -----------------------------------------------------------------------------
  
  #def plot(self, *args, **kwargs):
    #self.plot2d(*args, **kwargs)

  ## -----------------------------------------------------------------------------
  
  #def plot2d(self, **kwargs):
    #self.array = self.read(**kwargs)
    ##x1, x2 = self.proj.box[ :2]
    ##dx = self.proj.dx
    ##x = np.arange(x1, x2+dx, dx)
    ##x = np.arange(0, len(self.array))
    ##z = self.array[:,0,0]
    ##plt.plot(x, z)
    #return self.array

  ## -----------------------------------------------------------------------------
  #def plot3d(self, **kwargs):
    #from mpl_toolkits.mplot3d import Axes3D
    #ax = plt.gca(projection='3d')
    #ax.plot_surface(self.array)


# -------------------------------------------------------------------------------  


@traced
@logged
class TopoFile(SurfZFile):
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
class FsFile(SurfZFile, GridProjFile):
  """
  Free surface
  
  """
  def __init__(self, proj, path, **kwargs):
    suffix = 'FreeSurf'
    super().__init__(suffix, proj, path, **kwargs)

  def immerse(self, ibmcode='default', **kwargs):
    kwargs['truncate'] = kwargs.get('truncate', 1e4)
    if ibmcode == 'default':
      ibmcode = self.proj.ibmcode
    else:
      ibmcode = ibmcode

    ibmcode.run(\
      self.proj.name, \
      self.proj.i.path, **kwargs)
    
  def read(self, **kwargs):
    self.array = SurfZFile.read(self, **kwargs)
    self.array.extent = GridProjFile._extent(self, **kwargs)[ :-1] # SKIP Z-EXTENT
    return self.array
  
  def run(self, *args, **kwargs):
    self.immerse(*args, **kwargs)

  def run_old(self, **kwargs):
    """
    set log_lvl(n<=10) to see the output messages.
    """
    for fname in [self.proj.i.ghb.fname, self.proj.i.ght.fname]:
      self.__log.info('Deleting %s...' % fname)
      cmd = 'rm %s' % fname
      o, e = bash(cmd, **kwargs)

    exe = self.proj.exe['fsprep']
    path = exe[ :-len('fsprep')]
    self.__log.info('Making %s...' % exe)
    o, e = bash('make -C %s' % path)

    cmd = exe + " " + self.proj.name
    self.__log.info('Running %s...' % cmd)
    o, e = bash(cmd, path=self.proj.inp.path, **kwargs)
    self.__log.info(o) 
    if len(e) > 0:
      self.__log.warning(e)
  

# -------------------------------------------------------------------------------  


@traced
@logged
class ExtendedFsFile(SurfZFile, ExtenGridProjFile):
  """
  Extended-grid free surface
  
  """  
  def __init__(self, proj, path, suffix='FreeSurf_exten', **kwargs):
    super().__init__(suffix, proj, path, **kwargs)
  
  def read(self, **kwargs):
    self.array = SurfZFile.read(self, **kwargs) - self.proj.etop # COORECT THE Z-VALUE!
    self.array.extent = ExtenGridProjFile._extent(self, **kwargs)[ :-1] # SKIP Z-EXTENT
    return self.array


# -------------------------------------------------------------------------------


@traced
@logged
class InterpolFsFile(ExtendedFsFile):
  """
  Free surface interpolated on extended-grid.
  
  """    
  def __init__(self, proj, path, **kwargs):
    suffix = 'FreeSurf_exten_interp'
    super().__init__(proj, path, suffix, **kwargs)


# -------------------------------------------------------------------------------

