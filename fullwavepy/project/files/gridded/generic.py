"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import ArrayProjFile


@traced
@logged
class GridProjFile(ArrayProjFile):
  """
  It only sets the correct extent for imshow plots.
  
  The first cell should be centered at 1 
  (grid start with node no. 1)
  
  """
  def _extent(self, unit='m', **kwargs):
    """
    Nodes are numbered as 1, 2, ...
    
    Notes
    -----

    what's below is outdated:

    We CAN'T put in __init__ as these files are initalized 
    before params like proj.nx1 are known (they need input 
    files to be initialized already).
    
    VERY IMPORTANT
    x2, y2 and z2 seem NOT physical,
    e.g. if nx = 10, (x1,x2) one could think the correct 
    extent in X is (1,10), but instead here is (1,11).
    This is the ONLY way it is plot correctly both in 2D and 1D.
    
    """
    if unit == 'm':
      self.extent = np.array(self.proj.box).reshape(3,2)
    else:
      node1 = 1 # plot_image WILL SUBTRACT 0.5 TO CENTER AT INTEGERS AGAIN
      x1 = node1
      x2 = x1 + self.proj.nx1 - 1
      y1 = node1
      y2 = y1 + self.proj.nx2 - 1
      z1 = node1
      z2 = z1 + self.proj.nx3 - 1
      self.extent = np.array([[x1, x2], [y1, y2], [z1, z2]])
    return self.extent
  
  # -----------------------------------------------------------------------------  
  
  def read(self, *args, **kwargs):
    extent = self._extent(**kwargs)
    kwargs['shape'] = kw('shape', self.proj.dims, kwargs)
    self.array = super().read(*args, **kwargs)
    self.array.extent = extent
    self.__log.debug('self.array.extent %s' % str(self.array.extent))
    return self.array

  # -----------------------------------------------------------------------------

  def prep(self, arr, **kwargs):
    """
    Parameters
    ----------
    arr : arrau.a3d.Arr3d
        Or derived.
    """
    x1, x2, y1, y2, z1, z2 = self.proj.box 
    arr = arr.extract(extent=[[x1,x2],[y1,y2],[z1,z2]]).arr
    super().prep(arr, **kwargs)

  # def resize(self, *args, **kwargs):
  #   pass

  # -----------------------------------------------------------------------------
  
  def widg_plot_3slices(self):
    """
    Return widget-kwargs for @interact used in the notebook.

    """
    from ipywidgets import BoundedIntText as BIT
    
    p = self.proj
    nx, ny, nz = p.dims
    kwargs = dict(x=BIT(value=nx//2, min=0, max=nx-1),
                  y=BIT(value=ny//2, min=0, max=ny-1),
                  z=BIT(value=nz//2, min=0, max=nz-1))
    return kwargs
   
  def wp3s(self, *args, **kwargs):
    return self.widg_plot_3slices(*args, **kwargs)
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ExtenGridProjFile(GridProjFile):
  """
  """
  def _extent(self, unit='m', **kwargs):
    """
    """
    super()._extent(unit, **kwargs)
    [[x1, x2], [y1, y2], [z1, z2]] = self.extent

    if unit == 'm':
      # FIXME: test if it actually works as expected
      x1 = x1 - self.proj.elef * self.proj.dx
      y1 = y1 - self.proj.efro * self.proj.dx
      z1 = z1 - self.proj.etop * self.proj.dx
      
      x2 = x1 + (self.proj.enx1 - 1) * self.proj.dx
      y2 = y1 + (self.proj.enx2 - 1) * self.proj.dx
      z2 = z1 + (self.proj.enx3 - 1) * self.proj.dx
    else:
      x1 -= self.proj.elef
      y1 -= self.proj.efro
      z1 -= self.proj.etop
      
      x2 = x1 + self.proj.enx1
      y2 = y1 + self.proj.enx2
      z2 = z1 + self.proj.enx3
      
    if self.proj.nx2 == 1: # back to value from super()._extent
      y1, y2 = self.extent[1] 
      
    self.extent = np.array([[x1, x2], [y1, y2], [z1, z2]])
    return self.extent

  # -----------------------------------------------------------------------------  
  
  def read(self, *args, **kwargs):
    extent = self._extent(**kwargs)
    kwargs['shape'] = self.proj.edims
    self.array = super().read(*args, **kwargs)
    self.array.extent = extent
    return self.array

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


