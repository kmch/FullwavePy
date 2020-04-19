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
  def _extent(self, **kwargs):
    """
    Nodes are numbered as 1, 2, ...
    
    Notes
    -----
    We CAN'T put in __init__ as these files are initalized 
    before params like proj.nx1 are known (they need input 
    files to be initialized already).
    
    VERY IMPORTANT
    x2, y2 and z2 seem NOT physical,
    e.g. if nx = 10, (x1,x2) one could think the correct 
    extent in X is (1,10), but instead here is (1,11).
    This is the ONLY way it is plot correctly both in 2D and 1D.
    
    """
    node1 = 1 # plot_image WILL SUBTRACT 0.5 TO CENTER AT INTEGERS AGAIN
    x1 = node1
    x2 = x1 + self.proj.nx1
    y1 = node1
    y2 = y1 + self.proj.nx2 
    z1 = node1
    z2 = z1 + self.proj.nx3
    self.extent = np.array([[x1, x2], [y1, y2], [z1, z2]])
    return self.extent
  
  # -----------------------------------------------------------------------------  
  
  def read(self, *args, **kwargs):
    extent = self._extent(**kwargs)
    self.array = super().read(*args, **kwargs)
    self.array.extent = extent
    return self.array

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ExtenGridProjFile(GridProjFile):
  """
  """
  def _extent(self, **kwargs):
    super()._extent(**kwargs)
    [[x1, x2], [y1, y2], [z1, z2]] = self.extent
    
    x1 -= self.proj.elef
    y1 -= self.proj.efro
    z1 -= self.proj.etop
    
    x2 = x1 + self.proj.enx1
    y2 = y1 + self.proj.enx2
    z2 = z1 + self.proj.enx3
    
    if self.proj.nx2 == 1:
      y1 = 1
      y2 = 1
    
    self.extent = np.array([[x1, x2], [y1, y2], [z1, z2]])
    return self.extent


# -------------------------------------------------------------------------------


