"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists

from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.ioapi.fw3d import VtrFile
from fullwavepy.project.lists.basic import ShotFileList, TimestepFileList


# -------------------------------------------------------------------------------


@traced
@logged
class GridFile(ArrayProjFile):
  """
  File containing a property 
  defined on the model grid.
  
  Notes
  ----- 
  In particular shape of the grid 
  is assumed to be known a priori.
  
  """
  def create(self, array, **kwargs): # FIXME MOVE TO VTR
    from fullwavepy.ioapi.fw3d import save_vtr
    super().create(**kwargs)
    save_vtr(array , strip(self.fname) + '.vtr')

  # -----------------------------------------------------------------------------
  
  def dupl(self, source, cmd='cp',  *args, **kwargs): 
    """
    Standard file-duplication + resizing 
    to conform to the project-box.
    
    """
    super().dupl(source, cmd=cmd, **kwargs)
    
    if 'geom' in dir(self.proj):
      #self.resize(*args, **kwargs)
      self.__log.warn('Resize disabled until debugged')
    else:
      self.__log.warn('proj.geom not yet set.' + self.fname + 
                      ' will not be resized')

  # ----------------------------------------------------------------------------- 
  
  def read(self, **kwargs):
    self.array = super().read(**kwargs)
    self.array.extent = np.array([self.proj.box[ :2], 
                                  self.proj.box[2:4],
                                  self.proj.box[4: ]])
    return self.array
  
  # -----------------------------------------------------------------------------
  
  def resize(self, **kwargs):
    pass
  
  # -----------------------------------------------------------------------------
  
  def plot_3slices(self, *args, **kwargs):
    self.array = self.read(**kwargs)
    return self.array.plot_3slices(*args, **kwargs)
  
  # -----------------------------------------------------------------------------  
  
  def plot(self, **kwargs):
    kwargs['center_cmap'] = kw('center_cmap', False, kwargs)
    super().plot(**kwargs)

  # ----------------------------------------------------------------------------- 
  

# -------------------------------------------------------------------------------


@traced
@logged
class ExtendedGridFile(GridFile):
  """
  File defined on an extended grid 
  (i.e. with extra nodes).
  
  """
  def read(self, **kwargs):
    """
  #  #kwargs['scoord'] = None
  #  A = super().read(**kwargs)
  #  #if self.proj.dim == '3D':
  #  #  A = A[self.proj.eleft:-self.proj.eright, self.proj.efront:-self.proj.eback, self.proj.etop:-self.proj.ebot]
  #  #else:
  #  #  A = A[self.proj.eleft:-self.proj.eright, :, self.proj.etop:-self.proj.ebot]
  #  return A    
    """
    kwargs['shape'] = kw('shape', (self.proj.enx1, self.proj.enx2, self.proj.enx3), kwargs)
    A = super().read(**kwargs)
    return A


# -------------------------------------------------------------------------------  

