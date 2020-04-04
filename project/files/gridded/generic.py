"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists

from fullwavepy.ndat.arrays import Arr3d
from fullwavepy.ndat.points import Nodes
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.ioapi.fw3d import VtrFile
from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.project.lists.basic import ShotFileList, TimestepFileList


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
  def _shape(self, **kwargs):
    shape = kw('shape', (self.proj.nx1, self.proj.nx2, self.proj.nx3), kwargs)
    return shape 
 
  # ----------------------------------------------------------------------------- 
  
  #def _extent(self, **kwargs):
    #"""
    #"""
    #box = self.proj.box
    #extent = np.array([box[ :2], box[2:4], box[4: ]])
    ## CENTER VOXELS AT INTEGERS! ONLY IF unit='node'! FIXME
    #extent += 0.5 times dx possibly
    
    #self.__log.debug('extent: %s' % str(extent))    
    #return extent

  # -----------------------------------------------------------------------------
  
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
  
  #def read(self, **kwargs):
    #self.array = super().read(**kwargs)
    #self.__log.debug('self.array.extent %s' % str(self.array.extent))
    #self.array.extent = self._extent(**kwargs)
    
    #return self.array
 
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
  def _shape(self, **kwargs):
    shape = kw('shape', (self.proj.enx1, self.proj.enx2, self.proj.enx3), kwargs)
    return shape
  
  # -------------------------------------------------------------------------------

  def _extent(self, **kwargs):
    """
    
    """
    x1 = float(-self.proj.elef)
    x2 = float(x1 + self.proj.enx1)
    y1 = float(-self.proj.efro)
    y2 = float(y1 + self.proj.enx2)
    z1 = float(-self.proj.etop)
    z2 = float(z1 + self.proj.enx3)
    
    if self.proj.nx2 == 1:
      y1 = 1
      y2 = 1
    
    extent = np.array([[x1, x2], [y1, y2], [z1, z2]])
    # CENTER VOXELS AT INTEGERS! ONLY IF unit='node'! FIXME
    extent += 0.5
    self.__log.debug('extent: %s' % str(extent))
    return extent

  # -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------  


@traced
@logged
class InextFile(ExtendedGridFile):
  """
  Interior/exterior nodes of 
  extended grid.
  
  """
  def __init__(self, proj, path, **kwargs):
    self.name = proj.name + '-InextNodes.txt'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)
  
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    from fullwavepy.ioapi.generic import read_txt
    c = read_txt(self.fname)
    c = [[float(i) for i in j] for j in c]
    h = [int(i) for i in c[0]]
    d = np.array(c[1: ])
    self.array = Nodes(d.reshape(h + [4])[:,:,:,-1])
    self.array.extent = self._extent(**kwargs)
    return self.array
  
  # -----------------------------------------------------------------------------  
  
  
# -------------------------------------------------------------------------------


