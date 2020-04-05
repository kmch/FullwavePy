"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.project.files.gridded.generic import ExtenGridProjFile


@traced
@logged
class InextFile(ExtenGridProjFile):
  """
  Interior/exterior nodes of 
  extended grid.

  Notes
  -----
  M1
  super().__init__(proj, path, **kwargs)
  
  M2 (explicit call of the grand...parent
  from fullwavepy.project.files.generic import ProjFile,
  ProjFile.__init__(self, proj, path, **kwargs) 

  
  """
  def __init__(self, proj, path, **kwargs):
    self.name = proj.name + '-InextNodes.txt'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs) 
    
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    from fullwavepy.ioapi.generic import read_txt
    from fullwavepy.ndat.arrays import Grid
    c = read_txt(self.fname)
    c = [[float(i) for i in j] for j in c]
    h = [int(i) for i in c[0]]
    d = np.array(c[1: ])
    self.array = Grid(d.reshape(h + [4])[:,:,:,-1])
    self._extent()
    self.array.extent = self.extent
    self.__log.debug('self.array.extent %s' % str(self.array.extent))
    return self.array
  
  # -----------------------------------------------------------------------------  
  
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'magma_r', kwargs)
    super().plot(*args, **kwargs)

  # -----------------------------------------------------------------------------
  
  
# -------------------------------------------------------------------------------


