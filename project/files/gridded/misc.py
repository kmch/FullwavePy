"""
(c) 2019-2020 Kajetan Chrapkiewicz.
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
    
    self.in_flag = 0
    self.acc_flag = -1
    self.ext_flag = -666
    
    super().__init__(proj, path, **kwargs) 
    
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    in_flag = 0                ! FLAG FOR INTERIOR NODE
    integer, parameter         :: acc_flag = -1             ! FLAG FOR NODE LYING EXACTLY (WITHIN ACCURACY) ON FS
    integer, parameter         :: ext_flag = int(err_value)   
        #better plots? mapp = {-666.0 :  1,
    #          -1.0 :  0,
    #           0.0 : -1,
    #       }
    """
    
    from fullwavepy.ioapi.generic import read_txt
    from fullwavepy.ndat.arrays import Arr3d
    

            
    
    c = read_txt(self.fname)
    #c = [[float(i) for i in j] for j in c]
    h = [int(i) for i in c[0]]
    d = [[float(i) for i in j] for j in c[1: ]]
    #d = np.array(c[1: ])
    
    d = np.array(d)
    
    self.array = Arr3d(d.reshape(h + [4])[:,:,:,-1])
    self._extent()
    self.array.extent = self.extent
    self.array.in_flag = self.in_flag
    self.array.acc_flag = self.acc_flag
    self.array.ext_flag = self.ext_flag
    
    self.__log.debug('self.array.extent %s' % str(self.array.extent))
    
    return self.array
  
  # -----------------------------------------------------------------------------  
  
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'magma_r', kwargs)
    super().plot(*args, **kwargs)

  # -----------------------------------------------------------------------------
  
  
# -------------------------------------------------------------------------------


