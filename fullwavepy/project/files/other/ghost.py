"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import BinaryProjFile, TextProjFile


@traced
@logged
class SRsDataFileTxt(TextProjFile):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class SourcesDataFileTxt(SRsDataFileTxt):
  def __init__(self, proj, path, **kwargs):
    super().__init__(proj, path, **kwargs)
    self.name = self.pname + '-SourcesData.txt'
    self.fname = self.path + self.name


# -------------------------------------------------------------------------------


@traced
@logged
class ReceiversDataFileTxt(SRsDataFileTxt):
  def __init__(self, proj, path, **kwargs):
    super().__init__(proj, path, **kwargs)
    self.name = self.pname + '-ReceiversData.txt'
    self.fname = self.path + self.name


# -------------------------------------------------------------------------------

@traced
@logged
class GhostDataFileBin(BinaryProjFile):
  """
  """
  def __init__(self, proj, path, **kwargs):
    super().__init__(proj, path, **kwargs)
    self.name = self.pname + '-GhostData.bin'
    self.fname = self.path + self.name


# -------------------------------------------------------------------------------


@traced
@logged
class GhostDataFileTxt(TextProjFile):
  """
  """
  def __init__(self, proj, path, **kwargs):
    super().__init__(proj, path, **kwargs)
    self.name = self.pname + '-GhostData.txt'
    self.fname = self.path + self.name
    
  # -----------------------------------------------------------------------------   
  
  def read(self, **kwargs):
    """
    Read all the information contained 
    in a GhostData.txt file.
    
    Return
    -------
    ghosts : list
    intersects : list
    ficts : list
    auxs : list
    weights : list
    
    """
    # from fullwavepy.immerse.points import Ghosts
    
    ct = super().read(**kwargs)
    
    header = ct[0]
    data = ct[1: ]
    fict_no, auxs_no = [int(i) for i in header]
    
    ghosts_no, intersects_no = 1, 1
    ghosts, intersects, ficts, auxs, weights = [], [], [], [], []
    
    j = 0
    while j < len(data):
      ghosts.append(data[j])
      j += 1
      intersects.append(data[j])
      j += 1
      
      fcs = []
      for f in range(fict_no):
        fcs.append(data[j])
        j += 1
      ficts.append(fcs)
      
      fcs_auxs = []
      for f in range(fict_no):
        f_auxs = []
        for ax in range(auxs_no):
          f_auxs_x = []
          for ay in range(auxs_no):
            f_auxs_x.append(data[j])
            j += 1
          f_auxs.append(f_auxs_x)
        fcs_auxs.append(f_auxs)  
      auxs.append(fcs_auxs)
      
      fcs_weights = []
      for f in range(fict_no):
        f_weights = []
        for w in range(2 * auxs_no):
          f_weights.append(data[j])
          j += 1
        fcs_weights.append(f_weights)
      weights.append(fcs_weights)
     
    ghosts = [[int(i) for i in j] for j in ghosts]
    intersects = [[float(i) for i in j] for j in intersects]
    ficts = [[[float(i) for i in j] for j in k] for k  in ficts]
    auxs = [[[[[float(i) for i in j] for j in k] for k in l] for l in m] for m in auxs]  
    weights = [[[[float(i) for i in j] for j in k] for k in l] for l in weights]
    
    self.ghosts = np.array(ghosts) #Ghosts(ghosts)
    self.isects = intersects
    self.ficts = ficts
    self.auxs = auxs
    self.weights = weights
    
    #return ghosts, intersects, ficts, auxs, weights    

  # ----------------------------------------------------------------------------- 
  

# -------------------------------------------------------------------------------
