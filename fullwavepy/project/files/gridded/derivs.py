"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.gridded.generic import ExtenGridProjFile


@traced
@logged
class GradFile(ExtenGridProjFile):
  def __init__(self, suffix, proj, path, **kwargs):
    """
    """
    self.name = proj.name + '-' + suffix + '.vtr'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)  
  def plot(self, **kwargs):
    kwargs['cmap'] = kwargs.get('cmap', 'RdBu')
    kwargs['center_cmap'] = kwargs.get('center_cmap', True)
    super().plot(**kwargs) 


# ------------------------------------------------------------------------------- 


@traced
@logged
class PrecFile(ExtenGridProjFile):
  def __init__(self, suffix, proj, path, **kwargs):
    """
    """
    self.name = proj.name + '-' + suffix + '.vtr'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs) 
  def plot(self, **kwargs):
    kwargs['cmap'] = kwargs.get('cmap', 'hot')
    kwargs['center_cmap'] = kwargs.get('center_cmap', False)
    super().plot(**kwargs) 

# ------------------------------------------------------------------------------- 

