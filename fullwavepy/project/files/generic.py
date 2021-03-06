"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.ioapi.generic import File, TextFile, BinaryFile, ArrayFile


# =associate with project
@traced
@logged
class ProjFile(File):
  """
  Project file.
  
  """
  def __init__(self, proj, path, **kwargs):
    """
    
    Notes
    -----
    Name is set in a child class.
    
    """    
    if path[-1] != '/':
      path += '/'
    self.pname = proj.name    
    self.path = path
    self.proj = proj

  # -----------------------------------------------------------------------------
  
  def _is_compatible(self, **kwargs):
    raise NotImplementedError('Only in child classes')

  # -----------------------------------------------------------------------------
   
  def _make_compatible(self, **kwargs):
    raise NotImplementedError('Only in child classes')

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------  


# FIXME: MOVE TO text
@traced
@logged
class TextProjFile(ProjFile, TextFile):
  """
  File storing some meta data,
  not directlfullwavepy.plottable as an array.
  
  """
  def prep(self, *args, **kwargs):
    super().prep(*args, **kwargs)
    
    cat = kw('cat', True, kwargs)
    if cat:
      self.cat()
      
  # -----------------------------------------------------------------------------
  
  
# -------------------------------------------------------------------------------


# FIXME: DEL
@traced
@logged
class BinaryProjFile(ProjFile, BinaryFile):
  pass


# -------------------------------------------------------------------------------


# =create=save_vtr(array) 
@traced
@logged
class ArrayProjFile(ProjFile, ArrayFile):
  """
  We keep it as a generic interface because not all 
  project's array files are defined on a grid 
  (e.g. data files).
  
  """
  def create(self, array, **kwargs):
    from fullwavepy.ioapi.fw3d import save_vtr
    save_vtr(array, strip(self.fname) + '.vtr')

  # -----------------------------------------------------------------------------

  # def prep()


# -------------------------------------------------------------------------------

