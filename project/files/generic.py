"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, exten
from fullwavepy.generic.system import bash, exists
from fullwavepy.ioapi.generic import AsciiFile, BinaryFile, ArrayFile


# -------------------------------------------------------------------------------  


@traced
@logged
class ProjFile(object):
  """
  Project file.
  
  """
  
  # -----------------------------------------------------------------------------

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
  
  def prepare(self, *args, **kwargs):
    """
    Prepare means either duplicate 
    (as long as a source file is provided, dupl=...) 
    or create from scratch.
    
    Notes
    -----
    The flow is a bit convoluted:
    1. child.prepare calls this.
    2. this calls child.dupl/child.create!
    3. those can again call their parents
    
    """
    
    dupl = kw('dupl', None, kwargs)
    
    if dupl is not None:
      self.__log.debug('Using dupl='+dupl)
      del_kw('dupl', kwargs)
      self.dupl(dupl, *args, **kwargs)
    else:
      self.create(*args, **kwargs)
    
  def prep(self, *args, **kwargs):
    self.prepare(*args, **kwargs)
    
  # -----------------------------------------------------------------------------

  def create(self, *args, **kwargs):
    self.__log.info('Creating ' + self.fname + '...')    
    
  # -----------------------------------------------------------------------------

  def dupl(self, source, cmd='cp', **kwargs):
    """
    Duplicate (cp/mv/ln) a file.
    
    Notes
    -----
    
    """
    from fullwavepy.generic.parse import exten
    from fullwavepy.generic.system import duplicate
    
    destination = self.fname
    
    if exten(source) != exten(destination):
      raise IOError('Extension of source and destination must be the same.')    
    
    duplicate(source, destination, cmd)
    
  # -----------------------------------------------------------------------------
  
  def rm(self, cmd='trash', backup=True, **kwargs):
    """
    """
    if backup:
      bckp = self.fname + '_bckp'
      self.__log.debug('Creating a backup ' + bckp)
      o, e = bash('cp ' + self.fname + ' ' + bckp)
    
    self.__log.info('Removing ' + self.fname + '...')
    o, e = bash(cmd + ' ' + self.fname)
    
  # -----------------------------------------------------------------------------

  def compare(self, other_file, fig, gs=None, **kwargs):
    arr1 = self.read(**kwargs)
    arr2 = other_file.read(**kwargs)
    arr1.compare(arr2, fig, gs, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------  


@traced
@logged
class AsciiProjFile(ProjFile, AsciiFile):
  """
  File storing some meta data,
  not directlfullwavepy.plottable as an array.
  
  """
  
  # -----------------------------------------------------------------------------
  
  def prepare(self, **kwargs):
    super().prepare(**kwargs)
    
    cat = kw('cat', True, kwargs)
    if cat:
      self.cat()
      
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class BinaryProjFile(ProjFile, BinaryFile):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class ArrayProjFile(ProjFile, ArrayFile):
  pass


# -------------------------------------------------------------------------------

