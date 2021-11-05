"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import TextProjFile
from fullwavepy.seismic.srcrec import SRs, PointSR, Sources, Receivers


@traced
@logged
class SRFile(TextProjFile):
  """
  Sources/receivers files.
  
  Notes
  -----
  Creating a metaclass for pgy/geo seems 
  like an overkill.
  
  """  
  def create(self, dictio, **kwargs):
    """
    dictio : dictionary
      Dictionary with items:
        ID: [x, y, z]
    """
    from fullwavepy.ioapi.fw3d import save_pgy
    #from fullwavepy.ioapi.segy import save_geo
    
    if self.proj.io == 'fw3d':
      save_pgy(self.fname, dictio, **kwargs)
  
  # -----------------------------------------------------------------------------
  
  def read(self, extend=False, unit='node', **kwargs):
    """
    extend: if True, add extra nodes
    """
    #if not hasattr(self, 'd'): #IT'S A TINY FILE SO IT'S SAFER TO READ IT EVERYTIME
    io = self.proj.io
    if io == 'sgy':
      from fullwavepy.ioapi.segy import read_geo
      if 'dx' in dir(self.proj): kwargs['dx'] = self.proj.dx
      kwargs['x0'] = self.proj.box[0]
      kwargs['y0'] = self.proj.box[2]
      kwargs['z0'] = self.proj.box[4]
      sr = read_geo(self.fname, **kwargs)
      
    elif io == 'fw3d':
      from fullwavepy.ioapi.fw3d import read_pgy
      sr = read_pgy(self.fname, **kwargs)
    
    else:
      raise ValueError('Unknown io: ' + io)
    
    # self.__log.debug('self.d = SRs(sr)')
    # sr_SRs = SRs(sr)
    # print('sr_SRs', sr_SRs)
    #sr_SRs # it's not dict any more but we keep it for backward compatibility for now
    # return self.d
    return sr #SRs(sr)
  
  # -----------------------------------------------------------------------------  
  
  def spread_factors(self, *args, **kwargs):
    """
    """
    sr = self.read(**kwargs)
    sr.spread_factors(*args, **kwargs)

  # -----------------------------------------------------------------------------  
  
  def plot(self, *args, **kwargs): # GENERALIZE FOR ASCII PROJ FILE?
    return self.read(*args, **kwargs).plot(*args, **kwargs)
    
  # ----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class SourcesFile(SRFile):
  """

  """  
  def __init__(self, proj, path, **kwargs):
    """
    """
    super().__init__(proj, path, **kwargs)
    if proj.io == 'sgy':
      suffix = '-Sources.geo'
    elif proj.io == 'fw3d':
      suffix = '-PointSources.pgy'
    self.name = proj.name + suffix
    self.fname = self.path + self.name # self.path => trailing /  

  # -----------------------------------------------------------------------------

  def read(self, *args, **kwargs):
    return Sources(super().read(*args, **kwargs))
    
  
# -------------------------------------------------------------------------------


@traced
@logged
class ReceiversFile(SRFile):
  """

  """  
  def __init__(self, proj, path, **kwargs):
    """
    """
    super().__init__(proj, path, **kwargs)
    if proj.io == 'sgy':
      suffix = '-Receivers.geo'
    elif proj.io == 'fw3d':
      suffix = '-PointReceivers.pgy'
    
    self.name = proj.name + suffix
    self.fname = self.path + self.name # self.path => trailing /  
   
  # -----------------------------------------------------------------------------

  def read(self, *args, **kwargs):
    return Receivers(super().read(*args, **kwargs))
  
  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------

