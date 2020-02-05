"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.ioapi.fw3d import TtrFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.plot.generic import plot


# -------------------------------------------------------------------------------


@traced
@logged
class DataFile(ArrayProjFile):
  """
  Seismic-data file.
  
  Notes
  -----
  Data-processing methods can be applied 
  to the whole 'lump' (file before splitting).
  
  """
  
  # -----------------------------------------------------------------------------     
  
  def __init__(self, proj, path, **kwargs):
    """
    
    """
    # READ SHOT AND LINE ID IF PRESENT (IMPORTANT FEATURE USED BY compare ETC.)
    self.sid = kw('sid', None, kwargs)
    self.lid = kw('lid', None, kwargs)
    super().__init__(proj, path, **kwargs)

    self.name = proj.name + '-' + self.suffix + '.' + self.ext
    self.fname = path + self.name    
    self.__log.debug('self.fname: ' + self.fname)
  
  # -----------------------------------------------------------------------------   


# -------------------------------------------------------------------------------


@traced
@logged
class SynDataFile(DataFile):
  """
  Synthetic data.
  
  """
  
  # -----------------------------------------------------------------------------  
  
  def _get_first_breaks(self, fraction=0.01, **kwargs):
    """
    
    """
    from fullwavepy.signal.phase import first_breaks
    A = self.read(**kwargs)
    fb = first_breaks(A, fraction=fraction)
    fb = np.ravel(fb)
    
    fname = strip(self.fname) + '_firstbreaks.txt'
    with open(fname, 'w') as f:
      for pick in fb:
        f.write(str(pick) + '\n')
    
    return fb

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


@traced
@logged
class ObsDataFile(DataFile):
  """
  Observed data.
  
  """

  def __init__(self, proj, path, create=True, raw=True, filt=True, mute=True, **kwargs):
    """
    Creates subobjects of the same type but after some processing 
    (raw, filtering, muting).
    
    """
    super().__init__(proj, path, **kwargs)
    
    if create and raw:
      self.raw = self.__class__(proj, path, create=False, raw=True, filt=0, mute=0)
    elif raw: # THIS IS CALLED INSIDE self.raw.__init__!
      self.suffix += '_raw'

    if create and filt:
      self.filtered = self.__class__(proj, path, create=False, filt=True, raw=0, mute=0)
      self.fil = self.filtered
    elif filt:
      self.suffix += '_filtered'
      
    if create and mute:
      self.muted = self.__class__(proj, path, create=False, mute=True, raw=0, filt=0)
      self.mut = self.muted
    elif mute:
      self.suffix += '_muted'      

    self.name = proj.name + '-' + self.suffix + '.' + self.ext
    self.fname = path + self.name    
    self.__log.debug('self.fname: ' + self.fname)  
 
  # -----------------------------------------------------------------------------  
  
  def filt(self, **kwargs):
    self.__log.info('Creating a backup of the data:\n{}\nas:\n{}'.format(self.fname, self.raw.fname))
    self.raw.dupl(self.fname)
    #print(self.fil.fname)
    #kwargs['overwrite'] = kw('overwrite', True, kwargs)
    #try:
    #  self.filtered.dupl(self.raw.fname)
    #  self.filtered.filt(**kwargs)
    #  self.dupl(self.filtered.fname)
    #except AttributeError: # RECURSION
    #  super().filt(**kwargs)
  
  # -----------------------------------------------------------------------------  
  
  #def mute(self, **kwargs):

  # -----------------------------------------------------------------------------  
  
  def _get_first_breaks(self, **kwargs):
    """
    
    """
    if self.sid is None and self.lid is None:
      synth_file = self.proj.out.synth
    elif self.lid is None:
      synth_file = self.proj.out.synth.gather[self.sid]
    else:
      synth_file = self.proj.out.synth.line[self.sid][self.lid]

    fb = synth_file._get_first_breaks(**kwargs)
    return fb  

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------

