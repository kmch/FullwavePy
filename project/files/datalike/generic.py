"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, exten, strip, path_leave
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.ioapi.generic import CsvFile, read_txt
from fullwavepy.ioapi.fw3d import TtrFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.plot.generic import plot


# -------------------------------------------------------------------------------


@traced
@logged
class MetaDataFile(CsvFile):
  """
  
  csv is MUCH faster for Panda's read/write than json.
  Not to mention, we can choose columns.
  
  """
  def __init__(self, datafile, **kwargs):
    
    self.fname = strip(datafile.fname) + '_METADATA.csv'
    self.name = path_leave(self.fname)
    self.proj = datafile.proj
    self.path = datafile.path
    
  def read(self, overwrite=True, **kwargs):
    """
  
    Notes
    -----
    Overwrite=True by default because otherwise plots are not 
    updated even though they are supposed (e.g. you are passing 
    a different fname). They will be correct (updated) only
    if you delete self.array variable, e.g. by restarting the 
    notebook kernel.
    Disable overwrite only for PERFORMANCE (e.g. interactive plot)
    when the array remains unchanged unlike other (e.g. plotting)
    parameters.
    
    """
    if (not hasattr(self, 'df')) or overwrite:
      self.__log.warn('{}.df does not exist and will be read.'.format(type(self)))
      self.df = super().read(**kwargs)
    return self.df 


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

    self.meta = MetaDataFile(self, **kwargs)
  
  # -----------------------------------------------------------------------------   
  

# -------------------------------------------------------------------------------


@traced
@logged
class SynDataFile(DataFile):
  """
  Synthetic data.
  
  """
  
  # -----------------------------------------------------------------------------  
  
  def get_fbreaks(self, fraction=0.01, from_file=None, **kwargs):
    """
    
    """
    from fullwavepy.signal.phase import first_breaks
    A = self.read(**kwargs)
    self.fb = first_breaks(A, fraction=fraction)
    self.fb = np.ravel(self.fb)
    
    fname = strip(self.fname) + '_firstbreaks.txt'
    with open(fname, 'w') as f:
      for pick in self.fb:
        f.write(str(pick) + '\n')
    
    return self.fb

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
  
  #def filt(self, *args, **kwargs):
    #self.raw.dupl(self.fname)
    #self.__log.info('{} backuped as {}'.format(self.fname, self.raw.fname))
    #fname_fil = super().filt(*args, **kwargs)
    #self.__log.warn('Overwriting {} with {}'.format(fname_fil, fname_out))
    #o, e = bash('mv ' + fname_fil + ' ' + self.fname)
    #super().filt(*args, **kwargs)
    
    
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
  
  def get_fbreaks(self, syn_file, **kwargs):
    """
    
    """
    #if self.sid is None and self.lid is None:
    #  syn_file = self.proj.out.syn
    #elif self.lid is None:
    #  syn_file = self.proj.out.syn.gather[self.sid]
    #else:
    #  syn_file = self.proj.out.syn.line[self.sid][self.lid]

    self.fb = syn_file.get_fbreaks(**kwargs)
    return self.fb  

  # -----------------------------------------------------------------------------
    
  def read_first_breaks(self, fname, **kwargs):
    self.fb = read_txt(fname)
    self.fb = [float(i[0]) for i in self.fb]  
    return self.fb

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------

