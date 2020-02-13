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
from fullwavepy.project.files.generic import BinaryProjFile, ArrayProjFile
from fullwavepy.ioapi.fw3d import TtrFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.project.files.datalike.generic import DataFile, SynDataFile, ObsDataFile


# -------------------------------------------------------------------------------


@traced
@logged
class DataFileSgy(DataFile, SgyFile):
  """
  
  """

  # -----------------------------------------------------------------------------     
  
  def __init__(self, suffix, proj, path, **kwargs):
    """
    
    """  
    from fullwavepy.project.files.index import IndexFile
    from fullwavepy.project.files.templates import HedFile
    self.suffix = suffix
    self.ext = 'sgy'
    self.idx = IndexFile(suffix, proj, path, **kwargs) # SIMILAR TO TEMPLATE
    self.hed = HedFile(suffix, proj, path, **kwargs)
    super().__init__(proj, path, **kwargs)

  # -----------------------------------------------------------------------------   
  
  @timer
  def files_OLD(self, **kwargs):
    """
    Create hooks for 
    
    Notes
    -----
    Based on sources, receivers files.
    So if some files are missing that 
    will be raised as error when 
    trying to access such a file.
    (desired behaviour).
    
    we actually should store list of 
    shot ids as metadata...
    
    The following line:
    this_class = type(self)
    allows to create SynDataFile
    
    """
    from fullwavepy.ioapi.su import sugethw
    
    
    this_class = type(self) # CAN BE A CHILD CLASS
    self.__log.debug('this_class: ' + str(this_class))
    
    if not exists(self.fname):
      self.__log.warn(self.fname + ' not found. Returning.')
      return
    
    skey = self.proj.sgy.hw['sid']
    lkey = self.proj.sgy.hw['lid']
    
    sids = self.proj.env.csrefs_to_dump()
    if sids is None:
      sids = sugethw(self.fname, key=skey, unique_values=True, int_values=True, timer=True)
    
    if lkey is not None:
      fname = strip(self.fname) + '_' + str(skey) + str(sids[0]) + '.' + exten(self.fname)
      if exists(fname):
        self.__log.warn('Assuming each source gather contains the same lines')
        lids = sugethw(fname, key=lkey, unique_values=True, int_values=True, timer=True)
      else:
        lids = sugethw(self.fname, key=lkey, unique_values=True, int_values=True, timer=True)
    else:
      self.__log.warn('(Shot) line ID is not set in SEGY mapping. Files will not be split into lines.')
      lids = []
      
    self.gather = {} # IS IT A GOOD IDEA?
    self.line = {}
    self.__log.info('Creating instances of ' + str(this_class) + '...')
    for sid in sids:
      suffix = self.suffix + '_' + str(skey) + str(sid)
      self.gather[sid] = this_class(self.proj, self.path, 
                                    sid=sid, **kwargs)
      self.line[sid] = {}
      for lid in lids:
        suffix = str(self.suffix + '_' + str(skey) + str(sid) + '_' + 
                     str(lkey) + str(lid))

        self.line[sid][lid] = this_class(self.proj, self.path, 
                                         sid=sid, lid=lid, **kwargs)
 
  # -----------------------------------------------------------------------------   
  
  def split_OLD(self, sid=None, lid=None, **kwargs):
    """
    Split so that each file corresponds 
    to a single station and a single shot lines.
    
    Parameters
    ----------
    tracf - station
    ep - line
    
    Notes
    -----
    It collects file names
    using the list of project 
    shots and GetFiles().
    
    Only sgy.
    
    It overwrites all the old files. Otherwise 
    it would be mistakes-prone if old and new files
    were mixed.
    
    """
    from fullwavepy.ioapi.segy import split_sgy
    
    keys = kw('keys', [self.proj.sgy.hw['sid'],
                       self.proj.sgy.hw['lid']], kwargs)
    
    fnames = [self.fname]
    
    for key, value in zip(keys, [sid, lid]):
      for fname in fnames:
        nfnames = split_sgy(fname, key, value, **kwargs)
      
      fnames = list(nfnames)
      
    #self.files(**kwargs)
 
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# SPECIFIC FILES (UNIQUE IDs)
# -------------------------------------------------------------------------------


@traced
@logged
class SynDataFileSgy(DataFileSgy, SynDataFile):
  """
  Synthetic data.
  
  """

  # -----------------------------------------------------------------------------

  def __init__(self, proj, path, **kwargs):
    suffix = 'Synthetic'
    super().__init__(suffix, proj, path, **kwargs)
    self.hed = None # ACTUALLY IT IS NOT CREATED BY FULLWAVE
  
  # -----------------------------------------------------------------------------  
  
  def compare(self, **kwargs):
    """
    
    This assumefullwavepy.plotting lines.
    
    """
    outseis = self.proj.inp.outseis.line[self.sid][self.lid].fname
    kwargs['fname2'] = kw('fname2', outseis, kwargs)
    super().compare(**kwargs)
 
  # -----------------------------------------------------------------------------
  
  def compare_phase_all_shots(self, **kwargs):
    kwargs['cbar'] = kw('cbar', True, kwargs)
    gathers = self.gather
    for sid, synth in self.gather.items():
      for freq in [2, 3, 4]:
        outseis = self.proj.inp.outseis.gather[sid]

        title = str(sid) + '_' + str(freq) + 'Hz'
        plt.subplots(1,3, figsize=(20,10))
        plt.suptitle(title)
        
        plt.subplot(1,3,1)
        outseifullwavepy.plot_phase(freq=freq)
        plt.gca().set_aspect('equal')

        plt.subplot(1,3,2)
        syntfullwavepy.plot_phase(freq=freq)
        plt.gca().set_aspect('equal')

        plt.subplot(1,3,3)
        syntfullwavepy.plot_phase(freq=freq, subtract=outseis)
        plt.gca().set_aspect('equal')
        
        plt.savefig(self.path + self.pname + '-' + title)
        plt.close()
    
  # -----------------------------------------------------------------------------  
  
  def compare_all_lines(self, **kwargs):
    """
    
    """
    kwargs['cbar'] = kw('cbar', True, kwargs)
    
    lines = self.line
    for sid, val in lines.items():
      for lid in val:
        title = strip(lines[sid][lid].fname)
        plt.title(title)
        lines[sid][lid].compare(**kwargs)
        plt.savefig(title)
        plt.close()

  # -----------------------------------------------------------------------------  
  

# -------------------------------------------------------------------------------


@traced
@logged
class ObsDataFileSgy(DataFileSgy, ObsDataFile):
  """
  Observed (field) data.
  
  """

  # -----------------------------------------------------------------------------

  def __init__(self, proj, path, **kwargs):
    suffix = 'Observed'
    #self.suffix = suffix
    super().__init__(suffix, proj, path, **kwargs)  
 
  # -----------------------------------------------------------------------------    
  
  def _get_phase(self, freq, **kwargs):
    """
    
    This assumefullwavepy.plotting gathers.
    
    Notes
    -----
    First breaks must be extracted from synthetic data in all cases!
    Noise in observed.
    
    """
    file2 = self.proj.out.synth.gather[self.sid]
    phase = super()._get_phase(freq, file2, **kwargs)
    
    return phase

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------


@traced
@logged
class SignatureFileSgy(DataFileSgy):
  """
  
  """
  
  # -----------------------------------------------------------------------------     
  
  def __init__(self, proj, path, **kwargs):
    suffix = 'Signature'
    super().__init__(suffix, proj, path, **kwargs)

  # ----------------------------------------------------------------------------- 
  

# -------------------------------------------------------------------------------


@traced
@logged
class RawSignFile(DataFileSgy):
  """
  
  """  

  # ----------------------------------------------------------------------------- 

  def __init__(self, proj, path, **kwargs):
    """
    
    """
    suffix = 'RawSign'
    super().__init__(suffix, proj, path, **kwargs)
    self.idx = None # 
    self.hed = None ##
    
  # -----------------------------------------------------------------------------  

  def create(self, wavelet, **kwargs):
    """
    wavelet : str / array
      if string, it will use one of the 
      built-in wavelets (now only 'ricker'),
      if array, it will just save it to file.
    
    """
    from fullwavepy.ioapi.segy import array2sgy
    
    super().create(**kwargs)
    
    if type(wavelet) == type(np.array([])):
      array = wavelet
    
    elif isinstance(wavelet, str):
      if 'fpeak' in kwargs:
        fpeak = kwargs['fpeak']
      else:
        raise TypeError('You need to provide fpeak argument.')
      
      if wavelet == 'ricker':
        from fullwavepy.generic.math import ricker
        array = ricker(fpeak, self.proj.ns, self.proj.dt)
      
      elif wavelet == 'gauss':
        from fullwavepy.generic.math import gauss
        self.__log.warn('Actually fpeak is WRONG, it is tpeak')
        sigma = kw('sigma', 1, kwargs)
        self.__log.info('Assuming sigma=%s' % sigma)
        array = gauss(np.arange(self.proj.ns), mu=fpeak, sigma=sigma)
      
      elif wavelet == 'sine':
        t = np.arange(self.proj.ns) * self.proj.dt
        period = 1 / fpeak
        nodes = period / self.proj.dt
        self.__log.info('period: %s seconds, %s nodes' % (period, nodes))
        
        array = np.sin(2 * np.pi * fpeak * t)
        
      else:
        raise ValueError('Unknown wavelet type: ' + wavelet)
    
    else:
      raise TypeError('wavelet arg. must be either a string or an array')

    array2sgy(self.fname, array, self.proj.dt)
    
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    Convert into a 1D array.
    
    """
    from fullwavepy.generic.array import Arr1d
    a = super().read(**kwargs)
    
    #assert a.shape[:2] == (1,1) # WE ASSUME shape=(1,1,nsamp)
    if not isinstance(a, Arr1d):
      self.array = Arr1d(a[0,0,:])
    
    return self.array


# ------------------------------------------------------------------------------- 


