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
  
  def read(self, **kwargs):
    from fullwavepy.generic.array import WigglyData
    self.array = WigglyData(super().read(**kwargs))
    return self.array

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class SynDataFile(DataFile):
  """
  Synthetic data.
  
  """
  def __init__(self, *args, **kwargs):
    self.phase = {}
    super().__init__(*args, **kwargs) 
  
  # -----------------------------------------------------------------------------
  
  def get_fbreaks(self, fraction=0.01, overwrite=True, **kwargs):
    """
    
    """
    from fullwavepy.signal.phase import first_breaks
    
    fb_fname = strip(self.fname) + '_firstbreaks.txt'
    
    if (not exists(fb_fname)) or overwrite:
      A = self.read(**kwargs)
      self.fb = first_breaks(A, fraction=fraction)
      self.fb = np.ravel(self.fb)

      with open(fb_fname, 'w') as f:
        for pick in self.fb:
          f.write(str(pick) + '\n')

    else:
      self.fb = []
      with open(fb_fname, 'r') as f:
        for line in f:
          self.fb.append(float(line))
    
    return self.fb

  # -----------------------------------------------------------------------------
  
  def _get_phase(self, freq, overwrite=True, **kwargs):
    """
    
    Notes
    -----
    First breaks must be extracted from START-MOD
    synthetic data in all cases!
    Noise in observed.

    # I CHANGED BACK AFTER JO SAID THE ORIGINAL VERSION WAS OK
    # WE TAKE SYNTHETICS FROM THE START MOD FOR ALL ITERATIONS!
    #Bsyn, Bobs, Bdif = self.proj.out.dumpcomp.it[1][self.sid].read(**kwargs)
    #picks = first_breaks(Bsyn, **kwargs)
    
    Actually we should assume (ntraces, 1, nsamps) shape...
    
    """
    if (not freq in self.phase) or overwrite:
      from fullwavepy.signal.phase import first_breaks, extract_phase, wrap_phase
      from fullwavepy.generic.math import rms    
      self.__log.info('Getting phase info from ' + self.fname)
      
      self.read(**kwargs)
      self.read_header(**kwargs)
      
      if not hasattr(self, 'fb'):
        self.fb = first_breaks(self.array, **kwargs)
      
      ph_syn = np.ravel(extract_phase(self.array, self.fb, self.proj.dt, freq, **kwargs))
      #ph_obs = extract_phase(self.obs, self.fb, self.proj.dt, freq, **kwargs)
      #ph_dif = ph_syn - ph_obs
      
      #ph_dif = np.array([[wrap_phase(i) for i in j] for j in ph_dif])
      #ph_syn, ph_obs, ph_dif = [np.ravel(i) for i in [ph_syn, ph_obs, ph_dif]]
      
      #self.rms_value = rms(ph_dif)
      #self.__log.info('RMS of wrapped phase-differences: ' + 
                      #str(self.rms_value))
      
      self.head['phase syn (%s Hz)' % freq] = ph_syn
      #self.head['phase obs (%s Hz)' % freq] = ph_obs
      #self.head['phase dif (%s Hz)' % freq] = ph_dif
      
      self.phase[freq] = {'syn': ph_syn,}
                          #'obs': ph_obs,
                          #'dif': ph_dif}
    return self.phase[freq]

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
    (raw, filtering, muting). THIS IS NOT NECESSARY ACTUALLY, DELETE 
    AFTER YOU ARCHIVE THIS APPROACH FOR FUTURE USE
    
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
  
  def get_fbreaks(self, syn_file, **kwargs):
    """
    First breaks are found for corresponding synthetics.
    """
    self.__log.info('Getting first breaks from {}...'.format(syn_file.fname))
    self.fb = syn_file.get_fbreaks(**kwargs)
    return self.fb  

  # -----------------------------------------------------------------------------
  
  def filt(self, **kwargs):
    self.__log.info('Filtering {}...'.format(self.fname))
    super().filt(**kwargs)
    self.fil.dupl(self.fname) 

  # -----------------------------------------------------------------------------
  
  def mute(self, syn_file, **kwargs):
    fbreaks = self.get_fbreaks(syn_file, **kwargs)
    self.__log.info('Muting {}...'.format(self.fname))
    super().mute(fbreaks, **kwargs)

  # -----------------------------------------------------------------------------

  def process(self, filt_kwargs, mute_kwargs, **kwargs):
    self.__log.info('Processing {}...'.format(self.fname))
    self.filt(**filt_kwargs, **kwargs)
    self.mute(**mute_kwargs, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------

