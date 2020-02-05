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
class DataFileTtr(DataFile, TtrFile):
  """
  
  """

  # -----------------------------------------------------------------------------     
  
  def __init__(self, suffix, proj, path, **kwargs):
    """
    
    """  
    self.suffix = suffix
    self.ext = 'ttr'
    #self.name = proj.name + '-' + suffix + '.ttr'
    #self.fname = path + self.name
    super().__init__(proj, path, **kwargs)

  # ----------------------------------------------------------------------------- 

  def split(self, **kwargs):
    self.__log.warn('Not implemented. Skipping.')

  # -----------------------------------------------------------------------------
  
  def files(self, **kwargs):
    self.__log.warn('Not implemented. Skipping.')

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class SynDataFileTtr(DataFileTtr, SynDataFile):
  pass
  

# -------------------------------------------------------------------------------


@traced
@logged
class ObsDataFileTtr(DataFileTtr, ObsDataFile):
  """
  """
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class ObsRawDataFileTtr(DataFileTtr):
  """
  """
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class SignatureFileTtr(DataFileTtr):
  """
  
  """
  
  # -----------------------------------------------------------------------------     
  
  def __init__(self, proj, path, **kwargs):
    suffix = 'SourceSig-Time'
    super().__init__(suffix, proj, path, **kwargs)

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


@traced
@logged
class DumpDataFile(DataFileTtr): # NOT USED?
  """
  DUMPDAT format, i.e.
    1. wavelets
    2. observed traces
    3. modelled traces
    
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, suffix, proj, path, sid, **kwargs):
    self.sid = sid
    super().__init__(suffix, proj, path, **kwargs)
  
  # -----------------------------------------------------------------------------      
  
  def read(self, **kwargs):
    #from fullwavepy.ioapi.fw3d import read_ttr
    #A = read_ttr(self.fname)
    #A = A[:, 0, :]
    #return A
    A = super().read(**kwargs)
    
    # READ SOURCE IDs
    #nsrcs = len(self.proj.inp.s.read(unit='m'))
    #self.__log.info('Ommiting first ' + str(nsrcs) +
    #                 '  trace(s) as source-wavelet(s)')
    #A = A[nsrcs: ]
    self.__log.info('Ommiting first trace (the source wavelet)')
    A = A[1: ]
    
    imid = int(len(A) / 2)
    self.__log.info('Ommiting next ' + str(imid) +
                    '  trace(s) as observed data')    
    Aobs = A[ :imid]
    Asyn = A[imid: ]
    
    return Asyn

  # -----------------------------------------------------------------------------    


# -------------------------------------------------------------------------------


@traced
@logged
class DumpCompareFile(DataFileTtr):
  """
  DUMPCompare format, i.e.
    1. Traces of 1st kind
    3. Traces of 2nd kind
    
  """

  # ----------------------------------------------------------------------------- 

  def __init__(self, suffix, proj, path, it, sid, **kwargs):
    super().__init__(suffix, proj, path, **kwargs)
    self.it = it
    self.sid = sid
    self.syn = SynDataFileTtr(suffix+'_syn', proj, path)
    #self.obs = ObsDataFileTtr(suffix+'_obs', proj, path)
    #self.dif = SynDataFileTtr(suffix+'_dif', proj, path)
    
  # -----------------------------------------------------------------------------      
  
  def _get_rids(self, **kwargs): # FIXME: MERGE WITH BELOW
    from fullwavepy.ioapi.generic import read_txt
 
    srcs = self.proj.inp.s.read()
    recs = self.proj.inp.r.read()
    try:
      c = read_txt(self.proj.inp.path + self.proj.name + '-Observed.hed')
    except FileNotFoundError as err:
      self.__log.warn('Searching for Template.hed because: {}'.format(err))
      c = read_txt(self.proj.inp.path + self.proj.name + '-Template.hed')
    c = c[2: ]
    rids = []
    for l in c:
      if int(l[1]) == self.sid:
        rids.append(int(l[2]))
    return rids
  
  def _get_sr_coords(self, **kwargs):
    """
    Get HORIZONTAL (x,y) coordinates of sources and 
    receivers.
    
    """
    from fullwavepy.generic.parse import strip
    from fullwavepy.ioapi.generic import read_txt
    
    srcs = self.proj.inp.s.read()
    recs = self.proj.inp.r.read()
    
    try:
      c = read_txt(self.proj.inp.path + self.proj.name + '-Observed.hed')
    except FileNotFoundError as err:
      self.__log.warn('Searching for Template.hed because: {}'.format(err))
      c = read_txt(self.proj.inp.path + self.proj.name + '-Template.hed')
    c = c[2: ]
    hed = {} # IT'S A BIT REDUNDANT SINCE THERE'S ONLY ONE SOURCE PER FILE
    for sid in sorted(srcs.keys()):
      hed[sid] = []
      for l in c:
        if int(l[1]) == sid:
          rid = int(l[2])
          hed[sid].append(recs[rid])
    
    srcs = [srcs[self.sid]]
    recs = hed[self.sid] ## AS ABOVE
    #self.__log.info('Assuming that file ' + self.fname + ' contains all receivers')

    return srcs, recs
 
  # -----------------------------------------------------------------------------    
  
  def read(self, **kwargs):
    kwargs['scoord'] = None
    self.__log.info('Reading ' + self.fname + '...')
    A = super().read(**kwargs)
    isep  = int(len(A) / 3) # THIS GIVES EQUAL LENGTHS, CHECKED
    
    Asyn = A[ :isep]
    Aobs = A[isep:-isep]
    Adif = A[-isep: ] # SYN - OBS
    self.__log.debug('lengths of Asyn, Aobs, Adif', len(Asyn), len(Aobs), len(Adif))
    
    self.syn.array = Asyn
    #self.obs.array = Aobs
    
    return Asyn, Aobs, Adif
  
  # -----------------------------------------------------------------------------    
  
  def split(self, **kwargs):
    pass
  
  
  def plot(self, **kwargs):
    #from fullwavepy.plot.generic import plot, compare
    Asyn, Aobs, Adif = self.read(**kwargs)
    #plt.subplots(1,3, figsize=[15,5])
    #plt.subplot(1,3,1)
    #plot(Asyn)
    #plt.subplot(1,3,2)
    #plot(Aobs)
    #plt.subplot(1,3,3)
    #plot(Adif)
    #compare(Asyn, Aobs, **kwargs)
    #fig = plt.figure(figsize=(14,8))
    Aobs.plot_slice() #(Asyn, fig)

  # -----------------------------------------------------------------------------  
  
  def _get_first_breaks(self, *args, **kwargs):
    """
    """
    from fullwavepy.signal.phase import first_breaks    
    Asyn, Aobs, Adif = self.read(**kwargs)
    self.fb = np.ravel(first_breaks(Asyn, *args, **kwargs))
    return self.fb
  
  # -----------------------------------------------------------------------------  
  
  def _get_phase(self, freq, **kwargs):
    """
    
    Notes
    -----
    First breaks must be extracted from START-MOD
    synthetic data in all cases!
    Noise in observed.
    
    Actually we should assume (ntraces, 1, nsamps) shape...
    
    """
    from fullwavepy.signal.phase import first_breaks, extract_phase, wrap_phase
    from fullwavepy.generic.math import rms
    
    self.__log.info('Getting phase info from ' + self.fname)
    
    kwargs['scoord'] = None # DON'T SLICE
    Asyn, Aobs, Adif = self.read(**kwargs)
    
    # I CHANGED BACK AFTER JO SAID THE ORIGINAL VERSION WAS OK
    # WE TAKE SYNTHETICS FROM THE START MOD FOR ALL ITERATIONS!
    #Bsyn, Bobs, Bdif = self.proj.out.dumpcomp.it[1][self.sid].read(**kwargs)
    #picks = first_breaks(Bsyn, **kwargs)
    
    picks = first_breaks(Asyn, **kwargs) # NOTE
    ph_syn = extract_phase(Asyn, picks, self.proj.dt, freq, **kwargs)
    ph_obs = extract_phase(Aobs, picks, self.proj.dt, freq, **kwargs)
    ph_dif = ph_syn - ph_obs
    ph_dif = np.array([[wrap_phase(i) for i in j] for j in ph_dif])
    
    self.rms_value = rms(ph_dif)
    self.__log.info('RMS of wrapped phase-differences: ' + 
                    str(self.rms_value))
    
    self.phase = {}
    self.phase[freq] = [ph_syn, ph_obs, ph_dif]
    
    return ph_syn, ph_obs, ph_dif
    
  # -----------------------------------------------------------------------------   
  
  @timer
  def plot_phase(self, freq, **kwargs):
    """
    subtract : DataFile object
    
    Notes
    -----
    Cyclic colormaps ('hsv', 'twilight_shifted') would be preferred
    because so is the phase (periodic). But apparently jet, with 
    smoother colorscale (fewer colors) does better job...
    
    """
    from fullwavepy.plot.oned import plot_1d, slice_points
    from fullwavepy.signal.phase import wrap_phase 
    
    kwargs['scatt_cmap'] = kw('scatt_cmap', 'jet', kwargs)
    kwargs['scatt_cbar'] = kw('scatt_cbar', True, kwargs)

    try:
      ph_syn, ph_obs, ph_dif = self.phase[freq]
    except (AttributeError, KeyError):
      ph_syn, ph_obs, ph_dif = self._get_phase(freq, **kwargs)
    
    srcs, recs = self._get_sr_coords(**kwargs)

    
    if self.proj.dim.lower() == '2d':
      scoord = 'y'
    elif self.proj.dim.lower() == '3d':
      scoord = 'z'
    else:
      raise ValueError('proj.dim: ' + self.proj.dim)
    
    srcs = slice_points(srcs, scoord)
    recs = slice_points(recs, scoord)
    
    kwargs['scatt_vmin'] = -np.pi
    kwargs['scatt_vmax'] = np.pi
    
    
    plt.subplots(1,3, figsize=[20,8])
    
    #title = str('sid_' + str(self.sid) + '_it_' + str(self.it).rjust(3,'0') + 
                #'_freq_' + str(freq))
    plt.suptitle('Rms ' + str(self.rms_value))
    
    plt.subplot(1,3,1)
    plt.title('syn')
    plot_1d(scatts=[recs], scatt_ampl=[ph_syn.ravel()], **kwargs)
    plt.scatter(*srcs[0], s=20**2, marker='*', c='w', edgecolors='k')
    plt.xlabel('in-line node')
    plt.ylabel('cross-line node')
    #plt.gca().set_aspect('equal')
    #plt.gca().invert_yaxis() # IT IS ALREADY CORRECT
    
    plt.subplot(1,3,2)
    plt.title('obs')
    plot_1d(scatts=[recs], scatt_ampl=[ph_obs.ravel()], **kwargs)
    plt.scatter(*srcs[0], s=20**2, marker='*', c='w', edgecolors='k')
    #plt.gca().set_aspect('equal')
    #plt.gca().invert_yaxis() 
    plt.xlabel('in-line node')
    plt.ylabel('cross-line node')    
    
    plt.subplot(1,3,3)
    plt.title('syn-obs (wrapped)')
    plot_1d(scatts=[recs], scatt_ampl=[ph_dif.ravel()], **kwargs)
    plt.scatter(*srcs[0], s=20**2, marker='*', c='w', edgecolors='k')
    #plt.gca().set_aspect('equal')
    #plt.gca().invert_yaxis()
    plt.xlabel('in-line node')
    plt.ylabel('cross-line node')   

  # -----------------------------------------------------------------------------   
  
  def compare(self, **kwargs):
    from fullwavepy.plot.generic import compare
    Asyn, Aobs, Adif = self.read(**kwargs)
    compare(Asyn, Aobs)   
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------

