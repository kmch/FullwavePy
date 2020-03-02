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
from fullwavepy.generic.array import WigglyData
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
    super().__init__(proj, path, **kwargs)

  # ----------------------------------------------------------------------------- 

  def split(self, **kwargs):
    self.__log.warn('Not implemented. Skipping.')

  # -----------------------------------------------------------------------------
  
  def files(self, **kwargs):
    self.__log.warn('Not implemented. Skipping.')

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# SPECIFIC FILES (UNIQUE IDs)
# -------------------------------------------------------------------------------


@traced
@logged
class SynDataFileTtr(DataFileTtr, SynDataFile):
  def __init__(self, proj, path, **kwargs):
    suffix = 'Observed-Time'
    super().__init__(suffix, proj, path, **kwargs)
  

# -------------------------------------------------------------------------------


@traced
@logged
class ObsDataFileTtr(DataFileTtr, ObsDataFile):
  def __init__(self, proj, path, **kwargs):
    suffix = 'Observed-Time'
    super().__init__(suffix, proj, path, **kwargs)


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
  
  def read(self, **kwargs):
    try:
      self.array = super().read(**kwargs)
    except FileNotFoundError as err: #FIXME? MOVE TO A MORE GENERAL CLASS
      self.__log.warn('%s not found. Now looking for a txt version...' % err)
      fname = self.fname
      self.fname = strip(fname) + '.txt'
      self.array = super().read(**kwargs)
      self.fname = fname
    return self.array
      

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
    self.phase = {}
    #self.syn = SynDataFileTtr(suffix+'_syn', proj, path)
    #self.obs = ObsDataFileTtr(suffix+'_obs', proj, path)
    #self.dif = SynDataFileTtr(suffix+'_dif', proj, path)
    
  # -----------------------------------------------------------------------------      
  
  def read(self, **kwargs):
    """
    isep: 
    # SEPARATOR TO SPLIT ARRAY INTO 3 CHUNKS OF THE SAME LEN, CHECKED
    """
    from fullwavepy.generic.array import WigglyData    
    self.__log.info('Reading ' + self.fname + '...')
    self.array = super().read(**kwargs)
    self.read_header(**kwargs)
    
    isep  = int(len(self.array) / 3) 
    di = {'syn': WigglyData(self.array[      :isep]),
          'obs': WigglyData(self.array[+isep:-isep]),
          'dif': WigglyData(self.array[-isep:     ]),
         }
    
    for key in di.keys():
      setattr(self, key, di[key])
    
    return di

  # -----------------------------------------------------------------------------
  
  def read_header(self, overwrite=False, **kwargs):
    """
    Add reading csv, useful for heavy sgy files.
    
    Offset is measured by Fullwave3D in 3D, not 2D like SegyPrep!
    This was the source of bugs but is now fixed.
    
    """
    if (not hasattr(self, 'head')) or overwrite:
      df = self.proj.i.obs.read_header(overwrite, **kwargs)
      
      # SELECT THE CORRECT SHOT
      recipr = bool(self.proj.i.sp.read()['reciprocity'])
      if recipr:
        df = df[df.tracf == self.sid]
      else:
        df = df[df.fldr == self.sid]
      
      # CALCULATE THE OFFSET
      self.__log.warn('Taking -gelev as this is positive for OBS PROTEUS, double-check land stations!')
      df['offset3d'] = np.sqrt((df['sx'] - df['gx'])**2 + 
                               (df['sy'] - df['gy'])**2 + 
                               (df['selev'] + df['gelev'])**2)
      
      
      # SELECT THE OFFSET RANGE AS IN ITERATION INFO
      minoff = kw('minoff', 0, self.proj.i.rnf.iters[self.it])
      maxoff = kw('maxoff', 1e9, self.proj.i.rnf.iters[self.it])
      self.__log.debug('minoff, maxoff for this iteration: {}, {}'.format(minoff, maxoff))
      df = df[(df.offset3d >= minoff) & (df.offset3d <= maxoff)]
      # RESET INDEX, OTHERWISE BUGGY split() ETC.
      df.reset_index(drop=True, inplace=True)       
      
      self.head = df
      
    return self.head  
  
  # -----------------------------------------------------------------------------

  def split(self, overwrite=False, **kwargs):
    kwargs['overwrite'] = overwrite
    self.read(**kwargs)
    self.read_header(**kwargs)
    
    if not (hasattr(self.syn, 'lid')) or overwrite:
      # LINE IDs
      lid_hw = self.proj.sgy.hw['lid'] # IT SHOULD BE ep FOR PROTEUS
      lids = self.head[lid_hw].unique()
      
      for obj in [self.syn, self.obs, self.dif]:
        obj.lid = {}
        for lid in lids:
          subhead = self.head[:][self.head[lid_hw] == lid]
          obj.lid[lid] = WigglyData(np.take(obj, indices=subhead.index, axis=0))
  
  # -----------------------------------------------------------------------------
  
  def _get_first_breaks(self, *args, **kwargs):
    """
    Note: it is hard to merge with SynDataFile one.
    
    """
    from fullwavepy.signal.phase import first_breaks    
    Asyn, Aobs, Adif = self.read(**kwargs)
    self.fb = np.ravel(first_breaks(Asyn, *args, **kwargs))
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
    freq_max = self.proj.i.rnf.iters[self.it]['freq']
    if freq > freq_max:
      raise ValueError('You should not extract the phase at the freq. ' +
                       'above the iteration-block high-cut freq.: %s' % freq_max)
    
    if (not freq in self.phase) or overwrite:
      from fullwavepy.signal.phase import first_breaks, extract_phase, wrap_phase
      from fullwavepy.generic.math import rms    
      self.__log.info('Getting phase info from ' + self.fname)
      
      self.read(**kwargs)
      self.read_header(**kwargs)
      
      if not hasattr(self, 'fb'):
        self.fb = first_breaks(self.syn, **kwargs)
      
      ph_syn = extract_phase(self.syn, self.fb, self.proj.dt, freq, **kwargs)
      ph_obs = extract_phase(self.obs, self.fb, self.proj.dt, freq, **kwargs)
      ph_dif = ph_syn - ph_obs
      
      ph_dif = np.array([[wrap_phase(i) for i in j] for j in ph_dif])
      ph_syn, ph_obs, ph_dif = [np.ravel(i) for i in [ph_syn, ph_obs, ph_dif]]
      
      self.rms_value = rms(ph_dif)
      self.__log.info('RMS of wrapped phase-differences: ' + 
                      str(self.rms_value))
      
      self.head['phase syn (%s Hz)' % freq] = ph_syn
      self.head['phase obs (%s Hz)' % freq] = ph_obs
      self.head['phase dif (%s Hz)' % freq] = ph_dif
      
      self.phase[freq] = {'syn': ph_syn,
                          'obs': ph_obs,
                          'dif': ph_dif}
    return self.phase[freq]
      
  # -----------------------------------------------------------------------------   
  
  def exclude_bad_data(self, max_ph_dif, **kwargs):
    pass
  
  # -----------------------------------------------------------------------------  
  
  def plot(self, **kwargs):
    self.__log.info('Plotting whole content (syn, obs, dif) of ' + self.fname)
    self.read(**kwargs)
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs)
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    self.array.plot(**kwargs)
   
  # -----------------------------------------------------------------------------
  
  #def plot_line(self, **kwargs):
    #self.split(**kwargs)
    
   
  # -----------------------------------------------------------------------------
  
  #def compare(self, **kwargs):
    #from fullwavepy.plot.generic import compare
    #Asyn, Aobs, Adif = self.read(**kwargs)
    #compare(Asyn, Aobs)   
    
  # -----------------------------------------------------------------------------  
  
  def px_head(self, data_col, cmap='jet', **kwargs):
    import plotly.express as px
    fig = px.scatter(self.head, x='sx', y='sy', color=data_col, 
                     color_continuous_scale=cmap)
    return fig    

  # -----------------------------------------------------------------------------   
  
  def px_phase(self, freq, phase='dif', **kwargs):
    """
    Quick interactive plot. Can't do subplots with plotly :(
    
    Syntactic sugar: just a special case of px_head
    
    """
    import plotly.express as px
    
    self._get_phase(freq, **kwargs)
    
    data_col = 'phase %s (%s Hz)' % (phase, freq)
    fig = self.px_head(data_col, **kwargs)
    return fig
    
  # -----------------------------------------------------------------------------   
  
  def plot_phase(self, freq, annotate=False, **kwargs):
    """
    subtract : DataFile object
    
    Notes
    -----
    Cyclic colormaps ('hsv', 'twilight_shifted') would be preferred
    because so is the phase (periodic). But apparently jet, with 
    smoother colorscale (fewer colors) does better job...
    
    """
    cmap = kw('cmap', 'jet', kwargs)
    cbar = kw('cbar', True, kwargs)
    
    self._get_phase(freq, **kwargs)
    
    fig, ax = plt.subplots(1,3, figsize=[20,8])
    fig.suptitle('Rms ' + str(self.rms_value))
    
    for i, ph_type in enumerate(['syn', 'obs', 'dif']):
      ax = plt.subplot(1,3,i+1)
      ax.set_title(ph_type)
      ax.scatter(self.head.sx, self.head.sy, vmin=-np.pi, vmax=+np.pi, cmap=cmap,
                 c=self.head['phase %s (%s Hz)' % (ph_type, freq)])
      ax.scatter(self.head.gx[0], self.head.gy[0], s=20**2, 
                  marker='*', c='w', edgecolors='k')
      ax.set_xlabel('x, metres')
      ax.set_ylabel('y, metres')
      ax.set_aspect('equal')
      
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


