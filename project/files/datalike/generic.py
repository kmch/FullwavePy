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
  Formats other than .sgy not supported
  at the moment.
  
  Data-processing methods can be applied 
  to the whole 'lump' (file before splitting)
  or to gather files which are instances
  of a child class.
  
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

  def _read_picks(self, fname, kill_gaps=True, **kwargs):
    """
    fname 
      file containing picks in tlPicker format (UO).
    
    """
    from pandas import read_csv
    #big_number = 1e8
    
    if kill_gaps:
      self.__log.warning('Gaps in picks will be filled with (tracelen + 1) value ' +
                         'to effectively kill these traces by a top mute')
    
    # GET ALL PICKS FOR THIS STATION
    sid_times = read_csv(fname, header=None, sep=r"\s*", engine='python', 
                            usecols=[2,4], names=['sid','t'])
    
    picks = dict(zip(sid_times.sid, sid_times.t))
    
    # SELECT ONLY THOSE CHANNELS (SHOTS) PRESENT IN self
    file_chnls = self._gethw('fldr')
    
    file_times = []
    for fch in file_chnls:
      try:
        time = picks[fch]
      except KeyError:
        if kill_gaps:
          time = (self.proj.ns + 1) * self.proj.dt
        else:
          time = 0 # TOP MUTE WILL LEAVE IT INTACT
      file_times.append(time)
    
    return file_times
 
  # ----------------------------------------------------------------------------- 
 
  def plot(self, **kwargs):
    """
     fullwavepy.plot all of the split files.
    
    """
    from fullwavepy.plot.oned import colors
    
    spect = kw('spect', None, kwargs) # FIXME: (REPEATED) -> FUNCTION 
    if spect == 'ampl' or spect == 'power':
      kwargs['cmap'] = 'hot'
      kwargs['center_cmap'] = False
    else:
      kwargs['cmap'] = 'seismic'
      kwargs['center_cmap'] = True
    
    super().plot(**kwargs)

    # PICKS IF PRESENT
    pick_files = kw('pick_files', [], kwargs)    
    clrs = colors(len(pick_files), cmap='summer')
    for pick_file in pick_files:
      self.__log.info('Reading picks from: ' + pick_file)
      picks = self._read_picks(pick_file, **kwargs) 
      self.__log.debug('Converting picks from sec to samples.')
      picks = np.array(picks) / self.proj.dt
      plot(picks, 'o', c=next(clrs), fillstyle='none')
      #self.__log.info('Setting ylim to ')
      #plt.ylim(self.proj.ns, 0)
    
  # -----------------------------------------------------------------------------  
  
  def plot_all_lines(self, **kwargs):
    """
    Savfullwavepy.plots of all lines to png.
    
    """
    kwargs['cbar'] = kw('cbar', True, kwargs)
    
    lines = self.line
    for sid, val in lines.items():
      for lid in val:
        title = strip(lines[sid][lid].fname)
        plt.title(title)
        lines[sid][lid].plot(**kwargs)
        plt.savefig(title)
        plt.close()
        
  # -----------------------------------------------------------------------------
   
  #def compare(self, fname2, **kwargs):
    #"""
    
    #"""
    #pass
    #from fullwavepy.plot.generic import compare
    #
    #kwargs['norm'] = kw('norm', 'max', kwargs)
    #
    #spect = kw('spect', None, kwargs)
    #if spect == 'ampl' or spect == 'power':
    #  kwargs['cmap'] = 'hot'
    #  kwargs['center_cmap'] = False
    #else:
    #  kwargs['cmap'] = 'seismic'
    #
    #compare(self.fname, fname2, **kwargs)
    
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
    A = self.read(scoord=None)
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
      self.suffix += 'Raw'

    if create and filt:
      self.filtered = self.__class__(proj, path, create=False, filt=True, raw=0, mute=0)
    elif filt:
      self.suffix += 'Filt'
      
    if create and mute:
      self.muted = self.__class__(proj, path, create=False, mute=True, raw=0, filt=0)
    elif mute:
      self.suffix += 'Mute'      

    self.name = proj.name + '-' + self.suffix + '.' + self.ext
    self.fname = path + self.name    
    self.__log.debug('self.fname: ' + self.fname)  
  
  # -----------------------------------------------------------------------------  
  
  def filt(self, **kwargs):
    kwargs['overwrite'] = kw('overwrite', True, kwargs)
    try:
      self.filtered.dupl(self.raw.fname)
      self.filtered.filt(**kwargs)
      self.dupl(self.filtered.fname)
    except AttributeError: # RECURSION
      super().filt(**kwargs)
  
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
  
  def mute(self, fbreaks, ntaper=100, twin=1, **kwargs):
    #try:
      #self.muted.dupl(self.filtered.fname)
      #self.muted.mute(fbreaks, ntaper, twin, **kwargs) 
    #except AttributeError: # RECURSION
      #pass 
    
    from fullwavepy.signal.su import su_mute
    from fullwavepy.ioapi.generic import save_txt
    
    fbreaks = np.array(fbreaks)
    picks = fbreaks * self.proj.dt
    bpicks = picks + twin
    nmute = len(picks)
    xmute = range(1, nmute + 1)
    tmute = picks
    tmute2 = bpicks
    
    xmute = [str(i) for i in xmute]
    tmute = [str(i) for i in tmute]
    tmute2 = [str(i) for i in tmute2]
    
    for data, prefix in zip([xmute, tmute, tmute2], ['xmute', 'tmute', 'tmute2']):
      file_txt = self.proj.inp.path + prefix + '.txt'
      file_bin = self.proj.inp.path + prefix + '.bin'
      save_txt(file_txt, data)
      o, e = bash('a2b < {} n1=1 > {}'.format(file_txt, file_bin))
      
    # MAYBE IT NEEDS TO BE SPLIT IN TWO BITS WITH tmp.sg IN BETWEEN (AS WAS IN WORKING VERSION)
    cmd =  'segyread tape={} | '.format(self.filtered.fname)
    cmd += 'sumute key=tracr nmute={nmute} mode=0 ntaper={ntaper} xfile={xmute_bin} tfile={tmute_bin} | sumute key=tracr nmute={nmute} mode=1 ntaper={ntaper} xfile={xmute_bin} tfile={tmute2_bin} | '.format(nmute=nmute, ntaper=ntaper, xmute_bin='xmute.bin', tmute_bin='tmute.bin',
                tmute2_bin='tmute2.bin')
    cmd += 'segyhdrs | segywrite tape={fname_muted}'.format(fname_muted=self.muted.fname)
              
    o, e = bash(cmd)
    

# -------------------------------------------------------------------------------

