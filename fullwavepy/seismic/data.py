"""
Dealing with seismic data can be tricky for two reasons:
it is usually stored as heavy binaries and it has 
associated metadata.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, kwarg_over_attr
from fullwavepy.ndat.arrays import Arr3d
from fullwavepy.ioapi.generic import ArrayFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.plot.generic import *

# FIXME: MERGE WITH datalike!



def qc_datafile(datafile, ep, **kwargs):
  from fullwavepy.plot.generic import figure
  txlim = kw('txlim', None, kwargs)
  tylim = kw('tylim', None, kwargs)
  fxlim = kw('fxlim', None, kwargs)
  fylim = kw('fylim', None, kwargs)
  kwargs['win'] = dict(ep=[ep])
  datafile.read(**kwargs)
  figure(16,8)
  plt.suptitle(datafile.name + ', line ' + str(ep))    
  # time
  plt.subplot(121)
  datafile.array.plot(cmap='seismic', center_cmap=1, **kwargs)
  plt.xlim(txlim)
  plt.ylim(tylim)
  plt.xlabel('trace no.')
  plt.ylabel('sample')
  # frequency
  plt.subplot(122)
  datafile.array.plot(cmap='hot', center_cmap=0, spect='ampl', dt=datafile.dt, **kwargs)
  plt.xlim(fxlim)
  plt.ylim(fylim)
  plt.xlabel('trace no.')
  plt.ylabel('frequency [Hz]')
  plt.gca().set_aspect('auto') 


@traced
@logged
class DataFile(ArrayFile):
  """
  Generic seismic-data file. It is meant to 
  handle arbitrarily large files through `window` 
  and `split` methods that are I/O specific
  and need to be implement in child classes.

  It should generalize excellent implementation
  of split etc. from DumpCompareFile.
  
  no read is defined here!
  """
  def _split(self, *args, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
  
  def _window(self, *args, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
  
  def _read_shot_gather(self, sid, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
  
  def _read_shot_line(self, sid, lid, **kwargs):
    raise NotImplementedError('Overwrite in a child class')


  def cast(self, arr, **kwargs):
    self.array = Data(arr)
    return self.array
  
  def plotf(self, fig, tylim=None, fylim=None, *args, **kwargs): # FIXME: merge with plot.misc time_freq
    """
    both time and freq

    """
    if fig is None:
      figure(16,8)
    
    plt.subplot(121)
    self.plot(*args, **kwargs)
    plt.ylim(tylim)
    
    plt.subplot(122)
    kwargs['dt'] = kwarg_over_attr('dt', kwargs, self) # needed in DFT
    kwargs = dict(kwargs, spect='ampl', cmap='hot', center_cmap=False) # DFT, different cmap
    self.plot(*args, **kwargs)
    plt.ylim(fylim)
  

# -------------------------------------------------------------------------------


@traced
@logged
class DataFileSgy(DataFile, SgyFile):
  def read(self, **kwargs):
    arr = SgyFile.read(self, **kwargs)
    self.array = self.cast(arr, **kwargs)
    return self.array


# -------------------------------------------------------------------------------


@traced
@logged
class Data(Arr3d):
  """
  Seismic data.

  """
  def _xlabels_from_header(self, xlabels_hw='fldr', **kwargs):
    if hasattr(self, 'head'):
      self.__log.info('There is a header associated with this data.' +\
        '%s keyword will be used as xlabels.' % xlabels_hw)
      xlabels = self.head[xlabels_hw] # one per xtick
    # xlabel = xlabels_hw # one per xaxis
      return xlabels
    else:
      return None
  
  def interleave(self, othe, **kwargs):
    return super().interleave(othe, slice_at='y', node=0, **kwargs)
    
  def compare(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) #'twilight_shifted'
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    kwargs['xlabels'] = self._xlabels_from_header(**kwargs)    
    super().compare(*args, **kwargs)
  
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) #'twilight_shifted'
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    kwargs['xlabels'] = self._xlabels_from_header(**kwargs)
    
    super().plot(*args, **kwargs)
  
  def plot_phase_DEBUG(self, freq, ph_type, **kwargs):# sync it better with dumpcomp
    """
    ph_type : syn/obs/dif
    """
    if not hasattr(self, 'head'):
      self.__log.warning('self.head not found. Returning')
      return
  
    plt.scatter(self.head.sx, self.head.sy, vmin=-np.pi, vmax=+np.pi,
                c=self.head['phase %s (%s Hz)' % (ph_type, freq)])


# -------------------------------------------------------------------------------


@traced
@logged
class ObsData(Data):
  """
  Observed (recorded in the field) seismic data.
  
  """
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class SynData(Data):
  """
  Synthetic (calculated) seismic data.
  
  """
  pass


# -------------------------------------------------------------------------------
