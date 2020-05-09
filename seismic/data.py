"""
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


@traced
@logged
class DataFile(ArrayFile):
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

  def plot_metadata(self, **kwargs):

    return ax


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
  
  # def plot_phase(self, freq, ph_type, **kwargs):# sync it better with dumpcomp
  #   """
  #   ph_type : syn/obs/dif
  #   """
  #   if not hasattr(self, 'head'):
  #     self.__log.warn('self.head not found. Returning')
  #     return
    
  #   plt.scatter(self.head.sx, self.head.sy, vmin=-np.pi, vmax=+np.pi,
  #               c=self.head['phase %s (%s Hz)' % (ph_type, freq)])

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
