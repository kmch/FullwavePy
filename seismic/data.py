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
  
  def plotf(self, fig, tylim=None, fylim=None, *args, **kwargs):
    """
    both time and freq
    """
    if fig is None:
      figure(16,8)
    plt.subplot(121)
    self.plot(*args, **kwargs)
    if tylim is not None:
      plt.ylim(tylim)
    
    plt.subplot(122)
    self.__log.warn('REMEMBER ABOUT dt')
    kwargs['dt'] = kwarg_over_attr('dt', kwargs, self)
    kwargs = dict(kwargs, spect='ampl', cmap='hot', center_cmap=False)
    self.plot(*args, **kwargs)
    if tylim is not None:
      plt.ylim(fylim)
    # plt.gca().set_aspect('auto')

@traced
@logged
class DataFileSgy(DataFile, SgyFile):
  def read(self, **kwargs):
    arr = SgyFile.read(self, **kwargs)
    self.array = self.cast(arr, **kwargs)
    return self.array



@traced
@logged
class Data(Arr3d):
  """
  Seismic data.

  """
  def interleave(self, othe, **kwargs):
    return super().interleave(othe, slice_at='y', node=0, **kwargs)

  def compare(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) #'twilight_shifted'
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    super().compare(*args, **kwargs)
  
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) #'twilight_shifted'
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
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
