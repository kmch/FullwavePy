"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr3d
from fullwavepy.ioapi.generic import ArrayFile



@traced
@logged
class DataFileSgy(ArrayFile, SgyFile):
  pass


@traced
@logged
class DataFile(ArrayFile):
  pass


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
