"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.array import Arr3d


@traced
@logged
class Data(Arr3d):
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
