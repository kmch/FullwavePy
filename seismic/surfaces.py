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


@traced
@logged
class Surface(Arr3d):
  """

  """
  def _box2inds(self, box, **kwargs):
    """
    Interim fix. Should do sth more elegant. Why for Arr3d extent is not flat?
    
    """
    box = np.array(box)
    extent = np.array(self.extent)
    assert len(box.shape) == 1
    assert len(box) == len(extent.flatten())
    box = box.reshape(extent.shape)
    inds = np.zeros(box.shape)
    for axis, _ in enumerate(box[:-1]):
      self.__log.debug('axis %s' % axis)
      b0, b1 = box[axis]
      inds[axis][0] = self._metre2index(b0, axis)
      inds[axis][1] = self._metre2index(b1, axis) + 1 # NOTE: FOR np.arange(b1, b2) etc.
    
    inds[-1] = np.array([0, self.shape[-1]]) # FULL RANGE FOR DUMMY Z 

    return inds.astype(int)

  def plot(self, **kwargs):
    kwargs['slice_at'] = 'z'
    kwargs['node'] = 0
    super().plot(**kwargs)  

# -------------------------------------------------------------------------------


@traced
@logged
class BathyTopo(Surface):
  def plot(self, **kwargs):
    kwargs['cmap'] = kw('cmap', [], kwargs)
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    super().plot(**kwargs) 


# -------------------------------------------------------------------------------