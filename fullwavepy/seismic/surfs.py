"""
Various seismic surfaces: seabeds, 
free surfaces, model interfaces etc.

"""
import numpy as np
from autologging import logged #, traced

from arrau.a2d import Surf
from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr3d

class SeaLand(Surf):
  """
  Combined sea bed and land surface.
  """
  def plot(self, **kwargs):
    kwargs['mode'] = kwargs.get('mode', 'shade')
    kwargs['cmap'] = kwargs.get('cmap', ['cmo.phase', 'cmo.deep'])
    kwargs['label'] = kwargs.get('label', 'metres b.s.l.')
    kwargs['alpha'] = kwargs.get('alpha', 1)
    colors = kwargs.get('colors', 'w')
    super().plot(**kwargs)
    super().plot(mode='cr', colors=colors)

# -------------------------------------------------------------------------------
# Legacy classes
# -------------------------------------------------------------------------------
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
  def _2d(self):
    return self.slice(slice_at='z', node=0)
  def plot(self, **kwargs):
    slice_at = kw('slice_at', 'z', kwargs)
    if slice_at == 'z':
      return self.plot_map(**kwargs)
    else:
      return self._2d().plot(**kwargs)
  def plot_map(self, **kwargs):
    kwargs['slice_at'] = 'z'
    kwargs['node'] = 0
    return super().plot(**kwargs)
@logged
class FreeSurface(Surface):
  pass # NOTE it shouldn't be plotted the same way as BathyTopo!
@logged
class Seabed(Surface):
  pass
@logged
class BathyTopo(Surface):
  """
  """
  def extract_seabed(self, dx, **kwargs):
    self.sb = Seabed(np.clip(self, self.z_sea, None) / dx)
    self.sb.extent = self.extent
    return self.sb
  def extract_freesurf(self, add, dx, **kwargs):
    self.__log.debug('self.extent %s' % self.extent)
    self.fs = FreeSurface(np.clip(self, None, self.z_sea) / dx + add)
    self.fs.extent = self.extent
    return self.fs
  def plot(self, **kwargs):
    kwargs['cmap'] = kw('cmap', [], kwargs)
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    return super().plot(**kwargs)
  