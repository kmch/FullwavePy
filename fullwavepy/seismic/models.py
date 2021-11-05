"""
This module defines seismic models:
backgrounds and anomalies.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
from abc import ABC, abstractmethod
from autologging import logged#, traced
import numpy as np

from arrau.a3d import Arr3d
# -------------------------------------------------------------------------------
# to get rid of
# -------------------------------------------------------------------------------
from fullwavepy.generic.decor import timer #FIXME
from fullwavepy.generic.parse import kw, del_kw #FIXME
from fullwavepy.ndat.arrays import A3d #FIXME
import matplotlib.pyplot as plt #FIXME

class Mod:
  def __init__(self, arr, param='vp'):
    self.arr = arr
    self._init_param(param)
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kwargs.get('cmap', 'twilight')
    self.arr.plot(*args, **kwargs)
  def _init_param(self, param):
    if param == 'vp':
      self.param = Vp()
    else:
      raise ValueError()    
class Param(ABC):
  pass
class Vp(Param):
  pass
# -------------------------------------------------------------------------------
# Old-ish
# -------------------------------------------------------------------------------
@logged
class Model(A3d):
  """
  Generic seismic model storing an arbitrary,
  single parameter of the subsurface.

  """
  def __new__(cls, *args, **kwargs):
    """
    Parameters
    ----------
    k_fs : int (optional)
      Array index along Z-axis (depth) corresponding 
      to the origin of coordinates which is equivalent
      in SEGY to sea level. 
      In Fullwave3D for btop = 0 sea level is at node 0
      i.e. one node above the model grid. 
      See caldera.ipynb for details.
    
    Notes
    -----
    Run is_marine_or_land every time?

    """
    obj = super().__new__(cls, *args, **kwargs)
    obj.k_fs = kw('k_fs', None, kwargs) # see the docstring
    obj.vel_air = kw('vel_air', 0, kwargs)
    obj.marine_or_land = kw('marine_or_land', 'both', kwargs)
    return obj

  # -----------------------------------------------------------------------------

  def carve(self, *args, **kwargs):
    """
    #FIXME
    Two steps:
    1. Do the standard Arr.carve()
    2. Decide whether the extracted model is
       marine, land or both.
    """
    obj = super().carve(*args, **kwargs)
    # obj.marine_or_land = self.is_marine_or_land(**kwargs)
    return obj


  # -----------------------------------------------------------------------------

  def is_marine_or_land(self, **kwargs):
    """
    Decide whether a model is marine, land or both, 
    based on the topmost vertical slice and air-velocity.

    """
    raise NotImplementedError('Address chopped-topo and other special cases.')
    vel_air = self.vel_air
    top_slice = self[..., 0]

    if np.all(top_slice == vel_air):
      self.marine_or_land = 'land'
    elif np.any(top_slice == vel_air):
      self.marine_or_land = 'both'
    else:
      self.marine_or_land = 'marine'
    
    return self.marine_or_land
  
  # -----------------------------------------------------------------------------
    
  def topo_max(self, vel_air, **kwargs):
    """
    Find maximum elevation of land topography.
    
    """
    if self.marine_or_land == 'marine':
      self.__log.debug('Nothing to do here...')
    else:
      self.find_topo_max(vel_air, **kwargs)

  # -----------------------------------------------------------------------------

  def find_topo_max(self, vel_air, **kwargs):
    """
    FIXME: rename -> find... once the whole package
    is shippable and changes in the notebooks are easy
    and reliable.
    """
    for k in range(self.shape[-1]):
      if not np.all(self[...,k] == vel_air):
        self.__log.info('Maximum land elevation has array index (along Z axis): %s' % k)
        self.k_max_land = k
        return self.k_max_land
    raise ValueError('Maximum land elevation could not be found')
@logged
class Background(Model):
  pass
@logged
class Anomaly(Model):
  pass
@logged
class Chckrbrd(Anomaly):
  """
  FIXME: Rename -> ModelCheckerboard
  """
  pass
@logged
class ModelVp(Model):
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'twilight', kwargs)
    kwargs['center_cmap'] = kw('center_cmap', False, kwargs)
    return super().plot(*args, **kwargs) 

  def info(self, *args, **kwargs):
    self.__log.info('Property: velocity [m/s]')
    super().info(*args, **kwargs)
@logged
class StartVp(ModelVp):
  """
  FIXME: Rename to ModelVpStart
  """
  pass
# -------------------------------------------------------------------------------
# SOME USEFUL FUNCTIONS
# -------------------------------------------------------------------------------
def merge_model_and_fs(model, fs, vel_air, **kwargs):
  """
  Merge the model and an independently-derived
  free surface. 
  """
  pass
# -------------------------------------------------------------------------------
# GENERIC BACKGROUNDS AND ANOMALIES (MOVE AT THE EXPENSE OF COHESION?)
# -------------------------------------------------------------------------------
@logged
def grad(shape, vtop, vbot, **kwargs):
  """
  Create a vertical gradient model.
  
  """
  A = np.zeros(shape)
  dv = (vbot - vtop) / (shape[-1] - 1)
  grad._log.debug('velocity increment: {}'.format(dv))
  for iz in range(shape[-1]):
    A[..., iz] = vtop + iz * dv
  return A
@logged
def sphere(dims, center, radius, ampl=1, **kwargs):
  """
  Create a spherical anomaly.
  
  dist = np.linalg.norm(np.array((x,y,z)) - center)
  
  Notes
  -----
  See NumPy array boolean masking.
  
  """
  coords = np.indices(dims)
  distance = np.sqrt((coords[0] - center[0])**2 + (coords[1]-center[1])**2 + (coords[2]-center[2])**2) 
  A = (distance <= radius).astype(int) * ampl
  return A
@timer
@logged
def gauss(dims, centers, radius, **kwargs):
  """
  gaussian is a bottleneck, needs to be
  vectorized.
  
  """
  from fullwavepy.numeric.funcs import gauss as gaussian
  coords = np.indices(dims)
  A = np.zeros(dims)
  for center in centers: 
    distance = np.sqrt((coords[0] - center[0])**2 + (coords[1]-center[1])**2 + (coords[2]-center[2])**2) 
    a = gaussian(distance, 0, radius)
    A += a
  return A

