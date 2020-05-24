"""
This module defines seismic models:
backgrounds and anomalies.

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
class Model(Arr3d):
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

    """
    obj = super().__new__(cls, *args, **kwargs)
    obj.k_fs = kw('k_fs', None, kwargs) # see the docstring
    return obj


# -------------------------------------------------------------------------------


@traced
@logged
class Background(Model):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class Anomaly(Model):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class ModelVp(Model):
  def plot(self, **kwargs):
    kwargs['cmap'] = kw('cmap', 'twilight', kwargs)
    kwargs['center_cmap'] = kw('center_cmap', False, kwargs)
    return super().plot(**kwargs) 

  def info(self, *args, **kwargs):
    self.__log.info('Property: velocity [m/s]')
    super().info(*args, **kwargs)
    


# -------------------------------------------------------------------------------


@traced
@logged
class StartVp(ModelVp):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class LandModel(Model):
  def extract_freesurf(self, **kwargs):
    # self.__log.debug('self.extent %s' % self.extent)
    # self.fs = np.clip(self, None, self.z_sea)
    # self.fs.extent = self.extent
    # return self.fs  
    pass
  
  def topo_max(self, vel_air, **kwargs):
    for k in range(self.shape[-1]):
      if not np.all(self[...,k] == vel_air):
        self.__log.info('Maximum land elevation has array index (along Z axis): %s' % k)
        self.k_max_land = k
        return self.k_max_land
    raise ValueError('Maximum land elevation could not be found')


# -------------------------------------------------------------------------------


@traced
@logged
class MarineModel(Model):
  def extract_seabed(self, **kwargs):
    # self.sb = np.clip(self, self.z_sea, None)
    # self.sb.extent = self.extent
    # return self.sb
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class AmphibiousModel(LandModel, MarineModel):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class Chckr(Model):
  """
  Checkerboard model.
  """
  pass


# -------------------------------------------------------------------------------


@traced
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


# -------------------------------------------------------------------------------


@traced
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


# -------------------------------------------------------------------------------


@timer
@traced
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


# -------------------------------------------------------------------------------

