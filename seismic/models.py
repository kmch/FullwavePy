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
  pass


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
  pass



# -------------------------------------------------------------------------------


@traced
@logged
class StartVp(Model):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class LandModel(Model):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class MarineModel(Model):
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
  from fullwavepy.math.funcs import gauss as gaussian
  coords = np.indices(dims)
  A = np.zeros(dims)
  for center in centers: 
    distance = np.sqrt((coords[0] - center[0])**2 + (coords[1]-center[1])**2 + (coords[2]-center[2])**2) 
    a = gaussian(distance, 0, radius)
    A += a
  return A


# -------------------------------------------------------------------------------

