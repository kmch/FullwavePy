"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.generic.decor import timer


# -------------------------------------------------------------------------------


@traced
@logged
def grad(shape, vtop, vbot, **kwargs):
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
  from fullwavepy.generic.math import gauss as gaussian
  coords = np.indices(dims)
  A = np.zeros(dims)
  for center in centers: 
    distance = np.sqrt((coords[0] - center[0])**2 + (coords[1]-center[1])**2 + (coords[2]-center[2])**2) 
    a = gaussian(distance, 0, radius)
    A += a
  return A


# -------------------------------------------------------------------------------

