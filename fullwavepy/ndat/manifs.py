"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

About
-----
Definitions of multidimensional manifolds.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr3d, Arr2d


@traced
@logged
class Vol(Arr3d):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class Surf(Arr3d):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class SurfZ(Arr2d):
  """
  Surface of the form z = f(x,y).
  
  """
  def __new__(cls, source, **kwargs):
    return super().__new__(cls, source, **kwargs)
  #def read(self, **kwargs):
    #arr3d = Arr3d(super().read(**kwargs))
    #self = Arr2d(arr2d.slice(slice_at='z', node=0, **kwargs))
    #return self
  #def plot2d(self, **kwargs):
    #self.plot_slice(slice_at='z', **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class Plane(Surf):
  """
  Plane dipping along X.
  
  Parameters
  ----------
  z0 : float
    GENERALIZE TO ANY ANGLE
  tilt_angle : float
    Dip.
  tilt_azim : float 
    Angle 0-360 deg: 0 <=> dip along X.
  
  
  """
  def __new__(cls, dims, z0, tilt_angle, tilt_azim=0, **kwargs):
    """
    """
    from fullwavepy.numeric.funcs import linear
    assert len(dims) == 3
    if dims[-1] != 1:
      cls.__log.warning('Replacing nz=%s with nz=1' % dims[-1])
      dims = np.copy(dims) # OTHERWISE OVERWRITES proj.dims IF PASSED AS DIMS
      dims[-1] = 1
    
    source = np.zeros(dims)
    arr = super().__new__(cls, source, **kwargs)
    
    for y in range(dims[1]):
      arr[:, y, 0] = linear(np.arange(dims[0]), tilt_angle, z0)
    return arr
    

# -------------------------------------------------------------------------------


@traced
@logged
class SurfXYZ(Surf):
  """
  Surface in a parametric form
  X, Y, Z.

  E.g. a torus:
  angle = np.linspace(0, 2 * np.pi, 32)
  theta, phi = np.meshgrid(angle, angle)
  r, R = .25, 1.
  X = (R + r * np.cos(phi)) * np.cos(theta)
  Y = (R + r * np.cos(phi)) * np.sin(theta)
  Z = r * np.sin(phi)
  
  """
  pass


# -------------------------------------------------------------------------------



