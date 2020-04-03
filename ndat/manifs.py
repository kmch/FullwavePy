"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

About
-----
Definitions of multidimensional manifolds.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw#, exten, strip
from fullwavepy.generic.array import Arr3d


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
class SurfZ(Surf):
  """
  Surface of the form z = z(x,y).
  
  """
  def plot(self, **kwargs):
    kwargs['cmap'] = kw('cmap', [], kwargs)
    self.plot_slice(slice_at='z', **kwargs)
    #self.array = self.read(**kwargs)
    #shape = self.array.shape
    #self.array = 


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

