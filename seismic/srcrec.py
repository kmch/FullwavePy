"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.points import Point


@traced
@logged
class SrcRec(Point):
  def check_fs_pos(self, **kwargs):
    pass  
  def plot(self, **kwargs):
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class Monopole(SrcRec):
  """
  """
  def spread(self, r, **kwargs):
    from fullwavepy.generic.math import kaiser, sinc
    func = lambda x : kaiser(x, r) * sinc(x)
    funcs = [func for i in range(len(self))]
    return super().spread(r, funcs, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class Dipole(SrcRec):
  """
  axis : 0, 1 or 2
    corresponds to dipole along X, Y or Z axis respectively
  
  """
  def __new__(cls, xyz, axis, **kwargs):
    assert axis in [0, 1, 2]
    cls.axis = axis
    return super().__new__(cls, xyz, **kwargs)
  
  # -----------------------------------------------------------------------------
  
  def spread(self, r, **kwargs):
    from fullwavepy.generic.math import kaiser, sinc, dsinc_dx
    
    func1 = lambda x : kaiser(x, r) * sinc(x) 
    func2 = lambda x : kaiser(x, r) * dsinc_dx(x)
    
    funcs = [func1 for i in range(len(self))]
    funcs[self.axis] = func2
    
    return super().spread(r, funcs, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------







# -------------------------------------------------------------------------------


@traced
@logged
class Src(SrcRec):
  def spread_n_bounce(self, **kwargs):
    pass
  #def spread_factors(self, **kwargs):
    #self.find_neighs()
  def spread_bounce(self, **kwargs):
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class SuperSrc(Src):
  """
  """
  def check_fs_pos(self, **kwargs):
    pass
  
  def spread_factors(self, **kwargs):
    nsrcs = []
    while diverged:
      for src in srcs:
        nsrcs.append(src.spread_n_bounce())
      srcs = nsrcs
      self._check_convergence()
  
  def _check_convergence():
    pass
  
  def inject(self, wf, **kwargs):
    pass


# -------------------------------------------------------------------------------










@traced
@logged
def xyz2w(xyz, dims, **kwargs):
  """
  """
  x, y, z = xyz
  nx, ny, nz = dims
  return (x - 1) * ny * nz + (y - 1) * nz + z





  """
  """
  def check_fs_pos(self, **kwargs):
    pass
  
  def spread_factors(self, **kwargs):
    nsrcs = []
    while diverged:
      for src in srcs:
        nsrcs.append(src.spread_n_bounce())
      srcs = nsrcs
      self._check_convergence()
  
  def _check_convergence():
    pass
  
  def inject(self, wf, **kwargs):
    pass
  

# -------------------------------------------------------------------------------

