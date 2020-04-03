"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr


@traced
@logged
class Point(np.ndarray):
  """
  """
  def __new__(cls, xyz, **kwargs):
    return np.asarray(xyz).view(cls)

  # -----------------------------------------------------------------------------
  
  def __array_finalize__(self, obj):
    if obj is None: return

  # -----------------------------------------------------------------------------

  def find_neighs(self, r, **kwargs):
    """
    Find a cube of nodes surrounding the point.
    
    Notes
    -----
    Works in any dimension (checked 1-3), although in 1D 
    the resultant array should be flattened.
    
    Last axis of np.array(np.meshgrid(*ranges, indexing='ij')).T
    is a tuple of coordinates of a given point. e.g. (x,y,z) in 3D.
    
    We swap the first and the second (!) but last axis to have the Z-coordinate 
    to be the last and thus fastest-changing (np.arrays are row-major) index.
    
    """
    from fullwavepy.math.generic import neighs1d

    ranges = []
    for i in self:
      ranges.append(neighs1d(i, r))
    
    self.neighs = np.array(np.meshgrid(*ranges, indexing='ij')).T.swapaxes(0,-2)
    return self.neighs

  # -----------------------------------------------------------------------------    
  
  def spread(self, r, funcs, **kwargs):
    """
    Spread the point onto a cuboid using
    in general different functions along each coordinate 
    axis.
    
    Parameters
    ----------
    funcs : list
    
    Notes
    -----
    See Hicks 2002, Geophysics for details.
    
    We use the same r for neighbours and the window as 
    outside the window values are zero by definition.
    
    """
    assert len(funcs) == len(self)
    
    # CUBE OF (x,y,z) TUPLES
    cube = self.find_neighs(r, **kwargs)
    # CENTER THE COORDINATE SYSTEM AT self
    dists = cube - self
    # ND-DELTA IS A PRODUCT (!) OF 1D-ONES
    vol = np.ones(dists[...,0].shape)
    # APPLY ALONG TUPLE AXIS, I.E. TAKE POINTS COORDS AS AN ARGUMENT
    coord_axis = -1
    # DEAL WITH ONE COORDINATE AT A TIME
    for i, func in enumerate(funcs):
      # WRAPPER TO ACT ON A SINGLE COORDINATE OF AN ND-POINT
      func_of_xyz = lambda point : func(point[i])
      vol *= np.apply_along_axis(func_of_xyz, coord_axis, dists)
    
    self.vol = vol
    return self.vol

  # -----------------------------------------------------------------------------
  
  def interp(self, **kwargs):
    pass
  
  def interp_hicks(self, **kwargs):
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class Points(Arr):
  """
  """
  def info():
    raise NotImplementedError
  def compare():
    raise NotImplementedError
  def compare_subplots():
    raise NotImplementedError  


# -------------------------------------------------------------------------------


@traced
@logged
class Nodes(Points):
  pass # check if all int
  

# -------------------------------------------------------------------------------


