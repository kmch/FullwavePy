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
from fullwavepy.ndat.points import GenericPoint, Points3d
from fullwavepy.math.funcs import kaiser, sinc, dsinc_dx


@traced
@logged
class PointSR(GenericPoint):
  """
  
  """
  def spread(self, *args, **kwargs):
    vol = super().spread(*args, **kwargs)
    self.__log.debug('vol.extent' + str(vol.extent))
    self.vol = VolumeSR(vol)
    self.vol.extent = vol.extent
    self.vol.coords = vol.coords
    return self.vol
    
    
  # aka multipole
  def check_fs_pos(self, **kwargs):
    # not that trivial probably
    pass  
  def plot(self, **kwargs):
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class VolumeSR(Arr3d): # FIXME: WE NEED COORDINATES OF ARRAY ELEMENTS (JUST CORNERS?)
  """
  """
  def split(self, **kwargs):
    pass
# def check_fs_pos(self, **kwargs):
#   pass
# 
# def spread_factors(self, **kwargs):
#   nsrcs = []
#   while diverged:
#     for src in srcs:
#       nsrcs.append(src.spread_n_bounce())
#     srcs = nsrcs
#     self._check_convergence()
# 
# def _check_convergence():
#   pass
# 
# def inject(self, wf, **kwargs):
#   pass



# -------------------------------------------------------------------------------


@traced
@logged
class Monopole(PointSR):
  """
  """
  def spread(self, r, **kwargs):
    func = lambda x : kaiser(x, r) * sinc(x)
    funcs = [func for i in range(len(self))]
    return super().spread(r, funcs, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class Dipole(PointSR):
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
    func1 = lambda x : kaiser(x, r) * sinc(x) 
    func2 = lambda x : kaiser(x, r) * dsinc_dx(x)
    
    funcs = [func1 for i in range(len(self))]
    funcs[self.axis] = func2
    
    return super().spread(r, funcs, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class DipoleX(Dipole):
  def __new__(cls, xyz, **kwargs):
    return super().__new__(cls, xyz, 0, **kwargs)  


# -------------------------------------------------------------------------------


@traced
@logged
class DipoleY(Dipole):
  def __new__(cls, xyz, **kwargs):
    return super().__new__(cls, xyz, 1, **kwargs) 


# -------------------------------------------------------------------------------


@traced
@logged
class DipoleZ(Dipole):
  def __new__(cls, xyz, **kwargs):
    return super().__new__(cls, xyz, 2, **kwargs) 


# -------------------------------------------------------------------------------


@traced
@logged
class SRs(Points3d):
  def __new__(cls, dictio, **kwargs):
    for key, val in dictio.items():
      dictio[key] = PointSR(val)
    return super().__new__(cls, dictio, **kwargs)

  def set_type(self, srtype_ids, **kwargs):
    """
    it will be read from the file pgy / geo instead...?
    
    PROTEUS convention of naming data components.
    
    """
    mapp = {0 : Monopole,
            1 : DipoleZ, # NOT X!
            2 : DipoleY,
            3 : DipoleX,
           }
    
    assert len(srtype_ids) == len(self)
    for srtype_id, [k, v] in zip(srtype_ids, self.items()):
      self[k] = mapp[srtype_id](v)
  
# -------------------------------------------------------------------------------  


@traced
@logged
class Sources(SRs):
  def plot(self, *args, **kwargs):
    kwargs['marker'] = kw('marker', '*', kwargs)
    kwargs['markersize'] = kw('markersize', 10, kwargs)
    kwargs['markeredgecolor'] = kw('markeredgecolor', 'k', kwargs)
    kwargs['markerfacecolor'] = kw('markerfacecolor', 'w', kwargs)
    super().plot(*args, **kwargs)

  # -----------------------------------------------------------------------------

  def plotly(self, *args, **kwargs):
    kwargs['mode'] = kw('mode', 'markers', kwargs)
    kwargs['color'] = kw('color', 'black', kwargs)
    kwargs['size'] = kw('size', 2, kwargs)
    return super().plotly(*args, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class Receivers(SRs):
  def plot(self, **kwargs):
    kwargs['annotate'] = False
    kwargs['s'] = 1e-2
    kwargs['c'] = 'gray'
    kwargs['alpha'] = 1
    super().plot(**kwargs)
  
  # -----------------------------------------------------------------------------  

  def plotly(self, *args, **kwargs):
    kwargs['mode'] = kw('mode', 'markers', kwargs)
    kwargs['color'] = kw('color', 'grey', kwargs)
    kwargs['size'] = kw('size', 1, kwargs)
    return super().plotly(*args, **kwargs)

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------










#@traced
#@logged
#class Src(Multipole):
#  def spread_n_bounce(self, **kwargs):
#    pass
#  #def spread_factors(self, **kwargs):
#    #self.find_neighs()
#  def spread_bounce(self, **kwargs):
#    pass


# -------------------------------------------------------------------------------


#traced
#logged
#lass SuperSrc(Src):
# """
# """
# def check_fs_pos(self, **kwargs):
#   pass
# 
# def spread_factors(self, **kwargs):
#   nsrcs = []
#   while diverged:
#     for src in srcs:
#       nsrcs.append(src.spread_n_bounce())
#     srcs = nsrcs
#     self._check_convergence()
# 
# def _check_convergence():
#   pass
# 
# def inject(self, wf, **kwargs):
#   pass


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

