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
class SRs(Points3d):
  """
  """
  def __new__(cls, dictio, **kwargs):
    for key, val in dictio.items():
      dictio[key] = PointSR(val)
    return super().__new__(cls, dictio, **kwargs)

  # -----------------------------------------------------------------------------

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
 
  # ----------------------------------------------------------------------------- 
 
  def spread_factors(self, srtype_ids, **kwargs):
    self.set_type(srtype_ids, **kwargs)
    #self.sprd_fctrs = {}
    self.hyper = {}

    for srid, sr in self.items():
      #self.__log.info('Calculating spread_factors for SRID ' + str(srid)) 
      #self.sprd_fctrs[srid] = HyperPointSR(sr).spread_factors(**kwargs)
      self.hyper[srid] = HyperPointSR(sr)
    
  # -----------------------------------------------------------------------------


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


@traced
@logged
class HyperPointSR(object):
  """
  """
  def __init__(self, pointsr, **kwargs):
    self.rmax = 4
    self.pointsr = pointsr
    coords = self.pointsr.find_neighs(self.rmax, **kwargs)
    sh = np.array(coords.shape)
    sh[-1] += 1
    # THIS IS THE CUBOIDAL VOLUME OF SPREAD FACTORS 
    # IT WILL BE USED TO INJECT THE SOURCE INTO (OR INTERPOLATE RECEIVER)
    # THE WAVEVIELD
    self.vol = np.zeros(sh)
    self.vol[..., :3] = coords
    self.vol[..., 3] = 0.0 # this will store amplitude, i.e. sprd_fctrs
    
  # -----------------------------------------------------------------------------
  
  def spread_factors(self, **kwargs):
    #self.tospread = [
    
    
    return self.vol

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class PointSR(GenericPoint):
  """
  
  """
  def spread_factors(self, **kwargs):
    r = 3 # SAME AS IN FULLWAVE3D
    self.vol = self.spread(r)
    #self.bounce_off

  
  # -----------------------------------------------------------------------------
  
  def spread(self, *args, **kwargs):
    """
    """
    vol = super().spread(*args, **kwargs)
    self.__log.debug('vol.extent' + str(vol.extent))

    self.vol = VolumeSR(vol)
    self.vol.extent = vol.extent
    self.vol.coords = vol.coords
    return self.vol

  # -----------------------------------------------------------------------------    


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
class VolumeSR(Arr3d):
  """
  """
  def split(self, *args, **kwargs):
    """
    """
    self.flgs = self.flag(*args, **kwargs)
    
    attrs = ['ins', 'at', 'out']
    flags = [self.in_flag, self.acc_flag, self.ext_flag]

    for attr, flag in zip(attrs, flags):
      indices = self.flgs == flag
      coords = self.coords[indices]
      values = self[indices]
      arr = np.zeros(list(values.shape) + [4]) # must be neater syntax
      arr[..., 0:-1] = coords
      arr[..., -1] = values
      self.__log.info('Found %s %s-nodes' % (len(arr), attr))
      setattr(self, attr, arr)

  # -----------------------------------------------------------------------------    

  def flag(self, ine, *args, **kwargs):
    """
    ine: array of interior/accurately at FS/exterior nodes
    """
    self.indices = self._grid_coords_2_egrid_indices(*args, **kwargs)
    
    self.flgs = Arr3d(ine[self.indices])
    self.flgs.extent = self.extent
    
    self.in_flag = ine.in_flag
    self.acc_flag = ine.acc_flag
    self.ext_flag = ine.ext_flag
    
    return self.flgs
  
  # -----------------------------------------------------------------------------

  def _grid_coords_2_egrid_indices(self, elef, efro, etop, **kwargs):
    """
    Rationale
    ---------
    use vol.coords to pick subarray from inext-nodes array etc
    
    convert coords (x,y,z) of grid:
     x = 1, 2, ... nx1 
     
    (coords is an array as vol.coords laid out with 
    x  assumed to be the slowest and z the fastest index)
     
    to indices (!) of the array storing whole extended (!) grid 
     i = 0, ..., enx1
     ...
    order of indices as in coords  
     
    """
    grid = np.copy(self.coords.T.swapaxes(1,-1))
    grid[0] += elef
    grid[1] += efro
    grid[2] += etop
    indices = np.copy(grid-1) # subtract  1 to convert array indices!!!!!!
    indices = tuple(indices) # CRUCIAL, OTHERWISE SLICING WOULDN'T WORK
    return indices

  # -----------------------------------------------------------------------------

  def bounce_off(self, ghs, iss, elef, efro, etop, **kwargs):
    """
    ghs : list of ghosts
    iss : list of corresponding intersects (same len)
    """
    assert len(ghs) == len(iss)
    
    vout = np.copy(self.out)
    vout[:, 0] += elef
    vout[:, 1] += efro
    vout[:, 2] += etop
    # here we compare grid-coords, NOT grid-coords to array-indices
    # => we don't subtract 1

    aghs = np.array(ghs)

    reflected = []
    for vo in vout:
      x, y, z, A = vo
      G = aghs[((aghs[:,0] == x) & (aghs[:,1] == y) & (aghs[:,2] == z))][0]
      I = np.array(iss[ghs.index(list(G))])
      G = G[:3]
      R = np.zeros(4) # WE HAVE TO CREATE IT EVERY TIME
      R[:3] = 2 * I - G
      R[3] = -A
      #R = Monopole(
      reflected.append(R)
      #print('O', vo)
      #print('G', G)
      #print('I', I)
      #print('R', R)
      #print()   
    return np.array(reflected)

  # -----------------------------------------------------------------------------
  
  def plot(self, **kwargs):
    kwargs['slice_at'] = kw('slice_at', 'y', kwargs)
    kwargs['node'] = kw('node', self.shape[1]//2, kwargs)
    super().plot(**kwargs)

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

