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


clip_at = 1e-9 # spread-factors smaller than this will be neglected


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
      self[k] = mapp[srtype_id](v, **kwargs)
 
  # ----------------------------------------------------------------------------- 
 
  def spread_factors(self, srtype_ids, *args, **kwargs):
    """
    """
    self.set_type(srtype_ids, **kwargs)
    #self.sprd_fctrs = {}
    self.hyper = {}

    for srid, sr in self.items():
      self.__log.info('Calculating spread-factors for SRID ' + str(srid)) 
      self.hyper[srid] = HyperPointSR(sr, **kwargs)
      self.hyper[srid].spread_factors(*args, **kwargs)
    
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
    """
    """
    self.pointsr = pointsr
    self.r_hicks = kw('r_hicks', 3, kwargs)

  # -----------------------------------------------------------------------------
  
  def find_vol(self, **kwargs):
    """
    """
    rmax = kw('rmax', 10, kwargs)
    
    coords = self.pointsr.find_neighs(rmax, **kwargs)
    self.origin = coords[0,0,0] # COORDINATES OF THE FIRST NODE 
    
    
    # THIS IS THE CUBOIDAL VOLUME OF SPREAD FACTORS 
    # IT WILL BE USED TO INJECT THE SOURCE INTO (OR INTERPOLATE RECEIVER)
    # THE WAVEVIELD
    sh = np.array(coords.shape)
    sh[-1] += 1
    self.vol = np.zeros(sh)
    self.vol[..., :3] = coords
    self.vol[..., 3] = 0.0 # this will store amplitude, i.e. sprd_fctrs
    
  # -----------------------------------------------------------------------------
  
  def spread_factors(self, ghs, iss, ine, elef, efro, etop, **kwargs):
    """
    Iterative 
    """
    self.find_vol(**kwargs)
    self.vol[..., -1] = 0.0 # OTHERWISE GROWING IN IPYNB
    
    psrs_to_spread = [self.pointsr]
    i_max = 2
    i = 0
    self.snapshots = []
    while i < i_max:
      i += 1    
      n_to_spread = len(psrs_to_spread)
      self.__log.info('Iteration %s, no. of points to spread: %s' % (i, n_to_spread)) 

      if n_to_spread == 0:
        self.__log.info('No more points to spread. Returning.')
        return
      
      self.__log.debug('psrs_to_spread' + str(psrs_to_spread))
      new_psrs_to_spread = []
      for psr_to_spread in psrs_to_spread:
        self.__log.debug('psr_to_spread ' + str(psr_to_spread))
        psr_to_spread.spread_factors(r_hicks=self.r_hicks)
        psr_to_spread.vol.split(ine, elef, efro, etop)
        
        self._update_vol(np.copy(psr_to_spread.vol.ins))
        
        
        reflected = psr_to_spread.vol.bounce_off(ghs, iss, elef, efro, etop)
        if len(reflected) > 0:
          new_psrs_to_spread += reflected # append would create a nested list (undesired)
      
      self.snapshots.append(np.copy(self.vol))
      psrs_to_spread = list(new_psrs_to_spread)

  # -----------------------------------------------------------------------------
  
  def _update_vol(self, nodes_in, **kwargs):
    for node_in in nodes_in:
      #self.__log.info('node_in ' + str(node_in))
      xyz = node_in[:3]
      amp = node_in[3]
      if abs(amp) < clip_at: # skip tiny factors for speed and compactness of volume
        self.__log.debug('Skipping inside node %s with small ampl %s' % (str(xyz), str(amp)))
        continue      
      
      ijk = xyz - self.origin
      ijk = tuple((int(n) for n in ijk))
      try:
        prv = self.vol[ijk][-1]
      except IndexError:
        self.__log.debug('Skipping inside (medium) node %s which is outside the volume' % (str(xyz)))
        continue
      
      
      now = prv + amp
      self.vol[ijk][-1] = now
      #if abs(amp) > 1e-4:
      if True:
        self.__log.debug('Updated %s by %s from %s to %s ' % (str(ijk), 
                                                             str(amp), 
                                                             str(prv),
                                                             str(now)))

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class PointSR(GenericPoint):
  """
  ATTENTION 
  It is supposed to have coordinates in nodes 
  of the regular (not extended!) grid.
  NOTE
  Left-front-top corner of the grid is (1,1,1)
  Right-back-bottom corner of the grid is (nx1,nx2,nx3)
  
  Notes
  -----
  Is it really necessary?
  Monopole etc. could just inherit from GenericPoint directly...
  It is, provided VolumeSR is attribute of PointSR, not GenericPoint
  (even if we rename it to GenericVolume?). And we DO NEED VolumeSR 
  which has unique methods such as split and bounce_off.
  
  """
  def spread_factors(self, r_hicks, **kwargs):
    """
    This wrapper-class probably!  can't be avoided as it transmits
    info between children and parent. (needs different name than spread)
    """
    self.__log.info('r_hicks ' + str(r_hicks))
    self.vol = self.spread(r_hicks)

  # -----------------------------------------------------------------------------
  
  def spread(self, *args, **kwargs):
    """
    args are just passed from children to parent.
    this is especially useful here since the funcs
    used by the parent
    are defined in children under the hood.
    that's why we couldn't spread an instance
    of PointSR without providing them.
    
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
      self.__log.debug('Found %s %s-nodes' % (len(arr), attr))
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
    
    ATTENTION 
    all the calculation here is done in extended-grid coordinates.
    But PointSR.spread assumes regular (not extended) grid, so we 
    need to convert it back.
    
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
      x, y, z, amp = vo
      if abs(amp) < clip_at: # skip tiny factors for speed and compactness of volume
        self.__log.debug('Skipping outside node %s with small ampl %s' % (str([x,y,z]), str(amp)))
        continue
      
      Gs = aghs[((aghs[:,0] == x) & (aghs[:,1] == y) & (aghs[:,2] == z))]
      
      if len(Gs) == 1:
        G = Gs[0]
      elif len(Gs) > 1:
        raise ValueError('More than one ghost found: ' + str(Gs))
      else:
        self.__log.warn('No ghost found for the outside node %s. Skipping it.' % str([x,y,z]))
        continue
      
      G = aghs[((aghs[:,0] == x) & (aghs[:,1] == y) & (aghs[:,2] == z))][0]
      I = np.array(iss[ghs.index(list(G))])
      G = G[:3]
      
      #R[:3] = 2 * I - G
      #R[3] = -amp
      

      R = np.zeros(3) # WE HAVE TO CREATE IT EVERY TIME
      R = 2 * I - G # FIND A REFLECTION OF G WITH RESPECT TO I

      # NOTE return to not-extended-grid coordinates
      R[0] -= elef
      R[1] -= efro
      R[2] -= etop
      
      # NOTE my theory is we should always inject secondary (reflected) sources  
      # as monopoles regardless of the type of the primary source
      R = Monopole(R)
      R.value = -amp # flip the polarity as we reflect from a free surface (it will change for other boundaries!)
      
      reflected.append(R)
      #print('O', vo)
      #print('G', G)
      #print('I', I)
      #print('R', R)
      #print()   
    return reflected #np.array(reflected)

  # -----------------------------------------------------------------------------
  
  def plot(self, **kwargs):
    kwargs['slice_at'] = kw('slice_at', 'y', kwargs)
    kwargs['node'] = kw('node', self.shape[1]//2, kwargs)
    super().plot(**kwargs)

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
    func2 = lambda x : kaiser(x, r, dipole=True) * dsinc_dx(x)
    
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

