"""
(c) 2019-2021 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import time # to help find bottle-necks
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
from ipywidgets import Dropdown, IntSlider

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr3d
from fullwavepy.ndat.points import GenericPoint, Points3d
from fullwavepy.numeric.funcs import kaiser, sinc, dsinc_dx


# FIXME sync it somehow with hick_sources.f90
clip_at = 1e-9 # spread-factors smaller than this will be neglected


@traced
@logged
def vp2rho(vp, **kwargs):
  """
  Convert P-wave velocity into density using 
  Gardner relation of the same form as in Fullwave3D.
  
  Parameters
  ----------
  vp : float
    P-wave velocity in m/s.

  Returns
  -------
  rho : float
    Density in kg/m3.

  Notes
  -----
  Units seem OK. E.g. vp=3000 m/s gives 2294 kg/m3
  which looks like a physical density.

  """
  # Velocity range in which Gardner's relation is valid
  vp_min = 1600
  vp_max = 5000
  
  if (vp < vp_min):
    vp2rho._log.warning("vp=%s < vp_min=%s of Gardner's relation, setting water density")
    vp2rho._log.warning("Sync vp_min with Fullwave3D's Gardner cutoff etc.")
    rho = 1000
  elif (vp > vp_max):
    raise ValueError("Gardner relation not valid for this vp=%s > vp_max=%s" % (vp, vp_max))
  else:
    rho = 310 * (vp**0.25)
  
  return rho


# -------------------------------------------------------------------------------


@traced
@logged
class SRs(Points3d):
  """
  li : list
    List of [id, xyz] records, where xyz = (x,y,z)
    
  """
  def __init__(self, li, **kwargs):
    self.__log.debug('li: %s' % str(li))
    self.li = []
    for l in li:
      ID, xyz = l
      self.li.append(PointSR(xyz, ID=ID))
    
    # cls.__log.debug('type(new_li[0]) %s' % type(new_li[0]))
    # print(new_li)
    # return super().__new__(cls, new_li)
    # return new_li
    # super().__init__(new_li, **kwargs)

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
    
    assert len(srtype_ids) == len(self.li)
    # for srtype_id, sr[k, v] in zip(srtype_ids, self.items()):
    new_list = []
    for srtype_id, sr in zip(srtype_ids, self.li):
      clss = mapp[srtype_id]
      self.__log.debug('Appending instance of %s' % str(clss))
      new_list.append(clss(sr, ID=sr.ID, **kwargs))

    self.__log.debug('new_list: %s' % str(new_list))
    self.li = new_list
  
  # ----------------------------------------------------------------------------- 
 
  def spread_factors(self, srtype_ids, *args, **kwargs):
    """
    """
    self.set_type(srtype_ids, **kwargs)
    #self.sprd_fctrs = {}
    self.hyper = {}
    nsr = len(self.li)
    i = 1
    # for srid, sr in self.items():
    for sr in self.li:
      srid = sr.ID
      self.__log.info('ID %s (%s/%s)' % (srid, i, nsr))
      self.hyper[srid] = HyperPointSR(sr, ID=srid, **kwargs)
      self.hyper[srid].spread_factors(*args, **kwargs)
      i += 1
  
  # -----------------------------------------------------------------------------
  
  def qc_sprd_fctrs(self, **kwargs):
    pass
  
  def widgets_qc_sprd_fctrs(self, **kwargs):
    snap_max = 1
    raise ValueError(snap_max)
    widgets = dict(srid=Dropdown(options=[i.ID for i in self.li]),
                   snap=IntSlider(value=snap_max, min=0, max=snap_max, step=1),
                   slice_at=Dropdown(options=['y', 'x', 'z']),
                   node=IntSlider(value=rmax, min=0, max=2*rmax-1, step=1))
    return widgets


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
  FIXME: implement check if we're not out of bounds
  
  """
  def __init__(self, pointsr, **kwargs):
    """
    """
    self.pointsr = pointsr
    self.r_hicks = kw('r_hicks', 3, kwargs)

  # -----------------------------------------------------------------------------
  
  @timer
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
  
  
  # FIXME: should not depend on ghs, etc. in some cases
  @timer
  def spread_factors(self, ghs, iss, ine, elef, efro, etop, **kwargs):
    """
    Iterative 
    
    Notes
    -----
    We pass kwargs to allow for 

    """
    assert isinstance(ghs, np.ndarray)

    self.find_vol(**kwargs)
    self.vol[..., -1] = 0.0 # OTHERWISE GROWING IN IPYNB
    
    psrs_to_spread = [self.pointsr]
    i_max = 2
    i = 0
    self.snapshots = []
    while i <= i_max:
      n_to_spread = len(psrs_to_spread)
      self.__log.info('It. %s/%s: %s points to spread...' % (i, i_max, n_to_spread)) 

      if n_to_spread == 0:
        self.__log.info('No more points to spread. Returning.')
        return
      
      self.__log.debug('psrs_to_spread' + str(psrs_to_spread))
      new_psrs_to_spread = []
      for psr_to_spread in psrs_to_spread:
        self.__log.debug('psr_to_spread ' + str(psr_to_spread))
        
        del_kw('r_hicks', kwargs)
        psr_to_spread.spread_factors(r_hicks=self.r_hicks, **kwargs)
        psr_to_spread.vol.split(ine, elef, efro, etop, **kwargs)
        
        self._update_vol(np.copy(psr_to_spread.vol.ins), **kwargs)
        
        
        reflected = psr_to_spread.vol.bounce_off(ghs, iss, elef, efro, etop, **kwargs)
        if len(reflected) > 0:
          new_psrs_to_spread += reflected # append would create a nested list (undesired)
      
      self.snapshots.append(np.copy(self.vol))
      psrs_to_spread = list(new_psrs_to_spread)
      i += 1

  # -----------------------------------------------------------------------------
  
  @timer
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
  @timer
  def spread_factors(self, r_hicks, **kwargs):
    """
    This wrapper-class probably!  can't be avoided as it transmits
    info between children and parent. (needs different name than spread)
    """
    self.vol = self.spread(r_hicks)

  # -----------------------------------------------------------------------------
  
  @timer
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
  @timer
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
  
  @timer
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
  
  @timer
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
  
  @timer
  def bounce_off(self, ghs, iss, elef, efro, etop, **kwargs):
    """
    ghs : list of ghosts
    iss : list of corresponding intersects (same len)
    
    ATTENTION 
    all the calculation here is done in extended-grid coordinates.
    But PointSR.spread assumes regular (not extended) grid, so we 
    need to convert it back.
    
    """
    assert isinstance(ghs, np.ndarray) # lists are terribly slow and np.array(list)
    # for each point is even worse an idea. This was a 10x bottle-neck
    assert len(ghs) == len(iss)
    
    # t1 = time.perf_counter()
    vout = np.copy(self.out)
    # t2 = time.perf_counter()
    # self.__log.info('Copying self.out to vout took %s s' % "{:15.12f}".format(t2-t1))
    
    # t1 = time.perf_counter()
    vout[:, 0] += elef
    vout[:, 1] += efro
    vout[:, 2] += etop
    # t2 = time.perf_counter()
    # self.__log.info('Adding extra nodes to vout took %s s' % "{:15.12f}".format(t2-t1))    
    # here we compare grid-coords, NOT grid-coords to array-indices
    # => we don't subtract 1
    

    # THIS A BOTTLENECK!!!
    # t1 = time.perf_counter()
    # aghs = np.array(ghs)
    # t2 = time.perf_counter()
    # self.__log.info('aghs = np.array(ghs) took %s s' % "{:15.12f}".format(t2-t1))    

    reflected = []
    for vo in vout:
      x, y, z, amp = vo
      if abs(amp) < clip_at: # skip tiny factors for speed and compactness of volume
        self.__log.debug('Skipping outside node %s with small ampl %s' % (str([x,y,z]), str(amp)))
        continue
      
      # t1 = time.perf_counter()
      Gs = ghs[((ghs[:,0] == x) & (ghs[:,1] == y) & (ghs[:,2] == z))]
      # t2 = time.perf_counter()
      # self.__log.info('Locating the ghost took %s s' % "{:15.12f}".format(t2-t1))       
      
      if len(Gs) == 1:
        G = Gs[0]
      elif len(Gs) > 1:
        raise ValueError('More than one ghost found: ' + str(Gs))
      else:
        self.__log.warning('No ghost found for the outside node %s. Skipping it.' % str([x,y,z]))
        continue
      
      I = np.array(iss[np.where(np.all(ghs==G, axis=1))[0][0]])
      # I = np.array(iss[ghs.index(list(G))]) # this was the bottle-neck (see above)
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
    """
    Notes
    -----
    Particle velocity is proportional 
    to the NEGATIVE pressure gradient (v_i ~ -dp/dx_i),
    hence the minus sign in func2.

    """
    func1 = lambda x : kaiser(x, r) * sinc(x) 
    func2 = lambda x : kaiser(x, r, dipole=True) * dsinc_dx(x) * (-1)
    
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
def xyz2w(xyz, extended_dims, **kwargs):
  """
  """
  x, y, z = xyz
  nx, ny, nz = extended_dims
  return int((x - 1) * ny * nz + (y - 1) * nz + z)


# -------------------------------------------------------------------------------


@traced
@logged
def w2xyz(w, extended_dims, **kwargs):
  """
  For QC.
  
  FIXME: it's terribly slow so now only for 
  single-value checks. We could use Adrian's 
  formulas to vectorize it.
  
  """
  assert float(w).is_integer()
  assert w > 0
  enx, eny, enz = extended_dims
  
  if w > enx * eny * enz:
    raise ValueError('w > enx * eny * enz')


  i = 1
  for x in range(1, enx+1):
    for y in range(1, eny+1):
      for z in range(1, enz+1):
        if i == w:
          return x, y, z
        i += 1
  raise ValueError('Could not find x,y,z for w=%s' % w)


# -------------------------------------------------------------------------------

