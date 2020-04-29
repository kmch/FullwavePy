"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr3d


@traced
@logged
class GenericPoint(np.ndarray):
  """
  """
  def __new__(cls, xyz, value=1, **kwargs):
    cls.value = value
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
    
    Returns
    vol : Arr3d
      3D array of values.
      It is endowed with vol.coords and vol.extent attributes.
      The former stores coordinates.

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
    extent = [] # WHAT IS THIS FOR?! FIXME
    for i, func in enumerate(funcs):
      # WRAPPER TO ACT ON A SINGLE COORDINATE OF AN ND-POINT
      func_of_xyz = lambda point : func(point[i])
      vol *= np.apply_along_axis(func_of_xyz, coord_axis, dists)
    
    
    #nshape = np.array(cube.shape)
    #nshape[-1] += 1 # INCREASE THE LAST DIM TO INCLUDE VALUE -> (x,y,z,val)
    #
    #self.vol = np.zeros(nshape)
    #self.vol[..., :-1] = cube
    #self.vol[..., -1] = vol
    
    self.__log.debug('Scaling the cube with self.value=%s' % str(self.value))
    self.vol = Arr3d(vol) * self.value
    self.vol.extent = self.cube_extent(cube, **kwargs)
    #self.vol._set_coords()
    self.__log.debug('self.vol.extent' + str(self.vol.extent))
    self.vol.coords = cube
    return self.vol

  # -----------------------------------------------------------------------------
  
  def cube_extent(self, cube, **kwargs):
    if len(cube.shape) == 4:
      x1, y1, z1 = cube[0,0,0]
      self.__log.debug('Modified  extent of cube - tmp? DOUBLE-CHECK')
      x2, y2, z2 = cube[-1,-1,-1] # + 1
      extent = [[x1, x2], [y1, y2], [z1, z2]]

    else:
      raise ValueError('cube.shape ' + str(cube.shape))
    
    self.__log.debug('returning extent %s' % str(extent))    
    return extent
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class Points(dict):
  """
  dict stores ids (keys) and coords (values)
  ids can annotate plots
  
  __new__ is not necessary to redefine, unlike nd.array
  
  """  
  #def __new__(cls, source, ndims=None, **kwargs):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class Points3d(Points):
  """
  """    
  def slice(self, slice_at='y', **kwargs):
    if slice_at == 'x':
      i1, i2 = 1, 2
    elif slice_at == 'y':
      i1, i2 = 0, 2
    elif slice_at == 'z':
      i1, i2 = 0, 1
    else:
      raise ValueError('Wrong slice coord: %s' % slice_at)      
    for key, val in self.items():
      #assert len(val) == 3 # IT CAN HAVE METADATA
      self[key] = np.array([val[i1], val[i2]])

  # ----------------------------------------------------------------------------- 

  def plot_slice(self, ax=None, **kwargs):
    """
    """
    annotate = kw('annotate', False, kwargs)
    annoffset = kw('annoffset', 0, kwargs)
    alpha = kw('alpha', 0.7, kwargs)
    marker = kw('marker', '.', kwargs)
    markersize = kw('markersize', 5, kwargs)
    markeredgecolor = kw('markeredgecolor', 'k', kwargs)
    markerfacecolor = kw('markerfacecolor', 'none', kwargs) # EMPTY MARKERS
    if ax is None:
      ax = plt.gca()
    
    self.slice(**kwargs)
    
    if annotate: 
      for key, val in self.items():
        ax.annotate(key, (val[0]+annoffset, val[1]+annoffset), clip_on=True) # clip_on IS REQUIRED
    
    ax.plot([i[0] for i in self.values()], [i[1] for i in self.values()], 
            '.',
            alpha=alpha, 
            marker=marker, 
            markersize=markersize, 
            markeredgecolor=markeredgecolor,
            markerfacecolor=markerfacecolor,
           )
  
  # -----------------------------------------------------------------------------   
  
  def plot_3slices(self, fig, **kwargs): # LEGACY
    d = self.read(**dict(kwargs, unit='node'))
    
    s3 = kw('slice', 'y', kwargs) #FIXME: THIS MUST BE MERGED WITH arr3d
    s1, s2 = [i for i in ['x', 'y', 'z'] if i != s3]
    s = [s1, s2, s3]
    
    for i in range(3):
      self.plot_slice(s[i], fig.axes[i])

  # -----------------------------------------------------------------------------   

  def plot(self, *args, **kwargs):
    #if 'slice_at' in kwargs:
    self.plot_slice(*args, **kwargs)
    #else:
     
  # -----------------------------------------------------------------------------
  
  def plotly(self, fig=None, **kwargs): # LEGACY
    """
    """
    import plotly.graph_objects as go    
    color = kw('color', 'black', kwargs)
    mode = kw('mode', 'markers', kwargs)
    size = kw('size', 2, kwargs)
    
    di = self.read(unit='m')

    if fig is None:
      fig = go.Figure()
    fig.add_trace(go.Scatter(x=[i[0] for i in di.values()], 
                             y=[i[1] for i in di.values()], 
                             text=list(di.keys()), mode=mode,
                             marker=dict(color=color, size=size),
                             line=dict(color=color), showlegend=False))
    
    return fig       
     
  # -----------------------------------------------------------------------------   


# -------------------------------------------------------------------------------

