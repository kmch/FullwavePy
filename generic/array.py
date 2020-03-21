"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
from matplotlib.gridspec import GridSpec

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.parse import kw
from fullwavepy.plot.generic import new_figure
from fullwavepy.plot.twod import plot_image


# -------------------------------------------------------------------------------


@traced
@logged
class Arr(np.ndarray):
  """
  Wrapper around numpy's array.
  
  """

  # -----------------------------------------------------------------------------
  
  def __new__(cls, source, **kwargs):
    """
    Init by reading from source.
    
    Notes
    -----
    From https://docs.scipy.org/doc/numpy/user/basics.subclassing.html:
    Input array is an already formed ndarray instance

    """
    source = cls.read(source, **kwargs)
    
    obj = np.asarray(source).view(cls) # CAST THE TYPE
    obj = cls._set_extent(obj, **kwargs)
    
    return obj # NECESSARY!

  # -----------------------------------------------------------------------------   
  
  def _set_extent(obj, **kwargs):    
    if 'extent' in kwargs:
      obj.extent = kwargs['extent']
    else:
      obj.extent = []
      for dim in obj.shape:
        obj.extent.append([0, dim-1])
    return obj

  # -----------------------------------------------------------------------------

  def __array_finalize__(self, obj):
    if obj is None: return
  
  # -----------------------------------------------------------------------------  
  
  def read(source, **kwargs):
    """
    """
    if (type(source) == type(np.array([])) or 
        type(source) == Arr or
        type(source) == Arr1d or
        type(source) == Arr2d or
        type(source) == Arr3d or
        type(source) == WigglyData or
        type(source) == Surf or
        type(source) == SurfFunc or
        type(source) == np.memmap):
      A = source    

    elif isinstance(source, str):
      from fullwavepy.ioapi.generic import read_any
      A = read_any(source, **kwargs)

    else:
     raise TypeError('Arguments need to be either ' + 
                     'file-names or arrays or np.memmap, NOT: %s' %
                     type(source))
    return A

  # -----------------------------------------------------------------------------

  def info(self, **kwargs):
    self.__log.info('shape: {}'.format(self.shape))
    self.__log.info('min: {}, max: {}'.format(np.min(self), np.max(self)))

  # -----------------------------------------------------------------------------
  
  #@widgets()
  def compare(self, othe, mode='interleave', **kwargs): #fig, gs=None, widgets=False, 
    if mode == 'interleave':
      A = self.interleave(othe, **kwargs)
      A.plot(**kwargs)
    else:
      raise ValueError(mode)
    
   # -----------------------------------------------------------------------------

  # -----------------------------------------------------------------------------
  
  def compare_subplots(self, **kwargs):
      assert type(self) == type(othe)
      assert self.shape == othe.shape
      
      xlim = kw('xlim', None, kwargs)
      ylim = kw('ylim', None, kwargs)
      
      if widgets:
        figsize = (kw('figsize_x', 8, kwargs), kw('figsize_y', 8, kwargs))
        fig = plt.figure(figsize=figsize)
        kwargs['widgets'] = False
      
      if gs is None:
        gs = fig.add_gridspec(1,2)    
      
      
      ax1 = fig.add_subplot(gs[0,0])
      self.plot(**kwargs)
      ax2 = fig.add_subplot(gs[0,1])
      othe.plot(**kwargs)
      
      for ax in [ax1, ax2]:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
  

# -------------------------------------------------------------------------------


@traced
@logged
class Arr1d(Arr):
  """
  """
  
  # -----------------------------------------------------------------------------

  def plot(self, **kwargs):
    """
    """
    from fullwavepy.plot.oned import plot_1d
    assert len(self.shape) == 1 # FIXME: SHOULD PUT IT IN INIT
    plt.plot(self)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class Arr2d(Arr):
  """
  
  """
  def interleave(self, othe, **kwargs):
    self.interleaved = interleave_arrays(self, othe, **kwargs)
    return self.interleaved
  
  # -----------------------------------------------------------------------------

  def plot(self, wiggle=False, **kwargs):
    """
    """
    from fullwavepy.plot.twod import plot_image, plot_wiggl
    
    #kwargs['extent'] = np.ravel(self.extent) # ravel JUST IN CASE
    
    # IT SHOULDN'T BE APPLIED TWICE!
    self = modify_array(self, **kwargs)
    
    if wiggle:
      plot_wiggl(self, **kwargs)
    else:
      plot_image(self, **kwargs)     
  
  # -----------------------------------------------------------------------------

  #def compare(self, othe, **kwargs):
    #A = self.interleave(othe, **kwargs)
    #A.plot(**kwargs)
    

# -------------------------------------------------------------------------------


@traced
@logged
class Arr3d(Arr):
  """
  3D array.
  
  """
  @widgets('slice_at', 'node')
  def slice(self, slice_at='y', node=0, widgets=False, **kwargs):
    """
    """
    di = {'x': 0, 'y': 1, 'z': 2} # TRANSLATE slice_at INTO AXIS NO.
    axis = di[slice_at]
    A = Arr2d(np.take(self, indices=node, axis=axis))
    
    extent2d = [el for i, el in enumerate(self.extent) if i != di[slice_at]]
    self.__log.debug('extent2d: ' + str(extent2d))
    A.extent = np.ravel(extent2d)
    return A
  
  # -----------------------------------------------------------------------------
  
  #@widgets('chunk_size')
  def interleave(self, othe, *args, **kwargs):
    A1 = self.slice(*args, **kwargs)
    A2 = othe.slice(*args, **kwargs)
    A = Arr2d(interleave_arrays(A1, A2, **kwargs))
    return A

  # -----------------------------------------------------------------------------
  
  @widgets('cmap', 'slice_at', 'node')
  def plot_slice(self, slice_at='y', node=0, widgets=False, **kwargs):
    """
    """
    arr2d = self.slice(slice_at, node, widgets=False, **kwargs)
    arr2d.plot(**kwargs)
    if slice_at == 'z':
      plt.gca().invert_yaxis()
  
  # -----------------------------------------------------------------------------
  
  #@widgets('cmap', 'slice', 'x', 'y', 'z')
  def plot_3slices(self, fig, gs=None, widgets=False, **kwargs):
    """
    """
    from fullwavepy.plot.twod import plot_image
    
    # LABELS FOR EACH AXIS
    s2 = kw('slice', 'y', kwargs) # MAIN SLICE PLOTTED AT THE BOTTOM IN FULL WIDTH
    s0, s1 = [i for i in ['x', 'y', 'z'] if i != s2]
    s = [s0, s1, s2]
    # CONVERT THE LABELS INTO ARRAY DIMENSIONS (AXES)
    convert_s2a = {'x': 0, 'y': 1, 'z': 2} # TRANSLATE slice TO axis
 
    #if widgets: #FIXME BOILERPLATE
      #figsize = (kw('figsize_x', 8, kwargs), kw('figsize_y', 8, kwargs))
      #fig = plt.figure(figsize=figsize)
    
    if gs is None:
      gs = GridSpec(2,2)
      #gs = fig.add_gridspec(2,2)

    
    if widgets: #or fig is None:
      fig = new_figure(**kwargs)
      gs = fig.add_gridspec(2,2)
     

   
    axes = list(np.zeros(3))
    axes[0] = fig.add_subplot(gs[0,0])
    axes[1] = fig.add_subplot(gs[0,1])
    axes[2] = fig.add_subplot(gs[1,:]) 
    
    kwargs['vmin'] = kw('vmin', np.min(self), kwargs)
    kwargs['vmax'] = kw('vmax', np.max(self), kwargs)
    self.__log.debug('Setting vmin, vmax to: {}, {}'.format(kwargs['vmin'], 
                                                            kwargs['vmax']))
    kwargs['widgets'] = False
    self.__log.debug('Disabling widgets in inner functions.')
    
    
    for i, ax in enumerate(axes):
      plt.sca(ax)
      plot_image(np.take(self, kwargs[s[i]], convert_s2a[s[i]]), **kwargs)
      
      # PLOT SLICING LINES
      a, b = [j for j in ['x', 'y', 'z'] if j != s[i]]
      abcissae_horiz = range(self.shape[convert_s2a[a]])
      ordinate_horiz = np.full(len(abcissae_horiz), kwargs[b])
      ordinate_verti = range(self.shape[convert_s2a[b]])
      abcissae_verti = np.full(len(ordinate_verti), kwargs[a])
      
      if s[i] == 'z':
        abcissae_horiz, ordinate_horiz, abcissae_verti, ordinate_verti = abcissae_verti, ordinate_verti, abcissae_horiz, ordinate_horiz
        ax.invert_yaxis()
      plt.plot(abcissae_horiz, ordinate_horiz, '--', c='white')
      plt.plot(abcissae_verti, ordinate_verti, '--', c='white')
    
    #return ax1, ax2, ax3
    
  # -----------------------------------------------------------------------------  

  def plot(self, nslices=1, *args, **kwargs):
    """
    Note, it doesn't need to have @widgets!
    
    """
    if nslices == 1:
      self.plot_slice(*args, **kwargs)
    elif nslices == 3:
      self.plot_3slices(*args, **kwargs)
    else:
      raise ValueError('Wrong value of nslices: %s' %str(nslices))

  # -----------------------------------------------------------------------------
  
  def scroll(self, **kwargs):
    """
    
    """
    import matplotlib.pyplot as plt
    from fullwavepy.plot.events import IndexTracker
    
    A = self.read(scoord=None)
    
    fig, ax = plt.subplots(1, 1)
    tracker = IndexTracker(ax, A, **kwargs)
    return fig, ax, tracker

  def scrollall(self, fig, **kwargs):
    """
    
    """
    from fullwavepy.plot.events import IndexTrackerAll
    
    A = self.read(scoord=None)
    
    tracker = IndexTrackerAll(fig, A, **kwargs)
    return tracker
    #return tracker.onscroll

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class WigglyData(Arr3d):
  def interleave(self, othe, **kwargs):
    return super().interleave(othe, slice_at='y', node=0, **kwargs)

  def compare(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) #'twilight_shifted'
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    super().compare(*args, **kwargs)
  
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) #'twilight_shifted'
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    super().plot(*args, **kwargs)
  

# -------------------------------------------------------------------------------


@traced
@logged
class Surf(Arr3d):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class SurfFunc(Surf):
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
class SurfParam(Surf):
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


@traced
@logged
class Metadata(object):
  def plotly(self, **kwargs):
    pass







@traced
@logged
def interleave_arrays(A1, A2, **kwargs):
  """ 
  Create an array composed of 
  interleaved arrays Z1 & Z2.
  
  Parameters
  ----------
  
  A1, A2 : arrays
    2D arrays to interleave.
  
  **kwargs : keyword arguments (optional)  
   - chunk_size : int 
     No. of columns of 1 array
     before being proceeded 
     by 2nd array etc.
    
  Returns
  -------
  Z : array
    2D array.
  
  Notes
  -----
  
  """
  chunk_size = kw('chunk_size', 10, kwargs)    
  
  assert len(A1.shape) == 2
  
  if A1.shape != A2.shape:
    raise ValueError('Arrays must have same shapes.')
  
  A = np.array(A1)
  ncols = A.shape[0]
  
  if ncols < 2 * chunk_size:
    interleave_arrays._log.warn('No. of columns=' + str(ncols) + 
           ' < 2 * chunk_size! Outputting empty array')
    return []
  
  nchunks = ncols // chunk_size // 2
  
  for i, Ai in enumerate([A1, A2]):
    i_start = i * chunk_size
    for j in range(nchunks):
      i1 = i_start + j * 2 * chunk_size
      i2 = i_start + j * 2 * chunk_size + (chunk_size) # IT USED TO BE WRONG (-1)
      A[i1 : i2] = Ai[i1 : i2]

  A = np.array(A)

  return A


# -------------------------------------------------------------------------------


@traced
@logged
def slice_array(A, **kwargs):
  """
  
  Parameters
  ----------
  A : array 
   3D array, although other shapes
   should be handled too (work in progress).
   
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  2D or 1D array (special case).
  
  Notes
  -----
  
  """
  scoord = kw('scoord', 'y', kwargs)
  dims = A.shape

  # 1D ARRAY CONTAINING A SINGLE 1D TIME-SERIES
  if len(dims) == 1:
    slice_array._log.warn('1D array detected. Passing it on intact')
    return A
  
  # 3D ARRAY CONTAINING A SINGLE 1D TIME-SERIES
  if (dims[0] == 1) and (dims[1] == 1):
    return A[0][0]
  
  # SURFACES
  if (dims[-1] == 1):
    return A[..., 0]
  
  if scoord == 'x':
    svalue = kw('svalue', len(A)//2, kwargs)
    A = A[svalue]
  
  elif scoord == 'y':
    svalue = kw('svalue', len(A[0])//2, kwargs)
    A = [i[svalue] for i in A]
  
  elif scoord == 'z':
    svalue = kw('svalue', len(A[0][0])//2, kwargs)
    A = [[j[svalue] for j in i] for i in A]
  
  elif scoord is None:
    slice_array._log.warn('No slicing applied')
  
  else:
    raise ValueError('Wrong slice coord: ' + scoord)
  
  if scoord is not None:
    slice_array._log.info('Sliced at svalue=' + str(svalue) + ' node ' + 
                          'of scoord=' + scoord + ' axis')
  
  A = np.array(A)
  
  return A
  

# -------------------------------------------------------------------------------


@timer
@traced
@logged
def modify_array(A, *args, **kwargs):
  """
  Modify each trace (last dimension) 
  of a 1D/2D/3D array using a list 
  of functions.
  
  Parameters
  ----------
  A : array 
   1D/2D/3D array.
   
  **kwargs : keyword arguments (optional)
    Current capabilities:
    modifiers : list
      List of functions to apply subsequently 
      on each trace. The order is following:
      [func1, func2, ...] 
      first func1 will be applied and 
      followed by func2, etc. 
      Note that it is different from
      the composite function's notation:
      ...(func2(func1(trace))  
      Modifiers are set up in a separate function 
      for cleanliness.
      Modifiers are allowed to have *args and **kwargs
      so lambda functions are not recommended as 
      modifiers.

  Notes
  -----
  Always modify trace-wise where trace 
  is the last dimension of the array.
  
  """
  array_modifiers = _set_array_modifiers(**kwargs)
  tracewise_modifiers = _set_tracewise_modifiers(**kwargs)

  A = np.array(A)

  for func in array_modifiers:
    A = func(A, *args, **kwargs)

  for func in tracewise_modifiers:
    A = np.apply_along_axis(func, -1, A, *args, **kwargs)
  
  return A


# -------------------------------------------------------------------------------


@traced
@logged
def _set_array_modifiers(**kwargs):
  """
  
  """
  #from ..signal.su import su_process
  
  modifiers = kw('array_modifiers', [], kwargs)  
  
  clip = kw('clip', None, kwargs)
  clip_min = kw('clip_min', None, kwargs)  
  clip_max = kw('clip_max', None, kwargs)  
  func = kw('func', None, kwargs)
  
  if clip is not None or clip_min is not None or clip_max is not None:
    modifiers.append(clip_array)
    
  #if func is not None:
    #modifiers.append(su_process)
  
  return modifiers


# -------------------------------------------------------------------------------


@traced
@logged
def _set_tracewise_modifiers(**kwargs):
  """
  Set a list of functions to modify 
  a trace / an array of traces.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
    modifiers : list
      List of functions to apply subsequently 
      on each trace. The order is following:
      [func1, func2, ...] 
      first func1 will be applied and so on.
      Note that order of the elements is 
      opposite to the composite function's 
      notation:
      ...(func2(func1(trace))
  
  Returns
  -------
  modifiers : list 
    List of modifiers.
  
  Notes
  -----
  The order matters, they don't commute in general.
  
  We could use lambda functions, but we want to 
  pass **kwargs to modifiers, and it is bad to 
  define lambda functions with *args, **kwargs.
  
  Clipping is done before normalization.
  
  """
  from fullwavepy.generic.math import derivative, dft
  from fullwavepy.generic.math import normalize
  
  modifiers = kw('tracewise_modifiers', [], kwargs)
  norm = kw('norm', None, kwargs)
  spect = kw('spect', None, kwargs)
  deriv = kw('deriv', None, kwargs)
  
  # DERIVATIVE 
  if deriv is not None:
    modifiers.append(derivative)   
  
  # DISCRETE FOURIER TRANSFORM
  if spect is not None:
    modifiers.append(dft)
  
  # NORMALIZATION
  if norm is not None:
    modifiers.append(normalize)
  
  return modifiers


# -------------------------------------------------------------------------------


@traced
@logged
def clip_array(A, clip=None, **kwargs):
  """
  clip : float 
    Convenience to define both bounds 
    at once as [-clip, clip]
    
  
  """
  clip_min = kw('clip_min', None, kwargs)
  clip_max = kw('clip_max', None, kwargs)

  if clip is not None:
    clip_min = -clip
    clip_max = clip

  return np.clip(A, clip_min, clip_max)


# -------------------------------------------------------------------------------


@traced
@logged
def tseries2array(tseries, **kwargs):
  """
  Convert 1d time series into
  a 3d array (fw3d format).
  
  """
  A = np.zeros((1, 1, len(tseries)))
  A[0][0] = tseries
  return A


# -------------------------------------------------------------------------------


@traced
@logged
def list2str(li, **kwargs): #FIXME: MOVE SOMWHERE ELSE
  """
  For SU.
  
  """
  s = ''
  for i in range(len(li)):
    s += str(li[i]) + ','
    
  s = s[:-1] # REMOVE THE TRAILING COMMA
  return s


# ------------------------------------------------------------------------------

