"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.parse import kw
from ipywidgets import interact #, interactive, fixed, interact_manual


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
    
    Notes
    -----
    From https://docs.scipy.org/doc/numpy/user/basics.subclassing.html:
    Input array is an already formed ndarray instance

    """
    source = cls.read(source, **kwargs)
    
    # CAST THE TYPE
    obj = np.asarray(source).view(cls)
    
    return obj # NECESSARY!

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
  
  # -----------------------------------------------------------------------------

  def plot(self, wiggle=False, **kwargs):
    """
    
    Notes
    -----
    Function passed to interact must take kwargs only.

    We wrap plt.imshow because interact(plt.imshow, X=fixed(self.T)) 
    does not work for some reason.
    
    """
    from fullwavepy.plot.twod import plot_image, plot_wiggl
    
    if wiggle:
      plot_wiggl(self, **kwargs)
    else:
      plot_image(self, **kwargs)     
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class Arr3d(Arr):
  """
  """
  
  # -----------------------------------------------------------------------------
  
  @widgets
  def plot_3slices(self, fig, gs=None, widgets=False, **kwargs):
    from matplotlib.gridspec import GridSpec
    from fullwavepy.plot.twod import plot_image
    
    # LABELS FOR EACH AXIS
    s2 = kw('slice', 'y', kwargs) # MAIN SLICE PLOTTED AT THE BOTTOM IN FULL WIDTH
    s0, s1 = [i for i in ['x', 'y', 'z'] if i != s2]
    s = [s0, s1, s2]
    # CONVERT THE LABELS INTO ARRAY DIMENSIONS (AXES)
    convert_s2a = {'x': 0, 'y': 1, 'z': 2} # TRANSLATE slice TO axis
    
    if gs is None:
      gs = GridSpec(2,2)
    
    if widgets: #FIXME BOILERPLATE
      figsize = (kw('figsize_x', 8, kwargs), kw('figsize_y', 8, kwargs))
      fig = plt.figure(figsize=figsize)
    
    kwargs['vmin'] = kw('vmin', np.min(self), kwargs)
    kwargs['vmax'] = kw('vmax', np.max(self), kwargs)
    self.__log.debug('Setting vmin, vmax to: {}, {}'.format(kwargs['vmin'], 
                                                            kwargs['vmax']))
    
    kwargs['widgets'] = False
    self.__log.debug('Disabling widgets.')
    
    axes = list(np.zeros(3))
    axes[0] = fig.add_subplot(gs[0,0])
    axes[1] = fig.add_subplot(gs[0,1])
    axes[2] = fig.add_subplot(gs[1,:]) 
    
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

    
  def plot1d(self, **kwargs):
    pass
  
  def plot2d(self, axis, **kwargs):
    #axis = kw('axis', 0, kwargs)
    coor = kw('coor', 0, kwargs)
    def _plot2d(**kwargs):
      a = np.take(np.copy(self), indices=kwargs['coor'], axis=axis)
      plt.imshow(a.T)
    
    interact(_plot2d, coor=range(40))
    
  # -----------------------------------------------------------------------------

  def plot(self, svalue=0, **kwargs):
    from fullwavepy.plot.twod import plot_image
    self.__log.warn('scoord=y, svalue=0')
    plot_image(self[:, svalue, :], **kwargs)

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


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
      i2 = i_start + j * 2 * chunk_size + (chunk_size - 1) 
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

