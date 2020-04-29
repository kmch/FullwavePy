"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
from matplotlib.gridspec import GridSpec

from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.decor import widgets, timer
from fullwavepy.plot.generic import figure
from fullwavepy.plot.plt2d import plot_image

# NOTE
# Arr3d -> Arr2d etc. can only be achieved by their slice() methods
# Numpy's slice-index notation A[:,0,:] etc. works (i.e. reshapes)
# but doesn't convert the type



@traced
@logged
class Arr(np.ndarray):
  """
  Wrapper around numpy's array.
  
  """
  def __new__(cls, source, ndims=None, **kwargs):
    """
    Init by reading from source.
    
    Notes
    -----
    From https://docs.scipy.org/doc/numpy/user/basics.subclassing.html:
    Input array is an already formed ndarray instance

    """
    source = cls._read(source, **kwargs)
    
    obj = np.asarray(source).view(cls) # CAST THE TYPE
    obj = cls._set_extent(obj, **kwargs)
    #obj = cls._set_coords(obj, **kwargs)
    obj = cls._set_dx(obj, **kwargs)

    if ndims is not None:
      assert len(obj.shape) == ndims
    
    return obj # NECESSARY!

  # -----------------------------------------------------------------------------
 
  def _read(source, **kwargs):
    """
    """
    #from fullwavepy.seismic.data import Data
    #from fullwavepy.ndat.manifs import Surf, SurfZ, Plane
    #from fullwavepy.ndat.points import Points
    #
    #if (type(source) == type(np.array([])) or 
    #    type(source) == Arr or
    #    type(source) == Arr1d or
    #    type(source) == Arr2d or
    #    type(source) == Arr3d or
    #    type(source) == Data or
    #    type(source) == Surf or
    #    type(source) == SurfZ or
    #    type(source) == Plane or
    #    type(source) == Points or
    #    type(source) == np.memmap):
    #  A = source    

    if isinstance(source, str):
      from fullwavepy.ioapi.generic import read_any
      if hasattr(source, 'shape'): # FOR EFFICIENCY (SEE read_any)
        kwargs['shape'] = self.shape
      A = read_any(source, **kwargs)
    else:
      A = source
      
    #else:
    # raise TypeError('Arguments need to be either ' + 
    #                 'file-names or arrays or np.memmap, NOT: %s' %
    #                 type(source))
    return A
  
  # -----------------------------------------------------------------------------   
  
  def _set_extent(obj, func=None, **kwargs):    
    if 'extent' in kwargs:
      obj.__log.debug('Using extent from kwargs, even if it means overwriting')
      obj.extent = kwargs['extent']
    
    elif hasattr(obj, 'extent'):
      obj.__log.debug('obj.extent already set and not provided in kwargs')
      pass
    
    else:
      obj.__log.debug('Setting extent to default.')
      obj.extent = obj._default_extent(func, **kwargs)
    
    return obj
  
  # -----------------------------------------------------------------------------  
  
  def _default_extent(obj, func=None, **kwargs):
    """
    Redefined in child classes to account for vertical axis flipping 
    when plotting with imshow.
    
    """
    if func is None:
      func = lambda dim : [0, dim-1]   # outdated: # NOT dim-1; SEE GridProjFile ETC.
    extent = []
    for dim in obj.shape:
      extent.append(func(dim))
    
    # if len(obj.shape) == 1:
    #   extent = extent[0]
    
    return extent
  
  # -----------------------------------------------------------------------------

  def _set_dx(obj, **kwargs):
    """
    It is fully determined by extent and shape.
    In general, it is axis-dependent (dx != dy != dz != dx)
    """
    dx = []
    obj.__log.debug('obj.shape %s' % str(obj.shape))
    obj.__log.debug('obj.extent %s' % str(obj.extent))
    assert len(obj.shape) == len(obj.extent)
    for nx, (x1, x2) in zip(obj.shape, obj.extent):
      obj.__log.debug('nx=%s, x1=%s, x2=%s' % (nx, x1, x2))
      dx_1D = (x2 - x1) / (nx-1) if nx > 1 else None
      obj.__log.debug('dx_1D=%s' % dx_1D)
      dx.append(dx_1D)
    
    obj.dx = np.array(dx)
    return obj

  # -----------------------------------------------------------------------------
  
  def _set_coords(obj, **kwargs):
    obj.__log.debug('obj.extent' + str(obj.extent))
    obj.__log.debug('Setting coords to None. Fill it with actual code')
    obj.coords = None
    
    return obj
  
  #def _set_shape(obj, shape=None, **kwargs):
    #self.shape = None
   
  # -----------------------------------------------------------------------------

  def __array_finalize__(self, obj):
    if obj is None: return
  
  # -----------------------------------------------------------------------------  

  def _metre2index(self, m, axis, **kwargs):
    origin = self.extent[axis][0]
    i = (m - origin) / self.dx[axis]
    if not i.is_integer():
      raise ValueError('Index must be integer not %s' % i)
    return int(i)

  # -----------------------------------------------------------------------------  
  
  def _metre2gridnode(self, *args, **kwargs):
    return self._metre2index(*args, **kwargs) + 1  
  
  # -----------------------------------------------------------------------------

  def _box2inds(self, box, **kwargs):
    """
    Convert box into slicing-indices using extent.
    
    """
    box = np.array(box)
    extent = np.array(self.extent)
    assert len(box.shape) == 1
    assert len(box) == len(extent.flatten())
    box = box.reshape(extent.shape)
    inds = np.zeros(box.shape)
    for axis, _ in enumerate(box):
      b0, b1 = box[axis]
      inds[axis][0] = self._metre2index(b0, axis)
      inds[axis][1] = self._metre2index(b1, axis) + 1 # NOTE: FOR np.arange(b1, b2) etc.
      self.__log.debug('axis %s: i1=%s, i2=%s' % (axis, inds[axis][0], inds[axis][1]))    
    return inds.astype(int)

  # ----------------------------------------------------------------------------- 

  def carve(self, box, **kwargs):
    """
    Carve a box out of an array.

    """
    inds = self._box2inds(box, **kwargs)
    
    for axis in range(len(self.shape)):
      self = np.take(self, np.arange(*inds[axis]), axis=axis)
    
    self.extent = np.array(box).reshape(inds.shape)
    return self
  
  # -----------------------------------------------------------------------------  
 
  def save(self, fname, **kwargs):
    from fullwavepy.ioapi.fw3d import save_vtr
    save_vtr(self, fname)
  
  # -----------------------------------------------------------------------------

  def info(self, **kwargs):
    self.__log.info('shape: {}'.format(self.shape))
    self.__log.info('min: {}, max: {}'.format(np.min(self), np.max(self)))

  # -----------------------------------------------------------------------------
  
  ###@widgets()
  def compare(self, othe, mode='interleave', **kwargs): #fig, gs=None, widgets=False, 
    if mode == 'interleave':
      A = self.interleave(othe, **kwargs)
      A.plot(**kwargs)
    else:
      raise ValueError(mode)
    
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
  def __new__(cls, source, **kwargs):
    return super().__new__(cls, source, ndims=1, **kwargs)
    
  # -----------------------------------------------------------------------------

  def plot(self, **kwargs):
    """
    format of extent: [[x1,x2]] is for compatibility with 2d and 3d
    """
    c = kw('c', None, kwargs)

    assert np.array(self.extent).shape == (1,2)
    self.__log.debug('self.extent' + str(self.extent))
    
    x1, x2 = self.extent[0] 
    x = np.linspace(x1, x2, len(self))
    
    plt.plot(x, self, c=c)
    return plt.gca()

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class Arr2d(Arr):
  """
  
  """
  def __new__(cls, source, **kwargs):
    return super().__new__(cls, source, ndims=2, **kwargs)

  # -----------------------------------------------------------------------------

  ###@widgets('slice_at', 'node')
  def slice(self, slice_at='y', node=0, widgets=False, **kwargs):
    """
    """
    di = {'x': 0, 'y': 1} # TRANSLATE slice_at INTO AXIS NO.
    axis = di[slice_at]
    
    
    A = Arr1d(np.take(self, indices=node, axis=axis))
    
    assert len(self.extent) == 2
    extent1d = np.array([el for i, el in enumerate(self.extent) if i != di[slice_at]])
    self.__log.debug('extent1d %s' % str(extent1d))
    
    A.extent = extent1d
    # A = Arr1d(A)
    
    return A

  # -----------------------------------------------------------------------------
                    
  def interleave(self, othe, **kwargs):
    self.interleaved = interleave_arrays(self, othe, **kwargs)
    return self.interleaved
  
  # -----------------------------------------------------------------------------

  ###@widgets('cmap', 'slice_at', 'node')
  def plot_slice(self, slice_at='y', node=0, widgets=False, **kwargs):
    """
    """
    arr1d = self.slice(slice_at, node, widgets=False, **kwargs)
    ax = arr1d.plot(**kwargs)
    return ax

  # -----------------------------------------------------------------------------

  def plot_full(self, wiggle=False, **kwargs):
    """
    """
    kwargs['extent'] = np.ravel(self.extent) # ravel JUST IN CASE
    
    # IT SHOULDN'T BE APPLIED TWICE!
    self = modify_array(self, **kwargs)
    
    if wiggle:
      ax = plot_wiggl(self, **kwargs)
    else:
      ax = plot_image(self, **kwargs) 
    
    return ax
  
  # -----------------------------------------------------------------------------
  
  def plot(self, *args, **kwargs):
    """
    """
    if 'slice_at' in kwargs:
      ax = self.plot_slice(*args, **kwargs)
    else:
      ax = self.plot_full(*args, **kwargs)
    return ax
  
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
  def __new__(cls, source, **kwargs):
    return super().__new__(cls, source, ndims=3, **kwargs)
  
  # -----------------------------------------------------------------------------
  
  ###@widgets('slice_at', 'node')
  def slice(self, slice_at='y', node=0, widgets=False, **kwargs):
    """
    """
    di = {'x': 0, 'y': 1, 'z': 2} # TRANSLATE slice_at INTO AXIS NO.
    axis = di[slice_at]
    A = Arr2d(np.take(self, indices=node, axis=axis))
    
    assert len(self.extent) == 3
    # extent2d = np.ravel([el for i, el in enumerate(self.extent) if i != di[slice_at]])
    extent2d = np.array([el for i, el in enumerate(self.extent) if i != di[slice_at]])

    # if axis != 2:
    self.__log.debug('Setting extent2d so that no vertical-axis flipping is needed.')
    self.__log.debug('NOW ALSO FOR zslice (NOT TESTED BUT SEEMS TO HAVE FIXED THE BUG)')
    # extent2d[-2: ] = [extent2d[-1], extent2d[-2]]
    extent2d[-1] = extent2d[-1][::-1]
    self.__log.debug('extent2d: ' + str(extent2d))
    
    A.extent = extent2d
    
    return A
  
  # -----------------------------------------------------------------------------
  
  ###@widgets('chunk_size')
  def interleave(self, othe, *args, **kwargs):
    A1 = self.slice(*args, **kwargs)
    A2 = othe.slice(*args, **kwargs)
    A = Arr2d(interleave_arrays(A1, A2, **kwargs))
    return A

  # -----------------------------------------------------------------------------
  
  ##@widgets('cmap', 'slice_at', 'node')
  def plot_slice(self, slice_at='y', node=None, widgets=False, **kwargs):
    """
    """
    nx, ny, nz = self.shape
    if node is None:
      if slice_at == 'x':
        node = kw('node', nx//2, kwargs)
      elif slice_at == 'y':
        node = kw('node', ny//2, kwargs)
      elif slice_at == 'z':
        node = kw('node', nz//2, kwargs)      
      else:
        raise ValueError('Wrong slice_at: %s' % str(slice_at))

    arr2d = self.slice(slice_at, node, widgets=False, **kwargs)
    kwargs['title'] = 'slice at %s=%s' % (slice_at, node)
    del_kw('slice_at', kwargs) # JUST IN CASE
    
    ax = arr2d.plot(**kwargs)
  
    if slice_at == 'z': # DISABLE?
      ax.invert_yaxis()
    return ax
  
  # -----------------------------------------------------------------------------

  def plot_3slices(self, x, y, z, fig=None, gs=None, **kwargs):
    """
    """
    from fullwavepy.plot.plt2d import plot_image
    
    if fig is None:
      fig = figure(16,8)
    
    kwargs['x'] = x
    kwargs['y'] = y
    kwargs['z'] = z

    # LABELS FOR EACH AXIS
    s2 = kw('slice', 'y', kwargs) # MAIN SLICE PLOTTED AT THE BOTTOM IN FULL WIDTH
    s0, s1 = [i for i in ['x', 'y', 'z'] if i != s2]
    s = [s0, s1, s2]
    # CONVERT THE LABELS INTO ARRAY DIMENSIONS (AXES)
    convert_s2a = {'x': 0, 'y': 1, 'z': 2} # TRANSLATE slice TO axis
 
    if gs is None:
      gs = GridSpec(2,2)
   
    axes = list(np.zeros(3))
    axes[0] = fig.add_subplot(gs[0,0])
    axes[1] = fig.add_subplot(gs[0,1])
    axes[2] = fig.add_subplot(gs[1,:]) 
    
    kwargs['vmin'] = kw('vmin', np.min(self), kwargs)
    kwargs['vmax'] = kw('vmax', np.max(self), kwargs)
    self.__log.debug('Setting vmin, vmax to: {}, {}'.format(kwargs['vmin'], 
                                                            kwargs['vmax']))
    
    for i, ax in enumerate(axes):
      plt.sca(ax)
      aaxx = plot_image(np.take(self, kwargs[s[i]], convert_s2a[s[i]]), **kwargs)
      
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
  
  def p3s(self, *args, **kwargs):
    return self.plot_3slices(*args, **kwargs)
  
  # -----------------------------------------------------------------------------

  ###@widgets('cmap', 'slice', 'x', 'y', 'z')
  def plot_3slices_OLDish(self, fig=None, gs=None, widgets=False, **kwargs):
    """
    """
    from fullwavepy.plot.plt2d import plot_image
    
    if fig is None:
      fig = figure(16,8)
    
    kwargs['x'] = kw('x', 0, kwargs)
    kwargs['y'] = kw('y', 0, kwargs)
    kwargs['z'] = kw('z', 0, kwargs)

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
      fig = figure(**kwargs)
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

  def plot(self, *args, **kwargs):
    """
    Framework plotter.
    
    Notes
    -----
    This is a preferred function to call rather than
    plot_3slices directly. This is because plot 
    formatting is set in subclasses by overwriting
    plot method. This could be avoided by defining
    _format_plot() method or similar.

    Note, it doesn't need to have ##@widgets!
    
    """
    if not ('x' in kwargs or 'y' in kwargs or 'z' in kwargs):
      nslices = 1
    elif 'x' in kwargs and not ('y' in kwargs or 'z' in kwargs):
      nslices = 1
      kwargs['slice_at'] = 'x'
      kwargs['node'] = kwargs['x']
    elif 'y' in kwargs and not ('x' in kwargs or 'z' in kwargs):
      nslices = 1
      kwargs['slice_at'] = 'y'
      kwargs['node'] = kwargs['y']
    elif 'z' in kwargs and not ('x' in kwargs or 'y' in kwargs):
      nslices = 1
      kwargs['slice_at'] = 'z'
      kwargs['node'] = kwargs['z']  
    elif 'x' in kwargs and 'y' in kwargs and 'z' in kwargs:
      nslices = 3
    else:
      raise ValueError('Slicing arguments not understood.')
    
    if nslices == 1:
      self.plot_slice(*args, **kwargs)
    elif nslices == 3:
      self.plot_3slices(*args, **kwargs)
    else:
      raise ValueError('Wrong value of nslices: %s' %str(nslices))
    
    return plt.gca()

  # -----------------------------------------------------------------------------
  
  def scroll(self, **kwargs):
    """
    
    """
    import matplotlib.pyplot as plt
    from fullwavepy.plot.events import IndexTracker
    
    fig, ax = plt.subplots(1, 1)
    tracker = IndexTracker(ax, self, **kwargs)
    return fig, ax, tracker

  def scrollall(self, fig, **kwargs):
    """
    To make it work in a jupyter notebook:
     %matplotlib notebook
     %matplotlib notebook
     fig = plt.figure(figsize=(5,20))
     tracker = some_array.scrollall(fig, cmap='viridis')
     fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    
    """
    from fullwavepy.plot.events import IndexTrackerAll
    
    tracker = IndexTrackerAll(fig, self, **kwargs)
    return tracker
    #return tracker.onscroll

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


#@traced
#@logged
#class Grid(Arr3d):
  #pass
    #def _extent(self, **kwargs):
    #"""
    #"""
    #box = self.proj.box
    #extent = np.array([box[ :2], box[2:4], box[4: ]])
    ## CENTER VOXELS AT INTEGERS! ONLY IF unit='node'! FIXME
    #extent += 0.5 times dx possibly
    
    #self.__log.debug('extent: %s' % str(extent))    
    #return extent


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
  from fullwavepy.math.generic import normalize
  from fullwavepy.math.operators import derivative
  from fullwavepy.math.fourier import dft
  
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
  
