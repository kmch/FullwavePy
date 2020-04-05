"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
import cmocean
#from ipywidgets import interact, interactive, fixed, interact_manual


from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.decor import widgets


@traced
@logged
def compare_2d(A1, A2, **kwargs):
  """
  
  """
  mode = kw('mode', 'interleave', kwargs)
  
  if mode == 'interleave':
    from fullwavepy.ndat.arrays import interleave_arrays, Arr3d
    A = Arr3d(interleave_arrays(A1, A2, **kwargs))
    A.plot(**kwargs)
    #plt.grid(color='black', linestyle='-.', linewidth=0.7)
    #plot_2d(images=[A], **kwargs)
  else:
    raise ValueError('Unknown mode: ' + mode)


# -------------------------------------------------------------------------------  


@traced
@logged
@widgets('cmap')
def plot_image(image, widgets=False, center_cmap=False, cbar=True, **kwargs):
  """
  Wrapper around plt.imshow. 
  
  Parameters
  ----------
  image : 2D array
  
  Notes
  -----
  
  """
  from fullwavepy.plot.generic import new_figure
  from matplotlib.colors import LogNorm
  
  ax = kw('ax', plt.gca(), kwargs)
  title = kw('title', None, kwargs)
  cmap = kw('cmap', 'twilight', kwargs)
  ncolors = kw('ncolors', None, kwargs)  
  vmin = kw('vmin', np.min(image), kwargs)
  vmax = kw('vmax', np.max(image), kwargs)
  extent = kw('extent', None, kwargs)
  alpha = kw('alpha', 1, kwargs)
  lognorm = kw('lognorm', False, kwargs)
  norm = LogNorm() if lognorm else None
  if lognorm:
    plot_image._log.warn('lognorm on => setting vmin,vmax to None and center_cmap to False')
    vmin = None
    vmax = None
    center_cmap = False
    
  # THIS CENTERS CELLS AT INTEGERS AGAIN WHICH WAS OVERWRITTEN BY CUSTOM EXTENT
  if extent is not None:
    extent = np.array(extent) - .5
  
  if isinstance(cmap, list):
    cmap = _combine_2_cmaps(cmap)
  
  cmap = plt.cm.get_cmap(cmap, ncolors)  
  
  
  #if widgets:# or fig is None:
    #fig = new_figure(**kwargs)

  #if gs is None:
    #gs = fig.add_gridspec(1,1)
    
  if center_cmap:
    vmin, vmax = _center_around_zero(vmin, vmax)
  
  if title is not None:
    ax.set_title(title)
  
  #ax = fig.add_subplot()
  im = ax.imshow(image.T, cmap=cmap, extent=extent, 
                 vmin=vmin, vmax=vmax, norm=norm,
                 alpha=alpha)
  if cbar:
    colorbar(im, ax)


# -------------------------------------------------------------------------------


@traced
@logged
def plot_wiggl(image, **kwargs): #FIXME MOVE TO SEISMIC DATA
  """

  Parameters
  ----------
  image : 2D array  
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Notes
  -----
  It should be merged wit plot_trace.
  
  Taking the spectrum should be separated from plot_trace 
  or called from external.
  
  """
  from fullwavepy.plot.oned import plot_1d
  
  gap = kw('gap', 10, kwargs) # GAP BETWEEN TRACE
  
  t = np.arange(image.shape[-1])
  
  #print(image.shape)
  #raise NotImplementedError('Unfinished')
  
  for i, trace in enumerate(image):
    #trace = trace[0]
    
    #print(trace)
    trace += i * gap
    #plt.plot(trace)
    zero_axis = np.ones(len(t)) * i * gap
    plot_1d(lines=[trace], **kwargs)
    #plot_1d(t, zero_axis, trace, orient='verti',
            #c1='k', c2='w', c_line='k', lw=.1, **kwargs)
    
    #plt.plot(trace, c='k')
  #plt.gca().invert_yaxis() # DISABLED SINCE IT IS FLIPPED BY ANOTHER FUNCTION
  

@traced
@logged
def colorbar(imshow_object, ax, pos='right', size='3%', pad=0.2, **kwargs):
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  divider = make_axes_locatable(ax)
  cax = divider.append_axes(pos, size, pad)
  cbar = plt.colorbar(imshow_object, cax=cax) 
  plt.sca(ax)


# ------------------------------------------------------------------------------


@traced
@logged
def _center_around_zero(minn, maxx, **kwargs): #NOTE
  """
  Center a diverging colormap around zero.
  
  Parameters
  ----------
  minn, maxx : float
    Extreme value of the image tfullwavepy.plot.
  
  **kwargs : keyword arguments (optional)
    Current capabilities: 
  
  """  
  # SIGNED ZERO (PLATFORM DEPENDENT) - OTHERWISE WRONG BEHAVIOUR
  if minn == 0.0:
    maxx = -0.0 
          
  if abs(minn) > abs(maxx):
    a = abs(minn)
  else:
    a = abs(maxx)

  vmin = -a 
  vmax = a 
  
  return vmin, vmax


# ------------------------------------------------------------------------------


@traced
@logged
def _combine_2_cmaps(cmaps):
  """
  Combine 2 colormaps.
  
  Parameters
  ----------
  
  **kwargs : keyword arguments (optional)
    Current capabilities: 
  
  """    
  import cmocean
  from matplotlib.colors import LinearSegmentedColormap
  
  if len(cmaps) == 2:
    cmap1, cmap2 = cmaps
  elif len(cmaps) == 0:
   cmap1 = 'cmo.turbid'
   cmap2 = 'cmo.ice_r'    
  else:
    raise ValueError('len(cmaps): %s' % len(cmaps))
  
  cmap1 = plt.cm.get_cmap(cmap1)
  cmap2 = plt.cm.get_cmap(cmap2)
  
  colors1 = cmap1(np.linspace(0., 1, 128))
  colors2 = cmap2(np.linspace(0, 1, 128))
  colors = np.vstack((colors1, colors2))
  my_cmap = LinearSegmentedColormap.from_list('my_cmap', colors)

  return my_cmap


# ------------------------------------------------------------------------------
    

