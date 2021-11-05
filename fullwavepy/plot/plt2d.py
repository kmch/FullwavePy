"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
import cmocean
#from ipywidgets import interact, interactive, fixed, interact_manual


from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.decor import widgets
from fullwavepy.plot.generic import set_xlabels


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
@logged
##@widgets('cmap')
def plot_image(image, widgets=False, center_cmap=False, cbar=True, **kwargs):
  """
  Wrapper around plt.imshow. 
  
  Parameters
  ----------
  image : 2D array
  
  Notes
  -----
  
  """
  from fullwavepy.plot.generic import figure
  from matplotlib.colors import LogNorm, LightSource
  
  ax = kw('ax', plt.gca(), kwargs)
  aspect = kw('aspect', 'auto', kwargs)
  aspect = 'auto' if aspect == 'a' else aspect
  aspect = 'equal' if aspect == 'e' else aspect
  title = kw('title', None, kwargs)
  spect = kw('spect', None, kwargs)
  cmap = kw('cmap', 'Greys', kwargs)
  ncolors = kw('ncolors', None, kwargs)  
  vmin = kw('vmin', np.min(image), kwargs)
  vmax = kw('vmax', np.max(image), kwargs)
  extent = kw('extent', None, kwargs)
  noextent = kw('noextent', False, kwargs)
  if noextent: # useful for QC of extent
    extent = None  
  shade = kw('shade', False, kwargs)
  alpha = kw('alpha', 1, kwargs)
  lognorm = kw('lognorm', False, kwargs)
  norm = LogNorm() if lognorm else None
  if lognorm:
    plot_image._log.warning('lognorm on => setting vmin,vmax to None and center_cmap to False')
    vmin = None
    vmax = None
    center_cmap = False
    
  if spect is not None:
    from fullwavepy.numeric.fourier import dft_freqs
    ntraces, nsamps = image.shape
    plot_image._log.debug('nsamps=%s, ntraces=%s' % (nsamps, ntraces))
    y = dft_freqs(nsamps, which='positive', **kwargs)
    image = np.array(image[:, :len(y)])
    extent = [0, ntraces, y[-1], 0] 
  

  # THIS CENTERS CELLS AT INTEGERS AGAIN WHICH WAS OVERWRITTEN BY CUSTOM EXTENT
  if extent is not None:
    plot_image._log.debug('extent before: %s' % extent)
    extent = np.array(extent) - .5
    plot_image._log.debug('extent after: %s' % extent)
    # APPROACH BELOW IS WRONG BECAUSE IT STRETCHES THE CELLS
    #x1, x2, y1, y2 = extent
    #x1 -= .5
    #x2 += .5
    #y1 -= .5
    #y2 += .5 
    #extent = [x1, x2, y1, y2]

  if isinstance(cmap, list):
    cmap = _combine_2_cmaps(cmap)
  
  cmap = plt.cm.get_cmap(cmap, ncolors)  
  
  #if widgets:# or fig is None:
    #fig = figure(**kwargs)

  #if gs is None:
    #gs = fig.add_gridspec(1,1)
    
  if center_cmap:
    vmin, vmax = _center_around_zero(vmin, vmax)
  
  if title is not None:
    ax.set_title(title)
  
  if shade:
    raise NotImplementedError('Debug first!')
    ls = LightSource(azdeg=kw('azdeg', 45, kwargs), \
      altdeg=kw('altdeg', 45, kwargs))
    image = ls.shade(image, cmap=cmap, \
       blend_mode=kw('blend_mode', 'hsv', kwargs), 
       vert_exag=kw('vert_exag', 50, kwargs))  
  
  #ax = fig.add_subplot()
  # im = ax.imshow(image.T, cmap=cmap, extent=extent,
  im = ax.imshow(image.T, cmap=cmap, extent=extent, 
                 vmin=vmin, vmax=vmax, norm=norm,
                 alpha=alpha)
  if cbar:
    colorbar(im, ax)
  
  # xlabels = kw('xlabels', None, kwargs)
  # if 'xlabels' is not None:
  #   set_xlabels(xlabels, **kwargs)

  ax.set_aspect(aspect)
  return ax
@logged
def plot_wiggl(image, **kwargs):
  """
  Plot wiggly traces.

  Parameters
  ----------
  image : 2D array  
  
  Notes
  -----
  It should be merged wit plot_trace.
  Taking the spectrum should be separated from plot_trace 
  or called from external.
  """
  assert len(image.shape) == 2
  from fullwavepy.plot.plt1d import plot_1d
  gap = kwargs.get('gap', 10)
  t = np.arange(image.shape[-1])
  for i, trace in enumerate(image):
    trace += i * gap
    zero_axis = np.ones(len(t)) * i * gap
    plot_1d(lines=[trace], **kwargs)
  return plt.gca()
@logged
def colorbar(imshow_object, ax, pos='right', size='3%', pad=0.2, **kwargs):
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  divider = make_axes_locatable(ax)
  cax = divider.append_axes(pos, size, pad)
  cbar = plt.colorbar(imshow_object, cax=cax) 
  plt.sca(ax)
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
# @logged
def cat_cmaps(cmaps, vmin, vmax):
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
  
  n = 100
  ratio = np.abs(vmin)/np.abs(vmax)
  print('ratio',ratio)
  colors1 = cmap1(np.linspace(0., 1, int(ratio*n)))
  colors2 = cmap2(np.linspace(0, 1, n))
  colors = np.vstack((colors1, colors2))
  my_cmap = LinearSegmentedColormap.from_list('my_cmap', colors)

  return my_cmap
    

