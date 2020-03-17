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


# -------------------------------------------------------------------------------
# PLOT 2D
# -------------------------------------------------------------------------------


@traced
@logged
def compare_2d(A1, A2, **kwargs):
  """
  
  """
  mode = kw('mode', 'interleave', kwargs)
  
  if mode == 'interleave':
    from fullwavepy.generic.array import interleave_arrays, Arr3d
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
  
  ax = kw('ax', plt.gca(), kwargs)
  cmap = kw('cmap', 'twilight', kwargs)
  ncolors = kw('ncolors', None, kwargs)  
  vmin = kw('vmin', np.min(image), kwargs)
  vmax = kw('vmax', np.max(image), kwargs)
  extent = kw('extent', None, kwargs)
  
  if isinstance(cmap, list):
    cmap = _combine_2_cmaps(cmap)
  
  cmap = plt.cm.get_cmap(cmap, ncolors)  
  
  
  #if widgets:# or fig is None:
    #fig = new_figure(**kwargs)

  #if gs is None:
    #gs = fig.add_gridspec(1,1)
    
  if center_cmap:
    vmin, vmax = _center_around_zero(vmin, vmax)
  
  #ax = fig.add_subplot()
  im = ax.imshow(image.T, cmap=cmap, extent=extent, 
                  vmin=vmin, vmax=vmax)
  if cbar:
    colorbar(im, ax)


# -------------------------------------------------------------------------------


@traced
@logged
def plot_wiggl(image, **kwargs): #NOTE
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
  

# ------------------------------------------------------------------------------
# FORMAT 2D
# ------------------------------------------------------------------------------


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
    







@traced
@logged
def plot_image_OLD(image, **kwargs): # NOTE
  """
  
  Parameters
  ----------
  
  
  **kwargs : keyword arguments (optional)
    Current capabilities: 
    
    cmap : str 
      Add a trailing '_r' to get 
      a reversed version.
      Type plt.colormaps() to list all 
      registered ones.
      
      Sequential (perceptually uniform):
        'inferno', 'magma', 'plasma', 'viridis'
      Sequential (selected):
        'hot;, 'Greys', 'YlGnBu', 'GnBu', 'PuBuGn'
      Diverging:
        'seismic', 'bwr', 'RdBu', RdYlBu', 'Spectra'
      cmocean:
        'cmo.topo' etc.     
      Other:
        'ocean', 'gist_earth', 'cubehelix'
      
      
    shade : bool
      Default: False.
    
    Parsed only if shade=True:
    - shade_azim : float 
      Azimuth of the light source [deg].
      Default: 315.
    - shade_elev : float 
      Elevation of the light source [deg]. 
      Default: 45.
    - shade_mode : str
      Possible values: 'Hillshade', 'hsv', 'overlay', 'soft'.
      Default: 'overlay'.
    - shade_exag : float 
      Vertical exaggaration of shaded topography.
      Default: 10.
      
  Returns
  -------
  None
  
  Notes
  -----
  Mapping data onto colors using a colormap typically involves two steps: 
  a data array is first mapped onto the range 0-1 using a subclass of Normalize
  
  """  
  ncolors = kw('ncolors', None, kwargs)
  shade = kw('shade', False, kwargs)
  cbar = kw('cbar', False, kwargs)
  contour = kw('contour', False, kwargs)
  darken = kw('darken', False, kwargs)
  
  # DETECT CONFLICTS
  if shade and contour:
    raise ValueError('You cannot combine shading and a contoufullwavepy.plot.')
  if darken and contour:
    raise ValueError('You cannot combine darkening and a contoufullwavepy.plot.')  
  if darken and shade:
    raise ValueError('You cannot combine darkening and shading.')  
  
  kwargs['minn'], kwargs['maxx'] = np.min(image), np.max(image)
  imshow_kwargs = _set_cmap(**kwargs)
  imageT = image.T
  
  
  # SET CORRECT Y-AXIS FOR SPECTRAL PLOTS
  if 'spect' in kwargs:
    from fullwavepy.generic.math import dft_freqs
    nsamps, ntraces = imageT.shape
    y = dft_freqs(nsamps, which='positive', **kwargs)
    imageT = np.array(imageT[ :len(y)])
    imshow_kwargs['extent'] = [0, ntraces, y[-1], 0] 
    
    
  # THIS IS ONLY FOR COLORBAR TO DISPLAY VALUES BEFORE NORMALIZATION
  im_raw = plt.imshow(imageT, **imshow_kwargs)  
  ax = plt.gca()
  
  if shade:
    from matplotlib.colors import LightSource
    shade_azim = kw('shade_azim', 315, kwargs)
    shade_elev = kw('shade_elev', 45, kwargs)
    shade_mode = kw('shade_mode', 'overlay', kwargs)
    shade_exag = kw('shade_exag', 10, kwargs)
    ls = LightSource(azdeg=shade_azim, altdeg=shade_elev)
    imageT = ls.shade(image.T, **imshow_kwargs, 
                      blend_mode=shade_mode, vert_exag=shade_exag)
    
  if darken:
    from matplotlib.colors import Normalize
    dark_min = kw('dark_min', minn, kwargs)
    dark_max = kw('dark_max', maxx, kwargs)
    # CREATE AN ARRAY (nx,ny,nz,3) FULL OF VALUE 70 (GREY COLOR)
    # NOTE: '*' UNPACKS A LIST THAT FOLLOWS
    greys = np.full((*image.T.shape, 3), 70)
    alphas = np.array(image.T)
    # SET OPACITY TO 1 IF dark_min <= value <= dark_max 
    # AND TO 0.5 (HALF-TRANSPARENT) OTHERWISE 
    # IT'S DONE WEIRDLY TO MAKE SURE WE REPLACE VALUES CORRECTLY
    unique = 12345.12345
    alphas[np.where((alphas >= dark_min) & (alphas <= dark_max))] = unique
    alphas[alphas != unique] = 0.5
    alphas[alphas == unique] = 1.0 
    # THIS CONVERTS THE NP.ARRAY TO AN RGB IMAGE ADDING
    # AN EXTRA DIMENSION OF SIZE 4 (FOR R,G,B & ALPHA)
    cmap = imshow_kwargs['cmap']
    imageT = cmap(Normalize(imshow_kwargs['vmin'], imshow_kwargs['vmax'])(image.T))
    imageT[..., -1] = alphas # ELLPISIS IS EQUIVALENT TO :,:,: HERE
    plt.imshow(greys)
  
  # NOTE
  if 'extent' in kwargs:
    imshow_kwargs['extent'] = kwargs['extent']
    
  plt.imshow(imageT, **imshow_kwargs)  
  
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im_raw, cax=cax) 
    if ncolors:
      ticks = np.linspace(minn, maxx, ncolors+1)
      cbar.set_ticks(ticks)
      cbar.set_ticklabels(ticks)
    plt.sca(ax) # SET CURRENT AXIS BACK TO THE ACTUAL PLOT
    
  if contour:
    ncontours = kw('ncontours', 5, kwargs)
    contour_color = kw('contour_color', 'w', kwargs)
    contour_width = kw('contour_width', 1, kwargs)
    contour_style = kw('contour_style', 'dotted', kwargs)
    cntr = ax.contour(imageT, ncontours, colors=contour_color, 
                      linewidths=contour_width, linestyles=contour_style)  
    ax.clabel(cntr, cntr.levels)  
