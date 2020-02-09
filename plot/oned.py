"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from ipywidgets import interact, interactive, fixed, interact_manual

# -------------------------------------------------------------------------------
# PLOT 1D
# -------------------------------------------------------------------------------


@traced
@logged
def plot_1d(**kwargs):
  """
  Plot multiple lines and points.
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  None
  
  Notes
  -----
  Iterators allow for providing fewer labels
  than actual artists (lines/scatts). Only 
  the provided labels will shown in the legend.
  
  If only two lines are provided, the area 
  between them will be filled with color.
  
  """
  lines = kw('lines', [], kwargs)
  scatts = kw('scatts', [], kwargs)
  llabels = iter(kw('llabels', [], kwargs))
  slabels = iter(kw('slabels', [], kwargs))
  scatt_ampl = iter(kw('scatt_ampl', [], kwargs))
  
  line2 = lines[-1] if len(lines) == 2 else 0
  
  for line in lines:
    kwargs['label'] = next(llabels, None)
    plot_line(line, line2, **kwargs)    
  
  for scatt in scatts:
    kwargs['label'] = next(slabels, None)
    kwargs['scatt_ampl'] = next(scatt_ampl, None)
    plot_points(scatt, **kwargs)
    

# ------------------------------------------------------------------------------


@traced
@logged
def plot_line(line, line2=0, **kwargs):
  """
  Plot a line either as y(x) or x(y).
  
  Parameters
  ----------
  line : list 
    List of numbers.
  line2 : list 
    Second line to define a color filling.
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  None
  
  Notes
  -----
  
  """
  orient = kw('orient', 'h', kwargs)
  fill = kw('fill', False, kwargs)
  label = kw('label', None, kwargs)
  c = kw('c', None, kwargs)
  lw = kw('lw', None, kwargs)
  ls = kw('ls', None, kwargs)
  
  if type(line2) == type(0) and line2 == 0: #FIXME
    line2 = np.zeros(len(line))
  
  if 'spect' in kwargs:
    from fullwavepy.generic.math import dft_freqs
    x = dft_freqs(len(line), which='positive', **kwargs)
    line = np.array(line[ :len(x)])
    line2 = np.array(line2[ :len(x)])
  else:
    x = _set_xaxis_1d(line, **kwargs)
  
  if fill: # MUST BE BEFORE SWAPPING x, line
    _fill_area(plt.gca(), x, line, line2, **kwargs)
    
  if orient == 'h': # y(x)
    pass
  elif orient == 'v': # x(y)
    x, line = line, x
  else:
    raise ValueError('Wrong orient: ' + orient)
  
  plt.plot(x, line, label=label, c=c, lw=lw, ls=ls)


# -------------------------------------------------------------------------------


@traced
@logged
def slice_points(points, scoord, **kwargs):
  """
  Project points onto 1 of the planes.

  """
  if scoord == 'x':
    X1, X2 = [[i[1] for i in points], [i[2] for i in points]]
  elif scoord == 'y':
    X1, X2 = [[i[0] for i in points], [i[2] for i in points]]
  elif scoord == 'z':
    X1, X2 = [[i[0] for i in points], [i[1] for i in points]]
  else:
    raise ValueError('Wrong slice coord: %s' % scoord)
  
  points2d = list(zip(X1, X2))
  return points2d


# -------------------------------------------------------------------------------


@traced
@logged
def plot_points(scatt, **kwargs):
  """
  Plot points either as y(x) or x(y).
  
  Parameters
  ----------
  scatt : list 
    List of tuples (x,y).
  
  **kwargs : keyword arguments (optional)
    Current capabilities:
    - scatt_ampl : list 
      List of amplitudes, one per each element
      of scatt. The points will be colored.
  
  Returns
  -------
  None
  
  Notes
  -----
  
  """
  from fullwavepy.plot.twod import _center_around_zero, colorbar
  
  orient = kw('orient', 'h', kwargs)
  label = kw('label', None, kwargs)
  c = kw('scatt_ampl', None, kwargs)
  cmap = kw('scatt_cmap', 'viridis', kwargs)
  s = kw('scatt_size', 5, kwargs) # CAN BE A LIST!
  center_cmap = kw('center_cmap', False, kwargs)
  
  if type(s) == list:
    s = np.array(s)
  s = s ** 2 # SQUARED IS WHAT SCATTER EXPECTS
  vmin = kw('scatt_vmin', None, kwargs)
  vmax = kw('scatt_vmax', None, kwargs)
  cbar = kw('scatt_cbar', False, kwargs)
  
  if orient == 'h': # y(x)
    scatt_x = [i[0] for i in scatt]
    scatt_y = [i[1] for i in scatt]
  elif orient == 'v': # x(y)
    scatt_x = [i[1] for i in scatt]
    scatt_y = [i[0] for i in scatt]
  else:
    raise ValueError('Wrong orient: ' + orient)
  
  #scatter_kwargs = _set_cmap(**kwargs)
  #scatter_kwargs['vmin'] = vmin
  #scatter_kwargs['vmax'] = vmax
  
  #print(scatter_kwargs)

  if center_cmap:
    vmin, vmax = _center_around_zero(vmin, vmax)
  
  im = plt.scatter(scatt_x, scatt_y, label=label, 
                   s=s, c=c, vmin=vmin, vmax=vmax, cmap=cmap)
  
  if center_cmap:
    vmin, vmax = _center_around_zero(vmin, vmax)
  
  #ax = fig.add_subplot()
  #im = ax.imshow(image.T, cmap=cmap, extent=extent, 
                  #vmin=vmin, vmax=vmax)
  if cbar:
    ax = plt.gca()
    colorbar(im, ax)  
  
  #if cbar: # FIXME: MERGE WITfullwavepy.plot_image's
    #ax = plt.gca()
    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #cbar = plt.colorbar(im, cax=cax) 
    #plt.sca(ax) # SET CURRENT AXIS BACK TO THE ACTUAL PLOT
  

# -------------------------------------------------------------------------------
# FORMAT 1D
# -------------------------------------------------------------------------------


@traced
@logged
def _set_xaxis_1d(y, dt=1, **kwargs):
  """
  Define values of the X axis for 1fullwavepy.plots.
  
  Parameters
  ----------
  y : list 
    y = f(x).
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  x : list 
    X axis.
  
  Notes
  -----
  'extent' is for compatibility with 2D plots,
  y1, y2 are dummy here.
  
  """
  if 'extent' in kwargs:
    x1, x2, y1, y2 = kwargs['extent']
    x = np.linspace(x1, x2, len(y))
  else:  
    x = list(np.arange(len(y)) * dt)     
  
  return x


# ------------------------------------------------------------------------------


@traced
@logged
def _fill_area(ax, x, y1, y2=0, **kwargs):
  """
  Define values of the X axis for 1D plots.
  
  Parameters
  ----------
  y : list 
    y = f(x).
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  x : list 
    X axis.
  
  Notes
  -----
  It comes in handy because it allows for
  two different orientations and facecolors.
  Other than that it's just a built-in 
  'fill_between' function.
  
  """  
  orient = kw('orient', 'h', kwargs)
  if type(y2) == type(0):
    fc = kw('fc', 'bw', kwargs)
  else:
    fc = kw('fc', None, kwargs)
  
  if fc == 'bw':
    kws1 = {'facecolors' : 'black' }
    kws2 = {'facecolors' : 'white' }
  else:
    kws1 = {}
    kws2 = {}
    
  if orient == 'h':
    ax.fill_between(x, y1, y2, where=y1>=y2, **kws1)
    ax.fill_between(x, y1, y2, where=y1<y2, **kws2)
  elif orient == 'v':
    ax.fill_betweenx(x, y1, y2, where=y1>=y2, **kws1)
    ax.fill_betweenx(x, y1, y2, where=y1<y2, **kws2)
  else:
    raise ValueError('Wrong orient: ' + orient)  


# ------------------------------------------------------------------------------  


@traced
@logged
def colors(n, cmap='rainbow', **kwargs):
  """
  Create an iterator for rainbow colors.
  
  Parameters
  ----------
  n : int 
    Number of colors in the spectrum.
  
  
  Returns
  -------
  colors : iterator
    c = next(colors)
  
  Notes
  -----
  Usage plot(..., c=next(colors))
  
  """
  from matplotlib.cm import get_cmap
  
  cmap = get_cmap(cmap)
  cols = iter(cmap(np.linspace(0, 1, n)))
  
  return cols


# ------------------------------------------------------------------------------ 

