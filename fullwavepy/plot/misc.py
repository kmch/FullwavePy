"""
Misc plots. See fullwavepy.utils for more.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw

@traced
@logged
def time_freq(arg, kwargs_t={}, kwargs_f={}, **kwargs):
  """
  Plot file in both domains side by side.
  
  Should work both for gathers and single traces.
  
  """  
  from .generic import plot
  
  title = kw('title', None, kwargs)
  
  nx, ny = 1, 2
  plt.subplots(nx,ny)
  plt.suptitle(title)

  plt.subplot(nx,ny,1)
  kwargs_t['cmap'] = 'seismic' 
  plot(arg, **kwargs_t, **kwargs)
  xlim = kw('xlim', None, kwargs_t)
  ylim = kw('ylim', None, kwargs_t)
  plt.xlim(xlim)
  plt.ylim(ylim)

  plt.subplot(nx,ny,2)
  kwargs_f['cmap'] = kw('cmap', 'hot', kwargs_f)
  kwargs_f['norm'] = kw('norm', 'rms', kwargs_f)
  kwargs_f['center_cmap'] = kw('center_cmap', False, kwargs_f)
  plot(arg, spect='ampl', **kwargs_f, **kwargs)
  xlim = kw('xlim', None, kwargs_f)
  ylim = kw('ylim', None, kwargs_f)
  plt.xlim(xlim)
  plt.ylim(ylim)    
    

# ------------------------------------------------------------------------------

# alias
def plot_box(*args, **kwargs):
  plot_square(*args, **kwargs)

@traced
@logged
def plot_square(x1, x2, y1, y2, **kwargs):
  """
  Plot a 2D box (square) given its 2 vertices
  (along the diagonal).
  
  Parameters
  ----------
  x1 : float 
    Min. value of X-coord.
  x2 : float 
    Max. value of X-coord.
  y1 : float 
    Min. value of Y-coord.
  y2 : float 
    Max. value of Y-coord.
  **kwargs : keyword arguments (optional)
    Current capabilities:
    - c - color
    - ls - line style 
    - lw - line width
    - alpha - opacity
  
  Returns
  -------
  None
  
  """
  ls = kw('ls', '--', kwargs)
  c = kw('c', 'r', kwargs)
  lw = kw('lw', 2, kwargs)
  alpha = kw('alpha', 0.6, kwargs)
  label = kw('label', None, kwargs)
  
  plot_square._log.debug('(x1,x2,y1,y2) = ({},{},{},{})'.format(x1,x2,y1,y2))
  plt.plot([x1, x2], [y1, y1], ls, c=c, lw=lw, alpha=alpha, label=label)
  plt.plot([x1, x2], [y2, y2], ls, c=c, lw=lw, alpha=alpha)
  plt.plot([x1, x1], [y1, y2], ls, c=c, lw=lw, alpha=alpha)
  plt.plot([x2, x2], [y1, y2], ls, c=c, lw=lw, alpha=alpha)    


  
# ------------------------------------------------------------------------------
  
  
  
