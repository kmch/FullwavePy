"""
(c) 2019 Kajetan ChrapkiewicA.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw


def save(fname, **kwargs):
  """
  Save current figure.

  Notes
  -----
  save_as overwrites fname + suffix.
  
  """
  return NotImplementedError('Debug')
  # save_as = kw('save_as', None, kwargs)
  # suffix = kw('suffix', None, kwargs)
  
  # if suffix is not None:
  #   ext = exten(fname)
  #   fname = strip(fname) + '_' + suffix + '.' + ext
  #     # print(this_func, fname)
  
  # if not save_as:
  #   save_as = fname     
  
  # return save_fig(save_as, **kwargs)


# -------------------------------------------------------------------------------


def save_fig(fname, **kwargs):
  """
  
  """
  return NotImplementedError('Debug')
  # core = strip(fname)
  # plt.title(core)
  # fname = core + '.png'
  # plt.savefig(fname, format='png')
  # plt.close()
  
  # return fname

# -------------------------------------------------------------------------------

@traced
@logged
def compare(arg1, arg2, **kwargs):
  """
  Compare two arrays of data.
  
  Parameters
  ----------

  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  None
  
  Notes
  -----
  Data-pieces interact with each other. 
  They need to be known at the same time
  e.g. to color-fill the area between them.
  
  """
  from fullwavepy.ndat.arrays import Arr3d
  
  A1 = Arr3d(arg1)
  A2 = Arr3d(arg2)
  
  if A1.shape != A2.shape:
    raise ValueError('Arrays must have identical shapes.')
  
  ndims = len(A1.shape)  
  
  if ndims == 1:
    from fullwavepy.plot.plt1d import plot_1d
    plot_1d(lines=[A1, A2], **kwargs)
  
  elif ndims == 2:
    from fullwavepy.plot.plt2d import compare_2d
    compare_2d(A1, A2, **kwargs)
  
  else:
    raise ValueError('Wrong array shape: %s' % ndims)


# -------------------------------------------------------------------------------


@traced
@logged
def figure(figsize_x=6, figsize_y=6, **kwargs):
  """
  Apparently one has to create a new figure INSIDE 
  a function passed to interact. 
  This is the code that has to be put in every 
  function decorated with ##@widgets then.
  """
  figsize = (figsize_x, figsize_y)
  return plt.figure(figsize=figsize)
    

# -------------------------------------------------------------------------------


@traced
@logged
def figax(figsize_x=6, figsize_y=6, **kwargs):
  """
  """
  figsize = (figsize_x, figsize_y)
  fig, ax = plt.subplots(figsize=figsize)
  return fig, ax 


# -------------------------------------------------------------------------------


@traced
@logged
def flipy(ax=plt.gca(), **kwargs):
  ax.invert_yaxis()


# -------------------------------------------------------------------------------

@traced
@logged
def autect(ax=plt.gca(), **kwargs):
  ax.set_aspect('auto')


# -------------------------------------------------------------------------------


@traced
@logged
def aspeqt(ax=plt.gca(), **kwargs):
  ax.set_aspect('equal')


# -------------------------------------------------------------------------------


@traced
@logged    
def set_xlabels(labels, decim_xlabels=10, rotate_xlabels=None, **kwargs):
  """
  Decimate and rotate labels.

  labels : list 
    Labels before decimation
  """  
  locs = np.arange(len(labels))[::decim_xlabels]
  labels = labels[::decim_xlabels]
  rotation = lambda decim : np.clip(90 - 10 * (decim - 1), 0, 90)
  if rotate_xlabels is None:
    rotate_xlabels = rotation(decim_xlabels)
  set_xlabels._log.debug('Rotating xlabels %s degrees' % rotate_xlabels)
  locs, labels = plt.xticks(locs, labels, rotation=rotate_xlabels)
  return locs, labels


# ------------------------------------------------------------------------------
