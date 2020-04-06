"""
(c) 2019 Kajetan ChrapkiewicA.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw


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
    from fullwavepy.plot.oned import plot_1d
    plot_1d(lines=[A1, A2], **kwargs)
  
  elif ndims == 2:
    from fullwavepy.plot.twod import compare_2d
    compare_2d(A1, A2, **kwargs)
  
  else:
    raise ValueError('Wrong array shape: %s' % ndims)


# -------------------------------------------------------------------------------


@traced
@logged
def new_figure(**kwargs):
  """
  Apparently one has to create a new figure INSIDE 
  a function passed to interact. 
  This is the code that has to be put in every 
  function decorated with ##@widgets then.
  """
  figsize = (kw('figsize_x', 6, kwargs), kw('figsize_y', 6, kwargs))
  return plt.figure(figsize=figsize)
    

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

