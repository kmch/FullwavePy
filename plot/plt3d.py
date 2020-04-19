"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
#from ipywidgets import interact, interactive, fixed, interact_manual


from fullwavepy.generic.parse import kw, del_kw
#from fullwavepy.generic.decor import widgets


@traced
@logged
def figure3d(figsize=(10,10), **kwargs):
  """
  Create and set parameters of a 3d figure. 
  
  Parameters
  ----------
      
  Returns
  -------
  fig : figure 
    Current figure.
  ax : axes
    Current axes.

  Notes
  -----
  Sets axes' limits as extremes of each array.
  
  For some reason adding anything to ax
  invalidates zlimits...
  
  """
  from mpl_toolkits.mplot3d import Axes3D

  azim = kw('azim', None, kwargs)
  elev = kw('elev', None, kwargs)
  aspect = kw('aspect', 'auto', kwargs)
  
  fig = plt.figure(figsize=figsize)
  ax = fig.add_subplot(111, projection='3d')
  
  ax.view_init(azim=azim, elev=elev)  
  #ax.set_aspect(aspect)
  
  # SCALE THE AXES
  #ax.auto_scale_xyz([np.min(XX), np.max(XX)], 
  #                  [np.min(YY), np.max(YY)], 
  #                  [np.min(ZZ), np.max(ZZ)]) 
  #
  ## LABEL THE AXES
  #ax.set_xlabel(xlabel)
  #ax.set_ylabel(ylabel)
  #ax.set_zlabel(zlabel)
  #
  ## FLIP Z AXIS
  #if zflip:
  #  ax.invert_zaxis()
  return fig, ax


# ------------------------------------------------------------------------------
    

