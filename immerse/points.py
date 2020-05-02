"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import widgets
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.arrays import Arr3d, Arr


@traced
@logged
def Ghost(object):
  def __init__(self, x, y, z, w, lvl, flag, **kwargs):
    self.x = x
    self.y = y
    self.z = z
    self.w =  w
    self.lvl = lvl
    self.flag = flag


# -------------------------------------------------------------------------------
    
  
# see Points note about subclassing list
# @traced
# @logged
# class Ghosts(list):
#   """
#   for fancier subarray-ing check out np.ix_
#   """
#   def split_lvls(self, nlvls, **kwargs):
#     arr = np.array(self)
#     ilvl = 4 # FIXME MAKE IT GLOBAL IN class Ghost(Node)
#     self.lvl = {}
#     for lvl in range(1, nlvls+1):
#       self.lvl[lvl] = np.copy(arr[arr[:,ilvl] == lvl])
      
  #def plot(self, **kwargs):
    
  
# -------------------------------------------------------------------------------



