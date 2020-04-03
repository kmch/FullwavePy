"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import widgets
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.ndat.points import Nodes


@traced
@logged
def NodesIn(Nodes):
  pass


# -------------------------------------------------------------------------------
    
  
@traced
@logged
class Ghosts(Nodes):
  """
  for fancier subarray-ing check out np.ix_
  """
  def split_lvls(self, nlvls, **kwargs):
    ilvl = 4 # FIXME MAKE IT GLOBAL IN class Ghost(Node)
    self.lvl = {}
    for lvl in range(1, nlvls+1):
      self.lvl[lvl] = self[self[:,ilvl] == lvl]
      
  #def plot(self, **kwargs):
    
  
# -------------------------------------------------------------------------------



