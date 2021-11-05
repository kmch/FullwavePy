"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

This module contains handles for collections of similar files e.g. models
for different iterations. 

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.lists.basic import IterFileList, ShotFileList, JobFileList


@traced
@logged
class SchedFileList(IterFileList):
  """
  Files dumped by the scheduler, e.g. shared by all slaves
  not to mention timesteps.
  
  """
  def read(self, **kwargs):
    for f in self.it[1: ]: # SKIP ITER. 0
      try:
        f.read(**kwargs)
      except FileNotFoundError as err:
        self.__log.warning(err)    
  
  # -----------------------------------------------------------------------------

  def plot_all(self, **kwargs):
    """
    """
    suffix = kw('suffix', '', kwargs)
    
    for it in self.it:
      if it is None:
        continue
      title = it.fname
      self.__log.info('Plotting ' + it.fname)
      #plt.title(title)
      iplot(**kwargs)
      plt.savefig(strip(title)+suffix)
      plt.close()
  
  # ----------------------------------------------------------------------------- 

  ##@widgets('iter')
  # def plot(self, **kwargs):
  #   pass

  # -----------------------------------------------------------------------------
  
  #def plot(self, imin=0, imax=1, istep=1, **kwargs):
  #  """
  #  """
  #  its = self.it[imin:imax:istep]
  #  fnames = [i.fname for i in its]
  #  #return fullwavepy.ploti(fnames, **kwargs)
  #  raise NotImplementedError('syntax after sed')

  # -----------------------------------------------------------------------------
@traced
@logged
class SlaveFileList(SchedFileList, ShotFileList):
  """
  """
  def __init__(self, proj, **kwargs):
    """
    """
    super().__init__(proj, **kwargs)
    self.sids = self._read_sids(**kwargs)

  # -----------------------------------------------------------------------------    

  def read(self, **kwargs):
    for it in self.it[1: ]: # SKIP ITER. 0
      for f in it.values():
        try:
          f.read(**kwargs)
        except FileNotFoundError as err:
          self.__log.warning(err)    
  
  # -----------------------------------------------------------------------------    

  # def plot(self, **kwargs):
  #   for it in self.it[1: ]: # SKIP ITER. 0
  #     for f in it.values():
  #       try:
  #         f.plot(**kwargs)
  #         plt.figure()
  #       except FileNotFoundError as err:
  #         self.__log.warning(err)

  # -----------------------------------------------------------------------------    
@traced
@logged
class LogFileList(JobFileList):
  pass
@traced
@logged
class SubmitFileList(JobFileList):
  pass


