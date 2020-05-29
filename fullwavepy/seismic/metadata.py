"""
This module defines a seismic field-experiment as a whole.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced
import pandas as pd


from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, path_leave
from fullwavepy.generic.system import bash, get_files

from fullwavepy.ioapi.generic import *
from fullwavepy.ioapi.segy import *

from fullwavepy.seismic.data import DataFileSgy
from fullwavepy.seismic.models import *
from fullwavepy.seismic.surfaces import *


@traced
@logged
class Metadata(pd.DataFrame):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class MetadataFile(CsvFile, TextFile):
  def read(self, overwrite=True, **kwargs):
    """
  
    Notes
    -----
    Overwrite=True by default because otherwise plots are not 
    updated even though they are supposed (e.g. you are passing 
    a different fname). They will be correct (updated) only
    if you delete self.array variable, e.g. by restarting the 
    notebook kernel.
    Disable overwrite only for PERFORMANCE (e.g. interactive plot)
    when the array remains unchanged unlike other (e.g. plotting)
    parameters.
    
    """
    if hasattr(self, 'df') and overwrite == True:
       self.__log.warning('{}.df exists but it will be overwritten (overwrite=True).'.format(type(self)))    
    
    if not (hasattr(self, 'df') and overwrite == False):
      self.df = Metadata(super().read(**kwargs))

    return self.df 


  # -----------------------------------------------------------------------------  
 
 
# -------------------------------------------------------------------------------


@traced
@logged
class Experiment(object):
  """
  
  """
  def __init__(self, name, **kwargs):
    self.name = name
    #attrs = ['type', 
    #         'source', 
    #         'shot_lines']
    #
    #for attr in attrs:
    #  setattr(self, attr, kw(attr, '?', kwargs))

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------


@traced
@logged
class Dataset(dict):
  """
  """
  def __init__(self, path, pattern, experiment, **kwargs):
    self.path = path
    self.ex = experiment
    self.fnames = get_files(path, pattern, **kwargs)
    self.names = [path_leave(i) for i in self.fnames]
    self.ids = []
    self.fN = 1 / (2 * self.dt) # Nyquist
    self.ttime = self.ns * self.dt
    for name in self.names:
      id = self.get_station_id(name)
      self.ids.append(id)
      self[id] = DataFileSgy(name, path) # FIXME: MAKE IT GENERIC
      # they are already defined in child classes:
      self[id].dt = self.dt # FIXME it should rather be done for any SgyFile...
      self[id].ns = self.ns
    self.ids = sorted(self.ids)
  
  # -----------------------------------------------------------------------------

  def ls(self, **kwargs):
    o, e = bash('pwd', path=self.path)
    print('Content of ' + o)
    o, e = bash('ls -lth ' + self.path)
    print(o, e)

  # -----------------------------------------------------------------------------

  def get_station_id(self, **kwargs):
    raise NotImplementedError('Implement in a child-class. It is experiment-specific!')

  # -----------------------------------------------------------------------------
  
  def qc(self, **kwargs):
    pass

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


@traced
@logged
class ShotLine(object):
  """

  """
  def __init__(self, id, experiment, **kwargs):
    self.id = id
    self.ex = experiment

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------


