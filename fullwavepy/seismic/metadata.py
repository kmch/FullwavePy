"""
This module defines a seismic field-experiment as a whole.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
# from abc import ABC, abstractmethod
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
from fullwavepy.seismic.surfs import *

class Box3d:
  """
  Box = cuboidal domain.
  """
  def __init__(self, x1, x2, y1, y2, z1, z2):
    """
    Parameters
    ----------
    OBSOLETE
    extent : list
        Required format:
        [[x1, x2], [y1, y2], [z1, z2]]
    """
    self.extent = [[x1, x2], [y1, y2], [z1, z2]]
    self.x1, self.x2, self.y1, self.y2, self.z1, self.z2 = [x1,x2,y1,y2,z1,z2]
  def plot(self, **kwargs):
    from fullwavepy.plot.misc import plot_square, plot_box
    # kwargs['label'] = self.proj.name
    # plot_square(self.proj.box[0], self.proj.box[1], 
    #             self.proj.box[2], self.proj.box[3], **kwargs)       
  def plotly(self, fig=None, **kwargs):
    import plotly.graph_objects as go
    kwargs = {'line': dict(color=kw('color', 'red', kwargs), width=2, dash=None),
              'showlegend': False, 'name': self.proj.name # WILL APPEAR ON HOVER
             }
    
    X = np.arange(self.x1, self.x2+1)
    Y = np.arange(self.y1, self.y2+1)
    
    if fig is None:
      fig = go.Figure()    
    
    fig.add_trace(go.Scatter(x=X, y=[self.y1]*len(X), **kwargs))
    fig.add_trace(go.Scatter(x=X, y=[self.y2]*len(X), **kwargs))
    fig.add_trace(go.Scatter(x=[self.x1]*len(Y), y=Y, **kwargs))
    fig.add_trace(go.Scatter(x=[self.x2]*len(Y), y=Y, **kwargs))    
    
    return fig
  # -----------------------------------------------------------------------------
  def _box2dims(self, box, dx):
    x1, x2, y1, y2, z1, z2 = box
    assert x2 >= x1
    assert y2 >= y1
    assert z2 >= z1
    nx1 = int((x2 - x1) / dx) + 1 
    nx2 = int((y2 - y1) / dx) + 1  
    nx3 = int((z2 - z1) / dx) + 1     
    return nx1, nx2, nx3



@traced
@logged
class Metadata(pd.DataFrame):
  pass
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
@traced
@logged
class ShotLine(object):
  """

  """
  def __init__(self, id, experiment, **kwargs):
    self.id = id
    self.ex = experiment

  # -----------------------------------------------------------------------------  
