"""
This module defines a seismic field-experiment as a whole.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, path_leave
from fullwavepy.generic.system import bash, get_files

from fullwavepy.seismic.data import DataFileSgy


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


class ProteusDataset(Dataset):
  """
  """
  def __init__(self, path, pattern, experiment, **kwargs):
    self.cruise_id = 'MGL1521'
    self.sep = '_'
    self.pools = ['L', 'W', 'S']
    self.channels = ['1', '2', '3', '4']
    self.ext = 'sgy'
    self.len_station_id = 3
    super().__init__(path, pattern, experiment, **kwargs)

  # -----------------------------------------------------------------------------  
  
  def get_station_id(self, name):
    return name[len(self.cruise_id+self.sep+self.pools[0]):-len(self.sep+self.channels[0]+'.'+self.ext)]

  # -----------------------------------------------------------------------------

  def get_pool_id(self, name):  
    return name[len(self.cruise_id+self.sep):-(self.len_station_id+len(self.sep+self.channels[0]+'.'+self.ext))]

  # -----------------------------------------------------------------------------

  def get_channel_id(self, name):
    return name[len(self.cruise_id+self.sep+self.pools[0])+self.len_station_id+len(self.sep):-len('.'+self.ext)]
  
  # -----------------------------------------------------------------------------

  

# -------------------------------------------------------------------------------


class ProteusDatasetOBS(ProteusDataset):
    def __init__(self, *args, **kwargs):
      self.dt = 0.005 # s
      self.ns = 12000 # samples
      super().__init__(*args, **kwargs)


# -------------------------------------------------------------------------------


class ProteusDatasetLand(ProteusDataset):
    def __init__(self, *args, **kwargs):
      self.dt = 0.001 # s
      self.ns = 7000  # samples
      super().__init__(*args, **kwargs)


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
class ProteusExperiment(Experiment):
  """
  """
  def __init__(self, **kwargs):
    """
    """
    from fullwavepy.seismic.models import StartVp
    
    self.name = 'PROTEUS'
    self.cruise_id = 'MGL1521'
    self.type = 'marine-land'
    self.source = 'airgun'
    self.shot_lines = [ShotLine(i, self) for i in np.arange(1,61)]

    heavyphd = '/home/kmc3817/heavy_PhD/'
    self.path = {'start_mods': heavyphd + '/start_mods/',
                 'data': heavyphd + 'DATA/Santorini_2015/'}
    
    self.svp = {'bh_full': StartVp(self.path['start_mods']+'Ben_whole_model_18-04-24.sgy', shape=(2481,861,131)),
                'bh_clip': StartVp(self.path['start_mods']+'Ben_whole_model_18-04-24_sea-clipped.sgy', shape=(2481,861,101))}
    
    self.dataset = {'obshy': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*4.sgy', self),
                    'obsvx': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*3.sgy', self),
                    'obsvy': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*2.sgy', self),
                    'obsvz': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*1.sgy', self),
                    'lanvx': ProteusDatasetLand(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*3.sgy', self),
                    'lanvy': ProteusDatasetLand(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*2.sgy', self),
                    'lanvz': ProteusDatasetLand(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*1.sgy', self)}

    super().__init__(self.name, **kwargs)


# metadataobs = path_dataobs + 'metadata.csv'
# metadatalan = path_datalan + 'metadata.csv'
# metadata = path + 'metadata/prot
 
 
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


@traced
@logged
class Metadata(object): #FIXME?
  def __init__(self, dataset, **kwargs):
    pass

  # -----------------------------------------------------------------------------  
 

# -------------------------------------------------------------------------------

