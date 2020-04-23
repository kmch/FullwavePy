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
from fullwavepy.generic.system import get_files


@traced
@logged
class Dataset(dict):
  def __init__(self, path, pattern, experiment, **kwargs):
    self.path = path
    self.ex = experiment
    self.fnames = get_files(path, pattern, **kwargs)
    self.names = [path_leave(i) for i in self.fnames]
    for i in self.fnames:
      id = self.get_station_id(path_leave(i))
      self[id] = i

  # -----------------------------------------------------------------------------

  def get_station_id(self, **kwargs):
    raise NotImplementedError('It is experiment-specific!')

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
class Proteus(Experiment):
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
    
    self.dataset = {'obshy': ProteusDataset(self.path['data']+'seismic/OBS/segy_local_coords/', '*4.sgy', self),
                    'obsvx': ProteusDataset(self.path['data']+'seismic/OBS/segy_local_coords/', '*3.sgy', self),
                    'obsvy': ProteusDataset(self.path['data']+'seismic/OBS/segy_local_coords/', '*2.sgy', self),
                    'obsvz': ProteusDataset(self.path['data']+'seismic/OBS/segy_local_coords/', '*1.sgy', self),
                    'lanvx': ProteusDataset(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*3.sgy', self),
                    'lanvy': ProteusDataset(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*2.sgy', self),
                    'lanvz': ProteusDataset(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*1.sgy', self)}

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

