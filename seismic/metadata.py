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

from fullwavepy.ioapi.generic import *
from fullwavepy.ioapi.segy import *

from fullwavepy.seismic.data import DataFileSgy
from fullwavepy.seismic.models import *
from fullwavepy.seismic.surfaces import *


@traced
@logged
class ProteusBathyTopo(BathyTopo):
  def __new__(cls, *args, **kwargs):
    #cls.dx = [50, 50, 50]  # [dx, dy, dz] m
    cls.z_sea = 0.0
    cls.x1 = -8.0e4
    cls.x2 = +8.0e4
    cls.y1 = -4.0e4
    cls.y2 = +4.0e4 
    cls.z1 = 0
    cls.z2 = 0     
    cls.extent = [[cls.x1, cls.x2], [cls.y1, cls.y2], [cls.z1, cls.z2]]
    kwargs['shape'] = (3201, 1601, 1)
    
    obj = super().__new__(cls, *args, **kwargs)
    cls.dx = obj.dx    
    # obj = obj.slice(slice_at='z')
    # obj.extent = [[cls.x1, cls.x2], [cls.y1, cls.y2]]
    # cls.extent = obj.extent
    # obj.dx = obj.dx[:-1]
    # cls.dx = obj.dx # otherwise is not passed in .copy()

    return obj


# -------------------------------------------------------------------------------


@traced
@logged
class ProteusStartVp(StartVp, AmphibiousModel):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class BenStartVp(ProteusStartVp):
  def __new__(cls, *args, **kwargs):
    #cls.dx = [50, 50, 50]  # [dx, dy, dz] m
  
    cls.x1 = -6.0e4
    cls.x2 = +6.4e4
    cls.y1 = -1.4e4
    cls.y2 = +2.9e4 
    # Z-axis points downwards   
    cls.z1 = -1.5e3 
    cls.z2 = +5.0e3      
    cls.extent = [[cls.x1, cls.x2], [cls.y1, cls.y2], [cls.z1, cls.z2]]
    kwargs['shape'] = (2481, 861, 131)
    
    obj = super().__new__(cls, *args, **kwargs)
    cls.dx = obj.dx # otherwise is not passed in .copy()

    obj.k_peak = 23
    obj.z_peak = -obj.k_peak * obj.dx[-1]
    obj.k_sea = 30
    obj.z_sea = 0   # m it's important because all SR coords in SEGY data are relative to this
    return obj
  

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
                 'surfaces': heavyphd + '/surfaces/',  
                 'data': heavyphd + 'DATA/Santorini_2015/'}
    
    self.svp = {'bh': BenStartVp(self.path['start_mods']+'Ben_whole_model_18-04-24.sgy', shape=(2481,861,131)),}
                # 'bh_clip': BenStartVp(self.path['start_mods']+'Ben_whole_model_18-04-24_sea-clipped.sgy', shape=(2481,861,101))}
    
    self.bathytopo = ProteusBathyTopo(self.path['surfaces']+'bathy_x_-8e4_8e4_y_-4e4_4e4_cell_50.vtr')

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

