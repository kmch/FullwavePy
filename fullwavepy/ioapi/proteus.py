"""
This module defines all objects specific to the PROTEUS experiment.
It is useful to hard-wire some of the parameters as class attributes 
to exchange them between different jupyter notebooks. Keeping them
in a package seems like a better backup anyway.

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
from fullwavepy.seismic.surfs import *
from fullwavepy.seismic.metadata import * 
from fullwavepy.seismic.wavelets import *
from fullwavepy.project.types.deriv import *

@traced
@logged
class ProteusExperiment(Experiment):
  """
  """
  def __init__(self, **kwargs):
    """
    FIXME: I don't know why it creates an empty figure 
    when init. for the first time...    
    """
    from fullwavepy.seismic.models import StartVp

    self.sgyhw = {'sid': 'tracf',
                  'rid': 'fldr',
                  'lid': 'ep',
                 }

    self.name = 'PROTEUS'
    self.cruise_id = 'MGL1521'
    self.type = 'marine-land'
    self.source = 'airgun'
    self.shot_lines = [ShotLine(i, self) for i in np.arange(1,61)]
    # metadataobs = path_dataobs + 'metadata.csv'
    # metadatalan = path_datalan + 'metadata.csv'
    base_path = '/home/kmc3817/heavy_PhD/'
    self.path = {'data': base_path + 'DATA/Santorini_2015/',
                 'metadata': base_path + 'metadata/',
                 'start_mods': base_path + 'start_mods/',
                 'surfaces': base_path + 'surfaces/',  
                 'wavelets': base_path + 'wavelets/',
                 }
    
    self.md = MetadataFile('proteus_metadata.csv', self.path['metadata'])
    
    self.bathytopo = ProteusBathyTopo(self.path['surfaces'] + \
      'bathy_x_-8e4_8e4_y_-4e4_4e4_cell_50.vtr')

    self.dataset = {'obshy': ProteusDatasetOBS(self.path['data'] + \
      'seismic/OBS/segy_local_coords/', '*4.sgy', self),
                    'obsvx': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*3.sgy', self),
                    'obsvy': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*2.sgy', self),
                    'obsvz': ProteusDatasetOBS(self.path['data']+'seismic/OBS/segy_local_coords/', '*1.sgy', self),
                    'lanvx': ProteusDatasetLand(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*3.sgy', self),
                    'lanvy': ProteusDatasetLand(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*2.sgy', self),
                    'lanvz': ProteusDatasetLand(self.path['data']+'seismic/land/Santorini/segy_local_coords/', '*1.sgy', self)}

    self.wavelet = {'19-09-22': ProteusWavelet(self.path['wavelets'] + \
      'wavelet_19-09-22.sgy')}
    
    path = self.path['start_mods']
    # 'bh_clip': BenStartVp(self.path['start_mods']+'Ben_whole_model_18-04-24_sea-clipped.sgy', shape=(2481,861,101))}
    self.svp = {'bh': {'18-04-24': BenStartVp(path+'Ben_whole_model_18-04-24.sgy'),
                       '18-04-24_kameni1': StartVp(path+'bh_18-04-24_kameni1.vtr', k_fs=7, extent=[[-2500.,  2200.], [-2500.,  4300.], [ -400.,  3000.]], shape=(95, 137, 69)),
                       '19-04-25_v1': StartVp(path+'srModel_it5_velocity_19-04-25_-20000_20000_-5000_20000_-1000_5000.vtr', k_fs=20, extent=[[-2e4,2e4],[-5e3,2e4],[-1e3,5e3]], shape=(801, 501, 121)),
                      }
    }
                
    
    # Add associated free-surfaces to the models, if present
    self.svp['bh']['18-04-24_kameni1'].fs = FreeSurface(self.path['surfaces']+'fs_bh_18-04-24_kameni1.vtr', shape=(95, 137, 1))

    super().__init__(self.name, **kwargs)

   # -----------------------------------------------------------------------------  
# -------------------------------------------------------------------------------
# Datasets
# -------------------------------------------------------------------------------
@traced
@logged
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
@traced
@logged
class ProteusDatasetOBS(ProteusDataset):
    def __init__(self, *args, **kwargs):
      self.dt = 0.005 # s
      self.ns = 12000 # samples
      super().__init__(*args, **kwargs)
@traced
@logged
class ProteusDatasetLand(ProteusDataset):
    def __init__(self, *args, **kwargs):
      self.dt = 0.01  # s (value 10000 microseconds in SEGY header)
      self.ns = 7000  # samples

      super().__init__(*args, **kwargs)
# -------------------------------------------------------------------------------
# Wavelets
# -------------------------------------------------------------------------------
@traced
@logged
class ProteusWavelet(Wavelet):
    pass
# -------------------------------------------------------------------------------
# Arrays
# -------------------------------------------------------------------------------
@traced
@logged
class ProteusBathyTopo(BathyTopo):
  def __new__(cls, *args, **kwargs):
    #cls.dx = [50, 50, 50]  # [dx, dy, dz] m
    
    x1 = -8.0e4
    x2 = +8.0e4
    y1 = -4.0e4
    y2 = +4.0e4 
    z1 = 0
    z2 = 0     
    cls.extent = [[x1, x2], [y1, y2], [z1, z2]] # otherwise is not passed in .copy()
    kwargs['shape'] = (3201, 1601, 1)
    
    obj = super().__new__(cls, *args, **kwargs)
    cls.dx = obj.dx # otherwise is not passed in .copy()
    cls.z_sea = 0.0
    cls.k_sea = 29
    
    # obj = obj.slice(slice_at='z')
    # obj.extent = [[cls.x1, cls.x2], [cls.y1, cls.y2]]
    # cls.extent = obj.extent
    # obj.dx = obj.dx[:-1]
    # cls.dx = obj.dx # otherwise is not passed in .copy()

    return obj
@traced
@logged
class ProteusStartVp(StartVp):
  pass
@traced
@logged
class BenStartVp(ProteusStartVp):
  def __new__(cls, *args, **kwargs):
    #cls.dx = [50, 50, 50]  # [dx, dy, dz] m
  
    x1 = -6.0e4
    x2 = +6.4e4
    y1 = -1.4e4
    y2 = +2.9e4 
    # Z-axis points downwards   
    z1 = -1.5e3 
    z2 = +5.0e3      
    cls.extent = [[x1, x2], [y1, y2], [z1, z2]] # otherwise is not passed in .copy()
    kwargs['shape'] = (2481, 861, 131)
    
    obj = super().__new__(cls, *args, **kwargs)
    cls.dx = obj.dx # otherwise is not passed in .copy()

    # obj.k_peak = 23
    # obj.z_peak = -obj.k_peak * obj.dx[-1]
    cls.k_sea = 29
    cls.z_sea = 0   # m it's important because all SR coords in SEGY data are relative to this
    return obj
# -------------------------------------------------------------------------------
# Projects - NOT used at the moment
# -------------------------------------------------------------------------------
@traced
@logged
class ProteusProj(ProjExperiment):
  def __init__(self,  *args, **kwargs):
    # ex = ProteusExperiment() # disabled for speed, init once in ipynb and pass as kw
    # super().__init__(ex, *args, **kwargs) ##
    super().__init__(*args, **kwargs)
@traced
@logged
class ProteusProjSyn(ProteusProj, ProjSyn):
    pass
@traced
@logged
class ProteusProjInv(ProteusProj, ProjInv):
    pass
# -------------------------------------------------------------------------------
