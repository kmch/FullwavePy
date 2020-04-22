"""
This module defines objects that describe
a seismic field experiment as a whole.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw


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
    name = 'PROTEUS'
    self.type = 'marine-land'
    self.source = 'airgun'
    self.shot_lines = [ShotLine(i) for i in np.arange(1,61)]
    self.path = {'startmods': '/home/kmc3817/heavy_PhD/start_mods/'}
    
    path = self.path['startmods']
    self.svp = {'bh_full': StartVp(path+'Ben_whole_model_18-04-24.sgy', shape=(2481,861,131)),
                'bh_clip': StartVp(path+'Ben_whole_model_18-04-24_sea-clipped.sgy', shape=(2481,861,101))
    }

    super().__init__(name, **kwargs)

#     startvp_jm = path + 'start_mods/jm_inversecheck-StartVp.sgy'
# startvp_bh = path + 'start_mods/Ben_whole_model_18-04-24_sea-clipped.sgy'
# startvp_bm = path + 'start_mods/Brennah_whole_model_19-10-30_sea-clipped.sgy'

#  dataobs_hy_sgy = [SgyFile(path_leave(i), path=path_dataobs) for i in dataobs_hy]
# dataobs_vz_sgy = [SgyFile(path_leave(i), path=path_dataobs) for i in dataobs_vz]
# datalan_vz_sgy = [SgyFile(path_leave(i), path=path_datalan) for i in datalan_vz]
#  from fullwavepy.generic.system import get_files
# path = '/home/kmc3817/heavy_PhD/'
# path_dataobs = path + 'DATA/Santorini_2015/seismic/OBS/segy_local_coords/'
# path_datalan = path + 'DATA/Santorini_2015/seismic/land/Santorini/segy_local_coords/'
# dataobs_hy = get_files(path_dataobs, '*4.sgy')
# dataobs_vz = get_files(path_dataobs, '*1.sgy')
# datalan_vz = get_files(path_datalan, '*1.sgy')

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
  def __init__(self, id, **kwargs):
    self.id = id

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------


@traced
@logged
class Metadata(object): #FIXME?
  pass
 
  # -----------------------------------------------------------------------------  
 

# -------------------------------------------------------------------------------

