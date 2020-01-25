"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for peremoveission writing to k.chrapkiewicz17@imperial.ac.uk.

Notes
-----
Python devs usually don't write nested classes, neither do I.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.types.basic import ProjSyn, ProjInv


# -------------------------------------------------------------------------------


@traced
@logged
class ProjSynVsObs(ProjSyn): #FIXME?
  """
  Generate synthetics and compare 
  them to the observed data.

  Notes
  -----
  The main difference between this and ProjSyn 
  is handling OutSeis.sgy which is generated 
  ONLY for sgfullwavepy.ioapi, even though documentation
  claims to generate OutSeis.sgy (yes, .sgy)
  for 'fw3d' as well NOTE
  
  """
  
  # -----------------------------------------------------------------------------   
  
  def _init_input(self, **kwargs): 
    """
    
    """
    from fullwavepy.project.files.datalike.sgy import ObsDataFileSgy
    from fullwavepy.project.files.datalike.ttr import ObsDataFileTtr
    
    if self.io == 'sgy':
      ObsDataClass = ObsDataFileSgy
      suffix = 'OutSeis'
    elif self.io == 'fw3d':
      ObsDataClass = ObsDataFileTtr
      suffix = 'TO_CHECK_OUTSEIS_IN_FW3D'
      
    self.inp.outseis = ObsDataClass(suffix, self, self.inp.path, **kwargs)
    #self.inp.outseis_raw = ObsDataFile('OutSeis', self, self.inp.path, **kwargs)
    #self.inp.outseis_filt = ObsDataFile('OutSeisFilt', self, self.inp.path, **kwargs)
    #self.inp.outseis = ObsDataFile('OutSeisFiltMute', self, self.inp.path, **kwargs)
    
    #try:
    #  self.inp.outseis.files(timer=True)
    #except OSError as err_message: 
    #  self.__log.warning(str(err_message))
      
    #self.inp.outseis_filt.files(timer=True)
    #self.inp.outseis_raw.files(timer=True)
    
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------

@traced
@logged
class ProjInvSyn(ProjInv):
  """
  Inversion of synthetic data.
  
  """
  
  # -----------------------------------------------------------------------------  

  def __init__(self, name, syn_proj_name, **kwargs):
    """
    Notes
    -----
    FIXME this should be read from syn runfile.
    We don't pass kwargs (except for IO!) to ProjSyn 
    as they may differ (e.g. testing acoustic inv 
    of elastic data) and thus are dummy anyway.
    
    """
    super().__init__(name, **kwargs)
    del_kw('path', kwargs)
    
    self.__log.warn('Make sure ProjInvSyn kwargs match syn_proj ones!')
    
    if 'syn_proj_path' in kwargs:
      syn_proj_path = kwargs['syn_proj_path']
    else:
      self.__log.warn('syn_proj_path not provided. Assuming: ' + syn_proj_name)
      syn_proj_path = syn_proj_name
      
    self.__log.info('\n\nInitializing synthetic project to invert the data from...\n')
    self.syn = ProjSyn(syn_proj_name, path=syn_proj_path)

  # -----------------------------------------------------------------------------    

  def _prepare_input(self, **kwargs):
    """
    
    """
    #super()._prepare_input(**kwargs)
    
    #files = [self.inp.obs, self.inp.obs_idx,  
    
    if exists(self.inp.obs.fname):
      self.__log.warn(self.inp.obs.fname + ' already exists.')
    else:
      self.inp.obs.prepare(dupl=self.syn.out.syn.fname, **kwargs)

    if exists(self.inp.obs_idx.fname):
      self.__log.warn(self.inp.obs_idx.fname + ' already exists.')
    else:
      self.inp.obs_idx.prepare(dupl=self.syn.out.syn_idx.fname, **kwargs)

    #for mod in [self.inp.startvp]:
      #if exists(mod.fname):
      #  self.__log.warn(mod.fname + ' already exists.')
      #else:
      #  mod.prepare(dupl=self.syn.out.syn_idx.fname, **kwargs)
   
    if exists(self.inp.startvp.fname):
      self.__log.warn(self.inp.startvp.fname + ' already exists.')
    else:
      try:
        self.inp.startvp.prepare(dupl=self.syn.inp.truevp.bckgrnd.fname, **kwargs)   
      except AttributeError:
        self.__log.warn('syn.inp.truevp.bckg not implemented yet. Copying truevp')
        self.inp.startvp.prepare(dupl=self.syn.inp.truevp.fname, **kwargs)
    self.inp.rawsign.prepare(dupl=self.syn.inp.rawsign.fname, **kwargs)
    #self.inp.runfile.prepare(dupl=self.syn.inp.runfile.fname, **kwargs)
    self.__log.warn('You need to prepare the runfile manually!')
    
  # -----------------------------------------------------------------------------

  def plot_output(self, *args, **kwargs):
    super().plot_output(*args, **kwargs)

  # -----------------------------------------------------------------------------
    

# -------------------------------------------------------------------------------

