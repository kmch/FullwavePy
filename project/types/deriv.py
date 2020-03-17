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
class ProjInvSyn(ProjInv):
  """
  Inversion of synthetic data.
  
  Notes
  -----
  Input files of the synthetic project must be 
  available.
  
  """
  
  # -----------------------------------------------------------------------------  

  def __init__(self, name, syn_proj_name, **kwargs):
    """
    Same as ProjInv + book-keeping of an extra
    synthetic subproject.
    
    Notes
    -----
    We don't pass kwargs (except for IO!) to ProjSyn 
    as they may differ (e.g. testing acoustic inv 
    of elastic data) and thus are dummy anyway.
    
    """
    super().__init__(name, **kwargs)
    self.__log.warn('Make sure ProjInvSyn kwargs match syn_proj ones!')
    
    if 'syn_proj_path' in kwargs:
      syn_proj_path = kwargs['syn_proj_path']
    else:
      self.__log.warn('syn_proj_path not provided. Assuming: ' + syn_proj_name)
      syn_proj_path = syn_proj_name
      
    self.__log.info('\n\nInitializing synthetic sub-project to copy data from\n')
    self.psyn = ProjSyn(syn_proj_name, path=syn_proj_path)

  # -----------------------------------------------------------------------------    

  def prepare_input(self, file_z0, **kwargs):
    """
    Copy from the synthetic project.
    
    Notes
    -----
    synth-Synthetic.sgy -> inver-Observed.sgy 
    synth-Template.[idx/hed] -> inver-Observed.[idx/hed]
    
    Skeleton.key is essential in creating the Runfile.
    
    """
    self.__log.info('Duplicating relevant files from ' + self.psyn.name)
    
    for f_id in ['rawsign', 'sgn', 'sp', 's', 'r', 'skeleton']:
      finv = getattr(self.i, f_id)
      fsyn = getattr(self.psyn.i, f_id)
      finv.prep(dupl=fsyn.fname)
    
    for f_id in ['hed', 'idx']:
      # SIGNATURE
      finv = getattr(self.i.sgn, f_id) 
      fsyn = getattr(self.psyn.i.sgn, f_id) 
      finv.prep(dupl=fsyn.fname)
      # TEMPLATE
      finv = getattr(self.i.obs, f_id) 
      fsyn = getattr(self.psyn.i.tmpl, f_id) 
      finv.prep(dupl=fsyn.fname)
    
    self.i.obs.prep(dupl=self.psyn.o.syn.fname) #FOR 1ST ITER
    self.i.obs.raw.prep(dupl=self.psyn.o.syn.fname) #NOTE: raw!
    self.__log.warn(self.i.obs.fname + ' is meant to be a processed version of ' +
                    self.i.obs.raw.fname + '. Note that if the processing changes ' +
                    'no. of srcs/recs, SegyPrep has to be re-run!')
    
    try:
      self.i.svp.prep(dupl=self.psyn.i.tvp.bckgnd.fname, file_z0=file_z0)
    except AttributeError as err:
      self.__log.warn('Copying TrueVp because of: ' + str(err))
      self.i.svp.prep(dupl=self.psyn.i.tvp.fname, file_z0=file_z0)
    
    
    self.__log.warn('Runfile and PBS script need to be created manually!')
    
  # -----------------------------------------------------------------------------  
  
  def plot_input(self, **kwargs):
    super().plot_input(**kwargs)
    
  # -----------------------------------------------------------------------------

  def plot_output(self, *args, **kwargs):
    super().plot_output(*args, **kwargs)

  # -----------------------------------------------------------------------------
    

# -------------------------------------------------------------------------------

