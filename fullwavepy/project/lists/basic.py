"""
This module contains handles for collections of similar files e.g. models
for different iterations. 
It allows, among others, to plot all items of the collection at once.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.system import bash, exists
from fullwavepy.plot.generic import FileListPlotter


@traced
@logged
class ProjFileList(FileListPlotter):
  """
  A handle for a collection of similar files.
  
  Notes
  -----
  When extending __init__ in a child class, BEGIN it with:
  super().__init__(proj, **kwargs) 
  
  """
  def __init__(self, proj, **kwargs):
    """
    """
    self.proj = proj
  
  # -----------------------------------------------------------------------------  
  
  def _all_files(self, **kwargs):
    """
    Reminder what child classes must implement.
    """
    raise NotImplementedError()

  # -----------------------------------------------------------------------------
@traced
@logged
class IterFileList(ProjFileList):
  """ 
  A handle for a collection of files dumped at certain iteration(s).
  
  """
  def __init__(self, proj, **kwargs):
    """
    """
    super().__init__(proj, **kwargs)
    self._read_iters(**kwargs)
     
    self.init_err = None
    
    if hasattr(self, 'nits_total'):
      self.it = list(np.zeros(self.nits_total+1))
    else:
      #self.init_err = 'Cannot init {} because self.nits_total not set'.format(self.__class__)
      self.init_err = 'Cannot init {} because self.nits_total not set'.format('')
  
  # -----------------------------------------------------------------------------
      
  def _read_iters(self, **kwargs):
    """
    Notes
    -----
    We need it for all problems 
    (e.g. for synthetic problems to handle wavefield snapshots)
    
    This should actually be stored from the beginning in proj.inp.runfile 
    or somewhere.
    
    """
    if self.proj.problem == 'synthetic':
      self.__log.debug('synthetic problem => self.nits_total = 1 regardless' +
                       ' of the iteration blocks in the runfile')
      self.nits_total = 1
      return    
    
    try:
      self.nits_total = self.proj.inp.runfile.read_nits(**kwargs)
    except FileNotFoundError:
      self.__log.warning(self.proj.name + ' has no runfile yet, unable to set' + 
                      ' self.nits_total')
     
  # -----------------------------------------------------------------------------
@traced
@logged
class ShotFileList(ProjFileList):
  """ 
  A handle for a collection of files dumped for certain shot(s).
  
  """
  def _read_sids(self, **kwargs):
    """
    SLAVES_DUMPCSREFS=[list/ranges]
    Restrict dumps to those CSRefs in the list, which is a
    comma-separated set of values/ranges, e.g. "5,15-20,35-45,60"
    
    """
    proj = self.proj
    
    try:
      sids_all = [s.ID for s in self.proj.inp.s.read(dx=1).li] # dx IS DUMMY HERE
    except AttributeError as err:
      self.__log.warning(str(err))
      return []

    except FileNotFoundError as err:
      self.__log.warning('Returning [] due to FileNotFoundError: ' + str(err)) #FIXME: WHY DOESN'T err CONTAIN FileNotFoundError
      return []
    
    # .get IS A dict's METHOD
    key = 'SLAVES_DUMPCSREFS' # WILL BE USED SEVERAL TIMES
    csrefs = self.proj.env.var.get(key, None)
    
    if csrefs is None:
      self.__log.info(key + ' not defined. Assuming dumps for all sids')
      sids = sids_all
    else:
      sids = self._parse_sid_ranges(csrefs, **kwargs)
    
    # we don't sort any more since now they're read in the order they appear in geom files into a list
    # sids = sorted(sids) 
    return sids

  # -----------------------------------------------------------------------------      
    
  def _parse_sid_ranges(self, csrefs, **kwargs):
    """
    csrefs : str
      E.g. "5,15-20,35-45,60"
    
    """
    sids = []
    for r in csrefs.split(','):
      i = [int(i) for i in r.split('-')]
      try:
        i1, i2 = i
        i = list(np.arange(i1, i2+1))
      except ValueError:
        pass
      sids += i
    
    return sids

  # -----------------------------------------------------------------------------
@traced
@logged
class TimestepFileList(ProjFileList):
  """ 
  A handle for a collection of files dumped at certain timestep(s).
  
  """
  def _read_tsteps(self, **kwargs):
    raise NotImplementedError('Should be overwritten in a child-class')
@traced
@logged
class JobFileList(ProjFileList):
  """
  A handle for a collection of files associated with every job run either
  locally or on a cluster.
  
  """
  def __init__(self, proj, path, FileClass, **kwargs):
    """
    New approach to file-containers (static list).
    
    Notes
    -----
    Runs are counted from 0 now!
    
    """
    self.proj = proj  
    self.path = path
    self.fileclass = FileClass

    # CONTAINER OF ACTUAL FILES CORRESPONDING TO PARTICULAR RUNS
    self.max_no_runs = 100 
    self.__log.debug('Max. no. of runs set to {}'.format(self.max_no_runs))
    self.no = list(np.zeros(self.max_no_runs))
    
    for run_id in range(self.max_no_runs):
      self.no[run_id] = FileClass(proj, path, run_id, **kwargs)

    # ALSO, FOR BACKWARD COMPATIBILITY WITH OLDER PROJECTS, CREATE A HANDLE 
    # FOR A LOG FILE WITHOUT run_id
    self.single = FileClass(proj, path, run_id=None, **kwargs)

  
  # -----------------------------------------------------------------------------  
  
  def read(self, run_ids, **kwargs):
    """
    Concatenate all files of the list
    for which run_id is in run_ids.
    
    """
    from fullwavepy.ioapi.generic import read_txt_raw
    
    if None in run_ids and len(run_ids) > 1:
      raise TypeError('None in run_ids and len(run_ids) > 1.')

    c = ""
    for run_id in run_ids:
      if run_id is None:
        self.__log.debug('run_id is None. Assuming a single log file.')
        lines = read_txt_raw(self.single.fname, **kwargs)
      else:
        lines = read_txt_raw(self.no[run_id].fname, **kwargs)
      for line in lines:
        c += line
    return c

  # ----------------------------------------------------------------------------- 
#@traced
#@logged
#class DumpFileList(ProjFileList):
  #"""
  
  #"""
  
  ## ----------------------------------------------------------------------------- 
  
  #def _check_if_enabled(self, env_var, **kwargs):
    #"""
    #"""
    #if self.env.var[env_var] == 'yes':
      #enabled = True
    #else:
      #enabled = False
      
    #return enabled
  
  ## -----------------------------------------------------------------------------
