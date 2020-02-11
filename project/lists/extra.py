"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

This module contains handles for collections of similar files e.g. models
for different iterations. 

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.parse import kw
from fullwavepy.project.lists.basic import TimestepFileList
from fullwavepy.project.lists.deriv import SchedFileList, SlaveFileList

# THIS DOESN'T NEED BE SEPARATE FOR EACH DUMP!!!! #FIXME
# DERIVED TYPES SHOULD SUFFICE


# -------------------------------------------------------------------------------  


@traced
@logged
class CPFileList(SchedFileList):
  """
  Checkpoint file.
  
  Notes
  -----
  It doesn't seem productive to create
  child classes for all Grad, Vp, etc.
  That's why the file_id seems like 
  a sensible idea.
  
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, proj, FileClass, file_id, file_start, **kwargs):  
    """
    file_id : str 
      'Vp', 'RawGrad', etc.
      
    file_start : object / None
      file object at iteration 0; 
      can be None as for derivatives 
      (Grad, Prec etc.)
      
    Notes
    -----
    It assumes FileClass to have:
      def __init__(self, suffix, proj, path, **kwargs)    
    
    """
    super().__init__(proj, **kwargs)
    suffix = lambda file_id, it : 'CP' + str(it).rjust(5, '0') + '-' + file_id
    
    if self.init_err is None:
      self.it[0] = file_start
      for it in range(1, self.nits_total + 1):
        self.it[it] = FileClass(suffix(file_id, it), proj, proj.out.path, **kwargs)      
    else:
      self.__log.warn(self.init_err)
      
  # -----------------------------------------------------------------------------
  
  @widgets('it', 'slice_at', 'node', 'cmap')
  def plot(self, wdg=False, **kwargs):
    it = kw('it', 1, kwargs)
    self.it[it].plot(**kwargs)
    
  # -----------------------------------------------------------------------------

# -------------------------------------------------------------------------------  


#@traced
#@logged
#class DerivativeFileList(CheckpointFileList):
  #def __init__(self, proj, FileClass, file_id, file_start, **kwargs): 


# -------------------------------------------------------------------------------  


@traced
@logged
class DumpFileList(SlaveFileList):
  """
  A handle for dumped data-files (not wavefields!).
  
  """

  # ----------------------------------------------------------------------------- 
  
  def __init__(self, proj, FileClass, file_id, fwd=1, **kwargs):
    """
    file_id : str 
      Vp, RawGrad etc.
      
    Notes
    -----
    It assumes FileClass to have:
      def __init__(self, suffix, proj, path, **kwargs)
    
    fwd is not present in dump_adjoint # FIXME
    
    """
    super().__init__(proj, **kwargs)
    self.__log.warn('fwd set to %s' % fwd)
    self.__log.debug('file_id:' + file_id)
    
    # THIS COULD BE PUT SOMWHERE ELSE IF NEEDED
    suffix = lambda file_id, it, sid, fwd : str(file_id +
                                                '-csref' + str(sid).rjust(5,'0') + 
                                                '-iter' + str(it).rjust(5,'0')  + 
                                                'fwd' + str(fwd))
    
    if self.init_err is None:
      self.it[0] = None
      for it in range(1, self.nits_total + 1):
        self.it[it] = {}
        for sid in self.sids:
          sid = int(sid)
          self.it[it][sid] = FileClass(suffix(file_id, it, sid, fwd), 
                                       proj, proj.out.path, it=it, sid=sid, **kwargs)
    else:
      self.__log.warn(self.init_err)
    
  # -----------------------------------------------------------------------------
  
  @timer
  def load(self, **kwargs):
    """
    Load all the data for fast plotting
    
    """
    its = range(1, len(self.it))
    sids = sorted(self.it[1].keys())
    lids = sorted(self.proj.i.obs.read_header(overwrite=1)[self.proj.sgy.hw['lid']].unique())
    
    blocks = self.proj.i.rnf.blocks
    
    #freqs =[i['freq'] for i in self.proj.i.rnf.blocks]
    
    
    print(its)
    print(sids)
    print(lids)
    #print(freqs)
  

# -------------------------------------------------------------------------------


@traced
@logged
class WavefieldFileList(SlaveFileList, TimestepFileList):
  """
  A handle for dumped wavefield-snapshots.
  
  """
  
  # ----------------------------------------------------------------------------- 
  
  @timer
  def __init__(self, proj, FileClass, **kwargs):
    """
    
    Notes
    -----
    It assumes FileClass to have a certain sequence of args
    which is prone to bugs whenever FileClass changes.
    
    """
    super().__init__(proj, **kwargs)
    tsteps = self._read_tsteps(**kwargs)

    nfiles_max = 200
    o, e = bash('ls {} | wc -l'.format(self.proj.out.path))
    nfiles = int(o)    
    if nfiles > nfiles_max: # IT WOULD TAKE AAAGEEES OTHERWISE
      self.init_err = 'Cannot init {} because nfiles={} > nfiles_max={}'.format(self.__class__, 
                                                                                nfiles, nfiles_max)


    if self.init_err is None:
      self.it[0] = None
      for it in range(1, self.nits_total + 1):
        self.it[it] = {}
        for sid in self.sids:
          sid = int(sid)
          self.it[it][sid] = {}
          for ts in tsteps:
            ts = int(ts)
            tid = '00001' #'?????'
            self.it[it][sid][ts] = FileClass(proj, ts, sid, it, tid, **kwargs)
    else:
      self.__log.warn(self.init_err)
  
  # ----------------------------------------------------------------------------- 
  
  def _read_tsteps(self, **kwargs):
    """
    """
    proj = self.proj
    
    step = proj.env.var['SLAVES_WAVEFIELDSVTR']
    if step is not None:
      step = int(step)
      if step < 0:
        tsteps = np.arange(abs(step), proj.ns+1, abs(step))
        self.__log.info("proj.env.var['SLAVES_WAVEFIELDSVTR']=" + str(step) +
                        ' => wavefield dumped at timesteps: ' + str(tsteps))
      else:
        raise NotImplementedError('SLAVES_WAVEFIELDSVTR > 0')
        
    else:
      self.__log.warn('SLAVES_WAVEFIELDSVTR not set - returning empty list...')
      tsteps = []
    
    self.__log.debug('tsteps: %s' % str(tsteps))
    
    return tsteps

  # -----------------------------------------------------------------------------
  
  def plot(self, **kwargs):
    it = kw('it', 1, kwargs)
    sid = kw('sid', sorted(self.it[it].keys())[0], kwargs)
    tstep = kw('tstep', sorted(self.it[it][sid].keys())[0], kwargs)
    self.__log.debug('it: {}, sid: {}, tstep: {}'.format(it, sid, tstep))
    self.it[it][sid][tstep].plot(**kwargs)

  # -----------------------------------------------------------------------------
  
  def plot_all(self, **kwargs):
    """
    """
    self.__log.warn('assuming it[1: ] even for ProjSyn')
    for it in self.it[1: ]:
      for sid in it.values():
        for f in sid.values():
          try:
            f.plot(**kwargs)
            plt.figure()
          except FileNotFoundError as err:
            self.__log.warn(err)
  
  # -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------- 

