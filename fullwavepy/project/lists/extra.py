"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

This module contains handles for collections of similar files e.g. models
for different iterations. 

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.parse import kw, strip

from fullwavepy.project.lists.basic import TimestepFileList
from fullwavepy.project.lists.deriv import SchedFileList, SlaveFileList

from fullwavepy.project.files.gridded.wavefields import *

# THIS DOESN'T NEED BE SEPARATE FOR EACH DUMP!!!! #FIXME
# DERIVED TYPES SHOULD SUFFICE


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
      self.__log.warning(self.init_err)
      
  # -----------------------------------------------------------------------------
  
  ##@widgets('it', 'slice_at', 'node', 'cmap')
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
    self.__log.warning('fwd set to %s' % fwd)
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
      self.__log.warning(self.init_err)
    
  # -----------------------------------------------------------------------------
  
  @timer
  def load(self, **kwargs):
    """
    Load all the data for fast plotting
    
    """
    its = range(1, len(self.it))
    sids = sorted(self.it[1].keys())
    lids = sorted(self.proj.i.obs.read_header(overwrite=1)[self.proj.sgyhw['lid']].unique())
    
    blocks = self.proj.i.rnf.blocks
    
    #freqs =[i['freq'] for i in self.proj.i.rnf.blocks]
    
    
    print(its)
    print(sids)
    print(lids)
    #print(freqs)
  

# -------------------------------------------------------------------------------


# FIXME
@traced
@logged
class WavefieldPlotter(object):
  """

  """
  def _plot_one(self, **kwargs):
    """
    """
    it = kw('it', 1, kwargs)
    sid = kw('sid', sorted(self.it[it].keys())[0], kwargs)
    tstep = kw('tstep', sorted(self.it[it][sid].keys())[0], kwargs)
    self.__log.debug('it: {}, sid: {}, tstep: {}'.format(it, sid, tstep))
    # it's unlikely to have >= 1e3 iterations and >= 1e4 tsteps
    kwargs['title'] = 'it %s, sid %s, tstep %s' %\
      (str(it).rjust(3,'0'), str(sid).rjust(5,'0'), str(tstep).rjust(5,'0'))
    self.it[it][sid][tstep].plot(**kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class WavefieldFileList(SlaveFileList, TimestepFileList):
  """
  A handle for dumped wavefield-snapshots.
  
  """
  @timer
  def __init__(self, proj, FileClass, it_max, **kwargs):
    """
    
    Notes
    -----
    It assumes FileClass to have a certain sequence of args
    which is prone to bugs whenever FileClass changes.
    
    """
    super().__init__(proj, **kwargs)
    
    tsteps = self._read_tsteps(**kwargs)

    nfiles_max = 200
    o, e = bash('ls {}/*-?w-*.vtr | wc -l'.format(self.proj.out.path))
    nfiles = int(o)    
    if nfiles > nfiles_max: # IT WOULD TAKE AAAGEEES OTHERWISE
      self.init_err = 'Cannot init {} because fw/bw nfiles={} > nfiles_max={}'.format(self.__class__, 
                                                                                nfiles, nfiles_max)


    if self.init_err is None:
      self.it[0] = None
      # for it in range(1, self.nits_total + 1):
      for it in range(1, it_max + 1):
        self.it[it] = {}
        for sid in self.sids:
          sid = int(sid)
          self.it[it][sid] = {}
          for ts in tsteps:
            self.it[it][sid][ts] = FileClass(proj, ts, sid, it, **kwargs)
    else:
      self.__log.warning(self.init_err)
  
  # ----------------------------------------------------------------------------- 
  
  def _read_tsteps(self, **kwargs):
    """
    """
    proj = self.proj
    
    step = proj.env.var['SLAVES_WAVEFIELDSVTR']
    if step is not None:
      if step < 0:
        tsteps = np.arange(abs(step), proj.ns+1, abs(step))
        tsteps = [int(t) for t in tsteps]
        self.__log.info("proj.env.var['SLAVES_WAVEFIELDSVTR']=" + str(step) +
                        ' => wavefield dumped at timesteps: ' + str(tsteps))
      elif step > 0:
        tsteps = [step] 
        self.__log.info('SLAVES_WAVEFIELDSVTR=%s => a single snaphsot at %ss' % (step, step)) 
      
      else:
        raise ValueError('SLAVES_WAVEFIELDSVTR=%s' % step)
      
    else:
      self.__log.warning('SLAVES_WAVEFIELDSVTR not set - returning empty list...')
      tsteps = []
    
    self.__log.debug('tsteps: %s' % str(tsteps))
    
    self.tsteps = tsteps
    return self.tsteps

  # -----------------------------------------------------------------------------
  
  def _all_files(self, **kwargs):
    """
    FIXME Should share same abstraction with __init__. 
    """
    files = []
    for it in self.it[1: ]:
      for sid, sid_files in sorted(it.items()):
        for tstep, tstep_file in sorted(sid_files.items()):
          files.append(tstep_file)
    return files

  # -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------- 

@traced
@logged
class ForwardWavefieldFileList(WavefieldFileList):
  """
  A handle for dumped forward-wavefield snapshots.
  
  """
  def __init__(self, proj, *args, **kwargs):
    super().__init__(proj, ForwardWavefieldFile, *args, **kwargs)


# ------------------------------------------------------------------------------- 


@traced
@logged
class BackpropWavefieldFileList(WavefieldFileList):
  """
  A handle for dumped backpropagated-wavefield snapshots.
  
  """
  def __init__(self, proj, *args, **kwargs):
    super().__init__(proj, BackpropWavefieldFile, *args, **kwargs)


# -------------------------------------------------------------------------------

