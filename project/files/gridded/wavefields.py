"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, exten, strip, path_leave
from fullwavepy.ioapi.fw3d import VtrFile
from fullwavepy.project.files.gridded.generic import ExtenGridProjFile


@traced
@logged
class WavefieldFile(ExtenGridProjFile, VtrFile):
  """
  Shot-specific wavefield.
  
  """
  @timer
  def __init__(self, proj, file_id, ts, sid, it, **kwargs):
    """
    
    Notes
    -----
    E.g. p01-fw-004000-csref04160-iter00001-taskid00001.vtr
    
    Work-around to parse different taskids. 
    
    """
    #from fullwavepy.generic.system import get_files
    self.proj = proj
    self.path = proj.out.path
    self.ts = ts 
    self.sid = sid
    self.it = it
    
    if float(ts).is_integer():
      ts = str(ts).rjust(6,'0') # YES, 6 DIGITS
    else:
      ts = str(ts)
    
    tid = 1
    self.__log.debug('Assuming taskid = %s' % tid)
    
    self.name = str(proj.name + '-' + 
                    file_id + '-' + 
                    ts + 
                    '-csref' + str(sid).rjust(5,'0') + 
                    '-iter' + str(it).rjust(5,'0')  + 
                    '-taskid' + str(tid).rjust(5,'0') +
                    '.vtr')
    self.path = self.proj.out.path
    
    #fnames = get_files(path, pattern, **kwargs)
    #if len(fnames) > 1:
    # raise ValueError('Pattern {} matched by more than one file: {}'.format(pattern, fnames))
    #elif len(fnames) == 0:
    # self.__log.warn('Pattern {} matched by none of the files. Returning...'.format(pattern))
    #else:
    # suffix = strip(path_leave(fnames[0]))[len(proj.name + '-'): ]
    self.fname = self.path + self.name
    super().__init__(proj, self.path, **kwargs)
    
  # -----------------------------------------------------------------------------
  
  def plot(self, **kwargs):
    kwargs['cmap'] = kw('cmap', 'seismic', kwargs) # BrBG (PALER)
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    super().plot(**kwargs)

  # -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------- 


@traced
@logged
class ForwardWavefieldFile(WavefieldFile):
  """
  Forward wavefield.
  
  """
  def __init__(self, proj, ts, sid, it, **kwargs):
    file_id = 'fw'
    super().__init__(proj, file_id, ts, sid, it, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class BackproWavefieldFile(WavefieldFile):
  """
  Backropagated wavefield.
  
  """
  def __init__(self, proj, ts, sid, it, **kwargs):
    file_id = 'bw'
    super().__init__(proj, file_id, ts, sid, it, **kwargs)


# -------------------------------------------------------------------------------


