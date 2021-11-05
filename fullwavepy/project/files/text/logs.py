"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import TextProjFile
from fullwavepy.project.files.text.misc import JobFile


@traced
@logged
class LogFile(JobFile, TextProjFile):
  """
  Log files produced either by executables 
  or a cluster's batch system.

  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, path, suffix, run_id, **kwargs):
    """
    self.name = lambda run_id : proj.name + '-Out{run_id}.log'.format(run_id=run_id)
    
    """
    exten = 'log'
    super().__init__(proj, path, suffix, exten, run_id, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class OutLogFile(LogFile):
  """
  Standard output of
  Fullwave3D.exe.  
  
  """

  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, path, run_id, **kwargs):
    suffix = 'Out'
    super().__init__(proj, path, suffix, run_id, **kwargs)
  
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    Parse to extract resources info.
    
    """    
    c = super().read(**kwargs)
    return c

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ErrLogFile(LogFile):
  """
  Standard error of
  Fullwave3D.exe.  
  
  """

  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, path, run_id, **kwargs):
    suffix = 'Err'
    super().__init__(proj, path, suffix, run_id, **kwargs)
  
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class JobLogFile(LogFile):
  """
  Manage job logs of previous runs.
  
  Notes
  -----
  Job logs are created AFTER all commands from 
  the PBS script, that's why this script can't 
  manage it from the script itself.
  
  """
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class JobOutLogFile(JobLogFile):
  """
  Standard output of PBS job commands 
  defined in the PBS script plus extras
  from PBS system, if anything.
  
  Notes
  -----
  Example output:
  ============================================

        Job resource usage summary 

                  Memory (GB)    NCPUs
  Requested  :         2            40
  Used       :         0 (peak)   3.29 (ave)

  ============================================
  
  The more parallel code, the more peak values 
  approach the requested ones. But typically 
  they are much smaller.
  
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, path, run_id, **kwargs):
    suffix = 'JobOut'
    super().__init__(proj, path, suffix, run_id, **kwargs)
    
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    Parse to extract resources info.
    
    """
    c = super().read(**kwargs)

    resources = {'mem_req': [],
                 'mem_use': [],
                 'cpu_req': [],
                 'cpu_use': [],
                }
    
    found = False
    for i, line in enumerate(c):
      if line[0] == 'Requested':
        found = True
        resources['mem_req'].append(float(line[2]))
        resources['cpu_req'].append(float(line[3]))
        resources['mem_use'].append(float(c[i+1][2]))
        resources['cpu_use'].append(float(c[i+1][4]))
    
    if not found:
      raise IOError('Resources info not found in ' + self.fname)
        
    return resources  

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class JobErrLogFile(JobLogFile):
  """
  Standard error of PBS job commands 
  defined in the PBS script plus extras
  from PBS system if any.

  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, path, run_id, **kwargs):
    suffix = 'JobErr'
    super().__init__(proj, path, suffix, run_id, **kwargs) 
    
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    Parse to extract resources info.
    
    """
    c = super().read(**kwargs)
    return c

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------

