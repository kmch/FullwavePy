"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import BinaryProjFile, TextProjFile
from fullwavepy.ioapi.fw3d import TtrFile
from fullwavepy.ioapi.segy import SgyFile



@traced
@logged
class IndexFile(BinaryProjFile):
  """
  extension idx, are binary files that act as an index to their associated
  SEG-Y file, linking each SEG-Y trace number with the source and receiver numbers that are
  used within fullwave. These index files play a role that is analogous that played previously
  by ttm files.
  
  """
  def __init__(self, suffix, proj, path, **kwargs):
    self.name = proj.name + '-' + suffix + '.idx'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class IndexFileTtr(IndexFile, TtrFile):
  def __init__(self, suffix, proj, path, **kwargs):
    """
    Notes
    -----
    Note the extension! Still inherits from ttr though.
    
    """  
    self.suffix = suffix
    self.name = proj.name + '-' + suffix + '.ttm'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class SynIndexFileTtr(IndexFile, TtrFile):
  def __init__(self, proj, path, **kwargs):
    suffix = 'Observed-Time'
    super().__init__(suffix, proj, path, **kwargs)
    

# -------------------------------------------------------------------------------


@traced
@logged
class ObsIndexFileTtr(IndexFile, TtrFile):
  def __init__(self, proj, path, **kwargs):
    suffix = 'Observed-Time'
    super().__init__(suffix, proj, path, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class SynIndexFileSgy(IndexFile, SgyFile):
  def __init__(self, proj, path, **kwargs):
    suffix = 'Synthetic'
    super().__init__(suffix, proj, path, **kwargs)
    

# -------------------------------------------------------------------------------


@traced
@logged
class ObsIndexFileSgy(IndexFile, SgyFile):
  def __init__(self, proj, path, **kwargs):
    suffix = 'Observed'
    super().__init__(suffix, proj, path, **kwargs)


# -------------------------------------------------------------------------------


