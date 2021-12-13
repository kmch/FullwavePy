from abc import ABC, abstractmethod
from nsh.utils import extract_file_name, extract_path, get_file_names
from fullwavepy.ioapi.fw3d import VtrFile

class WaviefieldIO(ABC):
  def __init__(self):
    self._set_snapshot_fileclass()
    self._set_wavefield_pattern()
  @abstractmethod
  def _set_snapshot_fileclass(self):
    pass
  @abstractmethod
  def _set_wavefield_pattern(self):
    pass
class Fullwave3d(WaviefieldIO):
  def read(self, fname):
    pass
  def _extract_srcid(self, fname):
    srcid = int(fname.split('-csref')[1].split('-')[0].lstrip("0"))
    return srcid
  def _extract_tstep(self, fname):
    tstep = int(fname.split('-fw-')[1].split('-csref')[0].lstrip("0"))
    return tstep
  def _set_snapshot_fileclass(self):
    self.SnapshotFileCls = VtrFile
  def _set_wavefield_pattern(self):
    self.wavefield_pattern = '*-fw-*.vtr'
class Wavefield:
  def __init__(self, io='fullwave3d'):
    self.id = {}
    self._init_io(io)
  def get_files(self, proj_name, path):
    """
    Get all files present on the disk.
    Get = 'parse their names and init objects'.
    """
    fnames = self._get_file_names(proj_name, path)
    for fname in fnames:
      name = extract_file_name(fname)
      path = extract_path(fname)      
      srcid = self.io._extract_srcid(fname)
      tstep = self.io._extract_tstep(fname)      
      if srcid not in self.id:
        self.id[srcid] = {}
      self.id[srcid][tstep] = self.io.SnapshotFileCls(name, path)
  def _get_file_names(self, proj_name, path):
    pattern  = proj_name + self.io.wavefield_pattern
    self.fnames = get_file_names(path, pattern)
    return self.fnames
  def _get_src_ids(self, fnames):
    sids = set() # make them unique
    for fname in self.fnames:
      sid = self.io._extract_srcid(fname)
      sids.add(sid)
    return sorted(list(sids))
  def _init_io(self, io):
    if io == 'fullwave3d':
      self.io = Fullwave3d()
    else:
      raise NotImplementedError()
