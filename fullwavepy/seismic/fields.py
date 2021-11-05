from abc import ABC, abstractmethod
from nsh.utils import get_file_names

class WaviefieldIO(ABC):
  def __init__(self):
    self._set_wavefield_pattern()
  @abstractmethod
  def _set_wavefield_pattern(self):
    pass
class Fullwave3d(WaviefieldIO):
  def _extract_srcid(self, fname):
    srcid = int(fname.split('-csref')[1].split('-')[0].lstrip("0"))
    return srcid
  def _extract_tstep(self, fname):
    tstep = int(fname.split('-fw-')[1].split('-csref')[0].lstrip("0"))
    return tstep
  def _set_wavefield_pattern(self):
    self.wavefield_pattern = '*-fw-*.vtr'
class Wavefield:
  def __init__(self):
    self._init_io()
  def get_files(self):
    """
    Get all files present on the disk.
    Get = 'parse their names and init objects'.
    """
    fnames = self._get_file_names()
    # sids = self._get_src_ids(fnames)
    # self.list = []
    # self.dict = {}
    # self.shot = {}
    # self._get_file_names()
    for fname in fnames:
      srcid = self.io._extract_srcid(fname)
      tstep = self.io._extract_tstep(fname)
      if srcid in self.shot:
        self.shot[srcid] = WavefieldSnapshot(tstep)
  def _init_io(self):
    self.io = Fullwave3d()
  def _get_file_names(self):
    path = '/home/kmc3817/rds_home/my_ephemeral/PROJECTS/ch_kol/f01syn/out/'
    self.fnames = get_file_names(path, self.io.wavefield_pattern)
    return self.fnames
  def _get_src_ids(self, fnames):
    sids = set()
    for fname in fnames:
      sid = self.io._extract_srcid(fname)
      sids.add(sid)
    return sorted(list(sids))
class WavefieldSingleShot:
  pass
class WavefieldSnapshot:
  def __init__(self, tstep):
    pass
class WavefieldSingleShotSnapshot:
  pass
