"""
I/O interace for FWI codes.

Notes
-----
We have:
fwicode.io => preprocessor
generic_input --> preprocessor --> fwicode.io

One I/O may serve different FWI codes, and conversely,
one FWI code may be served by different pre-processors (less likely).
"""
# -------------------------------------------------------------------------------
# FWI i/o
# -------------------------------------------------------------------------------
class FwiIo(ABC):
  """
  Note that IO is something separate
  from an FWI code even if it is implemented
  only in that code at that moment. It can
  be transplanted to other codes in the 
  future.
  
  It should link project's objects with file-objects,
  or (better) read from those files.
  """    
  def __init__(self, prefix, shell='linux', path_inp='./inp/', path_out='./out/'):
    """
    prefix : str
      Typically, the project's name.
      
    """
    self.prefix = prefix
    self.shell = shell
    self.path_inp = path_inp
    self.path_out = path_out
  @abstractmethod
  def read_srcs_ids_coords_types(self):
    return ids, xyz, types
  @abstractmethod
  def read_recs(self):
    return ids, xyz, types    
class FwiIoFw3d(FwiIo):
  """
  One of the IOs implemented in the 
  Fullwave3D code.
  """
  def read_srcs(self):
    NIErr()
  def read_recs(self):
    NIErr()              
class FwiIoSegy(FwiIo):
  """
  One of the IOs implemented in the 
  Fullwave3D code.
  """    
  def read_recs(self):
    NIErr()       
  def read_srcs_ids_coords_types(self, **kwargs):
    ids, coords = self.read_srcs_ids_coords(**kwargs)
    types = self.read_srcs_types(**kwargs)
    return ids, coords, types
  def read_srcs_ids_coords(self, **kwargs):
    """
    List preserves the order from file. It is the 
    same as in obs file anyway.
    """
    from fwilight.fileio import FileGeo
    # We may endow IO with project-derived params for convenience...
    #         if hasattr(self, 'btop'):
    #         kwargs['btop'] = self.btop
    path = self.path_inp
    name = self.prefix + 'Sources.geo'
    fname = path + name
    #print('SEGY IO (Fullwave3D) => sources will be read from %s' % fname)
    f = FileGeo(name, path, self.shell)
    li = f.read(**kwargs)
    ids = [i[0] for i in li]
    coords = [i[1] for i in li]
    return ids, coords
  def read_srcs_types(self, **kwargs):
    path = self.path_inp
    return 'Nothing here yet.'
  def read_srcs_from_dobs(self):
    NIErr()
# -------------------------------------------------------------------------------
# FWI files
# -------------------------------------------------------------------------------
class Rnf(Runfile):
  pass
class Rnf728(Rnf):
  def create(self):
    self.dupl(self.proj.i.ske.fname)
  def prep(self, *args, **kwargs):
    self.create()
  def read(self):
    from fullwavepy.ioapi.generic import read_txt_raw
    content = read_txt_raw(self.fname)
    params = {}
    for line in content:
      line = line.lstrip().rstrip()
    if len(line) >= 3 and line[0][0] != '!':
      #print(line)
      line = line.split('!')[0] # do it first as there might be ':' in the trailing comment
      key, val = line.lower().split(':')
      key = key.rstrip().lstrip()
      val = val.split('!')[0].rstrip().lstrip()
      params[key] = val
    return sorted(params.items())
class Ske728(Skeleton,Rnf728):
  pass
