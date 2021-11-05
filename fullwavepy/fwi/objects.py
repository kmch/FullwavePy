"""
An interface between fwi projects (`fullwavepy.fwi.projects`) 
and seismic data-structures (`fullwavepy.seismic`).
It endowes seismic data-structures with project-specific parameters.
"""
class ProjMod(Mod):
  def __init__(self, proj):
    self.proj = proj
    self.dims = proj.dims
  def create(self, *args, **kwargs):
    return super().create(self.dims, *args, **kwargs)
class ProjFs(Fs):
  def __init__(self, proj):
    self.proj = proj
    self.dims = proj.dims
  def create(self, *args, **kwargs):
    return super().create(self.dims, *args, **kwargs)
class ProjWvl(Wvl):
  def __init__(self, proj):
    self.proj = proj
    self.dt = proj.dt
    self.ns = proj.ns
  def create(self, *args, **kwargs):
    return super().create(self.dt, self.ns, *args, **kwargs)  
class ProjSrcs(Srcs):
  def __init__(self, proj):
    self.proj = proj
    #     def read(self, *args, **kwargs):
    #         self.proj.fwicode.io.read_srcs(*args, **kwargs)
class ProjGrad:
  def __init__(self, proj):
    self.p = proj
class ProjGradGlob(ProjGrad):
  def plot(self, it, **kwargs):
    self.read(it)
    ax = self.arr.plot(**dict(kwargs, center_cmap=1))
    aspeqt(ax)
  def read(self, it):
    from fullwavepy.ndat.arrays import Arr3d
    from fullwavepy.ioapi.fw3d import read_vtr
    self._set_name(it)
    self.arr = Arr3d(read_vtr(self.fname))
  def _set_name(self, it):
    self.path = self.p.o.path
    pname = self.p.name
    it = str(it).rjust(5, '0')
    fid = 'RawGrad'
    self.name = '{}-CP{}-{}.vtr'.format(pname, it, fid)
    self.fname = self.path + self.name        
class ProjGradShot(ProjGrad):
  pass
