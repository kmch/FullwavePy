"""
Interface for different FWI software.
"""
from abc import ABC, abstractmethod
from autologging import logged, traced

class FwiCode(ABC):
  """
  A generic FWI code.
  """
  def __init__(self, proj, **kwargs):
    self.proj = proj
    self._init_exe()
    self._init_io(**kwargs)
    # self._init_pre()
    self._init_rnf()
    self._init_ske()
  @classmethod
  def create(cls, ID, *args, **kwargs):
    if ID == 'sofi3d':
      return Sofi3d(*args, **kwargs)
    elif ID == 'fw3d-688':
      return Fw3d688(*args, **kwargs)
    elif ID == 'fw3d-688-test':
      return Fw3d688test(*args, **kwargs)        
    elif ID == 'fw3d-690':
      return Fw3d690(*args, **kwargs)     
    elif ID == 'fw3d-690-p2v':
      return Fw3d690p2v(*args, **kwargs)          
    elif ID == 'fw3d-728':
      return Fw3d728(*args, **kwargs)
    else:
      raise ValueError
  @abstractmethod
  def _init_exe(self):
    pass
  @abstractmethod
  def _init_io(self):
    pass
  @abstractmethod
  def _init_pre(self):
    """
    Init the pre-processor.
    """
    pass    
  @abstractmethod
  def _init_rnf(self):
    pass
  @abstractmethod
  def _init_ske(self):
    pass    
  @abstractmethod
  def run(self):
    pass        
class Fw3d(FwiCode):
  """
  Generic Fullwave3d interface.

  """
  path_rds = '~/PhD/fullwave3D/'
  path_kmc = '/home/kmc3817/light_PhD/fullwave3D/'
  def _init_io(self, **kwargs):
    prefix = self.proj.name
    io = self.proj.io.lower()
    # note, proj.io is specific to Fullwave3D
    # so its parsing and further selection
    # should be done here, not in FwiIo.create
    fw3d = ['fw3d']
    segy = ['sgy', 'segy']
    args = [prefix]
    kwargs = {}
    if io in segy:
      self.io = FwiIoSegy(*args, **kwargs)
    elif io in fw3d:
      self.io = FwiIoFw3d(*args, **kwargs)
  def _set_fnames(self):
    pass
  def run(self, proj, no, machine='local', **kwargs):
    from fwipy.shell.generic.api import ShellFactory
    sh = ShellFactory.create('linux')        
    assert machine == 'local'
    if machine == 'local':
      if no is None:
        no = 0
        sh.remove(proj.i.bash.no[no].fname)
        sh.remove(proj.o.o.no[no].fname)
        sh.remove(proj.o.e.no[no].fname)
      proj.i.bash.no[no].prep(cat=0)
      proj.i.bash.no[no].run()
class Fw3d688(Fw3d):
  def _init_exe(self):
    self.proj.exe['fullwave'] = '~/PhD/fullwave3D/rev688/bin/fullwave3D.exe'
    self.proj.exe['fullwave_local'] = '/home/kmc3817/light_PhD/fullwave3D/rev688/bin/fullwave3D.exe'
  def _init_pre(self):
    self.pre = SegyPrep310()     
  def _init_rnf(self):
    pass
  def _init_ske(self):
    pass
class Fw3d688test(Fw3d):
  def _init_exe(self):
    self.proj.exe['fullwave'] = '~/PhD/fullwave3D/rev688_test/bin/fullwave3D.exe'
    self.proj.exe['fullwave_local'] = '/home/kmc3817/light_PhD/fullwave3D/rev688_test/bin/fullwave3D.exe'
    self.proj.exe['segyprep'] = '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16'
  def _init_pre(self):
    self.pre = SegyPrep316()     
  def _init_rnf(self):
    pass
  def _init_ske(self):
    pass
class Fw3d690(Fw3d688):
  def _init_exe(self):
    self.proj.exe['fullwave'] = '~/PhD/fullwave3D/rev690/bin/fullwave3D.exe'
    self.proj.exe['fullwave_local'] = '/home/kmc3817/light_PhD/fullwave3D/rev690/bin/fullwave3D.exe'
  def _init_pre(self):
    self.pre = SegyPrep310()   
    self.proj.exe['segyprep'] = '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16'
  def _init_rnf(self):
    pass
  def _init_ske(self):
    pass    
class Fw3d690p2v(Fw3d690):
  def _init_exe(self):
    self.proj.exe['fullwave'] = '~/PhD/fullwave3D/rev690_p2v/bin/fullwave3D.exe'
    self.proj.exe['fullwave_local'] = '/home/kmc3817/light_PhD/fullwave3D/rev690_p2v/bin/fullwave3D.exe'
class Fw3d728(Fw3d690):
  def _init_exe(self):
    self.proj.exe['fullwave'] = self.path_rds + '/rev728/bin/fullwave3D.exe'
    self.proj.exe['fullwave_local'] = self.path_kmc + 'rev728/bin/fullwave3D.exe'
  def _init_pre(self):
    self.pre = SegyPrep324()   
    self.proj.exe['segyprep'] = '/home/kmc3817/light_PhD/fullwave3D/rev728/bin/segyprep'
  def _init_rnf(self):
    self.proj.i.rnf = Rnf728(self.proj, self.proj.i.path)
  def _init_ske(self):
    self.proj.i.ske = Ske728(self.proj, self.proj.i.path)        
