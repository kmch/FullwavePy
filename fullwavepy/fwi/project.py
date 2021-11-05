"""
FWI projects.

Notes
-----
A project is defined by a single `input` to the FWI code
(in particular, a single set of runtime parameters).
This input can be fed into FWI code multiple times by 
running multiple `jobs` (e.g. when a single run has to be 
split up into smaller chunks, for computational reasons).

Two basic types of projects involve respectively:
- generation of synthetic data for QC purposes,
- inversion (either of synthetic or field data).
"""
from abc import ABC, abstractmethod
from autologging import logged, traced

from fullwavepy.fwi.throughput import Inp, Out
from nsh.generic import ShellFactory

class Project(ABC):
  """
  A generic FWI project.
  """
  def __init__(self, name, path='./', shell='linux'):
    self.args = [name, path]
    self.name = name
    self._init_shell(shell)
    self._init_path(path)
    self._init_directories()
  def init(self, **kwargs):
    self._init_params()
    self._init_input()
    self._create_aliases()
  def reinit(self):
    super().__init__(self.name, **self.kwargs)    
  def run(self, no=None, machine='local', **kwargs):
    self.fwicode.run(self, no, machine, **kwargs)
  # -----------------------------------------------------------------------------
  def _create_aliases(self):
    """
    Create shorthand aliases for faster typing.
    """
    self.i = self.inp
    # self.o = self.out  
  def _init_directories(self):
    pro_dir = self.path
    inp_dir = self.path + '/inp/'
    out_dir = self.path + '/out/'
    for d in [pro_dir, inp_dir, out_dir]:
      if not self.shell.exists(d):
        self.shell.mkdir(d)
  def _init_experiment(self, ID, **kwargs):
    self.exp = Experiment.create(ID, **kwargs)
  def _init_envvars(self, **kwargs):
    if 'snap' in kwargs and 'env' in kwargs:
      kwargs['env']['SLAVES_WAVEFIELDSVTR'] = kwargs['snap']        
    elif 'snap' in kwargs:
      kwargs['env'] = {'SLAVES_WAVEFIELDSVTR': kwargs['snap']}        
    return kwargs
  def _init_fwicode(self, ID, **kwargs):
    self.fwicode = FwiCode.create(ID, self, **kwargs)
  def _init_ibmcode(self, ID):
    if ID is not None:
      from fwipy.ibm.api import IbmCodeFactory
      self.ibmcode = IbmCodeFactory.create(ID)        
  def _init_input(self):
    self.inp = Inp(self)
    for param in self.params.all.values():
      param._init_input()
  def _init_params(self):
    self.params = ProjectParams(self)
  def _init_path(self, path):
    self.path = path + '/{name}/'.format(name=self.name)
    # self.__log.debug('Project path set to: ' + proj.path)
  def _init_precode(self, ID):
    self.precode = PreCodeFactory.create(ID, self)
  def _init_shell(self, shell):
    self.shell = ShellFactory.create(shell)
  # -----------------------------------------------------------------------------
class ProjectParams:
  """
  Generic parameters independent of 
  a project type, FWI code, etc.

  Examples: 
  - anisotropy ('none', 'vti', ...)
  - ibm (True, False)

  """
  def __init__(self, proj):
    self.all = {
      'aniso' : Anisotropy(proj)
    }
  def read(self, **kwargs):
    pass
  def _read_from_kwargs(self, **kwargs):
    pass
  def _read_from_file(self):
    pass
class Param(ABC):
  """
  API for creating new
  parameters.
  """
  def __init__(self, proj):
    self.proj = proj
  @abstractmethod
  def _init_input(self):
    pass
class Anisotropy(Param):
  def _init_input(self):
    self.proj.inp.true_delta = None
    self.proj.inp.true_epsil = None
    self.proj.inp.true_theta = None
class Density(Param):
  def init(self, **kwargs):
    rho = kwargs.get('rho', 'gardner')
class Equation(Param):
  def init(self, **kwargs):
    eq = kwargs.get('eq', 'ait')
    eq = eq.lower()
    if eq == 'ait':
      self.proj.models = ['vp', 'rho']
class Problem(Param):
  def _init_input(self, proj):
    self._init_models()
  def _init_models(self):
    self.proj.inp.models = None
    for m in self.proj.models:
      self.tvp = Model('vp')
      self.all.append(self.tvp)
