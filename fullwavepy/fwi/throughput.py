"""
Project throughput, i.e. 
input and output objects associated with 
a project.
"""
from abc import ABC, abstractmethod
from autologging import logged, traced

class Throughput(ABC):
  def __init__(self, proj):
    self.proj = proj
    # self._set_id()
    # self._set_attr_list()
    # self._copy_attrs()    
  @abstractmethod
  def _set_id(self):
    pass
  def _copy_attrs(self):
    """
    A. Martelli
    https://stackoverflow.com/questions/3818825/
    python-what-is-the-correct-way-to-copy-an-objects-
    attributes-over-to-another/3818861
    """
    if hasattr(self.proj, self.id):
      objfrom = self.proj.inp
      for n in self.attr_list:
        if hasattr(objfrom, n):
          v = getattr(objfrom, n)
          setattr(self, n, v)
class Inp(Throughput):
  # attr_list = ['path', 'fs', \
  #        'rsg', 'sgn', 'rse', 'sp', \
  #        's', 'r', 'ske', 'rnf', \
  #        'bash', 'pbs']
  # def __init__(self, *args, **kwargs):
  #   super().__init__(*args, **kwargs)
  # -----------------------------------------------------------------------------
  # @abstractmethod
  def create(self):
    pass
  def preprocess(self, **kwargs):
    """
    Notes
    -----
    The runfile keyword args are deep-copied,
    otherwise, blocks are modified by i.rnf.prep
    and 'nits' can be missing after that.
    """
    from copy import deepcopy
    kw_pre = kwargs.get('kw_pre', kwargs)
    kw_rnf = kwargs.get('kw_rnf', kwargs)
    self.proj.precode.prep_n_run(**kw_pre)
    self.proj.i.rnf.prep(**deepcopy(kw_rnf))  
    # if 'vacuum' in kw_rnf and kw_rnf['vacuum']:
    #   vacuum = 1
    #   kw_rnf = dict(kw_rnf, vacuum=0)
    # else:
    #   vacuum = 0
    # self.proj.i.rnf.prep(**rnf)
    # if 'ibfs' in rnf and float(rnf['ibfs']):
    #   self.proj.i.fs.immerse()
    # self.proj.i.rnf.modify(vacuum=vacuum)
  # -----------------------------------------------------------------------------
  # @abstractmethod
  def _set_attr_list(self):
    pass
  def _set_id(self):
    self.id = 'inp'
class Out(Throughput):
  attr_list = ['o', 'e'] # jo, je    
  def _set_id(self):
    self.id = 'out'

class InpSyn(Inp):
  def create(self, *args, **kwargs):
    if self.proj.exp is None:
      self._create_from_scratch(*args, **kwargs)
    else:
      self._create_from_experim(*args, **kwargs)
  def _create_from_experim(self, **kwargs):
    p = self.proj
    self.rsg.create(p.exp.wvl)
    self.tvp.create(p.exp.svp.carve(p.box))
  def _create_from_scratch(self, kw_tvp, kw_rsg, kw_fs, **kwargs):
    p = self.proj
    self.tvp.create(ProjMod(p).create(**kw_tvp))
    self.rsg.create(ProjWvl(p).create(**kw_rsg))
    if kw_fs is not None:
      self.fs.create(ProjFs(p).create(**kw_fs))
  def plot(self, **kwargs):
    fig = figure(15,5)
    gs = fig.add_gridspec(3,2, width_ratios=[3,1], height_ratios=[1,1,1])

    # SIGNATURES OF ALL SOURCES
    fig.add_subplot(gs[0,0])
    self.sgn.plot()
    # RAW SIGNATURE
    fig.add_subplot(gs[0,1])
    self.rsg.plot(**kwargs)
    # MODEL SUBGRIDSPEC
    gs_tvp = gs[1:, :].subgridspec(2,2)
    kwargs['fig'] = fig
    kwargs['gs'] = gs_tvp
    kwargs['nslices'] = 3
    self.tvp.plot(**kwargs)   
  def _set_attr_list(self):
    self.attr_list += ['tvp', 'ose']
class InpSynChecq(InpSyn):
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    from fullwavepy.project.files.gridded.models import ModelFileSgy
    p = self.proj
    m = self.tvp
    setattr(m, 'bckg', ModelFileSgy('BckgVp', p, m.path))
    setattr(m, 'anom', ModelFileSgy('AnomVp', p, m.path))
class InpInv(Inp):
  def create(self, projsyn=None):
    if projsyn is None:
      self._create_from_scratch()
    else:
      self._create_from_projsyn(projsyn)
  def _create_from_scratch(self):
    NIErr()
  def _create_from_projsyn(self, projsyn):
    p = self.proj
    self.rsg.dupl(projsyn.i.rsg.fname) ##### BUT NOT SIGNATURE!!!!!!!!!!!!
    self.svp.dupl(projsyn.i.tvp.fname)
    self.obs.dupl(projsyn.o.syn.fname)
    self.obs.idx.dupl(projsyn.o.syn.idx.fname)
    self.sp.dupl(projsyn.i.sp.fname)
    self.sp.modify(problem='tomography')
  def _set_attr_list(self):
    self.attr_list += ['svp', 'obs']# , 'rawgrad', 'rawprec', 'grad', 'prec']
class InpInvField(InpInv):
  def create(self, syn_proj, run=1, process=1, **kwargs):
    """
    Prepere inversion input based on the synthetic project from which 
    first breaks will be extracted.

    """
    cat = kwargs.get('cat', 1)

    self.rsg.dupl(syn_proj.i.rsg.fname)
    self.svp.dupl(syn_proj.i.tvp.fname)
    self.obs.dupl(syn_proj.i.ose.fname)
    self.obs.raw.dupl(syn_proj.i.ose.fname) # NOTE
    self.rse.prep(fnames=[self.obs.raw.name])
    self.sp.prep(**dict(syn_proj.i.sp.read(), \
              problem=self.proj.problem), cat=cat)
    if run:
      self.sp.run(cat=cat)
      self.rnf.prep(**kwargs)
    
    if process:
      self.obs.process(filt_kwargs=kwargs['filt_kwargs'], \
               mute_kwargs=dict(**kwargs['mute_kwargs'], 
               syn_file=syn_proj.o.syn))
class InpInvChecq(InpInv):
  def create(self, syn_proj, run=1, process=0, **kwargs):
    """
    Prepere inversion input based on the synthetic project from which 
    first breaks will be extracted.

    """
    cat = kwargs.get('cat', 1)

    self.rsg.dupl(syn_proj.i.rsg.fname)
    self.svp.dupl(syn_proj.i.tvp.bckg.fname) # NOTE
    self.obs.dupl(syn_proj.o.syn.fname) # NOTE
    self.obs.raw.dupl(syn_proj.o.syn.fname) # NOTE
    self.rse.prep(fnames=[self.obs.raw.name])
    self.sp.prep(**dict(syn_proj.i.sp.read(), \
              problem=self.proj.problem), cat=cat)
    if run:
      self.sp.run(cat=cat)
      self.rnf.prep(**kwargs)
    
    if process:
      raise NIErr('Implement muting only.')
      # self.obs.process(filt_kwargs=kwargs['filt_kwargs'], \
      #          mute_kwargs=dict(**kwargs['mute_kwargs'], 
      #          syn_file=syn_proj.o.syn))
class OutInv(Out):
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.grad = ProjGrad(self.proj)
