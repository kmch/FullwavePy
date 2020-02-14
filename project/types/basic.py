"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for peremoveission writing to k.chrapkiewicz17@imperial.ac.uk.

Notes
-----
Python devs usually don't write nested classes, neither do I.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.system import bash, exists
from fullwavepy.plot.generic import new_figure


# -------------------------------------------------------------------------------


@traced
@logged
class Proj(object):
  """
  fullwavepy.generic Fullwave project.
  
  """
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, name, **kwargs):
    """
    Initialize a Fullwave project.
    
    Parameters
    ----------
    name : str 
      Name of the project. It will be 
      a prefix of all the input/output files.
    kwargs : 
      Print help of the called functions to get info 
      about what kwargs are parsed.
    
    Returns
    -------
    None 
    
    Notes
    -----
    One usually initializes various projects in the 
    notebook many times.
    
    """
    from fullwavepy.project.io import ProjInput, ProjOutput
    from fullwavepy.project.files.misc import InfoFile, MetaDataProjFile
    from fullwavepy.project.files.runfiles import Runfile
    from fullwavepy.project.files.gridded.surfaces import TopographyFile
    from fullwavepy.project.aux import (ProjPath, ProjDirs, ProjDef, ProjGeometry, 
                                        ProjEnv, ProjSegyMapp, ProjBox)    
    self.name = name
    self.proj = self # USED IN wrapper_widgets (self.proj.dims)
    
    self.exe = kw('exe', {}, kwargs)
    #self.__log.debug('Paths to executables: ' + str(self.exe))
    #if len(self.exe) == 0:
    #  self.__log.warn('Empty paths to executables (exe dictionary)')
    
    
    ProjPath(self, **kwargs)
    del_kw('path', kwargs)
    ProjDirs(self, **kwargs)    
    
    self.info = InfoFile(self, **kwargs)
    ProjDef(self, **kwargs)
     
    self.sgy = ProjSegyMapp(self, **kwargs)
    self.env = ProjEnv(self, **kwargs)
    
    self.cluster = kw('cluster', 'cx1', kwargs)
    self.__log.info('PBS scripts will be prepared for the ' + 
                    str(self.cluster) + ' cluster')
    
    meta = kw('meta', None, kwargs)
    self.meta = MetaDataProjFile(self, self.path, **kwargs)
    
    self.inp = ProjInput(self, **kwargs)
    self.i = self.inp # ALIAS
    self.inp.init(**kwargs)

    self.geom = ProjGeometry(self, **kwargs) # OUTPUT USES IT
    self.pbox = ProjBox(self, **kwargs)
    
    self.out = ProjOutput(self, **kwargs)
    self.o = self.out # ALIAS
    self.out.init(**kwargs)
    
    topo = kw('topo', None, kwargs)
    self.topo = TopographyFile(self, self.inp.path, dupl=topo, **kwargs)  
    
  # -----------------------------------------------------------------------------
  
  def ls(self, **kwargs):
    """
    List content of project dirs.
    
    """
    self.inp.ls(**kwargs)
    self.out.ls(**kwargs)

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


@traced
@logged
class ProjSyn(Proj):
  """
  Generate synthetic data.
  
  """
  
  # -----------------------------------------------------------------------------  
  
  def __init__(self, name, **kwargs):
    """
    
    """
    kwargs['problem'] = 'synthetic'
    
    # REPEATS aux/ProjEnv...
    #env = kw('env', {}, kwargs)
    #kw('SLAVES_WAVEFIELDSVTR', None, env)
    
    key = 'SLAVES_WAVEFIELDSVTR'
    default = -1000
    
    if 'env' in kwargs:
      kwargs['env'][key] = kwargs['env'].get(key, default)
    else:
      kwargs['env'] = {key : default}
    
    super().__init__(name, **kwargs)
    
  # -----------------------------------------------------------------------------
  
  def init_input(self, **kwargs):
    """
    
    We don't probably need template idx.
    """
    from fullwavepy.project.files.gridded.models import ModelFileSgy, ModelFileVtr
    from fullwavepy.project.files.datalike.sgy import DataFileSgy
    from fullwavepy.project.files.templates import TemplateFileSgy, TemplateFileTtr
    
    if self.io == 'sgy':
      ModelClass = ModelFileSgy
      ObsDataClass = DataFileSgy
      TemplateClass = TemplateFileSgy
      
    elif self.io == 'fw3d':
      ModelClass = ModelFileVtr
      ObsDataClass = DataFileSgy # OutSeis FOR ioapi IS STILL .sgy!
      TemplateClass = TemplateFileTtr
    else:
      raise ValueError('Unknown io: ' + self.io)
    
    self.inp.truevp = ModelClass('TrueVp', self, self.inp.path, **kwargs)
    self.i.tvp = self.inp.truevp # ALIAS
    
    self.inp.outseis = ObsDataClass('OutSeis', self, self.inp.path, **kwargs)
    self.i.ose = self.inp.outseis
    
    self.inp.template = TemplateClass(self, self.inp.path, **kwargs)
    self.i.tmpl = self.i.template
    
    
    #self.inp.template_hed = TemplateHedClass(template_suffix, self, self.inp.path, 
                                      #**kwargs)
    
    #self.inp.template_idx = TemplateIndexClass(template_suffix, self, self.inp.path, 
                                      #**kwargs)
    
  # -----------------------------------------------------------------------------  
  
  def init_output(self, **kwargs): 
    """
    This should be called by
    proj.ouplot(...)

    Notes
    -----
    Leading underscore - convention 
    for 'private' methods.
    
    """
    from fullwavepy.project.files.datalike.sgy import SynDataFileSgy
    from fullwavepy.project.files.datalike.ttr import SynDataFileTtr
    from fullwavepy.project.files.index import SynIndexFileSgy, SynIndexFileTtr
    #, WavefieldFiles
    
    if self.io == 'sgy':
      SynDataClass = SynDataFileSgy
      SynIndexClass = SynIndexFileSgy
    elif self.io == 'fw3d':
      SynDataClass = SynDataFileTtr
      SynIndexClass = SynIndexFileTtr
    else:
      raise ValueError('Unknown io: ' + self.io)
    
    self.out.syn = SynDataClass(self, self.out.path, **kwargs)
    self.out.syn_idx = SynIndexClass(self, self.out.path, **kwargs)    
    #try:
    #  self.out.syn.files(timer=True)
    #except OSError as err_message: 
    #  self.__log.warning(str(err_message))
    
  # -----------------------------------------------------------------------------
  
  def prepare_input(self, pold, run=False, **kwargs):
    """
    Prepare the input based on another syn project.
    
    run is False by default because SP will operate on
    full receiver gathers => slow.
    
    """
    # THESE ARE ASSUMED TO BE IDENTICAL
    for a in ['rsg', 'tvp', 'rse']:
      sobj = getattr(self.i, a)
      oobj = getattr(pold.i, a)
      sobj.dupl(oobj.fname)

    # THIS CAN HAVE CHANGES, E.G. TOTALTIME
    self.i.sp.prep(reciprocity=bool(pold.i.sp.read()['reciprocity']))
    
    if run:
      self.i.sp.run()
      self.i.rnf.prep(**kwargs)
    else:
      self.__log.warn('You need to i.sp.run and i.rnf.prep!')

  # -----------------------------------------------------------------------------
      
  @widgets('cmap', 'x', 'y', 'z')  
  def plot_input(self, widgets=False, **kwargs):
    """
  
    Notes
    -----
    Reading files is slow only when generating the figure.
    After that each object has its array read and interaction
    should go then smoothly.
    
    """
    from matplotlib.gridspec import GridSpec
    

    
    gs = GridSpec(3,2, width_ratios=[3,1], height_ratios=[1,1,1])
    figsize = (kw('figsize_x', 15, kwargs), kw('figsize_y', 10, kwargs))
    fig = plt.figure(figsize=figsize)    
    
    # SIGNATURES OF ALL SOURCES
    fig.add_subplot(gs[0,0])
    self.i.sgn.plot()
    
    # RAW SIGNATURE
    fig.add_subplot(gs[0,1])
    self.i.rsg.plot(**kwargs)
    
    # BATHY      
    #fig.add_subplot(gs[1,:])
    
    # MODEL SUBGRIDSPEC
    gs_tvp = gs[1:, :].subgridspec(2,2)
    kwargs['fig'] = fig
    kwargs['gs'] = gs_tvp
    self.i.tvp.plot_3slices(**kwargs) 
    
    #if kwargs['sources']:
      #print('sources on')
      #self.i.s.plot_3slices(**kwargs)
    
  # -----------------------------------------------------------------------------
  
  def plot_output(self, figsize, layers, **kwargs):
    """
    """
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=figsize)

    gs = GridSpec(4, 2, height_ratios=[.5, 1, 1, 1])
#                  left=0, right=1, hspace=.1, wspace=.4) # width_ratios=[1, 2], )
    
    lays = []
    if 'tvp' in layers:
      lays.append(self.i.tvp)
    if 'fw' in layers:
      lays.append(self.o.fw)
    if 's' in layers:
      lays.append(self.i.s)
    if 'r' in layers:
      lays.append(self.i.r)    
    
    def annotate_2d():
      plt.gca().text(0.5, 0.5, "not available in 2D", va="center", ha="center")    
    
    ax1 = fig.add_subplot(gs[0, 0])
    self.i.rsg.plot()
    
    ax2 = fig.add_subplot(gs[0:2, 1])
    self.o.syn.plot()
    
    ax3 = fig.add_subplot(gs[1:3, 0])
    for a in lays:
      if self.dim == '3d':
        a.plot(scoord='z', **kwargs)
        plt.gca().set_aspect('equal')
      else:
        annotate_2d()
    ax4 = fig.add_subplot(gs[2, 1])
    for a in lays:
      if self.dim == '3d':
        a.plot(scoord='x', **kwargs)
        plt.gca().set_aspect('equal')
      else:
        annotate_2d()
    
    ax5 = fig.add_subplot(gs[3, :])
    for a in lays:
      a.plot(scoord='y', **kwargs)   
      plt.gca().set_aspect('equal')
   
  # -----------------------------------------------------------------------------     


# -------------------------------------------------------------------------------


@traced
@logged
class ProjInv(Proj):
  """
  Inversion.
  
  """

  # -----------------------------------------------------------------------------  

  def __init__(self, name, **kwargs):
    """
    
    """
    kwargs['problem'] = 'tomography'
    super().__init__(name, **kwargs)

  # -----------------------------------------------------------------------------  
  
  def init_input(self, **kwargs):
    """
    
    """
    from fullwavepy.project.files.gridded.models import ModelFileVtr, ModelFileSgy
    from fullwavepy.project.files.datalike.ttr import ObsDataFileTtr
    from fullwavepy.project.files.datalike.sgy import ObsDataFileSgy
    from fullwavepy.project.files.index import ObsIndexFileSgy, ObsIndexFileTtr
    from fullwavepy.project.files.templates import ObsHedFile
    
    if self.io == 'sgy':
      ModelClass = ModelFileSgy
      ObsDataClass = ObsDataFileSgy
      ObsIndexClass = ObsIndexFileSgy
    elif self.io == 'fw3d':
      ModelClass = ModelFileVtr
      ObsDataClass = ObsDataFileTtr
      ObsIndexClass = ObsIndexFileTtr
    else:
      raise ValueError('Unknown io: ' + self.io)
    
    self.inp.startvp = ModelClass('StartVp', self, self.inp.path, **kwargs)
    self.inp.svp = self.inp.startvp # ALIAS
    self.inp.obs = ObsDataClass(self, self.inp.path, **kwargs)    
    
  # -----------------------------------------------------------------------------

  def init_output(self, **kwargs):
    """
    
    Notes
    -----
    Not sure if for the multi-param FWI
    you get more gradients etc. and what their names are.
    Check it!
    
    # FIXME ADD SLAVE GRADS AND PRECS 
    
    """
    from fullwavepy.project.files.gridded.models import ModelFileVtr, ModelFileSgy
    from fullwavepy.project.files.gridded.derivs import GradFile, PrecFile
    from fullwavepy.project.files.datalike.ttr import DumpDataFile, DumpCompareFile
    from fullwavepy.project.lists.extra import CPFileList, DumpFileList
    from fullwavepy.project.files.misc import LastCheckpointFile
    from fullwavepy.project.qc import Functional

    self.out.lastcp = LastCheckpointFile(self, self.out.path, **kwargs)  
  
    self.out.fit = Functional(self, **kwargs)

    if self.io == 'sgy':
      ModelClass = ModelFileSgy
    elif self.io == 'fw3d':
      ModelClass = ModelFileVtr   
    else:
      raise ValueError('Unknown io: ' + self.io)

       
    cpnts = {'vp': ['Vp', self.inp.startvp, ModelClass],
             'grad': ['Grad', None, GradFile], #NOTE: GradFile doesn't imply 'Grad'
             'prec': ['Prec', None, PrecFile],
             'rawgrad': ['RawGrad', None, GradFile], #NOTE
             'rawprec': ['RawPrec', None, PrecFile]
             }
    for attr, [file_id, file_start, file_class] in cpnts.items():
      setattr(self.out, attr, 
              CPFileList(self, file_class, file_id, file_start, **kwargs))    
    
    
    dumps = {'dumpdat': ['SLAVES_DUMPDAT', DumpDataFile],
             'dumpcomp': ['SLAVES_DUMPCOMPARE', DumpCompareFile],
            }
    for attr, [file_id, file_class] in dumps.items():
      self.__log.debug('attr, file_id, file_class %s %s %s ' %
                       (attr, file_id, file_class))
      if (file_id in self.env.var) and (self.env.var[file_id] is not None):
        setattr(self.out, attr, DumpFileList(self, file_class, file_id, **kwargs))
        
    # ADD ALIASES MANUALLY
    self.out.dc = self.out.dumpcomp

  # -----------------------------------------------------------------------------
  
  def prepare_input(self, syn_proj, run=True, process=True, **kwargs):
    """
    Prepere inversion input based on the synthetic project from which 
    first breaks will be extracted.
    
    """
    if process:
      assert 'filt_kwargs' in kwargs
      assert 'mute_kwargs' in kwargs    
    
    self.i.rsg.dupl(syn_proj.i.rsg.fname)
    self.i.svp.dupl(syn_proj.i.tvp.fname)
    self.i.obs.dupl(syn_proj.i.ose.fname)
    self.i.obs.raw.dupl(syn_proj.i.ose.fname)
    self.i.rse.prep(fnames=[self.i.obs.raw.name])
    self.i.sp.prep(reciprocity=bool(syn_proj.i.sp.read()['reciprocity']))
    if run:
      self.i.sp.run()
      self.i.rnf.prep(**kwargs)
    else:
      self.__log.warn('You need to i.sp.run and i.rnf.prep!')
      
    if process:
      
      self.i.obs.process(filt_kwargs=kwargs['filt_kwargs'], 
                         mute_kwargs=dict(**kwargs['mute_kwargs'], 
                                          syn_file=syn_proj.o.syn))
    else:
      self.__log.warn('You need to i.obs.process!')
                      
  # -----------------------------------------------------------------------------
  
  def prepare_output(self, **kwargs):
    """
    Read all the array to make interactive plots
    run smoothly.
    
    """
    sources_dict = self.i.s.read(**kwargs)

  # -----------------------------------------------------------------------------

  def plot_input(self, **kwargs):
    print('Nothing yet.')
   
  # -----------------------------------------------------------------------------    
  
  @widgets('cmap', 'sids', 'run_ids', 'it')
  def plot_output(self, widgets=False, **kwargs):
    """
    """
    it = kw('it', 1, kwargs)
    #sid = kwargs['sid', kwargs[
    
    self.prepare_output(**kwargs)

    fig = new_figure(**kwargs)
    nrows = 3
    short = 1
    tall = 3
    height_ratios = [tall]*nrows
    height_ratios[0] = short
    gs = fig.add_gridspec(nrows, 1, height_ratios=height_ratios)
    
    fig.add_subplot(gs[0,0])
    self.o.fit.plot(**kwargs)
    
    fig.add_subplot(gs[1,0])
    self.o.vp.it[it].plot(**kwargs)
    #self.o.dc.it[it][4144].plot(**kwargs)
    
    #run_ids, it, sid, freq,
    #plt.figure()
    #self.o.fit.plot(run_ids)
    #plt.figure()
    #self.o.dc.it[it][sid].plot()
    #plt.figure()
    #self.o.dc.it[it][sid].plot_phase(freq)
    #plt.figure()
    #self.i.svp.compare(self.o.vp.it[it], plt.gcf(), **kwargs)


    #self.i.sgn.read(**kwargs)
    #self.i.rsg.read(**kwargs)
    #self.i.tvp.read(**kwargs)
    #
    #
    #figsize = (kw('figsize_x', 8, kwargs), kw('figsize_y', 8, kwargs))
    #fig = plt.figure(figsize=figsize)    
    #
    ## SIGNATURES OF ALL SOURCES
    #fig.add_subplot(gs[0,0])
    #self.i.sgn.plot(**kwargs)
    #
    ## RAW SIGNATURE
    #fig.add_subplot(gs[0,1])
    #self.i.rsg.plot(**kwargs)
    #
    ## BATHY      
    #fig.add_subplot(gs[1,:])
    #
    ## MODEL SUBGRIDSPEC
    #gs_tvp = gs[2:, :].subgridspec(2,2)
    #kwargs['fig'] = fig
    #kwargs['gs'] = gs_tvp
    #self.i.tvp.plot_3slices(**kwargs) 
    #
    #if kwargs['sources']:
    #  print('sources on')
      #self.i.s.plot_3slices(**kwargs)

  # -----------------------------------------------------------------------------    


# -------------------------------------------------------------------------------

