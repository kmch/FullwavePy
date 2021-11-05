"""
This module defines basic project types ProjSyn and ProjInv 
along their parent Proj that serves only 
as a library of common features and should not be called directly.

Other common features are implemented in ..generic.io.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer, widgets
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.system import bash, exists, current_dir
from fullwavepy.plot.generic import figure

# FIXME: tmp
from fullwavepy.project.generic.io import _InpPreparer 

# -------------------------------------------------------------------------------
# Mix-ins 
# -------------------------------------------------------------------------------
@traced
@logged
class _FwiRunner(object):
  def run(self, *args, **kwargs):
    proj = self
    self.fwicode.run(proj, *args, **kwargs)
  def run_old(self, no=0, runner='bash', **kwargs):
    self.__log.info('Running a job no. %s' % no)
    self.out.rm(ls=False)
    if runner == 'bash':
      runner = self.inp.bash.no[no]
    else:
      raise TypeError('Runner: %s' % str(runner))

    kwargs['cat'] = False
    runner.prep(**kwargs)
    runner.run(**kwargs)
    self.reinit()
@traced
@logged
class _ProjSyncer(object):
  def rsync(self, *args, **kwargs):
    """
    Sync with remote host.

    """
    self.__log.debug('Assuming inp and out rsync have the same (kw)args...')
    self.inp.rsync(*args, **kwargs)
    self.out.rsync(*args, **kwargs)
  
  # -----------------------------------------------------------------------------
@traced
@logged
class _ProjLister(object):
  def ls(self, **kwargs):
    """
    List content of project dirs.
    
    """
    self.inp.ls(**kwargs)
    self.out.ls(**kwargs)

  # -----------------------------------------------------------------------------
@traced
@logged
class _ProjMeta(object):
  def lids(self, overwrite=False, lid_hw='ep'):
    if self.problem == 'synthetic':
      f = self.i.ose
    elif self.problem == 'tomography':
      f = self.i.obs
    else:
      raise NotImplementedError
    head = f.read_header(overwrite=overwrite)
    return list(sorted(head[lid_hw].unique()))
# -------------------------------------------------------------------------------
# Generic abstract project
# -------------------------------------------------------------------------------
@traced
@logged
class Proj(_ProjMeta, _ProjSyncer, _ProjLister, _FwiRunner):
  """
  Abstract parent class of all full-waveform projects.
  It initializes all the project-related objects.

  """
  def __init__(self, name, **kwargs):
    """
    Initialize a Fullwave project.
    
    Parameters
    ----------
    name : str 
      Name of the project, 
      prefix of all the input/output files.

    Notes
    -----
    One usually initializes various projects in the 
    notebook many times => it needs to be fast.
    
    """
    from fullwavepy.project.generic.io import ProjInput, ProjOutput
    from fullwavepy.project.generic.au import (ProjPath, ProjDirs, ProjDef, ProjGeometry, 
                                               ProjEnv, ProjSgyMapp, ProjBox, ProjCluster,
                                               ProjBaseFiles)
    from fullwavepy.project.files.text.misc import InfoFile, MetaDataProjFile
    from fullwavepy.project.files.text.runfiles import Runfile
    
    self.kwargs = dict(kwargs) # crucial to inherit from project to project!
    # dict prevents deleting the path
    self.parent_dir = current_dir()
    self.name = name
    self.proj = self # USED IN wrapper_widgets (self.proj.dims)
    self.exe = kw('exe', {}, kwargs)
    #self.__log.debug('Paths to executables: ' + str(self.exe))
    #if len(self.exe) == 0:
    #  self.__log.warning('Empty paths to executables (exe dictionary)')
    
    # IMMERSED BOUNDARY
    #self.immerse = kw('immerse', False, kwargs)
    
    self.ibfs = kw('ibfs', False, kwargs)
    self.__log.debug('Is this project using immersed-boundary? Yes if one. %s' % str(self.ibfs))

    
    ProjPath(self, **kwargs)
    del_kw('path', kwargs)
    ProjDirs(self, **kwargs)    
    
    self.info = InfoFile(self, **kwargs)
    ProjDef(self, **kwargs)
    
    self.sgyhw = ProjSgyMapp(self, **kwargs)
    self.env = ProjEnv(self, **kwargs)
    self.cluster = ProjCluster(self, **kwargs)

    
    meta = kw('meta', None, kwargs)
    self.meta = MetaDataProjFile(self, self.path, **kwargs)
    
    self.inp = ProjInput(self, **kwargs)
    self.i = self.inp # ALIAS
    self.inp.init(**kwargs)

    self.geom = ProjGeometry(self, **kwargs) # OUTPUT USES IT
    self.pbox = ProjBox(self, **kwargs)
    self.base = ProjBaseFiles(self, **kwargs)    


    self.out = ProjOutput(self, **kwargs)
    self.o = self.out # ALIAS
    self.out.init(**kwargs)
  def reinit(self, path=None):
    """
    An interim fix to update some project objects
    automatically by the workflow.
    """
    if path is not None:
      self.kwargs['path'] = path
    self.__init__(self.name, **self.kwargs)
# -------------------------------------------------------------------------------
# Two basic projects to use in practice
# -------------------------------------------------------------------------------
@traced
@logged
class ProjSyn(_InpPreparer, Proj):
  """
  Generate synthetic data.
  
  """
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
      self.__log.warning('Setting env to default: %s' % str(kwargs['env']))
    
    super().__init__(name, **kwargs)
  def init_input(self, **kwargs):
    """
    
    We don't probably need template idx.
    """
    from fullwavepy.project.files.gridded.models import ModelFileSgy, ModelFileVtr
    from fullwavepy.project.files.datalike.sgy import DataFileSgy, TemplateFileSgy
    from fullwavepy.project.files.datalike.ttr import TemplateFileTtr
    
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
  def init_output(self, **kwargs): 
    """
    This should be called by
    proj.output(...)
    
    Notes
    -----
    Forward wavefield is common for both syn and inv => io

    """
    from fullwavepy.project.files.datalike.sgy import SynDataFileSgy
    from fullwavepy.project.files.datalike.ttr import SynDataFileTtr
    from fullwavepy.project.files.other.index import SynIndexFileSgy, SynIndexFileTtr
    from fullwavepy.project.lists.extra import ForwardWavefieldFileList
    
    if self.io == 'sgy':
      SynDataClass = SynDataFileSgy
      SynIndexClass = SynIndexFileSgy
    elif self.io == 'fw3d':
      SynDataClass = SynDataFileTtr
      SynIndexClass = SynIndexFileTtr
    else:
      raise ValueError('Unknown io: ' + self.io)
    
    self.out.syn = SynDataClass(self, self.out.path, **kwargs)
    self.out.fw = ForwardWavefieldFileList(self.proj, it_max=1, **kwargs) 
  def _TODELETE_prep_inp_from_scratch(self, run=False, **kwargs):  
    """
    """
    # PARSE kwargs
    # wavelet = kw('wavelet', True, kwargs)
    # rawseis = kw('rawseis', True, kwargs)
    # sp = kw('sp', True, kwargs)
    # runfile = kw('runfile', True, kwargs)
    rm = kw('rm', True, kwargs)
    run = kw('run', True, kwargs)
    #plot = kw('plot', True, kwargs)
    cat = kw('cat', True, kwargs)
    #anim = kw('anim', False, kwargs) #NOTE
    check = kw('check', True, kwargs)
    if rm:
      self.inp.rm(**kwargs)
      # self.out.rm(**kwargs)
    
    # PREPARE ALL THE OBJECTS, ONE AT A TIME
    obj_ids = ['rsg', 'tvp']
    for obj_id in obj_ids:
      obj = getattr(self.inp, obj_id)
      self.__log.info('Preparing %s...\n' % obj.name)
      
      # MAKE SURE THE BASE OBJECT EXISTS
      base_obj_id = 'base_%s' % obj_id
      assert hasattr(self, base_obj_id) # FIXME: p.base_tvp SHOULD BE MOVED
      
      # PREPARE THE OBJECT FROM THE BASE ONE
      arr = getattr(self, base_obj_id)
      # FIXME: TMP
      if obj_id == 'tvp':
        arr = arr.carve(self.box)
      obj.prep(arr)

    
    self.__log.warning('RETURNING FOR NOW TMP')
    return

    # -----------------------------------------------------------------------------
    # SOURCE WAVELET
    # -----------------------------------------------------------------------------
    if wavelet:
      assert hasattr(self, 'base_rsg')
      self.__log.info('\nPreparing the wavelet...\n') 
      self.inp.rsg.prep(self.base_rsg)
      #if plot > 0:
        #self.inp.wavelet.Plot(**kwargs)


    # -----------------------------------------------------------------------------
    # RAWSEIS TEXT FILE
    # -----------------------------------------------------------------------------    
    if rawseis:
      print(('\n\n' + 'Preparing RawSeis.txt...' + '\n\n'))      
      self.inp.rawseis.Create(**kwargs)
      if cat > 0:
        self.inp.rawseis.Cat(**kwargs)      

    # -----------------------------------------------------------------------------
    # SEGYPREP RUNFILE & RUN
    # -----------------------------------------------------------------------------    
    if sp:
      print('\n\n' + 'Preparing SegyPrep.key... (recipr=1)' + '\n\n')
      self.inp.sp.Create(reciprocity=1, **kwargs)
      if cat > 0:
        self.inp.sp.Cat(**kwargs)
      if run > 0:
        self.inp.sp.Run(**kwargs)

    
    # -----------------------------------------------------------------------------
    # SOURCES & RECEIVERS QC
    # -----------------------------------------------------------------------------          
    if cat > 0:
      self.inp.sr.Cat(**kwargs)
    #if plot > 0:
      #kwargs['dx'] = self.dx
      #self.inp.sr.Plot_All(**kwargs)
      #self.inp.sr.Show_Offsets(**kwargs)
    
    
    if runfile:
      print(('\n\n' + 'Preparing the runfile...' + '\n\n'))
      self.inp.runfile.Create(**kwargs)
      if cat > 0:
        self.inp.runfile.Cat(**kwargs)
     
    if verbos > 0:
      print('Note, PBS must be created manually')
    
    if check > 0:
      print('Checking all the input with Fullwave3D...')
      self.inp.Check()
  def _prep_inp_from_another(self, other_proj, **kwargs):
    """
    Prepare the input based on another syn project.
    
    run is False by default because SP will operate on
    full receiver gathers => slow.
    
    FIXME: THIS WILL FAIL AS SOON AS DIMS OF pold ARE DIFFERENT!
    (replace dupl with prep)

    """
    assert pold.dims == self.dims

    # THESE ARE ASSUMED TO BE IDENTICAL
    for a in ['rsg', 'tvp', 'rse']:
      sobj = getattr(self.i, a)
      oobj = getattr(pold.i, a)
      sobj.dupl(oobj.fname)

    # THIS CAN HAVE CHANGES, E.G. TOTALTIME
    self.i.sp.prep(reciprocity=bool(pold.i.sp.read()['reciprocity']))
    
    # if run:
    #   self.i.sp.run()
    #   self.i.rnf.prep(**kwargs)
    # else:
    #   self.__log.warning('You still need to call i.sp.run() and i.rnf.prep()!')
  ##@widgets('cmap', 'x', 'y', 'z')  
  def plot_input(self, widgets=False, **kwargs):
    """
  
    Notes
    -----
    Reading files is slow only when generating the figure.
    After that each object has its array read and interaction
    should go then smoothly.
    
    """
    from matplotlib.gridspec import GridSpec
    
    gs = GridSpec(3,2, \
      width_ratios=[3,1], height_ratios=[1,1,1])
    figsize = (kw('figsize_x', 15, kwargs), \
      kw('figsize_y', 10, kwargs))
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
    kwargs['nslices'] = 3
    self.i.tvp.plot(**kwargs) 
    
    #if kwargs['sources']:
      #print('sources on')
      #self.i.s.plot_3slices(**kwargs)
  def plot_output(self, **kwargs):
    self.out.syn.plot(**kwargs)
  def plot_output_FIXME(self, figsize, layers, **kwargs):
    """
    """
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=figsize)

    gs = GridSpec(4, 2, height_ratios=[.5, 1, 1, 1]) #left=0, right=1, hspace=.1, wspace=.4) # width_ratios=[1, 2], )
    
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
@traced
@logged
class ProjInv(Proj):
  """
  Inversion.
  
  """
  def __init__(self, name, **kwargs):
    """
    
    """
    kwargs['problem'] = 'tomography'
    super().__init__(name, **kwargs)
  def init_input(self, **kwargs):
    """
    
    """
    from fullwavepy.project.files.gridded.models import ModelFileVtr, ModelFileSgy
    from fullwavepy.project.files.datalike.ttr import ObsDataFileTtr
    from fullwavepy.project.files.datalike.sgy import ObsDataFileSgy
    from fullwavepy.project.files.other.index import ObsIndexFileSgy, ObsIndexFileTtr
    from fullwavepy.project.files.text.hed import ObsHedFile

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
    from fullwavepy.project.lists.extra import CPFileList, DumpFileList, \
      ForwardWavefieldFileList, BackpropWavefieldFileList
    from fullwavepy.project.files.text.misc import LastCheckpointFile
    from fullwavepy.project.generic.qc import Functional
    
    self.out.lastcp = LastCheckpointFile(self, self.out.path, **kwargs)  
    it_max = self.lastcp if self.lastcp > 0 else 1

    self.out.fw = ForwardWavefieldFileList(self.proj, it_max, **kwargs) 
    self.out.bw = BackpropWavefieldFileList(self.proj, it_max, **kwargs)

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
    
    #self.__log.warning('Disabled init of cp file for now (until debug)')
    #
    for attr, [file_id, file_start, file_class] in cpnts.items():
      self.__log.debug('attr=%s, file_id=%s, file_start=%s, file_class=%s' % \
        (attr, file_id, file_start, file_class))
      setattr(self.out, attr, 
              CPFileList(self, file_class, file_id, file_start, **kwargs))    
     # CPFileList init: def __init__(self, proj, FileClass, file_id, file_start, **kwargs)


    dumps = {'dumpcomp': ['SLAVES_DUMPCOMPARE', DumpCompareFile],
            #  'dumpdat': ['SLAVES_DUMPDAT', DumpDataFile],
            }
    for attr, [file_id, file_class] in dumps.items():
      self.__log.debug('attr, file_id, file_class %s %s %s ' %
                       (attr, file_id, file_class))
      if (file_id in self.env.var) and (self.env.var[file_id] is not None):
        setattr(self.out, attr, DumpFileList(self, file_class, file_id, **kwargs))
        #DumpFileList init: def __init__(self, proj, FileClass, file_id, fwd=1, **kwargs)

    # ADD ALIASES MANUALLY
    self.out.dc = self.out.dumpcomp
  def prepare_input(self, syn_proj, run=True, process=True, **kwargs):
    """
    Prepere inversion input based on the synthetic project from which 
    first breaks will be extracted.
    
    """
    if process:
      assert 'filt_kwargs' in kwargs
      assert 'mute_kwargs' in kwargs    
    
    cat = kw('cat', 1, kwargs)

    # Check if extra nodes changed
    run_fsprep = False
    extra_kw = ['e_abs', 'e_top']
    for ekw in extra_kw:
      if ekw in kwargs: # and kwargs[ekw] != syn_proj.i.rnf.read()[ekw]:
        # self.__log.info('%s found in kwargs => will run fsprep.' % ekw)
        # run_fsprep = True
        self.__log.warning('%s found in kwargs => you should run fsprep ' % str(ekw)+\
          'if it is different from synthetic project.' )
        break
    # If not, just copy the GhostData file
    if not run_fsprep:
      self.i.ghb.dupl(syn_proj.i.ghb.fname)  
      self.i.sdata.dupl(syn_proj.i.sdata.fname)
      self.i.rdata.dupl(syn_proj.i.rdata.fname)
    
    self.i.fs.dupl(syn_proj.i.fs.fname) # FIXME: ONLY IF IMMERSE
    self.i.rsg.dupl(syn_proj.i.rsg.fname)
    self.i.svp.dupl(syn_proj.i.tvp.fname)
    self.i.obs.dupl(syn_proj.i.ose.fname) # what for?
    self.i.obs.raw.dupl(syn_proj.i.ose.fname) # NOTE
    self.i.rse.prep(fnames=[self.i.obs.raw.name])
    self.i.sp.prep(**dict(syn_proj.i.sp.read(), problem=self.problem), cat=cat)
    if run:
      self.i.sp.run(cat=cat)
      self.i.rnf.prep(**kwargs)
      if run_fsprep:
        self.i.fs.run(**kwargs)
        # self.i.ght.read()
        # ghs = self.i.ght.ghosts
        # iss = self.i.ght.isects
        # ine = self.i.ine.read()
        # srcs = p.i.s.read()
        # recs = p.i.r.read()
        # r_hicks = 2
        # rmax = 5        
    else:
      self.__log.warning('You still need to:\n' + \
        '1. Run %s with p.i.sp.run()\n' % self.i.sp.name + \
        '2. Prepare %s with p.i.rnf.prep()\n' % self.i.rnf.name + \
        '3. Run FsPrep with p.i.fs.run() (if extra nodes changed).')
      
    if process:
      self.i.obs.process(filt_kwargs=kwargs['filt_kwargs'], 
                         mute_kwargs=dict(**kwargs['mute_kwargs'], 
                                          syn_file=syn_proj.o.syn))
    else:
      self.__log.warning('You need to i.obs.process!')
  def prepare_output(self, **kwargs):
    """
    Read all the array to make interactive plots
    run smoothly.
    
    """
    pass
  def plot_input(self, **kwargs):
    print('Nothing yet.')
  ##@widgets('cmap', 'sids', 'run_ids', 'it')
  def plot_output(self, widgets=False, **kwargs):
    """
    """
    it = kw('it', 1, kwargs)
    #sid = kwargs['sid', kwargs[
    
    self.prepare_output(**kwargs)

    fig = figure(**kwargs)
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

