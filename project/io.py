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

from fullwavepy.generic.decor import widgets
from fullwavepy.generic.system import bash
from fullwavepy.generic.parse import kw


@traced
@logged
class ProjThroughput(object):
  """
  Parent class of input and output.
  
  Notes
  -----
  Not much implemented - it's left 
  for subclasses to define their own
  versions of methods.
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    self.proj = proj
    
  # -----------------------------------------------------------------------------
  
  def ls(self, **kwargs):
    o, e = bash('pwd', path=self.path)
    print('Content of ' + o)
    o, e = bash('ls -lth ' + self.path)
    print(o, e)
    
  # -----------------------------------------------------------------------------
  
  def dupl(self, project=None, input_path=None, project_name=None, **kwargs):
    """
    Duplicate from another project.
    
    Used to be:
    def dupl(self, input_path, project_name, **kwargs):
    
    """
    from fullwavepy.generic.system import get_files, duplicate
    from fullwavepy.generic.parse import path_leave, exten
    
    if project is None:
      if input_path is None or project_name is None:
        raise TypeError('Input either a project-object or 2 strings: input_path and project_name')
    else:
      input_path = project.inp.path 
      project_name = project.name
      
    
    fnames = get_files(input_path, project_name+'-*')
    for fname in fnames:
      nfname = path_leave(fname)
      
      if exten(fname) == 'pbs': 
        self.__log.warn('Skipping PBS script ' + fname + 
                        ' prepare it yourself.')
        continue

      if 'Job' in fname or (('Out' in fname or 'Err' in fname) and exten(fname) == 'log'):
        self.__log.warn('Skipping ' + fname)
        continue
    
      nfname = self.path + self.proj.name + nfname[len(project_name): ]
      duplicate(fname, nfname, **kwargs)

  # -----------------------------------------------------------------------------    

  def rm(self, **kwargs):
    """
    Remove all the files associated with this object.
    
    Notes
    -----
    It actually renames files (by adding a 'rm_' prefix) 
    which need to be rm-ed manually (for safety).
    
    """
    from fullwavepy.generic.system import get_files
    from fullwavepy.generic.parse import path_leave
    
    fnames = get_files(self.path, '*')
    for fname in fnames:
      fname = path_leave(fname)
      o, e = bash('mv ' + fname + ' rm_' + fname, path=self.path)

    self.ls()

  # -----------------------------------------------------------------------------
  
  #@widgets#FIXME
  #def plot_3slices(self, fig, gs=None, widgets=False, **kwargs):
  #  #FIXME BOILERPLATE
  #  if widgets: 
  #    figsize = (kw('figsize_x', 8, kwargs), kw('figsize_y', 8, kwargs))
  #    fig = plt.figure(figsize=figsize)
  #  if gs is None:
  #    gs = fig.add_gridspec(2,2)  
  #  ax1 = fig.add_subplot(gs[0,0])
  #  ax2 = fig.add_subplot(gs[0,1])
  #  ax3 = fig.add_subplot(gs[1,:])        
      
      
    #if kwargs['true_vp']:
    #  self.tvp.plot_3slices(fig, gs, **dict(kwargs, widgets=False))
    #
    ##if kwargs['fw']:
    #  #self.tvp.plot_3slices(fig, gs, **dict(kwargs, widgets=False))
    #
    #if kwargs['receivers']: # ALWAYS MORE THAN SOURCES (HOPEFULLY)
    #  self.r.plot_3slices(fig, **kwargs)
    #  
    #if kwargs['sources']:
    #  self.s.plot_3slices(fig, **kwargs)


  #def plot(self, **kwargs): #FIXME
  #  """
  #  An elegant way of plotting all plottables.
  #  
  #  """
  #  for attr in [getattr(self, i) for i in vars(self)]:
  #    try:
  #      attr.plot(**kwargs)
  #      plt.figure()
  #    except (AttributeError, FileNotFoundError) as err:
  #      self.__log.debug(err)

  # -----------------------------------------------------------------------------

  #def qc(self, **kwargs): #FIXME
  #  self.ls(**kwargs)
  #  self.plot(**kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ProjInput(ProjThroughput):
  """
  
  Notes
  -----
  It contains some necessary work-arounds
  because proj.init_input etc. refer to
  proj.inp so they must be called 
  AFTER the constructor (__init__), not from it.  
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    super(ProjInput, self).__init__(proj, **kwargs)
    self.path = proj.path + '/inp/'
  
  # ----------------------------------------------------------------------------- 
  
  def init(self, **kwargs):
    """
    Common for all project types.
    
    Notes
    -----
    Also defined a couple of aliases to save typing-time.
    
    """
    from fullwavepy.project.files.datalike.sgy import RawSignFile, SignatureFileSgy
    from fullwavepy.project.files.datalike.ttr import SignatureFileTtr
    from fullwavepy.project.files.gridded.surfaces import FsFile, ExtendedFsFile, InterpolFsFile
    from fullwavepy.project.files.ghost import GhostDataFileBin, GhostDataFileTxt
    from fullwavepy.project.files.geom import SourcesFile, ReceiversFile
    from fullwavepy.project.files.misc import RawSeisTxtFile, JobInfoFile
    from fullwavepy.project.files.runfiles import SegyPrepFile, Runfile, Skeleton
    from fullwavepy.project.files.submit import BashFile
    from fullwavepy.project.lists.basic import JobFileList 
    from fullwavepy.ioapi.cx1 import PbsFileCx1
    
    self.__log.debug('Initializing project-type-specific  input...')
    self.proj.init_input(**kwargs)

    self.__log.debug('Initializing generic-project input...')
    self.rawsign = RawSignFile(self.proj, self.path, **kwargs)
    self.rsg = self.rawsign # ALIAS
    self.rawseis = RawSeisTxtFile(self.proj, self.path, **kwargs)
    self.rse = self.rawseis # NOTE: WE DON'T USE RawSeis.sgy, OTHERWISE AMBIGUITY
    self.sp = SegyPrepFile(self.proj, self.path, **kwargs)
     
    #SignatureFile(self.proj, self.path, **kwargs)[io] 
     
    if self.proj.io == 'sgy':
      self.signature = SignatureFileSgy(self.proj, self.path, **kwargs)
    elif self.proj.io == 'fw3d':
      self.signature = SignatureFileTtr(self.proj, self.path, **kwargs)                                            
    else:
      raise ValueError('Wrong io: ' + self.proj.io)    
    self.sgn = self.signature
    
    self.fs = FsFile(self.proj, self.path, **kwargs)
    self.fse = ExtendedFsFile(self.proj, self.path, **kwargs)
    self.fsi = InterpolFsFile(self.proj, self.path, **kwargs)
    self.ghostbin = GhostDataFileBin(self.proj, self.path, **kwargs)
    self.ghb = self.ghostbin
    self.ghosttxt = GhostDataFileTxt(self.proj, self.path, **kwargs)
    self.ght = self.ghosttxt
    
    self.s = SourcesFile(self.proj, self.path, **kwargs)
    self.r = ReceiversFile(self.proj, self.path, **kwargs)
    
    self.skeleton = Skeleton(self.proj, self.path, **kwargs)
    self.ske = self.skeleton
    self.runfile = Runfile(self.proj, self.path, **kwargs)
    self.rnf = self.runfile # ALIAS
    self.bash = JobFileList(self.proj, self.path, BashFile, **kwargs)
    
    if self.proj.cluster == 'cx1':
      self.pbs = JobFileList(self.proj, self.path, PbsFileCx1, **kwargs)
    else:
      raise NotImplementedError('Unknown cluster: ' + self.proj.cluster)    
    
    self.jobinfo = JobFileList(self.proj, self.path, JobInfoFile, **kwargs)
    self.jinfo = self.jobinfo # ALIAS
    
  # -----------------------------------------------------------------------------
  
  def prepare(self, *args, **kwargs):
    self.proj.prepare_input(*args, **kwargs)

  def prep(self, *args, **kwargs):
    self.prepare(*args, **kwargs)    
    
  # ----------------------------------------------------------------------------- 
  
  def plot(self, *args, **kwargs):
    self.proj.plot_input(*args, **kwargs)

  # -----------------------------------------------------------------------------   
  
  def check(self, by='both', **kwargs):
    if by == 'python' or by == 'both':
      self.check_by_python(**kwargs)
    if by == 'fullwave' or by == 'both':
      self.check_by_fullwave(**kwargs)

  # ----------------------------------------------------------------------------- 
  
  def check_stability(self, **kwargs):
    from fullwavepy.solver.checks import check_stability
    
    assert self.proj.equation == 'acoustic' # DOES ProjSynAIT MAKE SENSE?
    vp = self.proj.inp.truevp.read()
    
    check_stability(self.proj.dx, 
                    self.proj.dt, 
                    np.max(vp),
                    self.proj.kernel) 
  
  # -----------------------------------------------------------------------------
  
  def check_accuracy(self, **kwargs):
    from fullwavepy.solver.checks import check_accuracy
    from fullwavepy.signal.wavelet import find_passband
    
    # READ SOURCE'S AMP. SPECTRUM
    self.__log.debug('Assuming there is only 1 wavelet, RawSign.')
    rsg = self.rsg.read()
    f_min, f_max = find_passband(rsg, self.proj.dt)
    
    # READ THE MODEL
    vp = self.tvp.read()
    
    check_accuracy(self.proj.dx,
                   np.min(vp), 
                   f_max, 
                   self.proj.kernel)

  # -----------------------------------------------------------------------------
  
  def check_propag_dists(self, **kwargs): 
    from fullwavepy.solver.checks import check_propag_dists
    from fullwavepy.signal.wavelet import find_passband
    
    # READ SOURCE'S AMP. SPECTRUM
    self.__log.debug('Assuming there is only 1 wavelet, RawSign.')
    rsg = self.rsg.read()
    f_min, f_max = find_passband(rsg, self.proj.dt)    
  
    # READ THE MODEL
    vp = self.tvp.read()  
    
    check_propag_dists(self.proj.dims, 
                       self.proj.dx, 
                       self.proj.ttime, 
                       np.min(vp), np.max(vp), f_min, f_max)
    
  # -----------------------------------------------------------------------------
  
  def check_by_python(self, **kwargs):
    """
    """
    self.__log.info('\nPython is checking the input...')
    self.check_stability(**kwargs)
    self.check_accuracy(**kwargs)
    self.check_propag_dists(**kwargs)
    
  # ----------------------------------------------------------------------------- 

  def check_by_fullwave(self, **kwargs):
    """
    Check the project's input by running 
    fullwave with a -checkinput flag.
    
    """    
    exe = self.proj.exe['fullwave_local']
    text = '\n\n'
    text += ('Fullwave is checking the input...')
    text += str('\n(exe = %s)' % exe)
    self.__log.info(text)
    
    path = self.path
    o, e = bash(exe+ ' -checkinput ' + self.proj.name, path=path)
    
    mssg = 'Checked input & dumped canonical runfile'
    err = False
    if mssg in o:
      self.__log.info(mssg)
    else:
      err = True
    
    self.__log.warn('Note, accuracy criterion is not checked by Fullwave.')
    self.__log.warn('Make sure log_lvl<=20 in order to see the full log.')
    self.__log.info(o)
    if len(e) > 0:
      self.__log.info(e)
    if err:
      raise OSError('Error while checking the input by Fullwave')
    
  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


@traced
@logged
class ProjOutput(ProjThroughput):
  """
  As in ProjInput.
  
  """  
  
  # -----------------------------------------------------------------------------
  
  def __init__(self, proj, **kwargs):
    super(ProjOutput, self).__init__(proj, **kwargs) 
    self.path = proj.path + '/out/'
    
  # ----------------------------------------------------------------------------- 

  def init(self, **kwargs):
    """
    Initialize log files for all project types.
    
    Notes
    -----
    This is a necessary work-around 
    because proj.init_input refer to
    proj.inp so they must be called 
    AFTER the constructor, not from it.
    """    
    from fullwavepy.project.qc import JobStats
    from fullwavepy.project.files.logs import (OutLogFile, ErrLogFile, 
                                               JobOutLogFile, JobErrLogFile)
    
    from fullwavepy.project.lists.basic import JobFileList
    from fullwavepy.project.lists.extra import WavefieldFileList
    from fullwavepy.project.files.gridded.wavefields import ForwardWavefieldFile
     
    self.__log.debug('Initializing generic-project output...')
    self.jobout = JobFileList(self.proj, self.path, JobOutLogFile, **kwargs)
    self.joberr = JobFileList(self.proj, self.path, JobErrLogFile, **kwargs)
    self.jo = self.jobout
    self.je = self.joberr
    self.out = JobFileList(self.proj, self.path, OutLogFile, **kwargs)
    self.err = JobFileList(self.proj, self.path, ErrLogFile, **kwargs)
    self.o = self.out 
    self.e = self.err
    self.jobstats = JobStats(self.proj, **kwargs)
    self.jstat = self.jobstats
    self.fw = WavefieldFileList(self.proj, ForwardWavefieldFile, **kwargs) 
      
    self.__log.debug('Initializing project-type-specific output...')
    self.proj.init_output(**kwargs)

  # ----------------------------------------------------------------------------- 

  def log(self, run_id, **kwargs):
    """
    Cat all log files.
    
    """    
    self.joberr.no[run_id].cat(**kwargs)
    self.jobout.no[run_id].cat(**kwargs) 
    self.err.no[run_id].cat(**kwargs)
    self.out.no[run_id].cat(**kwargs)

  # -----------------------------------------------------------------------------     

  def prepare(self, *args, **kwargs):
    """
    Actions to take once the job finished running.
    
    """
    for f in [self.out, self.err, self.jobout, self.joberr]:
      f.prepare(**kwargs)
    
    self.proj.prepare_output(*args, **kwargs)
    
  def prep(self, *args, **kwargs):
    self.prepare(*args, **kwargs)
    
  # -----------------------------------------------------------------------------

  def plot(self, *args, **kwargs):
    self.proj.plot_output(*args, **kwargs)

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------

