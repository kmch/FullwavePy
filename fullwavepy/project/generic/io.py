"""
This module defines generic input and output of an FWI project 
(shared by all project-types). 

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import widgets
from fullwavepy.generic.system import bash, current_dir
from fullwavepy.generic.parse import kw

from fullwavepy.ioapi.clusters.cx1 import *
from fullwavepy.ioapi.clusters.archer import *
from fullwavepy.ioapi.clusters.thomas import *



@traced
@logged
class FileManager(object):
  """
  """
  def ls(self, **kwargs):
    o, e = bash('pwd', path=self.path)
    print('Content of ' + o)
    o, e = bash('ls -lth ' + self.path)
    print(o, e)

  # -----------------------------------------------------------------------------

  def rm(self, **kwargs):
    """
    Remove all the files from the self.path (!) NOTE
    
    Notes
    -----
    It actually renames files (by adding a 'rm_' prefix) 
    which need to be rm-ed manually (for safety).
    
    """
    from fullwavepy.generic.system import get_files
    from fullwavepy.generic.parse import path_leave
    
    ls = kwargs.get('ls', True)

    fnames = get_files(self.path, '*')
    for fname in fnames:
      fname = path_leave(fname)
      o, e = bash('mv ' + fname + ' rm_' + fname, path=self.path)
    
    if ls: 
      self.ls()

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


# class _PreprocSynInp():
#   attrs = ['rawsigns', 'rawseis', 'truemods', 'params']  

# class _PreprocSynOut():
#   attrs = ['wavelets', 'rawseis', 'truemods', 'params']  

# ENSURE PreprocOut = ProjInp

# class _ProjSynInp():
#   attrs = ['wavelets', 'rawseis', 'truemods', 'params']


# -------------------------------------------------------------------------------


@traced
@logged
class _InpPreparer(object):
  """
  Mix-in providing input-preparing functionality.

  """
  def prepare_input(self, other_proj=None, run=False, **kwargs):
    """
    Interface.
    
    Notes
    -----
    If we want to prepare some files from scratch and some 
    from another project, we need to do it manually for all
    of them or just overwrite selected ones. We need to 

    To pass run=True, we need to pass first arg too.

    """
    rm = kwargs.get('rm', True)
    if rm:
      self.inp.rm(**kwargs)
    
    if other_proj is None:
      self._prep_inp_from_scratch(**kwargs)
    else:
      self._prep_inp_from_another(other_proj, **kwargs)
    
    self._preprocess(**kwargs)
    self._postprocess(**kwargs)

  # -----------------------------------------------------------------------------

  def _preprocess(self, *args, **kwargs):
    self.inp.sp.run(*args, **kwargs)

  # -----------------------------------------------------------------------------

  def _postprocess(self, *args, **kwargs):
    rnf_kwargs = kwargs.get('rnf_kwargs', {})
    self.inp.rnf.prep(*args, **rnf_kwargs)

  # -----------------------------------------------------------------------------

  def _prep_inp_from_scratch(self, **kwargs):
    skip = kwargs.get('skip', [])
    obj_ids = [i for i in ['rsg', 'tvp', 'rse', 'sp'] if not i in skip]
    objects = [getattr(self.inp, obj_id) for obj_id in obj_ids]
    
    for obj in objects:
      self.__log.info('Preparing %s...\n' % obj.name)
      obj.prep(obj.base, **kwargs)

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


@traced
@logged
class ProjThroughput(FileManager):
  """
  Parent class of input and output.
  
  Notes
  -----
  Not much implemented - it's left 
  for subclasses to define their own
  versions of methods.
  
  """  
  def __init__(self, proj, **kwargs):
    self.proj = proj
    
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
        self.__log.warning('Skipping PBS script ' + fname + 
                        ' prepare it yourself.')
        continue

      if 'Job' in fname or (('Out' in fname or 'Err' in fname) and exten(fname) == 'log'):
        self.__log.warning('Skipping ' + fname)
        continue
    
      nfname = self.path + self.proj.name + nfname[len(project_name): ]
      duplicate(fname, nfname, **kwargs)

  # -----------------------------------------------------------------------------    

  def _configure_remote(self, host_nick, **kwargs):
    if host_nick == 'cx1':
      from fullwavepy.ioapi.clusters.cx1 import proj_path_cx1
      host = "kmc3817@login.cx1.hpc.ic.ac.uk"
      remote_path = host + ':' + proj_path_cx1
    elif host_nick == 'eph':
      from fullwavepy.ioapi.clusters.cx1 import proj_path_eph
      host = "kmc3817@login.cx1.hpc.ic.ac.uk"
      remote_path = host + ':' + proj_path_eph  
    elif host_nick == 'my_eph':
      remote_path = '/home/kmc3817/rds_home/my_ephemeral/PROJECTS/'
 
    else:
      raise NotImplementedError('host nick: %s' % host_nick)
    
    self.remote_path = remote_path + self.proj.parent_dir + self.proj.name

  # -----------------------------------------------------------------------------

  def rsync(self, thr, host_nick, only_logs=False, **kwargs):
    """
    thr : str 
      'inp' or 'out'
    host_nick : str
      Probably 'my_eph'

    Notes
    -----
    We could actually *sync* local and remote dirs so that 
    they are the same all the time (what Dropbox does) but 
    instead we push inp and pull out.

    """
    self._configure_remote(host_nick, **kwargs)
    local_path = self.proj.path
    remot_path = self.remote_path
    
    if thr == 'inp':
      self.__log.debug('Pushing inp from local to remote.')
      source = local_path + '/inp/'
      destin = remot_path + '/inp' # no slash
    elif thr == 'out':
      self.__log.debug('Pulling out from remote to local.')
      source = remot_path + '/out/' # NOTE: remote is now source!
      if only_logs:
        source += '*.log'
        self.__log.info('Syncing only %s' % source)
      destin = local_path + '/out'  # no slash either
    else:
      raise ValueError('thr %s' %s)
    
    cmd = 'rsync -azP --info=progress2 %s %s' % (source, destin)
    self.__log.debug('cmd %s' % str(cmd))
    
    o, e = bash(cmd)
    self.__log.info(str(o))
    if len(e) < 0:
      self.__log.warning(str(e))


  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ProjInput(ProjThroughput):
  """
  Project's input.
  Common for all Fullwave3D project types.
  FIXME: Should be generic (NOT ONLY FULLWAVE)

  Notes
  -----
  It contains some necessary work-arounds
  because proj.init_input etc. refer to
  proj.inp so they must be called 
  AFTER the constructor (__init__), not from it.  
  
  """  
  def __init__(self, proj, **kwargs):
    super(ProjInput, self).__init__(proj, **kwargs)
    self.path = proj.path + '/inp/'
  
  # ----------------------------------------------------------------------------- 
  
  def init(self, **kwargs):
    """
    Defines a couple of aliases to save typing-time.
    
    """
    from fullwavepy.project.files.datalike.sgy import RawSignFile, SignatureFileSgy
    from fullwavepy.project.files.datalike.ttr import SignatureFileTtr
    from fullwavepy.project.files.gridded.misc import InextFile
    from fullwavepy.project.files.gridded.surfaces import TopoFile, FsFile, ExtendedFsFile, InterpolFsFile
    from fullwavepy.project.files.other.ghost import GhostDataFileBin, GhostDataFileTxt, \
      SourcesDataFileTxt, ReceiversDataFileTxt
    from fullwavepy.project.files.text.srcrec import SourcesFile, ReceiversFile
    from fullwavepy.project.files.text.misc import RawSeisTxtFile, JobInfoFile
    from fullwavepy.project.files.text.runfiles import SegyPrepFile, Runfile, Skeleton
    from fullwavepy.project.files.text.submit import BashFile
    from fullwavepy.project.lists.basic import JobFileList 
    
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
    
    
    topo = kw('topo', None, kwargs)
    self.topo = TopoFile(self.proj, self.path, dupl=topo, **kwargs)  
    self.fs = FsFile(self.proj, self.path, **kwargs)
    self.fse = ExtendedFsFile(self.proj, self.path, **kwargs)
    self.fsi = InterpolFsFile(self.proj, self.path, **kwargs)
    self.ghostbin = GhostDataFileBin(self.proj, self.path, **kwargs)
    self.ghb = self.ghostbin
    self.ghosttxt = GhostDataFileTxt(self.proj, self.path, **kwargs)
    self.ght = self.ghosttxt
    self.ine = InextFile(self.proj, self.path, **kwargs)
    self.sdata = SourcesDataFileTxt(self.proj, self.path, **kwargs)
    self.rdata = ReceiversDataFileTxt(self.proj, self.path, **kwargs)
    
    self.s = SourcesFile(self.proj, self.path, **kwargs)
    self.r = ReceiversFile(self.proj, self.path, **kwargs)

    self.skeleton = Skeleton(self.proj, self.path, **kwargs)
    self.ske = self.skeleton
    self.runfile = Runfile(self.proj, self.path, **kwargs)
    self.rnf = self.runfile # ALIAS
    self.bash = JobFileList(self.proj, self.path, BashFile, **kwargs)
    
    if self.proj.cluster.name == 'cx1':
      self.pbs = JobFileList(self.proj, self.path, PbsFileCx1, **kwargs)
    elif self.proj.cluster.name == 'thomas':
      self.pbs = JobFileList(self.proj, self.path, PbsFileThomas, **kwargs)    
    elif self.proj.cluster.name == 'archer':
      self.pbs = JobFileList(self.proj, self.path, PbsFileArcher, **kwargs)    
    else:
      raise NotImplementedError('Unknown cluster: ' + self.proj.cluster.name)    
    
    self.jobinfo = JobFileList(self.proj, self.path, JobInfoFile, **kwargs)
    self.jinfo = self.jobinfo # ALIAS
  
  # ----------------------------------------------------------------------------- 

  def rsync(self, *args, **kwargs):
    super().rsync('inp', *args, **kwargs)

  # -----------------------------------------------------------------------------
  
  def prep(self, *args, **kwargs):
    self.proj.prepare_input(*args, **kwargs)

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
    from fullwavepy.fd.checks import check_stability
    
    assert self.proj.equation == 'acoustic' # DOES ProjSynAIT MAKE SENSE?
    vp = self.proj.inp.truevp.read()
    
    check_stability(self.proj.dx, 
                    self.proj.dt, 
                    np.max(vp),
                    self.proj.kernel) 
  
  # -----------------------------------------------------------------------------
  
  def check_accuracy(self, **kwargs):
    from fullwavepy.fd.checks import check_accuracy
    from fullwavepy.dsp.wavelet import find_passband
    
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
    from fullwavepy.fd.checks import check_propag_dists
    from fullwavepy.dsp.wavelet import find_passband
    
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
    
    self.__log.warning('Note, accuracy criterion is not checked by Fullwave.')
    self.__log.warning('Make sure log_lvl<=20 in order to see the full log.')
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
  Project's output.

  Notes
  -----  
  As in ProjInput.
  
  """  
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
    from fullwavepy.project.generic.qc import JobStats
    from fullwavepy.project.files.text.logs import (OutLogFile, ErrLogFile, 
                                               JobOutLogFile, JobErrLogFile)
    
    from fullwavepy.project.lists.basic import JobFileList 

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

      
    self.__log.debug('Initializing project-type-specific output...')
    self.proj.init_output(**kwargs)

  # ----------------------------------------------------------------------------- 

  def rsync(self, *args, **kwargs):
    super().rsync('out', *args, **kwargs)

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

  def prep(self, *args, **kwargs):
    """
    Actions to take once the job finished running.
    
    """
    for f in [self.out, self.err, self.jobout, self.joberr]:
      f.prepare(**kwargs)
    self.proj.prepare_output(*args, **kwargs)

  # -----------------------------------------------------------------------------

  def plot(self, *args, **kwargs):
    self.proj.plot_output(*args, **kwargs)

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------

