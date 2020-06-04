"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.ioapi.generic import CsvFile
from fullwavepy.project.files.generic import TextProjFile



@traced
@logged
class MetaDataProjFile(CsvFile, TextProjFile): # IS IT USED AT ALL?
  """
  
  csv is MUCH faster for Panda's read/write than json.
  Not to mention, we can choose columns.
  
  """
  def __init__(self, proj, path, **kwargs):
    self.suffix = 'MetaData'
    self.ext = 'csv'
    self.name = '{}-{}.{}'.format(proj.name, self.suffix, self.ext)
    self.fname = path + self.name

  # -----------------------------------------------------------------------------

  def read(self, overwrite=True, **kwargs):
    """
  
    Notes
    -----
    Overwrite=True by default because otherwise plots are not 
    updated even though they are supposed (e.g. you are passing 
    a different fname). They will be correct (updated) only
    if you delete self.array variable, e.g. by restarting the 
    notebook kernel.
    Disable overwrite only for PERFORMANCE (e.g. interactive plot)
    when the array remains unchanged unlike other (e.g. plotting)
    parameters.
    
    """
    if (not hasattr(self, 'df')) or overwrite:
      self.__log.warning('{}.df does not exist and will be read.'.format(type(self)))
      self.df = super().read(**kwargs)
    return self.df  

  # -----------------------------------------------------------------------------
  
  def plotly(self, fig=None, df=None, **kwargs):
    tracf0 = 4104
    fldr0 = 9882
    stride = 10
    import plotly.graph_objects as go
    if fig is None:
      fig = go.Figure() 
      
    if df is not None:
      import plotly.graph_objects as go
      data = df[df.tracf==tracf0][::stride]
      fig.add_trace(go.Scatter(x=data['sx'], y=data['sy']))
      data = df[df.fldr==fldr0]
      fig.add_trace(go.Scatter(x=data['gx'], y=data['gy'], mode='markers'))
      #fig.show()
  
    return fig
  
  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------


@traced
@logged
class JobFile(object):
  """
  File bound to a specific job 
  (code-run).
  
  """
  def __init__(self, proj, path, suffix, exten, run_id, **kwargs):
    """
    self.name = lambda run_id : proj.name + '-Out{run_id}.log'.format(run_id=run_id)
    
    """
    self.proj = proj
    self.path = path
    self.run_id = run_id
    self.suffix = suffix
    self.exten = exten
    if self.run_id is None:
      self.name = proj.name + '-' + suffix + '.' + self.exten
    else:
      self.name = proj.name + '-' + suffix + str(run_id) + '.' + self.exten
    self.fname = path + self.name
  
  # -----------------------------------------------------------------------------

  def _create_verbosity_triggers(self, **kwargs):
    """
    Create empty files for each MPI process. 
    Their presence makes Fullwave emit verbose 
    output if relevant env vars are set to yes.
    
    """
    fname = self.proj.inp.path + 'fullwave3d-verbose-scheduler'
    self.__log.debug('Creating ' + fname)
    with open(fname, 'w'):
      pass
    
    self.__log.debug('Creating fullwave3d-verbose-slave-? for each mpiproc')
    for i in range(1, self.mpiprocs):
      fname = self.proj.inp.path + 'fullwave3d-verbose-slave-' + str(i)
      self.__log.debug('Creating ' + fname)
      with open(fname, 'w'):
        pass    

  # -----------------------------------------------------------------------------  
  
  
# -------------------------------------------------------------------------------


@traced
@logged
class LastCheckpointFile(TextProjFile):
  """
  Stores a single number.
  
  """
  # ----------------------------------------------------------------------------- 
  
  def __init__(self, proj, path, **kwargs):
    self.name = proj.name + '-LastCheckpoint.txt'
    self.fname = path + self.name
    try:
      proj.lastcp = self.read()
      self.__log.info('Last checkpoint: ' + str(proj.lastcp))
    except FileNotFoundError:
      self.__log.warning(self.fname + ' not found. Setting proj.lastcp to 1')
      proj.lastcp = 1
    
    super().__init__(proj, path, **kwargs)

  # -----------------------------------------------------------------------------   
  
  def read(self, **kwargs):
    try:
      cp = super().read(**kwargs)
      cp = int(cp[0][0])
    except FileNotFoundError:
      self.__log.warning(self.fname + ' not found. Setting last checkpoint to 0')
      cp = 0
    return cp

  # ----------------------------------------------------------------------------- 
  

# -------------------------------------------------------------------------------


@traced
@logged
class InfoFile(TextProjFile):
  """
  A txt file with description
  of the PROJECT.
  
  Notes
  -----
  Use self.cat() to read it.
  
  """
  def __init__(self, proj, info=None, **kwargs):
    """
    Init. and save the file, if info is not None.
    
    """
    cat = kwargs.get('cat', True)
    self.name = proj.name + '-Info.txt'
    self.fname = proj.path + self.name
    del_kw('path', kwargs)
    super().__init__(proj, proj.path, **kwargs) 
    self._write(info, **kwargs)
    if cat:
      self.cat()

  # ----------------------------------------------------------------------------- 

  def _write(self, info=None, **kwargs):
    """
    Write a file with info.
    
    """
    if info is not None:
      with open(self.fname, 'w') as f:
        f.write('\n') # better looks
        f.write(info)
        f.write('\n') # better looks

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------


@traced
@logged
class JobInfoFile(JobFile, TextProjFile):
  """
  A txt file with the job description
  output by the qsub command.
  
  """
  
  # ----------------------------------------------------------------------------- 
  
  def __init__(self, proj, path, run_id, **kwargs):
    """
    
    """
    suffix = 'JobInfo'
    exten = 'log'
    super().__init__(proj, path, suffix, exten, run_id, **kwargs)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class RawSeisTxtFile(TextProjFile):
  """
  A txt file with paths to 
  data files used by SegyPrep.
  
  Notes
  -----
  We don't use RawSeis.sgy at the moment.
  
  """ 
  def __init__(self, proj, path, **kwargs): # OK
    """
    
    """
    self.name = proj.name + '-RawSeis.txt'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)
    
  # -----------------------------------------------------------------------------  
  
  def create(self, fnames, **kwargs):
    """
    Notes
    -----
    List only these OBSes that are contained in
    the model - speed gain in SegyPrep!
    
    """
    from fullwavepy.ioapi.generic import save_txt
    super().create()
    fnames = self._select_retained(fnames, **kwargs)
    n = len(fnames)
    self.__log.info('No. of fnames selected: ' + str(n))
    save_txt(self.fname, fnames, **kwargs)
    
  # -----------------------------------------------------------------------------
  
  @timer
  def _select_retained(self, fnames, bad_IDs=[], **kwargs):
    """
    fnames : 
      sgy
      
    bad_IDs : list
      List of IDs (ID is a string contained in, 
      and uniquely defining a station name) of 
      stations with bad data (to exclude).
      CAUTION: '152' appears in all PROTEUS file-names 
      as a part of the cruise ID!
      
    Assumes receiver gathers
    
    """
    from fullwavepy.ioapi.segy import json_header_suffix
    from pandas import read_json
    
    self.__log.debug('Selecting the stations contained in the' +
                    ' model box.')
    
    
    if len(bad_IDs) == 0:
      self.__log.warning('len(bad_IDs)=0, selecting all the data' +
                      ' contained in the box...')
    
    nfnames = []
    for fname in fnames:
      bad = False
      for bid in bad_IDs:
        bid = str(bid) + '_' # OTHERWISE SKIPS ALL
        if bid in fname:
          self.__log.info('Skipping bad-data file: ' + fname)
          bad = True
          break
      if bad:
        continue
      
      json = strip(fname) + json_header_suffix
      try:
        df = read_json(json)
      except ValueError:
        self.__log.warning('File ' + json + ' not found. Appending ' + 
                           fname + ' to the list of retained receivers.')
        nfnames.append(fname)
        continue
      
      x = float(df.gx.to_list()[0])
      y = float(df.gy.to_list()[0])
      if (x > self.proj.box[0]) and (x < self.proj.box[1]):
        if (y > self.proj.box[2]) and (y < self.proj.box[3]):
          nfnames.append(fname)
    
    return nfnames


# -------------------------------------------------------------------------------

