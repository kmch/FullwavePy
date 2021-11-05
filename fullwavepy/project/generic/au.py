"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.decor import timer
from fullwavepy.generic.system import bash, exists

from fullwavepy.ioapi.segy import SgyMapp


@traced
@logged
class ProjBaseFiles(object):
  """
  """
  def __init__(self, proj, **kwargs):
    self.proj = proj
    base = kwargs.get('base', {})
    for key, val in base.items():
      setattr(self, key, val)
      # REDUNDANCY
      obj = getattr(self.proj.inp, key)
      setattr(obj, 'base', val) 


# -------------------------------------------------------------------------------


@traced
@logged
class ProjBox(object):
  """
  """
  def __init__(self, proj, **kwargs):
    self.proj = proj
    self.x1, self.x2, self.y1, self.y2, self.z1, self.z2 = proj.box

  # -----------------------------------------------------------------------------
  
  def plot(self, **kwargs):
    from fullwavepy.plot.misc import plot_square
    kwargs['label'] = self.proj.name
    plot_square(self.proj.box[0], self.proj.box[1], 
                self.proj.box[2], self.proj.box[3], **kwargs)     
    
  # -----------------------------------------------------------------------------  
  
  def plotly(self, fig=None, **kwargs):
    import plotly.graph_objects as go
    kwargs = {'line': dict(color=kw('color', 'red', kwargs), width=2, dash=None),
              'showlegend': False, 'name': self.proj.name # WILL APPEAR ON HOVER
             }
    
    X = np.arange(self.x1, self.x2+1)
    Y = np.arange(self.y1, self.y2+1)
    
    if fig is None:
      fig = go.Figure()    
    
    fig.add_trace(go.Scatter(x=X, y=[self.y1]*len(X), **kwargs))
    fig.add_trace(go.Scatter(x=X, y=[self.y2]*len(X), **kwargs))
    fig.add_trace(go.Scatter(x=[self.x1]*len(Y), y=Y, **kwargs))
    fig.add_trace(go.Scatter(x=[self.x2]*len(Y), y=Y, **kwargs))    
    
    return fig

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ProjGeometry(object):
  """
  Class storing information on the geometry 
  of the project i.e.:
  - grid geometry (dimensions, interval)
  - acquisition geometry (sources and receivers positions
    of the experiment as a whole)
  
  """  
  def __init__(self, proj, reciprocity=True, **kwargs):  
    """
    
    """
    self.proj = proj
    
    # FOLD OUT geom IN kwargs IF IT'S PROVIDED AS A SUBDIR OF KWARGS
    if 'geom' in kwargs:
      geom = kwargs['geom']
      assert isinstance(geom, dict)
      del kwargs['geom']
      kwargs = {**kwargs, **geom}
    
    
    for f in [proj.inp.sp, proj.inp.runfile]:
      if exists(f.fname):
        # THIS SHOULD CREATE f.params DICTIONARY
        f.read()
      else:
        self.__log.warning(f.fname + ' not found.')
    
    self._set_discret(**kwargs)
    self._extend_grid(**kwargs)
    self.__log.debug('Double-check _retain_rreceivers')
    #self._retain_receivers(**kwargs)
    
    if self.proj.dims[1] == 1:
      self.proj.dim = '2d'
    
    # PRINT INFO
    box_attrs = ['x1', 'x2', 'y1', 'y2', 'z1', 'z2']
    box_values = self.proj.box
    attrs = ['dx', 'nx1', 'nx2', 'nx3', 'nn', 'ttime', 'dt', 'ns']
    values = [getattr(self.proj, attr) for attr in attrs]
    attrs = box_attrs + attrs
    values = box_values + values
    # WIDTH = 41 = 20 + 20 + 1 (colon)
    self.__log.info('\n' + '{:-^41}'.format('') + '\nProject geometry' + 
                    '\n' + '{:-^41}'.format('') + 
                    ''.join(['\n' + '{:<20}'.format(key) + ':' + '{:>20}'.format(value) 
                    for key, value in zip(attrs, values)]))     

  # -----------------------------------------------------------------------------
  
  def _set_discret(self, **kwargs):
    """
    Set discretization params for both 
    space and time.
    
    Notes
    -----
    Time discretisation can be read only
    from data files and/or segyprep.key.
    Parsing data files is not added yet.

    In SI units (metres, seconds).    
    
    """  
    from functools import reduce
    
    proj = self.proj
    
    sources = [type('Kwargs', (object,), 
                    {'fname' : 'kwargs', 'params': kwargs}),
               proj.inp.sp, proj.inp.runfile]
    
    self.__log.debug('Searching for grid dimensions in ' +
                    ', '.join([i.fname for i in sources]))
    
    # SEARCH FOR PARAMS IN VARIOUS PLACES
    for param in ['dx', 'box', 'dims', 'nx1', 'nx2', 'nx3', 'dt', 'dtms', 'ttime', 'ns']:
      val = None
      for source in sources:
        if val is None: # OVERWRITE ONLY IF NOT FOUND PREVIOUSLY
          val = kw(param, None, source.params)
          setattr(proj, param, val)
    
    # DEAL WITH TIME-DISCR. READ FROM SP.key
    if getattr(proj, 'dt') is None:
      if getattr(proj, 'dtms') is None:
        raise TypeError('Neither dt nor dtms found in any of: ' + 
                      ', '.join([i.fname for i in sources]))
      else:
        setattr(proj, 'dt', float(getattr(proj, 'dtms')))

    if getattr(proj, 'ns') is None:
      dt = getattr(proj, 'dt')
      ttime = getattr(proj, 'ttime')
      if dt is not None and ttime is not None:
        ns = int(float(ttime) / (float(dt) * 1000))
        setattr(proj, 'ns', ns)
    
    # SET THE PRIMARY PARAMS
    for param in ['dx', 'dt', 'ns']:
      if getattr(proj, param) is None:
        raise TypeError(param+str(' not found in any of: ' + 
                      ', '.join([i.fname for i in sources])))
      else:
        setattr(proj, param, float(getattr(proj, param)))
    
    param = 'ns'
    setattr(proj, param, int(getattr(proj, param)))
    
    if proj.box is not None:
      self.__log.debug('Setting proj.dims based on the box provided.')
      proj.dims = self.box2dims(proj.box, proj.dx)
    
    elif proj.dims is not None:
      self.__log.debug('Setting proj.box based on the dims provided.')
      proj.box = self.dims2box(proj.dims, proj.dx)
    
    elif ((proj.nx1 is not None) and (proj.nx2 is not None) and 
          (proj.nx3 is not None)):
      self.__log.debug('Setting proj.box and proj.dims based on' + 
                       ' the nx1, nx2, nx3 provided.')
      proj.dims = [int(float(i)) for i in [proj.nx1, proj.nx2, proj.nx3]]
      proj.box = self.dims2box(proj.dims, proj.dx)
    
    else:
      raise IOError('Grid dimensions not specified correctly. Provide' +
                    ' a box, dims or all nx1, nx2, nx3 in any of: ' + 
                    ', '.join([i.fname for i in sources]))
      
    proj.nx1, proj.nx2, proj.nx3 = proj.dims
    proj.nn = reduce((lambda x, y: x * y), proj.dims) # TOTAL NO. OF NODES
    proj.ttime = proj.ns * proj.dt
    
    if proj.nx2 == 1:
      proj.dim = '2d' # FIXME HERE?
    
  # ----------------------------------------------------------------------------- 
 
  def _extend_grid(self, **kwargs):
    """
    (by e nodes, see the runfile).
    
    Parameters
    ----------
    Returns
    -------
    etop : float 
      No. of e nodes added to the top edge.
    eleft  : float 
      No. of e nodes added to the left edge.
    efront  : float 
      No. of e nodes added to the front edge.
      
    Notes
    -----
    It must be consistent with what Fullwave3D does!
    It should also be synced with io_mod.f90/Set_Aux_Variables.
    
    """
    runfile = self.proj.inp.runfile
    
    if not exists(runfile.fname):
      self.__log.warning(runfile.fname + ' not found. Cannot read extra nodes.')
      return
    
    err = False
    params = ['btop', 'bbot', 'bleft', 'bright',
              'etop', 'ebot', 'elef', 'erig',]
    
    if self.proj.dim == '3d':
      params += ['bfront', 'bback', 
                 'efro', 'ebac']
    
    for param in params:    
      try:
        val = int(runfile.params[param])
        setattr(self.proj, param, val)
      except KeyError:
        err = True
        self.__log.warning(param + ' missing from ' + runfile.fname)
        
    if err:    
      self.__log.warning('Failed to extend_grid. Some of the boundaries info missing.')
    else:
      if (self.proj.btop == 0) and (self.proj.etop < 2):
        self.proj.etop += 2
      if (self.proj.bbot == 0) and (self.proj.ebot < 2):
        self.proj.ebot += 2
      if (self.proj.bleft == 0) and (self.proj.elef < 2):
        self.proj.elef += 2
      if (self.proj.bright == 0) and (self.proj.erig < 2):
        self.proj.erig += 2
      
      self.proj.enx1 = self.proj.nx1 + self.proj.elef + self.proj.erig
      self.proj.enx3 = self.proj.nx3 + self.proj.etop + self.proj.ebot
      
      if self.proj.dim == '3d':
        if (self.proj.bfront == 0) and (self.proj.efro < 2):
          self.proj.efro += 2
        if (self.proj.bback == 0) and (self.proj.ebac < 2):
          self.proj.ebac += 2    
      
        self.proj.enx2 = self.proj.nx2 + self.proj.efro + self.proj.ebac
      else:
        self.proj.enx2 = self.proj.nx2
        self.proj.efro = 0
        self.proj.ebac = 0
      
      self.proj.edims = (self.proj.enx1, self.proj.enx2, self.proj.enx3)
    
    
    
    
    if self.proj.dim == '2d': # FIXME AD-HOC
      self.proj.enx2 = 1
    

    #if verbos > 1:
      #print(this_func, 'etop, eleft, efront', etop, eleft, efront)  
      #print(this_func, 'ebot, eright, eback', ebot, eright, eback)  
 
  # -----------------------------------------------------------------------------
  
  def dims2box(self, dims, dx, **kwargs):
    """
    
    """
    nx1, nx2, nx3 = [int(i) for i in dims]
    x1 = 0
    x2 = (nx1 - 1) * dx
    y1 = 0
    y2 = (nx2 - 1) * dx
    z1 = 0
    z2 = (nx3 - 1) * dx
    box = [x1, x2, y1, y2, z1, z2]    
    return box

  # -----------------------------------------------------------------------------
  
  def box2dims(self, box, dx, **kwargs):
    """
    Add:
    Should be checking if it's really int.
    """
    x1, x2, y1, y2, z1, z2 = box
    assert x2 >= x1
    assert y2 >= y1
    assert z2 >= z1
    nx1 = int((x2 - x1) / dx) + 1 
    nx2 = int((y2 - y1) / dx) + 1  
    nx3 = int((z2 - z1) / dx) + 1     
    dims = (nx1, nx2, nx3)
    return dims

  # -----------------------------------------------------------------------------
    
  def _retain_receivers(self, receivers_csv=None, **kwargs): # FIXME:DEL
    """
    Find the receivers that are contained 
    in the box.
    
    Notes
    -----
    It makes RawSeis.txt contain as few data 
    files as possible. This is essential 
    for SegyPrep runs not to be 
    ridiculuously slow (as they are 
    for ~100 huge receiver files)
    
    It takes advantage of pre-prepared 
    lists of all sources and receivers 
    (within the box used by Ben as of 15.05.2019)
    
    """
    from pandas import read_csv
    from fullwavepy.ndat.arrays import list2str
    
    if receivers_csv is None:
      self.__log.warning('No receivers_csv file provided. Returning none.')
      return
    
    df = read_csv(receivers_csv)
    
    retained = []
    
    for i in df.index: 
      x = float(df['x'][i])
      y = float(df['y'][i])
      
      if (x > self.proj.box[0]) and (x < self.proj.box[1]):
        if (y > self.proj.box[2]) and (y < self.proj.box[3]):
          retained.append(df['id'][i])
    
    n = len(retained)
    self.__log.info('List of ' + str(n) + ' retained receivers: ' + 
                    list2str(retained))    
    
    if n == 0:
      raise ValueError('No receivers contained in the box!')
    
    self.proj.r_retained = retained
    
  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------


@traced
@logged
class ProjPath(object):
  def __init__(self, proj, **kwargs):
    """
      path : str
        A path where proj_name directory resides 
        (or will reside). 
        Default: ./ 
    """
    path = kw('path', './', kwargs)
    proj.path = path + '/{pname}/'.format(pname=proj.name)
    self.__log.debug('Project path set to: ' + proj.path)


# -------------------------------------------------------------------------------


@traced
@logged
class ProjDirs(object):
  def __init__(self, proj, **kwargs):
    """
    Make project's directory
    and subdirectories.
    
    """
    from fullwavepy.generic.system import exists
    
    pro_dir = proj.path 
    inp_dir = proj.path + '/inp/'
    out_dir = proj.path + '/out/'
    
    for d in [pro_dir, inp_dir, out_dir]:
      if not exists(d, **kwargs):
        o, e = bash('mkdir ' + d)


# -------------------------------------------------------------------------------


@traced
@logged
class ProjDef(object):
  def __init__(self, proj, **kwargs):
    """
    Set essential parameters defining 
    a problem-type to run by Fullwave.
    
    Notes
    -----
    The only way these params are already set 
    as project-attributes is via 
    proj.__init__:runfile.set_proj_params()
    
    qp, qs are not runfile params, they are used 
    by FullwavePy only. Attenuation is 
    switched on in Fullwave automatically if 
    relevant model files are found.
    
    In anisotropy, 'none' is a string directly
    writeable to a runfile.
    
    """
    
    # LIST OF PROBLEM-DEFINING proj-ATTRIBUTES 
    attrs = ['problem', 'domain', 'dim', 'equation', 'anisotropy', 'kernel',
             'io', 'units',
             'qp', 'qs', # NON-RUNFILE PARAMS, YET PROBLEM-DEFINING
            ]
    
    # DEFAULT VALUES OF attrs:
    defaults = ['synthetic', 'time', '3d', 'acoustic', 'none', 'low',
                'sgy', 'metric',
                False, False]
    
    for attr in attrs: 
      if attr in kwargs:
        setattr(proj, attr, kwargs[attr])
      else:
        i = attrs.index(attr)
        setattr(proj, attr, defaults[i])
    
    values = [getattr(proj, attr) for attr in attrs]
    
    # WIDTH CALCULATION 20 + 20 + 1 (colon) = 41
    self.__log.info('\n' + '{:-^41}'.format('') + '\nProject definition' + 
                    '\n' + '{:-^41}'.format('') + 
                    ''.join(['\n' + '{:<20}'.format(key) + ':' + '{:>20}'.format(value) 
                    for key, value in zip(attrs, values)]))
 
  
# -------------------------------------------------------------------------------


@traced
@logged
class ProjExe(object):
  def __init__(self, proj, **kwargs):
    """
      paths : dict
        A dictionary of paths (str) to executables 
        used by different functions. The functions
        will raise an error if any path they need 
        is not defined. 
        Default: {} (empty).
        Example:
        paths = {'segyprep': '~/segyprep/bin/segyprep',
                 'fullwave': '~/rev690/bin/fullwave3D.exe'}
    
    """  
    proj.exe = kw('exe', {}, kwargs)


# -------------------------------------------------------------------------------  


@traced
@logged
class ProjEnv(object):
  """
  Default values of environment variables.
    
  """
  def __init__(self, proj, dump_all=False, **kwargs): # TO FINISH
    """
    Set default values (mostly Fullwave's defaults,
    soome of them are overwritten in projtypes, the 
    rest can be overwritten anywhere).
    This serves as a lisiting of all env vars.
    
    Notes
    -----
    All defaults may as well be set up in projtypes.
    
    Up-to-date with rev690 version. Run fullwave with 
    env flag for more info.

    Fixme: update for 700-NewCodebase version after 
    consortium meeting.

    """
    max_ram_gb = kw('max_ram_gb', 40, kwargs) # HERE?
    
    self.proj = proj
    
    env = kw('env', {}, kwargs)
    self.var = {}
    self.var['SCHEDULER_SHOWLEVEL'] = kw('SCHEDULER_SHOWLEVEL', 2, env)
    self.var['SCHEDULER_CHECKVERBOSEFILE'] = kw('SCHEDULER_CHECKVERBOSEFILE', None, env)
    self.var['SCHEDULER_NTHREADS'] = kw('SCHEDULER_NTHREADS', None, env)
    self.var['SCHEDULER_CONVERGENCE'] = kw('SCHEDULER_CONVERGENCE', None, env)
    self.var['SCHEDULER_WRITEUPDATE'] = kw('SCHEDULER_WRITEUPDATE', None, env)
    self.var['SCHEDULER_GLOBPRECSTAB'] = kw('SCHEDULER_GLOBPRECSTAB', None, env)
    self.var['SCHEDULER_SPREADSRCS'] = kw('SCHEDULER_SPREADSRCS', None, env)
    self.var['SCHEDULER_SPREADRCVRS'] = kw('SCHEDULER_SPREADRCVRS', None, env)
    self.var['SCHEDULER_DUMPMODELPROPS'] = kw('SCHEDULER_DUMPMODELPROPS', None, env)
    self.var['SCHEDULER_DUMPGRAD'] = kw('SCHEDULER_DUMPGRAD', None, env)
    self.var['SCHEDULER_DUMPPREC'] = kw('SCHEDULER_DUMPPREC', None, env)
    self.var['SCHEDULER_DUMPRAWGRAD'] = kw('SCHEDULER_DUMPRAWGRAD', None, env) 
    self.var['SCHEDULER_DUMPRAWPREC'] = kw('SCHEDULER_DUMPRAWPREC', None, env)
    self.var['SCHEDULER_DUMPTTIPARAMS'] = kw('SCHEDULER_DUMPTTIPARAMS', None, env)
    self.var['SCHEDULER_SLAVESHAVEDATA'] = kw('SCHEDULER_SLAVESHAVEDATA', 'yes', env)
    self.var['SCHEDULER_SLAVESHAVEMODEL'] = kw('SCHEDULER_SLAVESHAVEMODEL', 'yes', env)
    self.var['SCHEDULER_SLAVESHAVEUPDATE'] = kw('SCHEDULER_SLAVESHAVEUPDATE', None, env)
    self.var['SCHEDULER_SLAVESTIMEOUT'] = kw('SCHEDULER_SLAVESTIMEOUT', None, env)
    self.var['SCHEDULER_SLAVESMAXRECOVER'] = kw('SCHEDULER_SLAVESMAXRECOVER', None, env)
    self.var['SCHEDULER_ERRSTREAM'] = kw('SCHEDULER_ERRSTREAM', 'stderr', env)
    self.var['SCHEDULER_SHOWTIMESTAMP'] = kw('SCHEDULER_SHOWTIMESTAMP', 'yes', env)
    self.var['SCHEDULER_SHOWDATESTAMP'] = kw('SCHEDULER_SHOWDATESTAMP', None, env)
    self.var['SCHEDULER_SHOWHOSTNAME'] = kw('SCHEDULER_SHOWHOSTNAME', None, env)
    self.var['SCHEDULER_FORCESLAVESGETMODEL'] = kw('SCHEDULER_FORCESLAVESGETMODEL', None, env)
    self.var['SCHEDULER_FORCESLAVESGETDATA'] = kw('SCHEDULER_FORCESLAVESGETDATA', None, env)

    self.var['SLAVES_SHOWLEVEL'] = kw('SLAVES_SHOWLEVEL', 2, env)
    self.var['SLAVES_CHECKVERBOSEFILE'] = kw('SLAVES_CHECKVERBOSEFILE', None, env)
    self.var['SLAVES_TIMEOUT'] = kw('SLAVES_TIMEOUT', None, env)
    self.var['SLAVES_SOURCEMULT'] = kw('SLAVES_SOURCEMULT', None, env)
    self.var['SLAVES_FORCESOURCESMOOTH'] = kw('SLAVES_FORCESOURCESMOOTH', None, env)
    self.var['SLAVES_WAVEFIELDSVTR'] = kw('SLAVES_WAVEFIELDSVTR', None, env)
    self.var['SLAVES_DUMPGRAD'] = kw('SLAVES_DUMPGRAD', None, env)
    self.var['SLAVES_DUMPPREC'] = kw('SLAVES_DUMPPREC', None, env)
    self.var['SLAVES_DUMPDAT'] = kw('SLAVES_DUMPDAT', None, env)
    self.var['SLAVES_DUMPLOWPASS'] = kw('SLAVES_DUMPLOWPASS', None, env) #proj.name+'-SLAVES_DUMPLOWPASS', env)
    self.var['SLAVES_DUMPAGC'] = kw('SLAVES_DUMPAGC', None, env) #proj.name+'-SLAVES_DUMPAGC', env)
    self.var['SLAVES_DUMPCOMPARE'] = kw('SLAVES_DUMPCOMPARE', proj.name+'-SLAVES_DUMPCOMPARE', env)
    self.var['SLAVES_DUMPRESIDS'] = kw('SLAVES_DUMPRESIDS', None, env)
    self.var['SLAVES_DUMPADJOINT'] = kw('SLAVES_DUMPADJOINT', None, env)
    self.var['SLAVES_DUMPDENTAB'] = kw('SLAVES_DUMPDENTAB', None, env)
    self.var['SLAVES_DUMPDENSITY'] = kw('SLAVES_DUMPDENSITY', None, env)
    self.var['SLAVES_DUMPMODELPROPS'] = kw('SLAVES_DUMPMODELPROPS', None, env)
    self.var['SLAVES_DUMPCSREFS'] = kw('SLAVES_DUMPCSREFS', None, env) # NOTE
    #   Restrict dumps to those CSRefs in the list, which is a
    #   comma-separated set of values/ranges, e.g. "5,15-20,35-45,60"    
    self.var['SLAVES_GARDNERPOWER'] = kw('SLAVES_GARDNERPOWER', None, env)
    self.var['SLAVES_GARDNERFACTOR'] = kw('SLAVES_GARDNERFACTOR', None, env)
    self.var['SLAVES_SRCPRECSTAB'] = kw('SLAVES_SRCPRECSTAB', None, env)
    self.var['SLAVES_RCVRPRECSTAB'] = kw('SLAVES_RCVRPRECSTAB', None, env)
    self.var['SLAVES_LOCALSTORE'] = kw('SLAVES_LOCALSTORE', '$work_dir', env)
    self.var['SLAVES_LOCALSTORETRIGGER'] = kw('SLAVES_LOCALSTORETRIGGER', str(max_ram_gb), env)
    self.var['SLAVES_ERRSTREAM'] = kw('SLAVES_ERRSTREAM', 'stderr', env)
    self.var['SLAVES_FORCEGRIDDECI'] = kw('SLAVES_FORCEGRIDDECI', None, env)
    self.var['SLAVES_FORCETIMEDECI'] = kw('SLAVES_FORCETIMEDECI', None, env)
    self.var['SLAVES_SHOWTIMESTAMP'] = kw('SLAVES_SHOWTIMESTAMP', 'yes', env)
    self.var['SLAVES_SHOWDATESTAMP'] = kw('SLAVES_SHOWDATESTAMP', None, env)
    self.var['SLAVES_SHOWHOSTNAME'] = kw('SLAVES_SHOWHOSTNAME', None, env)
    
    if dump_all:
      self.set_all_dumps_true(**kwargs)

    self._true2prefix(**kwargs)
    
  # -----------------------------------------------------------------------------
  
  def csrefs_to_dump(self, **kwargs): #FIXME IT IS IMPLEMENTED SOMEWHERE!!!
    """
    Get the list of IDs of composite-source  to dump.
    
    Notes
    -----
    It is used by external functions.
    
    """
    csrefs = self.var['SLAVES_DUMPCSREFS']

    if csrefs is None:
      pass
    else:
      if '-' in csrefs:
        raise NotImplementedError('Ranges not implemented. Provide a list of values.')
      csrefs = [int(i) for i in csrefs.split(',')]
    
    return csrefs
  
  # -----------------------------------------------------------------------------
    
  def decimate_shot_dumps(self, factor, **kwargs):
    """
    Dump info for every {factor}th shot.

    factor : int 
      Take every {factor}th shot.
    """
    pass
  
  # -----------------------------------------------------------------------------

  def set_all_dumps_true(self, **kwargs):
    """
    """
    for key, val in self.var.items():
      self.var[key] = val
      if ('DUMP' in key) and (key != 'SLAVES_DUMPCSREFS'):
          self.var[key] = True

  # -----------------------------------------------------------------------------

  def _true2prefix(self, **kwargs):
    """
    For SLAVES_DUMP* env vars convert 'yes', True, etc.
    into a required prefix.
    
    Notes
    -----
    The prefix can be arbitrary. Here:
    'proj_name-env_var_name' convention is used.
    
    """
    for key, val in self.var.items():
      if 'SLAVES_DUMP' not in key:
        continue
      
      if key == 'SLAVES_DUMPCSREFS': # THIS HAS DIFFERENT FORMAT
        continue
      
      convert = False
      if isinstance(val, str):
        s = val.lower()
        if s == 'yes' or s == 'true' or s == '1':
          convert = True
      elif val:
        convert = True
      
      if convert:
        self.var[key] = self.proj.name + '-' + key #NOTE

  # -----------------------------------------------------------------------------

  def info(self, **kwargs):
    """
    Print info about env vars by running fullwave3d executable
    with -env flag.
    
    """
    command = self.proj.exe['fullwave_local'] + ' -env ' + self.proj.name
    o, e = bash(command, path=self.proj.path)
    self.__log.info('Info output by a command: ' + command + '\n' + o + e)

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ProjSgyMapp(SgyMapp):
  """
  Mapping between physical quantities 
  and SEG-Y header-words.
  
  Notes
  -----
  The default values are what Fullwave3D
  actually produces.
  
  Values for the PROTEUS experiment 
  using reciprocity:
   sid = 'tracf'
   rid = 'fldr'
   lid = 'ep'
  
  """
  def __init__(self, proj, **kwargs):
    """
    """
    sgyhw = kw('sgyhw', {}, kwargs)
    if len(sgyhw) == 0:
      self.__log.info("Setting SEG-Y mapping (sgyhw) to Fullwave3D's default.")
    
    self['sid'] = kw('sid', 'fldr', sgyhw)
    self['rid'] = kw('rid', 'tracf', sgyhw)
    self['lid'] = kw('lid', 'ep', sgyhw)
    self['xmod'] = kw('xmod', 'sx', sgyhw)
    self['ymod'] = kw('ymod', 'sy', sgyhw)
    
  
# -------------------------------------------------------------------------------


@traced
@logged
class ProjCluster(object):
  """
  Computer cluster for which proj.inp.pbs files
  will be prepared.

  """
  def __init__(self, proj, **kwargs):
    """
    """
    self.proj = proj
    self.name = kw('cluster', 'cx1', kwargs)
    assert isinstance(self.name, str)
    self.__log.info('Setting cluster to:  %s.' % self.name)


# -------------------------------------------------------------------------------

