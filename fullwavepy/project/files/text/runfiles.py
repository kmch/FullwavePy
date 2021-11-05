"""
This module contains definitions of runfiles
(parameter files read at run-time by various programs).

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.decor import timer
from fullwavepy.project.files.generic import TextProjFile


@traced
@logged
class ParameterFile(TextProjFile):
  """
  
  """
  
  # -----------------------------------------------------------------------------

  def __init__(self, proj, path, **kwargs):
    super().__init__(proj, path, **kwargs)
    self.params = {} # NOT SURE IT'S A GOOD PRACTICE BUT 
    # ProjGeometry._set_grid USES IT

  # -----------------------------------------------------------------------------  

  def read(self, **kwargs):
    """
    Return dict of params 
    and assign it to a 'params'
    attribute as well.
    
    """
    from fullwavepy.ioapi.generic import read_dict

    d = read_dict(self.fname)
    
    # CONVERT ALL TO LOWER CASE
    d = {k.lower() : v for k, v in d.items()} 
    
    self.params = d
    
    return d 

  # -----------------------------------------------------------------------------
  
  def modify(self, **kwargs):
    """
    Replace some portion of the file
    with kwargs.
    
    """
    from fullwavepy.ioapi.generic import read_txt_raw

    fname = self.fname
    content = read_txt_raw(fname)
    
    if len(kwargs) == 0:
      self.__log.warning('Empty kwargs. Returning.')
      return
    
    f = open(self.fname, 'w')

    found_keys = []
    for line in content:
      nline = line
      line = line.split(None)
      
      if len(line) == 0:
        f.write(nline)
        continue      
      
      # FIRST CHECK IF THEY'RE PRESENT IN THE RUNFILE
      for key, value in kwargs.items():
        key = key.lower()
        if isinstance(value, str): value = value.lower()
        #print('key, line[0]',key, line[0])
        if key == line[0].lower():
          found_keys.append(key)
          if value == 'disable' or value is None:
            nline = '    ' + '{:<13}'.format('!' + key)  + '\n'
          elif value == 'enable':
            pass 
          else:
            nline = '     ' + '{:<13}'.format(key) + ' : ' + str(value) + '\n'
          break 
        
      f.write(nline)
      
    for key, value in kwargs.items():
      if key in ['cat', 'blocks', 'timer']:
        continue
      key = key.lower()
      if isinstance(value, str): value = value.lower()      
      if key not in found_keys:
        self.__log.warning('Parameter ' + key + 
                        ' not found and its value will be appended, not modified')
        f.write('     ' + '{:<13}'.format(key) + ' : ' + str(value) + '\n')
    
    f.close()  

  # -----------------------------------------------------------------------------
  
  def rename_key(self, old, new, **kwargs):
    """
    
    """
    cmd = 'sed -i -e "s/' + old + '/' + new + '/g" ' + self.fname
    o, e = bash(cmd)
    
    self.__log.debug(cmd)
    self.__log.debug(o)
    if len(e) > 0:
      self.__log.warning(e)
    
    
  # -----------------------------------------------------------------------------
@traced
@logged
class SegyPrepFile(ParameterFile):
  """
  A SegyPrep's runfile.
  
  """  
  def __init__(self, proj, path, **kwargs):
    """
    
    """
    self.name = proj.name + '-SegyPrep.key'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)

  # -----------------------------------------------------------------------------    
  
  def create(self, dummy=None, **kwargs):
    """
    Create a SegyPrep's runfile. 
    
    Notes
    -----
    ztype=ibfs => taking sources and receivers ', &
                  'depths as -selev, and +gelev respectively.

    """
    from fullwavepy.ioapi.generic import save_dict
    
    cat = kw('cat', True, kwargs)
    self.kwargs = kwargs
    proj = self.proj

    nx1, nx2, nx3  = proj.dims 
    dx = proj.dx 
    dt = proj.dt
    dt_ms = dt * 1000 # ms
    ns = proj.ns 
    ttime = proj.ttime 
    ttime_ms = ttime * 1000 # ms
    
    segyprep = {}
    
    segyprep['addtodepth'] = kw('addtodepth', 0, kwargs)
    
    ztype_default = 'ibfs'
    if 'ztype' not in kwargs:
      self.__log.warning('Setting z type to (kmc custom) default: %s' % str(ztype_default))
    segyprep['ztype'] = kw('ztype', ztype_default, kwargs)
   
    segyprep['problem'] = kw('problem', proj.problem, kwargs)
    segyprep['io'] = kw('io', proj.io, kwargs)  
    segyprep['reciprocity'] = kw('reciprocity', 0, kwargs)
    self.__log.warning('Reciprocity: ' + str(bool(segyprep['reciprocity'])))
    segyprep['fixedarray'] = kw('fixedarray', 'yes', kwargs)
    segyprep['unique'] = kw('unique', 'yes', kwargs)
    
    segyprep['minoffset'] = kw('minoffset', 0, kwargs) # IN METRES
    segyprep['maxoffset'] = kw('maxoffset', 1e6, kwargs) # IN METRES
    
    segyprep['FFID'] = kw('FFID', 'yes', kwargs) 
    segyprep['text'] = kw('text', 'yes', kwargs)
    segyprep['debug'] = kw('debug', 'yes', kwargs)
    segyprep['retain'] = kw('retain', 'yes', kwargs)
    segyprep['outseis'] = kw('outseis', 'yes', kwargs)
    segyprep['outsource'] = kw('outsource', 'yes', kwargs)
    
    segyprep['nx1'] = nx1
    segyprep['nx2'] = nx2
    segyprep['nx3'] = nx3
    segyprep['dx'] = dx
    segyprep['ttime'] = ttime_ms # YES, _ms!
    segyprep['dtms'] = dt_ms
    
    # GEOMETRY
    segyprep['geometry'] = kw('geometry', 'sgy', kwargs)
    if (segyprep['geometry'] == 'sgy') or (segyprep['geometry'] == 'segy'):
      segyprep = self.set_geom_segy(segyprep, **kwargs)
    elif segyprep['geometry'] == 'regular':
      segyprep = self.set_geom_regular(segyprep, **kwargs)
    else:
      raise ValueError('Unknown geometry: ' + str(geometry))
      
    save_dict(self.fname, segyprep)  
    self.di = segyprep
    # WARN ABOUT ANOMALOUS BEHAVIOUR
    if segyprep['text'] != 'yes':
      self.__log.warning("'text'=" + segyprep['text'] + ". No .hed files " +
                      "will be generated which will affect e.g. the phasubplots.")
                      
  # -----------------------------------------------------------------------------

  def set_geom_segy(self, segyprep, **kwargs):
    segyprep['geometry'] = 'segy'
    segyprep['xorigin'] = self.proj.box[0]
    segyprep['xshift'] = 0
    segyprep['yorigin'] = self.proj.box[2]
    segyprep['yshift'] = 0
    return segyprep
  
  # -----------------------------------------------------------------------------  
  
  def set_geom_regular(self, segyprep, **kwargs):
    geometry_in_nodes = kw('geometry_in_nodes', True, kwargs)
    
    if geometry_in_nodes:
      self.__log.warning('geometry_in_nodes=True => ' + 
                         'expecting regular geometry specified in nodes')
      a = -1
      b = self.proj.dx
      az = 0
      bz = b
    else:
      self.__log.info('Expecting regular geometry specified in metres')
      a = 0
      b = 1
      az = a
      bz = b
      
    self.__log.debug(kwargs)
    segyprep['souz']  = (float(kwargs['souz'] ) + az) * bz
    segyprep['recz']  = (float(kwargs['recz'] ) + az) * bz       
    segyprep['soux0'] = (float(kwargs['soux0']) + a) * b
    segyprep['soudx'] = (float(kwargs['soudx'])) * b
    segyprep['sounx'] = float(kwargs['sounx'])
    segyprep['recx0'] = (float(kwargs['recx0']) + a) * b
    segyprep['recnx'] = float(kwargs['recnx'])
    segyprep['recdx'] = (float(kwargs['recdx'])) * b
    if self.proj.dims[1] > 1:
      segyprep['souy0'] = (float(kwargs['souy0']) + a) * b
      segyprep['soudy'] = (float(kwargs['soudy']) + a) * b
      segyprep['souny'] = float(kwargs['souny'])
      segyprep['recy0'] = (float(kwargs['recy0']) + a) * b
      segyprep['recny'] = float(kwargs['recny'])
      segyprep['recdy'] = (float(kwargs['recdy']) + a) * b
    else:
      segyprep['souy0'] = 0
      segyprep['soudy'] = 0
      segyprep['souny'] = 1        
      segyprep['recy0'] = 0
      segyprep['recny'] = 0
      segyprep['recdy'] = 1  
    return segyprep
    
  # -----------------------------------------------------------------------------    
    
  @timer
  def run(self, overwrite=False, **kwargs):
    """
    overwrite : bool 
    
    """
    cat = kw('cat', True, kwargs)
    
    if 'segyprep' not in self.proj.exe:
      raise AttributeError("You need to define proj.exe['segyprep']")
    
    if overwrite:
      answer = str(1)
    else:
      answer = str(0)
    
    cmd = str('printf "' + answer + '\n" | ' + 
              self.proj.exe['segyprep'] + ' ' + self.proj.name)
     
    self.__log.info('Running SegyPrep...')
    o, e = bash(cmd, path=self.path)
    if cat:
      print(o, e)
    
    errors = ['Unable to decode SEG-Y headers']
    raise_error = False
    errors_raised = 'Errors raised: \n'
    for err in errors:
      if err in o:
        raise_error = True
        errors_raised += ' %s' % err
    
    if raise_error:
      raise ValueError(str(o) + errors_raised)
    
  # -----------------------------------------------------------------------------  
  
  @timer
  def postprocess(self, **kwargs):
    pass
  
  

  # -----------------------------------------------------------------------------

  def log(self, **kwargs):
    """
    Show the content of the log file.
    
    """
    o, e = bash('cat ' + self.fname[ :-len('key')] + 'log')
    print(o, e)

  # -----------------------------------------------------------------------------
@traced
@logged
class Runfile(ParameterFile):
  """
  Fullwave3D runfile.

  """  
  def __init__(self, proj, path, **kwargs):
    """
    Suffix is for reading CP runfiles too.
    
    """
    suffix = kw('suffix', 'Runfile', kwargs)
    self.name = proj.name + '-' + suffix + '.key'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)
    
    if suffix == 'Runfile' and exists(self.fname): # NOT SKELETON
      self.blocks = self.read_blocks(**kwargs)
      self.blocks2iters(**kwargs)
    
  # -----------------------------------------------------------------------------  

  def set_proj_params(self, **kwargs):  
    """
    
    """
    d = self.read(**kwargs)
    
    try:
      self.proj.problem = d['problem'].lower()
    except KeyError:
      pass
    
    try:
      self.proj.domain = d['domain'].lower()
      self.__log.debug('Set proj.domain to value from runfile: ' + 
                       str(self.proj.domain))
    except KeyError:
      pass
    
    try:
      self.proj.dim = d['dim'].lower()
      self.__log.debug('Set proj.dim to value from runfile: ' + 
                       str(self.proj.dim))      
    except KeyError:
      pass
    
    try:
      self.proj.equation = d['equation'].lower()
      self.__log.debug('Set proj.equation to value from runfile: ' + 
                       str(self.proj.equation))      
    except KeyError:
      pass
    
    try:
      self.proj.anisotropy = d['anisotropy']
      self.__log.debug('Set proj.anisotropy to value from runfile: ' + 
                       str(self.proj.anisotropy))      
    except KeyError:
      pass
    
    try:
      self.proj.kernel = d['kernel'].lower()
      self.__log.debug('Set proj.kernel to value from runfile: ' + 
                       str(self.proj.kernel))      
    except KeyError:
      pass
    
    try:
      self.proj.nx1 = int(d['nx1'])
      self.__log.debug('Set proj.nx1 to value from runfile: ' + 
                       str(self.proj.nx1))      
    except KeyError:
      pass
    
    try:
      self.proj.nx2 = int(d['nx2'])
      self.__log.debug('Set proj.nx2 to value from runfile: ' + 
                       str(self.proj.nx2))      
    except KeyError:
      pass
    
    try:
      self.proj.nx3 = int(d['nx3'])
      self.__log.debug('Set proj.nx3 to value from runfile: ' + 
                       str(self.proj.nx3))      
    except KeyError:
      pass
    
    try:
      self.proj.dx = float(d['dx'])
      self.__log.debug('Set proj.dx to value from runfile: ' + 
                       str(self.proj.dx))      
    except KeyError:
      pass

    try:
      self.proj.io = d['io'].lower()
      if self.proj.io == 'segy':
        self.proj.io = 'sgy'
      self.__log.debug('Set proj.io to value from runfile: ' + 
                       str(self.proj.domain))      
    except KeyError:
      pass
    
    try:
      self.proj.btop       = int(d['btop'])  
      self.proj.bbot       = int(d['bbot'])  
      self.proj.bleft      = int(d['bleft'])  
      self.proj.bright     = int(d['bright'])  
      self.proj.bfront     = int(d['bfront'])  
      self.proj.bback      = int(d['bback'])  
      self.proj.etop   = int(d['etop'])  
      self.proj.ebot   = int(d['ebot'])  
      self.proj.eleft  = int(d['elef'])  
      self.proj.eright = int(d['erig'])  
      self.proj.efront = int(d['efro'])  
      self.proj.eback  = int(d['ebac']) 
    except KeyError as err:
      self.__log.warning('Description of boundaries in the Runfile is not complete: ' + str(err))

  # -----------------------------------------------------------------------------
  
  #FIXME: dummy
  def prep(self, dummy=None, **kwargs):
    """
    
    """
    cat = kw('cat', True, kwargs)
    self.__log.debug('Setting cat=0')
    kwargs['cat'] = 0 # WILL CAT IN THE END
    super().prep(**kwargs)
    
    self.read(**kwargs)
    params2overwrite = {}
    for param in ['nx1', 'nx2', 'nx3', 'dx']:
      p_self = float(self.params[param])
      p_proj = float(getattr(self.proj, param))
      if p_self != p_proj:
        self.__log.warning('Found inconsistency for ' + param + '\n' + 
                         self.fname + '.' + param + ' = ' + str(p_self) + 
                         '\n' +'proj.' + param + ' = ' + str(p_proj))
        params2overwrite[param] = p_proj
    
    if len(params2overwrite) != 0:
      self.modify(**params2overwrite)
    
    if cat:
      self.cat()
    
  # -----------------------------------------------------------------------------

  def create(self, dummy=None, **kwargs):
    """
    1. Copy a template
    2. Update with values from Skeleton
    3. Update with extended grid.
    4. Set iteration blocks.
    
    
    """
    from fullwavepy.ioapi.generic import read_dict
    skelet_fname = self.path + self.pname + '-Skeleton.key'
    self.__log.info('Using some values from %s. To modify that file, '\
                    'update SegyPrep.key and re-run SegyPrep' % skelet_fname)
                    
    btop = kw('btop', 0, kwargs)
    etop = kw('etop', 0, kwargs)  
    del_kw('btop', kwargs)
    del_kw('etop', kwargs)
    # b_abs = kw('b_abs', 40, kwargs)
    # e_abs = kw('e_abs', 50, kwargs)
    # del_kw('b_abs', kwargs)
    # del_kw('e_abs', kwargs)
    
    b = {}
    for k in ['bbot', 'bleft', 'bright', 'bfront', 'bback']:
      if k in kwargs and ('b_abs' not in kwargs):
        b[k] = kwargs[k]
      elif 'b_abs' in kwargs:
        b[k] = kwargs['b_abs']
      else:
        b[k] = 40 # default 
      
      del_kw(k, kwargs)

    e = {}
    for k in ['ebot', 'elef', 'erig', 'efro', 'ebac']:
      if k in kwargs and ('e_abs' not in kwargs):
        e[k] = kwargs[k]
      elif 'e_abs' in kwargs:
        e[k] = kwargs['e_abs']
      else:
        e[k] = 50 # default 
      
      del_kw(k, kwargs)

    
    # READ THE TEMPLATE
    tmplt = self._get_template(**kwargs)
    with open(self.fname, 'w') as f:
      f.write(tmplt)
    
    # MERGE WITH SKELETON (CREATED BY SEGYPREP) IF PRESENT
    if exists(skelet_fname):
      skeleton = read_dict(skelet_fname)
      
      for key in ['ttime', 'end time', 'max time']:
        try:
          etime = skeleton[key]
        except KeyError:
          continue

      try:
        self.modify(io=skeleton['io'])
      except KeyError:
        pass
      self.modify(io=skeleton['io'], 
                  etime=etime,
                  NX1=skeleton['nx1'], 
                  NX2=skeleton['nx2'], 
                  NX3=skeleton['nx3'], 
                  DX=skeleton['dx'], 
                  NRECS=skeleton['nrecs'], 
                  MAXRC=skeleton['maxrc'], 
                  NCOMP=skeleton['ncomp'], 
                  NSHOTS=skeleton['nshot'])    
    else:
      self.__log.warning('File ' + skelet_fname + ' not found. ' + 
                      'No merge with template will be applied.')
    
    # APPLY THE REMAINDER 
    self.__log.debug('Setting other parameters, including boundary conditions...')
    self.modify(problem=self.proj.problem, 
                anisotropy=self.proj.anisotropy, 
                kernel=self.proj.kernel, 
                dim=self.proj.dim,
                domain=self.proj.domain,
                equation=self.proj.equation,
                units=self.proj.units,
                btop=str(btop),
                bbot=str(b['bbot']), 
                bleft=str(b['bleft']), 
                bright=str(b['bright']), 
                bfront=str(b['bfront']), 
                bback=str(b['bback']), 
                etop=str(etop), 
                ebot=str(e['ebot']), 
                elef=str(e['elef']), 
                erig=str(e['erig']), 
                efro=str(e['efro']), 
                ebac=str(e['ebac']),
                **kwargs)        
    
    self._set_iters(**kwargs)

  # -----------------------------------------------------------------------------  
    
  def _get_template(self, **kwargs):    
    """
    Get the template runfile.
    
    Notes
    -----
    Apparently bleft and elef 
    (etc.) is the only format
    allowed.
    
    To add:
    ! X. FREE SURFACE
     ibfs          : 1
     seaLevel      : 0
     maxGhostIter  : 2
     vacuum        : 0
     accuracy      : 0.001
     minGhostFS    : 0.001
     minFictFS     : 0.5
     interpMode    : 1
     srcType       : 0
     recType       : 0
    
    """
    tmplt = """
    ! May 2019, K. Chrapkiewicz 
    ! THIS RUNFILE CONFORMS TO
    ! FULLWAVE REV. 688 STANDARD
    
    
    ! A. PROBLEM DEFINITION
     problem       : synthetic
     domain        : time
     dim           : 3d
     equation      : acoustic
     units         : metric
     anisotropy    : none
     kernel        : low
    
     
    ! B. MODEL DEFINITION
     nx1           : 30
     nx2           : 20
     nx3           : 40
     dx            : 50.0
    
     
    ! C. DATA DEFINITION
     ncomp         : 1
     nshots        : 1
     nrecs         : 22
     maxrc         : 22
     maxps         : 1
     io            : fw3d 
    
     
    ! D. BOUNDARY CONDITIONS
     ibfs          : 0
     multisurf     : 0
     nosprdfctrs   : 0     
     seaLevel      : 0
     maxGhostIter  : 2
     vacuum        : 0
     accuracy      : 0.001
     minGhostFS    : 0.001
     minFictFS     : 0.5
     interpMode    : 1     
     
     btop          : 0
     bbot          : 2
     bleft         : 2
     bright        : 2
     bfront        : 2
     bback         : 2
     
     etop          : 0
     ebot          : 20
     elef          : 20
     erig          : 20
     efro          : 20
     ebac          : 20
    
     
    ! F. MUTES & DATA SELECTION
     etime         : 10
    
    
    ! G. INVERSION PARAMETERS
     iprop         : vp       ! INVERT PROPERTIES
     func          : twonorm  ! FUNCTIONAL
     agc           : off      ! DEFAULT: off
     amplitude     : no       ! DEFAULT: no
     normalise     : yes      ! DEFAULT: yes
     spatial       : yes      ! SPATIAL PRECONDITIONER; DEFAULT: yes
     slowness      : yes      ! DEFAULT: yes
     conj          : no       ! CONJUGATE GRADIENT; DEFAULT: no
     gwidth        : 0.0      ! GAUSSIAN WIDTH; DEFAULT: 0.0
     
     !OMITTED: huber
    
     
    ! H. SMOOTHING & WINDOWING
     usewin        : no
     win1          : 1
     win2          : 1
     win3          : 10000
     win4          : 10000
     
     sx1           : 2.0
     sx2           : 2.0
     sx3           : 1.0
     
     !OMITTED: smooth,
     !         asx1, asx2, asx3, fins, finsx1, finsx2, finsx3,
     !         finscut, finsramp
    
    
    ! I. P-WAVE VELOCITY MODEL
     velcut        : 1550     ! [m/s] DO NOT UPDATE FOR ALL x : ANY_PARAM(x) < velcut  
     
     !OMITTED: minvc, maxvc
    
     
    ! J. DENSITY MODEL
     wvel          : 1550     ! [m/s]   WATER VEL. (SETS wden FOR ALL x : Vp(x) <= wvel)
     wden          : 1000     ! [kg/m3] WATER DEN.
     gvel          : 1600     ! [m/s]   GARDNER CUTOFF VEL.
    
     !OMITTED: max water cell (gcell), min water cell
     
    
    ! K. ANISOTROPY MODEL
     !OMITTED: minec, maxec, mindc, maxdc, ttiu, angle
    
     
    ! L. ELASTIC MODEL
     !OMITTED: minvsc, maxvsc, vp to vs, vpzero, 
     
     
    ! M. EPICMods 
     !OMITTED: ALL
     
    
    ! N. GLOBAL INVERSION VIA qPSO
     !OMITTED: ALL 
    
    
    ! E. ITERATIONS & DATA SELECTION (THE LAST PART OF THE FILE)
     nblock         : 3
     nits           : 10
     freq           : 3.0     
     minoff         : 2000
     maxoff         : 20000
     
     nits           : 5
     freq           : 4.0
     minoff         : 2000
     maxoff         : 20000
     
     nits           : 5
     freq           : 5.0
     minoff         : 2000
     maxoff         : 20000
    """
    
    return tmplt
    
  # -----------------------------------------------------------------------------
  
  def _set_iters(self, **kwargs):
    """
    
    blocks : list
      List of blocks, each is a dictionary which 
      has to contain at least 
      - nits (no. of iterations)
      - freq (frequency in Hz)
      Order of this list obviously matters.
    
    1. Read the runfile
    2. Open a new file.
    3. Transcribe all the content except the iteration-blocks section
    4. Write new iteration-blocks section.
    
    
    """
    from fullwavepy.ioapi.generic import read_txt_raw
    
    if self.proj.problem == 'tomography':
      nits = kw('nits_per_block', 20, kwargs)
      minoff = kw('minoff', 5000, kwargs)
      blocks = kw('blocks', [{'nits': nits, 'freq': 3.0, 'minoff': minoff},
                             {'nits': nits, 'freq': 3.5, 'minoff': minoff},
                             {'nits': nits, 'freq': 4.0, 'minoff': minoff},
                             {'nits': nits, 'freq': 4.5, 'minoff': minoff},
                             {'nits': nits, 'freq': 5.0, 'minoff': minoff},
                             {'nits': nits, 'freq': 5.5, 'minoff': minoff},
                             {'nits': nits, 'freq': 6.0, 'minoff': minoff},
                             {'nits': nits, 'freq': 6.5, 'minoff': minoff},
                            ], kwargs)
      self.blocks = blocks # TO BE USED BY DUMPCOMPARE.read
    elif self.proj.problem == 'synthetic':
      self.blocks = kw('blocks', [], kwargs)
    else:
      raise ValueError('Unknown problem type: %s' %str(self.proj.problem))
    
    self.blocks2iters(**kwargs)
    
    fname = self.fname
    content = read_txt_raw(fname)
    
    
    with open(fname, 'w') as f:
      # REWRITE ALL SECTION BEFORE SECTION E (ITERATION BLOCKS)
      for line in content:
        nline = line
        line = line.split(None)
        if len(line) > 0 and line[0] == 'nblock':
          break
        else:
          f.write(nline) # NOTE
      
      # ADD CUSTOM SECTION E
      f.write('     ' + '{:<13}'.format('nblock') + ' : ' + str(len(self.blocks)) + '\n\n')
      
      for i, block in enumerate(self.blocks):
        if 'freq' not in block:
          raise ValueError('Specify freq for block no. ' + str(i+1))
        if 'nits' not in block:
          raise ValueError('Specify nits for block no. ' + str(i+1))
        
        # START WITH nits, OTHERWISE COULD SORT DIFFERENTLY
        key = 'nits'
        val = block[key]
        f.write('     ' + '{:<13}'.format(key) + ' : ' + str(val) + '\n')
        
        del block[key]
        
        for key, val in block.items():
          f.write('     ' + '{:<13}'.format(key) + ' : ' + str(val) + '\n')
        
        f.write('\n')
        
  # -----------------------------------------------------------------------------

  def read_blocks(self, new_block_activator='nits', **kwargs):
    """
    Read iteration blocks.
    
    Notes
    -----
    It assumes that each block starts with a line setting a 'nits' parameter.
    
    If no blocks are present (like for synthetic projects), empty list 
    is naturally returned.
    
    """
    from fullwavepy.ioapi.generic import read_txt
    

    content = read_txt(self.fname)
    
    self.__log.info('new_block_activator=%s => each iteration block must begin with %s keyword.' % \
      (new_block_activator, new_block_activator))
    self.__log.debug(str(content))
    
    

    blocks = []
    for line in content:
      if len(line) < 3:
        continue
      
      key = line[0]
      colon = line[1]
      val = line[2]
      
      if key == new_block_activator:
        self.__log.debug('Line %s starts with key %s => a new iteration block detected' % \
          (line, new_block_activator))
        blocks.append({})
        current_block = blocks[-1] 
      
      if len(blocks) > 0: 
        try:
          val = float(val)
        except ValueError as err:
          self.__log.debug('Non-numeric parameter: ' + str(err))
        
        current_block[key] = val
        
    self.blocks = blocks
    return blocks
  
  # -----------------------------------------------------------------------------
  
  def blocks2iters(self, **kwargs):
    """
    Prepare a list of all iterations, each
    having info about the block to which it belongs
    """
    self.iters =  [None]    
    for b in self.blocks:
      for i_b in np.arange(b['nits']):
        self.iters.append(b)

  # ----------------------------------------------------------------------------- 
 
  def read_nits(self, **kwargs):
    """
    Read total no. of iterations. 

    Notes
    -----
    We don't parse keys with spaces.
    That's why 'nits' is OK, but 'no of iterations' is not.
    
    """
    from fullwavepy.ioapi.generic import read_txt
    
    content = read_txt(self.fname)
    
    nits_total = 0
    
    for line in content:
      if len(line) > 0 and line[0] == 'nits':
        nits_total += int(float(line[2]))
      
    self.nits_total = nits_total
    return nits_total

  # -----------------------------------------------------------------------------
@traced
@logged
class Skeleton(Runfile):
  """
  
  """
  def __init__(self, proj, path, **kwargs):
    suffix = 'Skeleton'
    super().__init__(proj, path, suffix=suffix, **kwargs)

  # ----------------------------------------------------------------------------- 

