"""
This module serves to read/write files in formats 
defined by Fullwave3D's segy IO API, 
in particular SEG-Y files.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged #, traced

from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, strip, exten, \
  path_leave, path_extract
from fullwavepy.ioapi.generic import ArrayFile, CsvFile

json_header_suffix = '_HEAD.json'
filt_suffix = '_filt.sgy'
filt_mute_suffix = '_filt_mute.sgy'

# -------------------------------------------------------------------------------
# SGY
# -------------------------------------------------------------------------------
@logged
class SgyMapp(dict):
  pass
@logged
class SgyFile(ArrayFile):
  """
  SEG-Y file.
  
  Notes
  -----
  - it should be used as a second-parent class (multi-inheritance).
  - it just adds some methods for dealing with 
  SEGY files, usually through calling SeismicUnix.
  - it is meant to have no overlap with other classes 
  to avoid multi-inheritance problems.
  - .csv seems the best choice for storing the trace headers 
  since SEGY files are very structered (so json is an overkill)
  plus Pandas reads .csv more efficiently (only selected columns).
  
  """
  @classmethod
  def new(cls, name, path, **kwargs):
    return cls(name, path, **kwargs)
  def extract(self, suffix, *args, path=None, **kwargs):
    if path is None:
      fname = self.fname
    else:
      fname = path + self.name
    nfname = '%s_%s.%s' %\
      (strip(fname), suffix, exten(self.fname))
    kwargs['nfname'] = nfname
    d = self.read(*args, **kwargs)
    return self.new(path_leave(nfname), path_extract(nfname))
  def filt(self, pad, overwrite=False, **kwargs):
    """
    CAUTION
    It overwrites the content of the file.
    
    """
    from fullwavepy.ioapi.su import sushw
    from fullwavepy.dsp.su import su_filter_full
    
    self.__log.info('Using ' + str(pad) + ' samples of padding')
    self.__log.info('Setting dt in the header of: ' + self.fname)
    
    if hasattr(self, 'proj'):
      self.__log.debug('Taking dt = self.proj.dt')
      dt = self.proj.dt
    else:
      dt = kwargs['dt']

    sushw(self.fname, 'dt', (dt*1e6), **kwargs)
    # THIS PREVENTS OVERWRITING AND OUTPUTS INTERMEDIATE STEPS TOO
    fname_out = su_filter_full(self.fname, pad, **kwargs)
    self.__log.info('Filtered data output to ' + fname_out)
    self.__log.warning('Overwriting ' + self.fname)
    o, e = bash('mv {} {}'.format(fname_out, self.fname))
  def mute(self, fbreaks, ntaper=100, twin=1, **kwargs):
    """
    """
    from fullwavepy.dsp.su import su_mute # not used actually!
    from fullwavepy.ioapi.generic import save_txt
    
    fbreaks = np.array(fbreaks)
    picks = fbreaks * self.proj.dt
    bpicks = picks + twin
    nmute = len(picks)
    xmute = range(1, nmute + 1)
    tmute = picks
    tmute2 = bpicks
    
    # THIS LOOP ENSURES WE MUTE THE WHOLE TRACE 
    # IF THE FIRST ARRIVAL IS LATER THEN TOTALTIME (AS INDICATE BY pick==0)
    # WITHOUT THIS, IT WAS BUGGY!
    for i, pick in enumerate(picks):
      if pick == 0:
        bpicks[i] = 0    
    
    xmute = [str(i) for i in xmute]
    tmute = [str(i) for i in tmute]
    tmute2 = [str(i) for i in tmute2]
    
    for data, prefix in zip([xmute, tmute, tmute2], ['xmute', 'tmute', 'tmute2']):
      file_txt = prefix + '.txt'
      file_bin = prefix + '.bin'
      save_txt(file_txt, data)
      o, e = bash('a2b < {} n1=1 > {}'.format(file_txt, file_bin))
    
    fname_out = strip(self.fname) + '_tmp' + exten(self.fname)
    
    cmd =  'segyread tape={} | '.format(self.fname)
    cmd += 'sumute key=tracr nmute={nmute} mode=0 ntaper={ntaper} xfile={xmute_bin} tfile={tmute_bin} | sumute key=tracr nmute={nmute} mode=1 ntaper={ntaper} xfile={xmute_bin} tfile={tmute2_bin}'.format(nmute=nmute, ntaper=ntaper, xmute_bin='xmute.bin', tmute_bin='tmute.bin', tmute2_bin='tmute2.bin')
    cmd += ' | segyhdrs | segywrite tape={}'.format(fname_out)
    
    o, e = bash(cmd)
    o, e = bash('mv ' + fname_out + ' ' + self.fname)
  def read(self, overwrite=True, **kwargs):
    """
    To read and plot! windowed data in one line.
    """
    win = kw('win', None, kwargs)
    selfwin = getattr(self, 'win', None)

    if (win != selfwin) and (overwrite == False):
      overwrite = True
      self.__log.warning('Changed overwrite to {} because provided win {} is different from the previous {}'.format(overwrite, win ,selfwin))
    
    if (not hasattr(self, 'array')) or overwrite:
      if win is None:
        self.array = super().read(**kwargs)
      else:
        kwargs['win'] = win # OTHERWISE MULTIPLE VALUES
        self.array = self.window(**kwargs)
    
    self.__log.debug('self.array.extent %s' % str(self.array.extent))
    return self.array
  def read_header(self, overwrite=True, **kwargs):
    """
    Add reading csv, useful for heavy sgy files.
    """
    if (not hasattr(self, 'head')) or overwrite:
      self.head = header2csv(self.fname, keys='all', suffix='_HEAD', **kwargs)
    return self.head
  def rh(self, *args, **kwargs):
    return self.read_header(*args, **kwargs)
  def resize(self, box, **kwargs): #FIXME?
    """
    Cut the model to fit the proj.box.
    
    Parameters
    ----------
    file_z0 : integer
      physical coordinate in m of the 
      first (0th) sample in the file.
      Can be < 0.
      
    kwargs
    
    Returns
    -------
    None
    
    Notes
    -----
    First make sure the proj.box is correct!
    
    Although dt is defined in SEGY
    headers in microseconds, tmin and tmax are in
    seconds (these are not header words, and can be 
    floats).
    dt is expected to be in metres. If z1 and z2 are 
    also in metres, we have to divide them by 1e6,
    because suwind tmin=... takes values in microseconds!
    
    """    
    self.__log.info('Assuming integer box coords, as required by SEGY')
    x1, x2, y1, y2, z1, z2 = box
    
    #self.__log.warning('\n\n DISABLED BUGGY file_z0 CONVERSION!!!\n\n')
    #z1 -= file_z0
    #z2 -= file_z0
    
    self.__log.debug('z1={}, z2={}'.format(z1, z2))
    
    dt = int(self._gethw('dt', unique_values=True, timer=True, **kwargs)[0])
    self.__log.debug('Converting z1,z2 into microsec as required by suwind')
    if dt > 1000:
      self.__log.info('Header dt > 1000. Assuming miliseconds or milimetres')
      z1 /= 1e3
      z2 /= 1e3
    elif dt <= 1000:
      self.__log.info('Header dt <= 1000. Assuming seconds or metres')
      z1 /= 1e6
      z2 /= 1e6

    self.__log.debug('z1={}, z2={}'.format(z1, z2))
    
    scalco = int(self._gethw('scalco', unique_values=True, timer=True, **kwargs)[0])
    if scalco < 0: 
      scalco = abs(scalco)
    elif scalco == 0:
      scalco = 1
    else: # scalco > 0 MEANS IT IS USED AS A MULTIPLIER NOT A DIVISOR IN THE HEADER
      scalco = 1 / abs(scalco) # WE NEED TO DO THE OPPOSITE THING TO OUR BOX
      
    self.__log.debug('scalco' + str(scalco))
    x1, x2, y1, y2 = np.array([x1, x2, y1, y2]) * scalco
    
    key_x = self.proj.sgyhw['xmod']
    key_y = self.proj.sgyhw['ymod'] 
    
    tmp_fname = 'tmp.sgy'
    cmd = str('segyread tape=' + self.fname + ' | ' +
              'suwind key=' + key_x + 
              ' min=' + str(x1) + 
              ' max=' + str(x2) + ' | ' + 
              'suwind key=' + key_y + 
              ' min=' + str(y1) + 
              ' max=' + str(y2) + ' | ' +  
              'suwind tmin=' + str(z1) + ' tmax=' + str(z2) + ' | ' + 
              'segyhdrs | ' +
              'segywrite tape=' + tmp_fname)
    
    self.__log.debug(cmd)
    o, e = bash(cmd)
    o, e = bash('mv ' + tmp_fname + ' ' + self.fname) 
  def split(self, key, **kwargs):
    from fullwavepy.ioapi.su import suwind
    from fullwavepy.generic.parse import extend_fname
    
    setattr(self, key, {})
    selfkey = getattr(self, key)
    
    values = self._gethw(key, unique_values=True, **kwargs)
    
    for value in values:
      nfname = extend_fname(self.fname, [[key, value]])
      suwind(self.fname, nfname, {key: [value]}, **kwargs)
      selfkey[value] = self.__class__(path_leave(nfname), self.path)
  def surange(self, **kwargs):
    if not exists(self.fname):
      raise FileNotFoundError(self.fname)    
    
    o, e = bash('segyread tape=' + self.fname + ' | ' +
                'surange', **kwargs)
    print(e + '\n' + o) # NOT log.info IN ORDER TO BE ALWAYS PRINTED
  def suximage(self, **kwargs):
    if not exists(self.fname):
      raise FileNotFoundError(self.fname)    
    #FIXME
    self.__log.warning('suximage does not work here for some reason but it used to')
    o, e = bash('segyread tape=' + self.fname + ' | ' +
                'suximage', **kwargs)
    self.__log.info(e + '\n' + o)
  def window(self, win, nfname=None, **kwargs):
    """
    Extract ('window') a subset of data
    from the SEGY file using SU's suwind
    utility.

    Parameters
    ----------
    win : dict
      Dictionary in which keys are
      SEGY-header keywords used for windowing and
      each value may have 1 of 2 formats:
      (1) key: {'min': 10, 'max': 40} # NOTE a subdict!
      (2) key: [value1,value2,...] # list also for singles [value1]
    nfname : str
      If None, ...?    
    """
    from fullwavepy.ioapi.su import suwind
    
    self.win = win
    
    fname = self.fname
    if nfname is None:
      nfname = strip(fname) + '_windowed.' + exten(fname)
      path = self.path
    else:
      nfname = nfname 
      path = path_extract(nfname)
    
    suwind(fname, nfname, win, **kwargs)
    
    obj = ArrayFile(path_leave(nfname), path)
    self.array = obj.read(**kwargs)
    
    return self.array
  # -----------------------------------------------------------------------------  
  def _gethw(self, key, unique_values=False, **kwargs):
    """
    Get header-word values.
    
    Returns
    -------
    hw_values : list
      List of floats: either one per trace
      in order of traces (if unique_values=False)
      or as many as unique values, sorted.
    
    """
    from fullwavepy.ioapi.su import sugethw
    values = sugethw(self.fname, keys=[key], **kwargs)[key]
    if unique_values:
      values = sorted(list(set(values)))
    return values
  def _get_sr_coords(self, datafile=None, **kwargs):
    """
    Get HORIZONTAL (x,y) coordinates of sources and 
    receivers.
    
    datafile : DataFile 
      Different sgy file to 
      parse a header from.
    
    """
    self.__log.debug('Getting sx,sy,gx,gy from header')
    if datafile is None:
      datafile = self
      
    sx = datafile._gethw('sx', **kwargs)
    sy = datafile._gethw('sy', **kwargs)
    gx = datafile._gethw('gx', **kwargs)
    gy = datafile._gethw('gy', **kwargs)
    
    #plt.plot(sx, sy, '.')
    #plt.plot(gx, gy, '.')
    s = list(zip(sx, sy))
    r = list(zip(gx, gy))
    self.s = s
    self.r = r
    
    return s, r
@logged
class SgyHeader(CsvFile):
  sid = 'fldr'
  rid = 'tracf'
  lid = 'ep'
  def source_data(self):
    """
    It assumes each receiver has the same source data (same shots).
    """
    return self[self.tracf==1104] # whattt???
@timer
@logged
def split_sgy(fname, key, value=None, **kwargs):
  """
  Split a sgy file into smaller 
  files, one file for each value of 
  the key.
  
  if no value is provided, it will take 
  all values present in the file.
  
  Returns
  -------
  nfnames : list 
    List of files resulted from 
    the splitting. Each element
    can be feed into another splitting.
  
  Notes
  -----
  Useful for storing shot lines 
  separately for displaying.
  
  """
  from .su import suwind, sugethw
  from fullwavepy.generic.parse import extend_fname
  
  if not exists(fname):
    raise FileNotFoundError(fname)  
  
  if value is None:
    values = sugethw(fname, key, unique_values=True, **kwargs)
  else:
    values = [value]
  
  nfnames = []

  for value in values:
    nfname = extend_fname(fname, [[key, value]])
    nfnames.append(nfname)
    split_sgy._log.info('Output file: ' + nfname)
    
    suwind(fname, nfname, key, value, value, **kwargs)

  return nfnames 
@logged
def array2sgy(fname, A, dt, **kwargs):
  """
  Export an array to a .sgy file.
  
  Parameters
  ----------
  fname : str 
    SEGY file to be created.

  Notes
  -----
  It doesn't preserve a 3D structure.
  
  """
  from .fw3d import save_vtr
  
  # HANDLE VARIOUS DIMENSIONALITIES (-> func?)
  shape = A.shape
  if len(shape) == 3:
    An = A
  elif len(shape) == 2:
    array2sgy._log.warning('Array 2D detected, adding Y-dimension of length 1')
    An = np.zeros((shape[0], 1, shape[1]))
    An[:,0,:] = A
  elif len(shape) == 1:
    array2sgy._log.warning('Array 1D detected, adding X- and Y- dimensions of length 1')
    An = np.zeros((1, 1, shape[0]))
    An[0,0,:] = A
  else:
    raise IOError('Array has wrong no. of dimensions: ' + str(len(shape)))
  
  # ACTUAL JOB
  fvtr = strip(fname) + '.vtr'
  save_vtr(An, fvtr)
  vtr2sgy(fvtr, dt)
@logged
def vtr2sgy(fname, dt, **kwargs): 
  """
  fname : str 
    Name of the vtr file.
  
  dt : float 
    In seconds.
  
  Notes
  -----
  It is more tricky to convert the model # FIXME
  preserving nx,ny,nz.
  
  CAUTION
  fortran vtr2sgy fails to handle 1-letter file-names
  like d.vtr - it outputs d.vtr.sgy instead.
  
  """
  from .su import sushw
  
  if not exists(fname):
    raise FileNotFoundError(fname)  
  
  cmd = str('printf "yes\n' + fname + '\n\n\nyes\n" | ' + 
            'vtr2sgy')
  vtr2sgy._log.debug('cmd: ' + cmd)
  o, e  = bash(cmd)  
  #print(o, e)
  
  fsgy = strip(fname) + '.sgy'
  dtmicros = int(dt * 1e6)
  sushw(fsgy, 'dt', dtmicros)
@logged
def sgy2vtr(fname, nx=None, **kwargs):
  """
  Convert sgy to vtr using Fullwave3D's 
  sgy2vtr Fortran utility.
  
  Parameters
  ----------
  nx : int
    No. of inline nodes. This is needed 
    to split the SEG-Y in a way that retaining its 3D
    structure (as we want especially for 
    models).
    Default: None => set as a total no. of traces.
  
  Returns
  -------
  None
  
  Notes
  -----
  We treat the data as 2D here, i.e.
  having a shape: (ntraces, nsamples)
  by passing no. of all traces (ntraces)
  as no. of in-line traces (nx).
  
  """
  from .su import get_ntraces
  
  if not exists(fname):
    raise FileNotFoundError(fname)  
  
  if nx is None:
    sgy2vtr._log.debug('nx not provided, assuming nx=ntraces')
    nx = get_ntraces(fname, **kwargs)
  
  cmd = str('printf "yes\n' + fname + '\n' + str(nx) + '\n\nyes\n" | ' + 
            'sgy2vtr ' + fname) 
  o, e = bash(cmd)
  
  if len(e) > 0:
    sgy2vtr._log.warning(e)
  
  fname_vtr = strip(fname) + '.vtr'
  if not exists(fname_vtr):
    print(o, e)
    raise OSError('Conversion failed, most likely due to incorrect nx: %s' % nx)
@timer
@logged
def read_sgy(fname, overwrite=True, **kwargs):
  """
  Read a .sgy file through 
  conversion to .vtr.
  
  Parameters
  ----------
  fname : str    
    File name. It should include  extension. 
    It should include path if needed.
  **kwargs : keyword arguments (optional)
    Current capabilities:
      
  Returns
  -------
  data : array
    3D data array as read by read_vtr.
    Typically it will be (ntraces, 1, nsamples),
    see sgy2vtr and read_vtr for details.
  
  Notes
  -----
  
  """
  from fullwavepy.generic.system import exists
  from .fw3d import read_vtr
  
  if not exists(fname):
    raise FileNotFoundError(fname + ' does not exist. Will NOT look ' +
                           'for an associated vtr file')  
  
  if 'nx' in kwargs:
    if 'shape' in kwargs and kwargs['shape'][0] != kwargs['nx']:
      raise ValueError('nx ({}) inconsistent with shape {}'.format(nx, shape))
  elif 'shape' in kwargs:
    kwargs['nx'] = kwargs['shape'][0]
  
  else:
    read_sgy._log.debug('No nx nor shape provided. You need it to preserve ' + 
                       '3D array structure ' + 
                       '(nx, ny, nz). Otherwise it is gonna be (ntraces, 1, nsamps).')
  
  read_sgy._log.debug('File to read: ' + fname)
  fname_vtr = fname[:-len('sgy')] + 'vtr'
  
  convert = True
  if exists(fname_vtr):
    read_sgy._log.debug(fname_vtr + ' already exists.')
    if overwrite:
      read_sgy._log.debug('Overwriting ' + fname_vtr)
    else:
      read_sgy._log.debug('Skipping sgy2vtr because overwrite=0.')
      convert = False
  
  if convert:
    sgy2vtr(fname, **kwargs)
  
  A = read_vtr(fname_vtr, **kwargs)
  
  return A
@timer
@logged
def header2csv(fname, keys='all', suffix='_HEAD', overwrite=True, **kwargs):
  """
  Read trace-header values.
  
  Parameters
  ----------
  fname : str
    Name of the SEGY file.
  
  Returns
  -------
  header : pd.DataFrame
  
  Notes
  -----
  SEGY -> SU -> dict. It is quite slow.
  
  """
  from pandas import DataFrame, read_csv
  from fullwavepy.ioapi.su import sugethw
  
  nfname = strip(fname) + suffix + '.csv'
  
  if (not exists(fname)) or overwrite:
    header = DataFrame(sugethw(fname, keys, **kwargs))
    header2csv._log.debug('Overwriting %s' % nfname)
    header.to_csv(nfname, index=False)
  else:
    header2csv._log.debug('Reading %s' % nfname)
    header = read_csv(nfname)
  
  return header
@timer
@logged
def header2json(fname, **kwargs):
  """
  Wrapper around read_header!
  
  SEGY -> read_header -> dict -> JSON
  
  Export the SEG-Y file's header 
  to a JSON file.
  
  Notes
  -----
  Actually not used in the code.
  Handy in ipynb though.
  
  """
  from fullwavepy.generic.parse import strip
  from .json import save_json
  
  nfname = strip(fname) + json_header_suffix
  
  header2json._log.info('Exporting the header of ' + 
                        fname + ' to ' + nfname)
  
  h = read_header(fname, **kwargs)
  save_json(nfname, h, **kwargs)
  return nfname
# -------------------------------------------------------------------------------
# GEO
# -------------------------------------------------------------------------------
@logged
def read_geo(fname, unit='node', **kwargs):
  """
  Read .geo files (Fullwave's format 
  of sources/receivers files).
  
  Parameters
  ----------
  fname : str    
    File name. It should include  extension. 
    It can include path if needed.
  dx : float 
    Size of the grid cell in metres.

  Returns
  -------
  records : dict 
    records[id] = [x, y, z]
  
  Notes
  -----
  CAUTION Metres are just nodes converted with dx, 
  it has nothing to do with experiment coordinates 
  frame (xorigin etc.)
  
  x AND z ARE NOT SWAPPED IN THIS FORMAT 
    (IN CONTRAST TO .pgy)
  '1 +' ENSURES CONSISTENCY WITH .pgy FILES  

  For both btop=0 and btop=-99 the coordinates frame 
  is centered at the free surface which corresponds to
  grid node 0 and 1 respectively.
  That's why if we use btop=0 (as we do)
  free

    
  """
  from .generic import read_txt
  
  read_geo._log.debug('Assuming btop=0, not -99')
  
  content = read_txt(fname, **kwargs)
  header = content[0]
  data = content[1: ]

  # NOTE  
  # A list, not dict as we want to have a data-structure preserving the order of the records!
  # (this will matter)
  records = [] 
  for row in data:
    if unit == 'm':
      try:
        x0 = kwargs['x0']
        y0 = kwargs['y0']
        z0 = kwargs['z0']
      except KeyError as err:
        raise KeyError('For unit={} you need to provide origin of the coordinate frame.'.format(unit), err)
      
      dx = kwargs['dx']

      xyz = [float(row[1])+x0, 
             float(row[2])+y0, 
             float(row[3]) - dx +z0]    # this is for btop=0 only!   
    
    elif unit == 'node':
      dx = kwargs['dx']
      xyz = [1 + float(row[1]) / dx, 
             1 + float(row[2]) / dx, 
             1 + (float(row[3]) - dx) / dx]  # this is for btop=0 only!
    else:
      raise ValueError('Unknown unit: ' + str(unit))  
    
    records.append([int(row[0]), xyz])
  
  read_geo._log.debug('Assuming all keys are integer numbers')
  
  return records
