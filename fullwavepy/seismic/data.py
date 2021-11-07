"""
Handling seismic data which can be tricky for two reasons:
- it is usually stored as heavy binaries,
- it has associated metadata that need to be pass around.

Copywright (c) 2019- Kajetan Chrapkiewicz.
Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.
"""
from abc import ABC, abstractmethod
from autologging import logged #, traced
import numpy as np
import pandas as pd # only in _get_stations_within_extent

from arrau.a2d import Arr2d
from arrau.api.io import extent2str
from fullwavepy.ioapi.generic import save_txt  # for mute
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.ioapi.su import suwind
from nsh.generic import ShellFactory
from nsh.utils import core, exists, extract_file_name, \
  get_file_names, extract_ext, extract_path, strip
# ----------------------------------------------------------------------------
# To get rid of.
# ----------------------------------------------------------------------------
import matplotlib.pyplot as plt # FIXME
from fullwavepy.generic.parse import kwarg_over_attr, path_leave # FIXME
from fullwavepy.generic.system import get_files # FIXME
from fullwavepy.ndat.arrays import Arr3d # FIXME
from fullwavepy.ioapi.generic import ArrayFile # FIXME
from fullwavepy.plot.generic import * # FIXME
# ----------------------------------------------------------------------------

class Dat:
  def __init__(self, file_object, dt):
    self.dt = dt
    self.file = file_object
  def plot(self, *args, **kwargs):
    if not hasattr(self, 'arr'):
      self.read(*args, **kwargs)
    kwargs['aspect'] = kwargs.get('aspect', 'auto')
    kwargs['cmap'] = kwargs.get('cmap', 'seismic')
    kwargs['center_cmap'] = kwargs.get('center_cmap', 'True')
    kwargs['label'] = kwargs.get('label', 'amplitude')
    self.arr.plot(*args, **kwargs)
    ax = plt.gca()
    ttl = ''
    if 'station' in kwargs:
      ttl += 'Station %s. ' % kwargs['station']
    if 'line' in kwargs:
      ttl += 'Line %s. ' % kwargs['line']
    title = kwargs.get('title', ttl)
    ax.set_title(title)
    ax.set_xlabel('trace no.')
    ax.set_ylabel('time, s')
    ax.invert_yaxis()
  def read(self, station=None, line=None, normalise='max', **kwargs):
    if station is None and line is None:
      kws = {}
    elif line is None:
      kws = dict(win=dict(tracf=[station]))
    else:
      kws = dict(win=dict(tracf=[station], ep=[line]))
    arr = self.file.read(**kws)
    if len(arr.shape) == 3:
      assert arr.shape[1] == 1
      arr = arr[:,0,:]
    extent = self._get_extent(arr)
    self.arr = Arr2d(arr, extent=extent)
    self.arr.normalise(normalise)
    return self.arr
  def _get_extent(self, arr):
    ntr, ns = arr.shape
    tmax = self.dt * ns
    return [[1,ntr],[0,tmax]]
class DataIOFactory:
  subclasses = {}
  @classmethod
  def create(cls, ID, *args, **kwargs):
    import fullwavepy.seismic.proteus
    import fullwavepy.seismic.rainbow
    if ID not in cls.subclasses:
      raise ValueError('Wrong ID {}'.format(ID))
    return cls.subclasses[ID](*args, **kwargs)
  @classmethod
  def register_subclass(cls, ID):
    def decorator(subclass):
      cls.subclasses[ID] = subclass
      return subclass
    return decorator
class DataIO(ABC):
  def __init__(self):
    self._set_file_class()
    self._set_pattern()
  @abstractmethod
  def _get_srcids_within_extent(self, extent):
    pass  
  @abstractmethod
  def _set_file_class(self):
    pass    
  @abstractmethod
  def _set_pattern(self):
    pass
class DataMuter(ABC):
  @abstractmethod
  def mute(self):
    pass
class DataMuterSUSGY(DataMuter):
  def mute(self, fname, tmute, twin=1, ntaper=100, \
    fname_out=None, shell='linux'):
    """
    fname : str
        Name of the file with data mute.
    fname_out : str
    tmute : list
        List of the top-mute time-picks for each trace,
        in seconds. 
    """
    shell = ShellFactory.create(shell)
    nmute, xmute, tmute2 = self._prep(tmute, twin)
    self._save_mute_files(xmute, tmute, tmute2, shell)
    self._set_fname_out(fname, fname_out)
    cmd =  'segyread tape={} | '.format(fname)
    # top mute
    cmd += 'sumute key=tracr nmute={nmute} mode=0 '.format(nmute=nmute)
    cmd += 'ntaper={ntaper} '.format(ntaper=ntaper)
    cmd += 'xfile={xb} tfile={tb} | '.format(xb='xmute.bin',tb='tmute.bin')
    # bottom mute
    cmd += 'sumute key=tracr nmute={nmute} mode=1 '.format(nmute=nmute)
    cmd += 'ntaper={ntaper} '.format(ntaper=ntaper)
    cmd += 'xfile={xb} tfile={tb} | '.format(xb='xmute.bin',tb='tmute2.bin')
    cmd += 'segyhdrs | segywrite tape={fnout}'.format(fnout=self.fname_out)
    self.cmd = cmd
    o, e = shell.run(cmd)
  def _correct_tmute2_for_bad_picks(self, tmute, tmute2):
    # THIS LOOP ENSURES WE MUTE THE WHOLE TRACE 
    # IF THE FIRST ARRIVAL IS LATER THEN TOTALTIME (AS INDICATE BY pick==0)
    # WITHOUT THIS, IT WAS BUGGY!    
    for i, pick in enumerate(tmute):
      if pick == 0:
        tmute2[i] = 0 
    return tmute2 
  def _prep(self, tmute, twin):
    nmute = len(tmute)
    xmute = range(1, nmute+1)
    tmute2 = np.array(tmute) + twin # List of the bottom-mute time-picks
    tmute2 = self._correct_tmute2_for_bad_picks(tmute, tmute2)
    return nmute, xmute, tmute2
  def _set_fname_out(self, fname, fname_out):
    assert fname != fname_out
    if fname_out is None:
      fname_out = '%s_muted.%s' % (strip(fname), extract_ext(fname))
    self.fname_out = fname_out
  def _save_mute_files(self, xmute, tmute, tmute2, shell):
    xmute = [str(i) for i in xmute]
    tmute = [str(i) for i in tmute]
    tmute2 = [str(i) for i in tmute2]
    datas = [xmute, tmute, tmute2]
    prefixes = ['xmute', 'tmute', 'tmute2']
    for data, prefix in zip(datas, prefixes):
      file_txt = prefix + '.txt'
      file_bin = prefix + '.bin'
      save_txt(file_txt, data)
      cmd = 'a2b < {} n1=1 > {}'.format(file_txt, file_bin)
      o, e = shell.run(cmd)
class DataSet:
  def __init__(self, path, io='proteus_hy', regex=None, \
    shell='linux', get_files=True):
    """
    Parameters
    ----------
    path : str
        Path to directory with data files.
    io : str
    regex : str
        Regular-expression pattern matching
        the data file names. E.g. '*.sgy'.
        If None, it will be read from self.io
    shell : str
        By default 'linux'.
    get_files : bool
        Run _get_files if True. This may be slow
        hence optional.
    """
    self.files = {}
    self.io = DataIOFactory.create(io)
    self.path = path
    self._init_shell(shell) # used in `extract`
    self._set_regex(regex)
    self._get_file_names()
    self._get_sids()
    if get_files: 
      self._get_files()
  def extract(self, extent, name='extracted',\
     path='./', stations=None, exclude=[]):
    """
    Extract a subset of data recorded
    by the `stations` constrained by the `extent` of 
    the spatial domain.

    Parameters
    ----------
    extent : list / array
        Spatial domain (2d) constraining the 
        data subset: [[x1,x2],[y1,y2]]
    name : str
        Name of the extracted subset. It will 
        be used as a name of directory storing
        the files.
    path : str
        Path to save the files with extracted
        data to.
    stations : See `get_files`.
    exclude : See `get_files`.

    Returns
    -------
    DataSet
    """
    # prep for saving
    extent_str = extent2str(extent)
    regex = '*%s*' % extent_str
    [[x1, x2], [y1, y2]] = extent
    window =  dict(sx={'min': x1, 'max': x2},
                   sy={'min': y1, 'max': y2})    
    save_dir = self._create_save_dir(path, name)
    # select station-gathers
    files = self._get_files_within_extent(extent, stations, exclude)
    # iterate over gathers
    for f in files:
      # prep the file name
      core = f.name[:-len('.sgy')]
      new_fname = '{save_dir}/{core}_{extent_str}.sgy'.format(
        save_dir=save_dir, core=core, extent_str=extent_str)
      # extract traces
      suwind(f.fname, new_fname, window)
    return DataSet(save_dir, regex=regex , shell=self.shell)
  def _create_save_dir(self, path, name):
    """
    Notes
    -----
    This is actually taken care of by shell.mkdir
    # if not exists(full_path):
    #   if not exists(path):
    #     raise FileNotFoundError(path)       
    #   print('Creating %s directory in %s') % (name, path)
    """
    save_dir = '%s/%s' % (path, name)
    self.shell.mkdir(save_dir)
    return save_dir
  def _get_file_names(self):
    """
    Get file names with/without path (fnames/names)
    in the path matching the regex pattern.
    """
    self.fnames = get_file_names(self.path, self.regex)
    self.names = [extract_file_name(i) for i in self.fnames]
  def _get_sids(self):  
    self.sids =  [self.io._extract_srcid(i) for i in self.names]
    self.sids = sorted(self.sids)
    return self.sids
  def _get_stations(self, stations):
    if stations is None:
      stations = self.sids
    else:
      pass
    return stations  
  def _get_stations_within_extent(self, extent, stations=None):
    stations_in = self.io._get_srcids_within_extent(extent)
    if stations is None:
      stations = stations_in
    else:
      stations = [i for i in stations if i in stations_in]
    return stations
  def _get_files(self, stations=None, exclude=[]):
    """
    Get file-instances of `stations` excluding
    the ones in `exclude` list.

    Parameters
    ----------
    stations : list
         by default None.
    exclude : list, optional
    
    Returns
    -------
    dict

    """    
    self.files = {}
    stations = self._get_stations(stations)
    for name in self.names:
      sid = self.io._extract_srcid(name)
      if (sid not in stations) or (sid in exclude):
        continue
      self.files[sid] = self.io.DataFile(name, self.path)
    self.files = dict(sorted(self.files.items()))
    self._set_files_alias()
    return self.files
  def _get_files_within_extent(self, extent, stations=None, exclude=[]):
    """
    Returns
    -------
    List
        List of file objects.
    """
    stations = self._get_stations_within_extent(extent, stations)
    files = self._get_files(stations, exclude).values()
    return files
  def _init_shell(self, shell):
    if isinstance(shell, str):
      self.shell = ShellFactory.create(shell)
    else:
      self.shell = shell
  def _set_files_alias(self):
    self.id = self.files
  def _set_regex(self, regex):
    if regex is None:
      self.regex = self.io.pattern
    else:
      self.regex = regex
# ----------------------------------------------------------------------------
# Old-ish
# ----------------------------------------------------------------------------
@logged
class DataFile(ArrayFile):
  """
  Generic seismic-data file. It is meant to 
  handle arbitrarily large files through `window` 
  and `split` methods that are I/O specific
  and need to be implement in child classes.

  It should generalize excellent implementation
  of split etc. from DumpCompareFile.
  
  no read is defined here!
  """
  def plotf(self, fig, tylim=None, fylim=None, *args, **kwargs): # FIXME: merge with plot.misc time_freq
    """
    both time and freq

    """
    if fig is None:
      figure(16,8)
    
    plt.subplot(121)
    self.plot(*args, **kwargs)
    plt.ylim(tylim)
    
    plt.subplot(122)
    kwargs['dt'] = kwarg_over_attr('dt', kwargs, self) # needed in DFT
    kwargs = dict(kwargs, spect='ampl', cmap='hot', center_cmap=False) # DFT, different cmap
    self.plot(*args, **kwargs)
    plt.ylim(fylim)
  def _cast(self, arr, **kwargs):
    self.array = Data(arr)
    return self.array  
  def _read_shot_gather(self, sid, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
  def _read_shot_line(self, sid, lid, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
  def _split(self, *args, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
  def _window(self, *args, **kwargs):
    raise NotImplementedError('Overwrite in a child class')
class DataFileSgy(DataFile, SgyFile):
  def read(self, **kwargs):
    arr = SgyFile.read(self, **kwargs)
    self.array = self._cast(arr, **kwargs)
    return self.array
@logged
class Data(Arr3d):
  """
  Seismic data.

  """
  def compare(self, *args, **kwargs):
    kwargs['cmap'] = kwargs.get('cmap', 'seismic') #'twilight_shifted'
    kwargs['center_cmap'] = kwargs.get('center_cmap', True)
    kwargs['xlabels'] = self._xlabels_from_header(**kwargs)    
    super().compare(*args, **kwargs)
  def interleave(self, othe, **kwargs):
    return super().interleave(othe, slice_at='y', node=0, **kwargs)
  def plot(self, *args, **kwargs):
    kwargs['cmap'] = kwargs.get('cmap', 'seismic') #'twilight_shifted'
    kwargs['center_cmap'] = kwargs.get('center_cmap', True)
    kwargs['xlabels'] = self._xlabels_from_header(**kwargs)
    super().plot(*args, **kwargs)
  def plot_phase_DEBUG(self, freq, ph_type, **kwargs):# sync it better with dumpcomp
    """
    ph_type : syn/obs/dif
    """
    if not hasattr(self, 'head'):
      self.__log.warning('self.head not found. Returning')
      return
  
    plt.scatter(self.head.sx, self.head.sy, vmin=-np.pi, vmax=+np.pi,
                c=self.head['phase %s (%s Hz)' % (ph_type, freq)])
  def _xlabels_from_header(self, xlabels_hw='fldr', **kwargs):
    if hasattr(self, 'head'):
      self.__log.info('There is a header associated with this data.' +\
        '%s keyword will be used as xlabels.' % xlabels_hw)
      xlabels = self.head[xlabels_hw] # one per xtick
    # xlabel = xlabels_hw # one per xaxis
      return xlabels
    else:
      return None
class ObsData(Data):
  """
  Observed (recorded in the field) seismic data.
  """
  pass
class SynData(Data):
  """
  Synthetic (calculated) seismic data.
  
  """
  pass
