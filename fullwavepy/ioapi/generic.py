"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import exten, strip, kw, del_kw
from fullwavepy.generic.system import exists, bash

from fullwavepy.plot.generic import FilePlotter

@traced
@logged
class FilePreparer(object):
  """
  """
  def prep(self, *args, **kwargs):
    """
    Note, create is a method of child classes only!

    FIXME: add more cases: but maybe only as _ProjFilePreparer(?)
    """
    dupl = kw('dupl', None, kwargs)
    
    if dupl is None:
      self.create(*args, **kwargs)
    else:
      del_kw('dupl', kwargs)
      self.dupl(dupl, *args, **kwargs)

  # -----------------------------------------------------------------------------

  def _from_other_file(self, *args, **kwargs): # FIXME DUPL
    pass

  # -----------------------------------------------------------------------------
  
  def _from_scratch(self, *args, **kwargs): # FIXME: AKA CREATE
    pass

  # -----------------------------------------------------------------------------
  
  def _get(self, *args, **kwargs):
    pass

  # -----------------------------------------------------------------------------  
  
  # FIXME: -> _create
  def create(self, *args, **kwargs):
    #raise NotImplementedError('Define in a child class.') # with this, a lot fails
    pass

  # -----------------------------------------------------------------------------

  def dupl(self, source, cmd='cp', **kwargs):
    """
    Duplicate (cp/mv/ln) a file.
    
    """
    from fullwavepy.generic.parse import exten
    from fullwavepy.generic.system import duplicate
    
    destination = self.fname
    
    if exten(source) != exten(destination):
      raise IOError('File extensions of source and destination '\
                    'must be the same.')    
    
    self.__log.debug('Duplicating %s...' % source)    
    duplicate(source, destination, cmd)

  # -----------------------------------------------------------------------------
@traced
@logged
class File(FilePreparer):
  """
  Generic file.

  Notes
  -----
  Inherits from mix-ins to avoid violation of SRP 
  and god-class appearance, in particular.
  
  """
  def __init__(self, name, path, **kwargs):
    self.name = name
    self.path = path
    self.fname = path + '/' + name
  def rm(self, cmd='trash', backup=True, **kwargs):
    """
    
    """
    if backup:
      bckp = self.fname + '_bckp'
      self.__log.debug('Creating a backup ' + bckp)
      o, e = bash('cp ' + self.fname + ' ' + bckp)
    
    self.__log.info('Removing ' + self.fname + '...')
    o, e = bash(cmd + ' ' + self.fname)
@traced
@logged
class ArrayFile(FilePlotter, File):
  """
  File storing data (model, seismograms, etc.) 
  as a custom Arr3d type of array which is a 
  wrapper around np.ndarray.
  
  """
  def read(self, overwrite=True, **kwargs):
    """
  
    Notes
    -----
    Overwrite=True by default because otherwise plots are not 
    updated even though they are supposed to. 
    They will be correct (updated) only
    if you delete self.array variable, e.g. by restarting the 
    notebook kernel.
    Disable overwrite only for PERFORMANCE (e.g. interactive plot)
    when the array remains unchanged unlike other (e.g. plotting)
    parameters.
    
    """
    from fullwavepy.ndat.arrays import Arr3d
    
    if (not hasattr(self, 'array')) or overwrite:
      self.__log.debug('{}.array does not exist and will be read.'.format(type(self)))
      self.array = Arr3d(self.fname, **kwargs)
    
    # PASS extent (ESSENTIAL FOR PLOTS) FOR PLOTTING FROM FILE TO ARRAY
    if hasattr(self, 'extent'):
      self.__log.debug('Setting extent to %s' % str(self.extent))
      self.array.extent = self.extent 
    return self.array
  # def _plot(self, **kwargs):
  #   """
  #   We SHOULD self.read everytime, 
  #   see its docstring for rationale.
    
  #   """
  #   array = self.read(**kwargs)
  #   return array.plot(**kwargs)
  
  # ----------------------------------------------------------------------------

  def plot_3slices(self, *args, **kwargs):
    array = self.read(**kwargs)
    return array.plot_3slices(*args, **kwargs)

  # ----------------------------------------------------------------------------

  # def plotf(self, *args, **kwargs):
  #   self.array.plotf(*args, **kwargs)

  # ----------------------------------------------------------------------------  

  def scroll(self, **kwargs): # MOVE TO Arr3d?
    """
    
    """
    import matplotlib.pyplot as plt
    from fullwavepy.plot.events import IndexTracker
    
    A = self.read(scoord=None)
    
    fig, ax = plt.subplots(1, 1)
    tracker = IndexTracker(ax, A, **kwargs)
    return fig, ax, tracker

  def scrollall(self, fig, **kwargs):
    """
    
    """
    from fullwavepy.plot.events import IndexTrackerAll
    
    A = self.read(scoord=None)
    
    tracker = IndexTrackerAll(fig, A, **kwargs)
    return tracker
    #return tracker.onscroll

  # -----------------------------------------------------------------------------

  def compare(self, other_file, *args, **kwargs): #gs=None, **kwargs):
    arr1 = self.read(**kwargs)
    arr2 = other_file.read(**kwargs)
    # arr1.compare(arr2, fig, gs, **kwargs)
    arr1.compare(arr2, *args, **kwargs)

  # ----------------------------------------------------------------------------
@traced
@logged
class BinaryFile(File):
  pass
@traced
@logged
class TextFile(File):
  """
  File storing some meta data,
  not directlfullwavepy.plottable as an array.
  
  """
  def read(self, **kwargs):
    """
    
    """
    from fullwavepy.ioapi.generic import read_txt
    return read_txt(self.fname, **kwargs)
  
  # ----------------------------------------------------------------------------

  def cat(self, **kwargs):
    """
    Cat a content of a file (it assumes 
    it is possible to do so, that is 
    it is a an ASCII file.
    
    Notes
    -----
    Should be overwritten by NotIMplementedError 
    in not-ASCII files.
    
    """
    fname = self.fname
    o, e = bash('cat ' + fname)
    print('Content of ', fname, ': ')
    print(o, e)  

  # ----------------------------------------------------------------------------
@traced
@logged
class CsvFile(TextFile):
  @timer
  def read(self, **kwargs):
    from pandas import read_csv
    usecols = kw('usecols', None, kwargs)
    df = read_csv(self.fname, usecols=usecols)
    return df
# ------------------------------------------------------------------------------- 
# Functions
# -------------------------------------------------------------------------------
@traced
@logged
def read_any(fname, overwrite_mmp=False, **kwargs):
  """
  Read file of any format, possibly as a 
  memory-mapped file if both the file is 
  present on disk and the shape of the array 
  is provided. FILE MUST CONTAIN AN ARRAY.
  
  shape : tuple
    Shape of the array stored in fname.
    (this is necessary to read .mmp).
  
  """
  from .memmap import read_mmp, save_mmp
  
  if overwrite_mmp:
    read_any._log.info('Set overwrite_mmp=False for faster i/o!')
  else:
    read_any._log.info('If the array looks corrupted try overwrite_mmp=True.')

  shape = kw('shape', None, kwargs)
  
  fname_mmp = strip(fname) + '.mmp'
  
  if not exists(fname_mmp) or overwrite_mmp:
    read_any._log.debug(fname_mmp+' does not exist. Reading ' + 
                       fname + ' instead...')
    A = read_any_format(fname, **kwargs)
    read_any._log.info('Saving ' + fname_mmp + '...')
    save_mmp(A, fname_mmp)

  elif shape is None:
    read_any._log.debug('File ' + fname_mmp + ' exists, but you' +
                       ' need to provide its shape to read it. Reading ' + 
                       fname + ' instead...')
    A = read_any_format(fname, **kwargs) 
    
  else:
    read_any._log.debug(fname_mmp+' exists and its shape is provided: ' +
                       str(shape))
    A = read_mmp(fname_mmp, **kwargs)    
  
  return A
@traced
@logged
def read_any_format(fname, **kwargs):
  """
  Read an array from a file of
  any format.
  
  """
  from fullwavepy.ioapi.fw3d import read_vtr, read_ttr
  from fullwavepy.ioapi.segy import read_sgy
  from fullwavepy.ioapi.memmap import read_mmp
  
  read_any_format._log.debug('Reading ' + fname + '...')
  ext = exten(fname, **kwargs)
  
  ext = ext.lower() # convert to lower case

  if ext == 'vtr':
    A = read_vtr(fname, **kwargs)
  elif ext == 'ttr':
    A = read_ttr(fname, **kwargs)    
  elif ext == 'sgy' or ext == 'segy':
    A = read_sgy(fname, **kwargs)
  elif ext == 'mmp':
    A = read_mmp(fname, **kwargs)
  elif ext == 'txt':
    c = read_txt(fname, **kwargs)
    A = np.zeros((1,1,len(c)))
    A[0,0,:] = [float(i[0]) for i in c]    
  else:
    raise ValueError('Unknown extension: ' + ext)
  
  return A
@traced
@logged
def read_dict(fname, separator=':', comment_char='!', **kwargs):
  """
  Read a dictionary data-structure from the file
  of the format:
  
  key1 : value1
  key2 : value2
  .
  .
  .
  Lines that not adhere to this format will not be read!
  
  
  Parameters
  ----------
  fname : str 
    Full name of the file including 
    a path if needed.
  separator : str
    Character which separates the key
    from the value for each record (line).
    By default ':'.
  comment_char : str
    Character which begins a comment-line.
    By default '!'.
  Returns
  -------
  dictionary : dict 
    Dictionary to save.
    
  Notes
  -----
  Keys must be a single string without spaces.
    
  """
  read_dict._log.debug('File to read: ' + fname)
  read_dict._log.debug('Will convert all keys into lower case')
  
  # split into keys and values using the separator
  content = read_txt(fname, separator=separator)
  
  params = {}
  for record in content:
    # skip bad/empty? lines
    if len(record) < 2:
      continue
    # remove leading and trailing spaces and remove Caps
    key = record[0].strip().lower()
    val = record[1].strip().lower()
    # skip commented lines
    if key[0] == comment_char:
      continue    
    # ignore trailing comments and strip again (needed)
    val = val.split(comment_char)[0].strip()
    # THIS IS FOR PARSING SKELETON RUNFILE CREATED BY SEGYPREP
    # and should be moved somewhere else
    if key == 'max time':
      key = 'ttime'
      val = str(float(val) / 1000.) # convert ms to s
    # Finally, add to dict
    params[key] = val
  return dict(sorted(params.items()))
@traced
@logged
def read_dict_OLD(fname, **kwargs):
  """
  Read a dictionary data-structure from the file
  of the format:
  
  key1 : value1
  key2 : value2
  .
  .
  .
  Lines that not adhere to this format will not be read!
  
  
  Parameters
  ----------
  fname : str 
    Full name of the file including 
    a path if needed.
  
  Returns
  -------
  dictionary : dict 
    Dictionary to save.
    
  Notes
  -----
  Keys must be a single string without spaces.
    
  Examples
  --------
  
  """
  read_dict._log.debug('File to read: ' + fname)
  read_dict._log.debug('Will convert all keys into lower case')
  
  content = read_txt(fname)
  
  params = {}
  for record in content:
    if (len(record) >= 3) and (record[1] == ':') and (record[0][0] != '!'):
      key = record[0].lower()
      value = record[2]
      params[key] = value
    
    # THIS IS FOR PARSING SKELETON RUNFILE CREATED BY SEGYPREP
    if (record[0] == 'MAX') and (record[1] == 'TIME'):
      key = 'ttime'
      value = str(float(record[3]) / 1000.) # CONVERT ms TO s
      params[key] = value
  
  return params  
@traced
@logged
def save_dict(fname, dictionary, **kwargs):
  """
  Save a dictionary data-structure to a file.
  
  Parameters
  ----------
  fname : str 
    Full name of the file including 
    a path if needed.
  dictionary : dict 
    Dictionary to save.
  
  Returns
  -------
  Saves a file.
  
  Notes
  -----
  The key can be a number and probably
  any other object too.
  
  Examples
  --------
  >>> save_dict('test.txt', {'a' : 4, 'b' : 54})
  >>> !cat test.txt 
  a : 4
  b : 54  
  
  """
  with open(fname, 'w') as f:
    # NOTE WE SORT IT TO MAKE IT EASIER TO READ
    for key, value in sorted(dictionary.items()):
      f.write(str(key) + ' : ' + str(value) + '\n')
@traced
@logged
def read_txt(fname, separator=None, **kwargs):
  """
  Read a text-file's non-zero lines and
  split them using a blank-space separator.
  
  Parameters
  ----------
  fname : str 
    Name of the file including the path
    if needed.
  separator : str
    Character separating different words.
    By default None, which effectively is
    a space (or more spaces).
  
  Returns
  -------
  content : list
    List of split lines (lists of words).
  
  """
  content = []
  
  with open(fname, 'r') as f:   
    for line in f:
      line = line.split(separator)
      if len(line) != 0:
        content.append(line)

  return content
@traced
@logged
def read_txt_raw(fname, **kwargs):
  """
  Similar to read_txt but without
  splitting lines.
  
  Notes
  -----
  Useful in preserving formatting.
  
  """
  content = []
  with open(fname, 'r') as f:
    for line in f:
      content.append(line)

  return content
@traced
@logged
def save_txt(fname, lines, **kwargs):
  with open(fname, 'w') as f:
    for line in lines:
      f.write(line + '\n')

