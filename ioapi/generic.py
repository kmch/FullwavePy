"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import exten, strip, kw
from fullwavepy.generic.system import exists, bash


# -------------------------------------------------------------------------------
# CLASSES
# -------------------------------------------------------------------------------


@traced
@logged
class File(object):
  def __init__(self, name, path, **kwargs):
    self.name = name
    self.path = path
    self.fname = path + '/' + name
  #def __init__(self, fname, **kwargs):
  #  self.fname = fname


# -------------------------------------------------------------------------------


@traced
@logged
class AsciiFile(File):
  """
  File storing some meta data,
  not directlfullwavepy.plottable as an array.
  
  """
  
  # -----------------------------------------------------------------------------
  
  def read(self, **kwargs):
    """
    
    """
    from fullwavepy.ioapi.generic import read_txt
    return read_txt(self.fname, **kwargs)
  
  # -----------------------------------------------------------------------------

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

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class CsvFile(AsciiFile):
  @timer
  def read(self, **kwargs):
    from pandas import read_csv
    usecols = kw('usecols', None, kwargs)
    df = read_csv(self.fname, usecols=usecols)
    return df


# -------------------------------------------------------------------------------


@traced
@logged
class BinaryFile(File):
  pass
  

# -------------------------------------------------------------------------------


@traced
@logged
class ArrayFile(File):
  """
  File storing data (model, seismograms, etc.) 
  as a custom Arr3d type of array which is a 
  wrapper around np.ndarray.
  
  """
  
  # -----------------------------------------------------------------------------    

  def read(self, fname=None, overwrite=True, **kwargs):
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
    if (not hasattr(self, 'array')) or overwrite:
      from fullwavepy.generic.array import Arr3d
      self.__log.warn('{}.array does not exist and will be read.'.format(type(self)))
      kwargs['scoord'] = kw('scoord', None, kwargs)
      if fname is None:
        fname = self.fname
      self.array = Arr3d(fname, **kwargs)
    return self.array
  
  # -----------------------------------------------------------------------------

  def plot(self, **kwargs):
    """
    We SHOULD self.read everytime, see its docstring for rationale.
    
    """
    array = self.read(**kwargs)
    array.plot(**kwargs)
    
  # -----------------------------------------------------------------------------  

  def scroll(self, **kwargs):
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


# ------------------------------------------------------------------------------- 
# FUNCTIONS
# -------------------------------------------------------------------------------


@traced
@logged
def read_any(fname, overwrite_mmp=False, **kwargs):
  """
  Read file of any format, possibly as a 
  memory-mapped file if both the file is 
  present on disk and the shape of the array 
  is provided.
  
  shape : tuple
    Shape of the array stored in fname.
    (this is necessary to read .mmp).
  
  """
  from .memmap import read_mmp, save_mmp
  
  shape = kw('shape', None, kwargs)
  
  fname_mmp = strip(fname) + '.mmp'
  
  if not exists(fname_mmp) or overwrite_mmp:
    read_any._log.warn(fname_mmp+' does not exist. Reading ' + 
                       fname + ' instead...')
    A = read_any_format(fname, **kwargs)
    read_any._log.info('Saving ' + fname_mmp + '...')
    save_mmp(A, fname_mmp)

  elif shape is None:
    read_any._log.warn('File ' + fname_mmp + ' exists, but you' +
                       ' need to provide its shape to read it. Reading ' + 
                       fname + ' instead...')
    A = read_any_format(fname, **kwargs) 
    
  else:
    read_any._log.debug(fname_mmp+' exists and its shape is provided: ' +
                       str(shape))
    A = read_mmp(fname_mmp, **kwargs)    

  return A
  

# -------------------------------------------------------------------------------


@traced
@logged
def read_any_format(fname, **kwargs):
  """
  Read an array from a file of
  any format.
  
  """
  from .fw3d import read_vtr, read_ttr
  from .segy import read_sgy
  
  read_any_format._log.debug('Reading ' + fname + '...')
  ext = exten(fname, **kwargs)
  
  if ext == 'vtr':
    A = read_vtr(fname, **kwargs)
  elif ext == 'ttr':
    A = read_ttr(fname, **kwargs)    
  elif ext == 'sgy':
    A = read_sgy(fname, **kwargs)
  elif ext == 'txt':
    c = read_txt(fname, **kwargs)
    A = np.zeros((1,1,len(c)))
    A[0,0,:] = [float(i[0]) for i in c]    
  else:
    raise ValueError('Unknown extension: ' + ext)
  
  return A


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


@traced
@logged
def read_dict(fname, **kwargs):
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


# -------------------------------------------------------------------------------


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


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


@traced
@logged
def read_txt(fname, **kwargs):
  """
  Read a text-file's non-zero lines and
  split them using a blank-space separator.
  
  Parameters
  ----------
  fname : str 
    Name of the file including the path
    if needed.
  
  Returns
  -------
  content : list
    List of split lines (lists of words).
  
  """
  content = []
  
  with open(fname, 'r') as f:   
    for line in f:
      line = line.split(None)
      if len(line) != 0:
        content.append(line)

  return content


# -------------------------------------------------------------------------------


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
  

# -------------------------------------------------------------------------------


@traced
@logged
def save_txt(fname, lines, **kwargs):
  with open(fname, 'w') as f:
    for line in lines:
      f.write(line + '\n')


# -------------------------------------------------------------------------------

















@timer
@traced
@logged
def read_arrays(*args, **kwargs): #FIXME: DEL
  """
  Read data either from arrays
  or files. Both types can appear
  in *args.
  
  Parameters
  ----------
  *args : see below 
    List of string/arrays.
    If string, it is assumed to 
    stand for a file name (incl.
    the path if file is outside './'),
    otherwise it must be an array.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  arrays : list
    List of read arrays.
  
  Notes
  -----
  You can mix files and arrays defined ad hoc
  in the notebook which may come in handy.
  
  """
  from fullwavepy.generic.array import slice_array, modify_array
  
  arrays = []
  for arg in args:
    if isinstance(arg, str):
      A = read_any(arg, **kwargs)
    elif type(arg) == type(np.array([])) or type(arg) == np.memmap:
      A = arg
    else:
      raise TypeError('Arguments need to be either ' + 
                      'file-names or arrays or np.memmap.')
      
    A = slice_array(A, **kwargs)  
    A = modify_array(A, **kwargs)
    read_arrays._log.debug('Min of the array after modifs: ' + str(np.min(A)))
    read_arrays._log.debug('Max of the array after modifs: ' + str(np.max(A)))
    arrays.append(A)
  return arrays