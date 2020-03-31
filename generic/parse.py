"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced


# -------------------------------------------------------------------------------


#@traced
@logged
def exten(fname, **kwargs):
  """
  Extract the extension from the 
  file name.
  
  Parameters
  ----------
  fname : str 
    It can include path. It shouldn't include 
    more than one full stop that exactly separates
    the extension from the proceeding part of the 
    name.
  
  Returns
  -------
  ext : str 
    extension (without the full stop).
    
  Examples
  --------
  >>> from fullwavepy.generic.parse import exten
  >>> exten('aa..cc')  
  'cc'
  
  
  """
  split = [_f for _f in fname.split('.') if _f]
  
  if len(split) >= 2:
    ext = split[-1]
    if len(split) > 2:
      exten._log.warn(fname + ' contains more than 1 full stop ' +
                      'separated by other chars. ' + 
                      'Taking the last bit as an extension.')
  else:
    raise IOError('Failed to split ' + fname + 
                  ' using a full stop.')

  return ext


# -------------------------------------------------------------------------------


#@traced
@logged
def path_leave(path, **kwargs):
  """
  Extract a file name from a whole 
  path. 
  
  Parameters
  ----------
  path : str 
    Path including the file name.
    
  Returns
  -------
  file_name : str
    Extracted file name.
  
  Notes
  -----
  Understand return. 
  
  """
  import ntpath
  head, tail = ntpath.split(path)
  
  return tail or ntpath.basename(head)


# -------------------------------------------------------------------------------


#@traced
@logged
def strip(fname, **kwargs):
  """
  Strip the extension off the file name.
  
  Parameters
  ----------
  fname : str 
    String of the form:
    <string>.<extension>
      
  Returns
  -------
  nfname: str 
    Stripped string.
  
  Notes
  -----
  It assumes the string:
    <string_without_dots>.<extension>
  
  Examples
  --------
  >>> from fullwavepy.generic.parse import strip
  >>> strip('aa.bb.cc')
  WARNING: fullwavepy.generic.parse.exten: aa.bb.cc contains > 1 full stop. Taking last bit as an extension.
  'aa.bb'  

  >>> from fullwavepy.generic.parse import strip
  >>> strip('aa..cc')
  'aa..cc' 
  
  """
  ext = exten(fname)
  suffix = '.' + ext
  
  del_kw('suffix', kwargs)
  nfname = core(fname, ' ', suffix, **kwargs)
  
  return nfname


# -------------------------------------------------------------------------------


#@traced
@logged
def core(string, prefix, suffix, **kwargs): #FIXME: DEL??
  """
  
  
  Get a core of the string i.e. strip 
  both its prefix and suffix.
  
  Parameters
  ----------
  string : str 
    String to strip.
  prefix : str 
    Prefix contained in the string.
  suffix : str 
    Suffix contained in the string.
  
  Returns
  -------
  core : str 
    Stripped string.
  
  Notes
  -----
  For a blank prefix/suffix use space (' ', not '').

  Syntax 'list(filter...' is necessary in Python 3 since the 
  filter returns an iterator.
  
  """
  if string == '':
    raise ValueError("string == ''")
  if prefix == '':
    raise ValueError("prefix == ''")
  if suffix == '':
    raise ValueError("suffix == ''")  
  if len(set(string)) == 1:
    raise ValueError('Homogeneous string (a single char repeated).')
  if len(prefix) + len(suffix) > len(string):
    raise ValueError('len(prefix) + len(suffix) > len(string)')

    
  head = list(filter(None, string.split(suffix)))[0]
  core = list(filter(None, head.split(prefix)))[0]
  
  return core 


# -------------------------------------------------------------------------------


#@traced
@logged
def del_kw(kwarg_name, kwargs, warning=False):
  """
  Delete a key from a dictionary.
  
  Parameters
  ----------
  kwarg_name : str
    Name of the argument.
  kwargs : dict 
    Dictionary of all keyword 
    arguments.
  
  Returns
  -------
  kwargs : dict
  
  Notes
  -----
  Name of the function is as short as possible
  because it's frequently typed.
  
  We can't not include kwargs as an argument.
  
  """  
  try:
    del kwargs[kwarg_name]
  except KeyError:
    if warning:
      del_kw._log.warn('Value of ' + kwarg_name + ' was not specified.')
  
  return kwargs


# -------------------------------------------------------------------------------


#@traced
@logged
def kw(kwarg_name, default_value, kwargs):
  """
  FIXME? replace with 
  >>> kwargs.get(key, default_value)
  
  Read a value of a kwarg_name from
  kwargs dictionary or assign a default_value
  if it's missing.
  
  Parameters
  ----------
  kwarg_name : str
    Name of the argument.
  default_value : any 
    The value which will be assigned 
    if kwarg_name is not found 
    in kwargs.
  kwargs : dict 
    Dictionary of all keyword 
    arguments.
  
  Returns
  -------
  value : any 
    Value assigned to kwarg_name.
  
  Notes
  -----
  It should print the caller's name... #FIXME
  
  Don't decorate with @traced (too many messages).
  
  Name of the function is as short as possible
  because it's frequently typed.
  
  We can't not include kwargs as an argument.
  
  Examples
  --------
  >>> from fullwavepy.generic.parse import kw
  >>> kwargs = {'my_key' : 2}
  >>> kw('my_key', 1, kwargs, warning=1)  
  2
  >>> from fullwavepy.generic.parse import kw
  >>> kwargs = {}
  >>> kw('my_key', 1, kwargs, warning=1)  
  fullwavepy.generic.parse WARNING Value of my_key not specified. Setting default value 1
  
  """ 
  try:
    value = kwargs[kwarg_name]
  except KeyError:
    value = default_value 
    kw._log.debug('Set ' + kwarg_name + ' to default: ' + str(default_value))
  
  return value


# -------------------------------------------------------------------------------


#@traced
@logged
def extend_fname(fname, suffices, **kwargs): #FIXME: DEL??
  """
  Extend a name with suffices:
  _id_value.ext for each (id, pair)
  
  suffices : list
  
  """
  nfname = strip(fname)
  #for i, v in zip(ids, values):
  for suffix in suffices:
    nfname = nfname + '_' + suffix[0] + str(suffix[1])
    
  nfname = nfname + '.' + exten(fname)
  return nfname
    

# -------------------------------------------------------------------------------


#@traced
#@logged
#def parse_multiple_sources(params, sources, **kwargs):
#  """
#  
#  params : list 
#    List of string being names of the 
#    params to search for in sources.
#  sources : list
#    List of objects containing attributes
#    fname and params.
#  
#  """
#  parse_multiple_sources._log.info('Searching for params in ' +
#                  ', '.join([i.fname for i in sources]))
#  
#  # SEARCHING FOR PARAMS IN VARIOUS PLACES
#  for param in params:
#    val = None
#    for source in sources:
#      if val is None:
#        val = kw(param, None, source.params)
#        setattr(proj, param, val)
#



