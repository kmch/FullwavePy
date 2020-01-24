"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.parse import strip


# ------------------------------------------------------------------------------- 


@traced
@logged
def sgy2su(fname, **kwargs):
  """
  
  """
  if not exists(fname):
    raise FileNotFoundError(fname)  
  
  fsu = strip(fname) + '.su'
  cmd = str('segyread tape=' + fname +
            ' > ' + fsu)
  
  o, e = bash(cmd)
  

# -------------------------------------------------------------------------------


@traced
@logged
def array2su(fname, A, dt, **kwargs):
  """
  Useful for su_decon.
  
  """  
  from .segy import array2sgy
  
  fsgy = strip(fname) + '.sgy'
  array2sgy(fsgy, A, dt, **kwargs)
  sgy2su(fsgy, **kwargs)
 
 
# -------------------------------------------------------------------------------


@traced
@logged
def get_ntraces(fname, **kwargs):
  """
  Get no. of all traces present in the 
  SEG-Y file.
  
  """
  if not exists(fname):
    raise FileNotFoundError(fname) 
  
  cmd = str('segyread tape=' + fname + ' | ' + 
            "surange | awk 'NR==1 {print $1}'") # SINGLE ' MATTERS
  o, e = bash(cmd)
  ntraces = int(o)
  
  return ntraces


# -------------------------------------------------------------------------------


@traced
@logged
def get_keywords(fname, **kwargs):
  """
  Get all non-zero keywords present in the 
  SEG-Y file.
  
  """
  if not exists(fname):
    raise FileNotFoundError(fname) 
  
  cmd = str("segyread tape=" + fname + " | " + 
            "surange | " + 
            "sed '/^$/q' | " +
            "awk 'NR>1 {print $1}'")
  o, e = bash(cmd)
  keys = o.split(None)
  
  return keys


# -------------------------------------------------------------------------------


@timer
@traced
@logged
def sugethw(fname, key, unique_values=False, **kwargs):
  """
  
  unique_values : bool
    If false, values for each trace are output.
    If true, only non-repeating ones.
    
  """
  from fullwavepy.generic.parse import str2float 
  
  if not exists(fname):
    raise FileNotFoundError(fname) 
  
  sugethw._log.info('Getting ' + key + ' values from ' + fname + ' header')
  
  cmd = str("segyread tape=" + fname + " | " + 
            "sugethw " + key + " | " +
            "grep " + key + " | " +
            "sed 's@" + key + "=@@'")

  o, e = bash(cmd)
  hw_vals = o.split(None)
  hw_vals = str2float(hw_vals, **kwargs)
  
  if unique_values: # REMOVE DUPLICATES
    hw_vals = sorted(list(set(hw_vals)))
  
  # SEGY HEADER CAN'T STORE FLOATS SO IT'S ALWAYS INT
  hw_vals = [int(i) for i in hw_vals]
  
  return hw_vals  


# -------------------------------------------------------------------------------


@traced
@logged
def sushw(fname, key, value, **kwargs):
  """
  Set header value.
  
  Notes
  -----
  We write to a tmp file, otherwise
  the file gets corrupted!
  
  """
  if not exists(fname):
    raise FileNotFoundError(fname) 
  
  tmp_fname = fname + '.tmp'
  
  cmd = str('segyread tape=' + fname + ' | ' + 
            'sushw key=' + key + ' a=' + str(value) + ' | ' +
            'segyhdrs | segyclean | ' + 
            'segywrite tape=' + tmp_fname)
  
  o, e  = bash(cmd)
  #print(o,e)
  o, e = bash('mv ' + tmp_fname + ' ' + fname) 

  
# -------------------------------------------------------------------------------  
  

@traced
@logged
def suwind(fname, nfname, key, vmin, vmax, **kwargs):
  """
  
  nfname : str 
    Output file.
  vmin : float
    Min. value of key.
  vmax : float 
    Max.
  
  """
  if not exists(fname):
    raise FileNotFoundError(fname) 
  
  cmd = str('segyread tape=' + fname + ' | ' + 
            'suwind key=' + key + 
            ' min=' + str(vmin) + 
            ' max=' + str(vmax) + ' | ' 
            'segyhdrs | ' +
            'segywrite tape=' + nfname)
  
  o, e = bash(cmd)
  #print(o, e)


# -------------------------------------------------------------------------------

