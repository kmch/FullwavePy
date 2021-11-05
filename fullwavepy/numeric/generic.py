"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.numeric.const import *


@traced
@logged
def round_to_int(x):
  """
  Round to nearest int if it's 
  close enough.
  
  """
  threshold = 1e-6
  
  c = np.ceil(x)
  f = np.floor(x)
  
  return np.where(abs(x - f) < threshold, f, 
                  np.where(abs(x - c) < threshold, c, x))


# -------------------------------------------------------------------------------


@traced
@logged
def deg2rad(angle):
  """
  Convert degrees to radians.
  
  Parameters
  ----------
  angle : float 
    Angle in degrees.
    
  Returns
  -------
  angle_rad : float 
    Angle in radians.
    
  Notes
  -----
  
  """
  return angle * np.pi / 180


# -------------------------------------------------------------------------------


@traced
@logged
def rad2deg(angle):
  """
  Convert radians to degrees.
  
  Parameters
  ----------
  angle : float 
    Angle in radians.
    
  Returns
  -------
  angle_deg : float 
    Angle in degrees.
    
  Notes
  -----
  
  """
  return angle * 180 / np.pi 


# -------------------------------------------------------------------------------


@traced
@logged
def neighs_of_int(a, r):
  """
  Integers surrounding an integer 'a' within radius 'r'.
  
  Examples
  --------
  r=1, a=2   returns (1,3)
  r=2, a=2   returns (0,1,3,4)
  
  """
  return np.arange(a-r, a+r+1) # INCLUDES THE POINT ITSELF (SAME AS FULLWAVE3D's HICKS)
  #return np.array(list(np.arange(a-r, a)) + list(np.arange(a+1, a+r+1))) # EXCLUDES


# -------------------------------------------------------------------------------


@traced
@logged
def neighs_of_float(a, r):
  """
  Integers surrounding a float 'a' within radius 'r'.
  
  Examples
  --------
  r=1, a=2.5   returns (2,3)
  r=2, a=2.5   returns (1,2,3,4)
  
  """
  return np.arange(np.ceil(a)-r, np.ceil(a)+r)


# -------------------------------------------------------------------------------


@traced
@logged
def neighs1d(a, r):
  """
  Nearest nodes in 1D.
  
  a : float
  r : int
  
  Notes
  -----
  For integer 'a', 'a' is included itself.
  This is for ND discrete delta functions, 
  that needs non-zero values.
  
  This would be the version excluding int 'a':
  return np.where(float(a).is_integer(), neighs_of_int(a, r), 
                                         neighs_of_float(a, r))
  In this case we  can make use of np.where as the shape of the output 
  is the same in both situations.
  
  Examples
  --------
  r=1, a=2   returns (1,2,3)
  r=1, a=2.5 returns (2,3)
  
  """
  if float(a).is_integer():
    n = neighs_of_int(a, r)
  else:
    n = neighs_of_float(a, r)
  
  n = np.array([int(i) for i in n])
  neighs1d._log.debug('nodal-neighbours of %s within r=%s are: %s' % (a, r, n))
  return n
  #
  

# -------------------------------------------------------------------------------


@traced
@logged
def decimal(a):
  """
  Decimal (fractional) part of a float a.
  """
  return a - np.floor(a)


# -------------------------------------------------------------------------------


@traced
@logged
def divisors_of(n):
  """
  Return all divisors of n 
  in descending order.
  
  Returns
  -------
  a generator (sth you can
  directly iterate over)
  
  Notes
  -----
  
  OBSOLETE
  1 is excluded from the list
  as trivial and for purpose below.
  
  Used to optimal resources
  to request on a cluster.
  
  """
  yield n
  for i in np.arange(n/2+1, 0, -1):
    if n%i == 0: yield int(i)


# -------------------------------------------------------------------------------


@traced
@logged
def normalize(x, **kwargs):
  """
  Guard against small denominators.
  
  """
  norm = kw('norm', 'rms', kwargs)
  x = np.array(x)
  num_zero = 1e-10
  
  if norm == 'rms':
    norm = rms(x)
  elif norm == 'max':
    norm = np.max(np.abs(x)) # NOTE: abs
  else:
    raise ValueError('Unknown norm: ' + norm)
  
  normalize._log.debug('norm=%s' % norm)

  if norm > num_zero:
    x_normed = x / norm
  else:
    x_normed = np.zeros(x.shape)
  
  return x_normed  

def norm_bulk_max(A, **kwargs):
  return A / np.max(np.abs(A))


# -------------------------------------------------------------------------------


@traced
@logged
def rms(y):
  """
  Root mean square.
  
  Parameters
  ----------
  y : list / array
    A vector of numbers.
  
  Returns
  -------
  rms : float 
    RMS of the y.
  
  Notes
  -----
  rms(y) = (y_i * y_i / n)^(1/2)
  
  It uses only built-in numpy operations
  which ensures the best performance.
  
  """
  y = np.array(y)
  n = len(y)
  return np.sqrt(np.sum(y**2) / n)


# -------------------------------------------------------------------------------

