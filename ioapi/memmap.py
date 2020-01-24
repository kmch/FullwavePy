"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced


mmp_dtype = np.float32 # DATA TYPE FOR ARRAYS


# -------------------------------------------------------------------------------


@traced
@logged
def read_mmp(fname, shape, **kwargs):
  """
  Read a memory-mapped file.
  
  shape : tuple
    Shape of the array stored in fname.
    (this is necessary).
    
  Notes
  -----
  See top lines of this module:
  mmp_dtype = np.float32 
  
  """
  from fullwavepy.generic.system import exists
  
  if not exists(fname):
    raise FileNotFoundError(fname) 
  
  fA = np.memmap(fname, dtype=mmp_dtype, shape=shape)
  return fA


# -------------------------------------------------------------------------------


@traced
@logged
def save_mmp(A, fname, **kwargs):
  """
  Save an array to a memory-mapped file.
  
  fname : str 
    Including .mmp extension.
    
  Notes
  -----
  See top lines of this module:
  mmp_dtype = np.float32 
  
  """
  fA = np.memmap(fname, dtype=mmp_dtype, mode='w+', shape=A.shape)
  
  for x in range(A.shape[0]):
    for y in range(A.shape[1]):
        fA[x, y] = A[x][y]


# -------------------------------------------------------------------------------

