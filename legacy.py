"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for peremoveission writing to k.chrapkiewicz17@imperial.ac.uk.

Notes
-----
Legacy code likely to be chucked away for good.

"""


@timer
@traced
@logged
def read_arrays(*args, **kwargs):
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


# -------------------------------------------------------------------------------


@traced
@logged
def slice_array(A, **kwargs):
  """
  
  Parameters
  ----------
  A : array 
   3D array, although other shapes
   should be handled too (work in progress).
   
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  2D or 1D array (special case).
  
  Notes
  -----
  
  """
  scoord = kw('scoord', 'y', kwargs)
  dims = A.shape

  # 1D ARRAY CONTAINING A SINGLE 1D TIME-SERIES
  if len(dims) == 1:
    slice_array._log.warn('1D array detected. Passing it on intact')
    return A
  
  # 3D ARRAY CONTAINING A SINGLE 1D TIME-SERIES
  if (dims[0] == 1) and (dims[1] == 1):
    return A[0][0]
  
  # SURFACES
  if (dims[-1] == 1):
    return A[..., 0]
  
  if scoord == 'x':
    svalue = kw('svalue', len(A)//2, kwargs)
    A = A[svalue]
  
  elif scoord == 'y':
    svalue = kw('svalue', len(A[0])//2, kwargs)
    A = [i[svalue] for i in A]
  
  elif scoord == 'z':
    svalue = kw('svalue', len(A[0][0])//2, kwargs)
    A = [[j[svalue] for j in i] for i in A]
  
  elif scoord is None:
    slice_array._log.warn('No slicing applied')
  
  else:
    raise ValueError('Wrong slice coord: ' + scoord)
  
  if scoord is not None:
    slice_array._log.info('Sliced at svalue=' + str(svalue) + ' node ' + 
                          'of scoord=' + scoord + ' axis')
  
  A = np.array(A)
  
  return A
  

# -------------------------------------------------------------------------------
