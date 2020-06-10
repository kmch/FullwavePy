"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, strip, path_extract, path_leave
from fullwavepy.generic.system import exists
from fullwavepy.ioapi.generic import ArrayFile

# def read_sources_fw3d():
#   pass





# -------------------------------------------------------------------------------
# VTR 
# -------------------------------------------------------------------------------


@traced
@logged
class VtrFile(ArrayFile):
  """
  
  """
  def create(self, array, **kwargs):
    save_vtr(array , strip(self.fname) + '.vtr')

  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------


@traced
@logged
def save_vtr(data, fname, **kwargs):
  """
  Save data defined on a 3D 
  grid as a .vtr file.
  
  Parameters
  ----------
  data : array / list 
    Array of a shape (nx1,nx2,nx3).
  fname : str
    Name of the file including a path
    if not current directory.
    
  Returns
  -------
  Nothing.
  
  Notes
  -----
  The order of indices is the same as in Fullwave3D.   
  It might not be compatibile with Wavefor2D.
  
  """  
  from scipy.io import FortranFile
  
  save_vtr._log.debug('File to save: ' + fname)
    
  data = np.array(data)
  nx1, nx2, nx3 = data.shape
  
  f = FortranFile(fname, 'w')
  
  if nx2 == 1:
    ndims = 2
    f.write_record(np.array([1, ndims, 0], dtype=np.int32))
    f.write_record(np.array([nx3, nx1], dtype=np.int32))
  elif nx2 > 1:
    ndims = 3
    f.write_record(np.array([1, ndims, 0], dtype=np.int32))  
    f.write_record(np.array([nx3, nx2, nx1], dtype=np.int32)) #NOT COMPATIBLE WITH LLUIS CODE
  else:
    raise ValueError('Dimension nx2 = ' + str(nx2) + ' is < 1.')
    
  for x in range(int(nx1)):
    for y in range(int(nx2)):
      trace = data[x, y, :]
      f.write_record(np.array(trace, dtype=np.float32))
  
  f.close()  


# -------------------------------------------------------------------------------


@timer
@traced
@logged
def read_vtr(fname, **kwargs):
  """
  Read a .vtr file.
  
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
    Gridded 3D data structure of shape that should be 
    consistent with nx1, nx2, nx3 (see above).
  
  Notes
  -----
  .vtr is a Fullwave3D's native binary format. 
  
  """
  from scipy.io import FortranFile
  
  if not exists(fname):
    raise FileNotFoundError(fname)
  
  f = FortranFile(fname, 'r')
  
  #-------------------------------------------------------------
  # READ THE HEADER
  #-------------------------------------------------------------
  ncomponents, ndims, dummy_zero = f.read_ints(dtype = np.int32)
  
  if ndims == 3: 
    nx3, nx2, nx1 = f.read_ints(dtype = np.int32)
  elif ndims == 2:
    nx3, nx1 = f.read_ints(dtype = np.int32)
    nx2 = 1
  elif ndims == 1: #TEST IT!
    nx3, nx1 = f.read_ints(dtype = np.int32)
    nx2 = 1
  else:
    raise IOError('Wrong ndims: ' + str(ndims))

  #-------------------------------------------------------------
  # READ THE DATA
  #-------------------------------------------------------------
  data = []
  for x in range(int(nx1)):
    xslices = []
    for y in range(int(nx2)):
      trace = f.read_reals(dtype = np.float32)
      xslices.append(trace)
    data.append(xslices)    
     
  f.close()
  
  data = np.array(data)
  
  #-------------------------------------------------------------
  # DEBUG
  #-------------------------------------------------------------
  nx1_dat = data.shape[0]
  nx2_dat = data.shape[1]
  nx3_dat = data.shape[2]
  
  if (nx1_dat != nx1 or nx2_dat != nx2 or nx3_dat != nx3):
    read_vtr._log.warning('Shape of the data-array inconsistent with the header:')
    read_vtr._log.warning('nx1: data ' + str(nx1_dat) + ', header ' + str(nx1) + '\n')
    read_vtr._log.warning('nx2: data ' + str(nx2_dat) + ', header ' + str(nx2) + '\n')
    read_vtr._log.warning('nx3: data ' + str(nx3_dat) + ', header ' + str(nx3) + '\n')
  
  read_vtr._log.debug('min value: %s' % np.min(data))
  read_vtr._log.debug('max value: %s' % np.max(data))
  read_vtr._log.debug('nx1, nx2, nx3: %s %s %s' % (nx1, nx2, nx3))
  
  return data


# -------------------------------------------------------------------------------
# TTR 
# -------------------------------------------------------------------------------


@traced
@logged
class TtrFile(ArrayFile):
  """
  
  """
  def _get_sr_coords(self, **kwargs):
    raise NotImplementedError('It must be overwritten by a child class.')
    
  def plot(self, **kwargs):
    kwargs['center_cmap'] = kw('center_cmap', True, kwargs)
    return super().plot(**kwargs)


# -------------------------------------------------------------------------------


#save_ttr


# -------------------------------------------------------------------------------


@traced
@logged
def read_ttr(fname, **kwargs):
  """
  
  It is not as easy to implement as vtr 
  because FortranFile package fails to 
  read records of mixed data-types 
  (see ttr header).
  
  """
  from fullwavepy.generic.system import bash, exists
  from fullwavepy.generic.parse import strip
  
  if not exists(fname):
    raise FileNotFoundError(fname)
  

  
  # T
  name = path_leave(fname)
  path = path_extract(fname)

  len_max = 40
  if len(name) > len_max:
    read_ttr._log.warning('File name %s is longer than %s and may not be suitable for ttr2vtr code yet' % \
      (name, len_max))


  # FIXME rename
  cmd = str('convert_ttr2vtr_IMPROVED ' + name)

  o, e = bash(cmd, path=path)
  read_ttr._log.debug(o + e)
  
  fvtr = strip(fname) + '.vtr'
  
  if not exists(fvtr):
    raise IOError('Conversion ttr2vtr raised: ' + e)
  
  A = read_vtr(fvtr)

  return A


# -------------------------------------------------------------------------------
# PGY
# -------------------------------------------------------------------------------  


@traced
@logged
def save_pgy(fname, dictionary, dims=None, **kwargs):  #FIXME
  """
  Save Fullwave's .pgy geometry files.
  
  Parameters
  ----------
  fname : str    
    File name. It should include  extension. 
    It can include path if needed.

  Returns
  -------
  
  Notes
  -----
  header 
  id z y x
  
  """
  n = len(dictionary) # NO. OF POINTS
  sn   = '{:10}'.format(n)
  
  if dims is None:
    snx1 = ' notset '
    snx2 = ' notset '
    snx3 = ' notset '
  else:
    nx1, nx2, nx3 = dims
    snx3 = '{:15}'.format(nx3)
    snx2 = '{:15}'.format(nx2)
    snx1 = '{:15}'.format(nx1)
  
  with open(fname, 'w') as f:
    f.write(sn + snx3 + snx2 + snx1 + '\n')
    for key, value in sorted(dictionary.items()):
      xyz = [value[2], value[1], value[0]] # NOTE THE ORDER
      sz   = '{:15.8f}'.format(xyz[2])
      sy   = '{:15.8f}'.format(xyz[1])
      sx   = '{:15.8f}'.format(xyz[0])  
      ##print si, sz, sy, sx
      f.write(str(key) + sz + sy + sx + '\n')
      #i += 1


# -------------------------------------------------------------------------------


@traced
@logged
def read_pgy(fname, **kwargs):
  """
  Read Fullwave's .pgy geometry files.
  
  Parameters
  ----------
  fname : str    
    File name. It should include  extension. 
    It can include path if needed.

  Returns
  -------
  nx1, nx2, nx3 : int
    Dimension of the grid as in the geom.
    files.
  
  Notes
  -----
  Format of a pgy file (x AND z ARE SWAPPED compared to geo!)
  header line: Total no. items | nz | ny | nx
  data lines: item id | z [nodes] | y [nodes] | x [nodes]
  
  """
  from fullwavepy.ioapi.generic import read_txt
  
  if not exists(fname):
    raise FileNotFoundError(fname)  
  
  content = read_txt(fname)
  header = content[0]
  data = content[1: ]
  nx3, nx2, nx1 = [int(float(i)) for i in header[1: ]]
  
  records = [] # can't use dict as we want ordered, see read_geo
  for row in data:
    records.append([int(row[0]), np.array([float(row[3]), float(row[2]), float(row[1])])])
  
  return records
  

# -------------------------------------------------------------------------------  


