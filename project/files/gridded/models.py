"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.files.generic import ArrayProjFile
from fullwavepy.ioapi.segy import SgyFile
from fullwavepy.ioapi.fw3d import VtrFile
from fullwavepy.project.lists.basic import ShotFileList, TimestepFileList
from fullwavepy.project.files.gridded.generic import GridFile


# -------------------------------------------------------------------------------


@traced
@logged
class ModelFile(GridFile):
  """
  File containing a property 
  defined on the model grid.
  
  """
  def dupl(self, source, cmd='cp',  *args, **kwargs): 
    """
    Standard file-duplication + resizing 
    to conform to the project-box.
    
    """
    super().dupl(source, cmd=cmd, **kwargs)
    
    if 'geom' in dir(self.proj):
      #self.resize(*args, **kwargs)
      self.__log.warn('Resize disabled until debugged')
    else:
      self.__log.warn('proj.geom not yet set.' + self.fname + 
                      ' will not be resized')

  # -----------------------------------------------------------------------------  
  

# -------------------------------------------------------------------------------


@traced
@logged
class ModelFileSgy(ModelFile, SgyFile):
  """
  File containing a property defined 
  on the model grid.
  
  """

  # -----------------------------------------------------------------------------
  
  def __init__(self, suffix, proj, path, **kwargs):
    """
    """
    self.name = proj.name + '-' + suffix + '.sgy'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)

  # -----------------------------------------------------------------------------
  
  def create(self, *args, **kwargs):
    """
    
    """
    from fullwavepy.ioapi.segy import vtr2sgy
    super().create(*args, **kwargs)
    fvtr = strip(self.fname) + '.vtr'
    vtr2sgy(fvtr, self.proj.dx)

  # -----------------------------------------------------------------------------  

  def resize(self, file_z0=None, **kwargs):
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
    self.__log.warn('Resize disabled until debugged')
    #self.__log.info('Assuming integer box coords, as required by SEGY')
    #box = [int(i) for i in self.proj.box]
    #x1, x2, y1, y2, z1, z2 = box
    #
    #self.__log.warn('\n\n DISABLED BUGGY file_z0 CONVERSION!!!\n\n')
    ##z1 -= file_z0
    ##z2 -= file_z0
    #
    #self.__log.debug('z1={}, z2={}'.format(z1, z2))
    #
    #dt = int(self._gethw('dt', unique_values=True, timer=True, **kwargs)[0])
    #self.__log.debug('Converting z1,z2 into microsec as required by suwind')
    #if dt > 1000:
    #  self.__log.info('Header dt > 1000. Assuming miliseconds or milimetres')
    #  z1 /= 1e3
    #  z2 /= 1e3
    #elif dt <= 1000:
    #  self.__log.info('Header dt <= 1000. Assuming seconds or metres')
    #  z1 /= 1e6
    #  z2 /= 1e6
    #
    #self.__log.debug('z1={}, z2={}'.format(z1, z2))
    #
    #scalco = int(self._gethw('scalco', unique_values=True, timer=True, **kwargs)[0])
    #if scalco < 0: 
    #  scalco = abs(scalco)
    #elif scalco == 0:
    #  scalco = 1
    #else: # scalco > 0 MEANS IT IS USED AS A MULTIPLIER NOT A DIVISOR IN THE HEADER
    #  scalco = 1 / abs(scalco) # WE NEED TO DO THE OPPOSITE THING TO OUR BOX
    #  
    #self.__log.debug('scalco' + str(scalco))
    #x1, x2, y1, y2 = np.array([x1, x2, y1, y2]) * scalco
    #
    #key_x = self.proj.sgy.hw['xmod']
    #key_y = self.proj.sgy.hw['ymod'] 
    #
    #tmp_fname = 'tmp.sgy'
    #cmd = str('segyread tape=' + self.fname + ' | ' +
    #          'suwind key=' + key_x + 
    #          ' min=' + str(x1) + 
    #          ' max=' + str(x2) + ' | ' + 
    #          'suwind key=' + key_y + 
    #          ' min=' + str(y1) + 
    #          ' max=' + str(y2) + ' | ' +  
    #          'suwind tmin=' + str(z1) + ' tmax=' + str(z2) + ' | ' + 
    #          'segyhdrs | ' +
    #          'segywrite tape=' + tmp_fname)
    #
    #self.__log.debug(cmd)
    #o, e = bash(cmd)
    #o, e = bash('mv ' + tmp_fname + ' ' + self.fname) 
      
  # ----------------------------------------------------------------------------- 
  
  def read(self, **kwargs):
    """
    
    Notes
    -----
    Just adding nx info to retain 
    the 3D structure of the model while
    converting from SEG-Y.
    
    FIXME
    It's only needed when using sgy2vtr, I guess.
    
    """
    try:
      kwargs['nx'] = kw('nx', self.proj.dims[0], kwargs)
    except AttributeError:
      self.__log.warn('self.proj.dims not defined') 
    
    self.__log.debug("Set kwargs['nx'] to " + str(kwargs['nx']))
    A = super().read(**kwargs)
    return A
    
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class ModelFileVtr(ModelFile, VtrFile):
  """
  File containing a property defined 
  on the model grid.
  
  """

  # -----------------------------------------------------------------------------  
  
  def __init__(self, suffix, proj, path, **kwargs):
    """
    """
    self.name = proj.name + '-' + suffix + '.vtr'
    self.fname = path + self.name
    super().__init__(proj, path, **kwargs)

  # -----------------------------------------------------------------------------

  def resize(self, x1_file, y1_file, z1_file, dx_file, **kwargs):
      """
      
      """
      A = self.read(**kwargs)
      nx, ny, nz = A.shape
      raise NotImplementedError('resize vtr')
      #x2_file = x1_file + (nx - 1) * dx_file
      #y2_file = y1_file + (ny - 1) * dx_file
      #z2_file = z1_file + (nz - 1) * dx_file
      
      ## CONVERT METRES -> NODES
      #x1, x2, y1, y2, z1, z2 = self.proj.box
      #x1n = int((x1 - x1f) / self.dxs[0])
      #x2n = int((x2 - x1f) / self.dxs[0])
      #y1n = int((y1 - y1f) / self.dxs[1])
      #y2n = int((y2 - y1f) / self.dxs[1])
      #z1n = int((z1 - z1f) / self.dxs[2])
      #z2n = int((z2 - z1f) / self.dxs[2])
      
      ## SLICE THE ARRAY
      #model = np.array(model[x1n : x2n + 1, y1n : y2n + 1, z1n : z2n + 1]) # NOTE: DOUBLE-CHECK + 1!
      #print('model.shape after slicing', model.shape)
      
      #Save_vtr(self.fname[ :-len('.vtr')], model.shape[0], model.shape[1], model.shape[2], model, **kwargs)
      
      ## UPDATE FILE INFO TO PREVENT FURTHER RESIZING WHEN CALLED AGAIN
      #self.x_origin = x1
      #self.y_origin = y1    
      #self.z_origin = z1   
      
  # -----------------------------------------------------------------------------
      

# -------------------------------------------------------------------------------

