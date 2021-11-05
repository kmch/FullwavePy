"""
FWI pre-processing codes that transform generic 
data-structures  (e.g. an array of true Vp model) 
into specific input files for the given FWI code.

Motivation
----------
To generate synthetic data one needs: 
a model, a wavelet and an acquisition geometry.
However, different FWI codes require different input 
files. The role of a pre-processor is to prepare these
files based on the above data-structures.

Notes
-----
One pre-processor may serve different FWI codes, and conversely,
one FWI code may be served by different pre-processors.
E.g. Fullwave3D allows for different 'I/Os' (native and SEGY formats)
i.e. sets of input and output files.

"""
class PreCodeFactory:
  @classmethod
  def create(cls, ID, *args, **kwargs):
    if ID == 'sp-316':
      return SegyPrep316(*args, **kwargs)
    else:
      raise ValueError      
class PreCode(ABC):
  def __init__(self, proj):
    self.proj = proj
    self._init_exe()    
  @abstractmethod
  def _init_exe(self):
    pass    
  @abstractmethod
  def prep(self):
    pass
  @abstractmethod
  def prep_n_run(self):
    pass   
  @abstractmethod
  def run(self):
    pass    
class SegyPrep(PreCode):
  """
  Interface with SegyPrep -- a Fortran preprocessor
  for Fullwave3d written by Mike Warner et al.

  """
  def _create_acq(self, acq, **kwargs):
    """
    Create acquisition geometry in a form
    of self.sp_kwargs dictionary. 

    Parameters
    ----------
    acq : string
        Acquisition type.
    
    Returns
    -------
    None

    Raises
    ------
    ValueError
        If acq not in ['regular', 'single']
    """
    self.sp_kws = {}
    unit = kwargs.get('unit', 'nodes')
    self.sp_kws['geometry_in_nodes'] = 1 if unit != 'm' else 0 
    
    if acq == 'regular':
      self._create_acq_regular(**kwargs)
    elif acq == 'single':
      self._create_acq_single(**kwargs)
    else:
      raise ValueError
  def _create_acq_regular(self, **kwargs):
    """
    Very ad-hoc regular-geometry creator.

    Notes
    -----
    It assumes all sources are at z=sz
    and all receivers at z=rz.

    """
    proj = self.proj
    if 'src' in kwargs:
      src = kwargs['src']
      sz = src[-1]
      soux0 = src[0]
      soudx = 1
      sounx = 1
    else:
      assert 'sz' in kwargs
      assert 'sx' in kwargs
      sz = kwargs['sz']
      assert isinstance(sz, float) or isinstance(sz, int)
      sx = kwargs['sx']
      assert isinstance(sx, list) or isinstance(sx, np.ndarray)
      sdx = np.array(np.array(sx[1:] - sx[:-1]))
      assert np.all((sdx - sdx[::-1]) < 1e-6)
      sdx = sdx[0]      
      souz = sz
      soux0= sx[0]
      soudx = sdx
      sounx = len(sx)
      # else:
      #   src = np.array(proj.dims) / 2
      #   if 'sz' in kwargs:
      #     sz = kwargs['sz']
      #     assert isinstance(sz, float) or isinstance(sz, int)
      #     src[-1] = sz
    souy0 = self.proj.dims[1]/2
    soudy = 1
    souny = 1      

    # receivers
    rz = kwargs['rz']
    assert isinstance(rz, float) or isinstance(rz, int)
    rx = kwargs['rx']
    assert isinstance(rx, list) or isinstance(rx, np.ndarray)
    rdx = np.array(np.array(rx[1:] - rx[:-1]))
    assert np.all((rdx - rdx[::-1]) < 1e-6)
    rdx = rdx[0]

    self.sp_kws = dict(self.sp_kws, ztype='d', cat=0, geometry='regular',\
      souz=souz, soux0=soux0, soudx=soudx, sounx=sounx,\
      souy0=souy0, soudy=soudy, souny=souny,\
      recz=rz, recx0=rx[0], recdx=rdx, recnx=len(rx),\
      recy0=self.proj.dims[1]/2, recdy=1, recny=1) 
    
    if not self.sp_kws['geometry_in_nodes']:
      for k in ['souz', 'soux0', 'soudx', 'souy0', 'soudy',
            'recz', 'recx0', 'recdx', 'recy0', 'recdy']:
        self.sp_kws[k] *= self.proj.dx
  def _create_acq_single(self, src, rec, unit='m', **kwargs):
    kws = dict(geometry='regular',
           soudx=1, sounx=1, soudy=1, souny=1,
           recdx=1, recnx=1, recdy=1, recny=1,
           ztype='d', cat=0)
    kws = dict(kws, soux0=src[0], souy0=src[1], souz=src[2],
            recx0=rec[0], recy0=rec[1], recz=rec[2])
    self.sp_kws = dict(self.sp_kws, **kws)
  def prep(self, acq='regular', **kwargs):
    self._create_acq(acq, **kwargs)
    self.proj.i.sp.prep(**self.sp_kws)
  def prep_n_run(self, **kwargs):
    self.prep(**kwargs)
    self.run(cat=0)
  def run(self, **kwargs):
    self.proj.i.sp.run(**kwargs)
class SegyPrep310(SegyPrep):
  pass
class SegyPrep316(SegyPrep310):
  def _init_exe(self):
    self.path = '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16'
    self.proj.exe['segyprep'] = self.path
class SegyPrep324(SegyPrep316):
  pass
class SegyPrep325(SegyPrep324):
  pass
