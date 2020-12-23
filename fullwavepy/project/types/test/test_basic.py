"""
Unit tests of ..basic module.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import time
from autologging import logged, traced
from inspect import currentframe
from unittest import TestCase, skip

from fullwavepy.config.logging import log_lvl
from fullwavepy.project.types.basic import *

# 40 - suppress warnings and below
log_lvl(40)


@traced
@logged
class GenericSetUp(object):
  """
  """
  def generic_setup(self):
    self.startTime = time.time()
    # MOST GENERIC
    self._set_path()
    self._set_exe()
    self._set_sgy_hw()
    # EXPERIMENT-SPECIFIC, BASE-FILES CAN CHANGE OVER TIME
    self._set_experiment()
    self._set_base() # NOTE experiment -> base
    # PROJECT SPECIFIC
    self._set_geom()
    # WRAP-UP
    self.kwargs = dict(path=self.path, exe=self.exe, sgy_hw=self.sgy_hw,\
                       base=self.base, **self.geom, cat=0, ls=0)     

  # -----------------------------------------------------------------------------

  def _set_path(self):  
    self.path = "/home/kmc3817/projects_datadrive1/code_fullwavepy_tests/"

  # -----------------------------------------------------------------------------

  def _set_experiment(self):
    from fullwavepy.ioapi.proteus import ProteusExperiment
    self.exp = ProteusExperiment()    

  # -----------------------------------------------------------------------------

  def _set_sgy_hw(self):
    self.sgy_hw = {'sid': 'tracf',
                   'rid': 'fldr',
                   'lid': 'ep',
    }    

  # -----------------------------------------------------------------------------

  def _set_geom(self):
    # geom_thin18short
    dt = 0.0025 
    ns = 1500    
    dx = 50      
    x1 = -13400
    x2 = -10550
    y1 = 1500
    y2 = 3000
    z1 = 0 # CUT AT THE SEA SURFACE
    z2 = 500   
    box = [x1, x2, y1, y2, z1, z2]
    self.geom = {'box': box, 'dx': dx, 'ns': ns, 'dt': dt}  

  # -----------------------------------------------------------------------------

  def _set_exe(self):
    self.exe = {'fullwave':       '~/PhD/fullwave3D/rev690/bin/fullwave3D.exe',
       'fullwave_local': '/home/kmc3817/light_PhD/fullwave3D/rev690/bin/fullwave3D.exe',
       'segyprep':       '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16',
       'fsprep':         '/home/kmc3817/light_PhD/fsprep/fsprep',
       'modprep':        '/home/kmc3817/light_PhD/fullwave3D/modprep/modprep.exe',
    }      

  # -----------------------------------------------------------------------------

  def _set_base(self):
    self.base = dict(tvp=self.exp.svp['bh']['18-04-24'],\
                rsg=self.exp.wavelet['19-09-22'],\
                rse=['/home/kmc3817/projects_datadrive1/code_fullwavepy_tests/data/tmp-Synthetic.sgy'],\
                sp=None)    

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class GenericTearDown(object):
  """
  """
  def generic_teardown(self):
    self._print_runtime()

  # -----------------------------------------------------------------------------
  
  def _clean_dirs(self):
    #cmd = 'rm -r %s/%s' % (self.path, self.pname)
    pass

  # ----------------------------------------------------------------------------- 
  
  def _print_runtime(self):
    
    t = time.time() - self.startTime    
    print('Ran in %.3f s' % (t)) 

  # ----------------------------------------------------------------------------- 


# -------------------------------------------------------------------------------

@traced
@logged
class TestProjSynFromBaseLocalRun(GenericSetUp, GenericTearDown, TestCase):
  """
  """
  def setUp(self):
    self.generic_setup() 
    self.kwargs_inp = dict(reciprocity=1, rm=0, cat=0, ls=0,
                         rnf_kwargs=dict(b_abs=5, e_abs=10, cat=0))
    self.kwargs_run = dict()
    self.kwargs_out = dict()

  # -----------------------------------------------------------------------------

  def tearDown(self):
    p = ProjSyn('tmp', info=currentframe().f_code.co_name, **self.kwargs)
    p.inp.prep(**self.kwargs_inp)    
    p.run(**self.kwargs_run)
    p.out.plot(**self.kwargs_out)
    self.generic_teardown()  

  # -----------------------------------------------------------------------------
    
  #@skip('tmp')
  def test_segy_3d(self):
    self.kwargs['io'] = 'sgy'
    self.kwargs['dim'] = '3d'  

  # -----------------------------------------------------------------------------

  # @skip('tmp')  
  def test_fw3d_3d(self):
    self.kwargs['io'] = 'fw3d'
    self.kwargs['dim'] = '3d'  

  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------


@traced
@logged
class TestSyntheticInversion(TestCase): #FIXME: Fix and move to ProjInvSyn?
  @skip("tmp")
  def test_basic_3D_sgy_regular_geom_prep_and_run(self):
    p = ProjSyn(*self.default_args, **self.default_kwargs)
    p.i.rawsign.prep('ricker', fpeak=5)
    p.i.rsg.plot()
    
    p.i.tvp.prep(np.full(p.dims, 2000)) # dims are crucial!
    p.i.tvp.plot()
    cat = 0
    p.i.sp.prep(geometry='regular', souz=p.nx3//2, 
            soux0=p.nx1//2, soudx=1, sounx=1, 
            souy0=p.nx2//2, soudy=1, souny=1,
            recz=p.nx3//2, 
            recx0=p.nx1//2, recnx=1, recdx=1,
            recy0=p.nx2//2, recdy=1, recny=1, cat=cat)
    p.i.sp.read(cat=cat)   
    p.exe['segyprep'] = '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16'
    p.exe['fullwave_local'] = '~/my_phd/fullwave3D/rev688_kmc_LOCAL/bin/fullwave3D.exe'
    p.i.sp.run(cat=cat)
    
    b_abs = 2
    e_abs = 5

    p.i.runfile.prep(b_abs=b_abs, e_abs=e_abs, cat=cat)
    no = 0
    p.o.o.no[no].rm()
    p.o.e.no[no].rm()
    p.i.bash.no[no].prep(cat=cat)
    p.i.bash.no[no].run(cat=cat)

    p.o.syn.plot()
    
    psyn = p
    p = ProjInv('tmp_inv', exe=psyn.exe, **psyn.kwargs)
    # p.i.obs.raw.dupl(psyn.o.syn.fname) # FIXME: TMP
    psyn.i.ose.dupl(psyn.o.syn.fname) # FIXME: UGLY HACK
    kw_filt =  {'pad': 100, 'f1': 2, 'f2': 3, 'f3': 4.5, 'f4': 6.5, 'zerophase': False}
    p.i.prep(psyn, filt_kwargs=kw_filt, mute_kwargs={}, \
      b_abs=b_abs, e_abs=e_abs, \
      blocks=[{'nits': 1, 'freq': 3},
              {'nits': 1, 'freq': 4}])
    no = 0
    p.o.o.no[no].rm()
    p.o.e.no[no].rm()
    p.i.bash.no[no].prep(cat=cat)
    p.i.bash.no[no].run(cat=cat)              
    
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------












@traced
@logged
class TestProjSynOLD(object): #TestCase):
  """
  Priorities: 
  - sgy/vtr
  - 2d/3d
  - immerse
  - wavefields
  - management of files

  Less urgent:
  - anisotropy, Q, elastic
  """
  def setUp(self):
    """
    Common prerequisites for all the methods.
    It will be called as first for each and every test 
    function written in a test class.
    """
    from fullwavepy.config.logging import log_lvl
    log_lvl(40)
    self.startTime = time.time()

    self.pname = 'tmp'
    self.path = './'
    self.default_args = [self.pname]
    # for small ns, creating the wavelet with fpeak not large enough may fail!
    self.default_kwargs = dict(dt=0.001, dx=10, ns=400, dims=(21,11,11))

  # -----------------------------------------------------------------------------
    
  def tearDown(self):
    """
    Common 'post-processing' for all the methods.
    It will be called as first for each and every test 
    function written in a test class.    
    """
    from fullwavepy.generic.system import bash
    cmd = 'rm -r %s/%s' % (self.path, self.pname)
    t = time.time() - self.startTime
    
    # print('Ran %s in %.3f s' % (self.id(), t))
    print('\nRan in %.3f s\n' % (t))
    # bash(cmd)
    # self.__log.warning(cmd)

  # -----------------------------------------------------------------------------
  
  # def test_init_only_pname(self):

  # -----------------------------------------------------------------------------
  
  @skip("tmp")
  def test_basic_init(self):
    raise NotImplementedError('Useless')
    return
    p = ProjSyn(*self.default_args, **self.default_kwargs)

  # -----------------------------------------------------------------------------
  
  @skip("tmp")
  def test_project_attributes(self):
    p = ProjSyn(*self.default_args, **self.default_kwargs)
    attrs = [
        'anisotropy',
        'box',
        'cluster',
        'dim',
        'dims',
        'domain',
        'dt',
        'dtms',
        'dx',
        'env',
        'equation',
        'exe',
        'geom',
        'i',
        'ibfs',
        'info',
        'init_input',
        'init_output',
        'inp',
        'io',
        'kernel',
        'kwargs',
        'ls',
        'meta',
        'name',
        'nn',
        'ns',
        'nx1',
        'nx2',
        'nx3',
        'o',
        'out',
        'parent_dir',
        'path',
        'pbox',
        'plot_input',
        'plot_output',
        'prepare_input',
        'problem',
        'proj',
        'qp',
        'qs',
        'rsync',
        'sgyhw',
        'ttime',
        'units',
    ]
    for attr in attrs:
      answer = hasattr(p, attr)
      self.assertTrue(answer, msg='Missing attribute %s' % attr)

  # -----------------------------------------------------------------------------
  

# -------------------------------------------------------------------------------

