"""
Unit tests of ..basic module.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import time
from autologging import logged, traced
from unittest import TestCase

from ..basic import *


@traced
@logged
class TestProjSyn(TestCase):
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
    from fullwavepy.logging_config import log_lvl
    log_lvl(40)
    self.startTime = time.time()

    self.pname = 'tmp'
    self.path = './'
    self.default_args = [self.pname]
    self.default_kwargs = dict(dt=0.01, dx=10, ns=100, dims=(81, 41, 61))

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
    print('Ran in %.3f s' % (t))
    # bash(cmd)
    # self.__log.warning(cmd)

  # -----------------------------------------------------------------------------
  
  # def test_init_only_pname(self):

  # -----------------------------------------------------------------------------

  def test_basic_init(self):
    p = ProjSyn(*self.default_args, **self.default_kwargs)

  # -----------------------------------------------------------------------------

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

  def test_basic_prep_3D_sgy(self):
    p = ProjSyn(*self.default_args, **self.default_kwargs)
    p.i.rawsign.prep('ricker', fpeak=5)
    p.i.rsg.plot()
    
    p.i.tvp.prep(np.full(p.dims, 2000)) # dims are crucial!
    p.i.tvp.plot()
    
    p.i.sp.prep(geometry='regular', souz=p.nx3//2, 
            soux0=p.nx1//2, soudx=1, sounx=1, 
            souy0=p.nx2//2, soudy=1, souny=1,
            recz=p.nx3//2, 
            recx0=p.nx1//2, recnx=1, recdx=1,
            recy0=p.nx2//2, recdy=1, recny=1)
    p.i.sp.read()   
    p.exe['segyprep'] = '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16'
    p.i.sp.run()
    
  # -----------------------------------------------------------------------------


# -------------------------------------------------------------------------------
