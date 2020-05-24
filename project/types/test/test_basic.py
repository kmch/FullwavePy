"""
Unit tests of ..basic module.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
from autologging import logged, traced
from unittest import TestCase
from ..basic import *


# -------------------------------------------------------------------------------


@traced
@logged
class TestProjSyn(TestCase):
  """
  """
  def setUp(self):
    """
    Common prerequisites for all the methods.
    It will be called as first for each and every test 
    function written in a test class.
    """
    from fullwavepy.logging_config import log_lvl
    log_lvl(30)
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
    # bash(cmd)
    # self.__log.warning(cmd)

  # -----------------------------------------------------------------------------
  
  # def test_basic_init(self):
  #   p = ProjSyn(*self.default_args, **self.default_kwargs)

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


# -------------------------------------------------------------------------------