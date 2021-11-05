"""
Unit tests of ..arrays module.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
from unittest import TestCase
from fullwavepy.ndat.arrays import *


class TestArr(TestCase):
  """
  """
  def test_default_extent(self):
    a = Arr1d(np.zeros((2)))
    self.assertEqual(a.extent, [[0, 1]])
    a = Arr2d(np.zeros((2,3)))
    self.assertEqual(a.extent, [[0, 1], [0, 2]])
    a = Arr3d(np.zeros((2,3,4)))
    self.assertEqual(a.extent, [[0, 1], [0, 2], [0, 3]])

  def test_set_dx_1d(self):
    a = Arr(np.zeros(10))
    self.assertEqual(len(a.shape), 1)
    # s = 'foobar'
    # self.assertEqual(core(s, ' ', 'bar'), 'foo')
    # self.assertEqual(core(s, 'foo', ' '), 'bar')
    # self.assertRaises(ValueError, core, s, '', 'bar')
    # self.assertRaises(ValueError, core, s, 'foo', '')


# -------------------------------------------------------------------------------

class TestArr2d(TestCase):
  def setUp(self):
    pass
  
  # -----------------------------------------------------------------------------
  
  def tearDown(self):
    pass

  # -----------------------------------------------------------------------------
  
  def test_plot(self):
    Arr2d(np.random.rand(2,2)).plot()
    

# ------------------------------------------------------------------------------