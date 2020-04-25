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
  def test_set_dx_1d(self):
    a = Arr(np.zeros(nx))
    self.assertEqual(len(a.shape), 1)
    # s = 'foobar'
    # self.assertEqual(core(s, ' ', 'bar'), 'foo')
    # self.assertEqual(core(s, 'foo', ' '), 'bar')
    # self.assertRaises(ValueError, core, s, '', 'bar')
    # self.assertRaises(ValueError, core, s, 'foo', '')


# -------------------------------------------------------------------------------