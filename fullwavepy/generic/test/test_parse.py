"""
Unit tests of ..parse module.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
from unittest import TestCase
from fullwavepy.generic.parse import *


#class TestExten(TestCase)

#class TestPathLeave(TestCase)

#class TestStrip(TestCase)


# -------------------------------------------------------------------------------


class TestCore(TestCase):
  """
  """
  def test_onesided(self):
    s = 'foobar'
    self.assertEqual(core(s, ' ', 'bar'), 'foo')
    self.assertEqual(core(s, 'foo', ' '), 'bar')
    self.assertRaises(ValueError, core, s, '', 'bar')
    self.assertRaises(ValueError, core, s, 'foo', '')

  def test_empty_string(self):
    s = ''
    self.assertRaises(ValueError, core, s, ' ', ' ')
  
  def test_homog_strings(self):
    s = 'aaaaaaaaaaaaaaaa'
    self.assertRaises(ValueError, core, s, '', '')

  def test_strings_with_spaces(self):
    s = 'hello world!'
    self.assertEqual(core(s, 'h', '!'), 'ello world')
    
  def test_disjoint_strings(self):
    s = 'bcd'
    prefix = 'a'
    suffix = 'e'
    self.assertEqual(core(s, prefix, suffix), s)


# -------------------------------------------------------------------------------

#class TestDelKw(TestCase)    
