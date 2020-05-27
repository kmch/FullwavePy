"""
Unit tests of ..proteus module.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import time
from unittest import TestCase, skip

from ..proteus import *


class TestProteus(TestCase):
  def setUp(self):
    self.pe = ProteusExperiment()
  def tearDown(self):
    pass
  def test_plots(self):
    self.pe.svp['bh']['18-04-24'].plot()
    self.pe.bathytopo.plot()
