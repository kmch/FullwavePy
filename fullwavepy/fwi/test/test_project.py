import numpy as np
import os
from unittest import TestCase

import fullwavepy
from fullwavepy.fwi.project import Project

test_dir = os.path.dirname(fullwavepy.fwi.test.__file__)

class TestSynthetic(TestCase):
  def test_init(self):
    p = Project(name='test_init', path=test_dir)
    p.init()
    print(p.params.all)