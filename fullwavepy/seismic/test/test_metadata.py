from unittest import TestCase
from fullwavepy.seismic.metadata import Box3d

class TestBox3d(TestCase):
  def test_init(self):
    b = Box3d(0, 1, 0, 1, 0, 1)