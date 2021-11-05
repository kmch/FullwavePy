from unittest import TestCase
from fullwavepy.seismic.proteus import PROTEUS

class TestPROTEUS(TestCase):
  def test_read_metadata(self):
    xp = PROTEUS()
    xp.read_metadata()
    xp.read_bathy_topo()
    assert len(xp.pool) == 3