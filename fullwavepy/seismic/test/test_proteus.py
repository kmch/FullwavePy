from unittest import TestCase, skip
from fullwavepy.seismic.proteus import PROTEUS, CoordSystemPROTEUS
class TestPROTEUS(TestCase):
  def test_(self):
    xp = PROTEUS()
    # skip, as I moved the heavy model out of the package
    # assert xp.svp.arr.shape == (2481, 861, 131)
    assert xp.wvl.arr.shape == (2000,)
class TestCoordSystemPROTEUS(TestCase):
  def test_local2geogr(self):
    lon, lat, dep = CoordSystemPROTEUS.local2geogr(0,0,0)
    assert lon == 25.3971
    assert lat == 36.404199999999996
    assert dep == 0.0    
    # assert lon[0] == 25.3971
    # assert lon[1] == 'E'
    # assert lat[0] == 36.404199999999996
    # assert lat[1] == 'N'
    # assert dep == 0.0
  def test_rotate_anticlock(self):
    xy = [1,0]
    xy = CoordSystemPROTEUS.rotate_anticlock(xy, 90)
    assert abs(xy[0]) < 1e-8 # numeric zero
    assert xy[1] == 1
