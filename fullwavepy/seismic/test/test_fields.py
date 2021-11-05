from unittest import TestCase
from fullwavepy.seismic.fields import Wavefield

class TestWavefield(TestCase):
  def test_get_file_names(self):
    w = Wavefield()
    w._get_file_names()
    # for fn in w.fnames:
    #   print(fn)
  def test_io(self):
    w = Wavefield()
    fname = 'f01syn/out//f01syn-fw-000400-csref04125-iter00001-taskid00003.vtr'
    assert w.io._extract_srcid(fname) == 4125
    assert w.io._extract_tstep(fname) == 400
    w._get_file_names()
    print(w._get_src_ids(w.fnames))