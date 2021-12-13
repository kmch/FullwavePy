import os
from unittest import TestCase, skip
from fullwavepy.seismic.wavefields import Full
module_path = os.path.dirname(__file__)

class TestWavefield(TestCase):
  def test_get_file_names(self):
    w = Wavefield()
    proj_name = 'spike_fwhm1400_ampl050syn'
    path = '%s/example_wavefield/' % (module_path)    
    w._get_file_names(proj_name, path)
    assert len(w.fnames) == 2
  def test_get_files(self):
    w = Wavefield()
    proj_name = 'spike_fwhm1400_ampl050syn'
    path = '%s/example_wavefield/' % (module_path)
    w.get_files(proj_name, path)
    print(w.id)
  def test_io_extracts(self):
    w = Wavefield()
    fname = 'f01syn/out//f01syn-fw-000400-csref04125-iter00001-taskid00003.vtr'
    assert w.io._extract_srcid(fname) == 4125
    assert w.io._extract_tstep(fname) == 400
  