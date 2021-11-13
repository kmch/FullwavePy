import numpy as np
from unittest import TestCase, skip
from arrau.a2d import Arr2d
from fullwavepy.seismic.data import Dat, DataSet, DataFileSgy, DataMuterSUSGY

path = '/home/kmc3817/heavy_PhD/DATA/Santorini_2015/seismic/OBS/segy_local_coords/'

class TestDat(TestCase):
  def test_init(self):
    dat = Dat(.001)
  def test_read(self):
    dat = Dat(.001, arr=Arr2d(np.zeros((2,2))))
    dat.read()
  def test_interlace(self):
    d1 = Dat(.001, arr=Arr2d(np.zeros((2,2))))
    d2 = Dat(.001, arr=Arr2d(np.ones((2,2))))
    ic = d1.interlace(d2, chunk_size=1)
    assert np.all(ic == [[0., 0.], [1., 1.]])
class TestDataSet(TestCase):
  @skip
  def test_extract(self):
    ds = DataSet(path)
    extent = [[8000.0, 25000.0], [-3000.0, 15000.0]]
    exclude = [i for i in list(range(101,200)) if i != 139]
    pth = './'
    ds.extract(extent, path=pth, exclude=exclude)
  def test_get_files(self):
    ds = DataSet(path)
    # exclude all except 101
    exclude = list(range(102,200))
    files = ds._get_files(exclude=exclude)
    assert len(files) == 1
    assert isinstance(files[101], DataFileSgy)  
  def test_get_files_within_extent(self):
    ds = DataSet(path)
    extent = [[8000.0, 25000.0], [-3000.0, 15000.0]]
    exclude = [139]
    files = ds._get_files_within_extent(extent, exclude=exclude)
    assert len(files) == 15
  def test_get_stations_within_extent(self):
    ds = DataSet(path)
    extent = [[8000.0, 25000.0], [-3000.0, 15000.0]]
    sids = ds._get_stations_within_extent(extent)
    assert len(sids) == 16
class TestDataMuter(TestCase):
  def test_data_muter_susgy(self):
    dm = DataMuterSUSGY()
    dm.mute('fname_in.sgy', tmute=[1,2,4], twin=2, fname_out= 'fname_out.sgy')
    cmd = 'segyread tape=fname_in.sgy | sumute key=tracr nmute=3 mode=0 ntaper=100 xfile=xmute.bin tfile=tmute.bin | sumute key=tracr nmute=3 mode=1 ntaper=100 xfile=xmute.bin tfile=tmute2.bin | segyhdrs | segywrite tape=fname_out.sgy'
    assert dm.cmd == cmd
