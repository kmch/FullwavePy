from unittest import TestCase, skip
<<<<<<< HEAD
from fullwavepy.seismic.data import DataSet, DataFileSgy, DataMuterSUSGY
=======
from fullwavepy.seismic.data import DataSet, DataFileSgy
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08

path = '/home/kmc3817/heavy_PhD/DATA/Santorini_2015/seismic/OBS/segy_local_coords/'

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
<<<<<<< HEAD

class TestDataMuter(TestCase):
  def test_data_muter_susgy(self):
    dm = DataMuterSUSGY()
    dm.mute('fname_in.sgy', tmute=[1,2,4], twin=2, fname_out= 'fname_out.sgy')
    cmd = 'segyread tape=fname_in.sgy | sumute key=tracr nmute=3 mode=0 ntaper=100 xfile=xmute.bin tfile=tmute.bin | sumute key=tracr nmute=3 mode=1 ntaper=100 xfile=xmute.bin tfile=tmute2.bin | segyhdrs | segywrite tape=fname_out.sgy'
    assert dm.cmd == cmd
=======
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
