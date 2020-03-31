# -----------------------------------------------------------------------------
# SEG-Y MAPPING
# -----------------------------------------------------------------------------
sgy_hw = {'sid': 'tracf',
          'rid': 'fldr',
          'lid': 'ep',
          }

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------
from fullwavepy.generic.system import get_files
path = '/home/kmc3817/heavy_PhD/'
path_dataobs = path + 'DATA/Santorini_2015/seismic/OBS/segy_local_coords/'
path_datalan = path + 'DATA/Santorini_2015/seismic/land/Santorini/segy_local_coords/'
dataobs_hy = get_files(path_dataobs, '*4.sgy')
dataobs_vz = get_files(path_dataobs, '*1.sgy')
datalan_vz = get_files(path_datalan, '*1.sgy')

metadataobs = path_dataobs + 'metadata.csv'
metadatalan = path_datalan + 'metadata.csv'
metadata = path + 'metadata/proteus_metadata.csv' # CONCAT. OF ABOVE
md = pd.read_csv(metadata)

startvp_jm = path + 'start_mods/jm_inversecheck-StartVp.sgy'
startvp_bh = path + 'start_mods/Ben_whole_model_18-04-24_sea-clipped.sgy'
startvp_bm = path + 'start_mods/Brennah_whole_model_19-10-30_sea-clipped.sgy'

rawsign = path + 'wavelets/wavelet_19-09-22.sgy'
topography = path + 'surfaces/bathy_x_-8e4_8e4_y_-4e4_4e4_cell_50.vtr'

exe = {'fullwave':       '~/PhD/fullwave3D/rev690/bin/fullwave3D.exe',
       'fullwave_local': '/home/kmc3817/light_PhD/fullwave3D/rev690/bin/fullwave3D.exe',
       'segyprep':       '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16',
       'fsprep':         '/home/kmc3817/light_PhD/fsprep/fsprep',
       'modprep':        '/home/kmc3817/light_PhD/fullwave3D/modprep/modprep.exe',
      }


# -----------------------------------------------------------------------------
# ARRAYS
# -----------------------------------------------------------------------------
from fullwavepy.ioapi.generic import read_any
topo = Arr3d(read_any(topography, shape=(3201,1601,1)))
#topo = topo[...,0]