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
path = '/home/kmc3817/heavy_PhD/'
truevp = path + 'start_mods/jm_inversecheck-StartVp.sgy'
rawsign = path + 'wavelets/wavelet_19-09-22.sgy'
topography = path + 'surfaces/bathy_x_-8e4_8e4_y_-4e4_4e4_cell_50.vtr'
sources_csv = path + 'sources/all-Sources.csv'
receivers_csv = path + 'receivers/all-Receivers.csv'
data_path = path + 'DATA/Santorini_2015/seismic/OBS/segy_local_coords/'
from fullwavepy.generic.system import get_files
files_HY = get_files(data_path, '*4.sgy')
files_VZ = get_files(data_path, '*1.sgy')

exe = {'fullwave':       '~/PhD/fullwave3D/rev690/bin/fullwave3D.exe',
       'fullwave_local': '/home/kmc3817/light_PhD/fullwave3D/rev690/bin/fullwave3D.exe',
       'segyprep':       '/home/kmc3817/light_PhD/fullwave3D/segyprep_v3.16/bin/segyprep_v3.16',
       'fsprep':         '/home/kmc3817/light_PhD/fsprep/fsprep',
       'modprep':        '/home/kmc3817/light_PhD/fullwave3D/modprep/modprep.exe',
      }