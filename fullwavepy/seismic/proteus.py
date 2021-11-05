"""
This module defines all objects specific to the PROTEUS experiment.
The rationale is to back up and conveniently exchange (meta)data  
between different Jupyter notebooks.

(c) 2021- Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import os  # only for prot_path
import pandas as pd

<<<<<<< HEAD
from arrau.a1d import Arr1d
from arrau.a3d import Arr3d
import fullwavepy # only for prot_path
from fullwavepy.seismic.data import DataFileSgy, DataIO, DataIOFactory
from fullwavepy.seismic.exp import Experiment
from fullwavepy.seismic.misc import Box3d
from fullwavepy.seismic.models import Mod
from fullwavepy.seismic.plots import PlotExp
from fullwavepy.seismic.surfs import SeaLand
from fullwavepy.seismic.wavelets import SourceWavelet

from nsh.generic import ShellFactory
=======
import fullwavepy # only for prot_path
from fullwavepy.seismic.data import DataFileSgy, DataIO, DataIOFactory
from fullwavepy.seismic.exp import Experiment
from fullwavepy.seismic.plots import PlotExp
from fullwavepy.seismic.surfs import SeaLand

>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
path_prot = '%s/proteus/' % os.path.dirname(fullwavepy.seismic.__file__)

class PROTEUS(Experiment,PlotExp):
  """
  Handle for various (meta)data
  associated with the PROTEUS experiment
  (see Hooft et al., 2017).
  """
  def __init__(self, **kwargs):
    self.all_not_read = True
<<<<<<< HEAD
    self._init_boxes()
    self._init_start_vp()
    self._init_wavelet()
=======
    self._init_base_files(**kwargs)
    self._init_start_vp()
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
  def read_all(self):
    self.read_bathy_topo()
    self.read_metadata()
    self.all_not_read = False
  def read_bathy_topo(self):
    name = 'bathy_x-8e4_8e4_y_-4e4_4e4_res50_3201x1601.mmp'
    fname = '%s/%s' % (path_prot, name)
    arr = np.memmap(fname, dtype=np.float32, shape=(3201,1601))
    self.bt = SeaLand(arr, extent=[[-8e4,8e4],[-4e4,4e4]])
    return self.bt
  def read_metadata(self):
    self.recs = pd.read_csv('%s/recs.csv' % path_prot)
    self.srcs = pd.read_csv('%s/srcs.csv' % path_prot)
    self.pool = self._get_instrument_pools()
  # -----------------------------------------------------------------------------
  def _get_instrument_pools(self):
    df = self.recs
    sio = df[df.pool == 'SIO']
    who = df[df.pool == 'WHOI']
    lan = df[df.pool == 'land']
    pool = {'sio':sio,'who':who,'lan':lan}
    return pool
<<<<<<< HEAD
  def _init_boxes(self):
    self.box = {
      'kol1': Box3d(8000.0, 25000.0, -3000.0, 15000.0, 0, 4000.0)
    }
=======
  def _init_base_files(self, **kwargs):
    # from fullwavepy.ioapi.proteus import ProteusExperiment
    # self.pro = ProteusExperiment() # FIXME
    # self.met = kwargs.get('met_base', self.pro.md)
    # self.obs = kwargs.get('obs_base', self.pro.dataset)
    # self.wvl = kwargs.get('wvl_base', self.pro.wavelet['19-09-22'])
    # self.tpg = kwargs.get('tpg_base', self.pro.bathytopo)
    # self.svp = kwargs.get('vp_base', self.pro.svp['bh']['18-04-24'])
    pass
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
  def _init_paths(self):
    base_path = '/home/kmc3817/heavy_PhD/'
    self.path = {
      'data':       base_path + 'DATA/Santorini_2015/',
      'metadata':   base_path + 'metadata/',
      'start_mods': base_path + 'start_mods/',
      'surfaces':   base_path + 'surfaces/',
      'wavelets':   base_path + 'wavelets/',
    }
  def _init_start_vp(self):
    """
<<<<<<< HEAD
    Inits a `self.svp_all` attribute - a nested dictionary 
    containing starting models derived 
    by different authors and processed in 
    different ways, defined over different domains.
    Also sets a default - `self.svp`.
    """
    fname = 'vp_Heath2019_x-6e4_64e3_y-14e3_29e3_z-15e2_5e3_2481x861x131.mmp'
    fname = '%s/%s' % (path_prot, fname)
    vp = np.memmap(fname, dtype=np.float32, shape=(2481,861,131))
    x1 = -6.0e4
    x2 = +6.4e4
    y1 = -1.4e4
    y2 = +2.9e4 
    z1 = -1.5e3 
    z2 = +5.0e3      
    vp = Arr3d(vp, extent=[[x1, x2], [y1, y2], [z1, z2]])
    self.svp_all = {
      'heath19': Mod(vp, param='vp'),
    }
    self.svp = self.svp_all['heath19']
=======
    Inits a `self.svp` attribute which
    is a multi-nested dictionary 
    containing starting models derived 
    by different authors and processed in 
    different ways, constrained by different 
    boxes.
    """
    self.svp = {
      'bh': {
        '18-04-24': {
          'processed_v1': {
            'full': 'Ben_whole_model_18-04-24.sgy',
            'kol1': 'kol.sgy',
            'kam1': 'kam.sgy'
          },
          'processed_v2': {
            'full': None,
          }
        }
      },
      'bm': {}
    }    
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
  def _init_segy_mapping(self):
    self.sgyhw = {
      'sid': 'tracf',
      'rid': 'fldr',
      'lid': 'ep',
    }     
<<<<<<< HEAD
  def _init_wavelet(self):
    fname = 'wavelet_19-09-22_ns2000_dt25ms.mmp'
    fname = '%s/%s' % (path_prot, fname)
    ns = 2000
    dt = 0.0025 # s
    wvl = np.memmap(fname, dtype=np.float32, shape=(2000,))
    wvl = Arr1d(wvl, extent=[[0,ns*dt]])
    self.wvl_all = {
      '19-09-22': SourceWavelet(wvl)
    }
    self.wvl = self.wvl_all['19-09-22']
class CoordSystemPROTEUS:
  @classmethod
  def geogr2local(cls):
    """
    !cat {fname} | proj +proj=tmerc +lat_0={lat0} +lon_0={lon0} +ellps=WGS84 +units=km +datum=WGS84 
    """
    pass 
  @classmethod
  def rotation_matrix(cls, alpha):
    """
    Rotates points in the xy plane counterclockwise through an angle
    alpha (in degrees) with respect to the x axis about the
    origin of a two-dimensional Cartesian coordinate system.


    """
    alpha = alpha * np.pi / 180.
    R = [[+np.cos(alpha), -np.sin(alpha)],\
         [+np.sin(alpha), +np.cos(alpha)]]
    R = np.array(R)
    return R
  @classmethod
  def rotate_anticlock(cls, xy, alpha):
    """
    MP script uses this matrix with NEGATIVE alpha
    (i.e. -25.5 deg) to convert from latlon to local.
    So to convert from local to latlon, we will use POSITIVE alpha.
    """    
    xy = np.array(xy)
    return cls.rotation_matrix(alpha).dot(xy)
  @classmethod
  def local2geogr(cls, x, y, z, lon_0=25.3971, lat_0=36.4042, angle=25.5, \
    proj='tmerc', ellps='WGS84', shell='linux'):
    x, y = cls.rotate_anticlock([x,y], angle) # z is not affected
    shell = ShellFactory.create(shell)
    cmd =  'echo {x} {y} {z} | '.format(x=x,y=y,z=z)
    cmd += 'proj -I +proj=%s ' % proj
    cmd += '+lat_0=%s +lon_0=%s ' % (lat_0, lon_0)
    cmd += '+ellps=%s +units=m +datum=%s' % (ellps, ellps) 
    o, e = shell.run(cmd)
    lon, lat, dep = o.split(None)
    lon = cls._parse_proj_lonlat(lon)
    lat = cls._parse_proj_lonlat(lat)
    dep = float(dep)
    return lon, lat, dep
  @classmethod
  def _parse_proj_lonlat(cls, lonlat):
    deg, rest = lonlat.split('d')
    min, rest = rest.split("'")
    sec, NSWE = rest.split("\"")
    # print(deg, min, sec)
    deg = float(deg)
    min = float(min) / 60.
    sec = float(sec) / 3600.
    deg_float = deg + min + sec
    return deg_float # NOTE discard NSWE, it's always NE in our case
=======
class CoordSystem:
  def latlon2local(ix,iy,nhl,delim):
    """
    A function to automate coordinates transformation - 
    Michele Paulatto - 11/05/2017
    This is the transformation used by Ben for Santorini
    Must provide projection parameters (hard-wired), input and output files
     and column 
    indexes correspoding to longitude and latitude (assumed to be in decimal 
    degrees)
    Also provide number of header lines to skip (nhl) and the text delimiter.
    If nhl is set to -1 it will try to read variable names from the header. 
    This may not work if the names contain characters that are not allowed in 
    the Matlab 
    database.
    This function will convert the coordinates and add them as extra columns
    to the input.
    Usage example from the command ine:
    matlab -nodisplay -nojvm -r "latlon2local 'input.txt' 'output.txt' 1 2 1 '\t';
    exit"
    """
    lat0=36.4042
    lon0=25.3971
    rotation=25.5
    #    #%I am not sure if this is necessary, depends how the function is used
    # ix=str2num(ix);
    # iy=str2num(iy);
    # nhl=str2num(nhl);   
    #    #%Read data and extract coordinates using indexes provided
    # if nhl == -1
    #     data=readtable(inputfile,'Delimiter',delim,'MultipleDelimsAsOne',1, ...
    #     'ReadVariableNames',1,'FileType','text');
    # elseif nhl == 0
    #     data=readtable(inputfile,'Delimiter',delim,'MultipleDelimsAsOne',1, ...
    #     'HeaderLines',0,'FileType','text');
    # else
    #     data=readtable(inputfile,'Delimiter',delim,'MultipleDelimsAsOne',1, ...
    #     'HeaderLines',nhl,'FileType','text');
    # end
    # mapx=table2array(data(:,ix));
    # mapy=table2array(data(:,iy));

    # #    #%Set up a structure defining the origin and rotation
    # #    #%This unecessary, but is kept due to legacy code
    # srGeometry=struct('latitude',lat0,'longitude',lon0,'rotation',rotation);
  
    # #%Initialize projection
    # mstruct = defaultm('mercator');
    # #%Set projection origin x y
    # origin = [srGeometry.latitude srGeometry.longitude];
    # #%Add orientation?
    # mstruct.origin = [origin 0];
    # #%Not sure why we need this
    # mstruct = defaultm(mstruct);

    # #%Apply transformation from mapx mapy to dx dy using projection defined above.
    # [dx1,dy1] = mfwdtran( mstruct, mapy, mapx);

    # #%This seems to be using a spherical Earth. Implies that the transformation above
    # #%produced coordinates in radians?
    # dx2 = rad2km( dx1, 6371 )*1000;
    # dy2 = rad2km( dy1, 6371 )*1000;

    # #%Rotate coordinates
    # rota  = srGeometry.rotation*pi/180;
    # sinrota = sin(rota);
    # cosrota = cos(rota);
    # local_x =  cosrota*dx2 + sinrota*dy2;
    # local_y = -sinrota*dx2 + cosrota*dy2;

    # #%Output
    # data.local_x = local_x;
    # data.local_y = local_y;
    # if nhl == -1
    #     writetable(data,outputfile,'Delimiter',' ','WriteVariableNames',1,'FileType','text');
    # else
    #     writetable(data,outputfile,'Delimiter',' ','WriteVariableNames',0,'FileType','text');
    # end
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08

@DataIOFactory.register_subclass('proteus')
class DataIO_PROTEUS(DataIO):
  def _extract_chnnl(self, fname):
    return int(fname[13])  
  def _extract_srcid(self, fname):
    return int(fname[9:12])
  def _get_srcids_within_extent(self, extent):
    """
    Returns
    -------
    list of integer source-IDs.

<<<<<<< HEAD
    Note that we use < instead of <=.
    Apparently that's SP / fullwave3d are 
    even stricter and exclude stations / shots
    in some vicinity of the model boundaries.
    E.g. in the Kolumbo box from the thesis, 
    one extra shot (fldr) is excluded by SP,
    giving 1505 shots per OBS in total. 
    """
    [[x1, x2], [y1, y2]] = extent
    df = pd.read_csv('%s/recs.csv' % path_prot)
    srcids = list(df.loc[(df.gx > x1) & \
      (df.gx < x2) & (df.gy > y1) & (df.gy < y2)].id)
=======
    Note that we use <= instead of <
    Check if SU / fullwave3d does the same.
    """
    [[x1, x2], [y1, y2]] = extent
    df = pd.read_csv('%s/recs.csv' % path_prot)
    srcids = list(df.loc[(df.gx >= x1) & \
      (df.gx <= x2) & (df.gy >= y1) & (df.gy <= y2)].id)
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
    return srcids
  def _set_file_class(self):
    self.DataFile = DataFileSgy
  def _set_pattern(self):
    self.pattern = 'MGL1521_????_?.sgy'
@DataIOFactory.register_subclass('proteus_hy')
class DataIO_PROTEUS_HY(DataIO_PROTEUS):
  def _set_pattern(self):
    self.pattern = 'MGL1521_????_4.sgy'  
@DataIOFactory.register_subclass('proteus_vz')
class DataIO_PROTEUS_VZ(DataIO_PROTEUS):
  def _set_pattern(self):
    self.pattern = 'MGL1521_????_1.sgy'
@DataIOFactory.register_subclass('proteus_hyvz')
class DataIO_PROTEUS_HYVZ(DataIO_PROTEUS):
  def _set_pattern(self):
    self.pattern = 'MGL1521_????_[1,4].sgy'
<<<<<<< HEAD
=======

>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
