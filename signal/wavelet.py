"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw



# -------------------------------------------------------------------------------


@traced
@logged
def pad_traces(A, nsamp_pad, **kwargs):
  """
  Pad each trace with trailing zeroes.
  
  """
  nx, ny, nz = A.shape
  An = np.zeros((nx, ny, nz+nsamp_pad))
  for x in range(nx):
    for y in range(ny):
      An[x][y][ :-nsamp_pad] = A[x][y]

  return An


# -------------------------------------------------------------------------------


@traced
@logged
def read_sushaper(fname='shaper.txt', **kwargs):
  """
  return array
  """
  from fullwavepy.ioapi.generic import read_txt
  c = read_txt(fname)
  shaper = []
  for r in c:
    for i in r:
      shaper.append(float(i))
  
  print('len(shaper)', len(shaper))
  
  A = np.zeros((1,1,len(shaper)))
  A[0][0] = shaper

  return A


# -------------------------------------------------------------------------------


@traced
@logged
def shift_to_zero(A, fraction=0.001, **kwargs): # DEL BECAUSE fraction IS ARBITRARY?
  """
  """
  from fullwavepy.signal.phase import first_breaks
  
  Az = np.zeros(A.shape)
  nx, ny, nz = A.shape
  picks = first_breaks(A, fraction, **kwargs)
  
  for x in range(nx):
    for y in range(ny):
      pick = int(picks[x][y])
      assert pick > 0
      Az[x][y][ :-pick] = A[x][y][pick: ]
  
  return Az


# -------------------------------------------------------------------------------
  

@traced
@logged
def extract_wavelet(fname, proj_name, tend=1.5, ntaper=50, **kwargs):
  """
  fname contains all the good traces, raw.
  
  """
  from fullwavepy.ioapi.generic import read_any
  from fullwavepy.ioapi.su import sugethw
  from fullwavepy.signal.su import su_process, su_mute, su_filter_full
  from fullwavepy.project.t import ProjSynVsObs
  
  f1 = 2
  f2 = 3
  f3 = 4.5
  f4 = 6.5
  zerophase = False  
  
  dt = sugethw(fname, 'dt', unique_values=1, int_values=0, **kwargs)
  dt = dt[0] / 1e6 # microsec -> sec
  
  s, g = read_coords(fname)
  A = read_any(fname)
  
  # ALIGN TO THE FIRST TRACE AND GET ITS COORDS
  i_target = 0
  Aa = align(A, i_target)
  s_dobs = s[i_target]
  g_dobs = g[i_target]  
  
  dobs = stack(Aa)
  dobs_m = su_process(dobs, su_mute, dt, [tend], mode='bottom', ntaper=ntaper)
  
  kwargs =  {'f1': f1, 'f2': f2, 'f3': f3, 'f4': f4, 'zerophase': zerophase}
  dobs_mf = su_process(dobs_m, su_filter_full, dt, **kwargs)  
  
  margin = 5000 # m
  ns = 2000
  dx = 50      
  x1 = min(s_dobs[0], g_dobs[0]) - margin    
  x2 = max(s_dobs[0], g_dobs[0]) + margin     
  y1 = min(s_dobs[1], g_dobs[1]) - margin      
  y2 = max(s_dobs[1], g_dobs[1]) + margin   
  z1 = 0
  z2 = 2 * margin   
  box = [x1, x2, y1, y2, z1, z2]
  timespace = [box, dx, ns, dt]  
  
  proj = ProjSynVsObs(proj_name+'_1', paths=paths_kmc, io='sgy', timespace=timespace)
  export_dobs(dobs_mf, proj.inp.path+proj.name+'-RawSeis.sgy', dt, s_dobs, g_dobs)
  
  
  args_rawseis = [[proj.name+'-RawSeis.sgy']]
  kwargs_wavelet = {'source': proj.inp.path+proj.name+'-RawSeis.sgy'}
  kwargs_sp = {'reciprocity': False}
  args_truevp = [1500]
  proj.inp.prepare(args_rawseis=args_rawseis, 
                   kwargs_wavelet=kwargs_wavelet,
                   args_truevp=args_truevp,
                   kwargs_sp=kwargs_sp)  
  
  
# -------------------------------------------------------------------------------


@traced
@logged
def reverbs(fname, nmax, **kwargs):
  """
  Find arrival times of reverberation
  for each trace.
  
  nmax : int 
    Max. no. of multiples.
    E.g. nmax = 2 => 0, 1 
    (direct and first multiple)
  
  We assume shot is approx. at FS.
  
  """
  vel = 1500 # m/s
  srcs, recs = read_coords(fname)
  ntr = len(srcs)
  times = np.zeros((ntr, nmax))
  
  for i in range(ntr):
    s = srcs[i]
    r = recs[i]
    seabed = r[-1] # WE ASSUME RECEIVER TO BE AT THE SEABED
    off = offset(s, r)
    for n in range(nmax):
      nlegs = 1 + n * 2 # NO. OF PATHS BETWEEN SEABED AND SEASURFACE
      tleg = np.sqrt(seabed ** 2 + (off / nlegs) ** 2) / vel
      ttotal = nlegs * tleg
      times[i][n] = ttotal
  
  return times    


# -------------------------------------------------------------------------------


@traced
@logged
def offset(s, r):
  """
  It is a horizontal distance, thus
  we skip z coord.
  
  """
  return dist_l2(s[ :-1], r[ :-1])


# -------------------------------------------------------------------------------


@traced
@logged
def dist_l2(v, w):
  """
  Euclidean (L2) distance between
  two vectors.
  
  """
  s = 0
  for i in range(len(v)):
    s = s + (v[i] - w[i]) ** 2
  s = np.sqrt(s)
  return s


# -------------------------------------------------------------------------------

@traced
@logged
def export_dobs(dobs, fname, dt, s_dobs, g_dobs, **kwargs):
  """
  
  These fields are parsed by SP along with 
  a few others that should be correct by default.
  
  """
  from fullwavepy.ioapi.segy import array2sgy
  from fullwavepy.ioapi.su import sushw
  
  array2sgy(fname, dobs, dt)
  sushw(fname, 'sx', s_dobs[0])
  sushw(fname, 'sy', s_dobs[1])
  sushw(fname, 'sdepth', s_dobs[2])
  sushw(fname, 'selev', -s_dobs[2])
  sushw(fname, 'gx', g_dobs[0])
  sushw(fname, 'gy', g_dobs[1])
  sushw(fname, 'gelev', -g_dobs[2])


# -------------------------------------------------------------------------------


@traced
@logged
def read_coords(fname, **kwargs):
  """
  Read coords of source and receiver 
  for all traces in the file.
  
  """
  from fullwavepy.ioapi.segy import read_header
  
  h = read_header(fname, **kwargs)
  
  sz = h['sdepth']
  gz = list(-np.array(h['gelev'])) # ELEV. -> DEPTH
  
  s = zip(h['sx'], h['sy'], sz) # SHOTS
  g = zip(h['gx'], h['gy'], gz) # RECEIVERS
  
  # ITERATOR (PYTHON 3) -> LIST
  s = list(s)
  g = list(g)
  
  return (s, g)


# -------------------------------------------------------------------------------


@traced
@logged
def align(A, i_target, **kwargs):
  """
  A : array 
    (ntraces, 1, nsamps)
  
  i_target : int 
    Index of the trace in A
    that other traces will be 
    aligned to.
  
  Aa : array 
    Array of the same shape as A with 
    traces aligned to i_target-trace 
    of the array A.
  
  """
  from fullwavepy.signal.generic import xcorr
  
  if A.shape[1] != 1:
    raise IOError('A shape must have ny=1')
  
  Aa = np.zeros(A.shape) 
  
  target_trace = A[i_target][0]
  for i, data_trace in enumerate(A):
    data_trace = data_trace[0]
    Aa[i][0] = xcorr(data_trace, target_trace)
  
  return Aa


# -------------------------------------------------------------------------------


@traced
@logged
def stack(A):
  """
  Stack, i.e. calculate an 
  arithmetic mean of traces.
  
  A : array
    (ntr_inline, ntr_xline, nsamps)
  """
  nx, ny, nz = A.shape
  stck = np.ones((1,1,nz))

  for x in range(nx):
    for y in range(ny):
      stck += A[x][y]

  stck /= (nx + ny)
  
  return stck


# -------------------------------------------------------------------------------


def limit_offsets(fnames, offset_min, offset_max, output_path=None, **kwargs):
  """
  
  """
  from fullwavepy.generic.parse import path_leave, strip, extend_fname
  from fullwavepy.ioapi.su import suwind

  key = 'offset'
  vmin = 0
  vmax = 300

  for fname in fnames:
    print(fname)
    name = path_leave(fname)
    nfname = './data/' + extend_fname(name, [[key, str(vmin) + '-' + str(vmax)]])
    suwind(fname, nfname, key, vmin, vmax)


@traced
@logged
def derive_wavelet(**kwargs):
  """
  Using synthetic Green's function
  
  """
  
  M = kw('M', 12000, kwargs)
  mode = kw('mode', 'same', kwargs)
  xlim = kw('xlim', None, kwargs)
  ylim = kw('ylim', None, kwargs)
  plot = kw('plot', True, kwargs)
  err = kw('err', None, kwargs)
  
  
  G = np.zeros(M)
  G[M//2] = 1
  G = kw('G', G, kwargs)
  # plt.plot(G)


  S = read_any('S.vtr')
  S = S[0][0]
  S /= np.max(S)
  S = kw('S', S, kwargs)
  # plt.plot(S)


  DOBS = np.convolve(S, G, mode=mode)
  if type(err) != type(None):
    DOBS += err
  DOBS = kw('DOBS', DOBS, kwargs)
  #plt.plot(DOBS, label='DOBS')


  S0 = read_any('S0.vtr')
  S0 = S0[0][0]
  S0 /= np.max(S0)
  S0 = kw('S0', S0, kwargs)  
  #plt.plot(S0, label='S0')

  DSYN = np.convolve(S0, G, mode=mode)
  DSYN = kw('DSYN', DSYN, kwargs)
  #ifullwavepy.plot:
    #plt.figure()
    #plt.plot(DSYN, label='DSYN')
    #plt.legend()
  # plt.plot(DOBS)



  pnoise = kw('pnoise', 1e-4, kwargs)

  S1 = su_process(tseries2array(S0), su_decon, proj.dt, 
      d_inp=tseries2array(DSYN), 
      d_out=tseries2array(DOBS), 
      pnoise=pnoise)
  S1 = S1[0][0]
  # plt.plot(S1)




  DSYN2 = np.convolve(S1, G, mode=mode)
  #ifullwavepy.plot:
    #plt.figure()
    #plt.plot(DSYN2, label='DSYN2')
    #plt.plot(DOBS, label='DOBS')
    #if type(xlim) != type(None):
      #plt.xlim(xlim)
  ## plt.xlim(5e3, 7e3)
  #plt.legend()

  #ifullwavepy.plot:
    #plt.figure()
    #plt.plot(S, label='true wavelet')
    #plt.plot(S1, label='derived wavelet')
    #plt.legend()
  
  return G, S, DOBS, S0, DSYN, S1, DSYN2


# -------------------------------------------------------------------------------

