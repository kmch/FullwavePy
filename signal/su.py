"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.parse import kw, strip


# ------------------------------------------------------------------------------- 


@traced
@logged
def su_process(A, func, dt, *args, **kwargs):
  """
  Process an array. It wraps SU to
  handle arrays instead of files.
  
  Notes
  -----
  Calling SU instead of Python 
  makes sense for heavy files when SU
  is much more memory-efficient and thus
  faster (it operates on binaries).
  
  """
  from fullwavepy.ioapi.generic import read_any
  from fullwavepy.ioapi.segy import array2sgy
  
  #su_process._log.warn('shh %s' % A.shape)
  #su_process._log.warn('fffunc %s' % func)
  #su_process._log.warn('dddt %s' % dt)
  #su_process._log.warn('ars %s' % args)
  #su_process._log.warn('kars %s' % kwargs)
  
  fsgy = 'tmp.sgy'
  array2sgy(fsgy, A, dt, **kwargs)
  kwargs['dt'] = dt
  
  nfsgy = func(fsgy, *args, **kwargs)
  
  A = read_any(nfsgy)
  
  return A  


# -------------------------------------------------------------------------------


@traced
@logged
def su_decon(fname, d_inp, d_out, pnoise, **kwargs):
  """
  Wiener deconvolution using SU's
  sushape function.
  
  d_inp : 1d array
  
  Notes
  -----
  Is shaper a matrix?
  
  Bash complains about too long list of arguments
  => trough files.
  
  """
  from fullwavepy.ioapi.su import array2su
  from fullwavepy.generic.array import list2str
  from fullwavepy.ioapi.generic import read_txt
  
  if not exists(fname):
    raise IOError(fname + ' not found.')    
  
  nfname = strip(fname) + '_decon.sgy'
  #!sushape < s_inp.su wfile=d_syn.su dfile=d_obs.su pnoise=0.001 showshaper=0 > s_out.su
  
  #plt.plot(d_inp)
  
  wfile = 'wfile.su'
  dfile = 'dfile.su'
  
  array2su(wfile, d_inp, dt=1) 
  array2su(dfile, d_out, dt=1) 
  
  #return
  
  #w = list2str(d_inp)
  #d = list2str(d_out)
  
  #print(w)
  #print(d)
  
  cmd = str('segyread tape=' + fname + ' | ' +
            'sushape' + 
            ' wfile=' + wfile + 
            ' dfile=' + dfile + 
            ' pnoise=' + str(pnoise) + 
            ' showshaper=1 ' +
            ' verbose=0 2> shaper.txt' + ' | ' +
            'segyhdrs | ' + 
            'segywrite tape=' + nfname)  
  
  #c = read_txt('shaper.txt')
  #shaper = []
  #for r in c:
    #for i in r:
      #shaper.append(float(i))
      
  #print('len(shaper)', len(shaper))
  #plt.plot(shaper)
  #print(c[0])
  #print(c)
  
#  cmd = str('segyread tape=' + fname + ' | ' +
#            'sushape' + 
#            ' w=' + w + 
#            ' d=' + d + 
#            ' pnoise=' + str(pnoise) + 
#            ' showshaper=1 ' +
#            ' verbose=1' + ' | ' +
#            'segyhdrs | ' + 
#            'segywrite tape=' + nfname)
  # 2> shaper.txt
  #print(cmd)
  o, e = bash(cmd)
  
  
  #print('o,e', o,e)
  return nfname


# -------------------------------------------------------------------------------


#@traced
#@logged
#def su_trim(fname, pad_s, tracelen, **kwargs):
  #"""
  #opposit of su_pad
  #"""
  #pass


# -------------------------------------------------------------------------------


@traced
@logged
def su_taper(fname, taper_front, taper_back, **kwargs):
  """
  Taper all the traces.
  
  taper_...
  in samples
  
  """
  from fullwavepy.ioapi.su import get_ntraces, sugethw

  if not exists(fname):
    raise IOError(fname + ' not found.')  

  taper_type = kw('taper_type', 'cos', kwargs)
  if taper_type == 'linear':
    taper_type = '1'
  elif taper_type == 'sin':
    taper_type = '2'
  elif taper_type == 'cos':
    taper_type = '3'
  elif taper_type == 'gauss_wide':
    taper_type = '4'
  elif taper_type == 'gauss_narr':
    taper_type = '5' 
  else:
    raise ValueError('Wrong taper_type: ' + taper_type)
  
  nfname = strip(fname) + '_taper.sgy'
  
  ntr = get_ntraces(fname, **kwargs)
  
  # FIXME int_values=False NOT SUPPORTED
  dt_us = (sugethw(fname, ['dt'], unique_values=True, **kwargs))['dt'][0] # micro sec
  dt_ms = dt_us / 1000.  
  su_taper._log.info('dt_ms' + str(dt_ms))
  
  tbeg = taper_front * dt_ms
  tend = taper_back * dt_ms
  su_taper._log.info('tbeg %s, tend %s' % (tbeg, tend))
  

  cmd = str('segyread tape=' + fname + ' | ' +
            'sushw key=ntr a=' + str(ntr) + ' | ' +
            'sutaper' +
            ' tbeg=' + str(tbeg) +
            ' tend=' + str(tend) + 
            ' taper=' + str(taper_type) + ' | ' +
            'segyhdrs | ' + 
            'segywrite tape=' + nfname)
  
  o, e = bash(cmd)
  
  return nfname


# -------------------------------------------------------------------------------

   
@traced
@logged
def su_pad(fname, pad, trim=False, **kwargs):
  """
  Pad or trim at both ends of each trace.
  
  [pad]  = samples
  
  Notes
  -----
  Call su_shift?
  
  sign=-1 => shift towards later times.
  
  """
  from fullwavepy.ioapi.su import sugethw

  if not exists(fname):
    raise IOError(fname + ' not found.')  

  # READ NEEDED PARAMS  
  ns = (sugethw(fname, ['ns'], unique_values=True, **kwargs))['ns'][0]
  dt_us = (sugethw(fname, ['dt'], unique_values=True, **kwargs))['dt'][0] # micro sec
  dt_ms = dt_us / 1000.
  
  su_pad._log.debug('No. of samples read from header: ' + str(ns))
  su_pad._log.debug('dt (ms) read from header: ' + str(dt_ms))
  
  pad_ms = pad * dt_ms
  #pad_samples = pad_ms // dt_ms
  #su_pad._log.info('Padding: ' + str(pad_samples) + ' samples ')

  # BUILD A COMMAND FOR BASH
  cmd = str('segyread tape=' + fname + ' | ')
  
  if trim:
    nfname = strip(fname) + '_trim.sgy'
    cmd += str('sushw key=tstat a=' + str(pad_ms) + ' | ' +
               'sustatic hdrs=1 sign=1 | ' +
               'suvlength ns=' + str(ns - 2 * pad) + ' | ') 
  else:
    nfname = strip(fname) + '_pad.sgy'
    cmd += str('suvlength ns=' + str(ns + 2 * pad) + ' | ' +
               'sushw key=tstat a=' + str(pad_ms) + ' | ' +
               'sustatic hdrs=1 sign=-1 | ')
   
   
  cmd += str('segyhdrs | ' + 
             'segywrite tape=' + nfname)
  
  o, e = bash(cmd)

  return nfname


# -------------------------------------------------------------------------------  


@traced
@logged
def su_shift(fname, tshift_ms, up=False, **kwargs):
  """
  tshift_ms : float
    In mili seconds.
  
  """
  if not exists(fname):
    raise IOError(fname + ' not found.')  

  nfname = strip(fname) + '_shift.sgy'  
  
  if up:
    sign = 1
  else:
    sign = -1
  
  cmd = str('segyread tape=' + fname + ' | ' +
            'sushw key=tstat a=' + str(tshift_ms) + ' | ' +
            'sustatic hdrs=1 sign=' + str(sign) + ' | ' + 
            'segyhdrs | ' + 
            'segywrite tape=' + nfname)
  
  #print(cmd)
  o, e = bash(cmd)
  #print(o,e)
  return nfname 


# -------------------------------------------------------------------------------


@traced
@logged
def su_mute(fname, picks, mode, ntaper, **kwargs):
  """
   picks : list
     1D list of times [sec] 
    
  twin 
    seconds 
    
  ntaper 
    samples
    
  """
  from fullwavepy.generic.array import list2str
  
  if not exists(fname):
    raise IOError(fname + ' not found.')    
  
  key = kw('key', 'tracr', kwargs)
  nfname = strip(fname) + '_mute.sgy'
  picks = np.array(picks)
  
  # SELECT THE MODE
  if mode == 'top':
    mode_list = [0]
    picks_list = [picks]
  
  elif mode == 'bottom':
    mode_list = [1]
    picks_list = [picks]
  
  elif mode == 'both':
    try:
      twin = kwargs['twin']
    except KeyError:
      raise IOError('You need to provide twin for the mode: ' + mode)
    mode_list = [0, 1]
    picks_list = [picks, picks + twin]
  
  else:
    raise ValueError('Wrong mode: ' + str(mode))

  
  # START THE COMMAND
  cmd = str('segyread tape=' + fname + ' | ')
  
  # POSSIBLY MULTIPLE MUTES
  for mode, picks in zip(mode_list, picks_list):
    xmute = list2str(range(1, len(picks)+1))
    tmute = list2str(picks)
    cmd += str('sumute ' + 
               ' key=' + key + 
               ' nmute=' + str(len(picks)) +
               ' mode=' + str(mode) + 
               ' ntaper=' + str(ntaper) +                
               ' xmute=' + xmute + 
               ' tmute=' + tmute + ' | ')
  
  # FINISH THE COMMAND
  cmd += str('segyhdrs | ' + 
             'segywrite tape=' + nfname)    
             
  #print('koommadn', cmd)
  o, e = bash(cmd, **kwargs)
  ##print(o, e)
 
 
              #' twindow=' + twindow + ' | ' + 
              #' linvel=' + linvel +
              #' tm0=' + str(picks[0] / 2) + 
  #linvel = 'none'
  #twindow = 'none'
  #if mode == 3: # HYPERB.
    #linvel = '100'
  
  #if twindow: and mode == 4
  #  twin = ''
  #  for i in range(len(picks)):
  #    twin += str(twindow) + ','
  #  twindow = twin[:-1]
  #else:
  #  twindow = 'none' 
 
  return nfname  


# -------------------------------------------------------------------------------


@traced
@logged
def su_filter(fname, **kwargs):
  """
  Use Seismic Unix' Butterworth filter.
  
  Parameters
  ----------
  fname : str
  zerophase : bool 
  
  kwagrgs
  if f1 and f2 - highpass
  if f3 and f4 - lowpass
  if all - bandpass
  
  Returns
  -------
  nfname : str 
    Name of the result file.
  
  Notes
  -----
  Creates a new file.
  
  It is much faster than any Python filter 
  because it operates on binaries.
  
  To interface SU and Jupyter we always 
  need to pass a .sgy file.
  
  verbose=1 prints the information of 
  3db point and no. of poles that are 
  determined automatically for given
  fstop, fpass unless we choose the 
  2nd method of defining the filter:
   npoleselo=calculated     number of poles of the lo pass band
   npolesehi=calculated     number of poles of the lo pass band
   f3dblo=calculated	frequency of 3db cutoff frequency
   f3dbhi=calculated	frequency of 3db cutoff frequency
  
  """
  if not exists(fname):
    raise IOError(fname + ' not found.')    
  
  zerophase = kw('zerophase', False, kwargs)
  
  if zerophase:
    zerophase = str(1)
  else:
    zerophase = str(0)
  
  # NAME THE OUTPUT
  nfname = strip(fname) + '_filt.sgy'
  su_filter._log.debug('File name: ' + fname)
  su_filter._log.debug('New file name: ' + nfname)  
  
  su_filter._log.warning('Zerophase: ' + zerophase)
  
  cmd = str('segyread tape=' + fname + ' | ' +
             'subfilt')
  
  locut = '0'
  hicut = '0'
  
  if not (('f1' in kwargs and 'f2' in kwargs) or 
          ('f3' in kwargs and 'f4' in kwargs)):
    raise IOError('No filter frequencies provided.')
  
  # LOW-CUT BIT
  if ('f1' in kwargs) and ('f2' in kwargs):
    f1 = str(kwargs['f1'])
    f2 = str(kwargs['f2'])
    locut = '1'
    cmd += str(' fstoplo=' + f1 + 
               ' fpasslo=' + f2)
  # HIGH-CUT BIT
  if ('f3' in kwargs) and ('f4' in kwargs):
    f3 = str(kwargs['f3'])
    f4 = str(kwargs['f4'])   
    hicut = '1'
    cmd += str(' fpasshi=' + f3 + 
               ' fstophi=' + f4)
    
  cmd += str(' zerophase=' + zerophase +
             ' locut=' + locut  +
             ' hicut=' + hicut +
             ' verbose=1 | ' +
             'segyhdrs | ' + 
             'segywrite tape=' + nfname)
             
  #print(cmd)
  
  o, e = bash(cmd, **kwargs)
  #print(o, e)
  
  return nfname


# -------------------------------------------------------------------------------


@traced
@logged
def su_filter_full(fname, pad, **kwargs):
  """
  Includes tapering, padding, etc.
  
  [pad]  = samples
  
  """
  from fullwavepy.generic.parse import del_kw
  
  
  if not exists(fname):
    raise IOError(fname + ' not found.')  
  
  #tbeg_ms = kw('tbeg_ms', 1e4, kwargs)
  #del_kw('tbeg_ms', kwargs)
  #tend_ms = kw('tend_ms', tbeg_ms, kwargs)
  #del_kw('tend_ms', kwargs)  
  #pad_ms = kw('pad_ms', 100, kwargs)
  #del_kw('pad_ms', kwargs)
  
  
  
  #print(tbeg_ms, tend_ms, pad_ms)
  #print(kwargs)
  
  fname_taper = su_taper(fname, pad, pad, **kwargs) # tbeg_ms, tend_ms, **kwargs)  
  fname_taper_pad = su_pad(fname_taper, pad, **kwargs)
  fname_taper_pad_filt = su_filter(fname_taper_pad, **kwargs)
  fname_taper_pad_filt_trim = su_pad(fname_taper_pad_filt, pad, trim=True, **kwargs)
  
  return fname_taper_pad_filt_trim
 

# -------------------------------------------------------------------------------

