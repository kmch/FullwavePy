"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.numeric.const import *
from fullwavepy.numeric.generic import *


@traced
@logged
def dft(y, spect='ampl', **kwargs): 
  """
  Discrete Fourier transform.
  
  Parameters
  ----------
  y : list 
    Signal to transform.
  dt : float 
    Time interval [s].
  spectrum : str 
    Type of spectrum to return:
    'ampl'/'power'/'phase'.
    Default: 'ampl'.
    
  Returns
  -------
  Y : list 
    Amplitude/power/phase.
    
  Notes
  -----
  "The values in the result follow so-called “standard” order: 
  If A = fft(a, n), then A[0] contains the zero-frequency term 
  (the sum of the.dsp), which is always purely real for real inputs. 
  Then A[1:n/2] contains the positive-frequency terms, and A[n/2+1:] 
  contains the negative-frequency terms, in order of decreasingly 
  negative frequency. 
  For an even number of input points, A[n/2] represents 
  both positive and negative Nyquist frequency, and is 
  also purely real for real input. 
  For an odd number of input points, 
  A[(n-1)/2] contains the largest positive frequency, while A[(n+1)/2] 
  contains the largest negative frequency. The routine 
  >
  > np.fft.fftfreq(n)
  >
  returns an array giving the frequencies of corresponding elements in the 
  output."  
  
  (source: https://docs.scipy.org/doc/numpy/reference/routines.fft.html) 
  
  Phase is already wrapped within [-pi;pi] by np.angle.
  
  We get rid of negative frequencies output by fft.
  
  """
  # RETURN A ZERO-TRACE UNCHANGED
  if np.count_nonzero(y) == 0.0:
    return y
  
  if spect == 'ampl':
    Y = abs(np.fft.fft(y))
  elif spect == 'power':
    Y = abs(np.fft.fft(y))**2
  elif spect == 'phase':
    Y = np.angle(np.fft.fft(y))
  else:
    raise ValueError('Unknown spectrum: ' + spect)
  
  return Y


# -------------------------------------------------------------------------------


@traced
@logged
def _set_freq_axis(n, **kwargs): # DEL?
  """
  Get the list of frequencies, as output
  from np.fft.fft.
  
  Parameters
  ----------
  n : int 
    No. of samples of time series or its DFT.
  
  """
  x = _dft_freqs(n, **kwargs)
  
  # GET RID OF NEGATIVE FREQS
  n_positive = _dft_nfreqs_positive(len(x))  
  x = np.array(x[ :n_positive]) 

  return x


# -------------------------------------------------------------------------------


@traced
@logged
def dft_freqs(N, which='all', **kwargs):  # DEL?
  """
  
  Examples
  --------
  > n = 6
  > dt = 2
  > print(_dft_freqs(n, dt=dt))
  > print(np.fft.fftfreq(n, d=dt))
  
  """
  dt = kw('dt', None, kwargs)
  if not dt:
    dt = 1
    dft_freqs._log.warning('You should provide dt to scale freq. axis')  
  
  T = dt*N
  df = 1. / T
  
  freqs = np.array([df*n if n<N/2 else df*(n-N) for n in range(N)])
  
  if which == 'all':
    pass
  elif which == 'positive': # ACTUALLY NON-NEGATIVE
    freqs = freqs[:_dft_nfreqs_positive(N)]
  else:
    raise ValueError('which: ' + which)
  
  return freqs


# -------------------------------------------------------------------------------


@traced
@logged
def _dft_nfreqs_positive(N, **kwargs): # DEL?
  # NO. OF POSITIVE FREQS 
  return int(np.ceil(N/2))


# -------------------------------------------------------------------------------


@traced
@logged
def _dft_fmax_df(N, **kwargs): # DEL?
  """
  
  Notes
  -----
  It reaches Nyquist only for even 
  number of samples!
  
  Examples
  --------
  > n = 6
  > dt = 2
  > print(_dft_freqs(n, dt=dt))
  > print(np.fft.fftfreq(n, d=dt))
  
  """
  dt = kw('dt', None, kwargs)
  if not dt:
    dt = 1
    _dft_fmax_df._log.warning('You should provide dt to scale freq. axis')  
  
  T = dt*N
  df = 1. / T
  fmax = (np.ceil(N/2) - 1) * df
  
  return fmax, df


# -------------------------------------------------------------------------------

