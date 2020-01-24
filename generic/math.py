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
def decimal(a):
  """
  Decimal (fractional) part of a float a.
  """
  return a - np.floor(a)


# -------------------------------------------------------------------------------


@traced
@logged
def ricker(fpeak, ns, dt):
  """
  http://subsurfwiki.org/wiki/Ricker_wavelet
  """
  from ..signal.wavelet import shift_to_zero

  A  = np.zeros((1,1,ns))

  t = (np.arange(ns) - ns/2)* dt
  p = (np.pi * fpeak * t)**2
  A[0][0] = (1 - 2 * p) * np.exp(-p)

  #plt.plot(A[0][0])
  A = shift_to_zero(A)
  return A


# -------------------------------------------------------------------------------


@traced
@logged
def divisors_of(n):
  """
  Return all divisors of n 
  in descending order.
  
  Returns
  -------
  a generator (sth you can
  directly iterate over)
  
  Notes
  -----
  
  OBSOLETE
  1 is excluded from the list
  as trivial and for purpose below.
  
  Used to optimal resources
  to request on a cluster.
  
  """
  yield n
  for i in np.arange(n/2+1, 0, -1):
    if n%i == 0: yield int(i)


# -------------------------------------------------------------------------------


@traced
@logged
def normalize(x, **kwargs):
  """
  Guard against small denominators.
  
  """
  norm = kw('norm', 'rms', kwargs)
  x = np.array(x)
  num_zero = 1e-10
  
  if norm == 'rms':
    norm = rms(x)
  elif norm == 'max':
    norm = np.max(np.abs(x)) # NOTE: abs
  else:
    raise ValueError('Unknown norm: ' + norm)
  
  if norm > num_zero:
    x_normed = x / norm
  else:
    x_normed = np.zeros(x.shape)
  
  return x_normed  


# -------------------------------------------------------------------------------


@traced
@logged
def gauss(x, mu=0, sigma=1):
  """
  Unnormalized Gaussian distribution.
  
  Parameters
  ----------
  
  Returns
  -------
  y : type(x)
    Gaussian evaluated at x. 
  
  Notes
  -----
  Some people use alpha (1/e point)
  instead of the sigma (standard deviation)
  to define the width of the Gaussian. 
  They are related through: alpha = sigma * sqrt(2)
  
  """
  y = np.exp(-((x - mu)**2) / (2 * sigma**2))
  return y


# -------------------------------------------------------------------------------


@traced
@logged
def rms(y):
  """
  Root mean square.
  
  Parameters
  ----------
  y : list / array
    A vector of numbers.
  
  Returns
  -------
  rms : float 
    RMS of the y.
  
  Notes
  -----
  rms(y) = (y_i * y_i / n)^(1/2)
  
  It uses only built-in numpy operations
  which ensures the best performance.
  
  """
  y = np.array(y)
  n = len(y)
  return np.sqrt(np.sum(y**2) / n)


# -------------------------------------------------------------------------------


@traced
@logged
def derivative(y, **kwargs):
  """
  Calculate a derivative.
  
  Parameters
  ----------
  y : array 
    A vector of numbers.
  **kwargs : keyword arguments (optional)
    Current capabilities:
    - n : int 
      Order of derivative.
      Default: 1.
    - acc : int 
      Order of numerical accuracy.
      Default: 2.
    - dx : float
      Grid-cell size.
      Default: 1.
  
  Returns
  -------
  deriv : array 
    Derivative of y. It has 
    the same length as y.
  
  Notes
  -----
  We might as well use the scipy.misc.derivative or 
  numpy.gradient but an own implementation 
  is always more instructive.
  
  Benchmarked against numpy.gradient. Same results.
  
  Amplitude is the same as analytical derivative,
  no normalization needed - don't forget about 
  differentiation rules that are applied implicitly, 
  e.g.: d(sin(wx))/dx = w cos(wx)
  
  """
  n = kw('n', 1, kwargs) # ORDER OF DERIVATIVE
  acc = kw('acc', 2, kwargs) # ORDER OF NUMERICAL ACCURACY
  dx = kw('dx', 1, kwargs)
  
  deriv = np.zeros(len(y))
  
  if n == 1:
    if acc == 2:
      # INTERIOR
      deriv[1:-1] = [(-0.5*y[i-1] + 0.5*y[i+1]) / dx for i in range(1, len(y)-1)]

      # LEFT EDGE - FORWARD FD
      deriv[0] = ((-3/2.)*y[0] + 2*y[1] - 1/2.*y[2]) / dx 
      
      # RIGHT EDGE - BACKWARD FD (NOTE: OPPOSITE SIGNS)
      deriv[-1] = (3/2.*y[-1] - 2*y[-2] + 1/2.*y[-3]) / dx 
    
    else:
      raise NotImplementedError('Order of accuracy: %s' % acc)
      
  else:
    raise NotImplementedError('Order of derivative: %s' % n)      
      
  return deriv


# -------------------------------------------------------------------------------


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
  (the sum of the signal), which is always purely real for real inputs. 
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
  output." (https://docs.scipy.org/doc/numpy/reference/routines.fft.html)  
  
  Phase is already wrapped within [-pi;pi] by np.angle.
  
  We get rid of negative frequencies output by fft.
  
  """
  import math
  
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
    dft_freqs._log.warn('You should provide dt to scale freq. axis')  
  
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
  import math
  return math.ceil(N/2) 


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
  import math
  
  dt = kw('dt', None, kwargs)
  if not dt:
    dt = 1
    _dft_fmax_df._log.warn('You should provide dt to scale freq. axis')  
  
  T = dt*N
  df = 1. / T
  fmax = (math.ceil(N/2) - 1) * df
  
  return fmax, df


# -------------------------------------------------------------------------------

