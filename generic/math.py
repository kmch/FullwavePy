"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw


epsi = 1e-10 # epsilon (a small number)

# -------------------------------------------------------------------------------


@traced
@logged
def neighs_of_int(a, r):
  """
  Integers surrounding an integer 'a' within radius 'r'.
  
  Examples
  --------
  r=1, a=2   returns (1,3)
  r=2, a=2   returns (0,1,3,4)
  
  """
  return np.arange(a-r, a+r+1) # INCLUDES THE POINT ITSELF (SAME AS FULLWAVE3D's HICKS)
  #return np.array(list(np.arange(a-r, a)) + list(np.arange(a+1, a+r+1))) # EXCLUDES


# -------------------------------------------------------------------------------


@traced
@logged
def neighs_of_float(a, r):
  """
  Integers surrounding a float 'a' within radius 'r'.
  
  Examples
  --------
  r=1, a=2.5   returns (2,3)
  r=2, a=2.5   returns (1,2,3,4)
  
  """  
  return np.arange(np.ceil(a)-r, np.ceil(a)+r)


# -------------------------------------------------------------------------------


@traced
@logged
def neighs1d(a, r):
  """
  Nearest nodes in 1D.
  
  a : float
  r : int
  
  Notes
  -----
  For integer 'a', 'a' is included itself.
  This is for ND discrete delta functions, 
  that needs non-zero values.
  
  This would be the version excluding int 'a':
  return np.where(float(a).is_integer(), neighs_of_int(a, r), 
                                         neighs_of_float(a, r))
  In this case we  can make use of np.where as the shape of the output 
  is the same in both situations.
  
  Examples
  --------
  r=1, a=2   returns (1,2,3)
  r=1, a=2.5 returns (2,3)
  
  """
  if float(a).is_integer():
    n = neighs_of_int(a, r)
  else:
    n = neighs_of_float(a, r)
  
  neighs1d._log.debug('nodal-neighbours of %s within r=%s are: %s' % (a, r, n))
  return n
  #
  

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
def sine(t, freq_hz, phase=0, **kwargs):
  return np.sin(2*np.pi*freq_hz*t + phase)


# -------------------------------------------------------------------------------


@traced
@logged
def sinc(x, **kwargs):
  """
  Calculate value of the sinc function defined as:
         sinc(x) = sin(pi * x) / (pi * x).
  
  Parameters
  ----------
  x : float
    Argument of the sinc(x); Any real number.
  
  Returns
  -------
  value : float
    Value of the sinc(x).
    
  Notes
  -----
  Definition containing pi number is favorable 
  because in this case Sinc has its roots at integer 
  numbers (finite-difference nodes), not pi-multiples.
  
  Numerical stability is addressed.
  
  """
  from fullwavepy.generic.math import epsi
  
  pi = 3.141592654 # TO STAY CONSISTENT WITH FULLWAVE3D
  x = np.array(x)
  # ALLOWS VECTORIZATION WITH CONDITION INSIDE
  return np.where(abs(x) < epsi, 1.0, np.sin(pi * x) / (pi * x)) 


# -------------------------------------------------------------------------------


@traced
@logged
def dsinc_dx(x, **kwargs):
  """
  sinc'(x)
  
  Notes
  -----
  Unlike sinc:
  - max > 1
  - odd function
  
  See sinc docs too.
  
  """
  from fullwavepy.generic.math import epsi
  
  pi = 3.141592654 # TO STAY CONSISTENT WITH FULLWAVE3D
  x = np.array(x)
  # NOTE dsinc_dx(0) = 0
  return np.where(abs(x) < epsi, 0.0, (np.cos(pi * x) - sinc(x)) / x)


# -------------------------------------------------------------------------------


@traced
@logged
def ricker(fpeak, ns, dt):
  """
  http://subsurfwiki.org/wiki/Ricker_wavelet
  """
  from fullwavepy.signal.wavelet import shift_to_zero

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
def kaiser(x, r=3, dipole=False, **kwargs):
  """
  Value of the Kaiser-windowing function at point x.
  Bessel function can computed using different 
  implementations.
  
  Parameters
  ----------
  b: float
    Shape of window
  r: float
    Radius of window
  bessel : str 
    Implementation to choose: fw3D or python built-in.
  x : float
    Distance from the centre of sinc. 
    
  Returns
  -------
  value : float
    Value of the Kaiser-windowed sinc function at point x.
  
  Notes 
  -----
  See Hicks 2002 for details.
  
  It can be vectorized because 'if' statement was replaced 
  with np.where.
  
  """
  from scipy.special import i0
  from fullwavepy.generic.math import epsi
  
  assert isinstance(r, int)
  assert r <= 5 
  
  # OPTIMAL VALUES FOR r=0 (-> None), r=1, ... FROM Hicks2002/Table 1 & 2
  b_optim_mono = [None, 0.00, 1.84, 3.04, 4.14, 5.26] # kmax = 2/3 pi
  b_optim_dipo = [None, 0.00, 1.48, 3.25, 4.40, 5.44] # kmax = 2/3 pi 
  
  if dipole:
    b = kw('b', b_optim_dipo[r], kwargs)
  else:
    b = kw('b', b_optim_mono[r], kwargs)

  if r != 3:
    kaiser._log.warn('Using r=%s, not 3 as in Fullwave3D.' % r)
  
  if r == 3 and b != 4.14:
    kaiser._log.warn('Using b=%s, not 4.14 as in Fullwave3D.' % b)
  
  bessel = 'py'
  if bessel == 'py':
    bessel = lambda x : i0(x)
  elif bessel == 'fw3d':
    bessel = lambda x : bessel_fw3d(x)
  else:
    raise ValueError('Unknown implementation of the Bessel function: ' + bessel)
  
  frac = (x / r) ** 2
  return np.where(frac > 1, 0.0, bessel(b * np.sqrt(1 - frac)) / bessel(b))


# -------------------------------------------------------------------------------


@traced
@logged
def bessel_fw3d(x, **kwargs):
  """
  Fullwave3D's approximation of the modified zero-order Bessel's function of the first kind.
  
  Parameters
  ----------
  x : float
    Argument of the function.
    
  Returns
  -------
  s : float
    Value of the function. Named identical to fullwave3D.
  
  Notes
  -----
  From the original source code: 
  'accuracy is <1/4%, which is quite adequate for distributed source'.
  There are discrepancies compared to Python's built-in function
  but these are for x < 0 which does not matter in Kaiser window
  (we look only at x > 0). FIXME: d-c.
  
  """
  v = x * 0.5
  a = 1.0
  s = 1.0
  i = 0
  while a > 0.03:
    i = i + 1
    a = a * (v / i)
    s = s + a ** 2
  return s


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
  output."  
  
  (source: https://docs.scipy.org/doc/numpy/reference/routines.fft.html) 
  
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

