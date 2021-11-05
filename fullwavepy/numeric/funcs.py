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
def linear(x, alpha, y0, **kwargs):
  """
  Find y-values of a linear function y = ax + b
  given X coordinates, slope angle and
  y(x[0]).
  
  Parameters
  ----------
  x : scalar / array
    X value(s).
  alpha : float 
    Slope (angle) of the function [deg].
  y0 : float 
    Value at x[0].
    
  Returns
  -------
  y : scalar / array 
    Y value(s).
    
  Notes
  -----
  
  """
  x0 = x[0]
  a = np.tan(deg2rad(alpha))
  b = y0 - a * x0
  
  return a * x + b
  

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
  x = np.array(x)
  # NOTE dsinc_dx(0) = 0
  return np.where(abs(x) < epsi, 0.0, (np.cos(pi * x) - sinc(x)) / x)


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
  return np.exp(-((x - mu)**2) / (2 * sigma**2))


# -------------------------------------------------------------------------------


@traced
@logged
def ricker(fpeak, ns, dt):
  """
  http://subsurfwiki.org/wiki/Ricker_wavelet
  """
  from fullwavepy.dsp.wavelet import shift_to_zero

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
  
  assert isinstance(r, int)
  assert r <= 5 
  
  # OPTIMAL VALUES FOR r=0 (-> None), r=1, ... FROM Hicks2002/Table 1 & 2
  b_optim_mono = [None, 0.00, 1.84, 3.04, 4.14, 5.26] # kmax = 2/3 pi
  b_optim_dipo = [None, 0.00, 1.48, 3.25, 4.40, 5.44] # kmax = 2/3 pi 
  
  if dipole:
    b = kw('b', b_optim_dipo[r], kwargs)
  else:
    b = kw('b', b_optim_mono[r], kwargs)

  #if r != 3:
    #kaiser._log.debug('Using r=%s, not 3 as in Fullwave3D.' % r)
  
  #if r == 3 and b != 4.14:
    #kaiser._log.debug('Using b=%s, not 4.14 as in Fullwave3D.' % b)
  
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

