import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.numeric.const import *
from fullwavepy.numeric.generic import *


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
