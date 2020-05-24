"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw


# -------------------------------------------------------------------------------


@traced
@logged
def xcorr(data, target, **kwargs):
  """
  Cross-correlate two 1D.dsps.
  
  Notes
  -----
  Detrend / taper first?
  
  """
  from scipy.dsp import correlate

  shift = np.argmax(correlate(data, target)) - len(target) #+ 10
  datan = np.zeros(len(data))# - shift)
    
  if shift > 0:
    datan[ :-shift] = data[shift: ]
  elif shift < 0:
    shift = abs(shift)
    datan[shift: ] = data[ :-shift]
  else:
    datan = data
  
  return datan


# -------------------------------------------------------------------------------

