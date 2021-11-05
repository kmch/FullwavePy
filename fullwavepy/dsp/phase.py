"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw


@logged
def wrap_phase(phase, **kwargs):
  """
  To stay within -pi, pi.
  
  """
  if phase >= np.pi:
    phase = phase - 2 * np.pi
  elif phase <= -np.pi:
    phase = phase + 2 * np.pi
  
  return phase
@logged
def extract_phase(A, picks, dt, freq, **kwargs):
  """
  
  Parameters
  ----------
  A : array 
  
  picks : array
    Same shape as A.
  
  Returns
  -------
  phase : array 
    Same shape as A.
  
  Notes
  -----
  
  """
  from fullwavepy.numeric.fourier import dft
  
  A, _ = _window_data(A, picks, **kwargs)
  phase = np.zeros(picks.shape)
  #plot(A_win)
  
  # INDEX OF FREQ. 
  f_i = int(freq * dt * len(A[0][0]))
  
  i = 0
  nx, ny, nt = A.shape
  for x in range(nx):
    for y in range(ny):
      trace = A[x][y]
      phi_trace = dft(trace, dt=dt, spect='phase')
      phase[i] = phi_trace[f_i] 
      i += 1
  
  return phase
@logged
def _window_data(A, picks, **kwargs):
  """
  """
  A_win = np.zeros(A.shape)
  nx, ny, nt = A.shape
  ntraces = nx * ny
  if ntraces != len(picks):
    raise ValueError('No. of picks: ' + str(len(picks)) + 
                     ' must be equal to no. of traces: ' + str(ntraces))  
  i = 0
  wins = np.zeros(A.shape)
  for x in range(nx):
    for y in range(ny):
      trace = A[x][y] 
      pick = picks[i]
      win = _gauss_window(len(trace), pick, **kwargs)
      wins[x,y,:] = win
      A_win[x][y] = trace * win
      i += 1
  
  return A_win, wins
@logged
def _gauss_window(nsamp, center, gauss_alpha=80, **kwargs):
  """
  Compute Gaussian window.
  
  Parameters
  ----------
  
  Returns
  -------
  
  Notes
  -----
  
  """
  from fullwavepy.numeric.funcs import gauss
  
  sigma = gauss_alpha / np.sqrt(2)  
  win = np.arange(nsamp)
  
  win = gauss(win, mu=center, sigma=sigma)
  return win
@logged
def first_breaks(A, fraction=0.01, **kwargs):
  """
  Pick first breaks of the gather.
  
  Parameters
  ----------
  A : array
    Data array.
  fraction : float
    Fraction of max amplitude of each
    trace treated as a first break.
  
  Returns
  -------
  
  Notes
  -----
  WE NEED TO TAKE THEM FOR SYNTHETICS ONLY
  (OBS IS NOISY)
  
  """
  first_breaks._log.info('Picking first breaks as first sample above ' +
                         str(fraction) + ' of the max amplitude of the trace.')
  
  nx, ny, nt = A.shape
  picks = np.zeros((nx, ny))
  
  for x in range(nx):
    for y in range(ny):
      trace = A[x][y]
      maxx = np.max(np.abs(trace))
      for t in range(nt):
        if abs(trace[t]) > fraction * maxx:
          picks[x][y] = t
          break
        
  return picks
