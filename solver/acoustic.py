"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced


# -------------------------------------------------------------------------------


@traced
@logged
def solve2d(proj, sources, receivers):
  """
  Simple 2D solver.
  
  Parameters
  ----------

  Returns
  -------
  
  Notes
  -----
  Based on H. Igel's solver ("Computational Seismology")
  
  """
  dt = proj.dt
  dx = proj.dx
  dz = dx
  nt = proj.ns
  nx = proj.dims[0]
  nz = proj.dims[2]
  x = np.arange(0, nx) * dx
  z = np.arange(0, nz) * dz 
  v = proj.vp.true.array[:, 0, :].T # DOUBLE-CHECK TRANSPOSE
  
  # CHECK CFL CONDITION
  C_max = 0.5 # COURANT NUMBER
  vmax = np.max(v)
  C = vmax * dt / dx
  if C > C_max:
    raise ValueError('Stability criterion violated: C = %s > % s C_max' % (C, C_max))
  
  # INIT ARRAYS
  p = np.zeros((nz, nx))          # p at time n (now)
  pold = np.zeros((nz, nx))       # p at time n-1 (past)
  pnew = np.zeros((nz, nx))       # p at time n+1 (present)
  d2px = np.zeros((nz, nx))       # 2nd space derivative of p in x-direction
  d2pz = np.zeros((nz, nx))       # 2nd space derivative of p in z-direction
  wavefield = np.zeros((nt, nz, nx))  # WAVEFIELD
  for r in receivers:
    r.data = np.zeros(nt)
 
 
  for it in range(0, nt): # 5-point operator FD scheme

    # Space derivative in x-direction
    for i in range(2, nx - 2):
      d2px[i, :] = (- 1. / 12 * p[i + 2, :] + 4. / 3  * p[i + 1, :] - 5. / 2 * p[i, :] \
              + 4. / 3  * p[i - 1, :] - 1. / 12 * p[i - 2, :]) / dx ** 2

    # Space derivative in z-direction
    for j in range(2, nz - 2):
      d2pz[:, j] = (- 1. / 12 * p[:, j + 2] + 4. / 3  * p[:, j + 1] - 5. / 2 * p[:, j] \
              + 4. / 3  * p[:, j - 1] - 1. / 12 * p[:, j - 2]) / dz ** 2 

    # Time Extrapolation
    pnew = 2 * p - pold + dt ** 2 * v ** 2 * (d2px + d2pz)

    # Add Source Term at isx, isz
    for s in sources:
      #print('z,x: %s %s' % (s.z, s.x))
      pnew[s.z, s.x] = pnew[s.z, s.x] + s.wavelet[it] / (dx * dz) * (dt ** 2)

    wavefield[it, :, :] = pnew
    
    #FIXME
    #if it > 900 and it < 1400:
      #print(wavefield[it, s.x, s.z])

    # Remap Time Levels
    pold, p = p, pnew

    # Save Seismograms
    for r in receivers:
      r.data[it] = p[r.z, r.x]
    
  return wavefield, receivers

  
# -------------------------------------------------------------------------------