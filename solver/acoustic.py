"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.decor import timer
from fullwavepy.generic.array import Arr3d 

# TRY Devito EVENTUALLY


@traced
@logged
class Solver(object):
  """
  """
  def __init__(self, c, dx, ns, dt, rsg, srcs, recs, **kwargs):
    """
    c: vp model
    rsg: raw signature (source/wavelet)
    """
    self.dump = kw('dump', -100, kwargs)
    if self.dump < 0:
      self.dump = -self.dump
    else:
      raise ValueError('self.dump: ' + self.dump)
    
    self.c = c
    self.dims = c.shape
    self.dx = dx
    self.ns = ns
    self.dt = dt
    assert rsg.shape[ :2] == (1,1)
    self.rsg = rsg[0][0]
    assert len(self.rsg) == self.ns
    self.srcs = srcs
    self.recs = recs
    
    self.dtdx2 = (dt*dt)/(dx*dx)
    
    self.u_prv = np.zeros(self.dims)
    self.u_now = np.zeros(self.dims)
    self.u_nxt = np.zeros(self.dims)

    self.u = [self.u_prv, self.u_now, self.u_nxt]

    self.fw = np.zeros((self.ns // self.dump, *self.dims))
  
  def solve(self, **kwargs):
    """
    """
    for u in self.u:
      u[:] = 0


# -------------------------------------------------------------------------------


@traced
@logged
class Solver2d_2t2x(Solver):
  """
  Second-order in time, second-order in space.
  
  """  
  def __init__(self, c, *args, **kwargs):
    assert c.shape[1] == 1
    super().__init__(c, *args, **kwargs)
  
  def solve(self, **kwargs):
    """
    """
    super().solve(**kwargs)
    
    for i in range(self.ns):
      self.__log.info('Step %d (of %d)\r' % (i+1,self.ns)) 
      
      spatial = +self.u_now[ :-2,:,1:-1] + self.u_now[2:  ,:,1:-1] \
                +self.u_now[1:-1,:, :-2] + self.u_now[1:-1,:,2:  ] \
              -4*self.u_now[1:-1,:,1:-1]
            
      self.u_nxt[1:-1,:,1:-1] = 2 * self.u_now[1:-1,:,1:-1] \
                                   -self.u_prv[1:-1,:,1:-1] \
                      +self.dtdx2 * spatial * self.c[1:-1,:,1:-1]**2 
      
      if i + 1 < self.ns:
        for src in self.srcs:
          x, y, z = src
          self.u_nxt[x,y,z] = self.rsg[i+1]
      
      self.u_prv[:,:,:] = self.u_now[:,:,:]
      self.u_now[:,:,:] = self.u_nxt[:,:,:]      
      
      
      if i % self.dump == 0:
        self.fw[int((i+1)/self.dump - 1)] = Arr3d(self.u_now[:,:,:])


# -------------------------------------------------------------------------------





























@traced
@logged
def solve2d(proj, sources, receivers):
  """
  Simple 2D solver (Igel)
  
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