"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced


@traced
@logged
class Solver(object):
  """
  """
  def __init__(self, dims, dx, ns, dt, rsg, **kwargs):
    """
    rsg: raw signagure (source/wavelet)
    """
    self.dims = dims
    self.dx = dx
    self.ns = ns
    self.dt = dt
    self.rsg = rsg
    
    self.u_prv = np.zeros(dims)
    self.u_now = np.zeros(dims)
    self.u_nxt = np.zeros(dims)

    self.u = [self.u_prv, self.u_now, self.u_next]
  
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
  def solve(self, **kwargs):
    """
    """
    super().solve(self, **kwargs)
    spatial = (+u[:-2,1:-1]+u[2:,1:-1]-4*u[1:-1,1:-1]+u[1:-1,:-2]+u[1:-1,2:])
    self.u_nxt[1:-1,1:-1] = 2*u[1:-1,1:-1] - u_prv[1:-1,1:-1]  \
            + (c[1:-1,1:-1]**2) * dtdx2 * 



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



def propagate():
    
    # make sure everything starts off zero
    u[:] = 0.0
    u_prv[:] = 0.0

    u[sx,sz] = src[0]  # inject first source entry into current wavefield

    # begin time-stepping loop...

    for i in range(nt):

        if i%20==0:
            sys.stdout.write('Done step %d (of %d)\r' % (i+1,nt))

        # FILL IN CODE HERE TO CALCULATE THE NEW WAVEFIELD u_nxt (task 2)...
    
    #    NOTE: The commented code below uses loops, which is how we did it in the 1d code
    #    However, you'll find it's really slow doing it like that, so see quicker code afterwards... 
    #    for ix in range(1,nx-1):
    #        for iz in range(1,nz-1):
    #            u_nxt[ix,iz] = 2*u[ix,iz] - u_prv[ix,iz] + (c[ix,iz]**2) * dtdx2 *  \
    #                                (-4*u[ix,iz]+u[ix-1,iz]+u[ix+1,iz]+u[ix,iz-1]+u[ix,iz+1])

    #   This code is much quicker, since it works on (almost) the whole array at once
    #   (apart from the edges, which we want to leave alone - note the index bounds)
        u_nxt[1:-1,1:-1] = 2*u[1:-1,1:-1] - u_prv[1:-1,1:-1]  \
            + (c[1:-1,1:-1]**2) * dtdx2 * (+u[:-2,1:-1]+u[2:,1:-1]-4*u[1:-1,1:-1]+u[1:-1,:-2]+u[1:-1,2:]) 

        # inject source entry for this step at the source point
        if i+1<ns:
            u_nxt[sx,sz] = src[i+1]
    
        # shift wavefields for next time-step
        u_prv[:,:] = u[:,:]
        u[:,:] = u_nxt[:,:]
    

        # ADD CODE HERE TO PLACE WAVEFIELD VALUES ALONG RECEIVER LINE INTO ARRAY r (task 3)...
    
        r[:,i] = u[:,rz]


    
        if (i+1)%snapshot_gap == 0: # store the current wavefield u on every tenth step
            wavefield[int((i+1)/snapshot_gap-1)] = u[:,:]


# -------------------------------------------------------------------------------