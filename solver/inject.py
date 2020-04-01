"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

This library provides procedures for injection 
of a point source onto a finite-difference grid 
in the vicinity of an irregular free surface.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.generic.decor import timer


@traced
@logged
class Point(np.ndarray):
  """
  """
  def __new__(cls, xyz, **kwargs):
    return np.asarray(xyz).view(cls)

  # -----------------------------------------------------------------------------
  
  def __array_finalize__(self, obj):
    if obj is None: return

  # -----------------------------------------------------------------------------

  def find_neighs(self, r, **kwargs):
    """
    Find a cube of nodes surrounding the point.
    """
    from fullwavepy.generic.math import neighs1d
    x, y, z = self
    X = neighs1d(x, r)
    Y = neighs1d(y, r)
    Z = neighs1d(z, r)
    return np.array(np.meshgrid(X, Y, Z, indexing='ij')).T.swapaxes(0, 2) 

  # -----------------------------------------------------------------------------    
  
  def spread(self, r, funcx, funcy, funcz, **kwargs):
    """
    Spread the point onto a cuboid using
    possibly different discrete, band-limited 
    Dirac delta along each coordinate axis.
    
    Notes
    -----
    See Hicks 2002, Geophysics for details.
    
    We use the same r for neighbours and the window as 
    outside the window values are zero by definition.
    
    """
    # CUBE OF (x,y,z) TUPLES
    cube = self.find_neighs(r, **kwargs)
    # CENTER THE COORDINATE SYSTEM AT self
    dists = cube - self
    # START WITH ONES TO MULTIPLY BY NEW VALUES
    spread = np.ones(dists[...,0].shape)
    # APPLY ALONG TUPLE AXIS, I.E. TAKE POINTS COORDS AS AN ARGUMENT
    axis = 3
    # DEAL WITH ONE COORDINATE AT A TIME
    for i, func1d in enumerate([funcx, funcy, funcz]):
      func_of_xyz = lambda point : func1d(point[i])
      # ND-DELTA IS A PRODUCT OF 1D ONES
      spread *= np.apply_along_axis(func_of_xyz, axis, dists)
    return spread

  # -----------------------------------------------------------------------------
  
  def interp(self, **kwargs):
    pass
  
  def interp_hicks(self, **kwargs):
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class Monopole(Point):
  """
  """
  def spread(self, r, **kwargs):
    from fullwavepy.generic.math import kaiser, sinc
    func = lambda x : kaiser(x, r) * sinc(x)
    return super().spread(r, func, func, func, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class Dipole(Point):
  """
  axis : 0, 1 or 2
    corresponds to dipole along X, Y or Z axis respectively
  
  """
  def __new__(cls, xyz, axis, **kwargs):
    assert axis in [0, 1, 2]
    cls.axis = axis
    return super().__new__(cls, xyz, **kwargs)
  
  def spread(self, r, **kwargs):
    from fullwavepy.generic.math import kaiser, sinc, dsinc_dx
    
    func1 = lambda x : kaiser(x, r) * sinc(x) 
    func2 = lambda x : kaiser(x, r) * dsinc_dx(x)
    
    funcs = [func1, func1, func1]
    funcs[self.axis] = func2
    
    return super().spread(r, *funcs, **kwargs)
  

# -------------------------------------------------------------------------------















@traced
@logged
def xyz2w(xyz, dims, **kwargs):
  """
  """
  x, y, z = xyz
  nx, ny, nz = dims
  return (x - 1) * ny * nz + (y - 1) * nz + z


# -------------------------------------------------------------------------------


@traced
@logged
class SrcRec(Point):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class Src(SrcRec):
  def check_fs_pos(self, **kwargs):
    pass  
  def spread_n_bounce(self, **kwargs):
    pass
  #def spread_factors(self, **kwargs):
    #self.find_neighs()
  def spread_bounce(self, **kwargs):
    pass


# -------------------------------------------------------------------------------


@traced
@logged
class SuperSrc(Src):
  """
  """
  def check_fs_pos(self, **kwargs):
    pass
  
  def spread_factors(self, **kwargs):
    nsrcs = []
    while diverged:
      for src in srcs:
        nsrcs.append(src.spread_n_bounce())
      srcs = nsrcs
      self._check_convergence()
  
  def _check_convergence():
    pass
  
  def inject(self, wf, **kwargs):
    pass
  

# -------------------------------------------------------------------------------
















def Nearest_Neighbours(ndims, point, radius, include_point, **kwargs):
  """
  Find grid nodes which are the nearest neighbours
  of the point.
  
  Parameters
  ----------
  ndims : int 
    Dimensionality of the problem: 1, 2, or 3.
  point : array / list
    Vector of 3D coordinates [x, y, z].
    They do NOT need to be integers.
  radius : int 
    Radius of the ball (in the taxi-cab metric 
    i.e. along the grid lines; it's a cube) 
    centred on the point that contains only 
    the 'nearest neighbours' of the point.
  include_point : bool 
    If false, the neighbours cannot lie on 
    the same grid line as the given point. 
    This grid line is skipped and the 
    next one outwards is chosen for adjacent
    neighbours. This matters only if any of 
    point's coordinate is an integer 
    (<=> lies along a grid line).
    
  Returns
  -------
  coord_lists : list 
    List of lists of coordinates:
    [x_coord_list, y_coord_list, z_coord_list].
    They can differ in length, in particular 
    some of them can have unit-length.
   
  Notes
  -----
  2D => y-coord is fixed.
  1D => y-coord and z-coord are fixed.
  
  Before deleting the middle point of the list
  has always odd number of elements.
  
  FIXME; Check how Hicks is meant to work exactly.
  Is it alwasy 7x7x7?
  
  # 3 TO GET 7-POINT HICKS STENCIL # FIXME: IT GIVES 6 POINT STENCIL WITH False BELOW!!!
  
  """
  this_func = this_lib + 'Nearest_Neighbours: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  from lib_generic import Delete_List_Midpoint
  
  # CHECK THE INPUT
  if len(point) != 3:
    eprint(this_func + 'Error. Point has to have 3 coordinates.\n')
    eprint(this_func + 'Point: ' + str(point) + '\n')
    quit()
  
  # PRECISION OF TREATING COORDINATE AS INT
  precision = epsi # FIXME: TUNE IT.
  
  # ITERATE OVER COORDINATES
  coord_lists = []
  for i in range(len(point)):
    coord_list = []
    
    # SPECIAL TREATMENT OF Y-COORD IN 2D
    if (ndims == 2) and (i == 1):
      coord_lists.append([point[i]]) 
      continue
    
    # SPECIAL TREATMENT OF Y-COORD AND Z-COORD IN 1D
    if (ndims == 1) and (i > 0):
      coord_lists.append([point[i]]) 
      continue    
    
    # 'FLOOR-NEAREST' NODE OF THE POINT
    nint = int(point[i]) 
    
    # DEFINE THE RANGES OF COORDINATES
    rmin = nint - radius
    rmax = nint + radius
    if not Is_Int(point[i], precision):
      rmin += 1
   
    # GET THE LIST OF COORDINATES
    coord_list = np.arange(rmin, rmax + 1)
    
    # DELETE THE MIDDLE POINT IF NECESSARY
    if (Is_Int(point[i], precision)) and (not include_point):      
      coord_list = Delete_List_Midpoint(coord_list)
    
    coord_lists.append(list(coord_list))
  
  #print this_func + 'END' 
  return coord_lists


# -------------------------------------------------------------------------------





# -------------------------------------------------------------------------------




# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# PRECALCULATION (BASE)
# -------------------------------------------------------------------------------


def Source_Spread_Shallow_Iteratively(proj_name, orig_src, w0, **kwargs):
  """
  Spread a point source located arbitrarily close
  to the free surface.
  
  Parameters
  ----------
  
  
  Returns
  -------
  
  Notes
  -----
  The motivation is the need to reflect the part 
  of the source region that sticks out above the FS.
  
  One cannot use the immersed boundary method 
  because that assumes the knowledge of the wavefield
  inside the model at a current time step. In this 
  problem part of the wavefield is spread outside the model
  and one need to 'put it back' before even thinking about 
  immersed-boundary step.
  
  However, we utilize IBM in using flags (interior/exterior)
  and ghosts and associated intersections. In this way
  we know how to reflect the energy back inside the model.
  
  Nomenclature:
  hsrcs - point sources to spread, "Hicks sources"
  cube  - source region
  
  """
  this_func = this_lib + 'Source_Spread_Shallow_Iteratively: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  from lib_generic_CONST import hicks_radius
  from lib_fwi_fs import Read_GhostData_File_txt
  from lib_generic_PLOTT import Plot_Slices_XYZ
  
  # CHECK KEYWORD ARGUMENTS
  thresh = Kwarg('thresh', 1e-4, kwargs)
  i_max = Kwarg('i_max', 5, kwargs)
  ndims = len(w0.shape)
  
  # READ GHOST DATA ONCE
  ghost_data = Read_GhostData_File_txt(proj_name+'-GhostData.txt')
  flags = Flags_Read(proj_name, **kwargs)
  
  hsrcs = [orig_src] # START WITH ONLY 1 SOURCE
  i = 0
  while i < i_max:
    i += 1 
    print(this_func, 'Starting ', i, ' iteration,  no. of Hicks sources: ', len(hsrcs))
    
    # MAIN BIT
    new_hsrcs = []
    for hsrc in hsrcs: # SPREAD HICKS SOURCES 1-BY-1
      new_hsrcs, w0 = Source_Spread_Shallow(hsrc, new_hsrcs, ghost_data, flags, w0, ndims=ndims, **kwargs)    
    hsrcs = new_hsrcs
    
    #plt.figure()
    #Plot_Slices_XYZ(vols=[w0], cmap='seismic', **kwargs)
    
    # CONVERGENCE CHECK
    if Converged(hsrcs, thresh, orig_src[-1]):
      break
  
  #print this_func + 'END' 
  return hsrcs


# -------------------------------------------------------------------------------


def Source_Spread_Shallow(src, new_hsrcs, ghost_data, flags, w0, **kwargs):
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  new_hsrcs : list
    Point sources to spread in the next
    iteration.
  w0 : array
    Updated wavefield array.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Source_Spread_Shallow: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  ndims = Kwarg('ndims', 3, kwargs)
  if ndims == 2:
    two_d = True
  else:
    two_d = False
  
  ghosts, intersects, ficts, auxs, weights = ghost_data
  
  # SMEAR THE OFF-GRID POINT SOURCE ONTO A CUBE OF GRID-NODES
  cube = Source_Spread(src, ndims)
  
  # SPLIT THE CUBE INTO 2 PARTS ON OPPOSITE SIDES OF FS
  cube_below, cube_exact, cube_above = Cube_Split(cube, flags)

  # ADD THE INTERIOR PART TO THE WAVEFIELD
  w0 = Wavefield_Update(w0, cube_below, **kwargs)
  
  # FIND THE MIRROR-REFLECTIONS OF THE EXTERIOR PART
  cube_refl = Reflect_Beneath_FS(cube_above, ghosts, intersects, w0.shape, two_d)

  # ADD THEM TO THE EXISTING LIST (NOTE THE LOOP IN THE PARENT FUNCTION)
  new_hsrcs = Update_Sources_List(new_hsrcs, cube_refl)

  #print this_func + 'END'
  return new_hsrcs, w0


# -------------------------------------------------------------------------------


def Source_Spread(src, ndims, **kwargs):
  """
  Spread a point-source onto a finite grid.
  
  Parameters
  ----------
  ndims : int 
    Dimensionality of the problem: 1, 2, or 3.
  src_xyz : vector 
    Vector of 3D coordinates of the source to spread.
  src_amp : float 
    Amplitude of the source function at a given instant.
  radius : int 
    Radius of the source region. In FW3D it is 3 which 
    gives a 7x7x7 cube of the source region (in 3D).
  scaling : float
    Extra scaling factor for sinc to speed up the convergence.
  
  Returns
  -------
  source_nodes : list 
    List of nodes of the form [x, y, z, amplitude(x, y, z)] 
    onto which the point source has been spread.
  
  Notes
  -----
  Fullwave's Hicks corresponds to radius=3, include_point=T,
  i.e. in each direction, and b = 4.14.
  
  If all coords are off-grid one gets 6x6x6 = 216 source nodes.
  Otherwise more (6x6x7, 6x7x7, 7x7x7).
  
  """
  this_func = this_lib + 'Source_Spread: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START') 
  
  from lib_math_signal import Sinc, Windowed, Window_Kaiser
  
  scaling = Kwarg('scaling', 1.0, kwargs)
  radius = Kwarg('radius', hicks_radius, kwargs) 
  
  src_xyz, src_amp = np.array(src[ :3]), src[3]
  
  # COORDINATES OF THE SOURCE-REGION NODES
  include_point = True # OTHERWISE SIZE OF THE SOURCE REGION IS NOT CONSTANT
  x_list, y_list, z_list = Nearest_Neighbours(ndims, src_xyz, radius, include_point)
  
 
  source_nodes = []
  # ITERATE OVER SOURCE-REGION (CUBE) NODES
  # WE NEED TO DO IT 1 
  for x in x_list:
    for y in y_list:
      for z in z_list: 
        ampl = src_amp * scaling # NOTE: APPLIED ONLY ONCE # FIXME: D-C
        r = np.array([x, y, z])
        # ITERATE OVER ALL COORDINATES (AXES)
        for i in range(len(r)): # FIXME: SCALING SHOULD BE APPLIED ndims TIMES
          
          # DISTANCE ALONG THIS AXIS
          dr = abs(r[i] - src_xyz[i])
          #print this_func, 'dr', dr
          
          
          # HICKS AMPLITUDE
          ampl = ampl * Windowed(dr, Sinc, Window_Kaiser, **kwargs) #Sinc(dr) 
          #ampl = ampl * Sinc(dr) 
          # Kaiser_Sinc_1d(b_kaiser, r_kaiser, bessel, x)
        
        x, y, z = [int(i) for i in [x, y, z]] # NOTE: MAKES LIFE EASIER
        source_nodes.append([x, y, z, ampl]) 
  
  source_nodes = np.array(source_nodes)
  
  #print this_func + 'END' 
  return source_nodes


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


def Cube_Split(cube, flags, **kwargs):
  """
  Split the source cube into 2 parts:
  above and below FS respectively.
  
  Parameters
  ----------
  cube : list 
  
  flags : list
  
  Returns
  -------
  cube_below : list
  
  cube_above : list
  
  Notes
  -----
  We don't want to add to wavefield values 
  exactly at FS (because they should be 0)
  => we just forget about them.
  
  """
  this_func = this_lib + 'Cube_Split: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_generic_CONST import in_flag, acc_flag, ext_flag
  
  cube_below, cube_exact, cube_above = [], [], []
  
  for cn in cube:
    ix, iy, iz = [int(i) for i in np.array(cn[ :3]) - 1]
    
    # GET THE GHOST FLAG CATCHING A POSSIBLE OUT-OF-BOUNDS ERROR
    try:
      flag = flags[ix][iy][iz]
    except IndexError:
      print(this_func, 'cn', cn)
      print(this_func, 'ix, iy, iz', ix, iy, iz)
      raise ValueError('Node outside the grid or not int.' + '\n') 
    
    #print this_func, 'flag', flag
    # SPLIT INTO RELEVANT SUBSETS
    if flag == in_flag:
      cube_below.append(cn)
    
    elif flag == ext_flag:
      cube_above.append(cn)
    
    elif flag == acc_flag:
      #eprint(this_func + 'Cube-node is a FS-node ghost.\n')
      #print this_func, 'FS-node cn', cn
      cube_exact.append(cn)
      #quit()
    
    else:
      print(this_func, 'cn', cn)
      raise ValueError('Unknown flag: ' + flag + '\n')
    
  #Plot_3D(scatts=[cube_below, cube_exact, cube_above], zflip=1, **kwargs)
  #P_Slices_XYZ(scatts=[cube_below, cube_exact, cube_above], **kwargs)#, **kwargs)    
    
  #print this_func + 'END'
  return cube_below, cube_exact, cube_above


# -------------------------------------------------------------------------------


def Wavefield_Update(w, node_srcs, **kwargs):
  """
  Add amplitudes of nodal-sources
  to the wavefield.
  
  Parameters
  ----------
  w : array 
    Wavefield before the update. 
  srcs : list 
    Each is (x, y, z, amplitude).
  
  Returns
  -------
  w : array 
    Wavefield after the update.
  
  Notes
  -----
  
  """  
  this_func = this_lib + 'Wavefield_Update: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  qc = Kwarg('qc', 0, kwargs)

  dw = np.zeros(w.shape) # JUST FOR QC
  
  for src in node_srcs:
    ix, iy, iz = np.array(src[ :3]) - 1
    ix, iy, iz = int(ix), int(iy), int(iz)
    ampl = src[3]
    
    try:
      dw[ix][iy][iz] += ampl
    except IndexError:
      pass
      #eprint('Source outside the wavefield grid.')
  
  w_prev = np.copy(w) # QC 
  w += dw
  
  if qc > 0:
    print(this_func, 'min(w), max(w)', np.min(w), np.max(w))
    fig, ax, i = Subplots(1, 3)
    i = Subplot(fig, i)
    Plot(dw, cmap='seismic')
    i = Subplot(fig, i)
    Plot(w0_prev, cmap='seismic')
    i = Subplot(fig, i)
    Plot(w0, cmap='seismic')

  #print this_func + 'END'
  return w


# -------------------------------------------------------------------------------


def Reflect_Beneath_FS(cube_above, ghosts, intersects, extn, two_d, **kwargs):
  """
  Reflect beneath the FS this part of the source cube which
  sticks out above it.
  
  Parameters
  ----------
  two_d
  # FIXME: ONLY TMP
  
  Returns
  -------
  cube_refl : list
    List of vectors [x, y, z, A(x, y, z)]
 
  Notes
  -----
  
  """
  this_func = this_lib + 'Reflect_Beneath_FS: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_generic_CONST import normal_flag
  
  extn1, extn2, extn3 = extn
  
  i = 0
  cube_refl = []
  for cn in cube_above:
    ix, iy, iz = np.array(cn[:3]) - 1
    G = np.array(cn[:3]) # cube_above => cn IS ONE OF THE GHOSTS
    
    # FIND AN INDEX OF THE COLLOCATED GHOST AND ASSOCIATED INTRESECT. 
    found = False
    ig = -1
    for gh in ghosts: # BOTTLENECK - BUT HOW TO WORK AROUND? 
      ig += 1
      
      # SKIP GHOSTS WITH NO INTERSECT.
      gh_flag = gh[-1] 
      if gh_flag != normal_flag: 
        continue
      
      # PICK A GHOST WITH THE SAME COORDINATES
      if gh[0] == cn[0] and gh[1] == cn[1] and gh[2] == cn[2]:
        found = True
        #print this_func, 'found'
        break
    
    # ESSENTIAL # FIXME: DOUBLE-CHECK IT
    if not found:
      #eprint(this_func + 'No normal ghost found for this cube node.\n')
      #print this_func, 'cn', cn
      continue # => THIS EXTERIOR NODE IS NOT A NORMAL GHOST => HAS NO INTERSECT
    
    # PICK INTERSECT ASSOCIATED WITH THIS GHOST
    I = np.array(intersects[ig]) 
    
    # REFLECT THE GHOST USING THE INTERSECT
    R = 2 * I - G 
    
    if two_d: # FIXME: ONLY TMP
      R[1] = int(R[1])
    
    # DEBUGG.
    extn = [extn1, extn2, extn3]
    for i, q in enumerate(R):
      if (q < 1) or (q > extn[i]): # ADD UPPER TOO
        eprint(this_func + 'Error. R out of bounds\n')
        print(this_func, 'G', G)
        print(this_func, 'I', I)        
        print(this_func, 'R', R)
        quit()
        
    # OUTPUT
    new_ampl = -cn[3]
    cube_refl.append(list(R) + [new_ampl])
  
    i += 1
  
  if 0:
    fig, ax, i = Subplots(1, 2)
    i = Subplot(fig, i)
    P_Slice(scatts=[cube_above], ax=plt.gca(), yflip=1)
    i = Subplot(fig, i)
    P_Slice(scatts=[cube_refl], ax=plt.gca(), yflip=1)

  #print this_func + 'END'
  return cube_refl


# -------------------------------------------------------------------------------


def Update_Sources_List(old_srcs, new_srcs, **kwargs):
  """
  Append new sources and/or update values 
  of the old ones.
  
  Parameters
  ----------
  
  Returns
  -------
  old_srcs : list 
    Old sources updated based 
    on new ones.
  
  Notes
  -----
  
  """
  this_func = this_lib + 'Update_Sources_List: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  for nsrc in new_srcs:
    nx, ny, nz, nA = nsrc
    
    # CHECK IF THE nsrc IS ALREADY AMONG OLD ONES
    found = False
    for i, osrc in enumerate(old_srcs):
      ox, oy, oz, oA = osrc  
      if nx == ox and ny == oy and nz == oz:
        found = True
        break
    # IF SO, UPDATE AMPLITUDE OF THE EXISTING SOURCE
    if found: 
      old_srcs[i][-1] += nA
    # IF NOT, ADD A NEW SOURCE
    else:  
      old_srcs.append(nsrc)
      
  #print this_func + 'END'
  return old_srcs


# -------------------------------------------------------------------------------


def Converged(hsrcs, thresh, orig_src_amp, **kwargs):
  """
  Checks if the spreading the point source already 
  converged.
  
  Parameters
  ----------
  hsrcs : list 
    List of Hicks sources of the form: [x, y, z, A(x,y,z)].
  thresh : float 
    Maximum amplitude for which convergence is 
    regarded to be reached.
  orig_src_amp : float 
    Amplitude of the original point Hicks source to spread.
  
  Returns
  -------
  converged : bool 
    True if ALL the amplitudes of Hicks sources
    fall below the thresh. 
  
  Notes
  -----
  If any amplitude exceeds the original one - raise error.
  
  """
  this_func = this_lib + 'Converged: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')  
  
  from lib_generic_CONST import epsi
  
  converged = True 
  summ = 0
  maxx = epsi
  for hsrc in hsrcs:
    ampl = abs(hsrc[3])
    summ += ampl
    if ampl > maxx:
      maxx = ampl
    if ampl > thresh:
      converged = False
      if ampl > abs(orig_src_amp):
        eprint(this_func + 'Error. Ampl > orig_src_amp')
        quit()
  
  print(this_func, 'max ampl of hsrcs', maxx)
  print(this_func, 'sum of ampl of hsrcs', summ) 
  
  if converged:
    print('\n-----> ', this_func, 'Converged. \n')
  else:
    print(this_func, 'Not yet converged...')
  
  #print this_func + 'END'
  return converged
  

# -------------------------------------------------------------------------------
# UTILS
# -------------------------------------------------------------------------------


def Flags_Read(proj_name, **kwargs):
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Flags_Read: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  
  #print this_func, 'in_flag, acc_flag, ext_flag', in_flag, acc_flag, ext_flag
  
  fname = proj_name + '-InextNodes.txt'
  
  c = Read_File(fname)
  en1, en2, en3 = c[0]
  en1, en2, en3 = [int(i) for i in [en1, en2, en3]]
  data = c[1: ]
  
  flags = np.ones((en1, en2, en3)) * 33333 # FOR DEBUGG.
  #for record in c[1: ]:
    #print record
  i = 0
  for x in range(en3):
    for y in range(en2):
      for z in range(en1):
        #print this_func, 'data', data#[i]
        flag = int(data[i][-1])
        flags[x, y, z] = flag
        i += 1
  
  if verbos > 4:
    from lib_generic_PLOTT import Plot_Slices_XYZ
    Plot_Slices_XYZ(vols=[flags], minn=-2, maxx=2)
 
  #print this_func + 'END'
  return flags


# -------------------------------------------------------------------------------


def Plot_Points_4D(points, **kwargs): # FIXME: DEL/MOVE?
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Plot_Points_4D: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')

  plot_type = Kwarg('plot_type', 'scatter', kwargs)
  
  if plot_type == 'scatter':
    #for p in points:
    print()
      
    
    
  else:
    raise ValueError('Plot type ' + plot_type + ' not yet implemented\n')
  

  #print this_func + 'END'
  return 0


# -------------------------------------------------------------------------------


def Nodes_Flag(w0, ghosts, **kwargs): # FIXME: DEL?
  """
  
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  0
  
  Notes
  -----
  
  
  """
  this_func = this_lib + 'Nodes_Flag: '
  verbos = Kwarg('verbos', 1, kwargs)
  if verbos == verbos_func:
    print(this_func + 'START')
  
  from lib_generic_CONST import normal_flag, fs_node_flag, margin_flag # GHOSTS
  from lib_generic_CONST import in_flag, acc_flag, ext_flag # GRID NODES
  
  
  if len(w0.shape) != 3:
    raise ValueError('Ndims != 3 not implemented')  
  
  
  #in_flag, acc_flag, ext_flag = -100, 0, 100 #FIXME: TMP
  #if verbos > 2:
  #  print this_func, 'Ghosts flags: normal_flag, fs_node_flag, margin_flag', normal_flag, fs_node_flag, margin_flag
  #  print this_func, 'Node flags: in_flag, acc_flag, ext_flag', in_flag, acc_flag, ext_flag
  #
  #flags = np.ones(w0.shape) * (3333) # FIXME
  #
  #for g in ghosts:
  #  gx, gy, gz, gw, glvl, gflag = g
  #  gx, gy, gz = [int(i)-1 for i in [gx, gy, gz]]
  #  if glvl != 1:
  #    continue
  #  for z in range(w0.shape[-1]):
  #    if z <= gz:
  #      flag = ext_flag
  #    else:
  #      flag = in_flag
  #    
  #    print 'gx, gy, z', gx, gy, z
  #    flags[gx][gy][z] = flag
      
  
  #for x in range(len(w0)):
  #  for y in range(len(w0[0])):
  #    for z in range(5):#len(w0[0][0])):
  #      
  #      # FOR EACH GRID NODE SEARCH IN THE GHOSTS
  #      for g in ghosts:
  #        print this_func, 'x, y, z', x, y, z, 'g', g
  #        gx, gy, gz, gw, glvl, gflag = g
  #        gx, gy, gz = [int(i)-1 for i in [gx, gy, gz]]
  #        #print this_func, 'gflag', gflag
  #        if (x == gx) and (y == gy) and (z <= gz):
  #
  #          if gflag == fs_node_flag:
  #            flag = acc_flag
  #            print this_func, 'flag = acc_flag'
  #            return
  #          else:
  #            flag = ext_flag
  #            eprint(this_func + 'flag = ext_flag for ' + str(x) + ',' + str(y) + ',' + str(z) + '\n')
  #            #return
  #          #break
  #        
  #            #print this_func, 'flag = ext_flag'
  #          #if (x == gx) and (y == gy) and (z < gz):
  #          #  xyz = '(' + str(x) + ',' + str(y) + ',' + str(z) + ')' 
  #          #  eprint('Node ' + xyz + ' above ghost layers.\n')
  #
  #        else:
  #          flag = in_flag
  #          print this_func, 'flag = in_flag'
  #          #break
  #          
  #        return # FIXME  
  #        #else:
  #          #eprint('Node not covered laterally by ghosts\n')
  #          #flag = -666
  #          
  #      flags[x][y][z] = flag

  from lib_generic_PLOTT import Plot_Slices_XYZ          
  Plot_Slices_XYZ(vols=[flags], minn=in_flag, maxx=ext_flag, cmap='seismic')
  
  #print this_func + 'END'
  return flags


# -------------------------------------------------------------------------------


