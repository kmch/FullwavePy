"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced


# -------------------------------------------------------------------------------


@traced
@logged
def courant(dx, dt, v_max, **kwargs):
  """
  Calculate the Courant's number describing
  stability of the numerical scheme.
  
  Parameters
  ----------
  dx : float 
    Size of the spatial grid cell [m].
  dt : float 
    Time step [s].
  v_max : float
    Max. velocity of the model.
  
  Returns
  -------
  C : float 
    Courant's number used later in 
    CFL criterion C < C_max where 
    C_max is scheme-specific.
  
  Notes
  -----
  The smaller the more stable.
  
  If C > 1 the fastest wave in the model will 
  be covering more than 1 grid cell per time step.
  
  """
  # GRID VELOCITY (1 GRID-CELL PER TIME-STEP) IN PHYSICAL UNITS
  v_grid = dx / float(dt) 
  # RATIO OF MAX AND GRID VELOCITY 
  C = v_max / float(v_grid) 

  return C


# -------------------------------------------------------------------------------


@traced
@logged
def check_stability(dx, dt, v_max, kernel, **kwargs):
  """
  Check stability of the numerical scheme (kernel).
  
  Parameters
  ----------
  dx : float 
    Size of the spatial grid cell [m].
  dt : float 
    Time step [s].
  v_max : float
    Max. velocity of the model.
  kernel : str 
    Type of kernel. Current capabilities:
     - 'low' (fullwave3D)
     - 'high' (fullwave3D)
  
  Returns
  -------
  None
  
  Notes
  -----
  The smaller the more stable.

    Courant's number used later in 
    CFL criterion C < C_max where 
    C_max is scheme-specific.
  
  If C > 1 the fastest wave in the model will 
  be covering more than 1 grid cell per time step.
  
  """  
  C = courant(dx, dt, v_max)
  
  if kernel == 'high':
    C_max = 0.38
  elif kernel == 'low':
    C_max = 0.5
  else:
    raise NotImplementedError('Unknown kernel: %s' % kernel)
    
  if C < C_max:
    check_stability._log.info('\n\nAll good! Courant number C=%s < Cmax=%s' % (C, C_max))
  else:
    raise ValueError('The scheme is unstable for the current discretization\n' + 
                     '(Courant number C=%s >= Cmax=%s' % (C, C_max) + ')\n' + 
                     'To make it stable, do one of the following: \n' +
                     '1. Decrease v_max of the model below: ' + str(C_max * dx / dt) + ' m\n' +
                     '2. Increase grid cell above: ' + str(dt * v_max / C_max) + ' m\n' +
                     '3. Decrease time step below: ' + str(C_max * dx / v_max) + ' s\n')
  return C

# -------------------------------------------------------------------------------


@traced
@logged
def check_accuracy(dx, v_min, f_max, kernel, **kwargs):
  """
  Check accuracy of the scheme (kernel).
  
  Parameters
  ----------
  dx : float 
    Size of the spatial grid cell [m].
  v_min : float 
    Min. velocity of the model [m/s].
  f_max : float
    Max. frequency present in the source function.
  kernel : str 
    Type of kernel. Current capabilities:
     - 'low' (fullwave3D)
     - 'high' (fullwave3D)
  
  Returns
  -------
  None 
  
  Notes
  -----
  
  """    
  check_accuracy._log.debug('dx, v_min, f_max, kernel: %s, %s, %s, %s' % (dx, v_min, f_max, kernel))
  shortest_wavelength = v_min / f_max
  nodes_per_shortest_wavelength = shortest_wavelength / dx
  text  = '\n\n'
  text += str('Shortest wavelength for this model: %s m (%s nodes)' % ("{0:.1f}".format(shortest_wavelength), 
                                                                                    "{0:.1f}".format(nodes_per_shortest_wavelength)))
  # CHOOSE MIN NODES PER WAVELENGTH FOR A GIVEN SCHEME
  if kernel == 'high':
    min_nodes_per_wavelen = 3.6
  elif kernel == 'low':
    min_nodes_per_wavelen = 5
  else:
    raise NotImplementedError('Unknown kernel: %s' % kernel)
  
  f_allowed = v_min / (min_nodes_per_wavelen * dx)
  
  if f_max < f_allowed:
    text += str('\nAll good! f_max=%s < f_allowed=%s' % (f_max, "{0:.1f}".format(f_allowed)))
  else:
    good_f_max = v_min / (min_nodes_per_wavelen * dx)
    raise ValueError('The scheme is not accurate for the current discretization\n' + 
                     '(f_max=%s >= f_allowed=%s' % (f_max, f_allowed) + ')\n' + 
                     'To make it accurate, do one of the following: \n' +
                     '1. Decrease f_max below: ' + str(good_f_max) + ' Hz\n' +
                     '2. Decrease grid cell below: ' + 
                     "{0:.1f}".format(v_min / (min_nodes_per_wavelen * f_max)) + ' m\n' +
                     '3. Increase v_min of the model above: ' + 
                     "{0:.1f}".format(min_nodes_per_wavelen * f_max * dx) + ' m/s\n')   
   
  check_accuracy._log.info(text)
  return nodes_per_shortest_wavelength
  


# -------------------------------------------------------------------------------


@traced
@logged
def check_propag_dists(dims, dx, t, v_min, v_max, f_min, f_max):
  """
  Calculate propagation distance across 
  the model covered by the fastest waves.
  
  Parameters
  ----------
  dims : tuple
  dx : float 
    Size of the spatial grid cell [m].
  t : s

  f_min : float
    Min. frequency present in the source function.
  f_max : float
    Max. frequency present in the source function.
  v_min : float
    Min. vel. of the true model.
  v_max : float
    Max. vel. of the true model.   
  
  Returns
  -------
  dist : float 
    Distance in metres.
    
  Notes
  -----
  
  """    
  nx1, nx2, nx3 = dims

  assert f_min > 0 and f_max > 0
  shortest_wavelength = v_min / f_max
  longest_wavelength = v_max / f_min
  
  d1 = v_max * t # dist, m
  d2 = v_max * t / shortest_wavelength # dist_per_shortest_wavelength 
  d3 = v_max * t / longest_wavelength # dist_per_longest_wavelength
  d4 = v_min * t / dx # dist_in_nodes
  d5 = d4 / nx1 # dist_as_fraction_nx1
  d6 = d4 / nx2 # dist_as_fraction_nx2
  d7 = d4 / nx3 # dist_as_fraction_nx3
  
  text = '\n\n'
  text += 'Assuming t = ' + str(t) + ' s, the fastest wave will cover '
  text += "{0:.1f}".format(d1) + ' m, which corresponds to: \n'
  text += "{0:6.1f}".format(d2) + ' shortest wavelengths \n'
  text += "{0:6.1f}".format(d3) + ' longest wavelengths \n'
  text += "{0:6.1f}".format(d4) + ' nodes \n'
  text += "{0:6.1f}".format(d5) + ' model-sizes in X direction \n'
  text += "{0:6.1f}".format(d6) + ' model-sizes in Y direction \n'    
  text += "{0:6.1f}".format(d7) + ' model-sizes in Z direction \n'    
#.ljust(7, ' ')
  check_propag_dists._log.info(text)


# -------------------------------------------------------------------------------


