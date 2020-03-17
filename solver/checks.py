"""
(c) 2019 Kajetan Chrapkiewicz.
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
    print('All good! Courant number C=%s < Cmax=%s' % (C, C_max))
  else:
    raise ValueError('The scheme is unstable for the current discretization\n' + 
                     '(Courant number C=%s >= Cmax=%s' % (C, C_max) + ')\n' + 
                     'To make it stable, do one of the following: \n' +
                     '1. Decrease v_max of the model below: ' + str(C_max * dx / dt) + ' m\n' +
                     '2. Increase grid cell above: ' + str(dt * v_max / C_max) + ' m\n' +
                     '3. Decrease time step below: ' + str(C_max * dx / v_max) + ' s\n')


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
  check_accuracy._log.info('Shortest wavelength for this model: %s m (%s nodes)' % (shortest_wavelength, 
                                                                                    nodes_per_shortest_wavelength))
  check_accuracy._log.info('No. of nodes per (shortest) wavelength: %s' % nodes_per_shortest_wavelength)  
  
  # CHOOSE MIN NODES PER WAVELENGTH FOR A GIVEN SCHEME
  if kernel == 'high':
    min_nodes_per_wavelen = 3.6
  elif kernel == 'low':
    min_nodes_per_wavelen = 5
  else:
    raise NotImplementedError('Unknown kernel: %s' % kernel)
  
  f_allowed = v_min / (min_nodes_per_wavelen * dx)
  
  if f_max < f_allowed:
    check_stability._log.info('All good! f_max=%s < f_allowed=%s' % (f_max, f_allowed))
  else:
    good_f_max = v_min / (min_nodes_per_wavelen * dx)
    raise ValueError('The scheme is not accurate for the current discretization\n' + 
                     '(f_max=%s >= f_allowed=%s' % (f_max, f_allowed) + ')\n' + 
                     'To make it accurate, do one of the following: \n' +
                     '1. Decrease f_max below: ' + str(good_f_max) + ' Hz\n' +
                     '2. Decrease grid cell below: ' + 
                     str(v_min / (min_nodes_per_wavelen * f_max)) + ' m\n' +
                     '3. Increase v_min of the model above: ' + 
                     str(min_nodes_per_wavelen * f_max * dx) + ' m/s\n')   
    
    #eprint(this_func + '  (Decrease f_peak below: ' + str(nf_max / ricker_fmax2fpeak_ratio) + ' Hz)\n')
    
    #eprint(this_func + 'f_allowed = ' + str(f_allowed) + ', f_max = ' + str(f_max) + '\n')
    #eprint(this_func + 'Error! Inaccurate for ' + kernel + ' kernel\n')
    #eprint(this_func + 'Possible solutions which will work independently: \n')
    #nf_max = v_min / (min_nodes_per_wavelen * dx)
    #eprint(this_func + '1. Decrease f_max below: ' + str(nf_max) + ' Hz\n')
    #eprint(this_func + '  (Decrease f_peak below: ' + str(nf_max / ricker_fmax2fpeak_ratio) + ' Hz)\n')
    #eprint(this_func + '2. Decrease grid cell below: ' + str(v_min / (min_nodes_per_wavelen * f_max)) + ' m\n')
    #eprint(this_func + '3. Increase v_min of the model above: ' + str(min_nodes_per_wavelen * f_max * dx) + ' m/s\n')
    #quit()
  

# -------------------------------------------------------------------------------

