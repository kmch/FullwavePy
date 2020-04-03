"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for peremoveission writing to k.chrapkiewicz17@imperial.ac.uk.

Notes
-----
Legacy code, possibly to be chucked away for good.

"""
@timer
@traced
@logged
def read_arrays(*args, **kwargs):
  """
  Read data either from arrays
  or files. Both types can appear
  in *args.
  
  Parameters
  ----------
  *args : see below 
    List of string/arrays.
    If string, it is assumed to 
    stand for a file name (incl.
    the path if file is outside './'),
    otherwise it must be an array.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  arrays : list
    List of read arrays.
  
  Notes
  -----
  You can mix files and arrays defined ad hoc
  in the notebook which may come in handy.
  
  """
  from fullwavepy.ndat.arrays import slice_array, modify_array
  
  arrays = []
  for arg in args:
    if isinstance(arg, str):
      A = read_any(arg, **kwargs)
    elif type(arg) == type(np.array([])) or type(arg) == np.memmap:
      A = arg
    else:
      raise TypeError('Arguments need to be either ' + 
                      'file-names or arrays or np.memmap.')
      
    A = slice_array(A, **kwargs)  
    A = modify_array(A, **kwargs)
    read_arrays._log.debug('Min of the array after modifs: ' + str(np.min(A)))
    read_arrays._log.debug('Max of the array after modifs: ' + str(np.max(A)))
    arrays.append(A)
  return arrays


# -------------------------------------------------------------------------------


@traced
@logged
def slice_array(A, **kwargs):
  """
  
  Parameters
  ----------
  A : array 
   3D array, although other shapes
   should be handled too (work in progress).
   
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  2D or 1D array (special case).
  
  Notes
  -----
  
  """
  scoord = kw('scoord', 'y', kwargs)
  dims = A.shape

  # 1D ARRAY CONTAINING A SINGLE 1D TIME-SERIES
  if len(dims) == 1:
    slice_array._log.warn('1D array detected. Passing it on intact')
    return A
  
  # 3D ARRAY CONTAINING A SINGLE 1D TIME-SERIES
  if (dims[0] == 1) and (dims[1] == 1):
    return A[0][0]
  
  # SURFACES
  if (dims[-1] == 1):
    return A[..., 0]
  
  if scoord == 'x':
    svalue = kw('svalue', len(A)//2, kwargs)
    A = A[svalue]
  
  elif scoord == 'y':
    svalue = kw('svalue', len(A[0])//2, kwargs)
    A = [i[svalue] for i in A]
  
  elif scoord == 'z':
    svalue = kw('svalue', len(A[0][0])//2, kwargs)
    A = [[j[svalue] for j in i] for i in A]
  
  elif scoord is None:
    slice_array._log.warn('No slicing applied')
  
  else:
    raise ValueError('Wrong slice coord: ' + scoord)
  
  if scoord is not None:
    slice_array._log.info('Sliced at svalue=' + str(svalue) + ' node ' + 
                          'of scoord=' + scoord + ' axis')
  
  A = np.array(A)
  
  return A
  

# -------------------------------------------------------------------------------

@traced
@logged
def plot(*args, **kwargs):
  """
  A framework to plot (any number of) 
  arrays (possibly from files). 
  
  Parameters
  ----------
  *args : see below 
    List of string/arrays.
    If string, it is assumed to 
    stand for a file name (incl.
    the path if file is outside './'),
    otherwise it must be an array.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  None
  
  Notes
  -----
  You can mix files and arrays defined ad hoc
  in the notebook which may come in handy.
  
  There is no 'interaction' between args, one 
  arg is.plotted after another.
  
  """
  from fullwavepy.ioapi.generic import read_arrays
  
  arrays = read_arrays(*args, **kwargs)
  
  for A in arrays:
    plot_array(A, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
def plot_array(A, **kwargs):
  """
  Plot 1D/2D array.
  
  Parameters
  ----------
  A : array
    1D/2D array.

  **kwargs : keyword arguments (optional)
      
  Returns
  -------
  None
  
  Notes
  -----
  Just a simple framework.
  
  """    
  from fullwavepy.plot.oned import plot_1d
  from fullwavepy.plot.twod import plot_2d
  
  ndims = len(A.shape)
  
  if ndims == 1:
    plot_1d(lines=[A], **kwargs)

  elif ndims == 2:
    plot_2d(images=[A], **kwargs)
    
  else:
    raise ValueError('Wrong array shape: %s' % ndims)
    

# -------------------------------------------------------------------------------



@traced
@logged
def plot_image_OLD(image, **kwargs):
  """
  
  Parameters
  ----------
  
  
  **kwargs : keyword arguments (optional)
    Current capabilities: 
    
    cmap : str 
      Add a trailing '_r' to get 
      a reversed version.
      Type plt.colormaps() to list all 
      registered ones.
      
      Sequential (perceptually uniform):
        'inferno', 'magma', 'plasma', 'viridis'
      Sequential (selected):
        'hot;, 'Greys', 'YlGnBu', 'GnBu', 'PuBuGn'
      Diverging:
        'seismic', 'bwr', 'RdBu', RdYlBu', 'Spectra'
      cmocean:
        'cmo.topo' etc.     
      Other:
        'ocean', 'gist_earth', 'cubehelix'
      
      
    shade : bool
      Default: False.
    
    Parsed only if shade=True:
    - shade_azim : float 
      Azimuth of the light source [deg].
      Default: 315.
    - shade_elev : float 
      Elevation of the light source [deg]. 
      Default: 45.
    - shade_mode : str
      Possible values: 'Hillshade', 'hsv', 'overlay', 'soft'.
      Default: 'overlay'.
    - shade_exag : float 
      Vertical exaggaration of shaded topography.
      Default: 10.
      
  Returns
  -------
  None
  
  Notes
  -----
  Mapping data onto colors using a colormap typically involves two steps: 
  a data array is first mapped onto the range 0-1 using a subclass of Normalize
  
  """  
  ncolors = kw('ncolors', None, kwargs)
  shade = kw('shade', False, kwargs)
  cbar = kw('cbar', False, kwargs)
  contour = kw('contour', False, kwargs)
  darken = kw('darken', False, kwargs)
  
  # DETECT CONFLICTS
  if shade and contour:
    raise ValueError('You cannot combine shading and a contoufullwavepy.plot.')
  if darken and contour:
    raise ValueError('You cannot combine darkening and a contoufullwavepy.plot.')  
  if darken and shade:
    raise ValueError('You cannot combine darkening and shading.')  
  
  kwargs['minn'], kwargs['maxx'] = np.min(image), np.max(image)
  imshow_kwargs = _set_cmap(**kwargs)
  imageT = image.T
  
  
  # SET CORRECT Y-AXIS FOR SPECTRAL PLOTS
  if 'spect' in kwargs:
    from fullwavepy.import import dft_freqs
    nsamps, ntraces = imageT.shape
    y = dft_freqs(nsamps, which='positive', **kwargs)
    imageT = np.array(imageT[ :len(y)])
    imshow_kwargs['extent'] = [0, ntraces, y[-1], 0] 
    
    
  # THIS IS ONLY FOR COLORBAR TO DISPLAY VALUES BEFORE NORMALIZATION
  im_raw = plt.imshow(imageT, **imshow_kwargs)  
  ax = plt.gca()
  
  if shade:
    from matplotlib.colors import LightSource
    shade_azim = kw('shade_azim', 315, kwargs)
    shade_elev = kw('shade_elev', 45, kwargs)
    shade_mode = kw('shade_mode', 'overlay', kwargs)
    shade_exag = kw('shade_exag', 10, kwargs)
    ls = LightSource(azdeg=shade_azim, altdeg=shade_elev)
    imageT = ls.shade(image.T, **imshow_kwargs, 
                      blend_mode=shade_mode, vert_exag=shade_exag)
    
  if darken:
    from matplotlib.colors import Normalize
    dark_min = kw('dark_min', minn, kwargs)
    dark_max = kw('dark_max', maxx, kwargs)
    # CREATE AN ARRAY (nx,ny,nz,3) FULL OF VALUE 70 (GREY COLOR)
    # NOTE: '*' UNPACKS A LIST THAT FOLLOWS
    greys = np.full((*image.T.shape, 3), 70)
    alphas = np.array(image.T)
    # SET OPACITY TO 1 IF dark_min <= value <= dark_max 
    # AND TO 0.5 (HALF-TRANSPARENT) OTHERWISE 
    # IT'S DONE WEIRDLY TO MAKE SURE WE REPLACE VALUES CORRECTLY
    unique = 12345.12345
    alphas[np.where((alphas >= dark_min) & (alphas <= dark_max))] = unique
    alphas[alphas != unique] = 0.5
    alphas[alphas == unique] = 1.0 
    # THIS CONVERTS THE NP.ARRAY TO AN RGB IMAGE ADDING
    # AN EXTRA DIMENSION OF SIZE 4 (FOR R,G,B & ALPHA)
    cmap = imshow_kwargs['cmap']
    imageT = cmap(Normalize(imshow_kwargs['vmin'], imshow_kwargs['vmax'])(image.T))
    imageT[..., -1] = alphas # ELLPISIS IS EQUIVALENT TO :,:,: HERE
    plt.imshow(greys)
  
  # NOTE
  if 'extent' in kwargs:
    imshow_kwargs['extent'] = kwargs['extent']
    
  plt.imshow(imageT, **imshow_kwargs)  
  
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im_raw, cax=cax) 
    if ncolors:
      ticks = np.linspace(minn, maxx, ncolors+1)
      cbar.set_ticks(ticks)
      cbar.set_ticklabels(ticks)
    plt.sca(ax) # SET CURRENT AXIS BACK TO THE ACTUAL PLOT
    
  if contour:
    ncontours = kw('ncontours', 5, kwargs)
    contour_color = kw('contour_color', 'w', kwargs)
    contour_width = kw('contour_width', 1, kwargs)
    contour_style = kw('contour_style', 'dotted', kwargs)
    cntr = ax.contour(imageT, ncontours, colors=contour_color, 
                      linewidths=contour_width, linestyles=contour_style)  
    ax.clabel(cntr, cntr.levels)  


# -------------------------------------------------------------------------------


