"""
Generic plotters and formatters.

(c) 2019- Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
# from autologging import logged, traced

class Plotter(object):
  """
  Generic plotting mix-in. The plott method
  wraps the plot method of child classes.
  The rationale is that _initalize, _finalize and _save should be the same,
  regardless of the object properties to be plotted.
   
  If really needed in special cases, fine-grained customization 
  should done in the plot method of child classes. Note that they can 
  put the both 
  """
  def _initialize(self, **kwargs):
    """
    Set everything that needs to be set before the actual plotting
    function (imshow etc.) is called.
    """
    figure(**kwargs)
    return kwargs
  # -----------------------------------------------------------------------------
  def _finalize(self, grid={}, **kwargs):
    """
    Add final formatting.
    """
    assert isinstance(grid, dict)
    if len(grid) > 0:
      plt.grid(**grid)
  # -----------------------------------------------------------------------------
  def overlay(self, *layers, **kwargs):
    """
    Actually, we should develop this one as last, since it can 
    be achieved easily in the notebook.
    """
    self.plott(**kwargs)
    for layer in layers:
      layer.plot(**kwargs) # not plott, otherwise new figure
  # -----------------------------------------------------------------------------
  def plot(self, **kwargs):
    raise NotImplementedError('Overwritten in a child class.')
  # -----------------------------------------------------------------------------
  def plott(self, *args, **kwargs):
    """
    We use a different name (plott) to stay compatible with 
    plot methods of child classes.
    """
    save = kwargs.get('save', True)
    kwargs = self._initialize(**kwargs)
    ax = self.plot(*args, **kwargs)
    ax = self._finalize(**kwargs)
    if save:
      self._save(**kwargs)
  # -----------------------------------------------------------------------------
  def _save(self, fname, fmt='png', close=True, **kwargs):
    from fullwavepy.generic.parse import strip
    fname = strip(fname)
    fname = '%s.%s' % (fname, fmt)
    plt.savefig(fname)
    if close:
      plt.close()
  # -----------------------------------------------------------------------------
class FilePlotter(Plotter):
  """
  Generic file-plotting mix-in. Although mix-ins shouldn't inherit
  from any other class, it seems to be cleaner in this case, hence
  the exception. These two are the only options (either plot the data
  directly, or from the file).

  """
  def _initialize(self, **kwargs):
    """
    Here because currently we set title within imshow but it should
    change (-> finalize).
    """
    assert hasattr(self, 'name') # NOTE name, not fname (would be too long)
    kwargs['title'] = kwargs.get('title', self.name)
    return super()._initialize(**kwargs)

  # -----------------------------------------------------------------------------
  
  def plot(self, **kwargs):
    """
    This is a generic workflow for any file containing plottable data.
    We should self.read everytime, see its docstring for rationale.
    Note, it's plot not plott, so when plott is called it will always 
    lead to Plotter first (no clash).
    
    """
    data = self.read(**kwargs)
    return data.plot(**kwargs) 
  
  # -----------------------------------------------------------------------------  
  
  def _save(self, *args, **kwargs):
    """
    Save the plot using the file's original 
    name, just changing the extension.

    Parameters
    ----------
    *args 
      only to comply with the Liskov principle
    """
    assert hasattr(self, 'fname')
    super()._save(self.fname, **kwargs)

  # -----------------------------------------------------------------------------
class FileListPlotter(object):
  """
  """
  def plot(self, *args, **kwargs):
    self._plot('plot', *args, **kwargs)
  def plott(self, *args, **kwargs):
    self._plot('plott', *args, **kwargs)    
  def _plot(self, function, max_files_no=None, **kwargs):
    """
    """
    self.__log.debug('max no. of files allowed: %s' % max_files_no)
    
    for i, f in enumerate(self._all_files()):
      if (max_files_no is not None) and (i == max_files_no): 
        break # when put here, 0 value works as well
      
      self.__log.info('Plotting %s' % f.fname)
      getattr(f, function)(**kwargs)
  
  # -----------------------------------------------------------------------------

def compare(arg1, arg2, **kwargs):
  """
  Compare two arrays of data.
  
  Parameters
  ----------

  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  None
  
  Notes
  -----
  Data-pieces interact with each other. 
  They need to be known at the same time
  e.g. to color-fill the area between them.
  
  """
  from fullwavepy.ndat.arrays import Arr3d
  
  A1 = Arr3d(arg1)
  A2 = Arr3d(arg2)
  
  if A1.shape != A2.shape:
    raise ValueError('Arrays must have identical shapes.')
  
  ndims = len(A1.shape)  
  
  if ndims == 1:
    from fullwavepy.plot.plt1d import plot_1d
    plot_1d(lines=[A1, A2], **kwargs)
  
  elif ndims == 2:
    from fullwavepy.plot.plt2d import compare_2d
    compare_2d(A1, A2, **kwargs)
  
  else:
    raise ValueError('Wrong array shape: %s' % ndims)

def cm2inch(value, unit):
  inch = 2.54 # cm
  if unit == 'cm':
    # convert to inches
    value /= inch
  elif unit == 'inch':
    pass # because mpl expects inches
  elif unit != 'inch':
    raise ValueError('unit: %s'  % unit)  
  return value
def figure(figsize_x=6, figsize_y=6, unit='inch', **kwargs):
  """
  Apparently one has to create a new figure INSIDE 
  a function passed to interact. 
  This is the code that has to be put in every 
  function decorated with ##@widgets then.
  """
  figsize_x = cm2inch(figsize_x, unit)
  figsize_y = cm2inch(figsize_y, unit)  
  figsize = (figsize_x, figsize_y)
  return plt.figure(figsize=figsize)
def figax(figsize_x=15, figsize_y=5, unit='inch', **kwargs):
  """
  A4	210 x 297 mm,	8.3 x 11.7 inches
  """
  figsize_x = cm2inch(figsize_x, unit)
  figsize_y = cm2inch(figsize_y, unit)
  figsize = (figsize_x, figsize_y)
  fig, ax = plt.subplots(figsize=figsize)
  return fig, ax 

def flipy(ax=plt.gca(), **kwargs):
  ax.invert_yaxis()
def autect(ax=plt.gca(), **kwargs):
  ax.set_aspect('auto')
def aspeqt(ax=plt.gca(), **kwargs):
  ax.set_aspect('equal')

def set_xlabels(labels, decim_xlabels=10, rotate_xlabels=None, **kwargs):
  """
  Decimate and rotate labels.

  labels : list 
    Labels before decimation
  """  
  locs = np.arange(len(labels))[::decim_xlabels]
  labels = labels[::decim_xlabels]
  rotation = lambda decim : np.clip(90 - 10 * (decim - 1), 0, 90)
  if rotate_xlabels is None:
    rotate_xlabels = rotation(decim_xlabels)
  set_xlabels._log.debug('Rotating xlabels %s degrees' % rotate_xlabels)
  locs, labels = plt.xticks(locs, labels, rotation=rotate_xlabels)
  return locs, labels
