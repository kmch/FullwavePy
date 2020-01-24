"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from functools import wraps
from ipywidgets import interact, fixed

from fullwavepy.generic.parse import kw



# -------------------------------------------------------------------------------


def timer(func):
  """
  Print the time taken to run 
  the decorated function.
  
  """
  @wraps(func)
  def wrapper_timer(*args, **kwargs):
    import time
    timer = kw('timer', False, kwargs)
    
    start_time = time.perf_counter()
    value = func(*args, **kwargs)
    end_time = time.perf_counter()
    
    if timer:
      print('Function ' + func.__name__ + '() took ' + 
            str(end_time - start_time) + ' s to run')
    
    return value
  return wrapper_timer


# -------------------------------------------------------------------------------


def widgets(func): #, **widget_kwargs):
  """
  
  Notes
  -----
  It is meant to wrap class methods, hence
  the first argument is assumed to be self.
  
  """
  @wraps(func)
  def wrapper_widgets(*args, **kwargs):
    """
    Notes
    -----
    widgets=fixed(True) is passed to func 
    to force func to create a new figure, 
    otherwise it's not interactive.
    
    """
    from ipywidgets import (IntSlider, BoundedIntText, Dropdown, 
                            SelectMultiple, Checkbox,
                            Layout, TwoByTwoLayout) 
    widgets = kw('widgets', False, kwargs)
    
    def ifunc(**kwargs): # SKIP *args - THEY WOULD BREAK interact!
      return func(*args, **kwargs) 
    
    #try:
    #print('self.proj.dims', self.proj.dims)
    #except Eas err:
      #print(err)
    
    interact_kwargs = {
      #'figsize_x' : IntSlider(value=8, min=1, max=20, step=1, 
                              #layout=Layout(width='90%')),
      'figsize_y' : IntSlider(value=8, min=1, max=20, step=1),
      #'cmap'      : Dropdown(options=['twilight_r', cividis','seismic']+plt.colormaps()),
      #'slice'     : Dropdown(options=['y', 'x', 'z']),
      'x'         : BoundedIntText(value=0, min=0, max=100, step=5),
      'y'         : BoundedIntText(value=0, min=0, max=100, step=5),
      'z'         : BoundedIntText(value=0, min=0, max=100, step=5),
      'true_vp'   : Checkbox(True),
      #'bathy'     : Checkbox(True),
      #'freesurf'  : Checkbox(True),
      'sources'   : Checkbox(True),
      'receivers' : Checkbox(True),
    }

    #app = TwoByTwoLayout(top_left=Dropdown(options=['y']), 
                         #top_right=Dropdown(options=['y']),
                         #bottom_left=Dropdown(options=['y']),
                         #bottom_right=Dropdown(options=['y']))
      
    if widgets:
      interact(ifunc, widgets=fixed(True), **interact_kwargs)

    else: 
      # SET 'STATIC' VALUES FROM WIDGETS' DEFAULTS
      for key, widget in interact_kwargs.items():
        #NOTE: WE HAVE TO USE kw TO PASS VALUES FROM WIDGETS
        # TO INNER FUNCTIONS (LIKE plot_image)
        if isinstance(widget, (IntSlider, BoundedIntText, Checkbox)):
          kwargs[key] = kw(key, widget.value, kwargs)
        elif isinstance(widget, Dropdown):
          kwargs[key] = kw(key, widget.options[0], kwargs)
        else:
          raise ValueError(widget)
      
      return func(*args, **kwargs) # RETURN ONLY FOR widgets=False
  
  return wrapper_widgets


# -------------------------------------------------------------------------------

