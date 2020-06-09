"""
This module defines decorators.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from functools import wraps

from fullwavepy.generic.parse import kw


# def save(func):
#   """
#   Save a
  
#   """
#   @wraps(func)
#   def wrapper_save(*args, **kwargs):
#     save = kwargs.get('save', True)
#     value = func(*args, **kwargs)
    
#     if save:
#       print('Saving.')
#       plt.savefig(fname, format='png') # HOW TO ACCESS fname??? FIXME
#       plt.close()
#     return value
#   return wrapper_save


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
      t = end_time - start_time
      print('Exec. time: {} s <{}>'.format("{:15.12f}".format(t), func.__name__))    
    return value
  return wrapper_timer


# -------------------------------------------------------------------------------


@traced
@logged
def widgets(*widgets_args):
  """
  A wrapper around an actual decorator (see below)
  to make it take arguments other than function's ones.
  
  Notes
  -----
  This is apparently a canonical way of passing
  args to decorators.
  
  """
  from ipywidgets import interact, fixed
  @traced
  @logged
  def widgets_actual_decorator(func):
    """
    Rationale: cannot decorate class methods with @interact due to:
     > ValueError: cannot find widget or abbreviation for argument: 'self'
    
    """
    @traced
    @logged
    def wrapper_widgets(*args, **kwargs):
      """
      Notes
      -----
      widgets=fixed(True) is passed to func 
      to force func to create a new figure, 
      otherwise it's not interactive.
      Default values are set in the functions themselves.
      
      """
      from ipywidgets import (IntSlider, BoundedIntText, Dropdown, 
                              SelectMultiple, Checkbox,
                              Layout, TwoByTwoLayout) 
      widgets = kw('wdg', False, kwargs)
      
      try:
        proj = args[0].proj
        sids = sorted(list(proj.i.s.d.keys()))
        #it_max = proj.i.rnf.nits_total
      except AttributeError as err:
        wrapper_widgets._log.warning('Setting sids to [] because of %s' % str(err))
        sids = []
        
      ##print('wow', args[0]) USE THIS TO ACCESS PROJECT METADATA! FIXME
      
      # NOTE: WE ARE SKIPPING *args - THEY WOULD BREAK interact!
      def ifunc(**kwargs): 
        return func(*args, **kwargs) 
      
      interact_kwargs = {
        'figsize_x' : IntSlider(value=8, min=1, max=20, step=1), #layout=Layout(width='90%')),
        'figsize_y' : IntSlider(value=8, min=1, max=20, step=1),
        'cmap'      : Dropdown(options=['twilight','cividis','seismic']+plt.colormaps()),
        'slice_at'  : Dropdown(options=['y', 'x', 'z']),
        'node'      : BoundedIntText(value=0, min=0, max=100, step=5),
        'x'         : BoundedIntText(value=0, min=0, max=100, step=5),
        'y'         : BoundedIntText(value=0, min=0, max=100, step=5),
        'z'         : BoundedIntText(value=0, min=0, max=100, step=5),
        'true_vp'   : Checkbox(True),
        'bathy'     : Checkbox(True),
        'freesurf'  : Checkbox(True),
        'sources'   : Checkbox(True),
        'receivers' : Checkbox(True),
        'sids'      : SelectMultiple(options=sids, value=sids),
        'run_ids'   : SelectMultiple(options=range(20), value=[0]),
        'it'        : BoundedIntText(value=1, min=0, max=500, step=5)
      }
      
      # CHUCK AWAY ALL kwargs THAT ARE NOT LISTED IN widgets_args
      interact_kwargs = dict({(i, interact_kwargs[i]) for i in interact_kwargs.keys() if i in widgets_args})
      
      if widgets:
        interact(ifunc, widgets=fixed(True), **interact_kwargs)
      
      else: 
        return func(*args, **kwargs) # RETURN ONLY FOR widgets=False
    
    return wrapper_widgets
  return widgets_actual_decorator


# -------------------------------------------------------------------------------

