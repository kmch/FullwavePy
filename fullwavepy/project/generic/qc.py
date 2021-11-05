"""
This module contains tools for quality control of FWI runs.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, strip
from fullwavepy.generic.decor import widgets


@traced
@logged
class Functional(object):
  """
  """
  def __init__(self, proj, **kwargs):
    """
    
    """
    self.proj = proj
    self.path = self.proj.out.path
    #self.name = self.proj.name + '-Functional.png' # DUMP ASCII INSTEAD?
    #self.fname = self.path + self.proj.name + '-Functional.png'
  
  # ----------------------------------------------------------------------------- 

  def _prep_widgets(self, **kwargs):
    """
    Another take on robust, customized, automated
    interactive plotting.

    """
    from ipywidgets import IntSlider, BoundedIntText, Dropdown, \
                           SelectMultiple, Checkbox, Layout, TwoByTwoLayout
    widgets = {}
    p = self.proj
    
    sids = [s.ID for s in p.i.s.read().li]
    widgets['sids'] = SelectMultiple(options=sids, value=sids)
    
    self.widgets = widgets
    return widgets

  # -----------------------------------------------------------------------------    
  
  def read(self, run_ids, misfit=True, **kwargs):
    """
    Get value of fit for all sources and for 
    all iterations.
    
    Parameters
    ----------   
    
    Returns
    -------
    functional : dict 
      Dictionary with source-IDs as keys, 
      every key contains a list of fit (%) for 
      all iterations.
    
    Notes
    -----_
    PLOT FIT AS ~ SIZE, OR COLOR (TIM LIN, PHASE PLOTS)? 
    
    
    # timestamp TRUE IF EVERY LINE OF OUTPUT LOG STARTS WITH TIME (SLAVES_SHOWTIMESTAMP="yes")
    
    """
    from fullwavepy.ioapi.generic import read_txt_raw
    
    #timestamp = kw('timestamp', self.proj.env.var['slave_timestamp'],
                      #kwargs) # NOT SCHEDULER? FIXME
    
    #self.__log.debug('timestamp: ', timestamp)
    
    if ((self.proj.env.var['SLAVES_SHOWTIMESTAMP'] == 'yes') or 
        (self.proj.env.var['SLAVES_SHOWTIMESTAMP'] == 'on')):
      add_to_index = 1
    else: 
      add_to_index = 0
    
    functional = {}
    c = self.proj.out.out.read(run_ids) 
    for line in c.splitlines():
      split = line.split(None)
      # PARSE LINES CONTAINING FIT INFO
      if len(split) > 5 and split[1 + add_to_index] == 'calcResidsInfo:':
        #print split
        
        # SLAVE NO.
        first = split[0 + add_to_index]
        slave_no = first.split('(Slave')[1]
        slave_no = slave_no.split(';')[0]
        #print 'slave_no', slave_no 
        
        # SOURCE ID
        sid = first.split('CSRef')[1]
        sid = int(sid.split(';')[0])
        self.__log.debug('Converting shot ids into integers')
        #print 'sid', sid
        
        # FIT [%] #NOTE: D-C THAT'S THE CORRECT VALUE
        percent = float(split[4 + add_to_index][:-1]) # 4TH WORD OF LINE WITHOUT % SIGN
        
        if misfit:
          percent = 100 - percent
        
        # INITIALIZE A LIST FOR ALL ITERATIONS IF NOT DONE IT YET
        if not sid in functional:
          functional[sid] = []
        
        # APPEND THIS ITERATIONS (ACTUALLY TWICE EVERY ITERATION, SEE BELOW)
        functional[sid].append(percent)
    
    #NOTE: ESSENTIAL: Out.log CONTAINS 2 PIECES OF FIT INFO FOR EACH ITERATION
    # MOST PROBABLY THE SECOND ONE IS AFTER BACKPROPAGATION?! ANYWAY, WE GONNA GET 
    # RID OF IT FOR NOW BY DECIMATING LIST OF FITS (TAKING EVERY SECOND) SO THAT 
    # WE HAVE 1 NUMBER PER SOURCE PER ITERATIONS
    for key in functional:
      functional[key] = functional[key][::2]
      #print key, functional[key]
    
    # ADD TO BLOCK INFORMATION (FIXME: WE WOULD NEED TO PARSE MORE INFO, TO 
    # DISTINGUISH BETWEEN DIFFERENT BLOCKS, NO NEED FOR THIS, MAYBE IN THE FUTURE
    #for block in self.proj.inp.runfile.blocks:
    #  for iteration in range(block.niters):
    #    print 'a'
    #quit()  
    
    self.shot = functional 
    
    return self.shot
  
  # ----------------------------------------------------------------------------- 
  
  #@widgets('sids', 'run_ids')
  def plot(self, widgets=False, cmap='rainbow', **kwargs):
    """
    
    run_ids : list 
      Not every run has an Out.log parsable or meaningful.
      Choose which ones have.
    
    """
    from fullwavepy.plot.plt1d import colors
    
    misfit = kw('misfit', True, kwargs)   
    alpha = kw('alpha', 0.4, kwargs)
    ls = kw('ls', '.-', kwargs)
    lw = kw('lw', 3, kwargs)
    xlim = kw('xlim', None, kwargs)
    ylim = kw('ylim', (0,100), kwargs)
    
    srcs = self.proj.i.s.read(**kwargs)
    sids = kw('sids', [s.ID for s in srcs.li], kwargs)    
    kwargs['run_ids'] = kw('run_ids', [0], kwargs)
    misfit = kw('misfit', True, kwargs)
    
    functional = self.read(**dict(kwargs, misfit=misfit))
    
    # CONVERT KEYS TO INT TO SORT THEM
    functional = {int(k) : v for k, v in functional.items()}
    self.__log.debug('Converting shot IDs into integers')

    if sids is not None:
      functional = {sid : functional[sid] for sid in sids}

    clrs = colors(len(functional), cmap=cmap)
    
    for sid, fit in sorted(functional.items()):
      plt.plot(list(range(1, len(fit) + 1)), fit, ls, lw=lw, label=sid, c=next(clrs), 
               alpha=alpha)
    
    ax = plt.gca()
    ax.set_xlabel('Iteration')
    if misfit:
      ax.set_ylabel('Trace misfit [%]')
    else:
      ax.set_ylabel('Trace fit [%]')
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if len(functional) < 50:
      plt.legend(prop={'size': 6})
    
    return ax
    
  # -----------------------------------------------------------------------------    

    
# -------------------------------------------------------------------------------


#class StepLength(object):


# ------------------------------------------------------------------------------- 


@traced
@logged
class JobStats(object):
  """
  Statistics of code performance 
  in subsequent runs.
  
  """
  def __init__(self, proj, **kwargs):
    """
    """
    self.proj = proj
    #self.path = self.proj.out.path

  # -----------------------------------------------------------------------------
  
  def read(self, run_ids, **kwargs):
    """
    """
    stats = {}
    for run_id in run_ids:
      stats[run_id] = self.proj.out.jobout.no[run_id].read(**kwargs)
    
    return stats

  # -----------------------------------------------------------------------------

  def plot(self, run_ids, **kwargs):
    """
    See JobOutLog
    resources = {'mem_req': [],
                 'mem_use': [],
                 'cpu_req': [],
                 'cpu_use': [],
                }    
    """
    stats = self.read(run_ids, **kwargs)
    print('stats {}'.format(stats))    
    for key in stats[run_ids[0]].keys():
      plt.plot(run_ids, [stats[run_id][key] for run_id in run_ids]) 
     
  # -----------------------------------------------------------------------------  


# -------------------------------------------------------------------------------

