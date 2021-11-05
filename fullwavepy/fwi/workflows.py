"""
Ready-to-use FWI workflows.
"""
from abc import ABC, abstractmethod
from autologging import logged, traced

@logged
class Workflow(ABC):
  """
  A generic FWI workflow.
  
  Notes
  -----  
  The canonical sequence of steps is the following:
    0. Initialise synthetic and inversion projects
    1. Prepare input for synthetic calculation.
    2. Run the synthetic calculation.
    3. Plot the output of the synthetic calculation.
    4. Prepare input for inversion.
    5. Run the inversion.
    6. Plot the output of the inversion.  
  """
  def __init__(self, name, path):
    """
    Creates two projects: 
    - one for the synthetic calculation
    - one for the inversion,
    named after the workflow with 'syn' 
    and 'inv' suffixes, respectively.

    Parameters
    ----------
    name : str
        Name of the workflow.
    path : str
        Path where all project directories (!) 
        will reside in.
    
    Notes
    -----
    It does NOT initialise the projects yet.
    In this way, it allows to set all the necessary
    parameters as workflow attributes before initialising
    the projects which will use these attributes.
    
    """
    self.name = name
    self.path = path
    self.syn = self._SynClass(self.name + 'syn', self.path)
    self.inv = self._InvClass(self.name + 'inv', self.path)
  def run(self, steps, skip=[], **kwargs):
    """
    Run the workflow.

    Parameters
    ----------
    steps : list
        List of integers [0,1,...] denoting
        the subsequent steps of the workflow to run.
    skip : list, optional
        List of steps to skip, by default [].
    """
    kwargs = dict(kwargs, cat=0) # mute print-outs of text files
    self.steps = [i for i in steps if not i in skip]
    self._step00(**kwargs)
    self._step01(**kwargs)
    self._step02(**kwargs)
    self._step03(**kwargs)
    self._step04(**kwargs)
    self._step05(**kwargs)
    self._step06(**kwargs)  
  # -----------------------------------------------------------------------------
  @abstractmethod
  def _SynClass(self):
    """
    Define which class to use 
    for a synthetic project.

    Returns
    -------
    class
        Child-class of Pro.
    """
    return Syn
  @abstractmethod
  def _InvClass(self):
    """
    Define which class to use 
    for an inversion project.

    Returns
    -------
    class
        Child-class of Pro.
    """    
    return Invert
  # -----------------------------------------------------------------------------
  def _skip_step(self):
    if self.step_no in self.steps:
      return False 
    else:
      self.__log.debug('Skipping step %s' % self.step_no)
      return True   
  def _step_info(self, info):
    """
    Print information on what will be done
    during this step.

    Parameters
    ----------
    info : string
        Information to be printed.
    
    Notes
    -----
    It is an ordinary print, because log.info gets
    lost in project-logs etc.
    """
    print('Step %s: %s' % (str(self.step_no).rjust(3, '0'), info))  
  def _step00(self, **kwargs):
    self.step_no = 0
    if not self._skip_step():
      self._step_info('Initialising synthetic and inversion projects.')
      self._step00(**kwargs)
  def _step01(self, **kwargs):
    self.step_no = 1
    if not self._skip_step():
      self._step_info('Preparing input for synthetic calculation.')
      self._step01(**kwargs)
  def _step02(self, **kwargs):
    self.step_no = 2
    if not self._skip_step():
      self._step_info('Running synthetic calculation.')
      self._step02(**kwargs)
  def _step03(self, **kwargs):
    self.step_no = 3
    if not self._skip_step():
      self._step_info('Plotting output of synthetic calculation.')
      self._step03(**kwargs)
  def _step04(self, **kwargs):
    self.step_no = 4
    if not self._skip_step():
      self._step_info('Preparing input for inversion.')
      self._step04(**kwargs)
  def _step05(self, **kwargs):
    self.step_no = 5
    if not self._skip_step():
      self._step_info('Running inversion.')
      self._step05(**kwargs)
  def _step06(self, **kwargs):
    self.step_no = 6
    if not self._skip_step():
      self._step_info('Plotting output of inversion.')
      self._step06(**kwargs)
@logged
class GenerateSyn(Workflow):
  """
  Generate synthetic data.
  """
  pass
@logged
class InvertSyn(Workflow):
  """
  Generate and invert synthetic data.

  Notes
  -----
  The input of the inversion is the same as 
  the synthetic run except for relevant runfile parameters
  and the difference defined 'inv_vs_syn'
  parameter passed to _step04.
  """
  def _SynClass(self):
    return Syn
  def _InvClass(self):
    return Invert
  def _step00(self, **kwargs):
    kwargs['env'] = dict(
      SCHEDULER_DUMPRAWGRAD='yes',
      SCHEDULER_DUMPGRAD='yes',
      SCHEDULER_DUMPPREC='yes',
      # SLAVES_DUMPGRAD=1, 
      SLAVES_DUMPADJOINT=1, 
      SLAVES_DUMPRESIDS=1)
    self.syn.init(**kwargs)
    self.Invert.init(**kwargs)
  def _step01(self, **kwargs):
    self.syn.i.create(**kwargs)
    self.syn.i.preprocess(**kwargs)
  def _step02(self, **kwargs):
    self.syn.run(**kwargs)
  def _step03(self, **kwargs):
    pass
  def _step04(self, inv_vs_syn, kw_rnf, **kwargs):
    self.Invert.i.create(self.syn)
    
    if inv_vs_syn == 'homo':
      assert 'vp_homog' in kwargs
      self.Invert.i.svp.create(ProjMod(self.Invert).create('homo', val=kwargs['vp_homog']))
    elif inv_vs_syn == 'scaled_vp':
      if 'avp' in kwargs:
        avp = kwargs['avp']
        assert isinstance(avp, float)
        self.Invert.i.svp.create(avp * self.syn.i.tvp.read())
      else:
        raise ValueError('You need to provide a multiplicative factor avp such that '+\
            'Invert.svp = avp * syn.tvp')
    elif inv_vs_syn == 'anom':
      assert 'sphere_radius' in kwargs
      assert 'sphere_amplit' in kwargs
      assert 'sphere_backgr' in kwargs
      self.Invert.i.svp.create(ProjMod(self.Invert).create('sphere', r=kwargs['sphere_radius'], \
        ampl=kwargs['sphere_amplit']), bckd=kwargs['sphere_backgr'])
    else:
      raise ValueError('Unknown value of inv_vs_syn: %s' % inv_vs_syn)
      
    self.Invert.i.sp.run(cat=0)
    self.Invert.i.rnf.create(**kw_rnf)
  def _step05(self, **kwargs):
    self.Invert.run(**kwargs)
  def _step06(self, **kwargs):
    self.Invert.reinit() # necessary!
    figure(15,5)
    plt.subplot(121)
    self.Invert.o.rawgrad.it[1].plot(**kwargs)
    plt.subplot(122)
    self.Invert.o.vp.it[1].plot(**kwargs)    
@logged
class InvertField(Workflow):
  def _SynClass(self):
    return Syn
  def _InvClass(self):
    return InvertField  
  def _get_data_fnames(self, **kwargs):
    """
    Get a list of SEGY file names to extract 
    raw data from.

    Returns
    -------
    fnames : list
        List of file names.
    """
    p = self.syn
    self.dataset_id = kwargs['dataset_id']
    if isinstance(self.dataset_id, list):
      fnames = []
      for dataset_id in self.dataset_id:
        fnames += p.exp.obs[dataset_id].fnames
    else:
      fnames = p.exp.obs[self.dataset_id].fnames
    return fnames
  def _get_pre_params(self, **kwargs):
    """
    Get parameters for pre-processing 
    of the input.

    Returns
    -------
    pre_params : dict
      Parameters for the pre-processor.
    """
    cat = kwargs.get('cat', 0)
    if 'ztype' not in kwargs:
        raise TypeError('ztype has to be provided as a keyword argument.')
    if 'addtodepth' not in kwargs:
        raise TypeError('addtodepth has to be provided in as a keyword argument.\n' + \
          'It should be equal to MINUS z1 of the model box, in metres.')      
    if 'reciprocity' not in kwargs:
        print('Warning. Reciprocity set to default - True.')
        raise TypeError('reciprocity has to be provided in as a keyword argument.')
    pre_params = dict(ztype=kwargs['ztype'], addtodepth=kwargs['addtodepth'],\
      reciprocity=kwargs['reciprocity'], cat=cat)
    return pre_params
  def _step00(self, *args, **kwargs):
    self.syn.init(**kwargs)
    self.Invert.init(**kwargs)
  def _step01(self, **kwargs):
    p = self.syn
    # these two preceed i.create() as the latter is quite slow
    # and it's better to throw some errors before waiting for it
    p.i.rse.prep(self._get_data_fnames(**kwargs))
    p.i.sp.prep(self._get_pre_params(**kwargs))
    p.i.create() # true models and wavelet (what about fs? not yet)
    p.i.sp.run()
    p.i.rnf.prep(**kwargs)
    env_var = p.env.var
    p.env.var = dict(env_var,
      SCHEDULER_DUMPRAWGRAD='yes',
      SCHEDULER_DUMPGRAD='yes',
      SLAVES_DUMPGRAD=1, 
      SLAVES_DUMPADJOINT=1, 
      SLAVES_DUMPRESIDS=1)    
    self.syn.i.pbs.no[0].prep(q='debug')
    self.syn.i.pbs.no[1].prep(**kwargs) 
  def _step02(self, **kwargs):
    self.__log.info('Please schedule the synthetic calculation on your cluster.')
  def _step03(self, **kwargs):
    p = self.syn
    figure(15,10)
    plt.subplot(211)
    p.o.syn.plot(norm='max', overwrite=0, cmap='Greys')
    plt.xlim(0,500)
    plt.subplot(212)
    p.o.syn.plot(norm='max', overwrite=0, cmap='Greys')
  def _step04(self, **kwargs):
    self.__log.info('Preparing input for FWI using synthetic first breaks...\n')
    if 'kw_filt' not in kwargs:
      raise ValueError('You have to provide kw_filt dict as a keyword argument!')    
    self.twin = kwargs['twin']
    self.fmax = kwargs['fmax']
    p = self.Invert
    cat = kwargs.get('cat', 0)
    p.i.create(self.syn, b_abs=self.b_abs, e_abs=self.e_abs, 
          filt_kwargs=kwargs['kw_filt'], mute_kwargs={'twin': self.twin},
          blocks=[{'nits': 1, 'minoff': 0, 'freq': self.fmax}],
          cat=cat)
    qc_filt(p, self.syn, overwrite=1, overwrite_mmp=1)
    p.i.obs.plot(overwrite=1, overwrite_mmp=1)
    self.__log.info('Preparing PBS scripts for the synthetic run...\n')
    for p in [self.Invert]: 
      p.i.pbs.no[0].prep(q='debug')
      p.i.pbs.no[1].prep(**kwargs)     
  def _step05(self, **kwargs):
    self.__log.info('Please schedule the inversion on your cluster.')
  def _step06(self, **kwargs):
    freq = self.fmax - 1 # Hz
    overwrite = kwargs.get('overwrite', True)
    overwrite_mmp = overwrite
    p = self.Invert
    it = 1
    for sid in [s.ID for s in p.i.s.read().li]:
      for lid in sorted(p.i.obs.read_header(overwrite=0)['ep'].unique()):
        _ = plot_out_data(p, it=it, sid=sid, lid=lid, freq=freq, 
                interleave=1, 
                overwrite=overwrite, overwrite_mmp=overwrite_mmp)
        break # shortest line throws an error        
@logged
class InvertFieldCheq(InvertField):
  """
  Chequerboard test using field-data
  geometry.

  Notes
  -----
 
  """
  def _SynClass(self):
    return SynChecq  
  def _InvClass(self):
    return InvertFieldChecq    
  def _add_checker(self, kw_chqr, **kwargs):
    """
    Create three objects for each model parameter
    (currently only Vp): TrueVp, AnomVp and BckgVp.

    Parameters
    ----------
    kw_chqr : dict
        Keyword args passed to create_chqr function.
    """
    from fullwavepy.project.files.gridded.models import ModelFileSgy
    from fwilight.seis import create_chqr, Anom
    tvp = self.syn.i.tvp.read(overwrite=1, overwrite_mmp=1)
    anom = kw_chqr['anom']
    if anom == 'spike':
      a = Anom()
      x0, y0, z0 = kw_chqr['x0'], kw_chqr['y0'], kw_chqr['z0']
      ampl = kw_chqr['ampl']
      fwhm = kw_chqr['fwhm']
      chq = a.create(dims=tvp.shape, box=tvp.extent.ravel(), fwhm=fwhm, ampl=ampl,
        x0=x0, nx=1, dx=1, y0=y0, ny=1, dy=1, z0=z0, nz=1, dz=1)
    else:
      chq = create_chqr(tvp.shape, **kw_chqr)
    p = self.syn
    self.syn.i.tvp.bckg.create(tvp)
    self.syn.i.tvp.anom.create(chq)
    self.syn.i.tvp.create(tvp + tvp * chq)
  def _step01(self, data_fname, kw_chqr, **kwargs):
    """
    Prepare input for synthetic calculation.

    Parameters
    ----------
    data_fname : SGY file
        File (Observed, OutSeis etc.) to save time.
    """
    p = self.syn
    # these two preceed i.create() as the latter is quite slow
    # and it's better to throw some errors before waiting for it
    p.i.rse.prep([data_fname])
    p.i.sp.prep(**self._get_pre_params(**kwargs))
    p.i.create() # true models and wavelet (what about fs? not yet)
    self._add_checker(kw_chqr, **kwargs)
    p.i.sp.run()
    p.i.rnf.prep(**kwargs)
    env_var = p.env.var
    p.env.var = dict(env_var,
      SCHEDULER_DUMPRAWGRAD='yes',
      SCHEDULER_DUMPGRAD='yes',
      SLAVES_DUMPGRAD=1, 
      SLAVES_DUMPADJOINT=1, 
      SLAVES_DUMPRESIDS=1)    
    self.syn.i.pbs.no[0].prep(q='debug')
    self.syn.i.pbs.no[1].prep(**kwargs) 
  def _step04(self, blocks, **kwargs):
    self.__log.info('Preparing input for FWI using synthetic first breaks...\n')
    if 'kw_filt' not in kwargs:
      raise ValueError('You have to provide kw_filt dict as a keyword argument!')    
    self.twin = kwargs['twin']
    p = self.Invert
    cat = kwargs.get('cat', 0)
    def get_bounds(b_or_e):
      if hasattr(self, b_or_e):
        b_or_e = getattr(self, b_or_e)
      else:
        b_or_e = kwargs[b_or_e]
      return b_or_e
    b_abs = get_bounds('b_abs')
    e_abs = get_bounds('e_abs')
    p.i.create(self.syn, b_abs=b_abs, e_abs=e_abs,     
          filt_kwargs=kwargs['kw_filt'], mute_kwargs={'twin': self.twin},
          blocks=blocks,
          cat=cat)
    # qc_filt(p, self.syn, overwrite=1, overwrite_mmp=1)
    p.i.obs.plot(overwrite=1, overwrite_mmp=1)
    self.__log.info('Preparing PBS scripts for the synthetic run...\n')
    for p in [self.Invert]: 
      p.i.pbs.no[0].prep(q='debug')
      p.i.pbs.no[1].prep(**kwargs) 
