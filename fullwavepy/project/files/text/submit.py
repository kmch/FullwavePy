"""
Scripts submitting FWI jobs to local shell or
queuing systems on clusters.

(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.parse import kw, del_kw, exten, strip
from fullwavepy.generic.system import bash, exists
from fullwavepy.generic.decor import timer
from fullwavepy.project.files.generic import TextProjFile
from fullwavepy.project.files.text.misc import JobFile


@traced
@logged
class PbsFile(JobFile, TextProjFile):
  """
  Script to submit a job to a queueing system
  such as PBS.
  
  """  
  def __init__(self, proj, path, run_id, **kwargs):
    """
    Parameters
    ----------
    run_id : int 
      Number of runs so far + 1
      (track the performance of the code as 
      - a function of resources requested 
      - cluster setup variability (if any)
      - ...)
    """
    suffix = 'Run'
    exten = 'pbs'
    super().__init__(proj, path, suffix, exten, run_id, **kwargs)
  def create(self, **kwargs):
    self._prevent_overwriting(**kwargs)
    self._write_all(**kwargs)
  def dupl(self, source, **kwargs):
    raise NotImplementedError('Duplicating ' + self.fname + 
                              ' is not supported. Run prepare without dupl=')  
  # -----------------------------------------------------------------------------  
  def _write_all(self, **kwargs):
    with open(self.fname, 'w') as f:
      
      f = self._write_head(f, **kwargs)
      f = self._write_body(f, **kwargs)
      f = self._write_mpiexec(f, **kwargs)
      f = self._write_foot(f, **kwargs)
  def _write_body(self, f, **kwargs):
    """
    """
    if self.proj.problem == 'synthetic':
      f.write('\n#\n') 
      f.write('# SYNTHETIC RUN \n') 
      f.write('#\n') 
      
      f.write('\n# DELETE OUTPUT OF PREVIOUS SYNTH. RUNS, OTHERWISE FULLWAVE WILL TERMINATE\n')
      #f.write('rm $work_dir/' + self.proj.name + '-*.log\n') # WHY WAS IT ON?!
      f.write('rm $work_dir/' + self.proj.name + '-*fw*vtr\n')
      f.write('rm $work_dir/' + self.proj.name + '-Observed-Time.tt?\n')
      f.write('rm $work_dir/' + self.proj.name + '-Synthetic.*\n')
      f.write('rm $work_dir/' + self.proj.name + '-Observed.*\n')  

    
    elif self.proj.problem == 'tomography':
      f.write('\n#\n') 
      f.write('# INVERSION-RUN FLOW\n') 
      f.write('#\n') 
      #f.write('\nunset SLAVES_WAVEFIELDSVTR # OTHERWISE TOO MANY SNAPSHOTS \n')

    else:
      raise ValueError('self.proj.problem: {}'.format(self.proj.problem))

    return f
  def _write_foot(self, f, **kwargs):
    f.write('for f in fw*iter*vtr bw*iter*vtr\ndo\n  mv $f ' + self.proj.name + '-$f\ndone\n')
    f.write('\nend=`date +%s` \n')
    f.write('runtime=$((end-start))\n')
    f.write("\necho 'Runtime of the PBS script: '$runtime' s' \n")
    return f      
  def _write_head(self, f, **kwargs):
    raise NotImplementedError
  def _write_env(self, f, **kwargs):
    """
    Write environment variables set-up.
    
    """
    max_ram = kw('max_ram', 40, kwargs)
    
    for key, val in sorted(self.proj.env.var.items()):
      if val is not None:
        f.write('export ' + str(key) + '=' + str(val) + '\n')
    
    f.write('\n')   
    return f
  def _write_mpiexec(self, f, **kwargs):
    raise NotImplementedError
  def _prevent_overwriting(self, **kwargs):
    """
    """
    proj = self.proj
    
    for f in [proj.inp.jobinfo, 
              proj.out.out, proj.out.err, 
              proj.out.jobout, proj.out.joberr]:
      
      fname = f.no[self.run_id].fname
      if exists(fname):
        raise OSError("{} already exists!".format(fname))
  # -----------------------------------------------------------------------------
@traced
@logged
class BashFile(JobFile, TextProjFile):
  """
  Script launching local.submit 
  in a Bash shell (not on a cluster).
  
  """
  def __init__(self, proj, path, run_id, **kwargs):
    """
    
    """
    suffix = 'RunLocally'
    exten = 'sh'
    super().__init__(proj, path, suffix, exten, run_id, **kwargs) 
  def _prevent_overwriting(self, **kwargs):
    """
    """
    proj = self.proj
    
    for f in [proj.out.out, proj.out.err]:
      fname = f.no[self.run_id].fname
      if exists(fname):
        raise OSError("{} already exists!".format(fname))  
  def create(self, **kwargs):
    ompthreads = kw('ompthreads', 4, kwargs)
    path_fullwave = self.proj.exe['fullwave_local']
    path_make = self.proj.exe['fullwave_local'][ :-len('bin/fullwave3d.exe')]
    
    #outlog = '../../' + self.proj.out.out.no[self.run_id].fname
    #errlog = '../../' + self.proj.out.err.no[self.run_id].fname
    
    # CREATE EMPY FILES FOR VERBOSE SLAVE OUTPUT
    for i in range(1, ompthreads): # THIS EXCLUDES SCHEDULER
      with open(self.proj.inp.path + 'fullwave3d-verbose-slave-' + str(i), 'w'):
        pass
    
    with open(self.fname, 'w') as f:
      f.write('#!/bin/bash\n\n')
      
      #f.write('ompthreads=' + str(ompthreads) + '\n\n')
      
      #f.write('project=' + self.proj.name + '\n')
      #f.write('path_fullwave=' + path_fullwave + '\n')
      #f.write('work_dir=../out/\n')
      #f.write('ln -s * $work_dir\n')
      
      text = """
      ompthreads={ompthreads}
      
      project={pname}
      path_fullwave={path_fullwave}
      
      make -C {path_make}
      
      echo 'current dir: '
      pwd
      
      work_dir={outdir}
      # CAVEAT it is assumed $work_dir = proj/out/  
      rm {outdir}/*

      ln {inpdir}/* $work_dir 
      cd $work_dir # !!!
      """.format(ompthreads=ompthreads, pname=self.proj.name, path_fullwave=path_fullwave, path_make=path_make,
                 inpdir=self.proj.inp.path, outdir=self.proj.out.path, 
                 err=self.proj.o.e.no[self.run_id].name)
      
      f.write(text)
      
      #f.write('path_make=' + path_make + '\n')
      #f.write('outlog=' + outlog + '\n')
      #f.write('errlog=' + errlog + '\n')
      
      
      for key, val in sorted(self.proj.env.var.items()): 
        if val is not None:
          f.write('export ' + str(key) + '=' + str(val) + '\n')
      
      
      f.write('rm ' + self.proj.name + '-*fw*vtr\n')
      if self.proj.problem == 'synthetic': 
        f.write('rm ' + self.proj.name + '-Observed-Time.tt?\n')
        f.write('rm ' + self.proj.name + '-Synthetic.*\n')
        f.write('rm ' + self.proj.name + '-Observed.*\n')  
      
      #f.write('make -C ${path_make} -j fullwave3d\n')
      cmd = 'mpiexec -n {ompthreads} $path_fullwave {pname} 1> {out} 2> {err}'.format(
        ompthreads=ompthreads, pname=self.proj.name,
        out=self.proj.o.o.no[self.run_id].name, err=self.proj.o.e.no[self.run_id].name)
      
      f.write('\n\n' + cmd + '\n\n')
   
      f.write('for f in fw*iter*vtr bw*iter*vtr\ndo\n  mv $f ${project}-$f\ndone\n\n')
      
      #f.write('mv ' + self.proj.name + '-*Synthetic* ../out/\n')   
      #f.write('mv ' + self.proj.name + '-*Observed-Time* ../out/\n')   
      #f.write('mv ' + self.proj.name + '-*DUMP* ../out/\n')   
      #f.write('mv ' + self.proj.name + '-*iter* ../out/\n')   
      #f.write('mv ' + self.proj.name + '-*CP* ../out/\n')   
  @timer
  def run(self, **kwargs):
    """
    Run fullwave locally.
    Not advised for bigger projects.
    
    """
    self._prevent_overwriting(**kwargs)
    
    cat = kw('cat', True, kwargs)
    truncate = kw('truncate', 150, kwargs)
    
    self.__log.info('Running Fullwave3D...')
    
    o, e = bash(self.fname, truncate=truncate)
    
    if len(o) != 0 and cat:
      self.__log.info(o)
    if len(e) != 0:
      self.__log.warning(e)

