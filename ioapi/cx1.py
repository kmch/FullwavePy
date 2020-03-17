"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.project.files.submit import PbsFile


#FIXME REPLACE WITH CONTIGUOUS STRING

# -------------------------------------------------------------------------------


@traced
@logged
class PbsFileCx1(PbsFile):
  """
  
  Notes
  -----
    NOTE:(source: Imperial MPI.submit) 
    as soon as the walltime is reached the content
    of $TMPDIR is DELETED.
    
    pbsexec before mpiexec gives 20 min grace period,
    i.e.it terminates th pbsdsh2 code run with mpiexec 20 min 
    before end of the requested wall time. This allows 
    the rest of the .pbs script to be executed. 
    
    pbsdsh2 is necessary for throughputting on ALL nodes 
    of a multi-node job.
    
    Other Fullwave people don't use it at all (cd to dir). 
    But copying was recommended RCS (HPC) at Imperial.
    
    If * follows a prefix, it must be proceeded by \, e.g.:
    $default_work_dir/prefix\*.vtr
  
    The \ prior to * is necessary to escape the asterisk and ensure 
    it is not expanded by the shell when passed through pbsdsh2.
    
  """
  
  #prev_suffix = '_prev'
  
  # -----------------------------------------------------------------------------  
  
  def _set_resources(self, **kwargs):
    """
    Set computational resources to be requested for.
    
    """
    q = kw('q', 'pqmrwarn', kwargs)
    self.q = q # USED BY _write_mpiexec
    self.__log.info('Queue selected: ' + q)
    
    self._set_time(**kwargs)
      
    if q == 'pqmrwarn':
      select_min = 1
      select_max = 100 # TO CHECK ACTUAL VALUE

      node_type = kw('node_type', 'both', kwargs)
      if node_type == 'old' or node_type == 'both':
        ncpus_max =  40
      elif node_type == 'new':
        ncpus_max =  48
      ncpus_min = ncpus_max
      mem_min = 1
      mem_max = 128
      
    elif q == 'throughput':
      select_min = 1
      select_max = 1
      ncpus_min = 1
      ncpus_max = 8
      mem_min = 1
      mem_max = 96
    
    elif q == 'general':
      select_min = 1
      select_max = 16
      ncpus_min = 32
      ncpus_max = 32    
      mem_min = 1
      mem_max = kw('mem_max', 62, kwargs) # OR 124
    
    elif q == 'single':
      select_min = 1
      select_max = 1
      ncpus_min = 48
      ncpus_max = 48    
      mem_min = 1
      mem_max = 124 
   
    elif q == 'multi':
      raise ValueError('1 Feb 2020 -- multinode class retired -- minimum node count for large class is now reduced')
      select_min = 3
      select_max = 16
      ncpus_min = 12
      ncpus_max = 12    
      mem_min = 1
      mem_max = 46
      self.__log.warn('mem_max={}. Likely to exceed it and kill the job'.format(mem_max))
    
    elif q == 'debug':
      select_min = 1
      select_max = 1
      ncpus_min = 1
      ncpus_max = 8
      mem_min = 1
      mem_max = 96
    
    else:
      raise ValueError('Unknown queue type: ' + q)
    
    self.select = kw('select', select_min, kwargs)
    self.ncpus = kw('ncpus', ncpus_max, kwargs)
    self.mem = kw('mem', mem_max, kwargs)
    
    self._find_optimal_resources(**kwargs)
    
    if self.select < select_min or self.select > select_max:
      raise ValueError('select=' +  str(self.select) + ' incompatible with queue ' + q +
                       ': %s <= select <= %s' % (str(select_min), str(select_max)))   

    if self.ncpus < ncpus_min or self.ncpus > ncpus_max:
      raise ValueError('ncpus=' +  str(self.ncpus) + ' incompatible with queue ' + q +
                       ': %s <= ncpus <= %s' % (str(ncpus_min), str(ncpus_max)))

    if self.mem < mem_min or self.mem > mem_max:
      raise ValueError('mem=' +  str(self.mem) + ' incompatible with queue ' + q +
                       ': %s <= mem <= %s' % (str(mem_min), str(mem_max)))
    
    self._create_verbosity_triggers(**kwargs)
  
  # -----------------------------------------------------------------------------

  def _set_time(self, **kwargs): #FIXME
    """
    Calculate time request for a given queue type.
    
    Parameters
    ----------
    q : str 
      Queue.
    
    **kwargs : keyword arguments (optional)
      Current capabilities:
    
    FIXME: 
    if below 30 min. it is automatically a debug queue
    unless a pqmrwarn??
    
    """
    # MUST BE HH:00:00 OTHERWISE RAISES ERROR, AT LEAST FOR THE general QUEUE
    hours = kw('hours', 1, kwargs) 
    #err = False
    #try:
      #time = kwargs['time']
      #err = True
    #except:
      #pass
    
    #if err:
      #raise ValueError('Parameter <time> is obsolete. Use <hours> instead.')
    
    if isinstance(hours, str):
      time = hours
    else:
      time = str(hours)
      if hours < 10:
        time = '0' + time
      
      time += ':00:00'
    
    
    if self.q == 'debug':
      minutes = kw('minutes', 1, kwargs)
      assert isinstance(minutes, int) 
      assert minutes <= 30
      assert minutes >= 1
      if minutes < 10:
        minutes = '0' + str(minutes)
      else:
        minutes = str(minutes)
      time = '00:' + minutes + ':00'
  
    self.time = time

  # -----------------------------------------------------------------------------
  
  def _find_optimal_resources(self, optimize='idle', **kwargs):
    from fullwavepy.generic.math import divisors_of, decimal
    
    nshots = int(self.proj.inp.runfile.read()['ncomp'])
    self.__log.info('No. of shots in the runfile (ncomp): ' + str(nshots))
    
    #print('No. of shots = no. of slave processes to run: ', nshots)
    #mpi_to_run_total = nshots + 1
    #print('Total no. of tasks (slaves + scheduler) to run: ', mpi_to_run_total)
    #print('No. of nodes = ', self.select)
    #print()
    
    mpiprocs = kw('mpiprocs', None, kwargs)
    if 'mpiprocs' in kwargs:
      mpi = kwargs['mpiprocs']
      if self.ncpus % mpi != 0:
        raise ValueError('mpiprocs needs to be a divisor of ncpus: ' + str(self.ncpus))
      possible_mpis = [mpi]
    else:
      self.__log.info('mpiprocs not specified. It will be chosen to maximize a decimal ' + 
                      'part of nshots / (mpiprocs-1) for a given select=')
      possible_mpis = divisors_of(self.ncpus)
    
    mpis, omps, decims = [], [], []
    i = 0
    i_max = 0
    decim_max = 0
    ideals = []
    for mpi in possible_mpis:
      #print('No. of requested MPI processes per node: ', mpi)
      omp = self.ncpus // mpi
      #print('No. of requested OMP threads per node: ', omp)
      
      mpi_total = self.select * mpi
      #print('No. of all processes requested: ', mpi_total)
      if mpi_total == 1:
        continue
      a = nshots / (mpi_total - 1)
      decim = decimal(a)
      
      if decim == 0:
        ideals.append([mpi, omp])
      
      if decim > decim_max:
        decim_max = decim
        i_max = i
      
      mpis.append(mpi)
      omps.append(omp)
      decims.append(decim)
      
      #print('nshots % (mpi_total - 1): ', remainder)
      #print()
      i += 1
    
    self.__log.info('Max. decimal place ' + str(decim_max) +
                    ' is for mpiprocs: ' + str(mpis[i_max]) +
                    ' and ompthreads: ' + str(omps[i_max]))
    
    if len(ideals) > 0:
      self.__log.info('There are mpiprocs value(s) that give integer: ' +
                      'nshots / (mpiprocs-1)')
      
      mpi_max = 0
      best = ideals[0]
      for ideal in ideals:
        mpi, omp = ideal
        self.__log.info('mpiprocs: ' + str(mpi) + 
                        ', ompthreads: ' + str(omp))
        if mpi > mpi_max:
          best = [mpi, omp]
      self.__log.warn('Choosing the largest mpiprocs')
      self.mpiprocs = best[0]
      self.ompthreads = best[1]
    else:
      self.mpiprocs = mpis[i_max]
      self.ompthreads = omps[i_max]
    
  # -----------------------------------------------------------------------------  
  
  def _write_head(self, f, **kwargs):
    """
    Create a header (resources, modules etc.) 
    PBS script file for submitting a job to 
    a queueing system. 
    
    Parameters
    ----------
    f : file
      Already open file.
    proj_name : str
      Project name assumed to be provided 
      as a first argument of the FWI code
      being a prefix of all FWI input files
    **kwargs : keyword arguments (optional)
      
    Returns
    -------
    f : file 
      File with filled header.
    
    Notes
    -----
    Nomenclature:
    process = MPI rank = rank = task = job (SET AS mpiprocs AND FULLWAVE'S -K WHEN RUN LOCALLY)
    thread = ompthread  = virtual core (SET BOTH AS ompthreads FULLWAVE'S numthreads) 
    cpu = logic_core = logical core (SET AS ncpus)
    
    Unlike local runs,
    "it's not necessary to add "-n" or any other flag to specify the number of ranks."  
    (http://www.imperial.ac.uk/admin-services/ict/self-service/research-support/
    rcs/computing/high-throughput-computing/configuring-mpi.submit/)
    
    "default values of  mpiprocs==ncpus and ompthreads==1" 
    (-||-)
    
    """
    code_path = self.proj.exe['fullwave']
    #FIXME: cpnum DOESN'T MATTER (ONLY LastCheckpoint.txt), DOESN'T IT?
    cpnum = kw('cpnum', -1, kwargs) # CHECKPOINT NO. TO RESTART FROM (ARG. OF fullwave3D.exe)
    job_name = kw('job_name', self.proj.name, kwargs)
     
     
     
     
    self.jobout = '../out/' + self.proj.out.jobout.no[self.run_id].name
    self.joberr = '../out/' + self.proj.out.joberr.no[self.run_id].name
    
    # SHEBANG
    f.write('#!/bin/bash\n\n')
  
    f.write('##\n')
    f.write('# FULLWAVE3D WILL BE RUN FROM proj/out/\n')
    f.write('# For more explanation, see help(fullwavepy.project.files.runfiles.PbsFile).\n')
    f.write('#\n')
    f.write('##\n\n')

    f.write('#PBS -N ' + job_name + '\n')
    f.write('#PBS -o ' + self.jobout + '\n')
    f.write('#PBS -e ' + self.joberr + '\n')
    f.write('#PBS -l walltime=' + self.time + '\n')
    f.write('#PBS -l select=' + str(self.select) +
                 ':mpiprocs=' + str(self.mpiprocs) + 
               ':ompthreads=' + str(self.ompthreads) + 
                    ':ncpus=' + str(self.ncpus) + 
                      ':mem=' + str(self.mem) + 'gb\n')    
    f.write('#PBS -l place=scatter:excl\n') 

    f.write('\nstart=`date +%s`\n')
    
    # NOTE {{ ESCAPES {
    paths = """
    # PATHS
    code_path={code_path}
    echo 'code_path: '${{code_path}}

    work_dir=$PBS_O_WORKDIR/../out/ # $PBS_O_WORKDIR IS THE ONE CONTAINING THIS SCRIPT
    
    # HARD-LINK INPUT FILES TO OUTPUT DIR WHERE THE CODE WILL RUN
    ln $PBS_O_WORKDIR/* $work_dir
    
    # CHANGE DIRECTORY TO PROJECT OUTPUT
    cd $work_dir
    
    """.format(code_path=code_path, proj_name=self.proj.name)
    f.write(paths)
    #f.write('code_path=' + code_path + '\n')
    #f.write("echo 'code_path: '${code_path}\n")
    #f.write("echo 'PBS_O_WORKDIR: '$PBS_O_WORKDIR\n")
    #f.write('dir_with_pbs_script=$PBS_O_WORKDIR\n')
    #f.write('default_work_dir=$TMPDIR\n')
    #f.write('proj_inp_dir=$dir_with_pbs_script/../inp/\n') 
    #f.write('proj_out_dir=$dir_with_pbs_script/../out/\n') 
    #f.write('work_dir=$proj_out_dir # RUN FULLWAVE IN proj/out/!\n')
    #f.write('cd $work_dir \n')
    #f.write('ln -s $proj_inp_dir/* . \n')
    
    
    f.write("\n# DISABLE PINNING OF THE PROCESSES (MAKE ALL NODE CORES AVAILABLE TO ALL PROCESSES)\n")
    f.write('unset NCPUS\n')
    f.write('export I_MPI_PIN=no # (DOES NOT WORK IF RUNNING ON A WHOLE NODE) \n')

    f.write("\n# ARGUMENTS OF fullwave3D.exe\n")
    f.write('CPNUM=' + str(cpnum) + '\n')
    f.write('NTHREAD=' + str(self.ompthreads) + ' # EQUAL TO ompthreads\n')

    f.write("\n# FULLWAVE'S ENVIRONMENT VARIABLES\n")
    f = self._write_env(f, **kwargs)
    
    f.write('# LOAD MODULES\n')
    f.write('module unload mpi/intel-2019 \n')
    f.write('module load mpi/intel-2018\n')
    f.write('module load intel-suite\n')
    
    return f 
  
  # -----------------------------------------------------------------------------
  
  def _write_mpiexec(self, f, **kwargs):
    """

    Notes
    -----
    Double redirection >> guarantees _prev files are not overwritten
    if they already exist, and created if otherwise.     

    """
    f.write('\n# RUN FULLWAVE \n')
    if True:
      self.__log.warn('Apparently after 1Feb2020 pbsexec (grace period) does not work (regardless of the queue)')
      pbsexec = ''
    else: # NOW IT'S NOT RECOGNIZED BY ANY QUEUE
      pbsexec = 'pbsexec '
    f.write(pbsexec + 'mpiexec ${code_path} ' + self.proj.name + 
            ' ${CPNUM} ${NTHREAD}' + 
            ' 1>> $work_dir/' + self.proj.out.out.no[self.run_id].name +
            ' 2>> $work_dir/' + self.proj.out.err.no[self.run_id].name + 
            '\n')
    
    f.write('\nstat=$?\n')
    f.write('echo "Exit status: "$stat\n')
    
    return f  

  # -----------------------------------------------------------------------------  
  

# ------------------------------------------------------------------------------- 

