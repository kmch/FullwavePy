"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced

from fullwavepy.generic.parse import kw
from fullwavepy.project.files.text.submit import PbsFile


@traced
@logged
class PbsFileArcher(PbsFile):
  """
  Notes
  -----
  Hopefully, some methods of PbsFile
  can be reused here, at least _write_env
  seems promising.

  Others should be overwritten here,
  e.g. _write_head, for sure.

  Make sure the order they are written to the file
  specified in _write_all will work for Archer.
  If not, overwrite this method here as well.
  In general, it just writes directly to the file:
  with open(self.fname, 'w') as f:
    f.write(first_line_of_file)
    ...
    f.write(last_line_of_file)

  For reference, see:
  - fullwavepy.project.files.text.submit.PbsFile
  - same.as.above.BashFile - simpler case.

  The name of the file will stay as proj_name-Run.pbs.
  I don't think it matters for queueing systems.
  """
  def _write_all(self, **kwargs):
    with open(self.fname, 'w') as f:
      self._set_resources(**kwargs)
      f = self._write_head(f, **kwargs)
      f = self._write_body(f, **kwargs)
      f = self._write_srun(f, **kwargs)
      f = self._write_foot(f, **kwargs)
    
  def _find_optimal_resources(self, optimize='idle', **kwargs):
    from fullwavepy.numeric.generic import divisors_of, decimal
    
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
      self.__log.warning('Choosing the largest mpiprocs')
      self.mpiprocs = best[0]
      self.ompthreads = best[1]
    else:
      self.mpiprocs = mpis[i_max]
      self.ompthreads = omps[i_max]
  def _set_resources(self, **kwargs):
    """
    Set computational resources to be requested for.
    
    """
    q = kw('q', '', kwargs)
    self.q = q # USED BY _write_mpiexec
    self.__log.info('Queue selected: ' + q)
    
    self._set_time(**kwargs)
      
    if q == 'standard':
      select_min = 1
      select_max = 256 # TO CHECK ACTUAL VALUE
      node_type = kw('node_type', 'amd', kwargs)
      ncpus_max =  128
      ncpus_min = ncpus_max
      mem_min = 1
      mem_max = 256

    elif q == 'short':
      select_min = 1
      select_max = 8 # TO CHECK ACTUAL VALUE
      node_type = kw('node_type', 'amd', kwargs)
      ncpus_max =  128
      ncpus_min = ncpus_max
      mem_min = 1
      mem_max = 256

    elif q == 'long':
      select_min = 1
      select_max = 64 # TO CHECK ACTUAL VALUE
      node_type = kw('node_type', 'amd', kwargs)
      ncpus_max =  128
      ncpus_min = ncpus_max
      mem_min = 1
      mem_max = 256
   
    else:
      raise ValueError('Unknown queue type: ' + q)
    
    self.select = kw('select', select_min, kwargs)
    self.ncpus = kw('ncpus', ncpus_max, kwargs)
    self.mem = kw('mem', mem_max, kwargs)
    self.account = kw('account', 'n03-mp', kwargs)
    self.totalprocs = kw('totalprocs', '', kwargs)

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
  def _write_head(self, f, **kwargs):
    """
    Create a header (resources, modules etc.) 
    Slurm script file for submitting a job to 
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

    f.write('#SBATCH -J ' + job_name + '\n')
    f.write('#SBATCH -o ' + self.jobout + '\n')
    f.write('#SBATCH -e ' + self.joberr + '\n')
    f.write('#SBATCH -N ' + str(self.select) + '\n')
    f.write('#SBATCH -n ' + str(self.totalprocs) + '\n')
    f.write('#SBATCH --ntasks-per-node ' + str(self.mpiprocs) + '\n')
    f.write('#SBATCH -t ' + self.time + '\n')
    f.write('#SBATCH -c ' + str(self.ompthreads) + '\n')
    f.write('#SBATCH -A ' + self.account + '\n')
#    f.write('#SBATCH --mem=' + str(self.mem) + '\n') Archer2 does not want memorey specification
    f.write('#SBATCH --partition=standard\n')
    f.write('#SBATCH --qos='+ self.q +'\n')
    
    f.write('\nstart=`date +%s`\n')
    
    # NOTE {{ ESCAPES {
    paths = """
# PATHS
code_path={code_path}
echo 'code_path: '${{code_path}}

script_dir=`pwd`
work_dir=$script_dir/../out
    
# CAVEAT  
rm $work_dir/{proj_name}-Runfile.key   
rm $work_dir/{proj_name}-Ghost.*

# HARD-LINK INPUT FILES TO OUTPUT DIR WHERE THE CODE WILL RUN
ln $script_dir/* $work_dir
 
# CHANGE DIRECTORY TO PROJECT OUTPUT
cd $work_dir
   
    """.format(code_path=code_path, proj_name=self.proj.name)
    f.write(paths)

    f.write("\n# ARGUMENTS OF fullwave3D.exe\n")
    f.write('CPNUM=' + str(cpnum) + '\n')
    f.write('NTHREAD=' + str(self.ompthreads) + ' # EQUAL TO ompthreads\n')

    f.write("\n# FULLWAVE'S ENVIRONMENT VARIABLES\n")
    f = self._write_env(f, **kwargs)
    
    f.write('# LOAD MODULES\n')
    f.write('module load epcc-job-env \n')
    
    return f 
  def _write_srun(self, f, **kwargs):
    """

    Notes
    -----
    Double redirection >> guarantees _prev files are not overwritten
    if they already exist, and created if otherwise.     

    """
    f.write('\n# RUN FULLWAVE \n')

    jobmng = 'srun '
    jobopt = '--hint=nomultithread --distribution=block:block '
    
    f.write(jobmng + jobopt +' ${code_path} ' + self.proj.name + 
            ' ${CPNUM} ${NTHREAD}' + 
            ' 1>> $work_dir/' + self.proj.out.out.no[self.run_id].name +
            ' 2>> $work_dir/' + self.proj.out.err.no[self.run_id].name + 
            '\n')
    
    f.write('\nstat=$?\n')
    f.write('echo "Exit status: "$stat\n')
    
    return f  