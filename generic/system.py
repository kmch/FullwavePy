"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer


# -------------------------------------------------------------------------------


@traced
@logged
def duplicate(source, destination, cmd='cp', **kwargs):
  """
  Duplicate (cp/mv/ln) a file.
  
  Notes
  -----
  Hard (ln) vs. soft/symbolic (ln -s) links:
  "If we create a hard link to the file and then delete the file, 
  we can still access the file using hard link. 
  But if we create a soft link of the file and then delete the file, 
  we can’t access the file through soft link and soft link becomes dangling."
  (https://www.geeksforgeeks.org/soft-hard-links-unixlinux/)
  
  """
  duplicate._log.debug('cmd: ' + cmd)
  cmd = cmd + ' ' + source + ' ' + destination
  o, e  = bash(cmd, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
def bash(command, path='./', truncate=150, **kwargs):
  """
  Call a bash command/script.
  
  Parameters
  ----------
  command : str 
    Command as a 1 string.
  path : str 
    Where is it to run command.
  truncate : int 
    Max. no. of chars of the stdin/err to show.
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  
  Returns
  -------
  out : str 
    Standard output.
  err : str
    Standard error.
  
  Notes
  -----
  It respects wildcards.
  
  Examples
  --------
  >>> from fullwavepy.generic.system import bash
  >>> bash('ls', path='../')
  <content of the parent dir>
  
  """
  from subprocess import Popen, PIPE
  
  command = 'cd ' + path + ' && ' + command
  bash._log.debug('Command to run: ' + command)
  
  proc = Popen(command, shell=True,
               stdout=PIPE, stderr=PIPE)
  
  out, err = proc.communicate()
  
  # DECODE, OTHERWISE class 'bytes' OBJECT WITH \n PRINTED INSTEAD 
  # ACTUALLY BREAKING LINES
  
  try:
    out = out.decode('UTF-8') 
    err = err.decode('UTF-8') 
  except UnicodeDecodeError as error_message:
    bash._log.warning('Ignoring UnicodeDecodeError: ' + str(error_message))
  
  
  if len(out) > truncate:
    bash._log.debug('Truncated stdout (after truncate= ' + str(truncate) + 
                    ' chars):' + str(out[ :truncate]) + '...')
  elif len(out) > 0:
    bash._log.debug('stdout:' + out)
  
  if len(err) > 0:
    if len(err) > truncate:
      err = str('Truncated stderr (first ' + str(truncate) + 
                ' chars):' + err[ :truncate] + '...')
    bash._log.warn('stderr:' + err)
      
  
  return out, err


# -------------------------------------------------------------------------------


@traced
@logged
def exists(fname, **kwargs):
  """
  Tiny wrapper around os.path.exists
  (checks if the file exists).
  
  Parameters
  ----------
  **kwargs : keyword arguments (optional)
    Current capabilities:
  
  Returns
  -------
  exists : bool 
    True if the file exists.
  
  Notes
  -----
  Print a debug-message not a warning, because 
  lack of file is not necessarily an error.
  
  """
  import os

  found = False 
  if os.path.exists(fname):
    found = True  
  else:
    found = False
    exists._log.debug('File: ' + fname + ' does not exist.')
  
  return found


# -------------------------------------------------------------------------------


@timer
@traced
@logged
def get_files(path, pattern, **kwargs):
  """
  Get list of files located in path 
  defined by the REGEX pattern.
  
  Parameters
  ----------
  path : str 
    Full path to the directory supposed to
    contain the files.
  pattern : str 
    REGEX pattern defining names of the files;
    it can contain wildcards, e.g.
    pattern = 'project-*vtr'
  
  Returns
  -------
  files_list : list 
    Sorted alphabetically list of files names each of which includes
    the full path (!).  
  
  Notes
  -----
  Comes in handy in Jupyter notebooks.
  
  It gets reaaally slow if a no. of files in the path 
  is big (>10k)
  
  """
  from os import listdir as ls
  from fnmatch import filter as expand
  
  files_list = sorted([path + '/' + f for f in expand(ls(path), pattern)])
  
  if len(files_list) == 0:
    get_files._log.warn('No files matching ' + 
                         pattern + ' in ' + path)  
  
  
  return files_list


# -------------------------------------------------------------------------------


