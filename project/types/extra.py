"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

Notes
-----
Python devs usually don't write nested classes, neither do I.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw
from fullwavepy.generic.system import bash, exists
from fullwavepy.project.types.basic import *

# extra -> case_study?


@traced
@logged
class ProjExperiment(Proj):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class ProjSynSingleStation(ProjExperiment, ProjSyn):
  def __init__(self, name, experiment, dataset_id, station_id, *args, **kwargs):
    self.ex = experiment
    self.dataset = self.ex.dataset[dataset_id]
    self.sid = station_id
    self.shotgather = self.dataset[self.sid]
    self.sgh = self.shotgather
    # super().__init__(name, *args, **kwargs)


# -------------------------------------------------------------------------------


@traced
@logged
class ProjSynSingleLine(ProjExperiment, ProjSyn):
  pass


# -------------------------------------------------------------------------------


@traced
@logged
class ProjFsQC(ProjSyn):
  def _init_input(self, **kwargs):
    from fullwavepy.project.files.gridded.surfaces import FsFile
    super()._init_input(**kwargs)
    self.inp.fs = FsFile(self, self.inp.path, **kwargs)


# -------------------------------------------------------------------------------


#@traced
#@logged
#class ProjSourceQC(ProjSyn):
#  pass


# -------------------------------------------------------------------------------


#@traced
#@logged
#class ProjDiffractor(ProjInvSyn):


# -------------------------------------------------------------------------------

