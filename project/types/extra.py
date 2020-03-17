"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for peremoveission writing to k.chrapkiewicz17@imperial.ac.uk.

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
from fullwavepy.project.types.basic import ProjSyn


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

