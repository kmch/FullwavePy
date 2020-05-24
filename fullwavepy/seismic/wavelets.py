"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
import matplotlib.pyplot as plt
from autologging import logged, traced

from fullwavepy.generic.decor import timer
from fullwavepy.generic.parse import kw, del_kw

from fullwavepy.ndat.arrays import *

from fullwavepy.seismic.data import *


@traced
@logged
class Wavelet(Data):
  def plot(self, **kwargs):
    assert self.shape[0] == 1
    assert self.shape[1] == 1
    a = Arr1d(self[0,0]) # FIXME
    a.plot(**kwargs)


# -------------------------------------------------------------------------------

