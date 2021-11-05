"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
from autologging import logged
import numpy as np

from arrau.a1d import Arr1d
from fullwavepy.seismic.data import *

class SourceWavelet:
  def __init__(self, arr):
    self.arr = arr
  def plot(self, *args, **kwargs):
    self.arr.plot(*args, **kwargs)

# @traced
@logged
class Wavelet(Data):
  def plot(self, **kwargs):
    assert self.shape[0] == 1
    assert self.shape[1] == 1
    a = Arr1d(self[0,0]) # FIXME
    a.plot(**kwargs)

