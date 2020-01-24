import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ipywidgets import (interactive, interact, interact_manual, fixed,
                        IntSlider)
import pandas as pd

from fullwavepy.logging_config import *
from fullwavepy.generic.system import get_files
from fullwavepy.generic.parse import strip, exten
from fullwavepy.generic.array import Arr, tseries2array
from fullwavepy.plot.misc import time_freq
from fullwavepy.plot.generic import plot, compare
from fullwavepy.plot.oned import colors
from fullwavepy.plot.twod import plot_image
from fullwavepy.ioapi.generic import save_txt, read_txt, read_any
from fullwavepy.ioapi.fw3d import read_vtr, save_vtr
from fullwavepy.ioapi.json import save_json
from fullwavepy.ioapi.segy import header2json
from fullwavepy.project.types.basic import Proj, ProjSyn, ProjInv
from fullwavepy.project.types.deriv import ProjSynVsObs, ProjInvSyn
from fullwavepy.project.types.extra import ProjFsQC
from fullwavepy.signal.su import su_filter