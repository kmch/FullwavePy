import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 50)

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D

import plotly.express as px
import plotly.graph_objects as go
from ipywidgets import (interactive, interact, interact_manual, fixed,
                        IntSlider, BoundedIntText, Dropdown, 
                        SelectMultiple, Checkbox,
                        Layout, TwoByTwoLayout)


from fullwavepy.logging_config import *

from fullwavepy.generic.system import get_files
from fullwavepy.generic.parse import strip, exten
from fullwavepy.generic.array import Arr, Arr3d, Arr2d, Arr1d, tseries2array, WigglyData, SurfFunc, Surf, SurfParam

from fullwavepy.ioapi.generic import save_txt, read_txt, read_any
from fullwavepy.ioapi.fw3d import TtrFile, VtrFile, read_vtr, save_vtr
from fullwavepy.ioapi.segy import SgyFile

from fullwavepy.plot.generic import plot, compare, new_figure
from fullwavepy.plot.oned import colors
from fullwavepy.plot.twod import plot_image
from fullwavepy.plot.misc import time_freq

from fullwavepy.project.types.basic import Proj, ProjSyn, ProjInv
from fullwavepy.project.types.deriv import ProjInvSyn
from fullwavepy.project.types.extra import ProjFsQC

from fullwavepy.signal.su import su_filter