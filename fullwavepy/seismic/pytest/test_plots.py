"""
Unit tests of plotting functions 
using pytest's mpl plugin.

Notes
-----
First run as:
>>> pytest test_plots.py --mpl-generate-path=baseline
to save good figures in ./baseline directory

Then, after codebase modifications, run as:
>>> pytest test_plots.py --mpl
It will take AT LEAST a few sec to complete (can be slow).

NOTE It needs an X server running. 
(it works with Xming as of 4 Oct 2021)

`unittest`-style tests (Python standard-library) are understood
by pytest - it can run virtually all unittest-tests out of the box.

"""
import matplotlib.pyplot as plt
import pytest
from unittest import TestCase
from fullwavepy.seismic.proteus import PROTEUS

class Test(TestCase):
  """
  pytest-mpl tests comparing similarity of figures.
  
  Notes
  -----
  See the docstring of this module.
  """
  @pytest.mark.mpl_image_compare
  def test_PROTEUS_plot_acq(self):
    PROTEUS().plot_acq()
    plt.xlim(-1e4,1e4)
    plt.ylim(-8e3,1e4)
    return plt.gcf()
<<<<<<< HEAD
  @pytest.mark.mpl_image_compare
  def test_PROTEUS_svp_plot(self):    
    PROTEUS().svp.plot(0,0,0)
    return plt.gcf()
=======
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
